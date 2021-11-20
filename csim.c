/**
 * @author Jonathan Helland
 * @file csim.c
 * @brief A write-back write-allocate cache simulator that allows a
 * user-specified number of set, lines, and block size. Instructions are fed to
 * the cache via a trace files (cf. traces/csim/test.trace for syntax).
 *
 * The core data structure for the cache is essentially a 2D array, where the
 * rows corresponds to sets and the columns correspond to lines (row-major
 * ordering).
 * However, in order to implement a LRU policy, each row is treated as a doubly
 * linked list, meaning that there is some extra overhead required to keep track
 * of who is pointing at whom. Rather than implementing a traditional linked
 * list with pointers, I chose to implement a standard array, and then keep
 * track of indicies into the array (e.g. a parent index and child index for
 * each node). Of course, this is really just pointers under the hood, but my
 * approach makes it more obvious that memory is in continguous chunks for each
 * set of cache lines.
 *
 * The code should be pretty straightforward if you start from `main` and then
 * follow the function pointers, so to speak.
 */

#define _GNU_SOURCE // Need this for `getline`
#include <stdio.h>

// For `getopt`
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>

// misc
#include <assert.h>
#include <ctype.h> // isxdigit
#include <limits.h>
#include <stdbool.h>
#include <string.h>

#define ADDR_SIZE 64

#define MAX_SET_BITS 15
#define MAX_TAG_BITS ADDR_SIZE
#define MAX_BLOCK_BITS 16

#define DIRTY_BIT_MASK 1
#define VALID_BIT_MASK 2

#define NULL_IDX USHRT_MAX // for doubly linked list

/**
 * @brief This struct is just a simple container for user-passed args that
 * configure the cache simulation. An easy way to pass around cache metadata
 * needed by other functions.
 */
typedef struct {
    size_t nTagBits, nSetBits, nBlockBits, nLines;
    size_t cacheSize;
    char printCache;
    char verbose;
    char *traceFilePath;
} cache_cfg_t;

/**
 * @brief Cache line.
 *
 * @param  tag        Stores the tag to match against an instruction address.
 * @param  parentIdx  Doubly linked list overhead. Index of the parent cache
 *                    line within its set (note: not the same as the set index).
 * @param  childIdx   Index of the child cache line within its set.
 * @param  vs         The valid and dirty bits (the upper 2 bits of the char are
 *                    not used).
 */
typedef struct {
    size_t tag;

    // Doubly linked list overhead.
    u_int16_t parentIdx;
    u_int16_t childIdx;

    // valid and dirty bits
    char vd;
} cache_line_t;

/**
 * @brief Cache set containing an array of cache lines. The array is wrapped by
 * a doubly linked list implemented through indices that track the list head and
 * tail. Also note that each cache line contains indices to its parent and child
 * nodes in the list.
 *
 * @param  lines           Array of cache lines.
 * @param  headIdx         Index to the head of the doubly linked list. This
 *                         will always be the most recently used cache line.
 * @param  tailIdx         Index to the tail of the doubly linked list. This
 *                         will always be the least recently used cache line,
 *                         allowing for an LRU policy trivially.
 * @param  nextInvalidIdx  Index to the next cold block.
 */
typedef struct {
    cache_line_t *lines;

    // Linked list maintenance stuff.
    u_int16_t headIdx;
    u_int16_t tailIdx;
    u_int16_t nextInvalidIdx;
} cache_set_t;

/**
 * @brief The full cache struct -- essentially a 2D array with extra sugar.
 *
 * @param  sets    Array of cache sets, each entry of which contains its own
 * array of cache lines.
 * @param  nSets   Number of sets in the cache.
 * @param  nLines  Number of lines in each set. By putting this here, we assume
 * that the sets are homogeneous, although it would be easy to implement varying
 * size cache sets instead.
 */
typedef struct {
    cache_set_t *sets;
    u_int16_t nSets, nLines;
} cache_t;

/**
 * @brief Parse user arguments that determine the cache topology into a config
 * struct.
 *
 * @param[out]  cfg   The struct that contains the full cache configuration
 *                    data.
 * @param[in]   argc  The number of user-passed arguments.
 * @param[in]   argv  The user passed arguments in string form.
 */
void parse_args_(cache_cfg_t *cfg, const int argc, char *const argv[]) {
    // Variable declarations.
    const int numericBase = 10;
    int opt;
    const char *usageStr = "Usage: %s [-v verbose] [-s nSetBits] "
                           "[-E nLines] [-b nBlockBits] [-t]\n";
    size_t nSetBits, nBlockBits, nLines;
    nSetBits = nBlockBits = nLines = 0;

    // Parse user arguments (permutation invariant, as requested).
    cfg->verbose = false;
    cfg->nSetBits = cfg->nLines = cfg->nBlockBits = 0;
    cfg->traceFilePath = NULL;
    while ((opt = getopt(argc, argv, "pvs:E:b:t:")) != EOF) {
        switch (opt) {
        case 'v':
            cfg->verbose = true;
            break;
        case 's':
            nSetBits = strtoul(optarg, NULL, numericBase);
            break;
        case 'E':
            nLines = strtoul(optarg, NULL, numericBase);
            break;
        case 'b':
            nBlockBits = strtoul(optarg, NULL, numericBase);
            break;
        case 't':
            cfg->traceFilePath = optarg;
            break;
        default:
            // If the user has misspecified arguments, we should just halt right
            // now.
            fprintf(stderr, usageStr, argv[0]);
            exit(EXIT_FAILURE); // ok to direct exit because we haven't heap
                                // allocated yet
        }
    }

    if (nSetBits > MAX_SET_BITS) {
        fprintf(stderr, "Too many set bits! Expected <= %u, but got %lu\n",
                MAX_SET_BITS, nSetBits);
        exit(EXIT_FAILURE);
    }
    if (nBlockBits > MAX_BLOCK_BITS) {
        fprintf(stderr,
                "Too many block offset bits! Expected <= %u, but got %lu\n",
                MAX_BLOCK_BITS, nBlockBits);
        exit(EXIT_FAILURE);
    }
    if (nLines <= 0) {
        fprintf(stderr, "Not enough cache lines! Expected > 0, but got %lu\n",
                nLines);
        exit(EXIT_FAILURE);
    }
    if (nLines > USHRT_MAX) {
        fprintf(stderr, "Too many cache lines! Expected <= %u, but got %lu\n",
                USHRT_MAX, nLines);
        exit(EXIT_FAILURE);
    }

    cfg->nSetBits = nSetBits;
    cfg->nBlockBits = nBlockBits;
    cfg->nLines = nLines;

    size_t sum = cfg->nBlockBits + cfg->nSetBits;
    if (sum >= ADDR_SIZE) {
        fprintf(
            stderr,
            "b + s = %lu, which is too large relative to the address size %i\n",
            sum, ADDR_SIZE);
        exit(EXIT_FAILURE); // ok to direct exit because we haven't heap
                            // allocated yet
    }
    cfg->nTagBits = ADDR_SIZE - sum; // The tag will be the remaining bits.
    cfg->cacheSize =
        (1 << cfg->nBlockBits) * (1 << cfg->nSetBits) * cfg->nSetBits;
}

/**
 * @brief Check if a string is syntactically valid with respect to cache
 * instructions.
 *
 * This function will print to stderr with details if the syntax is invalid.
 *
 * @param[in]  lineBuf   A string containing syntax to check.
 * @param[in]  lineSize  The length of the string (excluding '\0').
 *
 * @return true if string contains a valid instruction, false otherwise.
 */
bool is_line_syntax_valid(const char *lineBuf, size_t lineSize) {
    char segment = 0;
    char delimCount = 0;
    for (size_t i = 0; i < lineSize; ++i) {
        switch (lineBuf[i]) {
        case ' ':
            if (delimCount) {
                fprintf(stderr, "Delimiter ' ' must come before ','\n");
                return false;
            }
            segment++;
            break;
        case ',':
            delimCount++;
            segment++;
            break;

        default:
            if (!isalpha(lineBuf[i]))
                break;

            if (segment == 0 && lineBuf[i] != 'L' && lineBuf[i] != 'S') {
                fprintf(stderr, "Not 'L' or 'S' '%c'\n", lineBuf[i]);
                return false;
            } else if (segment == 1 && !isxdigit(lineBuf[i])) {
                fprintf(stderr, "Not xdigit '%c'\n", lineBuf[i]);
                return false;
            } else if (segment == 2 && !isdigit(lineBuf[i])) {
                fprintf(stderr, "Not digit '%c'\n", lineBuf[i]);
                return false;
            } else if (segment > 2) {
                fprintf(stderr, "Too many delimiters\n");
                return false;
            }

            break;
        }
    }

    if (segment < 2) {
        fprintf(stderr, "Not enough delimiters\n");
        return false;
    }
    return true;
}

/**
 * @brief Return codes for trace file parsing functions.
 *
 * @param  SUCCESS     Means that parsing the current instruction was a success.
 * @param  BAD_SYNTAX  The current instruction contains invalid syntax.
 */
enum return_codes { SUCCESS, BAD_SYNTAX };

/**
 * @brief Codes for whether a cache instruction is a save or load operation.
 * Needed for the purposes of determining when to enact LRU, write-back, and
 * write allocate policies.
 */
enum command { LOAD, SAVE };

/**
 * @brief Struct that cache instructions will be parsed into. Makes code for
 * performing instructions more readable.
 */
typedef struct {
    enum command cmd;
    size_t addr;
    size_t bytes;
} cache_instr_t;

/**
 * @brief From an open file stream, parse the next line into a cache
 * instruction.
 *
 * @param[in]   traceFile  An open file stream.
 * @param[out]  instr      A `cache_instr_t` instance.
 *
 * @return EOF        (-1) if end of file stream has been reached.
 *         BAD_SYNTAX (1)  if the line could not be parsed into a valid cache
 *                         instruction.
 *         SUCCESS    (0)  if parsing was successful.
 */
int read_next_instr_(FILE *traceFile, cache_instr_t *instr) {
    char *lineBuf = NULL;
    size_t lineBufSize = 0;
    ssize_t lineSize = getline(&lineBuf, &lineBufSize, traceFile);

    // Early out if end of file.
    if (lineSize == EOF) {
        free(lineBuf);
        return EOF;
    }

    // Now check that the instruction syntax is valid in a first pass over the
    // string. This is slightly wasteful because we have to iterate over the
    // string twice now -- we could instead do parsing and syntax checking at
    // the same time. However, instruction lengths are pretty short, so this
    // shouldn't actually matter in practice. This way, we can avoid messier
    // code that results from simultaneous parsing/syntax checking.
    if (lineSize > 0 && !is_line_syntax_valid(lineBuf, lineSize)) {
        free(lineBuf);
        return BAD_SYNTAX;
    }

    // Parse isntruction assuming the format `instrType addr,numBytes`.
    // We can make this assumption because we know that the syntax of the line
    // is valid if we get to this point.
    char *lineToFree =
        lineBuf; // We must save the starting location of lineBuf for later
                 // freeing because strsep will mutate this pointer.
    char *instrType = strsep(&lineBuf, " ");
    char *addr = strsep(&lineBuf, ",");
    char *numBytes = strsep(&lineBuf, ",");

    instr->cmd = (instrType[0] == 'L') ? LOAD : SAVE;
    instr->addr = strtol(addr, NULL, 16);
    instr->bytes = strtol(numBytes, NULL, 10);

    free(lineToFree);
    return SUCCESS;
}

/**
 * @brief Moves a specified cache line to the head of the doubly-linked list of
 * lines maintained by each set. A helper function for enforcing LRU policy.
 *
 * @param[in]   idx  Index of the cache line to move to the front.
 * @param[out]  set  Pointer to a cache set containing the line associated with
 *                   `idx`.
 */
void move_line_to_head_(const u_int16_t idx, cache_set_t *set) {
    // Head node -- do nothing.
    u_int16_t parent = set->lines[idx].parentIdx;
    if (idx == set->headIdx) {
        return;
    }

    // Move tail to head.
    u_int16_t child = set->lines[idx].childIdx;
    if (idx == set->tailIdx) {
        // Update the head and tail pointers.
        u_int16_t head = set->headIdx;
        set->headIdx = idx;
        set->tailIdx = set->lines[idx].parentIdx;

        // Update the head parent/child pointers.
        set->lines[set->headIdx].parentIdx = NULL_IDX;
        set->lines[set->headIdx].childIdx = head;

        // Update the tail's parent parent/child pointers.
        set->lines[parent].childIdx = NULL_IDX;
        if (parent == head)
            set->lines[parent].parentIdx = idx;

        return;
    }

    // Move non-tail to head.
    set->lines[idx].childIdx = set->headIdx;
    set->lines[idx].parentIdx = NULL_IDX;
    set->headIdx = idx;
    // Patch up the hole.
    set->lines[parent].childIdx = child;
    set->lines[child].parentIdx = parent;
}

/**
 * @brief Heap-allocate the required amount of memory for the cache.
 *
 * @note Only uses one malloc call. Do something else if you want to avoid
 * malloc.
 *
 * @param[out]  cache  An uninitialized `cache_t` instance whose fields will be
 *                     updated.
 * @param[in]   cfg    A `cache_cfg_t` instance.
 */
void allocate_cache_(cache_t *cache, const cache_cfg_t *cfg) {
    cache->nSets = (1 << cfg->nSetBits);
    cache->nLines = cfg->nLines;

    // The data structure we're looking for here is analagous to a 2d array of
    // cache lines. However, if we just straight-up went for a 2d array of
    // `cache_line_t` structs, then we'd have to include a set id for every
    // single line, which wastes space. Instead, we can have a special
    // `cache_set_t` struct that contains an array of cache lines while only
    // needing a single bitfield.
    // Use calloc to zero-initialize everything because we'd have to do this
    // anyways.
    cache->sets = calloc(sizeof(cache_set_t), cache->nSets);
    if (cache->sets == NULL) {
        fprintf(stderr, "malloc returned NULL\n");
        exit(EXIT_FAILURE); // It's okay to direct exit here because we haven't
                            // heap allocated anything yet.
    }

    // Now we hand out blocks of the malloc memory as rows.
    // Also, initialize the fields of each line.
    for (size_t i = 0; i < cache->nSets; ++i) {
        cache_set_t *seti = &cache->sets[i];
        seti->lines = calloc(sizeof(cache_line_t), cache->nLines);

        // Now we have to go through and set up the doubly linked list
        // "pointers". Note that these are not actually pointers, I'm just doing
        // everything with indices into an array instead. There's no need to
        // initialize any of the struct fields because calloc already did the
        // work for us.
        for (size_t j = 0; j < cache->nLines; ++j) {
            // This code is a bit inelegant, but gets the job done.
            // It looks this stupid because of two edges cases: the current node
            // is the head or the tail. The real reason it looks stupid is
            // because I am not a clever man.
            if (j == 0) {
                seti->lines[j].parentIdx = -1;
                seti->lines[j].childIdx = j + 1;
                seti->headIdx = j;
            } else if (j == cache->nLines - 1) {
                seti->lines[j].parentIdx = j - 1;
                seti->lines[j].childIdx = -1;
                seti->tailIdx = j;
            } else {
                seti->lines[j].parentIdx = j - 1;
                seti->lines[j].childIdx = j + 1;
            }
        }
    }
}

/**
 * @brief Cleanup function for when we've finished the simulation.
 *
 * @param[out]  cache  Pointer to a cache that was allocated via the
 *                     `allocate_cache` function.
 */
void free_cache(cache_t *cache) {
    for (size_t i = 0; i < cache->nSets; ++i)
        free(cache->sets[i].lines);
    free(cache->sets);
}

/**
 * @brief Extract the set index from a 64-bit address.
 *
 * @param[in]  addr  64-bit address of a cache line, assumed to be little
 *                   endian. The structure should be
 *                          [tag | set index | block offset].
 * @param[in]  cfg   Cache config struct that stores the number of tag bits, set
 *                   bits, and block offset bits in the address.
 *
 * @returns The set index portion of the address.
 */
size_t get_set_bits(const size_t addr, const cache_cfg_t *cfg) {
    // We need to guard against the case where there are no set bits i.e. only
    // one set in cache.
    return cfg->nSetBits
               ? (addr << cfg->nTagBits) >> (cfg->nTagBits + cfg->nBlockBits)
               : 0;
}

/**
 * @brief Extract the tag from a 64-bit address.
 *
 * @param[in]  addr  64-bit address of a cache line, assumed to be little
 *                   endian. The structure should be
 *                          [tag | set index | block offset].
 * @param[in]  cfg   Cache config struct that stores the number of tag bits, set
 *                   bits, and block offset bits in the address.
 *
 * @returns The tag portion of the address.
 */
size_t get_tag_bits(const size_t addr, const cache_cfg_t *cfg) {
    // Guard against case when there are no tag bits.
    // GUARD CASE SHOULD NEVER HAPPEN
    return cfg->nTagBits ? addr >> (cfg->nSetBits + cfg->nBlockBits) : 0;
}

/**
 * @brief Print a representation of a cache instruction to stdout.
 *
 * @param[in]  instr  Pointer to a cache_instr_t struct.
 */
void print_instr(const cache_instr_t *instr) {
    printf("%c %lx,%zu ", instr->cmd ? 'S' : 'L', instr->addr, instr->bytes);
}

/**
 * @brief A verbosity function that will print a representation of the cache to
 * stdout.
 *
 * @param[in]  cache  Pointer to allocated cache.
 */
void print_cache(const cache_t *cache) {
    printf("\n");
    for (size_t i = 0; i < cache->nSets; ++i) {
        printf("set %lu | ", i);
        for (size_t j = 0; j < cache->nLines; ++j) {
            if (cache->sets[i].headIdx == j)
                printf("H ");
            else if (cache->sets[i].tailIdx == j)
                printf("T ");
            printf("vd(0x%1x) t(0x%8lx) p(%5u) c(%5u) | ",
                   cache->sets[i].lines[j].vd, cache->sets[i].lines[j].tag,
                   cache->sets[i].lines[j].parentIdx,
                   cache->sets[i].lines[j].childIdx);
        }
        printf("\n");
    }
}

/**
 * @brief The main workhorse of the cache simulator: simulates the effect of an
 * instruction on the cache, which is communicated via updates to the `stats`
 * struct.
 *
 * @param[in]   instr  Pointer to a cache_instr_t struct.
 * @param[out]  cache  Pointer to an allocated cache_t struct.
 * @param[out]  stats  Pointer to a csim_stats_t struct that will track how the
 *                     cache state mutates.
 * @param[in]   cfg    Pointer to a cache_cfg_t struct containing metadata about
 *                     the cache.
 */
void perform_cache_instr(const cache_instr_t *instr, cache_t *cache,
                         csim_stats_t *stats, const cache_cfg_t *cfg) {
    // Chunk up address into relevant bit fields that we can do stuff with.
    size_t tag = get_tag_bits(instr->addr, cfg);
    size_t set = get_set_bits(instr->addr, cfg);

    const size_t nCacheLineBytes = (1 << cfg->nBlockBits); // 2 ^ nBlockBits

    // Determine if we have a cache hit.
    cache_set_t *cacheSet = &cache->sets[set];
    bool hit = false;
    size_t lineIdx;
    for (lineIdx = 0; lineIdx < cache->nLines; ++lineIdx) {
        // Hit if block is valid and tags match.
        hit = ((tag == cacheSet->lines[lineIdx].tag) &&
               (VALID_BIT_MASK & cacheSet->lines[lineIdx].vd));

        // Early out to avoid unnecessary lineIdx increment.
        if (hit)
            break;
    }

    // We need to update the usage ordering in the case of a hit.
    if (hit) {
        if (cfg->verbose)
            printf("hit");
        stats->hits++;

        if (instr->cmd == SAVE) {
            if (cfg->verbose)
                printf(" write-back-flag(%lu)", nCacheLineBytes);

            // Write-back policy: flag as dirty but don't flush yet.
            cache_line_t *line = &cacheSet->lines[lineIdx];

            // Check the dirty bit.
            if (!(DIRTY_BIT_MASK & line->vd))
                stats->dirty_bytes += nCacheLineBytes;
            // Set dirty bit to 1 while preserving rest of bitfield.
            line->vd |= DIRTY_BIT_MASK;
        }

        move_line_to_head_(lineIdx, cacheSet);

    } else {
        if (cfg->verbose)
            printf("miss");
        stats->misses++;

        cache_line_t *line;

        // EVICTION MISS
        if (VALID_BIT_MASK & cacheSet->lines[cacheSet->tailIdx].vd) {
            if (cfg->verbose)
                printf(" eviction");
            stats->evictions++;

            line = &cacheSet->lines[cacheSet->tailIdx];

            // Write-back policy: update lower memory hierarchy.
            // Check the dirty bit.
            if (line->vd & 1) {
                if (cfg->verbose)
                    printf(" write-back-flush(-%lu)", nCacheLineBytes);

                line->vd &= 0xe; // set dirty bit to zero while preserving rest
                                 // of bitfield
                stats->dirty_evictions += nCacheLineBytes;
                stats->dirty_bytes -= nCacheLineBytes;
            }

            // LRU is now the most recent.
            move_line_to_head_(cacheSet->tailIdx, cacheSet);

            // We're done with this index, so set it to a non-index value.
            cacheSet->nextInvalidIdx = USHRT_MAX;

            // COLD MISS
        } else {
            if (cfg->verbose)
                printf("(cold)");

            // Need to mark the block as now valid.
            line = &cacheSet->lines[cacheSet->nextInvalidIdx];
            // set valid bit to one while preserving rest of bitfield
            line->vd |= VALID_BIT_MASK;

            // Update most recently used to this newly loaded line.
            move_line_to_head_(cacheSet->nextInvalidIdx, cacheSet);

            // Move to next invalid block for next cold miss.
            cacheSet->nextInvalidIdx++;
        }

        line->tag = tag;

        // Write-allocate policy: immediately update cache.
        if (instr->cmd == SAVE) {
            stats->dirty_bytes += nCacheLineBytes;
            // set dirty bit to one while preserving rest of bitfield
            line->vd |= DIRTY_BIT_MASK;
        }
    }
}

/**
 * @brief Store a summary of the cache simulation statistics.
 *
 * NOTE: Pasted here from another file for clarity and self-containment.
 *
 * @param[in] stats The simulation statistics to be stored
 */
void printSummary(const csim_stats_t *stats) {
    printf("hits:%ld misses:%ld evictions:%ld dirty_bytes_in_cache:%ld "
           "dirty_bytes_evicted:%ld\n",
           stats->hits, stats->misses, stats->evictions, stats->dirty_bytes,
           stats->dirty_evictions);

    FILE *output_fp = fopen(".csim_results", "w");
    if (output_fp == NULL) {
        fprintf(stderr, "Error: failed to open results file: %s\n",
                strerror(errno));
        return;
    }

    fprintf(output_fp, "%ld %ld %ld %ld %ld\n", stats->hits, stats->misses,
            stats->evictions, stats->dirty_bytes, stats->dirty_evictions);
    fclose(output_fp);
}

/**
 * @brief Main driver for the cache simulation.
 *
 * @param[in]  argc  Number of passed arguments.
 * @param[in]  argv  String representation of the passed arguments.
 */
int main(int argc, char *argv[]) {
    // Parse user arguments to get the cache configuration.
    // We should do this first to fail as quickly as possible if the user
    // specifies the target cache incorrectly.
    cache_cfg_t cfg;
    parse_args_(&cfg, argc, argv);
    if (cfg.verbose)
        printf("[VERBOSE MODE]\nRunning cache simulator with the following "
               "config:\n"
               "[-s] nSetBits %lu\n"
               "[-E] nLines %lu\n"
               "[-b] nBlockBits %lu\n"
               "     nTagBits %lu\n"
               "     cacheSize %lu\n"
               "[-t] traceFilePath %s\n"
               "--------------------------------------------------\n",
               cfg.nSetBits, cfg.nLines, cfg.nBlockBits, cfg.nTagBits,
               cfg.cacheSize, cfg.traceFilePath);

    // Now we need to open up the trace file to parse for instructions.
    // First, make sure that the damn thing actually exists.
    // We want to do this next due to the fail fast principle: if the user
    // has specified a cache trace file that does not exist, we don't want
    // to even bother allocating any resources for the cache because we're
    // just going to free those resources immediately anyways. As an
    // additional bonus, this means less case-specific cleanup is required.
    if (access(cfg.traceFilePath, F_OK)) {
        fprintf(stderr, "File `%s` does not exist!\n", cfg.traceFilePath);
        exit(EXIT_FAILURE); // ok to direct exit here because we haven't heap
                            // allocated yet
    }
    FILE *traceFile = fopen(cfg.traceFilePath, "r");

    // Allocate resources for the target cache we'll simulate.
    if (cfg.verbose)
        printf("Allocating cache\n");
    cache_t cache;
    allocate_cache_(&cache, &cfg);

    csim_stats_t stats;
    stats.hits = stats.misses = stats.evictions = stats.dirty_bytes =
        stats.dirty_evictions = 0;

    cache_instr_t instr;
    int opt;
    int count = 0;
    while ((opt = read_next_instr_(traceFile, &instr)) != EOF) {
        if (cfg.verbose) {
            printf("(%i) ", ++count);
            print_instr(&instr);
            print_cache(&cache);
        }

        // If the trace instruction is not valid, then we'll just go to the next
        // instruction.
        if (opt == BAD_SYNTAX) {
            if (cfg.verbose)
                fprintf(stderr, "Instruction gave BAD_SYNTAX error\n\n");
            continue;
        }

        perform_cache_instr(&instr, &cache, &stats, &cfg);

        if (cfg.verbose)
            printf("\n");
    }

    // This should always be called regardless of cfg->verbose state.
    printSummary(&stats);

    // Final cleanup.
    if (cfg.verbose)
        printf("Cleaning up...\n");

    free_cache(&cache);
    fclose(traceFile);
    exit(EXIT_SUCCESS);
}
