# CPU Cache simulator in C for CMU 15213
A write-back, write-allocate cache simulator that allows a user-specified number of set, lines, and block size. Includes a parser for trace files that specify CPU instructions.

I used this simulator to help write cache-friendly matrix transpose algorithms for a later assignment.
In particular, I used it to test my implementations of some blocking and cache-oblivious algorithms.

The cache implementation here formed the basis of my design of a thread-safe cache for proxy-lab later
in the course.

# Documentation
The code should be pretty straightforward if you start from [`csim.c::main(678)`](./csim.c) and then
follow the function pointers, so to speak.
Yep, my documentation for this is basically "read the code lol".

### **API**
```
Usage: [-v verbose] [-s nSetBits] [-E nLines] [-b nBlockBits] [-t traceFile]
```

### **Trace files**

Instructions are fed to the simulated CPU from trace files that use a simple language, where each line in the file has the form `[instruction] [address],[number of bytes]`.
- `[instruction]`:
    - `L` means read.
    - `S` means write.

Example trace files are provided in [`traces/`](./traces/).

For instance, [traces/test.trace](./traces/test.trace):
```
L 10,4
S 18,4
```

### **Data structures**
The core data structure for the cache is essentially a 2D array, where the
rows corresponds to sets and the columns correspond to lines (row-major
ordering).
However, in order to implement a LRU policy, each row is treated as a doubly
linked list, meaning that there is some extra overhead required to keep track
of who is pointing at whom. Rather than implementing a traditional linked
list with pointers, I chose to implement a standard array, and then keep
track of indicies into the array (e.g. a parent index and child index for
each node). Of course, this is really just pointers under the hood, but my
approach makes it more obvious that memory is in continguous chunks for each
set of cache lines.
