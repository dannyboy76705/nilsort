# nilsort

Sort files by [nilsimsa](https://en.wikipedia.org/wiki/Nilsimsa_Hash) similarity, so that similar files end up adjacent in the output list. Point it at a directory, get back an ordered file list on stdout — pipe that into an archiver (tar, DwarFS, etc.) and similar files land next to each other, which can meaningfully improve compression on redundant data (near-duplicate documents, versioned assets, video frames, and so on).

```
$ nilsort some_dir | tar -T- --no-recursion -czf out.tar.gz
```

## Why this exists

[`binsort`](https://github.com/tmueller/binsort) already does this using [simhash](https://en.wikipedia.org/wiki/SimHash), and [DwarFS](https://github.com/mhx/dwarfs) does something similar internally using nilsimsa as part of building its filesystem image. Neither exposes a standalone "just give me the nilsimsa-based sort order" tool — `binsort` uses a different hash entirely, and DwarFS's ordering step only runs as part of building a full filesystem image, not as an independent utility. `nilsort` fills that gap: a small, dependency-free C program that does nilsimsa hashing and similarity ordering on their own, with nothing else attached.

## How it works

1. **Hash.** Every file gets a 256-bit nilsimsa digest (multithreaded — one thread pool hashes all files in parallel, each digest independent of the others).
2. **Order.** Files are greedily chained: start anywhere, repeatedly append whichever unvisited file is most similar to the last one added.
3. **Refine.** [2-opt](https://en.wikipedia.org/wiki/2-opt) local-search passes then look for pairs of segments to reverse if doing so increases total adjacent similarity — the same idea as `binsort`'s `-o` optimization level. Runs until it stops finding improvements (or hits a sweep cap).

By default the **entire input is treated as one batch** — every file gets weighed against every other file for maximum accuracy. This is deliberately different from `binsort`, which has to manage an O(n²) *memory* footprint (a full pairwise distance matrix) and therefore documents real limits around very large file counts. `nilsort` only ever holds digests (32 bytes/file) and index arrays in memory — the O(n²) cost here is pure comparison *time*, not memory — so file count that would be memory-prohibitive for `binsort` is generally fine here. Batching (`-b`) is available if you want it, but it's opt-in, not the default.

## Build

```bash
gcc -O2 -Wall -pthread -march=native -o nilsort nilsort.c
```

`-march=native` builds for the exact CPU you're compiling on and picks up hardware POPCNT automatically (digest comparison is a Hamming-distance calculation — POPCNT measured ~8x faster than a lookup-table approach). Build separately on each machine you deploy to; don't copy the binary across machines with different CPUs.

If you're building once and running elsewhere on possibly-older hardware, use this instead:

```bash
gcc -O2 -Wall -pthread -mpopcnt -o nilsort nilsort.c
```

`-mpopcnt` enables just the POPCNT instruction (present on essentially all x86-64 CPUs since ~2008 Intel / ~2012 AMD) without tying the binary to your exact microarchitecture the way `-march=native` does.

## Usage

```
nilsort [options] <dir>
find /some/dir -type f | nilsort [options] -
```

| Option | Meaning |
|---|---|
| `-b N` | Batch size (files per similarity-sort batch). Default: unbounded — the whole input is one batch. |
| `-o N` | Max 2-opt refinement sweeps. Default: 50 (converges early once a sweep finds no improvement). `-o 0` disables refinement. |
| `-t N` | Number of hashing threads. Default: number of CPUs. |
| `-q` | Quiet — suppress progress bars and status messages. |
| `-v` | Print a memory/time estimate before running. |
| `-h` | Show help. |

Threading only affects the hashing phase's wall-clock time — each file's digest is written to a fixed array slot by index, independent of which thread computed it or in what order, so output is bit-for-bit identical regardless of `-t`. (Notably, `binsort` documents the opposite behavior for its own optimization stage.)

## Performance notes

Real numbers from testing, for a sense of scale:

- **~33,000 files:** hashing + greedy construction took ~14s with a table-based digest comparison, dropping to ~3s after switching to hardware POPCNT.
- **Memory:** ~11MB peak RSS at 33,000 files (digests + order arrays only) — for comparison, DwarFS's own sort phase used ~32MB at a similar scale, doing considerably more per-file bookkeeping (content-defined chunking, block-level dedup hashing, filesystem metadata) as part of the same pass.
- Both memory and time scale with file *count*, not file *size* — the comparison step only ever touches 32-byte digests, so this behaves the same whether you're sorting text files or multi-megabyte video frames.

Since comparison is O(n²), very large single batches (hundreds of thousands of files) will take a while, especially across multiple refinement sweeps — memory won't be the limiting factor, time will. `-v` gives you an upfront estimate; `-o` and `-b` are there if you want to trade some accuracy back for speed.

## License

MIT. See the header comment in `nilsort.c`.

The nilsimsa algorithm itself (the TRAN/POPC-derived tables, `tran3` mixing function, and digest construction) is a faithful port of the long-standing public reference implementation originally written by cmeclax, based on Damiani et al. 2004, "An Open Digest-based Technique for Spam Detection." Everything else — batching, threading, greedy + 2-opt ordering, progress reporting, and the POPCNT-based comparison — is original to this project.
