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
| `-m` | EBML/Matroska-aware hashing. For files that are actually `.mkv`/`.webm` (detected by magic bytes, not extension), hash only `Cluster` (payload) elements, skipping container framing — the EBML header, `SeekHead`, `Info`, `Tracks`, `Cues`, etc. Non-Matroska files in the same run are hashed normally either way. Implies `-M`. |
| `-M` | Use MinHash (literal shared-substring similarity) instead of nilsimsa (statistical-texture similarity) for the whole run, on any file type. On automatically with `-m`; pass alone to try it without container-stripping. |
| `-f PATH` | Write the sorted file list to `PATH` instead of stdout, and print a confirmation (`wrote N file path(s) to PATH`) once done. Handy so a forgotten `>` redirect doesn't silently lose the output. |
| `-h` | Show help. |

Threading only affects the hashing phase's wall-clock time — each file's digest is written to a fixed array slot by index, independent of which thread computed it or in what order, so output is bit-for-bit identical regardless of `-t`. (Notably, `binsort` documents the opposite behavior for its own optimization stage.)

## Performance notes

Real numbers from testing, for a sense of scale:

- **~33,000 files:** hashing + greedy construction took ~14s with a table-based digest comparison, dropping to ~3s after switching to hardware POPCNT.
- **Memory:** ~11MB peak RSS at 33,000 files (digests + order arrays only) — for comparison, DwarFS's own sort phase used ~32MB at a similar scale, doing considerably more per-file bookkeeping (content-defined chunking, block-level dedup hashing, filesystem metadata) as part of the same pass.
- Both memory and time scale with file *count*, not file *size* — the comparison step only ever touches 32-byte digests, so this behaves the same whether you're sorting text files or multi-megabyte video frames.

Since comparison is O(n²), very large single batches (hundreds of thousands of files) will take a while, especially across multiple refinement sweeps — memory won't be the limiting factor, time will. `-v` gives you an upfront estimate; `-o` and `-b` are there if you want to trade some accuracy back for speed.

## Duplicate detection

Alongside the similarity sort, `nilsort` also tracks exact duplicates — files whose content is byte-identical over the same bytes fed to nilsimsa (whole file normally, or just the Cluster payload when `-m` is active, so container metadata differences don't prevent a true duplicate frame from being flagged). This uses a 64-bit FNV-1a fingerprint plus a length check computed during the same read pass, at negligible extra cost. At the end of a run it prints a one-line summary to stderr (suppressed by `-q`, like other status output):

```
nilsort: 4 duplicate files found in 2 sets of identical content
```

This is a count only — it doesn't list which files or print anything to stdout, so it won't interfere with piping the sort order into an archiver.

## Similarity hasher: nilsimsa vs MinHash (-M / -m)

By default `nilsort` uses nilsimsa — a statistical-texture similarity hash (does this file have a similar overall byte-frequency histogram to that one). For some data types, that's a different question than what LZ-family compressors (zstd, xz) actually exploit: whether two files share literal, reusable byte substrings. `-M` switches the whole run to MinHash instead — sliding an 8-byte shingle window across each file and tracking, across 32 independent hash functions, the minimum value seen; the similarity score is how many of those 32 minimums match between two files (the standard MinHash/Jaccard estimator). `-m` (EBML-aware hashing) implies `-M` automatically, since MinHash measured better for the video-container use case `-m` targets; use `-M` alone to try it on non-Matroska data.

Unlike nilsimsa, MinHash here has no published reference implementation to validate against — the technique is standard, but this specific choice of shingle size and function count is this project's own. It's validated empirically instead: on synthetic files with a known, controlled fraction of content replaced, `minhash_compare` scores dropped monotonically from 32/32 (identical) through 29→23→18→9→4→2→0 as replacement went from 0% to 100%, and two fully independent random files scored 0/32 across repeated trials with no false-positive matches. Treat it as a validated heuristic, not a proven-correct standard the way nilsimsa is.

Duplicate detection (see above) is unaffected by which similarity hasher is active — it always uses the separate FNV-1a exact-match fingerprint.

## Container framing and -m

If you're sorting per-frame or per-GOP `.mkv` files from the same encoding pipeline, be aware that Matroska container overhead (EBML header, `SeekHead`, `Info`, `Tracks`, etc.) can be a large enough fraction of very small files — and near-identical across all of them — to distort similarity scoring. In one measured case, two single-frame FFV1 `.mkv` files with genuinely different frame content (red vs. blue) shared a 32% byte-identical prefix purely from container framing, and whole-file nilsimsa scores were compressed tightly enough that a different-content pair scored as *more* similar than a same-content pair. Passing `-m` fixed this: same-content frames scored the maximum possible 128, and different-content frames separated out cleanly. If your pipeline produces per-frame/per-GOP MKVs, `-m` is worth using by default; for full multi-GOP episode files the header is a much smaller fraction of total bytes and this matters far less.

## Benchmark results

Tested against three real video-frame corpora (three MPEG-2 episodes split into single-GOP files, one FFV1 episode split into single-frame files, and the same FFV1 episode's raw YUV frames), piped through `tar | zstd --long` and `xz`, comparing `nilsort`'s output against plain directory order as a baseline.

**Headline result: directory order won every single test.** Across every codec (MPEG-2, FFV1, raw YUV), every compressor and setting tried (zstd -3 through -13, xz default), every hasher (nilsimsa, MinHash), and both container modes, `nilsort`'s similarity-based ordering never beat simple directory order. The margins ranged from negligible (as low as +0.03%) to substantial (+10.15% for MinHash on raw YUV).

**Why, and what this means for whether to use `-m`/`-M`:** in every one of these corpora, files were named/traversed in true chronological (creation) order. Natural video is a physically continuous process — each frame differs from the last by a small increment — so chronological order is already very close to compression-optimal, for free. `nilsort`'s job in this benchmark was never "impose structure on unordered data," it was "beat ground truth using only an imperfect approximation of it inferred from content." **If your pipeline preserves true frame/creation order in filenames or directory traversal (as `mkvgopsplit` output does), sorting is unlikely to help and can measurably hurt.** `nilsort`'s more plausible use case — genuinely scrambled or unordered file sets with no reliable naming convention — was not tested here; a fair test of that would compare against a random shuffle, not chronological order.

**nilsimsa vs MinHash is not a fixed winner — it flips depending on whether the data is compressed or raw.** On entropy-coded formats (FFV1, MPEG-2), MinHash's literal-substring matching consistently beat nilsimsa's statistical-texture histogram, closing roughly 23% of nilsimsa's gap to directory order on the MPEG-2 set. But on **raw, uncompressed YUV, this reversed dramatically**: nilsimsa beat MinHash by a wide margin (compressed output +2.59% vs directory order, versus MinHash's +10.15%). The explanation traces to a measurable proxy: nilsimsa recovered 96.48% of true frame adjacency on raw YUV (vs. MinHash's 42.53%), because natural pixel noise breaks MinHash's exact-8-byte-shingle matching even between visually near-identical frames, while nilsimsa's histogram approach is far more tolerant of that noise. On compressed bitstreams the situation flips: literal substrings recur meaningfully (shared header structure, similar quantization patterns) in a way MinHash can exploit better than nilsimsa's histogram can. **Rule of thumb: MinHash (`-M`/`-m`) for compressed/entropy-coded formats, plain nilsimsa (no flags) for raw/uncompressed data.**

**`-m` (container-stripping) matters a lot for single-frame files, much less for multi-frame GOPs.** On FFV1 (one frame per file), `-m` dramatically changed the ordering (only 1.3% adjacency overlap between stripped and unstripped runs) and roughly 7x improved true-adjacency recovery. On MPEG-2 (many frames per GOP file), the effect was much smaller (35.7% adjacency overlap, duplicate counts identical with or without `-m`) since container overhead is a far smaller fraction of a whole GOP than of one frame.

**One genuine win: splitting didn't cost you anything on FFV1.** Despite the per-frame MKV container overhead bloating the split set to 7.1GiB (2.45x the 2.9GiB unsplit original), both directory-order and `-m`-sorted split-and-recompressed archives came in *smaller* than the original unsplit file (~0.5% smaller). So for FFV1 at least, the split-for-frame-access pipeline isn't a size tradeoff at all — you get per-frame access and don't lose anything on disk.

**Duplicate detection cross-validated cleanly.** The same 282-file/47-set duplicate count showed up independently on the FFV1 split, the raw YUV frames, and DwarFS's own duplicate detection — three different representations/tools agreeing is strong evidence these are real repeated frames in the source content (static shots, held frames), not an artifact of any one tool.

## Notes for future work

A few things worth knowing if you're picking this up again, extending it, or building a spinoff:

- **The real missing experiment from this benchmark run is nilsort vs. a random shuffle**, not vs. directory order. Every "default" baseline tested here was actually true chronological order, which turned out to be a very strong baseline for naturally continuous video. The scenario `nilsort` is more plausibly useful for — recovering good locality from a genuinely scrambled or unordered file set — was never directly tested. A quick way to test it: `shuf` the file list before piping to the archiver, and compare against that instead.
- **MinHash's shingle size (8 bytes) and function count (32) were chosen without tuning** — they work, and are validated to behave correctly on synthetic data (see the hasher section above), but nobody's swept those parameters against real compression outcomes. Given MinHash's clear edge over nilsimsa on compressed formats, a larger function count (more discriminative, but ~linearly more compute — see the throughput benchmarks above) or a different shingle size might do meaningfully better; untested.
- **Greedy nearest-neighbor + 2-opt optimizes pairwise similarity, not path smoothness.** This is likely part of why even a 96%-accurate similarity signal (nilsimsa on raw YUV) still lost to true order: the algorithm can and does take short "jumps" whenever a slightly-higher-scoring non-adjacent candidate exists, even when the true-sequence neighbor was already a very good match. An ordering algorithm that penalized non-local jumps, or that had access to a "prefer true-sequence order when the similarity signal doesn't clearly disagree" heuristic, might close more of the remaining gap to chronological order than pure greedy-similarity does.
- **The EBML parser (`-m`) only reads Cluster elements at the Segment's top level** — it doesn't descend into nested structures like `Tags` with embedded attachments, and treats "unknown size" elements (the EBML streaming convention) as running to the end of the enclosing element rather than fully implementing EBML's streaming semantics. Fine for finalized single-frame/single-GOP files from a known pipeline; would need hardening for arbitrary Matroska files in the wild.
- **Every performance and correctness claim in this file was empirically measured against real or controlled-synthetic data during development**, not just reasoned about — see the git history / prior conversation for the actual benchmark harnesses (POPCNT comparison speed, MinHash monotonicity on controlled content-overlap, hashing throughput, thread-order determinism, EBML parsing against real ffmpeg-generated files). Worth maintaining that standard for any new claims rather than trusting intuition alone — several of the findings above (MinHash losing badly on raw YUV, default beating everything) were genuinely counter to the initial hypothesis going in.

## License

MIT. See the header comment in `nilsort.c`.

The nilsimsa algorithm itself (the TRAN/POPC-derived tables, `tran3` mixing function, and digest construction) is a faithful port of the long-standing public reference implementation originally written by cmeclax, based on Damiani et al. 2004, "An Open Digest-based Technique for Spam Detection." Everything else — batching, threading, greedy + 2-opt ordering, progress reporting, and the POPCNT-based comparison — is original to this project.
