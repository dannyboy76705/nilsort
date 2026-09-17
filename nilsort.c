/*
 * nilsort - sort files by nilsimsa similarity
 *
 * MIT License
 *
 * Copyright (c) 2026 <Daniel Lee Witzel>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * The nilsimsa hashing algorithm implemented below (TRAN/POPC-derived
 * tables, tran3 mixing function, digest construction) is a faithful
 * port of the long-standing public reference implementation originally
 * written by cmeclax, based on Damiani et al. 2004, "An Open
 * Digest-based Technique for Spam Detection". The batching, threading,
 * greedy + 2-opt ordering, progress reporting, and POPCNT-based
 * comparison are original to this project.
 *
 * Scans a directory (or reads a file list from stdin), computes a
 * nilsimsa locality-sensitive hash digest for each file, then
 * greedily orders the files so that similar files end up adjacent
 * in the output list. This is the same general idea as `binsort`
 * (which uses simhash internally) and the file-ordering step in
 * DwarFS (which uses nilsimsa internally) -- except this tool just
 * prints the ordered file list to stdout, so you can pipe it into
 * an archiver the same way you would with binsort:
 *
 *   $ nilsort <dir> | tar -T- --no-recursion -czf out.tar.gz
 *
 * By default this tool treats its ENTIRE input as one batch, so every
 * file is weighed against every other file -- e.g. point it at one
 * episode's frame directory and all of that episode's frames are
 * compared against each other in a single sort, with nothing split
 * off into a separate, unrelated batch. This maximizes accuracy at
 * the cost of O(n^2) time. Memory is not the constraint here: digests
 * are only 32 bytes each, so even a million files is ~32 MB of digest
 * storage. If you deliberately want to cap a run (e.g. you're feeding
 * it multiple episodes' worth of files in one invocation and are fine
 * losing cross-batch accuracy for speed), pass -b to split into
 * batches; otherwise leave it unset for one full, maximally-accurate
 * sort per invocation. If you have multiple episodes, run nilsort once
 * per episode directory rather than once across all episodes combined
 * -- that keeps each episode's frames weighed as a single, self-
 * contained sort, which is usually what you want anyway.
 *
 * On top of the initial greedy nearest-neighbor construction, this
 * tool also runs 2-opt local-search refinement passes by default
 * (like binsort's -o optimization level, but here defaulting to
 * "run until it stops improving" since accuracy matters more than
 * speed). Each pass is another O(n^2) scan; use -o to cap or disable
 * it if a run is taking too long.
 *
 * Build:
 *   gcc -O2 -Wall -pthread -mpopcnt -o nilsort nilsort.c
 * (-mpopcnt tells gcc to emit the hardware POPCNT instruction for
 * digest comparison, ~8x faster than a lookup-table popcount in
 * testing. It's present on effectively all x86-64 CPUs since ~2008
 * (Intel) / ~2012 (AMD), so this is safe on any modern machine --
 * but the resulting binary will only run on CPUs that actually have
 * the instruction. If you need to build once and run on unknown or
 * older hardware, drop -mpopcnt: gcc then emits a portable software
 * popcount instead, which works everywhere but loses the speedup.)
 *
 * Usage:
 *   nilsort [options] <dir>
 *   find /some/dir -type f | nilsort [options] -
 *
 * Options:
 *   -b N   batch size (files per similarity-sort batch). Default:
 *          unbounded -- the whole input is treated as one batch.
 *   -o N   max 2-opt refinement sweeps after the initial greedy sort.
 *          Default: 50 (effectively "until convergence" for typical
 *          episode-sized batches). -o 0 disables refinement.
 *   -t N   number of hashing threads, default = number of CPUs
 *   -q     quiet (suppress progress messages on stderr)
 *   -v     print a memory/behavior estimate before running
 *   -h     show this help
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdatomic.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>

/* ---------------------------------------------------------------- */
/* Progress bar (stderr only, throttled to ~10 updates/sec)         */
/* ---------------------------------------------------------------- */

static void progress_bar(const char *label, size_t done, size_t total) {
    static struct timespec last_print;
    static int have_last = 0;
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    if (have_last && done != total) {
        double dt = (double)(now.tv_sec - last_print.tv_sec) +
                    (double)(now.tv_nsec - last_print.tv_nsec) / 1e9;
        if (dt < 0.1) return; /* throttle so huge n doesn't flood stderr */
    }
    last_print = now;
    have_last = 1;

    const int width = 30;
    double frac = total ? (double)done / (double)total : 1.0;
    if (frac > 1.0) frac = 1.0;
    int filled = (int)(frac * width);

    fprintf(stderr, "\r%-10s [", label);
    for (int i = 0; i < width; i++) fputc(i < filled ? '#' : '-', stderr);
    fprintf(stderr, "] %5.1f%% (%zu/%zu)", frac * 100.0, done, total);
    if (done >= total) fputc('\n', stderr);
    fflush(stderr);
}

/* ---------------------------------------------------------------- */
/* Nilsimsa core                                                    */
/* ---------------------------------------------------------------- */

static const unsigned char TRAN[256] = {
    0x02,0xD6,0x9E,0x6F,0xF9,0x1D,0x04,0xAB,0xD0,0x22,0x16,0x1F,0xD8,0x73,0xA1,0xAC,
    0x3B,0x70,0x62,0x96,0x1E,0x6E,0x8F,0x39,0x9D,0x05,0x14,0x4A,0xA6,0xBE,0xAE,0x0E,
    0xCF,0xB9,0x9C,0x9A,0xC7,0x68,0x13,0xE1,0x2D,0xA4,0xEB,0x51,0x8D,0x64,0x6B,0x50,
    0x23,0x80,0x03,0x41,0xEC,0xBB,0x71,0xCC,0x7A,0x86,0x7F,0x98,0xF2,0x36,0x5E,0xEE,
    0x8E,0xCE,0x4F,0xB8,0x32,0xB6,0x5F,0x59,0xDC,0x1B,0x31,0x4C,0x7B,0xF0,0x63,0x01,
    0x6C,0xBA,0x07,0xE8,0x12,0x77,0x49,0x3C,0xDA,0x46,0xFE,0x2F,0x79,0x1C,0x9B,0x30,
    0xE3,0x00,0x06,0x7E,0x2E,0x0F,0x38,0x33,0x21,0xAD,0xA5,0x54,0xCA,0xA7,0x29,0xFC,
    0x5A,0x47,0x69,0x7D,0xC5,0x95,0xB5,0xF4,0x0B,0x90,0xA3,0x81,0x6D,0x25,0x55,0x35,
    0xF5,0x75,0x74,0x0A,0x26,0xBF,0x19,0x5C,0x1A,0xC6,0xFF,0x99,0x5D,0x84,0xAA,0x66,
    0x3E,0xAF,0x78,0xB3,0x20,0x43,0xC1,0xED,0x24,0xEA,0xE6,0x3F,0x18,0xF3,0xA0,0x42,
    0x57,0x08,0x53,0x60,0xC3,0xC0,0x83,0x40,0x82,0xD7,0x09,0xBD,0x44,0x2A,0x67,0xA8,
    0x93,0xE0,0xC2,0x56,0x9F,0xD9,0xDD,0x85,0x15,0xB4,0x8A,0x27,0x28,0x92,0x76,0xDE,
    0xEF,0xF8,0xB2,0xB7,0xC9,0x3D,0x45,0x94,0x4B,0x11,0x0D,0x65,0xD5,0x34,0x8B,0x91,
    0x0C,0xFA,0x87,0xE9,0x7C,0x5B,0xB1,0x4D,0xE5,0xD4,0xCB,0x10,0xA2,0x17,0x89,0xBC,
    0xDB,0xB0,0xE2,0x97,0x88,0x52,0xF7,0x48,0xD3,0x61,0x2C,0x3A,0x2B,0xD1,0x8C,0xFB,
    0xF1,0xCD,0xE4,0x6A,0xE7,0xA9,0xFD,0xC4,0x37,0xC8,0xD2,0xF6,0xDF,0x58,0x72,0x4E
};

/* Similarity score via hardware POPCNT: compare the two 32-byte digests
 * as four 64-bit words, XOR, and count set bits with __builtin_popcountll
 * (compiles to a single POPCNT instruction with -mpopcnt / -march=native
 * on x86, or a portable software fallback otherwise -- correct either
 * way, just faster when the hardware instruction is available). This
 * replaces a per-byte lookup-table approach that measured ~8x slower
 * for this same computation on typical x86 hardware; note this affects
 * only *how fast* the comparison runs, not the digests or scores
 * themselves, so output is bit-for-bit identical either way.
 * 128 = identical digests, going down as they diverge. */
static inline int nilsimsa_compare(const unsigned char *d1, const unsigned char *d2) {
    uint64_t a[4], b[4];
    memcpy(a, d1, 32);
    memcpy(b, d2, 32);
    int bits = __builtin_popcountll(a[0] ^ b[0]) +
               __builtin_popcountll(a[1] ^ b[1]) +
               __builtin_popcountll(a[2] ^ b[2]) +
               __builtin_popcountll(a[3] ^ b[3]);
    return 128 - bits;
}

typedef struct {
    long count;         /* bytes seen so far */
    uint32_t acc[256];  /* trigram accumulators */
    int lastch[4];      /* last 4 bytes seen, -1 = not yet valid */
} nilsimsa_t;

static void nilsimsa_init(nilsimsa_t *n) {
    n->count = 0;
    memset(n->acc, 0, sizeof(n->acc));
    n->lastch[0] = n->lastch[1] = n->lastch[2] = n->lastch[3] = -1;
}

static inline unsigned char tran3(unsigned char a, unsigned char b, unsigned char c, int nn) {
    return (unsigned char)(((TRAN[(a + nn) & 255] ^ (TRAN[b] * (nn + nn + 1))) +
                             TRAN[c ^ TRAN[nn]]) & 255);
}

static void nilsimsa_update(nilsimsa_t *n, const unsigned char *data, size_t len) {
    int *l = n->lastch;
    for (size_t k = 0; k < len; k++) {
        unsigned char ch = data[k];
        n->count++;
        if (l[1] > -1) {
            n->acc[tran3(ch, (unsigned char)l[0], (unsigned char)l[1], 0)]++;
        }
        if (l[2] > -1) {
            n->acc[tran3(ch, (unsigned char)l[0], (unsigned char)l[2], 1)]++;
            n->acc[tran3(ch, (unsigned char)l[1], (unsigned char)l[2], 2)]++;
        }
        if (l[3] > -1) {
            n->acc[tran3(ch, (unsigned char)l[0], (unsigned char)l[3], 3)]++;
            n->acc[tran3(ch, (unsigned char)l[1], (unsigned char)l[3], 4)]++;
            n->acc[tran3(ch, (unsigned char)l[2], (unsigned char)l[3], 5)]++;
            n->acc[tran3((unsigned char)l[3], (unsigned char)l[0], ch, 6)]++;
            n->acc[tran3((unsigned char)l[3], (unsigned char)l[2], ch, 7)]++;
        }
        l[3] = l[2]; l[2] = l[1]; l[1] = l[0]; l[0] = ch;
    }
}

/* Produce the 32-byte (256-bit) digest from accumulated state. */
static void nilsimsa_digest(const nilsimsa_t *n, unsigned char out[32]) {
    long total;
    if (n->count == 3)      total = 1;
    else if (n->count == 4) total = 4;
    else if (n->count > 4)  total = 8 * n->count - 28;
    else                    total = 0;

    double threshold = total / 256.0;
    unsigned char code[32];
    memset(code, 0, sizeof(code));
    for (int i = 0; i < 256; i++) {
        if ((double)n->acc[i] > threshold) {
            code[i >> 3] |= (unsigned char)(1u << (i & 7));
        }
    }
    /* reference implementation reverses byte order in the result */
    for (int i = 0; i < 32; i++) out[i] = code[31 - i];
}

/* ---------------------------------------------------------------- */
/* File hashing                                                     */
/* ---------------------------------------------------------------- */

static int hash_file(const char *path, unsigned char digest_out[32]) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    nilsimsa_t n;
    nilsimsa_init(&n);
    unsigned char buf[1 << 16];
    size_t r;
    while ((r = fread(buf, 1, sizeof(buf), f)) > 0) {
        nilsimsa_update(&n, buf, r);
    }
    fclose(f);
    nilsimsa_digest(&n, digest_out);
    return 0;
}

typedef struct {
    char **paths;
    unsigned char (*digests)[32];
    size_t start, end;
    int quiet;
    atomic_size_t *done_counter;
} hash_job_t;

static void *hash_worker(void *arg) {
    hash_job_t *j = (hash_job_t *)arg;
    for (size_t i = j->start; i < j->end; i++) {
        if (hash_file(j->paths[i], j->digests[i]) != 0) {
            if (!j->quiet) {
                fprintf(stderr, "\nnilsort: warning: cannot read %s: %s\n",
                        j->paths[i], strerror(errno));
            }
            memset(j->digests[i], 0, 32);
        }
        atomic_fetch_add(j->done_counter, 1);
    }
    return NULL;
}

static void hash_all(char **paths, size_t n, unsigned char (*digests)[32],
                      int nthreads, int quiet) {
    if (n == 0) return;
    if (nthreads < 1) nthreads = 1;
    if ((size_t)nthreads > n) nthreads = (int)n;

    atomic_size_t done_counter = 0;
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)nthreads);
    hash_job_t *jobs = malloc(sizeof(hash_job_t) * (size_t)nthreads);
    size_t chunk = (n + (size_t)nthreads - 1) / (size_t)nthreads;
    int launched = 0;

    for (int t = 0; t < nthreads; t++) {
        size_t start = (size_t)t * chunk;
        if (start >= n) break;
        size_t end = start + chunk;
        if (end > n) end = n;
        jobs[t].paths = paths;
        jobs[t].digests = digests;
        jobs[t].start = start;
        jobs[t].end = end;
        jobs[t].quiet = quiet;
        jobs[t].done_counter = &done_counter;
        pthread_create(&th[t], NULL, hash_worker, &jobs[t]);
        launched++;
    }

    if (!quiet) {
        size_t seen;
        do {
            struct timespec nap = {0, 50 * 1000 * 1000}; /* 50ms */
            nanosleep(&nap, NULL);
            seen = atomic_load(&done_counter);
            progress_bar("hashing", seen, n);
        } while (seen < n);
    }

    for (int t = 0; t < launched; t++) pthread_join(th[t], NULL);
    free(th);
    free(jobs);
}

/* ---------------------------------------------------------------- */
/* Greedy similarity ordering (O(n^2) comparisons, O(n) memory)     */
/* ---------------------------------------------------------------- */

static size_t *greedy_order(unsigned char (*digests)[32], size_t n, int quiet) {
    if (n == 0) return NULL;
    size_t *order = malloc(sizeof(size_t) * n);
    unsigned char *visited = calloc(n, 1);

    order[0] = 0;
    visited[0] = 1;
    size_t cur = 0;

    for (size_t step = 1; step < n; step++) {
        int best_score = -1000;
        size_t best_idx = 0;
        int found = 0;
        for (size_t i = 0; i < n; i++) {
            if (visited[i]) continue;
            int s = nilsimsa_compare(digests[cur], digests[i]);
            if (!found || s > best_score) {
                best_score = s;
                best_idx = i;
                found = 1;
            }
        }
        order[step] = best_idx;
        visited[best_idx] = 1;
        cur = best_idx;
        if (!quiet) progress_bar("ordering", step, n - 1);
    }
    free(visited);
    return order;
}

/* One full 2-opt sweep over the current order. For every pair of
 * non-adjacent "edges" (order[i],order[i+1]) and (order[j],order[j+1]),
 * check whether reversing the segment between them increases the
 * total similarity of the two boundary edges; if so, apply it
 * immediately (first-improvement strategy) and keep scanning.
 * Returns 1 if any improving move was applied this sweep, else 0.
 */
static int twoopt_pass(unsigned char (*digests)[32], size_t *order, size_t n,
                        const char *label, int quiet) {
    int improved = 0;
    if (n < 4) return 0;
    size_t last_i = n - 4; /* i ranges 0 .. n-4 inclusive (i+3 < n) */
    for (size_t i = 0; i + 3 < n; i++) {
        int score_i = nilsimsa_compare(digests[order[i]], digests[order[i + 1]]);
        for (size_t j = i + 2; j + 1 < n; j++) {
            int old_score = score_i +
                nilsimsa_compare(digests[order[j]], digests[order[j + 1]]);
            int new_score =
                nilsimsa_compare(digests[order[i]], digests[order[j]]) +
                nilsimsa_compare(digests[order[i + 1]], digests[order[j + 1]]);
            if (new_score > old_score) {
                size_t lo = i + 1, hi = j;
                while (lo < hi) {
                    size_t tmp = order[lo];
                    order[lo] = order[hi];
                    order[hi] = tmp;
                    lo++; hi--;
                }
                improved = 1;
                score_i = nilsimsa_compare(digests[order[i]], digests[order[i + 1]]);
            }
        }
        if (!quiet) progress_bar(label, i, last_i);
    }
    return improved;
}

/* Run up to max_sweeps full 2-opt sweeps, stopping early once a sweep
 * makes no improving move (i.e. it has converged to a local optimum). */
static void optimize_order(unsigned char (*digests)[32], size_t *order,
                            size_t n, int max_sweeps, int quiet) {
    if (max_sweeps <= 0) return;
    char label[24];
    for (int s = 0; s < max_sweeps; s++) {
        snprintf(label, sizeof(label), "sweep %d", s + 1);
        int improved = twoopt_pass(digests, order, n, label, quiet);
        if (!improved) {
            if (!quiet) {
                fprintf(stderr, "nilsort: converged after %d sweep(s)\n", s + 1);
            }
            break;
        }
    }
}

/* ---------------------------------------------------------------- */
/* Path collection                                                  */
/* ---------------------------------------------------------------- */

typedef struct {
    char **paths;
    size_t count;
    size_t cap;
} pathlist_t;

static void pathlist_init(pathlist_t *pl) {
    pl->paths = NULL;
    pl->count = 0;
    pl->cap = 0;
}

static void pathlist_push(pathlist_t *pl, const char *p) {
    if (pl->count == pl->cap) {
        pl->cap = pl->cap ? pl->cap * 2 : 1024;
        pl->paths = realloc(pl->paths, pl->cap * sizeof(char *));
    }
    pl->paths[pl->count++] = strdup(p);
}

static void scan_dir(const char *dir, pathlist_t *pl, int quiet) {
    DIR *d = opendir(dir);
    if (!d) {
        if (!quiet) fprintf(stderr, "nilsort: warning: cannot open %s: %s\n",
                             dir, strerror(errno));
        return;
    }
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, "..")) continue;
        char full[PATH_MAX];
        snprintf(full, sizeof(full), "%s/%s", dir, ent->d_name);
        struct stat st;
        if (lstat(full, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) {
            scan_dir(full, pl, quiet);
        } else if (S_ISREG(st.st_mode)) {
            pathlist_push(pl, full);
        }
    }
    closedir(d);
}

static void read_stdin_list(pathlist_t *pl) {
    char line[PATH_MAX];
    while (fgets(line, sizeof(line), stdin)) {
        size_t l = strlen(line);
        while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) line[--l] = 0;
        if (l) pathlist_push(pl, line);
    }
}

/* ---------------------------------------------------------------- */
/* main                                                              */
/* ---------------------------------------------------------------- */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <dir>\n"
        "       find /some/dir -type f | %s [options] -\n\n"
        "Options:\n"
        "  -b N   batch size (files per similarity-sort batch). Default:\n"
        "         unbounded -- entire input is one batch (max accuracy).\n"
        "  -o N   max 2-opt refinement sweeps, default 50 (0 disables)\n"
        "  -t N   number of hashing threads, default = number of CPUs\n"
        "  -q     quiet (suppress progress messages on stderr)\n"
        "  -v     print a memory/behavior estimate before running\n"
        "  -h     show this help\n",
        prog, prog);
}

int main(int argc, char **argv) {
    size_t batch_size = 0;   /* 0 = unbounded; resolved after we know pl.count */
    int batch_explicit = 0;
    int max_sweeps = 50;
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1) nthreads = 1;
    int quiet = 0;
    int verbose = 0;
    int opt;

    while ((opt = getopt(argc, argv, "b:o:t:qvh")) != -1) {
        switch (opt) {
            case 'b':
                batch_size = (size_t)strtoull(optarg, NULL, 10);
                batch_explicit = 1;
                break;
            case 'o': max_sweeps = atoi(optarg); break;
            case 't': nthreads = atoi(optarg); break;
            case 'q': quiet = 1; break;
            case 'v': verbose = 1; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) {
        usage(argv[0]);
        return 1;
    }
    if (batch_explicit && batch_size < 1) batch_size = 1;

    const char *src = argv[optind];
    pathlist_t pl;
    pathlist_init(&pl);

    if (!strcmp(src, "-")) {
        read_stdin_list(&pl);
    } else {
        struct stat st;
        if (stat(src, &st) != 0) {
            fprintf(stderr, "nilsort: cannot stat %s: %s\n", src, strerror(errno));
            return 1;
        }
        if (S_ISDIR(st.st_mode)) {
            scan_dir(src, &pl, quiet);
        } else {
            pathlist_push(&pl, src);
        }
    }

    if (pl.count == 0) {
        if (!quiet) fprintf(stderr, "nilsort: no files found\n");
        return 0;
    }

    if (!batch_explicit) batch_size = pl.count; /* unbounded: one batch total */
    size_t nbatches = (pl.count + batch_size - 1) / batch_size;

    if (verbose) {
        /* Memory is O(batch_size): 32 bytes/digest plus a visited byte.
           Time is O(batch_size^2) per pass (greedy construction, and
           again for each 2-opt sweep) -- that's what actually limits
           how large a batch is practical, not memory. */
        double digest_mb = (double)(batch_size * 32) / (1024.0 * 1024.0);
        double comparisons_millions =
            ((double)batch_size * (double)batch_size) / 1e6;
        fprintf(stderr,
            "nilsort: %zu files, batch size %zu -> %zu batch(es), %d thread(s)\n"
            "nilsort: digest memory per batch: ~%.1f MB\n"
            "nilsort: greedy-sort comparisons per batch: ~%.0f million\n"
            "nilsort: up to %d refinement sweep(s), ~%.0f million comparisons each\n",
            pl.count, batch_size, nbatches, nthreads,
            digest_mb, comparisons_millions, max_sweeps, comparisons_millions);
    }

    unsigned char (*digests)[32] = malloc(batch_size * sizeof(*digests));
    if (!digests) {
        fprintf(stderr, "nilsort: out of memory allocating digest buffer\n");
        return 1;
    }

    size_t done = 0;
    size_t batch_num = 0;
    while (done < pl.count) {
        size_t n = pl.count - done;
        if (n > batch_size) n = batch_size;
        char **batch_paths = pl.paths + done;
        batch_num++;

        if (!quiet) {
            fprintf(stderr, "nilsort: batch %zu/%zu: hashing %zu files...\n",
                    batch_num, nbatches, n);
        }
        hash_all(batch_paths, n, digests, nthreads, quiet);

        if (!quiet) {
            fprintf(stderr, "nilsort: batch %zu/%zu: ordering...\n",
                    batch_num, nbatches);
        }
        size_t *order = greedy_order(digests, n, quiet);
        optimize_order(digests, order, n, max_sweeps, quiet);
        for (size_t i = 0; i < n; i++) {
            printf("%s\n", batch_paths[order[i]]);
        }
        fflush(stdout);
        free(order);

        done += n;
    }

    free(digests);
    for (size_t i = 0; i < pl.count; i++) free(pl.paths[i]);
    free(pl.paths);
    return 0;
}
