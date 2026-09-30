#!/usr/bin/env python3
"""
list_to_gensquashfs_sort.py - turn a y4mstore-style sort list (one path per
line, in the exact order to archive) into a gensquashfs --sort-file (one
"path priority" pair per line), so gensquashfs preserves that order exactly.

    ./list_to_gensquashfs_sort.py component-y4m-nilsimsa.txt yuv-elementary/ > component.sort
    gensquashfs out.sqfs --pack-dir yuv-elementary --sort-file component.sort

Unlike mksquashfs's -sort (a 16-bit priority field: only 65,536 distinct
values, so a list of hundreds of thousands of files has to be bucketed,
losing exact order within each bucket -- see list_to_sqfs_sort.py's header
for that whole story), gensquashfs's --sort-file priority is a 64-bit
number, confirmed directly in its own manual page: "sorted by priority,
with lower values preceding larger ones." That range comfortably holds one
unique priority per file at any size this project produces, so this script
just numbers the lines 0, 1, 2, ... in order -- no bucketing, no ties, no
direction to guess, and this one IS confirmed straight from gensquashfs's
own sort-file-format documentation, not inferred.

TWO THINGS THIS FORMAT REQUIRES THAT mksquashfs's DIDN'T
----------------------------------------------------------
1. Field order is PRIORITY FIRST, then filename -- the reverse of
   mksquashfs's "path priority". Confirmed directly from gensquashfs's own
   format documentation: "one entry per line, consisting of a numeric
   priority and a filename."

2. The filename is matched against the file's path INSIDE the image, not
   against your original input list's path -- confirmed the same way:
   "matched against the actual path of the file in the SquashFS file in
   the resulting image. It is not matched against the input path, which
   may differ." Since --pack-dir becomes the image's root, whatever prefix
   you pass as PACK_DIR here gets stripped from every line, the same way
   --pack-dir strips it when gensquashfs builds the image. Get this wrong
   and every line silently fails to match -- not an error, just every file
   quietly falling back to gensquashfs's undocumented default priority of
   0, with no warning that your ordering was never applied at all.

Still worth a quick sanity check before a real run:
    mkdir /tmp/gensqfs_dircheck && cd /tmp/gensqfs_dircheck
    touch aaa_first bbb_second
    printf '0 aaa_first\\n1 bbb_second\\n' > order.sort
    gensquashfs test.sqfs --pack-dir . --sort-file order.sort
    unsquashfs -ls test.sqfs
aaa_first (priority 0) should list before bbb_second (priority 1).
"""
import sys


def main():
    if len(sys.argv) not in (2, 3) or sys.argv[1] in ("-h", "--help"):
        print(f"usage: {sys.argv[0]} SORTLIST [PACK_DIR] > SORTFILE", file=sys.stderr)
        print("       SORTLIST: one path per line, in the order to archive", file=sys.stderr)
        print("       PACK_DIR: the prefix that will become --pack-dir's root (stripped", file=sys.stderr)
        print("                 from every line, matching what --pack-dir strips); omit", file=sys.stderr)
        print("                 only if SORTLIST's paths are already pack-dir-relative", file=sys.stderr)
        print("       SORTFILE: written to stdout, for gensquashfs --sort-file", file=sys.stderr)
        return 1
    prefix = sys.argv[2].rstrip("/") + "/" if len(sys.argv) == 3 else None

    with open(sys.argv[1], "r", encoding="utf-8", errors="strict") as f:
        paths = [line.rstrip("\n") for line in f if line.strip()]

    if not paths:
        print("list_to_gensquashfs_sort: input is empty; nothing to convert", file=sys.stderr)
        return 1

    stripped = 0
    for i, path in enumerate(paths):
        if "\n" in path:
            print(f"list_to_gensquashfs_sort: refusing: a path contains a newline: {path!r}", file=sys.stderr)
            return 1
        if " " in path:
            print(f"list_to_gensquashfs_sort: refusing: the sort-file format needs quoting for a "
                  f"path containing a space, not implemented here: {path!r}", file=sys.stderr)
            return 1
        rel = path
        if prefix is not None:
            if not path.startswith(prefix):
                print(f"list_to_gensquashfs_sort: refusing: {path!r} does not start with "
                      f"the given PACK_DIR prefix {prefix!r} -- check you passed the right one", file=sys.stderr)
                return 1
            rel = path[len(prefix):]
            stripped += 1
        print(f"{i} {rel}")

    print(f"list_to_gensquashfs_sort: {len(paths)} paths, each given its own priority "
          f"(0 to {len(paths) - 1}) -- exact order, no bucketing needed" +
          (f"; {prefix!r} prefix stripped from all {stripped}" if prefix is not None else
           "; NO prefix stripped -- make sure SORTLIST's paths already match what --pack-dir will see"),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
