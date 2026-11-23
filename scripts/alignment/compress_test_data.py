#!/usr/bin/env python3
"""Compresses the large oracle files of the test project (audit R14-20, B09 F299).

Every text oracle of tests/Seqeron/Seqeron.Genomics.Tests/TestData above 1 MiB is stored as `<name>.gz` (deterministic gzip: level 9, mtime 0, no file name) and embedded as the `.gz`; the tests open it by the
name of the PLAIN file through `TestDataResource.Open`, which decompresses with `GZipStream`. A generator writes the plain file; run this script afterwards:

  python3 scripts/alignment/compress_test_data.py            gzip every plain file above 1 MiB, delete it, point the csproj at the .gz
  python3 scripts/alignment/compress_test_data.py --check    exit 1 and list every plain text oracle above 1 MiB (what the contract test TestData_Size_Tests also asserts)
"""
import gzip
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
PROJECT = os.path.join(ROOT, "tests", "Seqeron", "Seqeron.Genomics.Tests")
DATA = os.path.join(PROJECT, "TestData")
CSPROJ = os.path.join(PROJECT, "Seqeron.Genomics.Tests.csproj")
LIMIT = 1 << 20
TEXT = (".tsv", ".csv", ".jsonl", ".json", ".txt", ".fa", ".fasta")


def plain_large():
    out = []
    for dirpath, _, files in os.walk(DATA):
        for f in files:
            p = os.path.join(dirpath, f)
            if f.lower().endswith(TEXT) and os.path.getsize(p) > LIMIT:
                out.append(p)
    return sorted(out)


def main(argv):
    big = plain_large()
    if "--check" in argv:
        for p in big:
            print("plain oracle above 1 MiB: %s (%.1f MB)" % (os.path.relpath(p, ROOT), os.path.getsize(p) / 1048576))
        return 1 if big else 0
    csproj = open(CSPROJ, encoding="utf-8").read()
    for p in big:
        with open(p, "rb") as src, open(p + ".gz", "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as gz:
                gz.write(src.read())
        os.remove(p)
        rel = os.path.relpath(p, PROJECT).replace("/", "\\")
        csproj = csproj.replace('Include="%s"' % rel, 'Include="%s.gz"' % rel)
        print("gzip %s -> %.1f MB" % (os.path.relpath(p, ROOT), os.path.getsize(p + ".gz") / 1048576))
    open(CSPROJ, "w", encoding="utf-8", newline="").write(csproj)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
