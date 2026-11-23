#!/usr/bin/env python3
"""Work package 1b.9: the bytes per input character of a SESA container, by section (the section table of the header page, layout of ADR-003 section 2 as decoded by analysis/em/sesa_sections.py).
usage: sections_1b9.py <container>...   one block per file: every section with its payload bytes and bytes per character, the file size, and
  Core  = TEXT + ALPHABET + COUNTS + SA + LCP + LCP_OVERFLOW + CHILD + BUILD_INFO (the container written with PSI off), payload bytes
  Full  = Core + PSI (what the default container holds)
  Full+ = Full-in-the-sense-of-ADR-003 section 5: Core + the navigation pack that the engine builds on first use in memory (ISA 4.000 B/char + the LCP pyramid 0.129 B/char of `nav`, nothing of it is written to the file)
The header page (4,096 bytes) and the alignment padding between the sections are in the file size and in no section."""
import struct
import sys

NAMES = {1: "TEXT", 2: "ALPHABET", 3: "SA", 4: "LCP", 5: "LCP_OVERFLOW", 6: "CHILD", 11: "BUILD_INFO", 12: "PSI", 13: "COUNTS"}
NAV_PACK_PER_CHAR = 4.0 + 0.129


def main(paths):
    for path in paths:
        with open(path, "rb") as f:
            page = f.read(4096)
            f.seek(0, 2)
            size = f.tell()
        magic, major, minor, probe, n, nodes, fsize, toff, count = struct.unpack_from("<8sHHIqqqqi", page, 0)
        assert magic == b"\x89SESA\r\n\x1a" and toff == 128 and fsize == size, "not a SESA container of the expected size"
        secs = {}
        for i in range(count):
            kind, flags, offset, length = struct.unpack_from("<iiqq", page, 128 + 64 * i)
            secs[kind] = length
        print(f"# {path.rsplit('/', 1)[-1]}  n={n}  file {size} bytes = {size / n:.4f} B/char")
        for kind, length in sorted(secs.items()):
            print(f"  {NAMES.get(kind, kind):<13}{length:>14,} B  {length / n:8.4f} B/char")
        core = sum(v for k, v in secs.items() if k != 12)
        full = sum(secs.values())
        print(f"  Core          {core:>14,} B  {core / n:8.4f} B/char   (sections) ; file without PSI about {(size - secs.get(12, 0)) / n:.4f}")
        print(f"  Full          {full:>14,} B  {full / n:8.4f} B/char   (sections, with PSI) ; file {size / n:.4f}")
        print(f"  Full+nav      {'':>14}    {core / n + NAV_PACK_PER_CHAR:8.4f} B/char   (Core + ISA + LCP pyramid built in memory on first use)")


main(sys.argv[1:])
