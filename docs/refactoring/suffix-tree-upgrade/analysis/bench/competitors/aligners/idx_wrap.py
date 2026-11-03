#!/usr/bin/env python3
"""Run an index-building command and print R lines: total bytes of the files matching the output prefix.
usage: idx_wrap.py <out-prefix> <n-chars> -- <cmd...>
benchlib.measure wraps this script, so wall time / peak RSS / CPU of the children are captured by wait4."""
import glob
import os
import subprocess
import sys

prefix, n = sys.argv[1], int(sys.argv[2])
cmd = sys.argv[4:]
for f in glob.glob(prefix + "*"):
    os.remove(f)
rc = subprocess.call(cmd, stdout=sys.stderr, stderr=sys.stderr)
if rc != 0:
    sys.exit(rc)
files = sorted(glob.glob(prefix + "*"))
tot = sum(os.path.getsize(f) for f in files)
print(f"R\tindex\tindex_bytes\t{tot}\tbytes\t{len(files)} files: " + ",".join(os.path.basename(f) for f in files)[:300])
print(f"R\tindex\tbytes_per_char\t{tot / n!r}\tbytes/char\tall index files / n")
for f in files:
    os.remove(f)
