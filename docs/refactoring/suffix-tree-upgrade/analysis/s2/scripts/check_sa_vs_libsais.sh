#!/usr/bin/env bash
# The checksum check of the S2 gate (ADR-005 sections 8, 9 S2.15): SHA-256 of the libsais suffix array of a text file (little-endian int32, SA[1..n]) against SHA-256 of SA[1..n] of the container built from it.
# libsais is compiled from its source tree (the clone of analysis/bench/run_1b9_gate.sh, commit 1aae6344b8a58f00626f001b1fe33bdbecd5e469) with sa_dump.c; it needs n x 4 bytes of RAM and of disk for the dump (2 GB at 5x10^8; 35 s).
# usage: check_sa_vs_libsais.sh <text file> <container> <libsais source dir>        prints two lines; they must carry the same hash
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
T=${1:?text}; C=${2:?container}; L=${3:?libsais source dir}
W=$OUT/libsais; mkdir -p "$W"
[ -x "$W/sa_dump" ] || gcc -O3 -I "$L/include" "$L/src/libsais.c" "$S2/scripts/sa_dump.c" -o "$W/sa_dump"
"$W/sa_dump" "$T" "$W/sa.bin"
echo "libsais_sa_sha256 $(sha256sum "$W/sa.bin" | cut -c1-64)"
rm -f "$W/sa.bin"
python3 - "$C" <<'PY'
import hashlib, mmap, struct, sys
f = open(sys.argv[1], 'rb'); m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
n = struct.unpack_from('<Q', m, 16)[0]; count = struct.unpack_from('<I', m, 48)[0]
for i in range(count):
    e = 128 + 64 * i
    kind, _ = struct.unpack_from('<II', m, e)
    off, _length = struct.unpack_from('<QQ', m, e + 8)
    if kind == 3: sa = off
h = hashlib.sha256(); pos = sa + 4; end = sa + 4 + 4 * n
while pos < end:
    q = min(end, pos + (1 << 26)); h.update(m[pos:q]); pos = q
print('container_sa_1_to_n_sha256', h.hexdigest(), 'n', n)
PY
