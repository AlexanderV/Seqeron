#!/usr/bin/env python3
"""Fault injection on REAL outputs of e6 (ADR-004 E8): each fault is written into the files in place, e8-verify runs, the file is restored byte for byte. A fault counts as detected when e8-verify exits non-zero.
usage: python3 -I inject_faults.py <e6 output prefix> <input text> <rows file> <case> [dotnet dll] [extra e8-verify args...]
The input text is COPIED (the original is never touched). Needs about 1 x n bytes of free disk for that copy. Rows: R<TAB>E8<TAB>case<TAB>inject<TAB>fault<TAB>detected|MISSED<TAB><TAB>failed checks."""
import os, random, shutil, struct, subprocess, sys, time

prefix, inp, rows, case = sys.argv[1:5]
dll = sys.argv[5] if len(sys.argv) > 5 else os.environ.get('E8_DLL', 'SuffixTree.EmProto.dll')
extra = sys.argv[6:]
work_input = prefix + '.inject-input.txt'
shutil.copyfile(inp, work_input)
rng = random.Random(8)


def read(path, off, n):
    with open(path, 'rb') as f:
        f.seek(off)
        return f.read(n)


def write(path, off, data):
    fd = os.open(path, os.O_WRONLY)
    try:
        os.pwrite(fd, data, off)
    finally:
        os.close(fd)


def u32(path, i):
    return struct.unpack('<I', read(path, 4 * i, 4))[0]


def put32(path, i, v):
    write(path, 4 * i, struct.pack('<I', v & 0xFFFFFFFF))


def verify(with_input=True):
    t = time.time()
    p = subprocess.run(['dotnet', dll, 'e8-verify', '--prefix', prefix] + (['--input', work_input] if with_input else []) + ['--quiet', '--no-hash'] + extra, capture_output=True, text=True)
    failed = [l.split('\t')[3] for l in p.stdout.splitlines() if l.split('\t')[0:1] == ['R'] and len(l.split('\t')) > 4 and l.split('\t')[4] == 'FAIL' and l.split('\t')[3] != 'verify_all']
    return p.returncode, failed, time.time() - t


def report(fault, rc, failed, secs):
    verdict = 'detected' if rc != 0 else 'MISSED'
    line = '\t'.join(['R', 'E8', case, 'inject', fault, verdict, '', 'failed checks: %s (exit code %d, %.1f s)' % (','.join(failed) or 'none', rc, secs)])
    print(line, flush=True)
    with open(rows, 'a') as f:
        f.write(line + '\n')
    return rc != 0


n = os.path.getsize(prefix + '.sa') // 4 - 1
width = os.path.getsize(prefix + '.lcp') // (n + 1)
ovf_pairs = os.path.getsize(prefix + '.ovf') // 8
results = []


def with_bytes(path, off, data, name):
    """overwrite bytes, verify, restore"""
    old = read(path, off, len(data))
    write(path, off, data)
    try:
        results.append(report(name, *verify()))
    finally:
        write(path, off, old)


def swap32(path, i, j, name):
    a, b = u32(path, i), u32(path, j)
    put32(path, i, b)
    put32(path, j, a)
    try:
        results.append(report(name, *verify()))
    finally:
        put32(path, i, a)
        put32(path, j, b)


sa = prefix + '.sa'
r1, r2 = rng.randrange(1, n // 2), rng.randrange(n // 2, n)
swap32(sa, r1, r2, 'swap two SA entries (ranks %d and %d)' % (r1, r2))
r = rng.randrange(1, n - 1)
swap32(sa, r, r + 1, 'swap two adjacent SA entries (ranks %d, %d)' % (r, r + 1))
r = rng.randrange(1, n - 10)
with_bytes(sa, 4 * r, struct.pack('<I', u32(sa, r + 5)), 'duplicate an SA entry (rank %d := rank %d)' % (r, r + 5))
r = rng.randrange(1, n)
with_bytes(sa, 4 * r, struct.pack('<I', n + 1), 'SA entry above n (rank %d := n + 1)' % r)
with_bytes(sa, 0, struct.pack('<I', n - 1), 'SA[0] is not n')

lcp = prefix + '.lcp'
for delta in (1, -1):
    while True:
        r = rng.randrange(2, n)
        v = read(lcp, width * r, 1)[0]
        if 5 < v < 200 and width == 1:
            break
        if width != 1:
            v = 7
            break
    with_bytes(lcp, width * r, bytes([v + delta]) + (b'' if width == 1 else read(lcp, width * r + 1, width - 1)), 'LCP value %+d (rank %d, %d -> %d)' % (delta, r, v, v + delta))

if ovf_pairs:
    ovf = prefix + '.ovf'
    k = rng.randrange(ovf_pairs)
    v = u32(ovf, 2 * k + 1)
    with_bytes(ovf, 8 * k + 4, struct.pack('<I', v + 1), 'overflow value +1 (pair %d)' % k)
    with_bytes(ovf, 8 * k + 4, struct.pack('<I', v - 1), 'overflow value -1 (pair %d)' % k)
    rk = u32(ovf, 2 * k)
    with_bytes(ovf, 8 * k, struct.pack('<I', rk + 1), 'overflow rank +1 (pair %d)' % k)

cld = prefix + '.cld'
for what in ('xor 1 of a set cell', 'a set cell := 0', 'an empty cell := 1', 'a set cell := 2^32 - 1'):
    while True:
        k = rng.randrange(0, n + 1)
        v = u32(cld, k)
        if (what == 'an empty cell := 1') == (v == 0):
            break
    new = {'xor 1 of a set cell': v ^ 1, 'a set cell := 0': 0, 'an empty cell := 1': 1, 'a set cell := 2^32 - 1': 0xFFFFFFFF}[what]
    with_bytes(cld, 4 * k, struct.pack('<I', new), 'child table: %s (cell %d)' % (what, k))
# the last cell and the first cell
with_bytes(cld, 0, struct.pack('<I', u32(cld, 0) + 1), 'child table: cell 0 (the root) + 1')

pk = prefix + '.pk'
off = rng.randrange(10, os.path.getsize(pk) - 100)
with_bytes(pk, off, bytes([read(pk, off, 1)[0] ^ 0x10]), 'packed text: one bit flipped (byte %d)' % off)

# truncation, deletion, extension
for ext in ('sa', 'lcp', 'cld', 'pk'):
    path = prefix + '.' + ext
    size = os.path.getsize(path)
    last = read(path, size - 1, 1)
    os.truncate(path, size - 1)
    try:
        results.append(report('truncate .%s by one byte' % ext, *verify()))
    finally:
        os.truncate(path, size)
        write(path, size - 1, last)
    with open(path, 'ab') as f:
        f.write(b'\0')
    try:
        results.append(report('extend .%s by one byte' % ext, *verify()))
    finally:
        os.truncate(path, size)
if ovf_pairs:
    path = prefix + '.ovf'
    size = os.path.getsize(path)
    last = read(path, size - 8, 8)
    os.truncate(path, size - 8)
    try:
        results.append(report('truncate .ovf by one pair', *verify()))
    finally:
        os.truncate(path, size)
        write(path, size - 8, last)
os.rename(prefix + '.lcp', prefix + '.lcp.bak')
try:
    results.append(report('delete .lcp', *verify()))
finally:
    os.rename(prefix + '.lcp.bak', prefix + '.lcp')

# off by one in n (metadata, statistics)
for ext, key in (('pk.meta', 'n'), ('stats', 'n'), ('stats', 'nodecount'), ('stats', 'lrsrank'), ('stats', 'internal'), ('stats', 'maxlcp')):
    path = prefix + '.' + ext
    old = open(path, 'rb').read()
    lines = old.decode('latin1').split('\n')
    for i, l in enumerate(lines):
        if l.startswith(key + '='):
            lines[i] = key + '=' + str(int(l.split('=')[1]) + 1)
    open(path, 'wb').write('\n'.join(lines).encode('latin1'))
    try:
        results.append(report('%s: %s + 1' % (ext, key), *verify()))
    finally:
        open(path, 'wb').write(old)

# the input text
with_bytes(work_input, rng.randrange(0, n), b'Z', 'input text: one symbol changed')

os.remove(work_input)
rc, failed, secs = verify(False)
ok = rc == 0
line = '\t'.join(['R', 'E8', case, 'inject', 'restored files verify again', 'pass' if ok else 'FAIL', '', 'after restoring every file: exit code %d (no input check: the input copy is gone)' % rc])
print(line)
with open(rows, 'a') as f:
    f.write(line + '\n')
line = '\t'.join(['R', 'E8', case, 'inject', 'summary', '%d of %d detected' % (sum(results), len(results)), '', 'MISSED = the verifier passed a corrupt set of files'])
print(line)
with open(rows, 'a') as f:
    f.write(line + '\n')
sys.exit(0 if all(results) and ok else 1)
