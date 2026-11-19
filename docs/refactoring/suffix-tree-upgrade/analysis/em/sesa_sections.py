#!/usr/bin/env python3
"""Outside the shipped code (work package 1b.7e, gate G2): reads a SESA container with a decoder written from the text of ADR-003 only (no code of the product),
checks the header and the section table (the header CRC in pure Python), and the CRC-32C of every section: with --crc-c <binary> by the independent C program crc32c_ref (a table
implementation of the definition, fast), otherwise in pure Python (slow: use it on small files only); --no-crc skips the section CRCs. Prints or extracts the payloads.
  usage: sesa_sections.py <file> hash [--no-crc]          one line per section: kind, name, offset, length, SHA-256 of the payload
         sesa_sections.py <file> extract <prefix> [--no-psi]   writes <prefix>.pk .sa .lcp .ovf .cld .psi (what exists), <prefix>.pk.meta and <prefix>.stats as the E8 verifier reads them
"""
import hashlib, struct, subprocess, sys

NAMES = {1: 'TEXT', 2: 'ALPHABET', 3: 'SA', 4: 'LCP', 5: 'LCP_OVERFLOW', 6: 'CHILD', 11: 'BUILD_INFO', 12: 'PSI', 13: 'COUNTS'}
EXT = {1: 'pk', 3: 'sa', 4: 'lcp', 5: 'ovf', 6: 'cld', 12: 'psi'}

def make_table():
    t = []
    for i in range(256):
        c = i
        for _ in range(8):
            c = (c >> 1) ^ 0x82F63B78 if c & 1 else c >> 1
        t.append(c)
    return t
TABLE = make_table()

def crc32c(data, crc=0):
    c = crc ^ 0xFFFFFFFF
    for b in data:
        c = TABLE[(c ^ b) & 0xFF] ^ (c >> 8)
    return c ^ 0xFFFFFFFF

def read_header(f):
    page = f.read(4096)
    if len(page) != 4096:
        raise SystemExit('shorter than the header page')
    magic, major, minor, probe, n, nodes, size, toff, count, esize, lrsd, lrsr, posw, sigma = struct.unpack_from('<8sHHIqqqqiiqqii', page, 0)
    sha = page[80:112]
    flags, reserved, hcrc = struct.unpack_from('<qII', page, 112)
    assert magic == b'\x89SESA\r\n\x1a', 'magic %r' % magic
    assert (major, minor, probe) == (1, 0, 0x01020304) and toff == 128 and esize == 64, 'version/probe/table'
    sections = []
    for i in range(count):
        e = struct.unpack_from('<iiqqqiiqqII', page, 128 + 64 * i)
        sections.append(dict(kind=e[0], flags=e[1], offset=e[2], length=e[3], count=e[4], elem=e[5], codec=e[6], p0=e[7], p1=e[8], crc=e[9]))
    expect = crc32c(page[128:128 + 64 * count], crc32c(page[:124]))
    assert expect == hcrc, 'header CRC %08x expected %08x' % (hcrc, expect)
    return dict(n=n, nodes=nodes, size=size, lrsd=lrsd, lrsr=lrsr, posw=posw, sigma=sigma, sha=sha, flags=flags), sections

def chunks(f, offset, length, step=1 << 20):
    f.seek(offset)
    left = length
    while left > 0:
        b = f.read(min(step, left))
        if not b:
            raise SystemExit('short read')
        left -= len(b)
        yield b

def main():
    path, cmd = sys.argv[1], sys.argv[2]
    with open(path, 'rb') as f:
        header, sections = read_header(f)
        f.seek(0, 2)
        assert f.tell() == header['size'], 'FILE_SIZE %d, file %d' % (header['size'], f.tell())
        if cmd == 'hash':
            check = '--no-crc' not in sys.argv
            cbin = sys.argv[sys.argv.index('--crc-c') + 1] if '--crc-c' in sys.argv else None
            for s in sections:
                h = hashlib.sha256()
                crc = 0
                for b in chunks(f, s['offset'], s['length']):
                    h.update(b)
                    if check and not cbin:
                        crc = crc32c(b, crc)
                if check and cbin:
                    crc = int(subprocess.check_output([cbin, path, str(s['offset']), str(s['length'])]).decode().strip(), 16)
                note = ''
                if check:
                    note = ' crc ok' if crc == s['crc'] else ' CRC MISMATCH %08x != %08x' % (crc, s['crc'])
                print('%d %s %d %d %s%s' % (s['kind'], NAMES.get(s['kind'], '?'), s['offset'], s['length'], h.hexdigest(), note))
            print('header n=%d nodecount=%d lrsdepth=%d lrsrank=%d sigma=%d textsha=%s' % (header['n'], header['nodes'], header['lrsd'], header['lrsr'], header['sigma'], header['sha'].hex()))
        elif cmd == 'extract':
            prefix = sys.argv[3]
            by = {s['kind']: s for s in sections}
            for kind, ext in EXT.items():
                if kind == 12 and '--no-psi' in sys.argv:
                    continue
                if kind not in by:
                    if kind == 5:
                        open('%s.%s' % (prefix, ext), 'wb').close()
                    continue
                s = by[kind]
                with open('%s.%s' % (prefix, ext), 'wb') as out:
                    for b in chunks(f, s['offset'], s['length']):
                        out.write(b)
            letters = struct.unpack('<%dH' % header['sigma'], bytes(b''.join(chunks(f, by[2]['offset'], by[2]['length']))))
            counts = struct.unpack('<%dQ' % header['sigma'], bytes(b''.join(chunks(f, by[13]['offset'], by[13]['length']))))
            bits = by[1]['p0']
            lcp = by[4]
            meta = 'n=%d\nbits=%d\nsigma=%d\nletters_hex=%s\ncounts=%s\nsha256=%s\nname=sesa\n' % (header['n'], bits, header['sigma'], bytes(letters).hex().upper() if all(l < 256 for l in letters) else '', ','.join(map(str, counts)), header['sha'].hex())
            open(prefix + '.pk.meta', 'w').write(meta)
            maxlcp = '?'
            stats = 'nodecount=%d\ninternal=%d\nlrsdepth=%d\nlrsrank=%d\nn=%d\nsigma=%d\nlcpwidth=%d\noverflow=%d\nmaxlcp=%s\n' % (header['nodes'], header['nodes'] - header['n'] - 2, header['lrsd'], max(header['lrsr'], 0), header['n'], header['sigma'], lcp['elem'], lcp['p0'], header['lrsd'])
            open(prefix + '.stats', 'w').write(stats)
        else:
            raise SystemExit('unknown command')

main()
