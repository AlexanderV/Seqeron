#!/usr/bin/env python3
"""Streaming, deterministic generators of the hand-off package (work package 1b.9f): a text of any length is written to a file in blocks, never held in memory, and its SHA-256 and letter counts are
computed on the way. The texts are the ones the benchmark and the earlier packages use, byte for byte (the test compares them with the in-memory generators):
  dna    random ACGT, the stream of gen_inputs.py (tag dna50): the first n symbols are dna_10m / dna_100m of the benchmark and the files of 1b.9d / 1b.9e
  rep    a repetitive text: a block of 100,000 random symbols and copies of it, each with 1 % of its symbols substituted (rep_10m of the benchmark is 100 such copies, rep_100m of 1b.9e 1,000)
  prot   a protein-like text (make_1b9e_inputs.py): the composition of UniProt, fresh runs of 100-2000 residues, in 12 % of the steps a copy of an earlier segment (30-3000 residues, 5 % substituted)
The genome-like text (E3-E8) is made by `e8-gen` of the program (the entry script calls it). Speed measured on the 4-vCPU x86 box, SHA-256 and letter counts included: dna 100 Mbp in 1.9 s, rep 10 Mbp
in 0.2 s, prot 10 Mbp in 1.2 s; ESTIMATE for 3x10^9 (linear, not measured): dna about 1 minute, rep 1-2 minutes, prot about 6 minutes, plus the write of the file.
usage: s1b_gen.py --kind dna|rep|prot --n N[k|m|g] --out FILE [--force] [--chunk-mb 16]"""
import hashlib
import importlib.util
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
AMINO = 'LAGVESIKRDTPNQFYMHCW'
WEIGHT = [9.9, 8.3, 7.1, 6.9, 6.8, 6.6, 5.9, 5.8, 5.5, 5.5, 5.3, 4.7, 4.1, 3.9, 3.9, 2.9, 2.4, 2.3, 1.4, 1.1]


def load_gi():
    for cand in (os.path.join(HERE, '..', 'bench', 'gen', 'gen_inputs.py'), os.path.join(HERE, 'gen_inputs.py')):
        if os.path.exists(cand):
            spec = importlib.util.spec_from_file_location('gi', cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    sys.exit('s1b_gen: gen_inputs.py of the benchmark was not found next to this script or in ../bench/gen')


def parse_n(text):
    unit = {'k': 10**3, 'm': 10**6, 'g': 10**9}.get(text[-1:].lower(), 1)
    return int(float(text[:-1]) * unit) if unit > 1 else int(text)


def dna_blocks(gi, n, chunk):
    r = gi.rng('dna50')
    done = 0
    chunk -= chunk % 4                       # a multiple of four bytes: getrandbits(8 * chunk) then consumes whole 32-bit words, so the stream does not depend on the chunk size
    while done < n:
        take = min(chunk, n - done)
        want = take + (-take) % 4            # the last block draws a whole number of words as well and is cut
        yield r.getrandbits(8 * want).to_bytes(want, 'little').translate(gi.DNA_TABLE)[:take]
        done += take


def rep_blocks(gi, n):
    rr = gi.rng('rep100')
    block = gi.rand_dna(rr, 100_000)
    done = 0
    while done < n:
        copy = gi.mutate(rr, block, 0.01, gi.DNA)
        yield copy[:n - done]
        done += len(copy)


def prot_blocks(n, out_path):
    """the protein-like text, streamed: the copies of earlier segments are read back from the file that is being written (pread), so memory stays at one segment"""
    r = random.Random('1b9e:prot_10m')
    fd = os.open(out_path, os.O_RDWR | os.O_CREAT | os.O_TRUNC, 0o644)
    total = 0
    try:
        while total < n:
            if total > 5000 and r.random() < 0.12:
                length = r.randint(30, 3000)
                start = r.randrange(0, total - length)
                seg = bytearray(os.pread(fd, length, start))
                for p in r.sample(range(length), max(1, length // 20)):
                    seg[p] = ord(r.choices(AMINO, WEIGHT)[0])
                piece = bytes(seg)
            else:
                piece = ''.join(r.choices(AMINO, WEIGHT, k=r.randint(100, 2000))).encode()
            os.pwrite(fd, piece, total)       # the file holds everything that the original keeps in memory, including what is cut at n at the end
            total += len(piece)
            yield piece[:max(0, n - (total - len(piece)))]
        os.ftruncate(fd, n)
    finally:
        os.close(fd)


def main(argv):
    opts = {}
    force = False
    i = 0
    while i < len(argv):
        k = argv[i]
        if k == '--force':
            force = True
            i += 1
        elif k.startswith('--') and i + 1 < len(argv):
            opts[k[2:]] = argv[i + 1]
            i += 2
        else:
            print(__doc__)
            return 2
    kind, out = opts.get('kind'), opts.get('out')
    if kind not in ('dna', 'rep', 'prot') or not out or 'n' not in opts:
        print(__doc__)
        return 2
    try:
        n = parse_n(opts['n'])
    except ValueError:
        print('s1b_gen: --n is not a number')
        return 2
    if n <= 0:
        print('s1b_gen: --n must be positive')
        return 2
    if os.path.exists(out) and not force:
        print(f"s1b_gen: '{out}' exists (use --force to replace it)")
        return 2
    chunk = int(float(opts.get('chunk-mb', 16)) * 1024 * 1024)
    part = out + '.part'
    sha = hashlib.sha256()
    counts = {}
    gi = load_gi() if kind in ('dna', 'rep') else None
    blocks = dna_blocks(gi, n, chunk) if kind == 'dna' else rep_blocks(gi, n) if kind == 'rep' else None
    try:
        if kind == 'prot':
            # the protein text writes through its own descriptor (it reads its earlier output back); the hash and the counts are taken from the finished file
            for _ in prot_blocks(n, part):
                pass
            with open(part, 'rb') as fh:
                while True:
                    b = fh.read(1 << 24)
                    if not b:
                        break
                    sha.update(b)
                    for c in set(b):
                        counts[c] = counts.get(c, 0) + b.count(bytes([c]))
        else:
            with open(part, 'wb') as fh:
                for b in blocks:
                    fh.write(b)
                    sha.update(b)
                    for c in set(b):
                        counts[c] = counts.get(c, 0) + b.count(bytes([c]))
                fh.flush()
                os.fsync(fh.fileno())
        os.replace(part, out)
    except BaseException:
        if os.path.exists(part):
            os.remove(part)
        raise
    letters = ''.join(chr(c) for c in sorted(counts))
    with open(out + '.sums', 'w') as fh:
        fh.write(f"format=s1b-gen-v1\nkind={kind}\ngenerator=s1b_gen.py ({'gen_inputs.py stream dna50' if kind == 'dna' else 'gen_inputs.py block rep100' if kind == 'rep' else 'make_1b9e_inputs.py protein_like'})\n"
                 f"n={n}\nalphabet={letters}\ncounts={','.join(str(counts[ord(c)]) for c in letters)}\nsha256={sha.hexdigest()}\n")
    print(f'R\tS1B\tgen\t0\tgen_sha256\t{sha.hexdigest()}\t\t{kind} n={n}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
