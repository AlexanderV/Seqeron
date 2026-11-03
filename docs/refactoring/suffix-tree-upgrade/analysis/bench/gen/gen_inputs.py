#!/usr/bin/env python3
"""Deterministic input + query generator shared by EVERY benchmark subject (Seqeron, Go, and the
hand-off tools).  Pure Python 3 stdlib.  Re-running with the same SEED yields byte-identical files
(verify with:  cd <outdir> && sha256sum -c MANIFEST.sha256).

Usage:  python3 gen_inputs.py <outdir>

Layout written under <outdir>:
  data/<name>.txt        raw text, one line, NO trailing newline  (what Seqeron/Go/sdsl read)
  data/<name>.fa         FASTA (80 columns)  (BWA / Bowtie2 / libsais wrappers can use it)
  queries/<name>.present20.txt   100 000 20-mers sampled from the text (all present), one per line
  queries/<name>.absent20.txt    100 000 random 20-mers (mostly/entirely absent), one per line
  queries/<name>.locate8.txt     1 000 random 8-mers (many hits on DNA), one per line
  queries/<name>.locate4.txt     (protein only) 1 000 random 4-mers (8-mers are ~never present at sigma=20)
  queries/<name>.memq100k.txt    one 100 000 bp query = window of the text with 1 % substitutions
  MANIFEST.tsv / MANIFEST.sha256

Texts:
  dna_1m, dna_10m, dna_20m, dna_50m : nested prefixes of ONE seeded random ACGT stream (50 Mbp)
  rep_10m : 100 kbp random ACGT block; 100 copies, each copy = the ORIGINAL block with exactly 1 % of
            positions substituted (to one of the 3 other bases); concatenated = 10 Mbp
  prot_1m : 1 Mbp uniform random text over the 20 amino-acid letters ACDEFGHIKLMNPQRSTVWY
"""
import hashlib
import os
import random
import sys

SEED = "seqeron-st-bench-2026-11-02"
DNA = b"ACGT"
PROT = "ACDEFGHIKLMNPQRSTVWY"
DNA_TABLE = bytes(DNA[b & 3] for b in range(256))


def rng(tag):
    # str seeding goes through SHA-512 -> identical on every CPython 3.x
    return random.Random(f"{SEED}:{tag}")


def rand_dna(r, n):
    return r.getrandbits(8 * n).to_bytes(n, "little").translate(DNA_TABLE)


def mutate(r, s, rate, alphabet):
    """exactly round(rate*len) distinct positions substituted by a DIFFERENT letter"""
    b = bytearray(s)
    k = round(rate * len(b))
    for p in r.sample(range(len(b)), k):
        c = b[p]
        while True:
            d = alphabet[r.randrange(len(alphabet))]
            if d != c:
                break
        b[p] = d
    return bytes(b)


def write(path, data):
    with open(path, "wb") as f:
        f.write(data)


def write_lines(path, items):
    with open(path, "wb") as f:
        f.write(b"\n".join(items))
        f.write(b"\n")


def fasta(path, name, data):
    with open(path, "wb") as f:
        f.write(b">" + name.encode() + b"\n")
        for i in range(0, len(data), 80):
            f.write(data[i:i + 80] + b"\n")


def make_queries(outdir, name, text, alphabet, protein=False):
    n = len(text)
    r = rng(name + ":q")
    pres = [text[p:p + 20] for p in (r.randrange(n - 19) for _ in range(100_000))]
    ab = alphabet
    absent = [bytes(ab[r.randrange(len(ab))] for _ in range(20)) for _ in range(100_000)]
    loc8 = [bytes(ab[r.randrange(len(ab))] for _ in range(8)) for _ in range(1000)]
    write_lines(f"{outdir}/queries/{name}.present20.txt", pres)
    write_lines(f"{outdir}/queries/{name}.absent20.txt", absent)
    write_lines(f"{outdir}/queries/{name}.locate8.txt", loc8)
    if protein:
        loc4 = [bytes(ab[r.randrange(len(ab))] for _ in range(4)) for _ in range(1000)]
        write_lines(f"{outdir}/queries/{name}.locate4.txt", loc4)
    off = r.randrange(n - 100_000)
    q = mutate(r, text[off:off + 100_000], 0.01, ab)
    write(f"{outdir}/queries/{name}.memq100k.txt", q + b"\n")


def main():
    out = sys.argv[1]
    os.makedirs(out + "/data", exist_ok=True)
    os.makedirs(out + "/queries", exist_ok=True)
    texts = {}
    d50 = rand_dna(rng("dna50"), 50_000_000)
    for nm, n in (("dna_1m", 1_000_000), ("dna_10m", 10_000_000), ("dna_20m", 20_000_000), ("dna_50m", 50_000_000)):
        texts[nm] = (d50[:n], DNA, False)
    del d50
    rr = rng("rep")
    block = rand_dna(rr, 100_000)
    texts["rep_10m"] = (b"".join(mutate(rr, block, 0.01, DNA) for _ in range(100)), DNA, False)
    rp = rng("prot")
    texts["prot_1m"] = ("".join(rp.choices(PROT, k=1_000_000)).encode(), PROT.encode(), True)

    rows = []
    for nm, (t, ab, prot) in texts.items():
        write(f"{out}/data/{nm}.txt", t)
        if not prot:
            fasta(f"{out}/data/{nm}.fa", nm, t)
        make_queries(out, nm, t, ab, prot)
        print("wrote", nm, len(t), flush=True)
        rows.append((nm, len(t), "protein20" if prot else "dna4"))
    with open(out + "/MANIFEST.tsv", "w") as f:
        f.write("name\tn\talphabet\n")
        for r_ in rows:
            f.write("\t".join(map(str, r_)) + "\n")
    with open(out + "/MANIFEST.sha256", "w") as f:
        for sub in ("data", "queries"):
            for fn in sorted(os.listdir(f"{out}/{sub}")):
                h = hashlib.sha256()
                with open(f"{out}/{sub}/{fn}", "rb") as g:
                    for blk in iter(lambda: g.read(1 << 22), b""):
                        h.update(blk)
                f.write(f"{h.hexdigest()}  {sub}/{fn}\n")


if __name__ == "__main__":
    main()
