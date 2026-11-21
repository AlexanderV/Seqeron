#!/usr/bin/env python3
"""Inputs of work package 1b.9e (the kill campaign): one byte per symbol, deterministic (str seeding goes through SHA-512, identical on every CPython 3.x).
  prot_10m   a protein-like text of 10^7 symbols: the composition of UniProt, fresh runs of 100-2000 residues, and in 12 % of the steps a copy of an earlier segment (30-3000 residues, 5 % substituted)
  dna_100m   the first 10^8 symbols of the stream of gen_inputs.py (tag dna50): the file of the benchmark, as 1b.9d used it (its first 5x10^7 symbols are dna_50m)
  rep_100m   a repetitive text of 10^8 symbols: 1000 copies of a block of 100,000 random symbols, each with 1 % of its symbols substituted (rep_10m is 100 such copies)
usage: make_1b9e_inputs.py <output directory> [prot_10m] [dna_100m] [rep_100m]   (no name = all three)"""
import importlib.util
import os
import random
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("gi", os.path.join(here, "..", "bench", "gen", "gen_inputs.py"))
gi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gi)

AMINO = "LAGVESIKRDTPNQFYMHCW"
WEIGHT = [9.9, 8.3, 7.1, 6.9, 6.8, 6.6, 5.9, 5.8, 5.5, 5.5, 5.3, 4.7, 4.1, 3.9, 3.9, 2.9, 2.4, 2.3, 1.4, 1.1]


def protein_like(n):
    r = random.Random("1b9e:prot_10m")
    out = bytearray()
    while len(out) < n:
        if len(out) > 5000 and r.random() < 0.12:
            length = r.randint(30, 3000)
            start = r.randrange(0, len(out) - length)
            seg = bytearray(out[start:start + length])
            for p in r.sample(range(length), max(1, length // 20)):
                seg[p] = ord(r.choices(AMINO, WEIGHT)[0])
            out += seg
        else:
            out += "".join(r.choices(AMINO, WEIGHT, k=r.randint(100, 2000))).encode()
    return bytes(out[:n])


def main():
    out = sys.argv[1]
    names = sys.argv[2:] or ["prot_10m", "dna_100m", "rep_100m"]
    os.makedirs(out, exist_ok=True)
    for name in names:
        path = os.path.join(out, name + ".txt")
        if name == "prot_10m":
            data = protein_like(10_000_000)
        elif name == "dna_100m":
            data = gi.rand_dna(gi.rng("dna50"), 100_000_000)
        elif name == "rep_100m":
            rr = gi.rng("rep100")
            block = gi.rand_dna(rr, 100_000)
            data = b"".join(gi.mutate(rr, block, 0.01, gi.DNA) for _ in range(1000))
        else:
            sys.exit("unknown input " + name)
        with open(path, "wb") as f:
            f.write(data)
        print("wrote", path, len(data), flush=True)


if __name__ == "__main__":
    main()
