#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/biopython_matrices_oracle.tsv.gz (audit R14-8, B09 F301) with the installed Biopython 1.88 `Bio.Align.substitution_matrices` and `Bio.Align.PairwiseAligner`:
section M: every cell of each of the 29 residue matrices (name, alphabet, cells in row order); section A: the optimal global / local score of seeded random pairs with each matrix and integer / fractional
affine gap scores (library convention GapOpen + k x GapExtend = Biopython open_gap_score GapOpen + GapExtend, extend_gap_score GapExtend).

  python3 scripts/alignment/generate_matrix_oracle.py
"""
import gzip
import os
import random

from Bio.Align import PairwiseAligner, substitution_matrices as sm

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
OUT = os.path.join(ROOT, "tests", "Seqeron", "Seqeron.Genomics.Tests", "TestData", "Alignment", "biopython_matrices_oracle.tsv.gz")
rnd = random.Random(20261124)
GAPS = [(-2.0, -1.0), (-1.5, -0.5), (-3.0, -0.25), (-0.4, -0.1), (-10.0, -0.5)]


def main():
    lines = ["kind\tname\ta\tb\tc\td\te\tf"]
    for name in sm.load():
        if name == "SCHNEIDER":
            continue
        m = sm.load(name)
        alphabet = "".join(m.alphabet)
        cells = [repr(float(m[a, b])) for a in alphabet for b in alphabet]
        lines.append("\t".join(["M", name, alphabet, " ".join(cells), "", "", "", ""]))
        letters = [c for c in alphabet if c.isalpha() and c not in "BZXJOUN" or name in ("TRANS", "HOXD70")] or list(alphabet)
        if name in ("BLASTN", "MEGABLAST", "NUC.4.4"):
            letters = list("ACGT")
        for _ in range(14):
            a = "".join(rnd.choice(letters) for _ in range(rnd.randint(4, 28)))
            b = "".join(rnd.choice(letters) for _ in range(rnd.randint(4, 28)))
            go, ge = rnd.choice(GAPS)
            for mode in ("global", "local"):
                al = PairwiseAligner()
                al.mode = mode
                al.substitution_matrix = m
                al.open_gap_score = go + ge
                al.extend_gap_score = ge
                score = al.score(a, b)
                lines.append("\t".join(["A", name, mode, repr(go), repr(ge), a, b, repr(float(score))]))
    with open(OUT, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as gz:
            gz.write(("\n".join(lines) + "\n").encode())
    print("rows", len(lines) - 1)


if __name__ == "__main__":
    main()
