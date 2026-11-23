#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/blast_significance_oracle.tsv.gz (audit R14-7, B09 F300) with the installed BLAST+ (blastp / blastn 2.12.0, `-subject`, `-outfmt 5`):
the score, bit score and E-value of every HSP of seeded related random pairs, for the standard protein matrices with several gap costs (blastp -comp_based_stats 0 -seg no) and for blastn reward / penalty /
gap-cost combinations (-dust no). Columns: kind, matrix_or_reward, penalty, gap_existence, gap_extension, query_length, subject_length, score, bit_score, evalue.

  python3 scripts/alignment/generate_significance_oracle.py
"""
import gzip
import os
import random
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
OUT = os.path.join(ROOT, "tests", "Seqeron", "Seqeron.Genomics.Tests", "TestData", "Alignment", "blast_significance_oracle.tsv.gz")
AA = "ACDEFGHIKLMNPQRSTVWY"
rnd = random.Random(20261123)


def rand(alpha, n):
    return "".join(rnd.choice(alpha) for _ in range(n))


def mutate(alpha, s, p):
    return "".join(rnd.choice(alpha) if rnd.random() < p else c for c in s)


def related_pair(alpha, mut):
    core = rand(alpha, rnd.randint(30, 90))
    q = rand(alpha, rnd.randint(10, 120)) + core + rand(alpha, rnd.randint(10, 120))
    s = rand(alpha, rnd.randint(10, 200)) + mutate(alpha, core, mut) + rand(alpha, rnd.randint(10, 200))
    return q, s


def run(cmd, q, s):
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "q.fa"), "w").write(">q\n" + q + "\n")
        open(os.path.join(d, "s.fa"), "w").write(">s\n" + s + "\n")
        out = subprocess.run(cmd + ["-query", os.path.join(d, "q.fa"), "-subject", os.path.join(d, "s.fa"), "-outfmt", "5", "-evalue", "10"], capture_output=True, text=True)
        if out.returncode != 0:
            return []
        root = ET.fromstring(out.stdout)
        rows = []
        for hsp in root.iter("Hsp"):
            rows.append((int(hsp.findtext("Hsp_score")), hsp.findtext("Hsp_bit-score"), hsp.findtext("Hsp_evalue")))
        return rows


def main():
    lines = ["kind\tmatrix_or_reward\tpenalty\tgap_existence\tgap_extension\tquery_length\tsubject_length\tscore\tbit_score\tevalue"]
    protein = [("BLOSUM62", 11, 1), ("BLOSUM62", 10, 1), ("BLOSUM62", 9, 2), ("BLOSUM62", 12, 1), ("BLOSUM45", 15, 2), ("BLOSUM45", 13, 3), ("BLOSUM50", 13, 2), ("BLOSUM50", 10, 3),
               ("BLOSUM90", 10, 1), ("BLOSUM90", 9, 2), ("PAM30", 9, 1), ("PAM30", 12, 1), ("PAM70", 10, 1), ("PAM70", 11, 1), ("PAM250", 14, 2), ("PAM250", 13, 2), ("BLOSUM80", 10, 1), ("BLOSUM80", 9, 2)]
    for matrix, go, ge in protein:
        for _ in range(14):
            q, s = related_pair(AA, rnd.choice([0.1, 0.25, 0.4]))
            cmd = ["blastp", "-comp_based_stats", "0", "-seg", "no", "-matrix", matrix, "-gapopen", str(go), "-gapextend", str(ge), "-soft_masking", "false", "-word_size", "2", "-threshold", "11"]
            for score, bits, ev in run(cmd, q, s):
                lines.append("\t".join(["protein", matrix, "0", str(go), str(ge), str(len(q)), str(len(s)), str(score), bits, ev]))
    nucl = [(2, -3, 5, 2), (2, -3, 4, 4), (2, -3, 6, 2), (1, -2, 0, 0), (1, -2, 2, 2), (1, -3, 2, 2), (1, -4, 1, 2), (1, -5, 3, 3), (1, -1, 3, 2), (4, -5, 6, 5), (5, -4, 10, 6), (3, -4, 6, 3), (2, -5, 2, 4), (2, -7, 4, 4), (3, -2, 5, 5)]
    for reward, penalty, go, ge in nucl:
        for _ in range(14):
            q, s = related_pair("ACGT", rnd.choice([0.05, 0.12, 0.2]))
            cmd = ["blastn", "-task", "blastn", "-dust", "no", "-soft_masking", "false", "-reward", str(reward), "-penalty", str(penalty), "-gapopen", str(go), "-gapextend", str(ge), "-word_size", "7"]
            for score, bits, ev in run(cmd, q, s):
                lines.append("\t".join(["nucleotide", str(reward), str(penalty), str(go), str(ge), str(len(q)), str(len(s)), str(score), bits, ev]))
    with open(OUT, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as gz:
            gz.write(("\n".join(lines) + "\n").encode())
    print("rows", len(lines) - 1, "->", OUT)


if __name__ == "__main__":
    main()
