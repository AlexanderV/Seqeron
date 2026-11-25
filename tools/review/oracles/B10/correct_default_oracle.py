"""R1-WP3 A9-1 oracle: musket_oracle.py two-sided stage to convergence (max_iters huge, num_errors=0) with a canonical vs a
forward-only spectrum, on the 8-read k=5, cut-off 2 rows of SequenceAssembler_ErrorCorrectReads_Tests.
Run: python3 -I correct_default_oracle.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import musket_oracle as m

ROWS = ["GCAAGGCCGGCTCCT,GGTCCCCCAAGGCCG,CATGGCCGGCTCCAC,TGGGGGACCGACGGA,CTTGGAGGACCGACG,TCGGTCCCCCAAGGC,CCAAGGCCGGCTCCA,GCTCTGGGGGACCGA",
        "ACAAAGCACCTCAAC,TACCAATCACTTCTA,TTGCGTTACCAACCA,GTGCTTGGGAACGCA,ACGTTACCAAGCACC,TGCTTGTTACCGCAA,CACCTCTACGGGTAT,CGTTACCAAGCACCT"]
if __name__ == "__main__":
    for r in ROWS:
        rd = r.split(',')
        for name, canon in (("CANONICAL (3-arg default)", True), ("FORWARD (TwoSidedForwardStrand)", False)):
            print(name, ",".join(m.correct_all(rd, 5, 2, canonical=canon, max_iters=10**9, num_errors=0)))
    rv = ["GAAGCCTGTAATC"] * 3 + ["GATTACGGGCTTC"]
    for name, canon in (("CANONICAL", True), ("FORWARD", False)):
        print("F32", name, ",".join(m.correct_all(rv, 3, 2, canonical=canon, max_iters=10**9, num_errors=0)))
