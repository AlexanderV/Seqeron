#!/usr/bin/env python3
"""Write runs/na_rows.tsv: explicit NA rows (capability not provided by the tool's own API) so that gaps in bench.tsv are visible
rather than silently missing.  Reads the raw files to learn (tool, version, input, n).  usage: make_na.py <bench-dir>"""
import glob
import sys

B = sys.argv[1]
OPS = [("lrs", "first_call_time", "s"), ("mem_q100k_min20", "time", "s"), ("mum_q100k_min20", "time", "s"), ("anchors_q100k_min20", "time", "s")]
REASON = {
    "go-suffixarray": "index/suffixarray offers Lookup/FindAllIndex/Write/Read only: no LRS, MEM, MUM or anchor API",
    "libsais": "SA construction library (libsais_lcp builds LCP) - no LRS/MEM/MUM/anchor API; would be user code on SA+LCP",
    "libdivsufsort": "SA construction library - no LRS/MEM/MUM/anchor API; would be user code on SA+LCP",
    "rust-bio": "bio 2.3 provides suffix array + FM-index; no LRS/MEM/MUM/anchor API (user code needed)",
    "pydivsufsort": "provides divsufsort/kasai/sa_search helpers; no LRS/MEM/MUM/anchor API (user code needed)",
    "sdsl-csa_wt_huff": "csa_* classes are indexes without LRS/MEM/MUM/anchor API",
    "sdsl-csa_wt_rrr": "csa_* classes are indexes without LRS/MEM/MUM/anchor API",
    "sdsl-csa_sada": "csa_* classes are indexes without LRS/MEM/MUM/anchor API",
    "sdsl-cst_sct3": "no MEM/MUM/anchor API (cst traversal would be user code); LRS measured only where an lrs row exists",
    "suffix-trees-py": "package offers find/find_all/lcs only: no LRS, MEM, MUM or anchor API",
}
seen = {}
for f in sorted(glob.glob(f"{B}/runs/raw*.tsv")):
    if f.endswith("raw_all.tsv"):
        continue
    with open(f) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10 or p[0] not in REASON:
                continue
            seen[(p[0], p[1], p[2], p[3])] = 1
have = set()
for f in sorted(glob.glob(f"{B}/runs/raw*.tsv")):
    if f.endswith("raw_all.tsv"):
        continue
    with open(f) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 10:
                have.add((p[0], p[2], p[4]))
with open(f"{B}/runs/na_rows.tsv", "w") as o:
    o.write("tool\tversion\tinput\tn\toperation\tmetric\tvalue\tunit\trepeats\tnotes\n")
    for (tool, ver, inp, n) in sorted(seen):
        if inp not in ("dna_1m", "dna_10m", "dna_1m_prefix_1000000"):
            continue
        for op, metric, unit in OPS:
            if (tool, inp, op) in have:
                continue
            o.write("\t".join([tool, ver, inp, n, op, metric, "NA", unit, "0", "NA: " + REASON[tool]]) + "\n")
print("na rows written")
