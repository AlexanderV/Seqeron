#!/usr/bin/env python3
"""B09 A5-10b3 oracle: iteration order of robin_hood::unordered_map (the robin_hood.h 3.9.1 bundled with Mash 2.3).

Downloads marbl/Mash v2.3 src/mash/robin_hood.h (raw.githubusercontent.com, User-Agent Seqeron-review/1.0), compiles
scripts/alignment/robin_hood_order_harness.cpp against it and runs 150 deterministic operation streams (splitmix64):
mode 0 = MinHashHeap::tryInsert on a bottom-s set, mode 1 = random inserts / erases.  Output: robin_hood_order_oracle.tsv
(bits mode seed n param size fnv1a-of-(key,value)-in-iteration-order).
"""
import os, random, subprocess, sys, tempfile, urllib.request
here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "..", "..", "tests", "Seqeron", "Seqeron.Genomics.Tests", "TestData", "Alignment", "robin_hood_order_oracle.tsv")
tmp = tempfile.mkdtemp()
req = urllib.request.Request("https://raw.githubusercontent.com/marbl/Mash/v2.3/src/mash/robin_hood.h", headers={"User-Agent": "Seqeron-review/1.0"})
open(os.path.join(tmp, "robin_hood.h"), "wb").write(urllib.request.urlopen(req).read())
exe = os.path.join(tmp, "harness")
subprocess.check_call(["g++", "-O2", "-std=c++14", "-I", tmp, "-o", exe, os.path.join(here, "robin_hood_order_harness.cpp")])
rng = random.Random(20261121)
cases = []
for i in range(150):
    bits = rng.choice([32, 64])
    mode = 0 if i % 3 != 2 else 1
    n = rng.choice([10, 100, 1000, 5000, 20000, 60000])
    param = rng.choice([1, 2, 7, 8, 9, 50, 300, 1000, 2000, 5000, 12000])
    cases.append((bits, mode, rng.getrandbits(62), n, param))
inp = "".join("%d %d %d %d %d\n" % c for c in cases)
res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
with open(out, "w") as f:
    for c, r in zip(cases, res):
        f.write("\t".join(map(str, c)) + "\t" + r.replace(" ", "\t") + "\n")
print(len(cases), "cases ->", out)
