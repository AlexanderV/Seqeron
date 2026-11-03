#!/usr/bin/env python3
"""Optional extra input dna_100m (used only to find where the Seqeron persistent build stops).
Same stream/seed as gen_inputs.py: first 50 Mbp are byte-identical to dna_50m (asserted).  usage: gen_extra_100m.py <outdir>"""
import importlib.util
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("gi", os.path.join(here, "gen_inputs.py"))
gi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gi)
out = sys.argv[1]
d = gi.rand_dna(gi.rng("dna50"), 100_000_000)
with open(f"{out}/data/dna_50m.txt", "rb") as f:
    assert f.read() == d[:50_000_000], "prefix property violated"
gi.write(f"{out}/data/dna_100m.txt", d)
gi.make_queries(out, "dna_100m", d, gi.DNA, False)
with open(f"{out}/MANIFEST.tsv", "a") as f:
    f.write("dna_100m\t100000000\tdna4\n")
print("wrote dna_100m")
