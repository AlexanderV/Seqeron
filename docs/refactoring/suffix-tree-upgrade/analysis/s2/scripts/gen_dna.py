#!/usr/bin/env python3
"""The benchmark DNA stream of analysis/bench/gen/gen_inputs.py (random ACGT, one symbol per byte) at any length: dna_10m, dna_50m and dna_100m are its prefixes (the generator's "prefix property").
usage: gen_dna.py <n> <out-file>      e.g. gen_dna.py 500000000 $S2_OUT/data/dna_500m.txt   (about 20 s and 1.5 GB of RAM at 5x10^8)"""
import importlib.util, os, sys
here = os.path.dirname(os.path.abspath(__file__))
gi_path = os.path.join(here, '../../bench/gen/gen_inputs.py')
spec = importlib.util.spec_from_file_location('gi', gi_path)
gi = importlib.util.module_from_spec(spec); spec.loader.exec_module(gi)
n, out = int(sys.argv[1]), sys.argv[2]
# Python's getrandbits cannot make more than 2^31 bits in one call; chunks of a multiple of 4 bytes consume whole 32-bit words of the generator, so the concatenation is the stream of one big call
r = gi.rng('dna50')
CH = 100_000_000
with open(out, 'wb') as f:
    left = n
    while left > 0:
        m = min(CH, left)
        f.write(r.getrandbits(8 * m).to_bytes(m, 'little').translate(gi.DNA_TABLE))
        left -= m
print(f'wrote {out}: {n} symbols')
