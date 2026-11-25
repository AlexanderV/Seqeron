# B10 R1-WP1 (A15-2): Biopython gc_fraction over range(0, len(seq), step) slices seq[i:i+w] (GC_skew windowing),
# for the three ambiguity modes. Prints C# expected arrays (percent, i.e. fraction * 100).
from Bio.SeqUtils import gc_fraction
seq = "ACGTNSWRYacgtnGGCC"
w, step = 5, 4
for mode in ("remove", "ignore", "weighted"):
    vals = [gc_fraction(seq[i:i + w], ambiguous=mode) * 100 for i in range(0, len(seq), step)]
    print(mode, [repr(v) for v in vals])
