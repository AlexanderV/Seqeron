"""Independent oracle for ASSEMBLY-COVER-001.
Placement: Bowtie2-style end-to-end ungapped (no soft clip), min Hamming, leftmost, forward kept on ties,
optional --norc off (both strands); MAPQ 0 when best score attained at >1 distinct locus (bowtie2 unique.h
BowtieMapq2: bestdiff==0 -> MAPQ<=1). Depth: real samtools depth (pysam) -a, with -Q 2 for excludeMultiMapped.
"""
import json, os, random, sys, tempfile
import pysam

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(s): return s.translate(COMP)[::-1]
def ham(a, b): return sum(x != y for x, y in zip(a, b))

def hits(ref, q):
    return [(ham(ref[p:p+len(q)], q), p) for p in range(len(ref) - len(q) + 1)]

def place(ref, read, minov, both):
    R, Q = ref.upper(), read.upper()
    if not R or not Q or len(Q) > len(R): return None
    f = hits(R, Q); bf = min(f)            # (dist, pos) -> leftmost among min
    best, strand = bf, 0
    r = hits(R, rc(Q)) if both else []
    if r:
        br = min(r)
        if br[0] < bf[0]: best, strand = br, 16
    if len(Q) - best[0] < minov: return None
    loci = {p for d, p in f if d == best[0]} | {p for d, p in r if d == best[0]}
    return best[1], strand, len(loci) > 1

def samtools_depth(ref, reads, minov, both, excl):
    d = tempfile.mkdtemp()
    sam = os.path.join(d, "a.sam")
    with open(sam, "w") as fh:
        fh.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:{len(ref)}\n")
        for i, read in enumerate(reads):
            pl = place(ref, read, minov, both)
            if pl is None: continue
            pos, flag, multi = pl
            seq = read.upper() if flag == 0 else rc(read.upper())
            fh.write(f"r{i}\t{flag}\tref\t{pos+1}\t{0 if multi else 42}\t{len(read)}M\t*\t0\t0\t{seq}\t*\n")
    bam = os.path.join(d, "a.bam")
    pysam.sort("-o", bam, sam); pysam.index(bam)
    args = ["-a"] + (["-Q", "2"] if excl else []) + [bam]
    out = pysam.depth(*args)
    depth = [0] * len(ref)
    for line in out.splitlines():
        _, p, c = line.split("\t"); depth[int(p) - 1] = int(c)
    return depth

if __name__ == "__main__":
    cases = []
    hand = [
        ("ACGTACGTAC", ["ACGTA", "GTACG"], 5, False, False),
        ("ACGTACGTAC", ["ACGTA", "GTACG"], 5, False, True),
        ("ACGGTCAATGCC", ["TTGACC"], 6, False, False),
        ("ACGGTCAATGCC", ["TTGACC"], 6, True, False),
        ("ACGGTCAATGCC", ["ttgacc", "GGTCAA"], 6, True, False),
        ("TTACGTTT", ["ACGT"], 4, True, True),
        ("AAAAAAAAAA", ["AAAAA"], 5, False, True),
        ("AAAAAAAAAA", ["AAAAA"], 5, False, False),
        ("GATTACAGATTACA", ["GATTAC"], 6, True, True),
        ("CCCCAGGTCTTTT", ["AGGTCT"], 5, True, False),
        ("CCCCAGGTCTTTT", ["AGCTCT"], 5, False, False),
        ("CCCCAGGTCTTTT", ["AGCTCT"], 6, False, False),
        ("ACGTTGCAAT", ["ACGTT", "TTGCA", "GCAAT"], 5, True, False),
        ("GGGAAACCC", ["TTTCC"], 4, True, True),
    ]
    for ref, reads, mo, b, e in hand:
        cases.append(dict(reference=ref, reads=reads, minOverlap=mo, rc=b, excl=e,
                          depth=samtools_depth(ref, reads, mo, b, e)))
    rng = random.Random(2026)
    for t in range(int(sys.argv[2]) if len(sys.argv) > 2 else 400):
        g = rng.randint(1, 24); ref = "".join(rng.choice("ACGT") for _ in range(g))
        reads = []
        for _ in range(rng.randint(0, 5)):
            L = rng.randint(1, g + 2)
            if rng.random() < 0.6 and L <= g:
                s = rng.randint(0, g - L); rd = ref[s:s+L]
                if rng.random() < 0.4: rd = rc(rd)
                rd = "".join(c if rng.random() > 0.1 else rng.choice("ACGT") for c in rd)
            else:
                rd = "".join(rng.choice("ACGT") for _ in range(L))
            if rng.random() < 0.2: rd = rd.lower()
            reads.append(rd)
        mo = rng.randint(1, 6); b = rng.random() < 0.5; e = rng.random() < 0.5
        cases.append(dict(reference=ref, reads=reads, minOverlap=mo, rc=b, excl=e,
                          depth=samtools_depth(ref, reads, mo, b, e)))
    json.dump(cases, open(sys.argv[1], "w"))
    for c in cases[:len(hand)]: print(c["reference"], c["reads"], c["minOverlap"], c["rc"], c["excl"], c["depth"])
