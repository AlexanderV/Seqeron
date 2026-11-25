# Evidence Artifact: ASSEMBLY-COVER-001

**Test Unit ID:** ASSEMBLY-COVER-001
**Algorithm:** Coverage (Depth) Calculation — per-base sequencing depth over a reference
**Date Collected:** 2026-06-13

---

## Online Sources

### Illumina — "Sequencing Coverage for NGS Experiments"

**URL:** https://sapac.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 2 (official vendor specification / standard usage)

**Key Extracted Points:**

1. **Coverage definition:** Retrieved text: "Next-generation sequencing (NGS) coverage describes the average number of reads that align to, or 'cover,' known reference bases." Coverage / sequencing depth = average number of reads covering each reference base.
2. **Average-coverage formula:** The page presents the Lander-Waterman formula **C = LN / G**, where **C** = coverage, **L** = read length, **N** = number of reads, **G** = haploid genome length. It is "a method for computing genome coverage."

### Daniel E. Cook — "Calculate Depth and Breadth of Coverage From a bam File"

**URL:** https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (reference-implementation / established workflow over samtools output)

**Key Extracted Points:**

1. **Depth of coverage:** Retrieved text states depth and coverage "refer to the same thing, the average number of reads aligned to an individual base." Computed as **Sum of Depths / genome size** (or contig length per chromosome).
2. **Per-base depth:** measured via `samtools depth` as "the number of reads aligned over a given base."
3. **Breadth of coverage:** "the percentage of the genome with aligned reads to it," computed as **Bases Mapped / genome size**, a fraction representing the proportion of bases with at least one aligned read.

### Metagenomics Wiki — "SAMtools: get breadth of coverage"

**URL:** https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (established bioinformatics tooling documentation)

**Key Extracted Points:**

1. **Per-base depth:** "the number of reads mapping to a specific reference position."
2. **Breadth of coverage:** "percentage/fraction of reference bases covered by at least one read" (example: 32,876 covered bases / 45,678 total = 0.719×).
3. **Average depth:** "total bases mapped / reference length" = "sum of per-base depths / genome size."

### Predicting the Number of Bases to Attain Sufficient Coverage (Daley et al.) — PMC7398442

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7398442/
**Accessed:** 2026-06-13 (fetched via WebFetch, after 301 redirect from ncbi.nlm.nih.gov)
**Authority rank:** 1 (peer-reviewed, restates the Lander-Waterman model)

**Key Extracted Points:**

1. **Lander-Waterman assumption:** Retrieved text: "the basic statistical assumption of the Lander/Waterman model is that reads are generated uniformly at random from the genome. Under this assumption, for each site in the genome, the number of reads that cover this site can be approximated by a Poisson distribution, where the Poisson rate λ is estimated by average sequencing depth."
2. **Sufficient coverage:** "a site as having sufficient coverage if the site is covered by at least r reads."

### Lander & Waterman (1988) Poisson model (via search-retrieved restatement)

**URL:** Search results for "Lander Waterman 1988 Genomics genomic mapping by fingerprinting random clones Poisson coverage" and "probability that a base is not covered e^-c gaps" (WebSearch; primary: Genomics 3:231-239).
**Accessed:** 2026-06-13 (WebSearch summary; primary not directly fetchable as full text in this session)
**Authority rank:** 1 (primary paper) / 4 (the retrieved restatement)

**Key Extracted Points:**

1. **Gap probability:** Retrieved text: "The probability that any given base in a sequence is not sequenced can be calculated by the equation **P = e^−m**, where m is the fold coverage." Worked numbers: 1× → e^−1 = 0.37 (≈37% uncovered); 5× → e^−5 = 0.0067.
2. **Breadth complement:** "The complementary probability (**1 − e^−c**) represents the probability that a base is covered at least once."

---

## Documented Corner Cases and Failure Modes

### From samtools / Metagenomics Wiki / Daniel Cook

1. **Position with zero overlapping reads:** depth at that reference position is 0 (it does not contribute to breadth-covered bases).
2. **Read extending past reference end:** only the overlapping portion of a read contributes to per-base depth (clipping at the reference boundary).
3. **Uncovered genome → breadth 0, average depth 0:** when no read aligns, every per-base depth is 0; average depth = 0 and breadth = 0.

### From Lander-Waterman / PMC7398442

1. **Highly variable coverage:** the uniform-Poisson expectation (mean ≈ C) holds only under uniform random read placement; observed per-base depth can deviate. This is a modeling caveat, not an arithmetic rule — the per-base depth array itself is exact regardless.

---

## Test Datasets

### Dataset: Hand-constructed exact-placement example (per-base depth)

**Source:** Definition from samtools / Metagenomics Wiki (per-base depth = number of reads covering each position); deterministic by construction.

Reference `ACGTTGCAAT` (length 10). Reads are distinct 5-mer substrings, each matching uniquely, so placement is unambiguous; per-base depth = count of reads spanning each position.

| Parameter | Value |
|-----------|-------|
| Reference | `ACGTTGCAAT` (len 10) |
| Read 1 | `ACGTT` aligns at pos 0 → covers [0,5) |
| Read 2 | `TTGCA` aligns at pos 3 → covers [3,8) |
| Read 3 | `GCAAT` aligns at pos 5 → covers [5,10) |
| Expected depth array | `[1,1,1,2,2,2,2,2,1,1]` |
| Total bases mapped (Σ depth) | 15 |
| Average depth (Σ/G) | 15/10 = 1.5 |
| Breadth (covered ≥1× / G) | 10/10 = 1.0 |

### Dataset: Lander-Waterman gap/breadth sanity (Poisson)

**Source:** Lander & Waterman (1988); P(uncovered) = e^−c; breadth = 1 − e^−c.

| Parameter | Value |
|-----------|-------|
| c = 1× | P(uncovered) = e^−1 ≈ 0.3679; breadth ≈ 0.6321 |
| c = 5× | P(uncovered) = e^−5 ≈ 0.006738 |

(Used as a property/derivation check, not asserted against the deterministic per-base depth array which is exact.)

---

## Assumptions

1. **ASSUMPTION: Mapping model for read placement** — The unit signature is `CalculateCoverage(reference, reads, minOverlap)`; the sources define depth given an *alignment*, but do not prescribe how an aligner places a read. The repository uses an ungapped best-match scan requiring ≥ `minOverlap` matching characters (`FindBestAlignment`). This placement rule is implementation-level (it determines *where* a read maps, not the depth-counting arithmetic). The depth-counting arithmetic (per-base = number of placed reads spanning the position; clip at reference end) is fully source-defined. Tests therefore use exact-match reads where placement is unambiguous, isolating the source-defined counting rule.

---

## Recommendations for Test Coverage

1. **MUST Test:** Per-base depth array equals the count of reads spanning each position for unambiguous exact-match placements. — Evidence: samtools / Metagenomics Wiki ("number of reads mapping to a specific reference position").
2. **MUST Test:** ~~Read extending past reference end contributes only its overlapping portion (boundary clipping).~~ Superseded 2026-11-25: placement is end-to-end (Bowtie 2 default `--end-to-end`, no soft clipping), so a read never overhangs; a read longer than the reference does not align (M4).
3. **MUST Test:** Unmatched read (below `minOverlap`) contributes 0 to all positions. — Evidence: a read that does not align does not add depth (samtools depth semantics).
4. **MUST Test:** Empty reads list → all-zero depth array of reference length. — Evidence: no aligned reads → every per-base depth 0 (Daniel Cook / Metagenomics Wiki).
5. **SHOULD Test:** Average depth = Σ(per-base depth) / reference length on the worked dataset (1.5). — Rationale: Daniel Cook "Sum of Depths / genome size."
6. **SHOULD Test:** Breadth = (# positions with depth ≥ 1) / reference length on the worked dataset (1.0). — Rationale: Metagenomics Wiki breadth definition.
7. **COULD Test:** Case-insensitive matching (lowercase read vs uppercase reference still maps). — Rationale: implementation normalizes case; verifies depth is placement-driven, not case-driven.

---

## 2026-11-25 review (B10, Stage A/B) — placement semantics

### Sources opened (fetched with a neutral User-Agent)

| Source | URL | Points used |
|---|---|---|
| Bowtie 2 MANUAL | https://raw.githubusercontent.com/BenLangmead/bowtie2/master/MANUAL.markdown | `--end-to-end` "requires that the entire read align from one end to the other, without any trimming (or "soft clipping")… `--end-to-end` is the default mode"; `--nofw/--norc` "Default: both strands enabled"; "Default mode: search for multiple alignments, report the best one"; equally good alignments → pseudo-random choice; MAPQ = uniqueness |
| Bowtie 2 `unique.h` | https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h | `BowtieMapq2`, end-to-end: `bestdiff == 0` (best AS = second-best AS) → MAPQ 1 or 0; `bestdiff > 0` → MAPQ ≥ 2 |
| samtools-depth(1) | https://raw.githubusercontent.com/samtools/samtools/develop/doc/samtools-depth.1 | `-a` all positions; `-Q` min MAPQ; default excludes UNMAP,SECONDARY,QCFAIL,DUP (one primary record per read counted); deletions ignored unless `-J` |
| bwa(1) | https://raw.githubusercontent.com/lh3/bwa/master/bwa.1 | samse: "Repetitive hits will be randomly chosen" (one placement per read) |

### Findings

1. **Reverse-strand reads never placed.** `CalculateCoverage("ACGGTCAATGCC", ["TTGACC"], 6)` → all 0, but `TTGACC` = revcomp(`GGTCAA`) at [2,8); Bowtie 2/BWA search both strands by default → samtools depth `[0,0,1,1,1,1,1,1,0,0,0,0]`. Fix: overload `CalculateCoverage(reference, reads, minOverlap, includeReverseComplement, excludeMultiMapped = false)`; 3-arg overload unchanged (forward only, backward compatible, as ASSEMBLY-DBG-001).
2. **No multi-mapping handling.** Equal-best reads were silently counted at the leftmost locus with no way to apply the standard MAPQ filter. Default unchanged (Bowtie 2 reports one alignment); opt-in `excludeMultiMapped` drops reads whose best Hamming distance occurs at > 1 distinct start (both strands pooled) = `samtools depth -Q 2` over Bowtie 2 MAPQ. `ACGTACGTAC`, [`ACGTA`,`GTACG`], 5 → `[0,0,1,1,1,1,1,0,0,0]`.
3. **minOverlap ≤ 0 placed reads with zero matching bases.** `("AAAA", ["GGGG"], 0)` → `[1,1,1,1]`. Now `ArgumentOutOfRangeException` for minOverlap < 1 (ASSEMBLY-OLC-001 F4 precedent); null read element → `ArgumentException` (was `NullReferenceException`).
4. **Doc/code mismatch (clipping).** XML doc and this Evidence claimed overhanging reads are clipped; the scan never placed an overhanging read (clip branch dead). Doc corrected to Bowtie 2 end-to-end semantics; dead branch removed.
5. **Duplication.** Private `FindBestAlignment` (char-loop max-matches scan) = `ApproximateMatcher.FindBestMatch` (leftmost min Hamming, case-insensitive) → removed; placement now delegates to `FindBestMatch`/`FindWithMismatches`; RC via Core `TryGetReverseComplement`. MCP `calculate_coverage` delegates (new optional `includeReverseComplement`, `excludeMultiMapped`).

### Numeric oracle

`pysam` 0.24.1 (pip) `samtools depth -a` (with `-Q 2` for `excludeMultiMapped`) over SAM records written by an independent Python placement oracle (flag 16 + reverse-complemented SEQ for reverse-strand hits, MAPQ 0 for equal-best multi-locus reads else 42, CIGAR `<L>M`). 14 hand cases + 2000 random cases (seed 2026; reads = mutated/reverse-complemented/lower-cased reference windows or random, minOverlap 1–6, both flags random): C# result equal on **2014/2014**, and the 3-arg overload equal on every forward-only, non-excluding case. Hand cases (locked in unit tests):

| reference | reads | minOverlap | RC | excl | samtools depth -a |
|---|---|---|---|---|---|
| ACGTACGTAC | ACGTA, GTACG | 5 | F | F | 1,1,2,2,2,1,1,0,0,0 |
| ACGTACGTAC | ACGTA, GTACG | 5 | F | T | 0,0,1,1,1,1,1,0,0,0 |
| ACGGTCAATGCC | TTGACC | 6 | F | F | 0×12 |
| ACGGTCAATGCC | TTGACC | 6 | T | F | 0,0,1,1,1,1,1,1,0,0,0,0 |
| ACGGTCAATGCC | ttgacc, GGTCAA | 6 | T | F | 0,0,2,2,2,2,2,2,0,0,0,0 |
| TTACGTTT | ACGT | 4 | T | T | 0,0,1,1,1,1,0,0 |
| AAAAAAAAAA | AAAAA | 5 | F | T | 0×10 |
| AAAAAAAAAA | AAAAA | 5 | F | F | 1,1,1,1,1,0,0,0,0,0 |
| GATTACAGATTACA | GATTAC | 6 | T | T | 0×14 |
| CCCCAGGTCTTTT | AGCTCT | 5 | F | F | 0,0,0,0,1,1,1,1,1,1,0,0,0 |
| CCCCAGGTCTTTT | AGCTCT | 6 | F | F | 0×13 |
| ACGTTGCAAT | ACGTT, TTGCA, GCAAT | 5 | T | F | 1,1,1,2,2,2,2,2,1,1 |
| GGGAAACCC | TTTCC | 4 | T | T | 0,1,1,1,1,1,0,0,0 |

Oracle code (`python3 -I cover_oracle.py cases.json 2000`):

```python
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
```

---

## References

1. Illumina, Inc. (n.d.). Sequencing Coverage for NGS Experiments. https://sapac.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html
2. Cook, D.E. (n.d.). Calculate Depth and Breadth of Coverage From a bam File. https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/
3. Metagenomics Wiki. (n.d.). SAMtools: get breadth of coverage. https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
4. Daley, T. et al. (2020). Predicting the Number of Bases to Attain Sufficient Coverage in High-Throughput Sequencing Experiments. PMC7398442. https://pmc.ncbi.nlm.nih.gov/articles/PMC7398442/
5. Lander, E.S., Waterman, M.S. (1988). Genomic mapping by fingerprinting random clones: a mathematical analysis. Genomics 3(2):231-239. https://doi.org/10.1016/0888-7543(88)90007-9 (primary; full text not fetched this session — Poisson gap/breadth formulas taken from retrieved restatements in refs 1 & 4 and WebSearch summaries).

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-11-25**: B10 review — Bowtie 2 / samtools / BWA sources, placement findings, samtools-depth oracle (2014/2014).
