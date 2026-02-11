# Evidence Artifact: SPLICE-DONOR-001

**Test Unit ID:** SPLICE-DONOR-001
**Algorithm:** Donor (5') Splice Site Detection
**Date Collected:** 2026-02-12

---

## Online Sources

### Wikipedia: RNA splicing

**URL:** https://en.wikipedia.org/wiki/RNA_splicing
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia with cited primaries)

**Key Extracted Points:**

1. **Donor site dinucleotide:** The splice donor site includes an almost invariant sequence GU (GT in DNA) at the 5' end of the intron, within a larger, less highly conserved region.
2. **Consensus sequence:** G-G-[cut]-G-U-R-A-G-U (donor site), where R = purine (A/G). This is the canonical 5' splice site consensus from Molecular Biology of the Cell (Alberts et al.).
3. **GU-AG rule:** The major spliceosome splices introns containing GU at the 5' splice site and AG at the 3' splice site. This accounts for >99% of splicing.
4. **Non-canonical splicing:** When intronic flanking sequences do not follow the GU-AG rule, noncanonical splicing occurs via the minor spliceosome.
5. **Extended consensus:** MAG|GURAGU, where M = A or C, the vertical bar indicates the exon-intron boundary, positions -3 to +6 relative to the splice junction.

### Wikipedia: Spliceosome

**URL:** https://en.wikipedia.org/wiki/Spliceosome
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia with cited primaries)

**Key Extracted Points:**

1. **U1 snRNP recognition:** U1 snRNP binds to the GU sequence at the 5' splice site of an intron during E complex formation.
2. **Intron elements:** Pre-mRNA introns contain specific sequence elements: 5' end splice site, branch point sequence, polypyrimidine tract, and 3' end splice site.
3. **Invariant dinucleotide:** Introns typically have a GU nucleotide sequence at the 5' end splice site.

### Shapiro & Senapathy (1987): Splice site analysis

**Citation:** Shapiro MB, Senapathy P (1987). RNA splice junctions of different classes of eukaryotes: sequence statistics and functional implications in gene expression. Nucleic Acids Res 15(17):7155-7174.
**Authority rank:** 1 (Peer-reviewed paper — foundational splice site consensus work)

**Key Extracted Points:**

1. **Position weight matrix:** Established nucleotide frequency distributions at positions -3 to +6 relative to the donor splice junction from a compilation of thousands of splice sites.
2. **Consensus:** The extended donor consensus is (C/A)AG|GU(A/G)AGU, with position 0 (G) and position +1 (U) being almost invariant (~100%).
3. **Position -1 (G):** Strongly conserved (~80% G), critical for U1 snRNA base-pairing.
4. **Position -2 (A):** Moderately conserved (~60% A).
5. **Position -3 (M):** A or C (approximately 35% each).
6. **Scoring approach:** Log-odds scoring of observed vs. background (uniform 0.25) frequency at each position is standard for PWM-based splice site prediction.

### Burge, Tuschl & Sharp (1999): Splicing precursors to mRNAs by the spliceosomes

**Citation:** Burge CB, Tuschl T, Sharp PA (1999). In: The RNA World, 2nd ed. Cold Spring Harbor Lab. Press, pp. 525-560.
**Authority rank:** 1 (Peer-reviewed textbook chapter)

**Key Extracted Points:**

1. **GC-AG introns:** A small fraction (~0.5–1%) of U2-type introns use GC instead of GT at the 5' splice site, still with AG at the 3' end.
2. **U12-type introns:** Minor spliceosome introns use AT-AC dinucleotides (AT at 5', AC at 3') instead of GT-AG.
3. **U12 donor consensus:** /ATATCC/ extended motif at the 5' end of U12-type introns.

### Yeo & Burge (2004): MaxEntScan

**Citation:** Yeo G, Burge CB (2004). Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals. J Comput Biol 11(2-3):377-394.
**Authority rank:** 1 (Peer-reviewed paper)

**Key Extracted Points:**

1. **MaxEntScan method:** Uses maximum entropy distributions to model short sequence motifs for splice site scoring.
2. **9-mer donor model:** Scores a 9-nucleotide window (positions -3 to +6) around the donor site.
3. **Score interpretation:** Higher scores indicate stronger splice sites; scores can be compared between candidate sites.
4. **Log-odds scoring:** PWM-based scoring uses log2(observed_frequency / background_frequency) at each position.

---

## Documented Corner Cases and Failure Modes

### From Shapiro & Senapathy (1987)

1. **Non-GT donor sites:** ~1% of introns use GC instead of GT at position 0,+1. These are valid U2-type introns processed by the major spliceosome.
2. **Cryptic splice sites:** Point mutations can activate cryptic splice sites that resemble the consensus but are not normally used.

### From Burge et al. (1999)

1. **U12-type AT donors:** Minor spliceosome introns with AT-AC boundaries are rare (~0.3% of introns) but biologically important.
2. **Context dependency:** Splice site strength depends not only on the dinucleotide but on the extended sequence context.

### From General Splice Site Biology

1. **Short sequences:** Sequences shorter than the PWM window (9 nt for donor) cannot be scored meaningfully.
2. **Empty/null input:** Should return no results, not throw exceptions.
3. **Case insensitivity:** DNA (T) and RNA (U) representations both valid; both uppercase and lowercase.

---

## Test Datasets

### Dataset 1: Canonical Human β-globin Intron 1 Donor Site

**Source:** Human β-globin gene (HBB, NCBI Gene ID: 3043), GenBank J00179

| Parameter | Value |
|-----------|-------|
| Exon-intron junction | ...CAG\|GTTGGT... (DNA) |
| RNA representation | ...CAG\|GUUGGU... |
| Position of GT | After CAG (exon), at intron start |
| Extended context (-3 to +6) | CAGGTTGGTG (DNA) / CAGGUUGGUG (RNA) |
| Expected classification | Strong canonical donor (GT at +1/+2) |

### Dataset 2: Splice Site Scoring Reference

**Source:** Derived from Shapiro & Senapathy (1987) consensus weights

| Motif (9-mer) | Type | Expected Relative Score |
|---------------|------|------------------------|
| CAGGUAAGU | Perfect consensus | Highest (near-optimal) |
| AAGGUAAAU | Weaker (-3 position A, +4 not G) | Moderate |
| UUUGUAAUU | Poor context | Low |
| CAGGCAAGU | GC donor (non-canonical) | Lower than GT, but valid |

### Dataset 3: Negative Controls

**Source:** Trivially correct — no GT/GU dinucleotide present

| Input | Expected Result |
|-------|-----------------|
| AAAAACCCCC | Empty (no GT/GU) |
| "" (empty string) | Empty |
| GTAA (< 6 nt with context) | Empty (too short) |

---

## Assumptions

1. **ASSUMPTION: PWM values in implementation** — The specific PWM values in the implementation are approximations of Shapiro & Senapathy statistics. Exact values not individually verified against the original paper's tables; the pattern of conservation (positions 0,+1 invariant; -1 strongly G; -2 moderately A) is consistent with published consensus.

2. **ASSUMPTION: Score normalization** — The implementation normalizes log-odds scores to [0, 1] range using an approximate linear transformation. The normalization formula is implementation-specific and not derived from MaxEntScan directly.

3. **ASSUMPTION: GC donor penalty** — The 0.7 multiplier for GC donor sites is an implementation heuristic. No specific penalty factor was found in authoritative sources; the general principle that GC donors are weaker than GT donors is well-established.

---

## Recommendations for Test Coverage

1. **MUST Test:** Canonical GT donor detection in well-known sequence — Evidence: Shapiro & Senapathy (1987), Wikipedia RNA splicing
2. **MUST Test:** Perfect consensus CAGGUAAGU scores highest — Evidence: Shapiro & Senapathy (1987) position frequency data
3. **MUST Test:** No GT/GU in sequence returns empty — Evidence: trivially correct
4. **MUST Test:** Strong context scores higher than weak context — Evidence: PWM model properties (Yeo & Burge, 2004)
5. **MUST Test:** Empty/short input returns empty — Evidence: trivially correct
6. **MUST Test:** GC non-canonical donor detected when includeNonCanonical=true — Evidence: Burge et al. (1999)
7. **MUST Test:** DNA T and RNA U input equivalence — Evidence: implementation converts T→U
8. **MUST Test:** Case insensitivity — Evidence: implementation calls ToUpperInvariant
9. **SHOULD Test:** Multiple GT sites in sequence all detected — Rationale: algorithm scans full sequence
10. **SHOULD Test:** Score and confidence are in [0, 1] range — Rationale: normalization invariant
11. **SHOULD Test:** Motif context string returned is non-empty — Rationale: caller needs context for display
12. **COULD Test:** U12-type AT donor detected with includeNonCanonical — Rationale: minor spliceosome support

---

## References

1. Shapiro MB, Senapathy P (1987). RNA splice junctions of different classes of eukaryotes: sequence statistics and functional implications in gene expression. Nucleic Acids Res 15(17):7155-7174.
2. Burge CB, Tuschl T, Sharp PA (1999). Splicing precursors to mRNAs by the spliceosomes. In: The RNA World, 2nd ed. CSHL Press, pp. 525-560.
3. Yeo G, Burge CB (2004). Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals. J Comput Biol 11(2-3):377-394.
4. Wikipedia: RNA splicing. https://en.wikipedia.org/wiki/RNA_splicing (accessed 2026-02-12).
5. Wikipedia: Spliceosome. https://en.wikipedia.org/wiki/Spliceosome (accessed 2026-02-12).
6. Alberts B et al. Molecular Biology of the Cell. Garland Science.

---

## Change History

- **2026-02-12**: Initial documentation.
