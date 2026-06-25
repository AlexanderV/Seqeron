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
4. **Log-odds scoring:** the maximum-entropy score is `log2( P_maxent(seq) / P_background(seq) )`.

**Retrieved this session (2026-06-25):** the reference `score5` factorisation and the precomputed
probability table were fetched VERBATIM from the MIT-licensed maxentpy port (kepbod/maxentpy):
- factorisation: `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/maxent.py` (WebFetch, `score5`)
- table: `https://raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/data/score5_matrix.txt` (curl, 16 384 records)
- licence: `https://github.com/kepbod/maxentpy/blob/master/LICENSE` (MIT; full text recorded in `Data/maxent_score3.LICENSE.md` and `Data/maxent_score5.LICENSE.md`)

5. **Factorisation (verbatim from maxentpy `score5`):** the 9-nt donor window (3 exon + 6 intron)
   has its conserved `GT` dinucleotide at 0-based positions 3..4 scored separately and removed; the
   remaining 7-nt "rest" sequence is `window[0:3] + window[5:9]`. The maximum-entropy probability of
   that rest sequence is looked up DIRECTLY in a single table keyed by the 7-mer string (4^7 = 16384
   entries) — unlike score3, score5 is single-matrix with no overlapping sub-windows. The score is:
   - GT term: `cons1_5[G] * cons2_5[T] / (bgd_5[G] * bgd_5[T])`
   - rest term: `matrix[rest]`
   - final score = `log2(GT_term * rest_term)`.
   - Probabilities: `bgd_5 = {A:0.27,C:0.23,G:0.23,T:0.27}`; `cons1_5 = {A:0.004,C:0.0032,G:0.9896,T:0.0032}`;
     `cons2_5 = {A:0.0034,C:0.0039,G:0.0042,T:0.9884}`.

6. **Worked examples (from the maxentpy `score5` docstring, the documented reference values):**
   - `score5('cagGTAAGT')` → **10.86** (reproduced 10.858313 this session — the canonical example)
   - `score5('gagGTAAGT')` → **11.08** (reproduced 11.078494)
   - `score5('taaATAAGT')` → **-0.12** (reproduced -0.116791; a non-GT donor)
   The canonical `cagGTAAGT → 10.86` value is the primary cross-check for the factorisation + table.

7. **PROVENANCE + LICENCE (flagged):** the embedded table (`Data/maxent_score5.txt`) and the
   factorisation come from **kepbod/maxentpy, which is MIT-licensed** (redistribution permitted — full
   MIT text recorded in `Data/maxent_score5.LICENSE.md`). The *original* Burge-lab Perl scripts/models
   carry academic terms (`http://genes.mit.edu/burgelab/maxent/download/READTHIS`); the artifact bundled
   here is the MIT-licensed port table, not the original.

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

### Dataset 4: MaxEntScan score5ss worked examples — Yeo & Burge (2004) / maxentpy

**Source:** the maxentpy `score5` docstring (documented reference values); reproduced this session.

| 9-nt window | Expected (2 dp) | Reproduced (full precision) |
|-------------|-----------------|------------------------------|
| `cagGTAAGT` | **10.86** | 10.858313 |
| `gagGTAAGT` | **11.08** | 11.078494 |
| `taaATAAGT` | **-0.12** | -0.116791 |

The 10.86 value is the canonical example from the MaxEntScan website / maxentpy README; it is the
primary cross-check for the score5 factorisation + embedded table (`ScoreDonorMaxEnt`).

---

## Assumptions

All previous assumptions have been eliminated:
- ~~A1 (PWM values)~~: Replaced with IUPAC consensus binary weights (1.0 = match, 0.0 = no match) derived from the universally verified MAG|GURAGU consensus.
- ~~A2 (Score normalization)~~: Replaced with simple consensus match fraction (matches / positions scored). No ad-hoc formula.
- ~~A3 (GC donor 0.7 penalty)~~: Removed. GC donors naturally score lower because position +1 (C) mismatches the invariant U consensus (max 8/9 ≈ 0.889 vs 9/9 = 1.0 for GT).

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
13. **MUST Test:** `ScoreDonorMaxEnt("cagGTAAGT")` == 10.86 bits — Evidence: maxentpy `score5` docstring (Dataset 4)
14. **MUST Test:** `ScoreDonorMaxEnt` second/third reference values (11.08 / -0.12) — Evidence: Dataset 4
15. **MUST Test:** `ScoreDonorMaxEnt` strong site ranks above weak site — Evidence: Dataset 4 ordering
16. **MUST Test:** `ScoreDonorMaxEnt` rejects wrong length / non-A/C/G/T(/U) / null — Evidence: implementation guard

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
- **2026-06-25 (MaxEntScan score5)**: IMPLEMENTED the Yeo & Burge (2004) MaxEntScan score5ss
  maximum-entropy 5' donor model as the opt-in `ScoreDonorMaxEnt`. Embedded the precomputed
  probability table (`Data/maxent_score5.txt`, 16 384 records) retrieved verbatim this session from
  the MIT-licensed maxentpy port; recorded the `score5` factorisation + provenance + MIT licence
  (Yeo & Burge source #6, dataset 4, `Data/maxent_score5.LICENSE.md`). Cross-checked the canonical
  `score5('cagGTAAGT') == 10.86` plus two further reference values exactly (11.08 / -0.12). The
  existing PWM/consensus donor scorers (`FindDonorSites`, `ScoreDonorSite`) and their defaults are
  unchanged; `ScoreDonorMaxEnt` is a new, additive, opt-in method.
