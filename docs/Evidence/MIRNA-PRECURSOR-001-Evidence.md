# Evidence Artifact: MIRNA-PRECURSOR-001

**Test Unit ID:** MIRNA-PRECURSOR-001
**Algorithm:** Pre-miRNA Hairpin Detection
**Date Collected:** 2026-02-10 (heuristic); 2026-06-24 (MFE-structure-based opt-in)

---

## Online Sources

### Wikipedia: MicroRNA

**URL:** https://en.wikipedia.org/wiki/MicroRNA
**Accessed:** 2026-02-10
**Authority rank:** 4 (Wikipedia with cited primaries)

**Key Extracted Points:**

1. **Pre-miRNA structure:** Hairpin loop structures composed of about 70 nucleotides each, produced from pri-miRNA by Drosha/DGCR8 cleavage.
2. **Mature miRNA length:** 21–23 nucleotides, single-stranded non-coding RNA.
3. **Biogenesis pathway:** pri-miRNA → Drosha/DGCR8 → pre-miRNA (~70 nt hairpin) → Exportin-5 → Dicer → mature miRNA:miRNA* duplex (~22 nt).
4. **3' overhang:** Drosha produces a 2-nucleotide overhang at 3' end of pre-miRNA.
5. **Dicer processing:** Cuts away terminal loop, yielding imperfect miRNA:miRNA* duplex. Loop size and hairpin length influence Dicer processing efficiency.
6. **Strand selection:** One strand (guide) is incorporated into RISC; the other (passenger/star) is degraded. Selection based on thermodynamic instability at 5' end.
7. **Nomenclature:** 5p and 3p arms denoted by position on hairpin.

### Bartel (2004): MicroRNAs: Genomics, Biogenesis, Mechanism, and Function

**URL:** doi:10.1016/S0092-8674(04)00045-5
**Accessed:** 2026-02-10
**Authority rank:** 1 (Peer-reviewed review, Cell 116:281-297)

**Key Extracted Points:**

1. **Pre-miRNA characteristics:** ~70 nt stem-loop hairpin intermediate.
2. **Stem structure:** Approximately 33 bp stem (imperfect, with mismatches and bulges).
3. **Terminal loop:** Variable size, typically 4–15 nt.
4. **Processing:** Drosha cleaves ~11 nt from hairpin base (one helical turn into stem).

### Ambros et al. (2003): A Uniform System for microRNA Annotation

**URL:** doi:10.1261/rna.2183803; PMC 1370393
**Accessed:** 2026-02-10
**Authority rank:** 1 (Peer-reviewed, RNA 9:277-279)

**Key Extracted Points:**

1. **Annotation criteria for miRNA:** (a) expression evidence (Northern blot, cloning); (b) biogenesis evidence (Dicer-dependent); (c) phylogenetic conservation.
2. **Structure criterion:** Must be predicted to fold into a hairpin precursor structure with the ~22-nt mature miRNA residing in one arm of the hairpin.

### Bartel (2009): MicroRNAs: Target Recognition and Regulatory Functions

**URL:** doi:10.1016/j.cell.2009.01.002; PMID 19167326
**Accessed:** 2026-02-10
**Authority rank:** 1 (Peer-reviewed, Cell 136:215-233)

**Key Extracted Points:**

1. **Mature miRNA:** ~22 nt, derived from one arm of ~70 nt hairpin precursor.
2. **Guide vs star strand:** Thermodynamic asymmetry determines which strand is loaded into RISC.

### Krol et al. (2004): Structural Features of miRNA Precursors

**URL:** doi:10.1074/jbc.M404931200; PMID 15292246
**Accessed:** 2026-02-10
**Authority rank:** 1 (Peer-reviewed, J Biol Chem 279:42230-42239)

**Key Extracted Points:**

1. **Stem length importance:** Stem length (typically 22-33 bp) is critical for Drosha and Dicer processing.
2. **Imperfect pairing:** Pre-miRNA hairpins feature internal mismatches, G:U wobble pairs, and bulges.
3. **Minimum stem:** Effective processing requires ~18-22 bp of stem (with allowed mismatches).

### miRBase (Griffiths-Jones et al., 2006)

**URL:** https://mirbase.org/
**Accessed:** 2026-02-10
**Authority rank:** 5 (Authoritative database)

**Key Extracted Points:**

1. **Canonical pre-miRNA lengths:** Database entries show pre-miRNA lengths ranging 55–120 nt.
2. **hsa-mir-21 (MI0000077):** pre-miRNA = 72 nt, mature 5p = 22 nt, sequence: `UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA`
3. **hsa-let-7a-1 (MI0000060):** pre-miRNA = 80 nt, mature 5p = 22 nt, sequence: `UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA`

---

## Documented Corner Cases and Failure Modes

### From Bartel (2004), Ambros (2003)

1. **Short sequences:** Sequences shorter than ~55 nt cannot form valid pre-miRNA hairpins.
2. **No complementarity:** Sequences without self-complementary regions yield no stem → no hairpin.
3. **Loop too small:** Loop < 3 nt is thermodynamically unfavorable and not observed in known pre-miRNAs.
4. **Loop too large:** Loop > ~25 nt — Drosha/Dicer may not process efficiently.
5. **Stem too short:** < 18 bp stem insufficient for enzyme recognition.
6. **T→U conversion:** DNA input must be converted to RNA (T→U) for analysis.

### From Krol et al. (2004)

1. **G:U wobble pairs in stem:** Must be counted as valid base pairs (common in pre-miRNAs).
2. **Multiple candidates in long sequences:** Scanning algorithm may find overlapping candidates.

---

## Test Datasets

### Dataset 1: Mixed-base Watson-Crick hairpin (57 nt)

**Source:** Constructed with biologically plausible mixed AUGC composition

| Parameter | Value |
|-----------|-------|
| 5' stem | `GCAUAGCUAGCUAGCUAGCUAGCUA` (25 nt) |
| Loop | `GAAAUUU` (7 nt) |
| 3' stem | `UAGCUAGCUAGCUAGCUAGCUAUGC` (25 nt, reverse complement) |
| Total length | 57 nt |
| Effective stem | 23 bp (capped by maxStem = n/2 - 5) |

### Dataset 2: Wobble-pair hairpin (57 nt)

**Source:** Constructed with 4 G:U wobble pairs per Krol (2004)

| Parameter | Value |
|-----------|-------|
| 5' stem | `GCAUAGCUAGCUAGCUAGCUAGCUA` (25 nt, same as Dataset 1) |
| Loop | `GAAAUUU` (7 nt) |
| 3' stem | `UAGCUAGUUAGCUAGUUAGUUAUGU` (25 nt, 4× C→U for wobble) |
| Total length | 57 nt |
| Wobble positions | Pairing positions 0, 5, 9, 17 (G:U instead of G:C) |

### Dataset 3: Real miRBase — hsa-mir-21 (MI0000077)

**Source:** miRBase v22 (authoritative database)

| Parameter | Value |
|-----------|-------|
| Sequence | `UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA` |
| Length | 72 nt |
| Expected result | **NOT detected** — known limitation of consecutive-pairing model |
| Reason | Internal mismatches and bulges yield only 16 consecutive pairs from ends (< 18 threshold) |

### Dataset 4: Real miRBase — hsa-let-7a-1 (MI0000060)

**Source:** miRBase v22 (authoritative database)

| Parameter | Value |
|-----------|-------|
| Sequence | `UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA` |
| Length | 80 nt |
| Expected result | **NOT detected** — known limitation of consecutive-pairing model |
| Reason | Only 5 consecutive pairs from ends (< 18 threshold) |

### Dataset 5: Short-stem hairpin (55 nt)

**Source:** Constructed to test stem-length rejection independently of n<55 early-exit

| Parameter | Value |
|-----------|-------|
| 5' stem | `GCAUAGCUAGCUAGC` (15 nt) |
| Loop | `AUGCAUGCAUGCAUGCAUGCAUGCA` (25 nt) |
| 3' stem | `GCUAGCUAGCUAUGC` (15 nt, reverse complement) |
| Total length | 55 nt (passes n<55 check) |
| Expected result | Rejected — stem 15 bp < 18 bp threshold |

---

## Design Limitations

1. **Consecutive stem pairing model** — The implementation requires consecutive base pairs scanning from the 5’/3’ ends inward, breaking at the first non-pairing position. Real pre-miRNAs have asymmetric internal loops (bulges) that offset the pairing alignment. A full RNA folding algorithm (Zuker/Nussinov) would be required to detect them. Tests M18/M19 document this limitation using real miRBase sequences.

---

## Recommendations for Test Coverage

1. **MUST Test:** Valid hairpin detection with known stem/loop structure — Evidence: Bartel (2004), Krol (2004)
2. **MUST Test:** Short sequence returns empty — Evidence: Pre-miRNA minimum ~55 nt (miRBase)
3. **MUST Test:** Null/empty input handling — Evidence: defensive coding
4. **MUST Test:** Stem length threshold (≥18 bp required) — Evidence: Krol (2004)
5. **MUST Test:** Loop size constraints (3–25 nt) — Evidence: Bartel (2004)
6. **MUST Test:** T→U conversion in input — Evidence: RNA biology standard
7. **MUST Test:** Mature and star sequence extraction from correct arms — Evidence: Bartel (2009)
8. **MUST Test:** Dot-bracket structure correctness — Evidence: Standard RNA structure notation
9. **MUST Test:** Free energy relative ordering (longer stem → more negative) — Evidence: Turner (2004) principles
10. **SHOULD Test:** Multiple hairpins in longer sequence — Evidence: miRNA clusters (Bartel 2004)
11. **SHOULD Test:** G:U wobble pairs accepted in stem — Evidence: Krol (2004)
12. **SHOULD Test:** maxHairpinLength parameter respected — Evidence: miRBase range
13. **MUST Test:** Real miRBase sequences NOT detected (known limitation) — Evidence: miRBase v22 (hsa-mir-21, hsa-let-7a-1)
14. **COULD Test:** Performance with long sequences — Rationale: O(n²) complexity

---

## References

1. Bartel DP (2004). MicroRNAs: Genomics, Biogenesis, Mechanism, and Function. Cell 116:281-297. doi:10.1016/S0092-8674(04)00045-5
2. Ambros V et al. (2003). A Uniform System for microRNA Annotation. RNA 9:277-279. doi:10.1261/rna.2183803
3. Bartel DP (2009). MicroRNAs: Target Recognition and Regulatory Functions. Cell 136:215-233. doi:10.1016/j.cell.2009.01.002
4. Krol J et al. (2004). Structural Features of miRNA Precursors and Their Relevance to miRNA Biogenesis. J Biol Chem 279:42230-42239. doi:10.1074/jbc.M404931200
5. Wikipedia: MicroRNA. https://en.wikipedia.org/wiki/MicroRNA (Accessed 2026-02-10)
6. miRBase. https://mirbase.org/ (Accessed 2026-02-10)
7. Griffiths-Jones S et al. (2006). miRBase: microRNA sequences, targets and gene nomenclature. Nucleic Acids Res 34:D140-144.
8. Turner DH, Mathews DH (2010). NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. Nucleic Acids Res 38:D280-282.

---

## MFE-Structure-Based Detection (opt-in) — Sources retrieved 2026-06-24

The opt-in `FindPreMiRnaHairpinsByMfe` / `AssessHairpinByMfe` methods fold each candidate with the
already-validated Zuker–Stiegler MFE engine (RNA-STRUCT-001,
`RnaSecondaryStructure.CalculateMfeStructure` / `CalculateMinimumFreeEnergy`, Turner 2004 NN model)
and derive the hairpin features (single terminal loop, paired-stem count, ΔG°, AMFE, MFEI) from the
ACTUAL MFE dot-bracket structure. No new thermodynamic parameters are introduced. The acceptance
thresholds below were retrieved this session from primary sources.

### Bonnet et al. (2004): microRNA precursors have lower folding free energies than random sequences

**Retrieved via:** WebFetch of https://pubmed.ncbi.nlm.nih.gov/15217813/ (2026-06-24).
**Authority rank:** 1 (Peer-reviewed, Bioinformatics 20(17):2911-2917, doi:10.1093/bioinformatics/bth374).

**Verbatim (abstract):** *"In contrast with transfer RNAs and ribosomal RNAs, the majority of the
microRNA sequences clearly exhibit a folding free energy that is considerably lower than that for
shuffled sequences, indicating a high tendency in the sequence towards a stable secondary
structure."* — justifies using a real MFE fold (a strongly negative ΔG° single hairpin) to recognise
pre-miRNAs, where the heuristic consecutive-pairing scan misses them.

### Zhang et al. (2006): Evidence that miRNAs are different from other RNAs (AMFE / MFEI)

**Retrieved via:** WebSearch + WebFetch of the PLOS One AMFE/MFEI definitions
(https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0113380, 2026-06-24) and the
cotton miRNA PLOS One paper quoting Zhang's criteria verbatim
(https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0033696, 2026-06-24), and the
miRkwood ab initio help (https://bioinfo.cristal.univ-lille.fr/mirkwood/abinitio/help.php).
**Authority rank:** 1 (Zhang BH, Pan XP, Cox SB, Cobb GP, Anderson TA (2006), Cell Mol Life Sci
63:246-254), supported by curated tools.

**Verbatim facts extracted:**

1. **AMFE** (adjusted minimal folding free energy): *"AMFE means the MFE of a RNA sequence with 100 nt
   in length, which is equal to MFE/(length of a potential pre-miRNA) * 100."* (PLOS One cotton paper;
   PLOS normalization paper: *"AMFE = 100·MFE/length … relate the index to a 100-nucleotides segment"*).
2. **MFEI** (minimal folding free energy index): *"MFEI is equal to MFE/(length of a potential
   pre-miRNA)/(The percentage of nucleotides G and C)"* = AMFE / (G+C)%. miRkwood: *"MFEI = [MFE /
   sequence length × 100] / (G+C%)."*
3. **MFEI threshold:** *"The value of MFEI for most miRNA precursors is greater than 0.85, which is
   remarkably higher than that of other non-coding small RNAs"* (Zhang et al. 2006, as quoted in the
   reviewed literature). Used as the default acceptance cutoff `minMfei = 0.85`.

**Sign convention used in the implementation:** Zhang computes MFEI from |MFE| (a positive
magnitude), so genuine pre-miRNAs yield MFEI > 0.85. The library's `CalculateMinimumFreeEnergy`
returns a **negative** ΔG°; `CalculateMfeIndex` therefore uses `Math.Abs(ΔG°)` so the published
"MFEI > 0.85" criterion applies directly. Documented so callers are not surprised.

### Ambros et al. (2003) — stem base-pairing criterion (re-applied to the MFE structure)

**Re-confirmed via:** WebFetch of the rnajournal abstract (https://rnajournal.cshlp.org/content/9/3/277.abstract)
and the Bartel-lab reprint search (2026-06-24). The fold-back hairpin must place the ~22 nt mature
miRNA in one arm base-paired to the opposite arm; the uniform-system criterion requires **≥16
complementary bases** to the opposite arm. Used as the minimum-stem cutoff
`MfeMinStemBasePairs = 16` on the count of base pairs in the MFE structure's single hairpin.

### Meyers/Bartel et al. (2008) — minimal bulges / single-arm duplex

**Retrieved via:** WebFetch of https://pmc.ncbi.nlm.nih.gov/articles/PMC2630443/ (2026-06-24).
**Verbatim:** *"base-pairing between the miRNA and the other arm of the hairpin, which includes the
miRNA*, is extensive such that there are typically four or fewer mismatched miRNA bases"* and
*"asymmetric bulges are minimal in size (one or two bases) and frequency"* — supports requiring a
**single dominant hairpin** (one terminal loop, nested stem with only small internal loops/bulges,
no multibranch) as the structural acceptance test.

### Reference values from the RNA-STRUCT-001 engine (reused, not re-derived)

Folded this session with `RnaSecondaryStructure.CalculateMfeStructure` (Turner 2004); ΔG° equals
`CalculateMinimumFreeEnergy` exactly. These are the exact expected values encoded in the tests:

| Candidate | Length | GC% | ΔG° (kcal/mol) | MFE dot-bracket (summary) | Stem bp | Loop | AMFE | MFEI | Verdict |
|-----------|--------|-----|----------------|----------------------------|---------|------|------|------|---------|
| `ValidHairpin57` (`GCAUAGCUAGCUAGCUAGCUAGCUA`+`GAAAUUU`+`UAGCUAGCUAGCUAGCUAGCUAUGC`) | 57 | 43.8596 | **−48.48** | 27×`(` 3×`.` 27×`)` | 27 | 3 | 85.052632 | **1.939200** | ACCEPT |
| `hsa-mir-21` (MI0000077) | 72 | 48.6111 | **−35.13** | single hairpin w/ internal loops, apical loop 3 | 32 | 3 | 48.791667 | **1.003714** | **ACCEPT** (heuristic rejects) |
| `hsa-let-7a-1` (MI0000060, 80 nt) | 80 | 42.5000 | **−34.31** | single hairpin w/ internal loops, apical loop 4 | 32 | 4 | 42.887500 | **1.009118** | **ACCEPT** (heuristic rejects) |
| `5S-rRNA-like` (120 nt, multibranch) | 120 | 64.1667 | **−47.04** | multibranch (`)`…`(`) | — | — | — | — | **REJECT** (not a single hairpin, despite strong ΔG°) |
| `NoComplementarity` (70 nt A/G) | 70 | 50.0 | **0.00** | all `.` | 0 | — | 0 | 0 | REJECT |

The `5S-rRNA-like` case proves the assessment rejects on **structure** (multibranch), not merely on a
weak ΔG°: its ΔG° (−47.04) is more negative than `ValidHairpin57`'s, yet it is correctly rejected
because the MFE fold is not a single dominant hairpin.

---

## Drosha/Dicer Cleavage-Site Prediction (opt-in) — Sources retrieved 2026-06-24

The opt-in `PredictDroshaDicerCleavage` method predicts the Drosha and Dicer cut coordinates and the
resulting mature-miRNA / miRNA* spans from a pri-miRNA hairpin, using the PUBLISHED measuring rules
(distances from the basal junction and from the Drosha-generated 5' end). No trained model is used.
Every distance/overhang below was retrieved verbatim this session.

### Han et al. (2006): Molecular Basis for the Recognition of Primary microRNAs by the Drosha-DGCR8 Complex

**Retrieved via:** WebFetch of https://pubmed.ncbi.nlm.nih.gov/16751099/ (2026-06-24).
**Authority rank:** 1 (Peer-reviewed, Cell 125(5):887-901, doi:10.1016/j.cell.2006.03.043).

**Verbatim (abstract):** *"A typical metazoan pri-miRNA consists of a stem of approximately 33 bp, with
a terminal loop and flanking segments. The terminal loop is unessential, whereas the flanking ssRNA
segments are critical for processing. **The cleavage site is determined mainly by the distance
(approximately 11 bp) from the stem-ssRNA junction.** Purified DGCR8, but not Drosha, interacts with
pri-miRNAs both directly and specifically, and the flanking ssRNA segments are vital for this binding
to occur. Thus, DGCR8 may function as the molecular anchor that measures the distance from the
dsRNA-ssRNA junction."*

1. **Drosha ruler distance:** ~**11 bp** from the basal stem-ssRNA (dsRNA-ssRNA) junction = the
   constant `DroshaCutBpFromBasalJunction = 11`.
2. **Stem length:** ~**33 bp** total stem.
3. **One helical turn:** ~11 bp ≈ one RNA helical turn (the "ruler" measured by DGCR8/Drosha) — also
   confirmed in Lee et al. (2003) and the Cell 2006 microprocessor review (search-result quote: *"count
   up ∼11 bp, one helical RNA turn, to the scissile phosphodiester bond"*).

### Lower-stem / upper-stem geometry (re-confirmed via search 2026-06-24)

**Retrieved via:** WebSearch (2026-06-24); corroborating reviews — *Drosha and Dicer: Slicers cut from
the same cloth*, Cell Research 2016 (https://www.nature.com/articles/cr201619) and the NAR lower-stem
paper (https://academic.oup.com/nar/article/48/5/2579/5709711).

**Verbatim facts extracted:**

1. *"The stem is divided into an upper stem where the mature miRNA sequence is embedded and a lower
   stem next to the flanking region."*
2. *"Drosha senses the basal junction between the lower stem and the single-stranded basal region,
   establishing a cleavage site 11 bp away."* — Drosha cut = basal junction + 11 bp.
3. *"Microprocessor cleaves at a distance of approximately 11-bp from the basal junction and at a
   distance of approximately 22-bp from the apical junction."* — the lower stem (~11 bp) is below the
   Drosha cut; the upper stem (~22 bp) carries the mature miRNA.
4. The Drosha cut produces the base of the pre-miRNA: the **5' end of the 5p mature** on the 5' arm and
   the **3' end of the 3p mature** on the 3' arm.

### Park et al. (2011): Dicer recognizes the 5' end of RNA for efficient and accurate processing

**Retrieved via:** WebFetch of https://pubmed.ncbi.nlm.nih.gov/21753850/ (2026-06-24).
**Authority rank:** 1 (Peer-reviewed, Nature 475(7355):201-205, doi:10.1038/nature10198).

**Verbatim (abstract):** *"… human Dicer anchors not only the 3' end but also the 5' end, with **the
cleavage site determined mainly by the distance (∼22 nucleotides) from the 5' end (5' counting
rule)**. … Dicer recognizes the 5'-phosphorylated end … a class of approximately 22-nucleotide RNAs."*

1. **Dicer ruler distance:** ~**22 nt** from the (Drosha-generated) 5' end = the constant
   `DicerCutNtFrom5PrimeEnd = 22`; this fixes the mature length at ~22 nt.
2. **Mature length:** ~**22 nt** (the miRNA:miRNA* duplex is ~22 nt; the 21–23 nt range is the
   accepted mature length, Wikipedia/Bartel).
3. Dicer's cut is *"determined largely by the Drosha cleavage site"* — i.e. Dicer measures from the
   Drosha-generated end, removing the terminal loop.

### RNase III 2-nt 3' overhang (Drosha and Dicer signature)

**Retrieved via:** WebSearch (2026-06-24); ScienceDirect RNase III overview and the EMBO Lee/Han
processing paper (https://www.embopress.org/doi/full/10.1038/sj.emboj.7600491).

**Verbatim facts extracted:**

1. *"Cleavage by RNase III domains results in 2-nt 3′-overhang end with a 5′-terminal monophosphate,
   and a 3′-hydroxyl in the product RNA. The two RNase III domains form an intramolecular heterodimer
   and make staggered cleavages in the two arms of a pri-miRNA or pre-miRNA."* — both Drosha and Dicer
   leave a **2-nt 3' overhang** = the constant `RNaseIII3PrimeOverhang = 2`.
2. *"Drosha … excises the upper part of this RNA hairpin to generate the precursor miRNA (pre-miRNA),
   which is ∼60 nt long with a 3′ 2 nt overhang."*
3. *"Dicer recognizes the 3′ 2 nt overhang of pre-miRNAs and then cuts ∼22 nt away to produce the
   miRNA:miRNA* duplex."*

### Auyeung et al. (2013): primary-sequence motifs that position Drosha (refinement)

**Retrieved via:** WebSearch + WebFetch of the PMC "Menu of features" review
(https://pmc.ncbi.nlm.nih.gov/articles/PMC4613790/, 2026-06-24); primary paper Cell 152(4):844-858,
doi:10.1016/j.cell.2013.01.031 (https://pubmed.ncbi.nlm.nih.gov/23415231/).
**Authority rank:** 1.

**Verbatim facts extracted (motifs that fine-position the Drosha cut):**

1. **Basal UG motif:** *"at the base of the pri-miRNA hairpin, at the junction between single-stranded
   and double-stranded RNA regions."* — anchors the basal junction.
2. **Apical UGU(G) motif:** *"a UGU/GUG motif in the apical loop"* (near the terminal loop).
3. **CNNC motif:** *"positioned 16–18 nt from the Drosha cut"* (downstream of the basal cleavage site)
   = the constants `CnncMinNtDownstreamOfDroshaCut = 16`, `CnncMaxNtDownstreamOfDroshaCut = 18`. Used
   only as an OPTIONAL confidence flag (reported, not required); when a CNNC (C-N-N-C) is present in
   the 16–18 nt window 3' of the Drosha basal cut, `HasCnncMotif = true`.

### Cross-check dataset: hsa-mir-21 (miRBase MI0000077)

**Retrieved via:** WebFetch of https://mirbase.org/hairpin/MI0000077 (2026-06-24).

**Verbatim from miRBase:**

- Precursor (stem-loop), 72 nt (lowercase = flanking/loop ssRNA as annotated by miRBase):
  `ugucgggUAGCUUAUCAGACUGAUGUUGAcuguugaaucucauggCAACACCAGUCGAUGGGCUGUcugaca`
- **hsa-miR-21-5p:** positions **8–29** (1-based) = `UAGCUUAUCAGACUGAUGUUGA` (22 nt).
- **hsa-miR-21-3p:** positions **46–66** (1-based) = `CAACACCAGUCGAUGGGCUGU` (21 nt).

**Derivation used as the expected values in the cross-check test.** The 7-nt basal ssRNA tail
`ugucggg` (positions 1–7) ends at the basal junction; the 5' stem (upper part of the lower stem) starts
at position 8. miRBase annotates the mature **5p starting exactly at the base of the stem (position 8,
0-based index 7)** — i.e. at the Drosha cut. With the model fed the basal-junction 0-based index 7 (the
first paired base of the 5' arm, after trimming the lower stem as miRBase does for this entry):

| Quantity | Rule | Value (0-based) |
|---|---|---|
| Drosha 5' cut (5' end of 5p) | basal junction (position 7) | index **7** |
| 5p mature span | Dicer 5' counting: 22 nt | indices **7–28** = `UAGCUUAUCAGACUGAUGUUGA` |
| 5p length | ~22 nt | **22** ✓ matches miRBase |
| 2-nt 3' overhang | RNase III | 5p 3' end protrudes 2 nt over the 3p 5' end |

This confirms the predicted 5p span and length match miRBase exactly. (The model's general path takes a
pri-miRNA with explicit flanking ssRNA and applies the +11 bp Drosha ruler; the miR-21 entry is the
already-trimmed special case where the annotated stem base = the Drosha cut, validating the +22 nt
Dicer ruler and the 22-nt mature length against the database.)

### Honest residual

The **trained natural-vs-background miRNA classifier** (miRDeep2-style: a probabilistic/ML score that
distinguishes genuine pre-miRNAs from random hairpins using read-stacking signatures and a fitted
model) is **NOT** implemented — it requires a trained model and labelled training data, which are
data-blocked. The cleavage-site prediction here is the published deterministic measuring-rule
("ruler") heuristic only.

> **Superseded 2026-06-25.** A trained **structure/sequence-feature** natural-vs-background classifier
> is now bundled (see the next section). The remaining honest residual is only the **read-stacking**
> (small-RNA-seq pileup) signal of miRDeep2, which needs the caller's sequencing reads and is out of
> scope.

---

## Trained natural-vs-background classifier (opt-in) — Sources retrieved 2026-06-25

A trained logistic-regression classifier that distinguishes genuine pre-miRNA hairpins from
length/composition-matched background was added. It is built ONLY from public-domain data and the
published method; **no GPL miRDeep2 code or weights were ported, copied, or linked** (miRDeep2's
published method was consulted as an offline reference only). Canonical methods:
`MiRnaAnalyzer.ExtractPreMiRnaFeatures`, `ClassifyPreMiRna`, `ScorePreMiRnaFeatures`,
`DinucleotideShuffle`.

### miRBase — licence (public domain) and training positives

**URL:** https://www.re3data.org/repository/r3d100010566 and https://mirbase.org/download/CURRENT/LICENSE/
**Accessed:** 2026-06-25
**Authority rank:** 5 (curated database) + its own licence file

**Key Extracted Points:**

1. **Licence (verbatim, from `mirbase.org/download/CURRENT/LICENSE/`):** *"miRBase is in the public
   domain. It is not copyrighted. You may freely modify, redistribute, or use it for any purpose."*
   re3data record r3d100010566 lists the data licence as **Public Domain** with **open** access.
   → miRBase sequence data is usable and redistributable. Confirmed.
2. **Positives (13 real human pre-miRNA hairpin precursors, retrieved verbatim this session from
   `mirbase.org/hairpin/<MI>`):** hsa-mir-21 (MI0000077, 72 nt), hsa-let-7a-1 (MI0000060, 82 nt),
   hsa-let-7b (MI0000063, 81 nt), hsa-mir-106a (MI0000113, 80 nt), hsa-mir-31 (MI0000089, 72 nt),
   hsa-let-7f-2 (MI0000068, 82 nt), hsa-mir-16-1 (MI0000070, 89 nt), hsa-mir-24-1 (MI0000080, 67 nt),
   hsa-mir-93 (MI0000095, 80 nt), hsa-mir-96 (MI0000098, 77 nt), hsa-mir-98 (MI0000100, 117 nt),
   hsa-mir-147a (MI0000262, 72 nt), hsa-mir-199a-2 (MI0000281, 109 nt). The two named test precursors
   (hsa-mir-21, hsa-let-7a-1) and their sequences are quoted in Datasets 3–4 above.

### Bonnet et al. (2004) — di-shuffle negative-set convention

**URL:** https://doi.org/10.1093/bioinformatics/bth374 (Bioinformatics 20(17):2911-2917)
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **Method (verbatim from the abstract retrieved this session):** *"the predicted minimum free energy
   of folding is compared with values obtained for structures inferred from randomly shuffling the
   original sequences."*
2. **Conclusion (verbatim):** *"the majority of the microRNA sequences clearly exhibit a folding free
   energy that is considerably lower than that for shuffled sequences, indicating a high tendency in
   the sequence towards a stable secondary structure."* → the canonical pre-miRNA-vs-shuffled MFE
   discrimination; shuffled, composition-matched sequences are the standard negative set.

### Xue et al. (2005) triplet-SVM / Batuwita & Palade (2009) microPred — features & negative convention

**URL:** https://academic.oup.com/bioinformatics/article/25/8/989/324698 (Bioinformatics 25(8):989-995);
Xue et al. BMC Bioinformatics 6:310 (doi:10.1186/1471-2105-6-310)
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **microPred features (retrieved this session):** sequential features (%C+G content, dinucleotide
   frequencies), **thermodynamic** measures (**MFE, MFEI** variants, normalised ensemble free energy),
   and **structural** descriptors (**base-pair composition / %paired**, stem characteristics, structural
   diversity). → MFE, MFEI, GC%, %paired, stem/loop geometry are the published feature family.
2. **Negative set:** microPred itself used genomic pseudo-hairpins; the **di-nucleotide-shuffle**
   background is the convention established by Bonnet (2004) and used throughout the pre-miRNA-classifier
   literature (Xue 2005 negatives are composition-controlled pseudo-hairpins). We adopt the di-shuffle
   negative (Bonnet 2004) because it is reproducible from the positives alone.

### Zhang et al. (2006) — MFEI / AMFE (already cited above)

MFEI = AMFE / (G+C)% with AMFE = 100·|MFE|/length — Zhang BH et al. (2006), Cell Mol Life Sci
63:246-254 (see "Zhang et al. (2006)" section above for the verbatim extract).

### Altschul & Erickson (1985) — dinucleotide-preserving shuffle algorithm

**URL:** uShuffle (Jiang et al. 2008, BMC Bioinformatics 9:192, doi:10.1186/1471-2105-9-192) describing
the Altschul-Erickson method; primary: Altschul SF, Erickson BW (1985) Mol Biol Evol 2(6):526-538.
**Accessed:** 2026-06-25
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points (retrieved this session):**

1. The Altschul-Erickson (1985) method is *"a method for random sequence permutation that preserves
   dinucleotide and codon usage"*; it performs an **Eulerian walk on a directed multigraph** built from
   the original sequence; *"the character labels on a random Euler tour of this graph (crossing each
   edge once) exactly preserves the dinucleotide frequencies."* The first and last characters are fixed.
   → `DinucleotideShuffle` implements exactly this (verified: identical dinucleotide counts, same
   first/last/length).

### Logistic regression — model & fit

**URL:** Hastie T, Tibshirani R, Friedman J, *The Elements of Statistical Learning* 2nd ed. §4.4.1,
doi:10.1007/978-0-387-84858-7; gradient-ascent on the log-likelihood (Stanford CS109 lecture notes,
web.stanford.edu/class/archive/cs/cs109).
**Accessed:** 2026-06-25
**Authority rank:** 1 (textbook)

**Key Extracted Points:** logistic model P(y=1|x) = σ(b + wᵀx), σ(z)=1/(1+e⁻ᶻ); fit by maximising the
log-likelihood Σ[ yᵢ log pᵢ + (1−yᵢ) log(1−pᵢ) ]; gradient ∂ℓ/∂wⱼ = Σᵢ (yᵢ − pᵢ) xᵢⱼ; iterate
w ← w + η·∇ℓ (batch gradient ascent). Features standardised by training-set mean/std. We add a small
L2 (ridge) penalty (λ=1e-3) — standard regularised logistic regression (ESL §4.4.4).

---

## Training procedure (reproducible) — `Training`

- **Positives:** the 13 miRBase precursors above (public domain).
- **Negatives:** `DinucleotideShuffle` (Altschul-Erickson) of each positive, **4 shuffles per
  positive** → 52 negatives. **52 + 13 = 65** examples; balanced enough (Bonnet-style) and strongly
  separable.
- **RNG seed:** fixed **20060101** (one shared `Random`, consumed in fixed positive→shuffle order).
- **Features (model inputs, in order):** `[FreeEnergy, AMFE, MFEI, GcContent, PairedFraction]`,
  standardised by the TRAIN-set mean/std. (`StemBasePairs`, `LoopSize`, `Length` are extracted and
  reported but excluded from the model — collinear / length-confounded.)
- **Split:** deterministic — first **70 %** by example index = **train (45)**, last **30 %** =
  **held-out test (20)**; the interleaved positive/shuffle order guarantees both classes appear in the
  held-out split (4 positives + 16 shuffles).
- **Fit:** batch gradient ascent, learning rate 0.1, 20 000 epochs, L2 λ=1e-3.
- **Bundled constants** (`MiRnaAnalyzer.cs`, regenerated by the offline trainer):
  - bias = `-4.340788257901692`
  - weights = `[-0.6299061497891402, 2.4689410267337104, 2.5770103435217253, -0.3315361036977469, 2.42751798356881]`
  - mean = `[-22.988666666666667, 29.297338207740115, 0.6157639504286563, 0.477134229110014, 0.6153156920227864]`
  - std = `[8.373410721510746, 10.37385177977526, 0.2084972456212268, 0.061249038538481466, 0.10661272255434924]`
- **Held-out metric (deterministic):** **accuracy = 1.0000, AUC = 1.0000** on the 20 held-out
  examples (real pre-miRNAs are well-separated from di-shuffled controls — Bonnet 2004; AUC ≳ 0.90
  expected and exceeded).

### Exact feature values (pins the extractor independent of the weights)

| Precursor | MFE (ΔG°) | AMFE | MFEI | GC | %paired | bp | loop | len |
|---|---|---|---|---|---|---|---|---|
| hsa-mir-21 (MI0000077) | −35.13 | 48.79166666666667 | 1.0037142857142858 | 0.48611111111111116 | 0.8888888888888888 | 32 | 3 | 72 |

(`ΔG°` is `RnaSecondaryStructure.CalculateMinimumFreeEnergy` verbatim, RNA-STRUCT-001 Turner 2004.)

### Discrimination spot-check

| Sequence | P(natural) | call |
|---|---|---|
| hsa-mir-21 (natural) | 0.9999949074161931 | natural ✓ |
| hsa-let-7a-1 (natural) | 0.9998877289367469 | natural ✓ |
| hsa-mir-21 di-shuffled (seed 999) | 0.0005220422779554753 | background ✓ |
