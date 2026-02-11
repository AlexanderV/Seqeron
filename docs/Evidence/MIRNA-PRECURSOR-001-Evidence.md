# Evidence Artifact: MIRNA-PRECURSOR-001

**Test Unit ID:** MIRNA-PRECURSOR-001
**Algorithm:** Pre-miRNA Hairpin Detection
**Date Collected:** 2026-02-10

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
2. **hsa-mir-21 (MI0000077):** pre-miRNA = 71 nt, mature 5p = 22 nt, sequence: `UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA`
3. **hsa-let-7a-1 (MI0000060):** pre-miRNA = 78 nt, mature 5p = 22 nt, sequence: `UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA`

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
| Length | 71 nt |
| Expected result | **NOT detected** — known limitation of consecutive-pairing model |
| Reason | Internal mismatches and bulges yield only ~16 consecutive pairs from ends (< 18 threshold) |

### Dataset 4: Real miRBase — hsa-let-7a-1 (MI0000060)

**Source:** miRBase v22 (authoritative database)

| Parameter | Value |
|-----------|-------|
| Sequence | `UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA` |
| Length | 78 nt |
| Expected result | **NOT detected** — known limitation of consecutive-pairing model |
| Reason | Only ~5 consecutive pairs from ends (< 18 threshold) |

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

## Assumptions

1. **ASSUMPTION: Energy model simplified** — The implementation uses a simplified energy model (-stemLength * 1.5 + loopSize * 0.5) rather than the Turner nearest-neighbor parameters. Tests should verify relative ordering (longer stem → more negative energy, larger loop → less negative energy) rather than absolute kcal/mol values.

2. **ASSUMPTION: Consecutive stem requirement** — The implementation requires consecutive base pairs from the ends inward (breaks at first non-pairing position). Real pre-miRNAs have bulges/mismatches. Tests should verify behavior against the implementation's definition, not against a full thermodynamic folding model.

3. **ASSUMPTION: ValidateHairpin is AnalyzeHairpin** — The checklist lists `ValidateHairpin(structure)` as a method, but the implementation has a private `AnalyzeHairpin(sequence, matureLength)`. The validation logic is tested indirectly via `FindPreMiRnaHairpins`.

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

---

## Change History

- **2026-02-11**: Replaced synthetic test datasets with real miRBase sequences (MI0000077, MI0000060) and biologically plausible mixed-base hairpins. Added known-limitation documentation for consecutive-pairing model.
- **2026-02-10**: Initial documentation.
