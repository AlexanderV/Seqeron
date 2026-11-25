# Donor (5') Splice Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Splicing |
| Test Unit ID | SPLICE-DONOR-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Shapiro & Senapathy 1987 score; MaxEntScan mode) |
| Last Reviewed | 2026-11-25 |

## 1. Overview

Donor splice-site detection identifies candidate 5' splice sites at the start of introns. In this repository, `SpliceSitePredictor.FindDonorSites` scans a sequence for canonical `GU` donor dinucleotides and, when requested, optional `GC` and U12-style `AU` donors. Candidates are scored with the Shapiro & Senapathy (1987) 5' splice-site score (default; Table 1 percentages, primate subtable) or the Yeo & Burge (2004) MaxEntScan `score5` (opt-in) and filtered by `minScore`. The method is deterministic, linear in sequence length, and returns a `SpliceSite` record for each retained candidate.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The donor splice site marks the 5' end of an intron. The original document records the extended consensus motif `MAG|GURAGU`, where the canonical intronic start is `GU` and non-canonical `GC` and U12-type `AU` donors are recognized in smaller subsets of introns.[1][2]

| Position | Consensus |
|----------|-----------|
| `-3` | `M` (`A` or `C`) |
| `-2` | `A` |
| `-1` | `G` |
| `0` | `G` |
| `+1` | `U` |
| `+2` | `R` (`A` or `G`) |
| `+3` | `A` |
| `+4` | `G` |
| `+5` | `U` |

### 2.2 Core Model

**Shapiro & Senapathy (1987) score (default).** Table 1 of the paper gives the percentage of each nucleotide at
5' ss positions −3..+6 for seven GENBANK categories and their total (primate n = 542, rodent 393, other mammal 65,
other vertebrate 209, invertebrate 80, plant 111, viral 46; total 1446). Methods (p. 7157–7158): every
subsequence of length 8 (positions −2..+6; −3 is tabulated but not scored) is scored as

`score = 100 (t − mint) / (maxt − mint)`

where `t` is the sum of the percentages of the observed nucleotides and `mint`/`maxt` the sums of the lowest /
highest percentages of the eight positions. Worked example (primate `AAGTGAGT`): t = 58+10+100+100+39+71+84+47 = 509,
mint = 47, maxt = 595 → 84.3. The library reports `score / 100` in [0, 1]. Because +1/+2 are 100 % G/T, a GC donor
loses the full +2 range (`CAG|GCAAGU` → 448/548 = 0.818 vs 1.0 for `CAG|GUAAGU`).

**MaxEntScan (opt-in, `DonorScoringModel.MaxEntScan`).** `ScoreDonorMaxEnt` (score5, bits) of the 9-nt window
−3..+6; see [SPLICE-MAXENT5-001]. Confidence = 2^s/(1+2^s) (posterior at prior odds 1 of the log2 likelihood ratio).

**U12-type donors.** Scored as the fraction of the 8 intron positions +1..+8 matching the U12 consensus
`/RTATCCTY` (Sharp & Burge 1997; Burge et al. 1999; Hall & Padgett 1994 `ATATCCTT`). Both AT-AC and GT-AG
U12 introns exist (Dietrich et al. 1997), so AT and GT starts are candidates.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Positions are independent (S&S position frequency score). | Dependencies between positions (captured by MaxEntScan) are ignored in the default model. |
| ASM-02 | `GC` and U12-style `AU` donors should be optional rather than always included. | Canonical-only scans can miss these alternative donor classes by design. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Shapiro-Senapathy and U12 `Score` values are in `[0, 1]` (MaxEntScan scores are bits). | (t − mint)/(maxt − mint) and matches/8. |
| INV-02 | Returned `Confidence` values are in `[0, 1]`. | `CalculateConfidence` clamps the result. |
| INV-03 | Canonical and `GC` donors are returned with `Type = Donor`; U12 `AU` donors use `Type = U12Donor`. | The public method assigns those enum values explicitly. |
| INV-04 | Higher `minScore` produces a subset of lower-threshold results. | Filtering is applied as `if (score >= minScore)`. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA or RNA sequence to scan. | Empty or `null` yields no results; a candidate whose window (S&S −2..+6, MaxEnt −3..+6, U12 +1..+8) runs off the sequence or holds a non-ACGTU character is not scored. |
| `minScore` | `double` | `0.5` | Minimum donor score required to emit a site. | Model units: [0,1] (S&S, U12) or bits (MaxEnt). |
| `includeNonCanonical` | `bool` | `false` | Whether to include `GC` donors and U12-type (`AU`/`GU` + /RUAUCCUY) donors. | Canonical `GU` scanning is always active. |
| `model` | `DonorScoringModel` | `ShapiroSenapathy` | Overload only: `ShapiroSenapathy` or `MaxEntScan`. | Unknown value → `ArgumentOutOfRangeException`. |
| `category` | `ShapiroSenapathyCategory` | `Primate` | Overload only: Table 1 subtable. | Unknown value → `ArgumentOutOfRangeException`. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | 0-based index of the G of the donor dinucleotide = first intron nucleotide (+1). |
| `Type` | `SpliceSiteType` | `Donor` for U2-type `GU`/`GC`, `U12Donor` for U12-type `AU`/`GU` candidates. |
| `Motif` | `string` | RNA-alphabet context −3..+6 (U12: −3..+8), truncated at the sequence ends. |
| `Score` | `double` | S&S score/100, MaxEntScan bits, or U12 match fraction. |
| `Confidence` | `double` | clamp((score − 0.5)/0.5) (S&S, U12) or 2^s/(1+2^s) (MaxEnt). |

### 3.3 Preconditions and Validation

The method uppercases the input and converts `T` to `U` before scanning. A candidate is scored only when its full
window lies inside the sequence and contains only A/C/G/T/U (no partial-window scoring: the paper scores whole
8-nt subsequences; MaxEntScan rejects other lengths / letters). An unknown `model` / `category` enum value throws
`ArgumentOutOfRangeException` eagerly.

## 4. Algorithm

### 4.1 High-Level Steps

1. Return no results for `null` or empty input.
2. Normalize the sequence to uppercase RNA notation.
3. For every dinucleotide `(i, i+1)`: if `GU` (or `GC` with `includeNonCanonical`), score the U2-type donor with the
   selected model; emit it as `Donor` if the window is scorable and `score ≥ minScore`.
4. With `includeNonCanonical` and the S&S model, if the dinucleotide is `AU` or `GU`, score the U12 consensus
   over +1..+8 and emit `U12Donor` when `score ≥ minScore`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Candidate class | Recognition rule | Scoring rule |
|-----------------|------------------|--------------|
| Canonical donor | `GU` at `i`, `i+1` | S&S (1987) score over −2..+6 (or MaxEntScan score5 over −3..+6) |
| `GC` donor | `GC` at `i`, `i+1` (non-canonical) | Same scorer as canonical sites |
| U12 donor | `AU` or `GU` at `i`, `i+1` (non-canonical, S&S model) | Fraction of +1..+8 matching `/RUAUCCUY` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindDonorSites` | `O(n)` | `O(1)` auxiliary | Single pass over the normalized sequence with streaming output. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SpliceSitePredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `SpliceSitePredictor.FindDonorSites(string, double, bool)` (S&S primate)
- `SpliceSitePredictor.FindDonorSites(string, double, bool, DonorScoringModel, ShapiroSenapathyCategory)`
- `SpliceSitePredictor.ScoreDonorSiteShapiroSenapathy` / `ScoreDonorSiteMaxEnt` / `ScoreU12DonorSite` (private helpers)
- `SpliceSitePredictor.ScoreDonorMaxEnt(string)` (opt-in Yeo & Burge 2004 MaxEntScan `score5ss` maximum-entropy 5' donor score, in bits; embedded probability table)
- `SpliceSitePredictor.CalculateMaxEntScore(string, SpliceSiteType.Donor)` (MCP `maxent_score`) delegates to `ScoreDonorMaxEnt` (9-nt window required)

### 5.2 Current Behavior

Canonical `GU` and optional `GC` donors are emitted with `Type = Donor`; the U2 / GC-AG distinction is made during
intron classification. With `includeNonCanonical`, a `GU` may also be emitted as a `U12Donor` candidate (GT-AG
U12 introns). Real human donors (Ensembl canonical transcripts of HBB, HBA1, TP53, CFTR; 40 donors) score
0.684–0.967 (mean 0.84) with the primate table, versus ~0.69 of all GT dinucleotides scoring ≥ 0.5.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical donor detection is anchored on `GU` and the extended donor consensus around the exon-intron junction.[1][2]
- Optional support exists for `GC` donors and U12-style donor recognition.[2]
- The opt-in `ScoreDonorMaxEnt` implements the Yeo & Burge (2004) **MaxEntScan `score5ss`** maximum-entropy 5' donor model on a 9-nt window (3 exon + 6 intron, conserved `GT` at 0-based positions 3–4): `log2(P_maxent/P_background)`, removing the `GT` (scored by a consensus/background model) and looking up the maximum-entropy probability of the 7 remaining positions (`window[0:3] + window[5:9]`) directly in a single embedded table. Unlike score3, score5 is single-matrix (no overlapping sub-windows). The factorisation and table are the MIT-licensed maxentpy port; the canonical documented example `ScoreDonorMaxEnt("cagGTAAGT")` reproduces `10.86` bits exactly (provenance + licence: `Data/maxent_score5.LICENSE.md`). Bit-identical to maxentpy `score5` on all 4^9 ACGT 9-mers and identical at 2 dp to score5.pl (log2 computed as ln/ln2, as the references do).[4]

- The default `FindDonorSites` scorer is the published Shapiro & Senapathy (1987) 5' ss score with all eight
  Table 1 subtables (B11 review 2026-11-25; reproduces 39/42 published Table 5 scores at 1 dp). `FindDonorSites`
  also offers the MaxEntScan score5 model.

**Intentionally simplified / not implemented:**

- U12-type donors are scored by consensus match fraction (/RTATCCTY); a trained U12 PWM classifier
  (SpliceRack, Sheth et al. 2006; intronIC, Moyer et al. 2020 — GPL-3, not bundleable) is not implemented.
- Each donor site is scored independently from local motif context only; exon structure and longer-range
  splicing signals are not considered at this stage.
- A distinct `SpliceSiteType` for `GC` donors; **users should rely on:** later intron classification in `DetermineIntronType(...)`.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No `GU` in canonical mode | Returns no donor sites. | No canonical donor candidate is available. |
| Window runs off the sequence (e.g. GU at index < 2 or within the last 5 nt) | Not scored, no site. | The S&S subsequence is −2..+6. |
| Non-ACGTU character in the scored window | Not scored, no site. | No Table 1 percentage. |
| DNA input with `T` | Same behavior as RNA input with `U`. | `T` is converted to `U`. |
| Lowercase input | Same behavior as uppercase input. | The sequence is uppercased before scanning. |
| Multiple donor motifs | Each candidate is scored independently. | The scan iterates over every candidate start position. |

### 6.2 Limitations

The default donor-site detection is the Shapiro & Senapathy (1987) position-frequency score (positions independent). It does not decide intron structure on its own. As an opt-in model it offers the Yeo & Burge (2004) MaxEntScan `score5ss` maximum-entropy 5' donor score (`ScoreDonorMaxEnt`, on a 9-nt window, embedded MIT-licensed table — see §5.3). The score5 model is human/mammalian-trained; other species may have different donor statistics.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SpliceSitePredictor_DonorSite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/SpliceSitePredictor_DonorSite_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [SPLICE-DONOR-001.md](../../../tests/TestSpecs/SPLICE-DONOR-001.md)
- Related algorithms: [Acceptor_Site_Detection.md](Acceptor_Site_Detection.md), [Gene_Structure_Prediction.md](Gene_Structure_Prediction.md)

## 8. References

1. Shapiro MB, Senapathy P. 1987. RNA splice junctions of different classes of eukaryotes. Nucleic Acids Research. doi:10.1093/nar/15.17.7155
2. Burge CB, Tuschl T, Sharp PA. 1999. Splicing precursors to mRNAs by the spliceosomes. The RNA World. N/A
3. Test specification: [SPLICE-DONOR-001.md](../../../tests/TestSpecs/SPLICE-DONOR-001.md)
4. Yeo G, Burge CB. 2004. Maximum entropy modeling of short sequence motifs. Journal of Computational Biology. doi:10.1089/106652704773135290
5. Wikipedia contributors. 2026. RNA splicing. Wikipedia. https://en.wikipedia.org/wiki/RNA_splicing
6. Wikipedia contributors. 2026. Spliceosome. Wikipedia. https://en.wikipedia.org/wiki/Spliceosome
