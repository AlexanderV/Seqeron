# Acceptor Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Splicing |
| Test Unit ID | SPLICE-ACCEPTOR-001 |
| Related Projects | N/A |
| Implementation Status | Complete (Shapiro & Senapathy 1987 3' ss score; MaxEntScan mode; yUnAy / UACUAAC branch points) |
| Last Reviewed | 2026-11-25 |

## 1. Overview

Acceptor splice-site detection identifies candidate 3' splice sites at the ends of introns. `SpliceSitePredictor.FindAcceptorSites` scans a sequence for `AG` acceptors and scores them with the Shapiro & Senapathy (1987) 3' splice-site score (default; Table 1 percentages, primate subtable) or the Yeo & Burge (2004) MaxEntScan `score3` (opt-in). With `includeNonCanonical`, U12-type acceptors (`AC` or `AG` terminus) are scored by their branch site. Candidates meeting `minScore` are returned as `SpliceSite` records.[1][3]

Branch points are detected by `FindBranchPoints` (scan; human `yUnAy` by default, yeast `UACUAAC` optional) and `FindAcceptorBranchPoint` (best `yUnAy` 18–40 nt upstream of a given `AG`, with the polypyrimidine-tract fraction).[9][10][11]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The 3' splice-site consensus is `(Y)nNCAG|G`: an almost invariant terminal `AG` preceded by a pyrimidine-rich region (polypyrimidine tract, PPT) and followed by the first exon nucleotide.[1][2] U12-type introns have no PPT and are recognised through a highly conserved branch site close to the 3' end.[12][13][14]

### 2.2 Core Model

**Shapiro & Senapathy (1987) 3' ss score (default).** Table 1 gives the percentage of each nucleotide at 3' ss positions −14..+1 for seven GENBANK categories and their total. Methods (p. 7158): "to the left of the NCAGG consensus, the best 8 of 10 locations were used … Thus a pyrimidine sequence containing two purines can score as well as one containing no purines … any sequence containing an AG in those 10 locations [is] eliminated from consideration":

```text
score = 100 · ((t1 − l1)/(h1 − l1) + (t2 − l2)/(h2 − l2)) / 2
t1 = sum of the best 8 of the 10 observed percentages at −14..−5
l1 = sum of the lowest 8 of the 10 per-position minima;  h1 = highest 8 of the 10 per-position maxima
t2 = sum of the observed percentages at −3, −2, −1, +1 (C A G | G)
l2 / h2 = sums of the per-position minima / maxima there
```

Position −4 (`N`) is not scored. The library reports `score / 100`. All eight subtables are transcribed (`ShapiroSenapathyAcceptorTables`); the formula reproduces the paper's published 3' ss scores of real introns (Table 5) at 1 dp for 40 of 41 sites (the remaining one differs by 0.1, table rounding).

**MaxEntScan (opt-in).** The Yeo & Burge (2004) `score3ss` of the 23-nt window −20..+3 (`ScoreAcceptorMaxEnt`, bit-identical to maxentpy), in bits.[3]

**U12-type acceptors (opt-in).** Fraction of the 9 positions matching the U12 branch-site consensus `UUCCUUAAC` (branch adenosine = second `A`; Brock, Dietrich & Padgett 2008), maximised over placements whose branch `A` lies 8–21 nt upstream of the last intron nucleotide (the branch-site-to-acceptor window of the Levine & Durbin 2001 human genome U12 scan; optimum 11–13 nt, Turunen et al. 2013). Terminus `AC` (AT-AC) or `AG` (GT-AG U12).[12][13][14]

**Branch points.** Human `yUnAy` (positions −3..+1, branch `A` at 0), scored as a conservation-weighted match fraction with the Gao et al. (2008) Table 4 frequencies (`y@−3` 0.470+0.320 = 0.790, `U@−2` 0.746, `A@0` 0.923, `y@+1` 0.331+0.420 = 0.751; `−1` uninformative):

```text
score = ( 0.790·[−3 ∈ {C,U}] + 0.746·[−2 = U] + 0.923·[0 = A] + 0.751·[+1 ∈ {C,U}] ) / 3.210
```

Yeast `UACUAAC` (strictly conserved in *S. cerevisiae*, branch `A` = sixth nucleotide; Langford & Gallwitz 1983, Gao et al. 2008) is scored as the fraction of matching positions. `FindAcceptorBranchPoint` searches the 18–40 nt window upstream of the `AG` (`DistanceFromAg = index(G of AG) − index(branch A)`; Gao: 83 % of branch points at −34..−21 from the 3' end, Mercer et al. 2015) and reports the pyrimidine fraction of the PPT +4..+24 downstream of the branch point (Gao et al. 2008).[9][10]

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The S&S table for the organism's GENBANK category applies (primate default). | Scores for other taxa should use their subtable (`category`). |
| ASM-02 | U12-type acceptors are scored by branch-site consensus match only. | No trained U12-vs-U2 classifier (see §6.2). |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | S&S and U12 `Score` ∈ [0, 1] (MaxEnt in bits). | Each half of the S&S score is a (t − l)/(h − l) fraction; U12 = matches/9. |
| INV-02 | `Confidence` ∈ [0, 1]. | clamp((s − 0.5)/0.5) or 2^s/(1+2^s). |
| INV-03 | `AG` sites → `Acceptor`; U12 candidates → `U12Acceptor`. | Assigned explicitly. |
| INV-04 | `Position` = index of the last intron nucleotide (`G` of `AG`, terminal `C`/`G` for U12). | `Position = i + 1`. |
| INV-05 | Branch-point `Score` ∈ [0, 1]; perfect `yUnAy` = 1.0. | Weighted fraction. |
| INV-06 | `FindAcceptorBranchPoint` reports only branch points 18–40 nt upstream of the `AG`. | Search window. |
| INV-07 | No site is scored from a partial window or a non-ACGTU character. | Full-window rule. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA or RNA sequence (case-insensitive; T = U). | Empty / `null` → no sites. |
| `minScore` | `double` | `0.5` | Minimum score to emit a site. | [0,1] for S&S / U12, bits for MaxEnt. |
| `includeNonCanonical` | `bool` | `false` | Also report U12-type acceptors (S&S model only). | |
| `model` | `AcceptorScoringModel` | `ShapiroSenapathy` | Overload `FindAcceptorSites(seq, minScore, nc, model, category)`. | |
| `category` | `ShapiroSenapathyCategory` | `Primate` | S&S Table 1 subtable. | Ignored by MaxEnt. |

`FindBranchPoints(seq, searchStart, searchEnd, minScore[, BranchPointModel])`: motif starts in [`searchStart`, `searchEnd`] (inclusive, −1 = to end) are scored.

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Index of the last intron nucleotide (acceptor) / branch `A` position (branch point). |
| `Type` | `SpliceSiteType` | `Acceptor`, `U12Acceptor` or `Branch`. |
| `Motif` | `string` | Acceptor: −20..+3 context (RNA alphabet, truncated at the ends); branch: the scored motif. |
| `Score` | `double` | S&S score/100, MaxEnt bits, U12 match fraction, or branch-point score. |
| `Confidence` | `double` | clamp((s − 0.5)/0.5) or 2^s/(1+2^s) (MaxEnt). |

### 3.3 Preconditions and Validation

A candidate is scored only when its whole window (S&S −14..+1, MaxEnt −20..+3, U12 branch site, branch motif) lies inside the sequence and is A/C/G/T/U; an `AG` inside the S&S pyrimidine region −14..−5 eliminates the candidate (S&S 1987). Unknown enum values throw `ArgumentOutOfRangeException`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case, T → U.
2. For every `AG` (A at `i`): score the S&S window `i−12..i+2` (or the MaxEnt window `i−18..i+4`); emit `Acceptor` at `i+1` when `score ≥ minScore`.
3. With `includeNonCanonical` (S&S): for every `AC` or `AG`, score the U12 branch site over distances 8–21; emit `U12Acceptor`.

### 4.2 Decision Rules

| Candidate class | Recognition | Scoring |
|-----------------|-------------|---------|
| U2 acceptor | `AG` | S&S 3' ss score (or MaxEnt score3) |
| U12 acceptor | `AC` / `AG` (non-canonical) | best `UUCCUUAAC` match fraction, branch A 8–21 nt upstream |

### 4.3 Complexity

`O(n)` for all scanners (constant-size windows).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SpliceSitePredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `FindAcceptorSites(string, double, bool)` (S&S primate) and `FindAcceptorSites(string, double, bool, AcceptorScoringModel, ShapiroSenapathyCategory)`
- `FindBranchPoints(string, int, int, double)` (human `yUnAy`) and `FindBranchPoints(string, int, int, double, BranchPointModel)`
- `FindAcceptorBranchPoint(string, int, double)` → `BranchPointResult`
- `ScoreAcceptorMaxEnt(string)`; `CalculateMaxEntScore(string, Acceptor)` (MCP `maxent_score`)
- private: `ScoreAcceptorSiteShapiroSenapathy`, `ScoreAcceptorSiteMaxEnt`, `ScoreU12AcceptorSite`, `ScoreBranchPointConsensus`, `ComputeBranchPointPptFraction`
- MCP: `find_acceptor_sites` (`model`, `category`), `find_branch_points` (`model` human/yeast)

### 5.2 Current Behavior

Real human acceptors (Ensembl canonical transcripts of HBB, HBA1, TP53, CFTR; 40 acceptors) score primate S&S 0.717–0.992 (mean 0.88; one eliminated by an upstream AG), MaxEnt 3.21–13.37 bits; a `yUnAy` ≥ 0.5 is found 18–40 nt upstream of all 40.

### 5.3 Conformance to Theory / Spec

- S&S 3' ss score: Table 1 + Methods verbatim; Table 5 reproduced 40/41 at 1 dp.
- MaxEnt: bit-identical to maxentpy `score3`.
- Independent Python oracle (all fields, all overloads/models) vs C# on 8,104 inputs (random/IUPAC/case/U/N + 4 real genes): 0 mismatches (Evidence Dataset 6).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | "AG in those 10 locations" read as both bases inside −14..−5. | Interpretation | An AG straddling −5/−4 is not eliminated. | accepted | Table 5 sites contain none. |
| 2 | U12 acceptor = branch-site consensus fraction. | Scope | No trained U12 PWM. | declared | §6.2. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior |
|------|-------------------|
| Empty / `null` input | No sites. |
| `AG` with < 12 upstream nt or no exon nt | No S&S site (window incomplete). |
| `N` inside the scored window | Candidate not scored. |
| `AG` inside −14..−5 | Candidate eliminated (S&S 1987). |
| DNA / lower case | Same as RNA / upper case. |

### 6.2 Limitations

U12-type acceptors are scored by the `UUCCUUAAC` branch-site consensus match fraction within the 8–21 nt window, not by a trained U12 PWM log-odds classifier (Sharp & Burge 1997; Levine & Durbin 2001 used 48-intron weight matrices; intronIC matrices are GPL-3, not bundleable). The branch-point weights and windows are human (Gao 2008); use the yeast model for *S. cerevisiae*.

## 7. Examples and Related Material

- Tests: [SpliceSitePredictor_AcceptorSite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/SpliceSitePredictor_AcceptorSite_Tests.cs)
- Test specification: [SPLICE-ACCEPTOR-001.md](../../../tests/TestSpecs/SPLICE-ACCEPTOR-001.md)
- Evidence: [SPLICE-ACCEPTOR-001-Evidence.md](../../Evidence/SPLICE-ACCEPTOR-001-Evidence.md)
- Related: [Donor_Site_Detection.md](Donor_Site_Detection.md), [Gene_Structure_Prediction.md](Gene_Structure_Prediction.md)

## 8. References

1. Shapiro MB, Senapathy P. 1987. RNA splice junctions of different classes of eukaryotes: sequence statistics and functional implications in gene expression. Nucleic Acids Res 15(17):7155–7174 (Table 1, Methods p. 7158, Table 5).
2. Burge CB, Tuschl T, Sharp PA. 1999. Splicing of precursors to mRNAs by the spliceosomes. The RNA World, 2nd ed. CSHL Press.
3. Yeo G, Burge CB. 2004. Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals. J Comput Biol 11(2–3):377–394.
4. Patel AA, Steitz JA. 2003. Nat Rev Mol Cell Biol 4:960–970.
5. Hall SL, Padgett RA. 1994. J Mol Biol 239:357–365.
6. Jackson IJ. 1991. Nucleic Acids Res 19:3795–3798.
7. Dietrich RC, Incorvaia R, Padgett RA. 1997. Mol Cell 1:151–160.
8. Test specification: [SPLICE-ACCEPTOR-001.md](../../../tests/TestSpecs/SPLICE-ACCEPTOR-001.md)
9. Gao K, Masuda A, Matsuura T, Ohno K. 2008. Human branch point consensus sequence is yUnAy. Nucleic Acids Res 36(7):2257–2267 (Table 4).
10. Mercer TR et al. 2015. Genome-wide discovery of human splicing branchpoints. Genome Res 25(2):290–303.
11. Langford CJ, Gallwitz D. 1983. Evidence for an intron-contained sequence required for the splicing of yeast RNA polymerase II transcripts. Cell 33:519–527.
12. Brock JE, Dietrich RC, Padgett RA. 2008. Mutational analysis of the U12-dependent branch site consensus sequence. RNA 14:2430–2439.
13. Levine A, Durbin R. 2001. A computational scan for U12-dependent introns in the human genome sequence. Nucleic Acids Res 29:4006–4013.
14. Turunen JJ, Niemelä EH, Verma B, Frilander MJ. 2013. The significant other: splicing by the minor spliceosome. WIREs RNA 4:61–76.
