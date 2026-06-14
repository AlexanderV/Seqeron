# Fusion Gene Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-FUSION-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

A gene fusion joins two genes via a breakpoint; in RNA-seq it is supported by two read classes: split
(junction-spanning) reads and discordant (spanning) mate-pairs [1][4]. This algorithm performs rule-based
candidate fusion calling: given per-class supporting-read counts for each candidate breakpoint, it applies
the STAR-Fusion minimum-support rule to decide which candidates are reported, computes total support, and
classifies the junction as in-frame or out-of-frame from codon phase [2][5]. It is a deterministic,
threshold-driven (heuristic) detector, not a probabilistic model. Use it as the first stage of a fusion
pipeline once chimeric/discordant evidence has been grouped per breakpoint.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A fusion transcript is produced when sequence from one gene is joined to sequence from another gene at a
breakpoint. RNA-seq evidence for a fusion is of two kinds: **split reads** — "reads with two segments
aligning in a noncontiguous fashion" — and **discordant mates** (a.k.a. spanning/bridge reads) —
"paired-end reads originating from the same fragment but with the mates aligning in a nonlinear way" [4].
The gene to which the longer segment of a split read aligns is its **anchor** [5].

### 2.2 Core Model

For a candidate with split-read counts anchored in the 5' and 3' partners (`split_reads1`, `split_reads2`)
and discordant mate-pairs (`discordant_mates`):

- **Total support** = `split_reads1 + split_reads2 + discordant_mates` [5].
- **Junction reads** = `split_reads1 + split_reads2` (the split/junction-spanning class) [2].
- **Detection rule** (STAR-Fusion defaults) [3]:
  - if `junction_reads ≥ 1`: report iff `junction_reads ≥ MIN_JUNCTION_READS (=1)` **and**
    `total_support ≥ MIN_SUM_FRAGS (=2)`;
  - if `junction_reads = 0`: report iff `discordant_mates ≥ MIN_SPANNING_FRAGS_ONLY (=5)`.
- **Reading frame** [6][7]: the junction is in-frame iff
  `(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0`, i.e. the 3' partner stays in codon phase
  across the breakpoint (codons are read in triplets; three reading frames exist, selected modulo 3).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Reported fusions have gene5p ≠ gene3p | A gene is not fused with itself (Registry invariant; same-gene candidates are skipped) |
| INV-02 | TotalSupport = split_reads1 + split_reads2 + discordant_mates | Arriba total-support definition [5] |
| INV-03 | Reported ⇒ (junction ≥ MIN_JUNCTION_READS ∧ total ≥ MIN_SUM_FRAGS) ∨ (junction = 0 ∧ discordant ≥ MIN_SPANNING_FRAGS_ONLY) | STAR-Fusion thresholds [3] |
| INV-04 | Results ordered by descending TotalSupport | STAR-Fusion scores by abundance of supporting reads [1][2] |
| INV-05 | InFrame ⇔ (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0 | Codon-phase / modulo-3 reading frame [6][7] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | This unit (rule-based count threshold) | Full STAR-Fusion / Arriba |
|--------|----------------------------------------|----------------------------|
| Input | Pre-grouped per-candidate read counts | Raw chimeric alignments from a BAM |
| FFPM/abundance normalization | Not applied (count thresholds only) | Applies min_FFPM and artifact filters |
| Frame call | Codon phase only | Phase + premature-stop-codon scan |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| candidates | IEnumerable&lt;FusionCandidate&gt; | required | Candidate breakpoints with per-class read counts | non-null; non-negative counts |
| thresholds | FusionDetectionThresholds? | STAR-Fusion defaults | min junction / sum / spanning-only thresholds | — |
| fivePrimeCodingBases | int | -1 (unknown) | Coding bases of 5' partner upstream of breakpoint | ≥ 0 for IsInFrame |
| threePrimeStartPhase | int | -1 (unknown) | Coding-start phase of 3' partner | 0, 1, or 2 for IsInFrame |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | IReadOnlyList&lt;FusionCall&gt; | Detected fusions, descending TotalSupport |
| FusionCall.JunctionReads | int | split_reads1 + split_reads2 |
| FusionCall.TotalSupport | int | split_reads1 + split_reads2 + discordant_mates |
| FusionCall.ReadingFrame | enum | InFrame / OutOfFrame / Unknown |

### 3.3 Preconditions and Validation

`DetectFusions(null, …)` throws `ArgumentNullException`. A candidate with any negative supporting-read
count throws `ArgumentException`. Empty input returns an empty list. Same-gene candidates (gene5p == gene3p,
case-insensitive) are silently skipped (not a fusion). `IsInFrame` throws `ArgumentOutOfRangeException` for a
negative base count or a phase outside {0,1,2}. When coding-phase fields are absent (-1), the reading frame
is reported as `Unknown` rather than guessed.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input (non-null; non-negative counts).
2. For each candidate: skip if gene5p == gene3p.
3. Compute junction reads and total support.
4. Apply the STAR-Fusion threshold rule (junction-present vs spanning-only branch).
5. For passing candidates, resolve reading frame from codon phase.
6. Return calls ordered by descending total support (gene pair as tie-break).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

| Constant | Value | Source |
|----------|-------|--------|
| MIN_JUNCTION_READS | 1 | STAR-Fusion source [3] |
| MIN_SUM_FRAGS | 2 | STAR-Fusion source [3] |
| MIN_SPANNING_FRAGS_ONLY | 5 | STAR-Fusion source [3] |
| CodonLength | 3 | Reading frame read in triplets [7] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectFusions | O(n log n) | O(n) | one pass over n candidates + final sort by support |
| IsInFrame / ComputeTotalSupport | O(1) | O(1) | arithmetic |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.DetectFusions(IEnumerable<FusionCandidate>, FusionDetectionThresholds?)`: rule-based detection.
- `OncologyAnalyzer.IsInFrame(int, int)`: codon-phase in-frame test.
- `OncologyAnalyzer.ComputeTotalSupport(FusionCandidate)`: split1 + split2 + discordant.

### 5.2 Current Behavior

Operates on candidate-level read counts (Arriba output schema), not raw BAM records: chimeric-read
extraction (`FindChimericReads`) and reference-genome validation (`ValidateFusion`) are separate
out-of-scope steps. No suffix tree is used: this unit performs no substring search or occurrence
enumeration — it applies arithmetic count thresholds and a single sort — so the repository suffix tree does
not apply. Reading-frame `Unknown` is returned when phase fields are unset, never guessed.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Total support = split_reads1 + split_reads2 + discordant_mates [5].
- STAR-Fusion thresholds: MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5 [3].
- In-frame iff codon phase preserved across the junction (modulo 3) [6][7].
- gene5p ≠ gene3p (a gene is not fused with itself).

**Intentionally simplified:**

- Frame call: codon phase only; **consequence:** a junction whose frame is numerically preserved but which
  has a premature stop codon (Arriba "stop-codon" value) is still reported InFrame here.
- Abundance normalization: STAR-Fusion's min_FFPM and artifact filters are not applied; **consequence:**
  count thresholds alone gate detection.

**Not implemented:**

- Chimeric-read extraction from BAM (`FindChimericReads`); **users should rely on:** an upstream aligner
  (STAR/Arriba) to produce per-candidate counts. (ONCO-FUSION-001 read-extraction row.)
- Known-fusion database lookup and breakpoint/protein analysis; **users should rely on:** ONCO-FUSION-002 /
  ONCO-FUSION-003 (future units).

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Input is per-candidate counts, not raw reads | Assumption | Limits unit to post-grouping stage | accepted | Mirrors Arriba output schema [5] |
| 2 | No premature-stop-codon detection | Assumption | InFrame may overstate translatability | accepted | Deferred to ONCO-FUSION-003 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Spanning-only, 4 discordant, 0 junction | Rejected | 4 < MIN_SPANNING_FRAGS_ONLY (5) [3] |
| Spanning-only, 5 discordant, 0 junction | Detected | 5 ≥ 5 [3] |
| Junction=1, total=1 | Rejected | total < MIN_SUM_FRAGS (2) [3] |
| gene5p == gene3p | Skipped | Not a fusion (INV-01) |
| Empty input | Empty list | trivial |
| Phase fields unset (-1) | ReadingFrame = Unknown | No guessing |

### 6.2 Limitations

Prokaryote/eukaryote agnostic on counts; does not detect read-through false positives beyond the
distinct-gene rule, does not normalize by library size (FFPM), and does not reconstruct the fusion
transcript (so premature stop codons are not detected). Reciprocal fusions are treated as distinct
candidates (different 5'/3' assignment).

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var candidates = new[]
{
    new OncologyAnalyzer.FusionCandidate("EML4", "ALK", 3, 2, 4, 300, 0), // junction=5, total=9, in-frame
    new OncologyAnalyzer.FusionCandidate("CD74", "ROS1", 0, 0, 5),         // spanning-only, total=5 → detected
    new OncologyAnalyzer.FusionCandidate("NCOA4", "RET", 0, 0, 4),         // spanning-only, 4 < 5 → rejected
};

IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);
// calls[0] = EML4-ALK, TotalSupport 9, ReadingFrame InFrame; CD74-ROS1 also present; NCOA4-RET absent.
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectFusions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectFusions_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ONCO-FUSION-001-Evidence.md](../../../docs/Evidence/ONCO-FUSION-001-Evidence.md)

## 8. References

1. Haas BJ, Dobin A, Li B, Stransky N, Pochet N, Regev A. 2019. Accuracy assessment of fusion transcript detection via read-mapping and de novo fusion transcript assembly-based methods. *Genome Biology* 20:213. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9
2. Haas BJ, Dobin A, Stransky N, et al. 2017. STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq. *bioRxiv* 120295. https://www.biorxiv.org/content/10.1101/120295
3. STAR-Fusion source (MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5). https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
4. Uhrig S, Ellermann J, Walther T, et al. 2021. Accurate and efficient detection of gene fusions from RNA sequencing data. *Genome Research* 31(3):448–460. https://pmc.ncbi.nlm.nih.gov/articles/PMC7919457/
5. Arriba output-file documentation (split_reads1/split_reads2/discordant_mates; reading_frame). https://github.com/suhrig/arriba/wiki/05-Output-files
6. Genomics England. 2021. Improving how we report gene fusion productivity. https://www.genomicsengland.co.uk/blog/gene-fusion-reporting
7. Wikipedia. Reading frame (citing Badger JH, Olsen GJ. 1999. *Mol Biol Evol* 16(4):512–524; Lodish, *Molecular Cell Biology* 6th ed., p.121). https://en.wikipedia.org/wiki/Reading_frame
