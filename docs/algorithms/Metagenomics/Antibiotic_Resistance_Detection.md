# Antibiotic Resistance Gene Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-RESIST-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Detects acquired antibiotic-resistance genes in assembled contigs by screening each contig
against a caller-supplied reference database of resistance genes, following the ResFinder
methodology [1]. For each reference gene the best ungapped alignment to the contig is located,
its BLAST-style percent identity [5] and its coverage of the reference gene length [1] are
computed, and the reference gene is reported only when both pass user-selectable thresholds.
The single best-matching reference gene per contig is returned, mirroring ResFinder's
"best-matching gene" output [1]. It is specification-driven (the detection rule is fixed by the
ResFinder definition); the gene catalogue is supplied by the caller because curated AMR databases
(ResFinder, CARD) are large tables that cannot be embedded verbatim [6].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Acquired antimicrobial-resistance (AMR) genes are horizontally transferred genes whose presence
in a genome/metagenome confers resistance to an antibiotic class. ResFinder identifies them by
sequence-aligning a curated database of known resistance genes against assembled contigs and
reporting database genes that match closely enough over enough of their length [1].

### 2.2 Core Model

For a contig *C* and a reference gene *R* of length *m*:

- **Percent identity** (BLAST definition [5]): the number of identical positions divided by the
  number of alignment columns. For a gapless (ungapped) alignment there are no gap columns, so the
  denominator equals the aligned window length *w*: `identity = matches / w` [5]. ResFinder defines
  %ID as "the percentage of nucleotides that are identical between the best-matching resistance
  gene in the database and the corresponding sequence in the genome" [1].
- **Coverage** (relative to the *reference* gene): `coverage = w / m`, the fraction of the reference
  gene length spanned by the alignment [1][3].
- **Reporting rule:** *R* is reported for *C* only if `identity ≥ idThreshold` and
  `coverage ≥ covThreshold` [1]. ResFinder requires a gene to "cover at least 2/5 of the length of
  the resistance gene in the database" [1]; the 60% coverage floor in later operation exists so
  "genes lying on the edge of a contig or spread over two contigs are not missed" [3].
- **Best-matching gene:** when several references pass, only the best match per contig is reported
  (highest identity, ties broken by greater coverage) [1][6].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ identity ≤ 1 and 0 ≤ coverage ≤ 1 | matches ≤ window ≤ m [5][1] |
| INV-02 | A hit is reported only if identity ≥ idThreshold AND coverage ≥ covThreshold | reporting rule [1] |
| INV-03 | At most one hit per contig (best match: max identity, tie → max coverage) | ResFinder best-matching gene [1][6] |
| INV-04 | Exact full-length match ⇒ identity = 1.0, coverage = 1.0 | CARD "Perfect" match [6]; identity formula [5] |
| INV-05 | Default thresholds = 0.90 identity, 0.60 coverage | ResFinder web-service / README defaults [2][3][4] |

### 2.5 Comparison with Related Methods

| Aspect | This detector | ResFinder (full) |
|--------|---------------|------------------|
| Alignment | gapless ungapped | gapped BLAST |
| Gene database | caller-supplied | bundled curated DB |
| Identity definition | matches / aligned columns [5] | same [1][5] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| contigs | `IEnumerable<(string ContigId, string Sequence)>` | required | Assembled contigs to screen | non-null; empty sequences skipped |
| referenceGenes | `IEnumerable<(string GeneId, string Sequence, string Name, string AntibioticClass)>` | required | Caller-supplied resistance-gene DB | non-null; empty sequences ignored |
| identityThreshold | `double` | 0.90 | Minimum percent identity | [0, 1] |
| coverageThreshold | `double` | 0.60 | Minimum reference coverage | [0, 1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| ContigId | string | Query contig identifier |
| ResistanceGene | string | Name of the matched reference gene |
| AntibioticClass | string | Antibiotic class of the matched gene |
| PercentIdentity | double | identical positions / gapless alignment length (0–1) |
| Coverage | double | aligned window / reference length (0–1) |

### 3.3 Preconditions and Validation

Null `contigs` or `referenceGenes` → `ArgumentNullException`. A threshold outside [0, 1] →
`ArgumentOutOfRangeException`. Empty contig sequences and empty-sequence reference genes are
skipped. Comparison is case-sensitive, nucleotide (caller normalizes case/alphabet); no T↔U
conversion is performed.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs and thresholds; materialize reference genes with non-empty sequences.
2. For each contig with a non-empty sequence, for each reference gene compute the best ungapped
   match (identity, coverage) via `BestUngappedMatch`.
3. Discard reference genes failing either threshold.
4. Keep the best-matching reference gene per contig (max identity, tie → max coverage).
5. Yield one `ResistanceHit` per contig that has a passing best match.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

`BestUngappedMatch` slides the reference across the contig at every offset from `-(m-1)` to
`n-1` (overhanging both ends so contig-edge / truncated genes are scored against the reference
length [3]); for each offset it counts identical positions over the overlapping window and keeps
the offset with the most matches (ties → **shorter**, i.e. higher-identity window, so the chosen
alignment is never padded with flanking mismatches; padding would dilute identity and could
spuriously fail the identity threshold even when a perfect HSP exists — mirroring BLAST reporting
the best-scoring HSP [5]). Thresholds (0.90 / 0.60) are named constants citing ResFinder [1][2][3][4].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindAntibioticResistanceGenes | O(c · d · n · m) | O(d) | c contigs, d reference genes; per pair the sliding ungapped scan is O(n·m) |
| BestUngappedMatch | O(n · m) | O(1) | n = contig length, m = reference length |

Per the checklist this unit is rated O(n × d); the inner ungapped scan adds the n·m alignment factor.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, referenceGenes, identityThreshold, coverageThreshold)`: best-match resistance-gene detector.
- `MetagenomicsAnalyzer.BestUngappedMatch(contig, reference)` (private): gapless identity/coverage computation.
- `MetagenomicsAnalyzer.ResistanceHit`: result record.
- `MetagenomicsAnalyzer.DefaultResistanceIdentityThreshold` / `DefaultResistanceCoverageThreshold`: constants.

### 5.2 Current Behavior

Uses a direct gapless sliding scan rather than the repository suffix tree: matching here is
approximate (mismatch-tolerant) and scoring-based (identity over a window), not exact-substring
occurrence enumeration, so the suffix tree's O(m) exact-match lookups do not fit. The legacy
`FindResistanceGenes` (motif-containment stub) is retained separately for the existing MCP tool
and is not part of this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- BLAST percent identity = matches / alignment columns (gapless case) [5].
- Coverage = aligned length / reference gene length [1][3].
- Dual-threshold reporting rule and best-matching-gene single output [1][6].
- ResFinder default thresholds 0.90 / 0.60 as named constants [2][3][4].

**Intentionally simplified:**

- Alignment model: gapless ungapped sliding match instead of full gapped BLAST; **consequence:**
  genes whose true alignment requires insertions/deletions may score lower identity/coverage than
  gapped BLAST would report. Substitution divergence and contig-edge truncation are scored exactly.

**Not implemented:**

- Bundled curated gene catalogue; **users should rely on:** supplying a ResFinder/CARD-derived
  reference set via `referenceGenes` (curated tables are not fabricated [6]).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Gapless alignment | Assumption | indel-containing matches under-scored vs gapped BLAST | accepted | See 5.3 "Intentionally simplified" |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty contig sequence | skipped (no hit) | nothing to align |
| Reference gene with empty sequence | ignored | no length to cover |
| No reference passes thresholds | contig produces no hit | reporting rule [1] |
| Multiple passing references | only best match returned | best-matching gene [1] |
| null contigs/referenceGenes | ArgumentNullException | contract |
| threshold ∉ [0,1] | ArgumentOutOfRangeException | contract |

### 6.2 Limitations

No gapped alignment (no indel handling); no SNP/protein homology models (CARD's homolog/variant
models are out of scope); detection quality depends entirely on the caller-supplied database;
nucleotide-only, case-sensitive comparison.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var contigs = new[] { ("c1", "AAACGTACGT") };
var db = new[] { ("blaX", "CGTACGT", "blaX-like", "beta-lactam") };
var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).ToList();
// hits[0]: ContigId="c1", ResistanceGene="blaX-like", PercentIdentity=1.0, Coverage=1.0
```

**Numerical walk-through:** reference `CGTACGT` (m=7) vs contig `CGTTCGT`: best ungapped offset
aligns all 7 positions, 6 identical (position 4 T vs A) → identity = 6/7 ≈ 0.857, coverage = 7/7 = 1.0.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Metagenomics/MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [META-RESIST-001-Evidence.md](../../../docs/Evidence/META-RESIST-001-Evidence.md)

## 8. References

1. Zankari E, Hasman H, Cosentino S, et al. 2012. Identification of acquired antimicrobial resistance genes. J Antimicrob Chemother 67(11):2640–2644. https://academic.oup.com/jac/article/67/11/2640/707208
2. genomicepidemiology/resfinder. ResFinder reference implementation (README threshold defaults). https://github.com/genomicepidemiology/resfinder
3. 2023. Pipeline validation for the identification of antimicrobial-resistant genes in carbapenem-resistant Klebsiella pneumoniae. Sci Rep 13. https://www.nature.com/articles/s41598-023-42154-6
4. 2016. Benchmarking of methods for identification of antimicrobial resistance genes in bacterial whole genome data. J Antimicrob Chemother 71(9):2484–2492. https://academic.oup.com/jac/article/71/9/2484/2238319
5. Li H. 2018. On the definition of sequence identity. https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
6. Alcock BP, et al. CARD Resistance Gene Identifier (RGI), McMaster University. https://card.mcmaster.ca/analyze/rgi
