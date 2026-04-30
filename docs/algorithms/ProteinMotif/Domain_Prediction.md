# Domain Prediction & Signal Peptide Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-DOMAIN-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

This document covers two related helpers in `ProteinMotifFinder`: domain prediction by signature matching and signal peptide prediction by rule-based cleavage scoring. The domain helper searches protein sequences for a fixed set of conserved signatures associated with Zinc Finger C2H2, WD40, SH3, PDZ, and Walker A / kinase ATP-binding motifs (PROSITE PS00028; PROSITE PS00017; Pfam PF00400, PF00018, PF00595). The signal peptide helper applies von Heijne's tripartite model to the N-terminus and returns the best cleavage candidate that satisfies the rule filters from the cited literature. Both helpers are deterministic, but both are simplified relative to broader profile- or training-based annotation systems.

## 2. Scientific / Formal Basis

> A = Domain prediction, B = Signal peptide prediction

### 2.A Domain Prediction

#### Domain Context

Protein domains are recurring structural and functional units whose conserved signatures can support functional annotation when they are detected in sequence data. Signature-based domain finding uses characteristic residue patterns as a lightweight alternative to richer family models such as Pfam profile HMMs (El-Gebali et al., 2019).

#### Core Model

The documented model is exact signature matching against a fixed library of domain-associated consensus patterns:

| Domain Type | Signature Pattern | Evidence |
|-------------|-------------------|----------|
| Zinc Finger C2H2 | `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` | PROSITE PS00028; Krishna et al. (2003) |
| WD40 Repeat | `[LIVMFYWC]-x(5,12)-[WF]-D` | Current repository pattern for Pfam PF00400 |
| SH3 | `[LIVMF]-x(2)-[GA]-W-[FYW]-x(5,8)-[LIVMF]` | Current repository pattern for Pfam PF00018 |
| PDZ | `[LIVMF]-[ST]-[LIVMF]-x(2)-G-[LIVMF]-x(3,4)-[LIVMF]-x(2)-[DEN]` | Current repository pattern for Pfam PF00595 |
| Protein Kinase ATP-binding | `[AG]-x(4)-G-K-[ST]` | PROSITE PS00017; Walker et al. (1982) |

This signature library is intentionally narrower than full family models and is best read as a heuristic domain screen rather than a complete domain annotation method.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-DOMAIN-01 | The conserved residues encoded in the listed signatures are sufficient to identify the intended domain family in the scanned sequence. | True domains can be missed when they diverge from the signature, and unrelated sequences can match by chance. |
| ASM-DOMAIN-02 | Short consensus patterns are acceptable surrogates for richer family models such as profile HMMs. | Sensitivity and specificity are lower than a Pfam-style profile search. |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-DOMAIN-01 | For a fixed sequence and fixed signature library, the set of signature hits is deterministic. | Exact pattern matching has no stochastic step. |
| INV-DOMAIN-02 | Every reported hit spans a contiguous subsequence that satisfies one complete signature from the library. | Each domain call is defined by one full consensus pattern match. |
| INV-DOMAIN-03 | Signature matching is heuristic and does not guarantee absence of false positives or false negatives. | Consensus patterns capture only part of domain family variation. |

#### Comparison with Related Methods

| Aspect | Signature Matching in This Document | Pfam Profile Search |
|--------|------------------------------------|---------------------|
| Model | Short residue patterns | Profile HMM family model |
| Coverage | Limited to the listed signatures | Broad family-level coverage |
| Output | Presence of a matching signature span | Family assignment from a learned profile |

### 2.B Signal Peptide Prediction

#### Domain Context

Signal peptides are short N-terminal targeting sequences that direct proteins to secretion or membrane-insertion pathways and are cleaved after translocation (von Heijne, 1986; Owji et al., 2018). Classical descriptions divide them into an N-region enriched in positive charge, an H-region enriched in hydrophobic residues, and a C-region that contains the cleavage neighborhood (von Heijne, 1985; von Heijne, 1986).

#### Core Model

The documented scoring model follows the tripartite signal peptide description and the `-1, -3` cleavage rule from the cited literature:

| Component | Definition | Evidence |
|-----------|------------|----------|
| Cleavage rule | Positions `-1` and `-3` relative to the cleavage site are small amino acids from `{A, G, S}` | von Heijne (1983) |
| N-region score | Positive-charge density from `K` and `R`, normalized so two positive residues map to `1.0` | von Heijne (1986) |
| H-region score | Fraction of hydrophobic residues from `{A, I, L, M, F, V, W}` | von Heijne (1985) |
| C-region score | Fraction of small or polar residues from `{A, G, S, T, N}` | von Heijne (1984) |
| Total score | `((n) + 2(h) + (c)) / 4` | von Heijne (1985); current repository formulation |

The doubled H-region weight reflects the current repository interpretation of von Heijne's statement that the hydrophobic core is central to membrane targeting.

#### Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-SIGNAL-01 | The target sequence follows the classical tripartite n/h/c organization. | Atypical or weak signal peptides may not be recognized. |
| ASM-SIGNAL-02 | The `-1, -3` small-residue rule remains informative for cleavage-site selection in the scanned sequence. | Valid cleavage sites outside this rule can be rejected. |

#### Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-SIGNAL-01 | Any accepted cleavage site satisfies the `-1, -3` rule with residues from `{A, G, S}`. | The rule is a mandatory precondition of the scoring model. |
| INV-SIGNAL-02 | Accepted predictions are evaluated through explicit N-, H-, and C-region components. | The tripartite model defines the scoring structure. |
| INV-SIGNAL-03 | When component scores are normalized to `[0, 1]`, the weighted total score also lies in `[0, 1]`. | The total is a convex combination of `n`, `h`, and `c`. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[DOMAIN] proteinSequence` | `string` | required | Protein sequence scanned for the five documented domain signatures | Null or empty input yields no results; matching is case-insensitive after uppercasing |
| `[SIGNAL] proteinSequence` | `string` | required | Protein sequence inspected for an N-terminal signal peptide | Null or empty input returns `null`; sequences shorter than `15` residues return `null` |
| `[SIGNAL] maxLength` | `int` | `70` | Maximum prefix length considered during signal peptide search | The implementation searches `min(maxLength, sequence.Length)` residues |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[DOMAIN] Name` | `string` | Domain label such as `Zinc Finger C2H2` or `PDZ` |
| `[DOMAIN] Accession` | `string` | Domain accession stored by the repository for the detected signature |
| `[DOMAIN] Start` | `int` | Inclusive 0-based start coordinate of the matched signature |
| `[DOMAIN] End` | `int` | Inclusive 0-based end coordinate of the matched signature |
| `[DOMAIN] Score` | `double` | Positive score inherited from the underlying motif match |
| `[DOMAIN] Description` | `string` | Repository description for the detected domain type |
| `[SIGNAL] CleavagePosition` | `int` | Predicted cleavage index in the scanned sequence |
| `[SIGNAL] NRegion` | `string` | Returned N-region substring |
| `[SIGNAL] HRegion` | `string` | Returned H-region substring |
| `[SIGNAL] CRegion` | `string` | Returned C-region substring |
| `[SIGNAL] Score` | `double` | Weighted signal peptide score |
| `[SIGNAL] Probability` | `double` | Repository probability field associated with the returned score |

### 3.3 Preconditions and Validation

`FindDomains(...)` returns an empty sequence for null or empty input. `PredictSignalPeptide(...)` returns `null` for null, empty, or shorter-than-`15` input. Both helpers uppercase the protein sequence before matching. Domain coordinates are reported as inclusive 0-based indexes. The repository does not perform full amino-acid alphabet validation beyond the pattern-matching and residue-set checks used internally by the two methods.

## 4. Algorithm

### 4.A Domain Prediction

#### High-Level Steps

1. Convert the input protein sequence to uppercase.
2. Evaluate each hard-coded domain signature against the full sequence.
3. Wrap each signature hit as a `ProteinDomain` with its repository metadata.
4. Return the union of all domain hits.

#### Decision Rules / Reference Tables

| Domain | Accession | Matching Rule |
|--------|-----------|---------------|
| Zinc Finger C2H2 | `PF00096` | `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` |
| WD40 Repeat | `PF00400` | `[LIVMFYWC].{5,12}[WF]D` |
| SH3 | `PF00018` | `[LIVMF].{2}[GA]W[FYW].{5,8}[LIVMF]` |
| PDZ | `PF00595` | `[LIVMF][ST][LIVMF].{2}G[LIVMF].{3,4}[LIVMF].{2}[DEN]` |
| Protein Kinase ATP-binding | `PF00069` | `[AG].{4}GK[ST]` |

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindDomains(...)` | `O(n x d)` | `O(1)` plus output | `n` = sequence length, `d = 5` hard-coded signatures |

### 4.B Signal Peptide Prediction

#### High-Level Steps

1. Consider cleavage candidates in the N-terminal search window.
2. Reject candidates that violate the `-1, -3` small-residue rule.
3. Derive N-, H-, and C-region substrings for each surviving candidate.
4. Score the three regions and combine them with the weighted formula `((n) + 2(h) + (c)) / 4`.
5. Return the highest-scoring candidate, or no prediction if no candidate passes the rule filters.

#### Decision Rules / Reference Tables

| Rule | Value |
|------|-------|
| `-1, -3` residues | `{A, G, S}` |
| N-region charge residues | `{K, R}` |
| H-region hydrophobic residues | `{A, I, L, M, F, V, W}` |
| C-region small or polar residues | `{A, G, S, T, N}` |
| Site ranking | Highest total score among candidates that satisfy the rule filters |

#### Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictSignalPeptide(...)` | `O(s)` | `O(1)` | `s` = number of cleavage candidates considered in the searched N-terminal window |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.FindDomains(string)`: scans the sequence against five hard-coded domain regexes and emits `ProteinDomain` values.
- `ProteinMotifFinder.PredictSignalPeptide(string, int)`: evaluates cleavage candidates in an N-terminal window and returns the best `SignalPeptide?` result.
- `ProteinMotifFinder.ScoreNRegion(string)`, `ProteinMotifFinder.ScoreHydrophobicRegion(string)`, and `ProteinMotifFinder.ScoreCRegion(string)`: compute the three signal peptide component scores used by `PredictSignalPeptide(...)`.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- `FindDomains(...)` delegates each domain signature to `FindMotifByPattern(...)`, so `ProteinDomain.Score` is the same information-content score used by motif matches.
- `FindDomains(...)` currently emits only `Zinc Finger C2H2`, `WD40 Repeat`, `SH3`, `PDZ`, and `Protein Kinase ATP-binding` hits, with accessions `PF00096`, `PF00400`, `PF00018`, `PF00595`, and `PF00069`.
- Returned domain coordinates are inclusive 0-based indexes inherited from `MotifMatch.Start` and `MotifMatch.End`.
- `PredictSignalPeptide(...)` rejects null, empty, and shorter-than-`15` input and limits the search to `Math.Min(maxLength, sequence.Length)` residues.
- Candidate cleavage sites are scanned from positions `15` through `35`; candidates must satisfy the `-1, -3` rule, have `HRegion.Length >= 7`, have positive N-region score, and have `hScore >= 0.5`.
- The highest-scoring candidate is returned when any site survives those filters. `Probability` is currently set equal to `Score`, and the source does not apply a separate total-score threshold before returning a prediction.

### 5.3 Conformance to Theory / Spec

#### 5.3.A Domain Prediction

**Implemented (verbatim from the cited theory/spec):**

- Signature matching for the documented Zinc Finger C2H2, WD40, SH3, PDZ, and Walker A / kinase ATP-binding patterns.
- Deterministic reporting of the full matched signature span for each detected pattern.

**Intentionally simplified:**

- Regex signatures are used instead of Pfam profile HMMs; **consequence:** remote homologs and family context outside the short signature are not modeled.
- Only five domain families are included in source; **consequence:** other domain classes are silently outside scope.
- Reported scores come from the repository's generic motif-scoring helper rather than a domain-specific confidence model; **consequence:** scores are useful as internal match-strength indicators, not calibrated domain probabilities.

**Not implemented:**

- Full Pfam or InterPro profile searches; **users should rely on:** external domain annotation tools.
- Domain-boundary refinement beyond the matched signature span; **users should rely on:** no current alternative in this repository.

#### 5.3.B Signal Peptide Prediction

**Implemented (verbatim from the cited theory/spec):**

- The `-1, -3` small-residue rule using `{A, G, S}`.
- Tripartite N-, H-, and C-region scoring with the weighted total `((n) + 2(h) + (c)) / 4`.
- Selection of the best-scoring cleavage candidate that satisfies the repository's rule filters.

**Intentionally simplified:**

- Search is restricted to a fixed N-terminal window and fixed region extraction rules; **consequence:** edge cases outside that window can differ from broader signal peptide predictors.
- `Probability` is set equal to `Score`; **consequence:** the reported probability is an internal quality field, not a separately calibrated estimate.
- The method requires positive N-region charge and `hScore >= 0.5`; **consequence:** weak or atypical signal peptides can be missed.

**Not implemented:**

- Organism-specific or machine-learned signal peptide models; **users should rely on:** external signal peptide predictors.
- Alternative ranked candidate output beyond the single best site; **users should rely on:** no current alternative in this repository.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Domain signatures are used as surrogates for richer family models | Assumption | Some true domains can be missed and some chance matches can be reported | accepted | This is the core simplification behind `FindDomains(...)` |
| 2 | Signal peptide search is limited to cleavage positions `15` through `35` within `maxLength` | Assumption | Informative cleavage sites outside that window are never returned | accepted | This bound is explicit in `PredictSignalPeptide(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `[DOMAIN]` Null or empty sequence | Returns no domains | `FindDomains(...)` yields no results for null or empty input |
| `[SIGNAL]` Null or empty sequence | Returns `null` | `PredictSignalPeptide(...)` rejects null or empty input before scanning |
| `[SIGNAL]` Sequence shorter than `15` residues | Returns `null` | The method enforces a minimum length before evaluating cleavage candidates |
| Lowercase or mixed-case protein sequence | Same result as uppercase input | Both helpers uppercase the sequence before matching |
| `[DOMAIN]` Sequence with multiple independent signatures | Returns one `ProteinDomain` per signature hit | Each signature is scanned independently |
| `[SIGNAL]` Sequence with no candidate satisfying the rule filters | Returns `null` | No cleavage candidate survives the rule checks |

### 6.2 Limitations

The domain helper is a compact signature scanner, not a full family annotation engine. The signal peptide helper is a rule-based predictor bound to a fixed N-terminal search strategy. Neither helper models broader sequence context, organism-specific effects, or external calibration data, so results should be interpreted as lightweight annotations rather than full proteome-scale annotation outputs.

## 8. References

1. von Heijne G. Signal sequences. The limits of variation. Journal of Molecular Biology. https://doi.org/10.1016/0022-2836(85)90046-4
2. von Heijne G. A new method for predicting signal sequence cleavage sites. Nucleic Acids Research. https://doi.org/10.1093/nar/14.11.4683
3. von Heijne G. How signal sequences maintain cleavage specificity. Journal of Molecular Biology. https://doi.org/10.1016/0022-2836(84)90192-X
4. von Heijne G. Patterns of amino acids near signal-sequence cleavage sites. European Journal of Biochemistry. https://doi.org/10.1111/j.1432-1033.1983.tb07624.x
5. Owji H, Nezafat N, Negahdaripour M, Hajiebrahimi A, Ghasemi Y. A comprehensive review of signal peptides. European Journal of Cell Biology. https://doi.org/10.1016/j.ejcb.2018.06.003
6. Walker JE, Saraste M, Runswick MJ, Gay NJ. Distantly related sequences in the alpha- and beta-subunits of ATP synthase, myosin, kinases and other ATP-requiring enzymes and a common nucleotide binding fold. The EMBO Journal. https://doi.org/10.1002/j.1460-2075.1982.tb01276.x
7. Krishna SS, Majumdar I, Grishin NV. Structural classification of zinc fingers. Nucleic Acids Research. https://doi.org/10.1093/nar/gkg161
8. Hulo N, Bairoch A, Bulliard V, et al. The PROSITE database. Nucleic Acids Research. https://doi.org/10.1093/nar/gkj063
9. El-Gebali S, Mistry J, Bateman A, et al. The Pfam protein families database in 2019. Nucleic Acids Research. https://doi.org/10.1093/nar/gky995
10. PROSITE PS00028. Zinc finger C2H2 type. https://prosite.expasy.org/PS00028
11. PROSITE PS00017. ATP/GTP-binding site motif A. https://prosite.expasy.org/PS00017
