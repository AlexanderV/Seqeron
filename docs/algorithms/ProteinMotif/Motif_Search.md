# Protein Motif Search

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-FIND-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Protein motif search in this repository finds short sequence patterns in protein strings by applying regular expressions to the input sequence. The repository exposes a direct single-pattern entry point and a multi-pattern helper that scans a fixed `CommonMotifs` dictionary containing both PROSITE-sourced motifs and several literature-based non-PROSITE motifs. The search is deterministic and case-insensitive, and the current implementation reports overlapping occurrences when a pattern can match at adjacent offsets. Match scores and E-values are repository-defined helpers rather than PROSITE or ScanProsite significance values.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Protein motifs are short conserved sequence patterns associated with functions such as post-translational modification, ligand binding, or localization. PROSITE is a standard curated collection of such patterns and defines a formal pattern language for motif descriptions (PROSITE Database; PROSITE User Manual; Hulo et al., 2007). Pattern-based motif search asks whether a sequence contains a contiguous subsequence satisfying a specified motif pattern.

### 2.2 Core Model

The formal model for this document is exact pattern occurrence search: given a protein sequence `S` and a motif pattern `P`, report every subsequence of `S` that satisfies `P`. In the repository, `P` is a regular expression supplied directly to `FindMotifByPattern(...)`, while `FindCommonMotifs(...)` repeats the same search across every stored `PrositePattern` in `CommonMotifs`. This model is appropriate for motifs whose defining evidence is a short residue pattern rather than a profile, HMM, or structural context model.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The biological site of interest is sufficiently characterized by a short local residue pattern. | Functional sites that depend on broader context can be missed or overcalled. |
| ASM-02 | A pattern match is treated as motif evidence without requiring structural or evolutionary confirmation. | False positives remain possible even when the sequence satisfies the pattern exactly. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | For a fixed sequence and fixed regex pattern, the set of reported matches is deterministic. | The scan contains no random or data-dependent sampling step. |
| INV-02 | Every reported match corresponds to a contiguous subsequence of the input sequence. | A motif hit is defined by one regex capture span. |
| INV-03 | Pattern search does not by itself prove biological activity. | Sequence-pattern matching is necessary evidence, not a full structural or functional model. |

### 2.5 Comparison with Related Methods

| Aspect | Regex Motif Search in This Document | Profile / HMM Search |
|--------|------------------------------------|----------------------|
| Pattern representation | Explicit residue pattern | Position-specific family model |
| Evidence required | Exact local motif occurrence | Statistical match to a trained profile |
| Scope | Best for short signature motifs | Better for broader or more variable families |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `proteinSequence` | `string` | required | Protein sequence scanned for motif hits | Null or empty input yields no matches |
| `regexPattern` | `string` | required for `FindMotifByPattern(...)` | Regular expression used for motif search | Null, empty, or invalid regex yields no matches |
| `motifName` | `string` | `Custom` | Name stored in each returned `MotifMatch` | Caller-supplied for direct scans; repository-supplied for `FindCommonMotifs(...)` |
| `patternId` | `string` | `""` | Identifier stored in `MotifMatch.Pattern` | Caller-supplied for direct scans; `FindCommonMotifs(...)` passes the dictionary accession |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Inclusive 0-based start index of the matched substring |
| `End` | `int` | Inclusive 0-based end index of the matched substring |
| `Sequence` | `string` | Uppercased matched substring |
| `MotifName` | `string` | Name associated with the motif entry |
| `Pattern` | `string` | Stored pattern identifier for the current call |
| `Score` | `double` | Repository motif score computed from regex information content |
| `EValue` | `double` | Repository E-value estimate derived from match length, sequence length, and score |

### 3.3 Preconditions and Validation

`FindMotifByPattern(...)` and `FindCommonMotifs(...)` return an empty sequence for null or empty `proteinSequence`. `FindMotifByPattern(...)` also returns an empty sequence if `regexPattern` is null, empty, or fails regex compilation. Matching is case-insensitive because the input is uppercased and the regex is compiled with `RegexOptions.IgnoreCase`. The repository reports inclusive 0-based coordinates and does not perform separate validation of the amino-acid alphabet beyond regex evaluation.

## 4. Algorithm

### 4.1 High-Level Steps

1. Uppercase the input protein sequence.
2. Compile the supplied regex pattern inside a lookahead wrapper so the scan can start at every position.
3. Enumerate every captured match span in the sequence.
4. For each hit, compute the repository score and E-value, then return a `MotifMatch`.
5. For `FindCommonMotifs(...)`, repeat the same scan for each stored motif entry in `CommonMotifs`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Item | Rule |
|------|------|
| Regex wrapper | `(?=(pattern))` so overlapping occurrences can be discovered |
| Score | Sum of per-position information content terms from the regex allowed-count model |
| E-value | `(N - L + 1) x 2^(-Score)` with `N` = sequence length and `L` = match length |
| Pattern library | `CommonMotifs` maps each accession or mnemonic ID to `Accession`, `Name`, `Pattern`, `RegexPattern`, and `Description` |
| Current motif families | PROSITE entries `PS00001`, `PS00004`, `PS00005`, `PS00006`, `PS00007`, `PS00008`, `PS00009`, `PS00016`, `PS00017`, `PS00018`, `PS00028`, `PS00029`, plus literature-based `NLS1`, `NES1`, `SIM1`, `WW1`, and `SH3_1` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindMotifByPattern(...)` | Regex scan over `n` residues | `O(k)` | `k` = number of returned matches; actual scan cost depends on the supplied pattern |
| `FindCommonMotifs(...)` | `O(m x scan)` | `O(k)` | `m` = number of entries in `CommonMotifs`; each entry delegates to `FindMotifByPattern(...)` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.FindMotifByPattern(string, string, string, string)`: scans one regex pattern and returns `MotifMatch` records.
- `ProteinMotifFinder.FindCommonMotifs(string)`: iterates the `CommonMotifs` dictionary and delegates each entry to `FindMotifByPattern(...)`.
- `ProteinMotifFinder.CalculateMotifScore(string, string)`, `ProteinMotifFinder.CalculateEValue(int, int, double)`, and `ProteinMotifFinder.ParseRegexAllowedCounts(string)`: implement the repository's score and E-value calculations.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- `FindMotifByPattern(...)` wraps the supplied regex in a lookahead, so overlapping occurrences are returned when the pattern allows them.
- Returned coordinates are inclusive 0-based indexes, and `Sequence` is the exact uppercased substring at `[Start..End]`.
- Invalid regex input is swallowed and yields no matches instead of throwing.
- `FindCommonMotifs(...)` scans the repository's fixed in-source `CommonMotifs` dictionary, which includes both PROSITE entries and the non-PROSITE entries `NLS1`, `NES1`, `SIM1`, `WW1`, and `SH3_1`.
- `Score` is the total information content derived from the regex pattern, and `EValue` is computed from that score under a uniform 20-amino-acid alphabet assumption.
- Current tests explicitly validate that the stored `PS00007` and `PS00018` definitions match the official PROSITE entries.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact motif occurrence search for supplied regex patterns.
- Deterministic reporting of match coordinates and matched substrings.
- Storage of official PROSITE pattern definitions for the current PROSITE-backed entries in `CommonMotifs`.

**Intentionally simplified:**

- `FindCommonMotifs(...)` scans only the fixed in-repository dictionary; **consequence:** the method is not a full PROSITE database search.
- Five non-PROSITE literature motifs are stored alongside PROSITE entries; **consequence:** not every returned hit corresponds to a PROSITE accession.
- `Score` and `EValue` are derived from regex allowed counts and a uniform amino-acid model; **consequence:** they are repository heuristics, not ScanProsite significance values.

**Not implemented:**

- Live synchronization with the full PROSITE catalog; **users should rely on:** external PROSITE or ScanProsite resources.
- Profile- or HMM-based motif scoring; **users should rely on:** no current alternative in this repository.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `CommonMotifs` mixes PROSITE and non-PROSITE entries | Assumption | Output from `FindCommonMotifs(...)` is broader than a pure PROSITE-only scan | accepted | Non-PROSITE entries are explicit in source and tests |
| 2 | Match scores and E-values assume a uniform 20-amino-acid background | Assumption | Statistical interpretation is limited to the repository's internal heuristic | accepted | Implemented in `CalculateMotifScore(...)` and `CalculateEValue(...)` |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty `proteinSequence` | Returns no matches | Both public search methods guard against null or empty sequence input |
| Null, empty, or invalid `regexPattern` | `FindMotifByPattern(...)` returns no matches | The method short-circuits or swallows regex compilation failure |
| Overlapping motif occurrences | All overlapping hits are reported when the pattern allows them | The regex wrapper is lookahead-based |
| No motif occurrence in the sequence | Returns an empty collection | Only explicit regex hits are emitted |
| Lowercase or mixed-case sequence | Same result as uppercase input | Matching is case-insensitive after uppercasing |

### 6.2 Limitations

This helper is a pattern scanner, not a family classifier. It cannot infer structural context, does not consult a live motif database, and does not calibrate its scores against empirical background distributions. `FindCommonMotifs(...)` is therefore best read as a convenience scan over the motifs currently encoded in source.

## 8. References

1. PROSITE Database. SIB Swiss Institute of Bioinformatics. https://prosite.expasy.org/
2. PROSITE User Manual. https://prosite.expasy.org/prosuser.html
3. Hulo N, Bairoch A, Bulliard V, et al. The 20 years of PROSITE. Nucleic Acids Research. https://doi.org/10.1093/nar/gkm977
4. De Castro E, Sigrist CJA, Gattiker A, et al. ScanProsite. Nucleic Acids Research. https://doi.org/10.1093/nar/gkl124
5. Schneider TD, Stephens RM. Sequence logos: a new way to display consensus sequences. Nucleic Acids Research. https://doi.org/10.1093/nar/18.20.6097
