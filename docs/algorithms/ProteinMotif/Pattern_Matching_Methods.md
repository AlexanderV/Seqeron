# Protein Pattern Matching Methods

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-PATTERN-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

This unit groups the four pattern-matching primitives of `ProteinMotifFinder`: enumerate occurrences of a regular-expression motif in a protein sequence (`FindMotifByPattern`), translate a PROSITE PA-line pattern into an equivalent .NET regular expression (`ConvertPrositeToRegex`), run a PROSITE pattern end-to-end against a sequence (`FindMotifByProsite`), and detect a small built-in set of domain signatures (`FindDomains`). The matching is specification-driven (exact, deterministic): a position matches iff it satisfies the regular expression derived from the PROSITE grammar. Each reported match carries an information-content score in bits and a model E-value. Use it to locate short functional motifs (glycosylation sites, phosphorylation sites, P-loops, zinc fingers) by pattern rather than by alignment.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

PROSITE describes short, biologically significant protein motifs as *patterns* — regular expressions over the 20-residue IUPAC amino-acid alphabet [1][2]. A pattern is a sequence of position elements separated by `-`; each element constrains the residue allowed at that position. Detecting a motif means finding substrings that satisfy every position constraint in order.

### 2.2 Core Model

**PROSITE PA-line grammar → regular expression** [1]:

| PROSITE element | Meaning | Regex |
|-----------------|---------|-------|
| `A` (letter) | exactly that residue (IUPAC one-letter code) | `A` |
| `x` | any residue | `.` |
| `[ALT]` | one of the listed residues ("Ala or Leu or Thr") | `[ALT]` |
| `{AM}` | any residue except those listed | `[^AM]` |
| `-` | element separator | (dropped) |
| `x(3)` | exactly 3 of the preceding | `.{3}` |
| `x(2,4)` | 2 to 4 of the preceding | `.{2,4}` |
| `A(3)` | exactly 3 of residue A (fixed counts on a letter are valid) | `A{3}` |
| `<` | N-terminus anchor | `^` |
| `>` | C-terminus anchor | `$` |
| trailing `.` | terminates the pattern | (ends parsing) |

Ranges `(n,m)` are valid only on `x`; a range on a residue letter (`A(2,4)`) is not a valid PROSITE element [1]. The Kleene star `*` (`<{C}*>`) belongs to the ScanProsite *query* extension, not the PA-line grammar [1].

**Information-content score** [3]. For an aligned position, Schneider & Stephens define the sequence conservation Rseq = log2 N − Σ pₙ log2 pₙ (bits), where N is the number of distinct symbols. For a pattern position that admits k of the 20 protein residues with uniform probability, this reduces to log2(20) − log2(k) = log2(20/k) bits. A whole-pattern score is the sum over positions:

> Score = Σᵢ log2(20 / kᵢ)

A fixed residue (k=1) contributes log2(20) ≈ 4.3219 bits; a wildcard `x` (k=20) contributes 0 bits; max per protein site is log2(20) ≈ 4.32 bits [3].

**E-value (model).** Expected number of random matches under a uniform i.i.d. amino-acid background: E = (N − L + 1) · 2^(−Score), where N is the sequence length and L the match length. This is the standard combinatorial expectation derived from the IC definition [3]; it is a model quantity, not ScanProsite's Swiss-Prot-frequency E-value.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Null/empty sequence or pattern yields no matches and no exception. | Guard clauses return early. |
| INV-02 | `Sequence == sequence.Substring(Start, End−Start+1)` (uppercased); `End = Start + len − 1`. | Match is sliced from the captured group. |
| INV-03 | Score = Σ log2(20/kᵢ) over pattern positions. | IC per position, Schneider & Stephens (1990) [3]. |
| INV-04 | `ConvertPrositeToRegex` maps each PA-line atom to its regex deterministically. | Direct grammar translation [1]. |
| INV-05 | Matching is case-insensitive; positions are 0-based. | Inputs upper-cased; `RegexOptions.IgnoreCase`. |
| INV-06 | Unsupported metacharacters (`*`,`?`,`+`) raise `FormatException`. | Reject-don't-drop policy [1]. |
| INV-07 | E-value ≥ 0 and equals (N−L+1)·2^(−Score). | Definition above [3]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| proteinSequence | string | required | Protein sequence | upper/lower-cased accepted; null/empty → empty result |
| regexPattern | string | required | .NET regex for the motif | null/empty → empty; invalid regex → empty |
| motifName | string | "Custom" | Label stored on each match | — |
| patternId | string | "" | Pattern identifier stored on each match | — |
| prositePattern | string | required (PROSITE entry points) | PROSITE PA-line | PA-line grammar only; `*`/`?`/`+` rejected |

### 3.2 Output / Return Value

`MotifMatch` record per occurrence: `Start` (0-based, inclusive), `End` (0-based, inclusive), `Sequence` (matched substring, uppercased), `MotifName`, `Pattern`, `Score` (bits), `EValue`. `FindDomains` returns `ProteinDomain` records (Name, Accession, Start, End, Score, Description).

### 3.3 Preconditions and Validation

Null/empty sequence or pattern → empty enumeration (no exception). Input is upper-cased; matching is case-insensitive. Coordinates are 0-based inclusive. An invalid .NET regex yields an empty enumeration; an unsupported PROSITE construct (`*`, `?`, `+`) throws `FormatException` from `ConvertPrositeToRegex` (and thus from `FindMotifByProsite`).

## 4. Algorithm

### 4.1 High-Level Steps

1. (PROSITE entry points) `ConvertPrositeToRegex` scans the PA line left-to-right, emitting the regex token for each atom; a trailing `.` ends parsing.
2. `FindMotifByPattern` wraps the regex in a zero-width lookahead `(?=(...))` so every start position (including overlapping ones) is enumerated.
3. For each capture, compute the IC score from the per-position allowed-residue counts, and the model E-value.
4. `FindDomains` runs a fixed set of signature regexes (e.g. P-loop `[AG].{4}GK[ST]`) through step 2–3 and wraps results as `ProteinDomain`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Allowed-residue count kᵢ per position drives the score: fixed letter → 1; `[ABC]` → class size; `[^ABC]` → 20 − size; `.` → 20; quantifier `{n}`/`{n,m}` repeats the count n times (minimum). Alphabet size 20 per the 20 standard amino acids [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindMotifByPattern | O(n·m) | O(n) | n = sequence length, m = pattern cost; lookahead scan over all starts |
| ConvertPrositeToRegex | O(p) | O(p) | p = pattern length, single left-to-right pass |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `ProteinMotifFinder.FindMotifByPattern(...)`: enumerate regex motif occurrences with IC score + E-value.
- `ProteinMotifFinder.ConvertPrositeToRegex(...)`: PROSITE PA-line → .NET regex.
- `ProteinMotifFinder.FindMotifByProsite(...)`: convert then match (delegate).
- `ProteinMotifFinder.FindDomains(...)`: run built-in signature patterns (delegate).

### 5.2 Current Behavior

Matching uses .NET `Regex` with `IgnoreCase` and a lookahead wrapper to recover overlapping occurrences. **Suffix tree was evaluated and not used:** the repository `SuffixTree` performs exact-substring search only; it cannot evaluate the character classes, negated classes, quantifiers, and anchors that PROSITE patterns require, so it does not fit this regex-based matcher. Unsupported PROSITE constructs (`*`, `?`, `+`) are rejected with `FormatException` rather than silently dropped.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Full PA-line atom translation (`x`, `[..]`, `{..}`, `x(n)`, `x(n,m)`, `A(n)`, `<`, `>`, trailing `.`) per the PROSITE syntax [1].
- Information-content score Σ log2(20/kᵢ) per Schneider & Stephens (1990) [3].

**Intentionally simplified:**

- E-value uses a uniform i.i.d. amino-acid background; **consequence:** the number differs from ScanProsite's Swiss-Prot-frequency E-value (a model quantity, see [3] for the IC basis).

**Not implemented:**

- PROSITE generalized *profiles* (weight matrices) [2]; **users should rely on:** dedicated profile tools (no in-repo alternative).
- ScanProsite extended query syntax (`*` Kleene star); **users should rely on:** rewriting the query as a PA-line pattern.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Overlapping enumeration via lookahead | Assumption | All overlapping starts are listed; PROSITE papers do not specify this | accepted | Repository contract; positions where the pattern can start are unchanged |
| 2 | Uniform-background E-value | Assumption | E-value is model-defined, not ScanProsite's | accepted | IC basis is [3] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null/empty sequence or pattern | empty enumeration | INV-01 |
| Invalid .NET regex | empty enumeration | guarded `try/catch` |
| Unsupported PROSITE `*`/`?`/`+` | `FormatException` | INV-06, reject-don't-drop [1] |
| Overlapping matches | all start positions listed | lookahead [5.2] |
| Mixed/lower case input | same matches as upper case | INV-05 |

### 6.2 Limitations

Patterns only — generalized profiles, position-specific scoring, and Swiss-Prot-frequency E-values are out of scope. `FindDomains` checks a small fixed signature set, not the full Pfam/PROSITE catalogue.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// PS00016 cell-attachment RGD: literal pattern, one match at 0-based position 3.
var m = ProteinMotifFinder.FindMotifByPattern("AAARGDAAA", "RGD", "RGD", "PS00016").Single();
// m.Start == 3, m.End == 5, m.Sequence == "RGD"
// m.Score == 3 * log2(20) == 12.965784284662089 bits  (3 fixed residues)
```

**Numerical walk-through:** RGD has three fixed positions, each allowing k=1 of 20 residues → IC = 3 × log2(20/1) = 3 × 4.321928094887363 = 12.965784284662089 bits. For `[ST].[RK]`, k = 2, 20, 2 → IC = log2(10) + 0 + log2(10) = 6.643856189774724 bits.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_FindMotifByPattern_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/ProteinMotifFinder_FindMotifByPattern_Tests.cs) — covers `INV-01`..`INV-07`
- Evidence: [PROTMOTIF-PATTERN-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md)
- Related algorithms: [PROSITE_Pattern_Matching](./PROSITE_Pattern_Matching.md), [Domain_Prediction](./Domain_Prediction.md), [Motif_Search](./Motif_Search.md)

## 8. References

1. ExPASy / PROSITE. PROSITE pattern syntax (ScanProsite documentation). https://prosite.expasy.org/scanprosite/scanprosite_doc.html (accessed 2026-06-14)
2. de Castro E, Sigrist CJA, Gattiker A, Bulliard V, Langendijk-Genevaux PS, Gasteiger E, Bairoch A, Hulo N. 2006. ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins. Nucleic Acids Research 34(Web Server):W362–W365. https://doi.org/10.1093/nar/gkl124
3. Schneider TD, Stephens RM. 1990. Sequence logos: a new way to display consensus sequences. Nucleic Acids Research 18(20):6097–6100. https://doi.org/10.1093/nar/18.20.6097
4. PROSITE entries PS00001, PS00005, PS00016, PS00017, PS00029. https://prosite.expasy.org/PS00001 (accessed 2026-06-14)
