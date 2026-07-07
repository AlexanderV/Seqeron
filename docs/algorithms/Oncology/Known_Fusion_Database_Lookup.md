# Known Fusion Database Lookup

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Fusion |
| Test Unit ID | ONCO-FUSION-002 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Given a detected gene fusion (a 5'/3' partner pair), this unit (a) formats its HGNC designation as
`gene5p::gene3p`, and (b) looks that designation up against a caller-supplied set of known fusions to attach
a clinical annotation. The designation format and the directional 5'→3' keying are specification-driven,
defined exactly by the HGNC gene-fusion nomenclature [1]. The known-fusion set itself is supplied by the
caller — the library bundles no curated database — so the unit is a *Framework* algorithm: the format and
matching rule are source-backed, the data is the caller's responsibility.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A gene fusion joins two genes; one contributes the 5' (upstream) portion and the other the 3' (downstream)
portion of the chimeric transcript. Curated cancer fusion resources (COSMIC, Mitelman, ChimerDB, the WHO
classification) record recurrent, clinically meaningful fusions. To query such a resource, a detected fusion
must first be written in a single, unambiguous, standardized form [1].

### 2.2 Core Model

HGNC designation rules [1]:

- The two partner symbols are joined by a **double colon** `::`, e.g. `BCR::ABL1` [1].
- The **5' partner is always listed first**, before the `::`, "irrespective of chromosomal location or the
  orientation of the gene" [1]. The designation is therefore **directional**: `A::B` and `B::A` denote two
  different fusions (e.g. reciprocal fusions).
- Partners are written with their HGNC approved gene symbols [1].
- Read-through transcripts keep the **hyphen** (`INS-IGF2`), not `::` [1]; this unit emits `::` for true
  fusions only.

Lookup model: build the directional key `gene5p::gene3p`; a fusion is *known* iff that exact directional key
is a member of the caller's set. The reciprocal key is not a match.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Designation = `gene5p` + `::` + `gene3p` | Direct construction; HGNC 5'-first + `::` rule [1] |
| INV-02 | Directional: designation(A,B) ≠ designation(B,A) for A≠B | 5' partner fixed first by role, not sorted [1] |
| INV-03 | A match requires the directional key `5'::3'`; the reciprocal key does not match | Directionality of the designation [1] |
| INV-04 | Symbol comparison is case-insensitive, order-preserving | Implementation choice (see 5.4); does not alter the `::`/order rule |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `gene5p` | `string` | required | 5' partner symbol | non-empty, non-whitespace |
| `gene3p` | `string` | required | 3' partner symbol | non-empty, non-whitespace |
| `fusion` | `FusionCall` | required | Detected fusion; its `Gene5Prime`/`Gene3Prime` form the key | partners non-empty |
| `knownFusions` | `IReadOnlyDictionary<string,string>` | required | Map from `5'::3'` designation to annotation | non-null; case-insensitive comparer recommended |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `GetFusionAnnotation` → | `string` | The `gene5p::gene3p` designation |
| `KnownFusionMatch.Designation` | `string` | The queried fusion's designation |
| `KnownFusionMatch.IsKnown` | `bool` | Whether the directional designation was in the supplied set |
| `KnownFusionMatch.Annotation` | `string?` | Caller-supplied annotation when matched, else `null` |

### 3.3 Preconditions and Validation

Null/empty/whitespace partner symbols throw `ArgumentException`. A null `knownFusions` map throws
`ArgumentNullException`. Symbol matching is case-insensitive (ordinal-ignore-case); the designation string
itself preserves the input case verbatim. No sequence alphabet is involved.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the partner symbols (non-empty).
2. Build the directional designation `gene5p + "::" + gene3p`.
3. Probe `knownFusions` by that key (using the dictionary's comparer, then an explicit case-insensitive scan).
4. Return `KnownFusionMatch(designation, IsKnown, annotation)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The single reference constant is the separator `::` ([1]). The known-fusion set is a caller-supplied
dictionary; no constants are bundled.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GetFusionAnnotation` | O(L) | O(L) | L = symbol length (string concat) |
| `MatchKnownFusions` (case-insensitive dict) | O(L) | O(1) | hash lookup |
| `MatchKnownFusions` (fallback scan) | O(k·L) | O(1) | k = set size; only when the dict comparer is case-sensitive |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.GetFusionAnnotation(gene5p, gene3p)`: formats the `5'::3'` HGNC designation.
- `OncologyAnalyzer.MatchKnownFusions(fusion, knownFusions)`: directional lookup against a caller-supplied set.
- `OncologyAnalyzer.KnownFusionMatch`: result record (designation, IsKnown, annotation).
- `OncologyAnalyzer.FusionDesignationSeparator`: the `::` constant.

### 5.2 Current Behavior

The designation is a verbatim concatenation, so it preserves the input case of the symbols. Matching is
case-insensitive: it first uses the supplied dictionary's comparer (so an `OrdinalIgnoreCase` dictionary
resolves in O(1)); if the supplied comparer is case-sensitive it falls back to a single linear case-insensitive
scan so callers are not silently order/case-trapped.

**Search reuse:** the repository suffix tree was evaluated and not used — this is an exact dictionary key
lookup over short symbol strings, not a substring/occurrence search over a long text, so a hash-map lookup is
the appropriate structure.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Double-colon `::` separator between partner symbols [1].
- 5' partner listed first; directional designation [1].
- `BCR::ABL1` worked example reproduced exactly [1].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Curated fusion-database content (Mitelman / COSMIC / ChimerDB lists); **users should rely on:** supplying
  their own `knownFusions` map — membership and annotations are out of scope per the unit mandate.
- Read-through (`INS-IGF2`, hyphen) designations [1]; **users should rely on:** this unit emits `::`
  designations for true fusions only.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Case-insensitive symbol comparison | Assumption | Inputs vary in case; matching tolerates this | accepted | INV-04; HGNC symbols are case-defined but the source does not address case folding |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Reciprocal `B::A` only in set, query `A::B` | not known | Directional designation [1] |
| Case-varied symbols (`eml4`/`alk`) | matches `EML4::ALK` | INV-04 |
| Null/empty partner | `ArgumentException` | Validation contract |
| Null known-fusion map | `ArgumentNullException` | Validation contract |

### 6.2 Limitations

No transcript/exon-level breakpoint detail (that is ONCO-FUSION-003); no read-through handling; the known set
and its annotations are supplied by the caller, so coverage/accuracy depend entirely on that data.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var known = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase)
{
    ["BCR::ABL1"] = "Chronic myeloid leukemia driver",
    ["EML4::ALK"] = "NSCLC driver, ALK TKI target",
};

string designation = OncologyAnalyzer.GetFusionAnnotation("BCR", "ABL1"); // "BCR::ABL1"

var call = new OncologyAnalyzer.FusionCall("EML4", "ALK", 5, 4, 9,
    OncologyAnalyzer.FusionReadingFrame.InFrame);
OncologyAnalyzer.KnownFusionMatch m = OncologyAnalyzer.MatchKnownFusions(call, known);
// m.IsKnown == true, m.Designation == "EML4::ALK", m.Annotation == "NSCLC driver, ALK TKI target"
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_MatchKnownFusions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_MatchKnownFusions_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-FUSION-002-Evidence.md](../../../docs/Evidence/ONCO-FUSION-002-Evidence.md)
- Related algorithms: [Fusion_Gene_Detection](Fusion_Gene_Detection.md)

## 8. References

1. Bruford EA, Antonescu CR, Carroll AJ, et al. 2021. HUGO Gene Nomenclature Committee (HGNC) recommendations for the designation of gene fusions. *Leukemia* 35(11):3040–3043. https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/ (DOI: 10.1038/s41375-021-01436-6)
2. 2021. Recommendations for future extensions to the HGNC gene fusion nomenclature. *Leukemia* 35(11):3044–3045. https://pmc.ncbi.nlm.nih.gov/articles/PMC8632684/ (DOI: 10.1038/s41375-021-01438-4)
