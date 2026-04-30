# PAM Site Detection Algorithm

| Field | Value |
|-------|-------|
| Algorithm Group | Molecular Tools |
| Test Unit ID | CRISPR-PAM-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

PAM (protospacer adjacent motif) site detection identifies DNA sequences adjacent to potential CRISPR targets that are required for Cas nuclease binding and cleavage. In this repository, the implementation searches both strands for system-specific PAM patterns, extracts the corresponding target region when it fits within sequence bounds, and returns site metadata describing the matched PAM, target, strand, and CRISPR system.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A PAM is a short DNA sequence adjacent to a CRISPR target site, and Cas nucleases depend on the presence of that PAM to bind and cleave. The original document states that PAM recognition helps distinguish self from non-self DNA in the CRISPR immune system and that Cas9 does not successfully bind or cleave a target that lacks the required PAM. Sources: Wikipedia (PAM, CRISPR), Jinek et al. (2012).

### 2.2 Core Model

The algorithm scans the input sequence for system-specific PAM patterns defined with IUPAC ambiguity codes. For PAM-after-target systems such as Cas9, the guide target lies immediately upstream of the PAM:

$$
targetStart = PAM_{pos} - guideLength, \quad targetEnd = PAM_{pos} - 1
$$

For PAM-before-target systems such as Cas12a, the guide target lies immediately downstream of the PAM:

$$
targetStart = PAM_{pos} + PAM_{length}, \quad targetEnd = targetStart + guideLength - 1
$$

The same search is repeated on the reverse complement to identify reverse-strand sites.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | A returned PAM site always satisfies the system's PAM pattern under IUPAC matching | `FindPamSitesCore(...)` yields only after `MatchesPam(...)` succeeds |
| INV-02 | A returned target sequence always fits within the scanned sequence bounds | The source checks `targetStart >= 0` and `targetEnd < seq.Length` before yielding |
| INV-03 | Both forward and reverse strands are searched | The implementation scans `seq` and `revComp` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to scan for PAMs | String overload normalizes to uppercase; empty string returns no results |
| `systemType` | `CrisprSystemType` | `SpCas9` | CRISPR system whose PAM definition is used | Must be one of the supported systems in `CrisprDesigner.GetSystem(...)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Start position of the PAM |
| `PamSequence` | `string` | Actual matched PAM sequence |
| `TargetSequence` | `string` | Guide-length target extracted from the scanned strand |
| `TargetStart` | `int` | Start coordinate recorded for the target sequence |
| `IsForwardStrand` | `bool` | `true` for forward-strand matches, `false` for reverse-strand matches |
| `System` | `CrisprSystem` | System metadata including name, PAM pattern, guide length, PAM orientation, and description |

### 3.3 Preconditions and Validation

The `DnaSequence` overload throws `ArgumentNullException` for null input. The raw-string overload returns an empty result for null or empty input and uppercases the sequence before scanning. If the sequence is shorter than the effective PAM-plus-target span or if a matched PAM would place the target outside sequence bounds, no `PamSite` is yielded.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the selected CRISPR system into its PAM pattern, guide length, and PAM orientation.
2. Scan the forward strand for PAM matches using IUPAC-aware pattern matching.
3. For each forward match, compute the target interval and yield a site only when the target fits within bounds.
4. Reverse-complement the input sequence and repeat the PAM scan.
5. Convert reverse-strand PAM positions back to forward-strand coordinates and yield the reverse-strand site metadata.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Supported CRISPR systems in the current implementation:

| System | PAM Sequence | PAM Location | Guide Length | Notes |
|--------|--------------|--------------|--------------|-------|
| SpCas9 | NGG | 3' of target | 20 bp | Canonical Cas9 from *Streptococcus pyogenes* |
| SpCas9-NAG | NAG | 3' of target | 20 bp | Lower-efficiency variant |
| SaCas9 | NNGRRT | 3' of target | 21 bp | *Staphylococcus aureus* Cas9 |
| Cas12a (Cpf1) | TTTV | 5' of target | 23 bp | `V = A, C, or G` |
| AsCas12a | TTTV | 5' of target | 23 bp | *Acidaminococcus sp.* |
| LbCas12a | TTTV | 5' of target | 24 bp | *Lachnospiraceae bacterium* |
| CasX | TTCN | 5' of target | 20 bp | Compact Cas protein |

IUPAC codes used by the PAM matcher in the original document:

| Code | Meaning | Nucleotides |
|------|---------|-------------|
| N | Any | A, C, G, T |
| R | Purine | A, G |
| V | Not T | A, C, G |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| PAM scan | `O(n)` | `O(k)` | `n` is sequence length and `k` is the number of yielded PAM sites |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CrisprDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs), [IupacHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/IupacHelper.cs)

- `CrisprDesigner.FindPamSites(DnaSequence, CrisprSystemType)`: Scans a validated DNA sequence for PAM sites.
- `CrisprDesigner.FindPamSites(string, CrisprSystemType)`: Raw-string overload that uppercases input and yields no results for empty input.
- `CrisprDesigner.GetSystem(CrisprSystemType)`: Maps the enum to PAM sequence, guide length, orientation, and description.
- `IupacHelper.MatchesIupac(char, char)`: Evaluates per-base IUPAC ambiguity matches.

### 5.2 Current Behavior

The current implementation searches both forward and reverse strands, uses `IupacHelper.MatchesIupac(...)` for ambiguity-code evaluation, and returns `PamSite` records carrying the matched PAM, target sequence, strand, and system metadata. For reverse-strand matches, `Position` is converted back to a forward-strand coordinate and `PamSequence` is reverse-complemented back to the forward-oriented PAM string. `TargetSequence` is extracted from the reverse-complement scan path, and `TargetStart` is recorded from that same reverse-strand traversal.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PAM recognition as a prerequisite for guide-target extraction.
- IUPAC-based ambiguity handling for PAM patterns.
- System-specific guide lengths and PAM-before/PAM-after orientation rules for the supported systems.

**Intentionally simplified:**

- The implementation is limited to the fixed set of CRISPR systems declared in `GetSystem(...)`; **consequence:** PAMs for other Cas proteins are not discoverable through this API.
- Site discovery is sequence-only; **consequence:** chromatin state, cleavage efficiency, and organism-specific activity effects are not reflected in the returned sites.

**Not implemented:**

- Guide efficacy or cleavage-efficiency prediction for the detected PAM sites; **users should rely on:** `Guide_RNA_Design.md` and downstream evaluation workflows rather than PAM detection alone.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns an empty collection | No scanable sequence content |
| Null `DnaSequence` input | Throws `ArgumentNullException` | Explicit guard in source |
| Sequence shorter than PAM plus guide | Returns no valid sites | Bounds checks reject incomplete target intervals |
| No PAM matches | Returns an empty collection | No valid PAM pattern exists in the sequence |
| PAM at sequence boundary | Included only if the target also fits within bounds | Target extraction is validated before yielding |
| Multiple overlapping PAMs | All distinct sites are returned | The scan yields every matching window |
| Lowercase input to the string overload | Case-insensitive matching | The string overload uppercases the input |

### 6.2 Limitations

The current implementation is a sequence-pattern detector rather than a full CRISPR activity model. It does not score cleavage efficiency, off-target risk, chromatin accessibility, or non-listed PAM systems, and reverse-strand target metadata follows the current scan-path representation used in source.

## 8. References

1. Wikipedia: Protospacer adjacent motif. https://en.wikipedia.org/wiki/Protospacer_adjacent_motif
2. Wikipedia: CRISPR. https://en.wikipedia.org/wiki/CRISPR
3. Jinek M, et al. (2012). "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity". Science 337(6096):816-821.
4. Zetsche B, et al. (2015). "Cpf1 is a single RNA-guided endonuclease of a class 2 CRISPR-Cas system". Cell 163(3):759-771.
5. Anders C, et al. (2014). "Structural basis of PAM-dependent target DNA recognition by the Cas9 endonuclease". Nature 513(7519):569-573.
