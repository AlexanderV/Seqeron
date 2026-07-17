---
type: concept
title: "Sequence validation (strict DNA/RNA alphabet membership)"
tags: [sequence-statistics, validation, algorithm]
mcp_tools:
  - is_valid_dna
  - is_valid_rna
sources:
  - docs/algorithms/Sequence_Composition/Sequence_Validation.md
source_commit: 600b8a4c4de9095474c221b2227966d0982588af
created: 2026-07-17
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: sequence-validation-spec
      evidence: "Test Unit ID: SEQ-VALID-001; Algorithm Group: Sequence Composition; Algorithm: Sequence Validation (strict DNA/RNA alphabet membership)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:fasta-parsing
      source: sequence-validation-spec
      evidence: "The FASTA parser's opt-in SequenceAlphabet overloads (StrictDna/IupacNucleotide/Rna/Protein) are the parse-time consumer of alphabet validation; this unit is the standalone character-membership primitive."
      confidence: medium
      status: current
---

# Sequence validation (strict DNA/RNA alphabet membership)

**Sequence validation** answers a single yes/no question: does a string contain **only** valid
nucleotide characters for a chosen alphabet? In this repository the documented surface is a
**strict character-membership scan** — DNA accepts only `{A, C, G, T}`, RNA only `{A, C, G, U}` —
run in linear time over the sequence. Validated as test unit **SEQ-VALID-001** (CLEAN, PASS/PASS
2026-06-24); [[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes
the artifact pattern. The implementation is [[research-grade-limitations|research-grade]].

## Scientific basis — the IUPAC alphabet

The authoritative basis is the **IUPAC/IUB nomenclature** for nucleic acids (IUPAC-IUB 1970;
NC-IUB 1984). The canonical *unambiguous* alphabets are `{A, C, G, T}` (DNA: Adenine, Cytosine,
Guanine, Thymine) and `{A, C, G, U}` (RNA, with Uracil replacing Thymine). The same standard also
defines **ambiguity codes** for incompletely specified bases — `R`(A/G), `Y`(C/T·U), `S`(G/C),
`W`(A/T·U), `K`(G/T·U), `M`(A/C), `B`(not A), `D`(not C), `H`(not G), `V`(not T·U), `N`(any) — plus
the gap `-`. Any validation rule derived from the standard reduces to a per-character set-membership
test against the selected valid set.

## Strict mode — the only documented mode

This library implements **strict validation only**: over the unambiguous alphabet, so IUPAC
**ambiguity codes and the gap `-` are rejected**, even though the standard defines them. There is no
IUPAC-permissive validation mode on this surface (listed as *not implemented*). Case is handled by
**uppercase normalization before comparison**, so validation is **case-insensitive** — a common
bioinformatics convenience the IUPAC standard itself does not prescribe.

| Aspect | IUPAC standard | This implementation | Reason |
|--------|----------------|---------------------|--------|
| Ambiguity codes (`N`, `R`, `Y`, …) | Defined | **Rejected** | Strict mode |
| Gap character `-` | Defined | **Rejected** | Strict mode |
| Case | Not specified | Case-insensitive | Common practice |

## Method surface and contract

Primary spec `docs/algorithms/Sequence_Composition/Sequence_Validation.md`; implementations in
`Seqeron.Genomics.Core` (`SequenceExtensions.cs`, `DnaSequence.cs`):

| Entry point | Input | Returns | Behavior |
|-------------|-------|---------|----------|
| `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)` | span | `bool` | `true` iff every char (uppercased) ∈ `{A,C,G,T}` |
| `SequenceExtensions.IsValidRna(ReadOnlySpan<char>)` | span | `bool` | `true` iff every char (uppercased) ∈ `{A,C,G,U}` |
| `DnaSequence.TryCreate(string, out DnaSequence?)` | string | `bool` + `DnaSequence?` | factory-style validate-and-materialize: `true` + instance on success; `false` + `null` when `DnaSequence` construction raises `ArgumentException` |

`IsValidDna`/`IsValidRna` are the direct predicates; `TryCreate` is the **exception-free
construction** wrapper that delegates validation to the `DnaSequence` constructor (so it enforces
whatever `DnaSequence` accepts). The MCP tools `is_valid_dna` / `is_valid_rna` are the thin wrappers
over the two predicates (see [[mcp-tool-catalog]]).

## Algorithm and complexity

A single left-to-right pass: uppercase each character and return `false` immediately on the first
symbol outside the allowed set; return `true` if the scan completes. **O(n)** time, **O(1)** space
for `IsValidDna`/`IsValidRna`; `TryCreate` is **O(n)** time / **O(1)** auxiliary, delegating to
construction. There is no search or indexing, so the repository suffix tree does not apply.

## Invariants

- **INV-01** DNA validation accepts only `{A, C, G, T}` (unambiguous alphabet).
- **INV-02** RNA validation accepts only `{A, C, G, U}` (unambiguous alphabet).
- **INV-03** validation fails (short-circuits) at the **first** symbol outside the allowed set.

## Edge cases

| Input | `IsValidDna` | `IsValidRna` | Rationale |
|-------|--------------|--------------|-----------|
| `""` (empty) | `true` | `true` | vacuous truth — no invalid character present |
| `"ACGT"` | `true` | `false` | `T` ∉ RNA alphabet |
| `"ACGU"` | `false` | `true` | `U` ∉ DNA alphabet |
| `"acgt"` | `true` | `false` | uppercase-normalized before comparison |
| `"ACGN"` | `false` | `false` | ambiguity code `N` rejected in strict mode |
| `"AC GT"` | `false` | `false` | whitespace is outside the alphabet |

## Relationship to neighbouring validation surfaces

- **FASTA parse-time validation** — [[fasta-parsing]] exposes opt-in `SequenceAlphabet` overloads
  that validate the *uppercased* record sequence against a **selectable** alphabet
  (`StrictDna` = A/C/G/T; `IupacNucleotide` = A/C/G/T/U + ambiguity codes + gap `-`; `Rna` = A/C/G/U;
  `Protein` = 20 residues + `B Z J X U O` + stop `*`), throwing `ArgumentException` on the first
  out-of-alphabet character. That is the multi-alphabet, IUPAC-permissive validator at the parser
  boundary; **this** unit is the standalone strict DNA/RNA-only character-membership primitive.
- **IUPAC degenerate matching** — [[iupac-degenerate-matching]] *interprets* ambiguity codes
  (a `V` matches A/C/G) rather than rejecting them; validation here is a membership test, not a
  degeneracy expansion.
- **Composition** — [[base-composition]] counts every character (routing non-canonical symbols into
  `CountN`/`CountOther`) instead of failing; validation is the go/no-go gate, composition the tally.

## References

IUPAC-IUB Commission on Biochemical Nomenclature (1970) *Biochemistry* 9(20):4022-4027
(doi:10.1021/bi00822a023); NC-IUB (1984) "Nomenclature for Incompletely Specified Bases in Nucleic
Acid Sequences" *Nucleic Acids Research* 13(9):3021-3030 (doi:10.1093/nar/13.9.3021).
