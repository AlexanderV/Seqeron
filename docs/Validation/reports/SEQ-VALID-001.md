# Validation Report: SEQ-VALID-001 — Sequence Validation

- **Validated:** 2026-06-24   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)`, `SequenceExtensions.IsValidRna(ReadOnlySpan<char>)`, `DnaSequence.TryCreate(string, out DnaSequence?)`, `DnaSequence(string)` ctor (+ `ValidateSequence`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened (this session) & what they confirm
- **Biopython `Bio/Data/IUPACData.py`** (raw master, fetched 2026-06-24) — exact verbatim values:
  `unambiguous_dna_letters = "GATC"`, `unambiguous_rna_letters = "GAUC"`,
  `ambiguous_dna_letters = "GATCRYWSMKHBVDN"`, `ambiguous_rna_letters = "GAUCRYWSMKHBVDN"`.
  The canonical methods validate against the **unambiguous** sets → {A,C,G,T} for DNA, {A,C,G,U}
  for RNA. This is the correct reference for "is this actual unambiguous sequence data".
- **bioinformatics.org/sms/iupac.html** (fetched 2026-06-24) — 4 standard bases (A, C, G, T-or-U);
  11 degenerate codes (R, Y, S, W, K, M, B, D, H, V, N); gap is `.` or `-`. Matches NC-IUB 1984.
  Ambiguity codes "represent positions where the exact nucleotide is unknown or variable" — i.e.
  uncertainty, not actual nucleotides — so they are correctly rejected by unambiguous validation.

### Formula / semantics check
The validation predicate is membership: every character (case-insensitive) ∈ the unambiguous
alphabet. DNA = {A,C,G,T}; RNA = {A,C,G,U}. This matches Biopython `unambiguous_*_letters` exactly.

### Edge-case semantics (all defined & source-backed)
| Edge case | Defined behaviour | Source |
|-----------|-------------------|--------|
| Empty sequence | **valid (true)** — vacuous truth | Biopython: "Zero-length sequences are always considered to be defined"; ∀ over ∅ = true |
| `ACGT` (DNA) / `ACGU` (RNA) | valid | Biopython `unambiguous_dna_letters="GATC"` / `unambiguous_rna_letters="GAUC"` |
| Lowercase / mixed case | valid | lowercase = same nucleotides (soft-mask convention) |
| `U` in DNA / `T` in RNA | invalid | U is RNA-only, T is DNA-only |
| N, R, Y, S, W, K, M, B, D, H, V | invalid | NC-IUB 1984 / bioinformatics.org: encode uncertainty, not a base |
| whitespace / digit / special / unicode / gap `-` | invalid | not part of unambiguous nucleotide notation |

### Independent cross-check (exact values)
`set("GATC") = {A,C,G,T}`, `set("GAUC") = {A,C,G,U}`. The membership rule
`all(c.upper() in "ACGT" for c in s)` reproduces every row: `""`→`all([])`→True; `"acgt"`→True;
`"ACGU"`→U∉set→False (DNA); `"ACGN"`→N∉set→False; `"ACGT"`→T∉{A,C,G,U}→False (RNA). The 11
ambiguity codes and gap `-` all lie outside both unambiguous sets → False.

### Findings / divergences
None. Stage A semantics are fully sourced and defensible. The one genuinely-arbitrary item
(empty → true/false) is resolved to **true**, the only choice consistent with Biopython and the
universal-quantifier reading of "all chars are valid".

## Stage B — Implementation

### Code path reviewed (current line numbers)
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:415` `IsValidDna` — loops,
  `char.ToUpperInvariant(c)`, rejects anything ≠ A/C/G/T. Empty span → loop body never runs → `true`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:430` `IsValidRna` — identical with A/C/G/U.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:22` ctor — null/empty → empty
  sequence (valid); else `ToUpperInvariant` then `ValidateSequence` (`:112`) throwing
  `ArgumentException` on first non-ACGT char.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:129` `TryCreate` — wraps ctor,
  catches `ArgumentException`, returns false + null.

### Formula realised correctly?
Yes. Membership test exactly implements "every char ∈ unambiguous alphabet, case-insensitively".
`ToUpperInvariant` gives the required case-insensitivity (INV-3/INV-4) using invariant culture
(avoids Turkish-I culture hazards). Ambiguity codes, U-in-DNA, T-in-RNA, whitespace, digits,
unicode, gap all fall through to the reject branch.

### Cross-verification table recomputed vs code (traced + tested)
| Input | IsValidDna | IsValidRna |
|-------|-----------|-----------|
| `""` | true | true |
| `"ACGT"` | true | false (T) |
| `"ACGU"` | false (U) | true |
| `"acgt"` / `"AcGt"` | true | — |
| `"ACGN"` | false | false |
| `"ACGX"` / `"X"` | false | false |
| `"AC GT"` / tab / newline | false | false |
| `"AC-GT"` (gap) | false | — |
| `"ACG日"` (unicode) | false | — |
All rows confirmed by the passing tests below.

### Variant/delegate consistency
- `TryCreate` ⟺ `IsValidDna` (INV-5): `TryCreate` succeeds exactly when `ValidateSequence` (same
  ACGT predicate) doesn't throw. Consistent.
- Ctor normalizes to uppercase (F5: `"acgt"`→`Sequence="ACGT"`), empty → empty sequence (F6),
  invalid → `ArgumentException` (F4). Consistent.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_SequenceValidation_Tests.cs` — assertions
are exact `Is.True`/`Is.False` (no no-throw tautologies), deterministic, and cover every Stage-A
edge case: empty, all-valid, lowercase, mixed, U-in-DNA, T-in-RNA, X, digit, whitespace, N, all 10
ambiguity codes R/Y/S/W/K/M/B/D/H/V, unicode, tab, newline, gap, boundary positions, all-same,
single base, long sequence, case-invariance. FsCheck properties `PureAcgt_IsValidDna` /
`PureAcgu_IsValidRna` (`Properties/SequenceCompositionProperties.cs:37,49`) lock INV-1/INV-2.

### Findings / defects
None.

## Verdict & follow-ups
Both stages PASS. Source-validated unambiguous-nucleotide semantics are faithfully realised for
every edge case. Source code and test file are unchanged since the prior (cb113ce) validation;
re-confirmed against freshly-fetched Biopython `IUPACData.py` and bioinformatics.org IUPAC table.
Filtered test run: **64 passed / 0 failed** (SequenceValidation + SequenceComposition properties).
No code or test changes required. No logged defects.

State: **CLEAN**.
