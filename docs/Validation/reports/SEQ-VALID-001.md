# Validation Report: SEQ-VALID-001 — Sequence Validation

- **Validated:** 2026-06-12   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)`, `SequenceExtensions.IsValidRna(ReadOnlySpan<char>)`, `DnaSequence.TryCreate(string, out DnaSequence?)`, `DnaSequence(string)` ctor, `SequenceBase.IsValid()`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Nucleic acid notation** (https://en.wikipedia.org/wiki/Nucleic_acid_notation) — Confirms the four canonical DNA bases are G, C, A, T and that U replaces T in RNA. Confirms degenerate/ambiguity symbols (N, R, Y, …) "represent uncertainty rather than actual nucleotides" — "each particular sequence will have in fact one of the regular bases"; they encode consensus/SNP/sequencing-error positions. Confirms lowercase letters are used in sequence files (soft-masking) and denote the *same* nucleotides as uppercase.
- **Bioinformatics.org IUPAC codes** (https://www.bioinformatics.org/sms/iupac.html) — Confirms standard bases A/C/G/T (T or U) and the full degenerate set R, Y, S, W, K, M, B, D, H, V, N, plus `.`/`-` as gap. Matches NC-IUB 1984.
- **Biopython `Bio.Data.IUPACData`** (raw master) — Exact values: `unambiguous_dna_letters = "GATC"`, `unambiguous_rna_letters = "GAUC"`, `ambiguous_dna_letters = "GATCRYWSMKHBVDN"`. The canonical methods validate against the *unambiguous* sets, which is the correct reference behaviour for "is this actual sequence data".

### Edge-case semantics check (all explicitly defined and source-defensible)
| Edge case | Defined behaviour | Source |
|-----------|-------------------|--------|
| Empty sequence | **valid (true)** — vacuous truth | Biopython: "Zero-length sequences are always considered to be defined"; logically, ∀-over-empty-set is true |
| All standard DNA `ACGT` | valid | IUPAC 1970 / Biopython `unambiguous_dna_letters="GATC"` |
| All standard RNA `ACGU` | valid | IUPAC 1970 / Biopython `unambiguous_rna_letters="GAUC"` |
| Lowercase / mixed case | valid (case-insensitive) | Wikipedia: lowercase = same nucleotides (soft-mask) |
| `U` in DNA | invalid | U is RNA-only |
| `T` in RNA | invalid | T is DNA-only |
| Ambiguity codes N,R,Y,S,W,K,M,B,D,H,V | invalid | NC-IUB 1984 / Wikipedia: encode uncertainty, not a base |
| Whitespace / numeric / special / unicode / gap `-` | invalid | not part of unambiguous nucleotide notation |

These are all *defined* (not implementation-defined) and match the chosen reference (Biopython unambiguous alphabets). The one genuinely-arbitrary checklist item ("Empty → true or false? define!") is resolved in favour of **true**, which is the only choice consistent with Biopython and with the universal-quantifier semantics of "all chars are valid".

### Independent cross-check (exact values)
Biopython `unambiguous_dna_letters` = `"GATC"` → `set("ACGT")` ; `unambiguous_rna_letters` = `"GAUC"` → `set("ACGU")`. A regex/alphabet membership check `all(c.upper() in "ACGT" for c in s)` reproduces every row above (e.g. `"ACGN"` → N ∉ set → False; `""` → `all([])` → True; `"acgt"` → True).

### Findings / divergences
None. Stage A semantics are fully sourced and defensible.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:287` `IsValidDna` — loops chars, `char.ToUpperInvariant(c)`, rejects anything ≠ A/C/G/T. Empty span → loop body never runs → returns `true`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:302` `IsValidRna` — identical with A/C/G/U.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:22` ctor — null/empty → empty sequence (valid); else `ToUpperInvariant` then `ValidateSequence` (`DnaSequence.cs:112`) throwing `ArgumentException` on first non-ACGT char.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:129` `TryCreate` — wraps ctor, catches `ArgumentException`, returns false + null.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/ISequence.cs:96` `SequenceBase.IsValid()` — membership against the type `Alphabet` set (sequence pre-uppercased in ctor `:84`).

### Formula realised correctly?
Yes. The membership test exactly implements "every char ∈ unambiguous alphabet, case-insensitively". `ToUpperInvariant` gives the case-insensitivity required by INV-3/INV-4 and avoids culture-dependent casing (e.g. Turkish-I) by using the invariant culture. Ambiguity codes, U-in-DNA, T-in-RNA, whitespace, digits, unicode, gap all fall through to the reject branch.

### Cross-verification table recomputed vs code (traced)
| Input | IsValidDna | IsValidRna | Matches Stage A |
|-------|-----------|-----------|-----------------|
| `""` | true | true | ✓ |
| `"ACGT"` | true | false (T) | ✓ |
| `"ACGU"` | false (U) | true | ✓ |
| `"acgt"` / `"AcGt"` | true | — | ✓ |
| `"ACGN"` | false | false | ✓ |
| `"ACGX"` / `"X"` | false | false | ✓ |
| `"AC GT"` / tab / newline | false | false | ✓ |
| `"ACG日"` | false | false | ✓ |
All 64 SEQ-VALID-001 tests (canonical + composition properties) pass, confirming each row.

### Variant/delegate consistency
- `TryCreate` ⟺ `IsValidDna` (INV-5): `TryCreate` succeeds exactly when ctor's `ValidateSequence` (same ACGT predicate) doesn't throw. Consistent.
- Ctor normalizes to uppercase (`F5`), empty → empty sequence (`F6`), invalid → `ArgumentException` (`F4`). Consistent.
- `SequenceBase.IsValid()` uses the alphabet set; for `DnaSequence` the canonical guarantee is enforced at construction so any constructed instance is valid by invariant.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_SequenceValidation_Tests.cs` — assertions are exact `Is.True`/`Is.False` (no "no-throw" tautologies), deterministic, and cover every Stage-A edge case: empty, all-valid, lowercase, mixed case, U-in-DNA, T-in-RNA, X, digit, whitespace, N, all 10 ambiguity codes (R,Y,S,W,K,M,B,D,H,V), unicode, tab, newline, gap, boundary positions (start/mid/end), all-same, single base, long sequence, and case-invariance. FsCheck properties `PureAcgt_IsValidDna` / `PureAcgu_IsValidRna` lock INV-1/INV-2.

### Findings / defects
None.

## Verdict & follow-ups
Both stages PASS. Implementation faithfully realises the source-validated unambiguous-nucleotide semantics for every edge case. No code or test changes required. No logged defects.

State: **CLEAN**.
