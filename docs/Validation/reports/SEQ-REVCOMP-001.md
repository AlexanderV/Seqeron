# Validation Report: SEQ-REVCOMP-001 — Reverse Complement

- **Validated:** 2026-06-12   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.TryGetReverseComplement(ReadOnlySpan<char>, Span<char>)`; helper `SequenceExtensions.GetComplementBase(char)`; variants `DnaSequence.ReverseComplement()`, `DnaSequence.GetReverseComplementString(string)`, `DnaSequence.TryWriteReverseComplement(Span<char>)`; complement-only `TryGetComplement`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Nucleic acid notation"** (https://en.wikipedia.org/wiki/Nucleic_acid_notation): full IUPAC ambiguity table with complements. Confirmed complement of every code matches the spec/implementation exactly: A↔T, U→A, G↔C, R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N.
- **Wikipedia — "Complementarity (molecular biology)"** (https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)): defines complementary sequence (A–T/U, G–C) and antiparallel pairing — confirms `ReverseComplement = Reverse(Complement(s))`.
- **Wikipedia — "Nucleic acid sequence"** (https://en.wikipedia.org/wiki/Nucleic_acid_sequence): quotes exactly *"the complementary sequence to TTAC is GTAA"* and *"…the base on each position in the complementary (i.e., A to T, C to G) and in the reverse order."*
- **Biopython 1.79 `Bio.Seq`** (https://biopython.org/docs/1.79/api/Bio.Seq.html): docstrings confirm `Seq("CCCCCGATAGNR").reverse_complement() → "YNCTATCGGGGG"` (with the note "R = G or A, its complement is Y"); `complement("ACTG-NH") → "TGAC-ND"`; gaps `-` are preserved in place; ambiguity codes complemented per IUPAC.

### Formula check
`ReverseComplement(s) = Reverse(Complement(s))`, complement base-by-base per Watson-Crick + IUPAC NC-IUB 1984. Matches cited Wikipedia definition and Biopython reference behaviour verbatim.

### Edge-case semantics check
- Empty → empty (`reverse_complement("") → ""`, returns true): sourced (Biopython).
- Single base → its complement: sourced (Watson-Crick).
- Palindrome (EcoRI `GAATTC`, BamHI `GGATCC`, HindIII `AAGCTT`) = own reverse complement: hand-verified below.
- Destination too small → false, no partial writes: API safety contract (not a biological claim).
- RNA `U` → `A` under DNA-centric `GetComplementBase` (output uppercase, DNA bases): consistent with IUPAC (U pairs with A); deviation D2 documented in spec — `GetRnaComplementBase` exists for RNA-output cases.
- Ambiguity codes + gaps: sourced (IUPAC table + Biopython gap pass-through).

### Independent cross-check (exact values, hand-computed against the IUPAC table)
- `CCCCCGATAGNR` → reverse `RNGATAGCCCCC` → complement `YNCTATCGGGGG` ✓ (Biopython)
- `ACTG-NH` revcomp → reverse `HN-GTCA` → `DN-CAGT` ✓ (Biopython)
- `ACTG-NH` complement → `TGAC-ND` ✓ (Biopython)
- `TTAC` revcomp → reverse `CATT` → `GTAA` ✓ (Wikipedia)
- `GAATTC` revcomp → `GAATTC` ✓ (palindrome)

### Findings / divergences
None affecting correctness. D1 (always-uppercase output) and D2 (DNA-centric complement of U) are explicitly documented and consistent with IUPAC convention.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:138` `GetComplementBase` — switch table; all 15 IUPAC codes (case-insensitive) + U; non-IUPAC (gaps) pass through via `_ => nucleotide`.
- `SequenceExtensions.cs:205` `TryGetReverseComplement` — length guard `destination.Length < sequence.Length → false` (checked before any write, so no partial writes); loop writes `GetComplementBase(sequence[len-1-i])`.
- `SequenceExtensions.cs:189` `TryGetComplement` — same guard, forward order.
- `DnaSequence.cs:68` `ReverseComplement()`, `:149` `GetReverseComplementString` (null/empty → input), `:191` `TryWriteReverseComplement` delegating to the span method.

### Formula realised correctly?
Yes. The index `sequence.Length - 1 - i` realises the reverse; `GetComplementBase` realises the complement. Guard placed before the write loop guarantees the documented "no partial writes" contract.

### Cross-verification table recomputed vs code
All five Stage-A values reproduced by the corresponding passing tests (`..._CCCCCGATAGNR`, `..._ACTG_NH`, `TryGetComplement_..._ACTG_NH`, `..._WikipediaExample_ReturnsGTAA`, palindrome TestCases). Hand-trace of `GetComplementBase` matches each.

### Variant/delegate consistency
`ReverseComplement()`, `TryWriteReverseComplement`, `GetReverseComplementString`, and `DnaSequence.TryGetReverseComplement` all route through `GetComplementBase` with identical reverse indexing — consistent with the canonical. (Note: `DnaSequence` ctor validates A/C/G/T only, so the instance variant never receives ambiguity codes; the span helper does and handles them correctly.)

### Test quality audit
`tests/.../SequenceExtensions_ReverseComplement_Tests.cs` (61 test runs) asserts exact sourced values, not tautologies. Edge cases covered: empty (buffer untouched), empty src+dst, all single bases incl. U, destination-too-small (asserts no partial write), empty destination, destination-larger, case-insensitivity, RNA U, all 11 ambiguity codes + involution + lowercase, gap pass-through, 3 Biopython cross-checks, 100-base sequences, static helper incl. null. Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
Both stages PASS. **STATE: CLEAN** — no defect found. Build succeeds; ReverseComplement filter = 61 passed; full suite = 4461 passed, 0 failed. No code or test changes required.
