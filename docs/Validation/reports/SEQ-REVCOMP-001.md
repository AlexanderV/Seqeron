# Validation Report: SEQ-REVCOMP-001 — Reverse Complement

- **Validated:** 2026-06-24   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.TryGetReverseComplement(ReadOnlySpan<char>, Span<char>)`; helper `SequenceExtensions.GetComplementBase(char)`; variants `DnaSequence.ReverseComplement()`, `DnaSequence.GetReverseComplementString(string)`, `DnaSequence.TryWriteReverseComplement(Span<char>)`, `DnaSequence.TryGetReverseComplement(...)`; complement-only `TryGetComplement`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Rosalind — Reverse Complement (REVC)** (https://rosalind.info/problems/revc/): defines the reverse complement as "the DNA string formed by reversing and complementing each symbol" (A↔T, C↔G). Sample dataset `AAAACCCGGT` → sample output `ACCGGGTTTT`. This confirms `ReverseComplement = Reverse(Complement(s))` with Watson-Crick complement.
- **Wikipedia — "Nucleic acid notation"** (https://en.wikipedia.org/wiki/Nucleic_acid_notation): full IUPAC NC-IUB 1984 ambiguity table with complements, retrieved this session. Confirms every code: A↔T, U→A, G↔C, R↔Y, Y↔R, S↔S, W↔W, M↔K, K↔M, B↔V, D↔H, H↔D, V↔B, N↔N — identical to the spec §3 table and to `GetComplementBase`.
- Prior session also cross-confirmed against Wikipedia "Complementarity (molecular biology)", Wikipedia "Nucleic acid sequence" (`TTAC` → `GTAA`), and Biopython `Bio.Seq` (`reverse_complement`/`complement` examples). Those values re-checked by hand below.

### Formula check
`ReverseComplement(s) = Reverse(Complement(s))`, complement base-by-base per Watson-Crick + IUPAC NC-IUB 1984. Matches Rosalind's definition and the Wikipedia IUPAC table verbatim.

### Edge-case semantics check
- Empty → empty (returns true): sourced (Biopython `reverse_complement("") → ""`).
- Single base → its complement: sourced (Watson-Crick / IUPAC).
- Biological palindromes (EcoRI `GAATTC`, BamHI `GGATCC`, HindIII `AAGCTT`) = own reverse complement: hand-verified.
- Destination too small → false, no partial writes: API safety contract (not a biological claim).
- RNA `U` → `A` under DNA-centric `GetComplementBase` (uppercase DNA output): consistent with IUPAC (U pairs with A); deviation D2 documented; `GetRnaComplementBase` exists for RNA-output cases.
- Ambiguity codes per IUPAC; gaps / non-IUPAC pass through unchanged (`_ => nucleotide`): sourced (IUPAC table + Biopython gap pass-through).

### Independent cross-check (exact values, hand-computed)
- Rosalind: `AAAACCCGGT` → reverse `TGGCCCAAAA` → complement `ACCGGGTTTT` ✓ (Rosalind sample output)
- `CCCCCGATAGNR` → reverse `RNGATAGCCCCC` → complement `YNCTATCGGGGG` ✓ (Biopython)
- `ACTG-NH` revcomp → reverse `HN-GTCA` → `DN-CAGT` ✓ (Biopython)
- `ACTG-NH` complement → `TGAC-ND` ✓ (Biopython)
- `TTAC` revcomp → reverse `CATT` → `GTAA` ✓ (Wikipedia)
- `GAATTC` revcomp → `GAATTC` ✓ (palindrome)

### Findings / divergences
None affecting correctness. D1 (always-uppercase output) and D2 (DNA-centric complement of U) are explicitly documented and consistent with IUPAC convention. Involution (MUST-02) is a genuine mathematical property: complement is its own inverse on the IUPAC table and reverse is self-inverse, so `RevComp∘RevComp = id`.

## Stage B — Implementation

### Code path reviewed (file:line)
- `SequenceExtensions.cs:251` `GetComplementBase` — switch table over all 16 IUPAC symbols (case-insensitive) incl. U; non-IUPAC (gaps) pass through via `_ => nucleotide`.
- `SequenceExtensions.cs:333` `TryGetReverseComplement` — guard `destination.Length < sequence.Length → false` checked before any write (no partial writes); loop writes `GetComplementBase(sequence[Length-1-i])`.
- `SequenceExtensions.cs:317` `TryGetComplement` — same guard, forward order.
- `DnaSequence.cs:68` `ReverseComplement()`, `:149` `GetReverseComplementString` (null/empty → input echoed; writes to `result[Length-1-i]` from forward read — equivalent reverse), `:191` `TryWriteReverseComplement`, `:200` static `TryGetReverseComplement` — all delegate to `GetComplementBase`.

### Formula realised correctly?
Yes. Index `sequence.Length - 1 - i` realises the reverse; `GetComplementBase` realises the complement. The length guard precedes the write loop, guaranteeing the documented "no partial writes" contract. `GetReverseComplementString` uses the mirror form (forward read, reverse-indexed write) which is equivalent.

### Cross-verification table recomputed vs code
All Stage-A values reproduced by passing tests: Rosalind logic matches `AAAACCCGGT` trace; `..._CCCCCGATAGNR`, `..._ACTG_NH` (revcomp), `TryGetComplement_..._ACTG_NH`, `..._WikipediaExample_ReturnsGTAA`, palindrome `[TestCase]`s. Hand-trace of `GetComplementBase` matches each character.

### Variant/delegate consistency
`ReverseComplement()`, `TryWriteReverseComplement`, `GetReverseComplementString`, static `TryGetReverseComplement`, and the span helper all route through `GetComplementBase` with reverse indexing — consistent. (Note: `DnaSequence` ctor validates A/C/G/T only, so its instance variant never sees ambiguity codes; the span helper does and handles them correctly.)

### Numerical robustness
N/A (no floating point). Span bounds guarded; no overflow.

### Test quality audit
`SequenceExtensions_ReverseComplement_Tests.cs` (32 runs) asserts exact sourced values, not tautologies. Edge cases covered: empty (buffer untouched), empty src+dst, all single bases incl. U, destination-too-small (asserts no partial write), empty destination, destination-larger (asserts no overwrite beyond source), case-insensitivity, RNA U×2, all 11 ambiguity codes + involution + lowercase, gap pass-through, 3 Biopython cross-checks, 100-base asymmetric sequences + involution, static helper incl. null. Deterministic.

### Findings / defects
None. Source file was touched since the prior validation (commits 6a1e8922, 6e900e92) but those changes added the RNA complement helper and Biopython GC-fraction modes; `GetComplementBase`, `TryGetReverseComplement`, and the `DnaSequence` revcomp variants are unchanged and remain correct.

## Verdict & follow-ups
Both stages PASS. **STATE: CLEAN** — no defect found. Build succeeds; ReverseComplement filter = 32 passed, 0 failed. No code or test changes required.
