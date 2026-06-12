# Validation Report: REP-PALIN-001 — DNA Palindrome Detection (reverse-complement palindromes / Rosalind REVP)

- **Validated:** 2026-06-12   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindPalindromes(DnaSequence, minLength, maxLength)`; overload `FindPalindromes(string, ...)`; alternate `GenomicAnalyzer.FindPalindromes(...)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Palindromic sequence** (https://en.wikipedia.org/wiki/Palindromic_sequence): confirms "a (single-stranded) nucleotide sequence is said to be a palindrome if it is equal to its reverse complement." A DNA palindrome is NOT a literal character palindrome. EcoRI example: 5'-GAATTC-3' / 3'-CTTAAG-5'; flipping the duplex gives the identical sequence.
- **Rosalind — REVP (Locating Restriction Sites)** (https://rosalind.info/problems/revp/): a reverse palindrome equals its reverse complement; report **every** such substring with **length between 4 and 12** inclusive; **1-based** position ("number of symbols to its left, including itself").
- **Wikipedia — Restriction enzyme**: Type II enzymes recognize palindromic sites, typically 4–8 bp (EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT, NotI GCGGCCGC).

### Reverse-complement definition (load-bearing)
A DNA palindrome P satisfies `P == ReverseComplement(P)`. Because each position i must complement position (len−1−i), and **no base is its own complement** (A↔T, G↔C), a contiguous reverse-complement palindrome must have **even length** (an odd-length one would require a central self-complementing base). This matches the spec's "Why Even Length Only" derivation. Note: Wikipedia mentions interrupted/spacer palindromes recognized by some enzymes can appear odd-length overall, but those are not contiguous reverse-complement palindromes and are out of scope for REVP; the even-length rule for the contiguous case is correct.

### Independent cross-check (hand-computed Rosalind REVP sample)
Input `TCAATGCATGCGGGTCTATATGCAT` (len 25). Independent Python computation (rev-comp equality, even lengths 4–12) reproduces Rosalind's published output **exactly**:

| 1-based pos | length | substring | 0-based pos (impl) |
|---|---|---|---|
| 4 | 6 | ATGCAT | 3 |
| 5 | 4 | TGCA | 4 |
| 6 | 6 | GCATGC | 5 |
| 7 | 4 | CATG | 6 |
| 17 | 4 | ATAT | 16 |
| 18 | 4 | TATA | 17 |
| 20 | 6 | ATGCAT | 19 |
| 21 | 4 | TGCA | 20 |

The spec's expected 0-based set `{(3,6),(4,4),(5,6),(6,4),(16,4),(17,4),(19,6),(20,4)}` is exactly the Rosalind output minus 1.

Standalone checks: `revcomp(GAATTC)=GAATTC` ✓ (palindrome), `revcomp(AT)=AT` ✓ (palindrome), `revcomp(GAATTA)=TAATTC ≠ GAATTA` ✓ (NOT a palindrome).

### Findings / divergences (notes)
- **Coordinate base:** Rosalind is 1-based; this implementation reports **0-based** positions. Documented in the spec and asserted by tests; consistent with the rest of the library. Note, not a defect.
- **minLength ≥ 4:** the implementation rejects `minLength < 4`. A 2-bp reverse-complement palindrome (AT, GC, CG, TA) is biologically real (AT verified above), but Rosalind REVP and the spec deliberately scope the minimum to 4 (smallest biologically relevant restriction site). This is a documented, sourced restriction — recorded as a PASS-WITH-NOTES, not a defect.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:483-538` — overloads + `FindPalindromesCore`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:127-144` — alternate.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:149-160` — `GetReverseComplementString`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:138-157` — `GetComplementBase` (correct A↔T, G↔C, plus full IUPAC set).

### Formula realised correctly
Core loop iterates `len` from minLength to maxLength in steps of **2** (even-length only), and for each window compares `candidate == GetReverseComplementString(candidate)` — i.e. true **reverse-complement equality**, NOT literal string reversal. Positions are the 0-based window start. Validation: `minLength` must be even and ≥ 4; `maxLength ≥ minLength`; null `DnaSequence` throws `ArgumentNullException`; empty string yields empty. The string overload normalizes via `ToUpperInvariant()` (case-insensitive). Both overloads share identical validation.

### Cross-verification recomputed vs code
The actual code path is exercised by test `FindPalindromes_RosalindSample_CorrectOutput`, which asserts the exact 8-pair 0-based set above (exact count + membership). Passes. EcoRI/HindIII/BamHI/NotI/GCGC/ATAT/12-bp cases all assert exact sequence+position+length and pass.

### Variant/delegate consistency
`GenomicAnalyzer.FindPalindromes` uses the same reverse-complement equality and even-length stepping (`len += 2`), 0-based positions, with `Math.Min(maxLength, seq.Length)` bound (equivalent outcome). Smoke-tested in `GenomicAnalyzerTests.cs`. `RepeatFinder` string vs DnaSequence overloads asserted equal by S1.

### Test quality audit
24 tests in `RepeatFinder_Palindrome_Tests.cs` plus property-based (reverse-complement invariant, even length, bounds) and a snapshot test. Assertions check exact sourced values (restriction-site sequences, Rosalind position/length pairs), edge cases (empty → empty, no-palindrome poly-A → empty, null → throw, odd/too-small/inverted length → throw), overlapping palindromes, entire-sequence palindrome, and the 0-based contract. Not tautological; deterministic.

### Findings / defects
None. Code faithfully realises the validated description.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (0-based reporting vs Rosalind 1-based; minLength floored at 4 — both documented/sourced).
- **Stage B:** PASS.
- **State:** CLEAN — no defect; no code changes.
- **Tests:** palindrome filter 41/41 pass; full `Seqeron.Genomics.Tests` 4461 passed, 0 failed (baseline preserved).
