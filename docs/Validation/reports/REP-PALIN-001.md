# Validation Report: REP-PALIN-001 — DNA Palindrome Detection (reverse-complement palindromes / Rosalind REVP)

- **Validated:** 2026-06-24   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindPalindromes(DnaSequence, minLength, maxLength)`; overload `FindPalindromes(string, ...)`; alternate `GenomicAnalyzer.FindPalindromes(...)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Palindromic sequence:** a (single-stranded) nucleotide sequence is a palindrome iff it equals its reverse complement. A DNA palindrome is NOT a literal character palindrome. EcoRI example 5'-GAATTC-3' / 3'-CTTAAG-5'.
- **Rosalind — REVP (Locating Restriction Sites):** a reverse palindrome equals its reverse complement; report **every** such substring with **length 4–12 inclusive**; **1-based** positions.
- **Wikipedia — Restriction enzyme:** Type II enzymes recognize palindromic sites, typically 4–8 bp (EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT, NotI GCGGCCGC).

### Reverse-complement definition (load-bearing)
A DNA palindrome P satisfies `P == ReverseComplement(P)`. Each position i must complement position (len−1−i); no base is its own complement (A↔T, G↔C), so a contiguous reverse-complement palindrome must have **even length**. Matches the spec's "Why Even Length Only" derivation. Interrupted/spacer palindromes (out of REVP scope) are not claimed here; the even-length rule for the contiguous case is correct.

### Independent cross-check (hand-computed Rosalind REVP sample)
Input `TCAATGCATGCGGGTCTATATGCAT` (len 25). Independent Python (rev-comp equality, even lengths 4–12) reproduces Rosalind's published output **exactly** (1-based position, length, substring):

| 1-based pos | length | substring | 0-based pos (impl) |
|---|---|---|---|
| 4 | 6 | ATGCAT | 3 |
| 5 | 4 | TGCA | 4 |
| 6 | 6 | GCATGC | 5 |
| 7 | 4 | CATG | 6 |
| 17 | 4 | TATA | 16 |
| 18 | 4 | ATAT | 17 |
| 20 | 6 | ATGCAT | 19 |
| 21 | 4 | TGCA | 20 |

The implementation's expected 0-based set `{(3,6),(4,4),(5,6),(6,4),(16,4),(17,4),(19,6),(20,4)}` is exactly Rosalind's output minus 1.

Restriction-site standalone checks (all `revcomp(s)==s` → True): GAATTC, GGATCC, AAGCTT, GCGGCCGC, GCGC, ATAT, CCCGGG, CAGCTG, GATATC. Negative: `revcomp(GAATTA)=TAATTC ≠ GAATTA` (not a palindrome).

### Findings / divergences (notes)
- **Coordinate base:** Rosalind is 1-based; this implementation reports **0-based**. Documented in spec and asserted by tests; consistent with library. Note, not a defect.
- **minLength ≥ 4:** implementation rejects `minLength < 4`. A 2-bp rev-comp palindrome (AT, GC, …) is real, but REVP and the spec scope the minimum to 4 (smallest biologically relevant site). Documented, sourced restriction — PASS-WITH-NOTES, not a defect.
- **Cosmetic substring-label swap (TestSpec + a test comment):** the TestSpec table (lines 146–147) labels 1-based pos 17 as `ATAT` and pos 18 as `TATA`, but the actual sequence has pos 17=`TATA`, pos 18=`ATAT` (verified char-by-char). Rosalind only scores (position, length) pairs, which are correct, so the asserted values are right; only the human-readable substring annotations are swapped. Minor doc nit, no functional impact.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:553-608` — overloads + `FindPalindromesCore`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:178-195` — alternate.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:149-160` — `GetReverseComplementString`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:251-270` — `GetComplementBase` (A↔T, G↔C + full IUPAC).

### Formula realised correctly
Core loop steps `len` from minLength to maxLength by **2** (even-length only); each window compares `candidate == GetReverseComplementString(candidate)` — true reverse-complement equality, NOT literal string reversal. Positions are the 0-based window start. Validation: `minLength` even and ≥ 4; `maxLength ≥ minLength`; null `DnaSequence` throws `ArgumentNullException`; empty/empty-string yields empty. String overload normalizes via `ToUpperInvariant()` (case-insensitive). Both `RepeatFinder` overloads share identical validation.

### Cross-verification recomputed vs code
Exercised by `FindPalindromes_RosalindSample_CorrectOutput` (RepeatFinder_Palindrome_Tests.cs:498), asserting the exact 8-pair 0-based set above (exact count + membership). Passes. EcoRI/HindIII/BamHI/NotI/GCGC/ATAT/12-bp cases assert exact sequence+position+length and pass.

### Variant/delegate consistency
`GenomicAnalyzer.FindPalindromes` uses the same rev-comp equality and even-length stepping (`len += 2`), 0-based positions, with `Math.Min(maxLength, seq.Length)` bound (equivalent outcome). It does NOT validate parameters (no minLength even/≥4 or maxLength check) — acceptable for a smoke-tested alternate; outputs agree with the canonical on valid inputs. String vs DnaSequence overloads asserted equal by S1.

### Test quality audit
103 tests under the Palindrome filter (24 named MUST/SHOULD/COULD + property-based + snapshot + related). Assertions check exact sourced values (restriction-site sequences, Rosalind pos/length pairs), edge cases (empty → empty, no-palindrome → empty, null → throw, odd/too-small/inverted-length → throw), overlapping palindromes, entire-sequence palindrome, 0-based contract. Not tautological; deterministic. (The Rosalind test's inline comment annotations at lines 510–511 carry the same cosmetic substring-label swap noted in Stage A; the asserted pairs are correct.)

### Findings / defects
None functional. Code faithfully realises the validated description.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (0-based vs Rosalind 1-based; minLength floored at 4; cosmetic substring-label swap in spec/test comment — all documented/sourced).
- **Stage B:** PASS.
- **State:** CLEAN — no defect; no code changes.
- **Tests:** Palindrome filter 103/103 pass; build succeeds. No code touched.
