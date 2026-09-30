# Evidence Artifact: REP-PALIN-001

**Test Unit ID:** REP-PALIN-001
**Algorithm:** Reverse-complement (restriction-site) palindrome detection — every window of even length in
[minLength, maxLength] equal to its reverse complement, A/C/G/T only (Rosalind REVP convention)
**Date Collected:** 2026-09-30 (campaign 2026-09, batch B04, findings F12–F13)

---

## Sources opened

### Rosalind REVP — "Locating Restriction Sites"

rosalind.info is blocked for curl; the statement was read from WebSearch snippets and from the problem
docstring reproduced in a public solution
(https://raw.githubusercontent.com/zonghui0228/rosalind-solutions/master/code/rosalind_revp.py, fetched):

- "A DNA string is a reverse palindrome if it is equal to its reverse complement. For instance, GCATGC is a
  reverse palindrome because its reverse complement is GCATGC."
- "Given: A DNA string of length at most 1 kbp in FASTA format."
- "Return: The position and length of **every** reverse palindrome in the string having length between 4 and 12.
  You may return these pairs in any order."
- Sample dataset `TCAATGCATGCGGGTCTATATGCAT` → sample output (1-based position, length):
  `4 6, 5 4, 6 6, 7 4, 17 4, 18 4, 20 6, 21 4` (listed by position, then length).

The fetched solution enumerates `i` then `j ∈ 4..12` (position-then-length order) and its complement function
maps only A/C/G/T (any other symbol is dropped, so a window with a non-ACGT symbol can never equal its reverse
complement).

### EMBOSS `palindrome.c` (raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/palindrome.c, fetched)

Different algorithm (maximal stems with a loop — used by REP-INV-001, not this unit). Relevant point for the
alphabet: the stem loop tracks `alln` and rejects a stem whose paired bases are all `n` (`!alln`): N–N symbolic
"complementarity" is explicitly not accepted as a real pairing. EMBOSS `einverted` scores a match only for
a/c/g/t (REP-INV-001 Evidence).

### IUPAC ambiguity codes (NC-IUB 1984, via the canonical `SequenceExtensions.GetComplementBase`)

N, S, W and R/Y pairs are symbolically self-/mutually complementary (N→N, S→S, W→W, R↔Y), so `NNNN`, `SSSS`,
`RYRY` equal their symbolic reverse complement. They denote *sets* of bases: `SS` ∈ {CC, CG, GC, GG}, of which
only CG/GC are palindromic, so symbolic equality does not make the underlying sequence a palindrome.
Non-IUPAC symbols pass through `GetComplementBase` unchanged (`A--T`, `1111` were "palindromes").

## Even length

For a window of odd length 2k+1 the middle base would have to equal its own complement; no base in {A,C,G,T} does
(A↔T, C↔G). Hence only even lengths exist; an odd `maxLength` is equivalent to `maxLength − 1`.

## Reference cross-check (2026-09-30)

Reference = Python, Biopython 1.88 `Seq(w).reverse_complement() == w` over every even window whose symbols are
all in ACGT (after upper-casing); order (position, length).

| Check | Cases | Mismatches |
|---|---|---|
| C# (new) vs Biopython reference — random 0–120 bp, alphabets ACGT/AT/GC/ACGTN/mixed case/ACGTNRYSWU-/ACGTS, minLength 4–10, maxLength up to min+200 | 6000 (2563 non-empty, 18 378 palindromes) | 0 |
| C# vs reference, 20 kb random ACGT (4–40) and 5 kb AT-only (4–100) | 2 (1656 / 2546 palindromes) | 0 |
| zonghui0228 REVP solution vs Biopython reference, random ACGT 0–200 bp | 2000 | 0 |
| REVP sample (reference, +1) | 1 | 0 — `4 6, 5 4, 6 6, 7 4, 17 4, 18 4, 20 6, 21 4` |

Old implementation (before F12/F13): `NNNNNN` → 4 hits (NNNN×3, NNNNNN), `ANNT`, `SSSS`, `RYRY`, `A--T`, `1111` →
1 hit each; `FindPalindromes("ACGT", 4, int.MaxValue)` → `ArgumentOutOfRangeException` from `Substring` after the
`len += 2` counter overflowed; output ordered by length then position. On pure-ACGT input the result *sets* were
already correct.

Timing (Release, random ACGT): 1 Mb, 4–12: 646 ms → 72 ms; 100 kb, 4–200: 3609 ms → 37 ms.
