# Validation Report: KMER-GENERATE-001 — K-mer Generation (enumerate all possible k-mers over an alphabet)

- **Validated:** 2026-06-16   **Area:** K-mer
- **Canonical method(s):** `KmerAnalyzer.GenerateAllKmers(int k, string alphabet = "ACGT")` (+ private `GenerateKmersRecursive`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent retrieval)

| # | Source | URL | What it confirms (verbatim) |
|---|--------|-----|------------------------------|
| 1 | Wikipedia — K-mer | https://en.wikipedia.org/wiki/K-mer | Definition: *"In bioinformatics, k-mers are substrings of length k contained within a biological sequence."* Universe size: *"There exist n^k total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA)."* (⇒ 4^k for DNA). Worked example AGAT → monomers A,G,A,T; 2-mers AG,GA,AT; 3-mers AGA,GAT; 4-mer AGAT (confirms a k-mer is a contiguous length-k string over {A,C,G,T}). |
| 2 | Python Std Library — itertools.product | https://docs.python.org/3/library/itertools.html | *"product(A, repeat=4) means the same as product(A, A, A, A)."* (k-fold Cartesian product = all length-k tuples). Ordering: *"The nested loops cycle like an odometer with the rightmost element advancing on every iteration. This pattern creates a lexicographic ordering so that if the input's iterables are sorted, the product tuples are emitted in sorted order."* Worked example: `product(range(2), repeat=3) → 000 001 010 011 100 101 110 111`. |
| 3 | Rosalind LEXF — Enumerating k-mers Lexicographically | https://rosalind.info/problems/lexf/ | Canonical "enumerate all length-n strings over an ordered alphabet, in lexicographic order" problem. Sample: alphabet `A C G T`, n=2 → exactly `AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT` (16 strings). Independently fixes the M2 expected list and lexicographic order. |

(Source 2 BioInfoLogics from the repo Evidence — "the possible combinations of k positions are computed as 4^k" — restates the same n^k fact; not re-fetched but the n^k claim is independently established by Sources 1 and 3.)

### Formula check
- **Universe size = n^k** (4^k for DNA): confirmed verbatim by Source 1; cross-checked by Source 3 (16 = 4^2 for the sample) and by an independent `itertools.product` run (counts 4,16,64,256,1024,4096 for k=1..6).
- **Enumeration = k-fold Cartesian product Σ^k**: confirmed by Source 2 (`product(A, repeat=k)`).
- **Order (sorted alphabet ⇒ lexicographic, rightmost position fastest / odometer)**: confirmed verbatim by Source 2 and demonstrated by Source 3's LEXF sample.

### Edge-case semantics
- **k ≤ 0**: a k-mer is a string *"of length k"* (Source 1); k≤0 has no valid k-mer length ⇒ error is the defined behaviour. (Convention; sourced from the definition.)
- **Single-letter alphabet**: 1^k = 1 (one homopolymer). Confirmed: `product('A', repeat=4) → ['AAAA']`.
- **Unsorted alphabet**: still yields all n^k strings, but in the alphabet's positional order — lexicographic guarantee holds *only if the alphabet is sorted* (Source 2). Confirmed: `product('TGCA', repeat=1) → ['T','G','C','A']`.
- **Null/empty alphabet**: no symbols ⇒ no k-mers can be formed ⇒ error (derived from contract; standard).

### Independent cross-check (exact numbers, retrieved/recomputed this session)
`python3 itertools.product` (reference implementation):
- counts k=1..6 over `ACGT` → 4, 16, 64, 256, 1024, 4096 ✓ (matches M4, Source 1's 4^k)
- k=2 list → `AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT` ✓ (matches M2 **and** Rosalind LEXF sample output)
- k=3 → first `AAA`, second `AAC`, last `TTT`; count 64 ✓ (matches M3)
- protein 20-letter alphabet, k=2 → 400 ✓ (matches M5)
- single-letter `A`, k=4 → `['AAAA']` ✓ (matches S1)
- unsorted `TGCA`, k=1 → `['T','G','C','A']` ✓ (matches C1)

### Findings / divergences
None. Every formula, ordering rule, and edge-case in the TestSpec/Evidence/doc traces to an authoritative source retrieved this session. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs:299–325`.

```
public static IEnumerable<string> GenerateAllKmers(int k, string alphabet = "ACGT")
{
    if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
    if (string.IsNullOrEmpty(alphabet)) throw new ArgumentException("Alphabet cannot be empty.", nameof(alphabet));
    return GenerateKmersRecursive("", k, alphabet);
}
private static IEnumerable<string> GenerateKmersRecursive(string prefix, int k, string alphabet)
{
    if (prefix.Length == k) { yield return prefix; yield break; }
    foreach (char c in alphabet)
        foreach (var kmer in GenerateKmersRecursive(prefix + c, k, alphabet))
            yield return kmer;
}
```

### Formula realised correctly? (evidence)
- The recursion extends a prefix one character per level, iterating the alphabet **in order** at each of the k positions — exactly the k-fold Cartesian product Σ^k. The leftmost position is the outermost loop and the rightmost (deepest) advances fastest ⇒ **odometer ordering**, matching Source 2. For a sorted alphabet this is lexicographic (INV-04); `ACGT` is sorted, so DNA output is `AAA, AAC, …, TTT`. ✓
- Count is n^k by construction (k nested n-way branches). ✓ Distinctness (INV-02): each path through the recursion is a distinct length-k tuple ⇒ no duplicates for a duplicate-free alphabet. ✓
- **Validation is eager**, not deferred: `GenerateAllKmers` uses `return` (not `yield`), so `k<=0` / null-or-empty checks fire on the *call*, before any enumeration. The S3/S4 tests wrap `GenerateAllKmers(...).ToList()` in the throwing lambda, so they pass either way; behaviour is correct and robust.

### Cross-verification table recomputed vs code (full suite run, actual code)
| Case | Sourced expected | Code (suite) |
|------|------------------|--------------|
| M1 k=1 ACGT | A,C,G,T (4) | ✓ |
| M2 k=2 ACGT | 16 two-mers AA..TT (lexicographic) — Rosalind LEXF | ✓ |
| M3 k=3 ACGT | 64; first AAA, second AAC, last TTT | ✓ |
| M4 k=1..6 | 4,16,64,256,1024,4096 | ✓ |
| M5 protein k=2 | 400 | ✓ |
| M6 k=4 | 256; distinct; = 4-fold Cartesian set | ✓ |
| S1 "A" k=4 | AAAA (1) | ✓ |
| S2 | every k-mer len 2, alphabet chars only | ✓ |
| S3 k=0/-1 | ArgumentOutOfRangeException | ✓ |
| S4 null/"" | ArgumentException | ✓ |
| C1 "TGCA" k=1 | T,G,C,A (alphabet order) | ✓ |

### Variant/delegate consistency
The only other public entry point, MCP `AnalysisTools.GenerateAllKmers` (`src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.cs:100`), is a thin pass-through to `KmerAnalyzer.GenerateAllKmers` with no added logic — fully covered by the canonical tests. No `*Fast`/instance variant exists.

### Test quality audit (against the hard gate)
- **Sourced, not code-echoes:** M4 derives the expected count from the 4^k formula (`Math.Pow(4,k)`), not from code output; M2/M3/S1/C1 assert exact ordered lists/elements taken from Source 2 / Rosalind LEXF; M6 builds an **independent** 4-fold Cartesian product in LINQ and asserts set-equality. M2 (exact order), M3 (first/second/last), M6 (set-equality), and C1 (alphabet-order) would each **fail** against a wrong-order or missing/extra-element implementation — none are green-washable.
- **No green-washing:** exact equalities everywhere a value is known; no Greater/AtLeast/range substituted for a known value; no widened tolerance; no skipped/ignored test.
- **Coverage:** all Stage-A branches exercised — happy path (M1–M6), single-letter degenerate (S1), INV-03 length/membership (S2), both error branches k≤0 (S3) and null/empty alphabet (S4), and the unsorted-order branch (C1). Both public k-input error paths and both alphabet error inputs (empty and null) are covered.
- **Honest green:** full unfiltered suite **Passed: 6607, Failed: 0** (the single MFE benchmark is a pre-existing explicit `[Explicit]` skip, unrelated). Build 0 errors; the 4 build warnings are pre-existing NUnit2007 warnings in the unrelated `ApproximateMatcher_EditDistance_Tests.cs` — the KMER test/impl files are warning-free.

**Test-quality gate: PASS.** No defect; no test change required.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS.** Description (TestSpec, Evidence, algorithm doc) is biologically/mathematically correct and fully sourced; n^k universe, Cartesian-product enumeration, sorted-alphabet lexicographic/odometer order all confirmed against Wikipedia, itertools.product, and Rosalind LEXF.
- **Stage B: PASS.** Code faithfully realises Σ^k with odometer ordering; all 11 tests assert exact sourced values and cover every Stage-A branch.
- **End-state: ✅ CLEAN.** No code or test change needed; algorithm fully functional. Full unfiltered suite green (6607/0).
