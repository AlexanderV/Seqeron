# Validation Report: KMER-BOTH-001 — Both-Strand K-mer Counting (forward + reverse complement)

- **Validated:** 2026-06-16   **Area:** K-mer
- **Canonical method(s):** `KmerAnalyzer.CountKmersBothStrands(string, int)`; delegate `CountKmersBothStrands(DnaSequence, int)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (and what they confirm)

1. **kPAL Methodology** — https://kpal.readthedocs.io/en/latest/method.html (WebFetch, 2026-06-16).
   Verbatim: *"kPAL can forcefully balance the k-mer profiles (if desired) by adding the values of
   each k-mer to its reverse complement."* Confirms the **additive ("balance") both-strand**
   operation: the both-strand value of a k-mer is the sum of its own count and that of its reverse
   complement.

2. **Shporer, Chor, Rosset, Horn (2016), BMC Genomics — Inversion symmetry of DNA k-mer counts** —
   https://pmc.ncbi.nlm.nih.gov/articles/PMC5006273/ (WebFetch, 2026-06-16).
   Verbatim: *"Inversion Symmetry (IS): the counts of a k-mer of nucleotides on a chromosomal strand
   are almost equal to those of its inverse (reverse-complement) string."* and *"the number of times
   a string of nucleotides of length k is observed on a strand, when read from 5' to 3', is almost
   equal to the number of times it is observed on the other strand when the latter is read from its
   5' end to 3' end."* This grounds the identity **count[w] = forward[w] + forward[RC(w)]**: counting
   w on the RC strand (read 5'→3') = counting RC(w) on the forward strand.

3. **Mash issue #45 (Ondov et al.)** — https://github.com/marbl/Mash/issues/45 (WebFetch, 2026-06-16).
   Verbatim: *"only the lexicographically smaller of the forward and reverse complement
   representations of a k-mer is hashed."* This is the **canonical-collapsing** convention (single
   lex-smaller key). It is explicitly NOT what KMER-BOTH-001 implements.

4. **KMC 2 / Jellyfish `-C` semantics** — WebSearch (2026-06-16), summarizing KMC2 (Oxford
   Bioinformatics 31(10):1569) and the Jellyfish manual: *"With the option -C, only the canonical
   representation of the mers are stored in the hash and the count value is the number of occurrences
   of both the mer and its reverse-complement."* Key cross-check: the canonical **value** per pair is
   forward[w]+forward[RC(w)] — identical to the additive value — the only difference is the **key
   set** (one lex-smaller key for canonical vs. both keys for additive).

### Formula check

- Model (doc §2.2, code XML-doc): `count[w] = forward[w] + forward[RC(w)]`, equivalently count k-mers
  of S and of RC(S) and sum per key. Matches Shporer 2016 (inversion symmetry) + kPAL balance exactly.
- Window count / grand total: a length-L sequence has L−k+1 overlapping k-mers per strand, so the
  both-strand total over all keys is 2·(L−k+1) (Marçais & Kingsford 2011 k-mer-occurrence definition).

### Definitions & conventions

- Coordinate base / overlap: overlapping windows i = 0..L−k (0-based), L−k+1 windows — standard.
- Case-insensitive (upper-cased internally); reverse complement via repository `GetComplementBase`
  (IUPAC-complete, verified table A↔T, G↔C, plus ambiguity codes). Standard.
- Additive vs canonical clearly distinguished in the doc §2.5 and the code XML-doc.

### Edge-case semantics (sourced/derived)

- k > L (L−k+1 ≤ 0) ⇒ no windows ⇒ empty (Marçais & Kingsford occurrence definition).
- Empty / null ⇒ empty (consistent with the window formula and sibling `CountKmers`; ASSUMPTION,
  boundary-only, non-correctness-affecting).
- k = L ⇒ one window per strand.
- Palindromic k-mer (RC(w)=w) ⇒ both contributions land on one key ⇒ count = 2·forward[w] (direct
  consequence of inversion symmetry).
- k ≤ 0 ⇒ ArgumentOutOfRangeException (API-shape assumption, matches sibling methods).

### Independent cross-check (hand-computed from the sourced identity, NOT from code)

Worked example S = ATGGC, k = 2. Forward 2-mers {AT,TG,GG,GC} each ×1. RC(ATGGC) = GCCAT, RC-strand
2-mers {GC,CC,CA,AT}. Per-key via inversion-symmetry identity count[w]=forward[w]+forward[RC(w)]:

| w | forward[w] | RC(w) | forward[RC(w)] | count[w] |
|---|-----------|-------|----------------|----------|
| AT | 1 | AT | 1 | **2** |
| TG | 1 | CA | 0 | **1** |
| GG | 1 | CC | 0 | **1** |
| GC | 1 | GC | 1 | **2** |
| CC | 0 | GG | 1 | **1** |
| CA | 0 | TG | 1 | **1** |

⇒ {AT:2, TG:1, GG:1, GC:2, CC:1, CA:1}, Σ = 8 = 2·(5−2+1). AT and GC are RC-palindromes, so they
exercise INV-04 (doubling on one key). Independently matches the spec's M1 dictionary and the doc
worked example.

Palindrome word check S = ACGT, k = 2: RC(ACGT)=ACGT ⇒ RC-strand 2-mers {AC,CG,GT} = forward ⇒ each
key doubles ⇒ {AC:2,CG:2,GT:2}. Non-palindrome S = AAA, k = 2: forward {AA:2}; RC(AAA)=TTT ⇒ {TT:2} ⇒
{AA:2,TT:2}. k=L S = ATGC, k=4: RC(ATGC)=GCAT ⇒ {ATGC:1, GCAT:1}. All match the spec.

### Findings / divergences (Stage A)

- The session prompt's parenthetical loosely equated this unit with the Jellyfish/KMC `--canonical`
  *collapsing* mode. The spec/doc/code deliberately and correctly implement the **non-collapsing
  additive (kPAL "balance")** variant instead. This is NOT an error: (a) both are legitimate
  both-strand semantics; (b) the per-k-mer **value** is identical to the canonical count
  (forward[w]+forward[RC(w)]) — only the key set differs (both keys vs one lex-smaller key); (c) the
  method name "BothStrands", the registry ("Forward + reverse complement"), and the doc §2.5 / §5.3 all
  consistently denote the additive view and explicitly note canonical collapsing as a separate,
  not-implemented variant. The description is mathematically and biologically correct. **PASS.**

## Stage B — Implementation

### Code path reviewed

- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs:476-490` —
  `CountKmersBothStrands(string, int)`: forward = `CountKmers(sequence, k)`; revComp =
  `CountKmers(GetReverseComplementString(sequence ?? ""), k)`; merge by summing per key.
- `:501-502` — `CountKmersBothStrands(DnaSequence, int)` delegates to the string overload.
- `:20-42` — `CountKmers`: null/empty ⇒ empty; k≤0 ⇒ ArgumentOutOfRangeException; k>L ⇒ empty;
  upper-cases; L−k+1 overlapping windows.
- `DnaSequence.cs:149-160` `GetReverseComplementString` (reverse + per-base complement);
  `SequenceExtensions.cs:138-157` `GetComplementBase` (IUPAC-complete, verified).

### Formula realised correctly?

Yes. Counting k-mers of RC(S) yields, for key u, the number of substrings of RC(S) equal to u, which
equals forward[RC(u)]. Merging forward[w] + (count in RC(S))[w] = forward[w] + forward[RC(w)] — exactly
the sourced inversion-symmetry identity. Null is coerced to "" before reverse-complement, so the null
path is safe. k≤0 throws via the inner `CountKmers` (forward pass runs first).

### Cross-verification table recomputed vs code (full suite executed)

| Case | Expected (sourced/hand) | Code result | Match |
|------|------------------------|-------------|-------|
| ATGGC k=2 | {AT:2,TG:1,GG:1,GC:2,CC:1,CA:1} | same (test M1) | ✅ |
| ACGT k=2 | {AC:2,CG:2,GT:2} | same (M2) | ✅ |
| AAA k=2 | {AA:2,TT:2} | same (M3) | ✅ |
| ATGGC total | 8 = 2·(5−2+1) | 8 (M4) | ✅ |
| ATGC k=4 | {ATGC:1,GCAT:1} | same (S2) | ✅ |
| empty / null / k>L | {} | {} (C1/C2/C3) | ✅ |
| k=0 | throws AOORE | throws (C4) | ✅ |

### Variant/delegate consistency

`DnaSequence` overload == string overload (test S3, `Is.EquivalentTo`). Case-insensitivity holds
(S1). ✅

### Numerical robustness

Integer counts; no division, overflow only at astronomically large sequence lengths (int counts), not
in scope. Null-safe via `?? string.Empty`. ✅

### Test quality audit (HARD gate)

Canonical test file: `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_CountKmersBothStrands_Tests.cs`
(13 tests).

- **Sourced, not code-echoes:** M1/M2/M3/S2 assert **exact** dictionaries (`Is.EquivalentTo`, no extras)
  hand-derived from the inversion-symmetry identity in the Evidence/this report — they would FAIL
  against a wrong implementation (e.g. canonical collapsing, missing RC strand, or wrong complement).
  M4 pins the sourced grand-total 2·(L−k+1). M5 (INV-01) and M6 (INV-03) add invariant cross-checks;
  M5 reuses `CountKmers`+`GetReverseComplementString` so it is a partial cross-check, but the absolute
  values are independently anchored by M1/M2/M3/S2, so the suite is not green-washable on values.
- **No green-washing:** exact equality everywhere a value is known; no weakened asserts, widened
  tolerances, skips, or expected-value-to-output adjustments.
- **Cover all logic:** both public overloads exercised (string + DnaSequence); all Stage-A branches —
  non-palindrome (M3), whole-word RC symmetry (M2), RC-palindrome doubling on one key (M1: AT, GC →
  INV-04), k=L boundary (S2), grand-total invariant (M4), strand symmetry (M6); and all edge/error
  cases (empty C1, null C2, k>L C3, k≤0 throws C4). INV-05 (positive integer counts) is implied by the
  exact-dictionary assertions.
- **Honest green:** full unfiltered suite `dotnet test` = **Failed: 0, Passed: 6607**; `dotnet build`
  0 errors (4 pre-existing warnings, all in unrelated `ApproximateMatcher_EditDistance_Tests.cs`, none
  in changed files — no files changed this session).

**Gate result: PASS** — tests validate the algorithm against external sources, not against the code.
No weak/code-echoing/green-washed/partial defects requiring a fix.

### Findings / defects (Stage B)

None. Implementation faithfully realises the validated additive both-strand identity; delegate
consistent; edge cases handled as sourced. **PASS.**

## Verdict & follow-ups

- **Stage A: PASS.** Description is mathematically and biologically correct; additive (kPAL balance)
  semantics confirmed by kPAL + Shporer 2016; correctly distinguished from canonical collapsing
  (Mash/Jellyfish/KMC). No defect.
- **Stage B: PASS.** Code computes count[w]=forward[w]+forward[RC(w)] exactly; tests lock exact
  sourced dictionaries and all edge cases; full suite green (6607/0).
- **End-state: ✅ CLEAN.** No defect found; no code or test changes required this session.
- No findings logged in FINDINGS_REGISTER (no defect).
