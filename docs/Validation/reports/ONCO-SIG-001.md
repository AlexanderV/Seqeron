# Validation Report: ONCO-SIG-001 — SBS-96 Single-Base-Substitution Trinucleotide Context Catalog

- **Validated:** 2026-06-16   **Area:** Oncology / Mutational Signatures
- **Canonical method(s):** `OncologyAnalyzer.ClassifySbsContext(char,char,char,char)`,
  `OncologyAnalyzer.EnumerateSbs96Channels()`, `OncologyAnalyzer.Build96ContextCatalog(...)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:2113-2261`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (test-coverage gap found and fixed in-session; no code/description defect)

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **COSMIC SBS96** (https://cancer.sanger.ac.uk/signatures/sbs/sbs96/) — WebFetched 2026-06-16.
   Verbatim: "the six substitution subtypes: C>A, C>G, C>T, T>A, T>C, and T>G, as well as the
   nucleotides immediately 5' and 3' to the mutation." Each substitution "is designated by the
   pyrimidine base of the mutated Watson-Crick pair." 96 = 6 substitutions × 4 5'-bases × 4 3'-bases.
2. **Wikipedia "Mutational signatures"** (https://en.wikipedia.org/wiki/Mutational_signatures) —
   WebFetched 2026-06-16. Verbatim: "The G>T substitution is considered equivalent to the C>A
   substitution because it is not possible to differentiate on which DNA strand (forward or reverse)
   the substitution initially occurred." Purine substitutions (G>A, G>C, G>T, A>T, A>G, A>C) are
   counted as their pyrimidine equivalents on the opposite strand via reverse complementation;
   96 = 4 × 6 × 4.
3. **WebSearch** ("SBS96 … pyrimidine reverse complement folding trinucleotide") surfaced COSMIC
   SBS96/SBS6 and confirmed the same six-subtype / 6×4×4 statements independently.

These corroborate the repo's cited Alexandrov 2013 / Bergstrom 2019 / COSMIC definitions without
relying on the repo's own Evidence/TestSpec.

### Formula check

- Six pyrimidine subtypes C>A, C>G, C>T, T>A, T>C, T>G — matches COSMIC/Wikipedia verbatim.
- Channel = `5'[REF>ALT]3'`, mutated base centred — matches the trinucleotide-context definition.
- Purine-reference fold = reverse-complement of context AND complement of substitution
  (complement map A↔T, C↔G) — matches the strand-equivalence rule (Wikipedia "G>T ≡ C>A").
- 96 channels = 6 × 4 × 4 — matches all sources.

### Edge-case semantics check

- Purine reference (A/G) → fold to pyrimidine strand: defined and sourced.
- ref == alt → not a substitution → reject: definitionally correct.
- Non-ACGT base → no defined context → reject: definitionally correct.
- Empty input → all 96 channels present, counts 0 (partition of the empty set): correct.

### Independent cross-check (hand-computed from the fold rule, traced to external sources)

| 5' | Ref | Alt | 3' | Fold step (revcomp context, complement substitution) | Expected channel |
|----|-----|-----|----|------------------------------------------------------|------------------|
| A | C | A | A | pyrimidine, none | A[C>A]A |
| T | C | T | G | pyrimidine, none | T[C>T]G |
| G | T | C | A | pyrimidine, none | G[T>C]A |
| T | G | T | A | revcomp(TGA)=TCA, G>T→C>A | T[C>A]A |
| C | A | G | T | revcomp(CAT)=ATG, A>G→T>C | A[T>C]G |
| G | G | C | C | revcomp(GGC)=GCC, G>C→C>G | G[C>G]C |
| A | A | T | A | revcomp(AAA)=TTT, A>T→T>A | T[T>A]T |
| A | G | A | C | revcomp(AGC)=GCT, G>A→C>T | G[C>T]T |
| G | A | C | T | revcomp(GAT)=ATC, A>C→T>G | A[T>G]C |

All nine rows hand-derived independently and reproduced exactly by the implementation. Rows 8–9
(the two purine folds the existing tests omitted) added as new locked tests.

### Findings / divergences

None. The description (algorithm doc + Evidence + TestSpec) is biologically and mathematically
correct and matches the independently retrieved external sources.

## Stage B — Implementation

### Code path reviewed

`OncologyAnalyzer.cs:2113-2261` — `Sbs96ChannelCount=96`, `PyrimidineSubstitutions` (the six
ordered pairs), `ClassifySbsContext` (normalise → reject ref==alt → fold if purine →
`5'[REF>ALT]3'`), `EnumerateSbs96Channels` (substitution-major × 5' × 3'), `Build96ContextCatalog`
(zero-init all 96, tally), `Complement` (A↔T, C↔G), `NormalizeBase` (upper-case, reject non-ACGT).

### Formula realised correctly?

Yes. The purine fold computes `foldedFive = Complement(threePrime)`, `foldedThree =
Complement(fivePrime)`, `reference = Complement(reference)`, `alternate = Complement(alternate)` —
exactly the reverse-complement of the trinucleotide plus the complemented substitution. Traced for
all nine cross-check rows; every channel matches the externally-sourced value.

### Variant/delegate consistency

`Build96ContextCatalog` classifies via `ClassifySbsContext` and keys on `EnumerateSbs96Channels`,
so the catalog keys are exactly the 96 canonical labels and every variant lands in the same channel
the standalone classifier returns (confirmed by M10/co-count test). Σ counts = #variants (partition).

### Test quality audit (HARD gate)

- **Sourced expectations:** every channel assertion is an exact `Is.EqualTo("…")` traced to the
  COSMIC/Wikipedia six-subtype set and the reverse-complement fold rule — not code echoes. M8
  asserts the full 6-subtype set and 96 distinct pyrimidine-ref channels (a wrong-strand or
  non-folding impl fails).
- **No green-washing:** no Greater/AtLeast/range, no widened tolerance, no skipped/ignored test.
- **Coverage gap found:** the original 16 tests folded only 4 of the 6 purine substitutions
  (G>T, A>G, G>C, A>T) — **G>A→C>T** and **A>C→T>G** were untested; and `ClassifySbsContext`
  validation tested only an invalid *reference* base, not an invalid *alternate* or *3'-flank* base.
  **Fixed in-session** by adding four exact-value tests (`PurineGtoA`→G[C>T]T, `PurineAtoC`→A[T>G]C,
  `InvalidAlternateBase`, `Invalid3PrimeBase`), all values hand-derived from the external fold rule.
  Fixture 16 → 20 [Test] methods.
- **Honest green:** full unfiltered `dotnet test` = **6637 passed, 0 failed**; `dotnet build` 0 errors.
  The changed test file builds warning-free (the 4 build warnings are pre-existing NUnit2007 in the
  unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects

No algorithm defect and no description defect. One test-coverage gap (two purine-fold substitution
branches + two validation branches untested), completely fixed this session. Logged as FINDINGS A45.

## Verdict & follow-ups

- **Stage A: PASS** — description independently confirmed correct against COSMIC SBS96 and Wikipedia
  Mutational signatures (six pyrimidine subtypes, reverse-complement strand folding, 96 = 6×4×4).
- **Stage B: PASS-WITH-NOTES** — implementation realises the validated formula exactly; the only
  issue was a test-coverage gap (now closed with sourced exact-value tests).
- **End-state: CLEAN** — gap fully fixed; build + full suite green (6637/0).
- **Test-quality gate: PASS** (after fix).
