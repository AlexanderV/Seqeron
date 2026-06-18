# Validation Report: RNA-INVERT-001 — RNA Inverted Repeats (potential stem regions)

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.FindInvertedRepeats(string sequence, int minLength = 4, int minSpacing = 3, int maxSpacing = 100)` → `IEnumerable<(int Start1, int End1, int Start2, int End2, int Length)>`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (retrieved, not trusted from citation labels)

1. **Wikipedia "Inverted repeat"** (cites Ussery, Wassenaar & Borini 2008) — https://en.wikipedia.org/wiki/Inverted_repeat
   - Verbatim: "An **inverted repeat** (or **IR**) is a single stranded sequence of nucleotides followed downstream by its reverse complement."
   - Worked example `5'---TTACGnnnnnnCGTAA---3'`; the page derives revcomp step-by-step: TTACG → complement AATGC → reverse **CGTAA**. Confirms right arm = *reverse complement* (antiparallel), not bare complement and not parallel direct repeat.
   - "The intervening sequence of nucleotides … can be any length **including zero**"; zero gap ⇒ palindrome.
2. **IUPACpal, Alamro et al. 2021, BMC Bioinformatics 22:51** — https://pmc.ncbi.nlm.nih.gov/articles/PMC7866733/
   - Verbatim: "A gapped IR is therefore a string that can be expressed in the form **WGW̄ᴿ** for some pair of strings W and G where |G|≥0." (left arm W, gap G, right arm = revcomp(W)).
   - Mismatch model: "gapped inverted repeat within k mismatches when it can be expressed in the form WGW̄ᴿ with **δH(W,W̄ᴿ)≤k**" (Hamming distance). Perfect IR = the k=0 case.
   - Complement basis: "A⟷T and C⟷G" over {A,C,G,T} (RNA: A⟷U, C⟷G).
3. **EMBOSS einverted manual** (Rice et al. 2000) — https://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted
   - Verbatim: "einverted finds inverted repeats (stem loops) … regions of local alignment of the input sequence and its reverse complement that exceed a threshold score." Confirms IR ≡ potential stem-loop; arms complementary; intervening region = loop. Default scoring match +3 / mismatch −4 / gap 12 / threshold 50 (the scored DP variant — explicitly out of scope here; this unit reports only perfect k=0 arms).

### Formula / definition check
The doc/TestSpec model `WGW̄ᴿ` with right arm = reverse complement of left, antiparallel pairing `complement(seq[s2+L−1−k]) == seq[s1+k]`, strict Watson-Crick (A⟷U, C⟷G) — matches IUPACpal and Wikipedia/Ussery verbatim. The k=0 restriction is a *documented special case* of the IUPACpal model, not an invented behavior. Loop bounds [minSpacing, maxSpacing] map to IUPACpal "maximum gap size" / einverted loop; min-arm-length maps to IUPACpal minimum length.

### Edge-case semantics
All sourced/defined: |G|≥0 (zero ⇒ palindrome; default RNA call requires loop ≥ 3, the biological minimal hairpin loop); arms below minLength not reported; parallel direct repeat is *not* an IR (right arm must be the reverse complement); strict WC, no G-U wobble (a documented, conservative scope restriction — over-reporting avoided).

### Independent cross-check (numbers)
Wrote an independent from-scratch reference (`/tmp/ir_ref.py`) with its **own** revcomp function and pure `WGW̄ᴿ` enumeration. Every exact-position expectation reproduced, and every reported right arm verified `== revcomp(left)` by the independent complement:
- M1 `UUACGAAAAAACGUAA` → `(0,4,11,15,5)`; right CGUAA = revcomp(UUACG) ✓ (matches Wikipedia example).
- M2 `GGCCAAAGGCC` → `(0,3,7,10,4)`; right GGCC = revcomp(GGCC) ✓.
- C1 `GGGGCAAAGCCCC` → `(0,4,8,12,5)`; right GCCCC = revcomp(GGGGC) ✓ (maximal arm, not truncated to 4).
- M3 `AAGGAAAAAGG` → ∅ (AAGG is parallel, revcomp(AAGG)=CCUU) ✓; M4 poly-A → ∅ ✓.
- A pure no-policy enumeration of *all* maximal perfect IRs confirms the canonical hits are the genuine maximal IRs (e.g. M1's `(0,3,12,15,4)` is a sub-arm correctly subsumed by the length-5 hit), so the expected values are definition-derived, not artifacts of the greedy reporting policy.

### Findings / divergences
None. Stage A **PASS**.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs:1803-1883` (`FindInvertedRepeats`).
- Complement: `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:175-194` (`GetRnaComplementBase`, strict WC + IUPAC degenerate codes; A⟷U, C⟷G).

### Formula realised correctly?
Yes. For each loop boundary `(p,q)` with loop `q−p−1 ∈ [minSpacing, maxSpacing]`, the stem extends **outward** while `GetComplement(seq[q+len]) == seq[p−len]` — exactly the antiparallel reverse-complement relation `complement(seq[s2+L−1−k]) == seq[s1+k]`. Maximal arm kept per boundary (≥ minLength), then greedy longest-first non-overlapping selection. Guard at :1809 returns empty for null/empty, `minLength<1`, `minSpacing<0`, `maxSpacing<minSpacing`; length guard `n < 2·minLength + minSpacing` at :1815. The historical bug (right-arm offset from `minLength` rather than matched arm length) is fixed — confirmed by C1 reporting length 5.

### Cross-verification table recomputed vs code
My independent reference (mirroring the documented reporting policy) and the actual test run agree on all 13 cases. The actual unit test run (`--filter ~RnaSecondaryStructure_FindInvertedRepeats`) is green at 13/13; values match the externally-sourced expectations above.

### Variant/delegate consistency
Single public method; no `*Fast`/delegate variants. `GetComplement` ≡ `GetRnaComplementBase` (RnaSecondaryStructure.cs:469). Consistent.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1/M2/C1 assert exact tuples via `Is.EquivalentTo` (a wrong implementation returning different positions fails); values trace to Wikipedia/IUPACpal, re-derived independently this session. Not weakenable to greater/contains.
- **No green-washing:** no weakened assertions, no widened tolerances, no skips. Exact tuples where exact values known; ranges (`Is.InRange(3,100)`, `GreaterThanOrEqualTo(4)`) appear **only** in C2, which is by design an INV-03/INV-04 *bound* property test (those invariants are bounds, not exact values), with the exact values already locked by M1/M2/C1.
- **Coverage gap found & fixed:** the documented parameter-guard branches `minLength<1`, `minSpacing<0`, `maxSpacing<minSpacing` (doc/TestSpec §3.1 "→ empty result") were untested. Added **S6** `FindInvertedRepeats_OutOfRangeParameters_ReturnEmpty` exercising all three on an input that *does* contain a real IR under defaults (isolating the guard). All other public branches (null/empty, too-short, parallel-reject, no-IR, minLength/minSpacing/maxSpacing loop filters, maximal extension, invariants) already covered.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6593` (was 6592 + 1 new); `dotnet build` 0 errors. The 4 build warnings are pre-existing NUnit2007 in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`; the changed file builds warning-free.

### Findings / defects
No behavioral defect. One test-coverage gap (untested guard branches) — fixed in-session by adding test S6. Stage B **PASS**.

## Verdict & follow-ups
- **Stage A: PASS.** Definition (`WGW̄ᴿ`, antiparallel reverse complement, k=0) matches IUPACpal, Wikipedia/Ussery, and EMBOSS einverted verbatim.
- **Stage B: PASS.** Code faithfully realises the validated description; all 13 tests green against externally re-derived expected values.
- **Test-quality gate: PASS** (after adding S6 to cover the documented out-of-range-parameter branches).
- **End-state: ✅ CLEAN** — no behavioral defect; the coverage gap was completely fixed; build + full suite green.
- No open follow-ups.
