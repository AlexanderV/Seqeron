# Validation Report: VARIANT-SNP-001 — SNP Detection

- **Validated:** 2026-06-15   **Area:** Variants
- **Canonical method(s):** `VariantCaller.FindSnpsDirect(string reference, string query)` (canonical, positional/Hamming-style); `VariantCaller.FindSnps(DnaSequence reference, DnaSequence query)` (delegate: aligns via `CallVariants`, filters to `VariantType.SNP`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened this session (independent of the repo artifacts)

1. **VCFv4.3 specification** — `https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex` (fetched 2026-06-15).
   Confirmed verbatim:
   - Simple SNP example record: `20  14370  rs6054257  G  A  29  PASS  …` described as "a good simple SNP" — a single-base REF `G` substituted by a single-base ALT `A`.
   - POS field: "The reference position, with the 1st base having position 1." → VCF POS is **1-based**.
   - REF field: "Each base must be one of A,C,G,T,N (case insensitive)." + "Tools processing VCF files are not required to preserve case in the allele Strings."
2. **Hamming distance** — `https://en.wikipedia.org/wiki/Hamming_distance` (fetched 2026-06-15).
   Confirmed verbatim: "The Hamming distance between two equal-length strings of symbols is the number of positions at which the corresponding symbols are different." Defined **only for equal-length strings** (the reference Python raises `ValueError` otherwise). Worked example: `karolin` vs `kathrin` → distance 3.
   This corroborates the repo's cited PMC5410656 definition ("the number of positions that two codewords of the same length differ").

### Formula check
- SNP = single-base substitution (single-base REF ≠ single-base ALT at one position); a position where REF == ALT is a match, not a variant. Matches VCFv4.3 §1.1. ✓
- Positional substitution set over equal-length inputs = the Hamming mismatch set; SNP count = Hamming distance. Matches the Hamming definition (equal-length only). ✓
- Position model: in-memory `Variant.Position` is 0-based; the 1-based VCF POS is produced only by `ToVcfLines` (`Position + 1`, verified in `VariantCaller.cs:367`). Internally consistent with the VCF spec and sibling VARIANT-CALL-001. ✓

### Edge-case semantics check
- Identical equal-length → 0 SNPs (Hamming 0). Sourced. ✓
- Empty / one-empty input → empty (nothing to compare). Reasonable, documented contract. ✓
- Unequal length → compare common prefix `min(|r|,|q|)` only; trailing region is indel territory, out of scope (ASSUMPTION-1). Consistent with "Hamming defined for equal length only". ✓
- Null inputs to `FindSnps` → `ArgumentNullException` (propagated from `CallVariants`). ✓

### Independent cross-check (hand computation against the fetched definitions)
| Input ref / query | Differing 0-based indices | Expected SNPs (ref→alt) | Source basis |
|---|---|---|---|
| `ATGC` / `ATGC` | {} | none (Hamming 0) | Hamming def |
| `ATGC` / `ATTC` | {2} | G→T @2 | Hamming + VCF SNP |
| `AAAA` / `TGTA` | {0,1,2} (idx3 A==A) | A→T, A→G, A→T | Hamming def |
| `ATGCAA` / `ATTC` | common prefix len 4: {2} (idx3 C==C) | G→T @2; trailing "AA" ignored | Hamming equal-length precondition |
| `ACGTACGT` / `TCATACGA` | {0,2,5,7} | count 4 | Hamming def |
| `ATGCATGC` / `ATGAATGC` | {3} | C→A @3 | VCF SNP (single substitution) |

All values match the repo's test datasets.

### Findings / divergences (Stage A)
- The Evidence/TestSpec phrasing "SNP comparison … must be case-insensitive" is broader than the source warrants: VCFv4.3's "case insensitive" governs allele-string *interpretation when parsing VCF*, not a mandate that a positional substitution detector normalize case. The algorithm doc (`SNP_Detection.md` §3.1) scopes this correctly ("A,C,G,T,N (case-insensitive **for classification**)"), and the only method that normalizes case (`ClassifyMutation`, `VariantCaller.cs:189-190`) is **not** a method under this unit. Recorded as a minor wording note; does not affect the Stage-A correctness of the detection formula.

**Stage A: PASS** — biology/maths of SNP detection is correct in the abstract and matches independently fetched authoritative sources.

## Stage B — Implementation

### Code path reviewed
- `FindSnpsDirect` — `VariantCaller.cs:125-144`: null/empty guard → `yield break`; `minLen = Math.Min(...)`; forward scan emitting `Variant(Position=i, REF=ref[i], ALT=query[i], Type=SNP, QueryPosition=i)` at each mismatch. Realises the validated positional/Hamming formula exactly. ✓
- `FindSnps` — `VariantCaller.cs:117-120`: `CallVariants(reference, query).Where(v => v.Type == VariantType.SNP)`. `CallVariants` (`:27-33`) does `ArgumentNullException.ThrowIfNull` on both args, then `CallVariantsCore` aligns (`SequenceAligner.GlobalAlign`) and walks columns (`CallVariantsFromAlignment :41-100`). SNP-only filter is correct; null propagation is correct. ✓

### Formula realised correctly?
- Equality test `reference[i] != query[i]` is the exact per-position substitution rule; no thresholds/heuristics. ✓
- For equal-length inputs the emitted-count equals the number of differing indices = Hamming distance (INV-05). ✓
- `Position`/`QueryPosition` are 0-based; VCF 1-based only via `ToVcfLines` (`:367`). ✓

### Cross-verification table recomputed vs code (all pass in the live suite)
| Test | Input | Asserted result | Matches hand/source value |
|---|---|---|---|
| M1 | `ATGC`/`ATGC` | empty | ✓ |
| M2 | `ATGC`/`ATTC` | 1 SNP pos2 G→T, QueryPos 2, Type SNP | ✓ |
| M3 | `AAAA`/`TGTA` | pos {0,1,2}, REF {A,A,A}, ALT {T,G,T} | ✓ |
| M4 | `AAAA`/`TGTA` | all Type==SNP, REF≠ALT | ✓ |
| M5 | `ATGCAA`/`ATTC` | 1 SNP pos2 G→T (prefix-only) | ✓ |
| C1 | `ACGTACGT`/`TCATACGA` | count == 4 (Hamming), all SNP, REF≠ALT | ✓ |
| M6 | `ATGCATGC`/`ATGAATGC` | 1 SNP pos3 C→A, all SNP | ✓ |
| S1 | `""`/`""` | empty | ✓ |
| S2 | `ATGC`/`""` | empty | ✓ |
| S3 | null ref | `ArgumentNullException` | ✓ |
| S4 | null query | `ArgumentNullException` | ✓ |
| S5 | identical DnaSequence | empty | ✓ |

### Variant/delegate consistency
- `FindSnps` ⊆ `CallVariants` filtered to SNP; `FindSnpsDirect` independent positional scan. For substitution-only equal-length inputs both agree (M6 vs M2/M5 logic). The two paths necessarily diverge only near indels (alignment-dependent positions), which is documented and out of this unit's scope. ✓

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M2/M3/M5/M6 lock exact 0-based positions and exact REF/ALT alleles derived by hand from the inputs and confirmed against the VCF/Hamming definitions; a deliberately-wrong implementation (off-by-one position, swapped allele, case-folding error) would fail them. ✓
- **No green-washing:** every detection assertion is exact (`Has.Count.EqualTo`, `Is.EqualTo` on positions/alleles/sequences); no ranges, no `Greater`, no widened tolerance, none skipped. ✓
- **Coverage:** both public methods exercised; branches covered = empty/null guard (S1,S2,S3,S4,S5), mismatch emit (M2–M6,C1), match-skip (M1,M3 idx3,M5 idx3), unequal-length prefix (M5), multi-SNP ordering (M3). Invariants INV-01..INV-06 all have a test. ✓
- **C1 note:** the property test computes the Hamming distance independently from the test inputs (`reference.Where((c,i) => c != query[i]).Count()`) rather than the implementation, so it is an honest cross-check, not a code echo; the differing positions {0,2,5,7} (count 4) are independently enumerable from the fetched definition. Acceptable.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6482`; SNP class `12/12`; `dotnet build` 0 errors (the 4 warnings are pre-existing NUnit2007 notices in unrelated files, not introduced here).

### Findings / defects (Stage B)
- **No code defect.** Both methods faithfully implement the validated description.
- **Minor coverage gap (not green-washed, intentionally not "fixed"):** there is no lowercase test for `FindSnpsDirect`. `FindSnpsDirect`/`FindSnps` perform **case-sensitive** char comparison, so `FindSnpsDirect("ATGC","atgc")` would report 4 SNPs. This is *not* a defect against the sources for the methods under test: VCF case-insensitivity governs VCF allele-string parsing, and the only case-normalizing method (`ClassifyMutation`) belongs to a different unit (VARIANT-CALL-001). Adding a test that asserts case-insensitive *detection* would assert behavior no retrieved source mandates; adding one that asserts case-*sensitive* detection would lock an undocumented convention. Left as a documented behavior note rather than a new test. The drives the PASS-WITH-NOTES.

## Verdict & follow-ups
- **Stage A: PASS.** **Stage B: PASS-WITH-NOTES** (case-handling of the detection path is undocumented at method level / case-sensitive by design; classification-level case-insensitivity is correctly out of unit scope).
- **End-state: CLEAN** — no defect found; description, implementation, and the 12 tests all agree with independently fetched authoritative sources; full suite green; no code/test change required.
- **Test-quality gate: PASS** — sourced exact expectations, no green-washing, all branches/invariants covered, honest full-suite green (6482/0).
