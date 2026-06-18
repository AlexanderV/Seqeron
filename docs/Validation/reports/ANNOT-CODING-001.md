# Validation Report: ANNOT-CODING-001 — Coding Potential Calculation (CPAT hexamer usage-bias score)

- **Validated:** 2026-06-15   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.CalculateCodingPotential(string, IReadOnlyDictionary<string,double>, IReadOnlyDictionary<string,double>, int, int)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS (after in-session fix)
- **End-state:** CLEAN (defect found and completely fixed)

## Stage A — Description

### Sources opened & what they confirm

| # | Source | Retrieved | Confirms |
|---|--------|-----------|----------|
| 1 | CPAT paper, Wang et al. (2013) NAR 41(6):e74, PMC3616698 | (Evidence file, 2026-06-13) | Hexamer score = log-likelihood ratio of in-frame hexamer frequencies coding vs noncoding; positive ⇒ coding, negative ⇒ noncoding |
| 2 | **Canonical CPAT** `liguowang/cpat` `src/cpmodule/FrameKmer.py` `kmer_ratio`/`word_generator` | **this session, 2026-06-15** | Exact branch logic incl. the both-zero `continue`; `math.log` = ln; word 6 / step 3, frame 0; full-length words only; `len(seq)<word_size→0`; empty-count→ except −1 |
| 3 | lncScore copy `WGLab/lncScore` `tools/cpmodule/FrameKmer.py` | **this session, 2026-06-15** | Identical `kmer_ratio` branch logic (same both-zero `continue`), confirming source 2 is not repo-specific |

### Formula check

Per-hexamer contribution (canonical `kmer_ratio`, frame-0 loop, lines 88-101 of `liguowang/cpat`):

- `coding>0 && noncoding>0`  → `+= ln(coding/noncoding)` (natural log — `math.log`); count++
- `coding>0 && noncoding==0` → `+= 1`; count++
- `coding==0 && noncoding==0` → **`continue` (NOT counted)**
- `coding==0 && noncoding>0` → `-= 1`; count++
- `else` (e.g. negatives) → `continue`
- hexamer missing from either table → `continue` (not counted)
- `len(seq) < word_size` → return 0
- return `sum/count`; if `count==0` the `ZeroDivisionError` is caught → reference returns −1

Score = (Σ contributions) / (count of scored hexamers). Matches the description's §2.2 **except** the
both-zero branch, which the original description/Evidence omitted (see findings).

### Edge-case semantics check
All sourced: too-short→0; missing-key→skip; coding-only→+1; noncoding-only→−1; both-zero→skip-not-counted;
no-scorable-hexamer→reference −1 (port returns 0, documented ASSUMPTION-1).

### Independent cross-check (numbers, hand-computed from canonical `kmer_ratio`)

| Input | Hexamers scored | Computation | Expected |
|-------|-----------------|-------------|----------|
| `ATGAAACCC`, C{ATGAAA:8,AAACCC:2} N{ATGAAA:2,AAACCC:4} | ATGAAA, AAACCC | (ln(8/2)+ln(2/4))/2 = (1.3862943611198906−0.6931471805599453)/2 | **0.34657359027997264** |
| `ATGAAA`, C{ATGAAA:4} N{ATGAAA:1} | ATGAAA | ln(4)/1 | **1.3862943611198906** |
| coding-only ATGAAA C=5 N=0 | ATGAAA | +1/1 | **1.0** |
| noncoding-only ATGAAA C=0 N=5 | ATGAAA | −1/1 | **−1.0** |
| both-zero ATGAAA(0,0)+AAACCC(4,1) on `ATGAAACCC` | **AAACCC only** (ATGAAA skipped, NOT counted) | ln(4)/1 | **1.3862943611198906** |

### Findings / divergences (Stage A)
- **DEFECT (description):** §2.2 of the algorithm doc, the Evidence `kmer_ratio` transcription, and TestSpec C1
  all omitted the canonical `elif coding[k]==0 and noncoding[k]==0: continue` branch and described the
  both-zero hexamer as "contributes 0 but is counted." The canonical CPAT and lncScore reference both
  `continue` (NOT counted). Description corrected this session in all three docs.
- ASSUMPTION-1 (no-scorable-hexamer returns 0 vs reference −1) is a defensible documented port choice; both are non-scores. Retained.

Stage A verdict: **PASS-WITH-NOTES** — biology/maths correct; one transcription error in the both-zero branch, fixed.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs:571-618`.

### Defect found
At line 610 the original code had **no** both-zero branch — it fell through and unconditionally executed
`scoredHexamers++`, so a both-zero in-both hexamer was **counted as a scored 0**, contradicting the
canonical reference (which `continue`s and does not count it). This inflates the denominator and dilutes
the mean whenever a both-zero hexamer occurs.

**Fix applied:** added `else if (coding == 0 && noncoding == 0) continue;` and a trailing `else continue;`
before the `scoredHexamers++`, verbatim to the canonical `kmer_ratio`.

### Formula realised correctly? (after fix)
Yes. In-frame extraction (offset 0, step `stepSize`, full-length `wordSize` words), `Math.Log` (= ln),
the four contribution branches, skip-if-missing, mean over scored count, `len<wordSize→0`,
`count==0→0` (ASSUMPTION-1). Matches canonical reference.

### Cross-verification table recomputed vs code
All five rows above reproduced by the test suite (6561 passed). C1 now asserts the sourced
**1.3862943611198906** (was the code-echoing 0.6931471805599453).

### Variant/delegate consistency
Single canonical method; MCP `AnnotationTools.cs:144` calls through with no behavioural dependency on the fixed branch.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** C1 previously asserted the implementation's wrong behaviour
  (both-zero counted). Rewritten from the canonical sourced value (both-zero skipped → score = ln 4).
  All exact values (M1–M6, M8, C1) trace to the canonical `kmer_ratio` + hand derivation, not to code output.
- **No green-washing:** no assertion weakened; exact `Is.EqualTo(...).Within(1e-10)` everywhere; no skips/ignores; expected values corrected to the source, never to actual output.
- **Coverage:** all Stage-A branches exercised — both-positive (M1,M2,M5,M6,M8), coding-only +1 (M3),
  noncoding-only −1 (M4), both-zero skip (C1), missing-key skip (M8), in-frame stepping (M9),
  too-short/empty→0 (M7,S2), no-scorable→0 (M15), case-insensitivity (S1), all four validation throws
  (M10–M14). The `else: continue` (negative values) branch is out of contract (values constrained non-negative) and is left untested by design.
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6561**; `dotnet build` 0 errors,
  no new warnings in the changed files (4 pre-existing warnings unrelated to this unit).

### Findings / defects
- FIND: both-zero hexamer counted instead of skipped (code + test C1 + 3 docs). FIXED this session.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS (post-fix). **End-state: CLEAN.**
- Defect logged in FINDINGS_REGISTER; code, test, Evidence, TestSpec and algorithm doc all corrected.
- Test-quality gate: **PASS** (no green-washing; sourced exact values; branches covered; honest full-suite green).
