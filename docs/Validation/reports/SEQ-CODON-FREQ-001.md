# Validation Report: SEQ-CODON-FREQ-001 — Codon Frequencies

- **Validated:** 2026-06-16   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateCodonFrequencies(string dnaSequence, int readingFrame = 0)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Duplicate / sibling check

The prompt flagged sibling codon units (CODON-USAGE-001, CODON-RSCU-001). They are **distinct
algorithms, not duplicates** of this unit:

| Unit | Method | Metric |
|------|--------|--------|
| SEQ-CODON-FREQ-001 (this) | `SequenceStatistics.CalculateCodonFrequencies` | frequency = count(codon) / total counted codons (count/total over all codons) |
| CODON-USAGE-001 | `CodonOptimizer.CalculateCodonUsage` | raw integer counts per codon |
| CODON-RSCU-001 | `CodonUsageAnalyzer.CalculateRscu` / `GenomeAnnotator` | RSCU = n·x_j / Σ_k x_k (denominator is the *synonymous family*, not all codons) |

The denominator differs in each case (all codons vs none vs synonymous family), so this is a genuine
distinct metric. **Full validation performed.**

## Stage A — Description

### Sources opened this session (retrieved, not trusted from the repo)
- **EMBOSS `cusp` documentation** — https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html
  (WebFetch, 2026-06-16). Confirmed the five output columns and definitions:
  - *Frequency* = "the expected number of codons, given the input sequence(s), per 1000 bases".
  - *Fraction* = "the proportion of usage of the codon among its redundant set" (per-amino-acid; a
    **different** metric, correctly excluded from scope).
  - *Number* = observed count. Sum of the Number column = **386**.
  - Published rows: CGC Number=22 Frequency=56.995; GGC=23/59.585; GCC=18/46.632; TGA(*)=1/2.591.
- **Kazusa CUTG README** — https://www.kazusa.or.jp/codon/readme_codon.html (WebFetch, 2026-06-16).
  Confirmed frequency is reported **per thousand**, and "**Codons containing ambiguous base were
  excluded from count**" — exactly the non-ACGT exclusion rule used here.

### Formula check
frequency(x) = count(x) / total, over non-overlapping in-frame triplets, ambiguous codons excluded.
This is the count/total convention of CUTG; the per-thousand value equals this fraction × 1000.
Matches the cited sources exactly.

### Edge-case semantics
- Non-ACGT triplet excluded from count **and** total — sourced (Kazusa "ambiguous excluded").
- Trailing 1–2 bases ignored — direct consequence of non-overlapping triplets.
- Frame offset changes the codon multiset — non-overlapping triplets from frame start.
- Zero valid codons (total=0) → empty table. This is the single Assumption (Evidence §Assumptions);
  no source defines a non-empty result, and empty is the only value consistent with count/total.
  Accepted as PASS (conservative, non-fabricating, non-correctness-affecting for any real CDS).

### Independent cross-check (numbers I recomputed this session)
| Codon | Number | count/total × 1000 (recomputed) | cusp published | Match |
|-------|--------|----------------------------------|----------------|-------|
| CGC | 22 | 56.99481865… | 56.995 | ✓ |
| GGC | 23 | 59.58549222… | 59.585 | ✓ |
| GCC | 18 | 46.63212435… | 46.632 | ✓ |
| TGA | 1  | 2.59067357…  | 2.591  | ✓ |

Number column reconstructed in the M5 test sums to **386** (verified by independent Python sum of all
64 entries). Hand-computed small sequences also confirmed: `ATGATGAAA` f0 → ATG 2/3, AAA 1/3;
f1 → TGA 1.0; f2 → GAT 1/2, GAA 1/2; `ATGNNNAAA` f0 → ATG 1/2, AAA 1/2; `ATGAA` f0 → ATG 1.0.

### Findings / divergences
None. Description is biologically and mathematically correct and every non-trivial value traces to a
source retrieved this session.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs:688-723`.

### Formula realised correctly?
Yes. Guard: null/empty/length<3 → empty. Upper-cases input (INV-04). Loops `i` from `readingFrame` to
`length-3` in steps of 3 (non-overlapping). Each triplet of only `ATGC` increments count and total;
others (ambiguous, including RNA `U`) excluded (INV-03). Final pass divides each count by total;
total=0 leaves the map empty (no div-by-zero). Exactly the validated count/total definition.

### Cross-verification table recomputed vs code
The M5 test reconstructs the full 64-codon cusp dataset and asserts CGC=22/386, GGC=23/386 and their
×1000 values against the published 56.995 / 59.585. Full suite run confirms the code reproduces these.

### Variant/delegate consistency
Single public method; no `*Fast`/instance variants. MCP wrapper `AnalysisTools` delegates to it
unchanged. Sibling `CalculateCodonUsage` / `CalculateRscu` are separate metrics (see duplicate check).

### Test quality audit (HARD gate)
- **Sourced, exact values, not code echoes** — every expected value is an exact rational or the
  published cusp per-thousand figure (Within 1e-10 / 1e-3), traceable to Kazusa/cusp, not to code
  output. A deliberately-wrong implementation would fail these.
- **No green-washing** — no Greater/AtLeast/Contains where an exact value is known; no widened
  tolerances; no skips. `Is.All.InRange(double.Epsilon, 1.0)` for INV-01 correctly encodes (0,1].
- **Coverage** — all 5 invariants and all MUST(M1–M5)/SHOULD(S1–S4)/COULD(C1) cases covered.
  **Gap found & fixed this session:** the `readingFrame=2` branch was untested (only frames 0 and 1).
  Added `CalculateCodonFrequencies_Frame2_ReturnsExactFrequencies` (`ATGATGAAA` f2 → GAT 1/2, GAA 1/2,
  hand-computed from the count/total definition), exercising the remaining frame path.
- **Honest green** — FULL unfiltered suite: **Failed: 0, Passed: 6618**; `dotnet build` 0 errors
  (4 pre-existing NUnit2007 warnings in an unrelated `ApproximateMatcher` file; no new warnings).
- Consolidation per TestSpec §5 is in place: codon-frequency tests live only in the canonical file;
  legacy weak tests and the SEQ-DINUC-001 duplicates were already removed.

### Findings / defects
No implementation defect. One test-coverage gap (frame-2 branch) found and fixed.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS.** End-state: **CLEAN** — no implementation defect; the single
  test-coverage gap was fully fixed in-session and locked with a sourced exact-value test.
- Test-quality gate: **PASS** (after adding the frame-2 test).
