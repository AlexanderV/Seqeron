# Validation Report: CODON-USAGE-001 — Codon Usage Analysis

- **Validated:** 2026-06-12   **Area:** Codon
- **Canonical method(s):** `CodonOptimizer.CalculateCodonUsage(string)` (codon counts), `CodonOptimizer.CompareCodonUsage(string, string)` (Total-Variation-Distance similarity)
- **Adjacent method validated for the task's RSCU emphasis:** `CodonUsageAnalyzer.CalculateRscu(...)` (RSCU)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

---

## Scope note (important)

The validation prompt for this session is written around **RSCU** (Relative Synonymous
Codon Usage). However, the **authoritative unit definition** for CODON-USAGE-001 — per
`tests/TestSpecs/CODON-USAGE-001.md`, `docs/Evidence/CODON-USAGE-001-Evidence.md`,
`ALGORITHMS_CHECKLIST_V2.md` (line 1291), and the named test file
`CodonOptimizer_CodonUsage_Tests.cs` — is the pair:

| Method | Class | Role |
|--------|-------|------|
| `CalculateCodonUsage(seq)` | `CodonOptimizer` | Canonical (codon counts) |
| `CompareCodonUsage(seq1, seq2)` | `CodonOptimizer` | Comparison (TVD similarity) |

RSCU is **not** part of CODON-USAGE-001. It is implemented in a *different* type,
`CodonUsageAnalyzer.CalculateRscu`, and the spec explicitly records this separation
(Decision **D1**: "Tests focus on `CodonOptimizer` methods, not `CodonUsageAnalyzer`
(separate test unit)"). To honour both the prompt and the real unit, this report
validates the CODON-USAGE-001 methods in full **and** independently verifies the RSCU
formula in `CodonUsageAnalyzer`. Both are correct.

---

## Stage A — Description

### Sources opened & what they confirm

1. **Wikipedia — "Codon usage bias"** (fetched). Confirms degeneracy of the genetic
   code (64 codons → 20 aa + 3 stop) and that codon-usage bias = differing frequencies
   of synonymous codons. The page lists CAI, Fop, RCA, Nc as bias measures; it does not
   itself give the RSCU equation (noted), so RSCU was cross-checked against primary/
   reference sources below.
2. **Kazusa Codon Usage Database** (per Evidence doc). Standard table format; codon
   frequencies reported per 1000 codons; trailing incomplete codons ignored. Used as the
   convention source for counting and the per-1000 frequency unit.
3. **Sharp & Li (1987), CAI paper** + **seqinr `uco` reference (R/CRAN)** (fetched) and
   the RSCU literature search. These give the canonical RSCU definition (below).

### Formula check

**Codon counts / frequency.** Count(c) = number of in-frame triplets equal to codon c;
frequency = count / total codons (optionally ×1000 for the Kazusa per-mille unit). Matches
sources.

**Comparison metric — Total Variation Distance similarity.** The unit uses
`Similarity = 1 − (Σ_c |f₁(c) − f₂(c)|) / 2`, where f_i(c) is the codon frequency in
sequence i and the sum runs over the union of observed codons. This is the standard TVD
between two discrete distributions (TVD = ½·L¹), giving a similarity in [0, 1]. Properties
(identity=1, symmetry, range [0,1], disjoint→0) are genuine mathematical facts of TVD.

**RSCU (canonical, Sharp & Li).**
RSCU_i = X_i / ( (1/n) Σ_j X_j ) = n · X_i / Σ_j X_j, where X_i is the observed count of
codon i, the sum Σ_j runs over the n synonymous codons of that amino-acid family, and n is
the family size. RSCU = 1 ⇒ used exactly as expected under uniform synonymous usage;
> 1 over-represented; < 1 under-represented. Confirmed verbatim in sense by seqinr `uco`
("number of times a codon is observed relative to uniform synonymous usage; absence of
bias ⇒ RSCU = 1.00") and the RSCU literature.

### Edge-case semantics

- Empty / null sequence → empty result; no codons. (Convention; sourced as "no data".)
- Incomplete trailing 1–2 nt → ignored (Kazusa / EMBOSS convention).
- DNA T treated as RNA U (biological equivalence).
- TVD: empty / one-empty → 0 similarity (no comparable data).
- RSCU: single-codon families (Met = AUG, Trp = UGG) ⇒ n=1, X_i = family total ⇒
  RSCU = 1. Absent amino acid: seqinr reports NA; Seqeron returns 0 — a defensible
  divergence (no observation ⇒ no meaningful ratio), not a scientific error.

### Independent cross-check (hand computation)

Amino acid with n = 4 synonymous codons, counts [10, 20, 30, 40], family total = 100.
Expected per codon = 100 / 4 = 25.
RSCU = [10/25, 20/25, 30/25, 40/25] = **[0.4, 0.8, 1.2, 1.6]**.
Sum = 0.4 + 0.8 + 1.2 + 1.6 = **4.0 = n**. The **sum-to-n** property holds
(Σ_i RSCU_i = Σ_i n·X_i/ΣX = n·(ΣX/ΣX) = n).

TVD worked example (spec M7): seq1 f(CUG)=3/4, f(CUA)=1/4; seq2 f(CUA)=3/4, f(CUG)=1/4 ⇒
Σ|Δ| = ½+½ = 1 ⇒ sim = 1 − ½ = **0.5**. Confirmed.

### Findings / divergences

PASS-WITH-NOTES: only because the session prompt's RSCU framing does not match the actual
CODON-USAGE-001 unit (counts + TVD). The unit's own description is mathematically correct;
RSCU (in the adjacent type) is also correct. No formula is wrong.

---

## Stage B — Implementation

### Code path reviewed

- `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`
  - `CalculateCodonUsage` (lines 634–652): upper-cases, T→U, splits into in-frame
    triplets via `SplitIntoCodons` (line 687), counts into a dictionary.
  - `CompareCodonUsage` (lines 657–681): computes per-codon frequencies over each
    sequence's own total and returns `1 − Σ|f₁−f₂|/2`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs`
  - `CalculateRscuCore` (lines 85–107): groups codons by amino acid, computes
    `expected = totalCount / numSynonymous`, then `rscu = observed / expected`
    (= n·X_i / ΣX_j) — exactly the canonical formula. Returns 0 when expected = 0.

### Formula realised correctly?

Yes. `CalculateCodonUsage` realises Count(c). `CompareCodonUsage` realises the TVD
similarity exactly (frequency = count/total per sequence; ½·L¹; 1 − TVD). `CalculateRscu`
realises RSCU with the **n factor** (`numSynonymous`) and **per-amino-acid-family**
normalisation (`totalCount` summed over the synonymous group only).

`SplitIntoCodons` loop `for (i=0; i+2 < len; i+=3)` is equivalent to `i ≤ len−3`, i.e.
only complete in-frame codons — trailing 1–2 nt correctly dropped (verified for len 3, 5, 9).

### Cross-verification table recomputed vs code

| Case | Input | Expected | Source |
|------|-------|----------|--------|
| Counts | "AUGGCUGCU" | {AUG:1, GCU:2} | M1 ✓ |
| Incomplete | "AUGGC" | {AUG:1} | M3 ✓ |
| All 64 | every codon ×1 | 64 keys, each 1 | M5 ✓ |
| TVD 0.5 | CUG/CUA mirror | 0.5 | M7 ✓ |
| TVD 0.75 | M9 3-codon | 0.75 | M9 ✓ |
| TVD 0.0 / 0.25 / 0.75 | M10 | 0.0 / 0.25 / 0.75 | M10 ✓ |
| TVD 2/3 | shared AUG | 0.666… | S6 ✓ |
| RSCU n=4 [10,20,30,40] | hand | [0.4,0.8,1.2,1.6], Σ=4 | this report ✓ |

All match (tests pass — see below).

### Variant/delegate consistency

`CompareCodonUsage` delegates to `CalculateCodonUsage`. `CalculateRscu(DnaSequence)` and
`CalculateRscu(string)` both delegate to `CalculateRscuCore`. Consistent.

### Test quality audit

`CodonOptimizer_CodonUsage_Tests.cs` (this unit) asserts **exact** values
(`Is.EqualTo(...).Within(1e-10)`), covering all M/S/edge cases, sum-to-total invariant,
symmetry, range, disjoint→0, empty/null. RSCU sum-to-n / unbiased / biased / Met-Trp /
empty edge cases are covered by `CodonUsageAnalyzerTests.cs` (the separate RSCU unit).
Tests are real (exact numbers, deterministic).

### Findings / defects

None. Code faithfully realises the validated formulas.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — unit description (counts + TVD) and the adjacent RSCU
  formula are both mathematically correct; the only "note" is that the session prompt's
  RSCU framing belongs to a *different* test unit (`CodonUsageAnalyzer`), per spec D1.
- **Stage B:** PASS — `CalculateCodonUsage`, `CompareCodonUsage`, and (cross-checked)
  `CalculateRscu` all match authoritative sources; RSCU carries the **n factor** with
  **per-family** normalisation and satisfies the **sum-to-n** property.
- **State:** CLEAN — no defect; no code change required. Build green; CODON-USAGE-001
  filter: 60 tests pass (includes the adjacent `CodonUsageAnalyzerTests`); full
  `Seqeron.Genomics.Tests` suite: 4484 passed, 0 failed (baseline preserved).

### Sources
- Wikipedia, "Codon usage bias": https://en.wikipedia.org/wiki/Codon_usage_bias
- seqinr `uco` (R/CRAN reference implementation): https://search.r-project.org/CRAN/refmans/seqinr/html/uco.html
- Sharp PM, Li WH (1987), Nucleic Acids Research 15(3):1281–1295 (CAI / RSCU).
- Kazusa Codon Usage Database: https://www.kazusa.or.jp/codon/
