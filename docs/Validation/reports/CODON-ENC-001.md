# Validation Report: CODON-ENC-001 — Effective Number of Codons (ENC / Nc), Wright 1990

- **Validated:** 2026-06-15   **Area:** Codon
- **Canonical method(s):** `CodonUsageAnalyzer.CalculateEnc(string)` (core), `CalculateEnc(DnaSequence)` (delegate); private `CalculateEncCore`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened & what they confirm
- **Fuglsang A. (2004) "The 'effective number of codons' revisited", BBRC 317:957–964** — retrieved the PDF and extracted the full text with `pdftotext -layout`. Confirms verbatim:
  - **Eq. (1)** codon homozygosity `F̂ = (n·Σpᵢ² − 1)/(n − 1)`, "n is the total count for the amino acid in the gene, and pᵢ is the codon frequency for the ith synonymous codon".
  - **Eq. (3)** `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`, "F̂₂ is the average homozygosity for the amino acids having a degeneracy of two (histidine, glutamine, etc.) and so on".
  - **Eq. (4)** missing **member** of a class: "If … a gene does not contain threonine, then F̂₄ will be the average of the codon homozygosities of glycine, valine, alanine, and proline" — i.e. average only the *estimable* members of the same class.
  - **Eq. (5a)** isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2`.
  - **Cap rule:** "There is a chance that N̂c … will exceed 61. In that case, Wright recommends re-adjusting the result down to 61."
  - **Range:** "a number between 20 and 61 … approaches 20 codons for the extremely biased genes, and approaches 61 for the genes where all possible codons are used with no preference."
- **Peden J. (codonW thesis), §"Equation 2-7/2-8", Nc** — retrieved PDF, extracted with `pdftotext`. The de-facto **reference implementation**. Confirms the formula and, decisively, the **whole-class-empty rule**:
  > "If amino-acids are rare or missing, adjustments must be made. For absent amino-acids the numerator is decreased, to reciprocally average for only the amino-acids present. **For sequences where a synonymous family of amino-acids is empty (F̂ₙ = 0), Nc is not calculated**, as the gene is assumed to be either too short or to have extremely skewed amino-acid usage (Wright 1990). An exception is made … isoleucine 3-fold … absent … F̂₃ estimated as the average of F̂₂ and F̂₄." Cap to 61 confirmed.
- **Standard genetic-code degeneracy partition** (NCBI table 1, cross-checked against the in-code `CodonToAminoAcid` map): 9 doublets, 1 triplet (Ile = ATT/ATC/ATA), 5 quartets, 3 sextets (Leu, Ser, Arg), 2 singlets (Met, Trp). Matches Eq. (3) constants `2, 9, 1, 5, 3`.

### Formula check
All four formulae (Eq. 1, 3, 4, 5a), the cap at 61, and the degeneracy partition are correctly stated in the TestSpec/Evidence and match the cited primary source word-for-word.

### Edge-case semantics check
- `n ≤ 1` for an amino acid → F̂ undefined (denominator n−1) → that aa excluded: **sourced** (Fuglsang).
- Missing **member** of a populated class → class mean of estimable members (Eq. 4): **sourced**.
- Empty **3-fold (Ile)** class → Eq. 5a: **sourced**.
- Empty **2/4/6-fold class** → reference (codonW) says **"Nc is not calculated"**. The TestSpec/Evidence (as originally authored) instead prescribed "absent class contributes its full codon count" and even listed *that very thing* as failure-mode #2 — a self-contradiction. **Divergence (note A1).**
- Lower clamp at 20 → **not** a Wright/codonW instruction (codonW caps only the top at 61). Recorded by the spec as a defensive assumption. **Note A2.**

### Independent cross-check (numbers)
Independent Python re-implementation of the Wright/codonW algorithm (Eq. 1 + class averaging + Eq. 5a + cap), run this session:
- Fully-populated biased gene (M3): **Nc = 41.288461538461526**.
- Same gene, Ile removed, Eq. 5a active (M5): **Nc = 39.47394540942927**.
- 2:1 bias across all families (C1): F̂ = 1/3 each → **Nc = 56.0**.
- One codon per aa (M1): F̂ = 1 each → **Nc = 20**.
The C# implementation reproduces all four to full double precision (see Stage B table).

### Findings / divergences
- **Note A1 (whole-class-absent):** description's prescribed behaviour for an entirely-empty 2/4/6-fold class diverges from the reference (codonW "Nc not calculated"). Minor, because real coding sequences never trigger it. Spec updated this session to record the divergence honestly (§7).
- **Note A2:** lower clamp at 20 is a defensive bound, not a sourced rule (unchanged conclusion; documented).

Stage A is **PASS-WITH-NOTES** — the core biology/maths is correct and sourced; two documented edge-case divergences remain, neither affecting real input.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs:274–360`.

### Formula realised correctly? (evidence)
- Eq. (1) at `:312–318` — `F = (n*Σp² − 1)/(n − 1)`, `p = nᵢ/n`. Exact.
- `n ≤ 1` skipped (`:309`); singlets skipped (`:306`). Correct.
- Eq. (4) within-class averaging via `AverageOrNull` (`:330–355`) — averages only estimable members. Correct.
- Eq. (5a) Ile fallback (`:337–338`) — fires only when F₃ unestimable *and* F₂,F₄ present. Correct.
- Eq. (3) aggregation (`:343–347`) with constant `2`. Cap to 61 (`:349`). Correct.
- **`ClassContribution` (`:357–360`)**: when a whole class has no estimable F̂, returns the raw codon count (full count) instead of declining to compute. **Diverges from codonW (defect B1).** Reachable only on genes missing an entire 2/4/6-fold class.
- Lower clamp `Math.Max(20, …)` (`:349`): library convention, not sourced (note A2).

### Cross-verification table recomputed vs code

| Case | Input | Independent reference Nc | C# `CalculateEnc` | Match |
|------|-------|--------------------------|-------------------|-------|
| M1 one-codon-per-aa | each aa ×2, single codon | 20.0 | 20.0 | ✅ |
| M2 near-uniform | all 61 sense codons ×2 | >61 → cap 61.0 | 61.0 | ✅ |
| M3 fully-populated biased | all classes estimable | 41.288461538461526 | 41.288461538461526 | ✅ |
| M5 Ile-absent (Eq. 5a) | all classes but Ile | 39.47394540942927 | 39.47394540942927 | ✅ |
| C1 uniform 2:1 bias | F̂=1/3 every class | 56.0 | 56.0 | ✅ |
| M5b whole-class-absent | only Phe | (undefined per codonW) | 29.0 (library convention) | n/a — divergence pinned |

The core algorithm is correct to double precision on every fully-populated input.

### Variant/delegate consistency
`CalculateEnc(DnaSequence)` null-checks then delegates to the string core (`:274–277`); equality verified by test S3. `GetStatistics` reuses `CalculateEncCore` (`:393`). Consistent.

### Test quality audit (HARD gate)
- **Code-echo defects found & fixed:** the original M3 (expected 29.0) and M5 (expected 40.4) asserted values produced by the *unsourced* full-count fallback (defect B1) — they would pass against the very behaviour the spec itself names as a defect. **Rewritten** to sourced exact values on fully-populated genes (M3 = 41.288461538461526; M5 = 39.47394540942927), each traced to the independent reference implementation run this session, and each genuinely exercising Eq. 1 / Eq. 3 / Eq. 4 (and, for M5, Eq. 5a with all of F₂/F₄/F₆ estimable).
- **Divergence pinned honestly:** added test **M5b** that pins the whole-class-absent value (29.0) but is explicitly labelled LIBRARY-SPECIFIC / NOT-Wright in name, comment, and assertion message, so it can never be mistaken for a sourced result. This keeps the `ClassContribution` fallback path under coverage without green-washing.
- **No weakening:** no tolerances widened, no assertions softened, no tests skipped. M1/M2/C1 (already sourced-correct) retained. M4 range-invariant (property test) retained; M6 (null→throws), M7 (empty→0), S1 (case), S2 (invalid codons), S3 (delegate) retained.
- **Coverage:** both public overloads, all formula paths (Eq. 1/3/4/5a), cap-at-61, range invariant, null/empty/invalid/case edge cases, and the divergent absent-class path are all exercised.
- **Honest green:** full unfiltered `dotnet test` = **6527 passed, 0 failed**, 1 unrelated benchmark skipped; `dotnet build` 0 errors (4 pre-existing unrelated warnings). Gate **PASS**.

### Findings / defects
- **B1 (logged):** whole-class-absent fallback adds the full codon count instead of declining to compute (codonW: "Nc not calculated"). Severity low — unreachable on real coding sequences; core formula proven correct. A complete fix (return a sentinel/NaN for "Nc undefined") would change the `double`-returning contract and ripple into `GetStatistics` and the MCP `EncResult`, which is a product/contract decision out of scope for one validation session. Behaviour is now documented (spec §7, FINDINGS_REGISTER) and pinned by an explicitly-labelled test rather than masquerading as sourced.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (formulae correct & sourced; two documented edge-case divergences).
- **Stage B:** PASS-WITH-NOTES (core formula exact vs independent reference; one low-severity documented divergence on degenerate input; code-echo tests fixed to sourced values).
- **End-state:** ✅ CLEAN — every defect found was either fixed (test code-echoes → sourced exact values; added honest divergence test) or fully documented with a clear rationale; no half-fix. Algorithm is fully functional for all real coding sequences.
- **Follow-up (optional, not blocking):** FR — decide a contract for "Nc undefined" (whole class empty) to match codonW semantics across `CalculateEnc`/`GetStatistics`/MCP.
