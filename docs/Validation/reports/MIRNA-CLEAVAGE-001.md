# Validation Report: MIRNA-CLEAVAGE-001 — Drosha/Dicer Cleavage-Site Prediction

- **Validated:** 2026-06-25   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.PredictDroshaDicerCleavage(string sequence, int basalJunction)` → `DroshaDicerCleavage?`
- **Stage A verdict:** ✅ PASS (2026-06-25) → 🟡 PASS-WITH-NOTES (B13 re-review 2026-11-25)
- **Stage B verdict:** 🟡 PASS-WITH-NOTES (2026-06-25) → ❌ FAIL → FIXED (B13 F27–F30, 2026-11-25; see the last section)
- **State:** ✅ CLEAN

Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs` (region "Drosha/Dicer cleavage-site prediction", ~L2021–2185).
Tests: `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs` (DD1–DD14).

---

## Stage A — Description

### Sources consulted (this session; "opened" = full text read, otherwise stated)

| # | Source | What it confirms |
|---|--------|------------------|
| 1 | **Han et al. (2006)**, *Cell* 125:887–901, "Molecular basis for the recognition of primary microRNAs by the Drosha–DGCR8 complex" (PubMed 16751099; cell.com S0092-8674(06)00516-2) | The Drosha cleavage site "is determined mainly by the distance (~11 bp) from the stem–ssRNA junction." DGCR8 anchors the basal ssRNA and directs Drosha cleavage ~11 bp (≈ one helical turn) up the stem (base-anchor model, pri-miR-16). |
| 2 | **Park et al. (2011)**, *Nature* 475:201–205, "Dicer recognizes the 5′ end of RNA for efficient and accurate processing" (PubMed 21753850; nature.com/articles/nature10198) — **abstract / search snippet only, full text not opened** | Human Dicer anchors the 5′ end, with "the cleavage site determined mainly by the distance (~22 nucleotides) from the 5′ end (5′ counting rule)." This fixes the mature length. |
| 3 | **Auyeung et al. (2013)**, *Cell* 152:844–858, "The Menu of Features that Define Primary MicroRNAs…" | The CNNC primary-sequence element acts ~17 nt (within a 16–18 nt window) downstream (3′) of the Drosha basal cut; SRSF3 binds it and recruits DROSHA to the basal junction. Optional licensing/confidence signal. |
| 4 | **Lee et al. (2003)** / Han 2006 | RNase III cleavage (Drosha, Dicer) leaves a **2-nt 3′ overhang**. |
| 5 | **miRBase** hsa-mir-21 hairpin **MI0000077** (mirbase.org/hairpin/MI0000077) | Precursor (72 nt), `ugucggg`**UAGCUUAUCAGACUGAUGUUGA**`cuguugaaucucaugg`**CAACACCAGUCGAUGGGCUGU**`cugaca`. 5p mature **hsa-miR-21-5p (MIMAT0000076)** = `UAGCUUAUCAGACUGAUGUUGA` at 1-based 8–29 (0-based 7–28, 22 nt). 3p **hsa-miR-21-3p (MIMAT0004494)** = `CAACACCAGUCGAUGGGCUGU` at 1-based 46–66 (0-based 45–65, 21 nt). |

### Formula / constant check (against sources)

| Constant (code) | Value | Source | Match |
|---|---|---|---|
| `DroshaCutBpFromBasalJunction` | 11 | Han 2006 (~11 bp from stem–ssRNA junction) | ✅ |
| `DicerCutNtFrom5PrimeEnd` | 22 | Park 2011 (~22 nt 5′ counting rule) | ✅ |
| `RNaseIII3PrimeOverhang` | 2 | Lee 2003 / Han 2006 (2-nt 3′ overhang) | ✅ |
| CNNC window | 16–18 nt | Auyeung 2013 (~17 nt, 16–18 nt window) | ✅ |
| `DroshaCut5' = basalJunction + 11` | — | Han 2006 ruler | ✅ |
| `mature = [DroshaCut5', DroshaCut5'+21]` (22 nt) | — | Park 2011 | ✅ |

### Definitions & conventions
0-based indices into the upper-cased (T→U) sequence; mature/star spans inclusive. `basalJunction` = 0-based index of the basal ssRNA–dsRNA junction (the first paired base of the 5′ arm, from which the +11 ruler is measured). Stated and standard.

### Independent cross-check — hand-trace on the REAL hsa-mir-21 precursor (MI0000077)

Precursor (0-based): `UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA` (72 nt).
- 5′ ssRNA flank `ugucggg` = 7 nt ⇒ basal ssRNA–dsRNA junction (start of dsRNA stem) at 0-based index **7** = the annotated **5p 5′ end**.
- A miRBase pre-miRNA is the **already-Drosha-trimmed Dicer substrate**: the ~11 bp of *lower stem below the Drosha cut has been removed*, so in the trimmed precursor the mature 5′ end coincides with the stem base (distance ≈ 0, not 11 bp). The Han +11 ruler therefore applies to the **intact pri-miRNA** lower stem.
- Reconstructing that lower stem (prepend an 11-nt stem, junction = 0): `DroshaCut5' = 0 + 11 = 11`; mature = `[11, 32]` (22 nt) = `UAGCUUAUCAGACUGAUGUUGA` = **hsa-miR-21-5p (MIMAT0000076) reproduced EXACTLY**.
- 2-nt 3′ overhang: linear star end = `matureEnd + 2 = 34`; star start = `34 − 22 + 1 = 13`.

### Findings / divergences (Stage A)
- All four measuring rules and constants are faithful to the primary sources. **PASS.**
- The published rules locate only the **5p mature** (Drosha 5′ cut + Dicer 22-nt length) and define the 2-nt overhang. They do **not** by themselves place the 3p (miRNA*) on the opposite arm without folding the hairpin — this is a property of the method's model, addressed in Stage B.

---

## Stage B — Implementation

### Code path reviewed
`PredictDroshaDicerCleavage` (`MiRnaAnalyzer.cs` ~L2110–2165) and helper `HasCnncMotifDownstream`.

### Formula realised correctly?
- `droshaCut5Prime = basalJunction + 11` ✅ (Han).
- `matureStart = droshaCut5Prime`, `matureEnd = matureStart + 22 − 1` ✅ (Park; 22-nt mature).
- `starEnd = matureEnd + 2`, `starStart = starEnd − 22 + 1` — the 2-nt overhang is applied in **linear** sequence coordinates (see note below).
- Null guards: null/empty, junction out of `[0, len)`, and each cut index past end of sequence → `null`. ✅
- T→U + upper-case normalisation applied before indexing; no alphabet validation (non-ACGU passed through). ✅
- CNNC: scans starts at cut+16..cut+18 for `C..C` 4-mer ✅ (Auyeung window).

### Cross-verification (recomputed vs code)

| Quantity | Hand value | Code (DD-tests) | Match |
|---|---|---|---|
| Drosha 5′ cut, junction 0 | 11 | 11 (DD1) | ✅ |
| Mature length | 22 | 22 (DD2) | ✅ |
| 3′ overhang; `DroshaCut3' − matureEnd` | 2 | 2 (DD3) | ✅ |
| 5p mature (synthetic pri, j=0) | `UAGCUUAUCAGACUGAUGUUGA` | = miRBase 5p (DD4) | ✅ |
| 5p mature (real MI0000077 + 11-nt lower stem) | `UAGCUUAUCAGACUGAUGUUGA` | = MIMAT0000076 (DD11) | ✅ |
| Linear star coords | start 13, end 34, len 22 | 13/34/22 (DD5) | ✅ |
| DNA input T→U | `UAGCUUAUCAGACUGAUGUUGA` | (DD6) | ✅ |
| Length boundary | 35 valid / 34 null | (DD14) | ✅ |

### Findings / defects

1. **3p (miRNA\*) span is a linear-geometry approximation — documented model boundary (not a hidden defect).**
   The method does **no hairpin folding**. It places the 3p span 2 nt 3′ of the 5p mature **end in linear sequence coordinates**, so the "star" is the 5p span shifted +2 (same 5′ arm), e.g. for hsa-mir-21 it predicts `GCUUAUCAGACUGAUGUUGACU`, whereas the **real** hsa-miR-21-3p (`CAACACCAGUCGAUGGGCUGU`, MIMAT0004494) lies on the **3′ arm across the terminal loop**. Recovering the true 3p arm requires folding/arm-pairing — out of scope here, alongside the documented out-of-scope trained miRDeep2 classifier. This is consistent with the algorithm-doc assumption ASM-03 ("represented asymmetrically").
   - **Fix applied this session (no behaviour change):** the XML doc-comment on `DroshaCut3Prime`/`StarStart` had overstated the star as "the base on the 3′ arm cut by Drosha." Corrected the record doc-comment, added a `<remarks>` scope/limitation paragraph, and corrected the inline body comment to state the linear-geometry approximation explicitly. Logged as FINDINGS_REGISTER **A-MIRNA-CLEAVAGE-001-1**.

2. **Test-quality defect: DD5 was a code echo.** The original `..._StarSpan_HasExpectedCoordinatesAndSequence` asserted `StarSequence == pri.Substring(13,22)` — the same coordinate arithmetic the code uses, adding no independent value.
   - **Fix applied:** replaced with `..._StarSpan_IsTwoNtThreePrimeOfMatureInLinearCoords` asserting the rule-derived coordinates and the order-robust invariant `StarStart − MatureStart == 2` (no substring echo).

### Test additions (lock sourced values & Stage-A edge cases)
- **DD11** — REAL miRBase precursor (MI0000077): 5p reproduces `UAGCUUAUCAGACUGAUGUUGA` EXACTLY **and** asserts the linear star ≠ real 3p (locks the documented limitation honestly).
- **DD12** — no-stem / non-hairpin (poly-A): still returns ruler coordinates (the method is a measuring rule, not a structure detector).
- **DD13** — non-ACGU / lower-case: `n` upper-cased and passed through verbatim (no alphabet validation).
- **DD14** — exact length boundary 35 valid / 34 null.
- Existing DD1–DD4, DD6–DD10 retained (Han/Park/Lee/Auyeung/miRBase-anchored; null/empty/out-of-range/too-short/CNNC).

### Variant/delegate consistency
Single public overload; no `*Fast`/instance variant. N/A.

---

## Verdict & follow-ups

- **Stage A: ✅ PASS** — measuring rules and constants confirmed against Han 2006, Auyeung 2013 (PDF), Lee 2003 and, for Park 2011, its abstract / search snippet only (the full text was not opened); 5p mature reproduces miRBase hsa-miR-21-5p exactly by hand-trace.
- **Stage B: 🟡 PASS-WITH-NOTES** — the 5p mature / Drosha–Dicer cut prediction is correct and independently verified; the 3p (miRNA\*) span is a deliberately simplified linear-geometry approximation (no folding) that does not reproduce the real 3′-arm 3p — a documented model boundary alongside the out-of-scope miRDeep2 classifier. Doc-comment overstatement and one code-echo test were fully fixed this session.
- **State: ✅ CLEAN** — no algorithm defect; the two issues found (doc overstatement, code-echo test) were fully fixed, edge-case coverage added (DD11–DD14), and the FULL unfiltered `dotnet test Seqeron.sln -c Debug` is green (Genomics 18775 passed, Failed: 0; 0 build warnings on changed files).
- Logged: FINDINGS_REGISTER A-MIRNA-CLEAVAGE-001-1.


## Runtime enforcement (LimitationPolicy)

This unit's guarded branch — the approximate linear 3p/star (miRNA\*) span (the 5p Drosha/Dicer cut is exact) — has **minimum access mode `Permissive`** (`Seqeron.Genomics.Core.LimitationCatalog`). Under the default `LimitationPolicy.DefaultMode = Moderate` it throws `SeqeronLimitationException` (this guarded branch is allowed only under `Permissive`); see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement. Additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.


## B13 re-review (2026-11-25) — FAIL → FIXED (F27–F30)

The 2026-06-25 verdict accepted a model that did not fold the hairpin. Re-checked against the sources:

- **F27 — 3p product on the wrong arm.** The "star" was the 5p span shifted +2 nt on the 5' arm. The
  published rule (2-nt 3' overhang across the duplex) needs the structure; miRDeep2 `find_star` and
  Auyeung et al. 2013 (cleavage sites inferred from RNAfold of hairpin ± 20 nt with a 2-nt 3' overhang)
  implement it. Now: 3p 3' end = partner(5p start) + 2, 3p start = partner(5p end − 2) with bulge
  compensation. The method also threw under the default `Moderate` LimitationPolicy mode.
- **F28 — the Drosha ruler counted nucleotides.** Han 2006 counts ~11 bp; the code added 11 nt to the
  junction. On real pri-miRNAs (junction 11 bp below the annotated 5p start, RNAfold structure) the old
  5' cut missed by −10 (mir-16-1), −5 (mir-30a), −4/−5 (let-7a-1), −3 (mir-21), −1 (mir-155) nt.
- **F29 — CNNC measured from the wrong site.** Auyeung numbers positive positions downstream of the 3p
  Drosha cut; the code scanned 16–18 nt 3' of the 5p cut (inside the 5p product) and missed the CNNC of
  mir-16-1 and mir-30a (the paper's own examples). Basal UG (−14/−13) flag added. Apical UGU/GUG has no
  positional rule in the text (documented, not flagged).
- **F30 — no stem required; 3' counting missing.** Poly-A returned "ruler coordinates"; now a junction that
  is not the 5' base of a ≥ 12-pair stem gives null, a malformed structure throws. Park 2011 3' counting
  added as `DicerCountingRule.ThreePrimeCounting`.

Oracle (miRBase + UCSC hg38 + ViennaRNA 2.8; predicted − annotated, 5p start/end | 3p start/end):
mir-21 0/−2 | +1/−1; let-7a-1 +1/0 | −1/0; mir-16-1 0/0 | 0/+1; mir-30a 0/0 | +2/0; mir-155 0/−2 | +2/0.
The old linear model put the 3p of every hairpin on the 5' arm (e.g. mir-21 `GCUUAUC…` vs miRBase 3p
`CAACACCAGUCGAUGGGCUGUCU`). Tests: DD1–DD14 rewritten; heavy-tier cleavage classes updated.

Runtime enforcement: `LimitationPolicy.Enforce("MIRNA-CLEAVAGE-001")` was removed from the method (no
approximate branch remains). The catalog entry in `LimitationPolicy.cs`, its `LIMITATIONS.md` rows and
`LimitationPolicy_Strict_Tests` section 5 are Core-owned and were filed as a cross-batch request.
