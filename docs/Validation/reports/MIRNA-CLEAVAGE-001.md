# Validation Report: MIRNA-CLEAVAGE-001 — Drosha/Dicer Cleavage-Site Prediction

- **Validated:** 2026-06-25   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.PredictDroshaDicerCleavage(string sequence, int basalJunction)` → `DroshaDicerCleavage?`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** 🟡 PASS-WITH-NOTES
- **State:** ✅ CLEAN

Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs` (region "Drosha/Dicer cleavage-site prediction", ~L2021–2185).
Tests: `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs` (DD1–DD14).

---

## Stage A — Description

### Sources opened (this session)

| # | Source | What it confirms |
|---|--------|------------------|
| 1 | **Han et al. (2006)**, *Cell* 125:887–901, "Molecular basis for the recognition of primary microRNAs by the Drosha–DGCR8 complex" (PubMed 16751099; cell.com S0092-8674(06)00516-2) | The Drosha cleavage site "is determined mainly by the distance (~11 bp) from the stem–ssRNA junction." DGCR8 anchors the basal ssRNA and directs Drosha cleavage ~11 bp (≈ one helical turn) up the stem (base-anchor model, pri-miR-16). |
| 2 | **Park et al. (2011)**, *Nature* 475:201–205, "Dicer recognizes the 5′ end of RNA for efficient and accurate processing" (PubMed 21753850; nature.com/articles/nature10198) | Human Dicer anchors the 5′ end, with "the cleavage site determined mainly by the distance (~22 nucleotides) from the 5′ end (5′ counting rule)." This fixes the mature length. |
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

- **Stage A: ✅ PASS** — measuring rules and constants confirmed verbatim against Han 2006, Park 2011, Auyeung 2013, Lee 2003; 5p mature reproduces miRBase hsa-miR-21-5p exactly by hand-trace.
- **Stage B: 🟡 PASS-WITH-NOTES** — the 5p mature / Drosha–Dicer cut prediction is correct and independently verified; the 3p (miRNA\*) span is a deliberately simplified linear-geometry approximation (no folding) that does not reproduce the real 3′-arm 3p — a documented model boundary alongside the out-of-scope miRDeep2 classifier. Doc-comment overstatement and one code-echo test were fully fixed this session.
- **State: ✅ CLEAN** — no algorithm defect; the two issues found (doc overstatement, code-echo test) were fully fixed, edge-case coverage added (DD11–DD14), and the FULL unfiltered `dotnet test Seqeron.sln -c Debug` is green (Genomics 18775 passed, Failed: 0; 0 build warnings on changed files).
- Logged: FINDINGS_REGISTER A-MIRNA-CLEAVAGE-001-1.

## Runtime enforcement (LimitationPolicy)

Under the default `LimitationPolicy.DefaultMode = Strict`, the approximate linear 3p/star (miRNA\*) span (the 5p Drosha/Dicer cut is exact) throws `Seqeron.Genomics.Core.SeqeronLimitationException` (named limitation + workaround; see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement and `LimitationCatalog`). `Permissive` mode returns the historical best-effort value. This is an additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
