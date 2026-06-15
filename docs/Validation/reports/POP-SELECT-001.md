# Validation Report: POP-SELECT-001 — Selection Signature Detection (integrated Haplotype Score, iHS)

- **Validated:** 2026-06-15   **Area:** PopGen
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateEhh`, `.CalculateIHS`, `.StandardizeIHS`, `.ScanForSelection` (+ internal `IntegrateEhh` / `IntegrateDirection` / `Choose2`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened this session
1. **Voight, Kudaravalli, Wen & Pritchard (2006), PLoS Biology 4(3):e72** — fetched the PLoS article (Materials & Methods) live this session. Confirmed verbatim:
   - Unstandardized iHS: "iHS = ln(iHH_A/iHH_D)".
   - Integration: "The EHH values at successive SNPs are joined by straight lines, and then we compute the total area under each curve, between the nearest points to the left and right of the core SNP where the EHH drops below 0.05." (= trapezoidal rule, both directions, 0.05 truncation).
   - Standardization: "The expectation and standard deviation of ln(iHH_A/iHH_D) are estimated from the empirical distribution at SNPs whose derived allele frequency p matches the frequency at the core SNP."
   - Scan: "the proportion of SNPs with |iHS| > 2"; "standard window size of 50 SNPs".
2. **Szpiech & Hernandez (2014), selscan, MBE 31(10):2824** — downloaded the arXiv PDF (1403.6854) and extracted text with `pdftotext` this session. Confirmed verbatim:
   - Eq. (3): EHH_c(x_i) = Σ_{h∈H_c} C(n_h,2)/C(n_c,2), with C(x_1)={11,10,00,01}, H_1={11,10}, H_0={00,01}.
   - Eq. (4): iHH_c = trapezoidal quadrature ½Σ(EHH_c(x_{i-1})+EHH_c(x_i))·g(x_{i-1},x_i) over downstream + upstream markers.
   - Eq. (5): unstandardized iHS = ln(iHH_1/iHH_0), with the explicit note: "this definition differs slightly from that in Voight et al. (2006), where unstandardized iHS is defined with iHH_1 and iHH_0 swapped" — confirming Voight = ln(iHH_A/iHH_D), the reciprocal/opposite sign of selscan.
   - Truncation: "the sums … are truncated at x_i — the marker at which the EHH … < 0.05" (Eq. 6/7 region; also "once EHH_c(x_i) < 0.05" for per-allele iHS).

### Formula check
- EHH_c = Σ C(n_h,2)/C(n_c,2) — matches selscan Eq. 3 exactly; the rehh form (1/(n_a(n_a−1)))·Σ n_k(n_k−1) is algebraically identical (factor of 2 cancels). ✅
- iHH = trapezoidal area both directions, truncated at first EHH<0.05 — matches Voight M&M and selscan Eq. 4 + truncation note. ✅
- Unstandardized iHS = ln(iHH_A/iHH_D) (Voight) — sign convention correctly chosen over selscan's reciprocal; divergence is documented (it is a pure sign flip). 🟡 (documented divergence, resolved in favour of the primary peer-reviewed source)
- Standardization (x − E_p)/SD_p in derived-freq bins — matches Voight. The N−1 (sample) SD and the 0.05 bin width are explicit assumptions Voight does not pin down; bin width matches rehh `freqbin`. ✅ (assumptions documented, magnitude-only impact)
- Scan: proportion of |iHS|>2 per window, threshold 2.0, default window 50 — matches Voight M&M. ✅

### Edge-case semantics
EHH=1 for a single chromosome (trivially homozygous), EHH=0 for all-distinct, EHH=0 empty sample — boundaries of C(n,2); monomorphic core / non-{0,1} core / null / length-mismatch / coreIndex OOR all have defined sourced throws (Voight: iHS only for polymorphic SNPs with ancestral state). Balanced decay ⇒ ln(1)=0. All defined and sourced.

### Independent cross-checks (numbers retrieved this session)
- M3 constructed panel hand-trace: derived (3× identical) ⇒ EHH=1 at each flank ⇒ iHH_D = ½(1+1)·10 ×4 trapezoids = **40.0**; ancestral (3× distinct) ⇒ EHH=0 at first flank ⇒ iHH_A = ½(1+0)·10 ×2 = **10.0**; ln(10/40) = **−1.3862943611198906** (verified with Python `math.log(0.25)`).
- S1 rehh ratio (SNP F1205400, IHH_A=284429.9, IHH_D=2057107.4): `math.log(284429.9/2057107.4)` = **−1.9785692742315621** (verified with Python this session).
- EHH worked values (selscan Eq. 3): {11,11,11,10} ⇒ (C(3,2)+C(1,2))/C(4,2) = 3/6 = **0.5**; {00,00,01,01} ⇒ (1+1)/6 = **0.333…**.

### Findings / divergences
No biological or mathematical error. One documented sign divergence vs selscan (resolved to Voight). Two documented standardization assumptions (N−1 SD, 0.05 bin width) affecting only standardized magnitude, not the canonical unstandardized score. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`:
- `CalculateEhh` (972–992): hashes whole-window strings, sums Choose2(n_h)/Choose2(n_c); nc=0→0, nc=1→1. Matches selscan Eq. 3.
- `CalculateIHS` (1015–1063): validates inputs, partitions on core allele {0,1}, integrates both alleles, returns ln(iHH_A/iHH_D) with Voight sign (1058–1060).
- `IntegrateDirection` (1182–1216): starts EHH=1 at the core, walks outward, adds trapezoid ½(prev+ehh)·|Δpos| **then** breaks once ehh<0.05 — exactly "area up to the nearest point where EHH drops below 0.05" (Voight) / "truncated at x_i" (selscan). Both directions summed in `IntegrateEhh`.
- `StandardizeIHS` (1077–1123): bins by `(int)(p*binCount)` (p=1 folded into top bin), per-bin (x−mean)/sampleSD, sd=0⇒0.
- `ScanForSelection` (1135–1162): non-overlapping windows, ExtremeCount = #(|score|>2.0), proportion = extreme/snpCount.
- Constants: `EhhIntegrationCutoff=0.05`, `ExtremeIhsThreshold=2.0`, default window 50, default bins 20. All match sources.

### Formula realised correctly?
Yes. Hand-trace of the M3 panel through the code reproduces iHH_A=10, iHH_D=40, iHS=ln(0.25); the S1 asymmetric panel reproduces UnstandardizedIHS == ln(IhhAncestral/IhhDerived). EHH worked values reproduce 0.5 and 1/3. Sign convention verified by negative-mutation test (below).

### Cross-verification table recomputed vs code
| Case | Source value | Code value | Match |
|------|-------------|-----------|-------|
| EHH {11,11,11,10} | 0.5 | 0.5 | ✅ |
| EHH {00,00,01,01} | 0.333… | 0.333… | ✅ |
| EHH single / distinct / empty | 1 / 0 / 0 | 1 / 0 / 0 | ✅ |
| M3 iHH_A / iHH_D / iHS / freq | 10 / 40 / −1.38629436 / 0.5 | identical | ✅ |
| S1 rehh ln ratio (arithmetic anchor) | −1.97856927 | −1.97856927 | ✅ |
| Scan 2-of-4 extreme, win 4 | proportion 0.5, count 2 | 0.5, 2 | ✅ |

### Variant/delegate consistency
The pre-existing non-canonical `CalculateIHS(ehh0,ehh1,positions)` and region `ScanForSelection(...)` overloads (MCP layer) are out of scope and untouched; the canonical haplotype overloads are self-consistent (INV-04 sign symmetry property test passes).

### Test quality audit (against sources, not code)
Two defects found and **fixed this session** in `PopulationGeneticsAnalyzer_SelectionSignature_Tests.cs`:
1. **S1 was a tautology** — the old `CalculateIHS_RehhRatio_MatchesReference` computed `Math.Log(ihhA/ihhD)` *inside the test* and asserted it equalled the same constant; it never invoked `CalculateIHS`, so it would pass against any implementation (even one with the sign swapped). Rewritten as `CalculateIHS_VoightSignConvention_AncestralOverDerived`: drives `CalculateIHS` end-to-end on an asymmetric panel and asserts `UnstandardizedIHS == ln(IhhAncestral/IhhDerived)` **and** `!= ln(IhhDerived/IhhAncestral)`, locking the Voight sign; retains the rehh −1.97857 constant as a documented arithmetic anchor.
2. **M5 used a one-sided bound** (`Is.LessThan(0.0)`) where the exact value is known (same panel as M3). Tightened to `Is.EqualTo(ln(10/40)).Within(1e-10)`.

Negative-mutation check: flipping the implementation to `ln(iHH_D/iHH_A)` now fails 4 selection tests (M3, M5, S1×2/C1) — previously S1 would have passed the swap. All other tests check exact sourced values within 1e-10, cover every public method/overload and all Stage-A edge/error branches (null, empty, monomorphic, length mismatch, invalid allele, coreIndex OOR, bad binCount/windowSize, singleton bin, two-bin independence, sign-symmetry property).

Remaining notes (not defects): M6/M7/S2 lock the **N−1 sample-SD + 0.05-bin** assumptions rather than a Voight-pinned value (Voight underspecifies the SD estimator and bin width); these are explicitly documented assumptions affecting standardized magnitude only, not the canonical unstandardized score. This is why Stage B is PASS-WITH-NOTES rather than plain PASS.

### Findings / defects
- Code: none. Implementation faithfully realises the validated formulas and edge cases.
- Tests: two test-quality defects (tautological S1, one-sided M5) — both fixed and verified.

## Verdict & follow-ups
- **Stage A: PASS.** Formulas, sign convention, truncation, standardization, and scan criterion all confirmed against Voight (2006) and selscan (2014) retrieved this session; worked numbers reproduced independently.
- **Stage B: PASS-WITH-NOTES.** Code is correct; two test-quality defects fixed in-session. Remaining notes are the two documented standardization assumptions (N−1 SD, 0.05 bin width) that affect standardized magnitude only.
- **End-state: CLEAN.** `dotnet build` 0 errors; full unfiltered suite 6548 passed / 0 failed.
- **Test-quality gate: PASS** (after fixes; sourced exact values, no green-washing, all branches covered, honest green on the full suite).
