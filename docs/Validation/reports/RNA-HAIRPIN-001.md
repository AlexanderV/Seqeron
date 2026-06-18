# Validation Report: RNA-HAIRPIN-001 — Hairpin Loop and Stem Free-Energy Calculation (Turner 2004)

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CalculateHairpinLoopEnergy(string, char, char, bool)`, `RnaSecondaryStructure.CalculateStemEnergy(string, IReadOnlyList<BasePair>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

This unit computes ΔG°37 (kcal/mol) of the two elementary motifs of an RNA stem-loop using the
Turner 2004 nearest-neighbor model (Mathews et al. 2004, PNAS; NNDB): a hairpin loop (additive
initiation + terminal mismatch + sequence bonuses + all-C penalty, with special-loop overrides) and
a stem (nearest-neighbor stacking + per-AU/GU-end penalty). The "RNA hairpin / stem-loop detection"
phrasing in the task maps here to thermodynamic scoring of the stem (stem base-pairing/stacking),
the loop, and the ≥3-nt minimum-loop rule — all of which this unit encodes.

## Stage A — Description

### Sources opened & what they confirm (all retrieved THIS session)

| Source | Retrieval | Confirms |
|--------|-----------|----------|
| NNDB Turner 2004 `hairpin.html` | Wayback `20240709061712` via curl (HTTP 200) | Formula for >3 nt and =3 nt verbatim; loops <3 nt "prohibited"; special-GU "GU closing pair (not UG) preceded by two Gs"; all-C linear `An+B`; 3-nt all-C penalty; `init(n>9)=init(9)+1.75 RT ln(n/9)` |
| NNDB `loop.txt` | Wayback `20240709061712` via curl (HTTP 200) | Hairpin initiation 3→5.4, 4→5.6, 5→5.7, 6→5.4, 7→6.0, 8→5.5, 9→6.4, 10→6.5 … 30→7.7 (matches impl table exactly) |
| NNDB `hairpin-example-1.html` | Wayback `20211027121741` via curl (HTTP 200) | Worked: stacks −2.11, −2.24, −2.11; AU end +0.45; terminal mismatch (AU/AA) −0.8; init(6) 5.4; **total −1.4**; "for unimolecular structures the helical intermolecular initiation does not appear" |
| NNDB `hairpin-example-2.html` | Wayback `20211027135550` via curl (HTTP 200) | Worked: same 3 stacks + 0.45; terminal mismatch (AU/GG) −0.8; **GG first mismatch −0.8**; init(5) 5.7; **total −1.9** |
| NNDB `wc-parameters.html` | Wayback `20240715153500` via curl (HTTP 200) | **"Per AU end +0.45 ± 0.04"**; helix initiation +4.09; symmetry correction (excluded for unimolecular) |
| NNDB `gu-parameters.html` | Wayback `20240714175709` via curl (HTTP 200) | **"Per GU end +0.45 ± 0.04"**; note b: **5'GGUC3'/3'CUGG5' = −4.12 ± 0.54**; 5'GU3'/3'UG5' = +1.29 |
| NNDB `tstack.txt` (terminal mismatch matrix) | Wayback `20240709063234` via curl (HTTP 200) | Decoded: closing A-U / A·A = −0.8; A-U / G·G = −0.8; C-G / U·U = −1.2; G-C / C·C = −0.7; G-C / A·A = −1.1 |
| NNDB `triloop.txt`, `tloop.txt`, `hexaloop.txt` | Wayback snapshots via curl (HTTP 200) | CAACG=6.8, GUUAC=6.9; CCUCGG=2.5, CUACGG=2.8, CAACGG=5.5; ACAGUGUU=1.8, ACAGUACU=2.8, ACAGUGAU=3.6 |
| **ViennaRNA `rna_turner2004.par`** (independent reference impl) | GitHub raw (HTTP 200, 381 KB) | Special-loop totals identical (÷100); all-C `MultipleCA37=30`(0.3), `MultipleCB37=160`(1.6); `TerminalAU37=50`(0.5); `TripleC37=100` |
| **RNAstructure `miscloop.dat`** (Mathews-lab Turner 2004 data, `maxhwardg/advanced_multiloops`) | GitHub raw (HTTP 200) | terminal AU 0.5; **GGG/special-GU −2.2**; c-slope 0.3; c-intercept 1.6; **c-of-3 (C3) 1.5**; intermolecular init 4.1 |

### Formula check
- Hairpin (>3 nt) and (=3 nt) formulas in the doc/spec/impl are **verbatim** the NNDB `hairpin.html` equations.
- Stem = Σ nearest-neighbor stacks (P−1 terms) + 0.45 per AU/GU helix end; intermolecular initiation
  and the self-complementary symmetry correction excluded — explicitly correct for a **unimolecular**
  stem per the NNDB Example-1 note. Confirmed.
- n>9 extrapolation `init(9)+1.75 RT ln(n/9)` with R=1.987 cal/mol·K, T=310.15 K — matches NNDB.

### Edge-case semantics
- Loops <3 nt: source gives **no value** ("prohibited"). Implementation returns a defined sentinel
  100.0 (INV-02) — a documented, sourced design choice, not an unstated behavior.
- 3-nt loops: no first-mismatch term (only initiation + all-C penalty) — confirmed.
- Special tri/tetra/hexaloops override the additive model — confirmed against three independent files.
- special-GU −2.2 is G-U-only (not U-G) — confirmed verbatim.

### Independent cross-check (numbers)
- Example 1 reproduced: loop = 5.4 − 0.8 = **+4.6**; helix = −2.11 −2.24 −2.11 +0.45 = **−6.01**; total **−1.4**.
- Example 2 reproduced: loop = 5.7 − 0.8 − 0.8 = **+4.1**; total **−1.9**.
- ViennaRNA and RNAstructure **independently** confirm the special-loop totals, the all-C A/B, and C3=1.5.

### Findings / divergences (Stage A)
- **Terminal-AU/GU-end penalty value:** the RNAstructure `miscloop.dat` and ViennaRNA both carry **0.50**,
  while the NNDB **web parameters and both worked examples use +0.45**. This is a known internal
  inconsistency in the Turner 2004 publication set (rounded data-file vs. fitted web value). The
  implementation uses **0.45**, which is consistent with the NNDB worked examples it is validated against
  (Example-1 helix = −6.01 only with 0.45). This is a **sourced, internally-consistent choice**, recorded
  as a note rather than a defect. Stage A → **PASS**.

## Stage B — Implementation

- **Code path:** `src/.../RnaSecondaryStructure.cs` — `CalculateHairpinLoopEnergy` (lines 684–759),
  `CalculateStemEnergy` (590–631), parameter tables (122–254).

### Formula realised correctly?
Yes. Special-loop key `c5+loop+c3` lookup first (override); else initiation (table / log-extrapolation /
<3-nt sentinel); for n≥4 adds terminal mismatch + UU/GA (−0.9) + GG (−0.8) + special-GU (−2.2 G-U only);
all-C penalty (C3 1.5 for n=3, `0.3·n+1.6` for n>3). Stem sums P−1 stacks (with the GGUC/CUGG −4.12
special context) + 0.45 per AU/GU end. Every constant traces to a source retrieved this session.

### Cross-verification table recomputed vs code (all match)

| Case | Inputs | Sourced expected | Code |
|------|--------|------------------|------|
| M1 | "AAAAAA",A,U | 5.4−0.8 = 4.6 | 4.6 ✓ |
| M2 | "GAAAG",A,U | 5.7−0.8−0.8 = 4.1 | 4.1 ✓ |
| M3 | "AAC",C,G | CAACG = 6.8 | 6.8 ✓ |
| M4 | "CUCG",C,G | CCUCGG = 2.5 | 2.5 ✓ |
| M5 | "AAA",G,C | init(3) = 5.4 | 5.4 ✓ |
| M6 | "CCC",G,C | 5.4 + C3 1.5 = 6.9 | 6.9 ✓ |
| M7 | "AA",G,C | prohibited → 100.0 sentinel | 100.0 ✓ |
| S1 | "AAAA",G,U flag | −2.2 isolated | −2.2 ✓ |
| S2 | "AAAA",U,G flag | no −2.2 (U-G) | 0 diff ✓ |
| S3 | "UAAU",C,G | 5.6 + (CUUG −1.2) + UU −0.9 = 3.5 | 3.5 ✓ |
| S4 | "CCCC",G,C | 5.6 + (GCCC −0.7) + (0.3·4+1.6) = 7.7 | 7.7 ✓ |
| S6* | "A"×40,G,C | init(40) 8.01 + (GAAC −1.1) = 6.91 | 6.91 ✓ |
| C1 | "CAGUGU",A,U | ACAGUGUU = 1.8 | 1.8 ✓ |
| M8 | C-G,A-U,C-G,A-U | −2.11−2.24−2.11+0.45 = −6.01 | −6.01 ✓ |
| M9 | [] | 0 | 0 ✓ |
| S5 | A-U,G-C,A-U | −2.08−2.35+0.45+0.45 = −3.53 | −3.53 ✓ |
| S7* | G-C,G-U(wobble) | (GG/CU −1.53) + GU-end 0.45 = −1.08 | −1.08 ✓ |
| S8* | G-C,G-U,U-G,C-G | special GGUC/CUGG = −4.12 | −4.12 ✓ |

(* = tests added this session.)

### Variant/delegate consistency
Only the two canonical statics exist (no `*Fast`/instance variants). `GetTerminalMismatchEnergy` /
`GetDanglingEndEnergy` helpers are out of this unit's scope (terminal-mismatch dictionary itself is
exercised through the hairpin path). Determinism (INV-01) covered by C2.

### Test quality audit
Pre-existing fixture had 16 tests; all assert **exact sourced values** with `Within(1e-9)`. Gate actions
taken this session:
- **M7 tightened** from `Is.GreaterThanOrEqualTo(100.0)` (a range that would pass against a wrong large
  value) to `Is.EqualTo(100.0)` — the sentinel is an exact known value (INV-02).
- **Added S6** (n>30 log-extrapolation branch — previously untested formula path).
- **Added S7** (stem GU-end / `Wobble` branch of `IsTerminalAUorGU` — previously untested).
- **Added S8** (GGUC/CUGG special-3-stack context, −4.12 — previously untested branch).
All new expected values trace to sources retrieved this session (NNDB hairpin.html extrapolation,
tstack GAAC, gu-parameters GG/CU and note-b −4.12). No assertion weakened, no tolerance widened, no
test skipped.

### Findings / defects
- **No code defect.** All sourced values reproduce exactly.
- **Note (not a defect):** AU/GU-end penalty 0.45 vs the 0.5 in RNAstructure/ViennaRNA data files — see
  Stage A; the 0.45 matches the NNDB worked examples the unit is validated against. Recorded as a
  documented divergence → Stage B **PASS-WITH-NOTES**.

## Verdict & follow-ups
- **Stage A: PASS.** Biology/maths independently confirmed against NNDB + two reference implementations.
- **Stage B: PASS-WITH-NOTES** (AU/GU-end = 0.45 vs 0.5 data-file value; sourced & internally consistent).
- **End-state: CLEAN.** No defect; test-quality gate satisfied (M7 tightened, 3 branch tests added);
  `dotnet build` 0 errors; full unfiltered suite **6589 passed, 0 failed** (1 pre-existing `[Explicit]`
  benchmark skipped). HAIRPIN fixture now 20 tests, all green.
- **Test-quality gate: PASS.**
