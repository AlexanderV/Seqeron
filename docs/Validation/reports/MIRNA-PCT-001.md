# Validation Report: MIRNA-PCT-001 — TargetScan PCT (branch-length score → logistic)

- **Validated:** 2026-06-25   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.ComputeBranchLengthScore(PhyloNode, IReadOnlyCollection<string>)` and `MiRnaAnalyzer.PctFromBranchLength(double, PctSigmoidParameters)`; wired into `ScoreTargetSiteContextPlusPlus` via `PctContribution` when a `PctConservation` input is supplied (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`).
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Authoritative sources opened (this session)

1. **Friedman, Farh, Burge & Bartel (2009)** "Most mammalian mRNAs are conserved targets of microRNAs", *Genome Research* 19:92–105. Defines the **branch-length score (Bls)** of a conserved site = "the total branch length in the phylogenetic tree connecting the subset of species having the [seed] sequence perfectly aligned", and PCT = P_CT as a signal-vs-background estimate of the probability that conservation is due to selective maintenance of targeting.
2. **TargetScan 7.0 reference Perl `targetscan_70_BL_PCT.pl`** (nsoranzo/targetscan GitHub mirror, raw source retrieved this session). Subroutines `calculatePCTthisBL` (the Bls→PCT logistic) and `getBranchLength` (tree traversal), plus `readAllPCTdata` (parameter file parsing).
3. **`PCT_parameters/8mer_PCT_parameters.txt`** (same mirror, raw bytes retrieved). Per-miRNA-family per-site-type b0..b3 tables.

## Stage A — Description

### Formula check (verbatim vs reference Perl)
`calculatePCTthisBL`, retrieved verbatim:
```perl
my $pct = $b0 + ( $b1 / (1 + $eConstant ** ( (0 - $b2) * $BL + $b3)));
$pct = sprintf ("%.4f", $pct);
if ($pct < 0) { $pct = "0.0"; }
```
with `$eConstant = 2.71828182845904523536` (Euler's e). This is exactly
**PCT(Bls) = b0 + b1 / (1 + e^(−b2·Bls + b3))**, truncated at 0 — the logistic the unit claims.
The repo XML-doc and code comment quote this Perl line verbatim and `PctFromBranchLength` realises it.
(The Perl additionally `sprintf "%.4f"` rounds to 4 dp before the negativity test; the library keeps full
double precision — a *more* precise, monotone-preserving choice, documented below as a benign divergence.)

### Parameter handling (verbatim vs reference Perl)
`readAllPCTdata` reads three files (`8mer`, `7mer_m8`, `7mer_a1` → site types), and for each data row
assigns `$familyPlusType2coeff{family}{siteType}{$i-1} = $f[$i]` for `$i = 1..4`. Confirmed against raw bytes:
the data rows have **5 tab fields** (`family_ID`, then four values) — note the file *header* lists six labels
(`family_ID #kmer B_0 B_1 B_2 B_3`) but `#kmer` does **not** appear in the data rows, so the four values are
**b0=col2, b1=col3, b2=col4, b3=col5**. The library does not bundle these tables (citation-required, documented
boundary); the caller supplies `PctSigmoidParameters(B0,B1,B2,B3)` in exactly this b0..b3 order. ✅

### Branch-length score (verbatim vs reference Perl)
`getBranchLength` climbs from each aligned species toward the root accumulating `branch_length`, detecting the
lowest common ancestor with previously-visited species and summing `$ref_cumul_dist{LCA} + $cumul_dist`. This is
the standard **total branch length of the minimal (Steiner) subtree connecting the aligned species**; a single
species (`$#species == 0`) yields BL = 0. The C# `ComputeBranchLengthScore` computes the same quantity by a
post-order walk that adds an edge iff there is an aligned species on **both** sides of it. ✅

### Edge-case semantics (all sourced)
- **Bls = 0** for 0 or 1 conserved species (no connecting subtree) — Perl `$#species == 0 ⇒ 0`. ✅
- **Bls = 0 logistic floor**: PCT(0) = b0 + b1/(1+e^b3), used as-is unless negative. ✅
- **Negative logistic output truncated to 0** — Perl `if ($pct < 0) {0.0}`. ✅
- **Monotonicity**: for b1>0, b2>0, PCT strictly increases in Bls (more conservation ⇒ higher PCT), the unit's
  contract invariant M and Friedman's biological claim. ✅

### Independent cross-check (exact numbers, hand/Python reference)
Worked Newick tree `((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);`:

| Conserved set | Hand-derived Bls | Library | Match |
|---|---|---|---|
| {A,B} | 1.0+2.0 = **3.0** (internal (A,B)→root edge not counted) | 3.0 | ✅ |
| {A,C} | 1.0+0.5+4.0+1.5 = **7.0** | 7.0 | ✅ |
| {A,B,C,D} | 1+2+1.5+3+0.5+4 = **12.0** | 12.0 | ✅ |
| {A} | single species ⇒ **0.0** | 0.0 | ✅ |

PCT logistic (sigmoid (b0,b1,b2,b3)=(0,1,1,0), independent Python `b0+b1/(1+e^(-b2·BL+b3))`):

| Bls | PCT (reference) | Library | Match |
|---|---|---|---|
| 3.0 | 1/(1+e⁻³) = **0.9525741268224334** | 0.952574126822433 | ✅ |
| 7.0 | 1/(1+e⁻⁷) = **0.9990889488055994** | 0.999088948805599 | ✅ |
| 0.0 | 1/(1+e⁰) = **0.5** (untruncated floor) | 0.5 | ✅ |

Truncation: (b0,b1,b2,b3)=(−0.5,0.3,1,5), Bls=0 ⇒ raw = −0.4979921447… < 0 ⇒ **0.0**. ✅

Real published params (8mer, family **UCCCUUU** / miR-30): b0=−1.15890235160174, b1=1.61944269599563,
b2=11.8180287443514, b3=36.260388031391 → PCT(3.0)=**0.0** (truncated, steep curve), PCT(4.0)=**0.460513612719198**,
PCT(5.0)=**0.46054034419686185**; strictly non-decreasing. C# `Math.Exp` reproduces these within 1e-9. ✅

**Stage A findings:** formula, BLS definition, truncation, and parameter ordering all match Friedman 2009 and the
TargetScan reference Perl. PASS.

## Stage B — Implementation

- **Code path:** `MiRnaAnalyzer.cs` — `PctFromBranchLength` (L1506–1512): `B0 + B1/(1 + Math.Exp(-B2·Bls + B3))`,
  then `pct < 0 ? 0 : pct`; `ComputeBranchLengthScore` (L1435–1449) + `AccumulateConnectingBranchLength` (L1456–1486)
  + `CountPresentInSubtree` (L1488–1496); `PctContribution` (L1517–1527) selects the per-site-type PCT row (8mer /
  7mer-m8 / 7mer-A1 / 6mer coeff+min+max from `Agarwal_2015_parameters.txt`) and min-max scales the PCT.
- **Formula realised correctly:** yes — the logistic line and the truncation each carry the verbatim Perl source
  as a comment and match it byte-for-byte in meaning. The BLS post-order "edge counted iff aligned species exist
  on both sides" is mathematically the Steiner-subtree total used by `getBranchLength`; verified equal on the
  worked tree and on a nested topology `(((A,B),C),D)` ({A,B}=3.0, {A,C}=3.0). `CountPresentInSubtree` is O(n²)
  but correct (not a defect for the small species sets in this domain).
- **Cross-verification:** every number in the Stage-A tables was recomputed against the running code (tests below).
- **Variant/delegate consistency:** the public static `ComputeBranchLengthScore`/`PctFromBranchLength` and the
  `ScoreTargetSiteContextPlusPlus` integration (Bls→PCT→`PctContribution`) agree (CTX-PCT-005/006/007).

### Test-quality audit
Existing fixture (`MiRnaAnalyzer_TargetPrediction_Tests.cs`, region CTX-PCT) had: BLS hand-derived cases
(CTX-PCT-001/002), the logistic worked value (003), truncation (004), and 8mer/7mer-m8 contribution wiring
(005/006). All expected values trace to Friedman 2009 / the Perl / independent computation — **not** code echoes.

Gaps found and **fixed this session** (7 new tests added; full suite re-run green):
- **CTX-PCT-007** — adds the missing **7mer-A1** (coeff −0.048, max 0.449) and **6mer** (coeff 0.005, max 0.193)
  PCT-parameter branches of `PctContribution` (previously only 8mer/7mer-m8 were exercised).
- **CTX-PCT-008** — **monotonicity** invariant (contract M): PCT strictly increases over Bls ∈ {0,…,8}; floor PCT(0)=0.5.
- **CTX-PCT-009** — **Bls=0 untruncated logistic floor**: (0.1,0.8,2,4) ⇒ 0.1+0.8/(1+e⁴)=0.11438896796967325.
- **CTX-PCT-010** — cross-check vs the **real published** miR-30 8mer parameter row (3 BL values, independent ref).

**Stage B findings:** no implementation defect. Test suite strengthened to cover all four site-type parameter
paths, the headline monotonicity invariant, the untruncated floor, and a real published-parameter anchor. PASS.

## Documented boundary (NOT a limitation)
TargetScan's compiled per-family/per-site-type PCT tables and the multi-species alignment+tree are
**caller-supplied** (`PctConservation` carries the tree, the conserved-species set, and the b0..b3 sigmoid
parameters). The library bundles only the published *equation* (BLS + logistic) and the Agarwal-2015 PCT
coefficient. Per the protocol and the campaign brief this is an acceptable, documented boundary: the BLS
computation and the PCT logistic both verify against Friedman 2009 / the reference Perl.

## Verdict & follow-ups
**Stage A ✅ PASS · Stage B ✅ PASS · State ✅ CLEAN.** No defects. One benign, documented divergence: the library
keeps full double precision where the Perl rounds PCT to 4 dp via `sprintf "%.4f"` before the negativity check —
strictly more precise and monotonicity-preserving. Full unfiltered `dotnet test Seqeron.sln -c Debug` = 0 failed
(Seqeron.Genomics.Tests 18769 passed; PCT fixture 56 passed, +7 new), 0 warnings on the changed test file.
