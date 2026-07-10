---
type: source
title: "Validation report: MIRNA-PAIR-001 (miRNA-target pairing — CanPair / wobble / reverse complement / antiparallel duplex)"
tags: [validation, mirna, governance]
doc_path: docs/Validation/reports/MIRNA-PAIR-001.md
sources:
  - docs/Validation/reports/MIRNA-PAIR-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: MIRNA-PAIR-001

The two-stage **validation write-up** for test unit **MIRNA-PAIR-001** (MiRNA-Target
Pairing Analysis), validated 2026-06-15. This is the *report* artifact that feeds one row
of the [[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The base-pairing rule, reverse complement, and
antiparallel duplex are summarized in [[rna-base-pairing]]; the two-stage methodology is
the [[validation-protocol]]. Distinct from the pre-implementation
[[mirna-pair-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End-state: CLEAN.** A green-washed /
contract-violating test plus a matching code gap (`CanPair`/`IsWobblePair` did not honour
the documented DNA-T→U contract) were found and completely fixed in-session; the code was
changed to match the description, not vice versa. Full unfiltered suite **6543 passed / 0
failed**; `dotnet build` 0 errors; changed files add 0 new warnings. Gate PASS.

## Stage A — description (algorithm faithfulness)

- Canonical methods: `MiRnaAnalyzer.AlignMiRnaToTarget(string,string)`,
  `GetReverseComplement(string)`, `CanPair(char,char)`, `IsWobblePair(char,char)` (+
  private `CalculateDuplexEnergy`, `StackingEnergies`).
- Sources opened: **Wikipedia "Wobble base pair"** (Crick 1966) + **ScienceDirect** —
  Watson-Crick = A·U, G·C; principal RNA wobble = G·U (not WC, ~as stable as A-U),
  confirming `CanPair` {A-U,U-A,G-C,C-G,G-U,U-G} and `IsWobblePair` {G-U,U-G};
  **miRBase/RNAcentral** hsa-let-7a-5p `UGAGGUAGUAGGUUGUAUAGUU`, seed (2-8) `GAGGUAG`, WC
  reverse complement `CUACCUC` (hand-verified); **Xia et al. (1998) / Turner 2004** — full
  10-stack Watson-Crick nearest-neighbor ΔG°37 set retrieved; **NNDB Turner 2004 GU page**
  (tandem-wobble values confirmed via search snippets; HTML 404s to the fetcher).
- Formula check: pairing classification (WC / wobble / mismatch), antiparallel orientation
  (miRNA `i` pairs target `len−1−i`, verified in code), and simplified Turner-2004 stacking
  ΔG (sign/ordering only, magnitude intentionally simplified) all match sources. All 16 WC
  `StackingEnergies` entries match Xia 1998 exactly.
- Edge cases: empty/null → empty duplex / `""`; unequal lengths → overlap = min(len); DNA
  T → U; unknown base → `N`, never pairs.
- Independent cross-check reproduced every oracle (RC(`GAGGUAG`)=`CUACCUC`; AAAA/UUUU 4 WC
  `||||`; GGGG/UUUU 4 wobble `::::`; AAAA/AAAA 4 mismatch; AGGU/AUCG 2 WC/1 wobble/1
  mismatch ` |:|`; GC-stem ΔG −20.76 ≤ 0). **Stage A: PASS.**

## Stage B — implementation (code review + cross-check)

- Code path: `Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs` — `CanPair` L322,
  `IsWobblePair` L335, `GetReverseComplement` L294, `AlignMiRnaToTarget` L350,
  `CalculateDuplexEnergy` L415, `StackingEnergies` L478. Formula realised correctly
  (antiparallel index, WC/wobble/mismatch, `|`/`:`/space symbols, ungapped overlap,
  nearest-neighbor sum); all 16 WC stacking constants match Xia 1998.
- **Defect (fixed this session):** `CanPair`/`IsWobblePair` did **not** honour the
  documented "DNA T treated as U" contract (doc §3.1, §6.1; TestSpec M4 "A-T → true").
  `CanPair('A','T')` returned **false** while `AlignMiRnaToTarget`/`GetReverseComplement`
  (which `.Replace('T','U')`) treat T as U — an inconsistency the MCP-exposed `CanPair`
  tool inherited. The canonical test `CanPair_LowercaseAndDnaT_ReturnsTrue` was
  **green-washed**: name + spec promised DNA-T handling but the body silently omitted the
  `A-T` assertion.
- **Fix (code):** added a private `NormalizeBase` (uppercase + T→U) used by `CanPair` /
  `IsWobblePair` so they honour the contract; internal callers pass pre-normalised
  sequences (unchanged there) — the change only *adds* correct T-handling at the
  public/MCP boundary. MCP `AnnotationTools.CanPair`/`IsWobblePair` delegate and inherit
  the fix.
- **Fix (tests):** `CanPair_LowercaseAndDnaT_ReturnsTrue` now asserts `('A','T')`,
  `('T','A')`, `('G','T')`, `('a','t')` == true; added `IsWobblePair_DnaT_TreatedAsU`
  (M4b); rewrote M15 (`AlignMiRnaToTarget_CountInvariant…`) from a vacuous perfect-
  complement input to a genuinely mixed duplex `AGGU`/`AUCG` (2 WC / 1 wobble / 1
  mismatch, `" |:|"`).

## Findings

- **FIXED:** T-normalisation gap in `CanPair`/`IsWobblePair` + green-washed M4 test +
  vacuous M15 invariant test — all corrected and locked with sourced assertions
  (Crick wobble, Agarwal/Wikipedia WC, Lewis RC, Xia 1998 stacking). Logged in
  FINDINGS_REGISTER. Test-quality audit (HARD gate) PASS; no assertions weakened, no
  tolerances widened, no skips. Free-energy *magnitude* remains intentionally simplified
  (sign-only validated, NNDB 404 documented limitation), not a coverage gap. Honest
  green: full suite 6543/0.
