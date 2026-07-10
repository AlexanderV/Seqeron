---
type: source
title: "Validation report: MIRNA-TARGET-001 (miRNA target-site prediction — seed-match site detection + TA_3UTR abundance)"
tags: [validation, mirna, governance]
doc_path: docs/Validation/reports/MIRNA-TARGET-001.md
sources:
  - docs/Validation/reports/MIRNA-TARGET-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: MIRNA-TARGET-001

The two-stage **validation write-up** for test unit **MIRNA-TARGET-001** (microRNA
target-site prediction — seed-match site detection plus the newly-added `TA_3UTR`
target-site-abundance capability), a fresh re-validation dated 2026-06-26 that supersedes
the 2026-06-25 entry after the limitation-fix commit 5f2fbd40. This is the *report* artifact
that feeds one row of the [[validation-ledger]]; it records the validator's **verdict** on both
the algorithm description and the shipped code. The two-stage methodology is the
[[validation-protocol]]; the algorithm itself is summarized in
[[mirna-target-site-prediction]]. Distinct from the pre-implementation
[[mirna-target-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS · State: ✅ CLEAN.** One divergence from the literal
Garcia "non-overlapping" requirement was found and **completely fixed in-session** (over-counting
on periodic seed cores), with a lock test added. Target-prediction test class **72 passed**; full
unfiltered `dotnet test Seqeron.sln -c Debug` gives **Seqeron.Genomics.Tests = 18861 passed / 0
failed**, whole solution green, **0 warnings** on the changed project. The residual (PCT params,
SPS, Len_ORF, ORF-8mer caller-supplied; no bundled transcriptome) is an honest, declared scope
boundary, not a defect.

## Stage A — description (algorithm faithfulness)

- Canonical methods: `MiRnaAnalyzer.FindTargetSites` (pre-existing site detection) plus the new
  TA surface `ComputeTa3Utr(MiRna, IEnumerable<string>)` / `CountSeedSites3Utr(...)`, with
  supporting `GetSeedSequence`, `GetReverseComplement`, private `CountSeedSitesInUtr`, and the
  context++ wiring `ScoreTargetSiteContextPlusPlus` / `TaContribution`.
- Formula validated: **`TA_3UTR = log10(N)`**, where **N = total non-overlapping 8mer + 7mer-m8
  + 7mer-A1 sites** of the miRNA seed across the supplied 3′UTR set. Confirmed **verbatim** against
  Garcia et al. 2011 Online Methods (PMC3190056: TA = number of non-overlapping 3′UTR 8mer/7mer-m8/
  7mer-A1 sites; log10 scale; bare 6mer / offset-6mer excluded), Agarwal et al. 2015 Table 1
  (`TA_3UTR` context++ feature, min-max scaled), and the **TargetScan 7.0 reference perl**
  (`read_TA_SPS` / `getAgarwalContribution` — the file value is consumed as-is, i.e. it already is
  log10(count)). Grimson et al. 2007 confirms the canonical site-type definitions.
- Independent cross-check (Python reimplementation written from the Garcia definition, NOT from the
  C# code): CTX-TA-001 synthetic seed `ACGUACG` over three UTRs → **N = 2 ⇒ TA = log10(2) =
  0.301029995663981**; CTX-TA-002 let-7a seed `GAGGUAG` → **N = 5 ⇒ TA = log10(5) =
  0.698970004336019**; both TA→context++ contributions reproduced exactly.
- Edge cases (defined & sourced): empty/null-only UTRs or no sites ⇒ N = 0 ⇒ TA = 0 (log10 floor);
  DNA UTRs normalised T→U; seed < 7 nt / degenerate ⇒ 0; null enumerable ⇒ `ArgumentNullException`;
  overlapping sites counted **non-overlapping** per Garcia (see Stage B).
- **Stage A finding: PASS.** No biological or mathematical error. The lone inaccuracy was the code
  comment's claim that "non-overlapping" was automatically satisfied — resolved in Stage B.

## Stage B — implementation (code review + cross-check)

- Code path: `MiRnaAnalyzer.cs` — `ComputeTa3Utr` (~888), `CountSeedSites3Utr` (~906), private
  `CountSeedSitesInUtr` (~938), `GetSeedSequence` (95), `GetReverseComplement` (296),
  `TaContribution` (1482), `ScoreTargetSiteContextPlusPlus` wiring (806/835); `FindTargetSites` (156)
  reconfirmed unchanged by the fix.
- Realises the description: `CountSeedSites3Utr` derives the seed, takes `seedRC`, `pos8Rc`,
  `sixmerCore`, and counts a site iff the 6mer core matches AND (m8 upstream OR A1 downstream) —
  8mer / 7mer-m8 / 7mer-A1, bare 6mer excluded; `ComputeTa3Utr` returns `Math.Log10(total)` (0 when
  total = 0). `TaContribution` applies the bundled Agarwal per-site-type `(coeff,min,max)` via
  `coeff × (ta−min)/(max−min)` — exactly the perl `getAgarwalContribution`. ✓
- **DEFECT (found & fixed this session):** Garcia requires **non-overlapping** sites; the original
  `CountSeedSitesInUtr` incremented once per qualifying 6mer-core anchor with a plain `for (i++)`
  scan, and the comment's claim that per-anchor site-type mutual-exclusivity guarantees
  non-overlap was **false**. On self-similar / periodic cores the same physical region was
  over-counted (`AAAAAA` vs `AAAAAAAAAA` → 5, should be 1; `ACACAC` vs `ACACACACACA` → 3, should be
  1). **Fix:** greedy left-to-right scan advancing to `i+6` on a counted site so the next core cannot
  overlap; non-qualifying positions still advance by 1. Identical to the prior count on real
  non-self-overlapping seeds (all existing CTX-TA fixtures unchanged: 2, 5, 2, 1); correct on
  periodic cores. Misleading comment + XML `<remarks>` corrected to cite Garcia. **Lock test
  CTX-TA-009** added (periodic seed `GUGUGUU`, core `ACACAC`: `ACACACACACA` ⇒ 1; two cores ≥6 apart
  ⇒ 2).
- Pre-existing surface (`FindTargetSites` reverse-complement targeting, seed scan, nt8/A1
  classification, antiparallel geometry, 0-based-inclusive coords, monotone base scorer, context++
  wiring) re-confirmed green and unchanged.
- Test-quality audit: 9 CTX-TA tests assert exact hand-derived values (log10(2), log10(5), scaled
  contributions, the non-overlapping counts 1 and 2), not code echoes; cover both public methods,
  both overloads' null contract, and the edge cases.

## Findings

- **One genuine defect, fully fixed in-session** (over-counting on periodic seed cores) — greedy
  non-overlapping scan + corrected comments + lock test CTX-TA-009; all existing fixtures unchanged;
  full suite green. Defect logged in FINDINGS_REGISTER.md.
- **Honest residual (declared boundary, NOT a defect):** PCT per-family parameters, SPS, Len_ORF and
  ORF-8mer count remain **caller-supplied** (reported in `OmittedFeatures` when absent), and **no
  default human transcriptome / 3′UTR set is bundled** — the caller must supply the 3′UTR set over
  which abundance is counted (LIMITATIONS.md §3, trimmed by the fix).
- **Runtime enforcement (LimitationPolicy):** the guarded branch — a **partial** context++ score
  (`OmittedFeatures` non-empty) — has minimum access mode `Permissive`; under the default `Moderate`
  it throws `SeqeronLimitationException`. Additive policy layer; the ✅ CLEAN verdict is unchanged.

See the full report at `docs/Validation/reports/MIRNA-TARGET-001.md`.
