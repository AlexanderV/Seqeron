---
type: gotcha
title: "CRISPR guide and off-target scores are heuristics, not calibrated efficacy/risk estimates"
tags: [crispr, gotcha]
mcp_tools:
  - evaluate_guide_rna
  - design_guide_rnas
  - crispr_specificity_score
  - find_off_targets
sources:
  - docs/algorithms/MolTools/Guide_RNA_Design.md
  - docs/algorithms/MolTools/Off_Target_Analysis.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# CRISPR guide and off-target scores are heuristics, not calibrated efficacy/risk estimates

**The trap.** Both scoring surfaces here are **honest heuristics**, not experimentally-calibrated
models. `EvaluateGuideRna` (used by `design_guide_rnas`) is a **composition-penalty** score — the
docs state plainly that **callers must not read the heuristic score as an efficacy estimate**. The
off-target side (`FindOffTargets` / `CalculateSpecificityScore`) is a **position-weighted heuristic**
with test-unit status *Simplified* — not the MIT/CFD experimentally-derived model.

**Why it bites.** If you rank candidate guides by `EvaluateGuideRna` and pick the top one as "most
efficient", or treat `CalculateSpecificityScore` as a validated off-target-risk probability, you are
over-reading a composition/positional heuristic. It will disagree with Doench/Rule-Set-2 (on-target)
and MIT/CFD (off-target) tools. Note also there are **three distinct "seed" notions** on this surface
(10-nt design seed, 12-bp off-target seed, full 20-position weight vector) — don't conflate them.

**What to rely on instead.** Use the scores for **relative screening / filtering** of candidates, not
as absolute efficacy or risk. For decision-grade calls, run an experimentally-calibrated model
(Doench 2016 on-target; MIT/CFD off-target) externally. `find_pam_sites` / `design_guide_rnas`
correctly scan **both strands** for PAMs — that geometry is sound; only the *scores* are heuristic.
Full model: [[crispr-guide-rna-design]]. Research-grade, [[research-grade-limitations]].
