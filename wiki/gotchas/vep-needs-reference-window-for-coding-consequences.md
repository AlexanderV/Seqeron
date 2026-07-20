---
type: gotcha
title: "annotate_variants without a reference window returns only the coarse coding term — no missense/synonymous/stop"
tags: [variants, gotcha]
mcp_tools:
  - annotate_variants
  - predict_variant_effect
sources:
  - docs/algorithms/Variants/Variant_Annotation.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# annotate_variants without a reference window returns only the coarse coding term — no missense/synonymous/stop

**The trap.** The consequence engine is **Simplified**, and refining a coding variant into
`missense_variant` / `synonymous_variant` / `stop_gained` / `stop_lost` requires a **reference
sequence window** so it can translate the affected codon. **Without a reference window**,
`Annotate`'s coordinate-only routing **returns only the coarse coding term** (no
synonymous/missense/stop refinement).

**Why it bites.** Call `annotate_variants` with coordinates only and you'll get back a blunt
`coding_sequence_variant`-level term and may conclude "no missense here" when the tool simply had no
sequence to translate. It also follows VEP's **most-specific-term** rule (the lowest-rank applicable
consequence), so a coarse term is a *lack of input*, not a biological finding.

**What to rely on instead.** Pass the `referenceSequence` (and correct `sequenceStart`) so
`PredictFunctionalImpact` can translate the reference vs alternate codon; only then are
missense/synonymous/stop distinctions valid. This is a *Simplified* consequence caller — not full
Ensembl-VEP / SnpEff, and not ACMG pathogenicity. Full model: [[variant-effect-annotation-vep]].
Research-grade, [[research-grade-limitations]].
