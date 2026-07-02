---
name: seqeron-rna-structure
description: >-
  Predict and analyze RNA secondary structure with Seqeron (MCP tools OR the C#
  API). Use to fold an RNA / compute minimum free energy (MFE, Turner 2004
  nearest-neighbour, kcal/mol), predict the secondary structure (dot-bracket +
  base pairs), enumerate stem-loops / hairpins, detect pseudoknots (crossing
  base pairs), find antiparallel complementary arms (RNA inverted repeats),
  classify base pairs (Watson-Crick / G-U wobble) and test if two bases can
  pair, parse / validate dot-bracket (WUSS) notation, and compute individual
  loop / stem / dangling-end / terminal-mismatch free-energy terms. Triggers:
  "fold this RNA", "predict secondary structure", "minimum free energy of…",
  "MFE of this RNA", "find stem-loops / hairpins", "detect pseudoknots", "is
  this dot-bracket valid", "parse this dot-bracket", "can G and U pair", "what
  base-pair type", "hairpin/internal/bulge/multibranch loop energy". Server:
  analysis (all RnaSecondaryStructure.*).
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-rna-structure — fold RNA, MFE, stem-loops, pseudoknots, dot-bracket

Routing + orchestration skill for the RNA-secondary-structure family on the **analysis**
server (17 tools, one backing class: `RnaSecondaryStructure.*`). It picks the right tool for a
folding / MFE / stem-loop / pseudoknot / dot-bracket question and gives a **dual-mode** recipe
(MCP tool name **and** the equivalent `Seqeron.Genomics` C# `Method ID`).

- **Rigor is delegated.** Tool-only computation, provenance, envelope, cross-check, units, and the
  alpha / not-for-clinical caveat are owned by **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies
  by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server analysis`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/analysis/*.md`;
  algorithm invariants in `docs/algorithms/{RnaStructure,RNA_Structure,RNA_Secondary_Structure}/*.md`.

## Scope boundary — this family vs neighbours

- **This skill OWNS** RNA *secondary structure*: folding/MFE, stem-loops, pseudoknot detection,
  RNA inverted repeats, base-pairing rules, dot-bracket parse/validate, per-motif energy terms.
  **[`bio-annotation`](../bio-annotation/SKILL.md)** name-drops this family shallowly — route here.
- **Pre-miRNA hairpins, miRNA seed / target pairing, RNA-seq, splicing** stay in
  **[`bio-annotation`](../bio-annotation/SKILL.md)** (`find_pre_mirna_hairpins`, `find_mirna_target_sites`, …).
- **DNA inverted repeats / palindromes** → `find_inverted_repeats` (core/analysis) via
  `bio-assembly` / `bio-annotation`; `find_rna_inverted_repeats` here uses **RNA** complement rules.
- **No base-pair-probability / partition-function MCP tool exists.** The McCaskill partition
  function is a documented C#-only algorithm (`docs/algorithms/RNA_Structure/`), not an exposed
  tool. For pairing feasibility use `can_pair` / `base_pair_type`; for structure use `predict_rna_structure`.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Fold an RNA → dot-bracket + pairs + stem-loops + pseudoknots + total MFE | `predict_rna_structure` / `RnaSecondaryStructure.PredictStructure` |
| Just the MFE number (kcal/mol) | `minimum_free_energy` / `RnaSecondaryStructure.CalculateMinimumFreeEnergy` |
| Enumerate hairpin stem-loop candidates | `find_stem_loops` / `RnaSecondaryStructure.FindStemLoops` |
| Antiparallel complementary arms (potential stems) | `find_rna_inverted_repeats` / `RnaSecondaryStructure.FindInvertedRepeats` |
| Pseudoknots = crossing base pairs | `detect_pseudoknots` / `RnaSecondaryStructure.DetectPseudoknots` |
| Can two bases pair? | `can_pair` / `RnaSecondaryStructure.CanPair` |
| Base-pair type (WatsonCrick / Wobble / null) | `base_pair_type` / `RnaSecondaryStructure.GetBasePairType` |
| RNA complement of one base | `rna_complement_base` / `RnaSecondaryStructure.GetComplement` |
| Parse dot-bracket → base-pair index tuples | `parse_dot_bracket` / `RnaSecondaryStructure.ParseDotBracket` |
| Validate dot-bracket (balanced + family-matched) | `validate_dot_bracket` / `RnaSecondaryStructure.ValidateDotBracket` |
| Stem stacking energy | `stem_energy` / `RnaSecondaryStructure.CalculateStemEnergy` |
| Hairpin-loop energy | `hairpin_loop_energy` / `RnaSecondaryStructure.CalculateHairpinLoopEnergy` |
| Internal-loop energy | `internal_loop_energy` / `RnaSecondaryStructure.CalculateInternalLoopEnergy` |
| Bulge-loop energy | `bulge_loop_energy` / `RnaSecondaryStructure.CalculateBulgeLoopEnergy` |
| Multibranch-loop energy | `multibranch_loop_energy` / `RnaSecondaryStructure.CalculateMultibranchLoopEnergy` |
| Dangling-end energy | `dangling_end_energy` / `RnaSecondaryStructure.GetDanglingEndEnergy` |
| Terminal-mismatch energy | `terminal_mismatch_energy` / `RnaSecondaryStructure.GetTerminalMismatchEnergy` |

## Canonical dual-mode pipelines

### (a) Fold an RNA → structure + MFE (the main entry point)
1. **[MCP]** `predict_rna_structure`(rnaSequence, minStemLength?=3, minLoopSize?=3, maxLoopSize?=10) → echoed sequence, `dotBracket`, `basePairs`, `stemLoops`, `pseudoknots`, total MFE.
2. **[MCP]** cross-check the energy: `minimum_free_energy`(rnaSequence, minLoopSize?=3) → MFE (kcal/mol; 0 if no structure forms).
- **[C# API]** `RnaSecondaryStructure.PredictStructure(...)` → `.CalculateMinimumFreeEnergy(...)`.
- Note: prediction is **greedy non-overlapping stem-loops** (nested pairs only) — a Zuker-optimal MFE structure is *not* guaranteed; `minimum_free_energy` is the DP energy bound.
```
Provenance
1) predict_rna_structure(rnaSequence,minStemLength=3,minLoopSize=3,maxLoopSize=10) → dotBracket,basePairs,stemLoops,pseudoknots,mfe
2) minimum_free_energy(rnaSequence,minLoopSize=3) → mfe (Zuker DP, Turner 2004) — cross-check magnitude
Envelope: RNA-STRUCT-001 documented-only (see STOP rule). Caveat: alpha — validate before use.
```

### (b) Enumerate stem-loops / hairpins directly
1. **[MCP]** `find_stem_loops`(rnaSequence, minStemLength?=3, minLoopSize?=3, maxLoopSize?=10, allowWobble?=true) → per-candidate stem bounds/pairs/energy, loop type/bounds/size, dot-bracket, Turner total energy.
2. **[MCP]** corroborate arm placement independently: `find_rna_inverted_repeats`(sequence, minLength?=4, minSpacing?=3, maxSpacing?=100) → left/right arm spans + arm length.
- **[C# API]** `RnaSecondaryStructure.FindStemLoops(...)` · `.FindInvertedRepeats(...)`.

### (c) Pseudoknot detection from a set of base pairs
1. **[MCP]** get pairs (either `parse_dot_bracket` on a known structure, or the `basePairs` from `predict_rna_structure`).
2. **[MCP]** `detect_pseudoknots`(basePairs=[{position1,position2,base1,base2,type}]) → crossing pair-of-pairs (`i<k<j<l`), reported once each.
- **[C# API]** `RnaSecondaryStructure.ParseDotBracket(...)` → `.DetectPseudoknots(pairs)`.
- Envelope: only csr-PK-grammar knots are in scope — see the STOP rule (RNA-STRUCT-001).

### (d) Dot-bracket handling (validate → parse → pairs)
1. **[MCP]** `validate_dot_bracket`(dotBracket) → balanced & family-matched? (`(]` is rejected).
2. **[MCP]** `parse_dot_bracket`(dotBracket) → `(openingIndex, closingIndex)` tuples (multi-family: `()[]{}<>` + letter pairs).
- **[C# API]** `RnaSecondaryStructure.ValidateDotBracket(...)` → `.ParseDotBracket(...)`. **Validate before parse.**

### (e) Manual base-pairing / energy-term inspection
- **[MCP]** `can_pair`(base1, base2) → bool (A-U, G-C, G-U). `base_pair_type`(base1, base2) → WatsonCrick / Wobble / null. `rna_complement_base`(base) → complement.
- **[MCP]** per-motif energies (all Turner 2004, kcal/mol): `stem_energy` · `hairpin_loop_energy` · `internal_loop_energy` · `bulge_loop_energy` · `multibranch_loop_energy` · `dangling_end_energy` · `terminal_mismatch_energy`.
- **[C# API]** `RnaSecondaryStructure.CanPair` · `.GetBasePairType` · `.GetComplement` · `.CalculateStemEnergy` · `.CalculateHairpinLoopEnergy` · … (see tool-map).

## Envelope — STOP rule (guarded unit in scope)

- **RNA-STRUCT-001** is a **documented-only** limitation (no runtime guard — nothing throws; the
  result is *exact for the stated NN-energy model / csr-PK grammar*, and the shortfall is
  undetectable per call). Two irreducible gaps: (1) tertiary-stabilised knots (e.g. BWYV / PDB 437D)
  are **not** recoverable as the MFE structure by *any* nearest-neighbour model; (2) pseudoknot
  classes **outside the csr-PK grammar** (kissing hairpins / loop-loop, triple-crossing / chained,
  non-canonical bulged helices) cannot be faithfully detected — the loop-loop energy was never
  measured. If a task asks to *recover a tertiary/pseudoknotted native fold* or *detect
  non-csr-PK pseudoknots*, **STOP** and report the limitation (cite `docs/Validation/LIMITATIONS.md`
  RNA-STRUCT-001) — do **not** present the NN-model output as the true structure. The general
  folding / MFE / stem-loop tools are fine within their model; just flag the model floor.

## End-to-end grounded example (extends `docs/mcp/README.md`)

**Task.** Given one RNA sequence, fold it, then independently verify each layer:
(1) fold → dot-bracket + pairs + MFE, (2) confirm the dot-bracket is valid and re-parse it to pairs,
(3) re-derive MFE via the DP, (4) check whether any predicted pairs cross (pseudoknots), (5) flag the model floor.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `predict_rna_structure`(rnaSequence) → `dotBracket`, `basePairs`, `stemLoops`, `pseudoknots`, `mfe`. (`RnaSecondaryStructure.PredictStructure`)
2. `validate_dot_bracket`(dotBracket) → true; then `parse_dot_bracket`(dotBracket) → pairs should match step 1's `basePairs`. (`.ValidateDotBracket` → `.ParseDotBracket`)
3. `minimum_free_energy`(rnaSequence) → MFE; magnitude should track step 1's `mfe`. (`.CalculateMinimumFreeEnergy`)
4. `detect_pseudoknots`(basePairs from step 1) → prediction is nested-only, so expect none; a nonempty result contradicts step 1 and must be reconciled. (`.DetectPseudoknots`)
```
Provenance
1) predict_rna_structure(rnaSequence) → dotBracket,basePairs,stemLoops,pseudoknots,mfe
2) validate_dot_bracket(dotBracket)=true; parse_dot_bracket(dotBracket) ≡ step-1 basePairs (cross-check)
3) minimum_free_energy(rnaSequence) → mfe consistent with step-1 mfe (independent DP path)
4) detect_pseudoknots(step-1 basePairs) → expect ∅ (greedy prediction is nested-only)
Envelope: RNA-STRUCT-001 documented-only — NN-model / csr-PK result; not a guaranteed native fold.
Caveat: alpha software; not for clinical use — independently validate before relying on any call.
```

## Reference

- **This family's tool map (all 17 — curated index; NOT in domain-map.json, so there is NO
  generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`Minimum_Free_Energy.md`](../../../docs/algorithms/RnaStructure/Minimum_Free_Energy.md) ·
  [`Pseudoknot_Detection.md`](../../../docs/algorithms/RnaStructure/Pseudoknot_Detection.md) ·
  [`Dot_Bracket_Notation.md`](../../../docs/algorithms/RnaStructure/Dot_Bracket_Notation.md) ·
  [`RNA_Base_Pairing.md`](../../../docs/algorithms/RnaStructure/RNA_Base_Pairing.md) ·
  [`RNA_Secondary_Structure.md`](../../../docs/algorithms/RNA_Structure/RNA_Secondary_Structure.md) ·
  [`RNA_Stemloop.md`](../../../docs/algorithms/RNA_Structure/RNA_Stemloop.md) ·
  [`RNA_Free_Energy.md`](../../../docs/algorithms/RNA_Structure/RNA_Free_Energy.md) ·
  [`Inverted_Repeats.md`](../../../docs/algorithms/RNA_Secondary_Structure/Inverted_Repeats.md)
- **Operating envelope / guarded unit:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (RNA-STRUCT-001)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (overlap: it name-drops RNA structure; pre-miRNA / miRNA / RNA-seq / splicing live there).
