# Golden regression set — Seqeron Claude Code skills

A curated set of **hard, realistic biological tasks**, each paired with the **expected routing**
(which skill(s) should fire), the **expected tool/Method-ID pipeline** (dual-mode: MCP tool calls +
C# `Method ID`s), the **rigor checkpoints** `bio-rigor` must enforce, and an **expected-shape output**.

This is the Phase-3 acceptance aid from [`../STRATEGY.md`](../STRATEGY.md) §7 ("golden-набір складних
задач як регрес-тести скілів"). It documents what "good" looks like so a human — or a future automated
eval — can confirm the 10 skills still route and orchestrate correctly.

## What this is (and is NOT)

- **It IS** a **manual / eval regression aid**: a frozen reference of the expected behaviour for a
  representative slice of the skill catalog.
- **It is NOT** an automated test. The skills are **markdown routers**, not code — there is no unit
  test to run. Nothing here executes; the "pass" is a reviewer's judgement against the expected
  behaviour below.
- It does **not** duplicate schemas or algorithm formulas — those stay in the source of truth
  ([`docs/mcp/tools/`](../../../docs/mcp/tools), [`docs/algorithms/`](../../../docs/algorithms),
  [`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md)). Golden tasks only
  cite **verified** tool names + Method IDs and link to the per-tool doc.

## The 10 skills under test

Cross-cutting: [`bio-rigor`](../../../.claude/skills/bio-rigor/SKILL.md),
[`seqeron-discovery`](../../../.claude/skills/seqeron-discovery/SKILL.md).
Domain: [`bio-qc`](../../../.claude/skills/bio-qc/SKILL.md),
[`bio-alignment`](../../../.claude/skills/bio-alignment/SKILL.md),
[`bio-assembly`](../../../.claude/skills/bio-assembly/SKILL.md),
[`bio-annotation`](../../../.claude/skills/bio-annotation/SKILL.md),
[`bio-moldesign`](../../../.claude/skills/bio-moldesign/SKILL.md),
[`bio-phylo-popgen`](../../../.claude/skills/bio-phylo-popgen/SKILL.md),
[`bio-metagenomics`](../../../.claude/skills/bio-metagenomics/SKILL.md),
[`bio-chromosome`](../../../.claude/skills/bio-chromosome/SKILL.md).

## How to use it (reviewer workflow)

1. Pick a task from [`tasks.md`](tasks.md). Paste **Task** (verbatim) into a Claude Code session that
   has the Seqeron skills installed (and, for the MCP path, the relevant servers attached).
2. **Check routing.** The skill(s) named under **Expected skill(s)** should trigger — and no obviously
   wrong skill should. Cross-domain tasks should pull in **several** skills in the right order.
3. **Check the pipeline.** The produced plan / tool calls should match the **Expected pipeline**:
   the same tools, in a compatible order, with the same C# `Method ID` on the dev path. Minor
   reordering of independent steps is fine; a **different tool for the same step is a finding**.
4. **Check rigor.** Confirm the **Rigor checkpoints** hold: tool-only computation (no manual parsing /
   mental math), a provenance block, unit + coordinate-base declarations, cross-checks where noted,
   and the alpha / not-for-clinical-use caveat on decision-relevant output.
5. **Check guarded-unit behaviour.** For tasks marked ⚠ **guarded**, confirm the skill **STOPs and
   reports the limitation** (or, in C#, surfaces the `MinimumMode` / `SeqeronLimitationException`
   contract) instead of silently forcing a number.
6. **Check shape.** The answer's structure should match **Expected-shape output** (table columns,
   provenance block). Exact numbers are **not** graded — the point is routing + orchestration +
   rigor, not recomputing biology by hand.

## Pass criteria (per task)

A task **passes** when all of the following hold:

| # | Criterion |
|---|---|
| 1 | **Routing** — every expected skill triggers; no clearly-wrong skill dominates. |
| 2 | **Pipeline** — the tool chain matches the expected tools + order (dual-mode: MCP names on the tool path, the listed `Method ID`s on the C# path). |
| 3 | **Rigor** — `bio-rigor` is honored: tool-only, provenance block present, units + 0-based/1-based coordinates stated, listed cross-checks performed. |
| 4 | **Envelope** — guarded-unit tasks **STOP + report** (name the limitation, its `MinimumMode`, the alternative); the mode is **not** raised just to force output. |
| 5 | **Disclaimer** — decision-relevant results (variant pathogenicity, real-assay primers/guides, telomere/replicative-age, resistance-gene calls) carry the alpha / not-for-clinical caveat. |
| 6 | **Shape** — output structure matches the expected table/provenance format. |

A task **fails** if a wrong tool is substituted for a graded step, a guarded unit is silently forced,
manual parsing / mental math replaces a tool call, or provenance / caveat is missing where required.

## Coverage

12 tasks. Domains touched: sequence QC, pairwise + MSA alignment, variant calling + classification,
PCR primer design, CRISPR design, phylogeny, popgen (Fst / HWE / Tajima's D), metagenomics diversity,
assembly / N50, chromosome centromere + GC-skew/origin, plus **4 cross-domain** chains (G1, G6, G10,
G12) and **3 guarded-unit** tasks (G2 → `PARSE-FASTQ-001`, G8 → `META-BIN-001`, G11 →
`PROBE-DESIGN-001`), with a variant-pathogenicity clinical-caveat exercise in G4. See the summary table
at the top of [`tasks.md`](tasks.md).

## Provenance of this set

- Every cited tool name + Method ID was verified with `python3 scripts/skills/find-tool.py <kw>`
  against [`docs/mcp/tools/`](../../../docs/mcp/tools) (per-tool doc path is cited in each task).
- Guarded units + their `MinimumMode` verified against
  [`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md).
- Task style follows the real worked workflows in [`docs/mcp/README.md`](../../../docs/mcp/README.md)
  ("cloning insert QC", "PCR primer QC").
- One anomaly noted: `detect_rearrangements` exists on **both** the Analysis server
  (`ComparativeGenomics.DetectRearrangements`) and the Chromosome server
  (`ChromosomeAnalyzer.DetectRearrangements`). Chromosome-scale synteny tasks (G10) use the
  **Chromosome** one, matching `bio-chromosome`.
