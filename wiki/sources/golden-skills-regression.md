---
type: source
title: "Golden skills regression set (docs/skills/golden/README.md)"
tags: [skills, testing]
doc_path: docs/skills/golden/README.md
sources:
  - docs/skills/golden/README.md
  - docs/skills/golden/tasks.md
source_commit: bcc98f65e7084b4311e01aa4007940c808c035a3
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-18
---

# Golden skills regression set

A curated set of hard, realistic biological tasks used as **manual / eval regression
aids** for the skill catalog — the Phase-3 acceptance instrument from [[skills-strategy]]
§7. Each task is paired with expected routing (which skills should fire), the expected
dual-mode pipeline (MCP tool calls + C# `Method ID`s), the rigor checkpoints `bio-rigor`
must enforce, and an expected-shape output. The task list itself is in
`docs/skills/golden/tasks.md` (coverage-excluded, reference-only).

## What it is and is not

It **is** a frozen reference of expected behaviour for a representative slice of the
catalog. It **is not** an automated test: the skills are markdown routers, not code, so
nothing executes — a "pass" is a reviewer's judgement. It does not duplicate schemas or
formulas; those stay in the source of truth, and it cites only *verified* tool names +
Method IDs.

## The 10 skills under test

The set exercises **10 skills**: two cross-cutting (`bio-rigor`, `seqeron-discovery`)
and eight domain skills (`bio-qc`, `bio-alignment`, `bio-assembly`, `bio-annotation`,
`bio-moldesign`, `bio-phylo-popgen`, `bio-metagenomics`, `bio-chromosome`).

## Reviewer workflow and pass criteria

Paste a task verbatim into a Claude Code session with the skills installed, then check
six things: **routing** (expected skills trigger, no clearly-wrong skill dominates),
**pipeline** (tool chain matches expected tools + order; a different tool for a graded
step is a finding), **rigor** (tool-only computation, provenance block, units +
coordinate base stated, listed cross-checks done), **envelope** (guarded-unit tasks
STOP and report the limitation + its `MinimumMode` instead of forcing a number),
**disclaimer** (decision-relevant results carry the beta / not-for-clinical caveat), and
**shape** (output structure matches). Exact numbers are not graded — the point is
routing + orchestration + rigor.

## Coverage

12 tasks spanning sequence QC, pairwise + MSA alignment, variant calling +
classification, PCR primer design, CRISPR design, phylogeny, popgen (Fst / HWE /
Tajima's D), metagenomics diversity, assembly / N50, and chromosome centromere +
GC-skew/origin — including 4 cross-domain chains (G1, G6, G10, G12) and 3 guarded-unit
tasks (G2 → `PARSE-FASTQ-001`, G8 → `META-BIN-001`, G11 → `PROBE-DESIGN-001`) plus a
variant-pathogenicity clinical-caveat exercise (G4). Every cited tool + Method ID was
verified with `find-tool.py`; guarded units + their `MinimumMode` verified against
`LIMITATIONS.md`. One provenance anomaly is flagged: `detect_rearrangements` exists on
**both** the Analysis server (`ComparativeGenomics.DetectRearrangements`) and the
Chromosome server (`ChromosomeAnalyzer.DetectRearrangements`); the chromosome-scale
synteny task (G10) uses the Chromosome one, matching `bio-chromosome`.

## The 12 golden tasks

The task list (`docs/skills/golden/tasks.md`) fixes each G-task's exact title, the
skills it exercises, and whether it drives a guarded unit. Coordinates are **0-based**
unless a tool doc says otherwise; every tool name + `Method ID` was verified with
`find-tool.py` against `docs/mcp/tools/**` and cited inline.

| # | Title | Skills exercised | Guard / caveat |
|---|---|---|---|
| G1 | Cloning-insert QC → find restriction sites | bio-qc, bio-moldesign | — (cross-domain) |
| G2 | FASTQ quality stats, overlap-confined qualities | bio-qc | ⚠ PARSE-FASTQ-001 (Permissive) |
| G3 | Pairwise + MSA of an ortholog family → consensus | bio-alignment | — |
| G4 | Call + classify + score variants in a CDS | bio-annotation | clinical caveat (pathogenicity) |
| G5 | Design + QC a PCR primer pair | bio-moldesign | — |
| G6 | CRISPR guides for an ORF located by annotation | bio-annotation, bio-moldesign | — (cross-domain) |
| G7 | NJ tree + neutrality (Tajima's D) for a population | bio-phylo-popgen | — |
| G8 | Metagenome: classify → profile → diversity → bin | bio-metagenomics | ⚠ META-BIN-001 (Moderate) |
| G9 | Assemble reads → N50 → k-mer QC | bio-assembly | — |
| G10 | Chromosome centromere + GC-skew origin | bio-chromosome, bio-annotation | — (cross-domain) |
| G11 | Design an MGB / dual-quencher qPCR probe | bio-moldesign | ⚠ PROBE-DESIGN-001 (Moderate) |
| G12 | reads → assemble → annotate ORFs → design primers | bio-assembly, bio-annotation, bio-moldesign | — (cross-domain, 4 skills) |

Every task runs under `bio-rigor` (tool-only computation, provenance, envelope, units +
coordinate base, alpha caveat); `seqeron-discovery` fires only when a tool name is
unknown. Each task also carries **one graded independent cross-check** — a second,
different code path that must reproduce the primary result: G1 `find_restriction_sites`
vs `suffix_tree_find_all`; G3 `multiple_align` consensus vs `compute_consensus`; G4
`classify_variant` types vs `variant_statistics` types; G6/G12 ORF confirmed by
`coding_potential`; G7 `diversity_statistics` vs split-path `nucleotide_diversity` +
`tajimas_d` (which takes **k̂ = π·L**, not per-site π); G8 `taxonomic_profile` diversity
vs `alpha_diversity`; G9/G12 engine `n50` vs `assembly_stats` n50 (and `totalLength ==
Σ|contig|`); G10 `analyze_centromere` vs `find_heterochromatin_regions` overlap. The
guarded tasks pin the exact `MinimumMode` behaviour to assert: G2 encoding must return
`Ambiguous` (all quality chars in the ASCII 64–74 Phred+33/+64 overlap; blocked in Strict
& Moderate); G8 `bin_contigs` throws `SeqeronLimitationException` under Strict (STOP) and
under Moderate returns only **domain-level CheckM approximations** of completeness /
contamination; G11 must **STOP on the MGB-ΔTm demand** (empirical/proprietary, no closed
form) and not fabricate an MGB Tm, delivering only a salt-corrected Tm clearly labelled
*not* MGB-corrected. G10 also notes the guarded `CHROM-CENT-001` (SF1/SF2 assignment) is
**not** invoked — general `analyze_centromere` is unguarded.

## Where this fits

The acceptance instrument for the [[skills-strategy]] and the [[skill-layer]]; enforces
the [[scientific-rigor]] discipline; task style follows the worked workflows in
[[mcp-readme]].
