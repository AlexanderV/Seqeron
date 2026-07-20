---
name: seqeron-discovery
description: >-
  Find the right Seqeron bioinformatics tool or algorithm for a task WITHOUT
  loading all 427 MCP tool schemas. Use to answer "which Seqeron tool does X?",
  "is there a tool/algorithm for <biology thing>?", "find the tool for
  <task>", or before attaching MCP servers / writing a Seqeron C# API call and
  you need to pick the right tool. Covers melting temperature (Tm), primers,
  CRISPR guides, alignment (NW/SW/MSA), variant calling, k-mers/assembly,
  motifs/ORFs, phylogeny, popgen, metagenomics, restriction, codon
  optimization — across the 11 MCP servers and 247 algorithm docs.
allowed-tools: Bash, Read
---

# Seqeron Discovery

Seqeron exposes **427 MCP tools across 11 servers** plus **247 algorithm docs**.
Attaching every server (or reading every schema) blows up context. This skill
routes a plain-language need to the **right tool name + `Method ID` + doc path**
by grepping the docs — cheap, deterministic, no schemas loaded.

## When to use

- The user asks "is there a Seqeron tool/algorithm for X?" or "which tool does X?".
- **Before** attaching MCP servers or writing a `Seqeron.Genomics` C# call, to pick
  the correct tool without loading 427 schemas into context.
- To locate the full I/O schema (follow the doc path) or the C# entry point
  (`Method ID`, e.g. `PrimerDesigner.CalculateMeltingTemperature`).

## How to run

```bash
python3 scripts/skills/find-tool.py <keywords...> [--server <name>] [--algorithms] [--limit N]
```

- Keywords are **AND**, case-insensitive. Matches are ranked: tool **name** >
  `Method ID` > description. Output is a compact table: tool · server · Method ID · doc.
- `--server <name>` filters to one server dir: `core, sequence, parsers,
  alignment, analysis, annotation, chromosome, metagenomics, moltools,
  phylogenetics, population`.
- `--algorithms` also searches `docs/algorithms/**` (title + Test Unit ID + path).
- `--limit N` caps each section (default 20). `--help` for usage.

Prefer specific keywords (`melting temperature`, `off target`) over broad ones
(`variant` alone returns many hits — add a qualifier or `--server`).

## Examples (verified output)

The **WIKI** column names the LLM-Wiki concept (or `(!) gotcha`) that curates the
science and sharp edges behind each tool — open it for meaning, caveats, and the
validated envelope after you have the tool name. It is read live from the wiki's
own `mcp_tools:` frontmatter, so it stays in sync as ingests add tools.

**Melting temperature → Tm tools in sequence + moltools:**

```
$ python3 scripts/skills/find-tool.py melting temperature
TOOL                              SERVER        METHOD ID                                WIKI (concept / (!) gotcha)          DOC
melting_temperature               Sequence      SequenceStatistics.CalculateMeltingTempe -                                    docs/mcp/tools/sequence/melting_temperature.md
primer_melting_temperature        MolTools      PrimerDesigner.CalculateMeltingTemperatu primer-dimer-thermodynamics-tm       docs/mcp/tools/moltools/primer_melting_temperature.md
primer_melting_temperature_salt   MolTools      PrimerDesigner.CalculateMeltingTemperatu primer-dimer-thermodynamics-tm       docs/mcp/tools/moltools/primer_melting_temperature_salt.md
```

**A tool with a known trap flags its gotcha inline:**

```
$ python3 scripts/skills/find-tool.py predict genes --limit 1
TOOL                              SERVER        METHOD ID                                WIKI (concept / (!) gotcha)          DOC
predict_genes                     Annotation    GenomeAnnotator.PredictGenes             (!) predict-genes-emits-every-orf    docs/mcp/tools/annotation/predict_genes.md
```

**CRISPR → the MolTools guide-design family:**

```
$ python3 scripts/skills/find-tool.py crispr --server moltools --limit 4
TOOL                              SERVER        METHOD ID                                WIKI (concept / (!) gotcha)          DOC
crispr_specificity_score          MolTools      CrisprDesigner.CalculateSpecificityScore -                                    docs/mcp/tools/moltools/crispr_specificity_score.md
crispr_system_info                MolTools      CrisprDesigner.GetSystem                 crispr-guide-rna-design              docs/mcp/tools/moltools/crispr_system_info.md
design_guide_rnas                 MolTools      CrisprDesigner.DesignGuideRnas           crispr-guide-rna-design              docs/mcp/tools/moltools/design_guide_rnas.md
evaluate_guide_rna                MolTools      CrisprDesigner.EvaluateGuideRna          crispr-guide-rna-design              docs/mcp/tools/moltools/evaluate_guide_rna.md
... 2 more tool hits (use --limit to see more).
```

**Alignment, including algorithm docs:**

```
$ python3 scripts/skills/find-tool.py alignment --algorithms --limit 4
TOOL                              SERVER        METHOD ID                                WIKI (concept / (!) gotcha)          DOC
alignment_statistics              Alignment     SequenceAligner.CalculateStatistics      alignment-statistics                 docs/mcp/tools/alignment/alignment_statistics.md
format_alignment                  Alignment     SequenceAligner.FormatAlignment          alignment-statistics                 docs/mcp/tools/alignment/format_alignment.md
...
ALGORITHMS
TITLE                                              TEST UNIT ID          DOC
Global Alignment (Needleman-Wunsch)                ALIGN-GLOBAL-001      docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md
Local Alignment (Smith-Waterman)                   ALIGN-LOCAL-001       docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md
Multiple Sequence Alignment (MSA)                  ALIGN-MULTI-001       docs/algorithms/Alignment/Multiple_Sequence_Alignment.md
```

## Reading the result

- **Doc path** (repo-relative) → open it for the full **Input/Output schema**,
  errors, and examples. Do not guess the schema; read the doc.
- **Method ID** → the direct `Seqeron.Genomics` C# entry point for the API path.
- **Server** → which MCP server to attach when using the MCP path.
- **WIKI** → the LLM-Wiki concept behind the tool (`wiki/concepts/<slug>.md`) for
  the science, parameters, and corner cases; `(!) <slug>` points at a
  `wiki/gotchas/` page recording a known trap — read it before relying on the tool.
- **Test Unit ID** (algorithms) → cross-reference into `docs/Validation` /
  `LIMITATIONS.md` for the validated operating envelope.

Point, don't duplicate: this skill routes; the docs are the source of truth.

## Fallbacks

- **No hits?** Loosen keywords (drop qualifiers), try a synonym, or add
  `--algorithms` — the concept may be an algorithm doc, not a wrapped MCP tool.
- **Optional enrichment:** if `docs/skills/_generated/catalog.json` exists (built
  by `scripts/skills/gen-catalog.py`), `find-tool.py` uses it automatically to
  annotate results; the script works fully without it.
- **Browse a whole server:** `python3 scripts/skills/find-tool.py <broad-term> --server <name> --limit 100`,
  or `ls docs/mcp/tools/<server>/` to list the per-tool docs. Each server also ships a
  README with its full tool table at `src/Seqeron/Mcp/Seqeron.Mcp.<Server>/README.md`
  (Core is `src/SuffixTree/Mcp/SuffixTree.Mcp.Core/README.md`).
- Once you know the tool, hand off to the relevant domain skill (e.g.
  `bio-moldesign`, `bio-alignment`, `bio-annotation`) and to `bio-rigor` for the
  actual computation.
- **To run it, MCP need not be registered anywhere.** Call the tool over the shipped
  MCP server on demand from Python — [`seqeron-python-client`](../seqeron-python-client/SKILL.md) —
  or via the C# API using the tool's `Method ID` — [`seqeron-dev`](../seqeron-dev/SKILL.md).
