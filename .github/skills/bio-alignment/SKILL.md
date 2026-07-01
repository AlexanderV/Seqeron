---
name: bio-alignment
description: >-
  Align biological sequences and quantify how similar they are with Seqeron
  (MCP tools OR the C# API). Use for pairwise alignment (global/Needleman-Wunsch,
  local/Smith-Waterman, semi-global/fitting/overlap), multiple sequence alignment
  (MSA) + consensus/conservation, percent identity / similarity / gap %, edit
  (Levenshtein) and Hamming distance, k-mer similarity, longest common region,
  approximate matching, and prompts like "align these sequences", "how similar
  are X and Y", "find the conserved region", "does this query fit in that
  reference". Servers: alignment + core.
allowed-tools: Read, Bash, Grep, Glob
---

# bio-alignment — pairwise/MSA alignment, identity/similarity, distances

Routing + orchestration skill for the **Alignment** and **Core** servers (34 tools).
It picks the right tool for an alignment/similarity question and gives a **dual-mode**
recipe (MCP tool calls **and** the equivalent `Seqeron.Genomics` C# `Method ID`s).

- **Rigor is delegated.** Parse-with-a-tool, envelope, provenance, cross-check, units/0-based
  coordinates, and the alpha / not-for-clinical-use caveat are all owned by **`bio-rigor`** —
  it applies here by default; do not restate its rules.
- **Don't know the tool name?** Use **`seqeron-discovery`**
  (`python3 scripts/skills/find-tool.py <kw> --server alignment|core`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/{alignment,core}/*.md`;
  algorithm invariants in `docs/algorithms/Alignment/*.md`. This skill links, it does not copy.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Align the **whole** of two seqs end-to-end | `global_align` / `SequenceAligner.GlobalAlign` |
| Find the **best-scoring shared region** (ignore flanks) | `local_align` / `SequenceAligner.LocalAlign` |
| **Fit** a short query into a long reference / overlap / containment | `semi_global_align` / `SequenceAligner.SemiGlobalAlign` |
| **Identity % / similarity % / gap %** of an alignment | `alignment_statistics` / `SequenceAligner.CalculateStatistics` |
| Pretty-print an alignment (`\|` `:` legend) | `format_alignment` / `SequenceAligner.FormatAlignment` |
| Align **N** sequences → MSA + consensus | `multiple_align` / `SequenceAligner.MultipleAlign` |
| Consensus from **pre-aligned** equal-length reads | `compute_consensus` / `SequenceAssembler.ComputeConsensus` |
| Quick **% identity** of two equal-length seqs (no gaps) | `sequence_identity` / `SequenceAssembler.CalculateIdentity` |
| **Edit (Levenshtein) distance** (any lengths) | `edit_distance` / `ApproximateMatcher.EditDistance` |
| **Hamming distance** (equal length) | `hamming_distance` / `ApproximateMatcher.HammingDistance` |
| **k-mer Jaccard similarity** (0–1) | `calculate_similarity` / `GenomicAnalyzer.CalculateSimilarity` |
| **Longest common substring / region** (+ positions) | `find_longest_common_region` / `GenomicAnalyzer.FindLongestCommonRegion`; text-only: `suffix_tree_lcs` / `SuffixTree.LongestCommonSubstring` |
| **Approximate match** of a pattern (mismatches / edits) | `find_with_mismatches`, `find_with_edits`, `find_best_match` |

Rule of thumb: **global** for two homologous full-length seqs; **local** to locate a motif/domain
inside dissimilar seqs; **semi-global** for overlap/containment (reads → reference, primer-in-amplicon);
**MSA** for a family; **distances/similarity** for quick pairwise scoring without a full alignment.

## Canonical dual-mode pipelines

Scoring knobs are shared by all pairwise/MSA tools: `match` (+1), `mismatch` (−1), `gapOpen` (−2),
`gapExtend` (−1 per position, linear gap). Change them together, and report them in provenance.

### (a) Global-align two DNA seqs → identity/similarity/gap %
1. **[MCP]** `global_align`(sequence1, sequence2) → `alignedSequence1`, `alignedSequence2`, `score`
2. **[MCP]** `alignment_statistics`(alignedSequence1, alignedSequence2) → `identity`, `similarity`, `gapPercent`
3. **[MCP]** `format_alignment`(…) → human-readable block (optional).
- **[C# API]** `SequenceAligner.GlobalAlign(s1,s2,match,mismatch,gapOpen,gapExtend)` → `SequenceAligner.CalculateStatistics(aln1,aln2)` → `SequenceAligner.FormatAlignment(aln1,aln2,lineWidth)`.
- **Note:** `alignment_statistics` uses the EMBOSS-needle denominator = **full alignment length incl. gap columns**; `sequence_identity` uses a **gapless equal-length** denominator — pick deliberately.
```
Provenance
1) global_align(s1,s2,match=1,mismatch=-1,gapOpen=-2,gapExtend=-1) → score, aln1, aln2
2) alignment_statistics(aln1,aln2) → identity=…% , similarity=…% , gapPercent=…% (denominator=aln length incl. gaps)
Cross-check: identity consistent with matches/alignmentLength from the same call.
Envelope: none guarded. Caveat: alpha — validate before any decision use.
```

### (b) Local-align to locate a conserved region
1. **[MCP]** `local_align`(sequence1, sequence2) → aligned substrings + `startPosition1/2`, `endPosition1/2` (**0-based**).
2. **[MCP]** `alignment_statistics`(aligned pair) → identity % of the shared core.
- **[C# API]** `SequenceAligner.LocalAlign(...)` → `SequenceAligner.CalculateStatistics(...)`.
- **Cross-check** the located region with `find_longest_common_region` / `GenomicAnalyzer.FindLongestCommonRegion` (independent suffix-tree path) — positions should agree.

### (c) Semi-global for overlap / containment / read-fitting
1. **[MCP]** `semi_global_align`(sequence1=query, sequence2=reference) → query gets free end gaps; reference reproduced in full; `startPosition2/endPosition2` = where the query lands (0-based).
- **[C# API]** `SequenceAligner.SemiGlobalAlign(query, reference, …)`.
- Use for: does read A overlap read B? does a primer/insert fit inside an amplicon? (See pipeline e for exact primer location.)

### (d) MSA of N sequences → consensus / conservation
1. **[MCP]** `multiple_align`(sequences=[…]) → `alignedSequences`, `consensus`, `totalScore` (sum-of-pairs).
2. **[MCP]** (optional, independent consensus check) `compute_consensus`(alignedReads=alignedSequences) → majority-vote consensus (ties → `N`).
- **[C# API]** `SequenceAligner.MultipleAlign(seqs,…)` → `SequenceAssembler.ComputeConsensus(alignedSequences)`.
- **Cross-check:** `multiple_align` consensus vs `compute_consensus` on its own rows should match (different code paths).

### (e) Quick similarity / distance without a full alignment
- **k-mer similarity (0–1):** `calculate_similarity`(s1,s2,kmerSize?) / `GenomicAnalyzer.CalculateSimilarity`.
- **Edit distance (any length):** `edit_distance` / `ApproximateMatcher.EditDistance`.
- **Hamming (equal length):** `hamming_distance` / `ApproximateMatcher.HammingDistance`.
- **Longest common substring:** `find_longest_common_region` (DNA, with positions) or `suffix_tree_lcs` (raw text).
- **Approximate motif search:** `find_with_mismatches` (Hamming), `find_with_edits` (Levenshtein), `find_best_match` (single best window).

## End-to-end grounded example (extends `docs/mcp/README.md`)

**Task.** A cloning insert `>seq1 GCGCGAATTCATGGATCCATAT` is thought to be a variant of a reference
amplicon `>ref GCGCGAATTCATGGATCCTTAT`. (1) Confirm the insert fits/aligns to the reference, (2) report
identity/similarity/gap %, (3) render the alignment, (4) corroborate the shared core independently.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `global_align`(sequence1="GCGCGAATTCATGGATCCATAT", sequence2="GCGCGAATTCATGGATCCTTAT")
   → `score`, `alignedSequence1`, `alignedSequence2` (expect a mismatch near the 3′ end, no gaps).
   (`SequenceAligner.GlobalAlign`)
2. `alignment_statistics`(alignedSequence1, alignedSequence2)
   → `identity` ≈ 90.9 %, `gapPercent` = 0 (denominator = 22 aligned columns). (`SequenceAligner.CalculateStatistics`)
3. `format_alignment`(alignedSequence1, alignedSequence2, lineWidth=80) → three-line `|…| …` block. (`SequenceAligner.FormatAlignment`)
4. Cross-check: `find_longest_common_region`(seq1, ref) → common region `GCGCGAATTCATGGATCC` at positions
   (0,0), length 18 — consistent with the aligned identical prefix. (`GenomicAnalyzer.FindLongestCommonRegion`)

Expected-shape output (values illustrative; **compute them with the tools, do not eyeball**):
```
| id   | length | identity_% | similarity_% | gap_% | best_score |
|------|-------:|-----------:|-------------:|------:|-----------:|
| seq1 |     22 |       90.9 |         90.9 |   0.0 |         18 |

Provenance
1) global_align(s1,s2,match=1,mismatch=-1,gapOpen=-2,gapExtend=-1) → score, aln1, aln2
2) alignment_statistics(aln1,aln2) → identity, similarity, gapPercent (denom = alignment length incl. gaps, 0-based cols)
3) format_alignment(aln1,aln2,lineWidth=80) → rendered block
Cross-check: find_longest_common_region(s1,ref) → common core positions agree with the aligned identical prefix.
Envelope: none guarded (Alignment + Core tools within contract).
Caveat: alpha software; not for clinical use — independently validate before relying on any construct decision.
```

## Reference

- **Full domain tool index (all 34, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, run `seqeron-discovery`).
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Tool map (34 tools by sub-task, one-liners + Method ID):** [`reference/tool-map.md`](reference/tool-map.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md`](../../../docs/algorithms/Alignment/Global_Alignment_Needleman_Wunsch.md) ·
  [`Local_Alignment_Smith_Waterman.md`](../../../docs/algorithms/Alignment/Local_Alignment_Smith_Waterman.md) ·
  [`Semi_Global_Alignment.md`](../../../docs/algorithms/Alignment/Semi_Global_Alignment.md) ·
  [`Multiple_Sequence_Alignment.md`](../../../docs/algorithms/Alignment/Multiple_Sequence_Alignment.md) ·
  [`Alignment_Statistics.md`](../../../docs/algorithms/Alignment/Alignment_Statistics.md)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup).
