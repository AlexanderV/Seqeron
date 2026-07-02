---
name: seqeron-mirna
description: >-
  microRNA (miRNA) analysis with Seqeron (MCP tools OR the C# API). Use for
  canonical miRNA target-site prediction (8mer / 7mer-m8 / 7mer-A1 / 6mer,
  Bartel/TargetScan hierarchy) with AU/positional context and local
  structure accessibility, seed extraction / seed-family grouping / similar-miRNA
  search / single-nucleotide seed variants, miRNA→target antiparallel duplex
  alignment, pre-miRNA hairpin detection (mature/star arms + dot-bracket + Turner
  energy), and RNA base-pairing / G-U wobble tests. Triggers: "find miRNA target
  sites", "predict miRNA targets in this 3'UTR/mRNA", "what's the seed of this
  miRNA", "miRNA seed family", "group these miRNAs by seed", "find similar
  miRNAs", "generate seed variants", "align this miRNA to its target", "find
  pre-miRNA hairpins", "is this a G-U wobble pair in a miRNA duplex", "can these
  two miRNA/target bases pair", "reverse-complement a seed to get its target
  motif", "miRNA target accessibility / context score". A context-free
  base-pairing / G-U question with no miRNA framing → seqeron-rna-structure.
  Server: annotation (all MiRnaAnalyzer.*). THREE guarded units — see the
  Envelope section.
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-mirna — miRNA targets, seeds/families, pre-miRNA hairpins, wobble pairing

Routing + orchestration skill for the microRNA family on the **annotation** server (14 tools, one
backing class: `MiRnaAnalyzer.*`). It picks the right tool for a target-site / seed-family /
pre-miRNA / pairing question and gives a **dual-mode** recipe (MCP tool name **and** the equivalent
`Seqeron.Genomics` C# `Method ID`). This is the **densest guarded family in the repo — three guarded
units** — so the ⚠ Envelope section below is the point of the skill.

- **Rigor is delegated.** Tool-only computation, provenance, cross-check, units/0-based coords, and
  the alpha / not-for-clinical caveat are owned by **[`bio-rigor`](../bio-rigor/SKILL.md)** — applies
  by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server annotation`) — never guess.
- **Point, don't duplicate.** Full I/O schemas: `docs/mcp/tools/annotation/*.md`; algorithm
  invariants: `docs/algorithms/MiRNA/*.md`.

## Scope boundary — this family vs neighbours

- **This skill OWNS** the miRNA workflow: canonical **target-site prediction** + context +
  accessibility, **seed** extraction / families / variants / similarity, **miRNA→target duplex**
  alignment, **pre-miRNA hairpin** (biogenesis precursor) detection, and miRNA base-pairing / wobble.
  **[`bio-annotation`](../bio-annotation/SKILL.md)** name-drops the miRNA row shallowly — route here.
- **Pure RNA secondary-structure folding** (MFE, general stem-loops, pseudoknots, dot-bracket) →
  **[`seqeron-rna-structure`](../seqeron-rna-structure/SKILL.md)**. Disambiguation: a **structural**
  hairpin/stem-loop question ("fold this RNA", "find stem-loops") is rna-structure;
  a **miRNA-biogenesis precursor** hairpin ("find pre-miRNA hairpins", mature/star arms) is
  `find_pre_mirna_hairpins` **here**. `MiRnaAnalyzer.CanPair` / `IsWobblePair` are the **miRNA**
  copies (distinct Method IDs from `RnaSecondaryStructure.CanPair`).
- **mRNA 3'UTR / gene structure / ORFs / splicing** → `bio-annotation`.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Build a miRNA record (T→U, seed 2–8 extracted) | `create_mirna` / `MiRnaAnalyzer.CreateMiRna` |
| Just the seed (positions 2–8) | `mirna_seed_sequence` / `MiRnaAnalyzer.GetSeedSequence` |
| Scan an mRNA for canonical target sites (8mer/7mer/6mer) | `find_mirna_target_sites` / `MiRnaAnalyzer.FindTargetSites` ⚠ |
| AU / positional context around a site | `analyze_target_context` / `MiRnaAnalyzer.AnalyzeTargetContext` ⚠ |
| Local structure accessibility of a site | `site_accessibility` / `MiRnaAnalyzer.CalculateSiteAccessibility` |
| Antiparallel miRNA↔target duplex (matches/wobbles/ΔG) | `align_mirna_to_target` / `MiRnaAnalyzer.AlignMiRnaToTarget` |
| Group miRNAs by identical seed (families) | `group_by_seed_family` / `MiRnaAnalyzer.GroupBySeedFamily` |
| Pairwise seed comparison (same family?) | `compare_seed_regions` / `MiRnaAnalyzer.CompareSeedRegions` |
| Find miRNAs within a seed-Hamming budget | `find_similar_mirnas` / `MiRnaAnalyzer.FindSimilarMiRnas` |
| Enumerate single-substitution seed variants | `generate_seed_variants` / `MiRnaAnalyzer.GenerateSeedVariants` |
| Detect pre-miRNA hairpins (mature/star arms) | `find_pre_mirna_hairpins` / `MiRnaAnalyzer.FindPreMiRnaHairpins` ⚠ |
| Seed's target pattern (reverse complement) | `rna_reverse_complement` / `MiRnaAnalyzer.GetReverseComplement` |
| Can two bases pair (A-U, G-C, G-U)? | `can_pair` / `MiRnaAnalyzer.CanPair` |
| Is a pair a G-U wobble? | `is_wobble_pair` / `MiRnaAnalyzer.IsWobblePair` |

⚠ = touches a guarded unit — read the Envelope section before relying on the output.

## Canonical dual-mode pipelines

### (a) Predict miRNA target sites in an mRNA (main entry)
1. **[MCP]** `create_mirna`(name, sequence) → record (T→U, seed 2–8). (`MiRnaAnalyzer.CreateMiRna`)
2. **[MCP]** `find_mirna_target_sites`(mRnaSequence, miRna, minScore?=0.5) → `sites[{start,end,type,score,freeEnergy,…}]` (0-based). (`.FindTargetSites`)
3. **[MCP]** cross-check a hit's seed pattern independently: `rna_reverse_complement`(seed) should equal the `targetSequence` seed match. (`.GetReverseComplement`)
- ⚠ **MIRNA-TARGET-001:** the returned `score` is **base site-type + pairing only**, *not* the full context++ score. See Envelope.

### (b) Context + accessibility scoring of a found site
1. **[MCP]** `analyze_target_context`(mRnaSequence, targetStart, targetEnd, contextWindow?=30) → `auContent`, `nearStart/nearEnd`, `contextScore∈[0,1]`. (`.AnalyzeTargetContext`)
2. **[MCP]** `site_accessibility`(mRnaSequence, siteStart, siteEnd) → `accessibility∈[0,1]` (±50 nt window). (`.CalculateSiteAccessibility`)
- ⚠ **MIRNA-TARGET-001:** `contextScore` is a **partial** context proxy (AU + position), not TargetScan context++. See Envelope.

### (c) Seed families / similarity / variants
1. **[MCP]** `group_by_seed_family`(miRnas) → `families[{seedFamily,members}]`. (`.GroupBySeedFamily`)
2. **[MCP]** `compare_seed_regions`(miRna1, miRna2) → `matches`, `mismatches`, `isSameFamily`. (`.CompareSeedRegions`)
3. **[MCP]** `find_similar_mirnas`(query, database, maxMismatches?=1) → seed-Hamming neighbours. (`.FindSimilarMiRnas`)
4. **[MCP]** `generate_seed_variants`(seedSequence) → original + all `3L` single substitutions. (`.GenerateSeedVariants`)
- No guarded unit here — seed handling is exact string logic.

### (d) miRNA↔target duplex alignment
1. **[MCP]** `align_mirna_to_target`(miRnaSequence, targetSequence) → `alignmentString` (`|`/`:`/` `), `matches`, `guWobbles`, `freeEnergy` (Turner 2004, kcal/mol). (`.AlignMiRnaToTarget`)
2. **[MCP]** spot-check a column: `can_pair`(b1,b2) / `is_wobble_pair`(b1,b2). (`.CanPair` / `.IsWobblePair`)
- No guarded unit — antiparallel (miRNA `i` vs `target[len-1-i]`), ungapped.

### (e) Pre-miRNA hairpin detection
1. **[MCP]** `find_pre_mirna_hairpins`(sequence, minHairpinLength?=55, maxHairpinLength?=120, matureLength?=22) → `hairpins[{start,end,sequence,matureSequence,starSequence,structure,freeEnergy}]`. (`.FindPreMiRnaHairpins`)
- ⚠ **MIRNA-CLEAVAGE-001** (star-arm span) **and** **MIRNA-PRECURSOR-001** (no read-support score). See Envelope.

## ⚠ Envelope — THREE guarded units (STOP rules)

Authoritative source: [`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md).
Do not raise `LimitationPolicy` to force output (bio-rigor rule 2) — STOP and report + name the alternative.

- **MIRNA-TARGET-001** — MinimumMode **`Permissive`** (guarded; throws in Strict & Moderate).
  **Scope limit:** a full **context++ score is not produced out-of-the-box** — `find_mirna_target_sites`
  gives site type + pairing only; `analyze_target_context` gives a *partial* AU+position proxy.
  `TA_3UTR`, `PCT`, `SPS`, `Len_ORF`, `ORF8m` are all **caller-supplied** (no default transcriptome /
  TargetScan sigmoid tables are bundled). **STOP rule:** if asked for a *TargetScan-grade context++
  score*, either supply a 3'UTR set (then `ComputeTa3Utr` derives `TA_3UTR`) + alignment/tree/sigmoid
  params for `PCT` + `SPS`/ORF features, **or** report the limitation and name TargetScan — do not
  present the partial score as context++.
- **MIRNA-CLEAVAGE-001** — MinimumMode **`Permissive`** (guarded; throws in Strict & Moderate).
  **Scope limit:** the 5p Drosha/Dicer cut reproduces miRBase exactly, but the **miRNA\*(3p)/star-arm
  boundary** returned by `find_pre_mirna_hairpins`/cleavage is an *approximate* linear 2-nt-3′-overhang
  offset, not the miRBase-annotated span (real boundaries encode dominant sequencing-read cut sites).
  **STOP rule:** if the exact star-arm coordinates matter, **supply the miRBase mature-3p (MIMAT)
  coordinates or small-RNA-seq read pileups**, or report the limitation — do not treat the offset span
  as the true annotation.
- **MIRNA-PRECURSOR-001** — **documented-only, no runtime guard** (nothing throws; the shortfall is
  undetectable per call). **Scope limit:** the **read-stacking / read-support score** of miRDeep2 is
  **not implemented** (the closed-form log-odds is only in a gated 2008 supplement + GPL source).
  `find_pre_mirna_hairpins` returns a *structural* candidate (fold + energy), **not** a read-supported
  confidence. **STOP rule:** if asked for read-support / miRDeep2-style confidence, report that it is
  out of scope and name **miRDeep2 with the caller's small-RNA-seq reads** — do not imply the
  structural candidate is read-validated.

## End-to-end grounded example

**Task.** Given let-7a and a candidate 3'UTR, find its target sites, score context + accessibility of
the top hit, and independently corroborate the seed pattern — flagging every guarded layer.

1. `create_mirna`("let-7a","UGAGGUAGUAGGUUGUAUAGUU") → seed `GAGGUAG`. (`MiRnaAnalyzer.CreateMiRna`)
2. `find_mirna_target_sites`(mRnaSequence, miRna, minScore=0.5) → `sites` (0-based; `type`, `score`). (`.FindTargetSites`) — ⚠ score is site-type+pairing (MIRNA-TARGET-001).
3. `rna_reverse_complement`("GAGGUAG") → `CUACCUC`; confirm it appears in the top site's `targetSequence` (independent cross-check). (`.GetReverseComplement`)
4. `analyze_target_context`(mRnaSequence, top.start, top.end, 30) → `auContent`, `contextScore`. (`.AnalyzeTargetContext`) — ⚠ partial context (MIRNA-TARGET-001).
5. `site_accessibility`(mRnaSequence, top.start, top.end) → `accessibility`. (`.CalculateSiteAccessibility`)
```
Provenance
1) create_mirna(name=let-7a,sequence=UGAGGUAGUAGGUUGUAUAGUU) → seed=GAGGUAG
2) find_mirna_target_sites(mRnaSequence,miRna,minScore=0.5) → sites[{start,end,type,score,freeEnergy}] (0-based)
3) rna_reverse_complement(GAGGUAG)=CUACCUC ≡ top-site seed match (cross-check)
4) analyze_target_context(mRnaSequence,top.start,top.end,contextWindow=30) → auContent,contextScore
5) site_accessibility(mRnaSequence,top.start,top.end) → accessibility
Envelope: MIRNA-TARGET-001 — score/context are partial, NOT TargetScan context++ (supply 3'UTR set/PCT params or cite TargetScan).
Caveat: alpha software; not for clinical use — independently validate before relying on any call.
```

## Reference

- **This family's tool map (all 14 — curated index; NOT in domain-map.json, so there is NO
  generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`Seed_Sequence_Analysis.md`](../../../docs/algorithms/MiRNA/Seed_Sequence_Analysis.md) ·
  [`Target_Site_Prediction.md`](../../../docs/algorithms/MiRNA/Target_Site_Prediction.md) ·
  [`MiRNA_Target_Pairing.md`](../../../docs/algorithms/MiRNA/MiRNA_Target_Pairing.md) ·
  [`Pre_miRNA_Detection.md`](../../../docs/algorithms/MiRNA/Pre_miRNA_Detection.md)
- **Operating envelope / guarded units:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md)
  (MIRNA-TARGET-001, MIRNA-CLEAVAGE-001, MIRNA-PRECURSOR-001)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (parent — it name-drops the miRNA row; route here) ·
  [`seqeron-rna-structure`](../seqeron-rna-structure/SKILL.md) (pure RNA folding / stem-loops vs miRNA-biogenesis hairpins — disambiguate).
