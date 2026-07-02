# seqeron-mirna — fuller pipelines & parameter guidance

Dual-mode recipes for the microRNA family (annotation server, `MiRnaAnalyzer.*`). Rigor (tool-only,
provenance, cross-check, alpha caveat) is delegated to [`bio-rigor`](../../bio-rigor/SKILL.md).
Schemas: `docs/mcp/tools/annotation/<tool>.md`. Tool index: [`tool-map.md`](tool-map.md).

## Parameter defaults (verify per tool doc)

- `find_mirna_target_sites`: `minScore` = 0.5 (site type + pairing).
- `analyze_target_context`: `contextWindow` = 30 nt each side; `contextScore = 0.5·auContent (+0.3 if mid-transcript)`, capped 1.
- `site_accessibility`: fixed ±50 nt window; `accessibility = max(0, 1 − density·10)`.
- `find_similar_mirnas`: `maxMismatches` = 1 (seed Hamming).
- `find_pre_mirna_hairpins`: `minHairpinLength` = 55, `maxHairpinLength` = 120, `matureLength` = 22; loop 3–25 nt.
- Seed = positions 2–8 (0-based `Substring(1,7)`; `seedStart=1`, `seedEnd=7`). Site positions **0-based inclusive**.
- Energies: **kcal/mol at 37 °C**, Turner 2004 NN. DNA `T` → RNA `U`; inputs case-insensitive.

## 1. Predict miRNA target sites in an mRNA (main entry)

**Goal.** Turn a miRNA + an mRNA into scored canonical target sites.

1. **[MCP]** `create_mirna`(name, sequence) → record with normalised sequence + seed 2–8.
2. **[MCP]** `find_mirna_target_sites`(mRnaSequence, miRna, minScore=0.5) → `sites[{start,end,targetSequence,type,seedMatchLength,score,freeEnergy,alignment}]` (0-based).
3. Cross-check: **[MCP]** `rna_reverse_complement`(miRna.seedSequence) should be present in a hit's `targetSequence` (the seed match is the seed's reverse complement).

- **[C# API]** `MiRnaAnalyzer.CreateMiRna(...)` → `.FindTargetSites(...)` → `.GetReverseComplement(...)`.
- ⚠ **MIRNA-TARGET-001:** `score` is **site-type + pairing only**, not a TargetScan context++ score.

```
Provenance
1) create_mirna(name,sequence) → miRna{sequence,seedSequence,seedStart=1,seedEnd=7}
2) find_mirna_target_sites(mRnaSequence,miRna,minScore=0.5) → sites[{start,end,type,score,freeEnergy}] (0-based)
3) rna_reverse_complement(miRna.seedSequence) ≡ hit.targetSequence seed match (cross-check)
Envelope: MIRNA-TARGET-001 — score is site-type+pairing, NOT context++ (supply 3'UTR set / PCT params or cite TargetScan).
Caveat: alpha — validate before use.
```

## 2. Context + accessibility scoring of a site

**Goal.** Add favourability signals to a found site.

1. **[MCP]** `analyze_target_context`(mRnaSequence, targetStart, targetEnd, contextWindow=30) → `auContent`, `nearStart`, `nearEnd`, `contextScore∈[0,1]`.
2. **[MCP]** `site_accessibility`(mRnaSequence, siteStart, siteEnd) → `accessibility∈[0,1]` (±50 nt window; higher = more open).
3. Cross-check: an AU-rich (`auContent` high) mid-transcript site should also read as more accessible; a strongly structured region depresses `accessibility` independently of AU content.

- **[C# API]** `MiRnaAnalyzer.AnalyzeTargetContext(...)` · `.CalculateSiteAccessibility(...)`.
- ⚠ **MIRNA-TARGET-001:** `contextScore` is a **partial** proxy (AU + position), **not** context++.

```
Provenance
1) analyze_target_context(mRnaSequence,targetStart,targetEnd,contextWindow=30) → auContent,nearStart,nearEnd,contextScore
2) site_accessibility(mRnaSequence,siteStart,siteEnd) → accessibility (±50 nt window)
Envelope: MIRNA-TARGET-001 — contextScore is partial (AU+position), not TargetScan context++.
Caveat: alpha.
```

### STOP rule (MIRNA-TARGET-001)

If the task asks for a **TargetScan-grade context++ score** (or `TA_3UTR`, `PCT`, `SPS`, `Len_ORF`,
`ORF8m`), the library does **not** produce it out-of-the-box (no default transcriptome / TargetScan
sigmoid tables bundled; the guarded branch throws below `Permissive`). Either:
- supply a **3'UTR set** (then `ComputeTa3Utr` derives `TA_3UTR`), the **alignment + tree + published
  per-family sigmoid parameters** for `PCT`, and `SPS` / ORF features yourself; **or**
- **STOP** and report the limitation, naming **TargetScan** as the reference —
  do not present the partial score as context++. Do not raise `LimitationPolicy` to force a number.
Source: [`../../../../docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

## 3. Seed families / similarity / variants

**Goal.** Cluster or expand miRNAs by seed (families predicted to share targets).

1. **[MCP]** `group_by_seed_family`(miRnas) → `families[{seedFamily,members}]` (exact seed equality).
2. **[MCP]** `compare_seed_regions`(miRna1, miRna2) → `matches`, `mismatches` (Hamming + length diff), `isSameFamily`.
3. **[MCP]** `find_similar_mirnas`(query, database, maxMismatches=1) → seed-Hamming neighbours (query excluded by name).
4. **[MCP]** `generate_seed_variants`(seedSequence) → original + `3L` single substitutions (enumerate near-seed patterns).

- **[C# API]** `.GroupBySeedFamily` · `.CompareSeedRegions` · `.FindSimilarMiRnas` · `.GenerateSeedVariants`.
- No guarded unit — exact string logic. `generate_seed_variants` `includeWobble` is reserved (no effect yet).

```
Provenance
1) group_by_seed_family(miRnas) → families[{seedFamily,members}]
2) compare_seed_regions(miRna1,miRna2) → matches,mismatches,isSameFamily
3) find_similar_mirnas(query,database,maxMismatches=1) → seed-Hamming neighbours
Caveat: alpha. (No guarded unit — seed handling is exact.)
```

## 4. miRNA↔target duplex alignment

**Goal.** Inspect the antiparallel duplex of one miRNA against a target window.

1. **[MCP]** `align_mirna_to_target`(miRnaSequence, targetSequence) → `alignmentString` (`|` WC / `:` G-U / ` ` mismatch), `matches`, `mismatches`, `guWobbles`, `gaps` (always 0), `freeEnergy` (Turner 2004).
2. Spot-check a column independently: **[MCP]** `can_pair`(b1,b2) / `is_wobble_pair`(b1,b2).

- Orientation: miRNA position `i` (5'→3') pairs against `target[len−1−i]` (target read 3'→5'). Ungapped.
- **[C# API]** `.AlignMiRnaToTarget(...)` · `.CanPair(...)` · `.IsWobblePair(...)`.

```
Provenance
1) align_mirna_to_target(miRnaSequence,targetSequence) → alignmentString,matches,guWobbles,freeEnergy
2) can_pair / is_wobble_pair on a spot-checked column (independent verification)
Caveat: alpha; freeEnergy is Turner 2004 NN over consecutive paired positions.
```

## 5. Pre-miRNA hairpin detection

**Goal.** Find biogenesis-precursor hairpins (mature + star arms) in a sequence.

1. **[MCP]** `find_pre_mirna_hairpins`(sequence, minHairpinLength=55, maxHairpinLength=120, matureLength=22) → `hairpins[{start,end,sequence,matureSequence,starSequence,structure,freeEnergy}]`.

- **[C# API]** `MiRnaAnalyzer.FindPreMiRnaHairpins(...)`.
- Simplified **consecutive-pairing** model — real pre-miRNAs with internal bulges may not be detected.
- For a **structural** hairpin/stem-loop question (not miRNA biogenesis), route to
  [`seqeron-rna-structure`](../../seqeron-rna-structure/SKILL.md) (`find_stem_loops`, `predict_rna_structure`).

```
Provenance
1) find_pre_mirna_hairpins(sequence,minHairpinLength=55,maxHairpinLength=120,matureLength=22) → hairpins[{start,end,matureSequence,starSequence,structure,freeEnergy}]
Envelope: MIRNA-CLEAVAGE-001 (star-arm span approximate) + MIRNA-PRECURSOR-001 (no read-support score).
Caveat: alpha.
```

### STOP rules (MIRNA-CLEAVAGE-001, MIRNA-PRECURSOR-001)

- **MIRNA-CLEAVAGE-001** (guarded, MinimumMode `Permissive`): the 5p Drosha/Dicer cut reproduces
  miRBase exactly, but the **miRNA\*(3p) / star-arm boundary** is an approximate linear 2-nt-3′-overhang
  offset — miRBase mature boundaries encode dominant sequencing-read cut sites, which folding + a fixed
  overhang does not recover. If exact star-arm coordinates matter, **supply the miRBase mature-3p (MIMAT)
  coordinates or small-RNA-seq read pileups**, or STOP and report the limitation. Do not present the
  offset span as the miRBase annotation, and do not raise `LimitationPolicy` to force output.
- **MIRNA-PRECURSOR-001** (documented-only, no runtime guard — nothing throws): the **read-stacking /
  read-support score** of miRDeep2 is **not implemented** (closed-form log-odds only in a gated 2008
  supplement + GPL miRDeep2 source). `find_pre_mirna_hairpins` returns a *structural* candidate, not a
  read-supported confidence. If asked for read support / miRDeep2-style confidence, report it is out of
  scope and name **miRDeep2 with the caller's small-RNA-seq reads** — do not imply the candidate is
  read-validated. Source: [`../../../../docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

## Scope reminders

- **General RNA folding / MFE / stem-loops / pseudoknots / dot-bracket** → [`seqeron-rna-structure`](../../seqeron-rna-structure/SKILL.md).
  `MiRnaAnalyzer.CanPair` / `.IsWobblePair` here are the **miRNA** copies (distinct Method IDs).
- **mRNA 3'UTR / ORFs / gene structure / splicing** → [`bio-annotation`](../../bio-annotation/SKILL.md).
- **Full context++ / TargetScan tables / small-RNA-seq read support** → out of scope; supply the data
  or use the named external tool (see the STOP rules above).
