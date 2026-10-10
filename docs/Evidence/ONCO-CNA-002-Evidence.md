# Evidence Artifact: ONCO-CNA-002

**Test Unit ID:** ONCO-CNA-002
**Algorithm:** Focal Amplification Detection (GISTIC2 length-based focal/broad split + oncogene mapping)
**Date Collected:** 2026-06-14

---

## Online Sources

### Mermel et al. (2011) — GISTIC2.0 (Genome Biology)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
**Accessed:** 2026-06-14 (fetched via WebFetch of the PMC full-text page)
**Authority rank:** 1 (peer-reviewed paper, Genome Biology 12:R41)

**Key Extracted Points:**

1. **Length-based focal/arm-level split:** GISTIC2.0 separates a copy-number profile into arm-level and focal components by length. The fetched text states focal SCNAs are those with *"length < 98% of a chromosome arm"*, while events *"occupying more than 98% of a chromosome arm"* are classified as arm-level. The procedure removes *"all SCNAs occupying more than 98% of a chromosome arm, leaving only the focal events."*
2. **Length is the natural classifier:** *"This reproducible distribution provides a natural basis for classifying events as 'arm-level' and 'focal' based purely on length."* — i.e. GISTIC2.0 moves away from amplitude-based filtering toward length-based filtering.
3. **Amplitude thresholds (historical/low-level):** A low-amplitude threshold of log2 ratio ±0.1 *"eliminates only low-level artifactual segments"*; older high-amplitude thresholds (log2 0.848 amp / −0.737 del) were used in prior versions to exclude arm-level events. GISTIC2.0 itself relies on length, not a single default amplitude cutoff.

### Broad Institute — GISTIC2 documentation (parameter reference)

**URL:** https://broadinstitute.github.io/gistic2/
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (canonical project / reference-implementation documentation)

**Key Extracted Points:**

1. **`broad_len_cutoff` (default 0.98):** verbatim — *"Threshold used to distinguish broad from focal events, given in units of fraction of chromosome arm."* Default value 0.98. This is the configurable form of the 98% figure in the paper.
2. **`t_amp` (default 0.1):** verbatim — *"Threshold for copy number amplifications. Regions with a copy number gain above this positive value are considered amplified."*
3. **`t_del` (default 0.1):** verbatim — *"Threshold for copy number deletions. Regions with a copy number loss below the negative of this positive value are considered deletions."* (out of scope here; recorded for completeness — ONCO-CNA-003 covers deletions.)

### CNVkit — Calling copy number gains and losses (documentation)

**URL:** https://cnvkit.readthedocs.io/en/stable/calling.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (established reference implementation documentation)

**Key Extracted Points:**

1. **Single-copy-gain log2 value:** verbatim — *"a single-copy gain in a perfectly pure, homogeneous sample has a copy ratio of 3/2. In log2 scale, this is log2(3/2) = 0.585."* Establishes that any gain (CN > 2) has log2 > 0; an amplitude threshold of 0.1 is well below a single-copy gain, so it admits all gains as amplitude-positive.
2. **Focal vs non-focal interpretation:** CNVkit's calling chapter and gainloss discussion describe that non-focal amplified segments show surrounding segments with similar copy ratios and are filtered out as artifacts; focal high-amplitude events are the therapeutically actionable ones. Confirms the length/extent-based notion of "focal".

### NCBI Gene — oncogene cytogenetic locations

**URLs (each fetched via WebFetch on 2026-06-14):**
- ERBB2: https://www.ncbi.nlm.nih.gov/gene/2064 → Location "17q12"
- MYC: https://www.ncbi.nlm.nih.gov/gene/4609 → Location "8q24.21"
- EGFR: https://www.ncbi.nlm.nih.gov/gene/1956 → Location "7p11.2"
- CCND1: https://www.ncbi.nlm.nih.gov/gene/595 → Location "11q13.3"
- MDM2: https://www.ncbi.nlm.nih.gov/gene/4193 → Location "12q15"
- CDK4: https://www.ncbi.nlm.nih.gov/gene/1019 → Location "12q14.1"

**Authority rank:** 5 (curated database, NCBI Gene).

**Key Extracted Points:**

1. **Cytogenetic locations (chromosome + arm) for the registry oncogene panel:** ERBB2 = chr17 q-arm (17q12); MYC = chr8 q-arm (8q24.21); EGFR = chr7 p-arm (7p11.2); CCND1 = chr11 q-arm (11q13.3); MDM2 = chr12 q-arm (12q15); CDK4 = chr12 q-arm (12q14.1). The chromosome+arm prefix (e.g. "17q", "7p", "8q") is what a focal-amplification segment's arm label is matched against for oncogene mapping.

---

### GISTIC2 reference implementation (MATLAB source) — 2026-09 review

**URL:** https://github.com/broadinstitute/gistic2 (cloned; submodule https://github.com/broadinstitute/snputil)
**Accessed:** 2026-09-28
**Authority rank:** 2 (original tool's published source code)

**Key Extracted Points:**

1. `snputil/reconstruct_genomes.m` (focal genome): `cur_segs = intersect(find(cur_Q(:,8) < params.broad_len_cutoff), find(cur_Q(:,12) >= thresh))`, `thresh = t_amp` for amplification fields — the length test is strict `<`; the focal filter combines the arm-fraction cutoff with the `t_amp` amplitude gate (confirms Assumption 1 below is GISTIC2's own rule, not an integration choice).
2. `source/score_genome.m`: `focal_del_segs = find(Qs.del(:,8) < params.broad_len_cutoff & Qs.del(:,12) >= params.t_del)`; `source/Qs.m`: column 8 = "event length as fraction of chromosome arm", column 12 = amplitude.
3. `source/gene_calls.m`: thresholded gain call `A.dat(:,i) > t_amp` (strict). GISTIC2 docs "above" ⇒ strict; `reconstruct_genomes.m` uses `>=` (differs only at exact equality).
4. `source/gp_gistic2_from_seg.m`: `t_amp = numeric_arg(a,'ta',0.1,[0,Inf])`, `broad_len_cutoff = numeric_arg(a,'brlen',0.98,[0 2])`; `numeric_arg` throws on NaN and on `val < lo || val > hi` ⇒ valid ranges t_amp ∈ [0, ∞), broad_len_cutoff ∈ [0, 2] (inclusive).
5. `source/perform_deconstruction.m`, `atomic_zigg_deconstruction.m`, `normalize_by_arm_length.m`: the filter is applied to ziggurat-deconstructed events (amplitude relative to the underlying level, cohort-learned broad levels), arm fraction measured in markers by default — not implemented here (see algorithm doc §5.3).

### GISTIC2 MATLAB source — arm-fraction units and gene mapping (FIN-B24 WP21, F48/F49)

**Source:** `git clone --depth 1 https://github.com/broadinstitute/gistic2` (master `26c590bd3330aafa27618ef3eb4b8d9b301f06b7`), accessed 2026-10-10.
**Authority rank:** 3 (reference implementation source)

1. **Arm fraction in markers (F48).** `source/make_sample_B.m`: `[B(:,6) chrarms] = normalize_by_arm_length(D,B,cyto,1,2,chrarms)`. `source/normalize_by_arm_length.m`: "norm_type: 1 = by number of snps (default); 2 = by length in bp"; case 1 → `lengths = Q(:,3) - Q(:,2)+1; fract = lengths./norms'` with `norms = armlengths_by_snp` = `band.snp_length = length(find_snps(D,band.chrn,band.start,band.end,0))`, `band.start = cyto(idx(1)).start+1`, `band.end = cyto(idx(end)).end`; `find_snps.m` counts markers with `cpos>=st & cpos<=en`. `ref_length = 2`: a centromere-spanning segment gets p-fraction + q-fraction. Event fraction = `sum(Bt(:,6))` (`perform_deconstruction.m`).
2. **Gene mapping by locus overlap (F49).** `source/genetables.m` calls `genes_at(rg, chr, start, end, 1, partial_hits(k))`, `partial_hits` default `ones(...)`. `source/genes_at.m` (partial_hits = 1): `in_reg = find((rg(in_chr).start <= pos_end) & (rg(in_chr).end >= pos_start))` — closed-interval overlap, reference-gene order; `pos_end < pos_start` → `error`. (`partial_hits = 0` = containment; closest-gene fallback printed as `[GENE]` only when no gene overlaps.)
3. **GRCh38 gene loci.** NCBI eutils / Ensembl REST / ncbi.nlm.nih.gov were unreachable from the sandbox (proxy 403 / DNS) on 2026-10-10, so the default panel uses GISTIC2's own hg38 gene source `refgenes/Gencode.v22.170324/gencode_genes.tsv` (GENCODE v22, 1-based closed): ERBB2 chr17:39687914-39730426; MYC chr8:127735434-127741434; EGFR chr7:55019021-55256620; CCND1 chr11:69641087-69654474; MDM2 chr12:68808172-68850686; CDK4 chr12:57747727-57756013; TP53 chr17:7661779-7687550; RB1 chr13:48303751-48481986; CDKN2A chr9:21967753-21995301; PTEN chr10:87863113-87971930; BRCA1 chr17:43044295-43125483; BRCA2 chr13:32315474-32400266. Cross-check vs GISTIC2's older `refgenes/hg38.UCSC.add_mir.160920/refGene.txt` (RefSeq 2016, union of transcripts, txStart+1..txEnd): same chromosomes and overlapping spans (e.g. BRCA1 identical 43044295-43125483; TP53 7668402-7687550; ERBB2 39688084-39728662) — RefSeq spans are equal or narrower.

**Reference runs (GNU Octave 8, original .m files unmodified, `addpath source`):**

| Case (`normalize_by_arm_length(D,Q,cyto,1,2)`; chr1 p = 10, q = 40 markers) | GISTIC2 fract (markers) | bp (norm_type 2) |
|---|---|---|
| 1p markers 1..3 | 0.29999999999999999 | 0.20000001 |
| whole 1p (10/10) | 1 | 0.90000001 |
| 1q 39/40 | 0.97499999999999998 | 0.950000005 |
| 1q 40/40 | 1 | 0.975000005 |
| spans centromere (6/10 p + 20/40 q) | 1.1 | — |

| Region (`genes_at(rg,chr,st,en)` on the 12-gene panel) | Genes |
|---|---|
| 17:39730426-39800000 | ERBB2 |
| 17:39730427-39800000 | (none) |
| 17:39600000-39687914 | ERBB2 |
| 17:39600000-39687913 | (none) |
| 12:57756013-68808172 | MDM2 CDK4 |
| 12:57756014-68808171 | (none) |
| 17:7000000-50000000 | ERBB2 TP53 BRCA1 |
| 13:32315474-32315474 | BRCA2 |
| 8:127736000-127740000 | MYC |
| 9:1-21967752 | (none) |
| 9:1-21967753 | CDKN2A |
| 17:60000000-61000000 | (none) |

All values are locked in `OncologyAnalyzer_DetectFocalAmplifications_Tests` / `OncologyAnalyzer_DetectHomozygousDeletions_Tests` (F48/F49 regions).

### GISTIC2 ziggurat deconstruction — source and Octave runs (FIN-B24 WP22–WP24, F50–F52)

**Source:** broadinstitute/gistic2 master `26c590bd` (`source/`) + submodule broadinstitute/snputil `cf3172b8` (`git submodule update --init snputil`), accessed 2026-10-10. Authority rank 3 (reference implementation).

**Octave harness.** GNU Octave 8.4.0, `addpath source; addpath snputil`, original `.m` files unmodified. One shim, documented: `bar3.m` (no-op) — Octave 8 has no `bar3`, which is called only in the plotting branch of `generate_2d_hists` (`do_plot = 1` in `perform_deconstruction`); figures are created invisible (`set(0,'defaultfigurevisible','off')`). `modi.m` (progress print) comes from snputil (no shim). Test D/cyto structures are built by a driver (`mkD.m`): chromosome c with P p-arm and Q q-arm markers has markers at k·10⁵ bp (k = 1..P+Q), p band [0, P·10⁵ + 5·10⁴], q band [P·10⁵ + 5·10⁴, (P+Q+1)·10⁵] — so `find_snps` puts exactly the first P markers on p (`band.start = cyto.start + 1`, `find_snps` closed, returns global indices `in_chr(snps)`). Output printed with `%.17g` and parsed exactly in the tests.

**F50 — per-sample building blocks.** `make_sample_B.m`: breakpoints `find(diff(D.dat(:,idx)) ~= 0)` ∪ `chrnEnd`, rows `[chrn st en D.dat(bpt) sample]`, column 6 = `normalize_by_arm_length(D,B,cyto,1,2,chrarms)` (marker units; centromere-spanning rows: `Q(:,3) = armstart_by_snp − 1` on p, `QQ(:,2) = armstart_by_snp` on q, fractions summed — `armstart_by_snp` is a global marker index because `find_snps` returns `in_chr(snps)`, verified on chr 2 with P ≠ Q). `deconstruct_sample.m` → `deconstruct_chr.m` (row with `en == chr_bpt`; last row → all p; **first row (and more rows) → all rows q** — reference quirk kept; else split; levels subtracted) → `prepare_B.m` (`B.*(B>0)`, `−(B.*(B<0))`, `merge_adj_segs`) → `atomic_zigg_deconstruction.m` (first maximum, higher neighbour with left on ties `>=`, event `[chrn st en amp sample cn_st cn_en fract]` when the step is > 0, merge equal neighbours, residual level from 0) → `add_broad_levels_to_zigg.m` (cn_st, cn_en += level; deletion levels stay in the negated space + level, as in the reference).

| Case (layout chr1 10p+10q, chr2 8p+12q) | B rows (chrn st en amp sample fract) | ZA / ZD (chrn st en amp sample cn_st cn_en fract) |
|---|---|---|
| A: 1p 0.5\|1.5\|0.5, 1q 0, chr2 0.2; levels 0 | 1 1 3 0.5 1 0.3; 1 4 6 1.5 1 0.3; 1 7 10 0.5 1 0.4; 1 11 20 0 1 1; 2 21 40 0.2 1 2 | ZA: 1 4 6 1 1 0.5 1.5 0.3; 1 1 10 0.5 1 0 0.5 1; 2 21 40 0.2 1 0 0.2 2 — ZD: none |
| A with levels (1p 0.5, 1q 0, chr2 0.2), bpt (10, 40) | — | ZA: 1 4 6 1 1 0.5 1.5 0.3 — ZD: none |
| B: 1q 0.6\|−0.8\|0.6, chr2 −0.3; levels 0 | 1 1 10 0 1 1; 1 11 14 0.6 1 0.4; 1 15 16 −0.8 1 0.2; 1 17 20 0.6 1 0.4; 2 21 40 −0.3 1 2 | ZA: 1 11 14 0.6 …0.4 (×2, 11–14, 17–20) — ZD: 1 15 16 0.8 1 0 0.8 0.2; 2 21 40 0.3 1 0 0.3 2 |
| B with levels (0, 0.6, −0.3, −0.3), bpt (10, 40) (first-row quirk) | — | ZA: none — ZD: 1 15 16 1.3999999999999999 1 0.6 2 0.2; 1 1 10 0.6 1 0.6 1.2 1 |
| C: chr2 0\|1.0\|0 with the gain on local markers 6..12, sample 2 | …; 2 26 32 1 2 0.70833333333333326 (3/8 + 4/12); 2 33 40 0 2 0.66666666666666663 | ZA: 1 1 20 0.3 2 0 0.3 2; 2 26 32 1 2 0 1 0.70833333333333326 |
| D: chr2 0.4\|1.2\|−0.5, levels (0.1, 0.1, 0.7, 0.4), bpt (20, 24) | 2 21 24 0.4 1 0.5; 2 25 29 1.2 1 0.58333333333333337; 2 30 40 −0.5 1 0.91666666666666663 | ZA: 2 25 29 0.79999999999999993 1 0.4 1.2 …; ZD: 2 30 40 0.9 1 0.4 1.3 … |

(Abbreviated decimals in the table; the tests hold the full `%.17g` strings.) All rows equal bit-for-bit in `OncologyAnalyzer_ZigguratDeconstruction_Tests` (F50 region).

**F51 — likelihood table and arm broad level.** `generate_2d_hists.m` (called as `generate_2d_hists(QA,QD,[],[],.01,1)`): `xamp = -2:.08:2`, `ylen = 0:0.04:2` (Octave ranges = start + i·step for all 51 elements, checked element-wise), `hd = hist2d(Q(:,4),Q(:,8),xamp,ylen)` (`hist2d.m`: 1 if ≤ first edge, 51 if ≥ last, else `r(i) ≤ x < r(i+1)`; `crosstab_full.m` counts), `hd1 = hd + pseudocount/100·sum(sum(hd))`, `log_hd = log(hd1/sum(sum(hd1)))` (note: the `isempty(ylen)` test guarding the pseudocount default never fires, so `.01` is used). `score_ziggs_by_table.m`: `xidx = find(xamp < amp,1,'last')`, `yidx = find(ylen < fract,1,'last')` (1 when empty). `find_max_broad_level_by_table.m`: `unique(B(:,4))`; > 1 level → `ziggurat_on_extremes.m` (above min / below max, merged, atomic, levels shifted back, scored) then per level `iterative_ziggurat.m` (`max(cn_st, level)` / `min(…)`, amplitude recomputed, wrong-sign events dropped, changed rows rescored) + broad row `[chrn min(st) max(en) level sample 0 level arm_fract]`, score = `sum(Q(:,9)) + broad score`, broad row appended when level ≠ 0, first maximum kept; 1 level → one broad row (score column left 0 — the 8-column `max_Q` is padded when column 10 is set), score by table; 0 levels → level 0, score 0.

Octave run (3-sample cohort on the same layout, copy number − 2 units: S1 1p 0.5\|1.5\|0.5 + chr2 0.2; S2 1q 0.6\|−0.8\|0.6 + chr2 −0.3\|0.9\|−0.3 spanning the centromere; S3 chr1 −0.4, 2p 0.25\|2.4\|0.25): initial events give 11 occupied (amplitude, length) bins, log_hd = −2.714898452106028 (count 1) / −2.0223507320495973 (count 2) / background −9.4415314548696934. All 32 `find_max_broad_level_by_table` calls over every breakpoint and arm part, e.g. S1 1p rows 1..3 → level 0.5, score −12.15642990697572, events +1.0 (markers 4–6, score −2.714898452106028) and broad 0.5 (fract 1, score −9.4415314548696934); S1 1p rows 2..4 → level 0, score −21.597961361845414. All values bit-identical in C# (`Generate2dHistogram_CohortInitialEvents_MatchOctaveLogTable`, `FindMaxBroadLevelByTable_EveryArmSplitOfCohort_MatchesOctave`, which embed the full Octave output).

---

## Documented Corner Cases and Failure Modes

### From Mermel et al. (2011) / GISTIC2 docs

1. **Whole-arm event excluded as focal:** a segment occupying ≥ 98% of its arm is arm-level/broad, not focal — it must NOT be reported as a focal amplification even if highly amplified.
2. **Boundary at exactly 98%:** the paper's wording is "more than 98%" → arm-level, "less than 98%" → focal. The 98% point itself is the cutoff; the doc parameter is "fraction of chromosome arm" = 0.98. We treat a fraction strictly less than the cutoff as focal (length/arm < 0.98 → focal).

### From CNVkit docs

1. **Amplitude below threshold:** a segment whose log2 gain does not exceed `t_amp` (0.1) is not amplified at all and is excluded regardless of length.

---

## Test Datasets

### Dataset: Synthetic GISTIC2-rule worked segments

**Source:** Derivation from Mermel et al. (2011) 98%-of-arm rule and GISTIC2 `t_amp` = 0.1; oncogene arms from NCBI Gene.

| Segment | Arm | ArmLength (bp) | Start–End (bp) | SegLen | SegLen/Arm | log2 | Amplified? (>0.1) | Focal? (<0.98) | Focal amplification? |
|---------|-----|---------------|----------------|--------|-----------|------|-------------------|----------------|----------------------|
| A (ERBB2) | 17q | 1,000,000 | 100,000–600,000 | 500,000 | 0.50 | 1.0 | yes | yes | **yes** |
| B (whole arm) | 8q | 1,000,000 | 0–990,000 | 990,000 | 0.99 | 1.5 | yes | no | no (arm-level) |
| C (low amp) | 7p | 1,000,000 | 0–300,000 | 300,000 | 0.30 | 0.05 | no | yes | no (not amplified) |
| D (boundary) | 11q | 1,000,000 | 0–980,000 | 980,000 | 0.98 | 1.0 | yes | no (= cutoff) | no (not < 0.98) |

---

## Assumptions

1. **ASSUMPTION: amplitude test for "amplified".** GISTIC2.0 itself classifies focal-vs-broad purely by length; "amplification" (gain direction with positive amplitude) is taken from the GISTIC2 `t_amp` = 0.1 parameter (gain above +0.1). This is source-backed (GISTIC2 docs) but combining the length rule (paper) with the `t_amp` amplitude rule (docs) into a single `DetectFocalAmplifications` predicate is the integration choice of this unit; it is documented rather than invented.
2. **ASSUMPTION: arm fraction provided as input.** GISTIC2 derives chromosome-arm boundaries from the genome assembly/cytoband file. This unit does not bundle a cytoband table; the caller supplies each segment's arm label and the arm's length, and the algorithm computes the segment-length / arm-length fraction. The 0.98 cutoff and the amplitude rule are unchanged.

---

## Recommendations for Test Coverage

1. **MUST Test:** A high-amplitude segment spanning < 98% of its arm is reported as a focal amplification; a ≥ 98%-of-arm segment is not. — Evidence: Mermel et al. (2011) 98%-of-arm rule; GISTIC2 `broad_len_cutoff` = 0.98.
2. **MUST Test:** A gain whose log2 does not exceed `t_amp` (0.1) is not reported even if focal in length. — Evidence: GISTIC2 `t_amp` = 0.1.
3. **MUST Test:** Boundary at exactly 0.98 fraction is NOT focal (treated as arm-level). — Evidence: paper "more than 98% ⇒ arm-level".
4. **MUST Test:** `IdentifyAmplifiedOncogenes` maps a focal amplification on arm "17q" to ERBB2, "8q" to MYC, "7p" to EGFR, "11q" to CCND1, "12q" to MDM2 and CDK4. — Evidence: NCBI Gene locations.
5. **SHOULD Test:** A non-amplified (neutral/loss) segment is never an oncogene amplification. — Rationale: only focal amplifications feed the mapper.
6. **COULD Test:** Null / empty input handling. — Rationale: documented failure modes.

---

## References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. (2011). GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12:R41. https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
2. Broad Institute. GISTIC2 documentation (parameter reference: `broad_len_cutoff`, `t_amp`, `t_del`). https://broadinstitute.github.io/gistic2/
3. Talevich E, Shain AH, Botton T, Bastian BC. CNVkit — Calling copy number gains and losses. https://cnvkit.readthedocs.io/en/stable/calling.html
4. NCBI Gene database: ERBB2 (2064), MYC (4609), EGFR (1956), CCND1 (595), MDM2 (4193), CDK4 (1019). https://www.ncbi.nlm.nih.gov/gene/

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: Review B24 — GISTIC2 MATLAB source cross-check; threshold range validation (t_amp ∈ [0,∞), broad_len_cutoff ∈ [0,2], NaN rejected); ziggurat-deconstruction limitation documented.
- **2026-10-10**: FIN-B24 WP21 — F48 marker-unit arm fraction (GISTIC2 `normalize_by_arm_length` norm_type 1), F49 locus-overlap gene mapping (GISTIC2 `genes_at` partial_hits 1) with GRCh38 default panels; Octave reference runs recorded.
- **2026-10-10**: FIN-B24 WP22 — F50 GISTIC2 ziggurat per-sample building blocks (`make_sample_B`, `deconstruct_sample`) ported (internal); Octave runs of the original `.m` files recorded.
- **2026-10-10**: FIN-B24 WP23 — F51 `generate_2d_hists` + `find_max_broad_level_by_table` (with `ziggurat_on_extremes`, `iterative_ziggurat`, `score_ziggs_by_table`, `hist2d`) ported (internal); Octave cohort run recorded.
