# bio-annotation — pipelines, parameters, coordinates, envelope

Fuller recipes and the gotchas that the router SKILL.md keeps short. **Point, don't duplicate**:
per-tool I/O schemas live in `docs/mcp/tools/{annotation,analysis}/<tool>.md`; algorithm invariants
in `docs/algorithms/*`. Rigor (parse-with-a-tool, provenance, cross-check, disclaimer) is owned by
[`bio-rigor`](../../bio-rigor/SKILL.md).

## Coordinate & unit conventions (verify before reporting — bio-rigor rule 5)

- **Variant positions are 0-based on the reference.** `call_variants` / `annotate_variants` return
  `position` (0-based ref) plus `queryPosition`. `predict_variant_effect` takes `variantPosition` =
  **0-based within the CDS**. State the base in every output.
- **ORFs:** `find_orfs` returns `start` (0-based, forward strand), exclusive `end` **including the
  stop codon**, `frame` 1..3, `isReverseComplement`, and `proteinSequence` ending in `*`.
- **Interval tools** (repeats, LCRs, CpG islands, motif hits, splice sites) return 0-based
  half-open-style coordinates unless a tool doc states otherwise — confirm per tool.
- **Units:** GC%/methylation as fractions or %, energies in kcal/mol (Turner 2004), Ti/Tv as a ratio,
  TPM/FPKM as normalized expression, ANI as % identity.
- **Alphabet:** validate DNA/RNA/protein input before deriving anything (bio-rigor rule 1). RNA tools
  (`create_mirna`, RNA structure) T→U normalize; don't feed protein into DNA tools.

## Envelope — guarded / documented-limited units (STOP rule)

Before using these, confirm the task is inside the validated envelope. A guarded call throws
`SeqeronLimitationException` below its **MinimumMode** (`Strict` < `Moderate` < `Permissive`, default
`Moderate`). **STOP and report the limitation — do not raise the mode to force output** (bio-rigor
rule 2). Single source of truth: [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

| Unit | Tools affected | Guard / limit |
|---|---|---|
| **DISORDER-REGION-001** | `predict_disorder` (uncalibrated confidence). MinimumMode **Permissive** → throws in Strict & Moderate on the uncalibrated-confidence branch. | Use the validated-boundary path for reportable region calls; report the disorder *propensity* only when the region-confidence branch is blocked. |
| **MIRNA-TARGET-001** | `find_mirna_target_sites`, `analyze_target_context` (partial "context++" context). MinimumMode **Permissive**. | Seed/duplex results are fine; the full context-score branch is best-effort — say so, don't imply a validated context++ score. |
| **MIRNA-CLEAVAGE-001** | pre-miRNA / cleavage-span estimates (approximate 3p/star span). MinimumMode **Permissive**. | Report spans as approximate. |
| **RNA-STRUCT-001** | `predict_rna_structure`, `minimum_free_energy` and the Turner-2004 energy terms — **documented-only** (no runtime guard). | The greedy predictor is not a full partition-function MFE; state the method and its scope, cross-check against `find_stem_loops` / dot-bracket validation. |

**C# note:** in tests, guarded units need the **Permissive bootstrap**
(`LimitationPolicy` → `Permissive`) — see [`bio-rigor/reference/envelope.md`](../../bio-rigor/reference/envelope.md).
Never flip the mode in production just to get a number.

## Pipeline details

### Structural annotation (family 1)
- `find_orfs` params: `minLength` (aa, default 100 — lower to ~30–50 for short bacterial CDS),
  `searchBothStrands` (default true), `requireStartCodon` (default true; ATG/GTG/TTG). Below-minimum
  ORFs are dropped → empty list, not an error.
- **Order matters:** ORFs → genes → regulatory motifs. Feed the ORF/gene coordinates when scanning
  for `find_promoter_motifs` (bacterial −10/−35) or `find_ribosome_binding_sites` (Shine–Dalgarno is
  looked for **upstream of forward-strand ORFs**).
- Serialize with `to_gff3`; round-trip with `parse_gff3` to confirm the model.
- **Cross-check:** `longest_orfs_per_frame` vs the per-frame max of `find_orfs`; `coding_potential`
  (CPAT log-likelihood) to distinguish a real CDS from a spurious ORF.

### Variant workflow (family 2)
- `call_variants(reference, query)` globally aligns then reports diffs: mismatches → `SNP`, ref-gap →
  `Insertion`, query-gap → `Deletion`. For inputs **already aligned**, use
  `call_variants_from_alignment` (no re-alignment). For a quick gapless SNP scan of equal-length
  seqs, `find_snps_direct`.
- `annotate_variants(reference, query, isCodingSequence=true)` is the one-shot: call + effect +
  mutationType. When `isCodingSequence=false` (default) it calls without protein effects.
- **Effect vs classification vs pathogenicity are distinct steps:** `classify_variant`
  (type: SNV/Ins/Del/MNV/Indel/Complex) → `predict_variant_effect` (protein consequence in a CDS) →
  `predict_pathogenicity` (ACMG-like, combines annotation + frequency + conservation + ClinVar).
  **`predict_pathogenicity` is decision-relevant** → always attach the alpha / not-for-clinical
  caveat and the "independently validate" note (bio-rigor rule 6).
- Summaries: `variant_statistics` (totals, Ti/Tv, density) and `titv_ratio`. Export via
  `variants_to_vcf` / `format_vcf_info`; import via `parse_vcf_variant`.
- **Cross-check:** `find_snps` + `find_indels` should reconcile with `call_variants`;
  `classify_mutation` (Ti/Tv/other) should agree with the type breakdown in `variant_statistics`.

### Structural variants (family 3)
- Read-signal SV calling: `find_discordant_pairs` → `cluster_discordant_pairs`; and/or
  `find_split_reads` → `cluster_split_reads` → `assemble_breakpoint_sequence` +
  `find_microhomology`. Then `genotype_sv`, `filter_svs`, `merge_overlapping_svs`, `annotate_svs`.
- CNV path: `segment_copy_number` (probe log-ratios → plateaus) → `identify_cnvs` (segments →
  del/dup SVs).
- **Cross-check:** breakpoints from split-read clustering vs discordant-pair clustering should
  co-localize.

### Motif discovery → scan (family 4)
- Discovery: `discover_motifs` (de novo overrepresented k-mers) or `find_known_motifs` (a supplied
  set). Build a model with `create_pwm` from aligned instances, then `scan_with_pwm(sequence, pwm,
  threshold)` for genome-wide hits.
- Exact/degenerate confirmation: `find_exact_motif` / `find_motif` (suffix-tree) or
  `find_degenerate_motif` (IUPAC). Protein motifs: `find_protein_motifs` / `find_motif_by_prosite` /
  `find_protein_domains`; convert PROSITE→regex with `prosite_to_regex`.
- **Cross-check (bio-rigor rule 4):** re-scan the reverse complement; verify an exact hit position
  via `kmer_positions`. **Mask first** if the region is low-complexity (family 6) — otherwise motif
  hits in a homopolymer/STR are spurious.

### Repeat / low-complexity masking (family 5–6)
- Detect: `find_tandem_repeats` / `find_microsatellites` / `find_direct_repeats` /
  `find_inverted_repeats` / `find_palindromes`; aggregate with `tandem_repeat_summary`.
- Complexity: `find_low_complexity_regions` (entropy), `dust_score` (DUST), `windowed_complexity`,
  `compression_ratio`. Mask with `mask_low_complexity` and feed the masked sequence into motif /
  alignment / variant steps.
- Protein LCRs: `find_protein_low_complexity_regions` / `predict_low_complexity_seg` (SEG).
- **Cross-check:** a window flagged by `find_low_complexity_regions` should exceed the `dust_score`
  threshold.

### k-mer / composition (family 7)
- Counting: `count_kmers` / `count_kmers_both_strands`; normalize with `kmer_frequencies`; profile
  with `kmer_spectrum` (coverage/error structure) and `analyze_kmers`.
- Enrichment/positions: `most_frequent_kmers`, `kmers_with_min_count`, `find_clumps` (windowed —
  e.g. DnaA boxes), `kmer_positions`, `unique_kmers`.
- Distance: `kmer_distance` (alignment-free composition distance between two sequences).
- GC/skew: `analyze_gc_content`, `gc_content_profile`, `gc_skew` / `cumulative_gc_skew` /
  `windowed_gc_skew`, `at_skew`; `predict_replication_origin` from the cumulative-skew minimum.
- **Cross-check:** the skew-minimum origin should co-localize with a `find_clumps` DnaA-box signal.

### Splicing / epigenetics / transcriptome / miRNA / RNA-struct / protein / comparative
See the family tables in [`tool-map.md`](tool-map.md). Typical chains:
- **Splicing:** `find_donor_sites` + `find_acceptor_sites` → `predict_introns` →
  `predict_gene_structure`; score sites with `maxent_score`.
- **Methylation:** `simulate_bisulfite_conversion` → `methylation_from_bisulfite` →
  `methylation_profile`; regions via `find_dmrs`; islands via `find_cpg_islands`.
- **Transcriptome:** `calculate_tpm` → (`log2_transform` / `quantile_normalize`) →
  `differential_expression` → `over_representation_analysis` / `enrichment_score`; `perform_pca` /
  `cluster_genes_by_expression` for structure.
- **miRNA (⚠):** `create_mirna` → `mirna_seed_sequence` → `find_mirna_target_sites` /
  `align_mirna_to_target` — mind MIRNA-TARGET-001 / MIRNA-CLEAVAGE-001.
- **RNA structure (⚠):** `predict_rna_structure` / `minimum_free_energy`; validate with
  `validate_dot_bracket`; mind RNA-STRUCT-001 (greedy, not full-MFE).
- **Protein (⚠ disorder):** `predict_disorder` (see DISORDER-REGION-001), `predict_signal_peptide`,
  `predict_transmembrane_helices`, `find_protein_domains`.
- **Comparative:** `calculate_ani` / `find_orthologs` / `find_reciprocal_best_hits` /
  `find_syntenic_blocks` / `generate_dot_plot`.

## Provenance template (bio-rigor rule 3)

```
Provenance
1) <tool>(<params>) → <key outputs>
2) <tool>(<params>) → <key outputs>
...
Coordinates: 0-based reference (or CDS) positions. Units: <GC% / bp / kcal·mol⁻¹ / Ti:Tv / TPM / ANI%>.
Cross-check: <second independent tool or invariant> → <agreement>.
Envelope: <none guarded | DISORDER-REGION-001 | MIRNA-TARGET-001 | MIRNA-CLEAVAGE-001 | RNA-STRUCT-001>.
Caveat (if decision-relevant, esp. predict_pathogenicity): alpha software; not for clinical use — validate independently.
```
