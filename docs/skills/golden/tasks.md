# Golden tasks — Seqeron skills regression set

12 tasks. See [`README.md`](README.md) for how to use this and the pass criteria. Every tool name +
`Method ID` below was verified with `scripts/skills/find-tool.py` against `docs/mcp/tools/**`; the
per-tool doc path is cited inline. Coordinates are **0-based** unless a tool doc says otherwise.

Rigor throughout is owned by [`bio-rigor`](../../../.claude/skills/bio-rigor/SKILL.md): tool-only
computation, provenance block, envelope respect, cross-checks, units/coordinates, alpha caveat.

## Summary

| # | Title | Skills exercised | Guarded? |
|---|---|---|:--:|
| G1 | Cloning-insert QC → find restriction sites | bio-rigor, bio-qc, bio-moldesign | — |
| G2 | FASTQ quality stats with overlap-confined qualities | bio-rigor, bio-qc | ⚠ PARSE-FASTQ-001 |
| G3 | Pairwise + MSA of an ortholog family → consensus | bio-rigor, bio-alignment | — |
| G4 | Call + classify + score variants in a CDS | bio-rigor, bio-annotation | clinical caveat (pathogenicity) |
| G5 | Design + QC a PCR primer pair | bio-rigor, bio-moldesign | — |
| G6 | CRISPR guides for an ORF located by annotation | bio-rigor, bio-annotation, bio-moldesign | — |
| G7 | NJ tree + neutrality (Tajima's D) for a population | bio-rigor, bio-phylo-popgen | — |
| G8 | Metagenome: classify → profile → diversity → bin | bio-rigor, bio-metagenomics | ⚠ META-BIN-001 |
| G9 | Assemble reads → N50 → k-mer QC | bio-rigor, bio-assembly | — |
| G10 | Chromosome centromere + GC-skew origin | bio-rigor, bio-chromosome, bio-annotation | — |
| G11 | Design an MGB / dual-quencher qPCR probe | bio-rigor, bio-moldesign | ⚠ PROBE-DESIGN-001 |
| G12 | Full chain: reads → assemble → annotate ORFs → primers | bio-rigor, bio-assembly, bio-annotation, bio-moldesign | — |

Fst / HWE are covered in G7's cross-checks; popgen Fst has its own note there. Cross-domain chains: **G1, G6, G10, G12**. Guarded-unit STOP exercises: **G2, G8, G11**.

---

## G1 — Cloning-insert QC → restriction sites *(cross-domain)*

**Task.**
> Use tools only; no manual parsing or calculations. I have a cloning insert in FASTA below. Read the
> sequence with a tool (do not interpret FASTA by hand), confirm it is valid DNA, report its length and
> GC% (2 decimals), and tell me the 0-based positions of every EcoRI (GAATTC) and BamHI (GGATCC) site.
> Return one Markdown table.
>
> ```
> >insert1
> GCGCGAATTCATGGATCCATATGAATTCGGG
> ```

**Expected skill(s).** `bio-qc` (parse/validate/GC%) + `bio-moldesign` (restriction sites), under
`bio-rigor`. (`seqeron-discovery` only if a tool name is unknown.)

**Expected pipeline** (MCP → C#):
1. `fasta_parse`(content) → sequence — `FastaParser.Parse` · [doc](../../../docs/mcp/tools/parsers/fasta_parse.md)
2. `dna_validate`(sequence) → valid, length — `DnaSequence.TryCreate` · [doc](../../../docs/mcp/tools/sequence/dna_validate.md)
3. `gc_content`(sequence) → gcContent, gcCount, totalCount — `SequenceExtensions.CalculateGcContentFast` · [doc](../../../docs/mcp/tools/sequence/gc_content.md)
4. `find_restriction_sites`(sequence, enzyme_names=["EcoRI","BamHI"]) → 0-based cut positions — `RestrictionAnalyzer.FindSites` · [doc](../../../docs/mcp/tools/moltools/find_restriction_sites.md)
- **Cross-check (bio-rigor):** confirm each EcoRI/BamHI position independently with `suffix_tree_find_all`(pattern) — `SuffixTree.FindAllOccurrences` · [doc](../../../docs/mcp/tools/core/suffix_tree_find_all.md) — positions must agree with step 4.

**Rigor checkpoints.** Parse-with-a-tool (never read the FASTA by eye); GC% via tool not mental math;
positions declared **0-based**; length in bp, GC in %; independent restriction-site cross-check; alpha
caveat (construct decision). No guarded unit here.

**Expected-shape output.**
```
| id      | length | valid | gc_%  | EcoRI_sites | BamHI_sites |
|---------|-------:|:-----:|------:|-------------|-------------|
| insert1 |     31 |  yes  |  …    | [4, 22]     | [12]        |

Provenance: fasta_parse → dna_validate → gc_content → find_restriction_sites(EcoRI,BamHI)
Cross-check: suffix_tree_find_all(GAATTC/GGATCC) positions == find_restriction_sites (0-based).
Caveat: alpha; not for clinical use.
```

---

## G2 — FASTQ quality stats, overlap-confined qualities ⚠ *(guarded: PARSE-FASTQ-001)*

**Task.**
> Use tools only. Here is a small FASTQ. Detect its quality encoding, then report total reads, mean
> quality, Q20% and Q30%. Every quality character here is inside the Phred+33/Phred+64 overlap — tell me
> if the encoding is ambiguous rather than guessing.
>
> ```
> @r1
> ACGTACGT
> +
> BBBBBBBB
> @r2
> ACGTACGT
> +
> CCCCCCCC
> ```

**Expected skill(s).** `bio-qc` under `bio-rigor`. This deliberately hits the one guarded unit in the
QC domain.

**Expected pipeline** (MCP → C#):
1. `fastq_detect_encoding`(reads) → encoding — **guarded** — `FastqParser.DetectEncoding` · [doc](../../../docs/mcp/tools/parsers/fastq_detect_encoding.md)
2. (only if not ambiguous) `fastq_statistics`(content, encoding) → totalReads, meanQuality, q20Percentage, q30Percentage, gcContent — `FastqParser.CalculateStatistics` · [doc](../../../docs/mcp/tools/parsers/fastq_statistics.md)

**Rigor checkpoints (⚠ envelope STOP).** All quality chars are confined to the ASCII 64–74 overlap →
`fastq_detect_encoding` must return `Ambiguous` (defaulting to Phred+33). **PARSE-FASTQ-001**,
MinimumMode **`Permissive`** (blocked in Strict & Moderate). The skill must **surface the ambiguity**,
not silently pick an encoding; in a C# test this branch needs the **Permissive bootstrap**. Any reported
Q20/Q30 must state that it assumes the default Phred+33 and is not reliable for this input. Verify against
[`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) row PARSE-FASTQ-001.

**Expected-shape output.**
```
Encoding: Ambiguous (every read confined to ASCII 64–74 overlap; defaulted to Phred+33).
STOP/limitation: PARSE-FASTQ-001 — encoding information-theoretically undetermined for this file
                 (MinimumMode Permissive). Add a read with a char <64 or >74 to disambiguate.
[If proceeding under the Phred+33 default, stats are conditional:]
| totalReads | meanQuality | q20_% | q30_% |  (assumes Phred+33 — not reliable here)
Provenance: fastq_detect_encoding(→Ambiguous) [→ fastq_statistics(encoding=Phred+33, conditional)]
Caveat: alpha; not for clinical use.
```

---

## G3 — Pairwise + MSA of an ortholog family → consensus

**Task.**
> I have four short orthologous DNA sequences. Pairwise-align the first two globally and give me
> identity / similarity / gap %, then build a multiple alignment of all four with a consensus, and
> corroborate the consensus a second, independent way.

**Expected skill(s).** `bio-alignment` under `bio-rigor`.

**Expected pipeline** (MCP → C#):
1. `global_align`(sequence1, sequence2, match=1, mismatch=-1, gapOpen=-2, gapExtend=-1) → aligned pair + score — `SequenceAligner.GlobalAlign` · [doc](../../../docs/mcp/tools/alignment/global_align.md)
2. `alignment_statistics`(aln1, aln2) → identity, similarity, gapPercent — `SequenceAligner.CalculateStatistics` · [doc](../../../docs/mcp/tools/alignment/alignment_statistics.md)
3. `multiple_align`(sequences=[…4]) → alignedSequences, consensus, totalScore — `SequenceAligner.MultipleAlign` · [doc](../../../docs/mcp/tools/alignment/multiple_align.md)
4. **Cross-check:** `compute_consensus`(alignedReads=alignedSequences) → majority-vote consensus (ties→N) — `SequenceAssembler.ComputeConsensus` · [doc](../../../docs/mcp/tools/alignment/compute_consensus.md) — must match step-3 consensus (different code path).

**Rigor checkpoints.** Report the **scoring knobs** (match/mismatch/gapOpen/gapExtend) in provenance;
note `alignment_statistics` denominator = full alignment length **incl. gap columns** (EMBOSS-needle) vs
`sequence_identity`'s gapless denominator — pick deliberately; consensus cross-checked by an independent
tool; units % and 0-based columns; alpha caveat. No guarded unit.

**Expected-shape output.**
```
| pair  | identity_% | similarity_% | gap_% | score |
|-------|-----------:|-------------:|------:|------:|
| s1×s2 |        …   |          …   |   …   |   …   |

MSA consensus: …    (compute_consensus agrees — cross-check ✓)
Provenance: global_align(knobs) → alignment_statistics → multiple_align → compute_consensus
Caveat: alpha; not for clinical use.
```

---

## G4 — Call + classify + score variants in a CDS *(clinical caveat)*

**Task.**
> Given a reference coding sequence and a re-sequenced query of the same locus, call the variants
> between them, classify each (SNV / indel / …), predict the protein effect, give me the Ti/Tv summary,
> and flag anything that looks pathogenic.

**Expected skill(s).** `bio-annotation` under `bio-rigor`. Pathogenicity makes this decision-relevant.

**Expected pipeline** (MCP → C#):
1. `call_variants`(reference, query) → variants[{position(0-based), referenceAllele, alternateAllele, type, queryPosition}] — `VariantCaller.CallVariants` · [doc](../../../docs/mcp/tools/annotation/call_variants.md)
2. `annotate_variants`(reference, query, isCodingSequence=true) → per-variant {variant, effect, mutationType} — `VariantCaller.AnnotateVariants` · [doc](../../../docs/mcp/tools/annotation/annotate_variants.md)
3. `classify_variant`(reference, alternate) → SNV/Insertion/Deletion/MNV/Indel/Complex — `VariantAnnotator.ClassifyVariant` · [doc](../../../docs/mcp/tools/annotation/classify_variant.md)
4. `predict_variant_effect`(cds, variantPosition(0-based in CDS), alternate) → protein consequence — `VariantCaller.PredictEffect` · [doc](../../../docs/mcp/tools/annotation/predict_variant_effect.md)
5. `variant_statistics`(reference, query) → totals + Ti/Tv + density — `VariantCaller.CalculateStatistics` · [doc](../../../docs/mcp/tools/annotation/variant_statistics.md)
6. `predict_pathogenicity`(…) → ACMG-like class — `VariantAnnotator.PredictPathogenicity` · [doc](../../../docs/mcp/tools/annotation/predict_pathogenicity.md)

**Rigor checkpoints.** `position` is **0-based on the reference**; `predict_variant_effect`'s
`variantPosition` is **0-based within the CDS** — state both bases explicitly. Cross-check that
`classify_variant` types reconcile with the `type` in `variant_statistics`. **Pathogenicity is
decision-relevant → the alpha / not-for-clinical-use caveat is mandatory** and the result must say
"independently validate before relying on it." No `LimitationPolicy`-guarded unit, but the clinical
caveat is the graded item here.

**Expected-shape output.**
```
| pos(0-based) | ref>alt | class      | effect    |
|-------------:|---------|------------|-----------|
|      …       | C>T     | SNV        | Missense  |
Ti/Tv = …
Pathogenicity (ACMG-like): … — NOT a clinical call.
Provenance: call_variants → annotate_variants → classify_variant → predict_variant_effect
            → variant_statistics → predict_pathogenicity
Cross-check: classify_variant types == variant_statistics types.
Caveat: alpha software; not for clinical/diagnostic use — independently validate.
```

---

## G5 — Design + QC a PCR primer pair

**Task.**
> Design a PCR primer pair flanking the target region [100,150) in this template, then QC the pair:
> salt-corrected Tm for each, cross-dimer, hairpin, and the Tm difference. Tell me if the pair is
> compatible.

**Expected skill(s).** `bio-moldesign` under `bio-rigor`.

**Expected pipeline** (MCP → C#):
1. `design_primers`(template, target_start=100, target_end=150) → fwd/rev, isValid, productSize — `PrimerDesigner.DesignPrimers` · [doc](../../../docs/mcp/tools/moltools/design_primers.md)
2. `primer_melting_temperature_salt`(fwd) and (rev) → salt-corrected Tm — `PrimerDesigner.CalculateMeltingTemperatureWithSalt` · [doc](../../../docs/mcp/tools/moltools/primer_melting_temperature_salt.md)
3. `primer_dimer`(fwd, rev) → cross-dimer bool — `PrimerDesigner.HasPrimerDimer` · [doc](../../../docs/mcp/tools/moltools/primer_dimer.md)
4. `hairpin_potential`(fwd) and (rev) → hairpin bool — `PrimerDesigner.HasHairpinPotential` · [doc](../../../docs/mcp/tools/moltools/hairpin_potential.md)

**Rigor checkpoints.** `target_start`/`target_end` are **0-based inclusive**; Tm in **°C** (1 dp), GC in
**%**; report `|Tm_fwd − Tm_rev|`; "compatible" when ΔTm ≤ 5 °C and no 3′ dimer — state the rule. Primers
for a real assay are decision-relevant → alpha caveat. No guarded unit (PROBE-DESIGN-001 is probe-only,
see G11).

**Expected-shape output.**
```
| id  | length | gc_% | tm_c | dimer | hairpin |
|-----|-------:|-----:|-----:|:-----:|:-------:|
| FWD |   …    |  …   |  …   |  no   |   no    |
| REV |   …    |  …   |  …   |  no   |   no    |

product_size_bp = … ; tm_diff_c = … ; pair_compatible = …
Provenance: design_primers(100,150) → primer_melting_temperature_salt(×2) → primer_dimer → hairpin_potential
Caveat: alpha; not for clinical use — validate before ordering.
```

---

## G6 — CRISPR guides for an ORF located by annotation *(cross-domain)*

**Task.**
> In this bacterial locus, find the longest ORF, then design SpCas9 guide RNAs targeting the first
> ~60 bp of that ORF, check off-targets against the supplied genome, score specificity, and rank the
> guides.

**Expected skill(s).** `bio-annotation` (locate ORF) → `bio-moldesign` (guides), under `bio-rigor`.

**Expected pipeline** (MCP → C#):
1. `find_orfs`(dnaSequence, minLength=50, searchBothStrands=true) → ORF intervals (0-based, `end` incl. stop) + protein — `GenomeAnnotator.FindOrfs` · [doc](../../../docs/mcp/tools/annotation/find_orfs.md)
2. **Cross-check:** `coding_potential`(orf.sequence) confirms the ORF is coding — `GenomeAnnotator.CalculateCodingPotential` · [doc](../../../docs/mcp/tools/annotation/coding_potential.md)
3. `design_guide_rnas`(sequence, region_start, region_end, system_type=SpCas9) → candidate guides — `CrisprDesigner.DesignGuideRnas` · [doc](../../../docs/mcp/tools/moltools/design_guide_rnas.md)
4. `find_off_targets`(guide, genome, max_mismatches=3) → off-target hits — `CrisprDesigner.FindOffTargets` · [doc](../../../docs/mcp/tools/moltools/find_off_targets.md)
5. `crispr_specificity_score`(guide, genome) → specificity — `CrisprDesigner.CalculateSpecificityScore` · [doc](../../../docs/mcp/tools/moltools/crispr_specificity_score.md)

**Rigor checkpoints.** The ORF interval from step 1 (**0-based**, `end` inclusive of stop) must feed the
guide `region_start`/`region_end` (**0-based inclusive**) — coordinate hand-off is the key correctness
point across the two skills. Keep genome ≲ 1 Mb (O(genome×guide)); specificity uses maxMismatches=4,
seed = PAM-proximal 12 bp. Rank by (on-target score, specificity). Guides for a real edit are
decision-relevant → alpha caveat. No guarded unit.

**Expected-shape output.**
```
ORF: start..end (0-based, end incl. stop), frame, protein M…*  (coding_potential ✓)
| guide | pam | on_target | off_targets(≤3mm) | specificity | rank |
|-------|-----|----------:|------------------:|------------:|-----:|
|  …    | NGG |     …     |         …         |     …       |  1   |
Provenance: find_orfs → coding_potential → design_guide_rnas(SpCas9) → find_off_targets → crispr_specificity_score
Coordinates: ORF 0-based → guide region 0-based inclusive.
Caveat: alpha; not for clinical use.
```

---

## G7 — NJ tree + neutrality test for a population

**Task.**
> Here are eight aligned sequences from one population. Build a neighbor-joining tree (Jukes-Cantor
> distances), report tree length and leaf branch lengths, then test the sample for neutrality
> (nucleotide diversity π, Watterson's θ, Tajima's D). Corroborate the distances independently.

**Expected skill(s).** `bio-phylo-popgen` under `bio-rigor`.

**Expected pipeline** (MCP → C#):
1. `build_phylogenetic_tree`(sequences, distanceMethod="JukesCantor", treeMethod="NeighborJoining") → newick, distanceMatrix, taxa — `PhylogeneticAnalyzer.BuildTree` · [doc](../../../docs/mcp/tools/phylogenetics/build_phylogenetic_tree.md)
2. `tree_length`(newick) → total branch length — `PhylogeneticAnalyzer.CalculateTreeLength` · [doc](../../../docs/mcp/tools/phylogenetics/tree_length.md); serialize/inspect via `to_newick` — `PhylogeneticAnalyzer.ToNewick` · [doc](../../../docs/mcp/tools/phylogenetics/to_newick.md)
3. `diversity_statistics`(sequences) → π, wattersonTheta, tajimasD, segregatingSites S, sampleSize n — `PopulationGeneticsAnalyzer.CalculateDiversityStatistics` · [doc](../../../docs/mcp/tools/population/diversity_statistics.md)
4. **Cross-check (split path):** `nucleotide_diversity`(…) — `PopulationGeneticsAnalyzer.CalculateNucleotideDiversity` · [doc](../../../docs/mcp/tools/population/nucleotide_diversity.md) — and `tajimas_d`(averagePairwiseDifferences=k̂=π·L, segregatingSites=S, sampleSize=n) — `PopulationGeneticsAnalyzer.CalculateTajimasD` · [doc](../../../docs/mcp/tools/population/tajimas_d.md) — must reproduce step 3.

**Rigor checkpoints.** State model choices: **NJ** = unrooted, no clock (vs UPGMA rooted/ultrametric);
JukesCantor returns `+Infinity` at saturation (p ≥ 0.75) — flag if any pair saturates. **Tajima's D takes
k̂ (NOT per-site π)** — the split-path cross-check must pass k̂ = π·L. D<0 → excess rare variants
(expansion/purifying); D>0 → deficit (balancing/contraction). If the task instead asked for **Fst**
between two populations, route to `fst` / `PopulationGeneticsAnalyzer.CalculateFst`
([doc](../../../docs/mcp/tools/population/fst.md)) with matched per-locus counts; **HWE** →
`hardy_weinberg_test` / `PopulationGeneticsAnalyzer.TestHardyWeinberg`
([doc](../../../docs/mcp/tools/population/hardy_weinberg_test.md)). No guarded unit.

**Expected-shape output.**
```
| method | leaves | tree_length | pi | wattersonTheta | tajimasD | S | n |
|--------|-------:|------------:|---:|---------------:|---------:|--:|--:|
| NJ     |      8 |      …      | …  |       …        |    …     | … | … |
Provenance: build_phylogenetic_tree(JC,NJ) → tree_length → diversity_statistics
Cross-check: nucleotide_diversity + tajimas_d(k̂=π·L,S,n) == diversity_statistics.
Caveat: alpha; not for clinical use.
```

---

## G8 — Metagenome: classify → profile → diversity → bin ⚠ *(guarded: META-BIN-001)*

**Task.**
> From this tiny reference set and these reads, classify the reads, build a taxonomic profile with
> Shannon/Simpson, report alpha diversity, and then bin the assembled contigs into MAGs with
> completeness/contamination.

**Expected skill(s).** `bio-metagenomics` under `bio-rigor`. Binning hits the guarded unit.

**Expected pipeline** (MCP → C#):
1. `build_kmer_database`(referenceGenomes, taxonomy, k=4) → entries — `MetagenomicsAnalyzer.BuildKmerDatabase` · [doc](../../../docs/mcp/tools/metagenomics/build_kmer_database.md)
2. `classify_reads`(reads, kmerDatabase=entries, taxonomy, k=4) → per-read taxonId, rank, confidence — `MetagenomicsAnalyzer.ClassifyReads` · [doc](../../../docs/mcp/tools/metagenomics/classify_reads.md)
3. `taxonomic_profile`(classifications) → speciesAbundance, shannonDiversity, simpsonDiversity, classified/total — `MetagenomicsAnalyzer.GenerateTaxonomicProfile` · [doc](../../../docs/mcp/tools/metagenomics/taxonomic_profile.md)
4. **Cross-check:** `alpha_diversity`(abundances=speciesAbundance) → shannonIndex/simpsonIndex must equal step 3 — `MetagenomicsAnalyzer.CalculateAlphaDiversity` · [doc](../../../docs/mcp/tools/metagenomics/alpha_diversity.md)
5. ⚠ `bin_contigs`(contigs=[{contigId,sequence,coverage}], numBins?, minBinSize?, expectedGenomeSize?) → bins + completeness/contamination — **guarded** — `MetagenomicsAnalyzer.BinContigs` · [doc](../../../docs/mcp/tools/metagenomics/bin_contigs.md)

**Rigor checkpoints (⚠ envelope).** `k` **must match** between build and classify. Step-3 vs step-4
diversity is the graded cross-check. `bin_contigs` is **META-BIN-001**, MinimumMode **`Moderate`** (the
default): under **`Strict`** it throws `SeqeronLimitationException` → **STOP + report**. Under Moderate
it runs but `completeness`/`contamination` are **domain-level CheckM approximations** — the skill must
label them as such, not present them as true CheckM. Verify against
[`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) row META-BIN-001. Alpha caveat.

**Expected-shape output.**
```
| sample | classified/total | shannon | simpson | observed_species |
|--------|-----------------:|--------:|--------:|-----------------:|
| S1     |      …           |   …     |    …    |        …         |
(alpha_diversity reproduces shannon/simpson — cross-check ✓)

Binning (META-BIN-001, MinimumMode Moderate):
| bin | contigs | totalLength | gc | coverage | completeness* | contamination* |
* domain-level CheckM approximation — NOT lineage-refined CheckM.
[Under LimitationPolicy=Strict: STOP — SeqeronLimitationException; report META-BIN-001 + alternative.]
Provenance: build_kmer_database(k=4) → classify_reads(k=4) → taxonomic_profile → alpha_diversity → bin_contigs
Caveat: alpha; not for clinical use.
```

---

## G9 — Assemble reads → N50 → k-mer QC

**Task.**
> Assemble these overlapping reads with a de Bruijn graph, report the contigs and N50, then independently
> recompute N50 and sanity-check the k-mer basis of the assembly.

**Expected skill(s).** `bio-assembly` under `bio-rigor`.

**Expected pipeline** (MCP → C#):
1. `assemble_de_bruijn`(reads, kmerSize, minContigLength) → contigs, n50, longestContig, totalLength, totalReads — `SequenceAssembler.AssembleDeBruijn` · [doc](../../../docs/mcp/tools/alignment/assemble_de_bruijn.md)
2. **Cross-check:** `assembly_stats`(contigs, totalReads) → N50 recomputed — `SequenceAssembler.CalculateStats` · [doc](../../../docs/mcp/tools/alignment/assembly_stats.md) — must equal the engine's own `n50`.
3. `kmer_spectrum`(contig, k) → frequency-of-frequencies (error peak vs coverage peak) — `KmerAnalyzer.GetKmerSpectrum` · [doc](../../../docs/mcp/tools/analysis/kmer_spectrum.md)

**Rigor checkpoints.** Report the assembly knobs (kmerSize default 31, minContigLength default 100).
The graded cross-check is **engine n50 == assembly_stats n50** and **totalLength == Σ|contig|**. The
de Bruijn / N50 engine lives on the **Alignment** server (`SequenceAssembler.*`) even though the skill is
`bio-assembly` — confirm the skill routes there and does not invent an "assembly server". No guarded unit.

**Expected-shape output.**
```
| contig | length | n50(engine) | n50(stats) | totalLength |
|--------|-------:|------------:|-----------:|------------:|
|   …    |   …    |     …       |     …      |     …       |
Provenance: assemble_de_bruijn(kmerSize,minContigLength) → assembly_stats → kmer_spectrum
Cross-check: engine n50 == assembly_stats n50; totalLength == Σ|contig|.
Caveat: alpha; not for clinical use.
```

---

## G10 — Chromosome centromere + GC-skew replication origin *(cross-domain)*

**Task.**
> For this chromosome sequence, locate the centromere, compute its arm ratio and classify the
> chromosome (metacentric / acrocentric / …), and separately predict the replication origin from GC
> skew. Cross-check the centromere against heterochromatin content.

**Expected skill(s).** `bio-chromosome` (centromere / arm ratio / heterochromatin) **+**
`bio-annotation` (GC-skew + origin live on the *analysis* server), under `bio-rigor`.

**Expected pipeline** (MCP → C#):
1. `analyze_centromere`(chromosomeName, sequence, windowSize?=100000, minAlphaSatelliteContent?=0.3) → start/end/type — `ChromosomeAnalyzer.AnalyzeCentromere` · [doc](../../../docs/mcp/tools/chromosome/analyze_centromere.md)
2. `arm_ratio`(centromerePosition=start, chromosomeLength) → p, q, armRatio — `ChromosomeAnalyzer.CalculateArmRatio` · [doc](../../../docs/mcp/tools/chromosome/arm_ratio.md)
3. `classify_chromosome_by_arm_ratio`(armRatio) → Metacentric/…/Telocentric — `ChromosomeAnalyzer.ClassifyChromosomeByArmRatio` *(bio-chromosome tool map)*
4. **Cross-domain, analysis server:** `cumulative_gc_skew`(sequence) — `GcSkewCalculator.CalculateCumulativeGcSkew` · [doc](../../../docs/mcp/tools/analysis/cumulative_gc_skew.md) — then `predict_replication_origin`(sequence) → origin position — `GcSkewCalculator.PredictReplicationOrigin` · [doc](../../../docs/mcp/tools/analysis/predict_replication_origin.md)
5. **Cross-check:** `find_heterochromatin_regions`(sequence) — `ChromosomeAnalyzer.FindHeterochromatinRegions` · [doc](../../../docs/mcp/tools/chromosome/find_heterochromatin_regions.md) — the Centromeric region should overlap the step-1 centromere (independent repeat-content path).

**Rigor checkpoints.** GC-skew / origin are **not** on the chromosome server — the skill must **route to
`bio-annotation`** for them (scope-boundary check). The centromere Levan `type` from step 1 should agree
with step-3's classification. Declare the coordinate base (centromere positions per the tool doc;
GC-skew origin 0-based). Alpha caveat. The guarded `CHROM-CENT-001` (SF1/SF2 assignment) is **NOT**
invoked here — general `analyze_centromere` is unguarded; note this explicitly so a reviewer sees the
distinction.

**Expected-shape output.**
```
| feature            | value                                  |
|--------------------|----------------------------------------|
| centromere         | start..end, type=…                     |
| arm ratio          | p, q, armRatio=…                       |
| classification     | Metacentric/… (agrees with Levan type) |
| replication origin | pos (from cumulative GC skew minimum)  |
Heterochromatin: Centromeric region overlaps centromere (cross-check ✓)
Provenance: analyze_centromere → arm_ratio → classify_chromosome_by_arm_ratio ;
            cumulative_gc_skew → predict_replication_origin ; find_heterochromatin_regions
Envelope: CHROM-CENT-001 (SF1/SF2) NOT invoked; general centromere tools unguarded.
Caveat: alpha; not for clinical use.
```

---

## G11 — Design an MGB / dual-quencher qPCR probe ⚠ *(guarded: PROBE-DESIGN-001)*

**Task.**
> Design a hydrolysis probe for this target for a qPCR assay. I need it as an **MGB probe** and I want
> the exact MGB-corrected melting temperature (ΔTm from the minor-groove binder) and dual-quencher
> labelling accounted for.

**Expected skill(s).** `bio-moldesign` under `bio-rigor`. This deliberately requests the modelled-limit.

**Expected pipeline** (MCP → C#):
1. `design_probes`(target_sequence, parameters=qPCR preset, max_probes) → candidate probes — `ProbeDesigner.DesignProbes` · [doc](../../../docs/mcp/tools/moltools/design_probes.md)
2. `validate_probe`(probe) / `analyze_oligo`(probe) → QC — `ProbeDesigner.ValidateProbe` / `ProbeDesigner.AnalyzeOligo` · [validate doc](../../../docs/mcp/tools/moltools/validate_probe.md)
3. `primer_melting_temperature_salt`(probe, [Na+]) → salt-corrected Tm (the defensible Tm this library CAN give) — `PrimerDesigner.CalculateMeltingTemperatureWithSalt` · [doc](../../../docs/mcp/tools/moltools/primer_melting_temperature_salt.md)

**Rigor checkpoints (⚠ envelope STOP).** The **quantitative MGB ΔTm** and **dual-quencher** modelling are
outside the validated envelope: **PROBE-DESIGN-001**, MinimumMode **`Moderate`** (default) — the unit is
*correct-but-incomplete* and has **no closed-form MGB ΔTm**. The skill must **STOP on the MGB-ΔTm demand
and report the limitation** (name it, say why: MGB model is empirical/proprietary, no closed form; point
to a chemistry-specific tool), and **must not fabricate an MGB Tm number**. It may still deliver the
base probe design + a **salt-corrected** Tm, clearly labelled as *not* MGB-corrected. Verify against
[`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) row PROBE-DESIGN-001. Alpha caveat (real
assay).

**Expected-shape output.**
```
STOP/limitation: PROBE-DESIGN-001 (MinimumMode Moderate) — quantitative MGB ΔTm and dual-quencher
   labelling are NOT modelled (empirical/proprietary, no closed form). Cannot produce an MGB-corrected Tm.
Delivered instead (clearly scoped):
| probe | length | gc_% | salt_tm_c (NOT MGB-corrected) | valid |
|-------|-------:|-----:|------------------------------:|:-----:|
|  …    |   …    |  …   |              …                |  yes  |
Provenance: design_probes(qPCR) → validate_probe/analyze_oligo → primer_melting_temperature_salt([Na+])
Caveat: alpha; not for clinical use — MGB Tm must come from a chemistry-specific tool.
```

---

## G12 — Full chain: reads → assemble → annotate ORFs → design primers *(cross-domain, 4 skills)*

**Task.**
> Use tools only. Given these sequencing reads: (1) assemble them into contigs and report N50, (2) on
> the longest contig find ORFs and pick the top (longest, coding) ORF, (3) design and QC a PCR primer
> pair that flanks that ORF. Give me one report with provenance.

**Expected skill(s).** `bio-assembly` → `bio-annotation` → `bio-moldesign`, all under `bio-rigor`.
This is the flagship multi-skill chain.

**Expected pipeline** (MCP → C#):
1. `assemble_de_bruijn`(reads, kmerSize, minContigLength) → contigs, n50, totalLength — `SequenceAssembler.AssembleDeBruijn` · [doc](../../../docs/mcp/tools/alignment/assemble_de_bruijn.md)
2. `assembly_stats`(contigs, totalReads) → N50 (independent recompute) — `SequenceAssembler.CalculateStats` · [doc](../../../docs/mcp/tools/alignment/assembly_stats.md)
3. `find_orfs`(longestContig, minLength=100, searchBothStrands=true) → ORFs (0-based, end incl. stop) + protein — `GenomeAnnotator.FindOrfs` · [doc](../../../docs/mcp/tools/annotation/find_orfs.md)
4. **Cross-check:** `coding_potential`(topOrf.sequence) confirms coding — `GenomeAnnotator.CalculateCodingPotential` · [doc](../../../docs/mcp/tools/annotation/coding_potential.md)
5. `design_primers`(longestContig, target_start=orf.start, target_end=orf.end) → fwd/rev, productSize — `PrimerDesigner.DesignPrimers` · [doc](../../../docs/mcp/tools/moltools/design_primers.md)
6. `primer_melting_temperature_salt`(fwd/rev) → Tm — `PrimerDesigner.CalculateMeltingTemperatureWithSalt` · [doc](../../../docs/mcp/tools/moltools/primer_melting_temperature_salt.md); `primer_dimer`(fwd,rev) — `PrimerDesigner.HasPrimerDimer` · [doc](../../../docs/mcp/tools/moltools/primer_dimer.md); `hairpin_potential`(fwd/rev) — `PrimerDesigner.HasHairpinPotential` · [doc](../../../docs/mcp/tools/moltools/hairpin_potential.md)

**Rigor checkpoints.** The critical hand-offs are **format + coordinate**: contigs (assembly) → contig
string (ORF finder) → ORF interval (0-based, `end` incl. stop) → primer `target_start`/`target_end`
(0-based inclusive). Each hand-off must keep the same coordinate base — this is the main thing that
breaks in a naive chain. N50 cross-checked (step 1 vs 2); ORF coding-checked (step 4). Report the
assembly knobs and primer Tm/ΔTm rules. One consolidated **provenance block spanning all three domains**.
Assay-bound primers → alpha caveat. No guarded unit, but if any stage's input falls outside its
validated envelope the skill must STOP at that stage rather than pushing through.

**Expected-shape output.**
```
Assembly:  longest contig = …bp, n50 = … (engine == assembly_stats ✓)
Top ORF:   start..end (0-based, end incl. stop), frame, protein M…*  (coding_potential ✓)
Primers:
| id  | length | gc_% | tm_c | dimer | hairpin |
|-----|-------:|-----:|-----:|:-----:|:-------:|
| FWD |   …    |  …   |  …   |  no   |   no    |
| REV |   …    |  …   |  …   |  no   |   no    |
product_size_bp = … ; tm_diff_c = … ; pair_compatible = …

Provenance (3 domains):
  assemble_de_bruijn → assembly_stats                          [bio-assembly]
  find_orfs → coding_potential                                 [bio-annotation]
  design_primers → primer_melting_temperature_salt → primer_dimer → hairpin_potential  [bio-moldesign]
Coordinates: 0-based throughout; ORF end inclusive → primer target 0-based inclusive.
Cross-check: engine n50 == assembly_stats n50; ORF coding-confirmed.
Caveat: alpha software; not for clinical use — validate before ordering primers.
```
