---
type: concept
title: "Algorithm validation evidence artifacts"
tags: [validation, testing]
sources:
  - docs/Evidence/ALIGN-GLOBAL-001-Evidence.md
  - docs/Evidence/ALIGN-MULTI-001-Evidence.md
  - docs/Evidence/ALIGN-SEMI-001-Evidence.md
  - docs/Evidence/ALIGN-STATS-001-Evidence.md
  - docs/Evidence/ANNOT-CODING-001-Evidence.md
  - docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
  - docs/Evidence/ANNOT-REPEAT-001-Evidence.md
  - docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md
  - docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
  - docs/Evidence/ASSEMBLY-COVER-001-Evidence.md
  - docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
  - docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md
  - docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
  - docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
  - docs/Evidence/ASSEMBLY-STATS-001-Evidence.md
  - docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
  - docs/Evidence/CHROM-ANEU-001-Evidence.md
  - docs/Evidence/CHROM-CENT-001-Evidence.md
  - docs/Evidence/CHROM-KARYO-001-Evidence.md
  - docs/Evidence/CHROM-SYNT-001-Evidence.md
  - docs/Evidence/CHROM-TELO-001-Evidence.md
  - docs/Evidence/CODON-CAI-001-Evidence.md
  - docs/Evidence/CODON-ENC-001-Evidence.md
  - docs/Evidence/CODON-OPT-001-Evidence.md
  - docs/Evidence/CODON-RARE-001-Evidence.md
  - docs/Evidence/CODON-RSCU-001-Evidence.md
  - docs/Evidence/CODON-STATS-001-Evidence.md
  - docs/Evidence/CODON-USAGE-001-Evidence.md
  - docs/Evidence/COMPGEN-ANI-001-Evidence.md
  - docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md
  - docs/Evidence/COMPGEN-COMPARE-001-Evidence.md
  - docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md
  - docs/Evidence/COMPGEN-ORTHO-001-Evidence.md
  - docs/Evidence/COMPGEN-RBH-001-Evidence.md
  - docs/Evidence/COMPGEN-REARR-001-Evidence.md
  - docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md
  - docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md
  - docs/Evidence/DISORDER-LC-001-Evidence.md
  - docs/Evidence/DISORDER-MORF-001-Evidence.md
  - docs/Evidence/DISORDER-PRED-001-Evidence.md
  - docs/Evidence/EPIGEN-AGE-001-Evidence.md
  - docs/Evidence/EPIGEN-BISULF-001-Evidence.md
  - docs/Evidence/EPIGEN-CHROM-001-Evidence.md
  - docs/Evidence/EPIGEN-CPG-001-Evidence.md
  - docs/Evidence/EPIGEN-DMR-001-Evidence.md
  - docs/Evidence/EPIGEN-METHYL-001-Evidence.md
  - docs/Evidence/GENOMIC-COMMON-001-Evidence.md
  - docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md
  - docs/Evidence/GENOMIC-ORF-001-Evidence.md
  - docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
  - docs/Evidence/KMER-BOTH-001-Evidence.md
  - docs/Evidence/KMER-DIST-001-Evidence.md
  - docs/Evidence/KMER-POSITIONS-001-Evidence.md
  - docs/Evidence/KMER-STATS-001-Evidence.md
  - docs/Evidence/KMER-UNIQUE-001-Evidence.md
  - docs/Evidence/META-ALPHA-001-Evidence.md
  - docs/Evidence/META-FUNC-001-Evidence.md
  - docs/Evidence/META-PATHWAY-001-Evidence.md
  - docs/Evidence/META-PROF-001-Evidence.md
  - docs/Evidence/META-RESIST-001-Evidence.md
  - docs/Evidence/META-TAXA-001-Evidence.md
  - docs/Evidence/MIRNA-PAIR-001-Evidence.md
  - docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md
  - docs/Evidence/MIRNA-SEED-001-Evidence.md
  - docs/Evidence/MIRNA-TARGET-001-Evidence.md
  - docs/Evidence/MOTIF-CONS-001-Evidence.md
  - docs/Evidence/MOTIF-DISCOVER-001-Evidence.md
  - docs/Evidence/MOTIF-SHARED-001-Evidence.md
  - docs/Evidence/ONCO-ACTION-001-Evidence.md
  - docs/Evidence/ONCO-ANNOT-001-Evidence.md
  - docs/Evidence/ONCO-ARTIFACT-001-Evidence.md
  - docs/Evidence/ONCO-ASCAT-001-Evidence.md
  - docs/Evidence/ONCO-CCF-001-Evidence.md
  - docs/Evidence/ONCO-CHIP-001-Evidence.md
  - docs/Evidence/ONCO-CLONAL-001-Evidence.md
  - docs/Evidence/ONCO-CNA-001-Evidence.md
source_commit: f794f2454365bdf09e3440ac95c4f59db1ff0553
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-global-001-evidence
      evidence: "Evidence Artifact: ALIGN-GLOBAL-001 ... Online Sources ... Test Dataset ... Deviations and Assumptions ... References"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:definition-of-done
      source: align-global-001-evidence
      evidence: "Deviations and Assumptions: None ... the implementation follows the standard Needleman–Wunsch linear gap penalty model exactly as given in the Wikipedia pseudocode"
      confidence: medium
      status: current
---

# Algorithm validation evidence artifacts

Each algorithm [[test-unit-registry|test unit]] has a per-unit **Evidence artifact** under
`docs/Evidence/<UnitID>-Evidence.md`. These are the literature-traced source record behind
the [[definition-of-done]]'s "Evidence documented" criterion and the
[[validation-and-testing]] campaign: they pin exactly which external references and worked
examples the implementation and its tests are validated against.

## Templated structure

Every evidence file follows the same shape:

1. **Header** — Test Unit ID, algorithm name, date collected.
2. **Online sources** — Wikipedia / primary-literature URLs with access dates and the key
   extracted points (definitions, recurrences, complexity, worked examples).
3. **Test dataset** — the canonical worked example(s) with exact parameters and expected
   outputs, used as the oracle / differential test fixture.
4. **Deviations and assumptions** — where the implementation departs from (or exactly
   follows) the reference, plus API-contract behaviours (null/empty handling) that sit
   outside the algorithm spec.
5. **References** — primary literature and encyclopedic citations.

Because these files are near-templated across the ~213 documented units, the wiki keeps
**one** shared page for the pattern (this page) plus a concise per-file source summary
(e.g. [[align-global-001-evidence]], [[align-multi-001-evidence]],
[[align-semi-001-evidence]], [[align-stats-001-evidence]], [[annot-coding-001-evidence]],
[[annot-codonusage-001-evidence]], [[annot-repeat-001-evidence]],
[[assembly-consensus-001-evidence]], [[assembly-correct-001-evidence]],
[[assembly-cover-001-evidence]], [[assembly-dbg-001-evidence]],
[[assembly-merge-001-evidence]], [[assembly-olc-001-evidence]],
[[assembly-scaffold-001-evidence]], [[assembly-stats-001-evidence]],
[[assembly-trim-001-evidence]], [[chrom-aneu-001-evidence]],
[[chrom-cent-001-evidence]], [[chrom-karyo-001-evidence]], [[chrom-synt-001-evidence]],
[[chrom-telo-001-evidence]], [[codon-cai-001-evidence]],
[[codon-enc-001-evidence]], [[codon-opt-001-evidence]], [[codon-rare-001-evidence]],
[[codon-rscu-001-evidence]], [[codon-stats-001-evidence]],
[[codon-usage-001-evidence]], [[compgen-ani-001-evidence]],
[[compgen-cluster-001-evidence]], [[compgen-compare-001-evidence]],
[[compgen-dotplot-001-evidence]], [[compgen-ortho-001-evidence]],
[[compgen-rbh-001-evidence]], [[compgen-rearr-001-evidence]],
[[compgen-reversal-001-evidence]], [[compgen-synteny-001-evidence]],
[[disorder-lc-001-evidence]], [[disorder-morf-001-evidence]],
[[disorder-pred-001-evidence]], [[epigen-age-001-evidence]],
[[epigen-bisulf-001-evidence]], [[epigen-chrom-001-evidence]],
[[epigen-cpg-001-evidence]], [[epigen-dmr-001-evidence]],
[[epigen-methyl-001-evidence]], [[genomic-common-001-evidence]],
[[genomic-motifs-001-evidence]], [[genomic-orf-001-evidence]],
[[genomic-repeat-001-evidence]], [[kmer-async-001-evidence]],
[[kmer-both-001-evidence]], [[kmer-dist-001-evidence]],
[[kmer-positions-001-evidence]], [[kmer-stats-001-evidence]],
[[kmer-unique-001-evidence]], [[meta-alpha-001-evidence]],
[[meta-func-001-evidence]], [[meta-pathway-001-evidence]],
[[meta-prof-001-evidence]], [[meta-resist-001-evidence]],
[[meta-taxa-001-evidence]], [[mirna-pair-001-evidence]],
[[mirna-precursor-001-evidence]], [[mirna-seed-001-evidence]],
[[mirna-target-001-evidence]], [[motif-cons-001-evidence]],
[[motif-discover-001-evidence]], [[motif-shared-001-evidence]],
[[onco-action-001-evidence]], [[onco-annot-001-evidence]],
[[onco-artifact-001-evidence]], [[onco-ascat-001-evidence]],
[[onco-ccf-001-evidence]], [[onco-chip-001-evidence]],
[[onco-clonal-001-evidence]], [[onco-cna-001-evidence]]). An
individual algorithm gets its own concept page only when it is itself distinct and wiki-worthy
— for example [[global-alignment-needleman-wunsch]], [[multiple-sequence-alignment]],
[[semi-global-alignment-fitting]], [[alignment-statistics]],
[[coding-potential-hexamer-score]], [[relative-synonymous-codon-usage]],
[[codon-adaptation-index]] (the CAI index in the codon-usage family),
[[effective-number-of-codons]] (the reference-free ENC/Nc measure in the codon-usage family),
[[codon-optimization]] (the CDS-rewriting operation of the codon-usage family),
[[rare-codon-analysis]] (the thresholded-frequency + %MinMax/Sherlocc cluster-detection unit of the codon-usage family),
[[codon-usage-comparison]] (the raw codon-count table + TVD distribution-comparison unit of the codon-usage family),
[[repetitive-element-detection]] (the anchor for the repeats/tandem family),
[[consensus-sequence]] (the anchor for the assembly CONSENSUS family),
[[kmer-spectrum-error-correction]] (the anchor for the assembly CORRECT family),
[[coverage-depth-calculation]] (the anchor for the assembly COVER family),
[[de-bruijn-graph-assembly]] (the anchor for the assembly DBG family),
[[contig-merge-overlap-collapse]] (the anchor for the assembly MERGE family),
[[overlap-layout-consensus-assembly]] (the anchor for the assembly OLC family), or
[[scaffolding]] (the anchor for the assembly SCAFFOLD family),
[[assembly-statistics]] (the anchor for the assembly STATS family), or
[[quality-trimming-running-sum]] (the anchor for the assembly TRIM family), or
[[aneuploidy-detection]] (the anchor for the chromosome-analysis copy-number/ploidy family), or
[[centromere-analysis]] (the anchor for the chromosome centromere / alpha-satellite family), or
[[karyotype-analysis]] (the anchor for the chromosome karyotyping / ploidy-detection family), or
[[synteny-and-rearrangement-detection]] (the shared anchor for the chromosome + comparative-genomics synteny/rearrangement family), or
[[telomere-analysis]] (the anchor for the chromosome telomere family), or
[[average-nucleotide-identity]] (the anchor for the comparative-genomics ANI genome-similarity family), or
[[conserved-gene-clusters-common-intervals]] (the comparative-genomics common-interval / conserved-cluster unit), or
[[genome-comparison-core-dispensable]] (the comparative-genomics end-to-end genome-comparison pipeline — core/dispensable partition + syntenic fraction), or
[[dot-plot-word-match]] (the comparative-genomics word-match / k-tuple dot-matrix visual sequence-comparison unit), or
[[ortholog-detection-reciprocal-best-hits]] (the comparative-genomics RBH ortholog + within-genome in-paralog detection unit, shared anchor for COMPGEN-RBH-001), or
[[genome-rearrangement-breakpoint-distance]] (the comparative-genomics signed-permutation / breakpoint-distance rearrangement-detection unit, the alternative formulation to the block-signal [[synteny-and-rearrangement-detection]]), or
[[protein-low-complexity-seg]] (the anchor for the protein disorder / features family — SEG low-complexity region detection, the first ingested DISORDER-* unit), or
[[morf-prediction-dip-in-disorder]] (the MoRF dip-in-disorder prediction unit of the protein-disorder family, the second ingested DISORDER-* unit), or
[[intrinsic-disorder-prediction-top-idp]] (the TOP-IDP `PredictDisorder` sliding-window intrinsic-disorder anchor of the protein-disorder family, the third ingested DISORDER-* unit that MoRF + region detection sit on), or
[[epigenetic-age-horvath-clock]] (the Horvath DNAm methylation-clock epigenetic-age estimator, the first ingested unit of the Epigenetics family), or
[[bisulfite-methylation-calling]] (the bisulfite-sequencing conversion/methylation-calling/profile unit of the Epigenetics family, the second ingested EPIGEN-* unit — distinct from the age clock), or
[[chromatin-state-prediction]] (the ChromHMM-style histone-mark chromatin-state annotation unit of the Epigenetics family, the third ingested EPIGEN-* unit — operates on histone ChIP-seq marks, not DNA methylation), or
[[cpg-island-detection]] (the sequence-only CpG site enumeration / O-E ratio / Gardiner-Garden & Frommer CpG-island detection unit of the Epigenetics family, the fourth ingested EPIGEN-* unit — touches no methylation state), or
[[differentially-methylated-regions]] (the methylKit tiling-window + Fisher's-exact-test two-sample DMR-detection unit of the Epigenetics family, the fifth ingested EPIGEN-* unit — compares methylation between two samples, consuming the bisulfite β-values), or
[[methylation-context-classification]] (the sequence-only CpG/CHG/CHH trinucleotide context classifier + methylation-site enumerator of the Epigenetics family, the sixth and final ingested EPIGEN-* unit — the IUPAC H = "not G" classification that partitions non-CpG cytosines, sharing the bisulfite unit's weighted-profile aggregator), or
[[longest-common-substring]] (the Sequence-Comparison longest-common-substring / common-region-detection unit — the contiguous maximal shared substring of two strings via a generalized suffix tree, with a documented deterministic tie-break), or
[[known-motif-search]] (the Motif-Analysis known-motif-search unit — multi-pattern **exact** substring matching of a set of known query motifs with all-overlapping-occurrences reporting, the exact-equality baseline of the motif family), or
[[open-reading-frame-detection]] (the Analysis ORF-detection unit — six-frame ATG→first-in-frame-stop enumeration `GenomicAnalyzer.FindOpenReadingFrames`, ATG-only / nucleotide-`minLength`, nested-ORFs-sharing-a-stop reported; distinct from the annotation-layer `GenomeAnnotator.FindOrfs` ANNOT-ORF-001), or
[[longest-repeated-substring]] (the Repeat-Analysis LRS + all-repeats-enumeration unit — `GenomicAnalyzer.FindLongestRepeat` / `FindRepeats` via the deepest-internal-node-with-≥2-leaves of a *single-string* suffix tree, the one-string sibling of [[longest-common-substring]] and distinct from the tandem/inverted [[repetitive-element-detection]] anchor), or
[[asynchronous-kmer-counting]] (the K-mer family's cancelable/progress-reporting `Task.Run` count wrapper — first ingested K-mer unit, numeric result identical to the synchronous count), or
[[both-strand-kmer-counting]] (the K-mer family's additive / kPAL-"balance" strand-aware count — `count[w] = forward[w] + forward[RC(w)]` summed over both strands, distinct from the not-implemented canonical-collapsing convention of Jellyfish `-C` / Mash), or
[[k-mer-euclidean-distance]] (the K-mer family's alignment-free word-frequency **distance** — Euclidean L2 over normalized k-mer frequency vectors `count/(L−k+1)`, an `alternative_to` the presence/absence [[kmer-jaccard-similarity]], the third ingested K-mer unit), or
[[k-mer-positions]] (the K-mer family's *occurrence-index* unit — `KmerAnalyzer.FindKmerPositions` returns the ascending 0-based positions where a given k-mer occurs in a sequence, *where* not *how many*; the single-pattern sibling of the multi-pattern exact-matcher [[known-motif-search]], distinct from the counting siblings), or
[[k-mer-statistics]] (the K-mer family's *summary-statistics* unit — `KmerAnalyzer.AnalyzeKmers` reduces the count profile to Total/Unique(=distinct)/Max/Min/Average multiplicity + Shannon **k-entropy** in bits, a companion summary layer over the shared `CountKmers`, distinct from the counting/generation/distance siblings), or
[[unique-and-mincount-kmers]] (the K-mer family's *frequency-filtering* unit — `KmerAnalyzer.FindUniqueKmers` returns the `Count == 1` singleton set and `FindKmersWithMinCount` the `Count ≥ t` recurrent k-mers ordered by count; the "unique"=singleton notion is distinct from [[k-mer-statistics]]'s `UniqueKmers`=distinct-count field), or
[[alpha-diversity]] (the anchor for the Metagenomics diversity family, the first ingested META-* unit — `MetagenomicsAnalyzer.CalculateAlphaDiversity` computes the six within-sample indices Shannon/Simpson/inverse-Simpson/Pielou-evenness/observed-richness/Chao1 over one taxon→abundance map; between-sample [[beta-diversity]] is a future sibling), or
[[functional-prediction]] (the Metagenomics family's functional-prediction unit, META-FUNC-001 — homology-based annotation transfer `MetagenomicsAnalyzer.PredictFunctions` scored with BLAST bit-score / E-value Karlin-Altschul statistics (ungapped BLOSUM62 λ=0.3176/K=0.134, best-hit = lowest E-value) plus the hypergeometric pathway ORA `FindPathwayEnrichment`; the ORA half is owned by [[pathway-enrichment-ora]]), or
[[pathway-enrichment-ora]] (the Metagenomics family's pathway-enrichment / over-representation-analysis unit, META-PATHWAY-001 — the GO::TermFinder / clusterProfiler hypergeometric right-tail test `FindPathwayEnrichment` / `HypergeometricUpperTail`; the dedicated owner of the ORA machinery that [[functional-prediction]] exercises as component B), or
[[taxonomic-profile]] (the Metagenomics family's community-abundance profiling unit, META-PROF-001 — `MetagenomicsAnalyzer.GenerateTaxonomicProfile` aggregates the per-read output of [[taxonomic-classification]] into normalized per-taxon relative abundances at four ranks plus inline Shannon/Simpson diversity and read counts; the aggregation step that the classification unit deliberately deferred, sibling of [[alpha-diversity]] / [[beta-diversity]]), or
[[antibiotic-resistance-gene-detection]] (the Metagenomics family's AMR gene-detection unit, META-RESIST-001 — the ResFinder-style `MetagenomicsAnalyzer.FindAntibioticResistanceGenes` screen of contigs against a caller-supplied resistance-gene DB, reporting the best-matching gene per contig by BLAST percent identity + reference coverage under a dual threshold (default 0.90 ID / 0.60 cov); a BLAST-style homology screen sharing machinery with [[functional-prediction]] but scoring nucleotide identity/coverage rather than a BLOSUM62 bit-score), or
[[significant-taxa-detection]] (the Metagenomics family's differential-abundance unit, META-TAXA-001 — the per-taxon two-group **Mann–Whitney U / Wilcoxon rank-sum** test `MetagenomicsAnalyzer.MannWhitneyU` / `FindSignificantTaxa` with midrank tie correction, optional 0.5 continuity correction, and the asymptotic normal-approximation two-tailed p-value `2·(1−Φ(z))`; consumes the per-sample [[taxonomic-profile]] abundance vectors, distinct from the hypergeometric [[pathway-enrichment-ora]] and the Fisher's-exact [[differentially-methylated-regions]] by its rank-sum test), or
[[rna-base-pairing]] (the MiRNA family anchor + shared RNA base-pairing primitive, MIRNA-PAIR-001 — the Watson-Crick {A-U,G-C} + standard G-U wobble pairing rule (`CanPair`/`IsWobblePair`), the antiparallel seed reverse complement (`GetReverseComplement`), and the ungapped antiparallel miRNA-target duplex `AlignMiRnaToTarget` with a sign-only-reliable simplified Turner-2004 stacking ΔG; the reusable pairing rule that RNA secondary-structure classification also builds on, distinct from target-site efficacy prediction), or
[[pre-mirna-hairpin-detection]] (the MiRNA family's precursor-hairpin unit, MIRNA-PRECURSOR-001 — `MiRnaAnalyzer.FindPreMiRnaHairpins` stem-loop detection (consecutive-pairing stem ≥18 bp, loop 3-25 nt, mature/star arm extraction, dot-bracket + Turner ΔG) plus three opt-in production paths: real-MFE-fold assessment `AssessHairpinByMfe` (RNA-STRUCT-001 Zuker–Stiegler, single-hairpin + stem ≥16 + MFEI ≥0.85 — detects the natural miRBase precursors the heuristic rejects), Drosha/Dicer cleavage-ruler prediction `PredictDroshaDicerCleavage`, and a trained logistic-regression natural-vs-background classifier `ClassifyPreMiRna`; depends on the [[rna-base-pairing]] pairing primitive), or
[[seed-sequence-analysis]] (the MiRNA family's seed-extraction unit, MIRNA-SEED-001 — `MiRnaAnalyzer.GetSeedSequence`/`CreateMiRna`/`CompareSeedRegions`: string-level extraction of the 7-nt positions-2-8 seed, the normalised `MiRna` record, and exact-seed **family equality** via a Hamming comparison; the seed determines animal targeting and feeds the [[rna-base-pairing]] `GetReverseComplement` seed→target motif and the target-site predictor, distinct from base-pairing and from site-type classification), or
[[mirna-target-site-prediction]] (the MiRNA family's target-site-prediction unit, MIRNA-TARGET-001 — the **completing** unit — `MiRnaAnalyzer.FindTargetSites` two-pass antiparallel seed-complement scan classifying the Bartel/TargetScan hierarchy (8mer/7mer-m8/7mer-A1/6mer + offset 6mer, higher classes suppress overlapping offset-6mer) with a heuristic site score, plus the opt-in fully-source-traced but partial TargetScan **context++** regression scorer `ScoreTargetSiteContextPlusPlus` (per-site-type coefficients, min-max-scaled continuous + raw indicator features; computed Local_AU/3P_score/Min_dist/Len_3UTR/Off6m + `ComputeTa3Utr` TA=log10 N + McCaskill-partition SA + Friedman-Bls PCT; SPS/Len_ORF/ORF8m/PCT-sigmoid caller-supplied); depends on [[seed-sequence-analysis]] (the seed determines targeting) and [[rna-base-pairing]] (`GetReverseComplement` + `AlignMiRnaToTarget` duplex)), or
[[consensus-from-alignment]] (the Motif-Analysis consensus-from-a-multiple-alignment unit, MOTIF-CONS-001 — `MotifFinder.CreateConsensusFromAlignment`, the **pure most-frequent** (plurality) column-wise consensus over equal-length aligned strings with a deterministic **alphabetical** tie-break (A<C<G<T) and **no threshold**; the `alternative_to` the assembly [[consensus-sequence]] (which uses a plurality threshold + tie→ambiguous `dumb_consensus` rule), and the motif-family sibling of the exact [[known-motif-search]]), or
[[overrepresented-kmer-discovery]] (the Motif-Discovery de-novo motif-discovery unit, MOTIF-DISCOVER-001 — `MotifFinder.DiscoverMotifs`, enumerate-count-rank of *unknown* over-represented k-mers by the observed/expected enrichment `Count / ((N−k+1)/4^k)` under a zero-order uniform background; `alternative_to` the [[known-motif-search]] matcher — discovery outputs motifs, known-search takes them as input — and distinct from the [[consensus-from-alignment]] plurality consensus), or
[[shared-motifs]] (the Motif-Discovery cross-sequence shared-motif unit, MOTIF-SHARED-001 — `FindSharedMotifs`, the van Helden / RSAT oligo-analysis **"matching sequences"** quorum: enumerate every fixed-`k` exact word across a *set* of sequences and report each word present in ≥ `minSequences` of them, each carrying its `SequenceIndices` set + `Prevalence`; the cross-sequence sibling of the single-sequence O/E [[overrepresented-kmer-discovery]] and `alternative_to` the fixed-k-quorum-vs-single-longest [[longest-common-substring]] LCSM framing), or
[[clinical-actionability-oncokb-levels]] (the anchor for the Oncology family, the first ingested ONCO-* unit — `Clinical Actionability Assessment` by the OncoKB Therapeutic Levels of Evidence: a **pure level-ranking** of caller-supplied leveled drug associations under the fixed combined order `R1 > 1 > 2 > 3A > 3B > 4 > R2` (with separate sensitive `1 > 2 > 3A > 3B > 4` and resistance `R1 > R2` axes), returning the maximum level or `NotActionable` when a variant carries no leveled association; the knowledgebase is a caller input, not embedded), or
[[cancer-variant-tier-classification-amp-asco-cap]] (the second ingested ONCO-* unit — `Cancer-Specific Variant Annotation` by the **AMP/ASCO/CAP 2017 four-tier** clinical-significance classification: a decision rule over caller-supplied evidence level (A–D) + population MAF + cancer-association flag → Tier I (strong, Level A/B) / II (potential, Level C/D) / III (unknown, rare + assoc) / IV (benign, MAF ≥ 1% or no assoc), with `GetCOSMICAnnotation` a null-on-miss lookup against a caller-supplied COSMIC catalog; distinct from but consistent with the OncoKB therapeutic-level ranking [[clinical-actionability-oncokb-levels]]), or
[[sequencing-artifact-detection]] (the third ingested ONCO-* unit — `FilterArtifacts`, rule-based **technical-artifact** filtering that removes false-positive somatic calls before clinical interpretation: OxoG oxidation (G>T/C>A) vs FFPE cytosine-deamination (C>T/G>A) substitution-class classification, the Chen 2017 / Damage-estimator **GIV** read-pair imbalance ratio (`R1 G>T / R2 G>T`, neutral 1 / damaged > 1.5), and the GATK **FisherStrand** Phred-scaled two-sided Fisher-exact strand-bias FS over the `[ref_fwd, ref_rev, alt_fwd, alt_rev]` 2×2 table; the QC sibling of the two clinical-significance ONCO units), or
[[allele-specific-copy-number-ascat]] (the fourth ingested ONCO-* unit — the upstream tumor copy-number layer: allele-specific integer copy number nA/nB derived from per-locus **logR/BAF** by a joint **purity ρ / ploidy ψ** grid fit (**ASCAT**, Van Loo 2010, integer-closeness objective with BAF=0.5 down-weighted ×0.05, γ=1 for sequencing), the **ASPCF** penalised-least-squares joint segmentation front-end (Nilsen 2012 DP + Ross 2021 BAF-mirroring), a **subclonal** two-state bracketing-integer mixture (Battenberg / Nik-Zainal 2012), and the downstream **multiplicity / CCF** VAF-inversion (McGranahan 2016 / PICTograph / DeCiFering); the allele-specific/clonal counterpart to the total-CN [[aneuploidy-detection]] and the copy-number substrate beneath the clinical-interpretation ONCO siblings), or
[[cancer-cell-fraction-clonal-clustering]] (the fifth ingested ONCO-* unit — the downstream clonal-structure layer: the standalone **`EstimateCCF`** point estimator `CCF = f·(ρ·N_T + 2(1−ρ))/(ρ·m)` (McGranahan 2016 / PICTograph / Tarabichi 2021, corroborated three ways) with a reported-value **[0,1] cap** exposing the uncapped raw (CNAqc's 1.06 noise case), plus the genuinely distinct **`ClusterCCFValues`** — a **deterministic 1D Lloyd k-means** (Lloyd 1982, quantile seeding, no RNG) that deconvolutes a CCF vector into clones/subclones with the **highest-centroid = clonal** rule (Tarabichi 2021); consumes the purity ρ / copy number N_T / multiplicity m produced by the upstream [[allele-specific-copy-number-ascat]], whose §4 carries the same CCF closed form), or
[[clonal-hematopoiesis-cfdna-filtering]] (the sixth ingested ONCO-* unit — the pre-interpretation **biological-origin filter** of the liquid-biopsy pipeline: removes **clonal-hematopoiesis (CHIP)** false positives — the dominant cfDNA confounder (Razavi 2019: 81.6% controls / 53.2% cancer) — via `IdentifyCHIPVariants` (driver-gene {DNMT3A,TET2,ASXL1,TP53,JAK2,SF3B1,SRSF2,PPM1D} + VAF ≥ 2%, Steensma 2015 / Genovese 2014), `FilterCHIP` (matched-WBC subtraction as the definitive origin test + a conservative labelled gene+VAF heuristic fallback), and the strict Bolton-2020 `CallVariantOrigin` (WBC VAF ≥ 2% AND ≥ 10 reads AND ≥ φ× tumour VAF, φ=2.0 / 1.5 lymph node); the **biological-origin sibling** of the technical-artifact QC filter [[sequencing-artifact-detection]], upstream of every clinical-interpretation and clonal-structure ONCO unit), or
[[clonal-subclonal-classification-ccf-posterior]] (the seventh ingested ONCO-* unit — the **probabilistic** clonal-structure classifier: a per-mutation **Bayesian CCF posterior** `P(c) ∝ Binom(a|N,f(c))` on a 100-point grid `c∈[0.01,1]` with a uniform prior, expected allele fraction `f(c) = αMc/(2(1−α)+αq)` (Landau 2013 ABSOLUTE-style; multiplicity M generalisation per Satas 2021 DeCiFering Eq. 1), classified **clonal iff `P(CCF>0.95) > 0.5`** — the two thresholds distinct (0.95 CCF cut, 0.5 posterior-mass cut) — plus the point-estimate `IdentifyClonalMutations` (strict CCF > 0.95, boundary excluded → indices {0,2,4}); the **`alternative_to`** the point-estimate + Lloyd-k-means [[cancer-cell-fraction-clonal-clustering]] answering the same clonal/subclonal question without clustering, consuming the purity α / local copy number q / multiplicity M produced by the upstream [[allele-specific-copy-number-ascat]]; one API-shape assumption (per-variant local copy number q over a genome-wide ploidy scalar), no source contradictions), or
[[copy-number-alteration-classification]] (the eighth ingested ONCO-* unit — the **total-copy-number classification** layer: a single **log2 copy ratio** → absolute integer CN `n = 2·2^log2` (CNVkit `_log2_ratio_to_absolute_pure`, diploid ref_copies=2) → discrete **CNA state** via CNVkit's `absolute_threshold` hard-threshold caller (default `−1.1/−0.25/0.2/0.7` → DeepDeletion/Loss/Neutral/Gain/Amplification, boundary-inclusive `log2 <= thresh` → lower bin, above-last-threshold `ceil(2·2^log2)` open-ended amplification, NaN → neutral CN 2), corroborated by **GISTIC2** (Mermel 2011: ±0.1 noise band + +0.848/−0.737 high-amplitude cutoffs); the **total-CN** counterpart to the allele-specific [[allele-specific-copy-number-ascat]] and the per-segment oncology-state sibling of the whole-chromosome [[aneuploidy-detection]], both of which share the `n = 2·2^log2` conversion; one diploid-ploidy=2 assumption, no source contradictions).
