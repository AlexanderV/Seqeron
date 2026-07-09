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
  - docs/Evidence/ONCO-CNA-002-Evidence.md
  - docs/Evidence/ONCO-CNA-003-Evidence.md
  - docs/Evidence/ONCO-CTDNA-001-Evidence.md
  - docs/Evidence/ONCO-DRIVER-001-Evidence.md
  - docs/Evidence/ONCO-EXPR-001-Evidence.md
  - docs/Evidence/ONCO-FUSION-001-Evidence.md
  - docs/Evidence/ONCO-FUSION-002-Evidence.md
  - docs/Evidence/ONCO-FUSION-003-Evidence.md
  - docs/Evidence/ONCO-HETERO-001-Evidence.md
  - docs/Evidence/ONCO-HLA-001-Evidence.md
  - docs/Evidence/ONCO-HRD-001-Evidence.md
  - docs/Evidence/ONCO-IMMUNE-001-Evidence.md
  - docs/Evidence/ONCO-LOH-001-Evidence.md
  - docs/Evidence/ONCO-MHC-001-Evidence.md
  - docs/Evidence/ONCO-MRD-001-Evidence.md
  - docs/Evidence/ONCO-MSI-001-Evidence.md
  - docs/Evidence/ONCO-NEO-001-Evidence.md
source_commit: 643a974d3e132ca5be1beedd823ade8c4e535528
created: 2026-07-09
updated: 2026-07-10
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
[[onco-clonal-001-evidence]], [[onco-cna-001-evidence]],
[[onco-cna-002-evidence]], [[onco-cna-003-evidence]],
[[onco-ctdna-001-evidence]], [[onco-driver-001-evidence]],
[[onco-expr-001-evidence]], [[onco-fusion-001-evidence]],
[[onco-fusion-002-evidence]], [[onco-fusion-003-evidence]],
[[onco-hetero-001-evidence]], [[onco-hla-001-evidence]],
[[onco-hrd-001-evidence]], [[onco-immune-001-evidence]],
[[onco-loh-001-evidence]], [[onco-mhc-001-evidence]],
[[onco-mrd-001-evidence]], [[onco-msi-001-evidence]],
[[onco-neo-001-evidence]]). An
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
[[copy-number-alteration-classification]] (the eighth ingested ONCO-* unit — the **total-copy-number classification** layer: a single **log2 copy ratio** → absolute integer CN `n = 2·2^log2` (CNVkit `_log2_ratio_to_absolute_pure`, diploid ref_copies=2) → discrete **CNA state** via CNVkit's `absolute_threshold` hard-threshold caller (default `−1.1/−0.25/0.2/0.7` → DeepDeletion/Loss/Neutral/Gain/Amplification, boundary-inclusive `log2 <= thresh` → lower bin, above-last-threshold `ceil(2·2^log2)` open-ended amplification, NaN → neutral CN 2), corroborated by **GISTIC2** (Mermel 2011: ±0.1 noise band + +0.848/−0.737 high-amplitude cutoffs); the **total-CN** counterpart to the allele-specific [[allele-specific-copy-number-ascat]] and the per-segment oncology-state sibling of the whole-chromosome [[aneuploidy-detection]], both of which share the `n = 2·2^log2` conversion; one diploid-ploidy=2 assumption, no source contradictions), or
[[focal-amplification-detection]] (the ninth ingested ONCO-* unit — **focal amplification detection**: a two-part predicate `DetectFocalAmplifications` keeps segments that are both **amplified** (log2 gain > GISTIC2 `t_amp` 0.1) **and** **focal** (segment length < `broad_len_cutoff` 0.98 of its chromosome arm — Mermel 2011's length-based, not amplitude-based, focal/arm-level split; strict `< 0.98`, exactly-0.98 → arm-level), then `IdentifyAmplifiedOncogenes` maps each focal amplification's arm prefix to a built-in oncogene panel (17q→ERBB2 · 8q→MYC · 7p→EGFR · 11q→CCND1 · 12q→MDM2 **and** CDK4, NCBI Gene cytobands); orthogonal to the log2→5-state [[copy-number-alteration-classification]] (ONCO-CNA-001) — a **length** question, not an amplitude bin — sharing only the GISTIC2 `t_amp`=0.1 amplitude gate; caller supplies each segment's arm label + arm length (no bundled cytoband table); deletions out of scope (ONCO-CNA-003); no source contradictions), or
[[homozygous-deletion-detection]] (the tenth ingested ONCO-* unit — **homozygous / deep deletion detection**, the deletion mirror of ONCO-CNA-002: filter segments whose classified integer copy number is **exactly 0** (a homozygous / deep deletion — cBioPortal `−2` "Deep Deletion", Cheng et al. 2017 "zero copies of both alleles"; **reuses** ONCO-CNA-001's CN-0 = `DeepDeletion` call, **no new threshold**), then `IdentifyDeletedTumorSuppressors` maps each arm prefix to a built-in tumour-suppressor panel (17p→TP53 · 13q→RB1 **and** BRCA2 · 9p→CDKN2A · 10q→PTEN · 17q→BRCA1, NCBI Gene cytobands); a single-copy loss (−1, heterozygous, CN ≥ 1) is NOT homozygous; the loss-side counterpart of the amplification unit [[focal-amplification-detection]] and a consumer of [[copy-number-alteration-classification]]'s CN-0 state; two assumptions (CN-0-reuse, caller-fixed panel); no source contradictions), or
[[ctdna-detection-and-tumor-fraction]] (the eleventh ingested ONCO-* unit — the liquid-biopsy **quantification / limit-of-detection** layer: the Poisson ctDNA detection probability **`DetectionProbability` p = 1 − e^(−n·d·k)** (n genome equivalents, d mutant allele fraction, k reporters; Patent US 11,085,084 restating Avanzini 2020, corroborated by Pessoa 2023's λ = n·d = 15) with a detectability test (caller threshold default 0.95 **and** λ = n·d·k ≥ 1), **`CalculateTumorFraction` = 2 × mean clonal-heterozygous VAF** (copy-neutral diploid, CNAqc/Antonello 2024 v = π/2), `CalculateMeanVaf` across reporters (Newman 2014 CAPP-Seq), and a genome-equivalents helper (3.3 pg/haploid, Devonshire 2014 ⇒ ≈303 GE/ng, Alcaide 2020); the **quantification counterpart** on the same cfDNA input as the biological-origin filter [[clonal-hematopoiesis-cfdna-filtering]], sharing the diploid-heterozygous VAF inversion with [[allele-specific-copy-number-ascat]]/[[cancer-cell-fraction-clonal-clustering]]; one flagged detection-threshold assumption, no source contradictions), or
[[driver-gene-classification-20-20-rule]] (the twelfth ingested ONCO-* unit — the **gene-level** driver classifier: the Vogelstein 2013 **20/20 rule** (`IdentifyDriverMutations`/`MatchCancerHotspots`/`ScoreDriverPotential`) classifying a cancer gene **Oncogene** (> 20% missense at **recurrent** positions, ≥ 2× = activating), **TumorSuppressor** (> 20% **truncating/inactivating** — nonsense/frameshift/splice/stop-gain-loss), or **Ambiguous** (neither, or exact dual-pass tie); IDH1 codon-132 all-recurrent-missense → Oncogene (1.00), dispersed 8/10-truncating → TSG (0.80), exactly-0.20 → NOT TSG (strict `>`) oracles; `ScoreDriverPotential` returns max-of-two-fractions in [0,1] (CADD/SIFT/PolyPhen are caller-supplied, not implemented); a **heuristic not a statistical test** (20/20+/MutSigCV successors out of scope), orthogonal to the variant-level [[cancer-variant-tier-classification-amp-asco-cap]] / [[clinical-actionability-oncokb-levels]]; three assumptions (strict `>0.20`, dual-pass tie-break, ScoreDriverPotential proxy), no source contradictions), or
[[expression-outlier-zscore-signature-score]] (the thirteenth ingested ONCO-* unit, ONCO-EXPR-001, and the wiki's first expression/transcriptome method — the per-gene **outlier z-score** `z = (r−μ)/σ` against a caller-supplied reference cohort with **σ = sample SD (n−1)** per cBioPortal `NormalizeExpressionLevels.java`, **strict ±2** over/under-expression thresholds, plus the Lee et al. 2008 / GSVA **combined-z signature activity** `a = (Σzᵢ)/√k`; zero-SD reference **throws** (reference-code behaviour, a deviation from the prose "NA" spec); reference `{2,2,4,6,6}`→σ=2, x=10→3.0/x=8→2.0-not-outlier/x=−1→−2.5 and signature `{3,1,−1,1}`→2.0 oracles; orthogonal to the copy-number/clonal ONCO units; one behavioural deviation (throw not NA), two scope assumptions), or
[[gene-fusion-detection-read-evidence]] (the fourteenth ingested ONCO-* unit, ONCO-FUSION-001, and the wiki's first gene-fusion / read-evidence structural-rearrangement method — the STAR-Fusion / Arriba **split-read + discordant-pair + minimum-support** candidate-fusion caller: DETECTED iff (junction ≥ `MIN_JUNCTION_READS`=1 AND total ≥ `MIN_SUM_FRAGS`=2) OR (zero junction AND discordant ≥ `MIN_SPANNING_FRAGS_ONLY`=5), with **total support = split_reads1+split_reads2+discordant_mates** and the **gene5p ≠ gene3p** distinct-gene invariant, results ordered by descending support; plus the standalone **exon-phase in-frame** check `(5' coding bases − 3' start phase) mod 3 == 0`; EML4-ALK/CD74-ROS1 DETECTED · NCOA4-RET(0,0,4) / KIF5B-RET(sum<2) / ALK-ALK(same gene) REJECTED oracles; consumes already-grouped candidate breakpoint counts (Arriba schema, not raw BAM) and does not scan for premature stops (ONCO-FUSION-003); orthogonal to the copy-number/clonal and clinical-interpretation ONCO units, distinct from the gene-order [[genome-rearrangement-breakpoint-distance]]; two scope assumptions, no source contradictions), or
[[gene-fusion-nomenclature-known-fusion-lookup]] (the fifteenth ingested ONCO-* unit, ONCO-FUSION-002, and the wiki's fusion **annotation/naming** method — the **HGNC gene-fusion designation** `GetFusionAnnotation(5p,3p) = "5p::3p"` (Bruford et al. 2021: `::` double-colon separator, **5′-partner-first** directional order, HGNC approved symbols, read-throughs keep the hyphen `INS-IGF2`) plus the **directional known-fusion match** `MatchKnownFusions` against a **caller-supplied** set keyed by `5′::3′`, case-insensitive; a **Framework** algorithm — format/keying source-backed (Bruford 2021 + consortium extension), set contents caller-supplied (no bundled Mitelman/COSMIC/ChimerDB); `BCR::ABL1` worked example, `A::B ≠ B::A` reciprocal-direction / hyphen ≠ `::` corner cases; the naming layer downstream of the read-evidence caller [[gene-fusion-detection-read-evidence]] (annotation round-trips through a `DetectFusions` `FusionCall`), distinct from the ONCO-FUSION-003 premature-stop scope; two assumptions (caller-supplied set, case-insensitivity), no source contradictions), or
[[fusion-breakpoint-frame-and-protein-prediction]] (the sixteenth ingested ONCO-* unit, ONCO-FUSION-003, and the **protein-consequence** third member of the fusion trio — the **four-state** `BreakpointFrameStatus` (`InFrame`/`OutOfFrame`/`StopCodon`/`NotPredicted`, Arriba's two-way native-frame model, **not** AGFusion's three-way `in-frame (with mutation)` class → that maps to `OutOfFrame`) gated by breakpoint-site classification (CDS vs UTR/intron/intergenic → `NotPredicted`), reusing ONCO-FUSION-001's exon-phase rule `(5' coding bases − 3' start phase) mod 3 == 0`, plus **`PredictFusionProtein`** (AGFusion `model.py`: 5′ CDS prefix + 3′ CDS suffix → concatenate → translate transl_table=1 → **truncate at first stop**; out-of-frame trims to whole codons first); oracles `ATGAAA|GATGGT`→`MKDG`, `ATGAAA|GATTAAGGT`→`MKD` (premature stop), `ATGA|AAGGT` phase-0→`OutOfFrame` yet clean `MKG` (the Arriba-vs-AGFusion divergence point) / phase-1→in-frame; downstream of the read-evidence caller [[gene-fusion-detection-read-evidence]] (elaborates its binary in-frame check into a full consequence with transcript/peptide reconstruction) and orthogonal to the naming unit [[gene-fusion-nomenclature-known-fusion-lookup]]; one API-shape assumption (caller supplies CDS strings + junction offsets, no bundled GTF), no source contradictions), or
[[intratumor-heterogeneity-metrics]] (the seventeenth ingested ONCO-* unit, ONCO-HETERO-001, and the wiki's first intratumor-heterogeneity **scalar-metric** method — the summary layer that reduces the reconstructed clonal structure to comparable numbers: the **MATH score** `100·1.4826·median(|VAF−median VAF|)/median(VAF)` (Mroz & Rocco 2013 / Mroz 2015 / maftools `mathScore.R`, three-way algebraically identical — the ratio of the width to the centre of the mutant-allele-fraction distribution, computed straight from VAFs with no clustering), the **Shannon clonal diversity** `H = −Σ pᵢ ln pᵢ` (natural log / nats, Liu 2017 / Shannon 1948), the **subclone count** (Liu 2017 richness = number of occupied CCF clusters), and the **subclonal fraction** `#(CCF < 0.95)/n` (reusing the Landau 2013 clonal CCF threshold 0.95); MATH oracles 49.42 (odd VAFs) / 59.304 (even), Shannon 0 / ln 2 / ln 4, zero-median-VAF → throw and MAD=0 → MATH=0 corner cases; **depends on** the CCF clustering of [[cancer-cell-fraction-clonal-clustering]] (its subclone count + Shannon pᵢ consume the CCF clusters) and reuses the 0.95 threshold of [[clonal-subclonal-classification-ccf-posterior]] — a **summary/metric** layer, not per-mutation reconstruction, sharing the Shannon-index math but not the domain with the metagenomics [[alpha-diversity]]; two source-consistent assumptions (Shannon pᵢ = per-cluster mutation proportions, R even-count median), no source contradictions), or
[[hla-nomenclature-and-allele-specific-loh]] (the eighteenth ingested ONCO-* unit, ONCO-HLA-001, and the wiki's first HLA / immuno-oncology antigen-presentation method — two disjoint pieces sharing the HLA locus: (1) **HLA allele nomenclature parsing/validation** of the WHO IPD-IMGT/HLA name `HLA-[Gene]*[F1]:[F2][:F3][:F4][suffix]` (Marsh 2010 colon-delimited convention; Field1 type/allele-group / Field2 specific protein / Field3 synonymous / Field4 non-coding, two-field "four-digit" minimum, four-field maximum, `HLA-` prefix required, N/L/S/C/A/Q expression-suffix set — `HLA-A*24:02:01:02L` valid, `HLA-A*02`/`A*02:01`/five-fields/`…X` rejected), and (2) **allele-specific HLA LOH** by the **LOHHLA** rule (McGranahan 2017, `mskcc/lohhla`): call HLA LOH iff min per-allele copy number **< 0.5** AND allelic-imbalance paired Student's t-test **p < 0.01** — both strict `<`, the p-test the explicit "avoid over-calling LOH" guard; oracles (1.8, 0.30, 0.001)→LOH allele 2 / (1.60, 0.40, 0.05)→no (p≥0.01) / (1.50, 0.50, 0.001)→no (0.5 not <0.5) / (1.70, 0.40, 0.01)→no (0.01 not <0.01); the HLA-locus specialization of the genome-wide allele-specific copy-number / LOH derivation [[allele-specific-copy-number-ascat]] (same logR + BAF allelic-contrast idea restricted to the HLA genes), a mechanism of immune escape via lost neoantigen presentation; one assumption (both alleles < 0.5 → distinct `HomozygousLoss` label, thresholds unchanged), no source contradictions; research-grade, not for clinical use), or
[[homologous-recombination-deficiency-score]] (the nineteenth ingested ONCO-* unit, ONCO-HRD-001 — the **HRD composite genomic-scar score** `HRD = LOH + TAI + LST`, an **unweighted sum** of three large-scale copy-number scar counts with the **HRD-high cutoff ≥ 42** inclusive (Telli 2016 + Stewart 2022 review, corroborated independently); all three components **derived per segment** from the [[allele-specific-copy-number-ascat]] major/minor CN substrate in `DetectHRD(segments)` — **HRD-LOH** = LOH regions > 15 Mb and < whole chromosome, excluding whole-chromosome-LOH chromosomes (Abkevich 2012 / oncoscanR / scarHRD `calc.hrd.R`, no centromere table needed); **TAI** = imbalanced (major ≠ minor CN) segments reaching a sub-telomere but not crossing the centromere (Birkbak 2012 / scarHRD `calc.ai_new.R`, sub-1 Mb dropped, single-segment whole-chr imbalance excluded); **LST** = adjacent ≥ 10 Mb same-arm regions separated by < 3 Mb after iterative 3 Mb smoothing (Popova 2012 / scarHRD `calc.lst.R`, autosomes only); TAI + LST need the per-chromosome centromere `[start,end]` **acen** table embedded for GRCh38/GRCh37 (UCSC cytoBand cross-verified vs NCBI GRC modeled centromeres — resolving the prior "centromere table unretrievable" blocker); (14,14,14)→42 HRD-high / (14,13,14)→41 negative / (0,0,0)→0 near-diploid oracles; the genomic-scar aggregation layer above [[allele-specific-copy-number-ascat]], distinct from the total-CN [[aneuploidy-detection]]; one even-ploidy AI-path assumption (major ≠ minor, ASCAT ploidy column absent), no source contradictions), or
[[immune-infiltration-deconvolution]] (the twentieth ingested ONCO-* unit, ONCO-IMMUNE-001, and the second expression/transcriptome method — **immune infiltration estimation** by expression deconvolution + enrichment scoring: the CIBERSORT linear-mixture model `m = S·f` solved by **ν-SVR** `DeconvoluteImmuneCellsNuSvr` (Newman 2015 / Schölkopf 2000; z-standardize, sweep ν∈{0.25,0.5,0.75} by lowest RMSE, zero-clip + normalize Σf=1; cross-checked vs scikit-learn `NuSVR`) with an NNLS/LLSR baseline `DeconvoluteImmuneCells` (Abbas 2009), plus **ESTIMATE** ssGSEA immune/stromal scoring `EstimateInfiltration` and the opt-in Affymetrix-only cosine tumor-purity transform `EstimateTumorPurity` (`cos(0.6049872018 + 0.0001467884·ESTIMATE)`, negative → NaN); **LM22** (547×22) is caller-supplied — Stanford forbids redistribution, no exact-CIBERSORT parity claim — while the permissive **ABIS-Seq** matrix (Monaco 2019, CC BY 4.0, 1296×17) is bundled via `LoadBundledAbisSignatureMatrix` so deconvolution works out-of-the-box; simplified-ssGSEA + Affymetrix-domain assumptions, collinear-cell-type looser-proportion corner case; the antigen-agnostic quantitative sibling of the immuno-oncology [[hla-nomenclature-and-allele-specific-loh]], sharing the ssGSEA signature-scoring layer with [[expression-outlier-zscore-signature-score]]; research-grade, not for clinical use), or
[[loss-of-heterozygosity-detection]] (the twenty-first ingested ONCO-* unit, ONCO-LOH-001 — the standalone genome-wide **LOH caller** `DetectLOH` + per-chromosome `CalculateLOHFraction`: counts **HRD-LOH** regions where **minor CN 0 AND major CN ≠ 0** (allelic loss, not homozygous deletion), length **strictly > 15 Mb** (length = end − start), excluding any **whole-chromosome** global-LOH chromosome (`chrDel`) — Abkevich 2012 + scarHRD `calc.hrd.R`/`scar_score.R` + oncoscanR `score_loh` (adjacent/overlapping LOH merged ≤ 1 bp first; major CN capped to 1 so LOH state drives merging); synthetic 7-segment dataset → score 1 (only the 20 Mb chr1 segment qualifies; exactly-15 Mb / het-retained / homozygous-deletion / whole-chr-LOH all excluded), LOH-fraction oracles 0.333/0.0/1.0; this **is** the LOH term the composite [[homologous-recombination-deficiency-score]] sums with TAI + LST (centromere-free, unlike its siblings), reads its segments off [[allele-specific-copy-number-ascat]], and is the genome-wide counterpart of the HLA-locus [[hla-nomenclature-and-allele-specific-loh]] (LOHHLA); two assumptions (segment input shape; the length-weighted LOH-fraction is a definitional/API choice), no source contradictions), or
[[mhc-peptide-binding-prediction]] (the twenty-second ingested ONCO-* unit, ONCO-MHC-001, and the wiki's first **MHC / peptide-presentation prediction** method — the affinity gate of neoantigen candidate scoring, with three layers over one `Strong`/`Weak`/`NonBinder` output space: (1) **classification** — IC50 tiers `< 50` Strong / `< 500` Weak / `≥ 500` NonBinder (Sette 1994 ≈ 500 nM CTL threshold + IEDB tiers) and NetMHCpan-4.1 %Rank tiers (class I SB `< 0.5` / WB `< 2`, class II `< 2` / `< 10`, all strict `<`; Reynisson 2020), plus peptide-length validity (class I **8–14**, class II 13–25); (2) **matrix prediction** — the BIMAS **product rule** (T½ = finalConstant · ∏ coefficients, 180 coeffs + const, unlisted residue = 1.0; Parker 1994) and the SMM **additive rule** with `IC50 = 50000^(1−score)` (Peters & Sette 2005 / IEDB log50k; score 0→50000, 0.5→223.6068, 1→1 nM); (3) a ported **MHCflurry 2.0** pan-allele class-I NN (BLOSUM62 `left_pad_centered_right_pad` 45×21 peptide + 37-residue allele pseudosequence 37×21, `tanh` Dense stack with `feedforward`/`with-skip-connections` topologies, `to_ic50 = 50000^(1−x)`, geometric-mean ensemble; port matches `mhcflurry` 2.1.5 to < 0.03 %); **Framework/packaging** boundary — trained matrices caller-supplied (BIMAS defunct CGI / Parker paywalled / IEDB non-commercial) and the 80 MB ensemble not embedded (one member for CI parity, rest via `LoadWeightPack`), explicitly analogized in-file to the CIBERSORT LM22 of [[immune-infiltration-deconvolution]]; the affinity gate downstream of neoantigen-peptide generation [[neoantigen-peptide-generation]] (ONCO-NEO-001) and the presentation-machinery sibling of [[hla-nomenclature-and-allele-specific-loh]] (HLA LOH removes the platform for the neoantigens scored against a lost allele); class I 8–14 window resolved vs NetMHCpan-4.1, no source contradictions), or
[[tumor-informed-mrd-detection]] (the twenty-third ingested ONCO-* unit, ONCO-MRD-001 — the liquid-biopsy **MRD calling / statistical-detection** layer built on top of the Poisson quantification primitive [[ctdna-detection-and-tumor-fraction]] (ONCO-CTDNA-001, which it `depends_on` and reuses for the panel Poisson `p = 1 − e^(−nfm)`): the tumor-informed **Signatera** workflow — up to 16 WES-selected clonal somatic SNVs tracked in plasma, `DetectMRD` **MRD-positive iff ≥ 2 of the 16 tracked variants are detected** (Reinert 2019 / PMC9265001 Table 1; 0/1→negative, 2/3→positive; sensitivity scales with panel size, >8 markers needed at ≤0.1% VAF, LoD 0.01% / specificity >99.5%; post-surgery ctDNA+ → HR 7.2 relapse) + longitudinal `TrackVariantsOverTime` (per-timepoint status + first-positive flag), plus the **INVAR** generalised-likelihood-ratio detector (Wan 2020 + INVAR2 verbatim formulas): per-locus mixture `q = p·g + e(1−p)`, EM ctDNA-fraction estimate, `LR = logL(p̂) − logL(0)` (pure background → p̂≈0/LR≈0/not detected; AF/SNR-weighted, LR monotone in signal), `IntegratedMutantAlleleFractionV2` background-subtracted depth-weighted `weighted.mean(max(0, MEAN_AF−BACKGROUND_AF), TOTAL_DP)`, fragment-size-weighted GLRT `EstimateInvarSignalWithSize`/`FragmentSizeProfile`/`InvarMolecule` (short tumour fragments up-weighted; discrete `COUNT/TOTAL` default, opt-in Gaussian-KDE `FromKernelDensity` via Silverman `bw.nrd0` + analytic Φ bin integral; flat P1==P0 → reduces to no-size GLRT), patient-specific outlier suppression `SuppressOutlierLoci` (repolish one-sided binomial + Bonferroni) and control-derived background `EstimateLocusBackground`/`PassesBothStrandsFilter`; **distinct** from ONCO-CTDNA-001 (multi-variant verdict vs single-reporter probability); one flagged per-variant "detected" = ≥1 alt read assumption (tunable, panel ≥2 rule unaffected), no source contradictions), or
[[microsatellite-instability-detection]] (the twenty-fourth ingested ONCO-* unit, ONCO-MSI-001 — the **microsatellite-instability biomarker** unit: the **MSIsensor/MSIsensor2** continuous MSI **score** `unstable loci / valid loci` (as %) with **MSI-High cutoff ≥ 20%** inclusive (Niu-lab README; per-site chi-square tumor-vs-normal repeat-length test at FDR 0.05, "percentage of sites with a somatic indel", Niu 2014; tumor-only mode unchanged; the original cohort's 3.5% separation is dataset-specific not the production cutoff), **and** the categorical **Bethesda** marker-count rule over the validated 5-marker panel (BAT-25/BAT-26/D2S123/D5S346/D17S250) 0/5→MSS · 1/5→MSI-L · ≥2/5→MSI-H (Boland 1998, cross-checked vs the 2004 revised-Bethesda ≥40%/≥20%&<40% fraction form BJC 2014); the **key modelling choice** = two inputs two classifiers — the continuous fraction gets only the binary MSI-H (≥20%) vs not-High, the three-way MSS/MSI-L/MSI-H applies to the discrete marker count (no MSI-L band fabricated for the fraction); 5/25→20%→MSI-H / 4/25→16%→not-High / 0-valid→undefined-throws / unstable>valid→throws oracles; an independent immunotherapy/mismatch-repair biomarker orthogonal to the copy-number-scar [[homologous-recombination-deficiency-score]] and consumed by the clinical-interpretation units [[cancer-variant-tier-classification-amp-asco-cap]] / [[clinical-actionability-oncokb-levels]]; one assumption (no MSI-L band on the continuous score), no source contradictions), or
[[neoantigen-peptide-generation]] (the twenty-fifth ingested ONCO-* unit, ONCO-NEO-001, and the **upstream partner** of the affinity gate [[mhc-peptide-binding-prediction]] — **neoantigen candidate peptide window generation** `GenerateNeoantigenPeptides`: given a somatic **missense SNV** at 1-based position `p` in a protein of length `L`, enumerate **every length-`k` window of the mutant protein that spans the mutated residue** (valid 0-based start `s ∈ [max(0, p−1−k+1), min(p−1, L−k)]`, count `= last − first + 1`; interior mutation → exactly `k` windows per length), each paired with its matched **wild-type (germline)** peptide at the same coordinates + the mutation offset — the mutant/WT **agretope** on which the differential agretopic index (DAI) is computed (Wells 2020, TESLA); class I search lengths **8–11** (pVACtools Hundal 2020; a 21-mer with the substitution centred, ProGeo-neo Li 2020's ±10-flank construction, contains exactly every overlapping 8–11-mer; NetMHCpan 8–14 / 9-mer-dominant Jurtz 2017 motivates enumerating 8–11 not 9 alone); interior **Y5C** on `MKTAYIAKQRSTVWLNDEFGH` → 5 windows per k = 8/9/10/11 → **20** candidates, each WT = mutant with `C`→`Y` at the offset; **terminal M1V** k = 9 → exactly 1 truncated window `VKTAYIAKQ`/`MKTAYIAKQ` offset 0; corner cases non-coding/synonymous → no peptide, length outside 8–11 not enumerated, protein shorter than `k` skips that length; **missense-only** (frameshift/indel/fusion neopeptides out of scope, separate translation) and **binding affinity out of scope** (IC50/%Rank scoring is caller-supplied via the downstream [[mhc-peptide-binding-prediction]] ONCO-MHC-001 gate, whose resolved class I 8–14 window propagates to this generator's default), the presentation-platform sibling [[hla-nomenclature-and-allele-specific-loh]] (HLA LOH removes the allele a neoantigen is restricted to); two source-backed scoping assumptions (single-residue substitution, affinity out of scope), no source contradictions).
