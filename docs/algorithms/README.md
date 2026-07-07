# Algorithms Documentation Index

This directory contains algorithm and method documentation grouped by area. Each link points to a
folder of per-algorithm documents. Identifier and folder canonicalization is defined in the
[Canonical Algorithm Map](CANONICAL_MAP.md).

## Sequence composition & statistics

- [Sequence Composition](Sequence_Composition) — GC skew, linguistic complexity, Shannon entropy, validation.
- [Statistics](Statistics) — descriptive sequence statistics and summaries.
- [Complexity](Complexity) — sequence-complexity measures.
- [Quality](Quality) — read/sequence quality metrics.

## Pattern matching, k-mers & repeats

- [Pattern Matching](Pattern_Matching) — exact/approximate search, edit distance, IUPAC, PWM, suffix tree.
- [K-mer](K-mer) — k-mer counting, frequency analysis, search.
- [Repeat Analysis](Repeat_Analysis) — direct/inverted/tandem repeats, microsatellites, palindromes.
- [Motif Discovery](Motif_Discovery) — motif finding and scoring.

## Alignment & comparative genomics

- [Alignment](Alignment) — pairwise and multiple sequence alignment.
- [Comparative Genomics](Comparative_Genomics) — ANI, orthology, synteny, rearrangements, dot plots.
- [Pan-genome](PanGenome) — core/accessory gene analysis.

## Annotation, variants & transcriptome

- [Annotation](Annotation) — gene prediction, ORF detection, promoter detection, GFF3 I/O.
- [Variants](Variants) — variant calling, annotation, SNP/indel handling.
- [Structural Variants](StructuralVar) — SV and CNV detection, breakpoints.
- [Splicing](Splicing) — splice-site and isoform analysis.
- [Transcriptome](Transcriptome) — expression and differential analysis.
- [Translation](Translation) — six-frame translation and related tools.

## Assembly

- [Assembly](Assembly) — OLC/de Bruijn assembly, consensus, scaffolding, correction, trimming, stats.

## Phylogenetics & population genetics

- [Phylogenetics](Phylogenetics) — distances, tree building (UPGMA/NJ), tree statistics, bootstrap.
- [Population Genetics](Population_Genetics) — diversity, differentiation, Hardy–Weinberg.

## Metagenomics

- [Metagenomics](Metagenomics) — taxonomic classification (Kraken-style), profiling, diversity, function.

## RNA structure & non-coding RNA

- [RNA Structure](RnaStructure) — secondary-structure prediction.
- [miRNA](MiRNA) — miRNA precursor, seed, and target analysis.

## Protein analysis

- [Protein Motifs](ProteinMotif) — protein motif detection.
- [Protein Prediction](ProteinPred) — disorder and property prediction.

## Epigenetics

- [Epigenetics](Epigenetics) — methylation, CpG islands, DMRs, bisulfite, chromatin, age.

## Oncology

- [Oncology](Oncology) — copy-number, drivers, clonality (CCF), ctDNA, mutational signatures, and more.

## Chromosome analysis

- [Chromosome Analysis](Chromosome_Analysis) — karyotype, centromere/telomere, aneuploidy, synteny.

## Codon usage & optimization

- [Codon](Codon) — codon usage, CAI, RSCU, rare-codon analysis.
- [Codon Optimization](Codon_Optimization) — codon optimization strategies.

## Molecular tools

- [MolTools](MolTools) — primer/probe design, CRISPR guide & PAM design, on-/off-target scoring,
  restriction analysis, melting temperature.

## General genomic analysis

- [Analysis](Analysis) — general genomic-analysis utilities.
- [GC-Skew Analysis](Extended_GC_Skew_Analysis) — extended GC-skew tooling.

## File I/O

- [File I/O](FileIO) — format parsers/writers (FASTA, FASTQ, GenBank, GFF, VCF, BED, EMBL).

---

## Canonicalization

One concept has exactly one canonical ID and one canonical documentation bucket; aliases remain
searchable but point to the canonical owner. See the full [Canonical Algorithm Map](CANONICAL_MAP.md).

### Current ID aliases

| Alias ID | Canonical ID |
|---|---|
| `SEQ-COMPOSITION-001` | `SEQ-STATS-001` |
| `SEQ-TM-001` | `SEQ-THERMO-001` |
| `GENOMIC-TANDEM-001` | `REP-TANDEM-001` |

### Legacy / alias folders

Consolidation status by folder. Names shown **without a link** have already been merged into their
canonical bucket (the folder is gone and all references were repointed); names shown **as links**
still exist and are pending consolidation.

| Legacy folder | Canonical bucket |
|---|---|
| `Molecular_Tools` | [MolTools](MolTools) |
| `K-mer_Analysis` | [K-mer](K-mer) |
| `PopGen` | [Population Genetics](Population_Genetics) |
| `RNA_Structure`, `RNA_Secondary_Structure` | [RnaStructure](RnaStructure) |
| [`Motif_Analysis`](Motif_Analysis) | [Motif Discovery](Motif_Discovery) |
| [`Sequence_Comparison`](Sequence_Comparison) | [Comparative Genomics](Comparative_Genomics) |
| [`Genomic_Analysis`](Genomic_Analysis) | [Analysis](Analysis) |
| [`Extended_Annotation`](Extended_Annotation) | [Annotation](Annotation) |
| [`Extended_Assembly`](Extended_Assembly) | [Assembly](Assembly) |
