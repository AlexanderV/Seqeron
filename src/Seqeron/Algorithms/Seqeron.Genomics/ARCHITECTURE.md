# Seqeron.Genomics - Modular Architecture

## Package Hierarchy

```
Level 0 (No dependencies):
└── Seqeron.Genomics.Infrastructure     # Base types, statistics, constants

Level 1 (Infrastructure only):
└── Seqeron.Genomics.Core               # Sequence models (DNA, RNA, Protein)

Level 2 (Core + Infrastructure):
├── Seqeron.Genomics.IO                 # Format parsers
├── Seqeron.Genomics.Alignment          # Sequence alignment
└── Seqeron.Genomics.Analysis           # K-mer, motif, repeat analysis

Level 3 (Multiple dependencies):
├── Seqeron.Genomics.Annotation         # Genome annotation
├── Seqeron.Genomics.Phylogenetics      # Phylogenetic analysis
├── Seqeron.Genomics.Population         # Population genetics
├── Seqeron.Genomics.Metagenomics       # Metagenomic analysis
├── Seqeron.Genomics.MolTools           # Molecular tools
└── Seqeron.Genomics.Chromosome         # Chromosome analysis

Level 4 (Top-level):
├── Seqeron.Genomics.Reports            # Report generation
└── Seqeron.Genomics                    # Meta-package (aggregator)
```

## Packages

### Seqeron.Genomics.Infrastructure
Base types with no external dependencies.

| File | Description |
|------|-------------|
| `StatisticsHelper.cs` | Statistical calculations (mean, std, correlation) |
| `ThermoConstants.cs` | Thermodynamic constants for Tm calculation |
| `AlignmentTypes.cs` | Scoring matrices, alignment result types |

### Seqeron.Genomics.Core
Fundamental sequence models.

| File | Description |
|------|-------------|
| `ISequence.cs` | Base sequence interface |
| `DnaSequence.cs` | DNA sequence with suffix tree |
| `RnaSequence.cs` | RNA sequence |
| `ProteinSequence.cs` | Protein sequence |
| `SequenceExtensions.cs` | Extension methods |
| `GeneticCode.cs` | Codon tables |
| `Translator.cs` | DNA/RNA to protein translation |
| `IupacHelper.cs` | IUPAC ambiguity codes |

### Seqeron.Genomics.IO
Format parsers and writers.

| File | Description |
|------|-------------|
| `FastaParser.cs` | FASTA format |
| `FastqParser.cs` | FASTQ format with quality |
| `GenBankParser.cs` | GenBank format |
| `EmblParser.cs` | EMBL format |
| `GffParser.cs` | GFF3 format |
| `BedParser.cs` | BED format |
| `VcfParser.cs` | VCF variant format |
| `SequenceIO.cs` | Unified I/O facade |
| `QualityScoreAnalyzer.cs` | Quality score analysis |
| `FeatureLocationHelper.cs` | Feature location parsing |

### Seqeron.Genomics.Alignment
Sequence alignment algorithms.

| File | Description |
|------|-------------|
| `SequenceAligner.cs` | Global/local/semi-global alignment |
| `ApproximateMatcher.cs` | Fuzzy pattern matching |
| `SequenceAssembler.cs` | Sequence assembly |

### Seqeron.Genomics.Analysis
Sequence analysis algorithms.

| File | Description |
|------|-------------|
| `KmerAnalyzer.cs` | K-mer counting and analysis |
| `SequenceStatistics.cs` | Sequence statistics |
| `SequenceComplexity.cs` | Complexity measures |
| `GcSkewCalculator.cs` | GC skew analysis |
| `GenomicAnalyzer.cs` | General genomic analysis |
| `MotifFinder.cs` | DNA motif discovery |
| `ProteinMotifFinder.cs` | Protein motif discovery |
| `RepeatFinder.cs` | Repeat detection |
| `RnaSecondaryStructure.cs` | RNA structure prediction |
| `DisorderPredictor.cs` | Protein disorder prediction |
| `ComparativeGenomics.cs` | Comparative analysis |

### Seqeron.Genomics.Annotation
Genome annotation and variant analysis.

| File | Description |
|------|-------------|
| `GenomeAnnotator.cs` | Gene prediction and annotation |
| `VariantAnnotator.cs` | Variant effect prediction |
| `VariantCaller.cs` | Variant detection |
| `SpliceSitePredictor.cs` | Splice site prediction |
| `TranscriptomeAnalyzer.cs` | Transcriptome analysis |
| `StructuralVariantAnalyzer.cs` | SV detection |
| `EpigeneticsAnalyzer.cs` | Epigenetic analysis |
| `MiRnaAnalyzer.cs` | miRNA analysis |

### Seqeron.Genomics.Phylogenetics
Phylogenetic analysis.

| File | Description |
|------|-------------|
| `PhylogeneticAnalyzer.cs` | Tree construction, distance matrices |

### Seqeron.Genomics.Population
Population genetics.

| File | Description |
|------|-------------|
| `PopulationGeneticsAnalyzer.cs` | Hardy-Weinberg, Fst, diversity |

### Seqeron.Genomics.Metagenomics
Metagenomic analysis.

| File | Description |
|------|-------------|
| `MetagenomicsAnalyzer.cs` | Taxonomic classification |
| `PanGenomeAnalyzer.cs` | Pan-genome analysis |

### Seqeron.Genomics.MolTools
Molecular biology tools.

| File | Description |
|------|-------------|
| `PrimerDesigner.cs` | PCR primer design |
| `ProbeDesigner.cs` | Hybridization probe design |
| `CrisprDesigner.cs` | CRISPR guide RNA design |
| `CodonOptimizer.cs` | Codon optimization |
| `CodonUsageAnalyzer.cs` | Codon usage analysis |
| `RestrictionAnalyzer.cs` | Restriction site analysis |

### Seqeron.Genomics.Chromosome
Chromosome-level analysis.

| File | Description |
|------|-------------|
| `ChromosomeAnalyzer.cs` | Karyotyping, centromere/telomere |
| `GenomeAssemblyAnalyzer.cs` | Assembly quality (N50, etc.) |

### Seqeron.Genomics.Reports
Report generation.

| File | Description |
|------|-------------|
| `ReportGenerator.cs` | HTML/JSON/Markdown reports |

### Seqeron.Genomics (Meta-package)
Aggregates all packages. Provides high-level coordinating APIs.

| File | Description |
|------|-------------|

## Usage

### Single Package
```csharp
// Use only what you need
using Seqeron.Genomics.Core;
using Seqeron.Genomics.IO;

var sequences = FastaParser.Parse("genome.fasta");
```

### Full Package
```csharp
// Use everything via meta-package
using Seqeron.Genomics;

var dna = new DnaSequence("ATCGATCG");
var report = ReportGenerator.GenerateHtml(results);
```

## Design Principles

1. **Minimal dependencies** - each package depends only on what it needs
2. **File-scoped namespaces** - modern C# style
3. **Global usings via csproj** - no repetitive using statements
4. **Immutable results** - readonly record struct for results
5. **Span/Memory** - efficient memory handling
6. **CancellationToken** - all long operations support cancellation
