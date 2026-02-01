# Seqeron.Genomics - Implementation Plan

## Architectural Principles

### Clean Architecture
```
Seqeron.Genomics/
├── Core/                    # Domain models and interfaces
│   ├── Sequences/           # DnaSequence, RnaSequence, ProteinSequence
│   ├── Interfaces/          # ISequence, ISequenceAnalyzer, IAligner
│   └── Results/             # Result types for algorithms
├── Analysis/                # Analytical algorithms
│   ├── Repeats/             # Repeat search
│   ├── Motifs/              # Motif search
│   ├── Structure/           # Secondary structure
│   └── Statistics/          # Statistical analysis
├── Alignment/               # Sequence alignment
├── Comparison/              # Comparative genomics
├── Applications/            # Practical applications (primers, CRISPR)
└── IO/                      # Format parsers (FASTA, GenBank, etc.)
```

### Coding Standards
- **Immutable results**: all Result types - readonly struct
- **Fluent API**: for configuring algorithm parameters
- **Lazy evaluation**: IEnumerable for large results
- **Memory efficient**: Span/Memory for working with sequences
- **Cancellation support**: CancellationToken for long operations
- **Nullable annotations**: full support for nullable reference types

---

## Phase 1: Core Foundation (Week 1-2)

### 1.1 Basic Sequence Types

#### Tasks
- [ ] `ISequence` - base interface
- [ ] `RnaSequence` - RNA sequence (A, C, G, U)
- [ ] `ProteinSequence` - amino acid sequence (20 AA)
- [ ] `IupacDnaSequence` - DNA with IUPAC codes (N, R, Y, etc.)
- [ ] `QualitySequence` - sequence with quality (FASTQ)

#### Files
```
Core/Sequences/
├── ISequence.cs
├── SequenceBase.cs
├── DnaSequence.cs (refactor existing)
├── RnaSequence.cs
├── ProteinSequence.cs
├── IupacDnaSequence.cs
└── QualitySequence.cs
```

#### Tests (25 tests)
- [ ] RnaSequenceTests.cs (10 tests)
- [ ] ProteinSequenceTests.cs (10 tests)
- [ ] IupacDnaSequenceTests.cs (5 tests)

---

### 1.2 Genetic Code and Translation

#### Tasks
- [ ] `GeneticCode` - codon table (Standard, Mitochondrial, etc.)
- [ ] `Translator` - DNA/RNA → protein translation
- [ ] `CodonTable` - codon frequencies for organisms

#### Files
```
Core/Translation/
├── GeneticCode.cs
├── CodonTable.cs
├── Translator.cs
└── TranslationResult.cs
```

#### Tests (15 tests)
- [ ] GeneticCodeTests.cs (5 tests)
- [ ] TranslatorTests.cs (10 tests)

---

## Phase 2: DNA Analysis (Week 3-4)

### 2.1 Exact Pattern Matching (partially exists)

#### Tasks
- [ ] Refactor `FindMotif` → `PatternMatcher`
- [ ] Support IUPAC wildcard patterns
- [ ] Batch pattern matching (multiple search)

#### Files
```
Analysis/Motifs/
├── PatternMatcher.cs
├── IupacPatternMatcher.cs
├── PatternMatchResult.cs
└── BatchPatternMatcher.cs
```

#### Tests (20 tests)
- [ ] PatternMatcherTests.cs (10 tests)
- [ ] IupacPatternMatcherTests.cs (10 tests)

---

### 2.2 Approximate Matching (with mutations)

#### Tasks
- [ ] `ApproximateMatcher` - search with k errors
- [ ] Support different error types:
  - Substitutions
  - Insertions
  - Deletions
- [ ] Edit distance (Levenshtein)
- [ ] Hamming distance (substitutions only)

#### Files
```
Analysis/Motifs/
├── ApproximateMatcher.cs
├── EditDistance.cs
├── ApproximateMatchResult.cs
└── MismatchType.cs
```

#### Tests (25 tests)
- [ ] ApproximateMatcherTests.cs (15 tests)
- [ ] EditDistanceTests.cs (10 tests)

---

### 2.3 Repeat Analysis (extending existing)

#### Tasks
- [ ] `SupermaximalRepeatFinder` - maximal repeats
- [ ] `InvertedRepeatFinder` - inverted repeats
- [ ] `MicrosatelliteFinder` - STR (Short Tandem Repeats)
- [ ] `MinisatelliteFinder` - VNTR

#### Files
```
Analysis/Repeats/
├── RepeatFinder.cs (refactor existing)
├── SupermaximalRepeatFinder.cs
├── InvertedRepeatFinder.cs
├── MicrosatelliteFinder.cs
├── MinisatelliteFinder.cs
└── RepeatClassifier.cs
```

#### Tests (30 tests)
- [ ] SupermaximalRepeatFinderTests.cs (10 tests)
- [ ] InvertedRepeatFinderTests.cs (10 tests)
- [ ] MicrosatelliteFinderTests.cs (10 tests)

---

### 2.4 Unique Substring Analysis

#### Tasks
- [ ] `ShortestUniqueSubstring` - for primers
- [ ] `MinimalUniqueSubstrings` - all minimal unique
- [ ] `UniqueKmerFinder` - unique k-mers

#### Files
```
Analysis/Unique/
├── ShortestUniqueSubstringFinder.cs
├── MinimalUniqueSubstringsFinder.cs
├── UniqueKmerFinder.cs
└── UniquenessResult.cs
```

#### Tests (20 tests)
- [ ] ShortestUniqueSubstringTests.cs (10 tests)
- [ ] UniqueKmerFinderTests.cs (10 tests)

---

## Phase 3: RNA Analysis (Week 5-6)

### 3.1 RNA Secondary Structure

#### Tasks
- [ ] `StemLoopFinder` - stem-loop structure search
- [ ] `HairpinFinder` - hairpin search
- [ ] `PseudoknotDetector` - pseudoknot detection
- [ ] `FreeEnergyCalculator` - ΔG calculation

#### Files
```
Analysis/Structure/
├── Rna/
│   ├── StemLoopFinder.cs
│   ├── HairpinFinder.cs
│   ├── PseudoknotDetector.cs
│   ├── FreeEnergyCalculator.cs
│   ├── SecondaryStructure.cs
│   └── BasePair.cs
```

#### Tests (25 tests)
- [ ] StemLoopFinderTests.cs (10 tests)
- [ ] HairpinFinderTests.cs (10 tests)
- [ ] FreeEnergyCalculatorTests.cs (5 tests)

---

### 3.2 Splicing Analysis

#### Tasks
- [ ] `SpliceSiteFinder` - donor/acceptor sites
- [ ] `IntronExonPredictor` - intron/exon prediction
- [ ] `AlternativeSplicingDetector` - alternative splicing

#### Files
```
Analysis/Splicing/
├── SpliceSiteFinder.cs
├── SpliceSite.cs
├── IntronExonPredictor.cs
└── AlternativeSplicingDetector.cs
```

#### Tests (20 tests)
- [ ] SpliceSiteFinderTests.cs (10 tests)
- [ ] IntronExonPredictorTests.cs (10 tests)

---

### 3.3 miRNA Analysis

#### Tasks
- [ ] `MiRnaTargetPredictor` - miRNA target prediction
- [ ] `SeedMatcher` - seed region match search
- [ ] `MiRnaScorer` - interaction scoring

#### Files
```
Analysis/MiRna/
├── MiRnaTargetPredictor.cs
├── SeedMatcher.cs
├── MiRnaScorer.cs
├── TargetSite.cs
└── MiRnaDatabase.cs
```

#### Tests (15 tests)
- [ ] MiRnaTargetPredictorTests.cs (10 tests)
- [ ] SeedMatcherTests.cs (5 tests)

---

## Phase 4: Protein Analysis (Week 7-8)

### 4.1 Protein Motifs

#### Tasks
- [ ] `ProteinMotifFinder` - PROSITE-like patterns
- [ ] `PrositeParser` - PROSITE format parser
- [ ] `RegexMotifMatcher` - regex for proteins

#### Files
```
Analysis/Protein/
├── ProteinMotifFinder.cs
├── PrositeParser.cs
├── PrositePattern.cs
├── RegexMotifMatcher.cs
└── ProteinMotifResult.cs
```

#### Tests (20 tests)
- [ ] ProteinMotifFinderTests.cs (10 tests)
- [ ] PrositeParserTests.cs (10 tests)

---

### 4.2 Protein Structure Prediction

#### Tasks
- [ ] `SignalPeptidePredictor` - signal peptides
- [ ] `TransmembranePredictor` - TM domains
- [ ] `CoiledCoilPredictor` - coiled-coil structures
- [ ] `DisorderPredictor` - disordered regions

#### Files
```
Analysis/Protein/
├── SignalPeptidePredictor.cs
├── TransmembranePredictor.cs
├── CoiledCoilPredictor.cs
├── DisorderPredictor.cs
└── ProteinRegion.cs
```

#### Tests (25 tests)
- [ ] SignalPeptidePredictorTests.cs (7 tests)
- [ ] TransmembranePredictorTests.cs (8 tests)
- [ ] CoiledCoilPredictorTests.cs (5 tests)
- [ ] DisorderPredictorTests.cs (5 tests)

---

### 4.3 Protein Properties

#### Tasks
- [ ] `ProteinPropertiesCalculator`:
  - Molecular weight
  - Isoelectric point (pI)
  - GRAVY (hydropathicity)
  - Instability index
  - Aliphatic index
- [ ] `AminoAcidComposition` - amino acid composition

#### Files
```
Analysis/Protein/
├── ProteinPropertiesCalculator.cs
├── AminoAcidComposition.cs
├── AminoAcidProperties.cs
└── ProteinProperties.cs
```

#### Tests (15 tests)
- [ ] ProteinPropertiesCalculatorTests.cs (10 tests)
- [ ] AminoAcidCompositionTests.cs (5 tests)

---

## Phase 5: Sequence Alignment (Week 9-10)

### 5.1 Pairwise Alignment

#### Tasks
- [ ] `GlobalAligner` - Needleman-Wunsch
- [ ] `LocalAligner` - Smith-Waterman
- [ ] `SemiGlobalAligner` - overlap alignment
- [ ] `ScoringMatrix` - BLOSUM, PAM matrices

#### Files
```
Alignment/
├── Pairwise/
│   ├── GlobalAligner.cs
│   ├── LocalAligner.cs
│   ├── SemiGlobalAligner.cs
│   ├── AlignmentResult.cs
│   └── AlignedSequence.cs
├── Scoring/
│   ├── IScoringMatrix.cs
│   ├── NucleotideScoringMatrix.cs
│   ├── BlosumMatrix.cs
│   ├── PamMatrix.cs
│   └── ScoringMatrixLoader.cs
```

#### Tests (30 tests)
- [ ] GlobalAlignerTests.cs (10 tests)
- [ ] LocalAlignerTests.cs (10 tests)
- [ ] ScoringMatrixTests.cs (10 tests)

---

### 5.2 Multiple Sequence Alignment

#### Tasks
- [ ] `ProgressiveAligner` - progressive MSA
- [ ] `ConsensusBuilder` - consensus building
- [ ] `AlignmentScorer` - MSA quality scoring

#### Files
```
Alignment/
├── Multiple/
│   ├── ProgressiveAligner.cs
│   ├── GuideTree.cs
│   ├── ConsensusBuilder.cs
│   ├── AlignmentScorer.cs
│   └── MultipleAlignment.cs
```

#### Tests (20 tests)
- [ ] ProgressiveAlignerTests.cs (10 tests)
- [ ] ConsensusBuilderTests.cs (10 tests)

---

### 5.3 Suffix Tree-based Alignment

#### Tasks
- [ ] `MummerAligner` - MUM-based alignment
- [ ] `MaximalUniqueMatchFinder` - MUM finder
- [ ] `AnchorChainer` - chaining anchors

#### Files
```
Alignment/
├── SuffixTreeBased/
│   ├── MummerAligner.cs
│   ├── MaximalUniqueMatchFinder.cs
│   ├── AnchorChainer.cs
│   └── MummerResult.cs
```

#### Tests (20 tests)
- [ ] MummerAlignerTests.cs (10 tests)
- [ ] MaximalUniqueMatchFinderTests.cs (10 tests)

---

## Phase 6: Comparative Genomics (Week 11-12)

### 6.1 Genome Comparison

#### Tasks
- [ ] `SnpDetector` - SNP detection
- [ ] `IndelDetector` - insertion/deletion detection
- [ ] `SyntenyFinder` - synteny blocks
- [ ] `GeneDuplicationFinder` - gene duplications

#### Files
```
Comparison/
├── SnpDetector.cs
├── IndelDetector.cs
├── SyntenyFinder.cs
├── GeneDuplicationFinder.cs
├── Snp.cs
├── Indel.cs
└── SyntenyBlock.cs
```

#### Tests (25 tests)
- [ ] SnpDetectorTests.cs (8 tests)
- [ ] IndelDetectorTests.cs (7 tests)
- [ ] SyntenyFinderTests.cs (10 tests)

---

### 6.2 Pan-Genome Analysis

#### Tasks
- [ ] `CoreGenomeFinder` - core genome
- [ ] `PanGenomeBuilder` - pan genome
- [ ] `AccessoryGeneFinder` - accessory genes

#### Files
```
Comparison/
├── PanGenome/
│   ├── CoreGenomeFinder.cs
│   ├── PanGenomeBuilder.cs
│   ├── AccessoryGeneFinder.cs
│   └── PanGenomeResult.cs
```

#### Tests (15 tests)
- [ ] CoreGenomeFinderTests.cs (8 tests)
- [ ] PanGenomeBuilderTests.cs (7 tests)

---

## Phase 7: Statistics (Week 13)

### 7.1 Sequence Statistics

#### Tasks
- [ ] `KmerAnalyzer` - k-mer spectrum
- [ ] `SequenceComplexity` - linguistic complexity
- [ ] `EntropyCalculator` - Shannon entropy
- [ ] `GcSkewCalculator` - GC skew analysis

#### Files
```
Analysis/Statistics/
├── KmerAnalyzer.cs
├── KmerSpectrum.cs
├── SequenceComplexity.cs
├── EntropyCalculator.cs
├── GcSkewCalculator.cs
└── StatisticsResult.cs
```

#### Tests (25 tests)
- [ ] KmerAnalyzerTests.cs (10 tests)
- [ ] SequenceComplexityTests.cs (5 tests)
- [ ] EntropyCalculatorTests.cs (5 tests)
- [ ] GcSkewCalculatorTests.cs (5 tests)

---

### 7.2 Codon Analysis

#### Tasks
- [ ] `CodonUsageAnalyzer` - codon usage bias
- [ ] `CaiCalculator` - Codon Adaptation Index
- [ ] `RscuCalculator` - Relative Synonymous Codon Usage

#### Files
```
Analysis/Statistics/
├── CodonUsageAnalyzer.cs
├── CaiCalculator.cs
├── RscuCalculator.cs
└── CodonUsageTable.cs
```

#### Tests (15 tests)
- [ ] CodonUsageAnalyzerTests.cs (8 tests)
- [ ] CaiCalculatorTests.cs (7 tests)

---

## Phase 8: Practical Applications (Week 14-15)

### 8.1 Primer Design

#### Tasks
- [ ] `PrimerDesigner` - PCR primer design
- [ ] `PrimerValidator` - primer validation
- [ ] `MeltingTemperatureCalculator` - Tm calculation
- [ ] `PrimerDimerChecker` - dimer checking

#### Files
```
Applications/Primers/
├── PrimerDesigner.cs
├── PrimerValidator.cs
├── MeltingTemperatureCalculator.cs
├── PrimerDimerChecker.cs
├── Primer.cs
├── PrimerPair.cs
└── PrimerDesignOptions.cs
```

#### Tests (30 tests)
- [ ] PrimerDesignerTests.cs (10 tests)
- [ ] PrimerValidatorTests.cs (10 tests)
- [ ] MeltingTemperatureCalculatorTests.cs (10 tests)

---

### 8.2 CRISPR Guide RNA

#### Tasks
- [ ] `GuideRnaDesigner` - gRNA design
- [ ] `PamFinder` - PAM sequence search
- [ ] `OffTargetPredictor` - off-target prediction
- [ ] `GuideRnaScorer` - gRNA scoring

#### Files
```
Applications/Crispr/
├── GuideRnaDesigner.cs
├── PamFinder.cs
├── OffTargetPredictor.cs
├── GuideRnaScorer.cs
├── GuideRna.cs
├── PamSequence.cs
└── CrisprSystem.cs
```

#### Tests (25 tests)
- [ ] GuideRnaDesignerTests.cs (10 tests)
- [ ] PamFinderTests.cs (5 tests)
- [ ] OffTargetPredictorTests.cs (10 tests)

---

### 8.3 Restriction Analysis

#### Tasks
- [ ] `RestrictionMapper` - restriction map
- [ ] `RestrictionEnzymeDatabase` - enzyme database
- [ ] `DigestSimulator` - restriction simulation
- [ ] `FragmentAnalyzer` - fragment analysis

#### Files
```
Applications/Restriction/
├── RestrictionMapper.cs
├── RestrictionEnzymeDatabase.cs
├── RestrictionEnzyme.cs
├── DigestSimulator.cs
├── FragmentAnalyzer.cs
└── RestrictionSite.cs
```

#### Tests (20 tests)
- [ ] RestrictionMapperTests.cs (10 tests)
- [ ] DigestSimulatorTests.cs (10 tests)

---

### 8.4 Probe Design

#### Tasks
- [ ] `ProbeDesigner` - hybridization probe design
- [ ] `ProbeSpecificityChecker` - specificity checking
- [ ] `OligoAnalyzer` - oligonucleotide analysis

#### Files
```
Applications/Probes/
├── ProbeDesigner.cs
├── ProbeSpecificityChecker.cs
├── OligoAnalyzer.cs
├── Probe.cs
└── ProbeDesignOptions.cs
```

#### Tests (15 tests)
- [ ] ProbeDesignerTests.cs (10 tests)
- [ ] OligoAnalyzerTests.cs (5 tests)

---

## Phase 9: File Formats (Week 16)

### 9.1 Input Formats

#### Tasks
- [ ] `GenBankParser` - GenBank format
- [ ] `EmblParser` - EMBL format
- [ ] `FastqParser` - FASTQ with quality
- [ ] `GffParser` - GFF/GTF annotations
- [ ] `BedParser` - BED format
- [ ] `VcfParser` - VCF (variants)

#### Files
```
IO/
├── Parsers/
│   ├── GenBankParser.cs
│   ├── EmblParser.cs
│   ├── FastqParser.cs
│   ├── GffParser.cs
│   ├── BedParser.cs
│   └── VcfParser.cs
├── Models/
│   ├── GenBankRecord.cs
│   ├── FastqRecord.cs
│   ├── GffFeature.cs
│   ├── BedRegion.cs
│   └── VcfVariant.cs
```

#### Tests (35 tests)
- [ ] GenBankParserTests.cs (8 tests)
- [ ] FastqParserTests.cs (7 tests)
- [ ] GffParserTests.cs (7 tests)
- [ ] BedParserTests.cs (6 tests)
- [ ] VcfParserTests.cs (7 tests)

---

### 9.2 Output Formats

#### Tasks
- [ ] `FastaWriter` - FASTA output (refactor)
- [ ] `GenBankWriter` - GenBank output
- [ ] `GffWriter` - GFF output
- [ ] `ReportGenerator` - HTML/JSON reports

#### Files
```
IO/
├── Writers/
│   ├── FastaWriter.cs
│   ├── GenBankWriter.cs
│   ├── GffWriter.cs
│   └── ReportGenerator.cs
```

#### Tests (15 tests)
- [ ] WriterTests.cs (15 tests)

---

## Summary

### Overall Statistics

| Phase | Week | Algorithms | Tests |
|------|---------|-----------|-------|
| 1. Core Foundation | 1-2 | 8 | 40 |
| 2. DNA Analysis | 3-4 | 12 | 95 |
| 3. RNA Analysis | 5-6 | 10 | 60 |
| 4. Protein Analysis | 7-8 | 12 | 60 |
| 5. Alignment | 9-10 | 9 | 70 |
| 6. Comparative Genomics | 11-12 | 7 | 40 |
| 7. Statistics | 13 | 7 | 40 |
| 8. Applications | 14-15 | 14 | 90 |
| 9. File Formats | 16 | 10 | 50 |
| **TOTAL** | **16 weeks** | **89 algorithms** | **545 tests** |

### Priorities (MVP)

**High Priority (implement first):**
1. ✅ DnaSequence (done)
2. ✅ GenomicAnalyzer basics (done)
3. ✅ FastaParser (done)
4. RnaSequence, ProteinSequence
5. ApproximateMatcher
6. KmerAnalyzer
7. PrimerDesigner

**Medium Priority:**
8. RNA Secondary Structure
9. Pairwise Alignment
10. Restriction Analysis
11. CRISPR gRNA

**Lower Priority:**
12. Multiple Alignment
13. Pan-Genome
14. All file formats

---

## Commands to Execute

```bash
# Create structure
mkdir -p Seqeron.Genomics/Core/Sequences
mkdir -p Seqeron.Genomics/Core/Translation
mkdir -p Seqeron.Genomics/Analysis/Repeats
mkdir -p Seqeron.Genomics/Analysis/Motifs
# ... etc

# Run tests
dotnet test --filter "FullyQualifiedName~Genomics"

# Run with coverage
dotnet test --collect:"XPlat Code Coverage"
```

---

## Definition of Done (per algorithm)

- [ ] Implementation with XML documentation
- [ ] Unit tests (min. 5 per algorithm)
- [ ] Integration tests with real data
- [ ] Performance benchmarks
- [ ] README documentation
