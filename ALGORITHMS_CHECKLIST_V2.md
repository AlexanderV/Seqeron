# Algorithms Checklist v2.0

**Date:** 2026-01-23
**Version:** 2.1 (100% Coverage)
**Library:** Seqeron.Genomics

---

## Quick Reference

| Metric | Value |
|--------|-------|
| **Total Test Units** | 100 |
| **Completed** | 41 |
| **In Progress** | 0 |
| **Blocked** | 0 |
| **Not Started** | 59 |

---

## Processing Registry

| Status | Test Unit ID | Area | Methods | Evidence | TestSpec | Test File(s) |
|--------|--------------|------|---------|----------|----------|--------------|
| ☑ | SEQ-GC-001 | Composition | 5 | Wikipedia, Biopython | [SEQ-GC-001.md](TestSpecs/SEQ-GC-001.md) | SequenceExtensions_CalculateGcContent_Tests.cs |
| ☑ | SEQ-COMP-001 | Composition | 3 | Wikipedia, Biopython | [SEQ-COMP-001.md](TestSpecs/SEQ-COMP-001.md) | SequenceExtensions_Complement_Tests.cs |
| ☑ | SEQ-REVCOMP-001 | Composition | 4 | Wikipedia, Biopython | [SEQ-REVCOMP-001.md](TestSpecs/SEQ-REVCOMP-001.md) | SequenceExtensions_ReverseComplement_Tests.cs |
| ☑ | SEQ-VALID-001 | Composition | 4 | Wikipedia, IUPAC 1970, Bioinformatics.org | [SEQ-VALID-001.md](TestSpecs/SEQ-VALID-001.md) | SequenceExtensions_SequenceValidation_Tests.cs |
| ☑ | SEQ-COMPLEX-001 | Composition | 4 | Wikipedia, Troyanskaya (2002), Orlov (2004) | [SEQ-COMPLEX-001.md](TestSpecs/SEQ-COMPLEX-001.md) | SequenceComplexityTests.cs |
| ☑ | SEQ-ENTROPY-001 | Composition | 2 | Wikipedia (Entropy, Sequence logo, K-mer), Shannon (1948) | [SEQ-ENTROPY-001.md](TestSpecs/SEQ-ENTROPY-001.md) | SequenceComplexityTests.cs |
| ☑ | SEQ-GCSKEW-001 | Composition | 4 | Wikipedia, Lobry (1996), Grigoriev (1998) | [SEQ-GCSKEW-001.md](TestSpecs/SEQ-GCSKEW-001.md) | GcSkewCalculatorTests.cs |
| ☑ | PAT-EXACT-001 | Matching | 4 | Wikipedia, Gusfield (1997), Rosalind | [PAT-EXACT-001.md](TestSpecs/PAT-EXACT-001.md) | FindAllOccurrencesTests.cs, ContainsTests.cs, CountOccurrencesTests.cs |
| ☑ | PAT-APPROX-001 | Matching | 2 | Wikipedia (Hamming), Rosalind (HAMM), Gusfield (1997), Navarro (2001) | [PAT-APPROX-001.md](TestSpecs/PAT-APPROX-001.md) | ApproximateMatcher_HammingDistance_Tests.cs |
| ☑ | PAT-APPROX-002 | Matching | 2 | Wikipedia (Levenshtein, Edit Distance), Rosetta Code, Navarro (2001) | [PAT-APPROX-002.md](TestSpecs/PAT-APPROX-002.md) | ApproximateMatcher_EditDistance_Tests.cs |
| ☑ | PAT-IUPAC-001 | Matching | 2 | Wikipedia (Nucleic acid notation), IUPAC-IUB 1970, Bioinformatics.org | [PAT-IUPAC-001.md](TestSpecs/PAT-IUPAC-001.md) | IupacMotifMatchingTests.cs |
| ☑ | PAT-PWM-001 | Matching | 2 | Wikipedia (PWM), Kel et al. (2003), Rosalind (CONS), Nishida (2008) | [PAT-PWM-001.md](TestSpecs/PAT-PWM-001.md) | MotifFinder_PWM_Tests.cs |
| ☑ | REP-STR-001 | Repeats | 4 | Wikipedia (Microsatellite, Trinucleotide repeat disorder), Richard et al. (2008) | [REP-STR-001.md](TestSpecs/REP-STR-001.md) | RepeatFinder_Microsatellite_Tests.cs |
| ☑ | REP-TANDEM-001 | Repeats | 2 | Wikipedia (Tandem repeat, Microsatellite), Richard et al. (2008) | [REP-TANDEM-001.md](TestSpecs/REP-TANDEM-001.md) | GenomicAnalyzer_TandemRepeat_Tests.cs, RepeatFinderTests.cs |
| ☑ | REP-INV-001 | Repeats | 1 | Wikipedia (Inverted repeat, Stem-loop, Palindromic sequence), EMBOSS einverted, Pearson (1996), Bissler (1998) | [REP-INV-001.md](TestSpecs/REP-INV-001.md) | RepeatFinder_InvertedRepeat_Tests.cs |
| ☑ | REP-DIRECT-001 | Repeats | 1 | Wikipedia (Direct repeat, Repeated sequence), Ussery (2009), Richard (2021) | [REP-DIRECT-001.md](TestSpecs/REP-DIRECT-001.md) | RepeatFinder_DirectRepeat_Tests.cs |
| ☑ | REP-PALIN-001 | Repeats | 2 | Wikipedia (Palindromic sequence, Restriction enzyme), Rosalind REVP | [REP-PALIN-001.md](TestSpecs/REP-PALIN-001.md) | RepeatFinder_Palindrome_Tests.cs |
| ☑ | CRISPR-PAM-001 | MolTools | 2 | Wikipedia (Protospacer adjacent motif, CRISPR), Jinek et al. (2012), Zetsche et al. (2015) | [CRISPR-PAM-001.md](TestSpecs/CRISPR-PAM-001.md) | CrisprDesigner_PAM_Tests.cs |
| ☑ | CRISPR-GUIDE-001 | MolTools | 2 | Wikipedia (Guide RNA, CRISPR gene editing, PAM), Addgene CRISPR Guide | [CRISPR-GUIDE-001.md](TestSpecs/CRISPR-GUIDE-001.md) | CrisprDesigner_GuideRNA_Tests.cs |
| ☑ | CRISPR-OFF-001 | MolTools | 2 | Wikipedia (Off-target genome editing), Hsu et al. (2013), Fu et al. (2013) | [CRISPR-OFF-001.md](TestSpecs/CRISPR-OFF-001.md) | CrisprDesigner_OffTarget_Tests.cs |
| ☑ | PRIMER-TM-001 | MolTools | 2 | Wikipedia (Nucleic acid thermodynamics), Marmur & Doty (1962), SantaLucia (1998), Owczarzy (2004) | [PRIMER-TM-001.md](TestSpecs/PRIMER-TM-001.md) | PrimerDesigner_MeltingTemperature_Tests.cs |
| ☑ | PRIMER-DESIGN-001 | MolTools | 3 | Wikipedia, Addgene, Primer3 | [PRIMER-DESIGN-001.md](TestSpecs/PRIMER-DESIGN-001.md) | PrimerDesigner_PrimerDesign_Tests.cs |
| ☑ | PRIMER-STRUCT-001 | MolTools | 3 | Wikipedia (Primer, Primer dimer, Stem-loop, Nucleic acid thermodynamics), Primer3 Manual, SantaLucia (1998) | [PRIMER-STRUCT-001.md](TestSpecs/PRIMER-STRUCT-001.md) | PrimerDesigner_PrimerStructure_Tests.cs |
| ☑ | PROBE-DESIGN-001 | MolTools | 3 | Wikipedia (Nucleic acid thermodynamics, Hybridization probe, FISH, DNA microarray), SantaLucia (1998) | [PROBE-DESIGN-001.md](TestSpecs/PROBE-DESIGN-001.md) | ProbeDesigner_ProbeDesign_Tests.cs |
| ☑ | PROBE-VALID-001 | MolTools | 2 | Wikipedia (Hybridization probe, DNA microarray, Off-target genome editing, BLAST), Amann & Ludwig (2000) | [PROBE-VALID-001.md](TestSpecs/PROBE-VALID-001.md) | ProbeDesigner_ProbeValidation_Tests.cs |
| ☑ | RESTR-FIND-001 | MolTools | 2 | Wikipedia (Restriction enzyme, Restriction site, EcoRI), Roberts (1976), REBASE | [RESTR-FIND-001.md](TestSpecs/RESTR-FIND-001.md) | RestrictionAnalyzer_FindSites_Tests.cs |
| ☑ | RESTR-DIGEST-001 | MolTools | 2 | Wikipedia (Restriction digest, Restriction enzyme, Restriction map), Addgene Protocol, Roberts (1976), REBASE | [RESTR-DIGEST-001.md](TestSpecs/RESTR-DIGEST-001.md) | RestrictionAnalyzer_Digest_Tests.cs |
| ☑ | ANNOT-ORF-001 | Annotation | 3 | Wikipedia (ORF), Rosalind ORF, NCBI ORF Finder, Deonier (2005), Claverie (1997) | [ANNOT-ORF-001.md](TestSpecs/ANNOT-ORF-001.md) | GenomeAnnotator_ORF_Tests.cs |
| ☑ | ANNOT-GENE-001 | Annotation | 2 | Wikipedia (Gene prediction, Shine-Dalgarno, RBS), Shine & Dalgarno (1975), Chen (1994), Laursen (2005) | [ANNOT-GENE-001.md](TestSpecs/ANNOT-GENE-001.md) | GenomeAnnotator_Gene_Tests.cs |
| ☑ | ANNOT-PROM-001 | Annotation | 1 | Wikipedia (Promoter, Pribnow box), Pribnow (1975), Harley & Reynolds (1987) | [ANNOT-PROM-001.md](TestSpecs/ANNOT-PROM-001.md) | GenomeAnnotator_PromoterMotif_Tests.cs |
| ☑ | ANNOT-GFF-001 | Annotation | 2 | Sequence Ontology GFF3 Spec v1.26, Wikipedia (GFF), RFC 3986 | [ANNOT-GFF-001.md](TestSpecs/ANNOT-GFF-001.md) | GenomeAnnotator_GFF3_Tests.cs |
| ☑ | KMER-COUNT-001 | K-mer | 3 | Wikipedia (K-mer), Rosalind (KMER, BA1E) | [KMER-COUNT-001.md](TestSpecs/KMER-COUNT-001.md) | KmerAnalyzer_CountKmers_Tests.cs |
| ☑ | KMER-FREQ-001 | K-mer | 3 | Wikipedia (K-mer, Entropy), Shannon (1948), Rosalind KMER | [KMER-FREQ-001.md](TestSpecs/KMER-FREQ-001.md) | KmerAnalyzer_Frequency_Tests.cs |
| ☑ | KMER-FIND-001 | K-mer | 3 | Wikipedia (K-mer), Rosalind BA1B (frequent words), Rosalind BA1E (clump finding) | [KMER-FIND-001.md](TestSpecs/KMER-FIND-001.md) | KmerAnalyzer_Find_Tests.cs |
| ☑ | ALIGN-GLOBAL-001 | Alignment | 1 | Wikipedia (Needleman–Wunsch, Sequence alignment) | [ALIGN-GLOBAL-001.md](TestSpecs/ALIGN-GLOBAL-001.md) | SequenceAligner_GlobalAlign_Tests.cs, PerformanceExtensionsTests.cs |
| ☑ | ALIGN-LOCAL-001 | Alignment | 1 | Wikipedia (Smith–Waterman, Sequence alignment) | [ALIGN-LOCAL-001.md](TestSpecs/ALIGN-LOCAL-001.md) | SequenceAligner_LocalAlign_Tests.cs |
| ☑ | ALIGN-SEMI-001 | Alignment | 1 | Wikipedia (Sequence alignment, Needleman–Wunsch, Smith–Waterman) | [ALIGN-SEMI-001.md](TestSpecs/ALIGN-SEMI-001.md) | SequenceAligner_SemiGlobalAlign_Tests.cs |
| ☑ | ALIGN-MULTI-001 | Alignment | 1 | Wikipedia (Multiple sequence alignment, Clustal) | [ALIGN-MULTI-001.md](TestSpecs/ALIGN-MULTI-001.md) | SequenceAligner_MultipleAlign_Tests.cs |
| ☑ | PHYLO-DIST-001 | Phylogenetic | 2 | Wikipedia (Models of DNA evolution, Distance matrices in phylogeny, Jukes-Cantor, Kimura 2-parameter) | [PHYLO-DIST-001.md](TestSpecs/PHYLO-DIST-001.md) | PhylogeneticAnalyzer_DistanceMatrix_Tests.cs |
| ☑ | PHYLO-TREE-001 | Phylogenetic | 1 | Wikipedia (UPGMA, Neighbor joining, Phylogenetic tree), Saitou & Nei (1987), Sokal & Michener (1958) | [PHYLO-TREE-001.md](TestSpecs/PHYLO-TREE-001.md) | PhylogeneticAnalyzer_TreeConstruction_Tests.cs |
| ☑ | PHYLO-NEWICK-001 | Phylogenetic | 2 | Wikipedia (Newick format), PHYLIP (Felsenstein), Olsen (1990) | [PHYLO-NEWICK-001.md](TestSpecs/PHYLO-NEWICK-001.md) | PhylogeneticAnalyzer_NewickIO_Tests.cs |
| ☑ | PHYLO-COMP-001 | Phylogenetic | 3 | Wikipedia (Robinson-Foulds, MRCA, Phylogenetic tree), Robinson & Foulds (1981) | [PHYLO-COMP-001.md](TestSpecs/PHYLO-COMP-001.md) | PhylogeneticAnalyzer_TreeComparison_Tests.cs |
| ☑ | POP-FREQ-001 | PopGen | 3 | Wikipedia (Allele frequency, Minor allele frequency, Genotype frequency), Gillespie (2004) | [POP-FREQ-001.md](TestSpecs/POP-FREQ-001.md) | PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs |
| ☑ | POP-DIV-001 | PopGen | 4 | Wikipedia (Nucleotide diversity, Watterson estimator, Tajima's D), Nei & Li (1979), Watterson (1975), Tajima (1989) | [POP-DIV-001.md](TestSpecs/POP-DIV-001.md) | PopulationGeneticsAnalyzer_Diversity_Tests.cs |
| ☑ | POP-HW-001 | PopGen | 1 | Wikipedia (Hardy-Weinberg principle, Chi-squared test), Hardy (1908), Weinberg (1908), Emigh (1980) | [POP-HW-001.md](TestSpecs/POP-HW-001.md) | PopulationGeneticsAnalyzer_HardyWeinberg_Tests.cs |
| ☑ | POP-FST-001 | PopGen | 2 | Wikipedia (Fixation index, F-statistics), Wright (1965), Weir & Cockerham (1984) | [POP-FST-001.md](TestSpecs/POP-FST-001.md) | PopulationGeneticsAnalyzer_FStatistics_Tests.cs |
| ☑ | POP-LD-001 | PopGen | 2 | Wikipedia (Linkage disequilibrium, Haplotype block), Lewontin (1964), Hill & Robertson (1968), Gabriel (2002) | [POP-LD-001.md](TestSpecs/POP-LD-001.md) | PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs |
| ☑ | CHROM-TELO-001 | Chromosome | 2 | Wikipedia (Telomere), Meyne (1989), Cawthon (2002) | [CHROM-TELO-001.md](TestSpecs/CHROM-TELO-001.md) | ChromosomeAnalyzer_Telomere_Tests.cs |
| ☑ | CHROM-CENT-001 | Chromosome | 1 | Wikipedia (Centromere, Karyotype, Chromosome), Levan (1964) | [CHROM-CENT-001.md](TestSpecs/CHROM-CENT-001.md) | ChromosomeAnalyzer_Centromere_Tests.cs |
| ☑ | CHROM-KARYO-001 | Chromosome | 2 | Wikipedia (Karyotype, Ploidy, Aneuploidy) | [CHROM-KARYO-001.md](TestSpecs/CHROM-KARYO-001.md) | ChromosomeAnalyzer_Karyotype_Tests.cs |
| ☐ | CHROM-ANEU-001 | Chromosome | 2 | - | - | - |
| ☐ | CHROM-SYNT-001 | Chromosome | 2 | - | - | - |
| ☐ | META-CLASS-001 | Metagenomics | 2 | - | - | - |
| ☐ | META-PROF-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-ALPHA-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-BETA-001 | Metagenomics | 1 | - | - | - |
| ☐ | META-BIN-001 | Metagenomics | 1 | - | - | - |
| ☐ | CODON-OPT-001 | Codon | 1 | - | - | - |
| ☐ | CODON-CAI-001 | Codon | 1 | - | - | - |
| ☐ | CODON-RARE-001 | Codon | 1 | - | - | - |
| ☐ | CODON-USAGE-001 | Codon | 2 | - | - | - |
| ☐ | TRANS-CODON-001 | Translation | 3 | - | - | - |
| ☐ | TRANS-PROT-001 | Translation | 1 | - | - | - |
| ☐ | PARSE-FASTA-001 | FileIO | 4 | - | - | - |
| ☐ | PARSE-FASTQ-001 | FileIO | 4 | - | - | - |
| ☐ | PARSE-BED-001 | FileIO | 6 | - | - | - |
| ☐ | PARSE-VCF-001 | FileIO | 4 | - | - | - |
| ☐ | PARSE-GFF-001 | FileIO | 3 | - | - | - |
| ☐ | PARSE-GENBANK-001 | FileIO | 3 | - | - | - |
| ☐ | PARSE-EMBL-001 | FileIO | 2 | - | - | - |
| ☐ | RNA-STRUCT-001 | RnaStructure | 4 | - | - | - |
| ☐ | RNA-STEMLOOP-001 | RnaStructure | 3 | - | - | - |
| ☐ | RNA-ENERGY-001 | RnaStructure | 2 | - | - | - |
| ☐ | MIRNA-SEED-001 | MiRNA | 3 | - | - | - |
| ☐ | MIRNA-TARGET-001 | MiRNA | 2 | - | - | - |
| ☐ | MIRNA-PRECURSOR-001 | MiRNA | 2 | - | - | - |
| ☐ | SPLICE-DONOR-001 | Splicing | 2 | - | - | - |
| ☐ | SPLICE-ACCEPTOR-001 | Splicing | 2 | - | - | - |
| ☐ | SPLICE-PREDICT-001 | Splicing | 3 | - | - | - |
| ☐ | DISORDER-PRED-001 | ProteinPred | 3 | - | - | - |
| ☐ | DISORDER-REGION-001 | ProteinPred | 2 | - | - | - |
| ☐ | PROTMOTIF-FIND-001 | ProteinMotif | 3 | - | - | - |
| ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | 2 | - | - | - |
| ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | 2 | - | - | - |
| ☐ | EPIGEN-CPG-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-METHYL-001 | Epigenetics | 3 | - | - | - |
| ☐ | EPIGEN-DMR-001 | Epigenetics | 2 | - | - | - |
| ☐ | VARIANT-CALL-001 | Variants | 3 | - | - | - |
| ☐ | VARIANT-SNP-001 | Variants | 2 | - | - | - |
| ☐ | VARIANT-INDEL-001 | Variants | 2 | - | - | - |
| ☐ | VARIANT-ANNOT-001 | Variants | 2 | - | - | - |
| ☐ | SV-DETECT-001 | StructuralVar | 3 | - | - | - |
| ☐ | SV-BREAKPOINT-001 | StructuralVar | 2 | - | - | - |
| ☐ | SV-CNV-001 | StructuralVar | 2 | - | - | - |
| ☐ | ASSEMBLY-OLC-001 | Assembly | 2 | - | - | - |
| ☐ | ASSEMBLY-DBG-001 | Assembly | 2 | - | - | - |
| ☐ | ASSEMBLY-STATS-001 | Assembly | 4 | - | - | - |
| ☐ | TRANS-EXPR-001 | Transcriptome | 3 | - | - | - |
| ☐ | TRANS-DIFF-001 | Transcriptome | 2 | - | - | - |
| ☐ | TRANS-SPLICE-001 | Transcriptome | 2 | - | - | - |
| ☐ | COMPGEN-SYNTENY-001 | Comparative | 2 | - | - | - |
| ☐ | COMPGEN-ORTHO-001 | Comparative | 2 | - | - | - |
| ☐ | COMPGEN-REARR-001 | Comparative | 2 | - | - | - |
| ☐ | PANGEN-CORE-001 | PanGenome | 2 | - | - | - |
| ☐ | PANGEN-CLUSTER-001 | PanGenome | 2 | - | - | - |
| ☐ | QUALITY-PHRED-001 | Quality | 3 | - | - | - |
| ☐ | QUALITY-STATS-001 | Quality | 2 | - | - | - |

**Statuses:** ☐ Not Started | ⏳ In Progress | ☑ Complete | ⛔ Blocked

---

## Definition of Done (DoD)

### Required Criteria

| # | Criterion | Artifact |
|---|-----------|----------|
| 1 | TestSpec created | `TestSpecs/{TestUnitID}.md` |
| 2 | Tests written | `*.Tests/{Class}_{Method}_Tests.cs` |
| 3 | Branch coverage ≥ 80% | Coverage report |
| 4 | Edge cases covered | null, empty, boundary, error |
| 5 | Tests pass | CI green |
| 6 | Evidence documented | PR/commit link in Registry |

### Test Quality Criteria

- [ ] Tests are independent (order-independent)
- [ ] Tests are deterministic (no random without seed)
- [ ] Naming: `Method_Scenario_ExpectedResult`
- [ ] Structure: Arrange-Act-Assert
- [ ] One assert per logical check

### For O(n²) and Higher Algorithms

- [ ] Property-based test for invariant
- [ ] Performance baseline recorded

---

## Test Units by Area

### 1. Sequence Composition (4 units)

#### SEQ-GC-001: GC Content Calculation

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ result ≤ 100 (percentage) or 0 ≤ result ≤ 1 (fraction) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateGcContent(ReadOnlySpan<char>)` | SequenceExtensions | Canonical |
| `CalculateGcFraction(ReadOnlySpan<char>)` | SequenceExtensions | Variant (0-1) |
| `CalculateGcContentFast(string)` | SequenceExtensions | Delegate |
| `CalculateGcFractionFast(string)` | SequenceExtensions | Delegate |
| `GcContent` (property) | DnaSequence | Delegate |

**Edge Cases:**
- [ ] Empty sequence → 0
- [ ] All G/C → 100% / 1.0
- [ ] All A/T → 0% / 0.0
- [ ] Mixed case input
- [ ] Non-ACGT characters (N, etc.)

---

#### SEQ-COMP-001: DNA Complement

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.GetComplementBase(char)` |
| **Complexity** | O(n) for sequence |
| **Invariant** | Complement(Complement(x)) = x |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetComplementBase(char)` | SequenceExtensions | Canonical |
| `TryGetComplement(ReadOnlySpan, Span)` | SequenceExtensions | Span API |
| `Complement()` | DnaSequence | Instance |

**Edge Cases:**
- [ ] A ↔ T, G ↔ C (both directions)
- [ ] Case insensitivity (a → T)
- [ ] RNA support (U → A)
- [ ] Unknown base → unchanged
- [ ] Destination too small → false

---

#### SEQ-REVCOMP-001: Reverse Complement

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.TryGetReverseComplement(ReadOnlySpan, Span)` |
| **Complexity** | O(n) |
| **Invariant** | RevComp(RevComp(x)) = x |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `TryGetReverseComplement(ReadOnlySpan, Span)` | SequenceExtensions | Canonical |
| `ReverseComplement()` | DnaSequence | Instance |
| `GetReverseComplementString(string)` | DnaSequence | Static helper |
| `TryWriteReverseComplement(Span)` | DnaSequence | Span API |

**Edge Cases:**
- [ ] Empty sequence
- [ ] Single nucleotide
- [ ] Palindrome (self-complementary)
- [ ] Destination too small → false

---

#### SEQ-VALID-001: Sequence Validation

| Field | Value |
|-------|-------|
| **Canonical** | `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `IsValidDna(ReadOnlySpan<char>)` | SequenceExtensions | Canonical DNA |
| `IsValidRna(ReadOnlySpan<char>)` | SequenceExtensions | Canonical RNA |
| `TryCreate(string, out DnaSequence)` | DnaSequence | Factory |
| `IsValid` (property) | DnaSequence | Instance |

**Edge Cases:**
- [ ] Empty → true (or false? define!)
- [ ] All valid bases
- [ ] Single invalid character
- [ ] Lowercase valid
- [ ] Whitespace handling

---

### 2. Pattern Matching (5 units)

#### PAT-EXACT-001: Exact Pattern Search

| Field | Value |
|------|----------|
| **Canonical** | `SuffixTree.FindAllOccurrences(string)` |
| **Complexity** | O(m + k) where m=pattern length, k=occurrences |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindAllOccurrences(string)` | SuffixTree/DnaSequence | Canonical |
| `Contains(string)` | SuffixTree/DnaSequence | Existence check |
| `CountOccurrences(string)` | SuffixTree/DnaSequence | Count only |
| `FindMotif(DnaSequence, string)` | GenomicAnalyzer | Wrapper |
| `FindExactMotif(DnaSequence, string)` | MotifFinder | Wrapper |

**Edge Cases:**
- [ ] Pattern not found → empty
- [ ] Pattern = entire sequence
- [ ] Overlapping occurrences
- [ ] Empty pattern
- [ ] Pattern longer than sequence

---

#### PAT-APPROX-001: Approximate Matching (Hamming)

| Field | Value |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindWithMismatches(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | HammingDistance requires equal length |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindWithMismatches(seq, pattern, max)` | ApproximateMatcher | Canonical |
| `HammingDistance(s1, s2)` | ApproximateMatcher | Distance |
| `HammingDistance(span1, span2)` | SequenceExtensions | Span API |

**Edge Cases:**
- [ ] maxMismatches = 0 → exact match
- [ ] maxMismatches ≥ pattern length → all positions
- [ ] Unequal lengths for HammingDistance → exception

---

#### PAT-APPROX-002: Approximate Matching (Edit Distance)

| Field | Value |
|------|----------|
| **Canonical** | `ApproximateMatcher.FindWithEdits(...)` |
| **Complexity** | O(n × m²) |
| **Invariant** | EditDistance(a,b) = EditDistance(b,a) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindWithEdits(seq, pattern, maxEdits)` | ApproximateMatcher | Canonical |
| `EditDistance(s1, s2)` | ApproximateMatcher | Distance |

**Edge Cases:**
- [ ] Identical strings → 0
- [ ] One empty string → length of other
- [ ] Single character difference
- [ ] Insertion vs deletion vs substitution

---

#### PAT-IUPAC-001: IUPAC Degenerate Matching

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.FindDegenerateMotif(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDegenerateMotif(seq, motif)` | MotifFinder | Canonical |
| `MatchesIupac(nucleotide, iupacCode)` | IupacHelper | Helper |

**Edge Cases:**
- [ ] All IUPAC codes: R, Y, S, W, K, M, B, D, H, V, N
- [ ] Mixed standard + IUPAC
- [ ] N matches everything

---

#### PAT-PWM-001: Position Weight Matrix

| Field | Value |
|------|----------|
| **Canonical** | `MotifFinder.ScanWithPwm(...)` |
| **Complexity** | O(n × m) where m=PWM width |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CreatePwm(alignedSequences)` | MotifFinder | Construction |
| `ScanWithPwm(seq, pwm, threshold)` | MotifFinder | Scanning |

**Edge Cases:**
- [ ] Empty alignment
- [ ] Single sequence alignment
- [ ] Threshold at boundary
- [ ] All positions below threshold

---

### 3. Repeat Analysis (5 units)

#### REP-STR-001: Microsatellite Detection (STR)

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindMicrosatellites(...)` |
| **Complexity** | O(n × U × R) where U=maxUnitLength, R=maxRepeats |
| **Invariant** | Result positions are non-overlapping (or document overlap policy) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMicrosatellites(DnaSequence, ...)` | RepeatFinder | Canonical |
| `FindMicrosatellites(string, ...)` | RepeatFinder | Overload |
| `FindMicrosatellites(..., CancellationToken)` | RepeatFinder | Cancellable |
| `FindMicrosatellites(..., IProgress)` | RepeatFinder | With progress |

**Edge Cases:**
- [ ] No repeats found
- [ ] Entire sequence is one repeat
- [ ] minRepeats = 2 (minimum)
- [ ] Unit length 1-6 (mono to hexa)
- [ ] Cancellation mid-operation

---

#### REP-TANDEM-001: Tandem Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomicAnalyzer.FindTandemRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTandemRepeats(seq, minUnit, maxUnit, minReps)` | GenomicAnalyzer | Canonical |
| `GetTandemRepeatSummary(seq, minRepeats)` | RepeatFinder | Summary |

---

#### REP-INV-001: Inverted Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindInvertedRepeats(...)` |
| **Complexity** | O(n² × L) where L=maxLoopLength |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindInvertedRepeats(seq, minArm, maxLoop)` | RepeatFinder | Canonical |

**Edge Cases:**
- [ ] Perfect palindrome (loop = 0)
- [ ] Maximum loop length
- [ ] Arm length at boundaries

---

#### REP-DIRECT-001: Direct Repeat Detection

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindDirectRepeats(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDirectRepeats(seq, minLen, maxLen, minSpacing)` | RepeatFinder | Canonical |

---

#### REP-PALIN-001: Palindrome Detection

| Field | Value |
|------|----------|
| **Canonical** | `RepeatFinder.FindPalindromes(...)` |
| **Complexity** | O(n²) |
| **Invariant** | Palindrome = reverse complement equals self |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPalindromes(seq, minLen, maxLen)` | RepeatFinder | Canonical |
| `FindPalindromes(seq, minLen, maxLen)` | GenomicAnalyzer | Alternate |

---

### 4. Molecular Tools (8 units)

#### CRISPR-PAM-001: PAM Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `CrisprDesigner.FindPamSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPamSites(seq, systemType)` | CrisprDesigner | Canonical |
| `GetSystem(systemType)` | CrisprDesigner | System info |

**CRISPR Systems:**
- [ ] SpCas9 (NGG)
- [ ] SpCas9-NAG
- [ ] SaCas9 (NNGRRT)
- [ ] Cas12a (TTTV)
- [ ] AsCas12a, LbCas12a
- [ ] CasX

---

#### CRISPR-GUIDE-001: Guide RNA Design

| Field | Value |
|------|----------|
| **Canonical** | `CrisprDesigner.DesignGuideRnas(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DesignGuideRnas(seq, start, end, type)` | CrisprDesigner | Canonical |
| `EvaluateGuideRna(guide, type, params)` | CrisprDesigner | Scoring |

---

#### CRISPR-OFF-001: Off-Target Analysis

| Field | Value |
|------|----------|
| **Canonical** | `CrisprDesigner.FindOffTargets(...)` |
| **Complexity** | O(n × m) with maxMismatches; may be higher |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOffTargets(guide, genome, maxMismatches)` | CrisprDesigner | Canonical |
| `CalculateSpecificityScore(guide, genome, type)` | CrisprDesigner | Score |

---

#### PRIMER-TM-001: Melting Temperature

| Field | Value |
|------|----------|
| **Canonical** | `PrimerDesigner.CalculateMeltingTemperature(...)` |
| **Complexity** | O(n) |
| **Formula** | Wallace rule (<14bp), Marmur-Doty (≥14bp) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateMeltingTemperature(primer)` | PrimerDesigner | Canonical |
| `CalculateMeltingTemperatureWithSalt(primer, Na)` | PrimerDesigner | Salt corrected |

**Constants (ThermoConstants):**
- `WallaceMaxLength` = 14
- `CalculateWallaceTm(at, gc)` = 2×AT + 4×GC
- `CalculateMarmurDotyTm(gc, len)` = 64.9 + 41×(GC-16.4)/len
- `CalculateSaltCorrection(Na)` = 16.6 × log10(Na/1000)

**Edge Cases:**
- [ ] Empty primer → 0
- [ ] Short primer (<14) uses Wallace
- [ ] Long primer (≥14) uses Marmur-Doty
- [ ] Salt concentration variations

---

#### PRIMER-DESIGN-001: Primer Pair Design

| Field | Value |
|------|----------|
| **Canonical** | `PrimerDesigner.DesignPrimers(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DesignPrimers(template, start, end, params)` | PrimerDesigner | Canonical |
| `EvaluatePrimer(seq, pos, isForward, params)` | PrimerDesigner | Single primer |
| `GeneratePrimerCandidates(template, region)` | PrimerDesigner | All candidates |

**Parameters (PrimerParameters):**
- MinLength=18, MaxLength=25, OptimalLength=20
- MinGcContent=40, MaxGcContent=60
- MinTm=55, MaxTm=65, OptimalTm=60
- MaxHomopolymer=4, MaxDinucleotideRepeats=4

---

#### PRIMER-STRUCT-001: Primer Structure Analysis

| Field | Value |
|------|----------|
| **Canonical** | `PrimerDesigner.HasHairpinPotential(...)` |
| **Complexity** | O(m²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `HasHairpinPotential(seq, minStemLength)` | PrimerDesigner | Hairpin |
| `HasPrimerDimer(primer1, primer2, minComp)` | PrimerDesigner | Dimer |
| `Calculate3PrimeStability(seq)` | PrimerDesigner | ΔG calculation |
| `FindLongestHomopolymer(seq)` | PrimerDesigner | Structure |
| `FindLongestDinucleotideRepeat(seq)` | PrimerDesigner | Structure |

---

#### RESTR-FIND-001: Restriction Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.FindSites(...)` |
| **Complexity** | O(n × k) where k=enzymes |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSites(seq, enzymeNames)` | RestrictionAnalyzer | Canonical |
| `FindAllSites(seq)` | RestrictionAnalyzer | All 40+ enzymes |
| `GetEnzyme(name)` | RestrictionAnalyzer | Lookup |

**Enzyme Database:** 40+ enzymes (EcoRI, BamHI, HindIII, NotI, etc.)

---

#### RESTR-DIGEST-001: Digest Simulation

| Field | Value |
|------|----------|
| **Canonical** | `RestrictionAnalyzer.Digest(...)` |
| **Complexity** | O(n + k log k) where k=cut sites |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Digest(seq, enzymeNames)` | RestrictionAnalyzer | Canonical |
| `GetDigestSummary(seq, enzymeNames)` | RestrictionAnalyzer | Summary |
| `CreateMap(seq, enzymeNames)` | RestrictionAnalyzer | Full map |

---

### 5. Genome Annotation (4 units)

#### ANNOT-ORF-001: ORF Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindOrfs(...)` |
| **Complexity** | O(n) |
| **Invariant** | ORF starts with start codon, ends with stop codon |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOrfs(dna, minLen, bothStrands, requireStart)` | GenomeAnnotator | Canonical |
| `FindLongestOrfsPerFrame(dna, bothStrands)` | GenomeAnnotator | Per-frame |
| `FindOpenReadingFrames(seq, minLen)` | GenomicAnalyzer | Alternate |

**Edge Cases:**
- [ ] No ORF found
- [ ] ORF extends to sequence end (no stop)
- [ ] All 6 reading frames
- [ ] Overlapping ORFs

---

#### ANNOT-GENE-001: Gene Prediction

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.PredictGenes(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictGenes(dna, minOrfLen, prefix)` | GenomeAnnotator | Canonical |
| `FindRibosomeBindingSites(dna, window)` | GenomeAnnotator | RBS/SD |

---

#### ANNOT-PROM-001: Promoter Detection

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.FindPromoterMotifs(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPromoterMotifs(dna)` | GenomeAnnotator | Canonical |

**Motifs:** -35 box (TTGACA), -10 box (TATAAT)

---

#### ANNOT-GFF-001: GFF3 I/O

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAnnotator.ParseGff3(...)` / `ToGff3(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ParseGff3(lines)` | GenomeAnnotator | Parse |
| `ToGff3(annotations, seqId)` | GenomeAnnotator | Export |

---

### 6. K-mer Analysis (3 units)

#### KMER-COUNT-001: K-mer Counting

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.CountKmers(...)` |
| **Complexity** | O(n) |
| **Invariant** | Sum of counts = n - k + 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CountKmers(seq, k)` | KmerAnalyzer | Canonical |
| `CountKmersSpan(seq, k)` | SequenceExtensions | Span API |
| `CountKmersBothStrands(dnaSeq, k)` | KmerAnalyzer | Both strands |

---

#### KMER-FREQ-001: K-mer Frequency Analysis

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.GetKmerSpectrum(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetKmerSpectrum(seq, k)` | KmerAnalyzer | Spectrum |
| `GetKmerFrequencies(seq, k)` | KmerAnalyzer | Normalized |
| `CalculateKmerEntropy(seq, k)` | KmerAnalyzer | Entropy |

---

#### KMER-FIND-001: K-mer Search

| Field | Value |
|------|----------|
| **Canonical** | `KmerAnalyzer.FindMostFrequentKmers(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMostFrequentKmers(seq, k)` | KmerAnalyzer | Most frequent |
| `FindUniqueKmers(seq, k)` | KmerAnalyzer | Count = 1 |
| `FindClumps(seq, k, window, minOcc)` | KmerAnalyzer | Clumps |

---

### 7. Alignment (4 units)

#### ALIGN-GLOBAL-001: Global Alignment (Needleman-Wunsch)

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.GlobalAlign(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | Optimal global alignment score |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GlobalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-LOCAL-001: Local Alignment (Smith-Waterman)

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.LocalAlign(...)` |
| **Complexity** | O(n × m) |
| **Invariant** | Score ≥ 0, finds best local match |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `LocalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-SEMI-001: Semi-Global Alignment

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.SemiGlobalAlign(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `SemiGlobalAlign(seq1, seq2, scoring)` | SequenceAligner | Canonical |

---

#### ALIGN-MULTI-001: Multiple Sequence Alignment

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAligner.MultipleAlign(...)` |
| **Complexity** | O(n² × m) progressive alignment |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `MultipleAlign(sequences)` | SequenceAligner | Canonical |

---

### 8. Phylogenetics (4 units)

#### PHYLO-DIST-001: Distance Matrix

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.CalculateDistanceMatrix(...)` |
| **Complexity** | O(n² × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateDistanceMatrix(seqs, method)` | PhylogeneticAnalyzer | Canonical |
| `CalculatePairwiseDistance(s1, s2, method)` | PhylogeneticAnalyzer | Single pair |

**Distance Methods:** p-distance, Jukes-Cantor, Kimura 2-parameter, Hamming

---

#### PHYLO-TREE-001: Tree Construction

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.BuildTree(...)` |
| **Complexity** | O(n³) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `BuildTree(seqs, distMethod, treeMethod)` | PhylogeneticAnalyzer | Canonical |

**Tree Methods:** UPGMA, Neighbor-Joining

---

#### PHYLO-NEWICK-001: Newick I/O

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.ToNewick(...)` / `ParseNewick(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ToNewick(treeNode)` | PhylogeneticAnalyzer | Export |
| `ParseNewick(newickString)` | PhylogeneticAnalyzer | Parse |

---

#### PHYLO-COMP-001: Tree Comparison

| Field | Value |
|------|----------|
| **Canonical** | `PhylogeneticAnalyzer.RobinsonFouldsDistance(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `RobinsonFouldsDistance(tree1, tree2)` | PhylogeneticAnalyzer | Topology |
| `FindMRCA(root, taxon1, taxon2)` | PhylogeneticAnalyzer | MRCA |
| `PatristicDistance(root, t1, t2)` | PhylogeneticAnalyzer | Path length |

---

### 9. Population Genetics (5 units)

#### POP-FREQ-001: Allele Frequencies

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(...)` |
| **Complexity** | O(n) |
| **Invariant** | p + q = 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateAlleleFrequencies(hom_maj, het, hom_min)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateMAF(genotypes)` | PopulationGeneticsAnalyzer | MAF |
| `FilterByMAF(variants, min, max)` | PopulationGeneticsAnalyzer | Filter |

---

#### POP-DIV-001: Diversity Statistics

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateDiversityStatistics(...)` |
| **Complexity** | O(n² × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateNucleotideDiversity(seqs)` | PopulationGeneticsAnalyzer | π |
| `CalculateWattersonTheta(segSites, n)` | PopulationGeneticsAnalyzer | θ |
| `CalculateTajimasD(pi, theta, S)` | PopulationGeneticsAnalyzer | D |
| `CalculateDiversityStatistics(seqs)` | PopulationGeneticsAnalyzer | All |

---

#### POP-HW-001: Hardy-Weinberg Test

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.TestHardyWeinberg(...)` |
| **Complexity** | O(1) per variant |
| **Invariant** | Expected: p², 2pq, q² |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `TestHardyWeinberg(variantId, counts)` | PopulationGeneticsAnalyzer | Canonical |

---

#### POP-FST-001: F-Statistics

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateFst(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ Fst ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateFst(pop1, pop2)` | PopulationGeneticsAnalyzer | Canonical |
| `CalculateFStatistics(variantData)` | PopulationGeneticsAnalyzer | Fis, Fit, Fst |

---

#### POP-LD-001: Linkage Disequilibrium

| Field | Value |
|------|----------|
| **Canonical** | `PopulationGeneticsAnalyzer.CalculateLD(...)` |
| **Complexity** | O(n) per pair |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateLD(var1, var2, genotypes)` | PopulationGeneticsAnalyzer | D', r² |
| `FindHaplotypeBlocks(variants)` | PopulationGeneticsAnalyzer | Blocks |

---

### 10. Chromosome Analysis (5 units)

#### CHROM-TELO-001: Telomere Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeTelomeres(...)` |
| **Complexity** | O(n) |
| **Constant** | Human telomere repeat: TTAGGG |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeTelomeres(chrName, seq, repeat, ...)` | ChromosomeAnalyzer | Canonical |
| `EstimateTelomereLengthFromTSRatio(tsRatio)` | ChromosomeAnalyzer | qPCR estimate |

---

#### CHROM-CENT-001: Centromere Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeCentromere(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeCentromere(chrName, seq, windowSize)` | ChromosomeAnalyzer | Canonical |

---

#### CHROM-KARYO-001: Karyotype Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.AnalyzeKaryotype(...)` |
| **Complexity** | O(k) where k=chromosomes |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AnalyzeKaryotype(chromosomes, ploidy)` | ChromosomeAnalyzer | Canonical |
| `DetectPloidy(depths, expected)` | ChromosomeAnalyzer | Ploidy |

---

#### CHROM-ANEU-001: Aneuploidy Detection

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.DetectAneuploidy(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectAneuploidy(depthData, medianDepth)` | ChromosomeAnalyzer | Canonical |
| `IdentifyWholeChromosomeAneuploidy(cnStates)` | ChromosomeAnalyzer | Classification |

---

#### CHROM-SYNT-001: Synteny Analysis

| Field | Value |
|------|----------|
| **Canonical** | `ChromosomeAnalyzer.FindSyntenyBlocks(...)` |
| **Complexity** | O(n log n) — requires verification |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSyntenyBlocks(orthologPairs, minGenes)` | ChromosomeAnalyzer | Canonical |
| `DetectRearrangements(syntenyBlocks)` | ChromosomeAnalyzer | Rearrangements |

---

### 11. Metagenomics (5 units)

#### META-CLASS-001: Taxonomic Classification

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.ClassifyReads(...)` |
| **Complexity** | O(n × m) where n=reads, m=read length |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClassifyReads(reads, kmerDB, k)` | MetagenomicsAnalyzer | Canonical |
| `BuildKmerDatabase(refGenomes, k)` | MetagenomicsAnalyzer | DB construction |

---

#### META-PROF-001: Taxonomic Profile

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.GenerateTaxonomicProfile(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GenerateTaxonomicProfile(classifications)` | MetagenomicsAnalyzer | Canonical |

---

#### META-ALPHA-001: Alpha Diversity

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.CalculateAlphaDiversity(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateAlphaDiversity(abundances)` | MetagenomicsAnalyzer | Canonical |

**Indices:** Shannon, Simpson, Inverse Simpson, Chao1, Pielou's evenness

---

#### META-BETA-001: Beta Diversity

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.CalculateBetaDiversity(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateBetaDiversity(sample1, sample2)` | MetagenomicsAnalyzer | Canonical |

**Metrics:** Bray-Curtis, Jaccard

---

#### META-BIN-001: Genome Binning

| Field | Value |
|------|----------|
| **Canonical** | `MetagenomicsAnalyzer.BinContigs(...)` |
| **Complexity** | O(n × k × i) where k=bins, i=iterations — needs verification |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `BinContigs(contigs, numBins, minBinSize)` | MetagenomicsAnalyzer | Canonical |

---

### 12. Codon Optimization (4 units)

#### CODON-OPT-001: Sequence Optimization

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.OptimizeSequence(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `OptimizeSequence(seq, organism, strategy)` | CodonOptimizer | Canonical |

**Strategies:** MaximizeCAI, BalancedOptimization, HarmonizeExpression, MinimizeSecondary, AvoidRareCodons

---

#### CODON-CAI-001: CAI Calculation

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.CalculateCAI(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ CAI ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCAI(codingSeq, codonTable)` | CodonOptimizer | Canonical |

---

#### CODON-RARE-001: Rare Codon Detection

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.FindRareCodons(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindRareCodons(seq, codonTable, threshold)` | CodonOptimizer | Canonical |

---

#### CODON-USAGE-001: Codon Usage Analysis

| Field | Value |
|------|----------|
| **Canonical** | `CodonOptimizer.CalculateCodonUsage(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateCodonUsage(seq)` | CodonOptimizer | Canonical |
| `CompareCodonUsage(seq1, seq2)` | CodonOptimizer | Comparison |

**Organism Tables:** E. coli K12, S. cerevisiae, H. sapiens

---

### 13. Translation (2 units)

#### TRANS-CODON-001: Codon Translation

| Field | Value |
|------|----------|
| **Canonical** | `GeneticCode.Translate(...)` |
| **Complexity** | O(1) per codon |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Translate(codon)` | GeneticCode | Canonical |
| `IsStartCodon(codon)` | GeneticCode | Check |
| `IsStopCodon(codon)` | GeneticCode | Check |
| `GetCodonsForAminoAcid(aa)` | GeneticCode | Reverse lookup |

**Genetic Codes:** Standard (1), Vertebrate Mito (2), Yeast Mito (3), Bacterial (11)

---

#### TRANS-PROT-001: Protein Translation

| Field | Value |
|------|----------|
| **Canonical** | `Translator.Translate(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Translate(dna, geneticCode)` | Translator | Canonical |

---

### 14. File I/O Parsers (7 units)

#### PARSE-FASTA-001: FASTA Parsing

| Field | Value |
|------|----------|
| **Canonical** | `FastaParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | FastaParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(fastaContent)` | FastaParser | Canonical |
| `ParseFile(filePath)` | FastaParser | File |
| `ParseFileAsync(filePath)` | FastaParser | Async |
| `ToFasta(entries, lineWidth)` | FastaParser | Export |
| `WriteFile(filePath, entries)` | FastaParser | Write |

**Edge Cases:**
- [ ] Empty content
- [ ] Multi-line sequences
- [ ] Special characters in headers
- [ ] Missing sequence after header

---

#### PARSE-FASTQ-001: FASTQ Parsing

| Field | Value |
|------|----------|
| **Canonical** | `FastqParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | FastqParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content, encoding)` | FastqParser | Canonical |
| `ParseFile(filePath, encoding)` | FastqParser | File |
| `CalculateStatistics(records)` | FastqParser | Stats |
| `FilterByQuality(records, minQ)` | FastqParser | Filter |

**Edge Cases:**
- [ ] Phred+33 vs Phred+64 encoding
- [ ] Malformed quality strings
- [ ] Empty records
- [ ] Auto-detect encoding

---

#### PARSE-BED-001: BED File Parsing

| Field | Value |
|------|----------|
| **Canonical** | `BedParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | BedParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content, format)` | BedParser | Canonical |
| `ParseFile(filePath, format)` | BedParser | File |
| `FilterByChrom(records, chrom)` | BedParser | Filter |
| `FilterByRegion(records, ...)` | BedParser | Filter |
| `MergeOverlapping(records)` | BedParser | Merge |
| `Intersect(records1, records2)` | BedParser | Set op |

**Edge Cases:**
- [ ] BED3 vs BED6 vs BED12 formats
- [ ] Block structures (BED12)
- [ ] Invalid coordinates

---

#### PARSE-VCF-001: VCF Parsing

| Field | Value |
|------|----------|
| **Canonical** | `VcfParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | VcfParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | VcfParser | Canonical |
| `ParseFile(filePath)` | VcfParser | File |
| `ParseWithHeader(content)` | VcfParser | With header |
| `GetVariantType(record)` | VcfParser | Classification |

**Edge Cases:**
- [ ] Multi-allelic variants
- [ ] Missing values (.)
- [ ] Complex INFO fields
- [ ] Sample genotypes

---

#### PARSE-GFF-001: GFF/GTF Parsing

| Field | Value |
|------|----------|
| **Canonical** | `GffParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | GffParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | GffParser | Canonical |
| `ParseFile(filePath)` | GffParser | File |
| `ToGff3(features)` | GffParser | Export |

---

#### PARSE-GENBANK-001: GenBank Parsing

| Field | Value |
|------|----------|
| **Canonical** | `GenBankParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | GenBankParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | GenBankParser | Canonical |
| `ParseFile(filePath)` | GenBankParser | File |
| `ExtractFeatures(record)` | GenBankParser | Features |

---

#### PARSE-EMBL-001: EMBL Parsing

| Field | Value |
|------|----------|
| **Canonical** | `EmblParser.Parse(...)` |
| **Complexity** | O(n) |
| **Class** | EmblParser |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Parse(content)` | EmblParser | Canonical |
| `ParseFile(filePath)` | EmblParser | File |

---

### 15. Sequence Complexity (3 units)

#### SEQ-COMPLEX-001: Linguistic Complexity

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateLinguisticComplexity(...)` |
| **Complexity** | O(n × k) |
| **Invariant** | 0 ≤ result ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateLinguisticComplexity(seq, maxLen)` | SequenceComplexity | Canonical |
| `CalculateLinguisticComplexity(string, maxLen)` | SequenceComplexity | String |
| `FindLowComplexityRegions(seq, window, threshold)` | SequenceComplexity | Regions |
| `MaskLowComplexity(seq, threshold)` | SequenceComplexity | Masking |

---

#### SEQ-ENTROPY-001: Shannon Entropy

| Field | Value |
|------|----------|
| **Canonical** | `SequenceComplexity.CalculateShannonEntropy(...)` |
| **Complexity** | O(n) |
| **Invariant** | 0 ≤ result ≤ 2 for DNA |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateShannonEntropy(sequence)` | SequenceComplexity | Canonical |
| `CalculateKmerEntropy(seq, k)` | SequenceComplexity | K-mer based |

---

#### SEQ-GCSKEW-001: GC Skew

| Field | Value |
|------|----------|
| **Canonical** | `GcSkewCalculator.CalculateGcSkew(...)` |
| **Complexity** | O(n) |
| **Invariant** | -1 ≤ result ≤ 1 |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateGcSkew(sequence)` | GcSkewCalculator | Canonical |
| `CalculateWindowedGcSkew(seq, window, step)` | GcSkewCalculator | Windowed |
| `CalculateCumulativeGcSkew(sequence)` | GcSkewCalculator | Cumulative |
| `FindOriginOfReplication(sequence)` | GcSkewCalculator | Origin detection |

---

### 16. RNA Secondary Structure (3 units)

#### RNA-STRUCT-001: Secondary Structure Prediction

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.Predict(...)` |
| **Complexity** | O(n³) |
| **Class** | RnaSecondaryStructure |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Predict(sequence)` | RnaSecondaryStructure | Canonical |
| `PredictWithConstraints(seq, constraints)` | RnaSecondaryStructure | Constrained |
| `ToDotBracket(structure)` | RnaSecondaryStructure | Notation |
| `FromDotBracket(notation)` | RnaSecondaryStructure | Parse |

---

#### RNA-STEMLOOP-001: Stem-Loop Detection

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.FindStemLoops(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindStemLoops(sequence, minStem, maxLoop)` | RnaSecondaryStructure | Canonical |
| `FindHairpins(sequence, params)` | RnaSecondaryStructure | Hairpins |
| `FindPseudoknots(sequence)` | RnaSecondaryStructure | Pseudoknots |

---

#### RNA-ENERGY-001: Free Energy Calculation

| Field | Value |
|------|----------|
| **Canonical** | `RnaSecondaryStructure.CalculateFreeEnergy(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateFreeEnergy(structure)` | RnaSecondaryStructure | Canonical |
| `CalculateStackingEnergy(bp1, bp2)` | RnaSecondaryStructure | Stacking |

---

### 17. MicroRNA Analysis (3 units)

#### MIRNA-SEED-001: Seed Sequence Analysis

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.GetSeedSequence(...)` |
| **Complexity** | O(1) |
| **Class** | MiRnaAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `GetSeedSequence(miRnaSequence)` | MiRnaAnalyzer | Canonical |
| `CreateMiRna(name, sequence)` | MiRnaAnalyzer | Factory |
| `CompareSeedRegions(mirna1, mirna2)` | MiRnaAnalyzer | Compare |

---

#### MIRNA-TARGET-001: Target Site Prediction

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindTargetSites(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindTargetSites(mRna, miRna, minScore)` | MiRnaAnalyzer | Canonical |
| `ScoreTargetSite(site)` | MiRnaAnalyzer | Scoring |

**Site Types:**
- [ ] 8mer, 7mer-m8, 7mer-A1, 6mer
- [ ] Supplementary pairing
- [ ] Centered sites

---

#### MIRNA-PRECURSOR-001: Pre-miRNA Detection

| Field | Value |
|------|----------|
| **Canonical** | `MiRnaAnalyzer.FindPreMiRnas(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindPreMiRnas(sequence)` | MiRnaAnalyzer | Canonical |
| `ValidateHairpin(structure)` | MiRnaAnalyzer | Validation |

---

### 18. Splice Site Prediction (3 units)

#### SPLICE-DONOR-001: Donor Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `SpliceSitePredictor.FindDonorSites(...)` |
| **Complexity** | O(n) |
| **Class** | SpliceSitePredictor |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDonorSites(sequence, minScore)` | SpliceSitePredictor | Canonical |
| `ScoreDonorSite(context)` | SpliceSitePredictor | Scoring |

---

#### SPLICE-ACCEPTOR-001: Acceptor Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `SpliceSitePredictor.FindAcceptorSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindAcceptorSites(sequence, minScore)` | SpliceSitePredictor | Canonical |
| `ScoreAcceptorSite(context)` | SpliceSitePredictor | Scoring |

---

#### SPLICE-PREDICT-001: Gene Structure Prediction

| Field | Value |
|------|----------|
| **Canonical** | `SpliceSitePredictor.PredictGeneStructure(...)` |
| **Complexity** | O(n²) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictGeneStructure(sequence)` | SpliceSitePredictor | Canonical |
| `FindIntrons(sequence)` | SpliceSitePredictor | Introns |
| `FindExons(sequence)` | SpliceSitePredictor | Exons |

---

### 19. Protein Disorder Prediction (2 units)

#### DISORDER-PRED-001: Disorder Prediction

| Field | Value |
|------|----------|
| **Canonical** | `DisorderPredictor.PredictDisorder(...)` |
| **Complexity** | O(n) |
| **Class** | DisorderPredictor |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictDisorder(sequence, windowSize, threshold)` | DisorderPredictor | Canonical |
| `CalculateDisorderScore(window)` | DisorderPredictor | Score |
| `CalculateHydropathy(sequence)` | DisorderPredictor | Hydropathy |

---

#### DISORDER-REGION-001: Disordered Region Detection

| Field | Value |
|------|----------|
| **Canonical** | `DisorderPredictor.IdentifyDisorderedRegions(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `IdentifyDisorderedRegions(predictions, minLen)` | DisorderPredictor | Canonical |
| `ClassifyRegionType(region)` | DisorderPredictor | Classification |

---

### 20. Protein Motif Finding (3 units)

#### PROTMOTIF-FIND-001: Motif Search

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.FindMotifs(...)` |
| **Complexity** | O(n × m) |
| **Class** | ProteinMotifFinder |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMotifs(sequence, patterns)` | ProteinMotifFinder | Canonical |
| `FindAllKnownMotifs(sequence)` | ProteinMotifFinder | All patterns |
| `ScanForPattern(sequence, pattern)` | ProteinMotifFinder | Single |

---

#### PROTMOTIF-PROSITE-001: PROSITE Pattern Matching

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.MatchPrositePattern(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `MatchPrositePattern(sequence, pattern)` | ProteinMotifFinder | Canonical |
| `ParsePrositePattern(pattern)` | ProteinMotifFinder | Parse |

**Common Patterns:**
- [ ] N-glycosylation (PS00001)
- [ ] Phosphorylation sites
- [ ] Zinc fingers
- [ ] Signal peptides

---

#### PROTMOTIF-DOMAIN-001: Domain Prediction

| Field | Value |
|------|----------|
| **Canonical** | `ProteinMotifFinder.PredictDomains(...)` |
| **Complexity** | O(n × d) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `PredictDomains(sequence)` | ProteinMotifFinder | Canonical |
| `PredictSignalPeptide(sequence)` | ProteinMotifFinder | Signal |

---

### 21. Epigenetics Analysis (3 units)

#### EPIGEN-CPG-001: CpG Site Detection

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindCpGSites(...)` |
| **Complexity** | O(n) |
| **Class** | EpigeneticsAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindCpGSites(sequence)` | EpigeneticsAnalyzer | Canonical |
| `FindCpGIslands(sequence, params)` | EpigeneticsAnalyzer | Islands |
| `CalculateObservedExpectedCpG(sequence)` | EpigeneticsAnalyzer | O/E ratio |

---

#### EPIGEN-METHYL-001: Methylation Analysis

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindMethylationSites(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindMethylationSites(sequence)` | EpigeneticsAnalyzer | Canonical |
| `CalculateMethylationProfile(sites)` | EpigeneticsAnalyzer | Profile |
| `GetMethylationContext(site)` | EpigeneticsAnalyzer | Context (CpG/CHG/CHH) |

---

#### EPIGEN-DMR-001: Differentially Methylated Regions

| Field | Value |
|------|----------|
| **Canonical** | `EpigeneticsAnalyzer.FindDMRs(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDMRs(profile1, profile2, threshold)` | EpigeneticsAnalyzer | Canonical |
| `AnnotateDMRs(dmrs, annotations)` | EpigeneticsAnalyzer | Annotate |

---

### 22. Variant Calling (4 units)

#### VARIANT-CALL-001: Variant Detection

| Field | Value |
|------|----------|
| **Canonical** | `VariantCaller.CallVariants(...)` |
| **Complexity** | O(n × m) |
| **Class** | VariantCaller |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CallVariants(reference, query)` | VariantCaller | Canonical |
| `CallVariantsFromAlignment(aligned1, aligned2)` | VariantCaller | From alignment |
| `ClassifyVariant(variant)` | VariantCaller | Classification |

---

#### VARIANT-SNP-001: SNP Detection

| Field | Value |
|------|----------|
| **Canonical** | `VariantCaller.FindSnps(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSnps(reference, query)` | VariantCaller | Canonical |
| `FindSnpsDirect(ref, query)` | VariantCaller | Direct (no alignment) |

---

#### VARIANT-INDEL-001: Indel Detection

| Field | Value |
|------|----------|
| **Canonical** | `VariantCaller.FindInsertions/FindDeletions(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindInsertions(reference, query)` | VariantCaller | Insertions |
| `FindDeletions(reference, query)` | VariantCaller | Deletions |

---

#### VARIANT-ANNOT-001: Variant Annotation

| Field | Value |
|------|----------|
| **Canonical** | `VariantAnnotator.Annotate(...)` |
| **Complexity** | O(v × g) |
| **Class** | VariantAnnotator |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `Annotate(variants, annotations)` | VariantAnnotator | Canonical |
| `PredictFunctionalImpact(variant)` | VariantAnnotator | Impact |

---

### 23. Structural Variant Analysis (3 units)

#### SV-DETECT-001: SV Detection

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectSVs(...)` |
| **Complexity** | O(n log n) |
| **Class** | StructuralVariantAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectSVs(readPairs, splitReads)` | StructuralVariantAnalyzer | Canonical |
| `FindDiscordantPairs(readPairs, params)` | StructuralVariantAnalyzer | Discordant |
| `ClassifySV(sv)` | StructuralVariantAnalyzer | Classification |

**SV Types:**
- [ ] Deletion, Duplication, Inversion
- [ ] Insertion, Translocation
- [ ] Complex rearrangements

---

#### SV-BREAKPOINT-001: Breakpoint Detection

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.FindBreakpoints(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindBreakpoints(splitReads)` | StructuralVariantAnalyzer | Canonical |
| `RefineBreakpoint(region, reads)` | StructuralVariantAnalyzer | Refinement |

---

#### SV-CNV-001: Copy Number Variation

| Field | Value |
|------|----------|
| **Canonical** | `StructuralVariantAnalyzer.DetectCNV(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectCNV(depthData, window)` | StructuralVariantAnalyzer | Canonical |
| `SegmentCopyNumber(logRatios)` | StructuralVariantAnalyzer | Segmentation |

---

### 24. Sequence Assembly (3 units)

#### ASSEMBLY-OLC-001: Overlap-Layout-Consensus

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleOLC(...)` |
| **Complexity** | O(n² × m) |
| **Class** | SequenceAssembler |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AssembleOLC(reads, params)` | SequenceAssembler | Canonical |
| `FindAllOverlaps(reads, minOverlap, minId)` | SequenceAssembler | Overlaps |

---

#### ASSEMBLY-DBG-001: De Bruijn Graph Assembly

| Field | Value |
|------|----------|
| **Canonical** | `SequenceAssembler.AssembleDeBruijn(...)` |
| **Complexity** | O(n × k) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `AssembleDeBruijn(reads, params)` | SequenceAssembler | Canonical |
| `BuildDeBruijnGraph(reads, k)` | SequenceAssembler | Graph construction |

---

#### ASSEMBLY-STATS-001: Assembly Statistics

| Field | Value |
|------|----------|
| **Canonical** | `GenomeAssemblyAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n log n) |
| **Class** | GenomeAssemblyAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(contigs)` | GenomeAssemblyAnalyzer | Canonical |
| `CalculateN50(contigs)` | GenomeAssemblyAnalyzer | N50 |
| `CalculateNx(contigs, threshold)` | GenomeAssemblyAnalyzer | Nx/Lx |
| `FindGaps(sequence)` | GenomeAssemblyAnalyzer | Gap detection |

---

### 25. Transcriptome Analysis (3 units)

#### TRANS-EXPR-001: Expression Quantification

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.CalculateTPM(...)` |
| **Complexity** | O(n) |
| **Class** | TranscriptomeAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateTPM(geneCounts)` | TranscriptomeAnalyzer | TPM |
| `CalculateFPKM(count, length, total)` | TranscriptomeAnalyzer | FPKM |
| `QuantileNormalize(samples)` | TranscriptomeAnalyzer | Normalization |

---

#### TRANS-DIFF-001: Differential Expression

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.FindDifferentiallyExpressed(...)` |
| **Complexity** | O(g × s) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindDifferentiallyExpressed(cond1, cond2, alpha)` | TranscriptomeAnalyzer | Canonical |
| `CalculateFoldChange(expr1, expr2)` | TranscriptomeAnalyzer | Fold change |

---

#### TRANS-SPLICE-001: Alternative Splicing

| Field | Value |
|------|----------|
| **Canonical** | `TranscriptomeAnalyzer.DetectAlternativeSplicing(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectAlternativeSplicing(isoforms)` | TranscriptomeAnalyzer | Canonical |
| `CalculatePSI(event, reads)` | TranscriptomeAnalyzer | Percent spliced in |

---

### 26. Comparative Genomics (3 units)

#### COMPGEN-SYNTENY-001: Synteny Detection

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindSyntenicBlocks(...)` |
| **Complexity** | O(n²) |
| **Class** | ComparativeGenomics |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindSyntenicBlocks(genes1, genes2, orthologs)` | ComparativeGenomics | Canonical |
| `VisualizeSynteny(blocks)` | ComparativeGenomics | Visualization |

---

#### COMPGEN-ORTHO-001: Ortholog Identification

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.FindOrthologs(...)` |
| **Complexity** | O(n × m) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `FindOrthologs(genes1, genes2, minIdentity)` | ComparativeGenomics | Canonical |
| `FindParalogs(genes, minIdentity)` | ComparativeGenomics | Paralogs |

---

#### COMPGEN-REARR-001: Genome Rearrangements

| Field | Value |
|------|----------|
| **Canonical** | `ComparativeGenomics.DetectRearrangements(...)` |
| **Complexity** | O(n log n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DetectRearrangements(blocks)` | ComparativeGenomics | Canonical |
| `ClassifyRearrangement(event)` | ComparativeGenomics | Classification |

---

### 27. Pan-Genome Analysis (2 units)

#### PANGEN-CORE-001: Core/Accessory Genome

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.ConstructPanGenome(...)` |
| **Complexity** | O(g² × s) |
| **Class** | PanGenomeAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ConstructPanGenome(genomes, idThreshold)` | PanGenomeAnalyzer | Canonical |
| `IdentifyCoreGenes(clusters, threshold)` | PanGenomeAnalyzer | Core genes |

---

#### PANGEN-CLUSTER-001: Gene Clustering

| Field | Value |
|------|----------|
| **Canonical** | `PanGenomeAnalyzer.ClusterGenes(...)` |
| **Complexity** | O(g² × s) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ClusterGenes(genomes, idThreshold)` | PanGenomeAnalyzer | Canonical |
| `GeneratePresenceAbsenceMatrix(clusters)` | PanGenomeAnalyzer | Matrix |

---

### 28. Quality Score Analysis (2 units)

#### QUALITY-PHRED-001: Phred Score Handling

| Field | Value |
|------|----------|
| **Canonical** | `QualityScoreAnalyzer.ParseQualityString(...)` |
| **Complexity** | O(n) |
| **Class** | QualityScoreAnalyzer |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ParseQualityString(qualStr, encoding)` | QualityScoreAnalyzer | Canonical |
| `ToQualityString(scores, encoding)` | QualityScoreAnalyzer | Export |
| `ConvertEncoding(qualStr, from, to)` | QualityScoreAnalyzer | Convert |

---

#### QUALITY-STATS-001: Quality Statistics

| Field | Value |
|------|----------|
| **Canonical** | `QualityScoreAnalyzer.CalculateStatistics(...)` |
| **Complexity** | O(n) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `CalculateStatistics(scores)` | QualityScoreAnalyzer | Canonical |
| `CalculateQ30Percentage(scores)` | QualityScoreAnalyzer | Q30 |

---

### 29. Probe Design (2 units)

#### PROBE-DESIGN-001: Hybridization Probe Design

| Field | Value |
|------|----------|
| **Canonical** | `ProbeDesigner.DesignProbes(...)` |
| **Complexity** | O(n²) |
| **Class** | ProbeDesigner |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `DesignProbes(template, params)` | ProbeDesigner | Canonical |
| `DesignTilingProbes(template, overlap)` | ProbeDesigner | Tiling |
| `ScoreProbe(sequence, params)` | ProbeDesigner | Scoring |

---

#### PROBE-VALID-001: Probe Validation

| Field | Value |
|------|----------|
| **Canonical** | `ProbeDesigner.ValidateProbe(...)` |
| **Complexity** | O(n × g) |

**Methods:**
| Method | Class | Type |
|-------|-------|-----|
| `ValidateProbe(probe, genome)` | ProbeDesigner | Canonical |
| `CheckSpecificity(probe, database)` | ProbeDesigner | Specificity |

---

## Appendix A: Method Index

| Method | Test Unit ID |
|--------|--------------|
| `AssembleDeBruijn` | ASSEMBLY-DBG-001 |
| `AssembleOLC` | ASSEMBLY-OLC-001 |
| `CalculateAlphaDiversity` | META-ALPHA-001 |
| `CalculateAlleleFrequencies` | POP-FREQ-001 |
| `CalculateBetaDiversity` | META-BETA-001 |
| `CalculateCAI` | CODON-CAI-001 |
| `CalculateCodonUsage` | CODON-USAGE-001 |
| `CalculateCumulativeGcSkew` | SEQ-GCSKEW-001 |
| `CalculateDistanceMatrix` | PHYLO-DIST-001 |
| `CalculateDiversityStatistics` | POP-DIV-001 |
| `CalculateFoldChange` | TRANS-DIFF-001 |
| `CalculateFPKM` | TRANS-EXPR-001 |
| `CalculateFreeEnergy` | RNA-ENERGY-001 |
| `CalculateFst` | POP-FST-001 |
| `CalculateGcContent` | SEQ-GC-001 |
| `CalculateGcFraction` | SEQ-GC-001 |
| `CalculateGcSkew` | SEQ-GCSKEW-001 |
| `CalculateKmerEntropy` | KMER-FREQ-001 |
| `CalculateLD` | POP-LD-001 |
| `CalculateLinguisticComplexity` | SEQ-COMPLEX-001 |
| `CalculateMeltingTemperature` | PRIMER-TM-001 |
| `CalculateMethylationProfile` | EPIGEN-METHYL-001 |
| `CalculateN50` | ASSEMBLY-STATS-001 |
| `CalculateNucleotideDiversity` | POP-DIV-001 |
| `CalculateObservedExpectedCpG` | EPIGEN-CPG-001 |
| `CalculatePSI` | TRANS-SPLICE-001 |
| `CalculateQ30Percentage` | QUALITY-STATS-001 |
| `CalculateShannonEntropy` | SEQ-ENTROPY-001 |
| `CalculateSpecificityScore` | CRISPR-OFF-001 |
| `CalculateStatistics` | ASSEMBLY-STATS-001, QUALITY-STATS-001 |
| `CalculateTajimasD` | POP-DIV-001 |
| `CalculateTPM` | TRANS-EXPR-001 |
| `CalculateWattersonTheta` | POP-DIV-001 |
| `CalculateWindowedGcSkew` | SEQ-GCSKEW-001 |
| `CallVariants` | VARIANT-CALL-001 |
| `CallVariantsFromAlignment` | VARIANT-CALL-001 |
| `CheckSpecificity` | PROBE-VALID-001 |
| `ClassifyReads` | META-CLASS-001 |
| `ClusterGenes` | PANGEN-CLUSTER-001 |
| `Complement` | SEQ-COMP-001 |
| `ConstructPanGenome` | PANGEN-CORE-001 |
| `Contains` | PAT-EXACT-001 |
| `ConvertEncoding` | QUALITY-PHRED-001 |
| `CountKmers` | KMER-COUNT-001 |
| `CreateMiRna` | MIRNA-SEED-001 |
| `CreatePwm` | PAT-PWM-001 |
| `DesignGuideRnas` | CRISPR-GUIDE-001 |
| `DesignPrimers` | PRIMER-DESIGN-001 |
| `DesignProbes` | PROBE-DESIGN-001 |
| `DesignTilingProbes` | PROBE-DESIGN-001 |
| `DetectAlternativeSplicing` | TRANS-SPLICE-001 |
| `DetectAneuploidy` | CHROM-ANEU-001 |
| `DetectCNV` | SV-CNV-001 |
| `DetectPloidy` | CHROM-KARYO-001 |
| `DetectRearrangements` | COMPGEN-REARR-001 |
| `DetectSVs` | SV-DETECT-001 |
| `Digest` | RESTR-DIGEST-001 |
| `EditDistance` | PAT-APPROX-002 |
| `EvaluateGuideRna` | CRISPR-GUIDE-001 |
| `EvaluatePrimer` | PRIMER-DESIGN-001 |
| `FilterByQuality` | PARSE-FASTQ-001 |
| `FindAcceptorSites` | SPLICE-ACCEPTOR-001 |
| `FindAllOccurrences` | PAT-EXACT-001 |
| `FindAllOverlaps` | ASSEMBLY-OLC-001 |
| `FindBreakpoints` | SV-BREAKPOINT-001 |
| `FindClumps` | KMER-FIND-001 |
| `FindCpGIslands` | EPIGEN-CPG-001 |
| `FindCpGSites` | EPIGEN-CPG-001 |
| `FindDegenerateMotif` | PAT-IUPAC-001 |
| `FindDifferentiallyExpressed` | TRANS-DIFF-001 |
| `FindDirectRepeats` | REP-DIRECT-001 |
| `FindDiscordantPairs` | SV-DETECT-001 |
| `FindDMRs` | EPIGEN-DMR-001 |
| `FindDonorSites` | SPLICE-DONOR-001 |
| `FindExactMotif` | PAT-EXACT-001 |
| `FindExons` | SPLICE-PREDICT-001 |
| `FindHaplotypeBlocks` | POP-LD-001 |
| `FindIntrons` | SPLICE-PREDICT-001 |
| `FindInvertedRepeats` | REP-INV-001 |
| `FindLowComplexityRegions` | SEQ-COMPLEX-001 |
| `FindMethylationSites` | EPIGEN-METHYL-001 |
| `FindMicrosatellites` | REP-STR-001 |
| `FindMostFrequentKmers` | KMER-FIND-001 |
| `FindMotif` | PAT-EXACT-001 |
| `FindMotifs` | PROTMOTIF-FIND-001 |
| `FindOffTargets` | CRISPR-OFF-001 |
| `FindOrfs` | ANNOT-ORF-001 |
| `FindOriginOfReplication` | SEQ-GCSKEW-001 |
| `FindOrthologs` | COMPGEN-ORTHO-001 |
| `FindPalindromes` | REP-PALIN-001 |
| `FindPamSites` | CRISPR-PAM-001 |
| `FindParalogs` | COMPGEN-ORTHO-001 |
| `FindPreMiRnas` | MIRNA-PRECURSOR-001 |
| `FindPromoterMotifs` | ANNOT-PROM-001 |
| `FindRareCodons` | CODON-RARE-001 |
| `FindRibosomeBindingSites` | ANNOT-GENE-001 |
| `FindSites` | RESTR-FIND-001 |
| `FindStemLoops` | RNA-STEMLOOP-001 |
| `FindSyntenicBlocks` | COMPGEN-SYNTENY-001 |
| `FindSyntenyBlocks` | CHROM-SYNT-001 |
| `FindTandemRepeats` | REP-TANDEM-001 |
| `FindTargetSites` | MIRNA-TARGET-001 |
| `FindUniqueKmers` | KMER-FIND-001 |
| `FindWithEdits` | PAT-APPROX-002 |
| `FindWithMismatches` | PAT-APPROX-001 |
| `GeneratePresenceAbsenceMatrix` | PANGEN-CLUSTER-001 |
| `GenerateTaxonomicProfile` | META-PROF-001 |
| `GetComplementBase` | SEQ-COMP-001 |
| `GetKmerFrequencies` | KMER-FREQ-001 |
| `GetKmerSpectrum` | KMER-FREQ-001 |
| `GetSeedSequence` | MIRNA-SEED-001 |
| `GlobalAlign` | ALIGN-GLOBAL-001 |
| `HammingDistance` | PAT-APPROX-001 |
| `HasHairpinPotential` | PRIMER-STRUCT-001 |
| `HasPrimerDimer` | PRIMER-STRUCT-001 |
| `IdentifyCoreGenes` | PANGEN-CORE-001 |
| `IdentifyDisorderedRegions` | DISORDER-REGION-001 |
| `IsStartCodon` | TRANS-CODON-001 |
| `IsStopCodon` | TRANS-CODON-001 |
| `IsValidDna` | SEQ-VALID-001 |
| `IsValidRna` | SEQ-VALID-001 |
| `LocalAlign` | ALIGN-LOCAL-001 |
| `MaskLowComplexity` | SEQ-COMPLEX-001 |
| `MatchesIupac` | PAT-IUPAC-001 |
| `MatchPrositePattern` | PROTMOTIF-PROSITE-001 |
| `MergeOverlapping` | PARSE-BED-001 |
| `MultipleAlign` | ALIGN-MULTI-001 |
| `OptimizeSequence` | CODON-OPT-001 |
| `Parse` | PARSE-FASTA-001, PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001 |
| `ParseFile` | PARSE-FASTA-001, PARSE-FASTQ-001, PARSE-BED-001, PARSE-VCF-001 |
| `ParseGff3` | ANNOT-GFF-001 |
| `ParseNewick` | PHYLO-NEWICK-001 |
| `ParsePrositePattern` | PROTMOTIF-PROSITE-001 |
| `ParseQualityString` | QUALITY-PHRED-001 |
| `PatristicDistance` | PHYLO-COMP-001 |
| `Predict` | RNA-STRUCT-001 |
| `PredictDisorder` | DISORDER-PRED-001 |
| `PredictDomains` | PROTMOTIF-DOMAIN-001 |
| `PredictGenes` | ANNOT-GENE-001 |
| `PredictGeneStructure` | SPLICE-PREDICT-001 |
| `PredictSignalPeptide` | PROTMOTIF-DOMAIN-001 |
| `QuantileNormalize` | TRANS-EXPR-001 |
| `ReverseComplement` | SEQ-REVCOMP-001 |
| `RobinsonFouldsDistance` | PHYLO-COMP-001 |
| `ScanForPattern` | PROTMOTIF-FIND-001 |
| `ScanWithPwm` | PAT-PWM-001 |
| `ScoreDonorSite` | SPLICE-DONOR-001 |
| `ScoreAcceptorSite` | SPLICE-ACCEPTOR-001 |
| `ScoreProbe` | PROBE-DESIGN-001 |
| `ScoreTargetSite` | MIRNA-TARGET-001 |
| `SegmentCopyNumber` | SV-CNV-001 |
| `SemiGlobalAlign` | ALIGN-SEMI-001 |
| `TestHardyWeinberg` | POP-HW-001 |
| `ToDotBracket` | RNA-STRUCT-001 |
| `ToFasta` | PARSE-FASTA-001 |
| `ToGff3` | ANNOT-GFF-001 |
| `ToNewick` | PHYLO-NEWICK-001 |
| `ToQualityString` | QUALITY-PHRED-001 |
| `Translate` | TRANS-CODON-001, TRANS-PROT-001 |
| `TryGetComplement` | SEQ-COMP-001 |
| `TryGetReverseComplement` | SEQ-REVCOMP-001 |
| `ValidateProbe` | PROBE-VALID-001 |

---

## Appendix B: Complexity Notes

| Test Unit | Claimed | Verified | Notes |
|-----------|---------|----------|-------|
| REP-STR-001 | O(n²) | ⚠️ | Actually O(n × U × R), depends on parameters |
| REP-INV-001 | O(n²) | ⚠️ | O(n² × L) with maxLoopLength |
| CHROM-SYNT-001 | O(n log n) | ⚠️ | Has nested loops, needs verification |
| ALIGN-MULTI-001 | O(n² × m²) | ⚠️ | Progressive is typically O(n² × m) |
| META-BIN-001 | O(n) | ⚠️ | K-means is O(n × k × i), verify iterations |
| CRISPR-OFF-001 | O(n × m) | ⚠️ | May be exponential with high mismatches |
| RNA-STRUCT-001 | O(n³) | ✓ | Standard Nussinov/Zuker |
| ASSEMBLY-OLC-001 | O(n² × m) | ⚠️ | Depends on overlap detection method |
| PANGEN-CORE-001 | O(g² × s) | ⚠️ | All-vs-all comparison |
| SV-CNV-001 | O(n) | ⚠️ | Segmentation may be O(n log n) |

---

## Appendix C: Canonical Implementations

### GC Content (SEQ-GC-001)

```
Canonical: SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)
         ↑
    ┌────┴────┐
    │         │
CalculateGcContentFast(string)   CalculateGcFraction(ReadOnlySpan)
    │                                      │
    │                              CalculateGcFractionFast(string)
    │                                      │
    ├──────────────────────────────────────┤
    │                                      │
PrimerDesigner.CalculateGcContent   MetagenomicsAnalyzer.CalculateGcContent
ChromosomeAnalyzer (internal)       DnaSequence.GcContent
```

**Test Strategy:** Test canonical `CalculateGcContent(ReadOnlySpan)` thoroughly. Delegates need only smoke tests verifying they call canonical correctly.

---

### DNA Complement (SEQ-COMP-001 + SEQ-REVCOMP-001)

```
Canonical: SequenceExtensions.GetComplementBase(char)
         ↑
    ┌────┴────────────────────┐
    │                         │
TryGetComplement(Span)   TryGetReverseComplement(Span)
    │                         │
DnaSequence.Complement   DnaSequence.ReverseComplement
                              │
                    GetReverseComplementString(string)
```

**Test Strategy:** Test `GetComplementBase` for all nucleotides. Test Span methods for sequence operations. Instance methods are wrappers.

---

## Appendix D: Class Coverage Summary

| Class | Test Units | Status |
|-------|------------|--------|
| ApproximateMatcher | PAT-APPROX-001, PAT-APPROX-002 | ✓ |
| BedParser | PARSE-BED-001 | ✓ |
| ChromosomeAnalyzer | CHROM-TELO-001 to CHROM-SYNT-001 | ✓ |
| CodonOptimizer | CODON-OPT-001 to CODON-USAGE-001 | ✓ |
| ComparativeGenomics | COMPGEN-SYNTENY-001 to COMPGEN-REARR-001 | ✓ |
| CrisprDesigner | CRISPR-PAM-001 to CRISPR-OFF-001 | ✓ |
| DisorderPredictor | DISORDER-PRED-001, DISORDER-REGION-001 | ✓ |
| EmblParser | PARSE-EMBL-001 | ✓ |
| EpigeneticsAnalyzer | EPIGEN-CPG-001 to EPIGEN-DMR-001 | ✓ |
| FastaParser | PARSE-FASTA-001 | ✓ |
| FastqParser | PARSE-FASTQ-001 | ✓ |
| GcSkewCalculator | SEQ-GCSKEW-001 | ✓ |
| GenBankParser | PARSE-GENBANK-001 | ✓ |
| GenomeAnnotator | ANNOT-ORF-001 to ANNOT-GFF-001 | ✓ |
| GenomeAssemblyAnalyzer | ASSEMBLY-STATS-001 | ✓ |
| GenomicAnalyzer | REP-TANDEM-001, REP-PALIN-001 | ✓ |
| GffParser | PARSE-GFF-001 | ✓ |
| GeneticCode | TRANS-CODON-001 | ✓ |
| IupacHelper | PAT-IUPAC-001 | ✓ |
| KmerAnalyzer | KMER-COUNT-001 to KMER-FIND-001 | ✓ |
| MetagenomicsAnalyzer | META-CLASS-001 to META-BIN-001 | ✓ |
| MiRnaAnalyzer | MIRNA-SEED-001 to MIRNA-PRECURSOR-001 | ✓ |
| MotifFinder | PAT-PWM-001, PAT-IUPAC-001 | ✓ |
| PanGenomeAnalyzer | PANGEN-CORE-001, PANGEN-CLUSTER-001 | ✓ |
| PhylogeneticAnalyzer | PHYLO-DIST-001 to PHYLO-COMP-001 | ✓ |
| PopulationGeneticsAnalyzer | POP-FREQ-001 to POP-LD-001 | ✓ |
| PrimerDesigner | PRIMER-TM-001 to PRIMER-STRUCT-001 | ✓ |
| ProbeDesigner | PROBE-DESIGN-001, PROBE-VALID-001 | ✓ |
| ProteinMotifFinder | PROTMOTIF-FIND-001 to PROTMOTIF-DOMAIN-001 | ✓ |
| QualityScoreAnalyzer | QUALITY-PHRED-001, QUALITY-STATS-001 | ✓ |
| RepeatFinder | REP-STR-001 to REP-PALIN-001 | ✓ |
| RestrictionAnalyzer | RESTR-FIND-001, RESTR-DIGEST-001 | ✓ |
| RnaSecondaryStructure | RNA-STRUCT-001 to RNA-ENERGY-001 | ✓ |
| SequenceAligner | ALIGN-GLOBAL-001 to ALIGN-MULTI-001 | ✓ |
| SequenceAssembler | ASSEMBLY-OLC-001, ASSEMBLY-DBG-001 | ✓ |
| SequenceComplexity | SEQ-COMPLEX-001, SEQ-ENTROPY-001 | ✓ |
| SequenceExtensions | SEQ-GC-001 to SEQ-VALID-001 | ✓ |
| SpliceSitePredictor | SPLICE-DONOR-001 to SPLICE-PREDICT-001 | ✓ |
| StructuralVariantAnalyzer | SV-DETECT-001 to SV-CNV-001 | ✓ |
| TranscriptomeAnalyzer | TRANS-EXPR-001 to TRANS-SPLICE-001 | ✓ |
| Translator | TRANS-PROT-001 | ✓ |
| VariantAnnotator | VARIANT-ANNOT-001 | ✓ |
| VariantCaller | VARIANT-CALL-001 to VARIANT-INDEL-001 | ✓ |
| VcfParser | PARSE-VCF-001 | ✓ |

**Total Classes Covered: 42/42 (100%)**

---

*Generated: 2026-01-22*
*Checklist version: 2.1 (100% Coverage)*
