# Checklist 03: Fuzzing

**Priority:** P2  
**Date:** 2026-03-19  
**Total algorithms:** 234

---

## Description

Fuzzing подає випадкові, невалідні або граничні дані для виявлення крашів, зависань та необроблених винятків. Критично важливо для парсерів файлових форматів і точок валідації вхідних даних. У геноміці: сміттєві послідовності, обрізані файли, невалідні символи.

**Поточне покриття:** `SuffixTreeFuzzTests` (corruption headers). Нуль для Seqeron.Genomics.

**Стратегії фазингу:**
- **RB** = Random Bytes (випадкові байти)
- **TF** = Truncated Fields (обрізані поля)
- **MC** = Malformed Content (невалідний контент)
- **BE** = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty)
- **INJ** = Injection (спецсимволи, null bytes, unicode)
- **OVF** = Overflow (екстремальні довжини, вкладеності)

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Fuzz Strategy | Fuzz Targets |
|---|--------|-----------|------|:---:|-------------|
| 1 | ☑ | SEQ-GC-001 | Composition | BE, INJ | Empty string, single char, non-ACGT chars, null, unicode, extremely long |
| 2 | ☑ | SEQ-COMP-001 | Composition | BE, INJ | Non-DNA chars, empty, null, mixed case, unicode |
| 3 | ☑ | SEQ-REVCOMP-001 | Composition | BE, INJ | Non-DNA chars, empty, null, single char, unicode |
| 4 | ☑ | SEQ-VALID-001 | Composition | RB, INJ, BE | Non-ASCII, null bytes, mixed-case, unicode, extremely long, control chars |
| 5 | ☑ | SEQ-COMPLEX-001 | Composition | BE, RB | Empty, single char, all same nucleotide, very long, random bytes |
| 6 | ☑ | SEQ-ENTROPY-001 | Composition | BE, RB | Empty, single symbol, all same, very long, non-nucleotide chars |
| 7 | ☑ | SEQ-GCSKEW-001 | Composition | BE | Empty, single base, no G or C, alternating GC, extremely long |
| 8 | ☑ | PAT-EXACT-001 | Matching | BE, MC | Pattern > seq length, empty pattern, 1-char pattern, empty seq, pattern = seq |
| 9 | ☑ | PAT-APPROX-001 | Matching | BE, MC | Empty strings, unequal lengths, maxDist > len, non-DNA chars |
| 10 | ☑ | PAT-APPROX-002 | Matching | BE, MC | Empty strings, maxEdits negative, maxEdits > seq len, non-DNA chars |
| 11 | ☑ | PAT-IUPAC-001 | Matching | MC, INJ | Invalid IUPAC codes, mixed DNA/protein, numbers in pattern, unicode |
| 12 | ☑ | PAT-PWM-001 | Matching | BE, MC | Zero-length matrix, NaN weights, empty training set, single-seq training |
| 13 | ☑ | REP-STR-001 | Repeats | BE | minRepeats=0, minRepeats=1, maxUnitLen > seqLen, empty seq |
| 14 | ☑ | REP-TANDEM-001 | Repeats | BE | minReps=0, minUnitLen=0, maxUnitLen=1, empty seq, single char seq |
| 15 | ☑ | REP-INV-001 | Repeats | BE | minLen=0, minLen > seqLen/2, empty seq, no complement possibilities |
| 16 | ☑ | REP-DIRECT-001 | Repeats | BE | minLen=0, minLen=1, empty seq, all unique chars, all same char |
| 17 | ☑ | REP-PALIN-001 | Repeats | BE | minLen=0, minLen odd, maxLen > seqLen, empty seq |
| 18 | ☑ | CRISPR-PAM-001 | MolTools | MC, BE | Invalid PAM sequences, non-DNA characters, empty seq, PAM longer than seq |
| 19 | ☑ | CRISPR-GUIDE-001 | MolTools | MC, BE | No PAM sites in seq, guide length > seq, guide length 0, non-DNA |
| 20 | ☑ | CRISPR-OFF-001 | MolTools | BE, MC | Zero mismatch tolerance, empty guide, guide of all N's |
| 21 | ☑ | PRIMER-TM-001 | MolTools | BE, INJ | Empty seq, 1-bp, 100+ bp, all-N, non-DNA chars, null |
| 22 | ☑ | PRIMER-DESIGN-001 | MolTools | BE, MC | Seq shorter than min primer, GC range 0-0, Tm range inverted |
| 23 | ☑ | PRIMER-STRUCT-001 | MolTools | BE | Extremely short, extremely long, palindromic primer, all-G |
| 24 | ☑ | PROBE-DESIGN-001 | MolTools | BE, MC | Seq shorter than min probe, Tm range inverted, empty seq |
| 25 | ☑ | PROBE-VALID-001 | MolTools | MC, INJ | Non-DNA probe, extremely short, null, probe with N's |
| 26 | ☑ | RESTR-FIND-001 | MolTools | MC, BE | Empty enzyme list, unknown enzyme names, empty sequence, non-DNA |
| 27 | ☑ | RESTR-DIGEST-001 | MolTools | BE | No cut sites, circular with 0 sites, 100+ enzymes, empty seq |
| 28 | ☑ | ANNOT-ORF-001 | Annotation | BE, MC | No ATG in seq, all start no stop, minLen=0, non-DNA chars |
| 29 | ☑ | ANNOT-GENE-001 | Annotation | BE, MC | No RBS in seq, overlapping genes, empty seq, non-coding only |
| 30 | ☑ | ANNOT-PROM-001 | Annotation | BE | No -10 box, threshold=0, empty seq, seq shorter than motif |
| 31 | ☑ | ANNOT-GFF-001 | Annotation | TF, MC, INJ | Malformed GFF string, missing fields, invalid strand, tab injection |
| 32 | ☑ | KMER-COUNT-001 | K-mer | BE | k=0, k > seqLen, k=seqLen, empty seq, k=1 |
| 33 | ☑ | KMER-FREQ-001 | K-mer | BE | k=0, k > seqLen, empty seq, single char seq |
| 34 | ☑ | KMER-FIND-001 | K-mer | BE | k=0, minFreq=0, minFreq > possible, empty seq |
| 35 | ☑ | ALIGN-GLOBAL-001 | Alignment | BE, MC | Empty vs empty, empty vs seq, very long seqs, non-DNA, single char |
| 36 | ☑ | ALIGN-LOCAL-001 | Alignment | BE | Empty seqs, identical seqs, completely different seqs, 1-char seqs |
| 37 | ☑ | ALIGN-SEMI-001 | Alignment | BE | Empty seqs, no overlap, complete overlap, 1-char seqs |
| 38 | ☑ | ALIGN-MULTI-001 | Alignment | BE | Single seq, 2 seqs, 100+ seqs, empty list, seqs of length 1 |
| 39 | ☑ | PHYLO-DIST-001 | Phylogenetic | BE, MC | Identical seqs, empty seqs, single seq, non-DNA chars |
| 40 | ☑ | PHYLO-TREE-001 | Phylogenetic | BE | 0 seqs, 1 seq, 2 seqs, 100+ seqs, all identical seqs |
| 41 | ☑ | PHYLO-NEWICK-001 | Phylogenetic | TF, MC, INJ | Malformed Newick, unbalanced parens, missing semicolon, empty string |
| 42 | ☑ | PHYLO-COMP-001 | Phylogenetic | BE | Same tree vs same tree, different leaf sets, empty tree |
| 43 | ☑ | POP-FREQ-001 | PopGen | BE | 0 samples, 1 allele, all same allele, negative counts |
| 44 | ☑ | POP-DIV-001 | PopGen | BE | Single sequence, all identical, 0 seqs, very long seqs |
| 45 | ☑ | POP-HW-001 | PopGen | BE | 0 genotypes, all heterozygous, all homozygous, single genotype |
| 46 | ☑ | POP-FST-001 | PopGen | BE | Single population, identical populations, 0 samples |
| 47 | ☑ | POP-LD-001 | PopGen | BE | Single locus, monomorphic loci, 0 samples |
| 48 | ☑ | CHROM-TELO-001 | Chromosome | BE, MC | No TTAGGG, only TTAGGG, empty seq, non-DNA chars |
| 49 | ☑ | CHROM-CENT-001 | Chromosome | BE | Empty seq, all AT, all GC, extremely short |
| 50 | ☑ | CHROM-KARYO-001 | Chromosome | BE | 0 chromosomes, 100+ chromosomes, empty data |
| 51 | ☑ | CHROM-ANEU-001 | Chromosome | BE | 0 depth, negative depth, extremely high depth, empty data |
| 52 | ☑ | CHROM-SYNT-001 | Chromosome | BE | 0 ortholog pairs, self-comparison, empty genomes |
| 53 | ☑ | META-CLASS-001 | Metagenomics | MC, BE | Empty database, read with no match, non-DNA read, extremely short read |
| 54 | ☑ | META-PROF-001 | Metagenomics | BE | 0 reads, single read, all same taxon |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | BE | 0 species, 1 species, all equal abundance, single sample with 0 |
| 56 | ☑ | META-BETA-001 | Metagenomics | BE | Identical samples, empty samples, single species samples |
| 57 | ☑ | META-BIN-001 | Metagenomics | BE | 0 contigs, 1 contig, huge number of contigs, contigs with 0 length |
| 58 | ☑ | CODON-OPT-001 | Codon | MC, BE | Non-coding seq (no codons), seq len not divisible by 3, non-DNA chars |
| 59 | ☑ | CODON-CAI-001 | Codon | BE | Empty seq, seq with only stop codons, seq not divisible by 3 |
| 60 | ☑ | CODON-RARE-001 | Codon | BE | Empty seq, all rare, none rare, threshold=0, threshold=1 |
| 61 | ☑ | CODON-USAGE-001 | Codon | BE | Empty seq, seq len=1, seq len=2, very long seq |
| 62 | ☑ | TRANS-CODON-001 | Translation | MC | Invalid table ID, null code, non-standard table |
| 63 | ☑ | TRANS-PROT-001 | Translation | MC, BE | Non-DNA chars, empty, len=1, len=2, all stop codons, no stop codon |
| 64 | ☑ | PARSE-FASTA-001 | FileIO | RB, TF, MC, INJ | Missing >, empty seq, binary garbage, huge headers, null bytes, no newline |
| 65 | ☑ | PARSE-FASTQ-001 | FileIO | RB, TF, MC, INJ | Missing @/+, quality len mismatch, invalid quality chars, truncated |
| 66 | ☑ | PARSE-BED-001 | FileIO | TF, MC, BE, INJ | Negative coords, start > end, non-numeric fields, tab injection |
| 67 | ☑ | PARSE-VCF-001 | FileIO | TF, MC, INJ | Missing header, invalid genotypes, huge alleles, missing ##fileformat |
| 68 | ☑ | PARSE-GFF-001 | FileIO | TF, MC, BE, INJ | Invalid strand, negative pos, malformed attributes, url encoding |
| 69 | ☑ | PARSE-GENBANK-001 | FileIO | RB, TF, MC | Missing LOCUS, truncated features, invalid locations, no // terminator |
| 70 | ☑ | PARSE-EMBL-001 | FileIO | RB, TF, MC | Missing ID line, truncated seq, invalid feature keys, no // terminator |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | BE, MC | Empty, 1 base, all AU, non-RNA chars, extremely long |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | BE | Empty, no complement, minLoop=0, minStem > seqLen/2 |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | BE | Empty, no stacks, temperature=0, temperature negative |
| 74 | ☑ | MIRNA-SEED-001 | MiRNA | BE, MC | miRNA shorter than 8nt, empty, non-RNA, miRNA = target |
| 75 | ☑ | MIRNA-TARGET-001 | MiRNA | BE | No complementarity, empty 3'UTR, miRNA longer than target |
| 76 | ☑ | MIRNA-PRECURSOR-001 | MiRNA | BE, MC | Too short, no hairpin, all one base, non-RNA chars |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | BE, MC | No GT sites, seq < window, non-DNA, all GT |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | BE, MC | No AG sites, seq < window, non-DNA, all AG |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | BE | No donor/acceptor, single exon, all intronic, empty seq |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | BE, MC | Empty protein, 1 AA, all same AA, non-standard AAs, X amino acid |
| 81 | ☑ | DISORDER-REGION-001 | ProteinPred | BE | threshold=0, threshold=1, all ordered, all disordered, empty |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | MC, BE | Empty protein, non-amino acid chars, extremely short, all same char |
| 83 | ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | MC, INJ | Invalid PROSITE pattern syntax, regex injection, empty pattern |
| 84 | ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | BE | Protein shorter than min domain, empty, all X residues |
| 85 | ☑ | EPIGEN-CPG-001 | Epigenetics | BE | No CG dinucleotides, all CG, empty seq, windowSize=0, windowSize > seqLen |
| 86 | ☑ | ONCO-IMMUNE-001 | Oncology | BE, MC | 0 expression, all NaN, negative expression, empty gene set, unknown genes |
| 87 | ☑ | ONCO-SOMATIC-001 | Oncology | BE, MC | 0 depth, alt>depth, tumor=normal, empty reads, all-N bases |
| 88 | ☑ | ONCO-VAF-001 | Oncology | BE | ref=alt=0, alt>total, negative counts, huge depth |
| 89 | ☑ | ONCO-DRIVER-001 | Oncology | BE, MC | empty mutation list, all-passenger, duplicate hotspots, unknown gene |
| 90 | ☑ | ONCO-ARTIFACT-001 | Oncology | BE | zero depth, extreme strand bias, all-pass, all-fail |
| 91 | ☑ | ONCO-ANNOT-001 | Oncology | MC | out-of-bounds coords, ref≠genome, empty alt, unknown chrom |
| 92 | ☑ | ONCO-TMB-001 | Oncology | BE | zero mutations, panel size 0, negative size, huge counts |
| 93 | ☑ | ONCO-MSI-001 | Oncology | BE | zero loci, all-stable, all-unstable, single read |
| 94 | ☑ | ONCO-HRD-001 | Oncology | BE | no events, negative component, extreme counts |
| 95 | ☑ | ONCO-LOH-001 | Oncology | BE | BAF=0.5 everywhere, BAF=0/1, single SNP |
| 96 | ☑ | ONCO-SIG-001 | Oncology | MC | non-SNV variant, ambiguous base, no flanking context, empty set |
| 97 | ☑ | ONCO-SIG-002 | Oncology | BE | zero catalogue, singular signature matrix, negative counts |
| 98 | ☑ | ONCO-SIG-003 | Oncology | BE | 0 bootstrap reps, 1 mutation, fixed-seed extremes |
| 99 | ☑ | ONCO-SIG-004 | Oncology | BE | tied exposures, all-zero exposures |
| 100 | ☑ | ONCO-FUSION-001 | Oncology | MC | no chimeric reads, self-fusion, identical genes, empty reads |
| 101 | ☑ | ONCO-FUSION-002 | Oncology | BE | empty DB, exact match, near-miss partner |
| 102 | ☑ | ONCO-FUSION-003 | Oncology | BE | breakpoint at gene boundary, intronic, out-of-bounds |
| 103 | ☑ | ONCO-CNA-001 | Oncology | BE | log2=±∞, NaN ratio, single bin |
| 104 | ☑ | ONCO-CNA-002 | Oncology | BE | genome-wide amp, single-bin focal, threshold edge |
| 105 | ☑ | ONCO-CNA-003 | Oncology | BE | CN exactly 0, near-0, single-bin deletion |
| 106 | ☑ | ONCO-PURITY-001 | Oncology | BE | all-VAF=0, all-VAF=1, single variant |
| 107 | ☑ | ONCO-PLOIDY-001 | Oncology | BE | flat diploid, fully amplified, empty segments |
| 108 | ☑ | ONCO-CLONAL-001 | Oncology | BE | CCF at threshold, CCF>1, CCF=0 |
| 109 | ☑ | ONCO-NEO-001 | Oncology | MC | mutation at protein terminus, stop-gain, non-coding, length<8 |
| 110 | ☑ | ONCO-MHC-001 | Oncology | BE | IC50=0, IC50=∞, peptide too short/long, unknown allele |
| 111 | ☑ | ONCO-CTDNA-001 | Oncology | BE | zero tumor reads, 100% tumor, ultra-low depth |
| 112 | ☑ | ONCO-MRD-001 | Oncology | BE | no tracked variants, all-detected, single low-VAF read |
| 113 | ☑ | ONCO-CHIP-001 | Oncology | BE | empty gene list, VAF at band edges, all-CHIP |
| 114 | ☑ | ONCO-PHYLO-001 | Oncology | BE | single clone, identical clones, no shared mutations |
| 115 | ☑ | ONCO-CCF-001 | Oncology | BE | purity=0, CN=0, VAF>purity |
| 116 | ☑ | ONCO-HETERO-001 | Oncology | BE | single VAF, all-equal VAFs, VAF=0 |
| 117 | ☑ | ONCO-HLA-001 | Oncology | MC | ambiguous allele, homozygous locus, no coverage |
| 118 | ☑ | ONCO-ACTION-001 | Oncology | BE | no evidence, conflicting tiers, unknown drug |
| 119 | ☑ | ONCO-SV-001 | Oncology | BE | zero breakpoints, single breakpoint, genome-wide shattering |
| 120 | ☑ | ONCO-EXPR-001 | Oncology | BE | zero variance, single sample, all-equal expression, NaN |
| 121 | ☑ | SEQ-COMPOSITION-001 | Statistics | BE | empty, single base, all-N, lowercase, non-ACGT |
| 122 | ☑ | SEQ-DINUC-001 | Statistics | BE | length<2, single base, all-N, homopolymer |
| 123 | ☑ | SEQ-HYDRO-001 | Statistics | BE, MC | empty, non-amino-acid char, all-X, single residue |
| 124 | ☑ | SEQ-MW-001 | Statistics | BE, MC | empty, unknown residue, single base, very long |
| 125 | ☑ | SEQ-PI-001 | Statistics | BE | no charged residues, all-acidic, all-basic, empty |
| 126 | ☑ | SEQ-SECSTRUCT-001 | Statistics | BE, MC | empty, single residue, unknown residue |
| 127 | ☑ | SEQ-STATS-001 | Statistics | BE | empty, all-N, lowercase, mixed alphabet |
| 128 | ☑ | SEQ-SUMMARY-001 | Statistics | BE | empty, single base, very long, mixed case |
| 129 | ☑ | SEQ-THERMO-001 | Statistics | BE | empty, single base, all-AT, all-GC |
| 130 | ☑ | SEQ-TM-001 | Statistics | BE | empty, single base, all-AT, all-GC, non-ACGT |
| 131 | ☑ | COMPGEN-ANI-001 | Comparative | BE | identical genomes, no shared k-mers, empty genome, single base |
| 132 | ☑ | COMPGEN-CLUSTER-001 | Comparative | BE | single genome, no conserved genes, identical genomes |
| 133 | ☑ | COMPGEN-COMPARE-001 | Comparative | BE | empty genome, A vs A, disjoint genomes |
| 134 | ☑ | COMPGEN-DOTPLOT-001 | Comparative | BE | empty, single base, palindrome, repeat-rich |
| 135 | ☑ | COMPGEN-ORTHO-001 | Comparative | BE, MC | no homologs, all-identical, empty gene set |
| 136 | ☑ | COMPGEN-RBH-001 | Comparative | BE | no hits, ties, single gene each |
| 137 | ☑ | COMPGEN-REARR-001 | Comparative | BE | identical order, full reversal, single gene |
| 138 | ☑ | COMPGEN-REVERSAL-001 | Comparative | BE | identity permutation, full reversal, singleton |
| 139 | ☑ | COMPGEN-SYNTENY-001 | Comparative | BE | no synteny, whole-genome block, single anchor |
| 140 | ☑ | ASSEMBLY-CONSENSUS-001 | Assembly | BE | single read, conflicting reads, empty, all-N |
| 141 | ☑ | ASSEMBLY-CORRECT-001 | Assembly | BE | zero coverage, all-error reads, empty |
| 142 | ☑ | ASSEMBLY-COVER-001 | Assembly | BE | no reads, zero-length ref, single read |
| 143 | ☑ | ASSEMBLY-DBG-001 | Assembly | BE | k>read length, single read, all-identical reads |
| 144 | ☑ | ASSEMBLY-MERGE-001 | Assembly | BE | no overlap, full containment, identical contigs |
| 145 | ☑ | ASSEMBLY-OLC-001 | Assembly | BE | minOverlap>read length, single read, no overlaps |
| 146 | ☑ | ASSEMBLY-SCAFFOLD-001 | Assembly | BE | no links, conflicting links, single contig |
| 147 | ☑ | ASSEMBLY-STATS-001 | Assembly | BE | empty assembly, single contig, equal-length contigs |
| 148 | ☑ | ASSEMBLY-TRIM-001 | Assembly | BE | all-low-quality, all-high-quality, empty, quality cutoff 0 |
| 149 | ☑ | RNA-DOTBRACKET-001 | RnaStructure | MC | unbalanced brackets, illegal chars, empty, length mismatch |
| 150 | ☑ | RNA-HAIRPIN-001 | RnaStructure | BE | loop<minLoop, no stem, empty |
| 151 | ☑ | RNA-INVERT-001 | RnaStructure | BE | palindrome, no complementarity, single base |
| 152 | ☑ | RNA-MFE-001 | RnaStructure | BE | empty, single base, all-A (no pairs), homopolymer |
| 153 | ☑ | RNA-PAIR-001 | RnaStructure | BE, MC | non-RNA base, lowercase, gap char |
| 154 | ☑ | RNA-PARTITION-001 | RnaStructure | BE | empty, single base, all-unpaired |
| 155 | ☑ | RNA-PSEUDOKNOT-001 | RnaStructure | BE | nested-only, fully crossing, empty |
| 156 | ☑ | KMER-ASYNC-001 | K-mer | BE | empty, k>len, cancellation, huge input |
| 157 | ☑ | KMER-BOTH-001 | K-mer | BE | palindromic k-mer, k>len, empty |
| 158 | ☑ | KMER-DIST-001 | K-mer | BE | identical, disjoint, empty, different k |
| 159 | ☑ | KMER-GENERATE-001 | K-mer | BE | k=0, large k, non-DNA alphabet |
| 160 | ☑ | KMER-POSITIONS-001 | K-mer | BE | absent k-mer, overlapping, k>len |
| 161 | ☑ | KMER-STATS-001 | K-mer | BE | empty, k>len, homopolymer |
| 162 | ☑ | KMER-UNIQUE-001 | K-mer | BE | all-identical, all-distinct, empty |
| 163 | ☑ | PROTMOTIF-CC-001 | ProteinMotif | BE, MC | empty, non-amino-acid, single residue |
| 164 | ☑ | PROTMOTIF-COMMON-001 | ProteinMotif | BE | single sequence, no common motif, identical inputs |
| 165 | ☑ | PROTMOTIF-LC-001 | ProteinMotif | BE | homopolymer, high-complexity, empty |
| 166 | ☑ | PROTMOTIF-PATTERN-001 | ProteinMotif | MC | empty pattern, invalid regex, no match |
| 167 | ☑ | PROTMOTIF-SP-001 | ProteinMotif | BE | no signal, very short, all-hydrophobic |
| 168 | ☑ | PROTMOTIF-TM-001 | ProteinMotif | BE | all-hydrophilic, all-hydrophobic, short |
| 169 | ☑ | MOTIF-CONS-001 | Matching | BE, MC | unequal row lengths, empty alignment, single row |
| 170 | ☑ | MOTIF-DISCOVER-001 | Matching | BE | k=1, k>len, no recurrence |
| 171 | ☑ | MOTIF-GENERATE-001 | Matching | BE | empty counts, single column, ties |
| 172 | ☑ | MOTIF-REGULATORY-001 | Matching | BE | empty, no element, overlapping |
| 173 | ☑ | MOTIF-SHARED-001 | Matching | BE | single input, disjoint inputs, identical |
| 174 | ☑ | PAT-APPROX-003 | Matching | BE | pattern>text, empty, exact present |
| 175 | ☑ | GENOMIC-COMMON-001 | Analysis | BE | single input, disjoint, identical |
| 176 | ☑ | GENOMIC-MOTIFS-001 | Analysis | BE | empty, no motif, overlapping |
| 177 | ☑ | GENOMIC-ORF-001 | Analysis | BE | no ATG, no stop, nested ORFs |
| 178 | ☑ | GENOMIC-REPEAT-001 | Analysis | BE | no repeat, full repeat, minLen edge |
| 179 | ☑ | GENOMIC-SIMILARITY-001 | Analysis | BE | identical, disjoint, different lengths |
| 180 | ☑ | GENOMIC-TANDEM-001 | Analysis | BE | no tandem, full tandem, single unit |
| 181 | ☑ | EPIGEN-AGE-001 | Epigenetics | BE | no clock sites, all-methylated, empty |
| 182 | ☑ | EPIGEN-BISULF-001 | Epigenetics | BE | no C, all-C, all-methylated, empty |
| 183 | ☑ | EPIGEN-CHROM-001 | Epigenetics | BE | empty, single region, no marks |
| 184 | ☑ | EPIGEN-DMR-001 | Epigenetics | BE | identical methylomes, threshold edge, single site |
| 185 | ☑ | EPIGEN-METHYL-001 | Epigenetics | BE | no reads, all-methylated, zero coverage |
| 186 | ☑ | VARIANT-ANNOT-001 | Variants | MC | out-of-bounds, unknown consequence, empty |
| 187 | ☑ | VARIANT-CALL-001 | Variants | BE | zero depth, tumor=normal, all-N |
| 188 | ☑ | VARIANT-INDEL-001 | Variants | BE | length-0 indel, indel at edge, empty |
| 189 | ☑ | VARIANT-SNP-001 | Variants | BE | ref=alt, multi-allelic, zero depth |
| 190 | ☑ | PANGEN-CLUSTER-001 | PanGenome | BE | single gene, all-identical, identity edge |
| 191 | ☑ | PANGEN-CORE-001 | PanGenome | BE | single genome, disjoint genomes, empty |
| 192 | ☑ | PANGEN-HEAP-001 | PanGenome | BE | single genome, 2 genomes, identical |
| 193 | ☑ | PANGEN-MARKER-001 | PanGenome | BE | empty core, more markers than core |
| 194 | ☑ | META-FUNC-001 | Metagenomics | BE | empty, unknown genes, no DB hit |
| 195 | ☑ | META-PATHWAY-001 | Metagenomics | BE | no pathway genes, all-pathway, empty |
| 196 | ☑ | META-RESIST-001 | Metagenomics | BE | no resistance gene, empty DB, partial hit |
| 197 | ☑ | META-TAXA-001 | Metagenomics | BE | identical samples, single taxon, empty |
| 198 | ☑ | TRANS-DIFF-001 | Transcriptome | BE | A=B, zero counts, single replicate |
| 199 | ☑ | TRANS-EXPR-001 | Transcriptome | BE | zero reads, single transcript, all-multimapped |
| 200 | ☑ | TRANS-SPLICE-001 | Transcriptome | BE | single isoform, no junction reads, empty |
| 201 | ☑ | SV-BREAKPOINT-001 | StructuralVar | BE | identical, single breakpoint, no split reads |
| 202 | ☑ | SV-CNV-001 | StructuralVar | BE | ratio=1, zero coverage, single bin |
| 203 | ☑ | SV-DETECT-001 | StructuralVar | BE | identical genomes, overlapping SVs, empty |
| 204 | ☑ | DISORDER-LC-001 | ProteinPred | BE | homopolymer, high-complexity, empty |
| 205 | ☐ | DISORDER-MORF-001 | ProteinPred | BE | fully ordered, fully disordered, short |
| 206 | ☐ | DISORDER-PROPENSITY-001 | ProteinPred | BE, MC | empty, non-amino-acid, single residue |
| 207 | ☐ | POP-ANCESTRY-001 | PopGen | BE | single population, admixed 50/50, empty |
| 208 | ☐ | POP-ROH-001 | PopGen | BE | all-heterozygous, all-homozygous, minLen edge |
| 209 | ☐ | POP-SELECT-001 | PopGen | BE | neutral, fixed locus, single locus |
| 210 | ☐ | SEQ-ATSKEW-001 | Composition | BE | balanced AT, all-A, no AT, window edge |
| 211 | ☐ | SEQ-REPLICATION-001 | Composition | BE | flat skew, single minimum, circular wrap |
| 212 | ☐ | SEQ-RNACOMP-001 | Composition | BE, MC | non-RNA base, T instead of U, lowercase |
| 213 | ☐ | CODON-ENC-001 | Codon | BE | single codon, uniform usage, length not %3 |
| 214 | ☐ | CODON-RSCU-001 | Codon | BE | single codon, missing amino acid, empty |
| 215 | ☐ | CODON-STATS-001 | Codon | BE | length not %3, empty, single codon |
| 216 | ☐ | ANNOT-CODING-001 | Annotation | BE | random seq, perfect ORF, empty |
| 217 | ☐ | ANNOT-CODONUSAGE-001 | Annotation | BE | empty, single codon, length not %3 |
| 218 | ☐ | ANNOT-REPEAT-001 | Annotation | BE | no repeat, full repeat, minLen edge |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | BE, MC | empty, out-of-range char, wrong offset |
| 220 | ☐ | QUALITY-STATS-001 | Quality | BE | empty, single base, all-equal quality |
| 221 | ☐ | PHYLO-BOOT-001 | Phylogenetic | BE | single tree, 0 replicates, identical sequences |
| 222 | ☐ | PHYLO-STATS-001 | Phylogenetic | BE | single leaf, star tree, deep ladder |
| 223 | ☐ | TRANS-SIXFRAME-001 | Translation | BE | length not %3, empty, single base |
| 224 | ☐ | RESTR-FILTER-001 | MolTools | BE | no sites, all-pass, all-fail |
| 225 | ☐ | MIRNA-PAIR-001 | MiRNA | BE | no complementarity, perfect match, short miRNA |
| 226 | ☐ | ALIGN-STATS-001 | Alignment | BE | identical, no overlap, all-gap, empty |
| 227 | ☐ | SEQ-CODON-FREQ-001 | Statistics | BE | empty, len not multiple of 3, non-ACGT, lowercase, very long |
| 228 | ☐ | SEQ-COMPLEX-COMPRESS-001 | Complexity | BE | empty, single char, homopolymer, random, very long |
| 229 | ☐ | SEQ-COMPLEX-DUST-001 | Complexity | BE | empty, shorter than window, homopolymer, non-ACGT |
| 230 | ☐ | SEQ-COMPLEX-KMER-001 | Complexity | BE | empty, len < k, k=0, homopolymer, very long |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | BE | window > len, window=0, empty, single char |
| 232 | ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | BE | window > len, window=0, empty, single char |
| 233 | ☐ | SEQ-GC-ANALYSIS-001 | Composition | BE | empty, all-GC, all-AT, non-ACGT, very long |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | BE | window > len, window=0, empty, single char |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 234 |
| ☑ Complete | 204 |
| ☐ Not started | 30 |
| High-priority (parsers + validation) | 12 |
| Medium-priority (boundary inputs) | 45 |
| Lower-priority (algorithm-specific edge cases) | 29 |
