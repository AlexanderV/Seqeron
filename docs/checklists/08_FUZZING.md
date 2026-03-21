# Checklist 08: Fuzzing

**Priority:** P2  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Fuzzing подаёт случайные, невалидные или пограничные данные для обнаружения крашей, зависаний и необработанных исключений. Критически важно для парсеров файловых форматов и точек валидации входных данных. В геномике: мусорные последовательности, обрезанные файлы, невалидные символы.

**Текущее покрытие:** `SuffixTreeFuzzTests` (corruption headers). Ноль для Seqeron.Genomics.

**Стратегии фаззинга:**
- **RB** = Random Bytes (случайные байты)
- **TF** = Truncated Fields (обрезанные поля)
- **MC** = Malformed Content (невалидный контент)
- **BE** = Boundary Exploitation (пограничные значения: 0, -1, MaxInt, empty)
- **INJ** = Injection (спецсимволы, null bytes, unicode)
- **OVF** = Overflow (экстремальные длины, вложенности)

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Fuzz Strategy | Fuzz Targets |
|---|--------|-----------|------|:---:|-------------|
| 1 | ☐ | SEQ-GC-001 | Composition | BE, INJ | Empty string, single char, non-ACGT chars, null, unicode, extremely long |
| 2 | ☐ | SEQ-COMP-001 | Composition | BE, INJ | Non-DNA chars, empty, null, mixed case, unicode |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | BE, INJ | Non-DNA chars, empty, null, single char, unicode |
| 4 | ☐ | SEQ-VALID-001 | Composition | RB, INJ, BE | Non-ASCII, null bytes, mixed-case, unicode, extremely long, control chars |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | BE, RB | Empty, single char, all same nucleotide, very long, random bytes |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | BE, RB | Empty, single symbol, all same, very long, non-nucleotide chars |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | BE | Empty, single base, no G or C, alternating GC, extremely long |
| 8 | ☐ | PAT-EXACT-001 | Matching | BE, MC | Pattern > seq length, empty pattern, 1-char pattern, empty seq, pattern = seq |
| 9 | ☐ | PAT-APPROX-001 | Matching | BE, MC | Empty strings, unequal lengths, maxDist > len, non-DNA chars |
| 10 | ☐ | PAT-APPROX-002 | Matching | BE, MC | Empty strings, maxEdits negative, maxEdits > seq len, non-DNA chars |
| 11 | ☐ | PAT-IUPAC-001 | Matching | MC, INJ | Invalid IUPAC codes, mixed DNA/protein, numbers in pattern, unicode |
| 12 | ☐ | PAT-PWM-001 | Matching | BE, MC | Zero-length matrix, NaN weights, empty training set, single-seq training |
| 13 | ☐ | REP-STR-001 | Repeats | BE | minRepeats=0, minRepeats=1, maxUnitLen > seqLen, empty seq |
| 14 | ☐ | REP-TANDEM-001 | Repeats | BE | minReps=0, minUnitLen=0, maxUnitLen=1, empty seq, single char seq |
| 15 | ☐ | REP-INV-001 | Repeats | BE | minLen=0, minLen > seqLen/2, empty seq, no complement possibilities |
| 16 | ☐ | REP-DIRECT-001 | Repeats | BE | minLen=0, minLen=1, empty seq, all unique chars, all same char |
| 17 | ☐ | REP-PALIN-001 | Repeats | BE | minLen=0, minLen odd, maxLen > seqLen, empty seq |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | MC, BE | Invalid PAM sequences, non-DNA characters, empty seq, PAM longer than seq |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | MC, BE | No PAM sites in seq, guide length > seq, guide length 0, non-DNA |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | BE, MC | Zero mismatch tolerance, empty guide, guide of all N's |
| 21 | ☐ | PRIMER-TM-001 | MolTools | BE, INJ | Empty seq, 1-bp, 100+ bp, all-N, non-DNA chars, null |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | BE, MC | Seq shorter than min primer, GC range 0-0, Tm range inverted |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | BE | Extremely short, extremely long, palindromic primer, all-G |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | BE, MC | Seq shorter than min probe, Tm range inverted, empty seq |
| 25 | ☐ | PROBE-VALID-001 | MolTools | MC, INJ | Non-DNA probe, extremely short, null, probe with N's |
| 26 | ☐ | RESTR-FIND-001 | MolTools | MC, BE | Empty enzyme list, unknown enzyme names, empty sequence, non-DNA |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | BE | No cut sites, circular with 0 sites, 100+ enzymes, empty seq |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | BE, MC | No ATG in seq, all start no stop, minLen=0, non-DNA chars |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | BE, MC | No RBS in seq, overlapping genes, empty seq, non-coding only |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | BE | No -10 box, threshold=0, empty seq, seq shorter than motif |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | TF, MC, INJ | Malformed GFF string, missing fields, invalid strand, tab injection |
| 32 | ☐ | KMER-COUNT-001 | K-mer | BE | k=0, k > seqLen, k=seqLen, empty seq, k=1 |
| 33 | ☐ | KMER-FREQ-001 | K-mer | BE | k=0, k > seqLen, empty seq, single char seq |
| 34 | ☐ | KMER-FIND-001 | K-mer | BE | k=0, minFreq=0, minFreq > possible, empty seq |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | BE, MC | Empty vs empty, empty vs seq, very long seqs, non-DNA, single char |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | BE | Empty seqs, identical seqs, completely different seqs, 1-char seqs |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | BE | Empty seqs, no overlap, complete overlap, 1-char seqs |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | BE | Single seq, 2 seqs, 100+ seqs, empty list, seqs of length 1 |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | BE, MC | Identical seqs, empty seqs, single seq, non-DNA chars |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | BE | 0 seqs, 1 seq, 2 seqs, 100+ seqs, all identical seqs |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | TF, MC, INJ | Malformed Newick, unbalanced parens, missing semicolon, empty string |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | BE | Same tree vs same tree, different leaf sets, empty tree |
| 43 | ☐ | POP-FREQ-001 | PopGen | BE | 0 samples, 1 allele, all same allele, negative counts |
| 44 | ☐ | POP-DIV-001 | PopGen | BE | Single sequence, all identical, 0 seqs, very long seqs |
| 45 | ☐ | POP-HW-001 | PopGen | BE | 0 genotypes, all heterozygous, all homozygous, single genotype |
| 46 | ☐ | POP-FST-001 | PopGen | BE | Single population, identical populations, 0 samples |
| 47 | ☐ | POP-LD-001 | PopGen | BE | Single locus, monomorphic loci, 0 samples |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | BE, MC | No TTAGGG, only TTAGGG, empty seq, non-DNA chars |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | BE | Empty seq, all AT, all GC, extremely short |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | BE | 0 chromosomes, 100+ chromosomes, empty data |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | BE | 0 depth, negative depth, extremely high depth, empty data |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | BE | 0 ortholog pairs, self-comparison, empty genomes |
| 53 | ☐ | META-CLASS-001 | Metagenomics | MC, BE | Empty database, read with no match, non-DNA read, extremely short read |
| 54 | ☐ | META-PROF-001 | Metagenomics | BE | 0 reads, single read, all same taxon |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | BE | 0 species, 1 species, all equal abundance, single sample with 0 |
| 56 | ☐ | META-BETA-001 | Metagenomics | BE | Identical samples, empty samples, single species samples |
| 57 | ☐ | META-BIN-001 | Metagenomics | BE | 0 contigs, 1 contig, huge number of contigs, contigs with 0 length |
| 58 | ☐ | CODON-OPT-001 | Codon | MC, BE | Non-coding seq (no codons), seq len not divisible by 3, non-DNA chars |
| 59 | ☐ | CODON-CAI-001 | Codon | BE | Empty seq, seq with only stop codons, seq not divisible by 3 |
| 60 | ☐ | CODON-RARE-001 | Codon | BE | Empty seq, all rare, none rare, threshold=0, threshold=1 |
| 61 | ☐ | CODON-USAGE-001 | Codon | BE | Empty seq, seq len=1, seq len=2, very long seq |
| 62 | ☐ | TRANS-CODON-001 | Translation | MC | Invalid table ID, null code, non-standard table |
| 63 | ☐ | TRANS-PROT-001 | Translation | MC, BE | Non-DNA chars, empty, len=1, len=2, all stop codons, no stop codon |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | RB, TF, MC, INJ | Missing >, empty seq, binary garbage, huge headers, null bytes, no newline |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | RB, TF, MC, INJ | Missing @/+, quality len mismatch, invalid quality chars, truncated |
| 66 | ☐ | PARSE-BED-001 | FileIO | TF, MC, BE, INJ | Negative coords, start > end, non-numeric fields, tab injection |
| 67 | ☐ | PARSE-VCF-001 | FileIO | TF, MC, INJ | Missing header, invalid genotypes, huge alleles, missing ##fileformat |
| 68 | ☐ | PARSE-GFF-001 | FileIO | TF, MC, BE, INJ | Invalid strand, negative pos, malformed attributes, url encoding |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | RB, TF, MC | Missing LOCUS, truncated features, invalid locations, no // terminator |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | RB, TF, MC | Missing ID line, truncated seq, invalid feature keys, no // terminator |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | BE, MC | Empty, 1 base, all AU, non-RNA chars, extremely long |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | BE | Empty, no complement, minLoop=0, minStem > seqLen/2 |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | BE | Empty, no stacks, temperature=0, temperature negative |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | BE, MC | miRNA shorter than 8nt, empty, non-RNA, miRNA = target |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | BE | No complementarity, empty 3'UTR, miRNA longer than target |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | BE, MC | Too short, no hairpin, all one base, non-RNA chars |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | BE, MC | No GT sites, seq < window, non-DNA, all GT |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | BE, MC | No AG sites, seq < window, non-DNA, all AG |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | BE | No donor/acceptor, single exon, all intronic, empty seq |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | BE, MC | Empty protein, 1 AA, all same AA, non-standard AAs, X amino acid |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | BE | threshold=0, threshold=1, all ordered, all disordered, empty |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | MC, BE | Empty protein, non-amino acid chars, extremely short, all same char |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | MC, INJ | Invalid PROSITE pattern syntax, regex injection, empty pattern |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | BE | Protein shorter than min domain, empty, all X residues |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | BE | No CG dinucleotides, all CG, empty seq, windowSize=0, windowSize > seqLen |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | BE, MC | 0 expression, all NaN, negative expression, empty gene set, unknown genes |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 0 |
| ☐ Not started | 86 |
| High-priority (parsers + validation) | 12 |
| Medium-priority (boundary inputs) | 45 |
| Lower-priority (algorithm-specific edge cases) | 29 |
