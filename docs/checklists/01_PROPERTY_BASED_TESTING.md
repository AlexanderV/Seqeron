# Checklist 01: Property-Based Testing (FsCheck)

**Priority:** P0
**Framework:** FsCheck + FsCheck.NUnit
**Date:** 2026-03-19
**Total algorithms:** 86

---

## Description

Property-based testing генерирует сотни случайных входов и проверяет, что алгоритм удовлетворяет математическим инвариантам (range bounds, symmetry, idempotence, monotonicity). FsCheck интегрирован в проект через `FsCheck.NUnit`. Каждый геномный алгоритм имеет как минимум один выражаемый инвариант.

**Текущее покрытие:** 22 файла в `Properties/` — GcContent, GcSkew, Sequence, EditDistance, Hamming, FASTA, FileIO, Alignment, Codon, CRISPR, K-mer, MiRNA, PatternMatching, Phylogenetic, PopGen, PrimerProbe, ProteinMotif, RepeatFinder, Restriction, RnaStructure, SequenceComposition, Splicing.

**Типы инвариантов:**
- **R** = Range (результат в допустимом диапазоне)
- **S** = Symmetry (f(a,b) = f(b,a))
- **I** = Idempotence / Involution (f(f(x)) = x)
- **M** = Monotonicity (больше X → больше/меньше Y)
- **P** = Preservation (свойство сохраняется при трансформации)
- **RT** = Round-trip (parse(serialize(x)) = x)
- **D** = Determinism (одинаковый вход → одинаковый выход)

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Invariants to Verify | Property File |
|---|--------|-----------|------|---------------------|---------------|
| 1 | ☑ | SEQ-GC-001 | Composition | R: GC% ∈ [0,100]; P: complement preserves GC%; R: GcFraction ∈ [0,1]; P: GC%=Frac×100 | GcContentProperties.cs |
| 2 | ☑ | SEQ-COMP-001 | Composition | I: complement(complement(x))=x; P: length preserved | SequenceProperties.cs |
| 3 | ☑ | SEQ-REVCOMP-001 | Composition | I: revcomp(revcomp(x))=x; P: length preserved | SequenceProperties.cs |
| 4 | ☑ | SEQ-VALID-001 | Composition | R: result ∈ {true,false}; P: case-insensitive; D: deterministic; P: valid DNA ⊂ valid IUPAC | SequenceCompositionProperties.cs |
| 5 | ☑ | SEQ-COMPLEX-001 | Composition | R: complexity ∈ [0,1]; M: homopolymer → min; M: random long → high; P: permutation invariance | SequenceCompositionProperties.cs |
| 6 | ☑ | SEQ-ENTROPY-001 | Composition | R: entropy ≥ 0; P: permutation invariance; M: uniform dist → max entropy; D: deterministic | SequenceCompositionProperties.cs |
| 7 | ☑ | SEQ-GCSKEW-001 | Composition | R: skew ∈ [-1,1]; S: complement reverses sign; P: cumulative array length = seq length | GcSkewProperties.cs |
| 8 | ☑ | PAT-EXACT-001 | Matching | R: positions ∈ [0, len-patLen]; M: substring → count≥1; D: deterministic; P: total matches ≤ len-patLen+1 | PatternMatchingProperties.cs |
| 9 | ☑ | PAT-APPROX-001 | Matching | R: distance ≥ 0; S: hamming(a,b)=hamming(b,a); I: hamming(x,x)=0; R: distance ≤ len | HammingDistanceProperties.cs |
| 10 | ☑ | PAT-APPROX-002 | Matching | R: edit dist ≥ 0; S: editDist(a,b)=editDist(b,a); I: editDist(x,x)=0; triangle inequality | EditDistanceProperties.cs |
| 11 | ☑ | PAT-IUPAC-001 | Matching | M: more degenerate code → ≥ matches; P: N matches all 4 bases; D: deterministic | PatternMatchingProperties.cs |
| 12 | ☑ | PAT-PWM-001 | Matching | R: scores ∈ ℝ; M: lower threshold → ≥ matches; D: deterministic; P: consensus from PWM valid | PatternMatchingProperties.cs |
| 13 | ☑ | REP-STR-001 | Repeats | R: positions ≥ 0; M: lower minRepeats → ≥ results; R: repeat count ≥ minRepeats; P: unit len in range | RepeatFinderProperties.cs |
| 14 | ☑ | REP-TANDEM-001 | Repeats | R: repeat count ≥ minReps; M: wider unit range → ≥ results; R: positions valid; D: deterministic | RepeatFinderProperties.cs |
| 15 | ☑ | REP-INV-001 | Repeats | P: right arm = revcomp(left arm); R: positions valid; R: arm len ≥ minLen; D: deterministic | RepeatFinderProperties.cs |
| 16 | ☑ | REP-DIRECT-001 | Repeats | R: positions valid; M: lower minLen → ≥ results; P: two copies identical; D: deterministic | RepeatFinderProperties.cs |
| 17 | ☑ | REP-PALIN-001 | Repeats | P: palindrome = revcomp of self; R: len ∈ [minLen, maxLen]; R: positions valid | RepeatFinderProperties.cs |
| 18 | ☑ | CRISPR-PAM-001 | MolTools | R: positions valid; P: PAM motif at each site; M: longer seq → ≥ sites; D: deterministic | CrisprProperties.cs |
| 19 | ☑ | CRISPR-GUIDE-001 | MolTools | R: guide length = specified; P: target strand correct; R: score ∈ valid range | CrisprProperties.cs |
| 20 | ☑ | CRISPR-OFF-001 | MolTools | R: off-target score ∈ [0,1]; M: more mismatches → lower score; D: deterministic | CrisprProperties.cs |
| 21 | ☐ | PRIMER-TM-001 | MolTools | R: Tm > 0; M: longer GC-rich → higher Tm; D: deterministic; P: Tm in biologically valid range | PrimerProbeProperties.cs |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | R: primer len ∈ [min,max]; P: GC% in range; R: Tm in range; D: deterministic | PrimerProbeProperties.cs |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | R: hairpin ΔG ≤ 0; R: dimer score ≥ 0; D: deterministic | PrimerProbeProperties.cs |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | R: probe len ∈ [min,max]; M: stricter Tm → fewer probes; P: GC% in range | PrimerProbeProperties.cs |
| 25 | ☐ | PROBE-VALID-001 | MolTools | R: validation pass/fail; D: deterministic; P: Tm in expected range | PrimerProbeProperties.cs |
| 26 | ☑ | RESTR-FIND-001 | MolTools | R: positions valid; P: site sequence matches enzyme recognition; D: deterministic | RestrictionProperties.cs |
| 27 | ☑ | RESTR-DIGEST-001 | MolTools | P: sum(fragment lengths) = seq length; R: fragments ≥ 1; D: deterministic | RestrictionProperties.cs |
| 28 | ☑ | ANNOT-ORF-001 | Annotation | R: ORF start < end ≤ seqLen; P: starts with ATG; M: longer seq → ≥ ORFs; R: len divisible by 3 | AnnotationProperties.cs |
| 29 | ☑ | ANNOT-GENE-001 | Annotation | R: gene start < end; P: contains RBS motif upstream; D: deterministic | AnnotationProperties.cs |
| 30 | ☑ | ANNOT-PROM-001 | Annotation | R: position ≥ 0; P: contains -10/-35 box; M: lower score threshold → ≥ promoters | AnnotationProperties.cs |
| 31 | ☑ | ANNOT-GFF-001 | Annotation | RT: parse(serialize(features))=features; R: well-formed GFF3; P: coordinates 1-based | AnnotationProperties.cs |
| 32 | ☑ | KMER-COUNT-001 | K-mer | R: count > 0; P: sum(counts) = seqLen - k + 1; M: larger k → ≤ distinct k-mers | KmerProperties.cs |
| 33 | ☑ | KMER-FREQ-001 | K-mer | R: freq ∈ [0,1]; P: sum(freqs) = 1.0; D: deterministic | KmerProperties.cs |
| 34 | ☑ | KMER-FIND-001 | K-mer | R: positions valid; M: lower minFreq → ≥ k-mers returned; D: deterministic | KmerProperties.cs |
| 35 | ☑ | ALIGN-GLOBAL-001 | Alignment | S: score(a,b)=score(b,a); R: aligned len ≥ max(len1,len2); P: identity → max score; P: aligned1.len = aligned2.len | AlignmentProperties.cs |
| 36 | ☑ | ALIGN-LOCAL-001 | Alignment | R: score ≥ 0; P: aligned1.len = aligned2.len; M: identical substring → score ≥ matchScore×len | AlignmentProperties.cs |
| 37 | ☑ | ALIGN-SEMI-001 | Alignment | P: aligned1.len = aligned2.len; R: score ≥ 0; D: deterministic | AlignmentProperties.cs |
| 38 | ☑ | ALIGN-MULTI-001 | Alignment | P: all aligned sequences same length; R: score ≥ 0; D: deterministic | AlignmentProperties.cs |
| 39 | ☑ | PHYLO-DIST-001 | Phylogenetic | S: d(a,b)=d(b,a); I: d(x,x)=0; R: d ≥ 0; triangle inequality | PhylogeneticProperties.cs |
| 40 | ☑ | PHYLO-TREE-001 | Phylogenetic | R: N leaves = N input sequences; P: tree connected; R: all branch lengths ≥ 0 | PhylogeneticProperties.cs |
| 41 | ☑ | PHYLO-NEWICK-001 | Phylogenetic | RT: parse(serialize(tree))=tree; P: leaf labels preserved; D: deterministic | PhylogeneticProperties.cs |
| 42 | ☑ | PHYLO-COMP-001 | Phylogenetic | S: RF(a,b)=RF(b,a); I: RF(t,t)=0; R: RF ≥ 0; R: RF ≤ 2(n-3) | PhylogeneticProperties.cs |
| 43 | ☑ | POP-FREQ-001 | PopGen | R: freq ∈ [0,1]; P: sum(allele freqs) = 1.0; D: deterministic | PopulationGeneticsProperties.cs |
| 44 | ☑ | POP-DIV-001 | PopGen | R: π ≥ 0; R: θ ≥ 0; M: more diverse → higher π; D: deterministic | PopulationGeneticsProperties.cs |
| 45 | ☑ | POP-HW-001 | PopGen | R: p² + 2pq + q² = 1.0; R: chi² ≥ 0; D: deterministic | PopulationGeneticsProperties.cs |
| 46 | ☑ | POP-FST-001 | PopGen | R: Fst ∈ [0,1]; S: Fst(A,B)=Fst(B,A); I: Fst(A,A)=0; D: deterministic | PopulationGeneticsProperties.cs |
| 47 | ☑ | POP-LD-001 | PopGen | R: D' ∈ [-1,1]; R: r² ∈ [0,1]; S: LD(a,b)=LD(b,a); D: deterministic | PopulationGeneticsProperties.cs |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | R: positions valid; P: telomere seq contains TTAGGG; M: more repeats → longer region; D: deterministic | ChromosomeProperties.cs (new) |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | R: centromere index ∈ [0, seqLen]; D: deterministic; P: AT-rich region | ChromosomeProperties.cs (new) |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | R: chromosome count > 0; D: deterministic; P: each chrom has classification | ChromosomeProperties.cs (new) |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | R: copy number ≥ 0; M: higher depth → higher CN; D: deterministic | ChromosomeProperties.cs (new) |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | R: block positions valid; S: synteny(A,B) symmetric with B,A; D: deterministic | ChromosomeProperties.cs (new) |
| 53 | ☐ | META-CLASS-001 | Metagenomics | R: confidence ∈ [0,1]; P: assigned taxon in database; D: deterministic | MetagenomicsProperties.cs (new) |
| 54 | ☐ | META-PROF-001 | Metagenomics | R: abundances sum to 1.0; R: abundance ∈ [0,1]; D: deterministic | MetagenomicsProperties.cs (new) |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | R: Shannon ≥ 0; R: Simpson ∈ [0,1]; M: more species → higher diversity; D: deterministic | MetagenomicsProperties.cs |
| 56 | ☑ | META-BETA-001 | Metagenomics | S: dist(a,b)=dist(b,a); R: dist ∈ [0,1]; I: dist(x,x)=0; D: deterministic | MetagenomicsProperties.cs |
| 57 | ☐ | META-BIN-001 | Metagenomics | R: each contig in ≤ 1 bin; P: bin GC% consistent within bin; D: deterministic | MetagenomicsProperties.cs (new) |
| 58 | ☑ | CODON-OPT-001 | Codon | P: optimized translates to same protein; R: only valid codons; D: deterministic | CodonProperties.cs |
| 59 | ☑ | CODON-CAI-001 | Codon | R: CAI ∈ [0,1]; M: all optimal codons → CAI close to 1.0; D: deterministic | CodonProperties.cs |
| 60 | ☑ | CODON-RARE-001 | Codon | R: rare codon positions valid; M: lower threshold → more rare codons; D: deterministic | CodonProperties.cs |
| 61 | ☑ | CODON-USAGE-001 | Codon | R: usage freqs ≥ 0; P: sum per amino acid = 1.0; D: deterministic | CodonProperties.cs |
| 62 | ☑ | TRANS-CODON-001 | Translation | R: 64 codons mapped; P: start codons → M; P: stop codons → *; D: deterministic | CodonProperties.cs |
| 63 | ☑ | TRANS-PROT-001 | Translation | R: protein len ≤ seqLen/3; P: starts with M if starts with ATG; D: deterministic | CodonProperties.cs |
| 64 | ☑ | PARSE-FASTA-001 | FileIO | RT: write(parse(fasta))=fasta; P: header preserved; P: sequence preserved; D: deterministic | FastaRoundTripProperties.cs |
| 65 | ☑ | PARSE-FASTQ-001 | FileIO | RT: round-trip; P: quality len = seq len; R: quality scores in valid range | FileIOProperties.cs |
| 66 | ☑ | PARSE-BED-001 | FileIO | R: start < end; R: chrom non-empty; R: start ≥ 0; D: deterministic | FileIOProperties.cs |
| 67 | ☑ | PARSE-VCF-001 | FileIO | R: pos > 0; P: ref allele non-empty; R: qual ≥ 0 or missing; D: deterministic | FileIOProperties.cs |
| 68 | ☑ | PARSE-GFF-001 | FileIO | RT: round-trip; R: start ≤ end; R: strand ∈ {+,-,.}; D: deterministic | FileIOProperties.cs |
| 69 | ☑ | PARSE-GENBANK-001 | FileIO | RT: round-trip; P: locus line present; P: sequence preserved; D: deterministic | FileIOProperties.cs |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | RT: round-trip; P: ID line present; P: sequence preserved; D: deterministic | FileIOProperties.cs |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | R: pairs count ≤ len/2; P: no crossing pairs (Nussinov); P: paired bases complementary; D: deterministic | RnaStructureProperties.cs |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | R: stem len > 0; P: loop len ≥ minLoop; P: stem arms complementary; D: deterministic | RnaStructureProperties.cs |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | R: ΔG ≤ 0 for stable structures; M: more GC pairs → lower energy; D: deterministic | RnaStructureProperties.cs |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | R: seed len = 6–8; P: seed at 5' end (pos 2–8); D: deterministic | MiRnaProperties.cs |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | R: score ∈ [0,1]; M: perfect seed match → higher score; D: deterministic | MiRnaProperties.cs |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | R: precursor len > mature len; P: hairpin structure present; D: deterministic | MiRnaProperties.cs |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | R: score ∈ [0,1]; P: canonical GT at donor site; D: deterministic | SplicingProperties.cs |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | R: score ∈ [0,1]; P: canonical AG at acceptor site; D: deterministic | SplicingProperties.cs |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | R: exon start < end; P: introns flanked by GT…AG; D: deterministic | SplicingProperties.cs |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | R: score ∈ [0,1]; P: len(scores) = len(sequence); D: deterministic | DisorderProperties.cs |
| 81 | ☑ | DISORDER-REGION-001 | ProteinPred | R: region start < end ≤ seqLen; M: lower threshold → larger regions; D: deterministic | DisorderProperties.cs (new) |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | R: positions valid; M: broader pattern → ≥ matches; D: deterministic | ProteinMotifProperties.cs |
| 83 | ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | R: match positions valid; P: match conforms to PROSITE pattern regex; D: deterministic | ProteinMotifProperties.cs |
| 84 | ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | R: domain start < end; D: deterministic; P: domain score above threshold | ProteinMotifProperties.cs |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | R: CpG ratio ≥ 0; R: GC% ∈ [0,1]; M: more CG → higher ratio; D: deterministic | EpigeneticsProperties.cs (new) |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | R: infiltration score ∈ [0,1]; P: sum(cell fractions) ≤ 1.0; D: deterministic | OncologyProperties.cs (new) |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 66 |
| ☐ Not started | 20 |
| New property files needed | 4 (Chromosome, Epigenetics, Oncology) |
| Existing property files to extend | 15 |
