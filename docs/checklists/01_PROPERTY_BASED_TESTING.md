# Checklist 01: Property-Based Testing (FsCheck)

**Priority:** P0
**Framework:** FsCheck + FsCheck.NUnit
**Date:** 2026-03-19
**Total algorithms:** 234

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
| 21 | ☑ | PRIMER-TM-001 | MolTools | R: Tm > 0; M: longer GC-rich → higher Tm; D: deterministic; P: Tm in biologically valid range | PrimerProbeProperties.cs |
| 22 | ☑ | PRIMER-DESIGN-001 | MolTools | R: primer len ∈ [min,max]; P: GC% in range; R: Tm in range; D: deterministic | PrimerProbeProperties.cs |
| 23 | ☑ | PRIMER-STRUCT-001 | MolTools | R: hairpin ΔG ≤ 0; R: dimer score ≥ 0; D: deterministic | PrimerProbeProperties.cs |
| 24 | ☑ | PROBE-DESIGN-001 | MolTools | R: probe len ∈ [min,max]; M: stricter Tm → fewer probes; P: GC% in range | PrimerProbeProperties.cs |
| 25 | ☑ | PROBE-VALID-001 | MolTools | R: validation pass/fail; D: deterministic; P: Tm in expected range | PrimerProbeProperties.cs |
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
| 48 | ☑ | CHROM-TELO-001 | Chromosome | R: positions valid; P: telomere seq contains TTAGGG; M: more repeats → longer region; D: deterministic | ChromosomeProperties.cs (new) |
| 49 | ☑ | CHROM-CENT-001 | Chromosome | R: centromere index ∈ [0, seqLen]; D: deterministic; P: AT-rich region | ChromosomeProperties.cs (new) |
| 50 | ☑ | CHROM-KARYO-001 | Chromosome | R: chromosome count > 0; D: deterministic; P: each chrom has classification | ChromosomeProperties.cs (new) |
| 51 | ☑ | CHROM-ANEU-001 | Chromosome | R: copy number ≥ 0; M: higher depth → higher CN; D: deterministic | ChromosomeProperties.cs (new) |
| 52 | ☑ | CHROM-SYNT-001 | Chromosome | R: block positions valid; S: synteny(A,B) symmetric with B,A; D: deterministic | ChromosomeProperties.cs (new) |
| 53 | ☑ | META-CLASS-001 | Metagenomics | R: confidence ∈ [0,1]; P: assigned taxon in database; D: deterministic | MetagenomicsProperties.cs (new) |
| 54 | ☑ | META-PROF-001 | Metagenomics | R: abundances sum to 1.0; R: abundance ∈ [0,1]; D: deterministic | MetagenomicsProperties.cs (new) |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | R: Shannon ≥ 0; R: Simpson ∈ [0,1]; M: more species → higher diversity; D: deterministic | MetagenomicsProperties.cs |
| 56 | ☑ | META-BETA-001 | Metagenomics | S: dist(a,b)=dist(b,a); R: dist ∈ [0,1]; I: dist(x,x)=0; D: deterministic | MetagenomicsProperties.cs |
| 57 | ☑ | META-BIN-001 | Metagenomics | R: each contig in ≤ 1 bin; P: bin GC% consistent within bin; D: deterministic | MetagenomicsProperties.cs (new) |
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
| 70 | ☑ | PARSE-EMBL-001 | FileIO | RT: round-trip; P: ID line present; P: sequence preserved; D: deterministic | FileIOProperties.cs |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | R: pairs count ≤ len/2; P: no crossing pairs (Nussinov); P: paired bases complementary; D: deterministic | RnaStructureProperties.cs |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | R: stem len > 0; P: loop len ≥ minLoop; P: stem arms complementary; D: deterministic | RnaStructureProperties.cs |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | R: ΔG ≤ 0 for stable structures; M: more GC pairs → lower energy; D: deterministic | RnaStructureProperties.cs |
| 74 | ☑ | MIRNA-SEED-001 | MiRNA | R: seed len = 6–8; P: seed at 5' end (pos 2–8); D: deterministic | MiRnaProperties.cs |
| 75 | ☑ | MIRNA-TARGET-001 | MiRNA | R: score ∈ [0,1]; M: perfect seed match → higher score; D: deterministic | MiRnaProperties.cs |
| 76 | ☑ | MIRNA-PRECURSOR-001 | MiRNA | R: precursor len > mature len; P: hairpin structure present; D: deterministic | MiRnaProperties.cs |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | R: score ∈ [0,1]; P: canonical GT at donor site; D: deterministic | SplicingProperties.cs |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | R: score ∈ [0,1]; P: canonical AG at acceptor site; D: deterministic | SplicingProperties.cs |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | R: exon start < end; P: introns flanked by GT…AG; D: deterministic | SplicingProperties.cs |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | R: score ∈ [0,1]; P: len(scores) = len(sequence); D: deterministic | DisorderProperties.cs |
| 81 | ☑ | DISORDER-REGION-001 | ProteinPred | R: region start < end ≤ seqLen; M: lower threshold → larger regions; D: deterministic | DisorderProperties.cs (new) |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | R: positions valid; M: broader pattern → ≥ matches; D: deterministic | ProteinMotifProperties.cs |
| 83 | ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | R: match positions valid; P: match conforms to PROSITE pattern regex; D: deterministic | ProteinMotifProperties.cs |
| 84 | ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | R: domain start < end; D: deterministic; P: domain score above threshold | ProteinMotifProperties.cs |
| 85 | ☑ | EPIGEN-CPG-001 | Epigenetics | R: CpG ratio ≥ 0; R: GC% ∈ [0,1]; M: more CG → higher ratio; D: deterministic | EpigeneticsProperties.cs (new) |
| 86 | ☑ | ONCO-IMMUNE-001 | Oncology | R: infiltration score ∈ [0,1]; P: sum(cell fractions) ≤ 1.0; D: deterministic | OncologyProperties.cs (new) |
| 87 | ☑ | ONCO-SOMATIC-001 | Oncology | R: VAF ∈ [0,1]; P: somatic calls absent in matched normal; M: higher tumor depth → ≥ confident calls; D: deterministic | OncologyProperties.cs (new) |
| 88 | ☑ | ONCO-VAF-001 | Oncology | R: VAF = alt/(alt+ref) ∈ [0,1]; M: more alt reads → higher VAF; R: VAF=0 when no alt reads; D: deterministic | OncologyProperties.cs (new) |
| 89 | ☑ | ONCO-DRIVER-001 | Oncology | R: driver score ≥ 0; M: more recurrent/hotspot → higher driver likelihood; P: known oncogenes flagged; D: deterministic | OncologyProperties.cs |
| 90 | ☑ | ONCO-ARTIFACT-001 | Oncology | P: passing set ⊆ input; M: stricter thresholds → ≤ survivors; R: strand-bias ∈ [0,1]; D: deterministic | OncologyProperties.cs (new) |
| 91 | ☑ | ONCO-ANNOT-001 | Oncology | R: consequence ∈ SO vocabulary; P: annotation preserves variant coordinates; D: deterministic | OncologyProperties.cs (new) |
| 92 | ☑ | ONCO-TMB-001 | Oncology | R: TMB ≥ 0 mut/Mb; M: more nonsynonymous mutations → higher TMB; P: TMB = count / panel-Mb; D: deterministic | OncologyProperties.cs (new) |
| 93 | ☑ | ONCO-MSI-001 | Oncology | R: MSI fraction ∈ [0,1]; M: more unstable loci → higher score; P: MSI-H ⟺ score ≥ threshold; D: deterministic | OncologyProperties.cs (new) |
| 94 | ☑ | ONCO-HRD-001 | Oncology | R: HRD = LOH+TAI+LST ≥ 0; M: more genomic scars → higher HRD; P: additive of three components; D: deterministic | OncologyProperties.cs (new) |
| 95 | ☑ | ONCO-LOH-001 | Oncology | R: region positions valid; P: LOH ⟺ BAF → 0/1 (one allele lost); M: lower BAF-dev threshold → ≥ LOH; D: deterministic | OncologyProperties.cs (new) |
| 96 | ☑ | ONCO-SIG-001 | Oncology | R: exactly 96 SBS channels; P: each SNV → one pyrimidine-centred trinucleotide; P: channel counts sum = #SNVs; D: deterministic | OncologyProperties.cs (new) |
| 97 | ☑ | ONCO-SIG-002 | Oncology | R: exposures ≥ 0; P: Σ exposures = #mutations (or normalised 1.0); M: better fit → lower reconstruction error; D: deterministic | OncologyProperties.cs (new) |
| 98 | ☑ | ONCO-SIG-003 | Oncology | R: CI_lower ≤ point ≤ CI_upper; P: bootstrap mean ≈ point estimate; D: deterministic given seed | OncologyProperties.cs (new) |
| 99 | ☑ | ONCO-SIG-004 | Oncology | P: dominant process = argmax exposure; R: confidence ∈ [0,1]; D: deterministic | OncologyProperties.cs (new) |
| 100 | ☑ | ONCO-FUSION-001 | Oncology | R: breakpoint positions valid; P: fusion joins two distinct genes; M: more split/spanning reads → higher confidence; D: deterministic | OncologyProperties.cs (new) |
| 101 | ☑ | ONCO-FUSION-002 | Oncology | P: matched fusions ⊆ known DB; R: match ∈ {true,false}; D: deterministic | OncologyProperties.cs (new) |
| 102 | ☑ | ONCO-FUSION-003 | Oncology | P: reading frame (in/out) correctly derived; R: breakpoint within gene bounds; D: deterministic | OncologyProperties.cs (new) |
| 103 | ☑ | ONCO-CNA-001 | Oncology | R: copy number ≥ 0; M: higher log2-ratio → higher CN; P: CN=2 → neutral; D: deterministic | OncologyProperties.cs (new) |
| 104 | ☑ | ONCO-CNA-002 | Oncology | P: focal segment length < arm-level cutoff; M: higher CN → amplified; R: positions valid; D: deterministic | OncologyProperties.cs (new) |
| 105 | ☑ | ONCO-CNA-003 | Oncology | P: CN ≈ 0 over deletion; R: positions valid; M: higher CN threshold → ≤ deletions; D: deterministic | OncologyProperties.cs (new) |
| 106 | ☑ | ONCO-PURITY-001 | Oncology | R: purity ∈ [0,1]; M: higher clonal VAF → higher purity; D: deterministic | OncologyProperties.cs (new) |
| 107 | ☑ | ONCO-PLOIDY-001 | Oncology | R: ploidy > 0; M: more amplified genome → higher ploidy; D: deterministic | OncologyProperties.cs (new) |
| 108 | ☑ | ONCO-CLONAL-001 | Oncology | P: clonal ⟺ CCF ≈ 1, else subclonal; R: class ∈ enum; D: deterministic | OncologyProperties.cs (new) |
| 109 | ☑ | ONCO-NEO-001 | Oncology | R: peptide length ∈ [8,11]; P: mutated residue inside every peptide window; P: peptides tile the mutation; D: deterministic | OncologyProperties.cs (new) |
| 110 | ☑ | ONCO-MHC-001 | Oncology | R: %rank ∈ [0,100]; M: lower IC50 → stronger binding; P: strong binder ⟺ rank ≤ threshold; D: deterministic | OncologyProperties.cs (new) |
| 111 | ☑ | ONCO-CTDNA-001 | Oncology | R: ctDNA fraction ∈ [0,1]; M: more tumor-supporting reads → higher fraction; D: deterministic | OncologyProperties.cs (new) |
| 112 | ☑ | ONCO-MRD-001 | Oncology | R: detection ∈ {positive,negative}; M: more tracked variants observed → MRD-positive; R: sensitivity ∈ [0,1]; D: deterministic | OncologyProperties.cs (new) |
| 113 | ☑ | ONCO-CHIP-001 | Oncology | P: CHIP flagged by gene list ∧ VAF band; P: survivors ⊆ input; D: deterministic | OncologyProperties.cs (new) |
| 114 | ☑ | ONCO-PHYLO-001 | Oncology | R: #leaves = #clones/samples; P: trunk mutations shared by all clones; R: branch lengths ≥ 0; D: deterministic | OncologyProperties.cs (new) |
| 115 | ☑ | ONCO-CCF-001 | Oncology | R: CCF ∈ [0,1]; P: CCF derived from VAF, CN, purity; M: higher VAF → higher CCF; D: deterministic | OncologyProperties.cs (new) |
| 116 | ☑ | ONCO-HETERO-001 | Oncology | R: heterogeneity (MATH) ≥ 0; M: wider VAF spread → higher heterogeneity; D: deterministic | OncologyProperties.cs (new) |
| 117 | ☑ | ONCO-HLA-001 | Oncology | P: alleles valid HLA nomenclature; R: ≤ 2 alleles per locus; D: deterministic | OncologyProperties.cs (new) |
| 118 | ☑ | ONCO-ACTION-001 | Oncology | R: evidence tier ∈ ordered levels (A>B>C>D); P: known actionable variant → tier assigned; D: deterministic | OncologyProperties.cs (new) |
| 119 | ☑ | ONCO-SV-001 | Oncology | P: chromothripsis/chromoplexy classified by breakpoint pattern; R: breakpoint count ≥ 0; D: deterministic | OncologyProperties.cs (new) |
| 120 | ☑ | ONCO-EXPR-001 | Oncology | P: outlier ⟺ |z-score| > threshold; M: lower threshold → ≥ outliers; R: z-scores finite; D: deterministic | OncologyProperties.cs (new) |
| 121 | ☑ | SEQ-COMPOSITION-001 | Statistics | R: each fraction ∈ [0,1]; P: Σ fractions = 1.0; P: counts sum = length; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 122 | ☑ | SEQ-DINUC-001 | Statistics | R: each frequency ≥ 0; P: Σ dinucleotide counts = length−1; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 123 | ☑ | SEQ-HYDRO-001 | Statistics | R: score finite within scale range (Kyte-Doolittle ∈ [−4.5,4.5]); M: more hydrophobic residues → higher mean; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 124 | ☑ | SEQ-MW-001 | Statistics | R: MW > 0 for non-empty; M: longer sequence → higher MW; P: MW additive over residues; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 125 | ☑ | SEQ-PI-001 | Statistics | R: pI ∈ [0,14]; P: net charge at pI ≈ 0; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 126 | ☑ | SEQ-SECSTRUCT-001 | Statistics | R: each propensity ≥ 0; P: every residue assigned H/E/C; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 127 | ☑ | SEQ-STATS-001 | Statistics | R: counts ≥ 0; P: Σ counts = length; R: GC% ∈ [0,100]; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 128 | ☑ | SEQ-SUMMARY-001 | Statistics | P: reported length = sequence length; R: GC% ∈ [0,100]; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 129 | ☑ | SEQ-THERMO-001 | Statistics | R: ΔG, ΔH, ΔS finite; M: more GC → more stable (lower ΔG); D: deterministic | SequenceStatisticsProperties.cs (new) |
| 130 | ☑ | SEQ-TM-001 | Statistics | R: Tm ≥ 0; M: more GC → higher Tm; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 131 | ☑ | COMPGEN-ANI-001 | Comparative | R: ANI ∈ [0,100]%; S: ANI(A,B)=ANI(B,A); I: ANI(A,A)=100; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 132 | ☑ | COMPGEN-CLUSTER-001 | Comparative | R: cluster size ≥ 1; P: cluster genes conserved across genomes; M: lower identity threshold → ≥ clusters; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 133 | ☑ | COMPGEN-COMPARE-001 | Comparative | S: compare(A,B) metrics symmetric; R: similarity ∈ [0,1]; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 134 | ☑ | COMPGEN-DOTPLOT-001 | Comparative | R: dot positions valid (i<lenA, j<lenB); P: main diagonal for identical seqs; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 135 | ☑ | COMPGEN-ORTHO-001 | Comparative | P: ortholog pairs bidirectional; R: positions valid; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 136 | ☑ | COMPGEN-RBH-001 | Comparative | S: RBH(A,B)=RBH(B,A); P: each gene in ≤ 1 RBH pair; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 137 | ☑ | COMPGEN-REARR-001 | Comparative | R: breakpoint count ≥ 0; P: rearrangement preserves gene content; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 138 | ☑ | COMPGEN-REVERSAL-001 | Comparative | R: distance ≥ 0; S: d(A,B)=d(B,A); I: d(A,A)=0; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 139 | ☑ | COMPGEN-SYNTENY-001 | Comparative | R: block positions valid; M: lower minBlockSize → ≥ blocks; P: blocks collinear; D: deterministic | ComparativeGenomicsProperties.cs (new) |
| 140 | ☑ | ASSEMBLY-CONSENSUS-001 | Assembly | R: consensus length ≥ 0; P: each position = majority base; D: deterministic | AssemblyProperties.cs (new) |
| 141 | ☑ | ASSEMBLY-CORRECT-001 | Assembly | P: corrected reads keep length; M: more coverage → fewer errors; D: deterministic | AssemblyProperties.cs (new) |
| 142 | ☑ | ASSEMBLY-COVER-001 | Assembly | R: coverage ≥ 0; P: mean coverage = total bases / ref length; D: deterministic | AssemblyProperties.cs (new) |
| 143 | ☑ | ASSEMBLY-DBG-001 | Assembly | R: contig positions valid; P: contigs are k-mer paths; M: larger k → ≤ branching; D: deterministic | AssemblyProperties.cs (new) |
| 144 | ☑ | ASSEMBLY-MERGE-001 | Assembly | P: merged length ≥ max input; P: overlap consistent; D: deterministic | AssemblyProperties.cs (new) |
| 145 | ☑ | ASSEMBLY-OLC-001 | Assembly | R: contig count ≥ 1; P: overlaps ≥ minOverlap; D: deterministic | AssemblyProperties.cs (new) |
| 146 | ☑ | ASSEMBLY-SCAFFOLD-001 | Assembly | P: scaffold order respects links; R: gap sizes ≥ 0; D: deterministic | AssemblyProperties.cs (new) |
| 147 | ☑ | ASSEMBLY-STATS-001 | Assembly | R: N50 ≥ 0; P: Σ contig lengths = total; M: longer contigs → higher N50; D: deterministic | AssemblyProperties.cs (new) |
| 148 | ☑ | ASSEMBLY-TRIM-001 | Assembly | P: trimmed length ≤ original; M: higher quality cutoff → shorter reads; D: deterministic | AssemblyProperties.cs (new) |
| 149 | ☑ | RNA-DOTBRACKET-001 | RnaStructure | RT: parse∘format = identity; P: balanced brackets → valid pairs; R: pair count ≤ len/2; D: deterministic | RnaStructureProperties.cs |
| 150 | ☑ | RNA-HAIRPIN-001 | RnaStructure | R: loop size ≥ minLoop; M: larger destabilising loop → higher energy; D: deterministic | RnaStructureProperties.cs |
| 151 | ☑ | RNA-INVERT-001 | RnaStructure | P: arms reverse-complementary; R: positions valid; D: deterministic | RnaStructureProperties.cs |
| 152 | ☑ | RNA-MFE-001 | RnaStructure | R: MFE ≤ 0; M: more GC pairs → lower energy; D: deterministic | RnaStructureProperties.cs |
| 153 | ☑ | RNA-PAIR-001 | RnaStructure | P: only A-U, G-C, G-U pair; S: canPair(a,b)=canPair(b,a); D: deterministic | RnaStructureProperties.cs |
| 154 | ☑ | RNA-PARTITION-001 | RnaStructure | R: Z > 0; R: base-pair probability ∈ [0,1]; D: deterministic | RnaStructureProperties.cs |
| 155 | ☑ | RNA-PSEUDOKNOT-001 | RnaStructure | P: detects crossing pairs; R: positions valid; D: deterministic | RnaStructureProperties.cs |
| 156 | ☑ | KMER-ASYNC-001 | K-mer | P: async result = sync result; D: deterministic | KmerProperties.cs |
| 157 | ☑ | KMER-BOTH-001 | K-mer | P: count = forward + reverse-complement; S: strand-symmetric; D: deterministic | KmerProperties.cs |
| 158 | ☑ | KMER-DIST-001 | K-mer | R: distance ≥ 0; S: d(a,b)=d(b,a); I: d(x,x)=0; D: deterministic | KmerProperties.cs |
| 159 | ☑ | KMER-GENERATE-001 | K-mer | R: count = 4^k (DNA); P: all distinct; D: deterministic | KmerProperties.cs |
| 160 | ☑ | KMER-POSITIONS-001 | K-mer | R: positions ∈ [0, len−k]; P: seq[pos..pos+k] = kmer; D: deterministic | KmerProperties.cs |
| 161 | ☑ | KMER-STATS-001 | K-mer | R: counts ≥ 0; P: Σ counts = len−k+1; D: deterministic | KmerProperties.cs |
| 162 | ☑ | KMER-UNIQUE-001 | K-mer | P: unique k-mers have count 1; R: minCount ≥ 1; D: deterministic | KmerProperties.cs |
| 163 | ☑ | PROTMOTIF-CC-001 | ProteinMotif | R: score ∈ [0,1]; P: heptad periodicity detected; D: deterministic | ProteinMotifProperties.cs |
| 164 | ☑ | PROTMOTIF-COMMON-001 | ProteinMotif | R: positions valid; M: more sequences sharing → ≥ support; D: deterministic | ProteinMotifProperties.cs |
| 165 | ☑ | PROTMOTIF-LC-001 | ProteinMotif | R: region start < end; M: lower complexity threshold → ≥ regions; D: deterministic | ProteinMotifProperties.cs |
| 166 | ☑ | PROTMOTIF-PATTERN-001 | ProteinMotif | P: match conforms to pattern; R: positions valid; D: deterministic | ProteinMotifProperties.cs |
| 167 | ☑ | PROTMOTIF-SP-001 | ProteinMotif | R: cleavage site ∈ [1, len]; P: N-terminal hydrophobic signal; D: deterministic | ProteinMotifProperties.cs |
| 168 | ☑ | PROTMOTIF-TM-001 | ProteinMotif | R: helix length in range; P: hydrophobic stretch; M: lower threshold → ≥ helices; D: deterministic | ProteinMotifProperties.cs |
| 169 | ☑ | MOTIF-CONS-001 | Matching | P: consensus length = alignment width; P: each column = majority residue; D: deterministic | PatternMatchingProperties.cs |
| 170 | ☑ | MOTIF-DISCOVER-001 | Matching | R: motif length = k; M: lower support → ≥ motifs; D: deterministic | PatternMatchingProperties.cs |
| 171 | ☑ | MOTIF-GENERATE-001 | Matching | P: consensus from counts; R: length = motif width; D: deterministic | PatternMatchingProperties.cs |
| 172 | ☑ | MOTIF-REGULATORY-001 | Matching | R: positions valid; P: match conforms to known element; D: deterministic | PatternMatchingProperties.cs |
| 173 | ☑ | MOTIF-SHARED-001 | Matching | P: shared motif present in all inputs; R: positions valid; D: deterministic | PatternMatchingProperties.cs |
| 174 | ☑ | PAT-APPROX-003 | Matching | R: best distance ≥ 0; P: exact match → distance 0; M: best ≤ any match distance; D: deterministic | PatternMatchingProperties.cs |
| 175 | ☑ | GENOMIC-COMMON-001 | Analysis | P: common region ⊆ all inputs; R: positions valid; D: deterministic | GenomicAnalyzerProperties.cs (new) |
| 176 | ☑ | GENOMIC-MOTIFS-001 | Analysis | R: positions valid; P: motif matches known set; D: deterministic | GenomicAnalyzerProperties.cs (new) |
| 177 | ☑ | GENOMIC-ORF-001 | Analysis | R: start < end; P: starts ATG, ends stop; R: len % 3 = 0; D: deterministic | GenomicAnalyzerProperties.cs (new) |
| 178 | ☑ | GENOMIC-REPEAT-001 | Analysis | R: positions valid; M: lower minLen → ≥ repeats; D: deterministic | GenomicAnalyzerProperties.cs (new) |
| 179 | ☑ | GENOMIC-SIMILARITY-001 | Analysis | R: similarity ∈ [0,1]; S: sim(a,b)=sim(b,a); I: sim(x,x)=1; D: deterministic | GenomicAnalyzerProperties.cs (new) |
| 180 | ☑ | GENOMIC-TANDEM-001 | Analysis | R: repeat count ≥ minReps; P: unit repeated contiguously; D: deterministic | GenomicAnalyzerProperties.cs (new) |
| 181 | ☑ | EPIGEN-AGE-001 | Epigenetics | R: age ≥ 0; M: more methylation at clock sites → higher age; D: deterministic | EpigeneticsProperties.cs (new) |
| 182 | ☑ | EPIGEN-BISULF-001 | Epigenetics | P: unmethylated C→T, methylated C preserved; P: length preserved; D: deterministic | EpigeneticsProperties.cs (new) |
| 183 | ☑ | EPIGEN-CHROM-001 | Epigenetics | P: each region assigned a state; R: positions valid; D: deterministic | EpigeneticsProperties.cs (new) |
| 184 | ☑ | EPIGEN-DMR-001 | Epigenetics | R: region start < end; M: lower threshold → ≥ DMRs; P: |Δmethylation| ≥ threshold; D: deterministic | EpigeneticsProperties.cs (new) |
| 185 | ☑ | EPIGEN-METHYL-001 | Epigenetics | R: methylation level ∈ [0,1]; P: = methylated/total; D: deterministic | EpigeneticsProperties.cs (new) |
| 186 | ☑ | VARIANT-ANNOT-001 | Variants | R: impact ∈ ordered categories; P: annotation preserves position; D: deterministic | VariantProperties.cs (new) |
| 187 | ☑ | VARIANT-CALL-001 | Variants | R: positions valid; P: called only where pileup differs from ref; M: higher depth → ≥ confidence; D: deterministic | VariantProperties.cs (new) |
| 188 | ☑ | VARIANT-INDEL-001 | Variants | R: indel length > 0; P: ref and alt lengths differ; R: positions valid; D: deterministic | VariantProperties.cs (new) |
| 189 | ☑ | VARIANT-SNP-001 | Variants | P: ref/alt single bases, ref≠alt; R: positions valid; D: deterministic | VariantProperties.cs (new) |
| 190 | ☑ | PANGEN-CLUSTER-001 | PanGenome | P: each gene in exactly one cluster; M: lower identity → ≥ merging; D: deterministic | PanGenomeProperties.cs (new) |
| 191 | ☑ | PANGEN-CORE-001 | PanGenome | P: core ⊆ every genome; P: core + accessory = pan; M: more genomes → ≤ core size; D: deterministic | PanGenomeProperties.cs (new) |
| 192 | ☑ | PANGEN-HEAP-001 | PanGenome | R: Heaps' α > 0; P: open vs closed by α threshold; D: deterministic | PanGenomeProperties.cs (new) |
| 193 | ☑ | PANGEN-MARKER-001 | PanGenome | P: markers ⊆ core genes; R: marker count ≤ requested; D: deterministic | PanGenomeProperties.cs (new) |
| 194 | ☑ | META-FUNC-001 | Metagenomics | R: function scores ≥ 0; P: assigned function in DB; D: deterministic | MetagenomicsProperties.cs |
| 195 | ☑ | META-PATHWAY-001 | Metagenomics | R: p-value ∈ [0,1]; M: more pathway genes → higher enrichment; D: deterministic | MetagenomicsProperties.cs |
| 196 | ☑ | META-RESIST-001 | Metagenomics | P: hit matches resistance DB; R: positions valid; D: deterministic | MetagenomicsProperties.cs |
| 197 | ☑ | META-TAXA-001 | Metagenomics | R: p-value ∈ [0,1]; P: significant ⟺ p ≤ α; D: deterministic | MetagenomicsProperties.cs |
| 198 | ☑ | TRANS-DIFF-001 | Transcriptome | R: |log2FC| ≥ 0; S: FC(A,B) = −FC(B,A); R: p-value ∈ [0,1]; D: deterministic | TranscriptomeProperties.cs (new) |
| 199 | ☑ | TRANS-EXPR-001 | Transcriptome | R: expression ≥ 0; P: Σ TPM = 1e6; D: deterministic | TranscriptomeProperties.cs (new) |
| 200 | ☑ | TRANS-SPLICE-001 | Transcriptome | R: PSI ∈ [0,1]; P: exon coordinates valid; D: deterministic | TranscriptomeProperties.cs (new) |
| 201 | ☑ | SV-BREAKPOINT-001 | StructuralVar | R: breakpoint positions valid; M: more split reads → ≥ confidence; D: deterministic | StructuralVariantProperties.cs (new) |
| 202 | ☑ | SV-CNV-001 | StructuralVar | R: copy number ≥ 0; M: higher coverage ratio → higher CN; D: deterministic | StructuralVariantProperties.cs (new) |
| 203 | ☑ | SV-DETECT-001 | StructuralVar | R: SV type ∈ enum; P: positions valid; D: deterministic | StructuralVariantProperties.cs (new) |
| 204 | ☑ | DISORDER-LC-001 | ProteinPred | R: region start < end; M: lower complexity threshold → ≥ regions; D: deterministic | DisorderProperties.cs |
| 205 | ☑ | DISORDER-MORF-001 | ProteinPred | P: MoRF lies within a disordered region; R: positions valid; D: deterministic | DisorderProperties.cs |
| 206 | ☑ | DISORDER-PROPENSITY-001 | ProteinPred | R: propensity ∈ [0,1]; P: len(scores) = len(sequence); D: deterministic | DisorderProperties.cs |
| 207 | ☑ | POP-ANCESTRY-001 | PopGen | R: each proportion ∈ [0,1]; P: Σ proportions = 1.0; D: deterministic | PopulationGeneticsProperties.cs |
| 208 | ☑ | POP-ROH-001 | PopGen | R: ROH start < end; M: lower minLen → ≥ ROH; P: homozygous within run; D: deterministic | PopulationGeneticsProperties.cs |
| 209 | ☑ | POP-SELECT-001 | PopGen | R: statistic finite; M: stronger differentiation → higher signal; D: deterministic | PopulationGeneticsProperties.cs |
| 210 | ☐ | SEQ-ATSKEW-001 | Composition | R: skew ∈ [−1,1]; S: complement reverses sign; D: deterministic | GcSkewProperties.cs |
| 211 | ☐ | SEQ-REPLICATION-001 | Composition | R: origin index ∈ [0, len]; P: at cumulative-skew extremum; D: deterministic | GcSkewProperties.cs |
| 212 | ☐ | SEQ-RNACOMP-001 | Composition | I: complement∘complement = identity; P: A↔U, G↔C; D: deterministic | SequenceProperties.cs |
| 213 | ☐ | CODON-ENC-001 | Codon | R: ENC ∈ [20,61]; M: more biased usage → lower ENC; D: deterministic | CodonProperties.cs |
| 214 | ☐ | CODON-RSCU-001 | Codon | R: RSCU ≥ 0; P: mean RSCU per amino acid = 1; D: deterministic | CodonProperties.cs |
| 215 | ☐ | CODON-STATS-001 | Codon | R: counts ≥ 0; P: Σ codon counts = len/3; D: deterministic | CodonProperties.cs |
| 216 | ☐ | ANNOT-CODING-001 | Annotation | R: coding score ∈ [0,1]; M: real ORF → higher score; D: deterministic | AnnotationProperties.cs |
| 217 | ☐ | ANNOT-CODONUSAGE-001 | Annotation | R: frequencies ≥ 0; P: sum per amino acid = 1.0; D: deterministic | AnnotationProperties.cs |
| 218 | ☐ | ANNOT-REPEAT-001 | Annotation | R: positions valid; M: lower minLen → ≥ elements; D: deterministic | AnnotationProperties.cs |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | R: Q ≥ 0; P: Q = ASCII − offset; RT: encode∘decode = identity; D: deterministic | FileIOProperties.cs |
| 220 | ☐ | QUALITY-STATS-001 | Quality | R: mean Q ≥ 0; P: len(scores) = len(sequence); D: deterministic | FileIOProperties.cs |
| 221 | ☐ | PHYLO-BOOT-001 | Phylogenetic | R: support ∈ [0,100]; M: more replicates → stable support; D: deterministic given seed | PhylogeneticProperties.cs |
| 222 | ☐ | PHYLO-STATS-001 | Phylogenetic | R: tree depth ≥ 0; P: leaf count consistent; D: deterministic | PhylogeneticProperties.cs |
| 223 | ☐ | TRANS-SIXFRAME-001 | Translation | R: exactly 6 frames; P: 3 forward + 3 reverse-complement; D: deterministic | CodonProperties.cs |
| 224 | ☐ | RESTR-FILTER-001 | MolTools | P: filtered ⊆ all sites; M: stricter criteria → ≤ sites; D: deterministic | RestrictionProperties.cs |
| 225 | ☐ | MIRNA-PAIR-001 | MiRNA | P: seed region paired; R: alignment score ≥ 0; D: deterministic | MiRnaProperties.cs |
| 226 | ☐ | ALIGN-STATS-001 | Alignment | R: identity ∈ [0,1]; P: matches+mismatches+gaps = alignment length; D: deterministic | AlignmentProperties.cs |
| 227 | ☐ | SEQ-CODON-FREQ-001 | Statistics | R: each freq ≥ 0; P: Σ codon counts = len/3; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 228 | ☐ | SEQ-COMPLEX-COMPRESS-001 | Complexity | R: ratio ∈ (0,1]; M: repetitive → lower ratio; D: deterministic | SequenceComplexityProperties.cs (new) |
| 229 | ☐ | SEQ-COMPLEX-DUST-001 | Complexity | R: DUST score ≥ 0; M: low-complexity → higher score; D: deterministic | SequenceComplexityProperties.cs (new) |
| 230 | ☐ | SEQ-COMPLEX-KMER-001 | Complexity | R: entropy ≥ 0; M: more distinct k-mers → higher entropy; P: homopolymer → 0; D: deterministic | SequenceComplexityProperties.cs (new) |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | R: each window score ∈ [0,1]; P: window count = len−w+1; D: deterministic | SequenceComplexityProperties.cs (new) |
| 232 | ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | R: each entropy ≥ 0; P: profile length = len−w+1; D: deterministic | SequenceStatisticsProperties.cs (new) |
| 233 | ☐ | SEQ-GC-ANALYSIS-001 | Composition | R: GC% ∈ [0,100]; P: windows tile sequence; D: deterministic | GcSkewProperties.cs |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | R: each GC% ∈ [0,100]; P: profile length = len−w+1; D: deterministic | SequenceStatisticsProperties.cs (new) |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 234 |
| ☑ Complete | 209 |
| ☐ Not started | 25 |
| New property files needed | 4 (Chromosome, Epigenetics, Oncology) |
| Existing property files to extend | 15 |
