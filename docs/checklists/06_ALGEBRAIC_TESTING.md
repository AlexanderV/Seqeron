# Checklist 06: Algebraic Testing

**Priority:** P1  
**Date:** 2026-03-19  
**Total algorithms:** 234

---

## Description

Algebraic testing перевіряє виконання алгебраїчних законів: identity, commutativity, associativity, involution, idempotence, round-trip (isomorphism), distributivity. Більш формальний підхід, ніж загальні property-тести. Багато геномних алгоритмів мають явні алгебраїчні закони: `complement(complement(x))=x`, `parse(serialize(x))=x`, `score(a,b)=score(b,a)`.

**Поточне покриття:** Деякі алгебраїчні властивості існують неявно в Property-файлах (involution complement у SequenceProperties, round-trip у FastaRoundTripProperties), але не систематизовані.

**Типи законів:**
- **ID** = Identity (f(x, e) = x або f(neutral) = 0)
- **COMM** = Commutativity (f(a,b) = f(b,a))
- **ASSOC** = Associativity (f(f(a,b),c) = f(a,f(b,c)))
- **INV** = Involution (f(f(x)) = x)
- **IDEMP** = Idempotence (f(f(x)) = f(x))
- **RT** = Round-trip / Isomorphism (g(f(x)) = x)
- **DIST** = Distributivity / Conservation law
- **TRI** = Triangle inequality

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete | ✗ Not Applicable

| # | Status | Test Unit | Area | Algebraic Laws |
|---|--------|-----------|------|---------------|
| 1 | ☑ | SEQ-GC-001 | Composition | ID: GC("")=0; IDEMP: GC(seq) same on recompute; DIST: GC(seq)=GC(complement(seq)) |
| 2 | ☑ | SEQ-COMP-001 | Composition | INV: complement(complement(x)) = x; ID: complement preserves length |
| 3 | ☑ | SEQ-REVCOMP-001 | Composition | INV: revcomp(revcomp(x)) = x; ID: revcomp preserves length |
| 4 | ✗ | SEQ-VALID-001 | Composition | IDEMP: validate(validate_result) is consistent; ID: empty string → defined result |
| 5 | ✗ | SEQ-COMPLEX-001 | Composition | ID: single nucleotide → complexity 0; IDEMP: same result on recompute |
| 6 | ✗ | SEQ-ENTROPY-001 | Composition | ID: single symbol → entropy 0; ID: uniform → max entropy = log2(|Σ|) |
| 7 | ☑ | SEQ-GCSKEW-001 | Composition | ID: empty → skew undefined or 0; DIST: skew(G-only)=1, skew(C-only)=-1 |
| 8 | ✗ | PAT-EXACT-001 | Matching | ID: pattern not in seq → count=0; ID: empty pattern → defined result |
| 9 | ☑ | PAT-APPROX-001 | Matching | ID: hamming(x,x)=0; COMM: hamming(a,b)=hamming(b,a); TRI: hamming(a,c) ≤ hamming(a,b)+hamming(b,c) |
| 10 | ☑ | PAT-APPROX-002 | Matching | ID: editDist(x,x)=0; COMM: editDist(a,b)=editDist(b,a); TRI: triangle inequality |
| 11 | ☑ | PAT-IUPAC-001 | Matching | ID: exact base match is reflexive; DIST: N matches all → N is top element |
| 12 | ✗ | PAT-PWM-001 | Matching | IDEMP: PWM from same sequences always the same; ID: all-zero PWM → no matches |
| 13 | ✗ | REP-STR-001 | Repeats | ID: no repeats → empty result; IDEMP: same seq → same microsatellites |
| 14 | ✗ | REP-TANDEM-001 | Repeats | ID: no tandem pattern → empty result; IDEMP: deterministic output |
| 15 | ☐ | REP-INV-001 | Repeats | ID: no complement in seq → no inverted repeats; DIST: revcomp(left)=right |
| 16 | ✗ | REP-DIRECT-001 | Repeats | ID: all unique chars → no direct repeats; IDEMP: deterministic |
| 17 | ☐ | REP-PALIN-001 | Repeats | ID: palindrome = revcomp(self); DIST: palindrome len is even (for DNA) |
| 18 | ✗ | CRISPR-PAM-001 | MolTools | ID: no PAM pattern in seq → 0 sites; IDEMP: same input → same sites |
| 19 | ✗ | CRISPR-GUIDE-001 | MolTools | ID: no valid PAM → no guides; IDEMP: deterministic guide generation |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | ID: 0 mismatches → max score; DIST: score decreases monotonically with mismatches |
| 21 | ✗ | PRIMER-TM-001 | MolTools | ID: Tm(empty)=0 or error; IDEMP: Tm(seq) deterministic |
| 22 | ✗ | PRIMER-DESIGN-001 | MolTools | ID: seq too short → no primers; IDEMP: same config → same primers |
| 23 | ✗ | PRIMER-STRUCT-001 | MolTools | ID: no self-complementary → dimer ΔG = 0; IDEMP: deterministic |
| 24 | ✗ | PROBE-DESIGN-001 | MolTools | ID: seq too short → no probes; IDEMP: deterministic |
| 25 | ✗ | PROBE-VALID-001 | MolTools | IDEMP: validate(seq) deterministic; ID: invalid probe fails consistently |
| 26 | ✗ | RESTR-FIND-001 | MolTools | ID: no recognition site → 0 sites; IDEMP: same seq + enzyme → same sites |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | DIST: sum(fragment.Length) = sequence.Length; ID: 0 cut sites → 1 fragment = full seq |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | ID: no ATG → no ORFs; DIST: ORF length divisible by 3 |
| 29 | ✗ | ANNOT-GENE-001 | Annotation | ID: no RBS + start codon → no genes; IDEMP: deterministic |
| 30 | ✗ | ANNOT-PROM-001 | Annotation | ID: no -10/-35 box → no promoters; IDEMP: deterministic |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | RT: parse(serialize(annotations)) = annotations; ID: empty → empty GFF |
| 32 | ☐ | KMER-COUNT-001 | K-mer | ID: sum(counts) = seqLen - k + 1; DIST: counting identity |
| 33 | ☐ | KMER-FREQ-001 | K-mer | ID: sum(frequencies) = 1.0; DIST: freq = count / total_kmers |
| 34 | ✗ | KMER-FIND-001 | K-mer | ID: k > seqLen → no results; IDEMP: deterministic |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | COMM: score(a,b)=score(b,a); ID: align(x,x) = perfect score; ID: align(x,"") = gap_penalty × len |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | ID: no common subsequence → score=0; COMM: local_score relationship (not always symmetric for alignments) |
| 37 | ✗ | ALIGN-SEMI-001 | Alignment | ID: identical overlap → max score region; IDEMP: deterministic |
| 38 | ✗ | ALIGN-MULTI-001 | Alignment | ID: single sequence → trivial alignment; IDEMP: deterministic output |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | ID: d(x,x)=0; COMM: d(a,b)=d(b,a); TRI: d(a,c)≤d(a,b)+d(b,c) |
| 40 | ✗ | PHYLO-TREE-001 | Phylogenetic | ID: 1 sequence → trivial tree; IDEMP: same input → same topology |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | RT: parse(toNewick(tree)) = tree; ID: leaf → simple Newick |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | ID: RF(t,t)=0; COMM: RF(a,b)=RF(b,a); TRI: RF(a,c)≤RF(a,b)+RF(b,c) |
| 43 | ☐ | POP-FREQ-001 | PopGen | ID: alleleFreq(all same)=1.0; DIST: sum(allele freqs)=1 |
| 44 | ✗ | POP-DIV-001 | PopGen | ID: identical sequences → π=0; ID: single seq → θ undefined or 0 |
| 45 | ☐ | POP-HW-001 | PopGen | DIST: p² + 2pq + q² = 1 (algebraic identity); ID: monomorphic → all one genotype |
| 46 | ☐ | POP-FST-001 | PopGen | ID: Fst(identical pops)=0; COMM: Fst(a,b)=Fst(b,a) |
| 47 | ☐ | POP-LD-001 | PopGen | ID: independent loci → D'≈0; COMM: LD(a,b)=LD(b,a) |
| 48 | ✗ | CHROM-TELO-001 | Chromosome | ID: no TTAGGG → no telomere; IDEMP: deterministic |
| 49 | ✗ | CHROM-CENT-001 | Chromosome | ID: no AT-rich region → no centromere signal; IDEMP: deterministic |
| 50 | ✗ | CHROM-KARYO-001 | Chromosome | ID: empty input → empty karyotype; IDEMP: deterministic |
| 51 | ✗ | CHROM-ANEU-001 | Chromosome | ID: normal diploid depth → CN=2; IDEMP: deterministic |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | COMM: synteny blocks(A,B) = reversed synteny blocks(B,A); ID: self-comparison → full synteny |
| 53 | ✗ | META-CLASS-001 | Metagenomics | ID: read from known reference → correct classification; IDEMP: deterministic |
| 54 | ☐ | META-PROF-001 | Metagenomics | DIST: sum(abundances)=1.0; ID: single species → abundance=1 |
| 55 | ✗ | META-ALPHA-001 | Metagenomics | ID: single species → Shannon=0; ID: single species → Simpson=0 |
| 56 | ☐ | META-BETA-001 | Metagenomics | ID: dist(x,x)=0; COMM: dist(a,b)=dist(b,a) |
| 57 | ✗ | META-BIN-001 | Metagenomics | ID: single contig → single bin; IDEMP: deterministic |
| 58 | ☐ | CODON-OPT-001 | Codon | RT: translate(optimize(dna)) = translate(dna); IDEMP: optimize(optimize(x)) = optimize(x) |
| 59 | ✗ | CODON-CAI-001 | Codon | ID: all optimal codons → CAI ≈ 1; IDEMP: CAI deterministic |
| 60 | ✗ | CODON-RARE-001 | Codon | ID: all common codons → empty result; IDEMP: deterministic |
| 61 | ☐ | CODON-USAGE-001 | Codon | DIST: sum per amino acid = 1.0; ID: single codon seq → usage=1 for that codon |
| 62 | ☐ | TRANS-CODON-001 | Translation | ID: genetic code is complete (64 entries); DIST: each AA mapped from ≥1 codon |
| 63 | ☐ | TRANS-PROT-001 | Translation | DIST: protein.len ≤ dna.len/3; ID: ATG → M; ID: stop codon → terminates |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | RT: parse(write(records)) = records; ID: empty file → 0 records |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | RT: parse(write(records)) = records; ID: empty → 0 records |
| 66 | ☐ | PARSE-BED-001 | FileIO | RT: parse(write(regions)) = regions; ID: empty → 0 regions |
| 67 | ☐ | PARSE-VCF-001 | FileIO | RT: parse(write(variants)) = variants; ID: header only → 0 variants |
| 68 | ☐ | PARSE-GFF-001 | FileIO | RT: parse(write(features)) = features; ID: empty → 0 features |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | RT: parse(write(record)) = record; ID: minimal record → valid |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | RT: parse(write(record)) = record; ID: minimal record → valid |
| 71 | ✗ | RNA-STRUCT-001 | RnaStructure | ID: empty sequence → no pairs; ID: single base → no pairs |
| 72 | ✗ | RNA-STEMLOOP-001 | RnaStructure | ID: no complementary region → no stem-loops; IDEMP: deterministic |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | ID: no structure → ΔG=0; DIST: energy is additive for independent stacks |
| 74 | ✗ | MIRNA-SEED-001 | MiRNA | IDEMP: seed extraction deterministic; ID: seed len always 6-8 |
| 75 | ✗ | MIRNA-TARGET-001 | MiRNA | ID: no seed match → score=0; IDEMP: deterministic scoring |
| 76 | ✗ | MIRNA-PRECURSOR-001 | MiRNA | ID: no hairpin → not precursor; IDEMP: deterministic |
| 77 | ✗ | SPLICE-DONOR-001 | Splicing | ID: non-GT site → score 0 or minimal; IDEMP: deterministic |
| 78 | ✗ | SPLICE-ACCEPTOR-001 | Splicing | ID: non-AG site → score 0 or minimal; IDEMP: deterministic |
| 79 | ✗ | SPLICE-PREDICT-001 | Splicing | ID: no GT/AG → no introns; IDEMP: deterministic |
| 80 | ✗ | DISORDER-PRED-001 | ProteinPred | ID: all ordered residues → all scores low; IDEMP: deterministic |
| 81 | ✗ | DISORDER-REGION-001 | ProteinPred | ID: fully ordered protein → no disordered regions; IDEMP: deterministic |
| 82 | ✗ | PROTMOTIF-FIND-001 | ProteinMotif | ID: no matching residues → no motifs; IDEMP: deterministic |
| 83 | ✗ | PROTMOTIF-PROSITE-001 | ProteinMotif | ID: non-matching seq → no pattern hit; IDEMP: deterministic |
| 84 | ✗ | PROTMOTIF-DOMAIN-001 | ProteinMotif | ID: seq too short → no domain; IDEMP: deterministic |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | ID: no CG dinucleotides → CpG ratio = 0; DIST: O/E formula verified |
| 86 | ✗ | ONCO-IMMUNE-001 | Oncology | ID: zero expression → zero infiltration; IDEMP: deterministic |
| 87 | ✗ | ONCO-SOMATIC-001 | Oncology | ID: empty tumor reads → no calls; IDEMP: re-calling the call set is stable |
| 88 | ☐ | ONCO-VAF-001 | Oncology | ID: alt=0 → VAF=0; HOMO: VAF(k·ref,k·alt)=VAF(ref,alt) (scale-invariant) |
| 89 | ✗ | ONCO-DRIVER-001 | Oncology | ID: empty mutation set → empty driver set; IDEMP: deterministic |
| 90 | ✗ | ONCO-ARTIFACT-001 | Oncology | IDEMP: filter(filter(x))=filter(x); ID: artifact-free input → unchanged |
| 91 | ✗ | ONCO-ANNOT-001 | Oncology | IDEMP: annotate∘annotate = annotate; ID: empty variants → empty |
| 92 | ☐ | ONCO-TMB-001 | Oncology | ID: zero mutations → TMB=0; HOMO: TMB scales inversely with panel-Mb |
| 93 | ✗ | ONCO-MSI-001 | Oncology | ID: all-stable loci → score 0; IDEMP: deterministic |
| 94 | ✗ | ONCO-HRD-001 | Oncology | ID: no scars → HRD=0; ADD: HRD = LOH + TAI + LST (additive decomposition) |
| 95 | ✗ | ONCO-LOH-001 | Oncology | ID: balanced BAF → no LOH; INVOL: A/B label swap applied twice = identity |
| 96 | ✗ | ONCO-SIG-001 | Oncology | ID: no SNVs → zero 96-vector; ADD: catalogue(A∪B)=catalogue(A)+catalogue(B) |
| 97 | ☐ | ONCO-SIG-002 | Oncology | ID: zero catalogue → zero exposures; HOMO: exposure(k·catalogue)=k·exposure |
| 98 | ✗ | ONCO-SIG-003 | Oncology | IDEMP: same seed → identical CI; ID: zero reps → degenerate CI |
| 99 | ✗ | ONCO-SIG-004 | Oncology | IDEMP: deterministic; ID: single signature → that process dominant |
| 100 | ✗ | ONCO-FUSION-001 | Oncology | ID: no chimeric reads → no fusions; IDEMP: deterministic |
| 101 | ✗ | ONCO-FUSION-002 | Oncology | IDEMP: match∘match = match; ID: empty DB → no matches |
| 102 | ✗ | ONCO-FUSION-003 | Oncology | IDEMP: deterministic; ID: intragenic breakpoint → no fusion frame |
| 103 | ✗ | ONCO-CNA-001 | Oncology | ID: log2=0 → CN=2 (neutral); IDEMP: deterministic |
| 104 | ✗ | ONCO-CNA-002 | Oncology | IDEMP: deterministic; ID: neutral genome → no focal amplifications |
| 105 | ✗ | ONCO-CNA-003 | Oncology | IDEMP: deterministic; ID: neutral genome → no homozygous deletions |
| 106 | ✗ | ONCO-PURITY-001 | Oncology | ID: all-germline VAFs → purity at boundary; IDEMP: deterministic |
| 107 | ✗ | ONCO-PLOIDY-001 | Oncology | ID: diploid genome → ploidy=2; IDEMP: deterministic |
| 108 | ✗ | ONCO-CLONAL-001 | Oncology | IDEMP: deterministic; ID: CCF≈1 → clonal |
| 109 | ✗ | ONCO-NEO-001 | Oncology | ID: silent mutation → no neoantigen; IDEMP: deterministic |
| 110 | ✗ | ONCO-MHC-001 | Oncology | IDEMP: deterministic; MONO: binding class order-preserving on IC50 |
| 111 | ☐ | ONCO-CTDNA-001 | Oncology | ID: no tumor reads → fraction 0; HOMO: scale-invariant to total depth |
| 112 | ✗ | ONCO-MRD-001 | Oncology | ID: no tracked variants → MRD negative; IDEMP: deterministic |
| 113 | ✗ | ONCO-CHIP-001 | Oncology | IDEMP: filter∘filter = filter; ID: no CHIP genes → unchanged |
| 114 | ✗ | ONCO-PHYLO-001 | Oncology | ID: single clone → trivial tree; IDEMP: deterministic |
| 115 | ☐ | ONCO-CCF-001 | Oncology | ID: VAF=0 → CCF=0; HOMO: CCF linear in VAF at fixed CN/purity |
| 116 | ✗ | ONCO-HETERO-001 | Oncology | ID: single VAF → heterogeneity 0; IDEMP: deterministic |
| 117 | ✗ | ONCO-HLA-001 | Oncology | IDEMP: deterministic; INVOL: normalise(normalise(allele))=normalise(allele) |
| 118 | ✗ | ONCO-ACTION-001 | Oncology | IDEMP: deterministic; ID: no evidence → lowest tier |
| 119 | ✗ | ONCO-SV-001 | Oncology | IDEMP: deterministic; ID: no breakpoints → no rearrangement |
| 120 | ✗ | ONCO-EXPR-001 | Oncology | ID: zero variance → no outliers; INVOL: z-score sign flips under expression negation |
| 121 | ✗ | SEQ-COMPOSITION-001 | Statistics | ID: empty → all-zero; ADD: composition(a+b)=composition(a)+composition(b) |
| 122 | ✗ | SEQ-DINUC-001 | Statistics | ID: length<2 → empty; ADD: counts additive on gap-free concatenation |
| 123 | ✗ | SEQ-HYDRO-001 | Statistics | ID: empty → 0; IDEMP: deterministic |
| 124 | ✗ | SEQ-MW-001 | Statistics | ID: empty → 0; ADD: MW additive (minus water per peptide bond) |
| 125 | ✗ | SEQ-PI-001 | Statistics | IDEMP: deterministic; INVOL: charge symmetry around pI |
| 126 | ✗ | SEQ-SECSTRUCT-001 | Statistics | IDEMP: deterministic; ID: empty → empty |
| 127 | ✗ | SEQ-STATS-001 | Statistics | ID: empty → zeros; ADD: counts additive on concatenation |
| 128 | ✗ | SEQ-SUMMARY-001 | Statistics | ID: empty → zero length; IDEMP: deterministic |
| 129 | ✗ | SEQ-THERMO-001 | Statistics | ID: empty → 0; IDEMP: deterministic |
| 130 | ☐ | SEQ-TM-001 | Statistics | ID: empty → 0; HOMO: homopolymer Tm scales linearly with length (Wallace) |
| 131 | ☐ | COMPGEN-ANI-001 | Comparative | ID: ANI(A,A)=100; COMM: ANI symmetric |
| 132 | ✗ | COMPGEN-CLUSTER-001 | Comparative | IDEMP: deterministic; ID: single genome → trivial clusters |
| 133 | ☐ | COMPGEN-COMPARE-001 | Comparative | COMM: symmetric; IDEMP: deterministic |
| 134 | ✗ | COMPGEN-DOTPLOT-001 | Comparative | IDEMP: deterministic; ID: A vs A → main diagonal |
| 135 | ☐ | COMPGEN-ORTHO-001 | Comparative | COMM: symmetric relation; IDEMP: deterministic |
| 136 | ☐ | COMPGEN-RBH-001 | Comparative | COMM: symmetric; IDEMP: RBH∘RBH = RBH |
| 137 | ✗ | COMPGEN-REARR-001 | Comparative | ID: identical genomes → empty; IDEMP: deterministic |
| 138 | ☐ | COMPGEN-REVERSAL-001 | Comparative | ID: d(A,A)=0; COMM: symmetric; TRIANGLE: d(A,C) ≤ d(A,B)+d(B,C) |
| 139 | ✗ | COMPGEN-SYNTENY-001 | Comparative | IDEMP: deterministic; ID: A vs A → whole-genome block |
| 140 | ✗ | ASSEMBLY-CONSENSUS-001 | Assembly | ID: single read → that read; IDEMP: consensus∘consensus = consensus |
| 141 | ✗ | ASSEMBLY-CORRECT-001 | Assembly | IDEMP: correct∘correct = correct; ID: clean reads unchanged |
| 142 | ✗ | ASSEMBLY-COVER-001 | Assembly | ID: no reads → 0; ADD: coverage additive |
| 143 | ✗ | ASSEMBLY-DBG-001 | Assembly | IDEMP: deterministic; ID: single read → single contig |
| 144 | ✗ | ASSEMBLY-MERGE-001 | Assembly | IDEMP: merge∘merge = merge; ID: non-overlapping contigs unchanged |
| 145 | ✗ | ASSEMBLY-OLC-001 | Assembly | IDEMP: deterministic; ID: single read → itself |
| 146 | ✗ | ASSEMBLY-SCAFFOLD-001 | Assembly | IDEMP: deterministic; ID: no links → unchanged contigs |
| 147 | ✗ | ASSEMBLY-STATS-001 | Assembly | ID: empty → 0; IDEMP: deterministic |
| 148 | ✗ | ASSEMBLY-TRIM-001 | Assembly | IDEMP: trim∘trim = trim; ID: high-quality read unchanged |
| 149 | ☐ | RNA-DOTBRACKET-001 | RnaStructure | RT: parse∘format = identity; ID: empty → no pairs |
| 150 | ✗ | RNA-HAIRPIN-001 | RnaStructure | IDEMP: deterministic; ID: loop < minLoop → no hairpin |
| 151 | ✗ | RNA-INVERT-001 | RnaStructure | INVOL: revcomp of arm twice = identity; ID: no complementarity → none |
| 152 | ✗ | RNA-MFE-001 | RnaStructure | ID: unpaired sequence → 0; IDEMP: deterministic |
| 153 | ☐ | RNA-PAIR-001 | RnaStructure | COMM: canPair symmetric; IDEMP: deterministic |
| 154 | ✗ | RNA-PARTITION-001 | RnaStructure | ID: single base → Z = 1; IDEMP: deterministic |
| 155 | ✗ | RNA-PSEUDOKNOT-001 | RnaStructure | ID: nested structure → none; IDEMP: deterministic |
| 156 | ✗ | KMER-ASYNC-001 | K-mer | ID: empty → empty; IDEMP: async = sync |
| 157 | ☐ | KMER-BOTH-001 | K-mer | ADD: counts additive; COMM: strand-symmetric |
| 158 | ☐ | KMER-DIST-001 | K-mer | ID: d(x,x)=0; COMM: symmetric; TRIANGLE: d(a,c) ≤ d(a,b)+d(b,c) |
| 159 | ✗ | KMER-GENERATE-001 | K-mer | ID: k=0 → empty; IDEMP: deterministic |
| 160 | ✗ | KMER-POSITIONS-001 | K-mer | ID: absent k-mer → empty; IDEMP: deterministic |
| 161 | ✗ | KMER-STATS-001 | K-mer | ID: empty → zeros; ADD: counts additive on concatenation |
| 162 | ✗ | KMER-UNIQUE-001 | K-mer | IDEMP: deterministic; ID: all-distinct → all unique |
| 163 | ✗ | PROTMOTIF-CC-001 | ProteinMotif | ID: empty → 0; IDEMP: deterministic |
| 164 | ✗ | PROTMOTIF-COMMON-001 | ProteinMotif | IDEMP: deterministic; ID: single sequence → trivial |
| 165 | ✗ | PROTMOTIF-LC-001 | ProteinMotif | IDEMP: deterministic; ID: high-complexity → none |
| 166 | ✗ | PROTMOTIF-PATTERN-001 | ProteinMotif | ID: empty pattern → none; IDEMP: deterministic |
| 167 | ✗ | PROTMOTIF-SP-001 | ProteinMotif | IDEMP: deterministic; ID: no hydrophobic N-term → none |
| 168 | ✗ | PROTMOTIF-TM-001 | ProteinMotif | IDEMP: deterministic; ID: hydrophilic → none |
| 169 | ✗ | MOTIF-CONS-001 | Matching | IDEMP: deterministic; ID: single row → that row |
| 170 | ✗ | MOTIF-DISCOVER-001 | Matching | IDEMP: deterministic; ID: no recurrence → none |
| 171 | ✗ | MOTIF-GENERATE-001 | Matching | IDEMP: deterministic; ID: single column → its base |
| 172 | ✗ | MOTIF-REGULATORY-001 | Matching | IDEMP: deterministic; ID: empty → none |
| 173 | ✗ | MOTIF-SHARED-001 | Matching | IDEMP: shared∘shared = shared; ID: single input → all its motifs |
| 174 | ✗ | PAT-APPROX-003 | Matching | ID: exact present → 0; IDEMP: deterministic |
| 175 | ✗ | GENOMIC-COMMON-001 | Analysis | IDEMP: common∘common = common; ID: single input → itself |
| 176 | ✗ | GENOMIC-MOTIFS-001 | Analysis | IDEMP: deterministic; ID: empty → none |
| 177 | ✗ | GENOMIC-ORF-001 | Analysis | ID: no ATG → none; IDEMP: deterministic |
| 178 | ✗ | GENOMIC-REPEAT-001 | Analysis | IDEMP: deterministic; ID: no repeat → none |
| 179 | ☐ | GENOMIC-SIMILARITY-001 | Analysis | ID: sim(x,x)=1; COMM: symmetric |
| 180 | ✗ | GENOMIC-TANDEM-001 | Analysis | IDEMP: deterministic; ID: no tandem → none |
| 181 | ✗ | EPIGEN-AGE-001 | Epigenetics | ID: no methylation → baseline; IDEMP: deterministic |
| 182 | ✗ | EPIGEN-BISULF-001 | Epigenetics | IDEMP: bisulfite∘bisulfite = bisulfite; ID: no C → unchanged |
| 183 | ✗ | EPIGEN-CHROM-001 | Epigenetics | IDEMP: deterministic; ID: empty → none |
| 184 | ☐ | EPIGEN-DMR-001 | Epigenetics | ID: identical methylomes → no DMR; COMM: |Δ| symmetric |
| 185 | ✗ | EPIGEN-METHYL-001 | Epigenetics | ID: no reads → 0; ADD: counts additive |
| 186 | ✗ | VARIANT-ANNOT-001 | Variants | IDEMP: annotate∘annotate = annotate; ID: empty → empty |
| 187 | ✗ | VARIANT-CALL-001 | Variants | ID: tumor = ref → no calls; IDEMP: deterministic |
| 188 | ✗ | VARIANT-INDEL-001 | Variants | ID: no indel → none; IDEMP: deterministic |
| 189 | ✗ | VARIANT-SNP-001 | Variants | ID: ref = alt → none; IDEMP: deterministic |
| 190 | ✗ | PANGEN-CLUSTER-001 | PanGenome | IDEMP: cluster∘cluster = cluster; ID: single gene → singleton |
| 191 | ✗ | PANGEN-CORE-001 | PanGenome | ID: single genome → core = that genome; IDEMP: deterministic |
| 192 | ✗ | PANGEN-HEAP-001 | PanGenome | IDEMP: deterministic; ID: single genome → degenerate fit |
| 193 | ✗ | PANGEN-MARKER-001 | PanGenome | SUB: markers ⊆ core; IDEMP: deterministic |
| 194 | ✗ | META-FUNC-001 | Metagenomics | IDEMP: deterministic; ID: empty → none |
| 195 | ✗ | META-PATHWAY-001 | Metagenomics | ID: no pathway genes → p=1; IDEMP: deterministic |
| 196 | ✗ | META-RESIST-001 | Metagenomics | IDEMP: deterministic; ID: empty → none |
| 197 | ✗ | META-TAXA-001 | Metagenomics | ID: no difference → p=1; IDEMP: deterministic |
| 198 | ✗ | TRANS-DIFF-001 | Transcriptome | ID: A = B → FC = 0; INVOL: FC(A,B) = −FC(B,A) |
| 199 | ☐ | TRANS-EXPR-001 | Transcriptome | HOMO: scaling depth → same TPM; ID: empty → 0 |
| 200 | ✗ | TRANS-SPLICE-001 | Transcriptome | ID: single isoform → PSI = 1; IDEMP: deterministic |
| 201 | ✗ | SV-BREAKPOINT-001 | StructuralVar | ID: identical → none; IDEMP: deterministic |
| 202 | ✗ | SV-CNV-001 | StructuralVar | ID: ratio = 1 → CN = 2; IDEMP: deterministic |
| 203 | ✗ | SV-DETECT-001 | StructuralVar | ID: identical genomes → none; IDEMP: deterministic |
| 204 | ✗ | DISORDER-LC-001 | ProteinPred | IDEMP: deterministic; ID: high-complexity → none |
| 205 | ✗ | DISORDER-MORF-001 | ProteinPred | IDEMP: deterministic; ID: fully ordered → none |
| 206 | ✗ | DISORDER-PROPENSITY-001 | ProteinPred | ID: empty → empty; IDEMP: deterministic |
| 207 | ✗ | POP-ANCESTRY-001 | PopGen | ID: single population → 100%; P: proportions sum to 1 |
| 208 | ✗ | POP-ROH-001 | PopGen | IDEMP: deterministic; ID: all-heterozygous → none |
| 209 | ✗ | POP-SELECT-001 | PopGen | ID: neutral → baseline; IDEMP: deterministic |
| 210 | ✗ | SEQ-ATSKEW-001 | Composition | ID: balanced AT → 0; INVOL: complement negates skew |
| 211 | ✗ | SEQ-REPLICATION-001 | Composition | IDEMP: deterministic; ID: flat skew → midpoint |
| 212 | ✗ | SEQ-RNACOMP-001 | Composition | INVOL: complement∘complement = identity; ID: defined per base |
| 213 | ✗ | CODON-ENC-001 | Codon | ID: uniform usage → ENC = 61; IDEMP: deterministic |
| 214 | ✗ | CODON-RSCU-001 | Codon | ID: uniform usage → RSCU = 1; P: per-AA mean = 1 |
| 215 | ✗ | CODON-STATS-001 | Codon | ID: empty → zeros; ADD: counts additive on concatenation |
| 216 | ✗ | ANNOT-CODING-001 | Annotation | IDEMP: deterministic; ID: random sequence → low score |
| 217 | ✗ | ANNOT-CODONUSAGE-001 | Annotation | ID: empty → empty; P: per-AA sum = 1 |
| 218 | ✗ | ANNOT-REPEAT-001 | Annotation | IDEMP: deterministic; ID: no repeat → none |
| 219 | ☐ | QUALITY-PHRED-001 | Quality | RT: encode∘decode = identity; ID: defined per char |
| 220 | ✗ | QUALITY-STATS-001 | Quality | ID: empty → 0; ADD: counts additive |
| 221 | ✗ | PHYLO-BOOT-001 | Phylogenetic | IDEMP: same seed → same support; ID: single tree → trivial |
| 222 | ✗ | PHYLO-STATS-001 | Phylogenetic | ID: single leaf → depth 0; IDEMP: deterministic |
| 223 | ✗ | TRANS-SIXFRAME-001 | Translation | INVOL: frames symmetric under revcomp; ID: defined for all 6 |
| 224 | ✗ | RESTR-FILTER-001 | MolTools | IDEMP: filter∘filter = filter; SUB: filtered ⊆ all |
| 225 | ✗ | MIRNA-PAIR-001 | MiRNA | IDEMP: deterministic; ID: no complementarity → no pairing |
| 226 | ☐ | ALIGN-STATS-001 | Alignment | ID: identity(x,x)=1; COMM: symmetric |
| 227 | ☐ | SEQ-CODON-FREQ-001 | Statistics | ID: freq("")=∅; IDEMP: deterministic; DIST: Σ counts = len/3 |
| 228 | ✗ | SEQ-COMPLEX-COMPRESS-001 | Complexity | IDEMP: deterministic; ID: homopolymer → minimal ratio |
| 229 | ✗ | SEQ-COMPLEX-DUST-001 | Complexity | IDEMP: deterministic; INVAR: DUST(complement(x))=DUST(x) |
| 230 | ✗ | SEQ-COMPLEX-KMER-001 | Complexity | ID: entropy(homopolymer)=0; IDEMP: deterministic; INVAR: reverse-invariant |
| 231 | ☐ | SEQ-COMPLEX-WINDOW-001 | Complexity | IDEMP: deterministic; DIST: window count = len−w+1 |
| 232 | ☐ | SEQ-ENTROPY-PROFILE-001 | Statistics | IDEMP: deterministic; INVAR: complement-invariant; DIST: length = len−w+1 |
| 233 | ☐ | SEQ-GC-ANALYSIS-001 | Composition | ID: GC("")=0; IDEMP: deterministic; DIST: GC(seq)=GC(complement(seq)) |
| 234 | ☐ | SEQ-GC-PROFILE-001 | Statistics | IDEMP: deterministic; INVAR: complement-invariant; DIST: length = len−w+1 |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 234 |
| ☑ Complete | 7 |
| ☐ Not started | 58 |
| ✗ Not applicable | 169 |
| Laws verified | ~172 (≈2 per algorithm) |
