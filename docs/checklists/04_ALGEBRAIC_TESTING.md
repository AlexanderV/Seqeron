# Checklist 04: Algebraic Testing

**Priority:** P1  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Algebraic testing проверяет выполнение алгебраических законов: identity, commutativity, associativity, involution, idempotence, round-trip (isomorphism), distributivity. Более формальный подход, чем общие property-тесты. Многие геномные алгоритмы имеют явные алгебраические законы: `complement(complement(x))=x`, `parse(serialize(x))=x`, `score(a,b)=score(b,a)`.

**Текущее покрытие:** Некоторые алгебраические свойства существуют неявно в Property-файлах (involution complement в SequenceProperties, round-trip в FastaRoundTripProperties) но не систематизированы.

**Типы законов:**
- **ID** = Identity (f(x, e) = x или f(neutral) = 0)
- **COMM** = Commutativity (f(a,b) = f(b,a))
- **ASSOC** = Associativity (f(f(a,b),c) = f(a,f(b,c)))
- **INV** = Involution (f(f(x)) = x)
- **IDEMP** = Idempotence (f(f(x)) = f(x))
- **RT** = Round-trip / Isomorphism (g(f(x)) = x)
- **DIST** = Distributivity / Conservation law
- **TRI** = Triangle inequality

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | Algebraic Laws |
|---|--------|-----------|------|---------------|
| 1 | ☐ | SEQ-GC-001 | Composition | ID: GC("")=0; IDEMP: GC(seq) same on recompute; DIST: GC(seq)=GC(complement(seq)) |
| 2 | ☐ | SEQ-COMP-001 | Composition | INV: complement(complement(x)) = x; ID: complement preserves length |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | INV: revcomp(revcomp(x)) = x; ID: revcomp preserves length |
| 4 | ☐ | SEQ-VALID-001 | Composition | IDEMP: validate(validate_result) is consistent; ID: empty string → defined result |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | ID: single nucleotide → complexity 0; IDEMP: same result on recompute |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | ID: single symbol → entropy 0; ID: uniform → max entropy = log2(|Σ|) |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | ID: empty → skew undefined or 0; DIST: skew(G-only)=1, skew(C-only)=-1 |
| 8 | ☐ | PAT-EXACT-001 | Matching | ID: pattern not in seq → count=0; ID: empty pattern → defined result |
| 9 | ☐ | PAT-APPROX-001 | Matching | ID: hamming(x,x)=0; COMM: hamming(a,b)=hamming(b,a); TRI: hamming(a,c) ≤ hamming(a,b)+hamming(b,c) |
| 10 | ☐ | PAT-APPROX-002 | Matching | ID: editDist(x,x)=0; COMM: editDist(a,b)=editDist(b,a); TRI: triangle inequality |
| 11 | ☐ | PAT-IUPAC-001 | Matching | ID: exact base match is reflexive; DIST: N matches all → N is top element |
| 12 | ☐ | PAT-PWM-001 | Matching | IDEMP: PWM from same sequences always the same; ID: all-zero PWM → no matches |
| 13 | ☐ | REP-STR-001 | Repeats | ID: no repeats → empty result; IDEMP: same seq → same microsatellites |
| 14 | ☐ | REP-TANDEM-001 | Repeats | ID: no tandem pattern → empty result; IDEMP: deterministic output |
| 15 | ☐ | REP-INV-001 | Repeats | ID: no complement in seq → no inverted repeats; DIST: revcomp(left)=right |
| 16 | ☐ | REP-DIRECT-001 | Repeats | ID: all unique chars → no direct repeats; IDEMP: deterministic |
| 17 | ☐ | REP-PALIN-001 | Repeats | ID: palindrome = revcomp(self); DIST: palindrome len is even (for DNA) |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | ID: no PAM pattern in seq → 0 sites; IDEMP: same input → same sites |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | ID: no valid PAM → no guides; IDEMP: deterministic guide generation |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | ID: 0 mismatches → max score; DIST: score decreases monotonically with mismatches |
| 21 | ☐ | PRIMER-TM-001 | MolTools | ID: Tm(empty)=0 or error; IDEMP: Tm(seq) deterministic |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | ID: seq too short → no primers; IDEMP: same config → same primers |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | ID: no self-complementary → dimer ΔG = 0; IDEMP: deterministic |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | ID: seq too short → no probes; IDEMP: deterministic |
| 25 | ☐ | PROBE-VALID-001 | MolTools | IDEMP: validate(seq) deterministic; ID: invalid probe fails consistently |
| 26 | ☐ | RESTR-FIND-001 | MolTools | ID: no recognition site → 0 sites; IDEMP: same seq + enzyme → same sites |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | DIST: sum(fragment.Length) = sequence.Length; ID: 0 cut sites → 1 fragment = full seq |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | ID: no ATG → no ORFs; DIST: ORF length divisible by 3 |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | ID: no RBS + start codon → no genes; IDEMP: deterministic |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | ID: no -10/-35 box → no promoters; IDEMP: deterministic |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | RT: parse(serialize(annotations)) = annotations; ID: empty → empty GFF |
| 32 | ☐ | KMER-COUNT-001 | K-mer | ID: sum(counts) = seqLen - k + 1; DIST: counting identity |
| 33 | ☐ | KMER-FREQ-001 | K-mer | ID: sum(frequencies) = 1.0; DIST: freq = count / total_kmers |
| 34 | ☐ | KMER-FIND-001 | K-mer | ID: k > seqLen → no results; IDEMP: deterministic |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | COMM: score(a,b)=score(b,a); ID: align(x,x) = perfect score; ID: align(x,"") = gap_penalty × len |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | ID: no common subsequence → score=0; COMM: local_score relationship (not always symmetric for alignments) |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | ID: identical overlap → max score region; IDEMP: deterministic |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | ID: single sequence → trivial alignment; IDEMP: deterministic output |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | ID: d(x,x)=0; COMM: d(a,b)=d(b,a); TRI: d(a,c)≤d(a,b)+d(b,c) |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | ID: 1 sequence → trivial tree; IDEMP: same input → same topology |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | RT: parse(toNewick(tree)) = tree; ID: leaf → simple Newick |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | ID: RF(t,t)=0; COMM: RF(a,b)=RF(b,a); TRI: RF(a,c)≤RF(a,b)+RF(b,c) |
| 43 | ☐ | POP-FREQ-001 | PopGen | ID: alleleFreq(all same)=1.0; DIST: sum(allele freqs)=1 |
| 44 | ☐ | POP-DIV-001 | PopGen | ID: identical sequences → π=0; ID: single seq → θ undefined or 0 |
| 45 | ☐ | POP-HW-001 | PopGen | DIST: p² + 2pq + q² = 1 (algebraic identity); ID: monomorphic → all one genotype |
| 46 | ☐ | POP-FST-001 | PopGen | ID: Fst(identical pops)=0; COMM: Fst(a,b)=Fst(b,a) |
| 47 | ☐ | POP-LD-001 | PopGen | ID: independent loci → D'≈0; COMM: LD(a,b)=LD(b,a) |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | ID: no TTAGGG → no telomere; IDEMP: deterministic |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | ID: no AT-rich region → no centromere signal; IDEMP: deterministic |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | ID: empty input → empty karyotype; IDEMP: deterministic |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | ID: normal diploid depth → CN=2; IDEMP: deterministic |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | COMM: synteny blocks(A,B) = reversed synteny blocks(B,A); ID: self-comparison → full synteny |
| 53 | ☐ | META-CLASS-001 | Metagenomics | ID: read from known reference → correct classification; IDEMP: deterministic |
| 54 | ☐ | META-PROF-001 | Metagenomics | DIST: sum(abundances)=1.0; ID: single species → abundance=1 |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | ID: single species → Shannon=0; ID: single species → Simpson=0 |
| 56 | ☐ | META-BETA-001 | Metagenomics | ID: dist(x,x)=0; COMM: dist(a,b)=dist(b,a) |
| 57 | ☐ | META-BIN-001 | Metagenomics | ID: single contig → single bin; IDEMP: deterministic |
| 58 | ☐ | CODON-OPT-001 | Codon | RT: translate(optimize(dna)) = translate(dna); IDEMP: optimize(optimize(x)) = optimize(x) |
| 59 | ☐ | CODON-CAI-001 | Codon | ID: all optimal codons → CAI ≈ 1; IDEMP: CAI deterministic |
| 60 | ☐ | CODON-RARE-001 | Codon | ID: all common codons → empty result; IDEMP: deterministic |
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
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | ID: empty sequence → no pairs; ID: single base → no pairs |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | ID: no complementary region → no stem-loops; IDEMP: deterministic |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | ID: no structure → ΔG=0; DIST: energy is additive for independent stacks |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | IDEMP: seed extraction deterministic; ID: seed len always 6-8 |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | ID: no seed match → score=0; IDEMP: deterministic scoring |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | ID: no hairpin → not precursor; IDEMP: deterministic |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | ID: non-GT site → score 0 or minimal; IDEMP: deterministic |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | ID: non-AG site → score 0 or minimal; IDEMP: deterministic |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | ID: no GT/AG → no introns; IDEMP: deterministic |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | ID: all ordered residues → all scores low; IDEMP: deterministic |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | ID: fully ordered protein → no disordered regions; IDEMP: deterministic |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | ID: no matching residues → no motifs; IDEMP: deterministic |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | ID: non-matching seq → no pattern hit; IDEMP: deterministic |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | ID: seq too short → no domain; IDEMP: deterministic |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | ID: no CG dinucleotides → CpG ratio = 0; DIST: O/E formula verified |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | ID: zero expression → zero infiltration; IDEMP: deterministic |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 0 |
| ☐ Not started | 86 |
| Laws verified | ~172 (≈2 per algorithm) |
