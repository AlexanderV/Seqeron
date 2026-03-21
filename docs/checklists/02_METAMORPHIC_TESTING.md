# Checklist 02: Metamorphic Testing

**Priority:** P0  
**Date:** 2026-03-19  
**Total algorithms:** 86

---

## Description

Metamorphic testing решает «проблему оракула» — когда для произвольного входа нет эталонного ответа. Вместо проверки конкретного результата проверяются **метаморфные отношения (MR)** — свойства, связывающие выходы нескольких запусков при трансформации входа. Идеально для биоинформатических алгоритмов поиска, выравнивания, скоринга и статистики.

**Текущее покрытие:** `MetamorphicTests.cs` — 7 юнитов: PAT-IUPAC-001, PAT-PWM-001, PAT-APPROX-002, REP-STR-001, REP-TANDEM-001, REP-INV-001, REP-DIRECT-001 (18+ MR).

**Типы отношений:**
- **SUB** = Subset/Superset (расширение параметров даёт надмножество)
- **MON** = Monotonicity (увеличение X → увеличение/уменьшение Y)
- **INV** = Invariance (трансформация входа не меняет выход)
- **SYM** = Symmetry (f(a,b) = f(b,a))
- **COMP** = Composition (f ⊆ g или f ∘ g = h)
- **SHIFT** = Positional shift (добавление фланка сдвигает позиции)

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | MR Relations |
|---|--------|-----------|------|-------------|
| 1 | ☐ | SEQ-GC-001 | Composition | INV: complement preserves GC%; INV: shuffle preserves GC%; INV: case-insensitive |
| 2 | ☐ | SEQ-COMP-001 | Composition | INV: complement(complement(x))=x (involution); INV: length unchanged |
| 3 | ☐ | SEQ-REVCOMP-001 | Composition | INV: revcomp(revcomp(x))=x; INV: length unchanged |
| 4 | ☐ | SEQ-VALID-001 | Composition | INV: case conversion preserves validity; COMP: valid DNA ⊂ valid IUPAC; INV: repeat seq → same result |
| 5 | ☐ | SEQ-COMPLEX-001 | Composition | INV: permutation preserves complexity; MON: homopolymer → min complexity; MON: random → higher |
| 6 | ☐ | SEQ-ENTROPY-001 | Composition | INV: permutation preserves Shannon entropy; MON: uniform → max entropy; MON: single symbol → 0 |
| 7 | ☐ | SEQ-GCSKEW-001 | Composition | SYM: complement flips skew sign; INV: reverse flips cumulative skew; INV: all-G → max positive |
| 8 | ☐ | PAT-EXACT-001 | Matching | SHIFT: prepend flank shifts positions by flank.len; COMP: exact ⊆ hamming(maxDist=0); INV: duplicate → doubled count |
| 9 | ☐ | PAT-APPROX-001 | Matching | SYM: hamming(a,b)=hamming(b,a); MON: higher maxDist → ≥ matches; COMP: exact ⊆ approx(d=0) |
| 10 | ☑ | PAT-APPROX-002 | Matching | MON: higher maxEdits → superset; COMP: exact ⊆ approximate; SYM: editDist(a,b)=editDist(b,a); INV: prefix extension; R: non-negativity |
| 11 | ☑ | PAT-IUPAC-001 | Matching | SUB: degeneracy hierarchy (N ⊇ all); INV: shift-invariance; MON: more degenerate → ≥ matches |
| 12 | ☑ | PAT-PWM-001 | Matching | INV: permuted training set → same PWM; INV: duplicated set → same consensus; SUB: lower threshold → ≥ matches |
| 13 | ☑ | REP-STR-001 | Repeats | MON: lower minRepeats → superset; MON: wider unit range → superset; MON: doubled repeat → increased count |
| 14 | ☑ | REP-TANDEM-001 | Repeats | SHIFT: prepend flank shifts positions; MON: higher minReps → subset; MON: higher minUnit → subset |
| 15 | ☑ | REP-INV-001 | Repeats | SYM: right arm = revcomp(left arm); INV: non-overlapping flank doesn't change core |
| 16 | ☑ | REP-DIRECT-001 | Repeats | INV: duplication → consistent results; MON: lower minLen → ≥ results; SHIFT: flank shifts positions |
| 17 | ☐ | REP-PALIN-001 | Repeats | INV: palindrome = revcomp of self; MON: wider len range → ≥ palindromes; SHIFT: flank shifts positions |
| 18 | ☐ | CRISPR-PAM-001 | MolTools | MON: longer sequence → ≥ PAM sites; INV: non-PAM region append → same count; SHIFT: flank shifts positions |
| 19 | ☐ | CRISPR-GUIDE-001 | MolTools | SUB: guide from PAM site ⊂ valid guides; MON: stricter scoring → ≤ guides; INV: downstream change → same guide |
| 20 | ☐ | CRISPR-OFF-001 | MolTools | MON: more mismatches → lower off-target score; MON: seed region mismatch penalized more; COMP: 0 mismatches → max score |
| 21 | ☐ | PRIMER-TM-001 | MolTools | MON: add GC → Tm increases; MON: add AT → Tm decreases; INV: same sequence → same Tm |
| 22 | ☐ | PRIMER-DESIGN-001 | MolTools | MON: wider Tm range → ≥ primers; SUB: stricter GC% → ⊆ results; INV: longer template → ≥ candidates |
| 23 | ☐ | PRIMER-STRUCT-001 | MolTools | MON: more self-complementary → higher dimer score; INV: non-complementary extension → same hairpin |
| 24 | ☐ | PROBE-DESIGN-001 | MolTools | MON: wider Tm → ≥ probes; SUB: stricter uniqueness → ⊆ results; INV: unrelated region append → same probes |
| 25 | ☐ | PROBE-VALID-001 | MolTools | MON: lower specificity threshold → more pass; INV: same input → same result |
| 26 | ☐ | RESTR-FIND-001 | MolTools | SHIFT: prepend flank shifts positions; MON: more enzymes → ≥ total sites; INV: non-site append → same sites |
| 27 | ☐ | RESTR-DIGEST-001 | MolTools | COMP: 0 sites → 1 fragment = full seq; MON: more enzymes → ≥ fragments; INV: fragment sum = seq length |
| 28 | ☐ | ANNOT-ORF-001 | Annotation | MON: lower minLen → ≥ ORFs; SHIFT: prepend shifts positions; INV: non-coding insert doesn't change upstream |
| 29 | ☐ | ANNOT-GENE-001 | Annotation | COMP: gene ⊃ ORF; INV: non-coding insertion doesn't affect upstream; MON: longer seq → ≥ genes |
| 30 | ☐ | ANNOT-PROM-001 | Annotation | MON: lower score threshold → ≥ promoters; SHIFT: prepend shifts; INV: downstream change → same promoter |
| 31 | ☐ | ANNOT-GFF-001 | Annotation | INV: re-serialization preserves line count; COMP: parse(write(x))=x; INV: attribute order irrelevant |
| 32 | ☐ | KMER-COUNT-001 | K-mer | INV: total k-mer instances = seqLen - k + 1; MON: k+1 → ≤ distinct k-mers; INV: reverse → same set (canonical) |
| 33 | ☐ | KMER-FREQ-001 | K-mer | INV: duplicate seq → same freqs; SYM: complement has related frequency profile; INV: sum(freq)=1 |
| 34 | ☐ | KMER-FIND-001 | K-mer | MON: lower minFreq → ≥ k-mers; SUB: top-1 ⊆ top-5; INV: repeat seq → high frequency for unit k-mer |
| 35 | ☐ | ALIGN-GLOBAL-001 | Alignment | SYM: score(a,b)=score(b,a); MON: more matches → higher score; COMP: identity → max score; INV: gap-only insert → score change = gap penalty |
| 36 | ☐ | ALIGN-LOCAL-001 | Alignment | SUB: local score ≥ 0; MON: extend matching region → ≥ score; COMP: local for identical = global; INV: distant flank → same local alignment |
| 37 | ☐ | ALIGN-SEMI-001 | Alignment | MON: more matching overlap → higher score; INV: extend non-overlapping part → same core alignment |
| 38 | ☐ | ALIGN-MULTI-001 | Alignment | INV: column permutation doesn't affect per-column scores; MON: add identical seq → score stays or increases |
| 39 | ☐ | PHYLO-DIST-001 | Phylogenetic | SYM: d(a,b)=d(b,a); MON: more mutations → higher distance; COMP: triangle inequality; INV: d(x,x)=0 |
| 40 | ☐ | PHYLO-TREE-001 | Phylogenetic | INV: UPGMA — permuting input order → same topology; MON: closer seqs → shorter branch lengths |
| 41 | ☐ | PHYLO-NEWICK-001 | Phylogenetic | COMP: parse(toNewick(tree))=tree; INV: whitespace doesn't affect parse |
| 42 | ☐ | PHYLO-COMP-001 | Phylogenetic | SYM: RF(a,b)=RF(b,a); COMP: RF(t,t)=0; MON: more rearrangements → higher RF |
| 43 | ☐ | POP-FREQ-001 | PopGen | INV: doubling all counts → same frequencies; COMP: sum(freq)=1.0; INV: reorder samples → same result |
| 44 | ☐ | POP-DIV-001 | PopGen | MON: more diverse sample → higher π; MON: more segregating sites → higher θ; INV: reorder → same |
| 45 | ☐ | POP-HW-001 | PopGen | INV: scaling sample size → same frequencies; COMP: p² + 2pq + q² = 1; MON: larger deviation → larger chi² |
| 46 | ☐ | POP-FST-001 | PopGen | SYM: Fst(A,B)=Fst(B,A); COMP: identical pops → Fst=0; MON: more differentiated → higher Fst |
| 47 | ☐ | POP-LD-001 | PopGen | SYM: LD(a,b)=LD(b,a); COMP: independent loci → D'≈0; MON: perfect linkage → D'=1 |
| 48 | ☐ | CHROM-TELO-001 | Chromosome | MON: more TTAGGG repeats → longer telomere; INV: flanking seq doesn't affect core; SHIFT: flank shifts position |
| 49 | ☐ | CHROM-CENT-001 | Chromosome | INV: non-centromeric flank append → same centromere position; MON: more AT-rich → higher score |
| 50 | ☐ | CHROM-KARYO-001 | Chromosome | COMP: N chromosomes → karyotype with N entries; INV: chromosome order doesn't affect classification |
| 51 | ☐ | CHROM-ANEU-001 | Chromosome | MON: doubled depth → doubled CN estimate; INV: neighbouring region doesn't affect local CN |
| 52 | ☐ | CHROM-SYNT-001 | Chromosome | SYM: synteny(A,B) blocks = synteny(B,A) reversed; INV: non-syntenic insert → same blocks |
| 53 | ☐ | META-CLASS-001 | Metagenomics | INV: duplicate read → same classification; MON: more reference genomes → ≥ classified reads |
| 54 | ☐ | META-PROF-001 | Metagenomics | INV: doubling all reads → same relative abundances; COMP: sum(abundances)=1 |
| 55 | ☐ | META-ALPHA-001 | Metagenomics | MON: remove species → diversity decreases; MON: equalize abundances → max Shannon; COMP: single species → Shannon=0 |
| 56 | ☐ | META-BETA-001 | Metagenomics | SYM: dist(a,b)=dist(b,a); COMP: identical samples → dist=0; MON: remove shared species → higher dist |
| 57 | ☐ | META-BIN-001 | Metagenomics | INV: adding non-overlapping contigs → existing bins unchanged; MON: more contigs → ≥ bins |
| 58 | ☐ | CODON-OPT-001 | Codon | INV: optimized → same protein; MON: more biased table → more codon changes; INV: already optimal → no change |
| 59 | ☐ | CODON-CAI-001 | Codon | MON: replace rare with optimal → CAI increases; COMP: all optimal codons → CAI=1; INV: same seq → same CAI |
| 60 | ☐ | CODON-RARE-001 | Codon | MON: lower usage threshold → more codons flagged; INV: no rare codons → empty result |
| 61 | ☐ | CODON-USAGE-001 | Codon | INV: duplicate sequence → same usage ratios; COMP: sum per AA = 1.0 |
| 62 | ☐ | TRANS-CODON-001 | Translation | COMP: standard code covers all 64 codons; INV: same code table always; INV: geneticCode idempotent |
| 63 | ☐ | TRANS-PROT-001 | Translation | INV: synonymous codon swap → same protein; COMP: stop codon → truncation; INV: frame shift → different protein |
| 64 | ☐ | PARSE-FASTA-001 | FileIO | COMP: write→parse→write = write; INV: adding empty lines → same records; INV: trailing newline irrelevant |
| 65 | ☐ | PARSE-FASTQ-001 | FileIO | COMP: round-trip; INV: quality encoding offset consistent; INV: interleaved order preserved |
| 66 | ☐ | PARSE-BED-001 | FileIO | INV: sorting doesn't change record content; COMP: round-trip; SHIFT: chromStart/End integrity |
| 67 | ☐ | PARSE-VCF-001 | FileIO | INV: comment lines don't affect variant records; COMP: round-trip; INV: INFO field order irrelevant |
| 68 | ☐ | PARSE-GFF-001 | FileIO | COMP: round-trip; INV: attribute order doesn't matter; INV: comment lines don't affect features |
| 69 | ☐ | PARSE-GENBANK-001 | FileIO | COMP: round-trip identity; INV: whitespace in sequence irrelevant |
| 70 | ☐ | PARSE-EMBL-001 | FileIO | COMP: round-trip identity; INV: whitespace in sequence irrelevant |
| 71 | ☐ | RNA-STRUCT-001 | RnaStructure | MON: adding complementary bases → ≥ pairs; INV: non-pairing insert near end → doesn't break existing pairs; COMP: empty → 0 pairs |
| 72 | ☐ | RNA-STEMLOOP-001 | RnaStructure | MON: longer complementary arms → longer stem; INV: loop content doesn't affect stem pairing; COMP: no complement → no stem |
| 73 | ☐ | RNA-ENERGY-001 | RnaStructure | MON: more GC stacks → lower ΔG; COMP: known energies additive; INV: same structure → same ΔG |
| 74 | ☐ | MIRNA-SEED-001 | MiRNA | INV: 3' end changes don't affect seed; COMP: seed match ⊂ full target prediction; INV: seed extraction deterministic |
| 75 | ☐ | MIRNA-TARGET-001 | MiRNA | MON: more seed complementarity → higher score; SUB: stringent → ⊂ lenient results; INV: distant 3'UTR change → same score |
| 76 | ☐ | MIRNA-PRECURSOR-001 | MiRNA | MON: extend stem → more stable precursor; INV: loop sequence doesn't affect structure classification |
| 77 | ☐ | SPLICE-DONOR-001 | Splicing | INV: downstream changes don't affect donor score; MON: consensus GT → higher score; COMP: non-GT → score=0 |
| 78 | ☐ | SPLICE-ACCEPTOR-001 | Splicing | INV: upstream far changes don't affect acceptor score; MON: consensus AG → higher score; COMP: non-AG → score=0 |
| 79 | ☐ | SPLICE-PREDICT-001 | Splicing | COMP: donor + acceptor → exon/intron boundary; MON: more consensus → higher confidence; INV: exonic mutations don't change boundaries |
| 80 | ☐ | DISORDER-PRED-001 | ProteinPred | INV: single residue change → local effect only; MON: more proline/charge → higher disorder; INV: same seq → same scores |
| 81 | ☐ | DISORDER-REGION-001 | ProteinPred | MON: lower threshold → larger/more regions; SUB: strict threshold ⊂ lenient regions; INV: ordered region insert → doesn't affect distant disorder |
| 82 | ☐ | PROTMOTIF-FIND-001 | ProteinMotif | INV: flanking change → same motif detected; MON: broader pattern → ≥ matches; SHIFT: prepend shifts positions |
| 83 | ☐ | PROTMOTIF-PROSITE-001 | ProteinMotif | SUB: specific pattern ⊂ generalized pattern matches; INV: non-matching flank doesn't affect detection |
| 84 | ☐ | PROTMOTIF-DOMAIN-001 | ProteinMotif | INV: domain intact after non-domain insertion; MON: longer domain seq → higher confidence |
| 85 | ☐ | EPIGEN-CPG-001 | Epigenetics | MON: more CG dinucleotides → higher CpG ratio; INV: non-CG flank doesn't change island detection; SHIFT: flank shifts positions |
| 86 | ☐ | ONCO-IMMUNE-001 | Oncology | INV: scaling expression → same relative infiltration; MON: higher marker expression → higher cell fraction; SYM: sample order independent |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 86 |
| ☑ Complete | 7 |
| ☐ Not started | 79 |
| MR relations defined | ~200+ |
