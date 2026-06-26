# Checklist 02: Metamorphic Testing

**Priority:** P0  
**Date:** 2026-03-19  
**Total algorithms:** 258

---

## Description

Metamorphic testing розв'язує «проблему оракула» — коли для довільного входу немає еталонної відповіді. Замість перевірки конкретного результату перевіряються **метаморфні відношення (MR)** — властивості, що пов'язують виходи кількох запусків при трансформації входу. Ідеально для біоінформатичних алгоритмів пошуку, вирівнювання, скорингу та статистики.

**Поточне покриття:** `MetamorphicTests.cs` — 7 юнітів: PAT-IUPAC-001, PAT-PWM-001, PAT-APPROX-002, REP-STR-001, REP-TANDEM-001, REP-INV-001, REP-DIRECT-001 (18+ MR).

**Типи відношень:**
- **SUB** = Subset/Superset (розширення параметрів дає надмножину)
- **MON** = Monotonicity (збільшення X → збільшення/зменшення Y)
- **INV** = Invariance (трансформація входу не змінює вихід)
- **SYM** = Symmetry (f(a,b) = f(b,a))
- **COMP** = Composition (f ⊆ g або f ∘ g = h)
- **SHIFT** = Positional shift (додавання фланку зсуває позиції)

---

## Checklist

### Statuses: ☐ Not Started | ⏳ In Progress | ☑ Complete

| # | Status | Test Unit | Area | MR Relations |
|---|--------|-----------|------|-------------|
| 1 | ☑ | SEQ-GC-001 | Composition | INV: complement preserves GC%; INV: shuffle preserves GC%; INV: case-insensitive |
| 2 | ☑ | SEQ-COMP-001 | Composition | INV: complement(complement(x))=x (involution); INV: length unchanged |
| 3 | ☑ | SEQ-REVCOMP-001 | Composition | INV: revcomp(revcomp(x))=x; INV: length unchanged |
| 4 | ☑ | SEQ-VALID-001 | Composition | INV: case conversion preserves validity; COMP: valid DNA ⊂ valid IUPAC; INV: repeat seq → same result |
| 5 | ☑ | SEQ-COMPLEX-001 | Composition | INV: permutation preserves complexity; MON: homopolymer → min complexity; MON: random → higher |
| 6 | ☑ | SEQ-ENTROPY-001 | Composition | INV: permutation preserves Shannon entropy; MON: uniform → max entropy; MON: single symbol → 0 |
| 7 | ☑ | SEQ-GCSKEW-001 | Composition | SYM: complement flips skew sign; INV: reverse flips cumulative skew; INV: all-G → max positive |
| 8 | ☑ | PAT-EXACT-001 | Matching | SHIFT: prepend flank shifts positions by flank.len; COMP: exact ⊆ hamming(maxDist=0); INV: duplicate → doubled count |
| 9 | ☑ | PAT-APPROX-001 | Matching | SYM: hamming(a,b)=hamming(b,a); MON: higher maxDist → ≥ matches; COMP: exact ⊆ approx(d=0) |
| 10 | ☑ | PAT-APPROX-002 | Matching | MON: higher maxEdits → superset; COMP: exact ⊆ approximate; SYM: editDist(a,b)=editDist(b,a); INV: prefix extension; R: non-negativity |
| 11 | ☑ | PAT-IUPAC-001 | Matching | SUB: degeneracy hierarchy (N ⊇ all); INV: shift-invariance; MON: more degenerate → ≥ matches |
| 12 | ☑ | PAT-PWM-001 | Matching | INV: permuted training set → same PWM; INV: duplicated set → same consensus; SUB: lower threshold → ≥ matches |
| 13 | ☑ | REP-STR-001 | Repeats | MON: lower minRepeats → superset; MON: wider unit range → superset; MON: doubled repeat → increased count |
| 14 | ☑ | REP-TANDEM-001 | Repeats | SHIFT: prepend flank shifts positions; MON: higher minReps → subset; MON: higher minUnit → subset |
| 15 | ☑ | REP-INV-001 | Repeats | SYM: right arm = revcomp(left arm); INV: non-overlapping flank doesn't change core |
| 16 | ☑ | REP-DIRECT-001 | Repeats | INV: duplication → consistent results; MON: lower minLen → ≥ results; SHIFT: flank shifts positions |
| 17 | ☑ | REP-PALIN-001 | Repeats | INV: palindrome = revcomp of self; MON: wider len range → ≥ palindromes; SHIFT: flank shifts positions |
| 18 | ☑ | CRISPR-PAM-001 | MolTools | MON: longer sequence → ≥ PAM sites; INV: non-PAM region append → same count; SHIFT: flank shifts positions |
| 19 | ☑ | CRISPR-GUIDE-001 | MolTools | SUB: guide from PAM site ⊂ valid guides; MON: stricter scoring → ≤ guides; INV: downstream change → same guide |
| 20 | ☑ | CRISPR-OFF-001 | MolTools | MON: more mismatches → lower off-target score; MON: seed region mismatch penalized more; COMP: 0 mismatches → max score |
| 21 | ☑ | PRIMER-TM-001 | MolTools | MON: add GC → Tm increases; MON: add AT → Tm decreases; INV: same sequence → same Tm |
| 22 | ☑ | PRIMER-DESIGN-001 | MolTools | MON: wider Tm range → ≥ primers; SUB: stricter GC% → ⊆ results; INV: longer template → ≥ candidates |
| 23 | ☑ | PRIMER-STRUCT-001 | MolTools | MON: more self-complementary → higher dimer score; INV: non-complementary extension → same hairpin |
| 24 | ☑ | PROBE-DESIGN-001 | MolTools | MON: wider Tm → ≥ probes; SUB: stricter uniqueness → ⊆ results; INV: unrelated region append → same probes |
| 25 | ☑ | PROBE-VALID-001 | MolTools | MON: lower specificity threshold → more pass; INV: same input → same result |
| 26 | ☑ | RESTR-FIND-001 | MolTools | SHIFT: prepend flank shifts positions; MON: more enzymes → ≥ total sites; INV: non-site append → same sites |
| 27 | ☑ | RESTR-DIGEST-001 | MolTools | COMP: 0 sites → 1 fragment = full seq; MON: more enzymes → ≥ fragments; INV: fragment sum = seq length |
| 28 | ☑ | ANNOT-ORF-001 | Annotation | MON: lower minLen → ≥ ORFs; SHIFT: prepend shifts positions; INV: non-coding insert doesn't change upstream |
| 29 | ☑ | ANNOT-GENE-001 | Annotation | COMP: gene ⊃ ORF; INV: non-coding insertion doesn't affect upstream; MON: longer seq → ≥ genes |
| 30 | ☑ | ANNOT-PROM-001 | Annotation | MON: lower score threshold → ≥ promoters; SHIFT: prepend shifts; INV: downstream change → same promoter |
| 31 | ☑ | ANNOT-GFF-001 | Annotation | INV: re-serialization preserves line count; COMP: parse(write(x))=x; INV: attribute order irrelevant |
| 32 | ☑ | KMER-COUNT-001 | K-mer | INV: total k-mer instances = seqLen - k + 1; MON: k+1 → ≤ distinct k-mers; INV: reverse → same set (canonical) |
| 33 | ☑ | KMER-FREQ-001 | K-mer | INV: duplicate seq → same freqs; SYM: complement has related frequency profile; INV: sum(freq)=1 |
| 34 | ☑ | KMER-FIND-001 | K-mer | MON: lower minFreq → ≥ k-mers; SUB: top-1 ⊆ top-5; INV: repeat seq → high frequency for unit k-mer |
| 35 | ☑ | ALIGN-GLOBAL-001 | Alignment | SYM: score(a,b)=score(b,a); MON: more matches → higher score; COMP: identity → max score; INV: gap-only insert → score change = gap penalty |
| 36 | ☑ | ALIGN-LOCAL-001 | Alignment | SUB: local score ≥ 0; MON: extend matching region → ≥ score; COMP: local for identical = global; INV: distant flank → same local alignment |
| 37 | ☑ | ALIGN-SEMI-001 | Alignment | MON: more matching overlap → higher score; INV: extend non-overlapping part → same core alignment |
| 38 | ☑ | ALIGN-MULTI-001 | Alignment | INV: column permutation doesn't affect per-column scores; MON: add identical seq → score stays or increases |
| 39 | ☑ | PHYLO-DIST-001 | Phylogenetic | SYM: d(a,b)=d(b,a); MON: more mutations → higher distance; COMP: triangle inequality; INV: d(x,x)=0 |
| 40 | ☑ | PHYLO-TREE-001 | Phylogenetic | INV: UPGMA — permuting input order → same topology; MON: closer seqs → shorter branch lengths |
| 41 | ☑ | PHYLO-NEWICK-001 | Phylogenetic | COMP: parse(toNewick(tree))=tree; INV: whitespace doesn't affect parse |
| 42 | ☑ | PHYLO-COMP-001 | Phylogenetic | SYM: RF(a,b)=RF(b,a); COMP: RF(t,t)=0; MON: more rearrangements → higher RF |
| 43 | ☑ | POP-FREQ-001 | PopGen | INV: doubling all counts → same frequencies; COMP: sum(freq)=1.0; INV: reorder samples → same result |
| 44 | ☑ | POP-DIV-001 | PopGen | MON: more diverse sample → higher π; MON: more segregating sites → higher θ; INV: reorder → same |
| 45 | ☑ | POP-HW-001 | PopGen | INV: scaling sample size → same frequencies; COMP: p² + 2pq + q² = 1; MON: larger deviation → larger chi² |
| 46 | ☑ | POP-FST-001 | PopGen | SYM: Fst(A,B)=Fst(B,A); COMP: identical pops → Fst=0; MON: more differentiated → higher Fst |
| 47 | ☑ | POP-LD-001 | PopGen | SYM: LD(a,b)=LD(b,a); COMP: independent loci → D'≈0; MON: perfect linkage → D'=1 |
| 48 | ☑ | CHROM-TELO-001 | Chromosome | MON: more TTAGGG repeats → longer telomere; INV: flanking seq doesn't affect core; SHIFT: flank shifts position |
| 49 | ☑ | CHROM-CENT-001 | Chromosome | INV: non-centromeric flank append → same centromere position; MON: more AT-rich → higher score |
| 50 | ☑ | CHROM-KARYO-001 | Chromosome | COMP: N chromosomes → karyotype with N entries; INV: chromosome order doesn't affect classification |
| 51 | ☑ | CHROM-ANEU-001 | Chromosome | MON: doubled depth → doubled CN estimate; INV: neighbouring region doesn't affect local CN |
| 52 | ☑ | CHROM-SYNT-001 | Chromosome | SYM: synteny(A,B) blocks = synteny(B,A) reversed; INV: non-syntenic insert → same blocks |
| 53 | ☑ | META-CLASS-001 | Metagenomics | INV: duplicate read → same classification; MON: more reference genomes → ≥ classified reads |
| 54 | ☑ | META-PROF-001 | Metagenomics | INV: doubling all reads → same relative abundances; COMP: sum(abundances)=1 |
| 55 | ☑ | META-ALPHA-001 | Metagenomics | MON: remove species → diversity decreases; MON: equalize abundances → max Shannon; COMP: single species → Shannon=0 |
| 56 | ☑ | META-BETA-001 | Metagenomics | SYM: dist(a,b)=dist(b,a); COMP: identical samples → dist=0; MON: remove shared species → higher dist |
| 57 | ☑ | META-BIN-001 | Metagenomics | INV: adding non-overlapping contigs → existing bins unchanged; MON: more contigs → ≥ bins |
| 58 | ☑ | CODON-OPT-001 | Codon | INV: optimized → same protein; MON: more biased table → more codon changes; INV: already optimal → no change |
| 59 | ☑ | CODON-CAI-001 | Codon | MON: replace rare with optimal → CAI increases; COMP: all optimal codons → CAI=1; INV: same seq → same CAI |
| 60 | ☑ | CODON-RARE-001 | Codon | MON: higher threshold → superset of flagged codons (f < τ monotone in τ); INV: no rare codons → empty result |
| 61 | ☑ | CODON-USAGE-001 | Codon | INV: duplicate sequence → same usage ratios; COMP: sum per AA = 1.0 |
| 62 | ☑ | TRANS-CODON-001 | Translation | COMP: standard code covers all 64 codons; INV: same code table always; INV: geneticCode idempotent |
| 63 | ☑ | TRANS-PROT-001 | Translation | INV: synonymous codon swap → same protein; COMP: stop codon → truncation; INV: frame shift → different protein |
| 64 | ☑ | PARSE-FASTA-001 | FileIO | COMP: write→parse→write = write; INV: adding empty lines → same records; INV: trailing newline irrelevant |
| 65 | ☑ | PARSE-FASTQ-001 | FileIO | COMP: round-trip; INV: quality encoding offset consistent; INV: interleaved order preserved |
| 66 | ☑ | PARSE-BED-001 | FileIO | INV: sorting doesn't change record content; COMP: round-trip; SHIFT: chromStart/End integrity |
| 67 | ☑ | PARSE-VCF-001 | FileIO | INV: comment lines don't affect variant records; COMP: round-trip; INV: INFO field order irrelevant |
| 68 | ☑ | PARSE-GFF-001 | FileIO | COMP: round-trip; INV: attribute order doesn't matter; INV: comment lines don't affect features |
| 69 | ☑ | PARSE-GENBANK-001 | FileIO | COMP: round-trip identity; INV: whitespace in sequence irrelevant |
| 70 | ☑ | PARSE-EMBL-001 | FileIO | COMP: round-trip identity; INV: whitespace in sequence irrelevant |
| 71 | ☑ | RNA-STRUCT-001 | RnaStructure | MON: adding complementary bases → ≥ pairs; INV: non-pairing insert near end → doesn't break existing pairs; COMP: empty → 0 pairs |
| 72 | ☑ | RNA-STEMLOOP-001 | RnaStructure | MON: longer complementary arms → longer stem; INV: loop content doesn't affect stem pairing; COMP: no complement → no stem |
| 73 | ☑ | RNA-ENERGY-001 | RnaStructure | MON: more GC stacks → lower ΔG; COMP: known energies additive; INV: same structure → same ΔG |
| 74 | ☑ | MIRNA-SEED-001 | MiRNA | INV: 3' end changes don't affect seed; COMP: seed match ⊂ full target prediction; INV: seed extraction deterministic |
| 75 | ☑ | MIRNA-TARGET-001 | MiRNA | MON: more seed complementarity → higher score; SUB: stringent → ⊂ lenient results; INV: distant 3'UTR change → same score |
| 76 | ☑ | MIRNA-PRECURSOR-001 | MiRNA | MON: extend stem → more stable precursor; INV: loop sequence doesn't affect structure classification |
| 77 | ☑ | SPLICE-DONOR-001 | Splicing | INV: downstream changes don't affect donor score; MON: consensus GT → higher score; COMP: non-GT → score=0 |
| 78 | ☑ | SPLICE-ACCEPTOR-001 | Splicing | INV: upstream far changes don't affect acceptor score; MON: consensus AG → higher score; COMP: non-AG → score=0 |
| 79 | ☑ | SPLICE-PREDICT-001 | Splicing | COMP: donor + acceptor → exon/intron boundary; MON: more consensus → higher confidence; INV: exonic mutations don't change boundaries |
| 80 | ☑ | DISORDER-PRED-001 | ProteinPred | INV: single residue change → local effect only; MON: more proline/charge → higher disorder; INV: same seq → same scores |
| 81 | ☑ | DISORDER-REGION-001 | ProteinPred | MON: lower threshold → larger/more regions; SUB: strict threshold ⊂ lenient regions; INV: ordered region insert → doesn't affect distant disorder |
| 82 | ☑ | PROTMOTIF-FIND-001 | ProteinMotif | INV: flanking change → same motif detected; MON: broader pattern → ≥ matches; SHIFT: prepend shifts positions |
| 83 | ☑ | PROTMOTIF-PROSITE-001 | ProteinMotif | SUB: specific pattern ⊂ generalized pattern matches; INV: non-matching flank doesn't affect detection |
| 84 | ☑ | PROTMOTIF-DOMAIN-001 | ProteinMotif | INV: domain intact after non-domain insertion; MON: more conserved signature → higher confidence (information content; not matched-substring length) |
| 85 | ☑ | EPIGEN-CPG-001 | Epigenetics | MON: more CG dinucleotides → higher CpG ratio; INV: non-CG flank doesn't change core-island detection; SHIFT: flank shifts positions |
| 86 | ☑ | ONCO-IMMUNE-001 | Oncology | INV: scaling expression → same relative infiltration; MON: higher marker expression → higher cell fraction; SYM: sample order independent |
| 87 | ☑ | ONCO-SOMATIC-001 | Oncology | MON: more tumor alt evidence → ≥ score / superset of calls; INV: adding pure-reference reads adds no somatic calls; SYM: read order independent |
| 88 | ☑ | ONCO-VAF-001 | Oncology | INV: scaling ref+alt depth equally → same VAF; MON: +k alt reads → higher VAF; INV: read order independent |
| 89 | ☑ | ONCO-DRIVER-001 | Oncology | MON: more samples sharing a mutation → ≥ driver score; INV: relabeling passenger genes preserves driver set |
| 90 | ☑ | ONCO-ARTIFACT-001 | Oncology | MON: stricter artifact evidence (higher GIV) → subset of survivors; INV: duplicating a passing variant keeps it passing |
| 91 | ☑ | ONCO-ANNOT-001 | Oncology | INV: uniform identity shift carries annotations equally (no coordinate field; identity is the analog); INV: variant order independent |
| 92 | ☑ | ONCO-TMB-001 | Oncology | INV: doubling panel-Mb and mutations → same TMB density; MON: +1 coding mutation → ≥ TMB; INV: order independent |
| 93 | ☑ | ONCO-MSI-001 | Oncology | MON: more unstable loci → ≥ MSI score; INV: locus order independent |
| 94 | ☑ | ONCO-HRD-001 | Oncology | MON: adding an LOH/TAI/LST event → ≥ HRD; INV: event order independent |
| 95 | ☑ | ONCO-LOH-001 | Oncology | INV: allele symmetry (retained-allele value irrelevant) + segment order; MON: looser size criterion / more regions → superset (CN-based, not BAF) |
| 96 | ☑ | ONCO-SIG-001 | Oncology | INV: reverse-complementing a variant maps to the same pyrimidine channel; INV: variant order independent |
| 97 | ☑ | ONCO-SIG-002 | Oncology | INV: scaling the catalogue by k scales exposures by k; MON: adding signature-consistent mutations → ≥ that exposure |
| 98 | ☑ | ONCO-SIG-003 | Oncology | INV: same seed → identical CI; INV: point estimate setting-independent; MON: wider confidence → non-narrower CI (percentile-bootstrap width converges, not monotone in reps) |
| 99 | ☑ | ONCO-SIG-004 | Oncology | INV: scaling all exposures preserves the dominant process; INV: signature order independent |
| 100 | ☑ | ONCO-FUSION-001 | Oncology | MON: more split reads → ≥ confidence (TotalSupport); INV: candidate order preserves fusion count (no breakpoint coordinate in API) |
| 101 | ☑ | ONCO-FUSION-002 | Oncology | SUB: matched ⊆ known DB; INV: 5'/3' orientation preserved + case-insensitive (no coordinate in API) |
| 102 | ☑ | ONCO-FUSION-003 | Oncology | INV: codon-multiple coordinate shift preserves in/out-of-frame classification (non-codon shift flips it) |
| 103 | ☑ | ONCO-CNA-001 | Oncology | MON: higher log2-ratio → ≥ CN class; INV: segment order independent |
| 104 | ☑ | ONCO-CNA-002 | Oncology | MON: higher CN keeps focal amplification; INV: prepend flank shifts focal coordinates |
| 105 | ☑ | ONCO-CNA-003 | Oncology | MON: lower CN keeps homozygous deletion; INV: segment order independent |
| 106 | ☑ | ONCO-PURITY-001 | Oncology | MON: scaling clonal VAFs up → ≥ purity; INV: variant order independent |
| 107 | ☑ | ONCO-PLOIDY-001 | Oncology | MON: amplifying more segments → ≥ ploidy; INV: segment order independent |
| 108 | ☑ | ONCO-CLONAL-001 | Oncology | MON: higher CCF keeps a clonal call clonal; INV: variant order independent |
| 109 | ☑ | ONCO-NEO-001 | Oncology | INV: flanking-context shift preserves the peptide set tiling the mutation |
| 110 | ☑ | ONCO-MHC-001 | Oncology | MON: lower IC50 → stronger-or-equal binding class; INV: peptide order independent |
| 111 | ☑ | ONCO-CTDNA-001 | Oncology | MON: spiking tumor signal → ≥ detection probability; INV: depends only on aggregate λ (read order/aggregation independent) |
| 112 | ☑ | ONCO-MRD-001 | Oncology | MON: observing more tracked variants keeps MRD positive; INV: variant order independent |
| 113 | ☑ | ONCO-CHIP-001 | Oncology | SUB: survivors ⊆ input; INV: duplicating a CHIP variant keeps it flagged |
| 114 | ☑ | ONCO-PHYLO-001 | Oncology | INV: sample-column permutation preserves topology; SYM: pairwise clone distance symmetric |
| 115 | ☑ | ONCO-CCF-001 | Oncology | MON: higher VAF → ≥ CCF at fixed CN/purity; INV: variant order independent |
| 116 | ☑ | ONCO-HETERO-001 | Oncology | INV: scaling all VAFs equally preserves MATH; MON: wider VAF spread → ≥ heterogeneity |
| 117 | ☑ | ONCO-HLA-001 | Oncology | INV: normalisation idempotent; INV: allele-string normalisation stable (case/whitespace) — API parses an allele string, not reads |
| 118 | ☑ | ONCO-ACTION-001 | Oncology | MON: stronger evidence → ≥ tier; INV: variant order independent |
| 119 | ☑ | ONCO-SV-001 | Oncology | INV: coordinate shift preserves rearrangement class; MON: more clustered breakpoints → chromothripsis |
| 120 | ☑ | ONCO-EXPR-001 | Oncology | INV: scaling all expression equally preserves z-scores/outliers; MON: lower threshold → superset |
| 121 | ☑ | SEQ-COMPOSITION-001 | Statistics | INV: permutation invariant; P: complement swaps A↔T and C↔G counts |
| 122 | ☑ | SEQ-DINUC-001 | Statistics | INV: reverse-complement maps each dinucleotide to its revcomp; SHIFT: prepend flank adds only boundary dinucleotides |
| 123 | ☑ | SEQ-HYDRO-001 | Statistics | INV: permutation changes profile but not mean; MON: adding a hydrophobic residue → ≥ mean |
| 124 | ☑ | SEQ-MW-001 | Statistics | ADD: MW(a+b) = MW(a)+MW(b) − water; INV: permutation invariant |
| 125 | ☑ | SEQ-PI-001 | Statistics | INV: permutation invariant; MON: more acidic residues → lower pI |
| 126 | ☑ | SEQ-SECSTRUCT-001 | Statistics | SHIFT: prepend flank shifts assignments; INV: deterministic |
| 127 | ☑ | SEQ-STATS-001 | Statistics | INV: permutation invariant; P: concatenation sums counts |
| 128 | ☑ | SEQ-SUMMARY-001 | Statistics | INV: permutation invariant for composition fields; SHIFT: length additive on concatenation |
| 129 | ☑ | SEQ-THERMO-001 | Statistics | MON: more GC pairs → lower ΔG; INV: permutation changes nearest-neighbour context only |
| 130 | ☑ | SEQ-TM-001 | Statistics | MON: more GC → higher Tm; INV: case-insensitive |
| 131 | ☑ | COMPGEN-ANI-001 | Comparative | SYM: ANI(A,B)=ANI(B,A); INV: ANI(A,A)=100; MON: more mutations → lower ANI |
| 132 | ☑ | COMPGEN-CLUSTER-001 | Comparative | MON: lower identity threshold → superset; INV: genome order independent |
| 133 | ☑ | COMPGEN-COMPARE-001 | Comparative | SYM: order independent; MON: more shared genes → higher similarity |
| 134 | ☑ | COMPGEN-DOTPLOT-001 | Comparative | INV: revcomp maps diagonal → anti-diagonal; SHIFT: prepend flank shifts dots |
| 135 | ☑ | COMPGEN-ORTHO-001 | Comparative | SYM: ortholog relation symmetric; INV: genome order independent |
| 136 | ☑ | COMPGEN-RBH-001 | Comparative | SYM: RBH symmetric; INV: input order independent |
| 137 | ☑ | COMPGEN-REARR-001 | Comparative | INV: identity → no rearrangements; SYM: (A,B) consistent with (B,A) |
| 138 | ☑ | COMPGEN-REVERSAL-001 | Comparative | SYM: symmetric; INV: identical permutation → 0; MON: more reversals applied → ≥ distance |
| 139 | ☑ | COMPGEN-SYNTENY-001 | Comparative | MON: lower minBlockSize → superset; INV: revcomp preserves block count |
| 140 | ☑ | ASSEMBLY-CONSENSUS-001 | Assembly | INV: read order independent; MON: adding a concordant read preserves consensus |
| 141 | ☑ | ASSEMBLY-CORRECT-001 | Assembly | INV: error-free reads unchanged; MON: more coverage → ≤ residual errors |
| 142 | ☑ | ASSEMBLY-COVER-001 | Assembly | INV: read order independent; ADD: coverage additive over reads |
| 143 | ☑ | ASSEMBLY-DBG-001 | Assembly | INV: read order independent; MON: larger k → ≤ spurious joins |
| 144 | ☑ | ASSEMBLY-MERGE-001 | Assembly | INV: merge order independent for compatible contigs |
| 145 | ☑ | ASSEMBLY-OLC-001 | Assembly | INV: read order independent; MON: higher minOverlap → ≤ joins |
| 146 | ☑ | ASSEMBLY-SCAFFOLD-001 | Assembly | INV: link order independent |
| 147 | ☑ | ASSEMBLY-STATS-001 | Assembly | INV: contig order independent; MON: splitting a contig → ≤ N50 |
| 148 | ☑ | ASSEMBLY-TRIM-001 | Assembly | MON: higher cutoff → subset of bases; INV: read order independent |
| 149 | ☑ | RNA-DOTBRACKET-001 | RnaStructure | RT: parse∘format identity; INV: pairing preserved under reparse |
| 150 | ☑ | RNA-HAIRPIN-001 | RnaStructure | MON: larger loop → higher (less stable) energy; INV: closing-pair context preserved |
| 151 | ☑ | RNA-INVERT-001 | RnaStructure | SYM: arms reverse-complementary; INV: revcomp preserves count |
| 152 | ☑ | RNA-MFE-001 | RnaStructure | MON: more GC pairs → lower MFE; INV: U/T case-insensitive |
| 153 | ☑ | RNA-PAIR-001 | RnaStructure | SYM: canPair(a,b)=canPair(b,a); INV: case-insensitive |
| 154 | ☑ | RNA-PARTITION-001 | RnaStructure | MON: more pairing options → higher Z; INV: deterministic |
| 155 | ☑ | RNA-PSEUDOKNOT-001 | RnaStructure | INV: nested structure → no pseudoknot; SHIFT: prepend flank shifts positions |
| 156 | ☑ | KMER-ASYNC-001 | K-mer | INV: async = sync; INV: read order independent |
| 157 | ☑ | KMER-BOTH-001 | K-mer | SYM: reverse-complement invariance; ADD: counts additive on concatenation |
| 158 | ☑ | KMER-DIST-001 | K-mer | SYM: d(a,b)=d(b,a); INV: d(x,x)=0 |
| 159 | ☑ | KMER-GENERATE-001 | K-mer | INV: order independent; P: set closed under all k-mers |
| 160 | ☑ | KMER-POSITIONS-001 | K-mer | SHIFT: prepend flank shifts positions; INV: order independent |
| 161 | ☑ | KMER-STATS-001 | K-mer | INV: permutation changes positions not counts; ADD: counts additive on concatenation |
| 162 | ☑ | KMER-UNIQUE-001 | K-mer | MON: duplicating a k-mer removes it from unique set; INV: order independent |
| 163 | ☑ | PROTMOTIF-CC-001 | ProteinMotif | INV: deterministic; SHIFT: prepend flank shifts positions |
| 164 | ☑ | PROTMOTIF-COMMON-001 | ProteinMotif | MON: more sequences sharing → ≥ support; INV: input order independent |
| 165 | ☑ | PROTMOTIF-LC-001 | ProteinMotif | MON: lower threshold → superset; SHIFT: prepend flank shifts regions |
| 166 | ☑ | PROTMOTIF-PATTERN-001 | ProteinMotif | SHIFT: prepend flank shifts matches; SUB: broader pattern → ≥ matches |
| 167 | ☑ | PROTMOTIF-SP-001 | ProteinMotif | INV: C-terminal extension doesn't change N-terminal signal |
| 168 | ☑ | PROTMOTIF-TM-001 | ProteinMotif | MON: lower threshold → superset; SHIFT: prepend flank shifts helices |
| 169 | ☑ | MOTIF-CONS-001 | Matching | INV: row order independent; INV: duplicating a row preserves consensus |
| 170 | ☑ | MOTIF-DISCOVER-001 | Matching | MON: lower support → superset; SHIFT: prepend flank shifts positions |
| 171 | ☑ | MOTIF-GENERATE-001 | Matching | INV: row order independent |
| 172 | ☑ | MOTIF-REGULATORY-001 | Matching | SHIFT: prepend flank shifts positions; SUB: broader set → ≥ matches |
| 173 | ☑ | MOTIF-SHARED-001 | Matching | INV: input order independent; SUB: fewer inputs → ⊇ shared set |
| 174 | ☑ | PAT-APPROX-003 | Matching | INV: exact match → 0; MON: best ≤ any candidate distance |
| 175 | ☑ | GENOMIC-COMMON-001 | Analysis | INV: input order independent; SUB: more inputs → ⊆ common |
| 176 | ☑ | GENOMIC-MOTIFS-001 | Analysis | SHIFT: prepend flank shifts positions; INV: deterministic |
| 177 | ☑ | GENOMIC-ORF-001 | Analysis | SHIFT: prepend in-frame flank shifts ORFs; INV: revcomp gives reverse-strand ORFs |
| 178 | ☑ | GENOMIC-REPEAT-001 | Analysis | MON: lower minLen → superset; SHIFT: prepend flank shifts positions |
| 179 | ☑ | GENOMIC-SIMILARITY-001 | Analysis | SYM: sim(a,b)=sim(b,a); INV: sim(x,x)=1 |
| 180 | ☑ | GENOMIC-TANDEM-001 | Analysis | MON: lower minReps → superset; SHIFT: prepend flank shifts positions |
| 181 | ☑ | EPIGEN-AGE-001 | Epigenetics | MON: more clock-site methylation → higher age; INV: site order independent |
| 182 | ☑ | EPIGEN-BISULF-001 | Epigenetics | INV: methylated-C set preserved; SHIFT: prepend flank shifts conversions |
| 183 | ☑ | EPIGEN-CHROM-001 | Epigenetics | INV: region order independent; SHIFT: prepend flank shifts states |
| 184 | ☑ | EPIGEN-DMR-001 | Epigenetics | MON: lower threshold → superset; SYM: DMR(A,B) consistent with (B,A) |
| 185 | ☑ | EPIGEN-METHYL-001 | Epigenetics | INV: read order independent; ADD: counts additive over reads |
| 186 | ☑ | VARIANT-ANNOT-001 | Variants | SHIFT: coordinate shift shifts annotations; INV: variant order independent |
| 187 | ☑ | VARIANT-CALL-001 | Variants | MON: deeper coverage → superset of confident calls; INV: read order independent |
| 188 | ☑ | VARIANT-INDEL-001 | Variants | SHIFT: prepend flank shifts indel positions; INV: read order independent |
| 189 | ☑ | VARIANT-SNP-001 | Variants | SHIFT: prepend flank shifts SNP positions; INV: read order independent |
| 190 | ☑ | PANGEN-CLUSTER-001 | PanGenome | MON: lower identity → coarser clusters; INV: gene order independent |
| 191 | ☑ | PANGEN-CORE-001 | PanGenome | MON: more genomes → ⊆ core; INV: genome order independent |
| 192 | ☑ | PANGEN-HEAP-001 | PanGenome | INV: genome order independent; MON: more genomes → better fit |
| 193 | ☑ | PANGEN-MARKER-001 | PanGenome | SUB: markers ⊆ core; INV: genome order independent |
| 194 | ☑ | META-FUNC-001 | Metagenomics | INV: read order independent; SUB: larger DB → ≥ assignments |
| 195 | ☑ | META-PATHWAY-001 | Metagenomics | MON: more pathway genes → higher enrichment; INV: gene order independent |
| 196 | ☑ | META-RESIST-001 | Metagenomics | INV: read order independent; SUB: larger DB → ≥ hits |
| 197 | ☑ | META-TAXA-001 | Metagenomics | INV: sample order independent; MON: larger effect → lower p-value |
| 198 | ☑ | TRANS-DIFF-001 | Transcriptome | SYM: FC(A,B) = −FC(B,A); INV: gene order independent |
| 199 | ☑ | TRANS-EXPR-001 | Transcriptome | INV: read order independent; HOMO: scaling depth preserves TPM |
| 200 | ☑ | TRANS-SPLICE-001 | Transcriptome | INV: read order independent; SHIFT: prepend flank shifts exon coords |
| 201 | ☑ | SV-BREAKPOINT-001 | StructuralVar | SHIFT: prepend flank shifts breakpoints; MON: more split reads → ≥ confidence |
| 202 | ☑ | SV-CNV-001 | StructuralVar | MON: higher coverage ratio → higher CN; INV: bin order independent |
| 203 | ☑ | SV-DETECT-001 | StructuralVar | INV: identical genomes → no SV; SHIFT: coordinate shift shifts SVs |
| 204 | ☑ | DISORDER-LC-001 | ProteinPred | MON: lower threshold → superset; SHIFT: prepend flank shifts regions |
| 205 | ☑ | DISORDER-MORF-001 | ProteinPred | INV: deterministic; SHIFT: prepend flank shifts MoRFs |
| 206 | ☑ | DISORDER-PROPENSITY-001 | ProteinPred | SHIFT: prepend flank shifts profile; INV: deterministic |
| 207 | ☑ | POP-ANCESTRY-001 | PopGen | INV: individual order independent; P: proportions sum to 1 |
| 208 | ☑ | POP-ROH-001 | PopGen | MON: lower minLen → superset; SHIFT: prepend flank shifts ROH |
| 209 | ☑ | POP-SELECT-001 | PopGen | INV: locus order independent; MON: stronger selection → higher signal |
| 210 | ☑ | SEQ-ATSKEW-001 | Composition | SYM: complement reverses sign; INV: cumulative length = seq length |
| 211 | ☑ | SEQ-REPLICATION-001 | Composition | INV: rotation shifts predicted origin; SYM: complement reflects origin |
| 212 | ☑ | SEQ-RNACOMP-001 | Composition | INV: complement∘complement = identity; P: A↔U, G↔C |
| 213 | ☑ | CODON-ENC-001 | Codon | MON: more biased usage → lower ENC; INV: codon order independent |
| 214 | ☑ | CODON-RSCU-001 | Codon | INV: codon order independent; P: per-AA RSCU mean = 1 |
| 215 | ☑ | CODON-STATS-001 | Codon | INV: order independent; ADD: counts additive on concatenation |
| 216 | ☑ | ANNOT-CODING-001 | Annotation | INV: deterministic; MON: real ORF → higher score |
| 217 | ☑ | ANNOT-CODONUSAGE-001 | Annotation | INV: codon order independent; P: per-AA sum = 1 |
| 218 | ☑ | ANNOT-REPEAT-001 | Annotation | MON: lower minLen → superset; SHIFT: prepend flank shifts elements |
| 219 | ☑ | QUALITY-PHRED-001 | Quality | RT: encode∘decode identity; INV: offset consistency |
| 220 | ☑ | QUALITY-STATS-001 | Quality | INV: order independent for mean; ADD: counts additive |
| 221 | ☑ | PHYLO-BOOT-001 | Phylogenetic | INV: same seed → same support; SYM: distance symmetric |
| 222 | ☑ | PHYLO-STATS-001 | Phylogenetic | INV: leaf relabeling preserves stats; INV: deterministic |
| 223 | ☑ | TRANS-SIXFRAME-001 | Translation | INV: frames 4–6 = translation of revcomp; P: exactly 6 frames |
| 224 | ☑ | RESTR-FILTER-001 | MolTools | SUB: filtered ⊆ all; MON: stricter criteria → subset |
| 225 | ☑ | MIRNA-PAIR-001 | MiRNA | SHIFT: prepend flank shifts alignment; INV: deterministic |
| 226 | ☑ | ALIGN-STATS-001 | Alignment | SYM: stats(a,b)=stats(b,a); P: identity(x,x)=1 |
| 227 | ☑ | SEQ-CODON-FREQ-001 | Statistics | INV: codon-preserving shuffle keeps frequencies; SCALE: triplicating seq preserves freqs |
| 228 | ☑ | SEQ-COMPLEX-COMPRESS-001 | Complexity | INV: case change preserves ratio; ORDER: concatenating repeats lowers ratio |
| 229 | ☑ | SEQ-COMPLEX-DUST-001 | Complexity | INV: complement preserves DUST; MONO: adding homopolymer run raises score |
| 230 | ☑ | SEQ-COMPLEX-KMER-001 | Complexity | INV: reverse preserves k-mer entropy; MONO: more distinct k-mers → higher entropy |
| 231 | ☑ | SEQ-COMPLEX-WINDOW-001 | Complexity | INV: complement preserves per-window score; SHIFT: prepend flank shifts profile |
| 232 | ☑ | SEQ-ENTROPY-PROFILE-001 | Statistics | INV: complement preserves profile; SHIFT: prepend flank shifts profile |
| 233 | ☑ | SEQ-GC-ANALYSIS-001 | Composition | INV: complement preserves GC%; INV: shuffle preserves GC% |
| 234 | ☑ | SEQ-GC-PROFILE-001 | Statistics | INV: complement preserves GC profile; SHIFT: prepend flank shifts profile |
| 235 | ☑ | ONCO-ASCAT-001 | Oncology | INV: constant logR shift preserves breakpoints; INV: A/B allele swap preserves total CN |
| 236 | ☑ | RNA-PKPREDICT-001 | Analysis | INV: known H-type knot recovered; INV: no spurious knot on a plain hairpin |
| 237 | ☑ | RNA-PKRECURSIVE-001 | Analysis | MON: recursive ΔG ≤ single-knot ΔG; INV: separable knots all recovered |
| 238 | ☑ | RNA-ACCESS-001 | RnaStructure | MON: extending the queried region cannot raise its unpaired probability; INV: sequence-independent constants reproduce analytic GAAAC value |
| 239 | ☑ | PROTMOTIF-HMM-001 | ProteinMotif | MON: appending random flank does not raise a true domain's bit score; SUB: stricter E-value → ⊆ hit set |
| 240 | ☑ | PRIMER-NNTM-001 | MolTools | MON: raising monovalent salt raises Tm; INV: reverse-complement has equal duplex Tm |
| 241 | ☑ | PRIMER-HAIRPIN-001 | MolTools | MON: lengthening a complementary stem lowers ΔG; INV: no stem possible → no hairpin |
| 242 | ☑ | PRIMER-DIMER-001 | MolTools | INV: self-dimer of S equals hetero-dimer(S,S); MON: extending WC alignment lowers ΔG |
| 243 | ☑ | PROBE-LNATM-001 | MolTools | MON: adding an LNA base → Tm ≥ unmodified Tm; INV: all-DNA input reduces to standard NN Tm |
| 244 | ☑ | PROBE-EVALUE-001 | MolTools | MON: increasing database size raises E for a fixed raw score; INV: λ solves Σ p_i p_j e^{λ s_ij}=1 |
| 245 | ☑ | MHC-NN-001 | Oncology | INV: BLOSUM-encoded peptide reproduces oracle within 0.03%; SUB: shorter peptide padded centred |
| 246 | ☑ | MHC-MATRIX-001 | Oncology | INV: IC50 = 50000^(1-score); MON: improving an anchor residue lowers IC50 |
| 247 | ☑ | IMMUNE-NUSVR-001 | Oncology | INV: mixing known fractions of pure profiles recovers those fractions; SUB: ν controls support-vector count |
| 248 | ☑ | META-CHECKM-001 | Metagenomics | MON: removing a marker lowers completeness; MON: duplicating a marker raises contamination |
| 249 | ☐ | META-TETRA-001 | Metagenomics | INV: reverse-complement-merged counts give identical z-vector; MON: identical sequences → correlation 1 |
| 250 | ☐ | SPLICE-MAXENT3-001 | Splicing | INV: canonical AG acceptor scores above shuffled background; SUB: window must contain the AG |
| 251 | ☐ | SPLICE-MAXENT5-001 | Splicing | INV: canonical GT donor scores above shuffled background; SUB: window must contain the GT |
| 252 | ☐ | MIRNA-CONTEXT-001 | MiRNA | MON: stronger local AU context → more negative score; INV: same site → same score |
| 253 | ☐ | MIRNA-PCT-001 | MiRNA | MON: deeper conservation (longer branch length) → higher PCT; INV: no conservation → PCT 0 |
| 254 | ☐ | MIRNA-CLASSIFY-001 | MiRNA | MON: more native-like MFEI → higher positive probability; INV: di-shuffled sequence → lower probability |
| 255 | ☐ | MIRNA-CLEAVAGE-001 | MiRNA | INV: hsa-miR-21-5p mature reproduced exactly; MON: shifting basal stem shifts Drosha site consistently |
| 256 | ☐ | REP-APPROX-001 | Repeats | MON: introducing a substitution lowers percent-matches; INV: perfect repeat → 100% matches, 0% indels |
| 257 | ☐ | CHROM-ALPHASAT-001 | Chromosome | INV: tandem 171-bp array detected as alpha-satellite; MON: more CENP-B boxes → stronger call |
| 258 | ☐ | CHROM-HOR-001 | Chromosome | MON: more HOR copies → stronger periodicity; INV: pure monomeric array → no HOR |

---

## Summary

| Metric | Value |
|--------|-------|
| Total algorithms | 258 |
| ☑ Complete | 248 |
| ☐ Not started | 10 |
| MR relations defined | ~200+ |
