# Validation Ledger — Phase 1 (Implemented Algorithms)

Independent re-validation of the 86 implemented (☑) test units, one fresh session per unit.
Protocol: [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md).
Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL · ⬜ pending.
State (end of session): ✅ CLEAN (fully functional) · 🔧 LIMITED (see report) · ⬜ pending.

**Progress:** 86 / 86 processed · 85 CLEAN · 1 LIMITED. ✅ PHASE 1 COMPLETE.

> The 86 per-unit validation reports were committed once for provenance and then
> consolidated into this ledger + [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md).
> Recover any report with: `git show cb113ce:docs/Validation/reports/{UNIT-ID}.md`.

**LIMITED units (not fully functional — see [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md)):**
- ~~META-CLASS-001~~ — **FIXED (C1, approved breaking change):** replaced the flat best-hit classifier with the faithful Kraken taxonomy-tree + k-mer-LCA + RTL (max-weight root-to-leaf path, LCA-of-leaves tie-break) classifier (Wood & Salzberg 2014). New public surface (`TaxonomyTree`, int-taxon `BuildKmerDatabase`/`ClassifyReads`); the 27 locked tests were rewritten to the new sourced semantics (28 tests on a hand-built taxonomy, mutation-checked). Now ✅ CLEAN. See [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md) §C1.
- ~~SPLICE-PREDICT-001~~ — **FIXED 2026-06-12**: `GenerateSplicedSequence` now concatenates the same reported-exon set as `DeriveExons` (consistent-filter); INV-3/INV-4 hold by construction; +3 strict tests; now ✅ CLEAN.

**Deferred BIG fixes (revisit later — too large for an in-session "perfect" autofix):**
- ~~META-CLASS-001: real taxonomy-DAG + LCA / weighted root-to-leaf classification~~ — **DONE (C1).** Taxonomy-tree + k-mer-LCA DB + RTL classification implemented (approved breaking API change).
- ~~Phylogenetics: full N-ary (multifurcating) tree model (PhyloNode is Left/Right binary)~~ — **DONE (C2, approved breaking change).** `PhyloNode` now stores `List<PhyloNode> Children` (Left/Right kept as convenience accessors); Newick parses/writes genuine multifurcations `(A,B,C);`, NJ produces its unrooted trifurcation (true round-trip), UPGMA stays binary, RF (rooted-clade + unrooted-bipartition) handles collapsed edges. (The unrooted-bipartition RF metric itself was the separate C3 fix.) See [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md) §C2.
- Alignment: guide-tree progressive MSA (current is star alignment).
- RESTR-DIGEST-001: circular-molecule digest (topology parameter + wrap-around fragments).
- Scientific scoring upgrades (only if desired): Doench/Azimuth on-target, MIT/Hsu off-target weights, SantaLucia NN Tm — currently honest heuristics, not defects.

| # | Unit ID | Area | Stage A | Stage B | State | Report |
|---|---------|------|:---:|:---:|:---:|--------|
| 1 | SEQ-GC-001 | Composition | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 2 | SEQ-COMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 3 | SEQ-REVCOMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 4 | SEQ-VALID-001 | Composition | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 5 | SEQ-COMPLEX-001 | Composition | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 6 | SEQ-ENTROPY-001 | Composition | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 7 | SEQ-GCSKEW-001 | Composition | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 8 | PAT-EXACT-001 | Matching | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 9 | PAT-APPROX-001 | Matching | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 10 | PAT-APPROX-002 | Matching | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 11 | PAT-IUPAC-001 | Matching | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 12 | PAT-PWM-001 | Matching | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 13 | REP-STR-001 | Repeats | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 14 | REP-TANDEM-001 | Repeats | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 15 | REP-INV-001 | Repeats | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 16 | REP-DIRECT-001 | Repeats | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 17 | REP-PALIN-001 | Repeats | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 18 | CRISPR-PAM-001 | MolTools | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 19 | CRISPR-GUIDE-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 20 | CRISPR-OFF-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 21 | PRIMER-TM-001 | MolTools | 🟡 | 🟡 | ✅ CLEAN | archived @cb113ce |
| 22 | PRIMER-DESIGN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 23 | PRIMER-STRUCT-001 | MolTools | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 24 | PROBE-DESIGN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 25 | PROBE-VALID-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 26 | RESTR-FIND-001 | MolTools | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 27 | RESTR-DIGEST-001 | MolTools | 🟡 | 🟡 | ✅ CLEAN | archived @cb113ce |
| 28 | ANNOT-ORF-001 | Annotation | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 29 | ANNOT-GENE-001 | Annotation | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 30 | ANNOT-PROM-001 | Annotation | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 31 | ANNOT-GFF-001 | Annotation | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 32 | KMER-COUNT-001 | K-mer | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 33 | KMER-FREQ-001 | K-mer | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 34 | KMER-FIND-001 | K-mer | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 35 | ALIGN-GLOBAL-001 | Alignment | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 36 | ALIGN-LOCAL-001 | Alignment | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 37 | ALIGN-SEMI-001 | Alignment | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 38 | ALIGN-MULTI-001 | Alignment | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 39 | PHYLO-DIST-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 40 | PHYLO-TREE-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 41 | PHYLO-NEWICK-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 42 | PHYLO-COMP-001 | Phylogenetic | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 43 | POP-FREQ-001 | PopGen | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 44 | POP-DIV-001 | PopGen | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 45 | POP-HW-001 | PopGen | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 46 | POP-FST-001 | PopGen | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 47 | POP-LD-001 | PopGen | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 48 | CHROM-TELO-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 49 | CHROM-CENT-001 | Chromosome | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 50 | CHROM-KARYO-001 | Chromosome | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 51 | CHROM-ANEU-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 52 | CHROM-SYNT-001 | Chromosome | 🟡 | 🟡 | ✅ CLEAN | archived @cb113ce |
| 53 | META-CLASS-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | C1 (Kraken tree+LCA+RTL) |
| 54 | META-PROF-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 55 | META-ALPHA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 56 | META-BETA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 57 | META-BIN-001 | Metagenomics | 🟡 | 🟡 | ✅ CLEAN | archived @cb113ce |
| 58 | CODON-OPT-001 | Codon | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 59 | CODON-CAI-001 | Codon | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 60 | CODON-RARE-001 | Codon | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 61 | CODON-USAGE-001 | Codon | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 62 | TRANS-CODON-001 | Translation | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 63 | TRANS-PROT-001 | Translation | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 64 | PARSE-FASTA-001 | FileIO | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 65 | PARSE-FASTQ-001 | FileIO | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 66 | PARSE-BED-001 | FileIO | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 67 | PARSE-VCF-001 | FileIO | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 68 | PARSE-GFF-001 | FileIO | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 69 | PARSE-GENBANK-001 | FileIO | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 70 | PARSE-EMBL-001 | FileIO | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 71 | RNA-STRUCT-001 | RnaStructure | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 72 | RNA-STEMLOOP-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 73 | RNA-ENERGY-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 74 | MIRNA-SEED-001 | MiRNA | 🟡 | 🟡 | ✅ CLEAN | archived @cb113ce |
| 75 | MIRNA-TARGET-001 | MiRNA | ✅ | 🟡 | ✅ CLEAN | archived @cb113ce |
| 76 | MIRNA-PRECURSOR-001 | MiRNA | 🟡 | 🟡 | ✅ CLEAN | archived @cb113ce |
| 77 | SPLICE-DONOR-001 | Splicing | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 78 | SPLICE-ACCEPTOR-001 | Splicing | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 79 | SPLICE-PREDICT-001 | Splicing | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 80 | DISORDER-PRED-001 | ProteinPred | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 81 | DISORDER-REGION-001 | ProteinPred | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 82 | PROTMOTIF-FIND-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 83 | PROTMOTIF-PROSITE-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 84 | PROTMOTIF-DOMAIN-001 | ProteinMotif | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |
| 85 | EPIGEN-CPG-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | archived @cb113ce |
| 86 | ONCO-IMMUNE-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | archived @cb113ce |

---

# Validation Ledger — Phase 2 (Phase-2 Registry units)

Independent re-validation of the 148 Phase-2 ☑ Registry units absent from Phase 1, one fresh
session per unit, same protocol. Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL · ⬜ pending.
State: ✅ CLEAN · 🔧 LIMITED · ↩︎ DUPLICATE-OF.

**Progress:** 148 / 148 processed · 147 ✅ CLEAN · 1 ↩︎ DUPLICATE-OF · 0 🔧 LIMITED. ✅ PHASE 2 COMPLETE.

> **Campaign result.** Every Phase-2 ☑ Registry unit independently re-validated against external
> first sources, one fresh session per unit. **13 genuine algorithm/behaviour defects** were found
> and completely fixed in-session (all ✅ CLEAN): SV-DETECT-001 (RF/everted→tandem-dup signature),
> SV-CNV-001 (half-integer CN rounding), TRANS-SPLICE-001 (A5SS/A3SS swapped), PANGEN-CORE-001
> (floor→fractional core rule), ASSEMBLY-TRIM-001 (missing cutadapt s<0 early break), MIRNA-PAIR-001
> (DNA-T→U pairing), META-RESIST-001 (ungapped tie-break false negatives), ANNOT-CODING-001 (CPAT
> both-zero hexamer), ANNOT-REPEAT-001 (ClassifyRepeat short-query), GENOMIC-REPEAT-001 (incomplete
> FindRepeats), PROTMOTIF-TM-001 (segment-End off-by-one), SEQ-COMPLEX-COMPRESS-001 (b<2 LZ base
> clamp), RNA-PARTITION-001 (McCaskill base-pair probabilities). Plus numerous green-washed/code-echoing
> tests rewritten to exact externally-sourced values and many untested-branch coverage gaps closed.
> Full suite grew 6484 → 6694 (Failed: 0); SEQ-COMPOSITION-001 is a documented duplicate of SEQ-STATS-001.

| # | Unit ID | Area | Stage A | Stage B | State | Commit |
|---|---------|------|:---:|:---:|:---:|--------|
| 1 | EPIGEN-METHYL-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | - |
| 2 | EPIGEN-DMR-001 | Epigenetics | 🟡 | ✅ | ✅ CLEAN | - |
| 3 | VARIANT-CALL-001 | Variants | ✅ | ✅ | ✅ CLEAN | - |
| 4 | VARIANT-SNP-001 | Variants | ✅ | 🟡 | ✅ CLEAN | - |
| 5 | VARIANT-INDEL-001 | Variants | ✅ | ✅ | ✅ CLEAN | - |
| 6 | VARIANT-ANNOT-001 | Variants | ✅ | ✅ | ✅ CLEAN | - |
| 7 | SV-DETECT-001 | StructuralVar | 🟡 | 🟡 | ✅ CLEAN | - |
| 8 | SV-BREAKPOINT-001 | StructuralVar | ✅ | 🟡 | ✅ CLEAN | - |
| 9 | SV-CNV-001 | StructuralVar | ✅ | 🟡 | ✅ CLEAN | - |
| 10 | ASSEMBLY-OLC-001 | Assembly | ✅ | 🟡 | ✅ CLEAN | - |
| 11 | ASSEMBLY-DBG-001 | Assembly | ✅ | 🟡 | ✅ CLEAN | - |
| 12 | ASSEMBLY-STATS-001 | Assembly | ✅ | ✅ | ✅ CLEAN | - |
| 13 | TRANS-EXPR-001 | Transcriptome | ✅ | 🟡 | ✅ CLEAN | - |
| 14 | TRANS-DIFF-001 | Transcriptome | ✅ | 🟡 | ✅ CLEAN | - |
| 15 | TRANS-SPLICE-001 | Transcriptome | 🟡 | ❌ | ✅ CLEAN | - |
| 16 | COMPGEN-SYNTENY-001 | Comparative Genomics | 🟡 | ✅ | ✅ CLEAN | - |
| 17 | COMPGEN-ORTHO-001 | Comparative Genomics | 🟡 | 🟡 | ✅ CLEAN | - |
| 18 | COMPGEN-REARR-001 | Comparative Genomics | 🟡 | 🟡 | ✅ CLEAN | - |
| 19 | PANGEN-CORE-001 | PanGenome | 🟡 | 🟡 | ✅ CLEAN | - |
| 20 | PANGEN-CLUSTER-001 | PanGenome | ✅ | ✅ | ✅ CLEAN | - |
| 21 | QUALITY-PHRED-001 | Quality | ✅ | ✅ | ✅ CLEAN | - |
| 22 | QUALITY-STATS-001 | Quality | 🟡 | ✅ | ✅ CLEAN | - |
| 23 | SEQ-STATS-001 | Statistics | 🟡 | ✅ | ✅ CLEAN | - |
| 24 | SEQ-MW-001 | Statistics | 🟡 | ✅ | ✅ CLEAN | - |
| 25 | SEQ-PI-001 | Statistics | ✅ | ✅ | ✅ CLEAN | - |
| 26 | SEQ-HYDRO-001 | Statistics | ✅ | 🟡 | ✅ CLEAN | - |
| 27 | SEQ-THERMO-001 | Statistics | ✅ | ✅ | ✅ CLEAN | - |
| 28 | SEQ-DINUC-001 | Statistics | ✅ | ✅ | ✅ CLEAN | - |
| 29 | SEQ-SECSTRUCT-001 | Statistics | 🟡 | ✅ | ✅ CLEAN | - |
| 30 | CODON-RSCU-001 | Codon | ✅ | 🟡 | ✅ CLEAN | - |
| 31 | CODON-ENC-001 | Codon | 🟡 | 🟡 | ✅ CLEAN | - |
| 32 | CODON-STATS-001 | Codon | 🟡 | ✅ | ✅ CLEAN | - |
| 33 | TRANS-SIXFRAME-001 | Translation | 🟡 | 🟡 | ✅ CLEAN | - |
| 34 | ASSEMBLY-MERGE-001 | Assembly | ✅ | ✅ | ✅ CLEAN | - |
| 35 | ASSEMBLY-SCAFFOLD-001 | Assembly | ✅ | ✅ | ✅ CLEAN | - |
| 36 | ASSEMBLY-COVER-001 | Assembly | 🟡 | ✅ | ✅ CLEAN | - |
| 37 | ASSEMBLY-CONSENSUS-001 | Assembly | 🟡 | ✅ | ✅ CLEAN | - |
| 38 | ASSEMBLY-TRIM-001 | Assembly | ✅ | ❌→✅ | ✅ CLEAN | - |
| 39 | ASSEMBLY-CORRECT-001 | Assembly | ✅ | ❌→✅ | ✅ CLEAN | - |
| 40 | PAT-APPROX-003 | Matching | ✅ | 🟡 | ✅ CLEAN | - |
| 41 | ALIGN-STATS-001 | Alignment | 🟡 | 🟡 | ✅ CLEAN | - |
| 42 | EPIGEN-BISULF-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | - |
| 43 | EPIGEN-CHROM-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | - |
| 44 | EPIGEN-AGE-001 | Epigenetics | ✅ | 🟡 | ✅ CLEAN | - |
| 45 | MIRNA-PAIR-001 | MiRNA | ✅ | 🟡 | ✅ CLEAN | - |
| 46 | PANGEN-HEAP-001 | PanGenome | ✅ | 🟡 | ✅ CLEAN | - |
| 47 | PANGEN-MARKER-001 | PanGenome | ✅ | ✅ | ✅ CLEAN | - |
| 48 | POP-SELECT-001 | PopGen | ✅ | 🟡 | ✅ CLEAN | - |
| 49 | POP-ANCESTRY-001 | PopGen | ✅ | ✅ | ✅ CLEAN | - |
| 50 | POP-ROH-001 | PopGen | ✅ | ✅ | ✅ CLEAN | - |
| 51 | META-FUNC-001 | Metagenomics | ✅ | 🟡 | ✅ CLEAN | - |
| 52 | META-RESIST-001 | Metagenomics | 🟡 | ❌→✅ | ✅ CLEAN | - |
| 53 | META-PATHWAY-001 | Metagenomics | 🟡 | 🟡 | ✅ CLEAN | - |
| 54 | META-TAXA-001 | Metagenomics | ✅ | 🟡 | ✅ CLEAN | - |
| 55 | PHYLO-BOOT-001 | Phylogenetics | ✅ | 🟡 | ✅ CLEAN | - |
| 56 | PHYLO-STATS-001 | Phylogenetics | ✅ | ✅ | ✅ CLEAN | - |
| 57 | ANNOT-CODING-001 | Annotation | 🟡 | ❌→✅ | ✅ CLEAN | - |
| 58 | ANNOT-REPEAT-001 | Annotation | ✅ | ❌→✅ | ✅ CLEAN | - |
| 59 | ANNOT-CODONUSAGE-001 | Annotation | ✅ | 🟡 | ✅ CLEAN | - |
| 60 | RESTR-FILTER-001 | MolTools | ✅ | 🟡 | ✅ CLEAN | - |
| 61 | KMER-DIST-001 | K-mer | 🟡 | ✅ | ✅ CLEAN | - |
| 62 | MOTIF-CONS-001 | Matching | 🟡 | ✅ | ✅ CLEAN | - |
| 63 | GENOMIC-REPEAT-001 | Analysis | ✅ | ❌→✅ | ✅ CLEAN | - |
| 64 | GENOMIC-COMMON-001 | Analysis | 🟡 | 🟡 | ✅ CLEAN | - |
| 65 | GENOMIC-MOTIFS-001 | Analysis | ✅ | ✅ | ✅ CLEAN | - |
| 66 | SEQ-RNACOMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | - |
| 67 | PROTMOTIF-PATTERN-001 | ProteinMotif | ✅ | 🟡 | ✅ CLEAN | - |
| 68 | PROTMOTIF-SP-001 | ProteinMotif | ✅ | 🟡 | ✅ CLEAN | - |
| 69 | PROTMOTIF-TM-001 | ProteinMotif | 🟡 | ❌→✅ | ✅ CLEAN | - |
| 70 | PROTMOTIF-CC-001 | ProteinMotif | 🟡 | 🟡 | ✅ CLEAN | - |
| 71 | PROTMOTIF-LC-001 | ProteinMotif | 🟡 | 🟡 | ✅ CLEAN | - |
| 72 | PROTMOTIF-COMMON-001 | ProteinMotif | ✅ | 🟡 | ✅ CLEAN | - |
| 73 | RNA-PAIR-001 | RnaStructure | ✅ | 🟡 | ✅ CLEAN | - |
| 74 | RNA-HAIRPIN-001 | RnaStructure | ✅ | 🟡 | ✅ CLEAN | - |
| 75 | RNA-MFE-001 | RnaStructure | ✅ | 🟡 | ✅ CLEAN | - |
| 76 | RNA-PSEUDOKNOT-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | - |
| 77 | RNA-DOTBRACKET-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | - |
| 78 | RNA-INVERT-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | - |
| 79 | RNA-PARTITION-001 | RnaStructure | 🟡 | ❌→✅ | ✅ CLEAN | - |
| 80 | SEQ-COMPLEX-KMER-001 | Complexity | ✅ | ✅ | ✅ CLEAN | - |
| 81 | SEQ-COMPLEX-WINDOW-001 | Complexity | 🟡 | ✅ | ✅ CLEAN | - |
| 82 | SEQ-COMPLEX-DUST-001 | Complexity | ✅ | ✅ | ✅ CLEAN | - |
| 83 | SEQ-COMPLEX-COMPRESS-001 | Complexity | 🟡 | ❌→✅ | ✅ CLEAN | - |
| 84 | COMPGEN-RBH-001 | Comparative | 🟡 | ✅ | ✅ CLEAN | - |
| 85 | COMPGEN-COMPARE-001 | Comparative | ✅ | ✅ | ✅ CLEAN | - |
| 86 | COMPGEN-REVERSAL-001 | Comparative | ✅ | ✅ | ✅ CLEAN | - |
| 87 | COMPGEN-CLUSTER-001 | Comparative | ✅ | 🟡 | ✅ CLEAN | - |
| 88 | COMPGEN-ANI-001 | Comparative | ✅ | 🟡 | ✅ CLEAN | - |
| 89 | COMPGEN-DOTPLOT-001 | Comparative | ✅ | 🟡 | ✅ CLEAN | - |
| 90 | MOTIF-DISCOVER-001 | Matching | ✅ | 🟡 | ✅ CLEAN | - |
| 91 | MOTIF-SHARED-001 | Matching | ✅ | ✅ | ✅ CLEAN | - |
| 92 | MOTIF-REGULATORY-001 | Matching | ✅ | ✅ | ✅ CLEAN | - |
| 93 | MOTIF-GENERATE-001 | Matching | 🟡 | ✅ | ✅ CLEAN | - |
| 94 | KMER-ASYNC-001 | K-mer Analysis | ✅ | 🟡 | ✅ CLEAN | - |
| 95 | KMER-UNIQUE-001 | K-mer | ✅ | ✅ | ✅ CLEAN | - |
| 96 | KMER-GENERATE-001 | K-mer | ✅ | ✅ | ✅ CLEAN | - |
| 97 | KMER-BOTH-001 | K-mer | ✅ | ✅ | ✅ CLEAN | - |
| 98 | KMER-STATS-001 | K-mer | ✅ | ✅ | ✅ CLEAN | - |
| 99 | KMER-POSITIONS-001 | K-mer | ✅ | ✅ | ✅ CLEAN | - |
| 100 | SEQ-ATSKEW-001 | Composition | ✅ | ✅ | ✅ CLEAN | - |
| 101 | SEQ-REPLICATION-001 | Composition | ✅ | ✅ | ✅ CLEAN | - |
| 102 | SEQ-GC-ANALYSIS-001 | Composition | ✅ | ✅ | ✅ CLEAN | - |
| 103 | DISORDER-MORF-001 | ProteinPred | 🟡 | ✅ | ✅ CLEAN | - |
| 104 | DISORDER-PROPENSITY-001 | ProteinPred | ✅ | ✅ | ✅ CLEAN | - |
| 105 | DISORDER-LC-001 | ProteinPred | 🟡 | ✅ | ✅ CLEAN | - |
| 106 | SEQ-COMPOSITION-001 | Statistics | — | — | ↩︎ DUPLICATE-OF SEQ-STATS-001 | - |
| 107 | SEQ-TM-001 | Statistics | ✅ | ✅ | ✅ CLEAN (DUPLICATE-OF SEQ-THERMO-001) | - |
| 108 | SEQ-ENTROPY-PROFILE-001 | Statistics | ✅ | ✅ | ✅ CLEAN | - |
| 109 | SEQ-GC-PROFILE-001 | Statistics | ✅ | 🟡 | ✅ CLEAN | - |
| 110 | SEQ-CODON-FREQ-001 | Statistics | ✅ | ✅ | ✅ CLEAN | - |
| 111 | SEQ-SUMMARY-001 | Statistics | 🟡 | ✅ | ✅ CLEAN | - |
| 112 | GENOMIC-TANDEM-001 | Analysis | ✅ | 🟡 | ✅ CLEAN (DUPLICATE-OF REP-TANDEM-001) | - |
| 113 | GENOMIC-SIMILARITY-001 | Analysis | ✅ | ✅ | ✅ CLEAN | - |
| 114 | GENOMIC-ORF-001 | Analysis | ✅ | ✅ | ✅ CLEAN | - |
| 115 | ONCO-SOMATIC-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 116 | ONCO-VAF-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 117 | ONCO-DRIVER-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 118 | ONCO-ARTIFACT-001 | Oncology | ✅ | 🟡 | ✅ CLEAN | - |
| 119 | ONCO-ANNOT-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | - |
| 120 | ONCO-TMB-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 121 | ONCO-MSI-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 122 | ONCO-HRD-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 123 | ONCO-LOH-001 | Oncology | ✅ | 🟡 | ✅ CLEAN | - |
| 124 | ONCO-SIG-001 | Oncology | ✅ | 🟡 | ✅ CLEAN | - |
| 125 | ONCO-SIG-002 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 126 | ONCO-SIG-003 | Oncology | 🟡 | 🟡 | ✅ CLEAN | - |
| 127 | ONCO-SIG-004 | Oncology | ✅ | 🟡 | ✅ CLEAN | - |
| 128 | ONCO-FUSION-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 129 | ONCO-FUSION-002 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 130 | ONCO-FUSION-003 | Oncology | 🟡 | 🟡 | ✅ CLEAN | - |
| 131 | ONCO-CNA-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 132 | ONCO-CNA-002 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 133 | ONCO-CNA-003 | Oncology | 🟡 | ✅ | ✅ CLEAN | - |
| 134 | ONCO-PURITY-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 135 | ONCO-PLOIDY-001 | Oncology | ✅ | 🟡 | ✅ CLEAN | - |
| 136 | ONCO-CLONAL-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 137 | ONCO-NEO-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 138 | ONCO-MHC-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | - |
| 139 | ONCO-CTDNA-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | - |
| 140 | ONCO-MRD-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | - |
| 141 | ONCO-CHIP-001 | Oncology | 🟡 | 🟡 | ✅ CLEAN | - |
| 142 | ONCO-PHYLO-001 | Oncology | ✅ | 🟡 | ✅ CLEAN | - |
| 143 | ONCO-CCF-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 144 | ONCO-HETERO-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 145 | ONCO-HLA-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 146 | ONCO-ACTION-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |
| 147 | ONCO-SV-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | - |
| 148 | ONCO-EXPR-001 | Oncology | ✅ | ✅ | ✅ CLEAN | - |

---

# Validation Ledger — Phase 3 (Independent re-validation of enhanced units)

After the open-questions program implemented 8 enhancements (Group 1 features + Group 2 breaking
changes), each touched unit is **independently re-validated in a fresh session** against external
first sources — the implementer and the validator are different contexts, per the protocol's premise
that code + tests authored together share blind spots. Scope = only the units whose code/tests changed.
Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL. State: ✅ CLEAN · 🔧 LIMITED.

| # | Unit ID | Enhancement re-validated | Stage A | Stage B | State | Commit |
|---|---------|--------------------------|:---:|:---:|:---:|--------|
| 1 | ANNOT-GENE-001 | Reverse-strand Shine-Dalgarno / RBS reporting (`FindRibosomeBindingSitesBothStrands`; `forwardPos = len − revPos − motifLen` mapping) | ✅ | ✅ | ✅ CLEAN | - |
| 2 | RESTR-DIGEST-001 | Circular-molecule digest (`Digest(DnaSequence, MoleculeTopology, params string[])`; `MoleculeTopology{Linear,Circular}`; circular k sites → k fragments, origin-spanning join len `(len−lastCut)+firstCut`; 0→1 uncut, 1→1 linearized) | ✅ | ✅ | ✅ CLEAN | - |
