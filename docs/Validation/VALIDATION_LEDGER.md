# Validation Ledger — Phase 1 (Implemented Algorithms)

Independent re-validation of the 86 implemented (☑) test units, one fresh session per unit.
Protocol: [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md).
Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL · ⬜ pending.
State (end of session): ✅ CLEAN (fully functional) · 🔧 LIMITED (see report) · ⬜ pending.

**Progress:** 86 / 86 processed · 85 CLEAN · 1 LIMITED. ✅ PHASE 1 COMPLETE.

**LIMITED units (not fully functional — see reports):**
- META-CLASS-001 — implements flat best-hit, not Kraken/LCA. **Overclaim now de-claimed** in code XML docs + spec + Evidence (2026-06-12). The remaining gap — a real taxonomy-DAG + LCA classifier — is NOT-POSSIBLE in-session (changes public API + 27 locked tests); see [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md) §C1.
- ~~SPLICE-PREDICT-001~~ — **FIXED 2026-06-12**: `GenerateSplicedSequence` now concatenates the same reported-exon set as `DeriveExons` (consistent-filter); INV-3/INV-4 hold by construction; +3 strict tests; now ✅ CLEAN.

**Deferred BIG fixes (revisit later — too large for an in-session "perfect" autofix):**
- META-CLASS-001: real taxonomy-DAG + LCA / weighted root-to-leaf classification (or rename method + de-overclaim docs).
- Phylogenetics: full N-ary (multifurcating) tree model (PhyloNode is Left/Right binary) — enables true NJ trifurcation round-trip; also an optional unrooted-bipartition RF metric alongside the current rooted-clade one.
- Alignment: guide-tree progressive MSA (current is star alignment).
- RESTR-DIGEST-001: circular-molecule digest (topology parameter + wrap-around fragments).
- Scientific scoring upgrades (only if desired): Doench/Azimuth on-target, MIT/Hsu off-target weights, SantaLucia NN Tm — currently honest heuristics, not defects.

| # | Unit ID | Area | Stage A | Stage B | State | Report |
|---|---------|------|:---:|:---:|:---:|--------|
| 1 | SEQ-GC-001 | Composition | 🟡 | ✅ | ✅ CLEAN | [report](reports/SEQ-GC-001.md) |
| 2 | SEQ-COMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | [report](reports/SEQ-COMP-001.md) |
| 3 | SEQ-REVCOMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | [report](reports/SEQ-REVCOMP-001.md) |
| 4 | SEQ-VALID-001 | Composition | ✅ | ✅ | ✅ CLEAN | [report](reports/SEQ-VALID-001.md) |
| 5 | SEQ-COMPLEX-001 | Composition | 🟡 | ✅ | ✅ CLEAN | [report](reports/SEQ-COMPLEX-001.md) |
| 6 | SEQ-ENTROPY-001 | Composition | 🟡 | ✅ | ✅ CLEAN | [report](reports/SEQ-ENTROPY-001.md) |
| 7 | SEQ-GCSKEW-001 | Composition | 🟡 | ✅ | ✅ CLEAN | [report](reports/SEQ-GCSKEW-001.md) |
| 8 | PAT-EXACT-001 | Matching | ✅ | ✅ | ✅ CLEAN | [report](reports/PAT-EXACT-001.md) |
| 9 | PAT-APPROX-001 | Matching | ✅ | ✅ | ✅ CLEAN | [report](reports/PAT-APPROX-001.md) |
| 10 | PAT-APPROX-002 | Matching | ✅ | ✅ | ✅ CLEAN | [report](reports/PAT-APPROX-002.md) |
| 11 | PAT-IUPAC-001 | Matching | ✅ | ✅ | ✅ CLEAN | [report](reports/PAT-IUPAC-001.md) |
| 12 | PAT-PWM-001 | Matching | ✅ | ✅ | ✅ CLEAN | [report](reports/PAT-PWM-001.md) |
| 13 | REP-STR-001 | Repeats | ✅ | 🟡 | ✅ CLEAN | [report](reports/REP-STR-001.md) |
| 14 | REP-TANDEM-001 | Repeats | ✅ | ✅ | ✅ CLEAN | [report](reports/REP-TANDEM-001.md) |
| 15 | REP-INV-001 | Repeats | ✅ | ✅ | ✅ CLEAN | [report](reports/REP-INV-001.md) |
| 16 | REP-DIRECT-001 | Repeats | ✅ | ✅ | ✅ CLEAN | [report](reports/REP-DIRECT-001.md) |
| 17 | REP-PALIN-001 | Repeats | 🟡 | ✅ | ✅ CLEAN | [report](reports/REP-PALIN-001.md) |
| 18 | CRISPR-PAM-001 | MolTools | ✅ | 🟡 | ✅ CLEAN | [report](reports/CRISPR-PAM-001.md) |
| 19 | CRISPR-GUIDE-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | [report](reports/CRISPR-GUIDE-001.md) |
| 20 | CRISPR-OFF-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | [report](reports/CRISPR-OFF-001.md) |
| 21 | PRIMER-TM-001 | MolTools | 🟡 | 🟡 | ✅ CLEAN | [report](reports/PRIMER-TM-001.md) |
| 22 | PRIMER-DESIGN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | [report](reports/PRIMER-DESIGN-001.md) |
| 23 | PRIMER-STRUCT-001 | MolTools | ✅ | ✅ | ✅ CLEAN | [report](reports/PRIMER-STRUCT-001.md) |
| 24 | PROBE-DESIGN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | [report](reports/PROBE-DESIGN-001.md) |
| 25 | PROBE-VALID-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | [report](reports/PROBE-VALID-001.md) |
| 26 | RESTR-FIND-001 | MolTools | ✅ | ✅ | ✅ CLEAN | [report](reports/RESTR-FIND-001.md) |
| 27 | RESTR-DIGEST-001 | MolTools | 🟡 | 🟡 | ✅ CLEAN | [report](reports/RESTR-DIGEST-001.md) |
| 28 | ANNOT-ORF-001 | Annotation | ✅ | 🟡 | ✅ CLEAN | [report](reports/ANNOT-ORF-001.md) |
| 29 | ANNOT-GENE-001 | Annotation | ✅ | 🟡 | ✅ CLEAN | [report](reports/ANNOT-GENE-001.md) |
| 30 | ANNOT-PROM-001 | Annotation | ✅ | ✅ | ✅ CLEAN | [report](reports/ANNOT-PROM-001.md) |
| 31 | ANNOT-GFF-001 | Annotation | 🟡 | ✅ | ✅ CLEAN | [report](reports/ANNOT-GFF-001.md) |
| 32 | KMER-COUNT-001 | K-mer | ✅ | ✅ | ✅ CLEAN | [report](reports/KMER-COUNT-001.md) |
| 33 | KMER-FREQ-001 | K-mer | ✅ | ✅ | ✅ CLEAN | [report](reports/KMER-FREQ-001.md) |
| 34 | KMER-FIND-001 | K-mer | ✅ | ✅ | ✅ CLEAN | [report](reports/KMER-FIND-001.md) |
| 35 | ALIGN-GLOBAL-001 | Alignment | ✅ | ✅ | ✅ CLEAN | [report](reports/ALIGN-GLOBAL-001.md) |
| 36 | ALIGN-LOCAL-001 | Alignment | ✅ | ✅ | ✅ CLEAN | [report](reports/ALIGN-LOCAL-001.md) |
| 37 | ALIGN-SEMI-001 | Alignment | ✅ | ✅ | ✅ CLEAN | [report](reports/ALIGN-SEMI-001.md) |
| 38 | ALIGN-MULTI-001 | Alignment | 🟡 | ✅ | ✅ CLEAN | [report](reports/ALIGN-MULTI-001.md) |
| 39 | PHYLO-DIST-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | [report](reports/PHYLO-DIST-001.md) |
| 40 | PHYLO-TREE-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | [report](reports/PHYLO-TREE-001.md) |
| 41 | PHYLO-NEWICK-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | [report](reports/PHYLO-NEWICK-001.md) |
| 42 | PHYLO-COMP-001 | Phylogenetic | 🟡 | ✅ | ✅ CLEAN | [report](reports/PHYLO-COMP-001.md) |
| 43 | POP-FREQ-001 | PopGen | ✅ | ✅ | ✅ CLEAN | [report](reports/POP-FREQ-001.md) |
| 44 | POP-DIV-001 | PopGen | ✅ | ✅ | ✅ CLEAN | [report](reports/POP-DIV-001.md) |
| 45 | POP-HW-001 | PopGen | ✅ | ✅ | ✅ CLEAN | [report](reports/POP-HW-001.md) |
| 46 | POP-FST-001 | PopGen | ✅ | ✅ | ✅ CLEAN | [report](reports/POP-FST-001.md) |
| 47 | POP-LD-001 | PopGen | ✅ | ✅ | ✅ CLEAN | [report](reports/POP-LD-001.md) |
| 48 | CHROM-TELO-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | [report](reports/CHROM-TELO-001.md) |
| 49 | CHROM-CENT-001 | Chromosome | 🟡 | ✅ | ✅ CLEAN | [report](reports/CHROM-CENT-001.md) |
| 50 | CHROM-KARYO-001 | Chromosome | ✅ | 🟡 | ✅ CLEAN | [report](reports/CHROM-KARYO-001.md) |
| 51 | CHROM-ANEU-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | [report](reports/CHROM-ANEU-001.md) |
| 52 | CHROM-SYNT-001 | Chromosome | 🟡 | 🟡 | ✅ CLEAN | [report](reports/CHROM-SYNT-001.md) |
| 53 | META-CLASS-001 | Metagenomics | 🟡 | 🟡 | 🔧 LIMITED | [report](reports/META-CLASS-001.md) |
| 54 | META-PROF-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | [report](reports/META-PROF-001.md) |
| 55 | META-ALPHA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | [report](reports/META-ALPHA-001.md) |
| 56 | META-BETA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | [report](reports/META-BETA-001.md) |
| 57 | META-BIN-001 | Metagenomics | 🟡 | 🟡 | ✅ CLEAN | [report](reports/META-BIN-001.md) |
| 58 | CODON-OPT-001 | Codon | ✅ | ✅ | ✅ CLEAN | [report](reports/CODON-OPT-001.md) |
| 59 | CODON-CAI-001 | Codon | 🟡 | ✅ | ✅ CLEAN | [report](reports/CODON-CAI-001.md) |
| 60 | CODON-RARE-001 | Codon | 🟡 | ✅ | ✅ CLEAN | [report](reports/CODON-RARE-001.md) |
| 61 | CODON-USAGE-001 | Codon | 🟡 | ✅ | ✅ CLEAN | [report](reports/CODON-USAGE-001.md) |
| 62 | TRANS-CODON-001 | Translation | ✅ | ✅ | ✅ CLEAN | [report](reports/TRANS-CODON-001.md) |
| 63 | TRANS-PROT-001 | Translation | ✅ | ✅ | ✅ CLEAN | [report](reports/TRANS-PROT-001.md) |
| 64 | PARSE-FASTA-001 | FileIO | ✅ | ✅ | ✅ CLEAN | [report](reports/PARSE-FASTA-001.md) |
| 65 | PARSE-FASTQ-001 | FileIO | ✅ | ✅ | ✅ CLEAN | [report](reports/PARSE-FASTQ-001.md) |
| 66 | PARSE-BED-001 | FileIO | ✅ | ✅ | ✅ CLEAN | [report](reports/PARSE-BED-001.md) |
| 67 | PARSE-VCF-001 | FileIO | 🟡 | ✅ | ✅ CLEAN | [report](reports/PARSE-VCF-001.md) |
| 68 | PARSE-GFF-001 | FileIO | ✅ | ✅ | ✅ CLEAN | [report](reports/PARSE-GFF-001.md) |
| 69 | PARSE-GENBANK-001 | FileIO | ✅ | ✅ | ✅ CLEAN | [report](reports/PARSE-GENBANK-001.md) |
| 70 | PARSE-EMBL-001 | FileIO | ✅ | ✅ | ✅ CLEAN | [report](reports/PARSE-EMBL-001.md) |
| 71 | RNA-STRUCT-001 | RnaStructure | ✅ | 🟡 | ✅ CLEAN | [report](reports/RNA-STRUCT-001.md) |
| 72 | RNA-STEMLOOP-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | [report](reports/RNA-STEMLOOP-001.md) |
| 73 | RNA-ENERGY-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | [report](reports/RNA-ENERGY-001.md) |
| 74 | MIRNA-SEED-001 | MiRNA | 🟡 | 🟡 | ✅ CLEAN | [report](reports/MIRNA-SEED-001.md) |
| 75 | MIRNA-TARGET-001 | MiRNA | ✅ | 🟡 | ✅ CLEAN | [report](reports/MIRNA-TARGET-001.md) |
| 76 | MIRNA-PRECURSOR-001 | MiRNA | 🟡 | 🟡 | ✅ CLEAN | [report](reports/MIRNA-PRECURSOR-001.md) |
| 77 | SPLICE-DONOR-001 | Splicing | ✅ | ✅ | ✅ CLEAN | [report](reports/SPLICE-DONOR-001.md) |
| 78 | SPLICE-ACCEPTOR-001 | Splicing | 🟡 | ✅ | ✅ CLEAN | [report](reports/SPLICE-ACCEPTOR-001.md) |
| 79 | SPLICE-PREDICT-001 | Splicing | ✅ | ✅ | ✅ CLEAN | [report](reports/SPLICE-PREDICT-001.md) |
| 80 | DISORDER-PRED-001 | ProteinPred | ✅ | ✅ | ✅ CLEAN | [report](reports/DISORDER-PRED-001.md) |
| 81 | DISORDER-REGION-001 | ProteinPred | 🟡 | ✅ | ✅ CLEAN | [report](reports/DISORDER-REGION-001.md) |
| 82 | PROTMOTIF-FIND-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | [report](reports/PROTMOTIF-FIND-001.md) |
| 83 | PROTMOTIF-PROSITE-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | [report](reports/PROTMOTIF-PROSITE-001.md) |
| 84 | PROTMOTIF-DOMAIN-001 | ProteinMotif | 🟡 | ✅ | ✅ CLEAN | [report](reports/PROTMOTIF-DOMAIN-001.md) |
| 85 | EPIGEN-CPG-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | [report](reports/EPIGEN-CPG-001.md) |
| 86 | ONCO-IMMUNE-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | [report](reports/ONCO-IMMUNE-001.md) |
