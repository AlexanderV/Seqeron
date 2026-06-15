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
| 53 | META-CLASS-001 | Metagenomics | 🟡 | 🟡 | 🔧 LIMITED | archived @cb113ce |
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

**Progress:** 6 / 148 processed.

| # | Unit ID | Area | Stage A | Stage B | State | Commit |
|---|---------|------|:---:|:---:|:---:|--------|
| 1 | EPIGEN-METHYL-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | - |
| 2 | EPIGEN-DMR-001 | Epigenetics | 🟡 | ✅ | ✅ CLEAN | - |
| 3 | VARIANT-CALL-001 | Variants | ✅ | ✅ | ✅ CLEAN | - |
| 4 | VARIANT-SNP-001 | Variants | ✅ | 🟡 | ✅ CLEAN | - |
| 5 | VARIANT-INDEL-001 | Variants | ✅ | ✅ | ✅ CLEAN | - |
| 6 | VARIANT-ANNOT-001 | Variants | ✅ | ✅ | ✅ CLEAN | - |
