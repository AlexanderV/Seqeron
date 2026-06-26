# Validation Ledger — Phase 1 (Implemented Algorithms)

> ⚠️ **FULL RE-VALIDATION RESET — 2026-06-24.** Every unit below was reset to **⬜ pending** for a fresh end-to-end re-verification (extensive code changes during the limitations-elimination campaign). Prior Stage/State results are **SUPERSEDED**, retained only as historical evidence; the per-unit reports are kept. Registry expanded to **237** (added ONCO-ASCAT-001, RNA-PKPREDICT-001, RNA-PKRECURSIVE-001 — pending first validation).
>
> ✅ **PHASE 1 RE-VALIDATION COMPLETE — 2026-06-24.** All **86 / 86** units independently re-validated, one fresh-context session per unit (Stage A description → Stage B implementation, authoritative external sources, hand/reference cross-checks). **Result: 86 / 86 ✅ CLEAN.** Stage-A: 63 PASS / 23 PASS-WITH-NOTES / 0 FAIL. Stage-B: 79 PASS / 7 PASS-WITH-NOTES / 0 FAIL. (Spec-alignment pass 2026-06-24 cleared 6 doc-only PASS-WITH-NOTES → PASS: SEQ-COMPLEX-001, REP-PALIN-001, PRIMER-DESIGN-001, CHROM-SYNT-001, MIRNA-SEED-001, PROTMOTIF-FIND-001 — stale/imprecise TestSpec wording corrected; code unchanged.) One real implementation defect found & fully fixed in-session (PARSE-GENBANK-001 multi-line qualifier reconstruction → spurious-space corruption of wrapped `/translation`); minor test/spec corrections elsewhere (SEQ-COMPLEX-001, TRANS-PROT-001 +3 lock tests, PROTMOTIF-FIND/DOMAIN spec fixes). Full unfiltered suite green at the end (18213 passed / 0 failed). The PASS-WITH-NOTES are documented by-design scope boundaries (declared heuristics, single-codon-AA CAI inclusion, stale-spec wording), not defects. RNA-STRUCT-001's non-csr-PK / tertiary-stabilised-knot residual and ONCO-MHC-001's trained-model gap remain the only LIMITATIONS.md entries (by-design, validated).
>
> 🔁 **POST-COMPLETION RE-RESET — 2026-06-25.** The limitation-elimination campaign continued after the 2026-06-24 sweep (Tiers G–N: McCaskill, HMMER E-value/null2/envelope, ntthal dimer/hairpin, MaxEntScan score3/score5, context++/PCT, ABIS, MHCflurry, CheckM markers, miRBase classifier, …). The following **19** previously-CLEAN units were touched again and are reset to **⬜ pending** re-validation: ANNOT-GFF-001, CHROM-CENT-001, CODON-CAI-001, CODON-RARE-001, DISORDER-REGION-001, META-BIN-001, MIRNA-PRECURSOR-001, MIRNA-TARGET-001, ONCO-IMMUNE-001, ONCO-MRD-001, PARSE-EMBL-001, PARSE-FASTA-001, PRIMER-TM-001, PROBE-DESIGN-001, PROBE-VALID-001, PROTMOTIF-DOMAIN-001, REP-STR-001, SPLICE-ACCEPTOR-001, SPLICE-DONOR-001 (ONCO-MHC-001 was already pending). Their per-unit reports are kept as historical evidence.

Independent re-validation of the 86 implemented (☑) test units, one fresh session per unit.
Protocol: [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md).
Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL · ⬜ pending.
State (end of session): ✅ CLEAN (fully functional) · 🔧 LIMITED (see report) · ⬜ pending.

**Progress:** 86 / 86 ✅ CLEAN — Phase 1 re-validation COMPLETE 2026-06-24 (0 FAIL, 0 LIMITED; 1 defect found & fixed: PARSE-GENBANK-001). — **NOTE 2026-06-25:** 19 of these units were changed again by later campaign tiers (G–N) and are re-reset to ⬜ pending below; first-time validation is also pending for the 21 new units in *New units (campaign)*. Re-validation count is therefore **no longer 86/86**.

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
- ~~Scientific scoring upgrades: Doench/Azimuth on-target, MIT/Hsu off-target weights~~ — **DONE (C7, fully resolved).** Doench 2014 "Rule Set 1" + Doench 2016 "Rule Set 2"/Azimuth on-target (`CalculateOnTargetDoench2014` / `CalculateOnTargetRuleSet2`) and MIT/Hsu 2013 + CFD 2016 off-target (`CalculateMitHitScore` / `CalculateCfdScore`) are all implemented from authoritative sources / trained models (Rule Set 2 reconstructed sklearn-free from the Azimuth pickles; reproducer `scripts/azimuth/extract_azimuth_model.py`). See [FINDINGS_REGISTER.md](FINDINGS_REGISTER.md) §C7 and `reports/CRISPR-GUIDE-001.md`. SantaLucia NN Tm is itself validated under SEQ-THERMO-001; only the primer-design *convenience* scoring remains an honest heuristic (not a defect).

| # | Unit ID | Area | Stage A | Stage B | State | Report |
|---|---------|------|:---:|:---:|:---:|--------|
| 1 | SEQ-GC-001 | Composition | 🟡 | ✅ | ✅ CLEAN | reports/SEQ-GC-001.md |
| 2 | SEQ-COMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | reports/SEQ-COMP-001.md |
| 3 | SEQ-REVCOMP-001 | Composition | ✅ | ✅ | ✅ CLEAN | reports/SEQ-REVCOMP-001.md |
| 4 | SEQ-VALID-001 | Composition | ✅ | ✅ | ✅ CLEAN | reports/SEQ-VALID-001.md |
| 5 | SEQ-COMPLEX-001 | Composition | ✅ | ✅ | ✅ CLEAN | reports/SEQ-COMPLEX-001.md |
| 6 | SEQ-ENTROPY-001 | Composition | ✅ | ✅ | ✅ CLEAN | reports/SEQ-ENTROPY-001.md |
| 7 | SEQ-GCSKEW-001 | Composition | 🟡 | ✅ | ✅ CLEAN | reports/SEQ-GCSKEW-001.md |
| 8 | PAT-EXACT-001 | Matching | ✅ | ✅ | ✅ CLEAN | reports/PAT-EXACT-001.md |
| 9 | PAT-APPROX-001 | Matching | ✅ | ✅ | ✅ CLEAN | reports/PAT-APPROX-001.md |
| 10 | PAT-APPROX-002 | Matching | ✅ | ✅ | ✅ CLEAN | reports/PAT-APPROX-002.md |
| 11 | PAT-IUPAC-001 | Matching | ✅ | ✅ | ✅ CLEAN | reports/PAT-IUPAC-001.md |
| 12 | PAT-PWM-001 | Matching | ✅ | ✅ | ✅ CLEAN | reports/PAT-PWM-001.md |
| 13 | REP-STR-001 | Repeats | ✅ | 🟡 | ✅ CLEAN | reports/REP-STR-001.md |
| 14 | REP-TANDEM-001 | Repeats | ✅ | ✅ | ✅ CLEAN | reports/REP-TANDEM-001.md |
| 15 | REP-INV-001 | Repeats | ✅ | ✅ | ✅ CLEAN | reports/REP-INV-001.md |
| 16 | REP-DIRECT-001 | Repeats | ✅ | ✅ | ✅ CLEAN | reports/REP-DIRECT-001.md |
| 17 | REP-PALIN-001 | Repeats | ✅ | ✅ | ✅ CLEAN | reports/REP-PALIN-001.md |
| 18 | CRISPR-PAM-001 | MolTools | ✅ | 🟡 | ✅ CLEAN | reports/CRISPR-PAM-001.md |
| 19 | CRISPR-GUIDE-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/CRISPR-GUIDE-001.md |
| 20 | CRISPR-OFF-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/CRISPR-OFF-001.md |
| 21 | PRIMER-TM-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | reports/PRIMER-TM-001.md |
| 22 | PRIMER-DESIGN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PRIMER-DESIGN-001.md |
| 23 | PRIMER-STRUCT-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PRIMER-STRUCT-001.md |
| 24 | PROBE-DESIGN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PROBE-DESIGN-001.md |
| 25 | PROBE-VALID-001 | MolTools | 🟡 | ✅ | ✅ CLEAN | reports/PROBE-VALID-001.md |
| 26 | RESTR-FIND-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/RESTR-FIND-001.md |
| 27 | RESTR-DIGEST-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/RESTR-DIGEST-001.md |
| 28 | ANNOT-ORF-001 | Annotation | ✅ | 🟡 | ✅ CLEAN | reports/ANNOT-ORF-001.md |
| 29 | ANNOT-GENE-001 | Annotation | ✅ | ✅ | ✅ CLEAN | reports/ANNOT-GENE-001.md |
| 30 | ANNOT-PROM-001 | Annotation | ✅ | ✅ | ✅ CLEAN | reports/ANNOT-PROM-001.md |
| 31 | ANNOT-GFF-001 | Annotation | ✅ | ✅ | ✅ CLEAN | reports/ANNOT-GFF-001.md |
| 32 | KMER-COUNT-001 | K-mer | ✅ | ✅ | ✅ CLEAN | reports/KMER-COUNT-001.md |
| 33 | KMER-FREQ-001 | K-mer | ✅ | ✅ | ✅ CLEAN | reports/KMER-FREQ-001.md |
| 34 | KMER-FIND-001 | K-mer | ✅ | ✅ | ✅ CLEAN | reports/KMER-FIND-001.md |
| 35 | ALIGN-GLOBAL-001 | Alignment | ✅ | ✅ | ✅ CLEAN | reports/ALIGN-GLOBAL-001.md |
| 36 | ALIGN-LOCAL-001 | Alignment | ✅ | ✅ | ✅ CLEAN | reports/ALIGN-LOCAL-001.md |
| 37 | ALIGN-SEMI-001 | Alignment | ✅ | ✅ | ✅ CLEAN | reports/ALIGN-SEMI-001.md |
| 38 | ALIGN-MULTI-001 | Alignment | ✅ | ✅ | ✅ CLEAN | reports/ALIGN-MULTI-001.md |
| 39 | PHYLO-DIST-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | reports/PHYLO-DIST-001.md |
| 40 | PHYLO-TREE-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | reports/PHYLO-TREE-001.md |
| 41 | PHYLO-NEWICK-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | reports/PHYLO-NEWICK-001.md |
| 42 | PHYLO-COMP-001 | Phylogenetic | ✅ | ✅ | ✅ CLEAN | reports/PHYLO-COMP-001.md |
| 43 | POP-FREQ-001 | PopGen | ✅ | ✅ | ✅ CLEAN | reports/POP-FREQ-001.md |
| 44 | POP-DIV-001 | PopGen | ✅ | ✅ | ✅ CLEAN | reports/POP-DIV-001.md |
| 45 | POP-HW-001 | PopGen | ✅ | ✅ | ✅ CLEAN | reports/POP-HW-001.md |
| 46 | POP-FST-001 | PopGen | ✅ | ✅ | ✅ CLEAN | reports/POP-FST-001.md |
| 47 | POP-LD-001 | PopGen | ✅ | ✅ | ✅ CLEAN | reports/POP-LD-001.md |
| 48 | CHROM-TELO-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-TELO-001.md |
| 49 | CHROM-CENT-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-CENT-001.md |
| 50 | CHROM-KARYO-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-KARYO-001.md |
| 51 | CHROM-ANEU-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-ANEU-001.md |
| 52 | CHROM-SYNT-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-SYNT-001.md |
| 53 | META-CLASS-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | reports/META-CLASS-001.md |
| 54 | META-PROF-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | reports/META-PROF-001.md |
| 55 | META-ALPHA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | reports/META-ALPHA-001.md |
| 56 | META-BETA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | reports/META-BETA-001.md |
| 57 | META-BIN-001 | Metagenomics | 🟡 | 🟡 | ✅ CLEAN | reports/META-BIN-001.md |
| 58 | CODON-OPT-001 | Codon | ✅ | ✅ | ✅ CLEAN | reports/CODON-OPT-001.md |
| 59 | CODON-CAI-001 | Codon | ✅ | ✅ | ✅ CLEAN | reports/CODON-CAI-001.md |
| 60 | CODON-RARE-001 | Codon | 🟡 | ✅ | ✅ CLEAN | reports/CODON-RARE-001.md |
| 61 | CODON-USAGE-001 | Codon | 🟡 | ✅ | ✅ CLEAN | reports/CODON-USAGE-001.md |
| 62 | TRANS-CODON-001 | Translation | ✅ | ✅ | ✅ CLEAN | reports/TRANS-CODON-001.md |
| 63 | TRANS-PROT-001 | Translation | ✅ | 🟡 | ✅ CLEAN | reports/TRANS-PROT-001.md |
| 64 | PARSE-FASTA-001 | FileIO | ✅ | ✅ | ✅ CLEAN | reports/PARSE-FASTA-001.md |
| 65 | PARSE-FASTQ-001 | FileIO | ✅ | 🟡 | ✅ CLEAN | reports/PARSE-FASTQ-001.md |
| 66 | PARSE-BED-001 | FileIO | ✅ | ✅ | ✅ CLEAN | reports/PARSE-BED-001.md |
| 67 | PARSE-VCF-001 | FileIO | 🟡 | ✅ | ✅ CLEAN | reports/PARSE-VCF-001.md |
| 68 | PARSE-GFF-001 | FileIO | ✅ | ✅ | ✅ CLEAN | reports/PARSE-GFF-001.md |
| 69 | PARSE-GENBANK-001 | FileIO | ✅ | 🟡 | ✅ CLEAN (fixed multi-line qualifier wrap) | reports/PARSE-GENBANK-001.md |
| 70 | PARSE-EMBL-001 | FileIO | ⬜ | ⬜ | ⬜ pending | reports/PARSE-EMBL-001.md |
| 71 | RNA-STRUCT-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | reports/RNA-STRUCT-001.md |
| 72 | RNA-STEMLOOP-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | reports/RNA-STEMLOOP-001.md |
| 73 | RNA-ENERGY-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | reports/RNA-ENERGY-001.md |
| 74 | MIRNA-SEED-001 | MiRNA | ✅ | ✅ | ✅ CLEAN | reports/MIRNA-SEED-001.md |
| 75 | MIRNA-TARGET-001 | MiRNA | ✅ | ✅ | ✅ CLEAN (re-validated 2026-06-26: `TA_3UTR = log10(N)`, N = **non-overlapping** 8mer+7mer-m8+7mer-A1 sites — confirmed verbatim vs Garcia 2011 Online Methods PMC3190056 + Agarwal 2015 Table 1 + TargetScan 7.0 perl. Found & fully fixed a non-overlapping over-count on periodic cores; +CTX-TA-009 lock; pre-existing site/scoring surface still holds; residual PCT/SPS/Len_ORF/ORF8m caller-supplied) | reports/MIRNA-TARGET-001.md |
| 76 | MIRNA-PRECURSOR-001 | MiRNA | ✅ | ✅ | ✅ CLEAN | reports/MIRNA-PRECURSOR-001.md |
| 77 | SPLICE-DONOR-001 | Splicing | ✅ | ✅ | ✅ CLEAN | reports/SPLICE-DONOR-001.md |
| 78 | SPLICE-ACCEPTOR-001 | Splicing | ✅ | ✅ | ✅ CLEAN | reports/SPLICE-ACCEPTOR-001.md |
| 79 | SPLICE-PREDICT-001 | Splicing | ✅ | ✅ | ✅ CLEAN | reports/SPLICE-PREDICT-001.md |
| 80 | DISORDER-PRED-001 | ProteinPred | ✅ | ✅ | ✅ CLEAN | reports/DISORDER-PRED-001.md |
| 81 | DISORDER-REGION-001 | ProteinPred | 🟡 | ✅ | ✅ CLEAN | reports/DISORDER-REGION-001.md |
| 82 | PROTMOTIF-FIND-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | reports/PROTMOTIF-FIND-001.md |
| 83 | PROTMOTIF-PROSITE-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | reports/PROTMOTIF-PROSITE-001.md |
| 84 | PROTMOTIF-DOMAIN-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | reports/PROTMOTIF-DOMAIN-001.md |
| 85 | EPIGEN-CPG-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | reports/EPIGEN-CPG-001.md |
| 86 | ONCO-IMMUNE-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | reports/ONCO-IMMUNE-001.md |

## New units (created during the limitations-elimination campaign)

First-time validation of the **24** units added while eliminating LIMITATIONS.md items (the original 3 below + **21** net-new algorithms from the later campaign tiers — McCaskill accessibility, Plan7 HMMER, ntthal dimer/hairpin Tm, MHCflurry NN, SMM/BIMAS, ν-SVR, CheckM, TETRA, MaxEntScan score3/score5, context++/PCT, pre-miRNA classifier, Drosha/Dicer, approximate TRF, α-satellite, HOR). Same protocol; all ⬜ pending first validation.

| Unit ID | Area | Stage A | Stage B | State | Report |
|---------|------|:---:|:---:|:---:|--------|
| ONCO-ASCAT-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-ASCAT-001.md |
| RNA-PKPREDICT-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | covered in reports/RNA-STRUCT-001.md (H-type csr-PK predictor) |
| RNA-PKRECURSIVE-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | covered in reports/RNA-STRUCT-001.md (recursive pknotsRG grammar) |
| RNA-ACCESS-001 | RnaStructure | ✅ | ✅ | ✅ CLEAN | reports/RNA-ACCESS-001.md |
| PROTMOTIF-HMM-001 | ProteinMotif | ✅ | ✅ | ✅ CLEAN | reports/PROTMOTIF-HMM-001.md |
| PRIMER-NNTM-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PRIMER-NNTM-001.md |
| PRIMER-HAIRPIN-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PRIMER-HAIRPIN-001.md |
| PRIMER-DIMER-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PRIMER-DIMER-001.md |
| PROBE-LNATM-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PROBE-LNATM-001.md |
| PROBE-EVALUE-001 | MolTools | ✅ | ✅ | ✅ CLEAN | reports/PROBE-EVALUE-001.md |
| MHC-NN-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/MHC-NN-001.md |
| MHC-MATRIX-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/MHC-MATRIX-001.md |
| IMMUNE-NUSVR-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/IMMUNE-NUSVR-001.md |
| META-CHECKM-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | reports/META-CHECKM-001.md |
| META-TETRA-001 | Metagenomics | ✅ | ✅ | ✅ CLEAN | reports/META-TETRA-001.md |
| SPLICE-MAXENT3-001 | Splicing | ✅ | ✅ | ✅ CLEAN | reports/SPLICE-MAXENT3-001.md |
| SPLICE-MAXENT5-001 | Splicing | ✅ | ✅ | ✅ CLEAN | reports/SPLICE-MAXENT5-001.md |
| MIRNA-CONTEXT-001 | MiRNA | ✅ | ✅ | ✅ CLEAN | reports/MIRNA-CONTEXT-001.md |
| MIRNA-PCT-001 | MiRNA | ✅ | ✅ | ✅ CLEAN | reports/MIRNA-PCT-001.md |
| MIRNA-CLASSIFY-001 | MiRNA | ✅ | ✅ | ✅ CLEAN | reports/MIRNA-CLASSIFY-001.md |
| MIRNA-CLEAVAGE-001 | MiRNA | ✅ | 🟡 | ✅ CLEAN | reports/MIRNA-CLEAVAGE-001.md |
| REP-APPROX-001 | Repeats | ✅ | ✅ | ✅ CLEAN | reports/REP-APPROX-001.md |
| CHROM-ALPHASAT-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-ALPHASAT-001.md |
| CHROM-HOR-001 | Chromosome | ✅ | ✅ | ✅ CLEAN | reports/CHROM-HOR-001.md |

---

# Validation Ledger — Phase 2 (Phase-2 Registry units)

Independent re-validation of the 148 Phase-2 ☑ Registry units absent from Phase 1, one fresh
session per unit, same protocol. Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL · ⬜ pending.
State: ✅ CLEAN · 🔧 LIMITED · ↩︎ DUPLICATE-OF.

**Progress:** 0 / 148 — RESET 2026-06-24, all units ⬜ pending re-validation.

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
| 1 | EPIGEN-METHYL-001 | Epigenetics | ⬜ | ⬜ | ⬜ pending | - |
| 2 | EPIGEN-DMR-001 | Epigenetics | ⬜ | ⬜ | ⬜ pending | - |
| 3 | VARIANT-CALL-001 | Variants | ✅ | ✅ | ✅ CLEAN | reports/VARIANT-CALL-001.md |
| 4 | VARIANT-SNP-001 | Variants | ⬜ | ⬜ | ⬜ pending | - |
| 5 | VARIANT-INDEL-001 | Variants | ⬜ | ⬜ | ⬜ pending | - |
| 6 | VARIANT-ANNOT-001 | Variants | ⬜ | ⬜ | ⬜ pending | - |
| 7 | SV-DETECT-001 | StructuralVar | ⬜ | ⬜ | ⬜ pending | - |
| 8 | SV-BREAKPOINT-001 | StructuralVar | ⬜ | ⬜ | ⬜ pending | - |
| 9 | SV-CNV-001 | StructuralVar | ⬜ | ⬜ | ⬜ pending | - |
| 10 | ASSEMBLY-OLC-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 11 | ASSEMBLY-DBG-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 12 | ASSEMBLY-STATS-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 13 | TRANS-EXPR-001 | Transcriptome | ⬜ | ⬜ | ⬜ pending | - |
| 14 | TRANS-DIFF-001 | Transcriptome | ⬜ | ⬜ | ⬜ pending | - |
| 15 | TRANS-SPLICE-001 | Transcriptome | ⬜ | ⬜ | ⬜ pending | - |
| 16 | COMPGEN-SYNTENY-001 | Comparative Genomics | ⬜ | ⬜ | ⬜ pending | - |
| 17 | COMPGEN-ORTHO-001 | Comparative Genomics | ⬜ | ⬜ | ⬜ pending | - |
| 18 | COMPGEN-REARR-001 | Comparative Genomics | ⬜ | ⬜ | ⬜ pending | - |
| 19 | PANGEN-CORE-001 | PanGenome | ⬜ | ⬜ | ⬜ pending | - |
| 20 | PANGEN-CLUSTER-001 | PanGenome | ⬜ | ⬜ | ⬜ pending | - |
| 21 | QUALITY-PHRED-001 | Quality | ⬜ | ⬜ | ⬜ pending | - |
| 22 | QUALITY-STATS-001 | Quality | ⬜ | ⬜ | ⬜ pending | - |
| 23 | SEQ-STATS-001 | Statistics | 🟡 | ✅ | ✅ CLEAN | reports/SEQ-STATS-001.md |
| 24 | SEQ-MW-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 25 | SEQ-PI-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 26 | SEQ-HYDRO-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 27 | SEQ-THERMO-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 28 | SEQ-DINUC-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 29 | SEQ-SECSTRUCT-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 30 | CODON-RSCU-001 | Codon | ⬜ | ⬜ | ⬜ pending | - |
| 31 | CODON-ENC-001 | Codon | ⬜ | ⬜ | ⬜ pending | - |
| 32 | CODON-STATS-001 | Codon | ⬜ | ⬜ | ⬜ pending | - |
| 33 | TRANS-SIXFRAME-001 | Translation | ⬜ | ⬜ | ⬜ pending | - |
| 34 | ASSEMBLY-MERGE-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 35 | ASSEMBLY-SCAFFOLD-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 36 | ASSEMBLY-COVER-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 37 | ASSEMBLY-CONSENSUS-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 38 | ASSEMBLY-TRIM-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 39 | ASSEMBLY-CORRECT-001 | Assembly | ⬜ | ⬜ | ⬜ pending | - |
| 40 | PAT-APPROX-003 | Matching | ⬜ | ⬜ | ⬜ pending | - |
| 41 | ALIGN-STATS-001 | Alignment | ⬜ | ⬜ | ⬜ pending | - |
| 42 | EPIGEN-BISULF-001 | Epigenetics | ⬜ | ⬜ | ⬜ pending | - |
| 43 | EPIGEN-CHROM-001 | Epigenetics | ⬜ | ⬜ | ⬜ pending | - |
| 44 | EPIGEN-AGE-001 | Epigenetics | ✅ | ✅ | ✅ CLEAN | reports/EPIGEN-AGE-001.md |
| 45 | MIRNA-PAIR-001 | MiRNA | ⬜ | ⬜ | ⬜ pending | - |
| 46 | PANGEN-HEAP-001 | PanGenome | ⬜ | ⬜ | ⬜ pending | - |
| 47 | PANGEN-MARKER-001 | PanGenome | ⬜ | ⬜ | ⬜ pending | - |
| 48 | POP-SELECT-001 | PopGen | ⬜ | ⬜ | ⬜ pending | - |
| 49 | POP-ANCESTRY-001 | PopGen | ⬜ | ⬜ | ⬜ pending | - |
| 50 | POP-ROH-001 | PopGen | ⬜ | ⬜ | ⬜ pending | - |
| 51 | META-FUNC-001 | Metagenomics | ⬜ | ⬜ | ⬜ pending | - |
| 52 | META-RESIST-001 | Metagenomics | ⬜ | ⬜ | ⬜ pending | - |
| 53 | META-PATHWAY-001 | Metagenomics | ⬜ | ⬜ | ⬜ pending | - |
| 54 | META-TAXA-001 | Metagenomics | ⬜ | ⬜ | ⬜ pending | - |
| 55 | PHYLO-BOOT-001 | Phylogenetics | ⬜ | ⬜ | ⬜ pending | - |
| 56 | PHYLO-STATS-001 | Phylogenetics | ⬜ | ⬜ | ⬜ pending | - |
| 57 | ANNOT-CODING-001 | Annotation | ⬜ | ⬜ | ⬜ pending | - |
| 58 | ANNOT-REPEAT-001 | Annotation | ⬜ | ⬜ | ⬜ pending | - |
| 59 | ANNOT-CODONUSAGE-001 | Annotation | ⬜ | ⬜ | ⬜ pending | - |
| 60 | RESTR-FILTER-001 | MolTools | ⬜ | ⬜ | ⬜ pending | - |
| 61 | KMER-DIST-001 | K-mer | ⬜ | ⬜ | ⬜ pending | - |
| 62 | MOTIF-CONS-001 | Matching | ⬜ | ⬜ | ⬜ pending | - |
| 63 | GENOMIC-REPEAT-001 | Analysis | ⬜ | ⬜ | ⬜ pending | - |
| 64 | GENOMIC-COMMON-001 | Analysis | ⬜ | ⬜ | ⬜ pending | - |
| 65 | GENOMIC-MOTIFS-001 | Analysis | ⬜ | ⬜ | ⬜ pending | - |
| 66 | SEQ-RNACOMP-001 | Composition | ⬜ | ⬜ | ⬜ pending | - |
| 67 | PROTMOTIF-PATTERN-001 | ProteinMotif | ⬜ | ⬜ | ⬜ pending | - |
| 68 | PROTMOTIF-SP-001 | ProteinMotif | ⬜ | ⬜ | ⬜ pending | - |
| 69 | PROTMOTIF-TM-001 | ProteinMotif | ⬜ | ⬜ | ⬜ pending | - |
| 70 | PROTMOTIF-CC-001 | ProteinMotif | ⬜ | ⬜ | ⬜ pending | - |
| 71 | PROTMOTIF-LC-001 | ProteinMotif | ⬜ | ⬜ | ⬜ pending | - |
| 72 | PROTMOTIF-COMMON-001 | ProteinMotif | ⬜ | ⬜ | ⬜ pending | - |
| 73 | RNA-PAIR-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 74 | RNA-HAIRPIN-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 75 | RNA-MFE-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 76 | RNA-PSEUDOKNOT-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 77 | RNA-DOTBRACKET-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 78 | RNA-INVERT-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 79 | RNA-PARTITION-001 | RnaStructure | ⬜ | ⬜ | ⬜ pending | - |
| 80 | SEQ-COMPLEX-KMER-001 | Complexity | ⬜ | ⬜ | ⬜ pending | - |
| 81 | SEQ-COMPLEX-WINDOW-001 | Complexity | ⬜ | ⬜ | ⬜ pending | - |
| 82 | SEQ-COMPLEX-DUST-001 | Complexity | ⬜ | ⬜ | ⬜ pending | - |
| 83 | SEQ-COMPLEX-COMPRESS-001 | Complexity | ⬜ | ⬜ | ⬜ pending | - |
| 84 | COMPGEN-RBH-001 | Comparative | ⬜ | ⬜ | ⬜ pending | - |
| 85 | COMPGEN-COMPARE-001 | Comparative | ⬜ | ⬜ | ⬜ pending | - |
| 86 | COMPGEN-REVERSAL-001 | Comparative | ⬜ | ⬜ | ⬜ pending | - |
| 87 | COMPGEN-CLUSTER-001 | Comparative | ⬜ | ⬜ | ⬜ pending | - |
| 88 | COMPGEN-ANI-001 | Comparative | ✅ | ✅ | ✅ CLEAN | reports/COMPGEN-ANI-001.md |
| 89 | COMPGEN-DOTPLOT-001 | Comparative | ⬜ | ⬜ | ⬜ pending | - |
| 90 | MOTIF-DISCOVER-001 | Matching | ⬜ | ⬜ | ⬜ pending | - |
| 91 | MOTIF-SHARED-001 | Matching | ⬜ | ⬜ | ⬜ pending | - |
| 92 | MOTIF-REGULATORY-001 | Matching | ⬜ | ⬜ | ⬜ pending | - |
| 93 | MOTIF-GENERATE-001 | Matching | ⬜ | ⬜ | ⬜ pending | - |
| 94 | KMER-ASYNC-001 | K-mer Analysis | ⬜ | ⬜ | ⬜ pending | - |
| 95 | KMER-UNIQUE-001 | K-mer | ⬜ | ⬜ | ⬜ pending | - |
| 96 | KMER-GENERATE-001 | K-mer | ⬜ | ⬜ | ⬜ pending | - |
| 97 | KMER-BOTH-001 | K-mer | ⬜ | ⬜ | ⬜ pending | - |
| 98 | KMER-STATS-001 | K-mer | ⬜ | ⬜ | ⬜ pending | - |
| 99 | KMER-POSITIONS-001 | K-mer | ⬜ | ⬜ | ⬜ pending | - |
| 100 | SEQ-ATSKEW-001 | Composition | ⬜ | ⬜ | ⬜ pending | - |
| 101 | SEQ-REPLICATION-001 | Composition | ⬜ | ⬜ | ⬜ pending | - |
| 102 | SEQ-GC-ANALYSIS-001 | Composition | ⬜ | ⬜ | ⬜ pending | - |
| 103 | DISORDER-MORF-001 | ProteinPred | ⬜ | ⬜ | ⬜ pending | - |
| 104 | DISORDER-PROPENSITY-001 | ProteinPred | ⬜ | ⬜ | ⬜ pending | - |
| 105 | DISORDER-LC-001 | ProteinPred | ⬜ | ⬜ | ⬜ pending | - |
| 106 | SEQ-COMPOSITION-001 | Statistics | — | — | ↩︎ DUPLICATE-OF SEQ-STATS-001 | - |
| 107 | SEQ-TM-001 | Statistics | ⬜ | ⬜ | ⬜ pending (was DUPLICATE-OF SEQ-THERMO-001) | - |
| 108 | SEQ-ENTROPY-PROFILE-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 109 | SEQ-GC-PROFILE-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 110 | SEQ-CODON-FREQ-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 111 | SEQ-SUMMARY-001 | Statistics | ⬜ | ⬜ | ⬜ pending | - |
| 112 | GENOMIC-TANDEM-001 | Analysis | ⬜ | ⬜ | ⬜ pending (was DUPLICATE-OF REP-TANDEM-001) | - |
| 113 | GENOMIC-SIMILARITY-001 | Analysis | ⬜ | ⬜ | ⬜ pending | - |
| 114 | GENOMIC-ORF-001 | Analysis | ⬜ | ⬜ | ⬜ pending | - |
| 115 | ONCO-SOMATIC-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 116 | ONCO-VAF-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 117 | ONCO-DRIVER-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 118 | ONCO-ARTIFACT-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 119 | ONCO-ANNOT-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 120 | ONCO-TMB-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 121 | ONCO-MSI-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 122 | ONCO-HRD-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-HRD-001.md |
| 123 | ONCO-LOH-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 124 | ONCO-SIG-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 125 | ONCO-SIG-002 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-SIG-002.md |
| 126 | ONCO-SIG-003 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-SIG-003.md |
| 127 | ONCO-SIG-004 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 128 | ONCO-FUSION-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 129 | ONCO-FUSION-002 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 130 | ONCO-FUSION-003 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 131 | ONCO-CNA-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 132 | ONCO-CNA-002 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 133 | ONCO-CNA-003 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 134 | ONCO-PURITY-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-PURITY-001.md |
| 135 | ONCO-PLOIDY-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-PLOIDY-001.md |
| 136 | ONCO-CLONAL-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 137 | ONCO-NEO-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 138 | ONCO-MHC-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-MHC-001.md (re-validated 2026-06-26 after the 8–11 → 8–14 widening, commit 66c24491. Stage A PASS: 8–14 = full NetMHCpan-4.1 class I window — Reynisson 2020 PMC7319546 "length range goes from 8 to 14 amino acids" — retrieved this session; MHCflurry `left_pad_centered_right_pad` 5–15 read from `mhcflurry` 2.1.5 source; prior 8–11 PASS-WITH-NOTES divergence resolved. Stage B PASS: lengths 12/13/14 IC50 reproduced vs the exact embedded net `PAN-CLASS1-1-3ed9fb2d2dcc9803` oracle to ≤ 9e-7 rel (25274.91 / 32389.13 / 32972.18 nM); 8/9/10-mer goldens byte-identical — `MhcflurryAffinityPredictor.cs` untouched by the fix; classification thresholds/chain unchanged. Full suite 18861/0. No defect; no change this session) |
| 139 | ONCO-CTDNA-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 140 | ONCO-MRD-001 | Oncology | 🟡 | ✅ | ✅ CLEAN | reports/ONCO-MRD-001.md |
| 141 | ONCO-CHIP-001 | Oncology | ✅ | ✅ | ✅ CLEAN | reports/ONCO-CHIP-001.md |
| 142 | ONCO-PHYLO-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 143 | ONCO-CCF-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 144 | ONCO-HETERO-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 145 | ONCO-HLA-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 146 | ONCO-ACTION-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 147 | ONCO-SV-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |
| 148 | ONCO-EXPR-001 | Oncology | ⬜ | ⬜ | ⬜ pending | - |

---

# Validation Ledger — Phase 3 (Independent re-validation of enhanced units)

After the open-questions program implemented 8 enhancements (Group 1 features + Group 2 breaking
changes), each touched unit is **independently re-validated in a fresh session** against external
first sources — the implementer and the validator are different contexts, per the protocol's premise
that code + tests authored together share blind spots. Scope = only the units whose code/tests changed.
Per-stage: ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL. State: ✅ CLEAN · 🔧 LIMITED.

**Progress:** 0 / 12 — RESET 2026-06-24, all enhancements ⬜ pending re-validation.

> Each enhanced unit was re-checked by an independent fresh validator against the **primary source**,
> under the rule **tests follow the source, code obeys the tests, tests are never bent to the code**.
> The independent pass found **1 real latent defect** the implementer + its tests had missed:
> PHYLO-NEWICK-001's parser silently accepted unbalanced parentheses (`((A,B);` → wrong topology) —
> now throws (650f18d). The published scoring tables were diffed element-by-element against freshly
> retrieved references (Doench-2014 70-coefficient table; MIT/Hsu 20-element W-vector + orientation;
> Kraken RTL/LCA/CQ re-derived in Python) — no transcription errors. Suite 6761 → 6780 (Failed: 0).

| # | Unit ID | Enhancement re-validated | Stage A | Stage B | State | Commit |
|---|---------|--------------------------|:---:|:---:|:---:|--------|
| 1 | ANNOT-GENE-001 | Reverse-strand Shine-Dalgarno / RBS reporting (`FindRibosomeBindingSitesBothStrands`; `forwardPos = len − revPos − motifLen` mapping) | ⬜ | ⬜ | ⬜ pending | - |
| 2 | RESTR-DIGEST-001 | Circular-molecule digest (`Digest(DnaSequence, MoleculeTopology, params string[])`; `MoleculeTopology{Linear,Circular}`; circular k sites → k fragments, origin-spanning join len `(len−lastCut)+firstCut`; 0→1 uncut, 1→1 linearized) | ⬜ | ⬜ | ⬜ pending | - |
| 3 | PARSE-EMBL-001 | INSDC location forms on `EmblParser.ParseLocation`: top-level remote ref `accession[.version]:descriptor` (`RemoteAccession`/`RemoteVersion`/`IsRemote`), remote ref NESTED in `complement(...)`/`join(...)` captured per-segment in `Location.RemoteParts` (`ExtractNestedRemoteReferences`), site-between `n^m` (`IsBetween`), deprecated single-period `n.m` (`IsSingleBaseFromRange`), `<..>` partials, `complement`/`join`/`order`, ordinary spans. Re-validated fresh 2026-06-25 verbatim to The DDBJ/ENA/GenBank FT Definition §3.4.2.1/§3.4.2.2/§3.4.3 (EBI `FT_current.txt`). All representative INSDC strings parse to hand-derived expected values with no version-digit leak (multi-digit `.12` & RefSeq `NC_000001.11` confirmed). `""`→`"` unescaping + shared GenBank path confirmed unaffected. N2: shared-helper splits `Parts` of `>`-partial spans into two single-base parts (Start/End/flag still correct) — pre-existing, out of scope. Boundary: remote entry sequence not fetched (needs network/DB) — documented, acceptable. **Enhancement 2026-06-26:** opt-in remote-aware FULL assembly `FeatureLocationHelper.ResolveLocationSequence` / `EmblParser.ResolveLocationSequence` with a **caller-supplied `RemoteSequenceResolver`** (library does no network I/O; caller supplies the fetch, library does the assembly) — segment order under `join`/`order`, reverse-complement under `complement(...)` (`complement(join(a,b))`=`join(complement(b),complement(a))`), 1-based inclusive slicing, `<`/`>` partials; INSDC §3.4/§3.5 re-retrieved this session (insdc.org + ddbj.nig.ac.jp); 17 hand-derived RESV tests. Pending re-validation. | ⬜ | ⬜ | ⬜ pending | - |
| 4 | ALIGN-MULTI-001 | Guide-tree **progressive** MSA (`MultipleAlignProgressive`) added beside the unchanged star `MultipleAlign`: all-pairs NW identity distances (d = 1 − identical-cols/aln-len) → UPGMA guide tree (proportional size-weighted averaging, lowest-index tie-break) → profile-profile NW over columns with sum-of-pairs profile scoring + "once a gap, always a gap". Sources this session: Feng & Doolittle 1987 (PubMed 3118049); Wikipedia "Multiple sequence alignment" §Progressive alignment (WebFetched — UPGMA/NJ guide tree, most-similar-first, "once aligned its alignment is not considered further"); Wikipedia UPGMA (WebFetched — proportional formula `d=(|A|·dAX+|B|·dBX)/(|A|+|B|)`, smallest-distance merge). Hand-derived end-to-end & re-confirmed by independent probe: ["ACGT","ACGT","AGT"] → dist d(0,1)=0/d(·,2)=0.25, merge (0,1) then add seq2 → rows ACGT/ACGT/**A-GT**, SP 8. Discriminating ["AAGAA","AACAA","GGTGG","GGTGG"]: progressive = gap-free len-5 SP −12 vs star gapped len-6 SP −13 (genuinely differ; progressive is the optimal). Profile-level once-a-gap independently derived ["ACGT","ACGT","AGT","AGT"] → ACGT/ACGT/A-GT/A-GT SP 15. Added strict M12 (profile-level once-a-gap, exact columns); mutation check (fill inserted gap column) kills 4 tests incl. M12. Star method byte-for-byte unchanged (S01). | ⬜ | ⬜ | ⬜ pending | - |
| 5 | CRISPR-GUIDE-001 | **Doench 2014 "Rule Set 1" on-target model** (`CalculateOnTargetDoench2014`, commit 129c2ca). Re-grounded **independently this session** by **re-downloading the reference `doenchScore.py`** (CRISPOR, github.com/maximilianh/crisporWebsite) to `/tmp` and reading it directly — NOT trusting the repo's copied constants. **30-mer layout** confirmed (4 up + 20 protospacer + 3 PAM + 3 down; GC over `seq[4:24]`, features by raw 0-based offset). **Intercept 0.59763615, gcLow −0.2026259 (≤10), gcHigh −0.1665878 (>10), GC term `abs(10−gc)·gcWeight`** all match. **Full tuple-by-tuple diff of the 70-entry table: EXACT MATCH (ordered), ∅ in ref-not-repo, ∅ in repo-not-ref — NO transcription error** (incl. the reference's own `(24,AG/CG/TG)`=`(24,A/C/T)` and `(26,GT)`=`(27,T)` quirks). Spot-checked weights vs reference file: (22,T)−0.8770074, (23,C)−0.8762358, (11,GG)−1.5169074, (20,GG)−0.7822076, (21,TC)−1.029693, (28,GG)−0.69774, (29,G)0.38634258, (6,C)−0.7411813, (17,G)−0.6780964. **Worked examples reproduced from scratch in Python from the coefficients** (not from the repo test literals): `TATAGCTGCGATCTGAGGTAGGGAGGGACC`→0.713089346955 vs published 0.713089368437 (Δ2.2e-08); `TCCGCACCTGTCACGGTCGGGGCTTGGCGC`→0.018983843193 vs 0.0189838463593 (Δ3.2e-09). M-003 all-A independently reproduced 4.4338168085440035. **Test gate real:** intercept +0.001 → score shift ~0.020 ≫ tol 1e-4, so a wrong intercept/weight fails M-001/M-002 (mutation-confirmed). Edge handling (wrong-length guard, lowercase, non-ACGT throw, 0–100 scaling) correct; added NGG-PAM guard is a documented SpCas9 stricter-contract divergence, not a defect. No code/test change required. | ⬜ | ⬜ | ⬜ pending | - |
| 6 | CRISPR-OFF-001 | **MIT/Hsu 2013 off-target model** (`CalculateMitHitScore` single-hit + `CalculateMitSpecificityScore` aggregate, commit 129c2ca) + the uncommitted W-orientation guard test. Primary source re-retrieved **independently this session**: **re-downloaded the reference `crispor.py`** (CRISPOR, github.com/maximilianh/crisporWebsite) — read `calcHitScore`/`calcMitGuideScore` and `hitScoreM` directly; docstring attributes the scheme to **"Scores of single hits" on crispr.mit.edu/about, "# The Patrick Hsu weighting scheme"** (Hsu et al. Nat Biotechnol 31:827, 2013, PMID 23873081). **W vector diffed element-by-element vs `hitScoreM` — EXACT MATCH, all 20:** `[0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]`. **Orientation pinned from the source** (not the code): `calcHitScore` indexes `hitScoreM[pos]`, `pos` 0..19 over the protospacer 5'→3'; index 0 = PAM-distal (5'), index 19 = PAM-proximal (3'/seed) — the high-weight tail (max W[13]=0.851) is the PAM-proximal seed, consistent with Wikipedia "Off-target genome editing" (PAM-proximal seed less mismatch-tolerant). **Formula re-derived from scratch in Python from the source** (NOT read off the C# array): `score1=Π(1−W[i])`, `score2=1/(((19−meanInterMmDist)/19)·4+1)` (≥2mm else 1), `score3=1/nmm²` (≥1mm else 1), `×100`; aggregate `100/(100+Σ)·100`. Reproduced every expected test value: perfect→100; W[0]→100; pos5(0.395)→60.5; pos19(0.583)→41.7; mm{5,15}→score1 0.10406·score2 0.345454·score3 0.25·100=0.8987; agg[60.5]→62.30529595; agg[60.5,41.7]→49.45598417. **Orientation guard test verified source-correct & is a genuine reversal-detector:** pos13(W=0.851)→14.9 (reversed vector would give W[6]=0.317→68.3), pos0(W=0)→100 (reversed→W[19]=0.583→41.7) — both re-derived independently incl. the reversed-vector counterfactual; KEPT as-is and committed. No transcription/orientation defect; no green-washing; no code/test change beyond committing the in-tree test. Off-target fixture 39 tests all green; full unfiltered suite 6773 passed / 0 failed. | ⬜ | ⬜ | ⬜ pending | - |
| 8 | PHYLO-NEWICK-001 | **N-ary (multifurcating) Newick I/O** (`ParseNewick`/`ToNewick`, `PhyloNode.Children`, commit c4f0190). Grammar re-sourced independently (WebFetched Wikipedia "Newick format" BNF `BranchSet→Branch \| Branch "," BranchSet` ⇒ polytomy, `Internal→"(" BranchSet ")" Name`, `Length→":" number`; WebFetched Felsenstein/Olsen PHYLIP `newicktree.html` "Trees can be multifurcating at any level", ex. `(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);`). **Multifurcation PASS:** hand-derived + probe-confirmed `(A,B,C);`→3 leaf children round-trips byte-for-byte; nested `(A,B,(C,D,E));` (5 leaves) exact; `(A:0.1,B:0.2,C:0.3);`→exact lengths→`(A:0.1000,B:0.2000,C:0.3000);`; 7-way `(A..G)`→7 children; Felsenstein ex.→3 top + nested clade, 5 leaves in order. **No-child-dropped pinned** (`Children.Count==3` + each name + ordered leaf set; OLD code threw on >2, regression-to-drop would be silent data loss). **CODE DEFECT FIXED:** closing `)` was not required — `(A,B`, `(A,B;`, `((A,B);` silently accepted as `(A,B)`/degenerate `((A,B))`, violating `Internal→"(" BranchSet ")"`; trailing-garbage guard only caught surplus *after* root. Added mandatory-`)` throw + 4 strict tests (extra-`)`, unclosed `(` w/ & w/o `;`, missing inner `)`); **mutation-verified** revert fails 2/2. Binary regression set + lenient missing-`;` `(A,B)` unchanged. Full unfiltered suite **6776 green**; build 0 errors. | ⬜ | ⬜ | ⬜ pending | - |
| 7 | META-CLASS-001 | **Kraken k-mer / LCA classifier** (`TaxonomyTree`+`Lca`, `BuildKmerDatabase`, `ClassifyReads`, commit 9839268). Independently re-validated against the **primary source re-retrieved this session** — Wood & Salzberg 2014, *Genome Biology* 15:R46 (open-access **PMC4053813**, WebFetched), + Kraken 2 manual for C/Q (WebFetched). **All five rules verbatim-confirmed:** DB build = LCA-collapse of a shared k-mer's owning taxa; classification = node-weighted-by-k-mer-count classification tree, RTL path scored as **sum of node weights along root→leaf**, leaf of max path assigned; tie → **LCA of tied leaves**; no hits → root/unclassified; confidence **C/Q** = clade-of-call k-mers / non-ambiguous k-mers. **Every test value hand-derived from the rules in a from-scratch Python reference** (NOT read off the code) and matched: LCA primitive (100,101→20; 100,200→10; 20,100→20; 100,100→100; 5,100→1; set→10; path [100,20,10,4,3,2,1]); SingleSpecies→100/RTL4/C4/Q4; SplitWithinGenus→20 (tie-LCA)/RTL2; SplitAcrossGenera→10 (tie-LCA)/RTL2; **RtlAncestorWeight** {100×1,101×2,20×1}→101 because RTL(101)=w(101)+w(20)=3 > RTL(100)=2 (ancestor weight summed into both paths)/C2/Q4/0.5; no-hit/empty/ambiguous→root/Q correct; **RC canonical** read AGGTT↔DB AACC→100/Q2; DB SharedKmer AGCT→Lca(100,101)=20, GCTA→100, GCTC↔GAGC→101. **Two out-of-suite adversarial cases** ({20×5,100×1}→species 100/RTL6; {100×2,200×2,10×1}→family 10/RTL3) re-derived & code matched — confirms RTL sums ancestors, not deepest-hit. **Mutation-checked this session:** tie-break→first-leaf fails 2 tests, DB-build→first-wins fails 1; restored, tree clean. **No code-echo, no defect, no code/test change.** Full unfiltered suite **6773 passed / 0 failed**; net8 core lib + genomics test project build 0 errors (pre-existing net9 MCP unbuildable on net8 SDK, unrelated). | ⬜ | ⬜ | ⬜ pending | - |
| 10 | PHYLO-COMP-001 | **Unrooted-bipartition RF (commit 4846294) + N-ary refactor of RF/MRCA (commit c4f0190)** (`CalculateUnrootedRobinsonFoulds`/`CalculateNormalizedUnrootedRobinsonFoulds` beside rooted-clade `RobinsonFouldsDistance`; `FindMRCA`/`PatristicDistance` traverse `PhyloNode.Children`). Re-grounded independently this session: **Robinson & Foulds 1981** (RF = symmetric difference of non-trivial bipartitions/splits; edge-contraction removes exactly that split) + **WebFetched Wikipedia "Robinson–Foulds metric"** (verbatim A+B definition; normalization 0–1). **Every RF/MRCA value hand-derived from the bipartition/clade definition** (NOT off the code): (a) identical→rooted 0/unrooted 0; (b) single NNI `((B,C),(A,(D,E)))`vs`((A,B),(C,(D,E)))` splits {BC\|ADE,DE\|ABC}vs{AB\|CDE,DE\|ABC}→**unrooted RF 2**; (c) **ROOT-INVARIANCE** X=`((A,B),(C,(D,E)))` Y=`(((A,B),C),(D,E))` SAME unrooted tree → both splits {AB\|CDE,ABC\|DE}→**unrooted RF 0**, rooted clades {CDE} only-X + {ABC} only-Y→**rooted RF 2≠0**; (d) **COLLAPSE** binary `(((A,B),C),(D,E))` vs multifurcating `((A,B,C),(D,E))` → rooted clades lose {AB}→**rooted RF 1**, unrooted splits lose {AB\|CDE}→**unrooted RF 1** (= splits/clades lost); max 5-taxon vs `(((A,C),E),(B,D))` disjoint→**RF 4 = 2(n−3)**, normalized 4/4=1, NNI 2/4=0.5, n=3 denom 2n−6=0→0; (e) **MRCA over 3-child polytomy** `((A,B,C),(D,E))`: MRCA(A,B)=(B,C)=(A,C)=polytomy node, (A,D)=(C,E)=root, (A,Z)=null. Code probe matched all. Unrooted reads rooted binary as unrooted via smaller-side `CanonicalSplitKey` (root-child dup dedups in HashSet, root side=total excluded); N-ary traversals iterate `Children`. **HARD mutation gate (run this session):** rooted-clade masquerade (`GetBipartitions`→`GetClades`) kills 5 tests incl. URF-ROOT; `Children.Take(2)` on all N-ary traversals kills all 4 multifurcation tests (RF-MULTI-ROOTED, RF-MULTI-UNROOTED, MRCA-POLYTOMY, STATS-POLYTOMY) — no test passes if multifurcation handling is wrong, no green-washing; source restored byte-for-byte. **No defect, no code/test change.** Full unfiltered suite **6779 passed / 0 failed**; build 0 errors. | ⬜ | ⬜ | ⬜ pending | - |
| 9 | PHYLO-TREE-001 | **N-ary tree construction: UPGMA (binary) + Neighbor-Joining (final trifurcation)** (`BuildUPGMA`, `BuildNeighborJoining`, `PhyloNode.Children`, commit c4f0190). Re-grounded independently this session: **Saitou & Nei 1987** (NJ — Q-matrix join, branch-length system, **final-step trifurcation** of the last 3 OTUs at the unrooted centre, `δ_i=(d_ij+d_ik−d_jk)/2`) and **Sokal & Michener 1958 / UPGMA** (proportional size-weighted averaging, height=d/2, strictly binary ultrametric). **NJ worked example (Wikipedia 5-taxon additive matrix) hand-run from the formulas (NOT off the code):** join (a,b) Q=−50 δ(a)=2 δ(b)=3 → join (u1,c) Q=−28 δ(c)=4 δ(u1)=3 → centre {u2,d,e} trifurcation δ(u2)=2 δ(d)=2 δ(e)=1; manual patristic reconstruction reproduces all 10 input distances exactly (additive recovery). **Code probe matched byte-for-byte:** root **Children.Count==3**, Newick `(((a:2.0000,b:3.0000):3.0000,c:4.0000):2.0000,d:2.0000,e:1.0000);` round-trips identically with the 3-child centre preserved on re-parse. **UPGMA worked example (Wikipedia 5S-rRNA) hand-run:** merges (a,b)@17→e@22→(c,d)@28→@33; δ(a)=δ(b)=8.5,(u,v)=2.5,e=11,(c)=(d)=14,(v,r)=5.5,(w,r)=2.5; ultrametric tip=16.5; code Newick `(((a:8.5000,b:8.5000):2.5000,e:11.0000):5.5000,(c:14.0000,d:14.0000):2.5000);` strictly binary (4 internal 2-child nodes). **Test-quality gap (test-only):** exact trifurcation branch lengths + exact Newick strings were not locked. **Added 3 strict sourced tests** (NJ exact trifurcation δ + leaf set; NJ exact-Newick round trip; UPGMA exact-Newick + strictly-binary), all expected values hand-derived in Stage A. **Mutation-checked:** NJ `>3`→`>2` (collapse to bifurcation) kills 3 tests; UPGMA proportional→naive average kills 3 tests; reverted, **production source byte-for-byte unchanged**. No code defect, no green-washing. Full unfiltered suite **6779 passed / 0 failed**; build 0 errors. | ⬜ | ⬜ | ⬜ pending | - |
| 12 | PHYLO-BOOT-001 | **N-ary refactor of bootstrap clade collection + support counting** (`Bootstrap`, `GetClades`/`CollectClades` over `PhyloNode.Children`; NJ-bootstrap test updated for the trifurcating root, commit c4f0190). Concern = correct N-ary clade collection during support counting (all children visited; clades over multifurcating replicate trees correct). Stage A re-confirmed vs Felsenstein 1985 (OSTI 6044842: "resample … with replacement … same size"; "keep all species while sampling characters"), Lemoine 2018 (PMC6030568: "pseudo-alignments of the same length"; support = "proportion of pseudo-trees containing that branch"; binary per-replicate), Biopython `Bio.Phylo.Consensus` (`find_clades(terminal=False)`, count/size). `CollectClades` (`PhylogeneticAnalyzer.cs:1011`) `foreach`-iterates **every** child — no first-two shortcut. **N-ary check hand-derived & probe-confirmed:** 4-taxon two-group → NJ root `((A,B),C,D)` genuine 3-child trifurcation (`root.Children.Count==3`), only rooted clade `{A,B}` support 1.0; **6-taxon three-pair → `((A,B),(C,D),(E,F))` whose 3rd root child is internal `(E,F)`** → clades `{A|B,C|D,E|F}` all 1.0. **Support invariants verified:** proportion in [0,1]; quantized to k/replicates; **1.0 on every-replicate-agreement** (invariant distance matrix → one topology each replicate); resample preserves length + all taxa; deterministic per seed. **HARD mutation gate:** `Children.Take(2)` drops the 3rd child → 6-taxon loses `{E,F}` + emits spurious `{A,B,C,D}` (`{A|B,A|B|C|D,C|D}`); 4-taxon emits spurious `{A,B,C}`. **GAP:** M7 caught the drop only via an *extra* key, not a *lost* clade → **added M8** (6-taxon trifurcation; keys ≡ `{A|B,C|D,E|F}`, each 1.0, no spurious union) as a lost-clade guard. All expected values hand-derived from the definition, not code-echoed; no assertion weakened, no skip. Bootstrap fixture 13→14. **No code/description defect.** Full unfiltered suite **6780 passed / 0 failed**; build 0 errors. | ⬜ | ⬜ | ⬜ pending | - |
| 11 | PHYLO-STATS-001 | **N-ary refactor of tree statistics + patristic distance** (`GetLeaves`, `CalculateTreeLength`, `GetTreeDepth`, `PatristicDistance` via `FindMRCAInternal`/`DistanceToTaxon`, `PhyloNode.Children`, commit c4f0190). Concern = correct N-ary traversal (all children visited, not just the first two). Stage A re-grounded **independently this session** against re-fetched reference docs: Biopython `Bio.Phylo.BaseTree` (`total_branch_length`="sum of all the branch lengths"; `count_terminals`="number of terminal (leaf) nodes"; `distance()`="sum of the branch lengths between two targets"=patristic distance through MRCA) + DendroPy `Tree.length()` ("sum of edge lengths … None→0"); tree height = longest root→leaf path in **edges** (Wikipedia graph-theory/ADT: single node 0, empty −1), distinct from Biopython branch-length `depths()`. **Polytomy hand-derivation** on `((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6);` (3-child node P over A,B,C): leaves=**5**; total length=0.1+0.2+0.3+0.0+0.4+0.5+0.6=**2.1** (incl. 3rd child C); depth=**2**; patristic A→B (MRCA=P)=**0.3**; patristic **A→C (3rd-child guard, MRCA=P)=0.4**. Code probe matched all four methods on the polytomy and on the binary M1–M10/S1/S2/C1 cases. All four methods use `foreach (… in node.Children)` — no `Take(2)`/Left-Right-only shortcut survived the refactor. **HARD mutation gate (run this session):** `Children.Take(2)` on `CalculateTreeLength` + `DistanceToTaxon` made length→**1.8** (drops C's 0.3) and `PatristicDistance(A,C)`→**NaN** (3rd child never reached), failing `TreeStatistics_OverMultifurcatingNode_TraverseAllChildren` exactly — the polytomy test is a genuine first-two-children regression guard, not green-washing; source restored byte-for-byte (working tree clean). No code defect, no test change required. Full unfiltered suite **6779 passed / 0 failed**; build 0 errors. | ⬜ | ⬜ | ⬜ pending | - |
