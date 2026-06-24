# Seqeron.Genomics — Current Limitations

**Library:** Seqeron.Genomics (mission-critical)   **Last reviewed:** 2026-06-24

This file lists **only what the library currently does NOT do** — genuine open limitations. Anything that has
been implemented (including opt-in alternatives, convention-parity modes, and resolved defects) is deliberately
**excluded**; those are recorded in the per-unit reports under `docs/Validation/reports/{UNIT}.md`, not here.

Every limitation below is `BY-DESIGN` and the unit is `✅ CLEAN` for its stated contract — these are honest scope
boundaries, never defects. They fall into three kinds:

- **Irreducible** — no algorithm or model can close it (physics / information theory).
- **Data-blocked** — needs a trained model / parameter set / database that is gated, non-redistributable, or has never been experimentally measured / published.
- **Scope** — a deliberate out-of-scope boundary; use the named reference tool instead.

---

## 1. Irreducible (cannot be closed by any implementation)

| Unit | Not done | Why irreducible |
|------|----------|-----------------|
| PARSE-FASTQ-001 | Auto-disambiguation of Phred+33 vs Phred+64 when a read carries only overlap-range (ASCII 64–73) qualities. | The encodings overlap; the input is genuinely ambiguous (information-theoretic). Callers must decode with an explicit offset. |
| RNA-STRUCT-001 | Recovering tertiary-stabilised knots (e.g. BWYV / PDB 437D) as the MFE structure. | Not representable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. |

## 2. Data-blocked (needs a trained model / parameter set / database that is gated, non-redistributable, or never measured)

| Unit | Not done | Blocking data |
|------|----------|---------------|
| ONCO-MHC-001 | Predicting peptide–MHC binding affinity / %Rank (the library only classifies a caller-supplied IC50 / %Rank). | No redistributable, cross-verifiable NetMHCpan-style trained model. |
| SPLICE-ACCEPTOR-001 | The Yeo & Burge (2004) MaxEntScan maximum-entropy 3' acceptor model. | Needs the Burge-lab trained score tables (`me2x3acc*`) — large data files, no clean redistribution licence. |
| ONCO-IMMUNE-001 | The CIBERSORT LM22 / ν-SVR immune deconvolution. | LM22 is gated behind academic registration (not plaintext-retrievable) and ν-SVR depends on it. |
| META-BIN-001 | Bin completeness / contamination from single-copy marker genes (CheckM). | CheckM's lineage-specific Pfam/TIGRFAM marker sets + reference genome tree (`checkm_data`) are a large trained DB with no clean plaintext source. |
| PROTMOTIF-DOMAIN-001 | Detecting SH3 (PROSITE PS50002) and PDZ (PS50106) domains, and any full Pfam HMM profile. | These are weight-matrix profiles / trained HMMs with no deterministic pattern to reproduce. |
| RNA-STRUCT-001 | Detecting pseudoknot classes outside the csr-PK grammar (kissing hairpins / loop–loop, triple-crossing / chained, non-canonical bulged helices). | The loop–loop interaction energy has never been experimentally measured (Sperschneider 2011 → heuristic estimate only); a faithful detector would require an unsourced energy constant. |
| DISORDER-REGION-001 | A calibrated per-residue / per-region disorder **confidence** value. | No disorder predictor publishes a calibrated confidence standard; the boundaries themselves follow the validated TOP-IDP threshold. |

## 3. Scope boundaries (deliberate — use the named reference tool)

| Unit | Not done | Use instead |
|------|----------|-------------|
| PRIMER-TM-001 | Mismatch / dangling-end / hairpin (secondary-structure) Tm terms. | Biopython `Tm_NN`, MELTING 5, UNAFold. |
| REP-STR-001 | TRF probabilistic k-tuple seeding + sum-of-Bernoulli significance + whole-genome-scale repeat discovery. | The reference Tandem Repeats Finder. |
| CHROM-CENT-001 | Higher-order repeat (HOR) structure / suprachromosomal-family classification (detection is monomer-level). | Dedicated centromere/HOR tooling. |
| PROBE-DESIGN-001 | MGB (minor-groove binder), LNA, and dual-quencher probe chemistries. | Chemistry-specific design tools. |
| PROBE-VALID-001 | A seeded genome-scale BLAST k-mer index + E-value model for off-target search. | BLAST / a genome-scale aligner. |
| MIRNA-TARGET-001 | The full-transcript context++ features (`3P_score`, `SPS`, `TA_3UTR`, `Min_dist`, `SA`, `PCT`, `Len_3UTR`, `Len_ORF`, `ORF8m`, `Off6m`); the computed score is a local-feature partial. | The TargetScan reference pipeline for the headline context++ score. |
| MIRNA-PRECURSOR-001 | Drosha/Dicer cleavage-site (mature/star coordinate) prediction and a trained natural-vs-background classifier. | miRDeep2 / a trained miRNA caller. |
| PARSE-EMBL-001 | Fetching the **sequence** of a remote entry referenced in a location (the reference itself is parsed; only local parts are extracted). | Retrieve the remote accession from the source database. |

---

## How to read this

- For **research / pipeline** use, every row is a normal scope boundary of a from-first-principles library — the algorithm is validated correct for its stated contract.
- For **clinical / decision-grade** use, the §1 (oncology) and §2 rows mark layers that require an external validated predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- **Irreducible** rows can never be closed; **data-blocked** rows reopen only if the gated / unmeasured data becomes available; **scope** rows point to the reference tool to use instead.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` (all `✅ CLEAN`).
