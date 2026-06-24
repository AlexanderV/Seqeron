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
| ONCO-MHC-001 | The **pan-allele NetMHCpan/MHCflurry neural** prediction, and a **bundled** trained coefficient matrix. (A matrix-based predictor is now implemented — opt-in `PredictBindingHalfLifeBimas` (BIMAS / Parker 1994 product rule) and `PredictIc50Smm` / `PredictAndClassifySmm` (SMM `IC50 = 50000^(1−score)`, Peters & Sette 2005) — chaining into the existing classifier; the published scoring rules are verified against the exact derived SMM anchors. The coefficient **matrix is caller-supplied** via `LoadScoringMatrix`.) | No redistributable, cross-verifiable trained HLA matrix was obtainable: the public BIMAS coefficient files are served only by a now-defunct dynamic CGI (unarchived), the Parker 1994 table is paywalled, and the IEDB SMM matrices carry a non-commercial / no-redistribution licence (like CIBERSORT LM22). The pan-allele neural model is a separate trained network. |
| SPLICE-ACCEPTOR-001 | The Yeo & Burge (2004) MaxEntScan **5' donor** `score5ss` maximum-entropy model. (The **3' acceptor** `score3ss` model is now bundled — opt-in `ScoreAcceptorMaxEnt`, with the precomputed probability tables embedded from the MIT-licensed maxentpy port; it reproduces the documented 2.89-bit reference exactly.) | The donor model needs the analogous Burge-lab 5' score tables (`me2x5` / `score5_matrix`) — a separate trained-data import not done this session. |
| ONCO-IMMUNE-001 | The CIBERSORT LM22 **signature matrix** (the ν-SVR solver + an LM22-format loader are now implemented — opt-in `DeconvoluteImmuneCellsNuSvr` / `LoadSignatureMatrix`, verified by planted-truth recovery and a scikit-learn/libsvm `NuSVR` cross-check). Bit-exact parity with the official CIBERSORT tool's per-sample fractions is not claimed (needs LM22 + the tool's quantile-normalisation/permutation pipeline). | LM22 is distributed by Stanford under a non-commercial, **no-redistribution** licence ("RECIPIENT shall not distribute the Program …") and is gated behind registration — so it is caller-supplied, not embedded. |
| META-BIN-001 | The **full CheckM lineage-specific marker database** (per-lineage Pfam/TIGRFAM collocated marker sets) + the **reference genome tree** for tree-based lineage placement. (CheckM-style completeness/contamination is **now implemented** — opt-in `EstimateBinQualityFromMarkers` / `EstimateBinQualityFromMarkerCounts`, the exact Parks 2015 per-collocated-marker-set formula, with marker detection via the Plan7 HMM engine; a small **CC0** universal single-copy ribosomal marker set (9 Pfam HMMs) is bundled and arbitrary marker HMMs are caller-supplied via `LoadMarkerHmms`. TNF/proxy metrics and `BinContigs` defaults unchanged.) | The `checkm_data` lineage-specific marker sets + reference tree are a large gated trained DB; lineage placement and operon-based marker-set collocation are not bundled — run CheckM itself, or supply lineage-specific marker HMMs. |
| PROTMOTIF-DOMAIN-001 | Exact HMMER `hmmsearch` bit-score / E-value parity, and Pfam coverage beyond the three bundled domains. (SH3 PF00018, PDZ PF00595 and WD40 PF00400 are **now detected** by the opt-in Plan7 profile-HMM engine `ProteinMotifFinder.FindDomainsByHmm` / `Plan7ProfileHmm`, scoring against bundled CC0 Pfam profiles; the exact-PROSITE-pattern `FindDomains` path is unchanged.) | The Plan7 Viterbi/Forward log-odds is verified exactly on a hand-built HMM and by correct true/false-positive ranking, but the full `hmmsearch` pipeline (MSV/bias filters, null2 biased-composition correction, Gumbel/exponential E-value calibration) and the full Pfam library are out of scope. |
| RNA-STRUCT-001 | Detecting pseudoknot classes outside the csr-PK grammar (kissing hairpins / loop–loop, triple-crossing / chained, non-canonical bulged helices). | The loop–loop interaction energy has never been experimentally measured (Sperschneider 2011 → heuristic estimate only); a faithful detector would require an unsourced energy constant. |
| DISORDER-REGION-001 | A calibrated per-residue / per-region disorder **confidence** value. | No disorder predictor publishes a calibrated confidence standard; the boundaries themselves follow the validated TOP-IDP threshold. |
| CHROM-CENT-001 | Suprachromosomal-family / specific α-satellite family (J1/J2/W/…) assignment for a detected HOR. (HOR *structure* — period, copy number, inter-/intra-HOR identity — is now detected by the opt-in `DetectHigherOrderRepeat`.) | Naming the family requires curated chromosome-specific reference HOR libraries (consensus HORs); these are external trained/curated data the library does not embed. |

## 3. Scope boundaries (deliberate — use the named reference tool)

| Unit | Not done | Use instead |
|------|----------|-------------|
| PRIMER-TM-001 | Hairpin / secondary-structure (folding-based) Tm. Internal single-mismatch and single dangling-end NN Tm are now implemented (opt-in `CalculateMeltingTemperatureNNMismatch`). | UNAFold, ViennaRNA, MELTING 5 for folding-based Tm. |
| REP-STR-001 | Whole-genome-scale seeded repeat discovery — TRF's probabilistic k-tuple seeding (the `R(d,k,pM)` 95% sum-of-heads percentile cut-off and the `W(d,pI)` random-walk band, whose values come from TRF's non-redistributable simulation tables). The per-repeat **Bernoulli statistical measures** (PM/PI, expected matches) are now computed by the opt-in `ComputeBernoulliStatistics`, and the deterministic exhaustive (start, period) scan already finds every candidate a seed would — so the residual is a performance index, not a correctness gap. | The reference Tandem Repeats Finder (for genome-scale seeded detection). |
| PROBE-DESIGN-001 | The **quantitative MGB (minor-groove binder) ΔTm**, and dual-quencher labelling. LNA-adjusted NN Tm is now implemented (opt-in `CalculateMeltingTemperatureNNLna`, McTigue 2004 NN increments), and the citable 3'-MGB *design rules* (3'-attachment, 12–20mer length window; Kutyavin 2000) are checked by `EvaluateMgbProbeDesign`. | The MGB ΔTm model (Kutyavin 2000 / MGB-Eclipse) is empirical/proprietary with no published closed-form formula — use a chemistry-specific design tool for an MGB Tm value; dual-quencher is a labelling note with no Tm impact. |
| PROBE-VALID-001 | Genome-scale **performance**: a precomputed seeded k-mer/seed index over a whole-genome database. The exhaustive sliding Smith–Waterman scan (`ScanOffTargetsGapped`) already finds every hit a seed would, and the off-target **Karlin–Altschul E-value / bit-score** statistics are now computed (opt-in `ComputeKarlinAltschul`), so the residual is a speed index, not a correctness gap. | BLAST / a genome-scale aligner (for whole-genome-scale speed only). |
| MIRNA-TARGET-001 | The remaining context++ features: `SA` (RNAplfold partition-function accessibility — not approximated by MFE), `PCT` (multi-species conservation), and `SPS` / `TA_3UTR` / `Len_ORF` / `ORF8m` (data-blocked: require the Garcia 2011 SPS table, a transcriptome, and the transcript ORF — accepted as optional `ContextPlusPlusInputs`). `3P_score`, `Min_dist`, `Len_3UTR`, and `Off6m` are now computed faithfully from the miRNA + 3'UTR; the score is still a partial context++. | The TargetScan reference pipeline for the headline context++ score (RNAplfold for SA, a multi-species alignment for PCT, the Garcia 2011 table for SPS, a transcriptome for TA, the ORF for Len_ORF/ORF8m). |
| MIRNA-PRECURSOR-001 | A **trained** natural-vs-background precursor classifier (miRDeep2-style fitted probabilistic model using read-stacking signatures) — data-blocked (needs a trained model + labelled data). Drosha/Dicer cleavage-site (mature/star coordinate) prediction is now implemented (opt-in `PredictDroshaDicerCleavage`, Han 2006 ~11 bp / Park 2011 ~22 nt 5'-counting rules). | miRDeep2 / a trained miRNA caller (for the probabilistic classifier only). |
| PARSE-EMBL-001 | Fetching the **sequence** of a remote entry referenced in a location (the reference itself is parsed; only local parts are extracted). | Retrieve the remote accession from the source database. |

---

## How to read this

- For **research / pipeline** use, every row is a normal scope boundary of a from-first-principles library — the algorithm is validated correct for its stated contract.
- For **clinical / decision-grade** use, the §1 (oncology) and §2 rows mark layers that require an external validated predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- **Irreducible** rows can never be closed; **data-blocked** rows reopen only if the gated / unmeasured data becomes available; **scope** rows point to the reference tool to use instead.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` (all `✅ CLEAN`).
