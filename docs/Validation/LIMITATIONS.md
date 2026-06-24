# Seqeron.Genomics — Current Limitations

**Library:** Seqeron.Genomics (mission-critical)   **Last reviewed:** 2026-06-25

This file lists **only what the library currently does NOT do** — genuine open limitations. Anything that has
been implemented (including opt-in alternatives, convention-parity modes, and resolved defects) is deliberately
**excluded**; those are recorded in the per-unit reports under `docs/Validation/reports/{UNIT}.md`, not here.

Every limitation below is `BY-DESIGN` and the unit is `✅ CLEAN` for its stated contract — these are honest scope
boundaries, never defects. They fall into three kinds:

- **Irreducible** — no algorithm or model can close it (physics / information theory).
- **Data-blocked** — needs a trained model / parameter set / database that is gated, non-redistributable, or has never been experimentally measured / published. (Where the *algorithm* is implemented, the gated data is **caller-supplied** via a loader — the limitation is that nothing usable ships out-of-the-box.)
- **Scope** — a deliberate out-of-scope boundary; use the named reference tool instead.

---

## 1. Irreducible (cannot be closed by any implementation)

| Unit | Not done | Why irreducible |
|------|----------|-----------------|
| PARSE-FASTQ-001 | Auto-disambiguation of Phred+33 vs Phred+64 when a read carries only overlap-range (ASCII 64–73) qualities. | The encodings overlap; the input is genuinely ambiguous (information-theoretic). Callers must decode with an explicit offset. |
| RNA-STRUCT-001 | Recovering tertiary-stabilised knots (e.g. BWYV / PDB 437D) as the MFE structure. | Not representable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. |
| RNA-STRUCT-001 | Detecting pseudoknot classes outside the csr-PK grammar (kissing hairpins / loop–loop, triple-crossing / chained, non-canonical bulged helices). | The loop–loop interaction energy has never been experimentally measured (Sperschneider 2011 → heuristic estimate only); a faithful detector would require an unsourced energy constant. |

## 2. Data-blocked (needs a trained model / matrix / database that is gated, non-redistributable, or never measured)

| Unit | Not done | Blocking data |
|------|----------|---------------|
| ONCO-MHC-001 | A **bundled** trained HLA coefficient matrix (the predictor matrix is caller-supplied), and the pan-allele NetMHCpan/MHCflurry neural model. | No redistributable, cross-verifiable trained HLA matrix is obtainable: the BIMAS coefficient files are served only by a now-defunct CGI (unarchived), the Parker 1994 table is paywalled, and the IEDB SMM matrices are non-commercial / no-redistribution. The pan-allele neural model is a separate trained network. |
| ONCO-IMMUNE-001 | A **bundled** CIBERSORT LM22 signature matrix (LM22 is caller-supplied). | LM22 is distributed by Stanford under a non-commercial **no-redistribution** licence and gated behind registration. |
| META-BIN-001 | The **full CheckM lineage-specific marker database** (per-lineage Pfam/TIGRFAM collocated sets) + the reference genome tree for lineage placement. | `checkm_data` is a large gated trained DB; only a small CC0 universal ribosomal marker set ships, the rest is caller-supplied. |
| PROTMOTIF-DOMAIN-001 | Exact HMMER `hmmsearch`-*reported* E-value **pipeline** parity, and Pfam coverage beyond the three bundled (CC0) / caller-supplied `.hmm` domains. | The Gumbel (MSV/Viterbi) and exponential-tail (Forward) P-value and `E = P·Z` formulas are now implemented from the profile's `STATS LOCAL` calibration (opt-in). What remains is the rest of the `hmmsearch` pipeline that HMMER applies before those formulas — the **null2 biased-composition correction** and the MSV/bias prefilters on a local-multihit bit score (this scorer is glocal) — plus the full Pfam library beyond the bundled profiles. |
| CHROM-CENT-001 | Suprachromosomal-family / specific α-satellite family (J1/J2/W/…) assignment for a detected HOR. | Naming the family needs curated chromosome-specific reference HOR (consensus) libraries — external curated data not embedded. |
| DISORDER-REGION-001 | A calibrated per-residue / per-region disorder **confidence** value. | No disorder predictor publishes a calibrated confidence standard (the region boundaries themselves follow the validated TOP-IDP threshold). |

## 3. Scope boundaries (deliberate — use the named reference tool)

| Unit | Not done | Use instead |
|------|----------|-------------|
| PRIMER-TM-001 | Hairpin / secondary-structure (folding-based) Tm. | UNAFold, ViennaRNA, MELTING 5. |
| PROBE-DESIGN-001 | The quantitative **MGB (minor-groove binder) ΔTm**, and dual-quencher labelling. | The MGB ΔTm model (Kutyavin 2000 / MGB-Eclipse) is empirical/proprietary with no published closed form — use a chemistry-specific tool; dual-quencher has no Tm impact. |
| MIRNA-TARGET-001 | The context++ feature `PCT` (multi-species conservation); `TA_3UTR` / `SPS` / `Len_ORF` / `ORF8m` are caller-supplied. The computed score is a partial context++. (`SA` is now computed — see note.) | The TargetScan reference pipeline for the headline context++ score (a multi-species alignment for PCT, a transcriptome for TA). `SA` now uses the library's own Turner-2004 McCaskill partition-function accessibility (`RnaSecondaryStructure.CalculateRegionUnpairedProbability`), computed for every site whose 14-nt window fits the 3'UTR. |
| MIRNA-PRECURSOR-001 | A **trained** natural-vs-background precursor classifier (read-stacking probabilistic model). | miRDeep2 / a trained miRNA caller. |
| PARSE-EMBL-001 | Fetching the **sequence** of a remote entry referenced in a location (the reference itself is parsed; only local parts are extracted). | Retrieve the remote accession from the source database. |

---

## How to read this

- This file lists capability/scope limitations only; pure **performance** optimisations (e.g. a genome-scale seeded index where an exhaustive scan already returns the same result) are not limitations and are not listed.
- For **research / pipeline** use, every row is a normal scope boundary of a from-first-principles library — the algorithm is validated correct for its stated contract.
- For **clinical / decision-grade** use, the §2 oncology rows mark layers that require an external validated predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- **Irreducible** rows can never be closed; **data-blocked** rows reopen only when the gated / unmeasured data becomes available (several already accept it via a caller-supplied loader); **scope** rows point to the reference tool to use instead.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` (all `✅ CLEAN`).
