# Seqeron.Genomics — Current Limitations

**Library:** Seqeron.Genomics (mission-critical)   **Last reviewed:** 2026-06-26

This file lists **only what the library currently does NOT do** — genuine open limitations. Anything already
implemented (including opt-in alternatives, convention-parity modes, and resolved defects) is deliberately
**excluded**; those are recorded in the per-unit reports under `docs/Validation/reports/{UNIT}.md`, not here.

Every limitation below is `BY-DESIGN` and the unit is `✅ CLEAN` for its stated contract — these are honest scope
boundaries, never defects. They fall into three kinds:

- **Irreducible** — no algorithm or model can close it (physics / information theory).
- **Data-blocked** — needs a trained model / matrix / database that is gated, non-redistributable, or has never been measured / published.
- **Scope** — a deliberate out-of-scope boundary; use the named reference tool, or supply the input yourself.

## Runtime enforcement (`LimitationPolicy`)

`Seqeron.Genomics.Core.LimitationPolicy` has **three modes**, least → most permissive — `Strict` <
`Moderate` < `Permissive`. A guarded call throws `SeqeronLimitationException` (naming the limitation,
what it is related to, and how to obtain the result another way) when the effective mode is *more
restrictive* than that limitation's **minimum access mode**. The default is **`Moderate`**. Set
`LimitationPolicy.DefaultMode`, or scope a region with `using (LimitationPolicy.Use(mode)) { … }`.

- **`Strict`** — only the ideal **and complete** result: throws on every guarded branch below.
- **`Moderate`** (default) — throws on the **non-ideal-output** branches; **allows** the
  **correct-but-incomplete / narrower-contract** branches.
- **`Permissive`** — allows everything (historical best-effort).

**Minimum access mode per guarded unit** (single source of truth: `LimitationCatalog`):

| Min mode | Units (guarded branch) |
|----------|------------------------|
| **`Permissive`** (non-ideal output — blocked in Strict & Moderate) | PARSE-FASTQ-001 (encoding undetermined), CHROM-CENT-001 (`Sf1OrSf2Dimeric`), DISORDER-REGION-001 (uncalibrated confidence — use `PredictDisorderRegions` for the validated boundaries), MIRNA-TARGET-001 (partial context++), MIRNA-CLEAVAGE-001 (approximate 3p/star span) |
| **`Moderate`** (correct-but-incomplete — blocked only in Strict) | ONCO-MHC-001 (SMM/BIMAS matrix score), ONCO-IMMUNE-001 (ABIS/caller-matrix deconvolution & ESTIMATE purity), META-BIN-001 (domain-level CheckM), PROBE-DESIGN-001 (qualitative MGB rules) |

Rows below with **no runtime guard** (documented only): the irreducible **RNA-STRUCT-001** pair (the
result is exact for the stated NN-energy model / csr-PK grammar, and the shortfall is undetectable per
call) and **MIRNA-PRECURSOR-001** (read-stacking is not implemented — nothing is returned to gate).

---

## 1. Irreducible (cannot be closed by any implementation)

| Unit | Not done | Why irreducible |
|------|----------|-----------------|
| PARSE-FASTQ-001 | Auto-disambiguation of Phred+33 vs Phred+64 for a FASTQ in which **every** read is confined to the overlap range (ASCII 64–74). | Such input is genuinely ambiguous (information-theoretic): both encodings decode it without error. `DetectEncoding(IEnumerable<string>)` resolves the normal case from the rest of the file (a single character &lt; 64 proves Phred+33; a character &gt; 74 with none below 64 infers Phred+64) and flags only the wholly overlap-confined file `Ambiguous`, defaulting to Phred+33 — that residual case cannot be resolved by any method. |
| RNA-STRUCT-001 | Recovering tertiary-stabilised knots (e.g. BWYV / PDB 437D) as the MFE structure. | Not representable by *any* nearest-neighbour thermodynamic model — an energy-model floor, not an algorithm gap. |
| RNA-STRUCT-001 | Detecting pseudoknot classes outside the csr-PK grammar (kissing hairpins / loop–loop, triple-crossing / chained, non-canonical bulged helices). | The loop–loop interaction energy has never been experimentally measured (Sperschneider 2011 → heuristic estimate only); a faithful detector would need an unsourced energy constant. |

## 2. Data-blocked (needs a trained model / matrix / database that is gated, non-redistributable, or never measured)

| Unit | Not done | Blocking data |
|------|----------|---------------|
| ONCO-MHC-001 | The *latest proprietary* NetMHCpan-4.1 ANN and the MHCflurry presentation / antigen-processing models; and a ready-to-use **matrix (SMM/BIMAS)** scorer. | NetMHCpan-4.1 is the vendor's trained model — use that tool. No redistributable, cross-verifiable matrix is obtainable: the BIMAS server is a now-defunct CGI, the Parker 1994 table is paywalled, IEDB SMM matrices are non-commercial. |
| ONCO-IMMUNE-001 | The **CIBERSORT-LM22-identical** signature matrix (exact-CIBERSORT parity). | LM22 is distributed by Stanford under a **no-redistribution** licence, gated behind registration; only the exact matrix is blocked. The ESTIMATE absolute-purity transform is calibrated on Affymetrix / ABSOLUTE data (Yoshihara 2013) — apply with that platform caveat off-array. |
| META-BIN-001 | The **per-lineage-specific** CheckM marker refinement + the **reference genome tree** for lineage placement. | The gated `checkm_data` (large trained DB) and the TIGRFAM-defined markers (CC BY-SA 4.0, not redistributable here) are caller-supplied; only domain-level completeness/contamination ships. |
| CHROM-CENT-001 | Separating **SF1 from SF2** (both dimeric with the identical A→B box pattern) and the chromosome-specific HOR identity (J1/J2/D1/D2/W1–W5/…). | Needs the SF-resolved consensus-monomer library, published only in non-machine-retrievable supplements and in no-licence HMM repos (HumAS-HMMER / logsdon-lab — not CC0/redistributable). Callers holding an SF-resolved reference pass it to `AssignSuprachromosomalFamily`. |
| DISORDER-REGION-001 | A calibrated per-residue / per-region disorder **confidence** value. | No disorder predictor publishes a calibrated confidence standard (the region boundaries themselves follow the validated TOP-IDP threshold). |

## 3. Scope boundaries (deliberate — use the named reference tool, or supply the input)

| Unit | Not done | Use instead |
|------|----------|-------------|
| PROBE-DESIGN-001 | The quantitative **MGB (minor-groove binder) ΔTm**, and dual-quencher labelling. | The MGB ΔTm model (Kutyavin 2000 / MGB-Eclipse) is empirical/proprietary with no published closed form — use a chemistry-specific tool; dual-quencher has no Tm impact. |
| MIRNA-TARGET-001 | A full context++ score out-of-the-box: `TA_3UTR` needs a caller-supplied 3'UTR set (no default transcriptome is bundled), `PCT` needs a caller-supplied alignment + tree + the published per-family sigmoid parameters (TargetScan's compiled tables are citation-required), and `SPS` / `Len_ORF` / `ORF8m` are caller-supplied. | Supply a 3'UTR set (`ComputeTa3Utr` then derives `TA_3UTR`); supply the alignment + tree + sigmoid parameters for `PCT`; supply `SPS` / ORF features. |
| MIRNA-CLEAVAGE-001 | The exact miRBase-annotated **miRNA\*(3p) / star** boundaries on the opposite hairpin arm (the 5p Drosha/Dicer cut reproduces miRBase exactly; the 3p span is a linear 2-nt-3′-overhang offset). | miRBase mature boundaries encode the dominant **sequencing-read** cut sites, not a deterministic fold + fixed-overhang rule (folding the precursor does not recover them — see report). Supply the miRBase mature-3p coordinates (MIMAT), or small-RNA-seq read pileups. |
| MIRNA-PRECURSOR-001 | The **read-stacking** (small-RNA-seq pileup) read-support score of miRDeep2. | Only the *qualitative* model is published (Friedländer 2008/2012); the closed-form log-odds exists only in a gated 2008 supplement and the GPL miRDeep2 source, so it is not re-implementable from retrievable sources. Use miRDeep2 with the caller's reads for the read-support signal. |

---

## How to read this

- This file lists capability/scope limitations only; pure **performance** optimisations (e.g. a genome-scale seeded index where an exhaustive scan already returns the same result) are not limitations and are not listed.
- For **research / pipeline** use, every row is a normal scope boundary of a from-first-principles library — the algorithm is validated correct for its stated contract.
- For **clinical / decision-grade** use, the §2 oncology rows mark layers that require an external validated predictor and clinical sign-off — the library computes the rule, not the trained model behind it.
- **Irreducible** rows can never be closed; **data-blocked** rows reopen only when the gated / unmeasured data becomes available (several already accept it via a caller-supplied loader); **scope** rows point to the reference tool, or accept the input from the caller.
- Each row traces to its per-unit validation report under `docs/Validation/reports/` (all `✅ CLEAN`). A unit that re-validated `PASS-WITH-NOTES` (🟡) but is **not** listed here carries only a by-design convention or an already-resolved doc/code note — recorded in its report, not a capability gap.
