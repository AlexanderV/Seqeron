---
type: source
title: "Evidence: TRANS-SPLICE-001 (Alternative splicing — event classification + Percent-Spliced-In / PSI)"
tags: [validation, transcriptome]
doc_path: docs/Evidence/TRANS-SPLICE-001-Evidence.md
sources:
  - docs/Evidence/TRANS-SPLICE-001-Evidence.md
source_commit: 82e3e03992f6e370559efdde3124a4b870a57893
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: TRANS-SPLICE-001

The validation-evidence artifact for test unit **TRANS-SPLICE-001** — **RNA-seq alternative /
differential splicing**: classify alternative-splicing events into the canonical five classes and
estimate **Percent-Spliced-In (PSI / Ψ)** from inclusion vs exclusion read support. A
**Transcriptome / RNA-seq family** unit and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm is synthesized in its own
concept, [[alternative-splicing-psi]]. [[test-unit-registry]] tracks the unit. Impl
`TranscriptomeAnalyzer.CalculatePSI` / `DetectAlternativeSplicing` (`Seqeron.Genomics.Annotation`).

**Scope note (disambiguation):** this is RNA-seq **read-quantification of splicing** (PSI from reads),
NOT the genomic **splice-site motif predictors** [[splice-donor-site-prediction]] /
[[splice-acceptor-site-prediction]] / [[gene-structure-prediction-intron-exon]], which score `GT`/`AG`
sequence signals. Same biology (splicing), orthogonal problem (quantify usage vs predict a site).

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Wang et al. (2008), *Nature* 456(7221):470–476** (PMID 18978772, rank 1) — the **five canonical
    AS event classes**: exon skipping (**SE**), intron retention (**RI**), alternative 5′ (**A5SS**)
    and 3′ (**A3SS**) splice sites, and mutually exclusive exons (**MXE**); PSI is the fraction of
    transcripts including a cassette exon, from inclusion-vs-exclusion read density; AS is pervasive
    (92–94% of multi-exon human genes).
  - **"Challenges in estimating percent inclusion…" BMC Bioinformatics 13(Suppl 6):S11** (PMC3330053,
    rank 1) — core definition **Ψ = inclusion / (inclusion + exclusion)**; Gaussian model
    `μ̃ = γᵢ/(γᵢ+γₑ)` (normalized inclusion/exclusion junction expression); read classification —
    reads on the alternative exon or its junctions with adjacent constitutive exons = **inclusion**,
    reads on the junction *between* the two adjacent constitutive exons (skipping the exon) =
    **exclusion**.
  - **Shen et al. (2014) — rMATS, *PNAS* 111(51):E5593–E5601** (PMC4280593, rank 1/3) —
    **length-normalized PSI** `ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ)` (I/S = inclusion/skipping counts, lᵢ/lₛ =
    effective isoform lengths = "number of unique isoform-specific read positions"); binomial model
    `I | ψ ∼ Binomial(n=I+S, p = lᵢψ/(lᵢψ + lₛ(1−ψ)))`. rMATS output event codes **SE / A5SS / A3SS /
    MXE / RI** match the Wang five-class taxonomy.
  - **rMATS project page + SUPPA2 (Trincado et al. 2018, PMID 29571299)** (rank 3, reference impls) —
    SUPPA **Ψ = inclusion/(inclusion+exclusion)** (inclusion = upstream-exon↔exon-of-interest
    junctions; exclusion = upstream↔downstream junction skipping it). rMATS-turbo coordinate columns
    confirm the **A5SS/A3SS convention**: on the **+ strand**, A5SS forms differ at the **downstream
    (3′/END) boundary** of the upstream exon (alternative **donor** = 5′ splice site), A3SS differ at
    the **upstream (5′/START) boundary** of the downstream exon (alternative **acceptor** = 3′ splice
    site).

- **Documented corner cases / failure modes:** **0/0 (no reads)** → PSI undefined (rMATS/Gaussian add
  a pseudo-count; here → NaN); **S = 0** → PSI = 1 (fully included); **I = 0** → PSI = 0 (fully
  excluded); **length normalization changes the estimate** — unnormalized Ψ = I/(I+S) and rMATS
  normalized differ whenever lᵢ ≠ lₛ; **an AS event requires ≥2 isoforms** of one gene (a single
  isoform defines no event).

- **Datasets (documented oracles):**
  - *Unnormalized PSI* — I=80, S=20 → `80/(80+20)` = **0.80** (PMC3330053; SUPPA).
  - *Length-normalized PSI (rMATS)* — I=80, S=20, lᵢ=200, lₛ=100 → I/lᵢ=0.40, S/lₛ=0.20 →
    `0.40/(0.40+0.20)` = **0.6666…** (Shen 2014).
  - *Event classification (two isoforms of one gene)* — SkippedExon (middle exon present in A,
    absent in B); RetainedIntron (intron retained as exon body); AlternativeFivePrimeSS (shared 5′
    start, different 3′ end ⇒ alternative donor); AlternativeThreePrimeSS (shared 3′ end, different
    5′ start ⇒ alternative acceptor); MutuallyExclusiveExons (exactly one of two alternative middle
    exons used).

- **Test-coverage recommendations:** MUST — `CalculatePSI` (80,20) → 0.80; with lengths
  (80,20,200,100) → 0.6666…; PSI=1 when S=0, PSI=0 when I=0; **0/0 → NaN**; `DetectAlternativeSplicing`
  classifies SkippedExon + RetainedIntron/A3SS/A5SS/MXE. SHOULD — invariant **0 ≤ PSI ≤ 1**; <2
  isoforms → no event. COULD — identical isoforms → no event.

## Deviations and assumptions

- **ASSUMPTION (length normalization is opt-in):** both the unnormalized Ψ = I/(I+S) (Wang 2008 /
  PMC3330053 / SUPPA) and the rMATS length-normalized Ψ are authoritative. `CalculatePSI` defaults to
  the **unnormalized read-count ratio** and switches to the rMATS form only when both effective
  lengths are supplied (> 0) — an API-shape choice; both numerical behaviors are source-backed.
- **ASSUMPTION (forward strand / ascending coordinates):** the event classifier compares exon
  coordinate lists assuming exons ordered 5′→3′ on one strand. Wang 2008 defines events per gene;
  strand handling is an input-normalization concern outside the formula.

No source contradictions — Wang 2008, PMC3330053, Shen 2014 (rMATS), and SUPPA2 are mutually
consistent on the PSI definition and the five-class taxonomy (they differ only in whether PSI is
length-normalized, which the API exposes as an option). Research-grade, not for clinical use.
