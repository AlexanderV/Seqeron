---
type: source
title: "Evidence: SV-DETECT-001 (discordant read-pair / PEM structural-variant detection)"
tags: [validation, structural-variant]
doc_path: docs/Evidence/SV-DETECT-001-Evidence.md
sources:
  - docs/Evidence/SV-DETECT-001-Evidence.md
source_commit: e525e3116b0a1c220283ce93c3fa751af524a7ae
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SV-DETECT-001

The validation-evidence artifact for test unit **SV-DETECT-001** — **Structural Variant Detection
from Paired-End Mapping (PEM) signatures** (discordant read pairs). The **discordant-read-pair
member of the germline structural-variant (SV) family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in [[discordant-pair-sv-detection]]. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the discordant-pair / PEM signature paradigm, corroborated
  across a peer-reviewed review, the BreakDancer reference implementation + protocol, the SAM
  proper-pair convention, and DELLY/SVXplorer):**
  - **Medvedev, Stanciu & Brudno 2009** (Nat Methods 6(11s):S13–S20, rank 1 review) — the PEM
    **signature catalogue**: a mate pair spanning a **deletion** maps with span **larger** than the
    insert size; an **insertion** gives a **smaller** span; an **inversion** flips one mate's
    orientation (the "basic inversion" signature); a **linking signature** connects arbitrarily
    distant regions or **different chromosomes** (translocation). Methods differ by "the signatures
    they can detect and … the way they cluster or window these signatures."
  - **BreakDancer README** (Chen et al. 2009, Nat Methods 6:677–681; rank 3 impl doc) — **anomaly
    cutoff `-c`** in units of s.d. (default **3**): `upper = mean + c·sd`, `lower = mean − c·sd`, a
    pair is span-discordant if `insertSize < lower OR > upper`; **minimum support `-r`** = min read
    pairs for a connection (default **2**); anomalies flagged by "unexpected separation distances or
    orientation"; SV type codes **DEL / INS / INV / ITX (intra-chromosomal translocation) / CTX
    (inter-chromosomal translocation)**.
  - **BreakDancer protocol** (Fan, Chen et al. 2014, Curr Protoc Bioinformatics, PMC3661775, rank 1)
    — read pairs "independently classified into six types" from (1) separation distance +
    orientation, (2) the user threshold, (3) the empirical insert-size distribution; thresholds of
    **3 s.d. or 4 s.d.** confirm the s.d.-cutoff convention.
  - **cureffi.org / BWA proper-pair** (rank 4 secondary, citing SAM FLAG `0x02`) — **concordant =
    FR** (upstream mate `+`, downstream `−`, pointing inward); **RF, FF, RR are abnormal**; BWA sets
    the proper-pair flag **only for FR** ⇒ FR is the sole formal concordant orientation, **RF is
    everted / discordant**.
  - **DELLY (Rausch 2012) + SVXplorer (Kumar 2020)** (rank 1 impls) — orientation → SV type:
    **FR** outliers at the far end of the insert distribution = **deletion**; **RF** (everted,
    changed relative order but default strands) = **tandem duplication**; **FF/RR** (one read's
    orientation change) = **inversion**. DELLY does "integrated paired-end and split-read analysis."

- **Documented corner cases / failure modes:**
  - **Insertion > insert size is invisible** — the small-span insertion signature does not appear
    once the insertion exceeds the fragment length, and the span **does not recover the inserted
    sequence** (Medvedev 2009).
  - **Below-support clusters** — clusters with fewer than the minimum supporting read pairs (default
    2) are not reported (BreakDancer).

- **Datasets (deterministic synthetic PEM signatures, derived from the cited signature rules;
  BreakDancer cutoff mean=400, sd=50, c=3 ⇒ bounds `[250, 550]`):** concordant FR span 400 same chr
  → not discordant; FR span 5000 same chr → **Deletion** (> 550); FR span 100 same chr →
  **Insertion** (< 250); **FF** same chr → **Inversion**; **RF** same chr → **Duplication**;
  chr1 ≠ chr2 → **Translocation**.

- **Coverage recommendations (8 items):** MUST — span `> μ+c·σ` (same chr, FR) → Deletion; span
  `< μ−c·σ` → Insertion; FF/RR intra-chromosomal → Inversion; different-chromosome pair →
  Translocation; concordant FR within bounds → no SV; cluster `< min-support` → no SV, `≥ min-support`
  → one SV. SHOULD — cutoff bounds exact (span at `μ±c·σ` concordant, one unit beyond discordant).
  COULD — empty input → empty output.

## Deviations and assumptions

- **ASSUMPTION — inter-chromosomal precedence over orientation.** When mates map to different
  chromosomes, the translocation (CTX) signature is reported **regardless of relative orientation**;
  chromosome difference is evaluated **first**. Justified because the sources define inter-chromosomal
  mapping as a linking/translocation signature (Medvedev 2009; BreakDancer CTX) and define inversion
  (INV) only for **intra-chromosomal** flipped pairs — a flipped orientation across chromosomes is
  undefined as an inversion.

No source contradictions — the Medvedev review, BreakDancer README + protocol, the SAM/BWA
proper-pair convention, and DELLY/SVXplorer agree on the discordant-pair span/orientation signatures,
the s.d. span cutoff, the FR-concordant / RF-everted convention, and the signature-then-cluster with
minimum-support paradigm.
</content>
