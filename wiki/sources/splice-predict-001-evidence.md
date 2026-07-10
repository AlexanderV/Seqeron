---
type: source
title: "Evidence: SPLICE-PREDICT-001 (gene structure prediction — intron/exon pairing, exon typing/phase, spliced sequence)"
tags: [validation, splicing]
doc_path: docs/Evidence/SPLICE-PREDICT-001-Evidence.md
sources:
  - docs/Evidence/SPLICE-PREDICT-001-Evidence.md
source_commit: ef8542055a460fa508da5902be1ed2256c1c6f83
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SPLICE-PREDICT-001

The validation-evidence artifact for test unit **SPLICE-PREDICT-001** — **Gene Structure
Prediction (Intron/Exon)**. The **composite / umbrella member of the splicing family**: it
assembles the two boundary detectors (SPLICE-DONOR-001 + SPLICE-ACCEPTOR-001) into a whole
split-gene model, and is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in [[gene-structure-prediction-intron-exon]]; [[test-unit-registry]] tracks the
unit.

## What this file records

- **Online sources (mutually consistent — the canonical GT-AG split-gene model):**
  - **Wikipedia Intron / Exon / Gene structure** (rank 4) — the **GT-AG rule** (5' donor
    `GT`/GU, 3' acceptor `AG`, citing Padgett 1986); intron length range **30 bp**
    (human *MST1L*, Piovesan 2015) to **>3.6 Mb** (*Drosophila DhDhc7*); four intron types
    (spliceosomal / tRNA / group I / group II); **branch point** near the 3' end forms the
    lariat; splicing ~99.999% accurate per intron (Hsu & Hertel 2009). Exon = any gene part
    surviving into mature RNA (Gilbert 1978); avg **5.48 exons/gene**, avg exon **~90–108 nt**
    (Sakharkar 2002); exon types **Initial / Internal / Terminal / Single**. Single-exon
    (intronless) genes exist (histones, most prokaryotic genes; Matera & Wang 2014).
  - **Breathnach & Chambon 1981** (rank 1) — **>99%** of introns begin `GT` end `AG`; donor
    consensus **`MAG|GURAGU`** (M=A/C, R=A/G); acceptor consensus **`(Y)nNCAG|G`**.
  - **Shapiro & Senapathy 1987** (rank 1) — the **PWM splice-site scoring** approach from
    >2,000 junctions (the scorer reused by the boundary units).
  - **Burge, Tuschl & Sharp 1999** (rank 1) — **U2** (GT-AG, ~99%) vs **U12** (AT-AC, ~0.5%)
    vs **GC-AG** variant (~0.5% of U2); and the pairing invariant: **each donor pairs with
    exactly one acceptor to define an intron**.

- **The composite pipeline:** find donors + acceptors → **pair** them into introns under
  GT-AG (bounded by `minIntronLength`/`maxIntronLength`, default `minIntronLength = 60`) →
  **greedy non-overlapping selection** by score → exons = the gaps between introns →
  **type** each exon (Initial/Internal/Terminal/Single) and compute **phase** =
  `(Σ preceding exon lengths) mod 3` → emit the **spliced sequence** (introns removed) and an
  **overall score** = mean of per-intron scores. Per-intron combined score =
  `(donor.Score + acceptor.Score) / 2`, normalized to **[0,1]**.

- **Documented corner cases / failure modes:** extremely short introns (30 bp) excluded by
  the default `minIntronLength=60`; overlapping donor/acceptor pairs need greedy
  non-overlapping selection; no splice sites found → single-exon structure; non-canonical
  `GC-AG` and `AT-AC` handled.

- **Datasets (deterministic):**
  - **D1 constructed two-exon gene** — Exon1 (35 nt) + donor `GUAAGU` + 60×A intron body +
    14-nt PPT `U…U` + acceptor `CAG` + Exon2 (35 nt) = **153 nt** → **1 intron, 2 exons**
    (Initial + Terminal), spliced length **70 nt**.
  - **D2 single-exon gene** (histone-like, 50 nt) → **0 introns, 1 exon** (Single),
    exon spans `[0, length−1]`.
  - **D3 empty/null** → 0 introns, 0 exons, empty spliced sequence, overall score **0**.

- **Coverage recommendations (16):** MUST — empty→empty structure; single-exon→Single;
  two-exon GT-AG→1 intron/2 exons; intron has GT donor + AG acceptor; `PredictIntrons`
  respects min/maxIntronLength; spliced sequence excludes intron; exon types assigned
  correctly; exon phase tracks frame; score ∈ [0,1]. SHOULD — non-overlapping selection;
  DNA T≡U; U2 (GT-AG) intron-type classification; short-sequence guard. COULD — multi-intron
  gene; case insensitivity.

## Deviations and assumptions

All four previously-listed assumptions are **RESOLVED**:

- **Greedy-by-score selection threshold** → reclassified a **design decision** (sources give
  only the "introns don't overlap" invariant, not the computational algorithm).
- **Exon phase calculation** → mathematically trivial (`Σ preceding lengths mod 3`; first
  exon = phase 0), standard convention per Alberts 2002.
- **Overall score = mean of intron scores** → reclassified a **definition** (no biological
  standard for a gene-structure quality metric exists).
- **Default branch-point score of 0.3** → **removed**; the combined score fell back to
  `(donor.Score + acceptor.Score) / 2`, eliminating an arbitrary magic constant.

No source contradictions — the encyclopedic (Wikipedia) and foundational (Breathnach &
Chambon, Shapiro & Senapathy, Burge et al.) sources are mutually consistent on the GT-AG
split-gene model.
