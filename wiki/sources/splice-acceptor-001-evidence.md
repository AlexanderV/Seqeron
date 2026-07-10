---
type: source
title: "Evidence: SPLICE-ACCEPTOR-001 (acceptor / 3' splice site detection — AG + PPT + branch point, MaxEntScan score3ss)"
tags: [validation, splicing]
doc_path: docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md
sources:
  - docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md
source_commit: 0176ee87cbff0a49e2b03a1d4b2f5e0d02034ab3
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SPLICE-ACCEPTOR-001

The validation-evidence artifact for test unit **SPLICE-ACCEPTOR-001** — **Acceptor Site
Detection** (3' splice site prediction). The **first ingested member of the splicing
family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in [[splice-acceptor-site-prediction]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the canonical GU-AG 3'ss model with three
  cis-elements):**
  - **Wikipedia RNA splicing / Polypyrimidine tract / Spliceosome** (rank 4) — the
    **GU-AG rule** (>99% of introns), 3'ss consensus `Y-rich-N-C-A-G|G`, PPT definition
    (C/U-rich, 15–20 nt, ~5–40 nt before intron 3' end, U2AF65 binding), branch point
    18–50 nt upstream, and the U12 minor-spliceosome AT-AC variant.
  - **Shapiro & Senapathy 1987** (rank 1, 3,700+ introns) — 3'ss nucleotide frequencies:
    position −2 A and −1 G ~100% conserved, −3 strong C (~65–70%), upstream −5..−15
    pyrimidine-enriched (C+U ≈ 70–80%), +1 first-exonic G most frequent (~50%). Source of
    the **AcceptorPwm** weights.
  - **Burge, Tuschl & Sharp 1999** (rank 1) — 3'ss structure `(Yn)NYAG|G`; PPT drives
    U2AF65 recruitment; **GC-AG** introns (~1%) keep AG acceptor; **U12** uses AC (~0.4%).
  - **Yeo & Burge 2004** (rank 1 model; rank 3 impl) — **MaxEntScan score3ss**: 23-nt
    window (20 intronic + 3 exonic, AG at 0-based 18–19), max-entropy model
    `log2(P_maxent/P_bgd)`. Factorisation + 82,560-record probability tables fetched
    **verbatim** this session from the **MIT-licensed maxentpy port** (kepbod/maxentpy),
    embedded as `Data/maxent_score3.txt`.
  - **Gao et al. 2008** (rank 1, PRIMARY) — human branch-point consensus **`yUnAy`**
    (positions −3..+1, branch A at 0); conservation y@−3 79.0%, U@−2 74.6%, **A@0 92.3%**,
    y@+1 75.1% (n = 181); median location −26 nt (150/181 at −34..−21); PPT spans 4–24 nt
    downstream of the branch point.
  - **Mercer et al. 2015** (rank 1, corroborating) — 59,359 high-confidence human
    branchpoints, predominantly adenosine, ~19–35 nt from the 3'ss — supports the 18–40 nt
    search envelope.
  - **Hall & Padgett 1994 / Jackson 1991** — U12 **YCCAC** 3'ss consensus (basis for the
    opt-in non-canonical scoring, replacing the old fixed 0.6 constant).

- **Documented corner cases / failure modes:** non-AG (U12 AC) acceptors; weak/interrupted
  PPT scores lower and may be skipped; cryptic intronic AG decoys (PPT quality + position
  disambiguate); sequence **< 20 nt** → empty (guard); no AG in range → empty.

- **Datasets (deterministic):**
  - **D1 canonical CAG + strong PPT** `UUUUUUUUUUUUUUUUCAGGG` (21 nt, AG@16–17) → ≥ 1 site.
  - **D2 no AG** `U×21` → empty.
  - **D3 strong vs weak PPT** — continuous-U > purine-interrupted context score.
  - **D4 branch-point yUnAy** — score `matched/maxScore`, `max = 0.790+0.746+0.923+0.751 =
    3.210`: perfect `CUUAC` @25 nt → 1.0; no-A@0 → 0.712 (< 0.8 ⇒ not found); purine@−3
    `AUUAC` → 0.754; window edges: 18 nt found / 17 nt not / 40 nt found / 41 nt not.
  - **D5 MaxEntScan score3ss** — `ttccaaacgaacttttgtAGgga` → **2.89** (2.886773),
    `tgtctttttctgtgtggcAGtgg` → **8.19** (8.190965), `ttctctcttcagacttatAGcaa` → **−0.08**
    (−0.080278). The 2.89 value is the canonical maxentpy/MaxEntScan cross-check.

- **Coverage recommendations (16):** MUST — canonical AG detected; no-AG/empty/null → empty;
  < 20 nt → empty; strong PPT > weak; score & confidence ∈ [0,1]; DNA T↔U equivalence;
  case-insensitive; multiple AG discovered. SHOULD — U12 YCCAC detected with
  `includeNonCanonical=true`, excluded by default; position reported after AG (i+1); motif
  contains AG context. COULD — threshold filtering; AG at position < 15 not detected.

## Deviations and assumptions

All previously-listed assumptions are **RESOLVED design decisions**, not open items:

- **PWM weights verified** against Shapiro & Senapathy 1987 (−3 C=0.70, −2 A=1.00,
  −1 G=1.00, +1 G=0.50; upstream C+U=0.80) — no assumption.
- **Normalization** `(score/(count+1)+2)/4` is a monotonic [0,1] linear mapping — an
  implementation **design choice**, not a biological claim.
- **U12 scoring** uses the published **YCCAC** consensus (Hall & Padgett 1994 / Jackson
  1991), replacing the earlier fixed 0.6.
- **Branch-point detection** — all constants (18–40 nt envelope, `yUnAy` weights, 4–24 nt
  PPT window) are source-traceable to Gao 2008 + Mercer 2015.
- **MaxEntScan** — implemented as opt-in `ScoreAcceptorMaxEnt` from the MIT-licensed
  maxentpy tables; reproduces all three reference values exactly. **Licence flag (not
  buried):** the bundled table is the MIT port; the original Burge-lab Perl models carry
  academic terms — recorded in `Data/maxent_score3.LICENSE.md`. A maintainer seeking
  belt-and-suspenders commercial clearance should review the upstream Burge-lab terms.

No source contradictions — the encyclopedic, statistical (Shapiro/Senapathy), structural
(Burge), max-entropy (Yeo/Burge), and branch-point (Gao/Mercer) sources are mutually
consistent.
