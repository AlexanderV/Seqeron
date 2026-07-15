---
type: source
title: "Evidence: SPLICE-DONOR-001 (donor / 5' splice site detection — GU/GT dinucleotide, MAG|GURAGU consensus, MaxEntScan score5ss)"
tags: [validation, splicing]
doc_path: docs/Evidence/SPLICE-DONOR-001-Evidence.md
sources:
  - docs/Evidence/SPLICE-DONOR-001-Evidence.md
source_commit: ce6f817f61151956d1e97909c1ccf5d70f0c333c
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SPLICE-DONOR-001

The validation-evidence artifact for test unit **SPLICE-DONOR-001** — **Donor (5') Splice
Site Detection**. The **donor member of the splicing family** (sibling of
SPLICE-ACCEPTOR-001) and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in [[splice-donor-site-prediction]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the canonical GU-AG 5'ss model):**
  - **Wikipedia RNA splicing / Spliceosome** (rank 4) — the near-invariant **GU** (GT in
    DNA) at the intron 5' end; donor consensus `G-G|GURAGU`; **GU-AG rule** (>99% of
    splicing); U1 snRNP binds the GU during E-complex formation; extended `MAG|GURAGU`
    (M = A/C), positions −3..+6; non-canonical minor-spliceosome variant.
  - **Shapiro & Senapathy 1987** (rank 1, foundational) — 5'ss position weight matrix,
    positions −3..+6: extended consensus `(C/A)AG|GU(A/G)AGU`; **position 0 (G) and +1 (U)
    ~100%**, −1 (G) ~80% (U1 base-pairing), −2 (A) ~60%, −3 (M) A/C ~35% each; log-odds
    scoring vs uniform 0.25 background.
  - **Burge, Tuschl & Sharp 1999** (rank 1) — **GC-AG** introns (~0.5–1% of U2-type use GC
    at the 5'ss); **U12-type** minor-spliceosome **AT-AC** introns (~0.3%), extended
    `/ATATCC/` donor motif; splice-site strength is extended-context-dependent.
  - **Yeo & Burge 2004** (rank 1) — **MaxEntScan score5ss**: max-entropy model over a
    **9-nt** window (3 exon + 6 intron), `log2(P_maxent/P_bgd)`. Factorisation + the
    16,384-record probability table fetched **verbatim** this session from the
    **MIT-licensed maxentpy port** (kepbod/maxentpy), embedded as `Data/maxent_score5.txt`.

- **MaxEntScan score5 factorisation (verbatim from maxentpy `score5`):** conserved `GT` at
  0-based positions 3..4 scored separately and removed; 7-nt "rest" = `window[0:3] +
  window[5:9]` looked up **directly in a single 4⁷ = 16,384-entry table** (contrast score3:
  overlapping sub-windows, 82,560 records). Score = `log2(GT_term · rest_term)`,
  `GT_term = cons1_5[G]·cons2_5[T]/(bgd_5[G]·bgd_5[T])`;
  `bgd_5 = {A:.27,C:.23,G:.23,T:.27}`, `cons1_5[G]=.9896`, `cons2_5[T]=.9884`.

- **Documented corner cases / failure modes:** non-GT (GC) donors valid U2-type but lower
  scoring; cryptic donor activation by point mutation; U12 AT donors rare but real;
  extended-context dependency; sequences shorter than the 9-nt window → empty; empty/null
  → empty (no throw); DNA T / RNA U + case equivalence.

- **Datasets (deterministic):**
  - **D1 canonical human β-globin (HBB, J00179) intron-1 donor** `...CAG|GTTGGT...`,
    extended −3..+6 `CAGGTTGGTG` → strong canonical GT donor.
  - **D2 relative scores** — `CAGGUAAGU` perfect consensus (highest) > `AAGGUAAAU`
    (weaker −3 A / +4 not G) > `UUUGUAAUU` (poor); `CAGGCAAGU` GC donor valid but lower.
  - **D3 negative controls** — `AAAAACCCCC` / `""` / short `GTAA` → empty.
  - **D4 MaxEntScan score5ss** — `cagGTAAGT` → **10.86** (10.858313, canonical primary
    cross-check), `gagGTAAGT` → **11.08** (11.078494), `taaATAAGT` → **−0.12**
    (−0.116791, a non-GT donor).

- **Coverage recommendations (16):** MUST — canonical GT detected; perfect `CAGGUAAGU`
  scores highest; no GT/GU → empty; strong > weak context; empty/short → empty; GC donor
  detected with `includeNonCanonical=true`; DNA T ↔ RNA U; case-insensitivity; the three
  `ScoreDonorMaxEnt` reference values (10.86 / 11.08 / −0.12) + strong-ranks-above-weak +
  wrong-length/non-ACGTU/null rejection. SHOULD — multiple GT sites all detected; score &
  confidence ∈ [0,1]; non-empty motif context. COULD — U12 AT donor with
  `includeNonCanonical`.

## Deviations and assumptions

All previously-listed assumptions are **ELIMINATED / RESOLVED design decisions**:

- **PWM values** replaced with **IUPAC consensus binary weights** (1.0 match / 0.0 no
  match) derived from the universally verified `MAG|GURAGU` consensus — no assumption.
- **Score normalization** replaced with a **simple consensus match fraction** (matches /
  positions scored) — no ad-hoc formula.
- **GC-donor 0.7 penalty removed** — GC donors naturally score lower because position +1
  (C) mismatches the invariant `U` consensus (max 8/9 ≈ 0.889 vs 9/9 = 1.0 for GT).
- **MaxEntScan** — opt-in `ScoreDonorMaxEnt` from the MIT-licensed maxentpy table;
  reproduces all three reference values exactly. **Licence flag (not buried):** the bundled
  table is the MIT port; the original Burge-lab Perl models carry academic terms — recorded
  in `Data/maxent_score5.LICENSE.md`. A maintainer seeking belt-and-suspenders commercial
  clearance should review the upstream Burge-lab terms.

No source contradictions — the encyclopedic, statistical (Shapiro/Senapathy), structural
(Burge), and max-entropy (Yeo/Burge) sources are mutually consistent.
