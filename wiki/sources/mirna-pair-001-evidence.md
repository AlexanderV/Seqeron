---
type: source
title: "Evidence: MIRNA-PAIR-001 (miRNA-target base pairing, Watson-Crick + G-U wobble)"
tags: [validation, mirna]
doc_path: docs/Evidence/MIRNA-PAIR-001-Evidence.md
sources:
  - docs/Evidence/MIRNA-PAIR-001-Evidence.md
source_commit: da06ef5590c4a92834125cfd650b942f61f9b586
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MIRNA-PAIR-001

The validation-evidence artifact for test unit **MIRNA-PAIR-001** — **MiRNA-Target Pairing
Analysis**, the base-pairing primitive that computes the antiparallel miRNA-mRNA duplex
(`MiRnaAnalyzer.AlignMiRnaToTarget` + the `CanPair` / `IsWobblePair` / `GetReverseComplement`
predicates). This is the **first ingested unit of the MiRNA family** and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing
concept is [[rna-base-pairing]]. [[test-unit-registry]] tracks the unit.

The unit validates the **shared RNA base-pairing rules** — Watson-Crick (A-U, G-C) plus the single
standard RNA **G-U wobble** — which the concept documents as a primitive reusable by RNA
secondary-structure work as well as miRNA targeting.

## What this file records

- **Online sources (all rank 1 except Wikipedia rank 4):**
  - **Bartel (2009)** *Cell* review — the miRNA **seed** = nucleotides **2-8** at the 5' end;
    silencing is directed "primarily through **Watson-Crick pairing** between the miRNA seed
    (nucleotides 2-7) and complementary sites within the 3'UTRs" of targets. (Cell full text 403'd
    the fetcher; facts taken from the PMC mirrors + search abstract.)
  - **PMC4870184** (transcriptome-wide non-canonical interactions review) — diagrams state "solid
    lines indicate **Watson-Crick base pairing**", explicitly distinguished from **wobble (G:U)**
    pairs; functional sites need "as few as 6-nt matches within the seed region (position 2-8)";
    "local short stretches (≥ 6 nt) of consecutive base-pairing significantly contribute to target
    recognition."
  - **Agarwal et al. (2015)** *eLife* (PMC4532895) — the canonical site-type ladder
    (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer) as Watson-Crick match patterns, the pairing
    chemistry "'A' pairs with 'U', and 'C' pairs with 'G'", and that the **A opposite position 1 is
    an Argonaute-recognised preference, independent of the miRNA nucleotide identity** — i.e. *not* a
    base-pairing event.
  - **Crick (1966)** "the wobble hypothesis" (via Wikipedia "Wobble base pair", rank 4 citing the
    primary *J. Mol. Biol.* 19(2):548-555) — standard pairs are A-U and G-C; the principal RNA
    **wobble pair is G-U** (other listed wobbles involve inosine, not standard RNA duplexes).
  - **Lewis et al. (2005)** *Cell* (PubMed 15652477) — targets are "mRNAs with conserved
    complementarity to the seed (nucleotides 2-7)", i.e. the target site is the **reverse complement
    of the miRNA seed read antiparallel**.
  - **NNDB Turner 2004** (Mathews/Turner; Xia et al. 1998) — the nearest-neighbor model: duplex free
    energy = Σ nearest-neighbor stacking over consecutive base pairs (+ initiation / terminal
    penalties); all Watson-Crick stacks at 37 °C are **negative (stabilising)**. Same `StackingEnergies`
    table reused from a prior committed unit (see Assumption).

- **Documented corner cases / failure modes:**
  - **G:U wobble is not Watson-Crick** — a duplex aligner must count wobbles separately from matches.
  - **A-opposite-position-1 is not a base-pairing event** — Argonaute-mediated, so not a "match" in a
    pure base-pairing aligner.
  - **Only G-U is a standard RNA wobble** — A-C, C-U, etc. do not pair and are mismatches.

- **Datasets (documented oracles):**
  - *Perfect Watson-Crick complement* — miRNA `AAAA` / target `UUUU` → 4 matches, 0 wobble,
    0 mismatch, alignment `||||`.
  - *G:U wobble duplex* — `GGGG` / `UUUU` → each G:U is a **wobble** → 0 Watson-Crick matches,
    4 wobbles, alignment `::::`.
  - *Non-pairing duplex* — `AAAA` / `AAAA` → A vs A cannot pair → 0 matches, 0 wobble, 4 mismatches.
  - *hsa-let-7a-5p seed reverse complement* — seed (pos 2-8) `GAGGUAG` →
    `GetReverseComplement` = `CUACCUC` (complement each base then reverse; trivially checkable from
    A-U/G-C).

- **Test-coverage recommendations:** MUST — `CanPair` true exactly for {A-U,U-A,G-C,C-G,G-U,U-G};
  `IsWobblePair` true only for G-U/U-G; `GetReverseComplement` antiparallel RNA reverse complement
  incl. the let-7a seed; `AlignMiRnaToTarget` perfect/wobble/non-pairing/empty duplex oracles.
  SHOULD — fully-paired duplex ΔG ≤ 0, all-mismatch not stabilising. COULD — DNA (T) normalised to
  RNA (U) before pairing.

## Deviations and assumptions

**One ASSUMPTION**, source-backed, not an open correctness gap:

- **Turner 2004 stacking numeric values not re-retrieved this session.** The exact per-stack
  free-energy values behind the ΔG estimate come from the NNDB Turner 2004 set already committed in a
  prior unit (`MiRnaAnalyzer.StackingEnergies`); the NNDB pages 404'd the fetcher this session, so the
  individual kcal/mol numbers were not independently re-opened. Consequence: `AlignMiRnaToTarget`
  tests assert only the **base-pairing structure** (match/mismatch/wobble counts, alignment symbols,
  antiparallel orientation) and the **sign** of the free energy (ΔG ≤ 0 for a fully paired duplex,
  ≥ 0 for an all-mismatch duplex), **not** exact ΔG magnitudes. The magnitude is documented
  "Intentionally simplified" (stacking-only; no loop/bulge/initiation/terminal-mismatch/AU-GU-end
  terms), so only its sign and relative ordering are reliable.

No source contradictions — Bartel 2009, PMC4870184, Agarwal 2015, Crick 1966, Lewis 2005, and Turner
2004 are mutually consistent on the pairing rules; the base-pairing classification is exact and
specification-driven while only the free-energy magnitude is a deliberate simplification.
