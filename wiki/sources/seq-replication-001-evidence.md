---
type: source
title: "Evidence: SEQ-REPLICATION-001 (Replication origin prediction — cumulative GC-skew minimum)"
tags: [validation, sequence-statistics, composition, chromosome]
doc_path: docs/Evidence/SEQ-REPLICATION-001-Evidence.md
sources:
  - docs/Evidence/SEQ-REPLICATION-001-Evidence.md
source_commit: c094b65e4a89b3c3c146d655c12489e6d28e8564
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-REPLICATION-001

The validation-evidence artifact for test unit **SEQ-REPLICATION-001** — **replication origin
prediction via the cumulative GC-skew minimum** (the classic ori-finding method). It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; [[test-unit-registry]] tracks the unit.

**A genuinely distinct algorithm → its own concept.** The scalar `(G−C)/(G+C)` skew and its
windowed profile already live on [[nucleotide-composition-skew]] and
[[windowed-gc-profile-and-variance]]; those pages *anticipated* the origin locator as a future
application. This unit realises it — the **cumulative integer skew** plus **argmin/argmax
locating** — so the method is synthesized on the concept
[[replication-origin-cumulative-skew]]. This page records only what the artifact itself pins.

## What this file records

- **Online sources:**
  - **Rosalind — Minimum Skew Problem (BA1F)** (rank 3, canonical worked example) — defines
    `Skew(Genome)` as *"the difference between the total number of occurrences of 'G' and 'C'"*,
    computed as a **per-position cumulative sum** (G:+1, C:−1, A/T:0, `Skew_0 = 0`) over prefix
    indices `i ∈ [0, |Genome|]` (so `|Genome|+1` values). Verbatim 100-nt sample input yields
    **min skew −4 at positions 53, 97** — re-derived independently in-session, reproducing the
    published output exactly.
  - **Grigoriev A (1998)** *Nucleic Acids Res* 26(10):2286 (rank 1, primary) — the cumulative
    skew diagram reaches its **global minimum at the replication origin** and its **global
    maximum near the terminus**, the two extrema separated by **≈ half the chromosome length**;
    leading strand is G-rich.
  - **Lobry JR (1996)** *Mol Biol Evol* 13(5):660 (PMID 8676740, rank 1, primary) — GC/AT skews
    **switch sign at origin and terminus**; the founding strand-asymmetry observation.
  - **Wikipedia "GC skew"** (rank 4, citing Lobry 1996 + Grigoriev 1998) — *"the maximum value
    of the cumulative skew corresponds to the terminal, and the minimum value corresponds to
    the origin of replication."*
- **Datasets:**
  - **BA1F sample** — 100-nt genome, min skew value **−4** at positions **53, 97**;
    `PredictedOrigin = 53` (first minimizing prefix index, this unit's tie-break).
  - **Tiny derived diagrams:** `CCGGGG` → `0,−1,−2,−1,0,+1,+2` (min −2 @ 2, max +2 @ 6);
    `GGGCCC` → `0,+1,+2,+3,+2,+1,0` (min 0 @ 0, max +3 @ 3); `AATT` → `0,0,0,0,0` (flat).
- **Corner cases / failure modes:**
  1. **Ties** — BA1F has two minimizers (`53 97`); a single-position API must fix a
     deterministic tie-break (this unit returns the **first/smallest** minimizing index).
  2. **Prefix indexing** — `i ∈ [0, |Genome|]`, `|Genome|+1` prefix values, `Skew_0 = 0`
     before any base.
  3. **No asymmetry → no signal** — a sequence with no net G/C strand bias gives a **flat**
     diagram (amplitude 0); origin/terminus not meaningfully resolved.

## Deviations and assumptions

**One documented assumption — `IsSignificant` semantics.** No authoritative source defines a
numeric "significance" cutoff for an origin call. The previous implementation used an invented
threshold (`amplitude > count × 0.01`); that untraceable constant was **removed**.
`IsSignificant` is redefined as the threshold-free, evidence-neutral predicate **`max > min`**
(the diagram has non-zero amplitude ⇒ a detectable strand-composition asymmetry exists per
Lobry 1996 / Grigoriev 1998). This is the weakest non-invented definition; callers needing a
quantitative confidence should inspect the skew amplitude directly.

Recommended coverage (from the artifact): MUST — BA1F sample ⇒ `PredictedOrigin = 53`;
per-nucleotide skew G:+1/C:−1/A,T:0 with `Skew_0 = 0` (`CCGGGG`, `GGGCCC`); `PredictedTerminus`
= global-maximum position; first-occurrence tie-break. SHOULD — flat diagram (no/balanced G/C)
⇒ origin = terminus = 0, `IsSignificant = false`; case-insensitive, A/T do not move the
diagram. COULD — empty/null ⇒ zero prediction (string) / throws (`DnaSequence`). No source
contradictions — Rosalind BA1F, Grigoriev 1998, Lobry 1996, and Wikipedia agree on the min =
origin / max = terminus locating rule and the skew counting convention.
