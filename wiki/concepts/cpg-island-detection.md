---
type: concept
title: "CpG site & CpG-island detection (Gardiner-Garden & Frommer)"
tags: [epigenetics, algorithm]
mcp_tools:
  - cpg_observed_expected
  - find_cpg_islands
  - find_cpg_sites
sources:
  - docs/Evidence/EPIGEN-CPG-001-Evidence.md
  - docs/algorithms/Epigenetics/CpG_Site_Detection.md
  - docs/Evidence/SEQ-DINUC-001-Evidence.md
source_commit: 4a7f3b50df393c2ccf0fe505da489d087d4f22f4
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: epigen-cpg-001-evidence
      evidence: "Test Unit ID: EPIGEN-CPG-001 ... Algorithm: CpG Site Detection ... Algorithm Group: Epigenetics"
      confidence: high
      status: current
---

# CpG site & CpG-island detection

Finding the **CpG dinucleotides** of a DNA sequence, scoring their **observed/expected density**, and
classifying **CpG islands** by the Gardiner-Garden & Frommer criteria. This is the **fourth ingested
unit of the Epigenetics family** and a **sequence-only** algorithm — distinct from every sibling:
[[bisulfite-methylation-calling]] *measures* the methylation state of these CpGs from reads,
[[epigenetic-age-horvath-clock]] scores age from already-measured β-values, and
[[chromatin-state-prediction]] works on histone-mark ChIP-seq signals. This unit touches **no
methylation state at all** — it is pure primary-sequence composition. Validated under test unit
**EPIGEN-CPG-001**; the record is [[epigen-cpg-001-evidence]], [[test-unit-registry]] tracks the
unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

A **CpG site** is a cytosine immediately followed by a guanine in the 5'→3' direction of one strand
(shorthand for 5'—C—phosphate—G—3'); it must not be confused with **GpC** (guanine-then-cytosine),
which is not a CpG. The unit exposes three deterministic operations on `EpigeneticsAnalyzer`, all
uppercasing the input before scanning:

1. `FindCpGSites(sequence)` — enumerate CpG positions.
2. `CalculateCpGObservedExpected(sequence)` — the canonical CpG O/E ratio.
3. `FindCpGIslands(sequence, minLength, minGc, minCpGRatio)` — sliding-window island detection.

## 1. CpG site enumeration

A single linear scan over adjacent dinucleotide windows; whenever the pair is `CG`, yield the
**0-based position of the C**. Exact, deterministic, `O(n)`. Adjacent CpGs are counted independently:
`CGCG` yields **two** sites (at 0 and 2) — each CpG is a distinct 2-nucleotide window, not overlapping.
A length-<2 input yields none. `FindCpGSites` is reused by [[bisulfite-methylation-calling]] to locate
the reference CpGs it calls methylation on (a single fixed two-character pattern in one pass is already
O(n)-optimal, so the repository suffix tree is deliberately not used).

## 2. CpG observed/expected ratio (Gardiner-Garden & Frommer 1987)

The canonical CpG-density statistic — how many CpGs are seen versus how many random base composition
would predict:

```
O/E = CpG_count / ((C_count × G_count) / L)
```

where `L` is sequence length. This is the formula used by the UCSC Genome Browser and most
bioinformatics tools. **Division-by-zero guard:** when the expected count is 0 (i.e. no C **or** no G,
or `L < 2`, or null/empty), the ratio returns **0.0**. A sequence with no CpG dinucleotides gives 0.

Worked oracles: `CGCGCG…` (20 bp, C=G=10) → Expected (10·10)/20 = 5.0, O/E = **2.0**; `ACGTCGACG`
(9 bp, C=G=3) → Expected 1.0, O/E = **3.0**; `ACGT` → Expected 0.25, O/E = **4.0**; AT-only → **0.0**.
Saxonov et al. (2006) give an **alternative** expected `((C+G)/2)²/L`; the repository uses the
Gardiner-Garden form. The CpG O/E ratio is a *dinucleotide* composition statistic — it is the
`CG`-specialized case of the general **Karlin dinucleotide odds ratio** `ρ_XY = f_XY/(f_X·f_Y)`
([[dinucleotide-relative-abundance]]); the two share the odds-ratio shape and differ only by the
dinucleotide-frequency normalization (this page divides the count by `L`/`N`, the general method by
`N − 1` — a documented modeling choice, `N/(N−1) → 1` for long sequences). Its single-base
compositional-asymmetry cousin is strand skew (`(A−T)/(A+T)`, `(G−C)/(G+C)`) —
see [[nucleotide-composition-skew]].

## 3. CpG-island detection (sliding window)

The default **Gardiner-Garden & Frommer** island criteria — a region qualifies when all three hold:

1. **length ≥ 200 bp** (`minLength`, default 200),
2. **GC content ≥ 50%** (`minGc`, default 0.5),
3. **CpG O/E ratio ≥ 0.6** (`minCpGRatio`, default 0.6).

A window of length `minLength` slides one base at a time; contiguous qualifying windows are **merged**
into one candidate island, then GC content and CpG ratio are **recomputed on the merged substring**
before it is yielded. Thresholds are compared **inclusively** (`gc >= minGc`, `cpgRatio >= minCpGRatio`)
— the `≥` operators are confirmed by Takai & Jones's explicit citation of the 1987 criteria. Islands
are `(Start, End, GcContent, CpGRatio)` tuples with **0-based inclusive `Start` / exclusive `End`** (the
stored boundary is `i + windowSize`). Oracle: 400 bp of `CGCG` repeats → one island (length ≥200, GC
100%, O/E 2.0). Complexity is `O(n × w)` (each window is rescanned, not a rolling count).

## Alternative criteria (not preset)

The stricter **Takai & Jones (2002)** criteria (**≥500 bp, GC% >55%, O/E >0.65**, which suppress
Alu-repeat false positives) and the **Saxonov (2006)** alternative expected formula are **not** offered
as named presets. A caller selects them by passing custom `minLength`/`minGc`/`minCpGRatio` arguments to
`FindCpGIslands`. This is an implementation scoping decision, not a spec deviation.

## Invariants and edge cases

- **INV:** every reported site marks a `C` whose next base is `G` (the defining CpG predicate);
  detection is deterministic for a given input.
- **INV:** the O/E ratio is uniquely determined by CpG/C/G counts and length; every reported island
  satisfies the configured length, GC, and O/E thresholds.
- **GpC ≠ CpG** — only the 5'→3' C-then-G order qualifies.
- **Case-insensitive:** mixed-case input is uppercased before scanning.
- Null / empty / length-<2 input → no sites, O/E `0.0`, no islands; a sequence shorter than
  `minLength` yields no islands; zero C or zero G → O/E `0.0`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for CpG enumeration, the O/E
ratio, and default-criteria island calling. **Sequence-only:** the island detector infers no
methylation state, chromatin context, or promoter status from the returned regions (those are the
sibling units' jobs). Not optimised with rolling counts (`O(n × w)` rather than incremental `O(n)`);
alternative literature criteria are supplied by argument, not preset. No source contradictions —
Gardiner-Garden & Frommer 1987, Takai & Jones 2002, Saxonov 2006, and Wikipedia are mutually
consistent; the Evidence file records **zero assumptions**.
