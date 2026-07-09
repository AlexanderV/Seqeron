---
type: concept
title: "Bisulfite sequencing methylation calling"
tags: [epigenetics, algorithm]
sources:
  - docs/Evidence/EPIGEN-BISULF-001-Evidence.md
  - docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md
source_commit: 3fdfb015426652fb931b21c72bd92c5b0a214a7e
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: epigen-bisulf-001-evidence
      evidence: "Test Unit ID: EPIGEN-BISULF-001 ... Algorithm: Bisulfite Sequencing Analysis ... Algorithm Group: Epigenetics"
      confidence: high
      status: current
---

# Bisulfite sequencing methylation calling

Turning **sodium-bisulfite sequencing** chemistry into base-level **DNA-methylation calls**. This is
the **second ingested unit of the Epigenetics family** and a **distinct algorithm** from the sibling
[[epigenetic-age-horvath-clock]] (which scores age from *already-measured* β-values — this concept is
how those per-CpG β-values are *produced* from bisulfite reads). Validated under test unit
**EPIGEN-BISULF-001**; the record is [[epigen-bisulf-001-evidence]], [[test-unit-registry]] tracks the
unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

Sodium bisulfite converts **unmethylated cytosine → uracil (read as thymine)** while leaving
**5-methylcytosine unchanged**, so methylation becomes a **C↔T difference** between a read and the
reference (Frommer et al. 1992 — the protocol "yields a positive display of 5-methylcytosine"). The
unit bundles three deterministic, exact operations (no probabilistic modelling) on
`EpigeneticsAnalyzer`:

1. `SimulateBisulfiteConversion(sequence, methylatedPositions)` — in-silico conversion of one strand.
2. `CalculateMethylationFromBisulfite(referenceSequence, bisulfiteReads)` — per-CpG methylation caller.
3. `GenerateMethylationProfile(sites)` — per-context weighted aggregation.

## 1. In-silico bisulfite conversion (Frommer 1992)

Per-base substitution on a **single strand**, iterating position `i`:

- **unprotected cytosine → thymine** (uracil, read as T);
- **protected (5-methyl) cytosine → cytosine** (unchanged, `methylatedPositions` are 0-based);
- **any non-cytosine → unchanged** (A, G, T never react with bisulfite).

Output length equals input length; **case is preserved** (`c → t`). Only the **supplied strand** is
converted — the complementary strand is a separate molecule and is *not* synthesised or merged
(Frommer's protocol is strand-specific). Oracle: `ACGTCGAA` with methylated `{1}` → `ACGTTGAA` (C@1
protected stays `C`; C@4 unmethylated → `T`; non-C bases unchanged).

## 2. Per-CpG methylation calling (Bismark rule)

Reference CpG cytosines are enumerated by reusing `FindCpGSites` (one linear scan; the repository
suffix tree is deliberately not used — a single fixed two-character pattern in one pass is already
O(n)-optimal). For each aligned read (`(ReadSequence, StartPosition)`, 0-based start), at every covered
reference CpG C the read base is tallied by the **Bismark C=methylated / T=unmethylated rule** (Krueger
& Andrews 2011). The methylation **level** is the Bismark percentage expressed as a fraction in [0, 1]:

```
level = methylated / (methylated + unmethylated)
```

Each emitted `MethylationSite` carries `Position` (0-based CpG C), `Type` = CpG, `Context`, `level`,
and `Coverage` (total valid calls). Oracle (reference `ACGTACGT`, CpG sites at index 1 and 5): one `C`
read + one `T` read at index 1 → level **0.5**, coverage 2; a single `T` read at index 5 → level
**0.0**, coverage 1.

Bismark's per-context call symbols are `z/Z` (CpG), `x/X` (CHG), `h/H` (CHH) — lowercase unmethylated,
uppercase methylated; this unit reports the numeric level rather than the symbol but uses the same rule.

## 3. Per-context weighted profile (Schultz 2012)

Sites are split by context and aggregated into a `MethylationProfile` (global / CpG / CHG / CHH
weighted levels + total/methylated CpG counts + per-position levels). The **weighted methylation
level** is the read-pooled fraction, *not* the unweighted mean:

```
weighted level = Σ(methylated reads) / Σ(methylated + unmethylated reads)
              = Σ(levelᵢ · coverageᵢ) / Σ(coverageᵢ)
```

Worked oracle: sites (level 1.0, coverage 10) and (level 0.0, coverage 30) →
`(1.0·10 + 0.0·30)/(10+30) = 10/40 = ` **0.25** — deliberately different from the unweighted mean
`(1.0+0.0)/2 = 0.5`. When total coverage is zero (e.g. sequence-only sites with no reads), the profile
falls back to the **unweighted mean** of per-site levels so such sites are not silently dropped.

## Invariants and edge cases

- **INV:** conversion output length = input length; protected C stays C, unprotected C→T (case
  preserved); non-C bases unchanged.
- **INV:** every reported level ∈ [0, 1], coverage ≥ 1; per-context profile = Σ(level·coverage)/Σcoverage.
- **Zero-coverage CpG sites are excluded** from calling output — the percentage is undefined when the
  denominator is 0.
- A CpG needs **both** the C and the following G, so the last reference base cannot start a CpG.
- A read base at a reference-C that is **neither C nor T** (e.g. an A/G mismatch) is **not a valid
  bisulfite call** and is ignored; read bases **past the reference end** are ignored.
- Null/empty `sequence` or `referenceSequence` → empty results; an empty `sites` enumerable → an
  all-zero `MethylationProfile`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the core conversion and
calling formulas, **not** a production BS-seq read aligner. **Not implemented:** complementary-strand
conversion/merging; bisulfite non-conversion-rate correction; sequencing-error filtering; statistical
significance testing; **CHG/CHH calling from reads** (the caller targets CpG sites — the profile
aggregator still buckets CHG/CHH sites if supplied). One accepted API-contract deviation: the registry
signature `CalculateMethylationFromBisulfite(bsSeq, refSeq)` is realised as `(referenceSequence,
bisulfiteReads)` because per-site coverage needs per-read multiplicity that a single converted string
cannot carry. For full-genome BS-seq, use dedicated pipelines (Bismark, methylKit). No source
contradictions — Frommer 1992, Krueger & Andrews 2011, the Bismark User Guide, and Schultz 2012 are
mutually consistent.
