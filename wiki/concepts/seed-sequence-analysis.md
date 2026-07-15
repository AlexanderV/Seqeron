---
type: concept
title: "miRNA seed-sequence analysis (seed extraction + seed-family equality)"
tags: [mirna, algorithm]
mcp_tools:
  - compare_seed_regions
  - group_by_seed_family
  - mirna_seed_sequence
sources:
  - docs/Evidence/MIRNA-SEED-001-Evidence.md
  - docs/algorithms/MiRNA/Seed_Sequence_Analysis.md
source_commit: 989c8a14c92271602bfd8ee10f709009b73178b3
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: mirna-seed-001-evidence
      evidence: "Test Unit ID: MIRNA-SEED-001 ... Algorithm Group: MiRNA ... Seed Sequence Analysis (GetSeedSequence / CreateMiRna / CompareSeedRegions)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-base-pairing
      source: mirna-seed-001-evidence
      evidence: "the stored canonical seed used by downstream target-site prediction ... the target site is the reverse complement of the seed read antiparallel (Lewis 2005) — GetReverseComplement lives in rna-base-pairing"
      confidence: high
      status: current
---

# miRNA seed-sequence analysis (seed extraction + seed-family equality)

Extracting the **5' seed region** of a mature miRNA and comparing seeds to test **miRNA-family
membership**. This is the **third ingested unit of the MiRNA family** (test unit **MIRNA-SEED-001**,
`MiRnaAnalyzer`); the record is [[mirna-seed-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern. The seed is the primary
determinant of animal target recognition (Bartel 2009; Lewis et al. 2005), so this unit produces the
normalised seed string that the sibling pairing unit [[rna-base-pairing]] turns into a target motif
(`GetReverseComplement`) and that [[mirna-target-site-prediction|target-site prediction]]
(`FindTargetSites`, MIRNA-TARGET-001) **depends on** to seek the seed's reverse complement on the
mRNA. It is **string-level analysis, not binding prediction**.

The unit exposes three operations on `MiRnaAnalyzer`:

1. `GetSeedSequence(miRnaSequence)` — extract the 7-nt seed.
2. `CreateMiRna(name, sequence)` — normalise a sequence and materialise the `MiRna` record.
3. `CompareSeedRegions(miRna1, miRna2)` — count seed matches/mismatches and flag same-family.

## 1. Seed extraction — `GetSeedSequence`

The seed used here is **positions 2-8** of the mature miRNA (1-based, from the 5' end) — 7 nt,
i.e. the **extended seed** including the m8 position. Implementation: uppercase the input and return
`Substring(1, 7)` (zero-based indices 1..7). `O(1)` on a fixed-width window.

- **Casing only, no T→U:** `GetSeedSequence` uppercases but does **not** convert DNA `T` to RNA `U`
  — a DNA string passed here yields a seed that still contains `T`. Only `CreateMiRna` normalises `T→U`.
- **Short input:** `null`, empty, or length `< 8` → returns `""` (positions 2-8 do not exist), never throws.

**Terminology nuance (not a contradiction):** the literature distinguishes the **canonical seed**
= positions 2-7 (6 nt) from the **extended seed** = 2-8 (7 nt, adds m8), and the TargetScan site-type
ladder (**8mer** = 2-8 + A opposite pos 1; **7mer-m8** = 2-8; **7mer-A1** = 2-7 + A opposite pos 1;
**6mer** = 2-7) rests on that 6-vs-7 distinction (Agarwal et al. 2015; TargetScan). This unit stores a
single 7-nt (2-8) region and calls it "the canonical seed," collapsing that distinction — matching a
seed to a target and assigning a site class is deliberately **out of scope** and belongs to
[[mirna-target-site-prediction|target-site prediction]] (MIRNA-TARGET-001).

## 2. The `MiRna` record — `CreateMiRna`

`CreateMiRna(name, sequence)` normalises the sequence (`ToUpperInvariant()` then `T→U`), extracts the
seed from the **normalised** sequence, and stores `Name`, `Sequence`, `SeedSequence`, and fixed
`SeedStart = 1` / `SeedEnd = 7`. `O(n)` (dominated by normalisation).

- **Coordinate note:** `SeedStart`/`SeedEnd` are **zero-based inclusive** indices (1..7) into
  `Sequence`, even though the biological description uses 1-based positions 2-8.
- `name` is stored verbatim (no validation); `sequence` is **not** null-guarded (expected non-null).

## 3. Seed comparison + family membership — `CompareSeedRegions`

Reads the two stored seed strings and returns a `SeedComparison`: `Matches` = count of equal
positions, `Mismatches` = count of unequal positions **plus any seed-length difference**, and
`IsSameFamily` = (seed1 == seed2) exactly. This is a **Hamming comparison** over the fixed-length seed.

- **miRNA family = identical seed 2-8:** members of one family share the same 2-8 sequence (Bartel
  2009), and share predicted targets. Oracle: **hsa-let-7a/-7b/-7c-5p** all have seed `GAGGUAG`
  → same family; **hsa-miR-21-5p** seed `AGCUUAU` differs → different family.
- **Self-comparison → 0 mismatches; completely different seeds → mismatches = seed length (7).**
- **Empty seed short-circuit:** if either stored seed is empty, the result is zeroed
  (`Matches = 0`, `Mismatches = 0`, `IsSameFamily = false`) rather than a partial comparison.

## Invariants and edge cases

- **INV:** `GetSeedSequence` returns `""` or a 7-char uppercase seed (`Substring(1,7)`).
- **INV:** `CreateMiRna(...).SeedSequence == GetSeedSequence(CreateMiRna(...).Sequence)` (seed computed
  from the normalised sequence); `SeedStart = 1`, `SeedEnd = 7`.
- **INV:** `IsSameFamily` ⟺ the two stored seed strings are exactly equal.
- **INV:** for two present canonical seeds, `Matches + Mismatches = 7`.
- `GetSeedSequence` uppercases only (T survives); `CreateMiRna` normalises `T→U`. Length `< 8` → `""`.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for **canonical seed extraction and
exact-seed family equality** — deterministic string operations, `O(1)`/`O(n)`/`O(k)`. **Intentionally
simplified:** family membership is exact equality of the stored 7-mer, so noncanonical seed classes,
offset/centered seeds, **isomiR-shifted seeds**, and broader curated family definitions are **not**
represented — it is an operational grouping inside this repository, not a complete biological family
taxonomy (use miRBase / TargetScan for curated families). Seed→target-motif reverse complement is
owned by [[rna-base-pairing]] (`GetReverseComplement`); site-type (8mer/7mer/6mer) classification and
binding-efficacy scoring belong to [[mirna-target-site-prediction|target-site prediction]] (MIRNA-TARGET-001). **No
source contradictions** — Bartel 2009, Lewis 2005, Agarwal 2015, TargetScan, and miRBase agree on the
seed definition and family concept; the only recorded item is the intentional exact-7-mer-equality
simplification, with the 2-7-vs-2-8 terminology collapse noted above.
