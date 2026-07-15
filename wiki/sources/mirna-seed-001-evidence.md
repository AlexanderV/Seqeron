---
type: source
title: "Evidence: MIRNA-SEED-001 (miRNA seed-sequence analysis — extraction + family equality)"
tags: [validation, mirna]
doc_path: docs/Evidence/MIRNA-SEED-001-Evidence.md
sources:
  - docs/Evidence/MIRNA-SEED-001-Evidence.md
source_commit: 989c8a14c92271602bfd8ee10f709009b73178b3
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MIRNA-SEED-001

The validation-evidence artifact for test unit **MIRNA-SEED-001** — **Seed Sequence Analysis**, the
string-level extraction and comparison of the mature-miRNA 5' seed
(`MiRnaAnalyzer.GetSeedSequence` / `CreateMiRna` / `CompareSeedRegions`). This is the **third
ingested unit of the MiRNA family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[seed-sequence-analysis]]. [[test-unit-registry]] tracks the unit.

The unit validates canonical **seed extraction** (positions 2-8, 7 nt), the normalised `MiRna`
record, and **exact-seed family equality** — the seed being the primary determinant of animal target
recognition and the input to the sibling pairing/target units.

## What this file records

- **Online sources (primary literature + reference DBs):**
  - **Wikipedia: MicroRNA** — seed = positions 2-7 of the miRNA, must be perfectly complementary for
    target recognition.
  - **TargetScan FAQ + 7mer docs** — canonical site classes over positions 2-8: **8mer** (2-8 + A
    opposite pos 1), **7mer-m8** (exact 2-8), **7mer-A1** (2-7 + A opposite pos 1), **6mer** (2-7).
  - **Lewis et al. (2005)** *Cell* (PMID 15652477) — "conserved seed pairing, often flanked by
    adenosines"; seed-based target prediction; target site = reverse complement of the seed.
  - **Bartel (2009)** *Cell* (PMID 19167326) — defines the seed region and the site hierarchy;
    **miRNA family = miRNAs with the same sequence at nucleotides 2-8**, sharing predicted targets.
  - **Grimson 2007 / Agarwal 2015 (TargetScan 7, PMID 26267216) / Friedman 2009** — context/PCT
    scoring beyond seed pairing (context of the site-type ladder; efficacy scoring is out of this
    unit's scope).
  - **miRBase** (Kozomara 2019) — authoritative miRNA sequences/nomenclature; source of the reference
    mature sequences and seeds.

- **Key definitions captured:** canonical seed = 2-7 (6 nt) vs extended seed = 2-8 (7 nt, adds m8);
  1-based indexing from the 5' end; family = identical 2-8 sequence.

- **Reference sequences / oracles (miRBase):** hsa-let-7a-5p / -7b / -7c-5p mature sequences all share
  seed **`GAGGUAG`** (pos 2-8) → **same family**; hsa-miR-21-5p → `AGCUUAU`; miR-155-5p → `UAAUGCU`;
  miR-1-3p → `GGAAUGU` (each a distinct family). `CreateMiRna("hsa-let-7a-5p",
  "UGAGGUAGUAGGUUGUAUAGUU")` → seed `GAGGUAG`; compared to let-7b → `IsSameFamily == true`.

- **Documented edge cases:**
  - miRNA shorter than 8 nt → cannot extract seed, returns `""` (miRNAs are ~23 nt; < 8 nt is invalid).
  - null/empty input → `""` (defensive).
  - DNA vs RNA: `T` handled equivalently to `U` **in `CreateMiRna`** (converts `T→U`); direct
    `GetSeedSequence` uppercases only, so a T-bearing seed can survive.
  - case-insensitivity: seed normalised to uppercase.
  - let-7 family → all members share the identical seed (validates family grouping).
  - self-comparison → 0 mismatches; completely different seeds → mismatches = seed length (7).

## Deviations and assumptions

**Intentionally simplified** (documented in the algorithm doc, not an open correctness gap):

- Family membership is reduced to **exact equality of the stored canonical 7-mer** (2-8). Consequence:
  noncanonical seed classes, shifted **isomiR** seeds, and broader/curated family definitions are not
  represented — an operational grouping inside the repository, not a complete biological taxonomy.
- Direct seed extraction returns only the canonical seed string; it interprets no pairing,
  conservation, or expression context and predicts no targeting strength.
- Coordinate note: the stored record uses **zero-based** `SeedStart = 1` / `SeedEnd = 7` even though
  the biological seed is described with 1-based positions 2-8. `CompareSeedRegions` also counts any
  seed-length difference as extra mismatches (canonical seeds are always 7 nt or empty).

**Terminology nuance flagged, not a contradiction:** the sources distinguish the canonical **2-7**
(6-nt) seed from the **2-8** (7-nt) extended seed, and the 8mer/7mer/6mer site ladder rests on that
distinction; this unit stores a single 7-nt (2-8) region and labels it "the canonical seed." Site-type
classification against a target is deferred to target-site prediction (MIRNA-TARGET-001). No source
contradictions — Bartel 2009, Lewis 2005, Agarwal 2015, TargetScan, and miRBase are mutually
consistent on the seed definition and the family concept.
