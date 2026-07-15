---
type: source
title: "Evidence: MOTIF-GENERATE-001 (IUPAC-degenerate consensus generation)"
tags: [validation, motif]
doc_path: docs/Evidence/MOTIF-GENERATE-001-Evidence.md
sources:
  - docs/Evidence/MOTIF-GENERATE-001-Evidence.md
source_commit: d36905351108ae77101357e168b9823952ca6dec
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MOTIF-GENERATE-001

The validation-evidence artifact for test unit **MOTIF-GENERATE-001** —
**IUPAC-Degenerate Consensus Generation** (`MotifFinder.GenerateConsensus`). It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the algorithm, its threshold→IUPAC decision rule, the NC-IUB symbol table,
worked oracles, and corner cases are synthesized in [[iupac-degenerate-consensus]]. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Cornish-Bowden / NC-IUB 1984** (rank 2, the official IUPAC/IUB nomenclature) — the
    bijective set→symbol mapping over the 11 non-trivial base subsets: two-base R/Y/S/W/K/M,
    three-base B/D/H/V (each "not-X"), N = any base. "Each symbol stands for the specific set
    of bases listed."
  - **UCSC Genome Browser IUPAC table** (rank 5) and **Wikipedia "Nucleic acid notation"**
    (rank 4, citing NC-IUB 1984) — both reproduce the same table verbatim, corroborating the
    standard.
  - **Bioconductor DECIPHER `ConsensusSequence`** (rank 3, reference implementation) — the
    **threshold-consensus family**: "removes the least frequent characters at each position,
    so long as they represent less than `threshold` fraction of the sequences," then encodes
    the survivors with IUPAC degeneracy codes; "degeneracy codes are always used in cases where
    multiple characters are equally abundant." DECIPHER's own default threshold is 0.05 — the
    threshold value is a tunable per-tool parameter, **not** a universal constant.
- **Datasets (documented oracles):**
  - *NC-IUB set→symbol table*: all 15 non-empty subsets of {A,C,G,T} → their IUPAC symbol.
  - *Threshold-inclusion worked examples* (this implementation, `threshold = n × 0.25`, base
    included iff `count > threshold`): two-base columns → R/Y/S/W/K/M; three-base → B/D/H/V;
    `A,A,G,G,C` (n=5, θ=1.25) → C dropped → **R**; `A,A,A,G` (n=4, θ=1.0) → G at ≤θ dropped →
    **A**; `A,C,G,T` (n=4, θ=1.0) → none exceed θ → **fallback to most-frequent** → **A**
    (tie broken alphabetically).

## Deviations and assumptions

Three documented **assumptions** — the *threshold-consensus family* and "encode the surviving
base set with an IUPAC code" rule are authoritative (DECIPHER/NC-IUB); only the specific
constants and the no-pass corner are implementation choices, all named and tested:

1. **25 % inclusion threshold, strict `>`.** A base is in a column's code iff
   `count > total × 0.25` (named constant `threshold = total * 0.25`). The 25 % cut and the
   strict `>` boundary are this implementation's design choice (DECIPHER's default is 0.05;
   tools vary). Tests pin the boundary; other inputs are chosen so the inclusion decision is
   unambiguous and the verified *symbol* is dictated solely by the authoritative NC-IUB table.
2. **Fallback when no base passes the threshold** (e.g. four bases each at exactly 25 %) →
   the single most-frequent base, ties broken alphabetically. No authoritative spec defines this
   corner; it is a verified implementation contract.
3. **Input normalisation:** column length taken from the first sequence; case-insensitive over
   {A,C,G,T}; non-ACGT characters ignored in the counts; empty collection → `""`; null →
   `ArgumentNullException`.

Recommended coverage (MUST): every two-base set → R/Y/S/W/K/M and every three-base set →
B/D/H/V; unanimous column → the input base; strict-`>` 25 % boundary (exactly-25 % excluded);
no-pass fallback → most-frequent. **No source contradictions** — NC-IUB, UCSC, and Wikipedia
agree on the symbol table; DECIPHER supplies the threshold-family framing.
