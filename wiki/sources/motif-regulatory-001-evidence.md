---
type: source
title: "Evidence: MOTIF-REGULATORY-001 (Regulatory-element consensus catalog scan)"
tags: [validation, motif]
doc_path: docs/Evidence/MOTIF-REGULATORY-001-Evidence.md
sources:
  - docs/Evidence/MOTIF-REGULATORY-001-Evidence.md
source_commit: 914ab57b7357e92a2f33c107af9b6fdabe63ac45
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MOTIF-REGULATORY-001

The validation-evidence artifact for test unit **MOTIF-REGULATORY-001** — scanning a DNA
sequence against a **curated catalog of 12 canonical regulatory consensus elements** and
reporting each occurrence with its name, pattern, and 0-based start position. The
algorithm, catalog, matching contract, corner cases, and the AP-1 corrected-consensus
defect are synthesized in the concept [[regulatory-element-detection]]. This is a
motif-family Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]] tracks
the unit.

## What this file records

- **Online sources (per-element primary literature, authority rank in parentheses):**
  - **Bucher (1990)**, *J Mol Biol* (rank 1) — weight-matrix descriptions of four Pol II
    promoter elements from 502 promoters; **TATA box** `TATAAA` (wider IUPAC
    `TATA(A/T)A(A/T)(A/G)`) and **CCAAT** pentanucleotide (core `CCAAT`, ~30% of promoters).
  - **Harley & Reynolds (1987)**, *Nucleic Acids Res* (rank 1) — prokaryotic σ70 promoter
    **−10** `TATAAT` and **−35** `TTGACA` hexamers, highly conserved, 17±1 bp spacing.
  - **Lundin, Nehlin & Ronne (1994)**, *Mol Cell Biol* (rank 1) — **GC box** `GGGCGG` (Sp1).
  - **Kozak (1987)**, *Nucleic Acids Res* (rank 1) — vertebrate **Kozak** context
    `gccRccAUGG`; most-preferred-base DNA string `GCCGCCACCATGG`; −3 purine / +4 G rules.
  - **Shine & Dalgarno (1974)**, *PNAS* (rank 1 via Wikipedia rank 4) — **SD/RBS** `AGGAGG`,
    ~8 bp upstream of AUG, complementary to 3' 16S rRNA `CCUCCU`.
  - **Proudfoot & Brownlee (1976)**, *Nature* (rank 1) — **poly(A) signal** `AATAAA`.
  - **Massari & Murre (2000)**, *Mol Cell Biol* (rank 1) — **E-box** `CANNTG` (canonical
    palindrome `CACGTG`); used as the IUPAC family pattern so any central dinucleotide matches.
  - **Lee, Mitchell & Tjian (1987)**, *Cell* (rank 1) — **AP-1 (TRE)** `TGACTCA` (palindrome
    `TGA(C/G)TCA`); the source of the corrected-defect regression.
  - **Sen & Baltimore (1986)**, *Cell* (rank 1) — **NF-κB** consensus `GGGRNWYYCC`; strong
    reference site scanned is `GGGACTTTCC`.
  - **Montminy et al. (1986)**, *PNAS* (rank 1) — **CREB (CRE)** palindrome `TGACGTCA`.

- **DEFECT FOUND (recorded in the file):** prior code used AP-1 `TGAGTCA` (G at position 4);
  corrected to `TGACTCA` per Lee, Mitchell & Tjian (1987). Also newly added the −10/−35
  prokaryotic promoter hexamers.

- **Documented corner cases / failure modes:** E-box `CANNTG` IUPAC-`N` degeneracy (any base
  in the two central positions); multiple occurrences of one element and of several elements
  in one sequence, each at its own 0-based start; substring-containment distinction between
  `TATAAT` (−10) and `TATAAA` (TATA); empty sequence → empty result (window `0 ≤ i ≤ n − m`);
  null sequence → `ArgumentNullException`.

- **Datasets (documented oracles):**
  - *Exact-match probes* — one per element, consensus copied verbatim from the primary
    sources, each expected at a specific 0-based position (e.g. TATA `GGGTATAAAGGG`→3,
    −10 `CCTATAATCC`→2, Kozak `TTGCCGCCACCATGGAA`→2, E-box `GGCACGTGGG`→2 via `CANNTG`,
    NF-κB `AAGGGACTTTCCAA`→2, CREB `CCTGACGTCAGG`→2).
  - *AP-1 negative control (defect regression)* — `AATGAGTCAGG` (contains old wrong
    `TGAGTCA`) → **0** AP-1 hits (correct consensus is `TGACTCA`).

- **Assumptions (both source-backed representative-site choices):** NF-κB scans the specific
  strong site `GGGACTTTCC` rather than expanding `GGGRNWYYCC`; Kozak scans the single
  most-preferred-base string `GCCGCCACCATGG` rather than expanding the −3/+4 degeneracy.

- **Test-coverage recommendations:** MUST — each of the 12 elements detected at the exact
  0-based position with correct Name/Pattern/Sequence; AP-1 `TGACTCA` matches while old
  `TGAGTCA` produces no hit; E-box degenerate `CACGTG` matches; null → `ArgumentNullException`,
  empty → no elements. SHOULD — multiple occurrences / multiple distinct elements. COULD —
  `KnownMotifs` constants equal their cited consensus strings.

## Deviations and contradictions

**No source contradictions** — the twelve consensus strings are each drawn from their
primary literature and are mutually consistent. The only deviations from a "fully degenerate"
reading are the two deliberate representative-site assumptions (NF-κB strong site, Kozak
exact context), both explicitly source-anchored in the file, neither an open correctness gap.
