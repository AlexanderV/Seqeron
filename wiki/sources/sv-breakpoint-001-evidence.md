---
type: source
title: "Evidence: SV-BREAKPOINT-001 (breakpoint detection from split / soft-clipped reads)"
tags: [validation, structural-variant]
doc_path: docs/Evidence/SV-BREAKPOINT-001-Evidence.md
sources:
  - docs/Evidence/SV-BREAKPOINT-001-Evidence.md
source_commit: e0ce1587e134c086efac45f05e5c8d110933190e
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SV-BREAKPOINT-001

The validation-evidence artifact for test unit **SV-BREAKPOINT-001** — **Breakpoint Detection
from Split (soft-clipped) Reads**. The **first ingested unit of the germline structural-variant
(SV) family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized
in [[breakpoint-detection-split-reads]]. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — the split-read / soft-clip-junction breakpoint
  paradigm, corroborated across the SAM/BAM format spec and three SV-caller papers):**
  - **SAM/BAM Format Specification** (samtools/hts-specs `SAMv1.tex`, rank 2, official
    file-format spec) — the **CIGAR consume table**: soft-clipped **`S`** bases consume the
    **query (read) but NOT the reference**; aligned `M/=/X` and reference-only `D/N` advance the
    reference position. **POS** = "1-based leftmost mapping position of the first CIGAR operation
    that consumes a reference base." **SEQ length** = sum of `M/I/S/=/X` lengths, so **clipped
    bases ARE present in SEQ** (the clipped sequence is recoverable from the read). **Clip
    placement**: `S` may only sit at the read ends → a **leading `S` = left clip**, a **trailing
    `S` = right clip**.
  - **Tattini, D'Aurizio & Magi 2015** (Front Bioeng Biotechnol 3:92, rank 1 review) — a
    **split read** has "one end anchored to the reference and the other end mapping imprecisely
    owing to a structural variant / indel breakpoint"; **split-read (SR) methods give single-base
    resolution**; SR callers "search for **clusters** of split reads" to call a breakpoint.
  - **ClipCrop** (Suzuki et al. 2011, BMC Bioinformatics 12(S14):S7, rank 1) — CIGAR `31S69M` =
    31 left-clipped + 69 matched bases; the **breakpoint = the marginal point between clipped and
    matched sequence** (the aligned/clipped junction coordinate); **L- vs R-breakpoints**
    distinguished by which side is clipped; candidate breakpoints are **sorted and clustered
    within 5-base differences**.
  - **SoftSearch** (Hart et al. 2013, PLoS ONE 8(12):e83356, rank 1) — a putative breakpoint is
    defined when **≥ x soft-clipped reads begin at position y** (default **x = 5**, configurable
    down to **2**); clipped reads are **combined only when their left/right orientation matches**;
    **minimum clip length** — reads with **≤ 5 unmapped (clipped) bases are discarded**.

- **Documented corner cases / failure modes:**
  - **Below-support clip stacks** — a clip position with fewer than the minimum supporting reads
    is not reported (SoftSearch, default 5 / configurable 2).
  - **Short clips** — reads with ≤ 5 clipped bases are too short to be a reliable signal and are
    dropped (SoftSearch).
  - **Position jitter** — mapping imprecision spreads a true breakpoint across nearby positions;
    clustering within a small tolerance (ClipCrop 5 b) merges them into one call.
  - **Clip side / no-signal reads** — a leading `S` and a trailing `S` define opposite-sided
    breakpoints; a read with **no `S` operation carries no breakpoint signal** (SAM spec).
  - **Per-chromosome** — a breakpoint position is chromosome-local (SAM POS is contig-local); reads
    on different chromosomes are not clustered into the same breakpoint.

- **Datasets (deterministic, derived from the cited clip/junction rules — no raw BAM table is
  openly downloadable, so split-read records are constructed directly from the rules). The
  repository `SplitRead` record stores the junction directly: `PrimaryPosition` = the read's
  anchored reference start (SAM POS), `SupplementaryPosition` = the breakpoint coordinate at the
  aligned/clipped junction, `ClipLength` = number of clipped bases. Parameters: cluster tolerance
  **5 bases**, default minimum support **2** (SoftSearch configurable minimum; baseline default 5),
  minimum clip length **> 5** bases.**
  - **3 split reads** at junction 5000 ± ≤ 5 b, same chromosome, support ≥ 2 → **one breakpoint at
    ~5000**.
  - **1 isolated** split read, min support 2 → **no breakpoint**.
  - **2 split reads at 5000 vs 5100** (gap > tolerance) → two separate clip groups, each below
    support 2 → **no breakpoint**.

- **Coverage recommendations (8 items):** MUST — junctions within tolerance + meeting min support
  → exactly one breakpoint at the junction; a clip stack below min support → no breakpoint; two
  groups > tolerance apart are not merged; the reported breakpoint carries support = number of
  clipped reads in the cluster; reads on different chromosomes are not clustered. SHOULD — tolerance
  boundary (`tolerance` apart cluster, `tolerance+1` apart do not); `RefineBreakpoint` narrows a
  region to the consensus junction supported by in-region reads. COULD — empty input → empty
  output, null input → throws.

## Deviations and assumptions

- **ASSUMPTION — cluster breakpoint-position estimator.** The cited sources fix the per-read
  breakpoint at the clip/aligned junction (ClipCrop) and cluster reads sharing a clip position
  within a tolerance (ClipCrop 5 b; SoftSearch "beginning at position y"), but they **do not
  prescribe a single summary statistic** (mean vs mode vs min) for the cluster's reported
  coordinate. This unit reports the **rounded mean** of the member junction coordinates, mirroring
  the sibling `ClusterSplitReads` in the same class. Justified because with reads inside a
  ≤ tolerance window the mean lies in the same single-base neighbourhood the sources define as one
  breakpoint; the choice changes only the sub-tolerance reported coordinate, not which reads form
  the cluster or the support count.

No source contradictions — the SAM spec, ClipCrop, SoftSearch and the Tattini/Magi review agree on
the soft-clip-junction = breakpoint model, positional clustering within a small tolerance, and a
minimum-support / minimum-clip-length filter.
