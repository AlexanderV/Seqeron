---
type: concept
title: "Breakpoint detection from split (soft-clipped) reads"
tags: [structural-variant, algorithm]
sources:
  - docs/Evidence/SV-BREAKPOINT-001-Evidence.md
source_commit: e0ce1587e134c086efac45f05e5c8d110933190e
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: sv-breakpoint-001-evidence
      evidence: "Test Unit ID: SV-BREAKPOINT-001, Algorithm: Breakpoint Detection from Split (soft-clipped) Reads"
      confidence: high
      status: current
---

# Breakpoint detection from split (soft-clipped) reads

The **first ingested germline structural-variant (SV) unit** (SV-BREAKPOINT-001) and the **anchor
of the germline SV family** in this wiki. It localizes an SV **breakpoint** to a single base from
**split (soft-clipped) reads** and clusters the per-read junctions into supported calls. Validated
under test unit **SV-BREAKPOINT-001** ([[sv-breakpoint-001-evidence]]); [[test-unit-registry]]
tracks the unit and [[algorithm-validation-evidence]] describes the artifact pattern.

The evidence model is the standard **soft-clip-junction** paradigm — the **SAM/BAM CIGAR consume
rules** (samtools/hts-specs) plus the **ClipCrop** (Suzuki 2011), **SoftSearch** (Hart 2013) and
**Tattini/Magi 2015** split-read caller literature — see [[sv-breakpoint-001-evidence]] for the
source-by-source trace.

## The soft-clip signature of a breakpoint

A **split read** has one end **anchored** to the reference and the other end **soft-clipped** —
mapped imprecisely because it crosses an SV / indel breakpoint. In CIGAR, a **soft clip `S`
consumes the read (SEQ) but not the reference**, and can only occur at a read end:

- a **leading `S`** (e.g. `31S69M`) = a **left / L-breakpoint**;
- a **trailing `S`** (e.g. `69M31S`) = a **right / R-breakpoint**.

The **breakpoint is the marginal point between the clipped and matched portions** — the
aligned/clipped **junction coordinate**, giving **single-base resolution**. A read with **no `S`
operation carries no breakpoint signal**. The repository `SplitRead` record stores the junction
directly: `PrimaryPosition` = anchored reference start (SAM POS), `SupplementaryPosition` = the
junction/breakpoint coordinate, `ClipLength` = number of clipped bases.

## Clustering rule (positional, orientation- and chromosome-consistent)

Individual split-read junctions are grouped into one supported breakpoint under three constraints,
all directly from the sources:

- **Positional tolerance** — junctions cluster only when within a small window (**ClipCrop
  "clustered within 5-base differences"**, default tolerance 5 b), which absorbs mapping jitter.
- **Minimum support** — a breakpoint is called only when the cluster has **≥ min supporting clipped
  reads** (SoftSearch default 5, configurable down to 2; this unit's default baseline 2). The
  reported support = **number of clipped reads in the cluster**.
- **Same side + same chromosome** — reads combine only when their **clip orientation matches**
  (left vs right, SoftSearch) and they are on the **same contig** (SAM POS is chromosome-local).

A **minimum clip length** filter drops reads with **≤ 5 clipped bases** (SoftSearch) as too short
to be reliable. `RefineBreakpoint` narrows a region to the consensus junction supported by the
reads inside it.

## Invariants and edge cases

- Junctions within tolerance and meeting support → **exactly one** breakpoint at that junction.
- A clip stack **below min support** → **no breakpoint** (e.g. one isolated split read).
- Two groups **more than tolerance apart** (5000 vs 5100) are **not merged** into one breakpoint.
- Tolerance boundary: junctions exactly `tolerance` apart cluster; `tolerance + 1` apart do not.
- **Input validation** — empty input → empty output; null input → throws.

Worked oracles (from [[sv-breakpoint-001-evidence]]): **3 reads @ 5000 ± ≤ 5 b, same chr, support
≥ 2** → one breakpoint at ~5000; **1 isolated read** (min support 2) → none; **2 reads @ 5000 vs
5100** (gap > tolerance) → two below-support groups → none.

## Scope, assumptions, and relation to other units

A [[research-grade-limitations|research-grade]] method, **not for clinical use**. It operates on
**already-extracted `SplitRead` records** (junction + clip length), **not** raw BAM alignment
parsing. One flagged assumption: the sources fix the per-read junction and cluster within a
tolerance but do **not** prescribe the cluster's reported coordinate statistic — this unit reports
the **rounded mean** of member junctions (mirroring the sibling `ClusterSplitReads`), a
sub-tolerance choice that does not affect cluster membership or support count.

As the germline SV-family anchor it is orthogonal to the **oncology** read-evidence rearrangement
units: the split-read + discordant-pair **gene-fusion caller**
[[gene-fusion-detection-read-evidence]] (which counts breakpoint-supporting reads to call a
transcript fusion) and the copy-number-pattern **chromothripsis** classifier
[[chromothripsis-inference]]. It is likewise distinct from the gene-order, signed-permutation
[[genome-rearrangement-breakpoint-distance]] (rearrangement *distance* between whole gene orders,
not read-evidence breakpoint localization) and the chromosome-scale
[[synteny-and-rearrangement-detection]]. Its germline-SV siblings are the read-depth
[[read-depth-cnv-segmentation]] (SV-CNV-001) — CNV deletion/duplication calls from **read depth of
coverage** rather than split-read junctions — and the discordant-read-pair
[[discordant-pair-sv-detection]] (SV-DETECT-001), which reads the **paired-end mapping (PEM)
span/orientation signature** to classify the SV type; all three are orthogonal read-evidence
channels of the same family (integrated callers such as DELLY combine split reads with discordant
pairs, and this split-read channel is what recovers the single-base breakpoint and the inserted
sequence that PEM span alone cannot). Further siblings — SV genotyping/merging — would enrich this
anchor. No source contradictions.
