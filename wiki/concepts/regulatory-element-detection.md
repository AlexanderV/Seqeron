---
type: concept
title: "Regulatory-element detection (canonical consensus catalog scan)"
tags: [motif, algorithm]
sources:
  - docs/Evidence/MOTIF-REGULATORY-001-Evidence.md
  - docs/algorithms/Motif_Discovery/Regulatory_Elements.md
source_commit: 914ab57b7357e92a2f33c107af9b6fdabe63ac45
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: motif-regulatory-001-evidence
      evidence: "Test Unit ID: MOTIF-REGULATORY-001 ... Algorithm: Regulatory Elements (scan a DNA sequence for known regulatory consensus motifs)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:known-motif-search
      source: motif-regulatory-001-evidence
      evidence: "Regulatory-element detection is the fixed-catalog specialization of the generic exact motif scan: each of the 12 named KnownMotifs is matched against the subject reporting every 0-based start position (MOTIF-REGULATORY-001 corner case: 'each occurrence is reported with its 0-based start position'), the same all-occurrences semantics known-motif-search applies to a caller-supplied set."
      confidence: high
      status: current
---

# Regulatory-element detection (canonical consensus catalog scan)

**Regulatory-element detection** scans a DNA sequence against a **fixed, curated catalog
of canonical regulatory sequence elements** — promoter boxes, translation-initiation
signals, and transcription-factor binding motifs — reporting each occurrence with its
name, matched pattern, and **0-based** start position. It is the *"you already know the
biology"* end of motif analysis: the query set is not caller-supplied literals but a
built-in constants table (`KnownMotifs`) whose entries are the **primary-literature
consensus strings** for each element. Validated under test unit **MOTIF-REGULATORY-001**;
the validation record is [[motif-regulatory-001-evidence]], [[test-unit-registry]] tracks
the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

## Why this is a distinct unit (not generic known-motif search)

The engine reuses the all-occurrences scanning of [[known-motif-search]]
(`GenomicAnalyzer.FindMotif`), but this unit owns something the generic finder does not: a
**named biological catalog**. The correctness of a regulatory-element scan is the
correctness of its *consensus constants* — each string must equal the value the source
literature assigns — so the validation is as much about the catalog contents (and their
citations) as about the scan. It also mixes two matching regimes: **exact** for most
elements and **IUPAC-degenerate** for the E-box (`CANNTG`, where `N` matches any base —
the [[iupac-degenerate-consensus|IUPAC ambiguity]] vocabulary applied to *matching*, not
generation). See **Scope** for the exact/degenerate boundary.

## The catalog (12 canonical elements, each source-anchored)

| Element | Consensus (5'→3') | Kind / role | Primary source |
|---|---|---|---|
| TATA box | `TATAAA` | Eukaryotic Pol II core promoter, −25..−35 of TSS | Bucher 1990 |
| −10 box (Pribnow) | `TATAAT` | Prokaryotic σ70 promoter −10 hexamer | Harley & Reynolds 1987 |
| −35 box | `TTGACA` | Prokaryotic σ70 promoter −35 hexamer (17±1 bp spacing) | Harley & Reynolds 1987 |
| CAAT box | `CCAAT` | Eukaryotic upstream promoter (NF-Y), ~30% of promoters | Bucher 1990 |
| GC box | `GGGCGG` | Sp1 zinc-finger binding site | Lundin et al. 1994 |
| Kozak | `GCCGCCACCATGG` | Vertebrate translation-initiation context (ATG at +1..+3) | Kozak 1987 |
| Shine-Dalgarno | `AGGAGG` | Prokaryotic ribosome-binding site (~8 bp upstream of AUG) | Shine & Dalgarno 1974 |
| Poly(A) signal | `AATAAA` | 3' cleavage/polyadenylation hexamer | Proudfoot & Brownlee 1976 |
| E-box | `CANNTG` | bHLH TF site (degenerate; canonical `CACGTG`) | Massari & Murre 2000 |
| AP-1 (TRE) | `TGACTCA` | AP-1/Jun-Fos binding palindrome | Lee, Mitchell & Tjian 1987 |
| NF-κB | `GGGACTTTCC` | κB site (strong reference of `GGGRNWYYCC`) | Sen & Baltimore 1986 |
| CREB (CRE) | `TGACGTCA` | cAMP-response-element palindrome | Montminy et al. 1986 |

The catalog spans **three functional layers** — transcription initiation (TATA / −10 /
−35 / CAAT / GC), translation initiation (Kozak / Shine-Dalgarno), 3' processing (poly(A)),
and TF binding (E-box / AP-1 / NF-κB / CREB) — deliberately covering both prokaryotic and
eukaryotic elements in one scan.

## The AP-1 corrected-consensus defect (regression anchor)

Prior repository code encoded AP-1 as `TGAGTCA` (G at position 4). Lee, Mitchell & Tjian
(1987) give the conserved recognition motif as **`TGACTCA`** (palindrome `TGA(C/G)TCA`);
the constant was corrected. The regression test is a **negative control**: the sequence
`AATGAGTCAGG` (containing the old wrong pattern) must yield **0** AP-1 hits. This is the
concrete failure mode the unit locks down — a catalog whose constants drift from the
literature silently mis-annotates.

## Matching contract and corner cases

| Aspect | Behaviour |
|---|---|
| Positions | **0-based** start index of each occurrence |
| Multiplicity | every occurrence of every element reported (multiple hits of one element and hits of several elements coexist) |
| Overlap / substring containment | `TATAAT` (−10) and `TATAAA` (TATA) differ only in the final base — both are distinct catalog entries and each exact hexamer present is reported independently |
| E-box degeneracy | `CANNTG` matches any central dinucleotide (`CACGTG`, `CAGCTG`, `CATATG`, …) via IUPAC `N` |
| Empty sequence | no occurrences → empty result (scan window `0 ≤ i ≤ n − m`) |
| Null sequence | `ArgumentNullException` |
| Per-hit payload | `Name`, `Pattern`, `Sequence` for each occurrence |

## Two source-backed representative-site assumptions

- **NF-κB scans the specific strong site `GGGACTTTCC`**, not the expanded degenerate
  consensus `GGGRNWYYCC` — the wild-type site from Sen & Baltimore (1986) / p50-binding
  studies. A source-backed representative-site choice, not a fabricated value.
- **Kozak is the single most-preferred-base string `GCCGCCACCATGG`**, not the degenerate
  `gccRccATGG`; the −3-purine and +4-G degeneracy is not expanded — only the canonical
  exact context is detected.

## Scope and siblings

This is the **canonical-catalog** motif scanner: a fixed, cited set of regulatory
consensus strings, mixing exact and one IUPAC-degenerate (`CANNTG`) match. It differs from
the generic exact [[known-motif-search]] (caller-supplied literal query set, no biology)
and from the *generation* direction of [[iupac-degenerate-consensus]] (which builds an
ambiguity code from an alignment rather than matching one against a subject). It is
unrelated to de novo [[overrepresented-kmer-discovery]] (finding **unknown** motifs) and to
alignment-derived [[consensus-from-alignment]]. The Kozak / Shine-Dalgarno translation
signals sit adjacent to [[open-reading-frame-detection]] (they mark the start context of a
coding ORF, though this scanner does not link them to ORFs) — the strand-aware,
spacing-scored SD-to-start ribosome-binding-site finder is
[[prokaryotic-gene-prediction-rbs]] (ANNOT-GENE-001), where this catalog only reports a
bare `AGGAGG` match anywhere. The distinct, scored −10/−35 promoter detector
[[promoter-detection]] (ANNOT-PROM-001, `GenomeAnnotator.FindPromoterMotifs`) is the
alternative treatment of the two prokaryotic promoter boxes: it reports the same `TATAAT` /
`TTGACA` boxes with **partial-variant score fractions** (prefix-5/suffix-5/prefix-4, score
∈ [0,1]), where this catalog scanner reports only the exact hexamer with no scoring. **No
source contradictions**; the only assumptions are the two representative-site choices above,
both source-anchored.
