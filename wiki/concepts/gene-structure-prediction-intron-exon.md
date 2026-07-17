---
type: concept
title: "Gene structure prediction (intron/exon): donor+acceptor pairing, exon typing/phase, spliced sequence"
tags: [splicing, algorithm]
mcp_tools:
  - predict_gene_structure
  - predict_introns
sources:
  - docs/algorithms/Splicing/Gene_Structure_Prediction.md
  - docs/Evidence/SPLICE-PREDICT-001-Evidence.md
source_commit: 2659a0374ef4d725033e34bd75a63dfb05662201
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: depends_on
      object: concept:splice-donor-site-prediction
      source: splice-predict-001-evidence
      evidence: "PredictIntrons pairs each 5' GT/GU donor with a 3' AG acceptor to define an intron; the donor end is detected by the SPLICE-DONOR-001 method. Burge et al. 1999: 'Each donor must pair with exactly one acceptor to define an intron.'"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:splice-acceptor-site-prediction
      source: splice-predict-001-evidence
      evidence: "PredictIntrons pairs each donor with a 3' AG acceptor (branch point + PPT + AG) detected by the SPLICE-ACCEPTOR-001 method to bound the intron."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: splice-predict-001-evidence
      evidence: "Test Unit ID: SPLICE-PREDICT-001 — Algorithm: Gene Structure Prediction (Intron/Exon)"
      confidence: high
      status: current
---

# Gene structure prediction (intron/exon)

Where [[splice-donor-site-prediction|donor]] and [[splice-acceptor-site-prediction|acceptor]]
detection each *score one boundary*, **gene structure prediction** is the **umbrella /
composite** that assembles those boundaries into a whole **eukaryotic split-gene model** —
an ordered list of **introns** and **exons**, a **spliced (mature) sequence**, and an
overall score. It is the third and integrative member of the **splicing family**
(validated under test unit **SPLICE-PREDICT-001** — record [[splice-predict-001-evidence]];
[[test-unit-registry]] tracks the unit, [[algorithm-validation-evidence]] the artifact
pattern). The distinct thing this unit adds over its two constituents is **pairing** and
**structure assembly**, not boundary scoring.

## From boundaries to a gene model

The pipeline is: find every candidate donor (a `GT/GU` 5' site) and acceptor (an `AG` 3'
site), then **pair** them into introns under the **GT-AG rule** (Breathnach & Chambon 1981:
>99% of spliceosomal introns begin `GT` and end `AG`; donor consensus `MAG|GURAGU`,
acceptor consensus `(Y)nNCAG|G`). Each **donor must pair with exactly one acceptor** to
define one intron (Burge et al. 1999). The exons are then the complement — the gaps
*between* the selected introns — and splicing them out yields the mature sequence
(Gilbert 1978: an exon is any part of a gene that survives into the mature RNA).

### Intron pairing constraints

- **`minIntronLength` / `maxIntronLength`** bound admissible donor→acceptor spans. The
  default `minIntronLength = 60` is biologically reasonable even though the shortest known
  metazoan intron is **30 bp** (human *MST1L*, Piovesan 2015) — a documented corner case
  where the default excludes real ultra-short introns.
- **Non-overlapping selection.** Candidate donor/acceptor pairs overlap; the implementation
  makes a **greedy-by-score** non-overlapping selection. The Evidence classifies this as a
  **design decision**, not a biological claim: external sources prescribe only the invariant
  (introns don't overlap), not the algorithm to enforce it.
- **Intron type.** Canonical `GT-AG` pairs are classified **U2** (major spliceosome);
  the doc also recognizes the `GC-AG` variant (~0.5% of U2 introns) and `AT-AC` **U12**
  (minor spliceosome, ~0.5%) — same non-canonical set as the two boundary units.

## Exon typing and phase

Once introns are fixed, each exon is assigned a **type** (Gilbert 1978 / Sakharkar 2002):

| Type | Position | Contents |
|------|----------|----------|
| **Initial** | first exon | 5' UTR + start codon |
| **Internal** | middle exons | coding only |
| **Terminal** | last exon | stop codon + 3' UTR |
| **Single** | whole gene, no introns | intronless gene (e.g. histones) |

Each exon also carries a **reading-frame phase** = `(sum of preceding exon lengths) mod 3`
(the first exon is phase 0 by definition; standard convention per Alberts 2002). This is
the same exon-phase arithmetic reused by the oncology fusion units
([[gene-fusion-detection-read-evidence]], [[fusion-breakpoint-frame-and-protein-prediction]]).

## Scoring and the resolved magic constants

The **overall gene-structure score** is the **mean of the per-intron scores** — a *defined*
quality metric, explicitly **not** a biological standard (no external source defines a gene
structure quality metric). Each intron's combined score is `(donor.Score + acceptor.Score) / 2`.
An earlier version added a **default branch-point score of 0.3** as a third term
(`(donor + acceptor + 0.3) / 3`); that arbitrary magic constant was **removed** — when no
branch point is found the score falls back to the plain two-term mean. Scores are normalized
to **[0, 1]**.

## Corner cases

- **Single-exon gene** (no `GT-AG`/`AT-AC` pair, or a high threshold) → **0 introns, 1
  exon of type Single** spanning the whole sequence. Prokaryotic genes and some eukaryotic
  genes (histones) are intronless (Matera & Wang 2014).
- **Empty / null input** → empty `GeneStructure`: 0 introns, 0 exons, empty spliced
  sequence, overall score 0.
- **DNA T-equivalence:** `T` is treated as `U` (implementation uses `ToUpperInvariant`,
  case-insensitive), so genomic DNA and mRNA inputs behave identically.

## Relation to other units

This is the **composite over** the two boundary detectors — it *consumes*
[[splice-donor-site-prediction]] (the `GT/GU` 5' site) and
[[splice-acceptor-site-prediction]] (the `AG` 3' site + branch point + polypyrimidine tract)
and adds intron/exon **pairing**, **exon typing/phase**, and the **spliced sequence** that
neither boundary unit produces on its own. It is distinct from — and downstream of — both:
the boundary units answer "is *this* site a real splice site?"; this unit answers "what is
the whole gene's exon-intron architecture?". The GT-AG worked model here is the same
canonical splicing model the two boundary pages document, so there are **no source
contradictions** across the three splicing units.

## Method contract (algorithm spec)

The canonical primary spec for this unit is
`docs/algorithms/Splicing/Gene_Structure_Prediction.md` (Test Unit **SPLICE-PREDICT-001**,
status *Simplified*). Two **public** entry points on
`SpliceSitePredictor` (`Seqeron.Genomics.Annotation`, `SpliceSitePredictor.cs`), plus three
private helpers:

- **`PredictGeneStructure(string sequence, int minExonLength = 30, int minIntronLength = 60, double minScore = 0.5)` → `GeneStructure`** — the orchestrator.
- **`PredictIntrons(string sequence, int minIntronLength = 60, int maxIntronLength = 100000, double minScore = 0.5)` → intron candidates** — the reusable pairing step. `maxIntronLength` is *public* here but **fixed to `100000`** when `PredictGeneStructure` calls it (callers can only widen the span via `PredictIntrons` directly).
- Private helpers: `SelectNonOverlappingIntrons(...)`, `DeriveExons(...)`, `GenerateSplicedSequence(...)`.

Both entry points first normalize the input to **uppercase RNA** (`T`→`U`); `null`/empty
returns an empty `GeneStructure` (0 introns, 0 exons, empty spliced sequence, score `0`).

### Pairing → scoring → selection pipeline

1. **Site detection.** `PredictIntrons` calls `FindDonorSites` and `FindAcceptorSites` with an
   **individual-site threshold of `minScore * 0.8`** and **`includeNonCanonical = true`** for
   both — i.e. the composite intentionally admits weaker and non-canonical (`GC-AG`, `AT-AC`
   U12) boundaries than a standalone boundary scan, and lets the *combined* score gate them.
2. **Pairing.** Every donor is paired with every downstream acceptor whose span satisfies the
   length window; **intron length = `acceptor.Position - donor.Position + 1`** (inclusive
   boundaries). `O(D · A)` pairwise.
3. **Branch-point contribution.** For each pair a branch point is searched in
   **`[acceptor.Position - 50, acceptor.Position - 18]`** with a **minimum branch-point score
   of `0.4`**. Combined intron score is **`(donor + acceptor + branch) / 3`** when a branch
   point is found, else **`(donor + acceptor) / 2`** (the old default-`0.3` third term is
   gone — see scoring section above).
4. **Candidate acceptance.** Keep introns with length in `[minIntronLength, maxIntronLength]`
   **and** combined score `≥ minScore`.
5. **Greedy non-overlap selection.** Sort candidates by **descending combined score** and take
   each intron only if it uses no already-occupied position (position-level overlap check).
6. **Exon derivation.** Exons are the gaps between selected introns; exons shorter than
   `minExonLength` are **omitted from the `Exons` list**. No introns selected ⇒ one `Single`
   exon over the whole sequence.
7. **Spliced sequence.** `GenerateSplicedSequence` removes the **selected intron intervals**
   directly from the normalized sequence.

### Output records and invariants

`GeneStructure { IReadOnlyList<Exon> Exons; IReadOnlyList<Intron> Introns; string
SplicedSequence; double OverallScore; }`; `Intron { int Start; int End; IntronType Type; }`
(inclusive boundaries from donor/acceptor positions); `Exon { ExonType Type; int? Phase; }`.

- **INV-01** — selected introns never overlap (`SelectNonOverlappingIntrons`).
- **INV-02** — `OverallScore` is the **mean of selected-intron scores**, or `0` when none are
  selected.
- **INV-03** — first exon phase `0`; later phases `(cumulative preceding exon length) mod 3`
  (`CalculatePhase`).
- **INV-04** — a sequence with no selected introns returns a single `Single`-typed exon
  covering the whole sequence.

### Spliced-sequence vs exon-record caveat

`GenerateSplicedSequence` is defined by **intron removal**, *not* by concatenating only the
retained `Exon` records — and the `minExonLength` filter is applied to the `Exons` list but
**not** to the spliced sequence. So a short inter-intron gap can be **dropped from `Exons`
yet still present in `SplicedSequence`**; the spec lists coordinate reconciliation between the
two as *not implemented* (caller-side interpretation required). This is a heuristic
splice-site-driven predictor: greedy (not globally optimal) selection, arithmetic-mean
(not probabilistic) scoring, and no HMM/GHMM gene model (no GenScan/Augustus-style DP).
