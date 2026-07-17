---
type: concept
title: "Discordant read-pair (PEM) structural-variant detection"
tags: [structural-variant, algorithm]
mcp_tools:
  - find_discordant_pairs
sources:
  - docs/algorithms/StructuralVar/SV_Detection.md
  - docs/Evidence/SV-DETECT-001-Evidence.md
  - docs/Validation/reports/SV-DETECT-001.md
source_commit: a408154d453bce5f8aa5336dd28c48ed6642fc5d
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: sv-detect-001-evidence
      evidence: "Test Unit ID: SV-DETECT-001, Algorithm: Structural Variant Detection from Paired-End Mapping (PEM) signatures"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:breakpoint-detection-split-reads
      source: sv-detect-001-evidence
      evidence: "SV-DETECT-001 is the discordant-read-pair (PEM) member of the germline SV family anchored by SV-BREAKPOINT-001; DELLY (ref 5) does 'integrated paired-end and split-read analysis' — the two are orthogonal read-evidence signatures of the same breakpoint."
      confidence: high
      status: current
---

# Discordant read-pair (PEM) structural-variant detection

The **discordant-read-pair member of the germline structural-variant (SV) family** (SV-DETECT-001).
Where the SV anchor [[breakpoint-detection-split-reads]] localizes a breakpoint from **within-read
soft-clip junctions** and [[read-depth-cnv-segmentation]] calls copy number from **aggregate depth**,
this unit reads the **paired-end mapping (PEM) signature** — how a mate pair's mapped **span** and
**orientation** deviate from the concordant expectation — and **classifies the SV type** it implies.
It is the third, orthogonal read-evidence channel of the same family (DELLY does "integrated
paired-end and split-read analysis"). Validated under test unit **SV-DETECT-001**
([[sv-detect-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. The two-stage validation
verdict is recorded in [[sv-detect-001-report]] — **Stage A / Stage B both PASS-WITH-NOTES,
State ✅ CLEAN** (the RF-everted pair was corrected from concordant to a discordant
tandem-Duplication signature in-session; full suite 6485/0).

The evidence is the standard **PEM / discordant-pair** paradigm — **Medvedev, Stanciu & Brudno 2009**
(the signature catalogue) plus the **BreakDancer** (Chen 2009) reference cutoffs and the
**DELLY / SVXplorer / LUMPY / Manta** orientation→type conventions; see [[sv-detect-001-evidence]]
for the source-by-source trace.

## Concordant vs discordant: the reference expectation

A short-insert Illumina library has one concordant orientation: **FR** — the upstream mate on the
`+` strand, the downstream mate on `−`, pointing toward each other, with a mapped span drawn from
the empirical **insert-size distribution** (mean μ, s.d. σ). This is the SAM/BAM **proper-pair flag
`0x02`** (BWA sets it only for FR). A pair is **discordant** when its span or its orientation
departs from that expectation. Note **RF (outward-facing / everted) is NOT concordant** — it is a
discordant signature, not a proper pair.

## The signature → SV-type classification

The **span cutoff** is expressed in standard deviations (BreakDancer `-c`, default **3**): a pair is
span-discordant iff `span < μ − c·σ` **OR** `span > μ + c·σ`. Combining the span and orientation
tests yields the six PEM classes (Medvedev 2009 signatures; BreakDancer DEL/INS/INV/ITX/CTX;
DELLY/SVXplorer orientation rules):

| Signature (same chromosome unless noted) | SV type |
|---|---|
| span within `[μ−c·σ, μ+c·σ]`, FR | **concordant** — no SV |
| span `> μ + c·σ`, FR | **Deletion** (mates span the deleted segment ⇒ larger span) |
| span `< μ − c·σ`, FR | **Insertion** (inserted bases push mates closer ⇒ smaller span) |
| same-orientation **FF / RR** | **Inversion** (one mate's orientation flipped) |
| everted **RF** | **tandem Duplication** |
| mates on **different chromosomes** | **Translocation** (inter-chromosomal linking signature) |

## Signature-then-cluster + minimum support

PEM callers "distinguish themselves by the signatures they can detect and the way they cluster …
these signatures." Discordant pairs supporting the **same event** are **clustered**, and a cluster
is reported as an SV only when it has **≥ min supporting read pairs** (BreakDancer `-r`, default
**2**). A single discordant pair below that threshold yields no call.

## Corner cases and failure modes

- **Insertion larger than the insert size is invisible** — the small-span insertion signature does
  not appear once the insertion exceeds the fragment length, and the PEM span **never recovers the
  inserted sequence itself** (Medvedev 2009). Localizing/assembling that sequence is the job of the
  split-read [[breakpoint-detection-split-reads]] channel, not PEM span.
- **Below-support clusters** — clusters with fewer than min-support pairs (default 2) are dropped
  (BreakDancer).
- **Cutoff boundary** — a span exactly at `μ ± c·σ` is concordant; one unit beyond is discordant.
- **Empty input → empty output.**

## Assumptions and scope

- **ASSUMPTION — inter-chromosomal precedence.** When mates map to **different chromosomes** the
  translocation (CTX) signature is reported **regardless of relative orientation**: chromosome
  difference is evaluated **first**. Justified because "inversion" (INV) is defined only for
  intra-chromosomal flipped pairs (Medvedev 2009; BreakDancer defines CTX for cross-chromosome
  mates), so a flipped orientation across chromosomes is undefined as an inversion and the event is
  by definition a translocation.

A [[research-grade-limitations|research-grade]] method, **not for clinical use.** It consumes
**already-mapped read-pair records** (span + orientation + chromosome per mate), not raw BAM
parsing.

## Method contract (algorithm spec)

The **primary spec** for this unit is `docs/algorithms/StructuralVar/SV_Detection.md`; the algorithm
lives in `StructuralVariantAnalyzer` (Seqeron.Genomics.Annotation, `StructuralVariantAnalyzer.cs`),
status **Simplified**. Three static entry points form the pipeline:

- **`FindDiscordantPairs(readPairs, expectedInsertSize=400, insertSizeStdDev=50, cutoffSd=3.0)`** →
  `IEnumerable<ReadPairSignature>` — streams the pairs flagged anomalous (interchromosomal, span
  outside `μ ± c·σ`, or non-FR orientation). O(n) time, O(1) streaming space.
- **`ClassifySV(signature)`** → `SVType` — maps one pair's PEM signature to a type. O(1).
- **`DetectSVs(readPairs, expectedInsertSize=400, insertSizeStdDev=50, cutoffSd=3.0, clusterDistance=500, minSupport=2)`**
  → `IEnumerable<StructuralVariant>` — the canonical entry point: find discordant → cluster → support
  gate → classify. O(n log n) time / O(n) space (dominated by the sort of discordant pairs; the
  internal `ClusterDiscordantPairs` is a linear adjacency sweep, not a windowing/CBS model).

**Input record.** `readPairs` is an `IEnumerable<(string ReadId, string Chr1, int Pos1, char Strand1,
string Chr2, int Pos2, char Strand2, int InsertSize)>` — positions **0-based**, strands `'+'`/`'−'`.
Chromosome names are compared by **raw string equality** (no normalization). Defaults: μ=400, σ=50,
c=3.0, clusterDistance=500, minSupport=2 (constraints: μ>0, σ≥0, c≥0, clusterDistance≥0, minSupport≥1).

**Output record.** `StructuralVariant` carries `Type, Start, End, Length, SupportingReads, Quality`
— one per qualifying cluster. `SVType` is `Deletion | Insertion | Inversion | Duplication |
Translocation | ComplexRearrangement`.

**Classification order (first match wins).** `ClassifySV` applies §2.2 in order: (1) `Chr1 ≠ Chr2` →
Translocation; (2) `Strand1 == Strand2` → Inversion; (3) RF (`Strand1=='−', Strand2=='+'`) →
Duplication; (4) FR and `span > μ+c·σ` → Deletion; (5) FR and `span < μ−c·σ` → Insertion; (6)
otherwise → **`ComplexRearrangement`** (the "everything else" bucket). `FindDiscordantPairs` also
flags any pair above a hard **`maxInsertSize` guard (default 10000)**; such a pair with no matching
basic signature falls through to `ComplexRearrangement`.

**Preconditions.** Null `readPairs` / `discordantPairs` → `ArgumentNullException`; empty input →
empty result. A span **exactly** at `μ ± c·σ` is concordant (bound inclusive; discordant iff strictly
outside). `ClassifySV` assumes its input pair is already anomalous — a concordant pair passed directly
is classified by the same rules but would not normally be produced by `FindDiscordantPairs`. Invariants
`INV-01`…`INV-06` (translocation precedence, same-strand⇒inversion, larger/smaller-span DEL/INS, span
cutoff, min-support gate) are covered by `StructuralVariantAnalyzer_DetectSVs_Tests.cs`. This unit runs
no substring/pattern search over a sequence, so the repository suffix tree is not applicable.

## Relation to the other germline-SV units

- Sibling read-evidence signature of the SV anchor [[breakpoint-detection-split-reads]]: PEM span/
  orientation gives **SV type + approximate locus from mate geometry**, while split reads give the
  **single-base breakpoint** at the clip junction — orthogonal, and integrated callers (DELLY) use
  both. PEM alone cannot recover an insertion's sequence; the split-read channel can.
- Complementary to the read-depth [[read-depth-cnv-segmentation]] (SV-CNV-001): read depth measures
  **copy number** (dosage) of del/dup, while a PEM del/dup signature is **breakpoint-spanning read
  geometry** — depth vs discordant-pair evidence for the same copy-number event.
- Distinct from the oncology read-evidence rearrangement units — the fusion caller
  [[gene-fusion-detection-read-evidence]] and the copy-number-pattern [[chromothripsis-inference]] —
  and from the gene-order [[genome-rearrangement-breakpoint-distance]]. No source contradictions.
