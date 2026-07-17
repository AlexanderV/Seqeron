---
type: concept
title: "Variant calling (SNP/indel from reference↔query alignment + Ti/Tv)"
tags: [annotation, algorithm]
mcp_tools:
  - find_indels
  - find_snps
  - titv_ratio
sources:
  - docs/algorithms/Variants/Variant_Detection.md
  - docs/algorithms/Variants/Indel_Detection.md
  - docs/algorithms/Variants/SNP_Detection.md
  - docs/Evidence/VARIANT-CALL-001-Evidence.md
  - docs/Evidence/VARIANT-INDEL-001-Evidence.md
  - docs/Evidence/VARIANT-SNP-001-Evidence.md
source_commit: 34caa3137ac3d5e8e3e9f380abb130154edf3da0
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: variant-call-001-evidence
      evidence: "Test Unit ID: VARIANT-CALL-001 — Variant Detection (SNP / Insertion / Deletion calling from a reference↔query comparison, with transition/transversion classification)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: variant-indel-001-evidence
      evidence: "Test Unit ID: VARIANT-INDEL-001 — Indel Detection (FindInsertions/FindDeletions filters over the aligned-column caller); the indel facet of this same calling unit"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: variant-snp-001-evidence
      evidence: "Test Unit ID: VARIANT-SNP-001 — SNP Detection (FindSnps alignment-based + FindSnpsDirect Hamming-style substitution enumeration); the SNP facet of this same calling unit"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:variant-effect-annotation-vep
      source: variant-call-001-evidence
      evidence: "Calling produces the variant that annotation interprets: Danecek 2011 — 'a variant is a difference from reference'; VARIANT-ANNOT-001 takes an already-called variant and predicts its consequence, so calling is the upstream step feeding annotation"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:somatic-variant-calling-tumor-normal
      source: variant-call-001-evidence
      evidence: "Both detect variants but by different evidence: this unit calls germline SNP/indel from a reference↔query global alignment (CallVariantsFromAlignment); ONCO-SOMATIC-001 classifies Somatic/Germline from tumor-vs-matched-normal VAF thresholds"
      confidence: high
      status: current
---

# Variant calling (SNP/indel from reference↔query alignment + Ti/Tv)

**Variant calling** is the **detection** step of the variant-analysis family: it *produces* the variants
every downstream unit interprets. This unit calls **SNPs, insertions, and deletions** (Danecek 2011 — a
variant is a *difference from reference*) by aligning a query sequence against a reference and reading out
the differing columns, then classifies each substitution as a **transition** or **transversion** and
reports the **Ti/Tv** ratio. `CallVariantsFromAlignment` runs `SequenceAligner.GlobalAlign(reference,
query)` and walks the aligned columns, emitting one `Variant` per difference. Validated as test unit
**VARIANT-CALL-001** ([[variant-call-001-evidence]]), with the **SNP facet** (`FindSnps` /
`FindSnpsDirect`) validated separately as **VARIANT-SNP-001** ([[variant-snp-001-evidence]]) and the
**indel facet** (`FindInsertions` / `FindDeletions`) as **VARIANT-INDEL-001**
([[variant-indel-001-evidence]]); see [[test-unit-registry]] for how units are tracked
and [[algorithm-validation-evidence]] for the evidence-artifact pattern. Research-grade
([[research-grade-limitations]]), **not for clinical or diagnostic use**.

## Three variant classes from aligned columns

Walking the reference/query alignment column by column yields exactly three difference types
(Danecek 2011 — the VCF variant classes):

| Column | Class | In-memory representation |
|--------|-------|--------------------------|
| both bases present, differ | **SNP** | ref base, alt base, 0-based `Position` |
| gap in reference, base in query | **Insertion** | ref = `"-"` gap sentinel, alt base |
| base in reference, gap in query | **Deletion** | ref base, alt = `"-"` gap sentinel |

Identical sequences yield **zero variants**. Mismatched aligned lengths throw `ArgumentException`; empty
input yields an empty call set.

## The umbrella `CallVariants` caller (Variant_Detection.md, VARIANT-CALL-001 spec)

`CallVariants` is the **parent** of the SNP and indel facets: it is the single alignment-based caller that
emits **all three classes at once**, and both the SNP filters (`FindSnps`) and the indel filters
(`FindInsertions` / `FindDeletions` / `FindIndels`) are just `.Where(v => v.Type == …)` projections over
its output. All entry points are static methods on `VariantCaller` (`Seqeron.Genomics.Annotation`,
`VariantCaller.cs`):

- **`CallVariants(DnaSequence reference, DnaSequence query) → IEnumerable<Variant>`** — the umbrella caller.
  Null `reference`/`query` throws `ArgumentNullException`; it runs `SequenceAligner.GlobalAlign` (the shared
  **Needleman–Wunsch** engine in `Seqeron.Genomics.Alignment`, `SimpleDna` matrix, match +1 / mismatch −1 /
  linear gap −1, uppercase-normalizing) and delegates to the aligned-column scan.
- **`CallVariantsFromAlignment(string alignedReference, string alignedQuery) → IEnumerable<Variant>`** — the
  **pre-aligned entry point**: it takes two already-gapped, equal-length strings and does *only* the column
  walk, bypassing the O(n×m) DP. Null/empty input → empty; **unequal aligned lengths → `ArgumentException`**.
  Use it when the alignment is supplied by an external aligner.
- **`CalculateStatistics(DnaSequence reference, DnaSequence query) → VariantStatistics`** — the summary
  aggregator: calls `CallVariants`, then counts by class and computes the Ti/Tv and variant density in one
  pass. Null input → `ArgumentNullException`.

**`VariantType` enumeration.** The enum declares **five** members — `SNP`, `Insertion`, `Deletion`, `MNP`
(multi-nucleotide polymorphism), `Complex` — but the column walk emits **only the first three**; `MNP` and
`Complex` are reserved and never produced (multi-base events surface as *adjacent* per-base SNP/indel
columns, not a combined record). The parallel `MutationType` enum (`Transition` / `Transversion` / `Other`)
is the Ti/Tv classification axis, orthogonal to `VariantType`.

**`VariantStatistics` record** (returned by `CalculateStatistics`) bundles the call-set summary:
`TotalVariants`, `Snps`, `Insertions`, `Deletions`, `TiTvRatio` (from `CalculateTiTvRatio`), plus
**`VariantDensity`** = `variants.Count / reference.Length × 1000` (**variants per kilobase**, `0` when the
reference is empty), `ReferenceLength`, and `QueryLength`. This is the only place the per-kb density is
defined.

**`VcfPosition` opt-in accessor.** Every `Variant` stores a **0-based** `Position` (array convention); the
record also exposes `VcfPosition => Position + 1` so a caller consuming the record directly gets the VCF
**1-based POS** (VCFv4.3 §1.4.1) without re-deriving the offset — the same +1 shift `ToVcfLines` applies
when serializing. The sibling `AnnotatedVariant` record (`Variant` + `VariantEffect` + `MutationType`) is
the hand-off shape into [[variant-effect-annotation-vep]].

## SNP detection: FindSnps / FindSnpsDirect and the Hamming-mismatch invariant

The SNP facet is validated separately as **VARIANT-SNP-001** ([[variant-snp-001-evidence]]). Two entry
points return **substitution columns only** (no indels):

- **`FindSnps`** — alignment-based; a **filter over the caller** that runs the same reference↔query
  global alignment and keeps only the columns where both bases are present and differ (`Type = SNP`).
- **`FindSnpsDirect`** — positional / **Hamming-style**; compares the two sequences index-by-index and
  reports each mismatched position as one SNP. Over two **equal-length** sequences this enumerates exactly
  the **Hamming-mismatch positions** (Acharya 2017, PMC5410656 — "the number of positions that two
  codewords of the same length differ"); the SNP count equals the **Hamming distance**. Each SNP carries
  `Position = i` (0-based mismatched index), `ReferenceAllele = ref[i]`, `AlternateAllele = query[i]`,
  with `ref[i] ≠ query[i]` — a position where **REF == ALT is not a variant**.

The load-bearing scope claim is the **equal-length precondition**: the Hamming distance is defined only
for equal lengths, so `FindSnpsDirect` on unequal-length inputs compares only the **common prefix**
(`min(reference.Length, query.Length)`) — the trailing region of the longer input is **indel territory**
(VARIANT-INDEL-001), never a SNP. Identical sequences → zero SNPs (Hamming distance 0); base comparison
and Ti/Tv classification are **case-insensitive** (VCFv4.3 REF alphabet A,C,G,T,N is case-insensitive);
null inputs to `FindSnps` throw `ArgumentNullException`, empty inputs to `FindSnpsDirect` yield empty.
Oracles: `ATGC`→`ATTC` = {2} (G→T); `AAAA`→`TGTA` = {0,1,2} (A→T, A→G, A→T); VCFv4.3 §1.1 simple SNP
`G→A` @POS 14370.

**Method contract (SNP_Detection.md, VARIANT-SNP-001 algorithm spec).** The two entry points live on
`VariantCaller` (`Seqeron.Genomics.Annotation`, `VariantCaller.cs`) and differ in **input type** as well
as strategy:

- **`FindSnpsDirect(string reference, string query) → IEnumerable<Variant>`** — the **canonical**
  positional Hamming-mismatch enumerator, operating on raw **`string`** bases (not `DnaSequence`). Single
  forward scan over `n = min(reference.Length, query.Length)`; emits one `Variant` per mismatched index.
- **`FindSnps(DnaSequence reference, DnaSequence query) → IEnumerable<Variant>`** — the **delegate**:
  `CallVariants` filtered to `VariantType.SNP` (drops insertion/deletion columns), so it inherits the
  Needleman–Wunsch global alignment.

Each emitted `Variant` carries a **0-based** reference `Position == i` (VCF POS is 1-based — the 1-based
form is produced only by `ToVcfLines`), single-base `ReferenceAllele = reference[i]` and
`AlternateAllele = query[i]` (INV-04), `Type == SNP` (INV-02), and a **0-based `QueryPosition`** (`== i`
for `FindSnpsDirect`). Contract edges: `FindSnps` throws `ArgumentNullException` on a null reference/query
(validated in `CallVariants`); `FindSnpsDirect` returns **empty** when either input is null or empty and
otherwise scans only the common prefix (INV-06). Cost: `FindSnpsDirect` is **O(n) time / O(1)** lazy
(O(k) materialized for `k` SNPs), a single equality-test scan; `FindSnps` is **O(n × m)** dominated by the
alignment DP, with an O(n) filter. The repository **suffix tree does not apply** — SNP detection is a
positional equality test between corresponding strings, not an occurrence/substring search. No scoring
tables or tunable constants: the only decision rule is `reference[i] != query[i]`.

## Indel detection: FindInsertions / FindDeletions and the directional length invariant

The indel facet is validated separately as **VARIANT-INDEL-001** ([[variant-indel-001-evidence]]).
`FindInsertions` and `FindDeletions` are **filters over `CallVariants`** — they run the same
alignment-column caller and return only the columns of their class (`FindInsertions` → insertions only,
`FindDeletions` → deletions only; no SNPs). A **multi-base** indel is reported as **consecutive per-base
indel columns**. The load-bearing correctness claim is the **directional length invariant**:

- **insertion ⇒ ALT longer than REF** — serialized VCF `C → CA`; in-memory, `ReferenceAllele = "-"` gap
  sentinel + a one-base `AlternateAllele`.
- **deletion ⇒ REF longer than ALT** — serialized VCF `TC → T`; in-memory, a one-base `ReferenceAllele`
  + `AlternateAllele = "-"`.

The VCFv4.3 §1.1 microsatellite record `GTC → G,GTCT` shows both at once (one 2-base deletion `TC`, one
1-base insertion `T`, each anchored at `G`). Reference implementations of the normalized form —
minimal_representation (Minikel) worked cases CFTR p.F508del `(7,117199646,CTT,-) → (7,117199644,ATCT,A)`
and BRCA2 `(13,32914438,T,-) → (13,32914437,GT,G)` — independently confirm the same length direction and
the left-anchor padding of empty alleles.

**Method contract (Indel_Detection.md, VARIANT-INDEL-001 algorithm spec).** All three indel entry points
live on `VariantCaller` (`Seqeron.Genomics.Annotation`) and take `(DnaSequence reference, DnaSequence
query)`, returning `IEnumerable<Variant>`:

- **`FindInsertions`** / **`FindDeletions`** — `CallVariants` filtered to `VariantType.Insertion` /
  `VariantType.Deletion` respectively (INV-02: every returned variant carries the matching `Type`).
- **`FindIndels`** — the **union** delegate (insertions ∪ deletions), the surface behind the `find_indels`
  MCP tool.

Each `Variant` carries a **0-based** reference `Position` (INV-06: `Position ∈ [0, reference.Length]`,
because `refPos` advances only on reference-consuming columns) plus a **`QueryPosition`** (0-based query
coordinate of the event) — the query-side coordinate is retained alongside the reference-side one. Bases
are **case-normalized to uppercase** by `SequenceAligner.GlobalAlign` (match +1 / mismatch −1 / linear gap
−1 `SimpleDna` matrix); a null `reference` or `query` throws `ArgumentNullException` (propagated from
`CallVariants`), empty sequences yield no variants. Worked oracle: `ATGCAT` → `ATGTCAT` (a `T` inserted
after index 2) gives one insertion `Type=Insertion, ReferenceAllele="-", AlternateAllele="T", Position=3`.

Cost is dominated by the Needleman–Wunsch DP: **O(n × m) time and space** over reference length `n` and
query length `m`; the column walk + class filter are only O(n + m). The repository **suffix tree
(exact-occurrence enumeration) does not apply** — indel detection is scoring-based optimal alignment with
gaps, not exact substring matching, so the shared `SequenceAligner` is reused instead.

## Transition / transversion classification

Every **SNP** is classified case-insensitively as a **transition** (a base change *within* a ring class —
purine↔purine `A↔G` or pyrimidine↔pyrimidine `C↔T`) or a **transversion** (*across* ring classes —
`A↔C`, `A↔T`, `G↔C`, `G↔T`). There are **2 possible transitions but 4 possible transversions** per base,
yet transitions occur *more often* — the **transitional bias** (Collins & Jukes 1994: transition rate
1.71×10⁻⁹ > transversion rate 1.22×10⁻⁹; ≈2/3 of SNPs are transitions). The summary statistic is the

```
Ti/Tv = (#transitions) / (#transversions)
```

**Convention:** the mathematically-undefined zero-transversion case (`#Tv = 0`) is mapped to **0**, not
`+∞` or an exception. Insertions/deletions are not SNPs, so they classify as `Other` (Ti/Tv is defined
only for substitutions). A genome-wide Ti/Tv well above 0.5 (the naive 2:4 expectation) is the standard
sanity check for a real vs artefactual SNP call set.

## Where it sits in the variant-analysis family

Calling is the **head** of the germline variant pipeline. Its output feeds
[[variant-effect-annotation-vep]] (VARIANT-ANNOT-001), which takes each *already-called* variant and
predicts its functional consequence (missense / stop_gained / frameshift …) plus Sequence-Ontology
IMPACT — annotation cannot run without a caller upstream. This unit is the **germline, reference↔query**
caller; the oncology sibling [[somatic-variant-calling-tumor-normal]] (ONCO-SOMATIC-001) is
**`alternative_to`** it — same goal (find variants), different evidence: tumor-vs-matched-normal **VAF**
thresholds rather than a pairwise alignment. Serialized output (`ToVcfLines`) targets the **VCFv4.3**
format; VEP-style annotation and the oncology tier/pathogenicity layers consume such calls.

## Scope and assumptions

- **Internal gap-sentinel indel representation** — the in-memory `Variant` uses `"-"` and a **0-based**
  `Position`; the VCF **padding base** and **1-based POS** conventions (VCFv4.3 field 4) apply only to the
  *serialized* `ToVcfLines` output, not the in-memory model. E.g. the VCFv4.3 §1.1 examples `GTC→G`
  (2-base deletion) and `GTC→GTCT` (1-base insertion) both anchor at the preceding `G` — a serialization
  concern, not a detection concern.
- **Indels are not left-aligned / parsimony-normalized** — per Tan 2015 the canonical representation
  requires **left-alignment + parsimony** (Algorithm 1: right-trim shared trailing nucleotides, re-padding
  left when an allele empties, then left-trim shared leading nucleotides while all alleles are length ≥2 —
  the suffix-then-prefix trimming that minimal_representation and PharmCAT implement). This caller reports
  the indel at the column `GlobalAlign` produces, with **no** normalization pass. This affects reported
  **position** in low-complexity / repeated regions only (the alignment there is non-unique — e.g.
  PharmCAT's tandem-repeat `AATGA→A` @97740414 left-shifts to `GATGA→G` @97740410, same count and type),
  **not** variant counts or types.
- **Alignment-based, not pileup-based** — this is a single-query-vs-reference comparison, not a
  read-pileup genotype caller; there is no depth model, base-quality weighting, or diploid genotype
  assignment.

Reference sources — **VCFv4.3** (samtools/hts-specs, the REF/ALT/POS + padding grammar), **Danecek 2011**
(the VCF variant classes), **Tan 2015** (left-align + parsimony normalization), **Collins & Jukes 1994**
(transitional bias), and Wikipedia Transition/Transversion (the classification table), and **Acharya 2017 / PMC5410656**
(Hamming distance = count of differing positions) — full record in [[variant-call-001-evidence]]; the
SNP-detection facet (`FindSnps`/`FindSnpsDirect`, Hamming-mismatch enumeration, equal-length/common-prefix
precondition) in [[variant-snp-001-evidence]]; the indel-detection facet (`FindInsertions`/`FindDeletions`,
directional length invariant, normalization theory) in [[variant-indel-001-evidence]]. No source
contradictions.
