---
type: concept
title: "Repetitive element detection and classification"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/ANNOT-REPEAT-001-Evidence.md
  - docs/Validation/reports/ANNOT-REPEAT-001.md
  - docs/algorithms/Annotation/Repetitive_Element_Detection.md
  - docs/Evidence/GENOMIC-TANDEM-001-Evidence.md
  - docs/algorithms/Genomic_Analysis/Tandem_Repeat_Detection.md
  - docs/Evidence/REP-STR-001-Evidence.md
  - docs/Evidence/RNA-INVERT-001-Evidence.md
source_commit: f11e8bc7feeba4d997051309d93f661db1f53382
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-repeat-001-evidence
      evidence: "Test Unit ID: ANNOT-REPEAT-001 ... Algorithm: Repetitive Element Detection and Classification (tandem repeats, inverted repeats, repeat-class assignment)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: genomic-tandem-001-evidence
      evidence: "Test Unit ID: GENOMIC-TANDEM-001 ... Algorithm: Tandem Repeat Detection (GenomicAnalyzer.FindTandemRepeats); duplicate Registry entry resolved by consolidation with REP-TANDEM-001"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rep-str-001-evidence
      evidence: "Test Unit ID: REP-STR-001 ... Algorithm: Microsatellite / Short Tandem Repeat (STR) detection — perfect (default) and approximate/imperfect/interrupted (opt-in, Tandem Repeats Finder model)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-invert-001-evidence
      evidence: "Test Unit ID: RNA-INVERT-001 ... Algorithm: RNA Inverted Repeats (potential stem regions) — the W G W̄ᴿ IUPACpal model in the RNA secondary-structure family"
      confidence: high
      status: current
---

# Repetitive element detection and classification

Finding and typing the repeated substrings that make up a large fraction of genomes. Seqeron's
repeat analyzer covers three distinct sub-problems, each with its own source-backed definition,
validated together as [[annot-repeat-001-evidence|ANNOT-REPEAT-001]]. This page is the shared
anchor for the whole **repeats / tandem family** — sibling units (GENOMIC-REPEAT,
GENOMIC-TANDEM, microsatellite/STR, low-complexity, etc.) should link here rather than
re-deriving the same definitions. See [[test-unit-registry]] for how units are tracked and
[[algorithm-validation-evidence]] for the evidence-artifact pattern.

## The three sub-problems

### 1. Tandem repeats (head-to-tail)

A **tandem repeat** is a motif whose copies are directly adjacent with no intervening sequence —
a *head-to-tail* arrangement (e.g. `ATTCG ATTCG ATTCG`). The minimum is **two** adjacent copies;
a motif occurring once is not a tandem repeat. By repeat-unit length:

- **Microsatellite / STR** — 1–6 bp unit (Simple_repeat).
- **Minisatellite** — 10–60 bp unit.
- **Satellite** — larger units.

**Primitive-unit rule:** the shortest period is the canonical unit. `AAAAAA` is the
mononucleotide `A` repeated, not the dinucleotide `AA` or trinucleotide `AAA` — reporting a
non-primitive unit double-counts. This is the **annotation `RepeatAnalyzer` convention**
([[annot-repeat-001-evidence|ANNOT-REPEAT-001]]).

**Two entry points, different period handling.** Seqeron detects tandem repeats through two
methods over the same **exact-copy** model. `GenomicAnalyzer.FindTandemRepeats`
([[genomic-tandem-001-evidence|GENOMIC-TANDEM-001]], a consolidated duplicate of
REP-TANDEM-001) is a brute-force detector that does **not** canonicalize competing periods:
a run like `AAAA` is reported once per unit-length interpretation meeting the threshold
(period 1 ×4 *and* period 2 ×2). The annotation `RepeatAnalyzer` path instead applies the
primitive-unit rule above. Both of these default paths are **exact** — they report only
head-to-tail copies with no substitutions or indels, and over exact repeats both match the
formal definition (period = unit length, copy number ≥ 2). The *approximate* tandem copies of
Benson's Tandem Repeats Finder (1999) — interrupted / imperfect tracts — are covered by a
separate **opt-in** detector; see *Approximate STR detection* below. The exact-only default
was historically a Framework/Simplified [[research-grade-limitations|limitation]]; the opt-in
approximate path closes that gap for the microsatellite/STR case.

#### Approximate STR detection (Benson TRF model)

The opt-in `FindApproximateTandemRepeats` + `ComputeBernoulliStatistics`
([[rep-str-001-evidence|REP-STR-001]]) implement the *approximate*-repeat model of Benson's
**Tandem Repeats Finder** (1999): "two or more contiguous, **approximate** copies of a
pattern". Where the exact detectors *fragment* an interrupted tract into short perfect runs,
this path reports it as **one** repeat with quantified imperfection. Per reported repeat it
records the TRF statistics — period size, copy number, consensus size (built by **majority
rule** over period-aligned columns, so `ConsensusSize == Period` in the subset), **percent
matches** and **percent indels** (as a fraction of total alignment columns), and an alignment
**score** using the recommended TRF weights match `+2`, mismatch/indel `−7`. Only repeats
scoring at least `Minscore` (default **50**, ≈25 aligned chars) are reported.

Worked oracle: `CAGCAGCAGTAGCAGCAG` (copy 4 `C→T`) → one period-3 repeat, `CAG`×6, 17 match /
1 mismatch, score `17·2 − 7 = 27`, 94.4% matches; the perfect detector reports only `CAG`×3.
A single-base deletion tract scoring 51 clears the default gate; a below-50 tract does not.

`ComputeBernoulliStatistics` adds Benson's statistical layer: alignment of two **adjacent**
copies is modelled as *n* independent Bernoulli trials (heads = match), yielding **PM**
(matching probability = average percent identity between adjacent copies) and **PI** (indel
probability), with defaults **PM = 0.80 / PI = 0.10**. These adjacent-copy statistics are
deliberately *distinct* from the consensus-based percent-matches (e.g. adjacent PM `13/15` vs
consensus `17/18` for the tract above), exactly as Benson specifies ("between adjacent copies
… not … the consensus pattern"). The expected-matches moment `E[heads] = PM·d` is reproduced;
the probabilistic **k-tuple seeding** (`R(d,k,pM)` / `W(d,pI)` 95% percentiles from
non-redistributable simulation tables) is **not** — a deterministic exhaustive (start, period)
scan replaces it, a documented [[research-grade-limitations|research-grade residual]] that
changes which windows are examined, not the statistics of a reported repeat.

### 2. Inverted repeats (reverse-complement arms)

An **inverted repeat** (IR) is a left arm `W` followed downstream by its reverse complement.
Following IUPACpal (Hampson et al. 2021), an IR has the form **W W̄ᴿ** (perfect, ungapped), or
**W G W̄ᴿ** for a gapped IR with a spacer/loop `G`, `|G| ≥ 0`. An **imperfect** IR allows up to
`k` mismatches: Hamming distance `δ_H(W, W̄ᴿ) ≤ k`.

- **Zero gap ⇒ palindrome.** When `|G| = 0` the composite is an even-length reverse-complement
  palindrome (e.g. `GAATTC` → arm `GAA`, revcomp arm `TTC`).
- Detection parameters: minimum arm length, maximum arm length, maximum gap, maximum mismatches.

The **RNA** flavour of this exact `W G W̄ᴿ` model — antiparallel complementary arms that form a
potential **stem-loop** — lives in the RNA secondary-structure family as test unit
RNA-INVERT-001 ([[rna-invert-001-evidence]]): same IUPACpal definition, but the arm complement is
the RNA base-pairing rule {A-U, G-C} of [[rna-base-pairing]] and a non-zero loop is required (the
object is then a stem-loop, cf. [[pre-mirna-hairpin-detection]] / [[rna-dot-bracket-notation]]).
That unit reports the **perfect, ungapped `k = 0`** case only.

### 3. Repeat-class assignment

Assigning a query to a repeat class from the RepeatMasker/Repbase vocabulary:
**SINE, LINE, LTR, DNA** (DNA transposons), **Satellite, Simple_repeat, Low_complexity,
Small RNA, Unclassified/Unknown**. RepeatMasker itself classifies by *homology* — best
Smith-Waterman-Gotoh match above a score threshold against the Repbase library — and returns
Unknown when nothing matches.

## Deviation: classification is exact-substring, not scored alignment

Seqeron's `ClassifyRepeat(sequence, repeatDb)` **does not** run Smith-Waterman against a curated
Repbase library (out of scope for one unit). Instead it screens the query for library elements
**exactly contained** within it (element ⊆ query, one-directional) and assigns the class of the
**longest** such match; with no match it falls back to motif-size Simple_repeat classification,
else Unknown. Only the *matching relaxation* (exact substring vs. scored homology) is assumed —
the class vocabulary is source-backed. The one-directional containment prevents a trivially short
query from being forced into a class just because a longer consensus happens to contain its
letters. Documented as a Framework/Simplified [[research-grade-limitations|limitation]], not an
invented constant. The one-directional fix was **not** cosmetic: the two-stage validation
([[annot-repeat-001-report|ANNOT-REPEAT-001 report]], ledger row 58) found and fixed a real
defect here — the original bidirectional `query.Contains(element) || element.Contains(query)`
misclassified a 1-bp query (`ClassifyRepeat("A", …)` → `"SINE/Alu"`) because a longer consensus
contained those letters; RepeatMasker instead screens the query *for* library elements.

## Structural invariants (good test oracles)

- Every reported tandem repeat's `sequence` equals `input[start..end]` and is an integer number
  of unit copies.
- IR arms are exact reverse complements (within `k` mismatches for imperfect IRs).
- Reported spans are half-open `[start, end)` (the worked example `ATTCGATTCGATTCG` → unit
  `ATTCG`, 3 copies, span `[0, 15)`).
