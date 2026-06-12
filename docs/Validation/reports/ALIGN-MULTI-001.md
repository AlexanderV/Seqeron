# Validation Report: ALIGN-MULTI-001 — Multiple Sequence Alignment

- **Validated:** 2026-06-12   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`
  (internal variant `MultipleAlignClassic`; helpers `SelectCenterSequence`, `ReconcileAlignments`,
  `BuildConsensus`, `ComputeSumOfPairsScore`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia "Multiple sequence alignment"** (https://en.wikipedia.org/wiki/Multiple_sequence_alignment) confirms verbatim:
  - Length invariant: "all conform to length L ≥ max{n_i | i = 1, ..., m} **and no values in the sequences of S of
    the same column consists of only gaps**".
  - Gap-only insertion: "A multiple sequence alignment is taken of this set of sequences S by **inserting any
    amount of gaps needed** into each of the S_i sequences" — i.e. only gaps inserted, residues unchanged/in order.
  - Reversibility: "To return from each particular sequence S'_i to S_i, **remove all gaps**".
  - Sum-of-pairs: "the so-called *sum of pair* score".
  - Progressive method (Feng & Doolittle 1987): "developed by Da-Fei Feng and Doolittle in 1987"; "a first
    stage in which the relationships between the sequences are represented as a phylogenetic tree, called a
    *guide tree*" and "a second step in which the MSA is built by adding the sequences sequentially to the
    growing MSA according to the guide tree", "beginning with the most similar pair and progressing to the most
    distantly related."
- **Wikipedia "Clustal"** (https://en.wikipedia.org/wiki/Clustal) confirms the three canonical progressive steps:
  (1) "computes a pairwise distance matrix between all pairs of sequences"; (2) guide tree via neighbor-joining
  (midpoint rooting) or UPGMA; (3) "the guide tree is used as an approximate template to generate a global
  alignment"; k-tuple similarity scoring; "The gap opening penalty and gap extension penalty parameters can be
  adjusted by the user."

### Edge-case semantics check

The spec/evidence give sourced, defined behaviour for: empty input → `Empty`; single sequence → itself, SP=0
(no pairs); identical sequences → no gaps, consensus = input; different lengths → padded to common L ≥ max{n_i};
no all-gap column; reversibility by gap removal. All are backed by the sources above (or .NET convention for null).

### Independent cross-check (numbers)

- M10 `["ATGC","ATGC","CTGC"]`, SimpleDna (Match=+1, Mismatch=-1), same length ⇒ no gaps. SP over C(3,2) pairs:
  (0,1)=+4, (0,2)=+2, (1,2)=+2 ⇒ **SP=8** — matches spec and `ComputeSumOfPairsScore`.
- M05 three identical 8-mers: 3 pairs × 8 matches × 1 = **24** — matches.
- S02 BlastDna Match=+2 on `["ATGC","ATGC"]`: 1 pair × 4 × 2 = **8** — matches.

### Findings / divergences (the PASS-WITH-NOTES)

The implemented algorithm is **star alignment** (one center sequence chosen by 4-mer cosine similarity; every
other sequence aligned to that center; pairwise alignments reconciled into one coordinate space), **not** the
guide-tree progressive method of Feng-Doolittle/ClustalW. The TestSpec header labels it "Anchor-based star
alignment (progressive alignment with suffix tree anchors)" and lists `O(k²×m)`. This is honestly disclosed:
the Evidence doc explicitly flags center-selection (cosine, not k-tuple), suffix-tree anchoring, gap
reconciliation, gap-gap=0 SP convention, and consensus gap/tie handling as **implementation choices, not from
external sources**, and lists "No guide tree: simplified star alignment lacks the phylogenetic guide tree used
by ClustalW" as a known limitation. Crucially, the tests assert **structural invariants** (equal length, count,
reversibility, no all-gap column, hand-computed SP, majority consensus) and do **not** overclaim
Clustal-identical columns. This is the correct scoping for a heuristic, tool-dependent algorithm. The note is
purely that "progressive" in the title is a loose label — the method is star alignment, accurately documented as
such in the body and Evidence.

## Stage B — Implementation

- **Code path reviewed:** `SequenceAligner.cs:644-1027` (`MultipleAlign`, `SelectCenterSequence`,
  `ReconcileAlignments`, `BuildConsensus`, `ComputeSumOfPairsScore`).
- **Method realised correctly?** Yes for the *documented* method: empty/single handled at 651-660;
  center via cosine similarity (738-795); pairwise to center via `AnchorBasedAligner.AlignWithAnchors`;
  gap reconciliation takes max gaps-before-each-center-position across alignments and re-projects every
  sequence into the merged space (825-958) — gaps only inserted, center column characters preserved in order;
  final pad to common L (681-686); majority consensus with gap-vs-nucleotide tie → nucleotide (965-989);
  column-based SP over all C(k,2) pairs with gap-gap=0, gap-nuc=GapExtend (997-1027).
- **Invariants verified independently** (ad-hoc probe over adversarial multi-length sets k=4..5, run against the
  compiled library): equal row length = true; reversibility (`Replace("-","")` == original) = true for every
  sequence (no residue dropped/reordered); no all-gap column = true; consensus ∈ {A,C,G,T,-} and
  `len == aligned len` = true. SP values finite and consistent (e.g. 54, -12, 26).
- **Cross-verification recomputed vs code:** M10 SP=8, M05 SP=24, S02 BlastDna SP=8, M04 SP=4 — all reproduced
  by the running tests.
- **Variant/delegate consistency:** `MultipleAlignClassic` (internal fallback) uses the same `BuildConsensus`
  and `ComputeSumOfPairsScore`, so consensus/SP semantics are identical to the canonical path.
- **Test quality audit:** 18 canonical + 2 property MSA tests use exact, hand-computed assertions (not "no
  throw"/range checks) for SP and consensus, and assert each Stage-A invariant; C01 uses a fixed seed with a
  diversity guard. Genuine and deterministic.
- **Findings / defects:** none.

## Verdict & follow-ups

Stage A **PASS-WITH-NOTES** (invariants and SP/consensus definitions are externally correct; the algorithm is
star alignment, not guide-tree progressive — honestly documented, tests scope accuracy as heuristic/invariant
rather than Clustal-identical). Stage B **PASS** (code faithfully realises the documented star-alignment method;
all structural invariants — equal length, residue order preserved, gaps-only insertion, reversibility, no
all-gap column, sum-of-pairs, majority consensus — verified by tests and an independent probe).

- **Build:** `Seqeron.Genomics.Tests` builds, 0 warnings.
- **Tests:** `~MultipleAlign` filter = 31 passed / 0 failed; full suite = **4461 passed / 0 failed**.
- **Code changed:** none.
- **End-state:** ✅ CLEAN.
