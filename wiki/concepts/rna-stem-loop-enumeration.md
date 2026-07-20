---
type: concept
title: "RNA stem-loop / hairpin enumeration (sequence scan, tetraloops)"
tags: [rna, algorithm]
mcp_tools:
  - find_stem_loops
  - hairpin_loop_energy
sources:
  - docs/algorithms/RnaStructure/RNA_Stemloop.md
  - docs/Evidence/RNA-STEMLOOP-001-Evidence.md
source_commit: 8945f802f509d383569c666345edbecea758083e
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-stemloop-001-evidence
      evidence: "Test Unit: RNA-STEMLOOP-001 ... Title: Stem-Loop Detection ... Area: RnaStructure"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-base-pairing
      source: rna-stemloop-001-evidence
      evidence: "Stem-loops form via Watson-Crick and wobble base pairing; Valid base pairs A-U, U-A, G-C, C-G, G-U, U-G; the enumerator checks CanPair (WC or, if allowWobble, wobble) when extending the stem"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-dot-bracket-notation
      source: rna-stemloop-001-evidence
      evidence: "Test category 7 (Must): Dot-bracket notation generation; GGGAAAACCC -> (((....))) (3 opening brackets, 4 dots, 3 closing brackets)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:rna-pseudoknot-detection
      source: rna-stemloop-001-evidence
      evidence: "FindPseudoknots detects crossing base pairs: for pairs (i,j) and (k,l), pseudoknot if i < k < j < l — the same crossing condition RNA-PSEUDOKNOT-001 scans"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:rna-minimum-free-energy-folding
      source: rna-stemloop-001-evidence
      evidence: "The enumerator scans the sequence for ALL stem-loops meeting size constraints (O(n^2*L)), whereas the MFE folder searches for the single optimal fold via O(n^3) DP; Known Limitation 4: 'Single structure: does not enumerate suboptimal structures'"
      confidence: medium
      status: current
---

# RNA stem-loop / hairpin enumeration (sequence scan, tetraloops)

The **enumeration layer** of the RNA secondary-structure family (test unit **RNA-STEMLOOP-001**, area
`RnaStructure`): given a sequence, **scan it for every stem-loop (hairpin) it can form** under
size constraints — a stem of antiparallel [[rna-base-pairing|Watson-Crick + G-U wobble]] pairs capped
by a single-stranded apical loop. The record is [[rna-stemloop-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern. Prior RNA
ingests flagged a "general RNA secondary-structure (stem-loop) concept" as a future neighbour of
[[pre-mirna-hairpin-detection]]; this is that page.

## 1. What distinguishes this unit

A stem-loop is a **stem** (double-stranded helix, antiparallel WC/wobble pairing) plus a
**loop** (single-stranded apex). This unit **enumerates** them from a raw sequence — which sets it
apart from three sibling units that touch hairpins but do a different job:

| Unit | Job | Input | Cost |
|------|-----|-------|------|
| **RNA-STEMLOOP-001** (this) | *enumerate* all stem-loops under size constraints | sequence | O(n² · L) |
| [[rna-hairpin-001-evidence\|RNA-HAIRPIN-001]] | *score* one given hairpin's ΔG° | one hairpin | — |
| [[rna-minimum-free-energy-folding\|RNA-MFE-001]] | *find the single optimal fold* (DP) | sequence | O(n³) |
| [[pre-mirna-hairpin-detection\|MIRNA-PRECURSOR-001]] | *detect miRNA precursor* stem-loops | sequence | — |

It is the **general, sequence-agnostic** hairpin finder — not miRNA-specific, not a single-answer
optimizer, and not an energy calculator (though it consumes Turner 2004 terms for the tetraloop bonus).

## 2. The enumeration algorithm

`FindStemLoops(sequence, minStem, minLoop, maxLoop, allowWobble)` (canonical; `FindHairpins` is a
parameter-struct variant):

1. Scan candidate **loop positions** across the sequence.
2. For each loop size in `[minLoopSize, maxLoopSize]`, **extend the stem outward** on both sides while
   the mirrored bases still pair — Watson-Crick always, G-U wobble only if `allowWobble`.
3. Collect stems meeting `minStemLength`.

**Complexity: O(n² · L)** (n = length, L = max loop size). The stem grows against the
[[rna-base-pairing]] `CanPair` predicate; the accepted hairpin is emitted as a balanced
[[rna-dot-bracket-notation|dot-bracket]] string — `GGGAAAACCC` → `(((....)))` (3 stem pairs, 4-nt
loop). **Defaults:** wobble allowed, `minLoopSize = 3` (steric), `minStemLength = 3` (stability).

## 3. Loop-size constraints (steric)

- **Minimum 3 nt** — loops < 3 nt are **sterically impossible** and do not form (Wikipedia).
- **Optimal 4-8 nt** (Wikipedia). Loops > ~30 nt carry a higher Turner 2004 energy penalty.

This is the biophysical basis for the family-wide floor that also appears in the
[[rna-minimum-free-energy-folding|MFE folder]] (a pair `(i,j)` with `j − i − 1 < 3` cannot close a
hairpin) and the [[rna-free-energy-turner-model|Turner energy model]] (no ΔG° defined for loops < 3).

## 4. Tetraloops — the distinctive biological surface

Four-nucleotide loops are the most common and best-characterised hairpin loops, and three families
are exceptionally stable:

- **GNRA** — G, any N, purine R, A: GAAA, GCAA, GGAA, GUAA.
- **UNCG** — UACG, UCCG, UGCG, **UUCG** (the single most stable tetraloop, Antao et al. 1991).
- **CUUG** — CUUG, CCUG.

**UUCG + GNRA make up ~70 % of all tetraloops in 16S rRNA** (Woese et al. 1990). A **~3.0 kcal/mol**
stability bonus applies to GNRA/UNCG loops; the closing GA first mismatch contributes −0.9 kcal/mol
via the standard [[rna-free-energy-turner-model|Turner 2004]] model. Worked oracle: `GGGCGAAAGCCC` →
4-bp stem + `GAAA` GNRA tetraloop (no NNDB special-loop entry exists for GNRA with a C-G closing pair,
so its stability comes from the GA mismatch bonus and tertiary interactions rather than a table
override — the NNDB special tetraloops are all UNCG-family with C-G closing).

## 5. Pseudoknot detection (shared primitive)

`FindPseudoknots(sequence)` flags **crossing** base pairs — pairs (i,j) and (k,l) cross when
`i < k < j < l` — exactly the [[rna-pseudoknot-detection]] (RNA-PSEUDOKNOT-001) crossing condition
(Antczak 2018; Rivas & Eddy 1999). Oracle: pairs (0,6)+(3,9) give `0<3<6<9` → crossing detected.
**Prediction of pseudoknots from sequence is out of scope** — NP-complete for the general case; the
family's energy-based crossing-helix predictor is [[rna-pseudoknot-prediction]].

## 6. Invariants and edge cases

- **INV-01 (non-overlap):** enumerated stem-loops must not overlap (simple structure).
- **INV-02 (one pair per base):** each base participates in at most one base pair.
- **INV-03 (contiguous loop):** the loop is a single unbroken single-stranded run.
- **INV-04 (antiparallel stem):** stem pairs read 5'→3' against 3'→5'.
- **No structure ⇒ empty:** a homopolymer (`AAAA…`, no complement) or a too-short sequence
  (`GC`, `GCAUC`) yields no stem-loops. Empty `""` → empty result; null → empty or exception.
- **Case-insensitive:** `gggaaaaccc` ≡ `GGGAAAACCC`. **Wobble-only stems** (all G-U) accepted iff
  `allowWobble = true`.

## 7. Scope, limitations, relationships

A [[scientific-rigor|research-grade]] enumerator. Documented limitations (from the artifact, all
accepted scope boundaries, not bugs): **no pseudoknot prediction from sequence** (detection only);
**simplified energy model** (less accurate than ViennaRNA); **no internal loops or bulges** (apical
hairpin loops only); **single structure** (does not enumerate suboptimal folds). It builds on
[[rna-base-pairing]], emits [[rna-dot-bracket-notation]], shares its crossing test with
[[rna-pseudoknot-detection]], and is the general-purpose complement to the single-optimal-fold
[[rna-minimum-free-energy-folding|MFE folder]], the [[rna-hairpin-001-evidence|hairpin energy unit]],
and the miRNA-specific [[pre-mirna-hairpin-detection|precursor detector]]. **No source
contradictions** — Wikipedia, Woese 1990, Antao 1991, Heus & Pardi 1991, Rivas & Eddy 1999, and NNDB
Turner 2004 are mutually consistent.

## 8. Implementation surface (RNA-STEMLOOP-001 primary spec)

The canonical algorithm spec is `docs/algorithms/RnaStructure/RNA_Stemloop.md`; the code lives in
`RnaSecondaryStructure` in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`
(Implementation Status **Simplified** — a local motif detector, *not* a global structure optimizer). The
spec pins the concrete surface behind the Evidence-level enumeration description:

- **`FindStemLoops(string rnaSequence, int minStemLength = 3, int minLoopSize = 3, int maxLoopSize = 10,
  bool allowWobble = true)`** — the greedy candidate generator. Enumerates loop-start positions, then loop
  sizes in `[minLoopSize, maxLoopSize]`, then **extends the stem outward** from the loop boundaries while
  `GetBasePairType` returns a valid pair (Watson-Crick always; G-U wobble only if `allowWobble`, otherwise a
  wobble pair *stops* extension — ASM-02). Accepted pairs are **reversed to 5′→3′ order** before the result
  is built. Candidates below `minStemLength` are discarded.
- **`DetectPseudoknots(IReadOnlyList<BasePair>)`** — the same post-hoc crossing-pair detector (`i<k<j<l`)
  described in §5; it consumes an *already-assembled* base-pair set, and does **not** predict pseudoknots
  from sequence.

**`StemLoop` output contract** — `Start` / `End` (inclusive 5′/3′ bounds of the detected span), `Stem`
(coordinates, base pairs, stem free energy), `Loop` (coordinates, `Size`, sequence), `TotalFreeEnergy`, and
`DotBracketNotation` — a **local** dot-bracket covering *only the detected span* (INV-04 below), not the
whole input sequence.

**Scoring** — `TotalFreeEnergy = ΔG_stem + ΔG_hairpin`, where `ΔG_stem` comes from `CalculateStemEnergy`
(stacked base pairs) and `ΔG_hairpin` from `CalculateHairpinLoopEnergy` (including the special `G-U`
closure-detection path), both drawing on the [[rna-free-energy-turner-model|Turner 2004]] terms.

**Preconditions** — the sequence is **uppercased internally**; `minLoopSize < 3` is **clamped up to 3**; a
sequence that is empty or shorter than `minStemLength * 2 + minLoopSize` yields **no candidates**.

**Spec invariants** — **INV-01** `Loop.Size >= 3` (min-loop clamp); **INV-02** `Stem.Length >=
minStemLength`; **INV-03** the returned loop type is **always `Hairpin`** (this API constructs only apical
hairpins); **INV-04** `DotBracketNotation.Length == End - Start + 1` (the helper allocates exactly the
detected span).

**Overlap policy — reconciliation.** `FindStemLoops` performs **no overlap suppression**: it yields **every**
qualifying candidate, overlapping or not (spec §6.1). The non-overlap property in §6 (INV-01 there) is *not*
a property of this raw enumerator — it is enforced downstream by
[[rna-secondary-structure-prediction|`PredictStructure`]]'s greedy position-based selection, which sorts
candidates by `TotalFreeEnergy` and marks each chosen span's positions used. So `FindStemLoops` is the
**unfiltered candidate source**; the single non-overlapping *structure* is assembled one layer up. This
refines — does not contradict — the Evidence view, which described the assembled simple-structure output.

**Complexity** — `O(n² · L)` time (`n` = length, `L` = scanned loop-size range), `O(k)` space (`k` =
returned stem-loops), matching §2.
</content>
