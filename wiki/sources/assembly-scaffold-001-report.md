---
type: source
title: "Validation report: ASSEMBLY-SCAFFOLD-001 (scaffolding — ordered contigs joined with N-gaps)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-SCAFFOLD-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-SCAFFOLD-001.md
source_commit: e187814f8c9c929bb1cbf71405eb9f99777a9c92
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-SCAFFOLD-001

The two-stage **validation write-up** for test unit **ASSEMBLY-SCAFFOLD-001** (scaffolding — ordering
already-built contigs along a path and concatenating them interspersed with `N`-runs whose length is
the estimated inter-contig gap), validated 2026-06-15 in a fresh context. This is the *report*
artifact that feeds one row of the [[validation-ledger]]; it records the validator's independent
**verdict** on both the algorithm description and the shipped code. The algorithm itself is summarized
in the concept [[scaffolding]] (anchor of the assembly SCAFFOLD family), and the wider campaign is
[[validation-and-testing]]. Distinct from [[assembly-scaffold-001-evidence]] (the pre-implementation
evidence artifact, sourced from `docs/Evidence/`) — this is the independent re-validation verdict.

Canonical method under test:
`SequenceAssembler.Scaffold(IReadOnlyList<string> contigs, IReadOnlyList<(int,int,int)> links, char gapCharacter='N')`
— consumes **pre-ordered** links (each a `(from, to, gapEstimate)` triple); it does not rank links,
resolve overlaps, or orient contigs.

## Verdict

**Stage A: ✅ PASS · Stage B: ✅ PASS · State: ✅ CLEAN.** No defect found; the implementation
faithfully realises the validated description and no production code was changed. The only gap was
missing branch coverage in the tests, closed with two sourced/invariant-derived tests. Full unfiltered
suite is **6531 passed, Failed: 0** (was 6529 + 2 new), `dotnet build` 0 errors (4 pre-existing
warnings in unrelated files, none introduced).

## Stage A — description (algorithm faithfulness)

Theory checked against primary sources opened this session, independent of repo artifacts:

- **Jackman et al. (2017), ABySS 2.0 scaffold paper** — the construction rule verbatim: *"The
  sequences of the vertices in a path are concatenated, interspersed with gaps represented by a run of
  the character N, whose length corresponds to the estimate of the distance between those two
  contigs."* Plus the negative-estimate rule: *"It is possible that the distance estimate is negative,
  indicating that the two contigs should in fact overlap. If such an overlap is indeed found in the
  contig overlap graph, the two contigs are merged."* The paper states **no** minimum gap size.
- **NCBI AGP Specification v2.1** — *"Gap lengths must be positive. Negative gaps and gap lines with
  zero length are not valid."* and *"For negative gaps, or gaps of unknown size, use U as the
  component_type and 100 as the gap size, since 100 is the GenBank/EMBL/DDBJ standard for gaps of
  unknown size."*
- **Sahlin et al. (2012)** — corroborates that a negative gap (overlap) is a common, expected input
  class for de Bruijn assemblers.
- **Pop et al. (2004), Bambus** — scaffold = a path of *distinct* contigs; the greedy scaffolder joins
  "most links first".

**Formula / invariant check.** Concatenate path contigs interspersing an `N`-run of length = distance
estimate → matches Jackman verbatim (INV-01/03/05). A non-positive (zero/negative) estimate → 100 gap
characters → matches the AGP unknown-gap length (INV-02). Each contig in exactly one scaffold (path of
distinct contigs) → matches Pop et al. (INV-04).

**Edge-case semantics (sourced).** Zero gap → invalid per AGP → unknown default 100 N; negative gap →
overlap, unresolved here → unknown default 100 N (the **constant 100 is source-backed** — AGP; only
the *scoping decision* to fall back rather than resolve the overlap is a documented assumption,
Assumption Register #1). Out-of-range / self link ignored; null inputs throw (consistent with the
sibling `MergeContigs`).

**Independent cross-check (hand-computed vs sources, not vs code).**

| Inputs | Gap handling | Result | Length | Source basis |
|--------|--------------|--------|--------|--------------|
| contigs `["ACGT","TTGG","CCAA"]`, links `[(0,1,3),(1,2,2)]` | positive → `g`×N | `ACGTNNNTTGGNNCCAA` | 17 = 4+3+4+2+4 | Jackman rule (M1) |
| contigs `["AAAA","TTTT"]`, links `[(0,1,-5)]` | negative → 100 N | `AAAA`+100N+`TTTT` | 108 | AGP unknown-size (M3) |
| gap `0` | zero → invalid → 100 N | 100-N gap | — | AGP (M4) |

**Findings (Stage A).** The only divergence from cited theory is documented and scoped: Bambus joins
"most links first" whereas this unit follows **input link order** (consumes pre-ordered links; no
link-count ranking, no overlap resolution, no contig orientation). All three are explicitly listed as
intentional simplifications. No biological or mathematical error. **Stage A = PASS.**

## Stage B — implementation (code review)

Code path `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:679–744`, plus the
`UnknownGapLength = 100` constant at line 651.

**Formula realised correctly.** The gap-length rule `gapLength = gapSize > 0 ? gapSize : UnknownGapLength`
(line 733) is exactly INV-01/INV-02. The path-follow loop appends `gapCharacter` × gapLength then
`contigs[nextIndex]` (lines 733–736), realising the concatenate-with-N-run model verbatim. A `used`
HashSet enforces INV-04 (each contig once); unreached contigs become single-contig scaffolds in
ascending index order (outer `for`, line 709). Null guards (lines 684–685); empty → empty (687).

**Cross-verification.** Traced M1, M3, M4, S1 and the two added branch cases by hand against the code —
all match the externally sourced values above.

**Test-quality audit (HARD gate: PASS).** M1/M2 lock exact strings/run-lengths from the Jackman rule;
M3/M4 lock exactly 100 N from the AGP standard — a wrong `Math.Max(1, gapSize)` implementation (the
pitfall called out in the TestSpec) would fail these. All assertions are exact `Is.EqualTo` — no
widened tolerances, nothing skipped/ignored. The original 13-test suite covered M1–M7, S1–S5, C1, but
left **two real logic branches uncovered**; both were closed with sourced/invariant-derived tests:

- `Scaffold_SuccessorAlreadyPlaced_SkippedKeepingContigsDistinct` — `["AA","CC","GG"]`,
  `[(0,1,2),(1,2,3),(2,1,4)]` → `"AANNCCNNNGG"` (1 scaffold, length 11), from INV-03/04/05 (skip an
  already-placed successor and stop the path).
- `Scaffold_MultipleForwardLinks_FollowsFirstDeclared` — `["AA","CC","GG"]`, `[(0,1,2),(0,2,3)]` →
  `["AANNCC","GG"]`, locking the documented input-order tie-break (first declared link wins) plus the
  unreached-contig length-1 scaffold.

## Findings

- **No algorithm defect. No production code changed. End-state ✅ CLEAN.** Description, formula,
  invariants (INV-01…INV-05) and every edge case are independently confirmed against ABySS 2.0, the
  NCBI AGP spec, Sahlin et al. and Bambus; tests assert exact source-traced values.
- **No contradictions** among the sources — all give the identical "ordered contigs + sized `N`-gaps"
  scaffold model; the AGP 100-N unknown-size default and the ABySS negative-gap = overlap rule are
  complementary (file-format side vs assembler side of the same non-positive-gap case).
- **Follow-ups:** none. The documented simplifications (link-count ranking, overlap resolution, contig
  orientation) are out of scope and clearly disclosed; they are not correctness defects for the stated
  contract.
