---
type: source
title: "Validation report: ASSEMBLY-OLC-001 (Overlap-Layout-Consensus assembly)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-OLC-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-OLC-001.md
source_commit: 288511090f725435564292b6bf3fc8fa05fb7b4c
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-OLC-001

The two-stage **validation write-up** for test unit **ASSEMBLY-OLC-001** (Overlap-Layout-Consensus —
the full OLC pipeline: overlap graph → layout → consensus), validated 2026-06-15 in a fresh context.
This is the *report* artifact that feeds one row of the [[validation-ledger]]; it records the
validator's **verdict** on both the algorithm description and the shipped code, within the wider
[[validation-and-testing]] campaign. The paradigm itself (three stages, overlap graph,
Hamiltonian-path layout, oracles and assumptions) is summarized in the concept
[[overlap-layout-consensus-assembly]], the anchor of the assembly OLC family. Distinct from
[[assembly-olc-001-evidence]] (the pre-implementation evidence artifact, sourced from
`docs/Evidence/`) — this is the independent re-validation verdict.

Canonical methods under test (`SequenceAssembler`):
`AssembleOLC(reads, parameters)`, `FindAllOverlaps(reads, minOverlap, minIdentity)` (+ cancellable
overload), and `FindOverlap(seq1, seq2, minOverlap, minIdentity)`.

## Verdict

**Stage A: ✅ PASS · Stage B: 🟡 PASS-WITH-NOTES · State: ✅ CLEAN.** No code defect; the
implementation faithfully realises the validated description. The single test-coverage gap was fixed
in session (see below). Full unfiltered suite **6494 passed, Failed: 0** (was 6493 + the new M5b);
`dotnet build` 0 warnings / 0 errors. The PASS-WITH-NOTES on Stage B reflects only the two documented,
sourced **intentional simplifications**, not a divergence.

## Stage A — description (algorithm faithfulness)

Theory checked against primary sources opened this session, independent of repo artifacts:

- **Compeau, Pevzner & Tesler (2011),** *Nat Biotechnol* 29:987–991 (DOI 10.1038/nbt.2023) — the
  overlap graph (one node per read; directed edge A→B when a suffix of A matches a prefix of B above
  a minimum threshold), layout = **Hamiltonian path**, and its **NP-completeness** (the motivation
  for heuristic layout).
- **Langmead (JHU) — "Overlap Layout Consensus assembly"** — three stages Overlap/Layout/Consensus
  (p.4); suffix-prefix overlap with minimum length `l` (p.5); "report only the longest suffix/prefix
  match" (p.10); transitive reduction + non-branching contig emission (p.21–25); consensus =
  per-column majority vote (p.28); all-pairs O(N²) vs suffix-tree O(N+a) (p.10/p.16).
- **Langmead (JHU) — "Assembly & Shortest Common Superstring"** — first law of assembly (p.16); the
  GTACGTACGAT 6-mer overlap graph (edges ≥ 4, weights 4 and 5, p.23–25); greedy-SCS and its
  sub-optimality (p.45/p.57); repeats foil assembly (p.58–62).

**Formula check.** Overlap = longest suffix of A equal to a prefix of B, length ≥ `l`, report only
the longest — matches Langmead OLC p.5/p.10 and SCS p.16/p.23. Coordinates are 0-indexed;
`Position2` is always 0 (prefix of seq2), internally consistent with the suffix-prefix definition.
Formulae are parameter-free (edge weight = overlap length; `MinOverlap`/`l` and `MinIdentity` are
caller-supplied) — no scoring tables or biological constants to source-check.

**Edge-case semantics.** Edgeless graph → singletons (Compeau overlap-graph definition); empty input
→ empty result (explicitly an ASSUMPTION, no source — an acceptable trivial identity).

**Independent cross-check (numbers, Python from the definition only — not repo code).** Every value
reproduced the Evidence table:

| Case | Re-derived value | Source basis |
|------|-----------------|--------------|
| `GTACGTACGAT` 6-mers @minOverlap 4 | exactly **12 directed edges**, weights {4,5} | Langmead SCS p.24–25 (exact set) |
| `FindOverlap("CTCTAGGCC","TAGGCCCTC")` | longest 6 (`TAGGCC`), pos1 3, pos2 0 | Langmead OLC p.5 |
| chain `AAAAACCCCC`·`CCCCCGGGGG`·`GGGGGTTTTT` @5-overlap | `AAAAACCCCCGGGGGTTTTT`, len 20 | superstring merge |
| identity 7/8 = 0.875 (suffix `ACGTACGT` vs prefix `ACGTACCT`) | length-8 overlap reported; len 9–12 windows all < 0.85 | identity gate |
| boundary `AAAACGTT`/`CGTTGGGG` | longest suffix-prefix 4 (`CGTT`); `ACGTACGT`/`CGTAAAAA` = 3 → no edge @minOverlap 4 | threshold boundary |

**Findings.** Stage A sound. The description is honest about the two intentional simplifications:
(1) **greedy best-successor chaining** instead of transitive-reduction + Hamiltonian search (justified
by NP-completeness), and (2) **consensus realized as chain concatenation** rather than column majority
vote (justified because for exact-overlap chains concatenation *is* the majority vote). Stage A =
**PASS**, no divergences.

## Stage B — implementation (code review)

Code path `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs`:

- `AssembleOLC` (L61–83) — overlap → `BuildContigsFromOverlaps` → `MinContigLength` filter (L77) → stats.
- `FindAllOverlaps` (L142–165) and cancellable overload (L177–212) — all ordered pairs i≠j, delegate to `FindOverlap`.
- `FindOverlap` (L221–242) — scans overlap lengths from `min(len1,len2)` down to `minOverlap`, returns the first (i.e. longest) window meeting `minIdentity`; pos1 = len1−L, pos2 = 0.
- `CalculateIdentity` (L247–260) — case-insensitive fraction of matching positions.
- `BuildContigsFromOverlaps` (L262–346) — greedy best-successor chaining; singleton fallback gated by `MinContigLength`.

**Formula realised correctly.** `FindOverlap`'s descending scan returns the longest qualifying window
(exact longest-suffix-prefix). `CalculateIdentity` is the case-insensitive matching fraction over the
window (matches doc §3.3). `FindAllOverlaps` excludes self-overlaps and reports one edge per ordered
pair. The greedy layout merges `successor[overlap:]` — exactly the superstring merge.

**Cross-verification (code vs external, test run).** Every case green: 12-edge GTACGTACGAT set (M1),
FindOverlap CTCTAGGCC/TAGGCCCTC len 6 pos1 3 pos2 0 (M3), chain merge (M4), 7/8 identity gate
accept@0.85 reject@0.95 (S1), minOverlap boundary accept@4 none@5 (S2), below-threshold no edges (S3),
and `MinContigLength` discard drop-5-base/keep-10-base (M5b). The cancellable `FindAllOverlaps`
(`CancellationToken.None`) yields the identical edge set (S5) and throws `OperationCanceledException`
on a pre-cancelled token (S6); both delegate to the same `FindOverlap`.

**Test-quality audit (HARD gate: PASS after adding M5b).** M1–M5, S1–S3, M3 assert exact
externally-derived values (12-edge set, exact contig strings, exact overlap length+positions, exact
accept/reject) — a deliberately-wrong implementation (shortest-overlap, off-by-one positions, wrong
merge) fails them. No green-washing: exact equality where exact values are known; the two property
tests (S4 bounds, C1 ≥-longest) are **inequalities by design** — they assert the documented
heuristic-limitation contract (INV-04 / ASM-02), the strongest honestly-assertable statement for a
greedy layout (Langmead SCS p.57: greedy is not optimal). All three public methods + the cancellable
overload + every Stage-A branch are exercised. **Gap found & fixed:** the `MinContigLength` discard
branch (doc §3.2) had no test; added **M5b** (`AssembleOLC_MinContigLength_DiscardsShortContigs`)
locking the contract (5-base read discarded, 10-base kept, `TotalLength` 10, `TotalReads` 2).

## Findings

- **No code defect. End-state ✅ CLEAN.** The single test-coverage gap (`MinContigLength` filter) was
  completely fixed in session with a sourced exact-value test; full suite green (6494/0).
- The Stage B **PASS-WITH-NOTES** is only the two documented, sourced intentional simplifications
  (greedy layout vs Hamiltonian-optimal; concatenation-consensus vs majority vote) — not a divergence.
- **No contradictions** among the sources — Compeau, Pevzner & Tesler (2011) and both Langmead notes
  give the identical overlap-graph / Hamiltonian-path / three-stage account; the re-derived 12-edge
  graph and all overlap/merge values match the source slides.
- **No follow-ups.**
