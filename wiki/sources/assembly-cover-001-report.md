---
type: source
title: "Validation report: ASSEMBLY-COVER-001 (coverage / depth calculation — per-base read depth)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-COVER-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-COVER-001.md
source_commit: 5a20c8c6a7047e13a84f9574bbe39793571da8f6
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-COVER-001

The two-stage **validation write-up** for test unit **ASSEMBLY-COVER-001** (Coverage /
Depth Calculation — per-base sequencing depth over a reference), validated 2026-06-15 in a
fresh context. This is the *report* artifact that feeds one row of the [[validation-ledger]];
it records the validator's **verdict** on both the algorithm description and the shipped code.
The algorithm itself is summarized in the concept [[coverage-depth-calculation]] (anchor of the
assembly COVER family), and the wider campaign is [[validation-and-testing]]. Distinct from
[[assembly-cover-001-evidence]] (the pre-implementation evidence artifact, sourced from
`docs/Evidence/`) — this is the independent re-validation verdict.

Canonical method under test:
`SequenceAssembler.CalculateCoverage(string reference, IReadOnlyList<string> reads, int minOverlap = 20)` → `int[]`
(per-base depth array), with internal helper `FindBestAlignment` tested indirectly.

## Verdict

**Stage A: 🟡 PASS-WITH-NOTES · Stage B: ✅ PASS · State: ✅ CLEAN.** No algorithm defect; the
depth-counting arithmetic is exactly source-faithful. The single PASS-WITH-NOTES flags a
documented **placement-model limitation** (below), not an arithmetic error. Full unfiltered
suite **6532 → 6533 passed, Failed: 0** (1 pre-existing skipped benchmark); `dotnet build` 0
errors (4 pre-existing warnings in unrelated files). One test-coverage gap (empty reference)
was found and closed in-session.

## Stage A — description (algorithm faithfulness)

Theory checked against sources opened this session, independent of repo artifacts:

- **Metagenomics Wiki — SAMtools breadth of coverage** — per-base depth = "number of reads
  mapping to a specific reference position"; breadth = covered bases / reference length (worked
  32876 / 45678 = 0.719).
- **Daniel Cook — depth/breadth from a BAM** — depth of coverage = `Sum of Depths / genome size`
  ("average number of reads aligned to an individual base"); breadth = `Bases Mapped / genome size`.
- **Illumina — Sequencing Coverage for NGS** — coverage = average number of reads covering known
  reference bases; Lander/Waterman **C = LN/G** (read length × read count / haploid genome length).

**Formula check.** `depth[i] = #{reads spanning i}`, `average = Σdepth / G`, `breadth =
#{i: depth[i] ≥ 1} / G` — each matches a source verbatim. `C = LN/G` cross-checked against the
per-base average on the worked set: L=5, N=3, G=10 → C = 15/10 = **1.5**, identical to Σdepth/G;
two independent methods agree. Edge semantics sourced: empty/no-aligned reads → all-zero depth,
average 0, breadth 0; empty reference → length-0 array; unaligned read (below `minOverlap`) → 0.

**Independent hand cross-check.** Reference `ACGTTGCAAT` (len 10); reads `ACGTT`@0 / `TTGCA`@3 /
`GCAAT`@5 (traced through `FindBestAlignment`). Depth `[1,1,1,2,2,2,2,2,1,1]`, Σ=15, average 1.5,
breadth 10/10 = 1.0 — matches spec M1/M7/M8 exactly.

**Note (why PASS-WITH-NOTES).** Sources describe a read overhanging the reference end
contributing only its overlapping portion (clip at the boundary). The repo's `FindBestAlignment`
scans only positions `0 … reference.Length − read.Length`, i.e. it places a read **only where it
fits entirely**. Consequently the `min(p + L, reference.Length)` clip in the depth loop is
genuinely **dead / defensive** code and the "overhang → partial contribution" case is
**unreachable** — an over-long read simply fails to place (depth 0, M4). This is an honest,
documented placement-model simplification (description §4.2, §5.3), **not** an arithmetic error;
the depth counting itself is exactly source-faithful.

## Stage B — implementation (code review)

Code path `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:772-819`.
Null `reference`/`reads` → `ArgumentNullException` (ThrowIfNull); allocates zeroed
`int[reference.Length]`; per read, `FindBestAlignment` returns the best position (≥ `minOverlap`
matches, leftmost on tie) or −1, then increments `[pos, min(pos+L, refLen))`. Over-long reads
return −1 (fit-only scan). `depth[i]` = count of placed reads spanning i, matching the validated
per-base definition; average/breadth are correctly left to the caller as Σ/len and covered/len.

**Cross-verification (test run vs sourced/hand oracles) — all ✅:** M1 depth
`[1,1,1,2,2,2,2,2,1,1]`; M7 Σ=15/avg 1.5; M8 breadth 1.0; M3 single read
`[0,0,0,0,0,1,1,1,1,1,0,0,0,0,0]`; M4 over-long → all 0; M5 unmatched → all 0; M6 empty reads →
all 0 len 10; S1 partial breadth/avg 0.5 → `[1,1,1,1,1,0,0,0,0,0]`; empty reference → `[]`
(new test).

**Variant/delegate consistency.** No public overload in `SequenceAssembler`; the MCP wrapper
`AlignmentTools.CalculateCoverage` delegates to it. `BedParser.CalculateCoverage` is a separate,
unrelated unit. No divergence.

**Test-quality audit (HARD gate: PASS).** M1/M7/M8/S1 assert exact arrays / `.Within(1e-10)`
values traced to the per-base definition + Cook/Metagenomics/Illumina formulas and hand
computation, not to code output; a deliberately-wrong impl (inclusive intervals, off-by-one,
wrong sum) would fail them. No Greater/AtLeast/Contains, no widened tolerances, no skips.
Coverage: all 8 MUST + 2 SHOULD + 1 COULD cases plus null-reference / null-reads validation.

## Findings

- **No algorithm defect.** State ✅ CLEAN.
- **One test-coverage gap fixed in-session** — the Stage-A-documented "empty reference → empty
  array" edge case (description §3.3) was untested; added
  `CalculateCoverage_EmptyReference_ReturnsEmptyArray`.
- **Carried note (not a defect)** — boundary-clip path is dead/defensive because
  `FindBestAlignment` places reads only where they fit entirely; the overhang partial-contribution
  case is unreachable. Intentional placement-model scope, recorded on [[coverage-depth-calculation]].
- **No follow-ups.**
