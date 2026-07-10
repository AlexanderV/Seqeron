---
type: source
title: "Validation report: ASSEMBLY-TRIM-001 (quality trimming — BWA/cutadapt running-sum)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-TRIM-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-TRIM-001.md
source_commit: 9e3afa723f48771feef69632da397e2992f74114
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-TRIM-001

The two-stage **validation write-up** for test unit **ASSEMBLY-TRIM-001** (quality trimming — the
BWA/cutadapt running-sum that removes low-quality bases from read ends before assembly), validated
2026-06-15 in a fresh context. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B). The running-sum core, the BWA argmax equivalence,
Phred+33 decoding and the published oracle are summarized in the concept
[[quality-trimming-running-sum]] (anchor of the assembly TRIM family), and the wider campaign is
[[validation-and-testing]]. Distinct from [[assembly-trim-001-evidence]] (the pre-implementation
evidence artifact, sourced from `docs/Evidence/`) — this is the independent re-validation verdict, and
it found and fixed a real defect the evidence pass did not surface.

Canonical methods under test:
`SequenceAssembler.QualityTrimReads(reads, minQuality, minLength)` plus the private `TrimEnd` /
`TrimStart` (`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:950`).

## Verdict

**Stage A: ✅ PASS · Stage B: ❌ FAIL → 🔧 FIXED · End-state: ✅ CLEAN.** A real algorithm defect was
found and completely fixed in-session; code now matches the cutadapt `quality_trim_index` reference
(**0 mismatches over 900k random reads**), tests lock the sourced behaviour, full unfiltered suite
**Passed: 6535, Failed: 0** (build 0 errors; the 4 warnings are pre-existing NUnit2007 in an unrelated
`ApproximateMatcher_EditDistance_Tests.cs`).

## Stage A — description (algorithm faithfulness): PASS

Theory re-checked against primary sources opened live this session:

- **cutadapt "Algorithm details"** — verbatim: subtract the cutoff from all qualities, compute partial
  sums from each index to the end, cut at the index where the sum is minimal; "repeat this for the
  other end." Worked example `42,40,26,27,8,7,11,4,2,3` @ threshold 10 → partial sums from the end
  `(70),(38),8,-8,-25,-23,-20,-21,-15,-7`, minimum −25 → **first four bases kept**.
- **cutadapt `quality_trim_index` (qualtrim.pyx)** — the *authoritative* implementation and the key
  reference for this pass: **two independent passes over the full read** — 5' (`s += cutoff - (q-base)`,
  `if s<0: break`, track `start=i+1` on a new max) and 3' (reversed, track `stop=i`) — then
  **`if start >= stop: start, stop = 0, 0`** (drop the read).
- **BWA `bwa_trim_read` (bwaseqio.c)** — the argmax loop, `qual - 33` (Phred+33) decode, the
  `trim_qual < 1 || qual == 0` no-trim guard, **the `if (s < 0) break;` early break**, and the
  `BWA_MIN_RDLEN` floor.
- **Cock et al. (2010) NAR 38(6):1767–1771** — Sanger FASTQ = ASCII offset 33 (Phred+33).

**Formula check.** The maths is correct and sourced: argmin of partial sums of `(q − cutoff)` from each
end **is equivalent** to BWA's argmax of `(cutoff − q)` (sign-flip). The cutadapt reference additionally
carries (a) the `s < 0` early break and (b) the `start ≥ stop ⇒ (0,0)` drop. Stage A is **PASS** — the
only issue was *description wording*: §4.1 said the 5' pass runs "over the surviving window" and §5.3
said the early break was "replaced by the global-minimum formulation." That text described the *buggy
code*, not the cited algorithm, and was corrected this session to match cutadapt `quality_trim_index`.

## Stage B — implementation (code review): FAIL → FIXED

**Defect found (real, fixed).** The original `TrimEnd`/`TrimStart` computed the **global minimum**
partial sum over the window with **no `s < 0` early break**, and chained the 5' pass onto the
3'-surviving window. This is NOT the cutadapt/BWA algorithm. Three consequences:

1. **Missing `s < 0` early break** — a high-quality base near an end should *protect* the interior from
   that end's pass (the documented "allow some good-quality bases among the bad" refinement). Without it
   the code finds the global minimum anywhere and over-trims.
2. **5' pass chained on the 3'-surviving window** instead of the full read `[0, n)`.
3. **No `start ≥ stop ⇒ drop` rule.**

**Minimal repro:** read qualities `20,0,8,40` (`"5!)I"`), cutoff 10. cutadapt: both passes break
immediately → **keep all 4 bases**. Old code: returned `(3,4)` → kept only the **last** base. A
randomized comparison of the old code vs the cutadapt reference showed **~19.5% of random reads
diverged** (38,989 / 200,000). Notably **every pre-existing test (M1–M5, S1–S3) agreed with cutadapt** —
all happy-path cases that never trigger the early break, so the suite shared the implementation's blind
spot (exactly the failure mode the [[validation-protocol]] warns about).

**Fix applied:** added the `s < 0` early break to both `TrimEnd`/`TrimStart` (sign-flipped `if (sum > 0)
break;`); both passes now run independently over the full read `[0, n)`; added the `start ≥ stop ⇒
(0,0)` drop rule in `QualityTrimReads`; updated `docs/algorithms/Assembly/Quality_Trimming.md` to match
the cutadapt reference. **Verification:** a Python port of the fixed algorithm vs cutadapt
`quality_trim_index` over **300,000 random reads × 3 cutoffs (1, 10, 20)**, Phred 0–93: **0 mismatches**.

**Cross-verification (recomputed vs cutadapt reference):**

| Read (Phred) | ASCII | cutoff | (start,stop) | kept |
|---|---|---|---|---|
| 42,40,26,27,8,7,11,4,2,3 | `KI;<)(,%#$` | 10 | (0,4) | first 4 |
| 40×6 | `IIIIII` | 20 | (0,6) | all |
| 0×6 | `!!!!!!` | 20 | (0,0) | dropped |
| 3,2,4,11,7,8,27,26,40,42 | `$#%,()<;IK` | 10 | (6,10) | last 4 |
| 20,0,8,40 | `5!)I` | 10 | (0,4) | **all 4 (refinement)** |
| 2,2,40,2,2 | `##I##` | 10 | (2,3) | only the good base |
| 20,5,5 | `5&&` | 20 | (0,0) | dropped (start≥stop) |

**Test-quality gate: PASS.** Three new tests lock the sourced behaviour and each fails against the old
(no-break) code: `GoodBaseAt3PrimeEnd_KeepsLowQualityInterior` (`5!)I`→`ACGT`),
`IsolatedGoodBaseInBadRead_KeepsOnlyGoodBase` (`##I##`→`G`), `WindowsCross_DropsRead` (`5&&`→dropped);
expected values traced to `quality_trim_index`, not to code output. All assertions exact, no skips,
honest green.

## Findings

- **One real defect (missing `s < 0` early break + 5'-on-surviving-window chaining + missing
  `start≥stop` drop), fully fixed. End-state ✅ CLEAN**; code conformant to BWA/cutadapt
  (0/900k mismatches).
- **No source contradictions** — cutadapt explicitly identifies its algorithm with BWA's; the argmax
  and argmin forms are algebraically equal; Cock et al. supply the Phred+33 encoding both rely on.
- **Concept reconciliation (flagged):** the corrected algorithm runs **both passes independently over
  the full read** (not a 5' pass on the 3'-surviving window) and adds the `start ≥ stop ⇒ drop` rule.
  The concept [[quality-trimming-running-sum]] originally described the pre-fix behaviour; it has been
  updated to match this report.
- **Follow-ups:** none.
