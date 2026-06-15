# Validation Report: ASSEMBLY-TRIM-001 — Quality Trimming (BWA / cutadapt running-sum)

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.QualityTrimReads(reads, minQuality, minLength)` (+ private `TrimEnd` / `TrimStart`)
- **Stage A verdict:** PASS
- **Stage B verdict:** FAIL → FIXED (defect found and completely fixed this session)
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session
- **Cutadapt — Algorithm details** (https://cutadapt.readthedocs.io/en/stable/algorithms.html), fetched 2026-06-15. Confirms verbatim: "Subtract the given cutoff from all qualities; compute partial sums from all indices to the end of the sequence; cut the sequence at the index at which the sum is minimal." Worked example qualities `42,40,26,27,8,7,11,4,2,3`, threshold `10` → after subtraction `32,30,16,17,-2,-3,1,-6,-8,-7`; partial sums from the end `(70),(38),8,-8,-25,-23,-20,-21,-15,-7`; minimum −25 → **first four bases kept**. "If both ends are to be trimmed, repeat this for the other end."
- **BWA `bwa_trim_read`** (https://raw.githubusercontent.com/lh3/bwa/master/bwaseqio.c), fetched 2026-06-15. Confirms the loop, the `qual[l] - 33` (Phred+33) decoding, the `if (trim_qual < 1 || p->qual == 0) return 0;` guard, **and the `if (s < 0) break;` early break** plus the `BWA_MIN_RDLEN` floor.
- **Cutadapt `quality_trim_index`** reference source (https://raw.githubusercontent.com/marcelm/cutadapt/main/src/cutadapt/qualtrim.pyx), fetched 2026-06-15. The authoritative implementation. Two independent passes over the **full** read: 5' (`s += cutoff_front - (qual[i]-base)`, `if s<0: break`, track `start=i+1` on new max) and 3' (reversed, track `stop=i`), then **`if start >= stop: start, stop = 0, 0`**.
- **Cock et al. (2010) NAR 38(6):1767–1771** — Sanger FASTQ = ASCII offset 33 (cited in repo Evidence; Phred+33 confirmed).

### Formula check
The repo description (§2.2, §4.1) quotes the cutadapt procedure and the BWA loop (which *includes* `if (s < 0) break;`). The maths is correct: argmin of partial sums of `(q − cutoff)` from each end **is equivalent** to BWA's argmax of `(cutoff − q)` — sign-flip. The cutadapt reference additionally has (a) the `s < 0` early break and (b) the `start >= stop ⇒ (0,0)` drop. The description text was *mostly* correct but its §4.1 step-3 said the 5' pass runs "over the surviving window" and §5.3 claimed the early break was "replaced by the global-minimum formulation" — that wording described the (buggy) code, not the cited algorithm. Corrected this session.

### Edge-case semantics
Threshold < 1 disables trimming (BWA guard) ✓. All-high → unchanged ✓. All-low → dropped ✓. Phred+33 ✓. The "good base among bad ones" refinement is precisely the `s < 0` early break.

### Independent cross-check (numbers)
Re-derived the cutadapt worked example by hand and via a Python port of `quality_trim_index`: (start, stop) = (0, 4), keep `42,40,26,27`. Matches the docs exactly.

### Findings / divergences
Stage A **PASS**. The biology/maths is correct and sourced. Minor description wording that mischaracterised the algorithm (omitting the early break, describing window-chaining) was corrected to match the cutadapt reference.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:950` (`QualityTrimReads`), `TrimEnd`, `TrimStart`.

### Defect found (DEFECT — FIXED)
The original `TrimEnd`/`TrimStart` computed the **global minimum** partial sum over the window with **no `s < 0` early break**, and chained the 5' pass onto the 3'-surviving window. This is NOT the cutadapt/BWA algorithm. Consequences:

1. **Missing `s < 0` early break.** A high-quality base near an end should *protect* the interior bases from that end's pass (the documented "allow some good-quality bases among the bad ones" refinement). Without the break the code finds the global minimum anywhere in the read and over-trims.
2. **5' pass chained on the 3'-surviving window** instead of the full read.
3. **No `start >= stop ⇒ drop` rule.**

**Minimal repro (hand-traced + cutadapt reference):** read qualities `20,0,8,40` (ASCII `"5!)I"`), cutoff 10. cutadapt: both passes break immediately (`10-40<0`, `10-20<0`) → keep all 4 bases. Old code: 3' pass kept all, 5' pass found min at index 2 → returned `(3,4)` → kept only the **last** base. A randomized comparison of the old code vs the cutadapt reference (`quality_trim_index`) showed **~19.5% of random reads diverged** (38,989 / 200,000).

Notably **every pre-existing test (M1–M5, S1–S3) agreed with cutadapt** — they are happy-path cases that never trigger the early break, so the suite shared the implementation's blind spot (exactly the failure mode the protocol warns about).

### Fix applied
- Added the `s < 0` early break to both `TrimEnd` and `TrimStart` (sign-flipped to `if (sum > 0) break;`).
- Both passes now run independently over the full read `[0, n)`.
- Added the `start >= stop ⇒ (0,0)` drop rule in `QualityTrimReads`.
- Updated `docs/algorithms/Assembly/Quality_Trimming.md` (§4.1, §5.2, §5.3, §2.5, Deviations) to match the cutadapt reference exactly.

### Verification of the fix
A Python port of the fixed algorithm was compared against the cutadapt `quality_trim_index` reference over **300,000 random reads × 3 cutoffs (1, 10, 20)** with Phred values spanning 0–93: **0 mismatches**.

### Cross-verification table (recomputed vs cutadapt reference)
| Read (Phred) | ASCII | cutoff | cutadapt (start,stop) | repo (fixed) | kept |
|---|---|---|---|---|---|
| 42,40,26,27,8,7,11,4,2,3 | `KI;<)(,%#$` | 10 | (0,4) | (0,4) | first 4 |
| 40×6 | `IIIIII` | 20 | (0,6) | (0,6) | all |
| 0×6 | `!!!!!!` | 20 | (0,0) | (0,0) | dropped |
| 3,2,4,11,7,8,27,26,40,42 | `$#%,()<;IK` | 10 | (6,10) | (6,10) | last 4 |
| 20,0,8,40 | `5!)I` | 10 | (0,4) | (0,4) | all 4 (refinement) |
| 2,2,40,2,2 | `##I##` | 10 | (2,3) | (2,3) | only good base |
| 20,5,5 | `5&&` | 20 | (0,0) | (0,0) | dropped (start≥stop) |

### Test quality audit / gate
- **Sourced expectations, not code echoes:** PASS. New tests `GoodBaseAt3PrimeEnd_KeepsLowQualityInterior` (`5!)I`→`ACGT`), `IsolatedGoodBaseInBadRead_KeepsOnlyGoodBase` (`##I##`→`G`), `WindowsCross_DropsRead` (`5&&`→dropped) each fail against the old (no-break) code and pass against the cutadapt-conformant code; expected values traced to `quality_trim_index`.
- **No green-washing:** PASS. All assertions are exact (`Is.EqualTo`, `Is.Empty`); no weakening, no skips.
- **Cover all logic:** PASS. M1–M5/S1–S3/C1–C3 (existing) + the three new refinement/break/cross cases now exercise the early-break branch, the `start>=stop` branch, the cutoff<1 guard, min-length boundaries, multiple reads, and null/empty/empty-seq edges.
- **Honest green:** PASS. Full unfiltered suite: **Failed: 0, Passed: 6535**. Build: 0 errors; changed files warning-free (the 4 build warnings are pre-existing NUnit2007 in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

## Verdict & follow-ups
- **Stage A: PASS** (description corrected for wording; algorithm sourced and correct).
- **Stage B: FAIL → FIXED.** Defect (missing `s < 0` early break, window chaining, missing `start>=stop` rule) was real and is now completely fixed; code matches cutadapt `quality_trim_index` (0/900k mismatches), tests lock the sourced behaviour, full suite green.
- **End-state: ✅ CLEAN.** Algorithm fully functional and conformant to BWA/cutadapt.
