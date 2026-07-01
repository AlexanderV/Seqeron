# Validation Report: SPLICE-MAXENT3-001 — MaxEntScan score3 (3' Acceptor)

- **Validated:** 2026-06-25   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.ScoreAcceptorMaxEnt(string window)`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
`ScoreAcceptorMaxEnt` (MaxEntScan score3ss — maximum-entropy 3' splice-acceptor model)

- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs` (lines ~1003–1227)
- **Embedded table:** `Data/maxent_score3.txt` (82 560 records, 9 matrices) + `Data/maxent_score3.LICENSE.md`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs` (region ME1–ME10)

## Authoritative sources opened (this session)
1. **Yeo G, Burge CB (2004)** "Maximum entropy modeling of short sequence motifs with applications
   to RNA splicing signals." *J Comput Biol* 11(2–3):377–394, DOI 10.1089/1066527041410418 — the
   original score3ss maximum-entropy model.
2. **maxentpy** (kepbod/maxentpy, MIT-licensed Python port) `maxent.py::score3` — reference
   factorisation, retrieved verbatim this session from
   `raw.githubusercontent.com/kepbod/maxentpy/master/maxentpy/maxent.py`.
3. Web confirmation of the documented conventions and reference value
   (23-mer = 20 intron + 3 exon; `score3('ttccaaacgaacttttgtAGgga') = 2.8867730651152104`).

## Stage A — Description

**Model conventions (confirmed):**
- Window length **23 nt = 20 intronic + 3 exonic**, with the conserved **AG at 0-based positions 18–19**.
- Score = **log2( P_maxent(window) / P_background(window) )** in bits; higher = stronger acceptor.
- The conserved AG dinucleotide is removed and scored by a consensus/background term;
  the remaining **21-nt "rest"** sequence's maximum-entropy probability is factorised over
  nine overlapping sub-sequences (inclusion–exclusion over 5 numerator + 4 denominator marginals).

**Factorisation (maxentpy `score3`, verbatim):**
```
key   = fa[18:20];  score = cons1_3[key0]*cons2_3[key1] / (bgd_3[key0]*bgd_3[key1])
rest  = fa[:18] + fa[20:]
rest_score *= m0[rest[0:7]]  *= m1[rest[7:14]] *= m2[rest[14:21]] *= m3[rest[4:11]] *= m4[rest[11:18]]
rest_score /= m5[rest[4:7]]  /= m6[rest[7:11]] /= m7[rest[11:14]] /= m8[rest[14:18]]
return log2(score * rest_score)
```
- `bgd_3 = {A:.27, C:.23, G:.23, T:.27}`; `cons1_3 = {A:.9903, C:.0032, G:.0034, T:.0030}`;
  `cons2_3 = {A:.0027, C:.0037, G:.9905, T:.0030}`. T==U; hashseq = base-4 (A0 C1 G2 T3).

**Edge-case semantics (sourced):**
- Window length ≠ 23 → invalid input (maxentpy `sys.exit('Wrong length of fa!')`).
- Non-A/C/G/T(/U) character → invalid (no hash code).
- Lowercase / DNA(T) vs RNA(U) → upper-cased, T≡U; identical scores.
- **No AG check:** maxentpy `score3` does **not** require AG at 18–19. A non-AG dinucleotide is
  *scored* (heavily penalised via the consensus term), not rejected. Confirmed by oracle:
  `…CCgga` → −13.220039, `…TTgga` → −14.078362 (vs +2.886773 for the AG window).

**Independent cross-check (numbers):** an independent Python reimplementation of maxentpy `score3`
(separate code path, reading the embedded table) reproduced the three documented worked examples
**exactly**:

| window | oracle score3 | round 2 dp |
|--------|---------------|-----------|
| `ttccaaacgaacttttgtAGgga` | 2.886773 | 2.89 |
| `tgtctttttctgtgtggcAGtgg` | 8.190965 | 8.19 |
| `ttctctcttcagacttatAGcaa` | −0.080278 | −0.08 |

**Stage A findings:** none. Description and embedded table are faithful to Yeo & Burge (2004) /
maxentpy. The repo's TestSpec/report stubs cited the right sources; the description is correct.

## Stage B — Implementation

**Code path reviewed:** `ScoreAcceptorMaxEnt` (SpliceSitePredictor.cs:1176–1217), `HashMaxEntSubsequence`
(1128–1146), `LoadMaxEntScore3Tables` (1094–1124), factor tables (1053–1068), AG/bgd dicts (1031–1046).

**Formula realised correctly?** Yes — line-by-line equivalent to maxentpy `score3`:
- AG term at `upper[18],upper[19]` = `cons1·cons2/(bgd·bgd)` ✓
- rest = `window[0..18) + window[20..23)` = 21 nt ✓
- numerator m0/rest[0:7], m1/[7:14], m2/[14:21], m3/[4:11], m4/[11:18] ✓
- denominator m5/[4:7], m6/[7:11], m7/[11:14], m8/[14:18] ✓
- `Math.Log2(consScore * restScore)` ✓; hashseq base-4 A0 C1 G2 T3, T≡U ✓.
- Embedded `maxent_score3.txt` = 82 560 records, matrix sizes {16384×5, 64, 256, 64, 256} match the
  documented factorisation exactly.

**Cross-verification table recomputed vs code** (differential test, C# `ScoreAcceptorMaxEnt`
vs the independent oracle, all to <5e-7):

| window | oracle | C# | match |
|--------|--------|----|-------|
| `ttccaaacgaacttttgtAGgga` | 2.886773 | 2.886773 | ✓ |
| `tgtctttttctgtgtggcAGtgg` | 8.190965 | 8.190965 | ✓ |
| `ttctctcttcagacttatAGcaa` | −0.080278 | −0.080278 | ✓ |
| `AAAAAAAAAAAAAAAAAAAGAAA` | −7.514732 | −7.514732 | ✓ |
| `cccccccccccccccccaAGttt` | 6.836558 | 6.836558 | ✓ |
| `TTTTTTTTTTTTTTTTTTAGGTC` | 13.774543 | 13.774543 | ✓ |
| `ttccaaacgaacttttgtCCgga` (non-AG) | −13.220039 | −13.220039 | ✓ |
| `ttccaaacgaacttttgtTTgga` (non-AG) | −14.078362 | −14.078362 | ✓ |

**Variant/delegate consistency:** single public entry point; no `*Fast`/instance variant for score3.
(score5ss is a separate region/unit.)

**Numerical robustness:** products of probabilities stay well within `double` range; the only
risk is a missing hash key (would `KeyNotFoundException`), but every A/C/G/T 7-/4-/3-mer is present
in the full 82 560-record table, so no lookup can miss for valid input. Invalid characters throw
before any lookup.

**Test quality audit:** ME1–ME9 pin exact sourced values (2.89 / 8.19 / −0.08 + their full-precision
forms, traced to maxentpy not to code), cover length-≠-23, non-ACGT, null, DNA≡RNA, case, and
ordering. Added **ME10** to lock the Stage-A "non-AG dinucleotide is scored not rejected" semantic
with oracle-derived values (−13.220039 / −14.078362). No green-washing; all expected values trace
to Yeo & Burge / maxentpy, not code echoes.

**Stage B findings:** none. Implementation is a faithful, exact reproduction of the reference.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- C# reproduces MaxEntScan score3 exactly on every tested window (8/8, <5e-7).
- Test added: ME10 (non-AG window scored, not rejected). Full suite: `Seqeron.Genomics.Tests`
  18762 passed / 0 failed; solution-wide Failed: 0.
- No defects logged.
