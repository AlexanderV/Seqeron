# Evidence Artifact: ASSEMBLY-TRIM-001

**Test Unit ID:** ASSEMBLY-TRIM-001
**Algorithm:** Quality Trimming (BWA / cutadapt running-sum)
**Date Collected:** 2026-06-13

---

## Online Sources

### Cutadapt — Algorithm details (Quality trimming)

**URL:** https://cutadapt.readthedocs.io/en/stable/algorithms.html
**Accessed:** 2026-06-13 (also cross-checked v1.18: https://cutadapt.readthedocs.io/en/v1.18/algorithms.html)
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **Algorithm identity:** "The trimming algorithm implemented in Cutadapt is the same as the one used by BWA", applied to both ends of the read in turn if requested. Retrieved via WebSearch ("cutadapt quality trimming algorithm BWA running sum threshold formula") then WebFetch of the algorithms page.
2. **Core procedure (verbatim):** "Subtract the given cutoff from all qualities; compute partial sums from all indices to the end of the sequence; cut the sequence at the index at which the sum is minimal."
3. **Both ends (verbatim):** "If both ends are to be trimmed, repeat this for the other end."
4. **Intent (verbatim):** "The basic idea is to remove all bases starting from the end of the read whose quality is smaller than the given threshold. This is refined a bit by allowing some good-quality bases among the bad-quality ones."
5. **Worked example (verbatim):** qualities `42, 40, 26, 27, 8, 7, 11, 4, 2, 3`, threshold `10`. After subtracting threshold: `32, 30, 16, 17, -2, -3, 1, -6, -8, -7`. Partial sums from the end: `(70), (38), 8, -8, -25, -23, -20, -21, -15, -7`. "The position of the minimum (-25) is used as the trimming position. Therefore, the read is trimmed to the first four bases, which have quality values 42, 40, 26, 27."

### BWA source — `bwa_trim_read` (bwaseqio.c)

**URL:** https://github.com/lh3/bwa/blob/master/bwaseqio.c (raw: https://raw.githubusercontent.com/lh3/bwa/master/bwaseqio.c)
**Accessed:** 2026-06-13
**Authority rank:** 3 (canonical reference implementation; Heng Li)

**Key Extracted Points:**

1. **Function (verbatim C):**
   ```c
   int bwa_trim_read(int trim_qual, bwa_seq_t *p)
   {
       int s = 0, l, max = 0, max_l = p->len;
       if (trim_qual < 1 || p->qual == 0) return 0;
       for (l = p->len - 1; l >= BWA_MIN_RDLEN; --l) {
           s += trim_qual - (p->qual[l] - 33);
           if (s < 0) break;
           if (s > max) max = s, max_l = l;
       }
       p->clip_len = p->len = max_l;
       return p->full_len - p->len;
   }
   ```
2. **Phred+33 decoding:** quality byte is decoded as `p->qual[l] - 33` (Sanger/Phred+33 ASCII offset).
3. **No-trim guard:** `if (trim_qual < 1 ...) return 0;` — a threshold below 1 disables trimming.
4. **Accumulator:** `s += trim_qual - (qual - 33)` from the 3' end; tracks the argmax position `max_l`; BWA additionally early-breaks when `s < 0` and floors at `BWA_MIN_RDLEN`. The argmax of accumulated `(threshold - q)` from the end is equivalent to cutadapt's argmin of partial sums of `(q - threshold)`.

### BWA source — `BWA_MIN_RDLEN` (bwtaln.h)

**URL:** https://github.com/lh3/bwa/blob/master/bwtaln.h
**Accessed:** 2026-06-13
**Authority rank:** 3

**Key Extracted Points:**

1. **Constant (verbatim):** `#define BWA_MIN_RDLEN 35 // for read trimming`. Retrieved via WebSearch then WebFetch. This is BWA's hard floor on trimmed read length, separate from the running-sum optimum.

### Cock et al. (2010) — Sanger FASTQ format (Phred+33)

**URL:** https://academic.oup.com/nar/article/38/6/1767/3112533 (also PMC/Strathprints PDF mirror)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper, Nucleic Acids Research)

**Key Extracted Points:**

1. **Sanger encoding (verbatim):** "Sanger FASTQ files use ASCII 33–126 to encode PHRED qualities from 0 to 93 (i.e. PHRED scores with an ASCII offset of 33)." Retrieved via WebFetch of the NAR article page.
2. **Phred definition:** PHRED quality score `Q` is defined in terms of the estimated error probability `P` as `Q = -10 · log10(P)` (the standard PHRED definition; the article presents the equation as a figure image). Used only as background; the trimming algorithm operates on integer Phred values directly.

---

## Documented Corner Cases and Failure Modes

### From Cutadapt / BWA

1. **Threshold disables trimming:** BWA returns 0 (no trim) when `trim_qual < 1`. A threshold of 0 means subtract 0 → all partial sums non-negative → minimum at the end → nothing trimmed.
2. **All-high-quality read:** every `q - threshold ≥ 0`, partial sums are all ≥ 0 and minimal at the last index → no bases removed.
3. **All-low-quality read:** every `q - threshold < 0`, partial sums strictly decreasing toward the front from the 3' side; trimming removes the entire 3' span and (for the 5' pass) the entire 5' span → read may be fully removed.
4. **Good base among bad ones:** a single high-quality base inside a low-quality tail is retained only if the cumulative sum does not reach a new minimum before it — exactly the "refinement" described by cutadapt.

---

## Test Datasets

### Dataset: Cutadapt worked example

**Source:** Cutadapt algorithm docs (verbatim, see Online Sources #5)

| Parameter | Value |
|-----------|-------|
| Qualities (Phred) | 42, 40, 26, 27, 8, 7, 11, 4, 2, 3 |
| Phred+33 ASCII string | `KI;<)(,%#$` |
| Threshold | 10 |
| Partial sums from end (q-cutoff) | (70), (38), 8, -8, -25, -23, -20, -21, -15, -7 |
| Minimum value | -25 (index 4, 0-based) |
| Expected 3' trim result | first 4 bases kept (qualities 42, 40, 26, 27) |

ASCII derivation (Phred+33): 42→`K`, 40→`I`, 26→`;`, 27→`<`, 8→`)`, 7→`(`, 11→`,`, 4→`%`, 2→`#`, 3→`$`.

---

## Assumptions

1. **ASSUMPTION: 5'/3' both-end pass order.** Cutadapt trims both ends "in turn"; the exact order (5' first vs 3' first) is not numerically significant because the two passes operate on disjoint ends of the surviving window. We implement 3'-end trim then 5'-end trim on the remaining window. This does not change which bases survive for the published example (a pure 3' case) and is consistent with cutadapt's "repeat for the other end" wording.
2. **ASSUMPTION: `minLength` filter semantics.** The Test Unit signature includes `minLength`; BWA/cutadapt's running-sum core does not itself define a minimum-length drop, but read filtering by minimum length after trimming is standard (cutadapt `--minimum-length`). Reads whose trimmed length `< minLength` are dropped from the output. This is a documented downstream filter, not part of the running-sum optimum.

---

## Recommendations for Test Coverage

1. **MUST Test:** Cutadapt worked example — qualities `KI;<)(,%#$` (Phred 42,40,26,27,8,7,11,4,2,3), threshold 10 → trimmed to first 4 bases. — Evidence: Cutadapt algorithm docs (Online Sources #5).
2. **MUST Test:** Phred+33 decoding — `q = ASCII - 33` per Cock et al. (2010). — Evidence: Cock et al. 2010.
3. **MUST Test:** All-high-quality read → unchanged. — Evidence: partial-sum minimum at end (cutadapt core).
4. **MUST Test:** All-low-quality read → dropped (length 0 < minLength). — Evidence: cutadapt core + min-length filter.
5. **SHOULD Test:** 5'-end trimming (low-quality prefix). — Rationale: cutadapt "repeat for the other end".
6. **SHOULD Test:** `minLength` filter drops short survivors; keeps survivors ≥ minLength. — Rationale: standard min-length filter.
7. **SHOULD Test:** threshold ≤ 0 disables trimming (BWA `trim_qual < 1`). — Rationale: BWA guard.
8. **COULD Test:** "good base among bad ones" retention behaviour. — Rationale: cutadapt refinement note.
9. **COULD Test:** empty / null inputs and empty quality string. — Rationale: defensive contract.

---

## References

1. Marcel, M. (Cutadapt project). 2024. Algorithm details — Cutadapt documentation. https://cutadapt.readthedocs.io/en/stable/algorithms.html
2. Li, H. BWA source: `bwa_trim_read` in bwaseqio.c. https://github.com/lh3/bwa/blob/master/bwaseqio.c ; `BWA_MIN_RDLEN` in bwtaln.h. https://github.com/lh3/bwa/blob/master/bwtaln.h
3. Cock, P.J.A., Fields, C.J., Goto, N., Heuer, M.L., Rice, P.M. (2010). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research 38(6):1767–1771. https://doi.org/10.1093/nar/gkp1137

---

## Change History

- **2026-06-13**: Initial documentation.
