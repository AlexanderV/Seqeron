# Evidence Artifact: ANNOT-CODING-001

**Test Unit ID:** ANNOT-CODING-001
**Algorithm:** Coding Potential Calculation (CPAT hexamer usage-bias score)
**Date Collected:** 2026-06-13

---

## Online Sources

### CPAT paper (Wang et al. 2013) — Nucleic Acids Research

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3616698/
**Accessed:** 2026-06-13 (fetched with WebFetch; prompt asked for the hexamer usage-bias formula, log base, window step, averaging, and citation details)
**Authority rank:** 1 (peer-reviewed paper); the reference implementation it describes is rank 3.

**Key Extracted Points:**

1. **Hexamer usage bias feature:** CPAT measures the *"log-likelihood ratio to measure differential hexamer usage between coding and noncoding sequences."* Per-hexamer score is the logarithm of the ratio of in-frame hexamer frequencies under the coding model `F(hi)` and the noncoding model `F'(hi)`.
2. **Interpretation:** *"Hexamer score determines the relative degree of hexamer usage bias in a particular sequence. Positive values indicate a coding sequence, whereas negative values indicate a noncoding sequence."*
3. **Sequence-level score:** *"For a given hexamer sequence S = H1, H2, … , Hm"* the sequence score aggregates the individual per-hexamer log-likelihoods; the exact aggregation equation in the paper is rendered as an image (um3.jpg) and is not extractable as text — the precise form (sum divided by number of hexamers) was therefore taken from the reference implementation below.
4. **Prior work cited:** the hexamer coding measure references Fickett JW, Tung C-S (1992) "Assessment of protein coding measures," Nucleic Acids Res 20(24):6441–6450, DOI 10.1093/nar/20.24.6441.

### CPAT/lncScore reference implementation — FrameKmer.py (`kmer_ratio`, `word_generator`, `kmer_freq_file`)

**URL:** https://raw.githubusercontent.com/WGLab/lncScore/master/tools/cpmodule/FrameKmer.py
**Accessed:** 2026-06-13 (fetched verbatim with WebFetch and re-fetched via `curl` to read `word_generator`, `kmer_ratio`, `kmer_freq_file`, `all_possible_kmer`)
**Authority rank:** 3 (reference implementation; this file is the CPAT `cpmodule.FrameKmer` module reused by lncScore).

**Key Extracted Points:**

1. **In-frame window (`word_generator`):**
   ```python
   def word_generator(seq,word_size,step_size,frame=0):
       for i in xrange(frame,len(seq),step_size):
           word =  seq[i:i+word_size]
           if len(word) == word_size:
               yield word
   ```
   Start at `frame` (0), advance by `step_size`, keep only full-length words. CPAT uses `word_size = 6`, `step_size = 3` (in-frame hexamers on codon boundaries).
2. **Per-hexamer score and averaging (`kmer_ratio`, frame 0 branch):**
   ```python
   def kmer_ratio(seq,word_size,step_size,coding,noncoding):
       if len(seq) < word_size:
           return 0
       sum_of_log_ratio_0 = 0.0
       frame0_count=0.0
       for k in word_generator(seq=seq, word_size=word_size, step_size=step_size, frame=0):
           if (not coding.has_key(k)) or (not noncoding.has_key(k)):
               continue
           if coding[k]>0 and noncoding[k] >0:
               sum_of_log_ratio_0  +=  math.log( coding[k] / noncoding[k])
           elif coding[k]>0 and noncoding[k] == 0:
               sum_of_log_ratio_0 += 1
           elif coding[k] == 0 and noncoding[k] == 0:
               continue
           elif coding[k] == 0 and noncoding[k] >0 :
               sum_of_log_ratio_0 -= 1
           else:
               continue
           frame0_count += 1
       try:
           return sum_of_log_ratio_0/frame0_count
       except:
           return -1
   ```
   **(Corrected 2026-06-15 during ANNOT-CODING-001 validation:** the initial transcription
   omitted the `elif coding[k] == 0 and noncoding[k] == 0: continue` branch and trailing
   `else: continue`. Verified verbatim against the **canonical CPAT** repo `liguowang/cpat`
   `src/cpmodule/FrameKmer.py` (lines 95-96, fetched 2026-06-15) and the lncScore copy
   `WGLab/lncScore` `tools/cpmodule/FrameKmer.py` (lines 90-91). A both-zero in-both hexamer is
   `continue`d and **does NOT increment `frame0_count`** — skipped, not counted as 0.)
3. **Log base:** `math.log` is Python's natural logarithm (base e). Imports at top: `import math`.
4. **Missing-key handling:** a hexamer absent from either table is skipped (`continue`) and does not count toward `frame0_count`. Hexamers containing `N` are not in the tables (see point 6) so they are skipped.
5. **Empty / too-short input:** `len(seq) < word_size` returns 0; if no hexamer is scored, division by zero is caught and returns -1 (frame-0-only path); the public CPAT path returns `sum/count` for frame 0.
6. **Hexamer tables are raw counts (`kmer_freq_file`):** `count_table = Counter(word_generator(...))`; entries with `'N' in kmer` are skipped (`if 'N' in kmer: continue`). Values stored in the returned dict are integer counts. Therefore `kmer_ratio` divides whatever values the table holds — `coding[k] / noncoding[k]` — and the score is the mean of `ln(coding[k]/noncoding[k])`.

### EMBOSS tcode manual — Fickett (1982) TESTCODE (context / alternative measure)

**URL:** https://www.bioinformatics.nl/cgi-bin/emboss/help/tcode
**Accessed:** 2026-06-13 (via WebSearch result summarizing the EMBOSS tcode manual)
**Authority rank:** 3 (EMBOSS reference tool documentation), citing Fickett (1982) rank 1.

**Key Extracted Points:**

1. **Alternative coding measure:** Fickett JW (1982) "Recognition of protein coding regions in DNA sequences," Nucleic Acids Res 10(17):5303–5318. The EMBOSS `tcode` program implements the Fickett TESTCODE statistic using four position values and four composition values, combined as `MAX(A1,A2,A3)/MIN(A1,A2,A3) + 1` per position.
2. **Relevance:** documents that "coding potential" is a defined family of measures; CPAT's hexamer score and Fickett TESTCODE are distinct, both authoritative. The canonical method title in the Registry ("hexamer frequency bias") selects the CPAT hexamer score; TESTCODE is recorded as a related, not-implemented alternative.

---

## Documented Corner Cases and Failure Modes

### From FrameKmer.py (`kmer_ratio`)

1. **Sequence shorter than word size:** `if len(seq) < word_size: return 0`.
2. **No scorable hexamer:** all hexamers missing from a table → `frame0_count == 0` → division by zero is caught → reference returns -1 in the frame-0-only helper. (Implementation note: the C# port returns 0 here, matching the "no information" semantics of the empty/too-short case; see Assumptions.)
3. **Hexamer present in only one table:** skipped via the `has_key` guard (does not count).
4. **Coding>0, noncoding==0:** contributes +1 (pseudo-score, avoids division by zero / log of infinity).
5. **Coding==0, noncoding>0:** contributes −1 (pseudo-score).
6. **Coding==0, noncoding==0:** `continue` — the hexamer is skipped and **does NOT increment `frame0_count`** (it is not counted as a scored 0). Verified against canonical CPAT `liguowang/cpat` and lncScore (2026-06-15).

### From CPAT paper

1. **Sign convention:** positive score = coding, negative score = non-coding (a documented, testable invariant).

---

## Test Datasets

### Dataset: Hand-derived worked example (from `kmer_ratio` formula)

**Source:** Derived from FrameKmer.py `kmer_ratio` (Wang et al. 2013 reference implementation), https://raw.githubusercontent.com/WGLab/lncScore/master/tools/cpmodule/FrameKmer.py

| Parameter | Value |
|-----------|-------|
| Sequence | `ATGAAACCC` (length 9) |
| wordSize / stepSize | 6 / 3 |
| In-frame hexamers (frame 0) | i=0 `ATGAAA`, i=3 `AAACCC` (i=6 `CCC` length 3 < 6 → not yielded) |
| coding table | `{ATGAAA: 8, AAACCC: 2}` |
| noncoding table | `{ATGAAA: 2, AAACCC: 4}` |
| `ATGAAA` term | ln(8/2) = ln 4 = 1.3862943611198906 |
| `AAACCC` term | ln(2/4) = ln 0.5 = −0.6931471805599453 |
| sum | 0.6931471805599453 |
| count | 2 |
| **Expected score** | 0.6931471805599453 / 2 = **0.34657359027997264** |

### Dataset: Pseudo-score cases (from `kmer_ratio` branches)

**Source:** FrameKmer.py `kmer_ratio` elif branches, same URL.

| Parameter | Value |
|-----------|-------|
| coding-only hexamer | `ATGAAA`: coding=5, noncoding=0 → contributes +1 |
| noncoding-only hexamer | `ATGAAA`: coding=0, noncoding=5 → contributes −1 |
| single hexamer, count=1 | score = contribution / 1 = the contribution itself |

---

## Assumptions

1. **ASSUMPTION: Empty-table / no-scorable-hexamer return value.** The reference `kmer_ratio` returns -1 when `frame0_count == 0` (a caught `ZeroDivisionError`), but only after the loop runs; this -1 is a degenerate signal, not a true score. The C# port returns **0** in this case, consistent with the documented `len(seq) < word_size → 0` branch (both mean "no hexamer information available"). Changing -1↔0 only affects inputs with zero scorable hexamers, where the score is meaningless either way. Recorded as an assumption because the two reference branches disagree on the sentinel; behavior is otherwise verbatim.
2. **ASSUMPTION: Table units (counts vs proportions).** `kmer_freq_file` stores raw counts; CPAT's distributed prebuilt tables are normalized proportions. `kmer_ratio` divides the stored values directly, so the implementation is agnostic to units as long as both tables use the **same** units (the difference is an additive constant `ln(Σcoding/Σnoncoding)` on the score). The contract documents this; tests use small explicit tables so the expected value is exact regardless.

---

## Recommendations for Test Coverage

1. **MUST Test:** Mean log-ratio on a multi-hexamer sequence equals the hand-derived value (0.34657359027997264 for the dataset above) — Evidence: FrameKmer.py `kmer_ratio` + worked derivation.
2. **MUST Test:** Coding-only hexamer contributes +1; noncoding-only contributes −1 — Evidence: `kmer_ratio` elif branches.
3. **MUST Test:** Sign invariant — coding-biased tables give a positive score, noncoding-biased give negative — Evidence: CPAT paper interpretation.
4. **MUST Test:** Sequence shorter than wordSize returns 0 — Evidence: `kmer_ratio` guard.
5. **MUST Test:** Hexamer missing from either table is skipped (not counted) — Evidence: `has_key` guard.
6. **SHOULD Test:** In-frame stepping (step 3) and full-length-word filter — only hexamers on codon boundaries are scored — Rationale: `word_generator` semantics.
7. **SHOULD Test:** Null sequence / null table → `ArgumentNullException`; non-positive wordSize/stepSize → `ArgumentOutOfRangeException` — Rationale: documented validation contract.
8. **COULD Test:** Score is unit-invariant — multiplying both tables by a constant within each table leaves per-hexamer ratios unchanged — Rationale: ASSUMPTION 2.

---

## References

1. Wang L, Park HJ, Dasari S, Wang S, Kocher J-P, Li W (2013). CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model. Nucleic Acids Research 41(6):e74. https://doi.org/10.1093/nar/gkt006 (open-access copy: https://pmc.ncbi.nlm.nih.gov/articles/PMC3616698/)
2. CPAT / lncScore reference implementation, `cpmodule/FrameKmer.py` (`kmer_ratio`, `word_generator`, `kmer_freq_file`). https://github.com/WGLab/lncScore/blob/master/tools/cpmodule/FrameKmer.py (raw: https://raw.githubusercontent.com/WGLab/lncScore/master/tools/cpmodule/FrameKmer.py)
3. Fickett JW, Tung C-S (1992). Assessment of protein coding measures. Nucleic Acids Research 20(24):6441–6450. https://doi.org/10.1093/nar/20.24.6441
4. Fickett JW (1982). Recognition of protein coding regions in DNA sequences. Nucleic Acids Research 10(17):5303–5318. https://doi.org/10.1093/nar/10.17.5303 (EMBOSS tcode manual: https://www.bioinformatics.nl/cgi-bin/emboss/help/tcode)

---

## Change History

- **2026-06-13**: Initial documentation.
