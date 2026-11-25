# Evidence Artifact: ASSEMBLY-MERGE-001

**Test Unit ID:** ASSEMBLY-MERGE-001
**Algorithm:** Contig Merging (suffix–prefix overlap merge / superstring collapse)
**Date Collected:** 2026-06-13

---

## Online Sources

### Langmead B., "Assembly & shortest common superstring" (JHU lecture notes)

**URL:** https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer/textbook-grade course notes; cites Gusfield, CLRS, Dasgupta-Papadimitriou-Vazirani)

**How retrieved:** WebFetch of the URL above returned the binary PDF (saved locally); text was
extracted with `pdftotext -layout` and read directly. Extracted points are copied from that text.

**Key Extracted Points:**

1. **Overlap definition:** "Overlap: length-l suffix of X matches length-l prefix of Y, where l is
   given." An overlap exists "when a suffix of X of length ≥ l exactly matches a prefix of Y."
2. **Merge primitive (`suffixPrefixMatch`):** the helper returns "length of longest suffix of x of
   length at least k that matches a prefix of y. Return 0 if there [is] no suffix/prefix match
   [of] length at least k." The returned value is the overlap length used to collapse the two strings.
3. **Superstring = concatenation with overlap removed:** "Without requirement of 'shortest', it's
   easy: just concatenate them." Example S = {BAA, AAB, BBA, ABA, ABB, BBB, AAA, BAB},
   plain concatenation = `BAAAABBBAABAABBBBBAAABAB` (length 24); collapsing overlaps yields a
   shorter superstring. The collapse keeps a single copy of the overlapping region.
4. **Greedy merge worked trace (l = 1):** input strings `ABA ABB AAA AAB BBB BBA BAB BAA`; the first
   round merges `BAA` + `AAB` (suffix `AA` = prefix `AA`, overlap length 2) into `BAAB` (length 4 =
   3 + 3 − 2). "Number in first column = length of overlap merged before that round."
5. **Greedy small example:** S = {AAA, AAB, ABB, BBB, BBA}, SCS = `AAABBBA`; the alignment
   AAA / AAB / ABB / BBB / BBA each overlap the next by 2 and collapse to the 7-char superstring.

### Langmead B., "Overlap Layout Consensus assembly" (JHU lecture notes)

**URL:** https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1

**How retrieved:** WebFetch of the URL above returned the binary PDF (saved locally); text extracted
with `pdftotext -layout` and read directly.

**Key Extracted Points:**

1. **Contig terminology:** "Fragments are contigs (short for contiguous)." OLC stages: "Overlap –
   Build overlap graph; Layout – Bundle stretches of the overlap graph into contigs; Consensus –
   Pick most likely nucleotide sequence for each contig."
2. **Overlap is a suffix/prefix match:** "an overlap is a suffix/prefix match … a directed edge is an
   overlap between suffix of source and prefix of sink."
3. **Longest match is the reported overlap:** "Assume for given string pair we report only the
   longest suffix/prefix match." (suffix-tree section) — i.e. when collapsing, the overlap length is
   the longest suffix-of-source/prefix-of-sink match.
4. **Suffix-tree applicability for overlaps:** overlaps "involving a prefix of x and a suffix of
   another string y" can be found by building a generalized suffix tree; "Time to build generalized
   suffix tree: O(N) … Overall: O(N + a)." This is the overlap *discovery* step, not the *merge* step.

### MIT 7.91J Foundations of Computational and Systems Biology, Lecture 6 (Spring 2014)

**URL:** https://ocw.mit.edu/courses/7-91j-foundations-of-computational-and-systems-biology-spring-2014/e885f0eb376ea6c2045eb9d8847f106f_MIT7_91JS14_Lecture6.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1 (university course, MIT OCW)

**How retrieved:** WebFetch of the URL returned the binary PDF (saved locally); text extracted with
`pdftotext -layout`. Used as an independent corroboration of the overlap/merge definition.

**Key Extracted Points:**

1. **Same overlap definition:** "an 'overlap' is when a suffix of X of length ≥ l exactly matches a
   prefix of Y." Overlap graph: "a directed edge is an overlap between suffix of source and prefix
   of sink."
2. **OLC consensus stage:** "Consensus – Pick most likely nucleotide sequence for each contig";
   contigs are built by bundling overlapping reads. Corroborates the suffix/prefix collapse model.

---

### Langmead B., comp-genomics-class notebook `CG_SCS.ipynb` (reference code)

**URL:** https://raw.githubusercontent.com/BenLangmead/comp-genomics-class/master/notebooks/CG_SCS.ipynb
**Accessed:** 2026-11-25 (B10 review; neutral User-Agent)
**How retrieved:** curl of the raw notebook; code cells extracted with `json`.

**Key Extracted Points (verbatim code):**

1. `overlap(a, b, min_length)` returns "length of longest suffix of 'a' matching a prefix of 'b' that is at
   least 'min_length' characters long. If no such overlap exists, return 0" — so the overlap length passed
   to the merge step is always in `[0, min(|a|,|b|)]`.
2. `scs`: `olen = overlap(ssperm[i], ssperm[i+1], min_length=1)`; `sup += ssperm[i+1][olen:]`
   ("add non-overlapping portion of B to superstring"). This is exactly `MergeContigs` default mode.
3. Out-of-range lengths are never produced by `overlap`; applied blindly, Python slicing would give
   `"BAA" + "AAB"[-2:] = "BAAAB"` and `"AC" + "GTAA"[5:] = "AC"` — not the old Seqeron fallback
   (`c1 + c2`), which therefore had no source (F-MERGE-1).

### EMBOSS `merger` — emboss/merger.c (reference implementation for overlap conflicts)

**URL:** https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/merger.c (552 lines, read in full)
**Accessed:** 2026-11-25 (neutral User-Agent). Authors G. Williams / A. Bleasby; EMBOSS: Rice, Longden &
Bleasby 2000, Trends Genet 16:276–277.

**Key Extracted Points:**

1. "Merge two overlapping nucleic acid sequences": global alignment (`embAlignPathCalc`/`embAlignWalkNWMatrix`),
   then `merger_Merge` walks the aligned columns. Both inputs are lower-cased first (`ajSeqFmtLower`).
2. Column rule: gap → base from the better-quality sequence (if not a gap); `p` is `n/N` → `q`; `q` is
   `n/N` → `p`; `p != q` → base of the sequence for which `merger_bestquality` holds; else `p`.
3. `merger_bestquality`: `merger_quality` with window 5, then 10, then 20 while equal; `qa >= qb` → first
   sequence ("If both sequences have the same quality, pick the first").
4. `merger_quality(seq, pos, window)`: sum over `i = pos .. pos+window-1` and `j = pos-1` down to
   `> pos-window` (and `>= 0`, `< strlen`): `aAcCgGtTuU` = 2 points, `mrwsykvhdbMRWSYKVHDB` = 1, else 0
   (N's and off-end = 0); result `(float)(value / (window*2.0 + 1.0))`. "THIS HEAVILY DISCRIMINATES
   AGAINST THE IFFY BITS AT THE END OF SEQUENCE READS".
5. Left/right flanks: the longer unaligned flank is copied; for a suffix/prefix overlap that is
   `contig1[..start]` and `contig2[l..]`.

**Oracle.** EMBOSS could not be installed in the sandbox (package install not permitted), so
`merger_Merge`/`merger_bestquality`/`merger_quality` were ported line by line to Python (ungapped
alignment `m = a[start1:]`, `n = b[:l]`, float32 quality) and run on the test inputs:

| contig1 | contig2 | l | EMBOSS-port result | deciding window (contig1 vs contig2) |
|---|---|---|---|---|
| AAAACCCCGGGG | CCCCGTGGTTTT | 8 | AAAACCCCGTGGTTTT | w5: 14/11 vs 18/11 → contig2 |
| ACGTACGTACGT | ACGAACGTACGTTT | 8 | ACGTACGTACGTACGTTT | w5: 18/11 vs 16/11 → contig1 |
| AGTACAC | ACGCGG | 4 | AGTACACGG | w5 tie; w10 0.667 vs 0.571 → contig1 |
| CTTTCA | TTTTCAA | 6 | TTTTCAA | w5 tie; w10 0.571 vs 0.667 → contig2 |
| CGGCTAGCAA | GGGCTAGCAAG | 10 | GGGCTAGCAAG | w5, w10 tie; w20 0.488 vs 0.537 → contig2 |
| GGGGGATAGCC | GGGGATAGCG | 10 | GGGGGATAGCC | w5, w10 tie; w20 0.537 vs 0.488 → contig1 |
| AGTATT | TTTTGA | 4 | AGTATTGA | tie at 5/10/20 → contig1 |
| ACGTACGTACGTNNAC | ACGTTTTTTTTT | 4 | ACGTACGTACGTACGTTTTTTTTT | N→contig2 base; then w5 8/11 vs 14/11 → contig2 |
| ACGTANGTAC | ANGTACCCGG | 6 | ACGTANGTACCCGG | N in both → N |
| ACGTACGG | ACNGTT | 4 | ACGTACGGTT | N in contig2 → contig1 base |
| BAA | AAB | 2 | BAAB | exact overlap = textbook collapse |
| acgtAC | GTAcgt | 4 | acgtacgt (EMBOSS lower-cases) | no conflict (case-insensitive) |

The tie-break inputs were found by a seeded random search (seed 7, 200 000 trials) over the port.

## Documented Corner Cases and Failure Modes

### From Langmead SCS notes

1. **No overlap → plain concatenation:** when no length-≥-l suffix/prefix match exists,
   `suffixPrefixMatch` returns 0; the strings are then simply concatenated ("just concatenate them").
   A merge with overlap length 0 is therefore equal to `X + Y`.
2. **Overlap bounded by the shorter string:** the overlap is a length-l suffix of X *and* a length-l
   prefix of Y, so l cannot exceed `min(|X|, |Y|)` (you cannot take a suffix/prefix longer than the
   string itself). `suffixPrefixMatch` guards `if len(x) < k or len(y) < k: return 0`.

### From Langmead OLC notes

1. **Longest match reported:** when several overlaps exist, only the longest suffix/prefix match is
   used; merging with that length removes exactly one copy of the overlapping region.

---

## Test Datasets

### Dataset: Greedy-SCS worked trace (l = 1)

**Source:** Langmead, "Assembly & shortest common superstring", greedy trace (point 4 above).

| Parameter | Value |
|-----------|-------|
| Inputs | `BAA`, `AAB` |
| Overlap length | 2 (suffix `AA` of `BAA` = prefix `AA` of `AAB`) |
| Merged result | `BAAB` |
| Merged length | 4 (= 3 + 3 − 2) |

### Dataset: Greedy-SCS chain {AAA, AAB, ABB, BBB, BBA}

**Source:** Langmead SCS notes (point 5 above).

| Parameter | Value |
|-----------|-------|
| `AAA` + `AAB`, overlap 2 | `AAAB` |
| `AAB` + `ABB`, overlap 2 | `AABB` |
| Full chain SCS | `AAABBBA` (length 7) |

### Dataset: Plain concatenation (no overlap)

**Source:** Langmead SCS notes ("just concatenate them").

| Parameter | Value |
|-----------|-------|
| Inputs | `BAA`, `AAB` |
| Overlap length | 0 |
| Merged result | `BAAAAB` (= `BAA` + `AAB`, length 6) |

---

## Assumptions

1. **ASSUMPTION: caller-supplied overlap length is trusted (not re-verified).**
   `MergeContigs(contig1, contig2, overlapLength)` is a low-level collapse primitive: it removes
   `overlapLength` characters from the front of `contig2` and appends the remainder to `contig1`. The
   sources define the *overlap* as a verified suffix/prefix match; in this library the verification is
   the separate responsibility of `FindOverlap` / `FindAllOverlaps` (which compute the length passed
   here). Merging with an unverified length does not change the merge arithmetic, only whether the
   collapsed region truly matched; this is a documented contract boundary, not a numeric/scoring value.
   It is therefore an API-contract assumption, not a correctness-affecting algorithm parameter.

2. ~~ASSUMPTION: out-of-range overlap (≤ 0 or > min length) → plain concatenation.~~ **Withdrawn
   2026-11-25 (F-MERGE-1).** No source defines a merge at a negative or over-long overlap (Langmead's
   `overlap()` only yields `0..min(|x|,|y|)`, and Python slicing of such values gives different strings),
   so these lengths now throw `ArgumentOutOfRangeException`, like the sibling OLC threshold checks (B10 F4).

3. **Default-mode mismatch handling.** The default keeps contig1's characters in the overlap (Langmead
   `sup += y[olen:]`); mismatch/`N` resolution is the opt-in EMBOSS `merger` rule (source above).

---

## Recommendations for Test Coverage

1. **MUST Test:** merge two strings with a published overlap length and assert the exact collapsed
   string and its length (`BAA`+`AAB`, overlap 2 → `BAAB`). — Evidence: Langmead SCS greedy trace.
2. **MUST Test:** chain three published overlaps to the full SCS (`AAABBBA`). — Evidence: Langmead SCS.
3. **MUST Test:** overlap length 0 → plain concatenation `X + Y`. — Evidence: Langmead SCS ("concatenate").
4. **MUST Test:** overlap = `min(|X|,|Y|)` boundary (full containment of the shorter prefix). — Evidence:
   overlap bounded by the shorter string (SCS `suffixPrefixMatch` guard).
5. **SHOULD Test:** overlapLength > min length → falls back to concatenation. — Rationale: documented bound.
6. **SHOULD Test:** negative overlap → concatenation. — Rationale: non-positive overlap = no overlap.
7. **SHOULD Test:** null contig1 / null contig2 → `ArgumentNullException`. — Rationale: sibling methods
   validate null inputs explicitly (`AssembleOLC`, `BuildDeBruijnGraph`).
8. **COULD Test:** empty + non-empty contig (overlap 0) returns the non-empty contig. — Rationale: identity.
9. **COULD Test:** result length invariant `|merge| = |c1| + |c2| − usedOverlap`. — Rationale: property check.

---

## References

1. Langmead, B. (n.d.). *Assembly & shortest common superstring* (JHU computational genomics lecture notes).
   https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf (accessed 2026-06-13).
2. Langmead, B. (n.d.). *Overlap Layout Consensus assembly* (JHU computational genomics lecture notes).
   https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf (accessed 2026-06-13).
3. MIT 7.91J / Burge C., Fraenkel E. (2014). *Foundations of Computational and Systems Biology,
   Lecture 6: Genome Assembly*. MIT OpenCourseWare.
   https://ocw.mit.edu/courses/7-91j-foundations-of-computational-and-systems-biology-spring-2014/e885f0eb376ea6c2045eb9d8847f106f_MIT7_91JS14_Lecture6.pdf (accessed 2026-06-13).
4. Compeau, P.E.C., Pevzner, P.A., Tesler, G. (2011). *How to apply de Bruijn graphs to genome
   assembly*. Nature Biotechnology 29:987–991. https://doi.org/10.1038/nbt.2023
   (background on overlap/de Bruijn assembly; cited in source 1).

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-11-25**: B10 review — Langmead `CG_SCS.ipynb` code and EMBOSS `merger.c` opened; out-of-range
  overlap now rejected (F-MERGE-1); opt-in EMBOSS conflict rule (F-MERGE-2) with port-oracle table.
