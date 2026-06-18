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

2. **ASSUMPTION: out-of-range overlap (≤ 0 or > min length) → plain concatenation.**
   The sources define overlap length 0 as "no overlap → concatenate" and bound a valid overlap by
   `min(|X|,|Y|)`. Treating a non-positive or oversized requested overlap as "no usable overlap, so
   concatenate `X + Y`" follows directly from those two source facts; it is the only behavior
   consistent with the suffix/prefix definition.

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
