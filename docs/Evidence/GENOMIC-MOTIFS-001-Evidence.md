# Evidence Artifact: GENOMIC-MOTIFS-001

**Test Unit ID:** GENOMIC-MOTIFS-001
**Algorithm:** Known Motif Search (multi-pattern exact substring matching)
**Date Collected:** 2026-06-13

---

## Online Sources

### Gusfield-based exact-matching course notes (Tufts COMP 150GEN)

**URL:** https://www.cs.tufts.edu/comp/150GEN/classpages/exact.html
**Accessed:** 2026-06-13
**Authority rank:** 1 (textbook-derived: Gusfield D., *Algorithms on Strings, Trees and Sequences*; course primary reading)
**Retrieved by:** WebSearch "exact string matching report all occurrences including overlapping positions definition" → WebFetch of the page.

**Key Extracted Points:**

1. **Exact matching problem (verbatim):** "find all occurrences of a pattern string P of length n in a text string T of length m". The result is the set of **all** start positions where P occurs in T.
2. **Overlapping occurrences (verbatim):** "If P = aaa and T = aaaaa, there are three (overlapping) occurrences". The problem counts every start position where P aligns, even when occurrences overlap; it is NOT a non-overlapping (greedy-skip) count.

### Biopython `Bio.Seq` reference implementation (master branch)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established bioinformatics library)
**Retrieved by:** WebSearch "Biopython Seq.search method returns all overlapping occurrences start positions subsequence" → WebFetch of the raw source file.

**Key Extracted Points:**

1. **Multi-pattern motif search (`Seq.search`):** "Search the substrings subs in self and yield the index and substring found." It accepts a list of multiple search terms (e.g. `dna.search(["CC", Seq("ATTG"), "ATTG", Seq("CCC")])`) and yields `(index, substring)` pairs — i.e., for a *set* of query motifs it reports where each motif is found in the subject. This is the multi-motif "known motif search" semantics of this unit.
2. **Overlapping occurrences counted (`count_overlap`):** docstring "Returns an integer, the number of occurrences of substring argument sub in the (sub)sequence". Example `Seq("AAAA").count_overlap("AA")` returns **3**, contrasted with non-overlapping `"AAAA".count("AA")` returning **2**. Confirms the biologically/algorithmically correct count for occurrence reporting includes overlaps.
3. **Case normalization:** the module upper-cases sequence data before processing (`sequence = sequence.upper()` in `_translate_str`); DNA work is done on an upper-cased alphabet.

### Wikipedia "Restriction site" (citing primary palindrome fact) — biological motif example

**URL:** https://en.wikipedia.org/wiki/Restriction_site
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia, used for a widely-documented sequence fact also held in REBASE/NEB rank-5 DBs)
**Retrieved by:** WebSearch "REBASE EcoRI recognition sequence GAATTC restriction enzyme" → WebFetch of the article.

**Key Extracted Points:**

1. **EcoRI recognition site (verbatim):** "the common restriction enzyme EcoRI recognizes the palindromic sequence GAATTC and cuts between the G and the A on both the top and bottom strands." Provides a real biological motif (`GAATTC`, length 6) for the test dataset.

---

## Documented Corner Cases and Failure Modes

### From Tufts COMP 150GEN / Gusfield

1. **Overlapping occurrences must all be reported:** `P=aaa`, `T=aaaaa` → 3 occurrences (positions 0,1,2 in 0-based indexing). A search that skips past a match (non-overlapping) under-reports and is a defect for this problem.
2. **Pattern longer than text / absent pattern:** "all occurrences" of a pattern that does not occur is the empty set.

### From Biopython `Bio.Seq`

1. **Multiple query motifs:** a set of motifs is searched; each motif independently maps to its set of positions in the subject (`search` yields per-motif hits).
2. **Empty query:** Biopython upper-cases and treats inputs as strings; an empty pattern is a degenerate query (handled by returning no meaningful occurrence — see Assumptions for repository policy).

---

## Test Datasets

### Dataset: Overlapping homopolymer (Gusfield worked example)

**Source:** Tufts COMP 150GEN exact.html (Gusfield-derived): "If P = aaa and T = aaaaa, there are three (overlapping) occurrences".

| Parameter | Value |
|-----------|-------|
| Text T | `AAAAA` (length 5) |
| Motif P | `AAA` (length 3) |
| Expected occurrences (0-based starts) | `0, 1, 2` (count = 3, overlapping) |

### Dataset: EcoRI motif in a synthetic sequence

**Source:** EcoRI recognition site `GAATTC` (Wikipedia "Restriction site"). Positions derived by the exact-matching definition (T[i..i+5] == "GAATTC"), 0-based.

| Parameter | Value |
|-----------|-------|
| Text T | `GAATTCAAAGAATTC` (length 15) |
| Motif P | `GAATTC` (length 6) |
| Expected occurrences (0-based starts) | `0, 9` (count = 2, non-overlapping by construction) |

### Dataset: Multi-motif set (Biopython `search` semantics)

**Source:** Biopython `Seq.search` (multiple query motifs yield per-motif hits). Positions derived by the exact-matching definition, 0-based.

| Parameter | Value |
|-----------|-------|
| Text T | `ACGTACGTAA` (length 10) |
| Motif set | `{ "ACGT", "AA", "TTT" }` |
| `ACGT` occurrences | `0, 4` |
| `AA` occurrences | `8` |
| `TTT` occurrences | (none — absent, not in result) |

---

## Assumptions

1. **ASSUMPTION: Empty motif handling** — Neither Gusfield nor Biopython defines a single canonical return for searching an *empty* pattern in the context of a known-motif set (suffix-tree `FindAllOccurrences("")` would return every position, which is not a meaningful "motif"). The repository convention (`GenomicAnalyzer.FindMotif` returns `Array.Empty<int>()` for null/empty motif) is reused: empty/whitespace motifs contribute no entry to the result. This is an API-shape policy, not an algorithm-correctness parameter (no authoritative source defines an empty-motif occurrence set as biologically meaningful).

2. **ASSUMPTION: Result key normalization** — Motifs are upper-cased before search (consistent with Biopython upper-casing DNA, and with `DnaSequence` upper-casing on construction); the result is keyed by the upper-cased motif. Source-supported by Biopython case handling (Point 3 above).

---

## Recommendations for Test Coverage

1. **MUST Test:** Overlapping occurrences fully reported (`AAA` in `AAAAA` → {0,1,2}). — Evidence: Tufts/Gusfield exact.html worked example.
2. **MUST Test:** Multi-motif set returns per-motif position lists; absent motifs are omitted. — Evidence: Biopython `Seq.search`.
3. **MUST Test:** Real biological motif `GAATTC` (EcoRI) located at all start positions. — Evidence: Wikipedia "Restriction site".
4. **MUST Test:** Returned positions are sorted ascending and 0-based (deterministic contract). — Evidence: exact-matching defines positions as a set; deterministic ordering required for a stable public contract.
5. **SHOULD Test:** Empty motif set → empty result; motif absent → omitted. — Rationale: documented degenerate inputs.
6. **SHOULD Test:** Lower-case / mixed-case motif normalized and matched; key is upper-cased. — Rationale: Biopython case handling.
7. **COULD Test:** Single-character motif counts all positions including consecutive (overlapping) ones. — Rationale: boundary of the overlap rule.

---

## References

1. Gusfield, D. (1997). *Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology*. Cambridge University Press. ISBN 0-521-58519-8. Exact-matching problem definition and overlapping-occurrence example as taught in Tufts COMP 150GEN: https://www.cs.tufts.edu/comp/150GEN/classpages/exact.html
2. Cock, P.J.A. et al. (Biopython). `Bio.Seq` module, `search` and `count_overlap` methods (master branch, accessed 2026-06-13): https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
3. Wikipedia contributors. "Restriction site" (accessed 2026-06-13). EcoRI recognizes GAATTC: https://en.wikipedia.org/wiki/Restriction_site

---

## Change History

- **2026-06-13**: Initial documentation.
