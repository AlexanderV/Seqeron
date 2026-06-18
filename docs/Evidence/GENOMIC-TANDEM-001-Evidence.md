# Evidence Artifact: GENOMIC-TANDEM-001

**Test Unit ID:** GENOMIC-TANDEM-001
**Algorithm:** Tandem Repeat Detection (`GenomicAnalyzer.FindTandemRepeats`)
**Date Collected:** 2026-06-14

---

## Online Sources

### Benson, G. (1999) — Tandem Repeats Finder, *Nucleic Acids Research* 27(2):573–580

**URL:** https://academic.oup.com/nar/article/27/2/573/1061099
**Accessed:** 2026-06-14 (fetched the article page in this session)
**Authority rank:** 1 (peer-reviewed primary paper; the foundational tandem-repeat detection paper)

**Key Extracted Points:**

1. **Formal definition:** "A tandem repeat in DNA is two or more contiguous … copies of a pattern of nucleotides." The repeated unit is the **pattern**; individual iterations are **copies**; copies are contiguous (directly adjacent).
2. **Period / copy number:** the **period** is the distance between matching positions in adjacent copies (the unit length for an exact repeat); **copy number** is how many copies are present. Reported per repeat: "period size; number of copies aligned with the consensus pattern; size of the consensus pattern."
3. **Minimum copies:** a tandem repeat requires **two or more** contiguous copies (k ≥ 2).
4. **DOI:** https://doi.org/10.1093/nar/27.2.573

### Wikipedia — "Tandem repeat" (using its cited primaries)

**URL:** https://en.wikipedia.org/wiki/Tandem_repeat
**Accessed:** 2026-06-14 (fetched the page, and re-confirmed via web search in this session)
**Authority rank:** 4 (Wikipedia; the definition and example below are quoted verbatim from the article)

**Key Extracted Points:**

1. **Definition (verbatim):** "tandem repeats occur in DNA when a pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other."
2. **Worked example (verbatim):** "ATTCG ATTCG ATTCG, in which the sequence ATTCG is repeated three times." → unit `ATTCG` (period 5), copy number 3, total length 15, contiguous.
3. **Classification by unit length:** dinucleotide (2 nt), trinucleotide (3 nt); microsatellites/STRs are shorter units; minisatellites = 10–60 nt repeated; macrosatellites ≈ 1,000 nt.
4. **Genome fraction / disease:** "Tandem repeats constitute about 8% of the human genome"; "implicated in more than 50 lethal human diseases."
5. **Detection:** "Tandem repeats in strings can be efficiently detected using suffix trees or suffix arrays."

---

## Documented Corner Cases and Failure Modes

### From Benson (1999) / Wikipedia (Tandem repeat)

1. **Exact vs approximate copies:** Benson's TRF detects *approximate* tandem copies; the repository implementation detects only *exact* contiguous copies (a documented simplification — see algorithm doc §5.3). Over exact repeats the definitions coincide.
2. **Minimum copy number:** fewer than two contiguous copies is not a tandem repeat (k ≥ 2). The implementation's `minRepetitions` default is 2, matching the floor.
3. **Multiple unit-length interpretations:** the same region (e.g. `AAAA`) satisfies the definition for several period sizes (period 1 ×4, period 2 ×2); the implementation reports each unit-length interpretation that meets the threshold (documented in algorithm doc §6).

---

## Test Datasets

### Dataset: Wikipedia worked example (Tandem repeat)

**Source:** Wikipedia "Tandem repeat" (accessed 2026-06-14), example quoted verbatim.

| Parameter | Value |
|-----------|-------|
| Sequence | `ATTCGATTCGATTCG` |
| Unit (pattern) | `ATTCG` |
| Period | 5 |
| Copy number (repetitions) | 3 |
| Total length | 15 |
| Start (0-based) | 0 |

### Dataset: Canonical fixture trinucleotide (REP-TANDEM-001 M1)

**Source:** definition (Benson 1999; Wikipedia trinucleotide class).

| Parameter | Value |
|-----------|-------|
| Sequence | `ATGATGATG` |
| Unit | `ATG` |
| Repetitions | 3 |
| Start (0-based) | 0 |

---

## Assumptions

1. **ASSUMPTION: exact-copy restriction** — The repository detector reports only exact contiguous copies, not the approximate copies of Benson's TRF. This is a documented simplification (algorithm doc §5.3); over exact repeats the output matches the formal definition, so no correctness-affecting parameter is invented.

---

## Recommendations for Test Coverage

1. **MUST Test:** exact contiguous trinucleotide unit returns unit/position/repetition count (e.g. `ATGATGATG` → `ATG`, pos 0, 3 copies) — Evidence: Wikipedia/Benson definition.
2. **MUST Test:** Wikipedia worked example `ATTCGATTCGATTCG` → `ATTCG`, pos 0, 3 copies, total length 15 — Evidence: Wikipedia worked example.
3. **SHOULD Test:** `minRepetitions` floor (k ≥ 2) and `minUnitLength` threshold honored — Rationale: definition requires two or more copies.
4. **COULD Test:** no-repeat and empty inputs return empty — Rationale: boundary behavior.

> All of the above are already implemented and verified under **REP-TANDEM-001** (canonical fixture `GenomicAnalyzer_TandemRepeat_Tests.cs`, 27 tests). GENOMIC-TANDEM-001 is a duplicate Registry entry for the identical method and is resolved by consolidation — see the TestSpec §7. No new or duplicate tests are created.

---

## References

1. Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research* 27(2):573–580. https://doi.org/10.1093/nar/27.2.573
2. Wikipedia. Tandem repeat. https://en.wikipedia.org/wiki/Tandem_repeat (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation. Sources re-retrieved this session to verify the existing `GenomicAnalyzer.FindTandemRepeats` implementation; unit resolved as a duplicate of REP-TANDEM-001 by consolidation.
