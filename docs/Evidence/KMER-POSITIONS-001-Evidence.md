# Evidence Artifact: KMER-POSITIONS-001

**Test Unit ID:** KMER-POSITIONS-001
**Algorithm:** K-mer Positions (find all start positions where a k-mer occurs in a sequence)
**Date Collected:** 2026-06-14

---

## Online Sources

### Rosalind — "Find All Occurrences of a Pattern in a String" (Problem BA1D)

**URL:** https://rosalind.info/problems/ba1d/
**Accessed:** 2026-06-14 (retrieved with WebFetch on the URL above)
**Authority rank:** 4 (educational problem set built on Compeau & Pevzner, *Bioinformatics Algorithms*; used here for its worked example with exact expected output)

**Key Extracted Points:**

1. **Problem:** "Find all occurrences of a pattern in a string." Input is two strings, *Pattern* and *Genome*.
2. **Output:** "All starting positions in *Genome* where *Pattern* appears as a substring. Use 0-based indexing." — positions are explicitly **0-based**.
3. **Overlapping:** Different occurrences of a substring can overlap and are all reported (overlapping occurrences are counted).
4. **Worked example (verbatim):** Pattern = `ATAT`, Genome = `GATATATGCATATACTT` → Output: `1 3 9`.

### Wikipedia — "k-mer" (using cited primaries; used for the formal definition)

**URL:** https://en.wikipedia.org/wiki/K-mer
**Accessed:** 2026-06-14 (retrieved with WebFetch on the URL above)
**Authority rank:** 4

**Key Extracted Points:**

1. **Definition:** A k-mer is a "substring of length k contained within a biological sequence."
2. **Count formula:** For a sequence of length L, the number of k-mers (overlapping, every start position) is **L − k + 1**.
3. **Overlapping example:** The sequence `AGAT` contains three 2-mers — `AG`, `GA`, `AT` — i.e. consecutive k-mers overlap by k−1 positions and start at every valid index.

### Compeau & Pevzner — *Bioinformatics Algorithms: An Active Learning Approach* (Pattern Matching Problem)

**URL:** https://gerdos.web.elte.hu/edu/bioinformatics_algorithms/week1.pdf (ELTE lecture notes built on the textbook)
**Accessed:** 2026-06-14 (retrieved with WebFetch; PDF is image/compressed and did not yield clean extractable text, so the Pattern Matching definition was taken from the Rosalind BA1D statement above, which restates the same textbook problem). Search query used: "Compeau Pevzner Bioinformatics Algorithms PatternMatching find all occurrences pattern 0-based positions" (WebSearch).
**Authority rank:** 1 (textbook)

**Key Extracted Points:**

1. **Pattern Matching Problem (search-result summary):** Find all occurrences of a pattern in a string; output all starting positions where the pattern appears as a substring. The textbook prose uses 1-based indexing in narrative form, but the associated Rosalind BA1D exercise (the canonical machine-checked form) specifies **0-based** output — recorded as the binding convention for this unit (see Assumptions).

### Rosalind — "Finding a Motif in DNA" (Problem SUBS) — added 2026-10-01 (B06 review)

**Opened:** rosalind.info is blocked for curl; statement and samples were read from GitHub mirrors via raw.githubusercontent.com: `breezedu/rosalind/master/FindingaMotifinDNA.java` (header: Sample Dataset `GATATATGCATATACTT` / `ATAT`, Sample Output `2 4 10`) and `mtarbit/Rosalind-Problems/master/e009-subs.py` (older statement: "The position of a symbol in a string is the total number of symbols found to its left, including itself"; sample `ACGTACGTACGTACGT` / `GTA` → `3 7 11`). Confirmed by a WebSearch snippet of rosalind.info/problems/subs/.

**Key point:** SUBS reports the same overlapping occurrence set as BA1D but **1-based**; 0-based = SUBS − 1 (`1 3 9`, `2 6 10`).

### Biopython / Python `re` reference implementations — added 2026-10-01

Biopython 1.88 `Bio.SeqUtils.nt_search(seq, subseq)` (source inspected: forward strand only, IUPAC regex, re-searches from `pos + 1`, so overlapping hits are returned, 0-based, upper-case input required) and `re.finditer("(?=P)", s)` gave identical results: GATATATGCATATACTT/ATAT → [1,3,9]; AAAA/AA → [0,1,2]; ATATATAT/ATAT → [0,2,4]; AAAAAAAAAA/AAA → [0..7]; ACGTACGTACGT/CGTA → [1,5]; GATATATGCATATACTT/GCAT → [7]; ACGTACGTACGTACGT/GTA → [2,6,10]. `Seq.count_overlap` equals the list length in every case. 5000 random cases (L ≤ 60, k ≤ 8, alphabets AC/ACGT/ACGTN/mixed-case/A) vs C#: 0 mismatches.

### Knuth–Morris–Pratt — added 2026-10-01

Knuth, Morris & Pratt (1977) *SIAM J. Comput.* 6:323–350 and CLRS §32.4 (KMP-MATCHER / COMPUTE-PREFIX-FUNCTION; recalled from the textbook, not re-opened): O(m) prefix function + O(n) matching; after a full match continue from π[m] to find overlapping occurrences. Failure-array reference code opened: TheAlgorithms/Python `strings/knuth_morris_pratt.py` (raw.githubusercontent.com; self-test text `ABABZABABYABABX` / `ABABX`). Used to replace the Θ((L−k+1)·k) window scan (measured 5.35 s for `A`^10⁶ vs `A`^99 999`C`; KMP 9–26 ms).

---

## Documented Corner Cases and Failure Modes

### From Rosalind BA1D / Wikipedia k-mer

1. **Overlapping occurrences:** when the pattern self-overlaps (e.g. `ATAT` in `...ATATAT...`, or `AA` in `AAAA`), every overlapping start position is reported, not just non-overlapping ones.
2. **Pattern longer than text:** if |kmer| > |sequence| there are 0 valid start positions (L − k + 1 ≤ 0) → empty result.
3. **No occurrence:** pattern absent → empty result.
4. **Pattern equals whole sequence:** exactly one occurrence at position 0.

---

## Test Datasets

### Dataset: Rosalind BA1D sample

**Source:** Rosalind BA1D, https://rosalind.info/problems/ba1d/

| Parameter | Value |
|-----------|-------|
| Pattern (kmer) | `ATAT` |
| Genome (sequence) | `GATATATGCATATACTT` |
| Expected positions (0-based) | `1 3 9` |

### Dataset: Wikipedia AGAT 2-mers

**Source:** Wikipedia "k-mer", https://en.wikipedia.org/wiki/K-mer

| Parameter | Value |
|-----------|-------|
| Sequence | `AGAT` |
| k | 2 |
| 2-mers in order | `AG`@0, `GA`@1, `AT`@2 |

### Dataset: Self-overlap derivation

**Source:** Direct application of the BA1D overlapping rule (overlapping occurrences counted, 0-based).

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAA` |
| kmer | `AA` |
| Expected positions | `0 1 2` (L − k + 1 = 3) |

---

## Assumptions

1. **ASSUMPTION: Indexing convention is 0-based.** The textbook prose narrates positions 1-based, but the canonical machine-checked exercise (Rosalind BA1D) and the repository's existing k-mer methods use 0-based indexing. 0-based is adopted; consistent with C# string indexing and with the repository SuffixTree (`FindAllOccurrences("ana")` → `[1,3]` for `banana`).
2. **ASSUMPTION: Case-insensitive matching.** No authoritative source mandates case-folding for k-mer position search. The repository convention (sibling `KmerAnalyzer` methods upper-case input) is adopted so that, e.g., `atat` matches `ATAT`. This is a repository interoperability choice, recorded as an assumption; it does not affect any all-uppercase evidence example.
3. **ASSUMPTION: Null/empty input returns empty.** No source defines behavior for null/empty `sequence` or `kmer`. The repository convention (sibling methods) of returning an empty result (no exception) is adopted.

---

## Recommendations for Test Coverage

1. **MUST Test:** Rosalind BA1D sample — `ATAT` in `GATATATGCATATACTT` → exactly `[1,3,9]` in order. — Evidence: Rosalind BA1D.
2. **MUST Test:** Overlapping self-occurrence — `AA` in `AAAA` → `[0,1,2]`. — Evidence: BA1D overlapping rule.
3. **MUST Test:** AGAT 2-mers / "ana" in "banana" → ascending start positions. — Evidence: Wikipedia k-mer; SuffixTree doc example.
4. **MUST Test:** Pattern absent → empty. — Evidence: BA1D (only matching starts reported).
5. **SHOULD Test:** Pattern longer than text → empty (L − k + 1 ≤ 0). — Rationale: count formula boundary.
6. **SHOULD Test:** Pattern equals whole sequence → `[0]`. — Rationale: single-occurrence boundary.
7. **SHOULD Test:** Case-insensitive match (`atat` ≡ `ATAT`). — Rationale: documents the case-folding assumption.
8. **COULD Test:** null/empty sequence or kmer → empty. — Rationale: documents null/empty handling assumption.

---

## References

1. Rosalind. (accessed 2026-06-14). Find All Occurrences of a Pattern in a String (Problem BA1D). https://rosalind.info/problems/ba1d/
2. Wikipedia contributors. (accessed 2026-06-14). k-mer. https://en.wikipedia.org/wiki/K-mer
3. Compeau, P., Pevzner, P. (2015). Bioinformatics Algorithms: An Active Learning Approach (Pattern Matching Problem). Active Learning Publishers. Lecture-notes mirror: https://gerdos.web.elte.hu/edu/bioinformatics_algorithms/week1.pdf
4. Rosalind. Finding a Motif in DNA (Problem SUBS). https://rosalind.info/problems/subs/ (via GitHub mirrors, see above).
5. Biopython 1.88, `Bio.SeqUtils.nt_search`.
6. Knuth, D.E., Morris, J.H., Pratt, V.R. (1977). Fast pattern matching in strings. SIAM J. Comput. 6(2):323–350.

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-10-01**: B06 review — Rosalind SUBS (1-based) samples, Biopython `nt_search` / regex cross-check, KMP source (implementation switched to KMP, O(L + k)).
