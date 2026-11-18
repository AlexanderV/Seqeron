# Evidence Artifact: GENOMIC-TANDEM-001

**Test Unit ID:** GENOMIC-TANDEM-001
**Algorithm:** Tandem Repeat Detection (`GenomicAnalyzer.FindTandemRepeats`)
**Date Collected:** 2026-06-14 (extended 2026-11-18, B09 F35–F38)

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
**Accessed:** 2026-06-14 (fetched the page, and re-confirmed via web search in this session) · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4 (Wikipedia; the definition and example below are quoted verbatim from the article)

**Key Extracted Points:**

1. **Definition (verbatim):** "tandem repeats occur in DNA when a pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other."
2. **Worked example (verbatim):** "ATTCG ATTCG ATTCG, in which the sequence ATTCG is repeated three times." → unit `ATTCG` (period 5), copy number 3, total length 15, contiguous.
3. **Classification by unit length:** dinucleotide (2 nt), trinucleotide (3 nt); microsatellites/STRs are shorter units; minisatellites = 10–60 nt repeated; macrosatellites ≈ 1,000 nt.
4. **Genome fraction / disease:** "Tandem repeats constitute about 8% of the human genome"; "implicated in more than 50 lethal human diseases."
5. **Detection:** "Tandem repeats in strings can be efficiently detected using suffix trees or suffix arrays."

### Bannai, I, Inenaga, Nakashima, Takeda, Tsuruta (2017) — The "Runs" Theorem, SIAM J. Comput. 46(5):1501–1514

**URL:** https://arxiv.org/abs/1406.0263 (abstract page and full text PDF opened 2026-11-18, neutral User-Agent)
**Authority rank:** 1

**Key Extracted Points (verbatim from the opened text):**

1. **Definition 1 (Runs):** "A triple r = (i, j, p) is a run of string w, if the smallest period p of w[i..j] satisfies |w[i..j]| ≥ 2p, and the periodicity cannot be extended to the left or right, i.e., i = 1 or w[i−1] ≠ w[i+p−1], and, j = n or w[j+1] ≠ w[j−p+1]. The rational number (j−i+1)/p is called the exponent of r." — the maximal-repetition definition (Kolpakov & Kucherov, FOCS 1999) used by this unit.
2. "Note that a Lyndon word w cannot have any period p < |w|" — a Lyndon root is primitive; hence a candidate (i, p) from a longest Lyndon word spanning ≥ 2p is a run with smallest period p.
3. **§3:** orders ≺0, ≺1 with a ≺0 b ⇔ b ≺1 a; ŵ = w$ with "$ ≺0 a (and thus a ≺1 $)"; l_ℓ(i) = longest Lyndon word of ŵ starting at i. **Lemma 7:** for a run r with ŵ[j+1] ≺ℓ ŵ[j+1−p], any L-root λ = [iλ..jλ] with iλ ≠ i equals l_ℓ(iλ). **Theorem 9:** ρ(n) < n.
4. **§4.1–4.2:** l_ℓ(i) by a stack pass using suffix comparisons (inverse suffix array); for a candidate [i..j] the run is (i′, j′, p) with p = j−i+1, j′ from the longest common prefix of w[i..n] and w[j+1..n], i′ from the longest common suffix of w[1..i−1] and w[1..j] (LCE queries on w and its reverse), kept iff |w[i′..j′]| ≥ 2p.

### Thiel et al. (2003) — MISA `misa.pl`

**URL:** https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl (opened 2026-11-18)
**Authority rank:** 2

1. Line 110: `#reject false type motifs [e.g. (TT)6 or (ACAC)5]` — non-primitive motifs are rejected; a motif is reported at its primitive length (the run is found at the shorter unit size).

### Kolpakov & Kucherov (1999) — Finding maximal repetitions in a word in linear time (FOCS, 596–604)

**URL:** https://dl.acm.org/doi/10.5555/795665.796470 — only search-engine snippets were available (page not opened): a repetition is a sub-word whose period is at most half its length; maximal repetitions cannot be extended while keeping the period; their number is linear in the word length (conjecture < n, proved by Bannai et al.).

### Benson (1999) — definition re-confirmed 2026-11-18

The publisher pages (academic.oup.com, doi.org) return HTTP 403 for scripted clients; the definition quoted below was re-confirmed from search-engine snippets of the indexed abstract/full text ("two or more contiguous, approximate copies of a pattern of nucleotides").

### Independent oracle (2026-11-18)

`scripts/analysis/generate_tandem_oracle.py` → `tests/.../TestData/Analysis/tandem_oracle.tsv` (1156 rows: `sequence minUnit minReps maxUnit(0 = none) count hits`; hits = `unit:pos:copies;…` ordered by unit length then start). Each row is computed by three independent Python methods that the script asserts equal:

* A: brute-force Kolpakov–Kucherov runs (every interval, its smallest period from the definition, maximality by extension, ⌊len/q⌋ copies).
* B: boolean-array maximal runs of `S[x] = S[x+p]` plus a primitivity test.
* C: regex lookahead `(?=((.{p})\2{k-1,}))` at left-maximal starts plus a primitivity test.

Numbers: both C# engines (per-unit-length scan = `RepeatFinder.FindMicrosatellites`, run enumeration) match all 1156 rows. The **pre-fix** O(n²·m) scan (Python port) differed on **276 / 1156 rows** (487 extra hits, every one a non-primitive unit; none missing, none non-maximal); on 20 000 further random strings it produced 26 747 extra non-primitive hits and 0 missing / non-maximal ones. Worked cases: `ATATATAT`, minUnit 1, minReps 2 → oracle `AT:0:4`; old code `AT:0:4;ATAT:0:2`. `AAAAA` (1,2) → `A:0:5` (old: also `AA:0:2`). `ACGACGACGACGTGTTAGCTTAGCTTAGC` (1,2) → `T:14:2;T:19:2;T:24:2;GT:11:2;ACG:0:4;TTAGC:14:3`; (1,2,max 3) → same without `TTAGC`; (5,2,5) → `TTAGC:14:3`. `GCAGCAGCAGCAGG` (3,2) → `GCA:0:4` (left-end phase), `TCAGCAGCAGCAGT` → `CAG:1:4`.

EMBOSS `etandem` 6.6.0 (installed) was run on 400 random cases (`-minrepeat = -maxrepeat = p`, `-threshold = (k−1)·p`) as an informational check: 297 agreed on (start, primitive consensus); the 103 differences are etandem's own diagonal-scoring model (threshold is a score over matching diagonal positions, not a copy count; different rotation phase; non-contiguous diagonals), so etandem is not an exact-run oracle for this contract. A `pip install pytrf` attempt was declined by the sandbox; the pytrf/MISA comparison of the canonical scan (3099/3132 and 2752/3132 with every difference one of two documented conventions) is B04 F5.

---

## Documented Corner Cases and Failure Modes

### From Benson (1999) / Wikipedia (Tandem repeat)

1. **Exact vs approximate copies:** Benson's TRF detects *approximate* tandem copies; `GenomicAnalyzer.FindTandemRepeats` detects only *exact* contiguous copies. Approximate repeats are implemented in-repo (`RepeatFinder.FindApproximateTandemRepeats`, TRF model; MCP `find_approximate_tandem_repeats`). Over exact repeats the definitions coincide.
2. **Minimum copy number:** fewer than two contiguous copies is not a tandem repeat (k ≥ 2). The implementation's `minRepetitions` default is 2, matching the floor.
3. **Multiple unit-length interpretations:** the same region (e.g. `AAAA`) satisfies the *formal* definition for several period sizes (period 1 ×4, period 2 ×2), but only the primitive unit is a tandem/microsatellite unit (MISA "reject false type motifs"; Kolpakov–Kucherov runs have the smallest period): `AAAA` → `A`×4 only. Runs of different, non-multiple periods that overlap are each reported.
4. **Maximality and phase:** one report per maximal run and unit length at its left end; rotations are not reported again; a trailing partial copy is not counted.

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

1. **Scope: exact copies** — the detector reports only exact contiguous copies; Benson's approximate copies are `RepeatFinder.FindApproximateTandemRepeats`. Over exact repeats the output matches the run definition, so no correctness-affecting parameter is invented.
2. **Convention: left-end phase** — the unit is the motif at the run's left end (MISA/pytrf), rotations are not canonicalised.

---

## Recommendations for Test Coverage

1. **MUST Test:** exact contiguous trinucleotide unit returns unit/position/repetition count (e.g. `ATGATGATG` → `ATG`, pos 0, 3 copies) — Evidence: Wikipedia/Benson definition.
2. **MUST Test:** Wikipedia worked example `ATTCGATTCGATTCG` → `ATTCG`, pos 0, 3 copies, total length 15 — Evidence: Wikipedia worked example.
3. **SHOULD Test:** `minRepetitions` floor (k ≥ 2) and `minUnitLength` threshold honored — Rationale: definition requires two or more copies.
4. **COULD Test:** no-repeat and empty inputs return empty — Rationale: boundary behavior.

> The sourced worked examples live in the canonical fixture `GenomicAnalyzer_TandemRepeat_Tests.cs` (27 tests); the maximal-run / primitive-unit contract, the oracle rows, both engines and the 1 Mbp behaviour are locked by `GenomicAnalyzer_TandemRepeat_Oracle_Tests.cs` (see the TestSpec §4).

---

## References

1. Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research* 27(2):573–580. https://doi.org/10.1093/nar/27.2.573
2. Wikipedia. Tandem repeat. https://en.wikipedia.org/wiki/Tandem_repeat (accessed 2026-06-14).
3. Kolpakov, R., Kucherov, G. (1999). Finding maximal repetitions in a word in linear time. *Proc. 40th IEEE FOCS*, 596–604.
4. Bannai, H., I, T., Inenaga, S., Nakashima, Y., Takeda, M., Tsuruta, K. (2017). The "Runs" Theorem. *SIAM J. Comput.* 46(5):1501–1514. https://arxiv.org/abs/1406.0263
5. Thiel, T. et al. (2003). *Theor. Appl. Genet.* 106:411–422; `misa.pl` mirror https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl

---

## Change History

- **2026-11-18**: B09 review (F35–F38): the "duplicate of REP-TANDEM-001" resolution hid a contract defect (non-primitive units, O(n²·m) scan, no unit-length bound); definition re-sourced (Kolpakov–Kucherov runs, Bannai et al. Definition 1, MISA), three-method Python oracle added, approximate repeats pointed to `RepeatFinder.FindApproximateTandemRepeats`.
- **2026-06-14**: Initial documentation. Sources re-retrieved this session to verify the existing `GenomicAnalyzer.FindTandemRepeats` implementation; unit resolved as a duplicate of REP-TANDEM-001 by consolidation.
