# Evidence Artifact: KMER-UNIQUE-001

**Test Unit ID:** KMER-UNIQUE-001
**Algorithm:** Unique K-mers / K-mers with Minimum Count (k-mer frequency filtering)
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia — K-mer

**URL:** https://en.wikipedia.org/wiki/K-mer
**Accessed:** 2026-06-14
**Authority rank:** 4 (encyclopedia article; used for the foundational definition and worked example, which are standard and uncontested)

**Retrieval:** WebSearch query `k-mer definition substring length k overlapping bioinformatics number of k-mers n-k+1`, then WebFetch of the article URL above.

**Key Extracted Points:**

1. **Definition:** "In bioinformatics, k-mers are substrings of length k contained within a biological sequence." (verbatim from fetched text)
2. **Total-count formula:** A sequence of length L contains **L − k + 1** total k-mers (overlapping, step size +1). The number of *possible* k-mers over an alphabet of n monomers is **n^k** (n = 4 for DNA).
3. **Worked example (AGAT):** monomers (k=1): A, G, A, T → 4; 2-mers: AG, GA, AT → 3; 3-mers: AGA, GAT → 2; 4-mers: AGAT → 1. Confirms overlapping extraction with single-position advancement.
4. **Overlap:** k-mers are extracted by a sliding window of step 1, so adjacent k-mers overlap by k−1 characters.

### BioInfoLogics — k-mer counting, part I: Introduction

**URL:** https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
**Accessed:** 2026-06-14
**Authority rank:** 3 (technical reference accompanying an established k-mer-counting toolchain; used for the precise total/distinct/unique terminology and a worked frequency table)

**Retrieval:** WebSearch query `singleton unique k-mers appearing exactly once sequence definition`, then WebFetch of the article URL above.

**Key Extracted Points:**

1. **Total k-mers:** the sum of all k-mers extracted from a sequence counting duplicates; for length L there are L − k + 1 of them.
2. **Distinct k-mers:** "_Distinct k-mers_ are counted only once, even if they appear more times." (verbatim) — i.e. the number of different k-mer strings present.
3. **Unique k-mers:** "_Unique k-mers_ are those that appear only once." (verbatim) — k-mers whose frequency is exactly 1. This is the canonical definition `FindUniqueKmers` must implement.
4. **Worked example (ATCGATCAC, k=3, non-canonical):** 7 total 3-mers, 6 distinct, 5 unique. The 5 unique 3-mers occurring exactly once are TCG, CGA, GAT, TCA, CAC; ATC occurs twice (positions 0 and 4) and is therefore NOT unique.

### Compeau & Pevzner — Bioinformatics Algorithms: An Active Learning Approach (k-mer counting / Frequent Words)

**URL:** https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637 (catalogue/description page surfaced for the textbook)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed textbook)

**Retrieval:** WebSearch query `Compeau Pevzner Bioinformatics Algorithms k-mer Count Pattern Text frequency definition L-k+1`.

**Key Extracted Points:**

1. **Count(Text, Pattern):** the number of times a k-mer `Pattern` appears as a substring of `Text` (overlapping occurrences counted). This is the per-k-mer frequency on which both "unique" (Count = 1) and "min-count" (Count ≥ t) filters operate.
2. **Most-frequent / recurrent k-mers:** a k-mer is a most frequent k-mer if it maximises Count(Text, Pattern); selecting k-mers whose Count ≥ a threshold t is the standard way to isolate recurrent/over-represented k-mers (the basis for `FindKmersWithMinCount`).

### Jellyfish source — dump / stats (B06 review, 2026-10-01)

**URL:** https://raw.githubusercontent.com/gmarcais/Jellyfish/master/sub_commands/dump_main.cc, `dump_main_cmdline.yaggo`, `stats_main.cc` (opened)
**Authority rank:** 2 (reference implementation; Marçais & Kingsford 2011, *Bioinformatics* 27:764)

**Key Extracted Points:**

1. `dump`: `if(it.val() < lower_count || it.val() > upper_count) continue;` — inclusive range filter; defaults `lower_count = 0`, `upper_count = numeric_limits<uint64_t>::max()`. yaggo: `-L` "Don't output k-mer with count < lower-count", `-U` "Don't output k-mer with count > upper-count".
2. Output order is the database iteration (hash) order — no sorted order is promised; this API therefore documents its own deterministic order (count desc, ties ordinal k-mer).
3. `stats` "Unique" = `uniq += val == 1` — same definition of unique as BioInfoLogics; Distinct is separate.

### KMC 3 CLI (B06 review, 2026-10-01)

**URL:** https://raw.githubusercontent.com/refresh-bio/KMC/master/kmc_CLI/kmc.cpp (opened)

1. `-ci<value> - exclude k-mers occurring less than <value> times (default: 2)`; `-cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)` — lower and upper count cut-offs, as in Jellyfish.

### Reference cross-check (Python Counter replica of `jellyfish dump`)

BA1B sample ACGTTGCATGTCGCATGATGCATGAGAGCT, k=4: `-L1 -U1` → 17 k-mers (ACGT … TTGC, = Jellyfish stats Unique 17); `-L2` → CATG:3 GCAT:3 ATGA:2 TGCA:2; `-L2 -U2` → ATGA:2 TGCA:2; `-L3 -U3` → CATG:3 GCAT:3. GTAGAGCTGT k=2 unique → CT GA GC TA TG; `-L2` → AG:2 GT:2. AAAACGTAAA k=2 `-L2` → AA:5, `-L2 -U2` → none. Full table: docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md §7.2.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia — K-mer

1. **k > L:** when the k-mer length exceeds the sequence length, L − k + 1 ≤ 0, so the sequence contains zero k-mers.

### From BioInfoLogics

1. **Repeated k-mers excluded from "unique":** a k-mer that appears two or more times (e.g. ATC in ATCGATCAC) is distinct but NOT unique; the unique set is strictly the frequency-1 subset.

---

## Test Datasets

### Dataset: ATCGATCAC (BioInfoLogics worked example)

**Source:** BioInfoLogics — k-mer counting, part I (worked table, non-canonical k=3)

| Parameter | Value |
|-----------|-------|
| Sequence | ATCGATCAC |
| k | 3 |
| Total 3-mers | 7 (ATC, TCG, CGA, GAT, ATC, TCA, CAC) |
| Distinct 3-mers | 6 |
| Unique 3-mers (Count = 1) | 5 → {TCG, CGA, GAT, TCA, CAC} |
| Non-unique 3-mers | ATC (Count = 2) |

### Dataset: ACGTACGT (derived from the k-mer definition for min-count filtering)

**Source:** k-mer definition (Wikipedia L−k+1 + Compeau & Pevzner Count); occurrences enumerated directly.

| Parameter | Value |
|-----------|-------|
| Sequence | ACGTACGT |
| k | 4 |
| Total 4-mers | 5 (ACGT, CGTA, GTAC, TACG, ACGT) |
| Counts | ACGT=2, CGTA=1, GTAC=1, TACG=1 |
| FindKmersWithMinCount(..,4,2) | {(ACGT, 2)} |
| FindKmersWithMinCount(..,4,1) | all 4 distinct, ordered by Count desc (ACGT first, Count 2) |
| FindUniqueKmers(..,4) | {CGTA, GTAC, TACG} (the Count=1 set) |

### Dataset: AGAT (Wikipedia worked example, all distinct)

**Source:** Wikipedia — K-mer (AGAT example)

| Parameter | Value |
|-----------|-------|
| Sequence | AGAT |
| k | 2 |
| 2-mers | AG, GA, AT (all Count = 1) |
| FindUniqueKmers(AGAT, 2) | {AG, GA, AT} (all three unique) |

---

## Assumptions

1. **ASSUMPTION: minCount ≤ 0 behaviour** — Authoritative sources define min-count filtering only for a meaningful threshold (t ≥ 1, recurrent k-mers). For minCount ≤ 1 the predicate `Count ≥ minCount` is satisfied by every observed k-mer; the implementation returns all distinct k-mers ordered by count. This is the mathematically consistent extension of `Count ≥ t`, not an invented value, so it is treated as defined behaviour rather than a correctness-affecting unknown.
2. **ASSUMPTION: case normalisation** — sources use upper-case DNA; the implementation upper-cases input (consistent with sibling `KmerAnalyzer` methods) so that case variants count as the same k-mer. No source contradicts this.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindUniqueKmers(ATCGATCAC, 3)` returns exactly {TCG, CGA, GAT, TCA, CAC} (Count = 1 set; ATC excluded). — Evidence: BioInfoLogics worked table.
2. **MUST Test:** `FindUniqueKmers(AGAT, 2)` returns all three 2-mers {AG, GA, AT}. — Evidence: Wikipedia AGAT example.
3. **MUST Test:** `FindKmersWithMinCount(ACGTACGT, 4, 2)` returns exactly {(ACGT, 2)}. — Evidence: definition + enumerated occurrences.
4. **MUST Test:** `FindKmersWithMinCount(ACGTACGT, 4, 1)` returns all 4 distinct 4-mers ordered by Count descending, ACGT (Count 2) first. — Evidence: Compeau & Pevzner Count ≥ t + ordering.
5. **MUST Test:** homopolymer (e.g. AAAAA, k=3) has zero unique k-mers (the single distinct 3-mer AAA has Count = 3). — Evidence: definition (Count > 1 ⇒ not unique).
6. **SHOULD Test:** empty sequence and k > length return empty for both methods. — Rationale: L − k + 1 ≤ 0 ⇒ no k-mers (Wikipedia).
7. **SHOULD Test:** k ≤ 0 throws ArgumentOutOfRangeException. — Rationale: k-mer length must be positive (definition: substrings of length k).
8. **COULD Test:** case-insensitivity (lower-case input yields same unique set). — Rationale: documented normalisation assumption.

---

## Audit round 1, WP3 (B06, F12) — Jellyfish `count -C` + `dump -L/-U` and `histo` options (executed reference)

- **Sources opened** (raw.githubusercontent.com, gmarcais/Jellyfish master): `sub_commands/dump_main.cc` (`if(it.val() < lower_count || it.val() > upper_count) continue;`), `sub_commands/histo_main.cc` (`compute_histo`: `base = inc >= low ? 0 : low - inc`, `ceil = high + inc`, `nb_buckets = (ceil + inc - base) / inc`, `< base → histo[0]`, `> ceil → histo[nb_buckets-1]`, else `(val - base) / inc`; rows printed when `histo[i] > 0 || full`), `sub_commands/histo_main_cmdline.yaggo` (`-l` default 1, `-h` default 10000, `-i` default 1, `-f` "Full histo. Don't skip count 0."; "The last bucket in the output behaves as a catchall").
- **Reference program executed:** Jellyfish 2.3.1 (`apt jellyfish 2.3.1-3build1`), `count -m k -s 10000 -t 1 [-C]`, then `dump -c -L l [-U u]` (20 runs) and `histo [-l] [-h] [-i] [-f]` (72 runs). A Python replica of `histo_main.cc` reproduces all 72 histo rows.
- **Numbers (dump, `-C`):** BA1B k=4 `-L 1 -U 1` → 16 k-mers (AACG … GCGA); `-L 2` → ATGC:4 CATG:3 ATGA:2 TGCA:2; `-L 2 -U 3` → CATG:3 ATGA:2 TGCA:2; `acgtNNacgtacgRtTTGCAnA` k=3 `-L 1 -U 1` → AAA CAA; ATGATG k=3 → ATC TCA / ATG:2; Rosalind KMER k=4 `-L 7` → AACT:10 ACTC:10 ACTG:9 AGTC:9 CTCA:9 AGAC:7 AGTA:7 CAGC:7 GTGA:7. Full table: docs/algorithms/K-mer/Unique_And_MinCount_Kmers.md §7.3.
- **Numbers (histo):** Rosalind k=4 `-C` defaults → 1 23 2 34 3 27 4 25 5 7 6 5 7 4 9 3 10 2; `-h 5` → … 5 7 6 14 (cap bin); `-i 2` → 0 23 2 61 4 32 6 9 8 3 10 2; `-l 3 -h 8 -i 2` → 1 57 3 52 5 12 7 4 9 5. Table: docs/algorithms/K-mer/K-mer_Frequency_Analysis.md §7.3.
- **Implementation:** `FindUniqueKmers` / `FindKmersWithMinCount` option overloads (one 5-argument implementation over the option-aware `CountKmers` + the shared `SelectByCountRange`); `GetKmerHistogram` (sequence and count-table overloads). C# equals every row.

---

## References

1. Wikipedia contributors. 2026. *K-mer*. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Bernardo Clavijo et al. 2018. *k-mer counting, part I: Introduction*. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach* (2nd ed.). Active Learning Publishers. https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-10-01**: B06 audit round 1 WP3 — option-aware overloads and Jellyfish histo cross-check (executed Jellyfish 2.3.1).
