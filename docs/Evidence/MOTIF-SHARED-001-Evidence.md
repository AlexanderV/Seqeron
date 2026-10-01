# Evidence Artifact: MOTIF-SHARED-001

**Test Unit ID:** MOTIF-SHARED-001
**Algorithm:** Shared Motifs via fixed-length word enumeration with matching-sequence quorum (oligo-analysis "matching sequences")
**Date Collected:** 2026-06-14

---

## Online Sources

### RSAT oligo-analysis manual — output statistics definitions (reference implementation)

**URL:** https://rsat.eead.csic.es/plants/help.oligo-analysis.html
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 3 (reference implementation: RSAT is the published implementation of the van Helden oligo-analysis method)

**Key Extracted Points:**

1. **Occurrences (occ):** Retrieved text: "a simple count of the number of occurrences of each oligonucleotide. Overlapping matches are detected and summed in the counting."
2. **Matching sequences (mseq):** Retrieved text — verbatim: "the number of sequences from the input set which contain at least one occurrence of the oligonucleotide." This is the central statistic for the shared-motif (quorum) decision: each input sequence is counted at most once per oligonucleotide.
3. **Fixed-length enumeration:** Retrieved text: the oligonucleotide size parameter analyses "with oligonuleotides of any size between 1 and 8"; a fixed oligo length (k) is enumerated and counted throughout the entire input sequence set.
4. **Per-sequence counting for the matching-sequence statistic:** Retrieved text: for matching-sequence probability calculations "only the first occurrence of each sequence is taken into consideration" — i.e. presence/absence per sequence, not multiplicity, drives the matching-sequence count.

### RSAT source code — `perl-scripts/oligo-analysis` (reference implementation, review 2026-09-29)

**URL:** https://raw.githubusercontent.com/rsa-tools/rsat-code/master/perl-scripts/oligo-analysis (opened with curl, 2026-09-29)

**Key Extracted Points:**

1. Matching sequences: `### matching sequences (only count first occurrence per sequence)` — per window `$current_mseq{$pattern_seq} = 1` (with `-2str`, a word whose reverse complement is already marked is not re-marked), then after each sequence `foreach $pattern_seq (keys %current_mseq) { $patterns{$pattern_seq}->{mseq} += 1; }`.
2. Overlapping windows are the default (`$noov="-ovlp"`); `-1str` counts the given strand only. `FindSharedMotifs` = mseq under `-1str -ovlp`.
3. Cross-check: independent Python port of this loop on ACGTACGTTAGC / TTACGTAGCAAC / GGTAGCACGTTT / CATTTTACG, k=4, q=2 → ACGT{0,1,2}, CGTA{0,1}, TACG{0,1,3}, CGTT{0,2}, TAGC{0,1,2}, TTAC{1,3}, GTAG{1,2}, AGCA{1,2}; Rosalind LCSM sample, k=2, q=3 → TA, AC, CA.

### A survey of DNA motif finding algorithms (Das & Dai, 2007) — word-enumeration family

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2099490/
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 1 (peer-reviewed survey, BMC Bioinformatics 8(Suppl 7):S21)

**Key Extracted Points:**

1. **Oligo-analysis (van Helden et al. 1998):** Retrieved text describes it as a word-based / enumerative method whose methodology involves "(1) constitution of regulatory families and (2) calculation of expected oligonucleotide frequencies."
2. **Exact words only:** Retrieved text — the stated limitation: "The greatest shortcomings of the algorithm of van Helden et al. is that there are no variations allowed within an oligonucleotide." This justifies exact (non-degenerate) k-mer matching in the shared-motif word-enumeration approach.
3. **Quorum across sequences:** The companion search summary (WebSearch, 2026-06-14, query "shared motif finding common k-mer across multiple DNA sequences bioinformatics algorithm", BMC survey) states one approach "records the number of sequences containing occurrences of each k-mer" — i.e. count input sequences per k-mer and report those over a threshold.

### van Helden, André & Collado-Vides (1998) — primary source identification

**URL:** WebSearch query "van Helden Andre Collado-Vides 1998 extracting regulatory sites oligonucleotide frequencies Journal Molecular Biology" (2026-06-14); result item links to https://www.sciencedirect.com/science/article/abs/pii/S0022283698919477 (HTTP 403 on direct fetch — abstract captured from the search result summary).
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary article, J Mol Biol 281(5):827–842)

**Key Extracted Points:**

1. **Method basis:** Retrieved search summary: the method "is based on the detection of over-represented oligonucleotides, with statistical significance defined based on tables of oligonucleotide frequencies observed in all non-coding sequences." Confirms fixed-length oligonucleotide (word) counting as the basis of the method that RSAT oligo-analysis implements.
2. **Citation:** van Helden, J., André, B. & Collado-Vides, J. (1998). "Extracting regulatory sites from the upstream region of yeast genes by computational analysis of oligonucleotide frequencies." J Mol Biol 281(5), 827–842. (Direct article body could not be opened — HTTP 403; cited as the named primary, with the RSAT manual and Das & Dai survey carrying the operative verbatim definitions.)

### ROSALIND "Finding a Shared Motif" (LCSM) — contrasted alternative framing

**URL:** https://rosalind.info/problems/lcsm/
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 4 (educational problem set; used only to delineate the *alternative* "longest common substring" framing that this unit does NOT implement)

**Key Extracted Points:**

1. **Common substring:** Retrieved text: "A substring contained in all strings from a collection."
2. **Longest common substring:** Retrieved text: "A common substring of a collection of maximum length."
3. **Non-uniqueness:** Retrieved text: "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".
4. **Example:** "CG" is a common substring of "ACGTACGT" and "AACCGTATA"; "CGTA" is a longer common substring of both.
5. **Scope note:** LCSM requires a substring present in *all* sequences with variable length. The unit under test instead fixes the length (k) and uses a *quorum* (≥ minSequences), so LCSM is documented as a related-but-distinct algorithm, not the contract here. LCSM itself is available since B05 F26 as `SuffixTree.FindLongestCommonSubstrings(texts, minSupport)` (MCP `suffix_tree_k_common_substrings`; Rosalind sample → {AC, CA, TA}).

---

## Documented Corner Cases and Failure Modes

### From RSAT oligo-analysis manual

1. **Each sequence counted once per word:** The matching-sequence statistic counts presence/absence per sequence ("at least one occurrence"; "only the first occurrence of each sequence is taken into consideration"), so a word repeated many times within one sequence still contributes exactly 1 to its matching-sequence count.
2. **Overlapping occurrences:** For the raw occurrence count, "Overlapping matches are detected and summed" — relevant to occurrence counting, but the matching-sequence count is unaffected (still 1 per sequence).

### From Das & Dai (2007)

1. **No variations within a word:** Exact word matching only in the plain overloads; degenerate (`-oneN` / `-onedeg`) matching-sequence statistics are available since B05 F33 via `FindSharedMotifs(…, OligoDegeneracy)` (see the degenerate-word section below).

---

## Test Datasets

### Dataset: Hand-traced k=3 quorum example (derived from RSAT matching-sequence definition)

**Source:** Definition from RSAT oligo-analysis manual (matching sequences = number of input sequences containing ≥ 1 occurrence). Trace computed by hand.

Sequences (0-based index):
- S0 = `ATGATG`
- S1 = `ATGCCC`
- S2 = `CCCGGG`

| Parameter | Value |
|-----------|-------|
| k | 3 |
| minSequences | 2 |
| Word `ATG` | in S0 (pos 0,3), in S1 (pos 0) → matching sequences = {0,1}, count 2 |
| Word `CCC` | in S1 (pos 3), in S2 (pos 0) → matching sequences = {1,2}, count 2 |
| `ATG` prevalence | 2/3 |
| Shared motifs (count ≥ 2) | {`ATG`, `CCC`}; `ATG` SequenceIndices = [0,1]; `CCC` = [1,2] |

### Dataset: Rosalind LCSM sample (used only to assert this unit is NOT LCSM)

**Source:** https://rosalind.info/problems/lcsm/

| Parameter | Value |
|-----------|-------|
| Sequences | GATTACA, TAGACCA, ATACA |
| LCSM answer | `AC` (length 2; an LCSM-specific output) |
| This unit at k=2, minSeq=3 | reports all 2-mers present in all three (e.g. `AC`, `TA`) — fixed-k quorum, not a single longest substring |

---

## Assumptions

1. **ASSUMPTION: Default parameters k=6, minSequences=2** — The defaults are API ergonomics, not biological constants. RSAT permits any oligo length 1–8 and any quorum; the unit's defaults sit inside that range but are not prescribed by a source. Changing them changes which words are reported, so they are documented but treated as caller-supplied parameters (the algorithm is correct for any valid k ≥ 1 and minSequences ≥ 1).
2. **ASSUMPTION: Prevalence = matchingSequences / totalSequences** — RSAT reports raw matching-sequence counts; expressing it as a fraction of total input sequences is a presentation convenience consistent with the definition (a value in (0,1]). Not a source formula, so flagged.

---

## Recommendations for Test Coverage

1. **MUST Test:** A word present in exactly the quorum number of distinct sequences is reported with the correct `SequenceIndices` set and matching-sequence count — Evidence: RSAT "matching sequences" definition.
2. **MUST Test:** A word repeated multiple times within a single sequence contributes 1 (not its occurrence count) to the matching-sequence count — Evidence: RSAT "at least one occurrence" / "only the first occurrence of each sequence."
3. **MUST Test:** A word present in fewer than minSequences sequences is excluded — Evidence: quorum criterion (Das & Dai word-enumeration family).
4. **MUST Test:** `Prevalence` = matchingSequences / totalSequences exactly — Evidence: definition + Assumption 2.
5. **SHOULD Test:** Exact-word semantics — a near-miss (one substitution) is not matched — Rationale: Das & Dai "no variations allowed within an oligonucleotide."
6. **SHOULD Test:** k longer than the shortest sequence yields no words from that sequence (window loop empty) — Rationale: boundary of window enumeration.
7. **COULD Test:** Empty input collection returns no motifs; k < 1 throws — Rationale: input validation contract.

---

## RSAT oligo-analysis run (review 2026-09 follow-up, 2026-09-30)

**Source opened:** git clone of rsa-tools/rsat-code master 10043f2 (2026-09-23): `perl-scripts/oligo-analysis` (v1.169),
`perl-scripts/lib/RSA.disco.lib` (`NbPossibleOligos`, `MultiTestCorrections`, `GroupRC`), `perl-scripts/lib/RSA.lib`
(`SumExpectedFrequencies`, `ReadPatternFrequencies`), `perl-scripts/lib/RSA.seq.lib` (`SmartRC`),
`perl-scripts/lib/RSAT/stats.pm` (`sum_of_binomials`, `binomial_boe`), `perl-scripts/lib/RSAT/MarkovModel.pm`
(`load_from_file_oligos`, `add_pseudo_freq`, `normalize_transition_frequencies`, `segment_proba`).

**How it was run:** `RSAT=<clone> perl -I perl-scripts/lib perl-scripts/oligo-analysis -i x.fa -l k ...` (an empty
`RSAT_config.props`; no other installation needed), plus a copy with one debug line after `MultiTestCorrections`
printing every pattern's exp_freq, occ_P, occ_E, occ_sig, exp_ms, ms_P, ms_E, ms_sig with `%.17g`.

**Source facts established (code, confirmed by the runs):**
- Trials n = `sum_occurrences` = all N − k + 1 windows (also with `-noov`, where overlaps are added back); `-2str` pairs
  are grouped before summing, so n is not doubled.
- occ_P = `sum_of_binomials(exp_freq, sum_occurrences, occ, sum_occurrences)` (right tail, exact summation).
- occ_E = occ_P × `$nb_tested_patterns` — the number of patterns reaching `CalcProba`, i.e. after `-lth occ`
  (the manual's "NPO" equals it only with `-zeroocc`); occ_sig = −log10 occ_E (350 when occ_E = 0).
- ms_P = `binomial_boe(exp_ms/S, S, mseq)`, exp_ms = S·(1 − (1 − exp_freq)^(nb_pos/S)), S = all sequences,
  nb_pos over sequences with L ≥ k; ms_E = ms_P × `nb_possible_oligos` (4^k; with `-2str` (4^k + 4^{k/2})/2 for even k, 4^k/2 for odd k).
- `-markov m` (slow mode): sub-word frequencies of all overlapping (m+1)- and m-mers, single strand;
  exp_freq = f(w[0..m]) ∏ f(w[o..o+m]) / f(w[o..o+m−1]); m ≤ k − 2 for m > 0.
- `-2str`: occ summed with the reverse complement (palindromes once); input Bernoulli pools complementary residues.
  The code copies the kept member's exp_freq to its partner (2·p(min)), contradicting the manual's
  exp_freq(W) + exp_freq(W') for asymmetric backgrounds (e.g. `-markov 2` gives exp_freq 0 for observed pairs);
  the documented sum is implemented.
- `-bgfile` probabilities are returned as `%5g` strings by `segment_proba` (6 significant digits).

**Independent cross-checks:** Python port of the above (scipy.stats.binom.sf) = RSAT to ≤ 2e-12 relative on 68 runs
(t2/long sequences, k = 3/4/6, equi/input/`-markov 0..4`/`-noov`, 1str/2str, mseq on 4 sequences) and to %5g on 18
`-bgfile` runs, tables of order 0–2 (strand-sensitive and -insensitive); C# = port to ≤ 3e-13 on 400 random configurations.
Binomial tail vs mpmath exact sums ≤ 1e-13 (scipy ≤ 5e-14).

**Locked values (see TestSpec §8):** ACGTACGTTAGC / TTACGTAGCAAC / GGTAGCACGTTT / CATTTTACG, k = 4, `-1str` input
Bernoulli: ACGT mseq 3, exp_ms 0.12225827585979898, ms_P 0.00011159466135309564, ms_E 0.028568233306392483,
ms_sig 1.5441166160753867; `-2str`: CGTA|TACG mseq 3, ms_P 0.00087319490302470882, ms_sig 0.9253498996301337.

---

## Degenerate-word matching-sequence statistics (B05 F33, 2026-10-01)

Overload `FindSharedMotifs(IEnumerable<DnaSequence>, int k, int minSequences, OligoBackgroundModel, OligoStrandMode, double pseudoFrequency, OligoDegeneracy degeneracy)` (`MotifFinder.SharedMotifsDegenerate.cs`) applies RSAT `-oneN` / `-onedeg` to mseq / ms_P: a degenerate word's matching sequences are the **union** of its words' sequence sets (RSAT 1.169 sums them — documented deviation), expected frequency = Σ over matching words; `OligoDegeneracy.None` = the plain overload exactly. Reference: RSAT `oligo-analysis` 1.169 run with the env-guarded mseq-union fix (B05.md F33). MCP: `shared_motifs_significance` optional `degenerate` ('none' | 'oneN' | 'onedeg'). Tests and spec: `Unit/Analysis/MotifFinder_SharedMotifsDegenerate_Tests.cs`; TestSpec [MOTIF-DISCOVER-001 §11](../../tests/TestSpecs/MOTIF-DISCOVER-001.md) (S1–S7, M1–M3).

---

## References

1. van Helden J, André B, Collado-Vides J. (1998). Extracting regulatory sites from the upstream region of yeast genes by computational analysis of oligonucleotide frequencies. J Mol Biol 281(5):827–842. https://www.sciencedirect.com/science/article/abs/pii/S0022283698919477
2. Das MK, Dai HK. (2007). A survey of DNA motif finding algorithms. BMC Bioinformatics 8(Suppl 7):S21. https://pmc.ncbi.nlm.nih.gov/articles/PMC2099490/
3. RSAT — oligonucleotide analysis (oligo-analysis) manual. Regulatory Sequence Analysis Tools. https://rsat.eead.csic.es/plants/help.oligo-analysis.html (accessed 2026-06-14)
4. ROSALIND. Finding a Shared Motif (LCSM). https://rosalind.info/problems/lcsm/ (accessed 2026-06-14; cited to delineate the alternative LCSM framing, which this unit does not implement — the generalized-suffix-tree LCS of k strings is implemented separately since 2026-10-01 (B05 F26): `SuffixTree.FindLongestCommonSubstrings(texts, minSupport)` / MCP `suffix_tree_k_common_substrings`)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-29**: Review 2026-09 — RSAT source-code mseq loop added with Python cross-check values.
- **2026-09-30**: RSAT oligo-analysis significance / Markov / `-2str` implemented; RSAT run + Python port + mpmath references (review 2026-09 follow-up, B05).
- **2026-10-01**: LCSM pointer → implemented (B05 F26, `SuffixTree.FindLongestCommonSubstrings` / `suffix_tree_k_common_substrings`); degenerate-word overload `FindSharedMotifs(…, OligoDegeneracy)` cross-referenced (B05 F33, spec in TestSpec MOTIF-DISCOVER-001 §11).
