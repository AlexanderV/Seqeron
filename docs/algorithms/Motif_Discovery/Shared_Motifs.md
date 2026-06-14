# Shared Motifs (fixed-length word enumeration with matching-sequence quorum)

| Field | Value |
|-------|-------|
| Algorithm Group | Motif Discovery / Matching |
| Test Unit ID | MOTIF-SHARED-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

`FindSharedMotifs` identifies fixed-length words (oligonucleotides / k-mers) that recur across a
collection of DNA sequences. Each length-k word is scored by its **matching-sequence count** — the
number of input sequences that contain at least one exact occurrence of it [3] — and any word whose
matching-sequence count reaches a caller-supplied quorum is reported. It belongs to the word-based /
enumerative family of motif-finding methods (oligo-analysis) [1][2], and is exact (no degenerate or
substituted matches) [2]. It is deterministic. It is distinct from the longest-common-substring
("Finding a Shared Motif", LCSM) problem [4], which seeks a single variable-length substring present
in *all* sequences rather than fixed-k words meeting a quorum.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Transcription-factor binding sites and other regulatory signals appear as short, recurrent words
shared by co-regulated upstream regions. Word-enumeration methods enumerate all oligonucleotides of a
fixed length over the input set and rank them by how widely/often they occur [1][2]. The RSAT
oligo-analysis implementation of the van Helden method reports, for each word, both its raw occurrence
count and its "matching sequences" count [1][3].

### 2.2 Core Model

Let `S = {s_0, …, s_{m-1}}` be the input sequences and `k` the word length. For each word `w` of
length `k`, define the matching-sequence set `M(w) = { i : s_i contains ≥ 1 exact occurrence of w }`.
"Matching sequences" is "the number of sequences from the input set which contain at least one
occurrence of the oligonucleotide" [3]; a word repeated several times within one sequence still
contributes 1 to `|M(w)|` ("only the first occurrence of each sequence is taken into consideration") [3].
A word `w` is a shared motif iff `|M(w)| ≥ q`, where `q` (`minSequences`) is the quorum — the
word-enumeration "records the number of sequences containing occurrences of each k-mer" and reports
those over a threshold [2]. Matching is exact: van Helden's method allows "no variations within an
oligonucleotide" [2]. Prevalence is reported as `|M(w)| / m`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported word has length exactly k | only length-k windows are enumerated [3] |
| INV-02 | Sequence indices are distinct; each sequence counted at most once per word | per-sequence presence/absence, "at least one occurrence" [3] |
| INV-03 | A word is reported iff `|M(w)| ≥ minSequences` | quorum criterion [2] |
| INV-04 | `Prevalence = |M(w)| / m ∈ (0, 1]` | matching count over total sequences [3] |
| INV-05 | Matching is exact (a 1-substitution variant is a different word) | "no variations allowed within an oligonucleotide" [2] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Shared Motifs (this) | LCSM (Rosalind) [4] |
|--------|----------------------|---------------------|
| Word length | fixed k | variable (maximal) |
| Membership | quorum (≥ minSequences) | present in all sequences |
| Output | all words meeting quorum | the longest common substring(s) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequences | `IEnumerable<DnaSequence>` | required | Input sequences (uppercase-normalized by `DnaSequence`) | non-null |
| k | `int` | 6 | Word length | ≥ 1 |
| minSequences | `int` | 2 | Quorum (matching-sequence threshold) | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Sequence | `string` | The shared word (length k) |
| SequenceIndices | `IReadOnlyList<int>` | Distinct 0-based indices of sequences containing the word |
| Prevalence | `double` | `SequenceIndices.Count / totalSequences` ∈ (0,1] |

### 3.3 Preconditions and Validation

`sequences` null → `ArgumentNullException`. `k < 1` or `minSequences < 1` → `ArgumentOutOfRangeException`.
Input is 0-based; `DnaSequence` normalizes to uppercase A/C/G/T. Empty collection → no results. A
sequence shorter than k yields no words (no length-k window). Matching is exact and case-normalized.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize the sequences; record total count `m`.
2. For each sequence index `i`, slide a length-k window; collect the **distinct** words seen in `s_i`
   (a per-sequence `HashSet`) so each word contributes 1 to that sequence.
3. For each distinct word seen in `s_i`, append `i` to the word's matching-sequence list.
4. Emit every word whose matching-sequence list size ≥ `minSequences`, with its indices and prevalence.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

- Matching-sequence count (RSAT mseq): per-sequence presence, not occurrence multiplicity [3].
- Quorum filter: `|M(w)| ≥ minSequences` [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindSharedMotifs | O(Σᵢ (nᵢ − k + 1) · k) | O(distinct words · k) | nᵢ = length of sequence i; per-sequence HashSet of length-k words |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.FindSharedMotifs(IEnumerable<DnaSequence>, int, int)`: enumerates fixed-length words and reports those meeting the matching-sequence quorum.

### 5.2 Current Behavior

Each sequence is scanned once with a per-sequence `HashSet<string>` so a word repeated within a
sequence is counted once (matching-sequence semantics) [3]. Words are accumulated in a dictionary
`word → list of sequence indices`; results are filtered by the quorum and yielded lazily.

**Suffix tree decision — not used.** The repository `SuffixTree` is a single-text generalized suffix
tree; its `LongestCommonSubstring(other)` computes a *two-string* longest common substring, which
solves neither the fixed-k enumeration nor the k-sequence matching-sequence count needed here. The
required operation is "for every length-k word, in how many distinct sequences does it occur," which a
linear per-sequence window scan with a HashSet computes directly in O(Σᵢ nᵢ · k); building one suffix
tree per sequence would add construction overhead without changing the result. Therefore a direct
scan is used; correctness is unchanged.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Matching-sequence count = number of input sequences containing ≥ 1 exact occurrence of the word [3].
- Per-sequence single count (repeats within a sequence count once) [3].
- Fixed-length word enumeration over the input set [1][3].
- Quorum reporting (≥ minSequences) [2].
- Exact (non-degenerate) matching [2].

**Intentionally simplified:**

- Statistical significance: van Helden / RSAT also rank words by an over-representation P-value against
  a genome-wide background table [1]; **consequence:** this method reports the matching-sequence quorum
  only, without the P-value ranking, so users get presence-based shared words, not significance-ranked ones.

**Not implemented:**

- Degenerate / substituted matching (e.g. dyad-analysis, spacer words); **users should rely on:**
  `FindDegenerateMotif` for IUPAC-degenerate single-pattern search, or external RSAT for significance ranking.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default k=6, minSequences=2 | Assumption | Changes which words are reported | accepted | API defaults inside RSAT's allowed range; caller-supplied (Evidence ASSUMPTION 1) |
| 2 | Prevalence as a fraction | Assumption | Presentation of matching count | accepted | `|M(w)|/m`; RSAT reports the raw count (Evidence ASSUMPTION 2) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty collection | no results | no words to enumerate |
| Sequence shorter than k | contributes no words | no length-k window |
| Word repeated within one sequence | counts once toward its quorum | matching-sequence definition [3] |
| 1-substitution near-miss | treated as a different word | exact matching [2] |
| k < 1 / minSequences < 1 | `ArgumentOutOfRangeException` | contract |
| null collection | `ArgumentNullException` | contract |

### 6.2 Limitations

Exact words only — no mismatches, gaps, or degeneracy; no statistical significance / background model;
fixed k (does not find variable-length shared substrings — see LCSM [4]); does not consider the
reverse-complement strand.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

```csharp
var seqs = new[] { new DnaSequence("ATGATG"), new DnaSequence("ATGCCC"), new DnaSequence("CCCGGG") };
var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();
// "ATG" -> SequenceIndices [0,1], Prevalence 2/3 (ATG repeats in seq0 but counts once)
// "CCC" -> SequenceIndices [1,2], Prevalence 2/3
// "GGG" only in seq2 -> excluded (matching sequences 1 < 2)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_FindSharedMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_FindSharedMotifs_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [MOTIF-SHARED-001-Evidence.md](../../../docs/Evidence/MOTIF-SHARED-001-Evidence.md)
- Related algorithms: [Overrepresented_Kmer_Discovery](./Overrepresented_Kmer_Discovery.md)

## 8. References

1. van Helden J, André B, Collado-Vides J. 1998. Extracting regulatory sites from the upstream region of yeast genes by computational analysis of oligonucleotide frequencies. J Mol Biol 281(5):827–842. https://www.sciencedirect.com/science/article/abs/pii/S0022283698919477
2. Das MK, Dai HK. 2007. A survey of DNA motif finding algorithms. BMC Bioinformatics 8(Suppl 7):S21. https://pmc.ncbi.nlm.nih.gov/articles/PMC2099490/
3. RSAT. oligo-analysis manual (Regulatory Sequence Analysis Tools). https://rsat.eead.csic.es/plants/help.oligo-analysis.html (accessed 2026-06-14)
4. ROSALIND. Finding a Shared Motif (LCSM). https://rosalind.info/problems/lcsm/ (accessed 2026-06-14)
