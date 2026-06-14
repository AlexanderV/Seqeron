# Test Specification: MOTIF-SHARED-001

**Test Unit ID:** MOTIF-SHARED-001
**Area:** Matching (ProteinMotif/Matching — `Seqeron.Genomics.Analysis`)
**Algorithm:** Shared Motifs via fixed-length word enumeration with matching-sequence quorum
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | RSAT oligo-analysis manual (reference implementation; defines "matching sequences") | 3 | https://rsat.eead.csic.es/plants/help.oligo-analysis.html | 2026-06-14 |
| 2 | Das & Dai (2007), A survey of DNA motif finding algorithms, BMC Bioinformatics 8(S7):S21 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2099490/ | 2026-06-14 |
| 3 | van Helden, André, Collado-Vides (1998), J Mol Biol 281(5):827–842 | 1 | https://www.sciencedirect.com/science/article/abs/pii/S0022283698919477 | 2026-06-14 (403; named primary) |
| 4 | ROSALIND LCSM (alternative framing, not implemented) | 4 | https://rosalind.info/problems/lcsm/ | 2026-06-14 |

### 1.2 Key Evidence Points

1. "Matching sequences" = "the number of sequences from the input set which contain at least one occurrence of the oligonucleotide." — RSAT oligo-analysis manual.
2. Each input sequence is counted at most once per word ("at least one occurrence"; "only the first occurrence of each sequence is taken into consideration"). — RSAT.
3. Fixed-length words (oligonucleotides) are enumerated across the entire input set; only exact (non-degenerate) words match — "no variations allowed within an oligonucleotide." — Das & Dai (2007).
4. The shared-motif decision is a quorum: report words appearing in ≥ a threshold number of sequences ("records the number of sequences containing occurrences of each k-mer"). — Das & Dai word-enumeration family.
5. LCSM (single longest common substring present in *all* sequences) is a distinct algorithm; this unit fixes k and uses a quorum, so it is not LCSM. — Rosalind LCSM.

### 1.3 Documented Corner Cases

- A word repeated many times in one sequence still adds 1 to its matching-sequence count (RSAT).
- Exact match only; a one-substitution near-miss is a different word (Das & Dai).

### 1.4 Known Failure Modes / Pitfalls

1. Counting total occurrences instead of distinct sequences would inflate the matching-sequence count — RSAT definition.
2. Allowing degenerate/substituted matches would over-report — Das & Dai "no variations allowed."

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindSharedMotifs(IEnumerable<DnaSequence>, int k = 6, int minSequences = 2)` | `MotifFinder` | **Canonical** | Fixed-length word enumeration; matching-sequence quorum |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned `SharedMotif.Sequence` has length exactly k | Yes | Fixed-length enumeration (RSAT oligo length) |
| INV-2 | `SequenceIndices` are distinct 0-based input indices; each sequence appears at most once per motif | Yes | RSAT "matching sequences" (at least one occurrence per sequence) |
| INV-3 | A motif is returned iff its matching-sequence count ≥ minSequences | Yes | Quorum criterion (Das & Dai) |
| INV-4 | `Prevalence` = (matching-sequence count) / (total input sequences), in (0,1] | Yes | RSAT definition + Evidence Assumption 2 |
| INV-5 | Matching is exact: a word differing by ≥1 substitution is a distinct word | Yes | Das & Dai "no variations allowed within an oligonucleotide" |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Quorum word reported with exact indices | S0=ATGATG, S1=ATGCCC, S2=CCCGGG; k=3, minSeq=2; word ATG | ATG returned, SequenceIndices = [0,1] | RSAT matching sequences |
| M2 | Second quorum word | same input; word CCC | CCC returned, SequenceIndices = [1,2] | RSAT matching sequences |
| M3 | Repeats in one sequence count once | S0=ATGATG (ATG twice) → still contributes 1 | ATG SequenceIndices contain 0 exactly once; count = 2 (with S1) | RSAT "at least one occurrence" |
| M4 | Below-quorum word excluded | same input; word GGG only in S2 | GGG NOT returned (count 1 < 2) | Quorum (Das & Dai) |
| M5 | Prevalence exact | ATG matching 2 of 3 sequences | Prevalence = 2.0/3.0 (.Within 1e-10) | RSAT def + INV-4 |
| M6 | Exact-word semantics | S0=ACGT, S1=ACTT (differ at pos 2); k=4, minSeq=2 | No length-4 motif shared (ACGT≠ACTT) | Das & Dai no variations |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | k longer than a sequence | S0=AT (len2), S1=ATGAAT; k=3, minSeq=2 | S0 yields no windows → no shared 3-mer | Window boundary i ≤ len−k |
| S2 | minSequences = full set | three seqs all containing ATG; minSeq=3 | ATG returned with all three indices | Quorum at full quorum |
| S3 | INV-1 length | all returned motifs have length k | Length == k for every result | Fixed-length enumeration |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty collection | no sequences | empty result | Input validation |
| C2 | k < 1 | k = 0 | throws ArgumentOutOfRangeException | Implementation contract |
| C3 | null collection | null | throws ArgumentNullException | Implementation contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/`. Existing `MotifFinder_DiscoverMotifs_Tests.cs` covers the sibling unit MOTIF-DISCOVER-001 only. No existing tests for `FindSharedMotifs`. New canonical file `MotifFinder_FindSharedMotifs_Tests.cs` created.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| C1 | ❌ Missing | new |
| C2 | ❌ Missing | new |
| C3 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_FindSharedMotifs_Tests.cs` — all cases above.
- **Remove:** none (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_FindSharedMotifs_Tests.cs` | Canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |
| 11 | C2 | ❌ Missing | Implemented | ✅ Done |
| 12 | C3 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact indices asserted |
| M2 | ✅ Covered | exact indices asserted |
| M3 | ✅ Covered | repeat-once verified |
| M4 | ✅ Covered | below-quorum excluded |
| M5 | ✅ Covered | exact prevalence |
| M6 | ✅ Covered | exact-word negative |
| S1 | ✅ Covered | k>len boundary |
| S2 | ✅ Covered | full quorum |
| S3 | ✅ Covered | INV-1 length |
| C1 | ✅ Covered | empty collection |
| C2 | ✅ Covered | k<1 throws |
| C3 | ✅ Covered | null throws |

**✅ count = 12 = total in-scope cases.**

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default parameters k=6, minSequences=2 are API ergonomics within RSAT's allowed ranges, not biological constants | Evidence Assumptions §1; defaults not asserted as required values |
| 2 | Prevalence = matchingSequences / totalSequences (presentation of the matching-sequence count) | INV-4, M5 |

---

## 7. Open Questions / Decisions

1. **Decision:** The unit implements the word-enumeration / matching-sequence quorum framing (van Helden / RSAT), NOT Rosalind LCSM. LCSM is documented as a related-but-distinct algorithm in the algorithm doc.
2. **Decision (suffix tree):** Suffix tree not used — see algorithm doc §5.2. The repo SuffixTree is single-text; its `LongestCommonSubstring` is a two-string LCS and does not compute fixed-k matching-sequence counts across k sequences. A per-sequence HashSet word scan is the correct O(Σ(nᵢ)·k) approach.
