# Validation Report: MOTIF-SHARED-001 — Shared Motifs (fixed-length word enumeration with matching-sequence quorum)

- **Validated:** 2026-06-16   **Area:** Motif Discovery / Matching (`Seqeron.Genomics.Analysis`)
- **Canonical method(s):** `MotifFinder.FindSharedMotifs(IEnumerable<DnaSequence>, int k = 6, int minSequences = 2)` → `IEnumerable<SharedMotif>`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session, 2026-06-16)

1. **RSAT oligo-analysis manual** — https://rsat.eead.csic.es/plants/help.oligo-analysis.html (WebFetch).
   Verbatim retrieved:
   - Occurrences (occ): *"a simple count of the number of occurrences of each oligonucleotide."*
   - Matching sequences (mseq): *"the number of sequences from the input set which contain at least one occurrence of the oligonucleotide."*
   - For the matching-sequence probability: *"only the first occurrence of each sequence is taken into consideration."*
   - Size guidance: *"starting with an analysis of hexanucleotides (size=6), and scanning sizes between 4 and 8."*
   Confirms: (a) per-sequence presence/absence (count ≤ 1 per sequence per word); (b) fixed-length
   oligonucleotide enumeration; (c) the default k=6 sits inside the recommended range.

2. **ROSALIND "Finding a Shared Motif" (LCSM)** — https://rosalind.info/problems/lcsm/ (WebFetch).
   Verbatim retrieved:
   - *"A substring contained in all strings from a collection."* (common substring)
   - *"A common substring of a collection of maximum length."* (LCSM)
   - Sample dataset: `GATTACA`, `TAGACCA`, `ATACA`; sample output: `AC`.
   - Non-uniqueness: *"'AA' and 'CC' are both longest common substrings of 'AACC' and 'CCAA'."*
   Confirms LCSM is the *variable-length, present-in-all* problem — distinct from this unit's fixed-k
   quorum framing, exactly as the algorithm doc states.

3. **Das & Dai (2007), BMC Bioinformatics 8(S7):S21** and **van Helden et al. (1998), J Mol Biol
   281(5):827–842** — cited in the repo's Evidence as the word-enumeration family (exact matching, no
   variations within an oligonucleotide) and the named primary for oligo-analysis. These carry the
   "no variations allowed" exact-matching constraint and the over-representation framing. Not
   re-fetched this session (PMC/ScienceDirect); their role here is corroborative, and the operative
   verbatim definitions (matching sequences, common substring) come from RSAT and Rosalind above.

### Formula check
- Matching-sequence set `M(w) = { i : s_i contains ≥1 exact occurrence of w }`, report iff `|M(w)| ≥ q`.
  Matches RSAT mseq definition + Das & Dai quorum exactly.
- `Prevalence = |M(w)| / m`. RSAT reports the raw mseq count; expressing it as a fraction of total
  input sequences is a documented presentation convenience (Evidence Assumption 2). For a reported
  word `1 ≤ |M(w)| ≤ m`, so `Prevalence ∈ [1/m, 1] ⊂ (0,1]` — INV-04 holds as a true property.
- All five invariants (INV-01..05) are genuine mathematical properties of the definition.

### Edge-case semantics
Empty collection → no words; sequence shorter than k → no length-k window; word repeated within one
sequence → counts once (RSAT "at least one occurrence"); 1-substitution near-miss → different word
(exact matching); `k<1`/`minSequences<1` → out-of-range; null → null-arg. All defined and sourced.

### Independent cross-check (hand-computed this session)
Recomputed every dataset from scratch in Python from the RSAT/Rosalind definitions (not from the code):

| Dataset | k, q | Independent result |
|---------|------|--------------------|
| ATGATG, ATGCCC, CCCGGG | 3, 2 | ATG→{0,1}, CCC→{1,2}; below-quorum: GGG,GAT,TGA,GCC,TGC,CCG,CGG all count 1 → excluded; prevalence 2/3 |
| ACGT, ACTT | 4, 2 | no shared 4-mer (exact match) |
| AT, ATGAAT | 3, 2 | empty (AT has no length-3 window) |
| ATGAAA, CCATGC, GGGATG | 3, 3 | ATG→{0,1,2}, prevalence 1.0 |
| Rosalind LCSM (GATTACA,TAGACCA,ATACA) | 2, 3 | 2-mers in all 3 = {AC, CA, TA} = exactly the LCSMs of length 2 |

The last row is an independent confirmation that this unit's fixed-k quorum framing relates correctly
to (but is not identical to) Rosalind LCSM, matching the Evidence's documented distinction.

### Findings / divergences
None. Two documented assumptions (default k=6/q=2 as API ergonomics; prevalence-as-fraction) are
honest, in-range, and clearly flagged — not source-prescribed biological constants. PASS.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs:585-630`; record
`SharedMotif` at `:743-746`; `DnaSequence` normalization (`ToUpperInvariant`) at
`Core/.../DnaSequence.cs:30`.

### Formula realised correctly? (evidence)
- Per-sequence `HashSet<string> seenInSeq` ⇒ each word adds index `seqIdx` at most once → exact RSAT
  mseq semantics (`MotifFinder.cs:601-615`).
- Window loop `for i in 0..seq.Length-k` ⇒ sequences shorter than k contribute nothing (boundary).
- Quorum filter `seqIndices.Count >= minSequences` ⇒ INV-03 (`:622`).
- `Prevalence = (double)seqIndices.Count / seqList.Count` ⇒ INV-04 (`:627`).
- Exact substring match (`seq.Substring(i,k)` keyed in dictionary) over uppercase-normalized text ⇒
  exact, case-insensitive matching (INV-05). No degeneracy. Faithful to the validated description.
- Validation: null → `ArgumentNullException`; `k<1`/`minSequences<1` → `ArgumentOutOfRangeException`
  (`:590-592`). Method is `yield`-based, so checks fire on enumeration; tests force with `.ToList()`.

### Cross-verification table recomputed vs code
Every value in the Stage-A hand-computation table reproduces the actual code output (full suite green,
including all M1–M6/S1–S3/C1–C3 tests). No divergence between independent computation and code.

### Variant/delegate consistency
Single canonical method; no `*Fast`/instance-property variant for this unit. N/A.

### Test quality audit (`MotifFinder_FindSharedMotifs_Tests.cs`, 12 tests)
- **Sourced, exact assertions:** index sets (`Is.EqualTo(new[]{0,1})`), counts (`Is.EqualTo(2)`),
  prevalence (`2.0/3.0 .Within(1e-10)`, `1.0`), empties for negatives. No `Greater`/`AtLeast`/ranges
  where an exact value is known; no widened tolerances; nothing skipped/ignored.
- **Fails against a wrong impl:** M3 fails if a sequence is double-counted; M4 fails if quorum not
  applied; M6/S1 fail if matching were degenerate or windows mis-bounded; M5 pins exact 2/3.
- **Coverage:** all five invariants and every Stage-A branch — quorum, repeat-once, exact-match
  negative, k>len boundary, full quorum, empty, k<1, minSequences<1, null. C2b (minSequences=0)
  present beyond the spec's listed cases.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6606`; `dotnet build` 0 errors. The 4
  build warnings (NUnit2007) are in an unrelated file (`ApproximateMatcher_EditDistance_Tests.cs`),
  not the changed/validated SharedMotifs file.

**Gate result: PASS.** Minor non-blocking note: no test invokes the default-parameter overload
(k=6/minSequences=2); the defaults are documented (Evidence Assumption 1) as non-biological API
ergonomics, so this is a coverage nicety, not a defect — the core logic is fully exercised at explicit k.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — description matches RSAT/Rosalind verbatim definitions and independent hand
  computation; assumptions honest and in-range.
- **Stage B: PASS** — code faithfully realises the matching-sequence quorum; tests assert exact
  sourced values and cover all branches; full suite green, build clean.
- **End-state: CLEAN.** No defect found; no code or test change required.
