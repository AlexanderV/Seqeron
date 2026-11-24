# Validation Report: GENOMIC-MOTIFS-001 — Known Motif Search

- **Validated:** 2026-06-16   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindKnownMotifs(DnaSequence, IEnumerable<string>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS-WITH-NOTES. **State:** FIXED (not 'no defect found': the 2026-06 verdict below is superseded).
- **Fixes (B09.md, F-numbers):** F18 `FindMotif` ascending order and delegation to `MotifFinder.FindExactMotif`, F19 IUPAC-degenerate motifs and the reverse strand (`FindKnownMotifHits`), F20 KMP scans below 128 queries / suffix tree from 128, F49 mismatch-tolerant search (`maxMismatches`), F50 N / IUPAC subjects (string overloads).
- **Audit round 2 (A2-18 / A2-32 / A2-35, F91 / F96 / F99):** fuzznuc's "overhanging hit" artefacts decided against the binary (300 runs; `FindKnownMotifHits` reports none), the fuzznuc pattern language (`FuzznucPattern`, `FindPatternHits`, pattern files; 3 300 random syntax-tree rows vs fuzznuc 6.6.0, 0 disagreements), string overloads; MCP `find_pattern_hits`.
- **Oracle numbers:** EMBOSS fuzznuc 6.6.0 = Python `re` on 2 448 pairs; F49 fuzznuc `-pmismatch` oracle (see the Evidence file).
- **Open items:** registered, not doable here: `MotifFinder.FindExactMotif` first-query suffix-tree policy (23 s / 1.7 GB at 10 Mbp; needs a 'tree built?' accessor on `DnaSequence`, QUEUE R19). The statement 'a non-ACGT motif character simply does not match - documented limitation' in the original pass is superseded by F19 / F49 / F50.
- The unit's TestSpec (`tests/TestSpecs/GENOMIC-MOTIFS-001.md`), Evidence (`docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## Summary

GENOMIC-MOTIFS-001 ("Known Motif Search") is the classical **exact set-matching** problem:
given a subject sequence T and a set of fixed query motifs, report for each motif the 0-based
start positions of **all** its occurrences, including overlapping ones. It is *not* de novo
motif discovery / the (l,d) motif problem — those concern enriched/approximate motifs and are
covered by other units (MOTIF-DISCOVER-001, MOTIF-SHARED-001, etc.). This unit is exact
substring matching of a curated motif set against one reusable suffix-tree index.

## Stage A — Description

### Sources opened & what they confirm

1. **Rosalind "SUBS — Finding a Motif in DNA"** (https://rosalind.info/problems/subs/, retrieved
   this session via WebFetch). The problem is "find all starting positions where a substring
   appears within a larger DNA string", reporting **every** occurrence including overlapping ones.
   Sample: motif `ATAT` in `GATATATGCATATACTT` → output `2 4 10` (1-based). Positions 2 and 4 are
   overlapping. This directly confirms INV-01 (exact match), INV-02 (overlaps all reported), and
   the "set of all start positions" model. (Rosalind uses 1-based; the repo uses 0-based — see
   convention note below.)
2. **Wikipedia "EcoRI"** (https://en.wikipedia.org/wiki/EcoRI, retrieved this session). Quote:
   "The nucleic acid recognition sequence where the enzyme cuts is G↓AATTC". Confirms the
   biological motif `GAATTC` (length 6) used in dataset M2. (The Evidence file cited the
   "Restriction site" article; I cross-checked the dedicated EcoRI article and it agrees.)
3. **Gusfield exact-matching definition** (cited in Evidence via Tufts COMP 150GEN). The
   exact-matching problem is "find all occurrences of a pattern string P ... in a text string T";
   `P=aaa`, `T=aaaaa` → three overlapping occurrences at 0-based 0,1,2. Independently corroborated
   by the Rosalind overlapping example above.

### Formula / model check

The model is `{ i : T[i..i+n-1] == P }`, 0-based, overlaps included, per-motif independent sets,
absent motifs omitted, keys upper-cased. Every clause matches the cited sources. No log base,
normalization, or scoring is involved (exact match, not scored/approximate), so there is no
formula subtlety to mis-state.

### Definitions & conventions

- **Coordinate base:** 0-based start positions. The doc/spec state this explicitly and
  consistently. Rosalind is 1-based, but that is a presentation choice; the underlying set
  {1,3,9} (0-based) = {2,4,10} (1-based) is identical. Internal consistency holds.
- **Case handling:** `DnaSequence` upper-cases on construction (src line 30) and validates ACGT;
  `FindKnownMotifs` upper-cases each motif and keys the result by the upper-cased form
  (Biopython processes DNA upper-cased). Consistent and sourced.
- **Overlap convention:** all occurrences reported (not greedy non-overlapping). Sourced.

### Edge-case semantics

- Empty motif set → empty dictionary (no motif to search). Defined.
- Absent motif → omitted (empty occurrence set). Defined, INV-04.
- Empty/whitespace motif → skipped. This is an **assumption** (documented): no authoritative
  source defines an empty-pattern occurrence set as meaningful; a suffix-tree `FindAllOccurrences("")`
  would return every position, which is not a motif. The repo policy (skip) is reasonable and
  matches `FindMotif`. Accepted as an API-shape policy, not an algorithmic claim.
- `motifs` null → `ArgumentNullException`. Standard input validation.

### Independent cross-check (numbers)

Reference implementation: a naive overlapping exact-match scan
(`[i for i in range(len(T)-len(P)+1) if T[i:i+len(P)]==P]`), run this session:

| Dataset | Motif | Reference (0-based) | Spec expected |
|---------|-------|---------------------|----------------|
| `AAAAA` | `AAA` | `[0,1,2]` | `[0,1,2]` ✓ |
| `GAATTCAAAGAATTC` | `GAATTC` | `[0,9]` | `[0,9]` ✓ |
| `ACGTACGTAA` | `ACGT` | `[0,4]` | `[0,4]` ✓ |
| `ACGTACGTAA` | `AA` | `[8]` | `[8]` ✓ |
| `ACGTACGTAA` | `TTT` | `[]` | omitted ✓ |
| `AAAAA` | `A` | `[0,1,2,3,4]` | `[0,1,2,3,4]` ✓ |
| `ACAA` | `AA` | `[2]` | `[2]` ✓ |
| `GAATTC` | `GAATTC` | `[0]` | `[0]` ✓ |
| `ACGT` | `TTTT` | `[]` | omitted ✓ |
| `GATATATGCATATACTT` (Rosalind) | `ATAT` | `[1,3,9]` = 1-based `2 4 10` | matches Rosalind ✓ |

All values trace to externally-retrieved sources (Rosalind sample output; EcoRI site) and an
independent reference implementation, not to the repo's own code.

### Findings / divergences

None. Description is biologically and mathematically correct. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed

- `GenomicAnalyzer.FindKnownMotifs` — `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:190-228`.
- `SuffixTree.FindAllOccurrences` — `src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs:81-93`
  (and span overload :148). Collects all leaves under the matched node via `CollectLeaves`
  (DFS), so all start positions (overlapping included) are returned, unsorted.
- `DnaSequence` ctor — `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:22-31`
  (upper-case + ACGT validation); `SuffixTree` property :48 (built once, cached).

### Formula realised correctly?

Yes. For each non-empty motif: upper-case → `tree.FindAllOccurrences(normalized)` (all occurrences,
overlaps included) → if non-empty, `OrderBy(p => p)` (ascending sort for deterministic contract,
INV-03) → store under upper-cased key. Null guard throws `ArgumentNullException`. Empty/whitespace
skip via `string.IsNullOrWhiteSpace` (also covers a null element). Duplicate upper-cased keys
collapse (`result.ContainsKey` guard). Every Stage-A clause is realised.

### Cross-verification table recomputed vs code

The 12 unit tests assert the exact-value table above and all pass (see below). The values were
sourced independently first, then confirmed equal to actual code output — not the reverse.

### Variant/delegate consistency

Only one public entry point (`FindKnownMotifs`); no `*Fast`/instance variant. It delegates to the
shared `SuffixTree` used across the library. Consistent.

### Numerical robustness

Integer positions only; no floating point, overflow, or div-by-zero. Empty subject → empty suffix
tree → no occurrences → empty result. Non-ACGT motif characters (2026-06 text: "would simply not match,
documented limitation" - superseded 2026-11: IUPAC-degenerate motifs, both strands, mismatches and N / IUPAC subjects
are implemented, B09 F19 / F49 / F50). Robust on the stated domain.

### Test quality audit (HARD gate)

File: `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindKnownMotifs_Tests.cs` (12 tests).

- **Sourced, not code-echoes:** every expected value matches the externally-sourced /
  independently-recomputed table above (Rosalind overlap definition, EcoRI site, naive reference).
  A deliberately-wrong implementation (e.g. greedy non-overlapping, or DFS-unsorted output) would
  FAIL M1/S4 (overlap) and M4 (ascending) respectively. Not code echoes.
- **No greenwashing:** all assertions use exact `Is.EqualTo(new[]{…})` or `Is.Empty` /
  `ContainsKey(...) Is.False`. No `Greater`/`AtLeast`/`Contains`-only/ranges. No widened
  tolerances, no skip/ignore. M4 pairs `Is.Ordered.Ascending` with an exact-value check.
- **Coverage of all branches:** overlapping (M1, S4), multi-motif per-motif sets (M3), absent
  omitted (M5 + TTT in M3), ascending sort (M4), empty motif set (S1), case normalization +
  upper-cased key (S2), empty/whitespace skip (S3), duplicate-key dedup (C1), null throws (Edge).
  This exercises every branch of the method body. INV-01..INV-05 each have a covering test.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6573`; FindKnownMotifs subset
  `Failed: 0, Passed: 12`; project builds with **0 errors** (4 pre-existing NUnit2007 warnings are
  in an unrelated file, ApproximateMatcher_EditDistance_Tests.cs — not touched this session).

**Test-quality gate: PASS.**

### Findings / defects

None. Code faithfully realises the validated description; tests are sourced, exact, and complete.
No code or test change was required this session.

## Verdict & follow-ups

- **Stage A: PASS** — description is correct against Rosalind SUBS, Wikipedia EcoRI, and Gusfield
  exact-matching definition; overlapping-occurrence and 0-based conventions are standard and
  internally consistent.
- **Stage B: PASS** — implementation computes the exact-matching set with overlaps via the suffix
  tree and sorts for determinism; all 12 tests pass; full suite green (6573/0); build clean.
- **End-state: CLEAN** — no defect found; algorithm fully functional. No follow-ups.

---

## B09 review 2026-11-18 (campaign 2026-09) — re-validation

- **Stage A: PASS-WITH-NOTES.** The 2026-06 description (exact set matching, overlaps, 0-based, absent motifs omitted) stands. Notes: (a) the unit has no catalogue of "known" motifs — the set is caller-supplied (the fixed TATA/Kozak/CAAT/GC-box library is `MotifFinder.FindRegulatoryElements`, B05); (b) "not implemented: IUPAC motifs, reverse strand" was a doable gap (EMBOSS fuzznuc does both) — implemented; (c) the algorithm doc declared a single engine although the suffix-tree index costs 23 s / 1.7 GB at 10 Mbp for one query (B06 measurement) — engine choice now documented and measured (1 Mbp: tree 924 ms vs KMP 6.9 ms per query; break-even 73-330 queries).
- **Stage B: PASS-WITH-NOTES → FIXED.** Results of `FindKnownMotifs` were correct (sorted, overlapping, deduplicated, empty skipped). Defects: F18 `FindMotif` returned suffix-tree order (Rosalind SUBS `ATAT`: `[1,9,3]`) and duplicated `MotifFinder.FindExactMotif` (null sequence → `NullReferenceException`); F19 no IUPAC / both-strand search; F20 one-off and few-motif searches paid the suffix-tree build.
- **Oracle:** EMBOSS fuzznuc 6.6.0 `-pmismatch 0 -complement Y` equals an independent Python `re` implementation on 2448 (sequence, pattern) pairs (16 exact + 20 IUPAC patterns × 68 sequences); `FindKnownMotifHits` equals both on every pair.
- **Audit r1 (A1-20):** mismatch-tolerant search (fuzznuc `-pmismatch k`) implemented as `FindKnownMotifHits(maxMismatches)` over the canonical `ApproximateMatcher` (additive IUPAC-pattern overload, "B09 for B05"), `KnownMotifHit.Mismatches`, MCP `mismatches` + `matchMismatches`/`reverseMismatches`; 6576-row fuzznuc oracle (B09.md F49).
- **Fixes / tests / docs:** see B09.md F18-F20; `GenomicAnalyzer_KnownMotifs_Oracle_Tests` (+10 tests; +2 in `Seqeron.Mcp.Analysis.Tests`), MCP `find_known_motifs` (`bothStrands`, `iupac`, `reverseMatches`), algorithm doc §2.4/§4/§5/§6, TestSpec INV-6..9 + O1-O7, Evidence (fuzznuc).
- **End-state: CLEAN.**
- **Audit r3 (A3-12, F140):** `fuzzpro` / `fuzztran` ported (`ProteinPattern`, `FindProteinPatternHits`, `FindTranslatedPatternHits`, MCP `fuzzpro` / `fuzztran`): 900 / 900 and 760 / 760 oracle rows equal; the "fuzznuc reads a sequence ambiguity code as a base set" claim corrected (letters compared literally; A3-12a registered with the measured rules).
- **Audit r3 (A3-12a, F146):** fuzznuc with IUPAC codes in the sequence: `FindPatternHits(..., sequenceAmbiguity: true)`, MCP `sequenceAmbiguity`; 1 350 / 1 350 oracle rows (subset rule, classes plain only, exclusions everything unlisted, U = T, reverse = IUPAC revcomp).

- **Audit A4-3 (F153):** the EMBOSS reader drops gap characters and every other non-letter (except `*`) from plain input sequences; one shared `EmbossSequenceReader` (Infrastructure) serves fuzznuc, fuzzpro, fuzztran: 240 oracle rows each equal (`emboss_gap_oracle.tsv`).

- **Audit A4-1 (F157):** `FindTranslatedPatternHits`, `FindProteinPatternHits` poll the ambient cancellation token (mid-run cancel tests); MCP tools take the request token.

- **Audit r7 (A5-8b, F184):** all 24 `-rformat` of the report (simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq; diffseq refused as the program does), `-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq`, several sequences per file and the sequence qualifiers through ONE writer and `EmbossReportPrograms`: 12 000 report files (fuzznuc / fuzzpro / fuzztran / etandem / equicktandem, 100 per format) equal the installed programs byte for byte; MCP `emboss_report_file`.
- **Audit r5 / r6 (A5-7, A5-8, A6-2, A6-4, A6-10; F172, F176, F177, F179, F183):** `fuzznuc` / `fuzzpro` / `fuzztran` read the sequence-input qualifiers (180 oracle rows each; fuzztran prints the program's own, negative reverse-frame coordinates) and write the `-rformat gff` / `excel` reports (`EmbossReports`, `EmbossReportWriter`: 1 200 whole reports equal for the five programs); the qualifier text resolves like the ACD library (`fuzznuc a.fa -pattern ACG -sbegin1 5` = `-sbegin 5`; `-se` is ambiguous in these one-sequence programs); `EmbossReports.Fuzztran` refuses a hit without `EmbossStart` / `EmbossEnd`, `EmbossReportWriter.Write` an undefined format; the MCP report mode runs the search once.
- **Audit A5-8c (F192):** NCBI-style FASTA ids (`ajSeqParseNcbi`), `seqSetName`, `EMBOSS_nnn` names and pattern files (`EmbossPatternList`, multi-pattern reports of fuzznuc / fuzzpro / fuzztran): 1 110 / 1 110 report files equal the installed programs.
- **WP47 (A5-8c2 / A5-8c3, F195 / F196):** `fuzzpro` / `fuzztran` `-pmismatch` (`FindProteinPatternHits(..., maxMismatches)`: Baeza-Yates-Perleberg = Hamming scan, `EmbossTubSearch` = `embPatTUBInit` / `embPatTUBSearch` (skip table of the letters A-Z only), the recursive matcher `patOUBrute` with one hit per matching unrolling; 750 fuzzpro + 250 fuzztran oracle rows equal, `fuzzpro_mismatch_oracle.tsv`); the same duplicate-hit order for `fuzznuc` with a range pattern and mismatches (`fuzznuc_dup_oracle.jsonl.gz`, 300 reports equal); `diffseq` (`EmbossDiffseqPrograms`, 336 pairs x all 24 report formats and both GFF3 feature files equal).
- **WP51 (A5-8c3b / A5-8c3c, F216 - F220):** `fuzznuc` / `fuzzpro` / `fuzztran` patterns with every element an optional range list the EMPTY hits (`end = start - 1`; the empty match at the end of an end-anchored pattern; the translation of an empty fuzztran hit): 400 / 400 whole reports equal; `diffseq` takes the sequence qualifiers (260 / 260) and GenBank / EMBL inputs with their feature tables (420 / 420: the features that overlap a difference are listed with it); every tag with a special validator of the library is checked as the program does (450 / 450 runs: bad values are notes, correctable ones are corrected, `rpt_unit` / `rpt_unit_seq` share one compiled pattern per run), and both feature files are written in every output feature format of the program (`-offormat4` / `-offormat5`: 500 / 500 runs, nucleotide and protein pairs); details in `docs/algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md` and B09.md F216 - F220.


- **Audit R14-3 (F290):** `FindMotif` / `FindExactMotif` do not build the suffix tree for a one-off query (`DnaSequence.HasSuffixTree`, `ReleaseSuffixTree`; 2 Mbp 4.37 s / 542 MB -> 66 ms); the tree route stays for a sequence whose tree is already built.

- **R14-18 (F360, WP78):** the mismatch-tolerant `FindKnownMotifHits` runs on the bit-parallel Hamming scan of `ApproximateMatcher` (2 Mbp, 50 motifs x 10 bp, both strands: mm = 2 4.6 s -> 1.3 s, IUPAC mm = 2 8.2 s -> 0.95 s), an exact A/C/G/T motif on a `DnaSequence` too; results identical (`ApproximateMatcher_PackedHamming_Tests`).
