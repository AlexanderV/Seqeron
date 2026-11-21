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
**Accessed:** 2026-06-13 · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4 (Wikipedia, used for a widely-documented sequence fact also held in REBASE/NEB rank-5 DBs)
**Retrieved by:** WebSearch "REBASE EcoRI recognition sequence GAATTC restriction enzyme" → WebFetch of the article.

**Key Extracted Points:**

1. **EcoRI recognition site (verbatim):** "the common restriction enzyme EcoRI recognizes the palindromic sequence GAATTC and cuts between the G and the A on both the top and bottom strands." Provides a real biological motif (`GAATTC`, length 6) for the test dataset.

### EMBOSS `fuzznuc` 6.6.0 (reference tool, run locally) — IUPAC patterns and both strands

**Source:** Debian `emboss` 6.6.0 binary (`fuzznuc -sequence s.fa -pattern P -pmismatch 0 -complement Y -rformat gff`), Rice, Longden & Bleasby (2000) Trends Genet 16:276. Accessed 2026-11-18 (B09 review).
**Authority rank:** 3 (reference implementation).

**Key Extracted Points (observed, then confirmed against an independent Python `re` lookahead implementation on 2448 sequence/pattern pairs, 816 with forward and 776 with reverse hits, `scripts/analysis/generate_known_motif_oracle.py`):**

1. Overlapping occurrences are all reported (AAA in AAAAA → 3 hits).
2. IUPAC codes in the pattern match their NC-IUB base sets (`TATAWAW`).
3. With `-complement` the reverse strand is searched: a '-' hit is a forward-strand window equal to the reverse complement of the pattern, reported with forward-strand coordinates; a palindromic pattern (GAATTC) is reported on both strands at the same coordinates (`GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC`: `+` 1-6, 10-15 and `-` 1-6, 10-15; TATAWAW `+` 18-24, `-` 17-23, 1-based).
4. Reference for the cheap-engine claim: Rosalind SUBS sample (`GATATATGCATATACTT`, `ATAT` → 2 4 10, 1-based) = `[1,3,9]` 0-based; the Seqeron suffix tree enumerates it 1, 9, 3 (B09 F18).
5. **Mismatches (A1-20, `-pmismatch k`, run 2026-11-18 on the same binary).** A hit is a window with at most k disagreeing positions (IUPAC pattern positions agree when the base is in the code's set); each hit line carries the note `*mismatch n` (absent = 0); with `-complement Y` the reverse-complement pattern is searched with the same bound and the report line is a forward-strand window (`AACCGT`, k = 1 on `AACCGTACGGTA`: `+` 1-6 exact, `-` 7-12 `*mismatch 1`; `GAATTC` k = 1 on `GAATTCAAGATTCGGATTTTGAATTTT`: `+`/`-` 1-6 exact and 21-26 `*mismatch 1`; `RAATTC` k = 2: `+` 0:0 7:1 13:2 14:2 20:1 21:2, `-` 0:0 7:2 13:1 14:1 20:0 21:2, 0-based). Equal to an independent Python Hamming scan (positions and counts) on 6576 (sequence, pattern, k) rows (20 patterns, k = 1..3, sequences incl. planted 0-3-substitution copies on both strands; `generate_known_motif_oracle.py --mismatches`, byte-reproducible). Two fuzznuc quirks, both outside the oracle: a hit overhanging the end of a sequence shorter than the pattern is reported (`AAA`, `-pmismatch 2`, sequence `AC` → 1 hit; 28 of 3330 probe cases, all with sequence shorter than pattern); k at or above the pattern length is refused for some patterns (`embPatCompileII: Mismatches (4) must be less than the real pattern length`, `RYRY`).

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

## 2026-11-19 audit round 2 (A2-32 / A2-35, B09 F96 / F99): the fuzznuc pattern language and string overloads

- **Source definitions.** EMBOSS 6.6.0 `fuzznuc -help` (the installed binary: the `-pattern` text quoted in `FuzznucPattern`), and the original program source `nucleus/embpat.c`, `nucleus/embpatlist.c`, `emboss/fuzznuc.c` read on `raw.githubusercontent.com/kimrutherford/EMBOSS` (User-Agent `Seqeron-review/1.0`): `embPatClassify` (grammar and errors), `embPatGetTypeII` (search class by mismatch / range / class flags), `patOUBrute` / `embPatBruteForce` (recursive matcher), `embPatMatchFindAllC` (regex all-matches, `nterm` break), `embPatternSeqSearch` (reverse strand = reverse-complemented sequence, same pattern).
- **Findings decided on the binary.** (1) Fixed-length patterns (classes, exclusions, repeats, anchors, any k < pattern length) report every window within k substitutions — 2000+ rows agree with a plain Hamming scan; anchors flip on the reverse strand. (2) A range reports every matching window length (k = 0: regex DFA all-matches). (3) With -pmismatch > 0 a range pattern is searched by `patOUBrute`: a mismatch is allowed only in the fixed elements and the first n repeats of `X(n,m)`; the extra repeats must match exactly, and a path that still has pattern elements when the scanned sequence ends fails (a trailing optional element chosen zero times cannot end on the last base) — 960 range-with-mismatch rows. (4) An exact range pattern starting with `<` returns only the longest match at the anchor (`nterm`) — on the reverse strand only the longest window ending at the last base. (5) The binary accepts malformed spellings in unspecified ways (`A<C` = `AC`, `{}` = any base, `A(0,0)` empty hits, `A.C`, `[A-C]`) and rejects ambiguity letters inside brackets; `-pattern @file` takes one pattern per line or FASTA `>name` + pattern (mixed files are not named consistently). (6) A sequence ambiguity code matches when the pattern symbol covers it (`R` vs `N` hits, vs `A` or `[AG]` does not); not reproduced — the canonical B05 contract (the sequence side is literal) is kept and documented.
- **Oracle numbers.** `scripts/analysis/generate_fuzznuc_pattern_oracle.py`: 3300 random syntax-tree (sequence, pattern, k) rows (1758 with ranges, 960 range + mismatches, 987 anchored, 2298 with reverse hits; k < the shortest hit), fuzznuc 6.6.0 `-pmismatch k -complement Y -rformat gff` == the independent Python reference with the two search classes above, 0 disagreements.
- **Implementation.** `FuzznucPattern` (parser, `Expand` / `ExpandWithExtras`, `ReverseComplement`, `ParsePatternFile`) and `GenomicAnalyzer.FindPatternHits` (+ named-pattern overload) compile to the canonical `ApproximateMatcher.FindWithMismatches(…, iupacPattern: true)` over each fixed-length expansion; no scanning loop was added. String overloads of `FindMotif`, `FindKnownMotifs`, `FindKnownMotifHits`: the sequence's non-A/C/G/T symbols are replaced by a character that equals nothing.
- **Tests.** `GenomicAnalyzer_PatternHits_Oracle_Tests` (14), `GenomicAnalyzer_StringOverloads_Tests` (3); MCP new `find_pattern_hits`, IUPAC input of `find_motif` / `find_known_motifs`.

## 2026-11 audit round 3 (A3-12, B09 F140): fuzzpro, fuzztran and the sequence-ambiguity rules of fuzznuc

Sources: `emboss/fuzzpro.c`, `emboss/fuzztran.c` (`fuzztran_SourceFeature`), `emboss/fuzznuc.c`, `nucleus/embpat.c` (EMBOSS 6.6.0, raw.githubusercontent.com/kimrutherford/EMBOSS) and the installed binaries (`fuzzpro -help`: only `-pattern`, no `-pmismatch`; `fuzztran -help`: `-frame`, `-table`, no `-pmismatch`).

**fuzzpro.** `scripts/analysis/generate_fuzzpro_oracle.py`: 900 random syntax-tree rows (protein 12-140 aa, with X, B, Z, U, J, O, * and lower case; letters, `x`, classes, exclusions, repeats, ranges, `-`, `<` / `>`, `.`; 57 rows the binary REJECTS - B / Z / X / J nested in `[]` / `{}`; 275 with hits; 456 with ranges, 263 anchored) = `FindProteinPatternHits` on **900 / 900**. Findings of the measurement (an independent Python reference first failed on 13 rows): (1) `B` = D or N, `Z` = E or Q, `J` = I or L as pattern letters and never the letters themselves; `O` and `U` are literal; `x` / `X` any; (2) a start-anchored range pattern returns only the longest match (as the nucleotide regex class); (3) a single trailing `*` of the sequence is read off (`ACDEF*` and `ACDEF` give the same hits, `ACDEF**` keeps one `*`); (4) `[]`, `{}` (matches everything), `A<C` / `A>C` and a repeat `(1,0)` are accepted or fail in unspecified ways by the binary and are rejected by the parser.

**fuzztran.** `scripts/analysis/generate_fuzztran_oracle.py`: 760 random rows (DNA 24-300 nt with N and lower case; patterns built around windows of a translated frame so that hits occur; frame in `1 2 3 F -1 -2 -3 R 6`; table 0, 1, 2, 3, 4, 5, 9, 11; 384 rows with hits, 57 rejected patterns): frame, protein start / end and translation of every hit equal `FindTranslatedPatternHits` on **760 / 760**, the printed nucleotide coordinates equal on every forward-frame hit. Findings: (1) fuzztran searches only complete codons - the dangling-base residue that transeq appends (`ATW`... at the end of frames 2, 3, -2, -3) is not searched; (2) a stop `*` of the translation is not read off (unlike fuzzpro's reader); (3) for the reverse frames the coordinates the binary prints are negative numbers (`start - frame - 3*S + 4`, `fuzztran_SourceFeature`: the program passes its range BEGIN as `start`, and the range END for frame -3 of `-frame R` / `6`; WP42 A5-7 found the rule: they are deterministic, not state-dependent) - since WP42 `EmbossStart` / `EmbossEnd` carry them (the earlier version left them null and reported only the forward-strand location `Position`, `Length`, which stay); (4) ambiguous codons: EMBOSS resolves them with its own table (`ajax/core/ajtranslate.c` `trnComplete`), the library with the Biopython-style `GeneticCode.Translate`; on all 3375 three-letter IUPAC codons they differ for 13 patterns (`ATW`, `TRA`, `YTA`, `YTG`, `YTR`, `MGA`, `MGG`, `MGR` = EMBOSS X where the codon is unique; `ATS`, `TRR`, `YTN`, `WSN`, `MGN` = EMBOSS resolves where the set of amino acids is not unique) - the existing B02 limitation, so the oracle uses only N as ambiguity code (the tool description says so).

**fuzznuc and a sequence ambiguity code (registered A3-12a).** The earlier statement "fuzznuc treats a sequence ambiguity code as a base set" is refuted. Measured on fuzznuc 6.6.0 (`-pmismatch 0`): sequence `ACGRTT`: pattern `ACGA` no, `ACGG` no, `ACGR` yes, `ACGN` yes, `ACGM` no, `ACG[AG]` no, `ACG{C}` yes, `ACG{AG}T` yes (sequence R is not in the exclusion list), `ACGAT` no; sequence `ACGNTT`: `ACGN` yes, `ACGA` no; `ACGBTT` / `ACGN` yes. Rule: letters are compared literally, a pattern IUPAC letter stands for its bases AND itself, `N` for every symbol, an exclusion for every symbol not listed, a class for the listed bases only. Porting needs an element kind (letter / class / exclusion) in `PatternElement` and a scan over letter sets (the expansion to IUPAC strings loses the distinction between `R` and `[AG]`); registered with these rules and a 1 000-row oracle plan.

## 2026-11 audit round 3 (A3-12a, B09 F146): fuzznuc with an IUPAC code in the sequence implemented

Rule measured and refined on fuzznuc 6.6.0 (EMBOSS `embpat.c`, which the earlier notes read, uses a bit-set per symbol): a pattern IUPAC letter L matches a sequence symbol S iff the base set of S lies within the base set of L (pattern `R` matches A, G, R; sequence `R` matches the letters R, V, D, N, not A, M, Y; the first note's "L stands for its bases AND itself" is the single-base special case); pattern `N` matches every symbol including a non-IUPAC `X`; a class matches the listed plain bases only; an exclusion matches every symbol not listed (so `{ACGT}` matches `R`, `N`, `X`); `U` reads as `T`; a `-` gap is removed by the sequence reader (outside the contract); the reverse strand is the IUPAC reverse complement of the sequence searched with the same pattern (R<->Y, K<->M, B<->V, D<->H). Everything else is the A2-32 behaviour (mismatch counting, repeats, ranges, anchors, the two search classes).

Oracle `fuzznuc_iupac_sequence_oracle.tsv` (`generate_fuzznuc_iupac_sequence_oracle.py`): random syntax-tree patterns x random sequences over A C G T plus the 11 ambiguity letters, U and X (14-40 symbols, planted instances on both strands), `-pmismatch 0-2 -complement Y`: **1 350 rows, fuzznuc == the independent Python reference, 0 disagreements** (an earlier hypothesis "pattern letter matches itself and its bases" failed 59 / 450, the subset rule 0 / 1 350).

Implementation: `PatternElementKind` (letter / class / exclusion) on `PatternElement`; `FuzznucPattern.ExpandPaths` (element per position; `ExpandWithExtras` is now derived from it: the plain IUPAC expansion is unchanged, 3 300 / 3 300 original rows still equal); `FuzznucPattern.Matches(element, symbol)`; `GenomicAnalyzer.FindPatternHits(string ..., sequenceAmbiguity: true)` scans the paths with the literal rule; MCP `find_pattern_hits` parameter `sequenceAmbiguity`. Tests: `GenomicAnalyzer_PatternHits_IupacSequence_Tests` (5), MCP `StrandPatternAndRepeatToolsTests` (+1).

## Audit A4-3 (B09 F153): EMBOSS reader - gap characters and other non-letters

Sources actually opened: EMBOSS `ajax/core/ajseqread.c` (`seqAppend`: `ajStrKeepSetAlphaC(.., "*.~?#+-")`), `ajseqtype.c` (`ajSeqTypeCheckIn`, `seqTypeFix`: an ungapped type calls `ajStrRemoveGap`; `seqType[]` table: every ACD type of the programs is ungapped, `gapany` for consambig / showalign), `ajstr.c` (`ajStrRemoveGap` keeps `A-Za-z*`), `ajseq.c` (`ajSeqIsNuc`, `ajSeqGapStandard`) from raw.githubusercontent.com/kimrutherford/EMBOSS (User-Agent `Seqeron-review/1.0`), the installed ACD files (`/usr/share/EMBOSS/acd`: `type:` of every sequence input) and the installed programs. Measured with `scripts/alignment/measure_emboss_reader_chars.py` (each of the characters 1..255 except newline / carriage return / `>` inserted into the middle of a sequence, both sequences, nucleotide and protein): for fuzznuc, fuzzpro, fuzztran the output equals the output without the character for exactly 199 characters (all except the 52 ASCII letters and `*`); `*` is kept (a protein stop; an error in a nucleotide-typed input).

Oracle `emboss_gap_oracle.tsv` (`generate_emboss_gap_oracle.py`, seed 20261126): 18 programs x 240 rows, fuzznuc (240), fuzzpro (240) and fuzztran (240; hit coordinates and the sequence length of the GFF header count the cleaned sequence): **equal on every row**. Audit example: `needle AC-GTTAGC.TT* ACGTTAGCTT` prints `# Length: 11`, `# Score: 58.0` (the port scored 43 on the uncleaned text).

## Audit A4-1 (B09 F157): cancellation

Mechanism: `AlignmentCancellation` (SequenceAligner.Cancellation.cs, A2-40); tests `SequenceAligner_Cancellation_Tests` (`PreCancelled_A4_1_...`, one `CancelledWhileRunning_...` test per class) and the MCP `CancellationTokenTests` / `EmbossProgramToolsTests`: `FindTranslatedPatternHits`, `FindProteinPatternHits` poll the ambient token (a mid-run cancel stops with `OperationCanceledException` inside the named poll) and the MCP tool takes the request token.

**A5-7 / A5-8 (WP42).** EMBOSS 6.6 `ajacd.c`, `ajseq.c`, `ajseqread.c`, `ajseqtype.c`, `ajmath.c`, `ajalign.c`, `fuzznuc.c`, `fuzzpro.c`, `fuzztran.c`, `seqmatchall.c` opened (raw.githubusercontent.com/kimrutherford/EMBOSS, User-Agent `Seqeron-review/1.0`); installed 6.6.0: `emboss_seqopts_oracle.tsv` (180 rows per program) and `emboss_report_oracle.tsv` (120 runs per program and format, gff + excel). fuzztran reverse-frame numbers `S - |frame| - 3 start + 4` (frame -3 of R / 6: `S` = the end) reproduced.

**A6-2 / A6-4 / A6-10 (WP43).** The qualifier text of the fuzznuc / fuzzpro / fuzztran tools is resolved as the ACD library resolves it (`fuzznuc a.fa -pattern ACG -sbegin1 5` = `-sbegin 5`, `-se` ambiguous in the one-sequence programs: installed 6.6.0 and `ajacd.c`; docs/Evidence/ALIGN-MULTI-001-Evidence.md A6-2); `EmbossReports.Fuzztran` refuses a hit without the program's coordinates.
