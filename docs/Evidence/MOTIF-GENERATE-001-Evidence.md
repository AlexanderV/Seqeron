# Evidence Artifact: MOTIF-GENERATE-001

**Test Unit ID:** MOTIF-GENERATE-001
**Algorithm:** IUPAC-Degenerate Consensus Generation (`MotifFinder.GenerateConsensus`)
**Date Collected:** 2026-06-14

---

## Online Sources

### Cornish-Bowden / NC-IUB (1984) — Nomenclature for incompletely specified bases in nucleic acid sequences

**URL:** https://academic.oup.com/nar/article/13/9/3021/2381659 (DOI 10.1093/nar/13.9.3021)
**Accessed:** 2026-06-14
**Authority rank:** 2 (official IUPAC/IUB nomenclature standard)
**Retrieved how:** WebSearch "Cornish-Bowden 1985 Nomenclature for incompletely specified bases…" located the NAR article; WebFetch of the Oxford Academic article page returned the symbol table.

**Key Extracted Points:**

1. **Single-letter degenerate symbols (verbatim mapping from the fetched table):** R = A or G (purine); Y = C or T/U (pyrimidine); M = A or C (amino); K = G or T/U (keto); S = G or C (strong); W = A or T/U (weak); B = C, G, or T/U (not-A); D = A, G, or T/U (not-C); H = A, C, or T/U (not-G); V = A, C, or G (not-T/U); N = any base (A, C, G, T/U).
2. **Meaning rule:** "each symbol stands for the specific set of bases listed" — a degenerate symbol denotes exactly the set of standard bases it abbreviates. The set→symbol mapping is bijective over the 11 non-trivial subsets.

### UCSC Genome Browser — IUPAC Nucleotide Code Table

**URL:** https://genome.ucsc.edu/goldenPath/help/iupac.html
**Accessed:** 2026-06-14
**Authority rank:** 5 (well-maintained genomics database documentation; corroborates the standard)
**Retrieved how:** WebSearch "IUPAC nucleotide ambiguity codes…"; WebFetch of the UCSC page returned the full table.

**Key Extracted Points:**

1. **Full table (verbatim):** G=G, A=A, T=T, C=C; R=G or A; Y=T or C; M=A or C; K=G or T; S=G or C; W=A or T; H=A or C or T (not-G); B=G or T or C (not-A); V=G or C or A (not-T); D=G or A or T (not-C); N=G or A or T or C (any).
2. **Use context:** single-character codes represent multiple observed alleles at a single position — i.e., the same set-of-bases semantics used when summarising an aligned column.

### Wikipedia — Nucleic acid notation (citing NC-IUB 1984)

**URL:** https://en.wikipedia.org/wiki/Nucleic_acid_notation
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing the NC-IUB 1984 primary)
**Retrieved how:** WebFetch of the article; Table 1 returned with its primary reference.

**Key Extracted Points:**

1. **Table 1 mapping** matches Cornish-Bowden/UCSC exactly (W,S,M,K,R,Y two-base; B,D,H,V three-base; N four-base).
2. **Primary reference cited:** "Nomenclature Committee of the International Union of Biochemistry (NC-IUB) (1984), Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences, Nucleic Acids Research" — confirming the table is the IUPAC/NC-IUB standard, not Wikipedia-original.

### Bioconductor DECIPHER — `ConsensusSequence` (reference implementation of threshold-based degenerate consensus)

**URL:** https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (established library reference implementation)
**Retrieved how:** WebSearch "degenerate consensus sequence IUPAC threshold…"; WebFetch returned the man page text.

**Key Extracted Points:**

1. **Threshold-consensus mechanism (verbatim):** "ConsensusSequence removes the least frequent characters at each position, so long as they represent less than `threshold` fraction of the sequences in total"; remaining characters are then "represented using IUPAC degeneracy codes."
2. **Tie / equal-abundance rule (verbatim):** "Degeneracy codes are always used in cases where multiple characters are equally abundant."
3. **Establishes the family:** a degenerate consensus is parameterised by a frequency threshold; bases above the threshold at a column are combined into the IUPAC symbol for that base set. (DECIPHER's default threshold is 0.05; the threshold value is a tunable parameter, not a fixed universal constant.)

---

## Documented Corner Cases and Failure Modes

### From DECIPHER `ConsensusSequence`

1. **Frequency threshold governs inclusion:** a minority base whose frequency is below the threshold is dropped and not encoded in the ambiguity code; the threshold value is a design parameter chosen per tool.
2. **Equal abundance → degeneracy code:** when ≥2 bases pass the inclusion rule, the column emits the IUPAC symbol for that set rather than picking one base arbitrarily.

### From NC-IUB 1984 / UCSC

1. **N is the four-base symbol:** N denotes the full set {A,C,G,T}; any single missing base yields a three-base not-X symbol (B/D/H/V) instead of N.

---

## Test Datasets

### Dataset: IUPAC set→symbol mapping (NC-IUB 1984 / UCSC)

**Source:** Cornish-Bowden NAR 13(9):3021 (DOI 10.1093/nar/13.9.3021); UCSC IUPAC table.

| Base set present at column | IUPAC symbol |
|----------------------------|--------------|
| {A} | A |
| {C} | C |
| {G} | G |
| {T} | T |
| {A,G} | R |
| {C,T} | Y |
| {C,G} | S |
| {A,T} | W |
| {G,T} | K |
| {A,C} | M |
| {C,G,T} | B |
| {A,G,T} | D |
| {A,C,T} | H |
| {A,C,G} | V |
| {A,C,G,T} | N |

### Dataset: Threshold-inclusion worked examples (this implementation, threshold = n × 0.25, base included iff count > threshold)

**Source:** derivation from the DECIPHER threshold-consensus family applied to this implementation's documented 25 % design constant; the resulting *symbol* per base set is dictated by the NC-IUB table above.

| Input column (n seqs) | n | threshold = n×0.25 | bases with count > threshold | symbol |
|-----------------------|---|--------------------|------------------------------|--------|
| A,G | 2 | 0.5 | {A,G} (each 1>0.5) | R |
| C,T | 2 | 0.5 | {C,T} | Y |
| C,G | 2 | 0.5 | {C,G} | S |
| A,T | 2 | 0.5 | {A,T} | W |
| G,T | 2 | 0.5 | {G,T} | K |
| A,C | 2 | 0.5 | {A,C} | M |
| C,G,T | 3 | 0.75 | {C,G,T} (each 1>0.75) | B |
| A,G,T | 3 | 0.75 | {A,G,T} | D |
| A,C,T | 3 | 0.75 | {A,C,T} | H |
| A,C,G | 3 | 0.75 | {A,C,G} | V |
| A,A,G,G,C | 5 | 1.25 | {A(2),G(2)}; C(1) dropped | R |
| A,A,A,G | 4 | 1.0 | {A(3)}; G(1) at ≤threshold dropped | A |
| A,C,G,T | 4 | 1.0 | none (each 1, not >1.0) → all four tied at max | N (F13, 2026-09; was A) |
| A,C,-,- | 4 | 1.0 | none → A,C tied at max | M (F13) |

---

## Assumptions

1. **ASSUMPTION: 25 % inclusion threshold is a documented design constant.** This implementation includes a base in a column's IUPAC code iff its count is strictly greater than 25 % of the number of sequences (`count > total × 0.25`). The *threshold-consensus family* and the "include bases above a frequency threshold, encode the set with an IUPAC code" rule are authoritative (DECIPHER); the specific 25 % cut and the strict `>` boundary are this implementation's design choice (DECIPHER's own default is 0.05 and tools vary). It is correctness-affecting but documented and named (`threshold = total * 0.25`), not invented-untraceable. Tests pin the boundary behaviour explicitly and otherwise use inputs where the inclusion decision is unambiguous so the verified *symbol* is dictated solely by the authoritative NC-IUB table.
2. **Fallback when no base passes the threshold (sourced since 2026-09-29, B05 F13a).** When no base exceeds the threshold, the column is the IUPAC code of **all bases tied at the maximum count** (DECIPHER `ConsensusSequence`: "degeneracy codes are always used in cases where multiple characters are equally abundant"; Biopython `degenerate_consensus` agrees): four equal bases ACGT → `N`, `A,C,-,-` → `M`; a column without any A/C/G/T → `N`. (Previously the single first most-frequent base was emitted, e.g. ACGT → `A` — superseded.)
3. **Input contract (since 2026-09-29, B05 F13b): all rows must have equal length; case-insensitive over {A,C,G,T}; non-ACGT characters at a position are ignored in the counts.** Ragged rows (and null rows) throw `ArgumentException` ("All sequences must have the same length."), as Biopython `MultipleSeqAlignment` / `motifs.create`; the column length is no longer silently taken from the first sequence. Inputs are upper-cased; only A/C/G/T are counted per the four-base alphabet.

---

## Recommendations for Test Coverage

1. **MUST Test:** every two-base set maps to the correct IUPAC code (R,Y,S,W,K,M) and every three-base set to (B,D,H,V) — Evidence: NC-IUB 1984 / UCSC table.
2. **MUST Test:** unanimous columns reproduce the input base (A/C/G/T) — Evidence: NC-IUB (singleton set → standard base).
3. **MUST Test:** strict-`>` 25 % boundary — a base at exactly 25 % is excluded; a base above 25 % is included — Evidence: implementation design constant (documented).
4. **MUST Test:** fallback to most-frequent base when no base passes the threshold (four bases each at 25 %) — Evidence: implementation contract.
5. **SHOULD Test:** minority base below threshold dropped (A,A,G,G,C → R, not a three-base code) — Rationale: confirms threshold filtering precedes IUPAC encoding.
6. **SHOULD Test:** case-insensitivity (lowercase input) and empty-collection → "" — Rationale: documented input normalisation / guard.
7. **COULD Test:** null input throws `ArgumentNullException` — Rationale: documented guard.

---

## References

1. Cornish-Bowden A. (1985). Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. Nucleic Acids Research 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021
2. UCSC Genome Browser. IUPAC ambiguity codes. https://genome.ucsc.edu/goldenPath/help/iupac.html (accessed 2026-06-14)
3. Wikipedia. Nucleic acid notation (Table 1, citing NC-IUB 1984). https://en.wikipedia.org/wiki/Nucleic_acid_notation (accessed 2026-06-14)
4. Wright E.S. DECIPHER `ConsensusSequence` (Bioconductor). https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-10-01**: Assumptions 2/3 synced to B05 F13 (tied-max IUPAC fallback, ACGT → N; ragged/null rows → `ArgumentException`).

---

## Review 2026-09 (B05, F13) — sources actually opened

- **Biopython 1.88** installed source `Bio/motifs/matrix.py` `GenericPositionMatrix.degenerate_consensus` — Cavener rules verbatim: `counts[0] > sum(counts[1:]) and counts[0] > 2*counts[1]` → single; `4*sum(counts[:2]) > 3*sum(counts)` → pair; `counts[3] == 0` → triple; else N ("The same rules are used by TRANSFAC").
- **Biopython Tutorial** `Doc/Tutorial/chapter_motifs.rst` (raw.githubusercontent.com): `m.degenerate_consensus` = `WACVC`; reverse complement `GBGTW`; slice `m[2:-1]` = `CV`, "constructed following the rules specified by Cavener".
- **Cavener 1987** NAR 15(4):1353 (PMID 3822832) — WebSearch snippet only: single base if frequency > 50 % and > twice the second; two bases if their sum > 75 %.
- **DECIPHER `ConsensusSequence`** — WebSearch snippet (rdrr.io/bioc manuals blocked for curl) confirms verbatim "Degeneracy codes are always used in cases where multiple characters are equally abundant." DECIPHER's `threshold` is cumulative (least-frequent characters removed while together < threshold, default 0.05) — the library's per-base 25 % cut is therefore a design constant, not DECIPHER's rule.
- Biopython `MultipleSeqAlignment` rejects unequal rows ("Sequences must all be the same length"); `motifs.create` likewise.

### Reference numbers (Biopython 1.88 `degenerate_consensus`)

| Rows | Result |
|---|---|
| A,C,G,T | N |
| A,A,A,C | A |
| A,A,G,G | R |
| A,C,G | V |
| A,A,C,G,T | N |
| A×5 C×3 G T | M |
| A×6 C×2 G T | A |
| AAAA, AAGT, AACT, AATT | AANT |
| tutorial 7 rows | WACVC |

Plus 12 random alignments (seed 20260930) locked in `MotifFinder_GenerateConsensus_Tests.CavenerBiopythonCases`.

## Configurable inclusion threshold (B05 follow-up, 2026-09-30)

`GenerateConsensus(sequences, θ)` generalises the 25 % design constant (declared limitation) to θ ∈ [0, 1]; θ = 0.25 shares the code path of the parameterless overload (bit-identical; property test on 500 random alignments). The rule (count > θ·n, tie fallback to the maximum-count bases, no A/C/G/T → N, NC-IUB map) was cross-checked against an independent Python implementation on 700 random alignments (1–12 rows, 1–40 columns, alphabets ACGT/ACGTN with gaps and lower case; θ ∈ {0, 0.1, 0.25, 0.3, 1/3, 0.5, 0.75, 1, random}) → 700/700 identical; 20 locked in `MotifFinder_AlignmentConsensus_Tests.ThresholdCases`. Not a published rule (DECIPHER's threshold is cumulative; Cavener = `GenerateCavenerConsensus`).

## DECIPHER `ConsensusSequence` (B05 audit group C, F28, 2026-10-01)

Opened (raw.githubusercontent.com/bioc/DECIPHER/devel — Bioconductor git mirror; bioconductor.org and codeload were proxy-blocked): `DESCRIPTION` (Version 3.9.4), `R/ConsensusSequence.R` (argument checks: threshold ∈ [0, 1), minInformation ∈ (0, 1], noConsensusChar ∈ DNA_/RNA_/AA_ALPHABET; `includeNonLetters` passed to the C argument `ignoreNonLetters`; RNA `T` → `U`; `?` → noConsensusChar), `src/ConsensusSequence.c` (`frontTerminalGaps`/`endTerminalGaps`(`AA`), `alphabetFrequency`(`AA`), `makeConsensus`, `makeConsensusAA`, `consensusSequence`(`AA`) — threshold passed as 1 − threshold), `man/ConsensusSequence.Rd` (examples with stated outputs), `src/DECIPHER.h`, `src/Biostrings_stubs.c`.

Reference build: the verbatim `ConsensusSequence.R` + `ConsensusSequence.c` compiled as an R package (`R CMD INSTALL`, package name DECIPHER, `useDynLib(DECIPHER)`) against Ubuntu noble `r-base-core` 4.3.3 and `r-bioc-biostrings` 2.70.2 (the full DECIPHER package is not packaged for Ubuntu and Bioconductor is unreachable; ConsensusSequence depends only on these two files). All `.Rd` examples reproduce their stated outputs (`W`, `A`, `+`, `N`; `W`/`W`/`X`/`X`/`J`/`J`; `ANSCT-`/`+NSCT-`; `ABZJX-`; `ANNNA`/`AAAAA`; `AWNDA`/`AAAAA`) plus `SWD`/`GTD`/`++D` for the majority example.

| Check | Cases | Result |
|---|---|---|
| `GenerateDecipherConsensus` vs the R/C build, random (seeds 20261001 ×4000, 7 ×6000): DNA/RNA/AA, pure/IUPAC/`BZJX*UO` cores, 0–60 % `-`/`.`/`+`, terminal gaps, 20 % ragged rows, lower case, 1–41 rows, threshold ∈ {0, .05, .1, .25, .3, .5, .75, .9, .99, random}, minInformation default or random, random noConsensusChar from the alphabet, random ambiguity/includeNonLetters/includeTerminalGaps | 10,000 | 10,000/10,000 identical |
| Locked in `MotifFinder_DecipherConsensus_Tests` | 90 stratified + 22 manual + hand-derived | — |

Behaviour notes (from the source, confirmed by the build): the tests run in source order (singletons, `Y K W S R M`, `B D H V`, `N`), each requiring every included base strictly above every excluded one, so A 0.9 / C 0.06 / G 0.04 gives `M`; protein `B`/`Z`/`J` require N/Q/I to be the (possibly D/E/L-tied) maximum; a column with nothing counted is `-` without terminal gaps and `noConsensusChar` with them; with `includeNonLetters = TRUE` gaps and masks are dropped (the R argument is passed to C `ignoreNonLetters`, consistent with the manual's examples although its prose says the opposite); with `ambiguity = FALSE` DNA gaps/masks are always counted.
