# Evidence Artifact: GENOMIC-ORF-001

**Test Unit ID:** GENOMIC-ORF-001
**Algorithm:** Open Reading Frame (ORF) Detection — `GenomicAnalyzer.FindOpenReadingFrames`
**Date Collected:** 2026-06-14

---

## Online Sources

### Rosalind — "Open Reading Frames" (ORF) problem

**URL:** https://rosalind.info/problems/orf/
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 4 (educational reference with an exact worked dataset; used here for the worked example, cross-checked against ranks 1–2)

**Key Extracted Points:**

1. **ORF definition:** A sequence in DNA/RNA potentially able to encode a protein that begins with a start codon and ends with a stop codon, containing no internal stop codons.
2. **Reading frames:** A DNA string generates six total reading frames — three from the forward strand and three from its reverse complement.
3. **Start codon:** AUG (codes for methionine, M), indicating the beginning of translation.
4. **Stop codons:** one of three RNA codons indicating termination (UAA/UAG/UGA → DNA TAA/TAG/TGA).
5. **Return value:** all **distinct** candidate protein strings derived by translating ORFs into amino acids until a stop codon is encountered. (Nested/overlapping ORFs that share a stop are both reported — see Test Dataset, e.g. `MGMTPRLGLESLLE` and `MTPRLGLESLLE`.)

### Wikipedia — "Open reading frame" (using its cited statements)

**URL:** https://en.wikipedia.org/wiki/Open_reading_frame
**Accessed:** 2026-06-14 (retrieved with WebFetch) · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4

**Key Extracted Points:**

1. **Definition:** "Reading frames are defined as spans of DNA sequence between the start and stop codons." An ORF may contain a start codon (usually AUG) and by definition cannot extend beyond a stop codon (usually UAA, UAG or UGA in RNA).
2. **Six-frame translation:** "a DNA strand has three distinct reading frames"; with double-stranded DNA "there are six possible frame translations" — three reading frames on each complementary strand.
3. **Minimum length:** "some authors say that an ORF should have a minimal length, e.g. 100 codons or 150 codons." Length divisible by three.

### NCBI ORFfinder (tool reference)

**URL:** https://www.ncbi.nlm.nih.gov/orffinder/
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 5 (NCBI tool)

**Key Extracted Points:**

1. **Start-codon options:** "ATG" only; "ATG" and alternative initiation codons; Any sense codon ("find all stop-to-stop ORFs").
2. **Minimal ORF length:** selectable nucleotide values 30 / 75 / 150 / 300 / 600; the tool restricts the search to ORFs with length **equal or more than** the selected value (nucleotide units).

### NCBI C++ Toolkit — `COrf::FindOrfs` (the ORFfinder engine) — opened 2026-11-18 (B09 review)

**URL:** https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/sequence/orf.cpp (and `include/algo/sequence/orf.hpp`)
**Authority rank:** 1 (tool source). ORFfinder's own parameters `-ml`, `-g`, `-s`, `-n` map onto it (ProteinIQ's wrapper documentation of the native ORFfinder 0.4.3 executable, https://proteiniq.io/app/orf-finder, opened the same day, lists them).

**Key Extracted Points:**

1. `FindForwardOrfs`: "Find all ORFs in forward orientation with length in *base pairs* >= min_length_bp ... Returned range does not include the stop codon" — the minimum length does **not** count the stop (`to + 1 >= min_length_bp + from`, `to` = last base of the codon before the stop); the reported interval is extended by the stop (`AddInterval`: `to += 3`). `min_length_bp < 3` is raised to 3.
2. Stops per frame from the genetic code (`CTrans_table::IsOrfStop`); starts from `allowable_starts` (`GetStartCodons(genetic_code, include_atg, include_alt)` = ATG and/or the table's alternative starts); `"STOP"` in the list = stop-to-stop ORFs.
3. `longest_orfs = true` (ORFfinder): each ORF begins at the most upstream allowed start after the previous stop (`from = *starts.begin()`); `longest_orfs = false`: all starts of the stop-free stretch each yield an ORF (`for (auto s: starts) AddInterval(...)`) — exactly the Rosalind "every ATG" semantics.
4. The region after the last stop (`stops[frame].push_back(seq.size())`) is reported with `to_fuzz = true` (3'-partial, no stop appended) when it has a start: open-ended ORFs end at the last complete codon.
5. Reverse strand: the reverse complement is scanned and coordinates converted `from' = size − to − 1`, `to' = size − from − 1` (forward-strand coordinates, strand `minus`).
6. ProteinIQ (ORFfinder 0.4.3 wrapper) documents: "The terminating stop codon does not count toward the threshold", start modes "ATG only / ATG + code-specific alternative starts / Any sense codon", "ORFfinder reports each ORF from the most upstream start codon allowed", "Ignore nested ORFs: drop ORFs placed completely within another ORF on the searched strands" (a cross-frame filter, not in-frame nesting), and the lacI example 1,092 nt = 363 aa + stop (the reported range includes the stop).

### NCBI ORFfinder 0.4.3 standalone binary - measured (B09 A1-18 / A1-19, 2026-11-18)

**Source:** the unmodified NCBI program `https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz` (8 386 303 bytes, Last-Modified 2021-05-04, sha256 `f30a0806b984931a967d8024066a53dbb811035ca14564f36f7f7045a6a509f2`; unpacked binary sha256 `933cb6ff...9a3`; `ORFfinder -version` = `ORFfinder: 0.4.3`, Build-Signature `GCC_730-Release64--x86_64-unknown-linux3.10.0-gnu2.17-coremake1`, Build-Date May 4 2021) plus `USAGE.txt` (sha256 `368d0fe3...b409`) and `CHANGELOG.txt`. The sandbox of the audit refuses to execute a downloaded binary, so the repository owner ran the hand-off package `docs/Validation/review-2026-09/handoff/A1-18-orffinder/` in Docker (`--platform linux/amd64`, Debian 12 bookworm, `uname -m` x86_64, 18 CPUs) and returned `returned/A1-18-orffinder-out.zip` (`out/ENV.txt`, `USAGE.txt`, `jobs.tsv`, `runs.tsv`, `runs/<id>/{cmd,stdout,stderr,rc,meta,orfs}`). Authority rank 1 (the reference program itself); the archive was treated as untrusted data (text parsed with `python3 -I`, nothing executed).

**Run matrix (1750 invocations; exit codes 1738 x 0, 4 x 1, 3 x 254, 5 x 255, all the non-zero ones being deliberate invalid-option runs):** defaults; the factorial `-g {1,2,4,11} x -s {0,1,2} x -n {false,true} x -strand {both,plus,minus}` (`-ml 30`) in `-outfmt` 0 and 3 (1 and 2 on subsets); an `-ml` sweep (-5 ... 1000) with and without `-n`; every genetic code 0..33; `-c true/false` on 122 origin-spanning genomes and on the designed set; a `-b/-e` grid (24 pairs incl. invalid) x strand x `-n` x `-s` x `-c`; 116 designed records one per process; 24 edge files (length 0..3, 29, 30, lower/mixed case, N, IUPAC, U, gaps, CRLF, protein, ...); 6 circular genomes without a stop; `-h -help -xmlhelp -version -version-full -dryrun -logfile -conffile -id` and invalid option values. Inputs: 629 FASTA records generated by `make_inputs.py` (244 designed incl. reverse complements, 170 random, 122 circular, 60 region, 33 edge).

**Converted oracle** (`scripts/analysis/generate_orffinder_oracle.py`, deterministic): 620 distinct sequences, 260 option sets (`-c` merged: identical outputs), **54 820 (option set, sequence) results = 54 460 ORF lists (232 162 ORFs) + 360 rejected regions**, plus 7 rejected genetic codes. The converter asserts while reading: `-outfmt 0`, 3 (45 240 sequence results), 1 (5 835) and 2 (3 350) list the same ORFs in the same order and numbering, the `-outfmt 2` fuzz flags agree with the `partial` header flag, and 4 134 (option set, sequence) pairs measured by several jobs agree (factorial / minlength / circular / one-sequence-per-process). 27 records the FASTA reader rejects ("invalid residue(s)": very short, very C-rich, U, gap, protein) have no result.

**Format facts** (observed, `outfmt 0`): `>lcl|ORF<n>_<seqid>:<from>:<to> unnamed protein product[, partial]`, `<from>:<to>` 0-based inclusive with the stop codon, `from > to` = minus strand; `-outfmt 1` headers `lcl|<seqid>:<from+1>-<to+1>` / `c<hi+1>-<lo+1>`; `-outfmt 3` feature table has no `<`/`>` partial markers; `-outfmt 2` `fuzz-from lim lt` / `fuzz-to lim gt` give the open end.

**Measured semantics - the library equals the binary on every row (strand, forward coordinates, order, 5'/3' flags, protein hash):**

| Item | Measured |
|------|----------|
| `-n true` | drops an ORF whose forward span lies completely inside another ORF of any frame / strand. 15 306 (`-n false`, `-n true`) row pairs: 67 373 ORFs dropped = 58 732 strictly inside + 3 373 sharing only the start + 3 539 sharing only the end + **1 729 identical spans (plus-strand ORF dropped, minus-strand kept)**; 23 231 of the dropped ORFs are contained only by an ORF of the other strand; the containers are the ORFs left by `-ml` / `-strand`. Strict containment (variant tested: 1 786 mismatching rows) and "keep both of an equal pair" (324 rows) are refuted. |
| `-ml` | `max(ml, 30)`, then whole codons: coding nt >= `3 * floor(ml / 3)` (stop not counted for terminated ORFs; open ORFs count all nt). The sweep -5, 0, 1, 10, 29, 30, 31, 45, 60, 74, 75, 76, 90, 150, 300, 1000 on 414 sequences: 9 552 rows equal. The plain "coding >= ml" rule is refuted (384 mismatching rows at `-ml` 31, 74, 76). |
| `-b/-e` | 1-based inclusive, `-b 0` = 1, `-e 0` / beyond the end = to the end; `-b -5` ("Requested interval start -5 is out of range") and `-b 250 -e 100`, `-b 5 -e 3` ("Interval Start value is expected to be smaller than Stop") are errors; **cut-then-scan** with input-space coordinates reproduces all 7 080 region rows (partial flags relative to the region). |
| `-c true` | **no effect**: `orfs.txt` and stderr byte-identical to `-c false` for 36 / 36 circular pairs, 20 / 20 region pairs, 12 / 12 no-stop genomes (0 timeouts); no wrap-around ORF in any output. |
| partial | 5'-partial <=> strand offset of the start < 3 (22 000+ ORFs, 0 exceptions, also with ATG); 3'-partial <=> no stop codon. |
| protein | table translation, stop dropped, first residue `M` unless 5'-partial (254 010 ORFs, 0 exceptions). |
| `-s`, `-g` | `-s 0/1/2` and all tables 1-6, 9-16, 21-33 equal the scanner (9 936 + 10 736 + 9 552 sequence rows); tables 27, 28, 31 accepted with the codons the table lists as stops (488 rows each); 0, 7, 8, 17-20 rejected ("Unable to find genetic code number N"). |
| order | plus strand then minus; each by frame of the strand-coordinate start, then start (50 250 sequence results, numbering `ORF1..n` consecutive). |
| edge | lower / mixed case accepted; `N`/IUPAC follow the `OrfScanner` rules (`TAR` a stop, `NTG`/`NNN` never a start, N read through); `U`, `-`, protein rejected by the reader; `-strand x/left`: warning, both strands; `BOTH`, `Plus` accepted; `-n`/`-c` spellings `t f yes no TRUE 1` accepted, `maybe` rc 1; `-s -1/3/99` rc 254; `-outfmt -1/4/5` rc 255. |

### EMBOSS getorf 6.6.0 (binary run locally, Debian `emboss`)

`getorf -find 3 -table T -minsize M` = nucleotide ORFs START→STOP, one per stop (nested in-frame starts do not open an ORF), start set = the table's starts (EGC files; EGC.3 and EGC.13 carry older start sets than NCBI gc.prt 4.6 and are excluded from the comparison), `-minsize` compares the length without the stop, open ORFs written to the last complete codon, header `[s - e]` 1-based forward coordinates (`s > e` on the reverse strand). Corner: a START that is the last complete codon of its frame is never written by getorf (the NCBI `COrf` reports it).

### NCBI Genetic Codes — Standard code (transl_table=1)

**URL:** https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 2 (official NCBI specification)

**Key Extracted Points:**

1. **Standard start codon:** ATG (Methionine, M).
2. **Stop codons:** TAA, TAG, TGA — all designated "*" (termination).

---

## Documented Corner Cases and Failure Modes

### From Rosalind / Wikipedia

1. **Nested ORFs sharing a stop:** A downstream ATG in the same frame opens a second, shorter ORF terminated by the same stop; both protein candidates are reported (Rosalind sample: `MGMTPRLGLESLLE` vs `MTPRLGLESLLE`).
2. **No in-frame stop:** A reading begun at an ATG with no subsequent in-frame stop is NOT a complete ORF and yields no protein candidate (Rosalind translates "until a stop codon is encountered").
3. **Both strands:** ORFs must be searched on the reverse complement as well as the forward strand (six-frame search).

### From NCBI ORFfinder

1. **Minimum length filter:** ORFs shorter than the selected nucleotide length are excluded (length ≥ threshold).

### Review finding (B09, F23)

The earlier unit documents claimed the default `minLength` "matches NCBI ORFfinder's nucleotide length filter". The NCBI source shows the filter excludes the stop codon; this method's default keeps `Length ≥ minLength` including it (a self-consistent contract for existing callers) and offers `MinLengthExcludesStop` for the NCBI/getorf convention. The Rosalind problem itself defines no minimum length.

### Oracle numbers (B09, 2026-11-18)

`scripts/analysis/generate_orf_oracle.py` (reproducible byte-identical): 1746 (sequence, options) rows; (1) definition-based Python finder (Biopython 1.88 `translate`, `CodonTable`, `reverse_complement`), (2) line-by-line port of `COrf::FindForwardOrfs`/`s_FindOrfs` (1515 rows; nested + any-sense is not an NCBI mode), (3) EMBOSS getorf 6.6.0 `-find 3` (377 rows) — all asserted equal inside the generator; the library equals the TSV on all 1746 rows. Differential tests: `Translator.FindOrfs` (240 sequences × tables 1, 2, 4, 11) and `GenomeAnnotator.FindOrfs` (150 sequences) agree with `FindOpenReadingFrames` under the matching options.

### The three ORF engines (DUP_MAP §3) — measured relation

| | `GenomicAnalyzer.FindOpenReadingFrames` (options) | `Translator.FindOrfs` | `GenomeAnnotator.FindOrfs` |
|--|--|--|--|
| Scan core | `OrfScanner` (Core) | `OrfScanner` (since B09 F25) | own loop (B11), same algorithm as `NestedOrfs = true` |
| Starts | ATG / table starts / any sense | table starts (TTG/CTG/ATG in table 1) | ATG/GTG/TTG (Prodigal set, not a table) |
| Nested starts | yes (default) / no | no | yes |
| Open-ended ORFs | option | reported | only for `requireStartCodon = false` stop-to-stop mode |
| Min length unit | nt incl. stop (option: excl.) | amino acids | amino acids |
| Coordinates | scanned strand + `GetForwardStrandRange` | scanned strand, frames −1..−3 | forward strand, half-open |

---

## Test Datasets

### Dataset: Rosalind ORF sample (`Rosalind_99`)

**Source:** Rosalind ORF problem, https://rosalind.info/problems/orf/ (retrieved 2026-06-14).

| Parameter | Value |
|-----------|-------|
| Input DNA | `AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG` |
| Expected distinct proteins (4) | `MLLGSFRLIPKETLIQVAGSSPCNLS`, `M`, `MGMTPRLGLESLLE`, `MTPRLGLESLLE` |
| Frames searched | 6 (forward + reverse complement, 3 each) |
| Start codon | ATG |
| Stop codons | TAA, TAG, TGA |

Independent re-derivation (Standard code, every-ATG-to-first-in-frame-stop, both strands, distinct proteins) reproduces exactly these 4 proteins.

### Dataset: Single canonical forward ORF (derivation)

**Source:** Standard genetic code (NCBI transl_table=1); construction.

| Parameter | Value |
|-----------|-------|
| Input DNA | `ATG` + `AAA` + `AAA` + `TAA` = `ATGAAAAAATAA` |
| ORF DNA span (incl. stop) | `ATGAAAAAATAA` (12 nt, divisible by 3) |
| Protein candidate (excl. stop) | `MKK` |
| Start position (0-based) | 0, frame 1 |

---

## Assumptions

1. **ASSUMPTION: ORF span includes the terminating stop codon in the reported DNA `Sequence`.** Sources define the ORF as start-codon to stop-codon; the reported nucleotide span includes the stop so that `Length % 3 == 0` and the boundaries are explicit (Wikipedia: "length divisible by three … bounded by stop codons"). The translated protein candidate excludes the stop (Rosalind translates "until a stop codon"). This affects the reported `Sequence`/`Length` only, not which ORFs are detected.
2. **ASSUMPTION: `minLength` is measured in nucleotides** (matching NCBI ORFfinder's nucleotide length options), inclusive lower bound (length ≥ minLength). Default 100 nt retained from the existing public API; it is a caller-supplied parameter and any value is honored.
3. **ASSUMPTION: Standard start codon ATG only by default** (the Rosalind / "ATG only" convention; the alternative initiation codons and any sense codon are `OrfSearchOptions.StartCodons`, implemented since F25 and measured on the ORFfinder 0.4.3 binary).

---

## Recommendations for Test Coverage

1. **MUST Test:** Rosalind sample — six-frame, both strands, returns exactly the 4 distinct proteins. — Evidence: Rosalind ORF dataset.
2. **MUST Test:** Single forward ORF `ATGAAAAAATAA` → one ORF, sequence `ATGAAAAAATAA`, position 0, frame 1, protein `MKK`. — Evidence: Standard code derivation.
3. **MUST Test:** Nested ORFs sharing a stop are both reported. — Evidence: Rosalind (`MGM…`/`MTP…`).
4. **MUST Test:** Reading begun at ATG with no in-frame stop yields no ORF. — Evidence: Rosalind / Wikipedia.
5. **MUST Test:** Reverse-complement-only ORF is detected with `IsReverseComplement=true`. — Evidence: Rosalind six-frame requirement.
6. **MUST Test:** minLength filtering excludes shorter ORFs and includes exactly-at-threshold ORFs. — Evidence: NCBI ORFfinder length filter.
7. **MUST Test:** Invariants — every ORF starts with ATG, ends with a stop codon, length divisible by 3. — Evidence: Wikipedia/NCBI definition.
8. **SHOULD Test:** Lowercase input handled (case-insensitive). — Rationale: real-data robustness.
9. **SHOULD Test:** Multiple distinct stop codons (TAA/TAG/TGA) all recognized. — Rationale: completeness of stop set.
10. **COULD Test:** Empty / too-short sequence → empty result. — Rationale: edge guard.

---

## References

1. Rosalind. "ORF — Open Reading Frames." https://rosalind.info/problems/orf/ (accessed 2026-06-14).
2. Wikipedia contributors. "Open reading frame." https://en.wikipedia.org/wiki/Open_reading_frame (accessed 2026-06-14).
3. NCBI. "ORFfinder (Open Reading Frame Finder)." https://www.ncbi.nlm.nih.gov/orffinder/ (accessed 2026-06-14).
4. NCBI. "The Genetic Codes (transl_table=1, Standard)." https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (accessed 2026-06-14).

---

## Change History

- **2026-11-18 (B09 A1-18 / A1-19, F57)**: NCBI ORFfinder 0.4.3 run by the user (Docker, 1750 invocations, 629 inputs); converter `generate_orffinder_oracle.py`, oracle `orffinder_cli_oracle.tsv` (54 820 results, 232 162 ORFs); implemented `IgnoreNestedOrfs`, `RegionStart/End`, `Strands`, `OrfFinderMinimumLength`, `DualCodingStopsTerminate`, `IsStartPartial`, `OrfFinderCompatible`; `-c` recorded as a no-op of the binary.
- **2026-11-18 (B09 A1-21)**: string overloads `FindOpenReadingFrames(string, …)` for IUPAC input (ambiguous codon = stop only if every concrete codon is a stop, never a start; = `Translator.FindOrfs(string)` / EMBOSS getorf rules); oracle `generate_orf_oracle.py --ambiguous` (2844 rows: Biopython definition; getorf 6.6.0 on the 184 A/C/G/T/N rows it resolves exactly — getorf mis-resolves OR-ed codons such as TRA); MCP `find_open_reading_frames` accepts IUPAC input.
- **2026-06-14**: Initial documentation.
