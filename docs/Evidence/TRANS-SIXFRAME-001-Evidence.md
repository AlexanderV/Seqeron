# Evidence Artifact: TRANS-SIXFRAME-001

**Test Unit ID:** TRANS-SIXFRAME-001
**Algorithm:** Six-Frame Translation and Open Reading Frame (ORF) finding
**Date Collected:** 2026-06-13

---

## Online Sources

### EMBOSS transeq — application documentation (reference implementation, frame numbering)

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established library — EMBOSS)

**Key Extracted Points:**

1. **Reverse-frame definition (phase-locked convention):** "The reverse frame '-1' is defined as the translation you get when you use the reverse-complement of the sequence with the same codon phase as the codon in frame '1'." Frames -2 and -3 use the phase of frames 2 and 3 respectively.
2. **Worked example (phase-locked):** For sequence `ACTGG`, frame -1 translation is `S` (phase of frame 1: reverse-complement read as `...AGT`), frame -2 is `QX`, frame -3 is `PV`.
3. **Documented alternative convention (the one used here):** The doc explicitly records a second, widely used convention: "The alternative way of generating the reverse translation frames used by some people is that frame -1 is made by taking the frame '1' of the reverse complement. There is no correspondence between the codons used in frame 1 and -1 …" — i.e., frame -1 = offset 0 of the reverse complement.
4. **Frame option values:** Permitted `-frame` values are `1, 2, 3, F` (forward three), `-1, -2, -3, R` (reverse three), and `6` (all six frames) — confirming exactly six frames.

### Biopython `Bio.SeqUtils.six_frame_translations` — reference implementation (governing convention)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — Biopython)

**Key Extracted Points:**

1. **Reverse strand computed once:** `anti = reverse_complement(seq)` (`reverse_complement_rna` for RNA).
2. **Reverse-frame offsets (independent / "alternative" convention):** the reverse frames are filled as `frames[-(i + 1)] = translate(anti[i : i + fragment_length], genetic_code)[::-1]` for `i in range(3)`. Thus frame -1 = translate(reverse-complement, offset 0), frame -2 = offset 1, frame -3 = offset 2. The trailing `[::-1]` only reverses the *string for positional display alignment* under the forward sequence; the residue multiset/content of each reverse-frame protein is the translation of the reverse complement read 5'→3' at that offset.
3. **Forward frames:** `frames[i + 1] = translate(seq[i : i + fragment_length], genetic_code)` for `i in range(3)`, i.e. forward frames are offsets 0, 1, 2 of the input.
4. **Six frames keyed ±1..±3:** the dictionary uses keys `+1,+2,+3` and `-1,-2,-3`.

### NCBI — The Genetic Codes, Standard Code (transl_table=1)

**URL:** https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
**Accessed:** 2026-06-13
**Authority rank:** 2 (official NCBI specification / standard)

**Key Extracted Points:**

1. **Codon table 1 rows (verbatim):**
   `AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG`
   `Starts = ---M------**--*----M---------------M----------------------------`
   `Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG`
   `Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG`
   `Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`
2. **Start codons (the three `M` positions in the Starts line):** TTG, CTG, ATG.
3. **Stop codons:** TAA, TAG, TGA.
4. **Initiator translated as Met:** "The initiator codon - whether it is AUG, CTG, TTG or something else, - is by default translated as methionine (Met, M)." `FindOrfs` reports the initiator as `M` (review 2026-09; getorf `-methionine`).

### EMBOSS getorf — application documentation (ORF definition)

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — EMBOSS)

**Key Extracted Points:**

1. **ORF definitions:** getorf recognises an ORF as either "a region that is free of STOP codons" or "a region that begins with a START codon and ends with a STOP codon."
2. **`-find` option values:** 0 = translation between STOP codons; 1 = translation between START and STOP codons; 2/3 = corresponding nucleic sequences; 4/5/6 = nucleotides flanking START / initial STOP / ending STOP codons. Option **1 (START→STOP)** is the legacy `FindOrfs` model; all seven modes are available through `FindOrfs(DnaSequence, OrfFindMode, …)` (2026-10-14, A1-12).
3. **Minimum size:** default `-minsize` is 30 bases (= 10 amino acids); it measures nucleotide length.
4. **Both strands & reverse numbering:** by default both forward and reverse strands are searched; reverse-strand ORFs are labelled `(REVERSE SENSE)`.

### Wikipedia — Reading frame (uses cited primaries: Lodish 2007; Pierce 2012; Badger & Olsen 1999; Anderson et al. 1981)

**URL:** https://en.wikipedia.org/wiki/Reading_frame
**Accessed:** 2026-06-13
**Authority rank:** 4 (encyclopedic, citing primary sources)

**Key Extracted Points:**

1. **Six frames:** "any given sequence of DNA can therefore be read in six different ways: Three reading frames in one direction (starting at different nucleotides) and three in the opposite direction."
2. **Reverse strand direction:** the three additional frames "may be read from the other, complementary strand in the 5′→3′ direction along this strand"; the 5′→3′ direction on the second strand corresponds to 3′→5′ on the first — i.e. the reverse-complement read forward.

### EMBOSS getorf — C source (opened 2026-09-28, B02 code review)

**URL:** https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/getorf.c and `emboss/acd/getorf.acd`
**Authority rank:** 3 (reference implementation source)

1. `-find 1` (P_START2STOP): an ORF opens at a START only when none is open in that frame (`codon == START && !ORF[frame]`) — nested starts are not reported.
2. `-methionine` (acd default `Y`, "Change initial START codons to Methionine"): the first residue is appended as `'M'`.
3. On a STOP: `getorf_WriteORF(..., start[frame], pos-1, ...)` — the printed range EXCLUDES the stop codon.
4. At the end without a STOP (`pos >= seqlen-5`): the last codon is appended and `WriteORF(..., pos+2)` — the ORF ends at the last base of the last complete codon.
5. `minsize` (default 30 nt) is divided by 3 and compared to the protein length.
6. EMBOSS ships genetic codes EGC.1 … but not EGC.25–33 (HTTP 404 for 25, 26, 27, 28, 31, 33), so getorf defines nothing for the dual-coding tables 27/28/31.

### EMBOSS getorf — `-find` 0–6 (getorf.c / getorf.acd, opened 2026-10-14, FIN-B02 WP11)

**URL:** https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/getorf.c and `emboss/acd/getorf.acd`
**Authority rank:** 3 (reference implementation source)

1. acd `-find` values: "0:Translation of regions between STOP codons; 1:Translation of regions between START and STOP codons; 2:Nucleic sequences between STOP codons; 3:Nucleic sequences between START and STOP codons; 4:Nucleotides flanking START codons; 5:Nucleotides flanking initial STOP codons; 6:Nucleotides flanking ending STOP codons". Defaults: `-minsize 30`, `-maxsize 1000000` ("Minimum/Maximum nucleotide size of ORF to report"), `-flanking 100`, `-reverse Y`, `-methionine Y`, `-circular N`.
2. `main`: `if(find == P_STOP2STOP || find == P_START2STOP || find == AROUND_START) { minsize /= 3; maxsize /= 3; }` — size bounds compared with the residue count for modes 0, 1, 4, with the nucleotide count for 2, 3, 5, 6.
3. `getorf_FindORFs` (non-circular): STOP→STOP modes start with `ORF[f] = ajTrue; start[f] = f` ("assume already in a ORF so we get ORFs at the start of the sequence"); `if(codon == STOP || pos >= seqlen-5)`: the last codon (if not a STOP) is appended, the region is written (`WriteORF(start, pos-1)` on a STOP, `pos+2` at the end), and `start[frame] = pos+3`. START modes: see `-find 1` above.
4. `getorf_AppORF`: modes 2, 3, 5, 6 store the codon's nucleotides; 0, 1, 4 the residue.
5. `getorf_WriteORF`: AROUND_START `codonpos = s` (the START); AROUND_INIT_STOP `codonpos = s-3` (the STOP before the region); AROUND_END_STOP `codonpos = e+1` (the terminating STOP); window `[codonpos − around, codonpos + around + 2]`, `return` (no output, no ORF number consumed) if `s < 1` or `e > seqlen`. Reverse sense: `s = len-s+1; e = len-e+1`, then shifted into range; header `[s - e] (REVERSE SENSE)`, "Around codon at <codonpos>" (codonpos in reverse-complement coordinates).
6. acd `-flanking` help: "If the region of flanking nucleotides crosses the start or end of the sequence, no output is given for this codon."

### EMBOSS 6.6.0 getorf binary — reference output (2026-10-14)

`getorf -sequence in.fa -outseq stdout -find N -table T -minsize M -maxsize X -flanking F -reverse R -methionine Y|N -auto`. EMBOSS EGC files equal the NCBI/Biopython tables for 1, 2, 4, 5, 6, 9, 10, 11, 12, 14, 15, 16, 21, 22, 23; EGC.3 lacks start GTG and EGC.13 lacks starts ATA/GTG/TTG (older NCBI start sets), so tables 3/13 were compared only in the START-independent modes 0/2/5/6. EMBOSS also writes empty regions when `-minsize 0` (adjacent STOPs).

Worked example `GGGATGAAACCCTAAGGG`, `-minsize 9 -flanking 2 -table 1`:
- `-find 0`: `[1 - 12] GMKP`, `[8 - 16] NPK`, `[3 - 17] DETLR`, `[11 - 3] (REVERSE SENSE) GFI`, `[16 - 2] (REVERSE SENSE) LRVSS`, `[18 - 1] (REVERSE SENSE) PLGFHP`.
- `-find 2`: same ranges, `GGGATGAAACCC`, `AACCCTAAG`, `GATGAAACCCTAAGG`, `GGTTTCATC`, `CTTAGGGTTTCATCC`, `CCCTTAGGGTTTCATCCC`.
- `-find 3`: `[4 - 12] ATGAAACCC`. `-find 4`: `[2 - 8] Around codon at 4. GGATGAA`. `-find 5`: `[3 - 9] Around codon at 5. GATGAAA`, `[16 - 10] (REVERSE SENSE) Around codon at 5. CTTAGGG`. `-find 6`: `[11 - 17] Around codon at 13. CCTAAGG`.
- `TAATAA -find 0 -minsize 0 -reverse N`: `[1 - 0]` (empty), `[2 - 4] N`, `[3 - 5] I`, `[4 - 3]` (empty).

Bulk comparison (C# `FindOrfs(dna, OrfFindMode, …)` mapped by `OrfFeature.ToGetOrfRange` vs the binary): 2 500 random sequences × 290 runs, 0 mismatching (sequence, run) pairs over 3 474 722 reference records; 25 084 records excluded as the documented one-codon open-ORF divergence (modes 1/3/4). Locked subset: `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Getorf/getorf_find_modes.tsv`.

### EMBOSS — ambiguous codons and `-circular` (ajtranslate.c, getorf.c; opened 2026-10-14, FIN-B02 WP12)

**URL:** https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajtranslate.c and `emboss/getorf.c`
**Authority rank:** 3 (reference implementation source)

1. `ajTrnCodonstrTypeC`: `*aa = trnObj->GC[tc1][tc2][tc3]; if(trnObj->Starts[tc1][tc2][tc3] == 'M') return 1; if(*aa == '*') return -1; return 0;` — `trnconv` maps A C G T/U R Y S W K M B D H V (and N / any other symbol to index 14).
2. `ajTrnReadFile` fills `Starts` only from the 64 Base1/Base2/Base3 codons (all other entries `'-'`), so an ambiguous codon is never a START.
3. `trnComplete`: ambiguous `GC` entries = third-base wobble ("Calculate wobble in base3 and populate") plus, per residue, one "most ambiguous" codon built by OR-ing the bases of its codons; everything else stays `'X'`. Dumped with `transeq` over all 3 375 IUPAC codons, the OR-ed entries disagree with the all-expansions rule for a few codons: table 1/11 `TRR` → `*` (TGG = W), `TRA` (TAA/TGA) → `X`, `YTN` → L (TTT/TTC = F), `MGN` → R, `ATS` → I, `WSN` → S; table 2 `WRR` → `*`, `YTN` → L, `WSN` → S.
4. `getorf.c main`: "if the sequence is circular, append it to itself to triple its length so can deal easily with wrapped ORFs, but don't update len"; the reverse pass is `ajSeqReverseForce` of the tripled sequence.
5. `getorf_FindORFs` circular branch: `ORF[f] = ajFalse` for every mode; frames with `LASTORF[frame]` are skipped ("don't want to find extra ORFs when already been round circ"); a START opens an ORF only `if(pos < len)`; a STOP / end codon at `pos >= len` sets `LASTORF` (STOP→STOP: `ORF[frame] = ajFalse`, no reopening); after the loop, for `-find 0/2`, a frame without any STOP (`!GOTSTOP[f]`) is written as `WriteORF(seq, …, f, seqlen-1, ajTrnSeqOrig(trnTable, seq, f+1))` — frame 3 tested with `ajSeqGetLen(pep) >= maxsize` (frames 1, 2: `<= maxsize`).
6. `getorf_WriteORF`: `while(s > len) { s -= len; e -= len; }` ("Move the reported positions back to start in the range 1..len"); flank windows checked against `seqlen` (3L); "(ORF crosses the breakpoint)" when `s > len || e > len`.
7. `ajTrnSeqOrig` → `ajTrnSeqDangleC`: two dangling bases translated as `GC[b1][b2][N]`, one as `X`.

### EMBOSS 6.6.0 getorf binary — IUPAC and `-circular` reference output (2026-10-14, FIN-B02 WP12)

- `-find 1 -table 1 -minsize 3 -reverse N`: `ATGNNNCCCTAA` → `[1 - 9] MXP`; `atgaaacccrtaa` → `[1 - 12] MKPX`; `AUGAAACCCUAAGG` → `[1 - 9] MKP` (U written as T; input case kept in nucleotide output); `YTGAAACCCTAA` → no ORF; `ATGAAACCCTRR` → `[1 - 9] MKP` (TRR treated as STOP); `ATGAAACCCTRA` → `[1 - 12] MKPX` (TRA read through).
- `CCCTAAGGGATGAAA -circular Y -reverse N -minsize 3`: `-find 1` `[10 - 18] (ORF crosses the breakpoint) MKP` (linear: `[10 - 15] MK`); `-find 3` `ATGAAACCC`; `-find 0` `[7 - 18] GMKP`, `[14 - 25] NPKG`; `-find 6 -flanking 2` `[17 - 23] Around codon at 19. CCTAAGG`, `[24 - 30] Around codon at 26. GATGAAA`.
- `GCGCGCCCGGG -find 0 -minsize 3 -circular Y -reverse N`: `[1 - 33] ARPGRARGAPG`, `[2 - 33] RARGAPGARPG`; `-find 2` prints the same two proteins; with `-maxsize 30`: only `[3 - 33] APGARPGRARX` (frame-3 `>= maxsize` test).

Bulk comparison (C# string overload vs the binary; 1 000 random sequences per suite, 1–400 nt, × 93 runs: modes 0–6 × tables 1/2/11 × 4 minsize/maxsize/flanking/reverse sets + `-methionine N`):

| Suite | Reference records | (sequence, run) pairs | exact | residue-only (getorf X / inconsistent EMBOSS residue) | EMBOSS STOP-table codon present | unexplained |
|---|---|---|---|---|---|---|
| IUPAC, linear | 384 927 | 93 000 | 86 052 | 4 653 | 2 295 | 0 |
| IUPAC, `-circular Y` | 344 889 | 93 000 | 83 790 | 6 516 | 2 694 | 0 |
| A/C/G/T, `-circular Y` | 447 526 | 93 000 | 93 000 | 0 | 0 | 0 |

In the circular suites the full-length records of STOP-free frames are compared with getorf's literal block (frame-3 `>= maxsize`, protein in mode 2), which is where this API diverges by design. Locked subset: `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Getorf/getorf_iupac_circular.tsv` (24 sequences × 86 runs, 6 143 records; sequences free of the inconsistent codons of item 3 and of STOP-free circular frames).

### Biopython 1.88 — six frames of IUPAC strings (2026-10-14, FIN-B02 WP12)

`Seq(s)[f:].translate(table=t)` and `Seq(s).reverse_complement()[f:].translate(table=t)` (complete codons) vs `Translator.TranslateSixFrames(string, GeneticCode)`: 20 000 random strings (0–200 nt, 0 / 5 / 20 / 50 / 100 % ambiguity codes, 30 % lower case, 30 % U for T; U → T before Biopython) over all 27 Biopython tables — 0 mismatching frames. Worked values: `ATGNNNTARGCNRAYTAAYTG` table 1 → `MX*AB*L`, `XXXXXN`, `XXXXLX` / reverse `QLXXLXH`, `X*XAXX`, `XXXXXX`; `TTGTANSARGAYMGR` table 11 → `LXZDR`, `CXXX`, `VXXX` / `XXXXQ`, `XSXY`, `XXXT`.

### EMBOSS 6.6.0 transeq — frame conventions (2026-10-14, FIN-B02 WP17, audit A2-10)

Sources read: `emboss/transeq.c`, `emboss/acd/transeq.acd`, `ajax/core/ajtranslate.c` (raw.githubusercontent.com/kimrutherford/EMBOSS/master); binary `/usr/bin/transeq` (EMBOSS:6.6.0.0).

1. **Frame dispatch (transeq.c):** frame i < 3 → `ajTrnSeqOrig(seq, i+1)`; reverse frames → `ajTrnSeqOrig(seq, 2-i)` (frames −1..−3, comment "frame -1 uses frame 1 codons") by default, `ajTrnSeqOrig(seq, -i-1)` (frames −4..−6, "frame -1 uses codons starting at end") with `-alternative Y`. acd: "The default definition of frame '-1' is the reverse-complement of the set of codons used in frame 1 … a common standard, used by the Staden package … If you prefer to define frame '-1' as using the set of codons starting with the last codon of the sequence, then set this to be true."
2. **Phasing (ajtranslate.c `ajTrnSeqFrameC`):** frame −k (k = 1..3) = `ajTrnSeqRevC(&seq[k-1], len-k+1)` — the full codons of forward frame k read backwards and complemented; frame −(k+3) = `ajTrnSeqAltRevC(seq, len-k+1)` — codons from base `len-k` downwards, i.e. `reverse_complement(seq)[k-1:]` (the Biopython phasing).
3. **Dangling bases (`ajTrnSeqDangleS`, called by `ajTrnSeqOrig`):** forward frame k: `(len-k+1) mod 3` trailing bases; alternative −k: `(len-k+1) mod 3`; default −k: always `k-1` (the leading forward bases). Two bases are translated as `GC[b1][b2][N]` (reverse: `comp(seq[1]) comp(seq[0]) N`) — a residue when the third base does not matter, else X; one base → `X`. For `len = 1`, frame −3 reads past the base and prints X.
4. **`-trim`** strips trailing X/`*`; **`-clean`** replaces `*` by X (both default N, post-processing only). Stops are printed `*`.
5. **Ambiguity:** transeq's `trnComplete` table (item 3 of the getorf section above) never emits B/Z/J and differs from the all-expansions rule on 28 / 15 / 28 of the 3 375 IUPAC codons for tables 1 / 2 / 11 (e.g. table 1 `RAY` → X vs B, `ATS` → I vs X, `ATW` → X vs I, `TRA` → X vs `*`, `TRR` → `*` vs X).
6. **Worked example** (`transeq -frame 6`, `ATGAAACCCGGGTT`): 1 `MKPGX`, 2 `*NPGX`, 3 `ETRV`, −1 `PGFH`, −2 `TRVSX`, −3 `NPGFX`; with `-alternative Y`: −1 `NPGFX`, −2 `TRVSX`, −3 `PGFH`. `A`: `X`, ``, ``, ``, `X`, `X` (default) / `X`, ``, ``, `X`, ``, `` (alternative).
7. **Binary comparison:** `Translator.TranslateSixFrames(string, code, SixFrameConvention.EmbossTranseq | EmbossTranseqAlternative)` vs `transeq -frame 6 -table T -alternative N|Y` for 2 000 random A/C/G/T sequences (1–150 nt; 200 of 1–8 nt; 1/7 lower case, 1/11 with U; seed 20261014) × tables 1/2/11 × N/Y: 72 000 frames, **0 mismatches**. 600 random IUPAC sequences (1–120 nt, ~45 % ambiguity codes) × the same 6 runs: 21 600 frames, every residue of both tools equals its own per-codon table applied to the modelled codon list (forward, alternative-reverse and phase-locked-reverse codons plus dangling `b1 b2 N`), 5 442 residues differ, all at codons of item 5. Locked: `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Transeq/transeq_six_frames.tsv` (80 A/C/G/T(U) sequences × tables 1/2/11 × N/Y = 480 rows).

### EMBOSS 6.6.0 transeq — `-regions`, `-methionine`, `-trim`, `-clean` (2026-10-14, FIN-B02 WP21, audit A3-3)

Sources read: `emboss/transeq.c`, `emboss/acd/transeq.acd`, `ajax/core/ajrange.c`, `ajax/core/ajtranslate.c` (raw.githubusercontent.com/kimrutherford/EMBOSS/master); binary `/usr/bin/transeq` (EMBOSS 6.6.0), data `EGC.1`/`EGC.2`/`EGC.11` (Starts: table 1 TTG CTG ATG; table 2 ATT ATC ATA ATG GTG; table 11 TTG CTG ATT ATC ATA ATG GTG).

1. **Order (transeq.c `main`):** `wholeseq = ajRangeIsWhole(regions)` (true only when `-regions` is empty); if not whole, `ajRangeSeqExtract(regions, seq)`; per frame `ajTrnSeqInit` when `methionine && !wholeseq`, else `ajTrnSeqOrig`; then `if(trim) transeq_Trim(pep)`; `if(clean) transeq_Clean(pep) /* clean after the trim */`.
2. **`-regions` (ajrange.c `ajRangeSeqExtract`):** `ajStrAppendSubS(&outstr, seq, st - 1, en - 1)` for each region in the given order — the regions are concatenated and the result replaces the sequence (all six frames translate the concatenation). Range parsing rejects `from < 1` and `from > to`; an end past the sequence end is clamped (binary: `-regions 5-40` on 18 nt = `5-18`); a start past the end reads the last base (binary: `1-3,20-25` on 18 nt → `ATG` + `G`) — rejected by the API instead.
3. **`-methionine` (acd: default `Y`, "Change initial START codons to Methionine"; ajtranslate.c):** `ajTrnSeqInitC` (forward frames): `if(lenmod3 > 2)` and `Starts[c1][c2][c3] == 'M'` → first residue `M`; `ajTrnSeqRevinitC` (phase-locked reverse frames, codon = complement of the last whole frame codon read backwards): `end = (len/3)*3-1; if(iend > 2)` — needs two whole codons; `ajTrnSeqAltRevinitC` (`-alternative` reverse frames): `iend = len-1; if(iend > 1)` — one codon. The Starts table is filled for the 64 A/C/G/T codons only (`ajTrnReadFile`), so an ambiguous codon is never replaced. Binary, table 11 `-regions 1-3` on `CAA`: frame −1 `L` (phase-locked) vs `M` (`-alternative Y`); `GGGCAA` `-regions 1-6`: −1 `MP` under both.
4. **`-trim` (`transeq_Trim`):** scan from the right end while the residue is `X` or `*`, truncate there (acd: "removes all 'X' and '*' characters from the right end"). **`-clean` (`transeq_Clean`):** `ajStrExchangeSetCC(&str, "*", "X")`.
5. **Worked example (binary, table 11, frame 1, `TTGAAATAGCCCTGA`):** whole `LK*P*`; `-regions 1-6,13-15,7-9` → `TTGAAA`+`TGA`+`TAG` → `MK**`; `-methionine N` → `LK**`; `-trim -clean` → `MK`; whole `-clean` → `LKXPX`.
6. **Binary comparison:** `Translator.TranslateSixFrames(string, code, EmbossTranseq | EmbossTranseqAlternative, TranseqOptions)` vs `transeq -frame 6 -table T -alternative A -trim TR -clean CL -methionine M [-regions R]` over 1 000 random A/C/G/T sequences (1–80 nt, 15 % start-codon insertions, 1/10 lower case, 1/10 with U; seed 2103) × tables 1/2/11 × N/Y × trim × clean × regions {none, `4-30`, `2-10,15-40`, `1-3,10-12,20-200`, `7-7,1-20`, `3-14,1-9`} (sequences at least as long as the largest region start) × methionine Y/N (with regions): 236 544 rows = 1 419 264 frames, **0 mismatches**. 60 IUPAC sequences (seed 2104, 14 544 rows): untrimmed frames have equal lengths and differ only at ambiguous codons resolved differently by transeq's table (item 5 of the conventions section: J/L/Z/B/R → X, X → L/*/R/I, …); with `-trim` the trimmed length follows those residues. Locked subset: `TestData/Transeq/transeq_options.tsv` (5 sequences covering 205 of the 252 methionine/trim effects seen in the sweep, all 264 parameter runs, 1 320 rows).

Generator (summary; run with `python3 -I`): for each `(table, alt, trim, clean, regions, methionine)` write the selected sequences as FASTA, run `transeq -sequence in.fa -outseq stdout -frame 6 -table T -alternative A -trim TR -clean CL -methionine M [-regions R] -auto`, split the output into records (six per input, frames 1, 2, 3, −1, −2, −3; empty records kept) and emit `table, alt, trim, clean, methionine, regions|-, sequence, f1..f-3`.

### INSDC Feature Table Definition — CDS (WebSearch snippet; insdc.org / ddbj blocked for fetch)

"CDS: sequence of nucleotides that corresponds with the sequence of amino acids in a protein (location includes stop codon)". This is the stop-inclusive end-coordinate convention used by `OrfResult.EndPosition`.

### Biopython 1.88 `Bio.Seq._translate_str` (installed source)

`cds=True`: first codon must be a start and "will be translated as methionine"; final codon must be a stop and is excluded. For tables whose stop codons are also in `forward_table` (27: TGA; 28: TAA, TAG, TGA; 31: TAA, TAG — i.e. every stop codon of these tables) `to_stop=True` raises `ValueError`.

---

## Documented Corner Cases and Failure Modes

### From EMBOSS transeq / Biopython

1. **Incomplete trailing codon:** translation stops at the last complete codon; trailing 1 or 2 nucleotides that cannot form a full codon are ignored (`fragment_length = 3 * ((length - i) // 3)` in Biopython).
2. **Reverse-frame numbering ambiguity:** two conventions exist (phase-locked vs. independent offset). The repository defaults to the Biopython "independent offset" convention (frame -1 = reverse-complement offset 0), the documented "alternative" in EMBOSS transeq; `SixFrameConvention.EmbossTranseq` gives transeq's phase-locked default and `EmbossTranseqAlternative` its `-alternative Y` (both with transeq's dangling-base translation).

### From EMBOSS getorf

1. **No START codon present:** under `-find 1` (START→STOP) no ORF is emitted for a region lacking a START codon.
2. **ORF running off the end of the sequence:** a region beginning at a START codon with no downstream STOP before the sequence end is an incomplete ORF.
3. **minsize filtering:** ORFs shorter than the minimum are discarded.

### From NCBI Standard Code

1. **IUPAC-ambiguous codons:** codons containing ambiguity codes (N, R, Y, …) are not in the 64-codon table.

---

## Test Datasets

### Dataset: Six-frame translation of `ACTGG` (EMBOSS transeq worked example)

**Source:** EMBOSS transeq documentation (worked example), cross-checked with Biopython convention.

| Frame | Phase-locked (EMBOSS default) | Independent offset (Biopython / this repo) |
|-------|-------------------------------|--------------------------------------------|
| +1 | T (ACT,GG→T) | T |
| -1 | S | P (reverse-complement `CCAGT` offset 0 → CCA=P) |

Reverse complement of `ACTGG` = `CCAGT`. This repo returns frame -1 = `P`; EMBOSS default returns `S`. The conflict is documented (different conventions).

### Dataset: Six-frame translation of a 39-nt sequence (derived from cited algorithm)

**Source:** Standard genetic code (NCBI table 1) applied per Biopython six-frame algorithm.

Input DNA: `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG` (length 39).
Reverse complement: `CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT`.

| Frame | Protein |
|-------|---------|
| +1 | `MAIVMGR*KGAR*` |
| +2 | `WPL*WAAERVPD` |
| +3 | `GHCNGPLKGCPI` |
| -1 | `LSGTLSAAHYNGH` |
| -2 | `YRAPFQRPITMA` |
| -3 | `IGHPFSGPLQWP` |

### Dataset: ORF finding (forward strand, START→STOP)

**Source:** EMBOSS getorf `-find 1` model; positions derived from the START→STOP definition.

Input DNA: `GGGATGAAACCCTAAGGG`. ATG begins at index 3; stop `TAA` occupies indices 12–14.

| Field | Value |
|-------|-------|
| StartPosition (0-based, start codon first base) | 3 |
| EndPosition (0-based, stop codon last base, inclusive; getorf prints 1-based [4 - 12], stop excluded) | 14 |
| Frame | 1 |
| Protein (start residue included, stop excluded) | `MKP` |
| AminoAcidLength | 3 |
| NucleotideLength (End − Start + 1) | 12 |

---

## Assumptions

1. **ASSUMPTION: Reverse-frame numbering convention.** Two documented conventions exist (EMBOSS phase-locked vs. Biopython independent-offset). The repository defaults to the **Biopython** convention (frame -k = reverse-complement offset k−1; EMBOSS conventions selectable via `SixFrameConvention` since FIN-B02 WP17), which is the dominant reference-implementation behaviour and is explicitly listed as an accepted alternative in the EMBOSS transeq documentation. This is a convention choice, not an invented value; both produce correct biology, only the −1/−2/−3 labels differ.
2. **ASSUMPTION: Stop codons rendered as `*`; ambiguous IUPAC codons resolved as in Biopython (all-same → that residue, all-stop → `*`, D/N → `B`, E/Q → `Z`, I/L → `J`, otherwise `X`; TRANS-CODON-001 review 2026-09).** The `*` for stop is universal (NCBI). The ambiguity resolution is Biopython's `Bio.Data.CodonTable` ambiguous-table rule, implemented and oracle-tested in `GeneticCode.Translate` (TRANS-CODON-001); it is not exercised as a six-frame-specific MUST.
3. **ASSUMPTION: ORF length filter is in amino acids.** getorf's `-minsize` is in nucleotides; the repository's `FindOrfs(minLength)` parameter counts amino acids (protein length). This is an API-shape choice documented in the contract; behaviour is well-defined for any value. The getorf-modes overload (`FindOrfs(DnaSequence, OrfFindMode, …)`) takes getorf's nucleotide `minSize`/`maxSize` unchanged.

---

## Recommendations for Test Coverage

1. **MUST Test:** `TranslateSixFrames` returns exactly six frames keyed +1,+2,+3,−1,−2,−3 with the exact protein strings of the 39-nt dataset — Evidence: Biopython six_frame_translations algorithm + NCBI table 1.
2. **MUST Test:** Forward frames +1/+2/+3 equal `Translate` at offsets 0/1/2 — Evidence: Biopython forward-frame loop.
3. **MUST Test:** Reverse frames −1/−2/−3 equal translation of the reverse complement at offsets 0/1/2 — Evidence: Biopython reverse-frame loop.
4. **MUST Test:** Incomplete trailing codon is ignored (length not a multiple of 3) — Evidence: `fragment_length` truncation in Biopython.
5. **MUST Test:** `FindOrfs` returns START→STOP ORF with exact StartPosition/EndPosition/Frame/Protein for the ORF dataset — Evidence: getorf `-find 1` model.
6. **MUST Test:** No START codon ⇒ no ORF; ORF below `minLength` filtered — Evidence: getorf `-find 1` + `-minsize`.
7. **MUST Test:** Null input throws `ArgumentNullException`; empty sequence yields six empty frames / no ORFs — Evidence: implementation contract; truncation rule.
8. **SHOULD Test:** Both-strand ORF search finds reverse-strand ORFs with negative frame labels — Rationale: getorf searches both strands by default.
9. **COULD Test:** Alternative start codon (TTG/CTG) initiates an ORF and is translated by its actual residue — Rationale: NCBI table 1 lists TTG/CTG as starts.

---

## References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite — transeq application documentation. https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
2. Rice P, Longden I, Bleasby A. 2000. EMBOSS — getorf application documentation. https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html
3. Cock PJA et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11):1422–1423. Source: Bio/SeqUtils/__init__.py (six_frame_translations). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
4. NCBI. The Genetic Codes (transl_table=1, The Standard Code). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Wikipedia contributors. Reading frame. (cites Lodish 2007; Pierce 2012; Badger & Olsen 1999; Anderson et al. 1981). https://en.wikipedia.org/wiki/Reading_frame

---

## Change History

- **2026-06-13**: Initial documentation.
- **2026-10-14**: getorf `-find` 0–6 sources and EMBOSS 6.6.0 binary reference output (FIN-B02 WP11, audit A1-12).
- **2026-10-14**: EMBOSS ambiguous-codon (`ajtranslate.c`) and `-circular` (`getorf.c`) sources; getorf 6.6.0 IUPAC / circular reference output; Biopython IUPAC six-frame comparison (FIN-B02 WP12, audits A1-13, A1-20).
- **2026-10-14**: EMBOSS transeq frame conventions (`transeq.c`, `transeq.acd`, `ajTrnSeqFrameC`/`ajTrnSeqDangleS`) and transeq 6.6.0 reference output (FIN-B02 WP17, audit A2-10).
- **2026-10-14**: transeq `-regions`/`-methionine`/`-trim`/`-clean` (`transeq.c`, `ajRangeSeqExtract`, `ajTrnSeqInit*`) and transeq 6.6.0 reference output (FIN-B02 WP21, audit A3-3).
