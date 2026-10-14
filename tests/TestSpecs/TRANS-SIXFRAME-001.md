# Test Specification: TRANS-SIXFRAME-001

**Test Unit ID:** TRANS-SIXFRAME-001
**Area:** Translation
**Algorithm:** Six-Frame Translation and ORF finding
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-10-14 (FIN-B02 WP12: IUPAC string overloads; getorf `-circular`)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | EMBOSS transeq documentation (frame numbering) | 3 | https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html | 2026-06-13 |
| 2 | Biopython `Bio.SeqUtils.six_frame_translations` source | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py | 2026-06-13 |
| 3 | NCBI The Genetic Codes — Standard Code (table 1) | 2 | https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi | 2026-06-13 |
| 4 | EMBOSS getorf documentation (ORF definition) | 3 | https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html | 2026-06-13 |
| 5 | Wikipedia — Reading frame (cites Lodish 2007; Pierce 2012) | 4 | https://en.wikipedia.org/wiki/Reading_frame | 2026-06-13 |
| 6 | EMBOSS getorf source `getorf.c` + `getorf.acd` (opened 2026-09-28) | 3 | https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/getorf.c | 2026-09-28 |
| 7 | INSDC Feature Table Definition — CDS "location includes stop codon" (WebSearch snippet; site blocked) | 2 | https://www.insdc.org/submitting-standards/feature-table/ | 2026-09-28 |
| 9 | EMBOSS 6.6.0 `getorf` binary (`/usr/bin/getorf`, EGC.1–23 genetic-code files) — reference output for `-find` 0–6 | 3 | EMBOSS 6.6.0 | 2026-10-14 |
| 8 | Biopython 1.88 `Bio.Seq._translate_str` (`cds=True` → initiator M; `to_stop` refused for dual-coding tables) | 3 | installed package | 2026-09-28 |
| 10 | EMBOSS `ajax/core/ajtranslate.c` (`ajTrnCodonstrTypeC`, `trnComplete`, `ajTrnSeqOrig`/`ajTrnSeqDangleC`) and `getorf.c` circular branch | 3 | https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajtranslate.c | 2026-10-14 |

### 1.2 Key Evidence Points

1. A double-stranded sequence has exactly six reading frames: three forward (offsets 0,1,2) and three reverse — Wikipedia Reading frame; EMBOSS transeq `-frame 6`.
2. Forward frames: `frames[i+1] = translate(seq[i:])`, i = 0,1,2 — Biopython six_frame_translations.
3. Reverse frames (this repo's convention): `frames[-(i+1)] = translate(reverse_complement(seq)[i:])`, i = 0,1,2 — Biopython six_frame_translations; documented as the "alternative" convention in EMBOSS transeq.
4. Standard code (table 1): start codons TTG, CTG, ATG; stop codons TAA, TAG, TGA — NCBI Genetic Codes Starts line `---M------**--*----M---------------M----…`.
5. ORF under getorf `-find 1` = region from a START codon to a STOP codon; both strands searched by default — EMBOSS getorf.
7. getorf `-find` 0–6 (getorf.acd / getorf.c): 0/2 translation / nucleotides of STOP-free regions (every frame open at the sequence start; a STOP or the last complete codon writes the region); 1/3 translation / nucleotides START→STOP; 4/5/6 `-flanking` nucleotides each side of the START / initial STOP / ending STOP, no output when the window crosses a sequence end. `-minsize`/`-maxsize` are nucleotides, divided by 3 (integer) for the protein-storing modes 0, 1, 4 — Sources 6, 9.
6. Incomplete trailing codon is dropped (`fragment_length = 3*((len-i)//3)`) — Biopython.
8. IUPAC input (A1-13): Biopython translates any IUPAC string per frame (`Seq(s)[f:].translate(table)`, `reverse_complement()` with R↔Y, K↔M, B↔V, D↔H) using its ambiguous codon tables. getorf classifies a codon with `ajTrnCodonstrTypeC`: START iff the Starts table entry is `M` — filled for the 64 A/C/G/T codons only, so an ambiguous codon is never a START; STOP iff its precomputed residue is `*`. EMBOSS precomputes ambiguous residues in `trnComplete` (third-base wobble + one "most ambiguous" codon per residue built by OR-ing bases), which is inconsistent for a few codons: table 1/11 TRR → `*` (TGG = W), TRA → X (all expansions stops), YTN → L, MGN → R, ATS → I, WSN → S; table 2 WRR → `*` — Sources 8, 9, 10.
9. getorf `-circular` (A1-20): the sequence (length L) is tripled and its reverse complement likewise; every frame starts closed; a START at pos ≥ L is ignored; the first STOP / end codon at pos ≥ L ends the frame (`LASTORF`); a region starting past L is shifted back by L (`while(s > len)`); a frame with no STOP at all is written after the scan as the full-length record `WriteORF(f, seqlen-1)` translated by `ajTrnSeqOrig` (dangling bases: 2 → codon b1 b2 N, 1 → X; frame 3 tested with `>= maxsize`, a typo) — Sources 9, 10.

### 1.3 Documented Corner Cases

- Trailing 1–2 nt that cannot form a codon are ignored (Biopython truncation).
- No START codon ⇒ no ORF emitted under START→STOP model (EMBOSS getorf).
- ORF may run to the end of the sequence without a STOP (incomplete ORF) (EMBOSS getorf).
- IUPAC-ambiguous codons are outside the 64-codon table (NCBI); the string overloads resolve them per Biopython (residue) and getorf (START / STOP role), see 1.2 #8.
- Circular ORFs crossing the origin end past L − 1 (getorf "(ORF crosses the breakpoint)").

### 1.4 Known Failure Modes / Pitfalls

1. Reverse-frame numbering convention ambiguity — EMBOSS transeq documents two conventions; mixing them mislabels −1/−2/−3 (Source 1).
2. Off-by-one in inclusive stop-codon end position — getorf prints the range WITHOUT the stop (`WriteORF(start, pos-1)`, Source 6); this API includes it (INSDC CDS convention, Source 7), so getorf end = EndPosition − 3.
3. Open ORF with trailing partial codon — end must be the last complete codon (`WriteORF(start, pos+2)`, Source 6), not Length−1.
4. Alternative start codons — initiator reported as `M` (getorf `-methionine` default Y, Source 6; Biopython `cds=True`, Source 8).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `TranslateSixFrames(DnaSequence, GeneticCode?)` | Translator | Canonical | Six-frame translation; deep evidence-based testing |
| `FindOrfs(DnaSequence, GeneticCode?, int minLength, bool searchBothStrands)` | Translator | Canonical | ORF finding (START→STOP model) |
| `FindOrfs(DnaSequence, OrfFindMode, GeneticCode?, int minSize, int maxSize, int flanking, bool searchBothStrands, bool methionine, bool circular = false)` → `OrfFeature` | Translator | Canonical | getorf `-find` 0–6 (A1-12) and `-circular` (A1-20); `OrfFeature.ToGetOrfRange` = getorf `[s - e]` header |
| `TranslateSixFrames(string, GeneticCode?)` | Translator | Canonical | IUPAC DNA/RNA string, Biopython semantics (A1-13) |
| `FindOrfs(string, GeneticCode?, int minLength, bool searchBothStrands)` | Translator | Canonical | IUPAC string, START→STOP (A1-13) |
| `FindOrfs(string, OrfFindMode, GeneticCode?, int, int, int, bool, bool, bool circular = false)` → `OrfFeature` | Translator | Canonical | IUPAC string, getorf modes 0–6 / `-circular` (A1-13, A1-20) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `TranslateSixFrames` returns exactly 6 entries with keys {+1,+2,+3,−1,−2,−3} | Yes | EMBOSS transeq `-frame 6`; Wikipedia Reading frame |
| INV-2 | Frames +1/+2/+3 equal `Translate` of the input at offsets 0/1/2 | Yes | Biopython forward-frame loop |
| INV-3 | Frames −1/−2/−3 equal translation of the reverse complement at offsets 0/1/2 | Yes | Biopython reverse-frame loop |
| INV-4 | Each frame length = floor((effectiveLength)/3); trailing partial codon ignored | Yes | Biopython `fragment_length` |
| INV-5 | Every `FindOrfs` result starts at a START codon, Protein[0] = `M`; if terminated EndPosition is the stop's last base (inclusive); if open, the last base of the last complete codon; Protein excludes the stop | Yes | EMBOSS getorf `-find 1` (Source 6); INSDC (Source 7) |
| INV-7 | getorf-modes overload: the records (as `ToGetOrfRange`, strand, sequence), in order, equal the EMBOSS getorf output for the same `-find/-table/-minsize/-maxsize/-flanking/-reverse/-methionine`, except the documented one-codon open-ORF divergence (algorithm doc §5.4 #5) | Yes | Sources 6, 9 |
| INV-6 | `OrfResult.NucleotideLength = EndPosition − StartPosition + 1` = 3·(aa+1) terminated / 3·aa open; `AminoAcidLength = Protein.Length` | Yes | Sources 6, 7 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | SixFrames_KeysAndCount | 39-nt input returns 6 frames, keys ±1..±3 | Count = 6; keys = {1,2,3,−1,−2,−3} | INV-1; EMBOSS transeq |
| M2 | SixFrames_ForwardProteins | Forward frames of 39-nt dataset | +1=`MAIVMGR*KGAR*`, +2=`WPL*WAAERVPD`, +3=`GHCNGPLKGCPI` | Biopython fwd loop + NCBI table 1 |
| M3 | SixFrames_ReverseProteins | Reverse frames of 39-nt dataset | −1=`LSGTLSAAHYNGH`, −2=`YRAPFQRPITMA`, −3=`IGHPFSGPLQWP` | Biopython rev loop + NCBI table 1 |
| M4 | SixFrames_ForwardEqualsTranslateOffsets | +1/+2/+3 equal `Translate(dna, frame:0/1/2)` | Equal strings | INV-2 |
| M5 | SixFrames_ReverseEqualsRevCompOffsets | −1/−2/−3 equal `Translate(revComp, frame:0/1/2)` | Equal strings | INV-3 |
| M6 | SixFrames_PartialCodonIgnored | `ATGAAATAGGC` (11 nt) frame +1 ignores trailing 2 nt | +1=`MK*` (3 aa) | INV-4; Biopython truncation |
| M7 | SixFrames_NullInput_Throws | null DnaSequence | `ArgumentNullException` | implementation contract |
| M8 | SixFrames_EmptySequence_SixEmptyFrames | empty input | 6 frames, all empty strings | INV-1 + INV-4 |
| M9 | FindOrfs_ForwardStartToStop_Positions | `GGGATGAAACCCTAAGGG`, minLength 1, fwd only | one ORF: Start=3, End=14, Frame=1, Protein=`MKP` | INV-5; EMBOSS getorf `-find 1` |
| M10 | FindOrfs_NoStartCodon_Empty | sequence with no ATG/TTG/CTG, fwd only | no ORFs | EMBOSS getorf START→STOP |
| M11 | FindOrfs_BelowMinLength_Filtered | ORF protein shorter than minLength | filtered out | EMBOSS getorf `-minsize` |
| M11b | FindOrfs_OrfRunsToSequenceEndWithoutStop | `ATGAAACCCGGG`, minLength 1, fwd only | open ORF: Start=0, End=11, Frame=1, Protein=`MKPG` | INV-5; EMBOSS getorf incomplete-ORF (doc §6.1) |
| M12 | FindOrfs_NullInput_Throws | null DnaSequence | `ArgumentNullException` | implementation contract |
| M13 | FindOrfs_OrfResult_LengthDerivations | check Nucleotide/AminoAcid length of M9 ORF | NucleotideLength=12, AminoAcidLength=3 | INV-6 |
| M11c | FindOrfs_OpenOrfWithTrailingPartialCodon_EndsAtLastCompleteCodon | `ATGAAACCCGG`, fwd only | Start=0, End=8, NucleotideLength=9, `MKP` (getorf [1 - 9]) | Source 6 (review 2026-09 fix) |
| M11d | FindOrfs_ReverseStrandOpenOrf_EndsAtLastCompleteCodon | `CCGGGTTTCAT`, both strands | Frame −1, Start 0, End 8, `MKP` | Source 6 |
| M14 | FindOrfs_Table12_MatchesGetorfReference | 35-nt, table 12, minLength 1 | (−2,16,33,`MHNMKI`), (2,13,33,`MQKIYFA`) | getorf port (Source 6) + Biopython table 12 |
| M15 | FindOrfs_DualCodingStopTable_Throws | tables 27, 28, 31 | `ArgumentException` | Source 8 (to_stop refusal); all stops dual-coding |
| M16 | FindOrfs_Find0_WorkedExample_MatchesGetorf | `GGGATGAAACCCTAAGGG`, `-find 0 -minsize 9` | 6 regions, e.g. (0,14,+1,terminated,`GMKP`) = [1 - 12]; (0,17,−1,`PLGFHP`) = [18 - 1] (REVERSE SENSE) | Source 9 |
| M17 | FindOrfs_Find2_WorkedExample_MatchesGetorf | same, `-find 2` | `GGGATGAAACCC`, `AACCCTAAG`, … `CCCTTAGGGTTTCATCCC` | Source 9 |
| M18 | FindOrfs_Find3_WorkedExample_MatchesGetorf | same, `-find 3` | (3,14,+1,`ATGAAACCC`) = [4 - 12] | Source 9 |
| M19 | FindOrfs_FlankingModes_WorkedExample_MatchGetorf | same, `-find 4/5/6 -flanking 2` | 4: [2 - 8] `GGATGAA`; 5: [3 - 9] `GATGAAA`, [16 - 10] rev `CTTAGGG`; 6: [11 - 17] `CCTAAGG` | Source 9 |
| M20 | FindOrfs_AllModes_MatchGetorfBinaryOracle (65 cases) | 30 sequences × 65 getorf runs (modes 0–6 × tables 1/2/11 × {defaults; minsize 0, flanking 0; minsize 10, maxsize 60, flanking 5, reverse N} + `-methionine N` + table 3) | equal to `TestData/Getorf/getorf_find_modes.tsv` | INV-7; Source 9 |
| M21 | FindOrfs_Find1_EqualsLegacyStartToStopOverload | oracle sequences, tables 1/2/11 | `-find 1 minSize 6` ORF set = legacy `minLength 2` | INV-5, INV-7 |
| M22 | FindOrfs_Find3_StartAsLastCodon_ReportedOneCodonOpenOrf | `ATG`, `-find 3 -minsize 0` | (0,2,+1,open,`ATG`); getorf writes nothing (oracle c28/s4) | algorithm doc §5.4 #5 |
| M23 | FindOrfs_Find0_MinSizeZero_ReportsEmptyRegions | `TAATAA`, `-find 0 -minsize 0 -reverse N` | [1 - 0] ``, [2 - 4] `N`, [3 - 5] `I`, [4 - 3] `` | Source 9 |
| M24 | FindOrfs_FlankingWindowCrossingSequenceEnd_IsNotReported | worked example, flanking 4 / 6 | no record | getorf.acd `-flanking` help (Source 6) |
| M25 | FindOrfs_Find1_MethionineOff_InitiatorKeepsItsResidue | `TTGAAAGGGTAA`, `-methionine N` | `LKG` (Y: `MKG`) | Sources 6, 9 |
| M26 | FindOrfs_GetOrfModes_InvalidArguments_Throw | null, mode 7, negative sizes/flanking, table 27 | `ArgumentNullException` / `ArgumentOutOfRangeException` / `ArgumentException` | contract |
| M27 | TranslateSixFrames_IupacString_MatchesBiopython | `ATGNNNTARGCNRAYTAAYTG` (t1), `TTGTANSARGAYMGR` (t11) | t1: `MX*AB*L`, `XXXXXN`, `XXXXLX`, `QLXXLXH`, `X*XAXX`, `XXXXXX` | Source 8 (Biopython 1.88 output) |
| M28 | TranslateSixFrames_StringOverload_CaseAndRnaInsensitive_EqualsDnaOverload | lower-case RNA spelling; A/C/G/T string, tables 1/2/11 | same frames as upper-case DNA / `DnaSequence` overload | 1.2 #8 |
| M29 | FindOrfs_IupacString_Find1_MatchesGetorf | `ATGNNNCCCTAA`, `atgaaacccrtaa`, `AUGAAACCCUAAGG`, `-find 1 -minsize 3 -reverse N` | [1 - 9] `MXP`; open [1 - 12] `MKPX`; [1 - 9] `MKP` | Source 9 |
| M30 | FindOrfs_IupacString_AmbiguousStartCodon_DoesNotOpenOrf | `YTGAAACCCTAA` | no ORF (getorf writes nothing); `TTG…` gives `MKP` | Sources 9, 10 |
| M31 | FindOrfs_IupacString_AllStopAmbiguousCodon_TerminatesOrf | `ATGAAACCCTRA` | (0,11,terminated,`MKP`) = [1 - 9]; getorf `[1 - 12] MKPX` (EMBOSS TRA → X; documented difference) | 1.2 #8 |
| M32 | FindOrfs_LegacyStringOverload_EqualsDnaOverload_AndReadsIupac | A/C/G/T string (lower case), tables 1/2/11; `ATGNNNCCCTAA` | equal to the `DnaSequence` overload; (0,11,+1,`MXP`) | 1.2 #8 |
| M33 | FindOrfs_Circular_Find1_OrfCrossesOrigin | `CCCTAAGGGATGAAA`, `-find 1/3 -minsize 3 -circular Y -reverse N` | (9,20,+1,`MKP`) = [10 - 18]; `-find 3` `ATGAAACCC`; linear [10 - 15] `MK` | Source 9 |
| M34 | FindOrfs_Circular_StopModes_MatchGetorf | same, `-find 0`, `-find 6 -flanking 2` | [7 - 18] `GMKP`, [14 - 25] `NPKG` (+ frame-3 full-length record [3 - 45], getorf omits it: typo); [17 - 23] `CCTAAGG`, [24 - 30] `GATGAAA` | Sources 9, 10 |
| M35 | FindOrfs_Circular_StopFreeFrame_FullLengthRecord | `GCGCGCCCGGG`, `-find 0/2 -minsize 3 -circular Y` | [1 - 33] `ARPGRARGAPG`, [2 - 33] `RARGAPGARPG`, [3 - 33] `APGARPGRARX`; mode 2: the frames' nucleotides | Sources 9, 10 |
| M36 | FindOrfs_Circular_DnaAndStringOverloadsAgree | A/C/G/T sequence, modes 0–6 | string = `DnaSequence` overload, circular and linear | contract |
| M37 | FindOrfs_IupacAndCircular_MatchGetorfBinaryOracle (86 cases) | 24 IUPAC sequences × 86 getorf runs (modes 0–6 × tables 1/2/11 × {defaults + flanking 20; minsize 9, maxsize 90, flanking 5, reverse N} × linear / `-circular Y`, + `-methionine N`) | equal to `TestData/Getorf/getorf_iupac_circular.tsv` (residues: equal or getorf X) | Sources 9, 10 |
| M38 | IupacCircularOracle_CoversAmbiguityAndOriginCrossing | the oracle | ≥ 20 ambiguous sequences, lower case, U, 43 circular runs, > 100 origin-crossing records | — |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindOrfs_BothStrands_FindsReverseOrf | ORF only present on reverse strand, searchBothStrands true | ORF with negative Frame found | getorf searches both strands |
| S2 | FindOrfs_ForwardOnly_SkipsReverse | same input, searchBothStrands false | no reverse-frame ORF | searchBothStrands semantics |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FindOrfs_AltStartCodon_Initiates | ORF beginning with TTG, minLength 1 | ORF found; Protein = `MKG` (initiator → M; corrected 2026-09 from `LKG`) | Sources 6, 8; NCBI |
| C1b | FindOrfs_AltStartCodonCtg | `CTGCTGTAA` | `ML` (initiator M, internal CTG = L) | Biopython `cds=True` |
| C2 | SixFrames_RunningStop_NoTermination | TranslateSixFrames does not stop at internal stop (renders `*`) | `*` present mid-protein | TranslateSixFrames uses toFirstStop=false |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/TranslatorTests.cs` belongs to a different unit (TRANS-PROT-001). It contains weak six-frame/ORF tests (count/existence checks, no exact reverse-frame values, permissive `Has.Count`). Those tests are out of scope for this unit's file and are left in place under TRANS-PROT-001.
- No `Translator_SixFrames_Tests.cs` exists. This is the canonical file for TRANS-SIXFRAME-001.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M13, S1–S2, C1–C2 | ❌ Missing | New canonical file; none of the exact-value cases exist for this unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/Translator_SixFrames_Tests.cs` — all TRANS-SIXFRAME-001 cases.
- **Remove:** nothing. Existing `TranslatorTests.cs` six-frame/ORF tests stay under TRANS-PROT-001 (different unit; out of scope to modify).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `Translator_SixFrames_Tests.cs` | Canonical for TRANS-SIXFRAME-001 | 17 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | S1 | ❌ Missing | Implemented | ✅ Done |
| 15 | S2 | ❌ Missing | Implemented | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented | ✅ Done |
| 17 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `TranslateSixFrames_Returns_SixFramesKeyedPlusMinus` |
| M2 | ✅ Covered | `TranslateSixFrames_ForwardFrames_MatchEvidenceProteins` |
| M3 | ✅ Covered | `TranslateSixFrames_ReverseFrames_MatchEvidenceProteins` |
| M4 | ✅ Covered | `TranslateSixFrames_ForwardFrames_EqualTranslateAtOffsets` |
| M5 | ✅ Covered | `TranslateSixFrames_ReverseFrames_EqualReverseComplementOffsets` |
| M6 | ✅ Covered | `TranslateSixFrames_PartialTrailingCodon_IsIgnored` |
| M7 | ✅ Covered | `TranslateSixFrames_NullInput_ThrowsArgumentNullException` |
| M8 | ✅ Covered | `TranslateSixFrames_EmptySequence_ReturnsSixEmptyFrames` |
| M9 | ✅ Covered | `FindOrfs_ForwardStartToStop_ReturnsExactPositionsAndProtein` |
| M10 | ✅ Covered | `FindOrfs_NoStartCodon_ReturnsEmpty` |
| M11 | ✅ Covered | `FindOrfs_OrfBelowMinLength_IsFiltered` |
| M12 | ✅ Covered | `FindOrfs_NullInput_ThrowsArgumentNullException` |
| M13 | ✅ Covered | `FindOrfs_OrfResult_LengthDerivations_AreCorrect` |
| S1 | ✅ Covered | `FindOrfs_BothStrands_FindsReverseStrandOrf` |
| S2 | ✅ Covered | `FindOrfs_ForwardOnly_DoesNotReturnReverseStrandOrf` |
| C1 | ✅ Covered | `FindOrfs_AlternativeStartCodonTtg_InitiatesOrfWithMethionine` |
| C1b, M11c, M11d, M14, M15 | ✅ Covered | added in the 2026-09 review (see §7) |
| M16–M26 | ✅ Covered | `Translator_GetOrfModes_Tests.cs` (FIN-B02 WP11, A1-12) |
| M27–M38 | ✅ Covered | `Translator_IupacCircular_Tests.cs` (FIN-B02 WP12, A1-13, A1-20) |
| C2 | ✅ Covered | `TranslateSixFrames_InternalStop_IsRenderedNotTerminated` |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Reverse-frame numbering follows the Biopython independent-offset convention (frame −k = revcomp offset k−1) | M3, M5, INV-3 |
| 2 | Stop = `*`, ambiguous IUPAC codon resolved per Biopython (e.g. `GCN`→`A`, `NNN`→`X`; inherited from `GeneticCode.Translate`) | C2 |
| 3 | `FindOrfs.minLength` counts amino acids (not nucleotides as in getorf); the getorf-modes overload's `minSize` is getorf's nucleotide value | M11, M20 |

---

## 7. Open Questions / Decisions

1. Decision: the repository's reverse-frame numbering follows Biopython (the EMBOSS-documented "alternative"); EMBOSS's phase-locked default would relabel −1/−2/−3 differently. Documented in Evidence and algorithm doc §5.4 as an accepted convention, not a defect.
2. Review 2026-09 (B02): cross-checked against a Python port of EMBOSS `getorf_FindORFs` (`-find 1`, `-methionine` Y) and Biopython six-frame translation on 30 000 random sequences (0–300 nt, 26 NCBI tables): six-frame 0 mismatches; ORFs 0 mismatches over 153 483 reference ORFs after the fixes (initiator → M; open-ORF end at last complete codon; dual-coding tables rejected). The only remaining divergence is getorf never reporting a 1-codon open ORF whose START is the frame's last complete codon (getorf branch-order artefact; only visible with `minLength ≤ 1`) — accepted, algorithm doc §5.4 #5.
3. FIN-B02 WP11 (2026-10-14, audit A1-12): getorf `-find` 0, 2, 3, 4, 5, 6 added as `FindOrfs(DnaSequence, OrfFindMode, …)` (port of `getorf_FindORFs`/`getorf_WriteORF`/`getorf_AppORF`, linear only — `-circular` not supported). Verified against the EMBOSS 6.6.0 `getorf` binary on 2 500 random sequences (1–400 nt; uniform / GC-rich / AT-rich) × 290 parameter sets (modes 0–6 × tables 1, 2, 4, 11, 12 [+ 3, 13 for the START-independent modes 0/2/5/6, whose EMBOSS start sets predate NCBI's] × 5 minsize/maxsize/flanking/reverse sets, + `-methionine N` for 1/3/4): 0 mismatching (sequence, run) pairs over 3 474 722 reference records (modes 0/2/3/4/5/6: 723 069 / 681 150 / 469 100 / 409 574 / 348 228 / 351 943); 25 084 one-codon open-ORF records excluded as the documented divergence (§5.4 #5 of the algorithm doc). A 30-sequence × 65-run subset is locked in `TestData/Getorf/getorf_find_modes.tsv` (M20).
4. FIN-B02 WP12 (2026-10-14, audits A1-13, A1-20): IUPAC string overloads `TranslateSixFrames(string, …)`, `FindOrfs(string, …)` (legacy and getorf-modes) and getorf `-circular` (`circular: true` on the getorf-modes overloads). Six frames vs Biopython 1.88 (`Seq(s)[f:].translate(t)` / `reverse_complement()`): 20 000 random IUPAC strings (0–200 nt, 0–100 % ambiguity codes, lower case / U) over all 27 Biopython tables, 0 mismatching frames. getorf string overload vs the EMBOSS 6.6.0 binary, 1 000 random IUPAC sequences × 93 runs (modes 0–6 × tables 1/2/11 × 4 parameter sets, + `-methionine N`): linear 384 927 reference records, circular 344 889 — 0 unexplained mismatches; (sequence, run) pairs that differ only in residues where getorf prints X (or one of EMBOSS's inconsistent wobble/OR entries): 4 653 linear / 6 516 circular of 93 000; pairs whose sequence contains a codon with a different STOP status in EMBOSS's table (TRA/TRR, WRR): 2 295 / 2 694. `-circular` on A/C/G/T: 1 000 random sequences × 93 runs, 447 526 reference records, 0 mismatches once getorf's full-length records (STOP-free frame) are taken literally; the two documented divergences there are the frame-3 `>= maxsize` typo and mode 2's protein output (algorithm doc §5.4 #7). A 24-sequence × 86-run subset is locked in `TestData/Getorf/getorf_iupac_circular.tsv` (M37).
