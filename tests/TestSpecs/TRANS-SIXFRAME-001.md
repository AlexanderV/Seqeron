# Test Specification: TRANS-SIXFRAME-001

**Test Unit ID:** TRANS-SIXFRAME-001
**Area:** Translation
**Algorithm:** Six-Frame Translation and ORF finding
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

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

### 1.2 Key Evidence Points

1. A double-stranded sequence has exactly six reading frames: three forward (offsets 0,1,2) and three reverse — Wikipedia Reading frame; EMBOSS transeq `-frame 6`.
2. Forward frames: `frames[i+1] = translate(seq[i:])`, i = 0,1,2 — Biopython six_frame_translations.
3. Reverse frames (this repo's convention): `frames[-(i+1)] = translate(reverse_complement(seq)[i:])`, i = 0,1,2 — Biopython six_frame_translations; documented as the "alternative" convention in EMBOSS transeq.
4. Standard code (table 1): start codons TTG, CTG, ATG; stop codons TAA, TAG, TGA — NCBI Genetic Codes Starts line `---M------**--*----M---------------M----…`.
5. ORF under getorf `-find 1` = region from a START codon to a STOP codon; both strands searched by default — EMBOSS getorf.
6. Incomplete trailing codon is dropped (`fragment_length = 3*((len-i)//3)`) — Biopython.

### 1.3 Documented Corner Cases

- Trailing 1–2 nt that cannot form a codon are ignored (Biopython truncation).
- No START codon ⇒ no ORF emitted under START→STOP model (EMBOSS getorf).
- ORF may run to the end of the sequence without a STOP (incomplete ORF) (EMBOSS getorf).
- IUPAC-ambiguous codons are outside the 64-codon table (NCBI).

### 1.4 Known Failure Modes / Pitfalls

1. Reverse-frame numbering convention ambiguity — EMBOSS transeq documents two conventions; mixing them mislabels −1/−2/−3 (Source 1).
2. Off-by-one in inclusive stop-codon end position — EMBOSS getorf positions include the stop codon (Source 4).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `TranslateSixFrames(DnaSequence, GeneticCode?)` | Translator | Canonical | Six-frame translation; deep evidence-based testing |
| `FindOrfs(DnaSequence, GeneticCode?, int minLength, bool searchBothStrands)` | Translator | Canonical | ORF finding (START→STOP model) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `TranslateSixFrames` returns exactly 6 entries with keys {+1,+2,+3,−1,−2,−3} | Yes | EMBOSS transeq `-frame 6`; Wikipedia Reading frame |
| INV-2 | Frames +1/+2/+3 equal `Translate` of the input at offsets 0/1/2 | Yes | Biopython forward-frame loop |
| INV-3 | Frames −1/−2/−3 equal translation of the reverse complement at offsets 0/1/2 | Yes | Biopython reverse-frame loop |
| INV-4 | Each frame length = floor((effectiveLength)/3); trailing partial codon ignored | Yes | Biopython `fragment_length` |
| INV-5 | Every `FindOrfs` result starts at a START codon and (if terminated) ends at a STOP codon; EndPosition is the stop's last base (inclusive); Protein excludes the stop | Yes | EMBOSS getorf `-find 1` |
| INV-6 | `OrfResult.NucleotideLength = EndPosition − StartPosition + 1`; `AminoAcidLength = Protein.Length` | Yes | implementation contract; getorf inclusive positions |

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

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindOrfs_BothStrands_FindsReverseOrf | ORF only present on reverse strand, searchBothStrands true | ORF with negative Frame found | getorf searches both strands |
| S2 | FindOrfs_ForwardOnly_SkipsReverse | same input, searchBothStrands false | no reverse-frame ORF | searchBothStrands semantics |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FindOrfs_AltStartCodon_Initiates | ORF beginning with TTG, minLength 1 | ORF found; first residue = `L` (TTG) | NCBI table 1 lists TTG as start |
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
| C1 | ✅ Covered | `FindOrfs_AlternativeStartCodonTtg_InitiatesOrf` |
| C2 | ✅ Covered | `TranslateSixFrames_InternalStop_IsRenderedNotTerminated` |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Reverse-frame numbering follows the Biopython independent-offset convention (frame −k = revcomp offset k−1) | M3, M5, INV-3 |
| 2 | Stop = `*`, ambiguous IUPAC codon = `X` (inherited from `GeneticCode.Translate`) | C2 |
| 3 | `FindOrfs.minLength` counts amino acids (not nucleotides as in getorf) | M11 |

---

## 7. Open Questions / Decisions

1. Decision: the repository's reverse-frame numbering follows Biopython (the EMBOSS-documented "alternative"); EMBOSS's phase-locked default would relabel −1/−2/−3 differently. Documented in Evidence and algorithm doc §5.4 as an accepted convention, not a defect.
