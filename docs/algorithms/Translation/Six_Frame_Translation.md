# Six-Frame Translation and ORF Finding

| Field | Value |
|-------|-------|
| Algorithm Group | Translation |
| Test Unit ID | TRANS-SIXFRAME-001 |
| Related Projects | Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-10-14 (WP12) |

## 1. Overview

Six-frame translation produces the six possible protein readings of a double-stranded
nucleotide sequence: three forward frames (reading offsets 0, 1, 2 of the given strand)
and three reverse frames (offsets 0, 1, 2 of its reverse complement) [1][5]. Open
Reading Frame (ORF) finding then locates coding-candidate regions that begin at a START
codon and end at a STOP codon [4]. Both operations are exact and specification-driven:
they apply a genetic-code table (NCBI translation tables) deterministically, with no
heuristic or probabilistic step. They are used in gene prediction, unannotated-sequence
exploration, and proteogenomics.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A codon is a triplet of nucleotides read 5'→3'; the genetic code maps each of the 64
codons to one amino acid or a stop signal [3]. Because reading can start at any of three
offsets and the antiparallel complementary strand provides three more, a duplex sequence
has exactly six reading frames [5].

### 2.2 Core Model

Let `s` be the input sequence (T→U normalised) and `G` the genetic-code map.
Forward frame `f ∈ {1,2,3}` is the concatenation
`G(s[o]s[o+1]s[o+2]) · G(s[o+3]…) · …` with offset `o = f-1`, consuming only complete
codons [2]. Reverse frame `-f` is the same construction applied to `revcomp(s)` with
offset `o = f-1` [2]; this is Biopython's `frames[-(i+1)] = translate(anti[i:])`
convention (see §5.4). An ORF (getorf `-find 1`) is a maximal region beginning at a
START codon and ending at the next in-frame STOP codon, both selected from the active
genetic code [4]. The initiator residue is reported as Met whatever the start codon
(getorf `-methionine`, default Y [6]; Biopython `translate(cds=True)` [2]; NCBI [3]). The Standard code (NCBI table 1) has start codons TTG, CTG, ATG and
stop codons TAA, TAG, TGA [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TranslateSixFrames` yields exactly 6 frames keyed {+1,+2,+3,−1,−2,−3} | three forward + three reverse frames of a duplex [1][5] |
| INV-02 | Forward frame `+f` = translation at offset `f−1` of the input | direct codon reading [2] |
| INV-03 | Reverse frame `−f` = translation at offset `f−1` of the reverse complement | Biopython reverse-frame loop [2] |
| INV-04 | Each frame length = ⌊(len−offset)/3⌋; trailing partial codon ignored | only complete codons consumed [2] |
| INV-05 | Every ORF starts at a START codon and its protein starts with `M`; if terminated, EndPosition is the STOP's last base (inclusive, INSDC CDS convention [7]) and the protein excludes the STOP; if open, EndPosition is the last base of the last complete codon | getorf START→STOP model [4][6]; INSDC [7] |
| INV-06 | `NucleotideLength = EndPosition − StartPosition + 1` = 3·(AminoAcidLength + 1) if terminated, 3·AminoAcidLength if open; `AminoAcidLength = Protein.Length` | inclusive coordinates [6][7] |

### 2.5 Comparison with Related Methods

| Aspect | This module (Biopython convention) | EMBOSS transeq default |
|--------|-----------------------------------|------------------------|
| Reverse-frame label | frame −1 = reverse complement at offset 0 | frame −1 = phase-locked to forward frame 1 |
| Trailing 1–2 bases | dropped | translated (`b1 b2 N` → residue or X; one base → X) |
| Selectable here | `SixFrameConvention.Biopython` (default) | `SixFrameConvention.EmbossTranseq`; `-alternative Y` = `EmbossTranseqAlternative` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| dna | DnaSequence | required | Sequence to translate / scan | non-null |
| geneticCode | GeneticCode? | Standard (table 1) | Codon→amino-acid table | NCBI translation table; `FindOrfs` rejects the dual-coding-stop tables 27, 28, 31 (`ArgumentException`) |
| minLength | int | 100 | Minimum ORF length, in **amino acids** | ≥ 0 |
| searchBothStrands | bool | true | Also scan the reverse complement for ORFs | — |
| convention | SixFrameConvention | `Biopython` (overloads without it) | `TranslateSixFrames(…, geneticCode, convention)`: `Biopython`, `EmbossTranseq` (transeq default, phase-locked reverse frames), `EmbossTranseqAlternative` (transeq `-alternative Y`) | defined enum value (`ArgumentOutOfRangeException`) |

IUPAC string overloads (A1-13): `TranslateSixFrames(string, GeneticCode?)`, `FindOrfs(string, GeneticCode?, int minLength, bool)` and `FindOrfs(string, OrfFindMode, …)` take the same parameters with `sequence` = case-insensitive IUPAC DNA or RNA (A C G T U R Y S W K M B D H V N; U read as T); null → `ArgumentNullException`, any other symbol → `ArgumentException`.

getorf-modes overload `FindOrfs(DnaSequence | string, OrfFindMode mode, GeneticCode?, int minSize = 30, int maxSize = 1000000, int flanking = 100, bool searchBothStrands = true, bool methionine = true, bool circular = false)` — parameters mirror getorf.acd [6]:

| Name | getorf | Default | Description |
|------|--------|---------|-------------|
| mode | `-find` | required | `OrfFindMode` 0–6 (values equal getorf's codes): 0 translation STOP→STOP, 1 translation START→STOP, 2 nucleotides STOP→STOP, 3 nucleotides START→STOP, 4 flanking START, 5 flanking initial STOP, 6 flanking ending STOP |
| minSize / maxSize | `-minsize` / `-maxsize` | 30 / 1 000 000 | **Nucleotides**, ≥ 0. Modes 0, 1, 4 compare `size / 3` (integer division) with the residue count; modes 2, 3, 5, 6 compare with the nucleotide count (STOP excluded). Flanking modes filter on the underlying ORF |
| flanking | `-flanking` | 100 | Bases reported each side of the codon (modes 4–6), ≥ 0; a window crossing a sequence end is not reported |
| searchBothStrands | `-reverse` | true | Also scan the reverse complement |
| methionine | `-methionine` | true | Mode 1: initiator residue reported as `M` |
| circular | `-circular` | false | Circular sequence: ORFs may cross the origin (§4.1 step 7) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (six frames) | IReadOnlyDictionary<int, ProteinSequence> | Keys +1,+2,+3,−1,−2,−3 → protein per frame |
| OrfResult.StartPosition | int | 0-based first base of the START codon (in the scanned strand's coordinates) |
| OrfResult.EndPosition | int | 0-based, inclusive: last base of the STOP codon (INSDC CDS convention [7]; getorf prints this − 3); for an open ORF the last base of the last complete codon |
| OrfResult.Frame | int | +1..+3 forward, −1..−3 reverse |
| OrfResult.Protein | ProteinSequence | Translated residues; initiator reported as `M` (getorf `-methionine` [6]); STOP excluded |
| OrfFeature.Mode | OrfFindMode | getorf `-find` mode that produced the record |
| OrfFeature.StartPosition / EndPosition | int | 0-based inclusive, scanned-strand coordinates. Modes 0–3: the ORF, End = last base of the terminating STOP (as `OrfResult`) or of the last complete codon when open. Modes 4–6: the flank window `[codon − flanking, codon + 2 + flanking]`. Circular: Start ∈ [0, L); End ≥ L when the region crosses the origin (base i is base i mod L) |
| OrfFeature.Frame / IsReverseStrand | int / bool | +1..+3 / −1..−3: (offset of the ORF's first base mod 3) + 1, reverse complement negative |
| OrfFeature.IsTerminated | bool | Underlying ORF ends with a STOP codon |
| OrfFeature.CodonPosition | int | Modes 4–6: 0-based first base of the flanked START / STOP codon; −1 otherwise |
| OrfFeature.Sequence | string | getorf output sequence: protein (0, 1), ORF nucleotides (2, 3), window nucleotides (4–6); STOP never included in 0–3 |
| OrfFeature.ToGetOrfRange(L) | (int, int) | getorf's 1-based `[s - e]` header: `s = Start + 1`, `e = End + 1` (ORF modes: `End − 3 + 1` when terminated); reverse strand `s = L − s + 1`, `e = L − e + 1` with getorf's range shifts |

`OrfFeature` records are yielded in getorf output order (forward strand, then reverse; within a
strand in the order ORFs close), so the k-th record is getorf's `name_k`.

### 3.3 Preconditions and Validation

Null `dna` throws `ArgumentNullException`. Input is upper-cased and T→U normalised
before codon lookup; indexing is 0-based; ORF EndPosition is inclusive. An empty
sequence yields six empty frames and no ORFs. IUPAC-ambiguous codons are resolved as in
Biopython (inherited from `GeneticCode.Translate`: e.g. `GCN`→`A`, `TAR`→`*`, `RAY`→`B`, `NNN`/`TAN`→`X`).
`DnaSequence` input is A/C/G/T only; ambiguous genomes go through the string overloads. In ORF finding an
ambiguous codon is a STOP only when every expansion is a stop (`GeneticCode.Translate` = `*`) and is never
a START (EMBOSS `ajTrnCodonstrTypeC` reads a Starts table filled for the 64 A/C/G/T codons only [8]).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; resolve the genetic code (default Standard).
2. Compute the reverse complement once.
3. For offsets 0,1,2: translate the forward strand (→ +1,+2,+3) and the reverse
   complement (→ −1,−2,−3), consuming only complete codons.
4. For ORFs: in each frame, on entering a START codon begin accumulating residues;
   (the first residue is written as `M`); on the next in-frame STOP, emit an ORF if its
   protein length ≥ `minLength`; if the strand ends mid-ORF, emit the open ORF (ending at
   the last complete codon) if long enough. Starts inside an open ORF are not reported
   separately (getorf `ORF[frame]` flag [6]).
5. If `searchBothStrands`, repeat ORF scanning on the reverse complement with negative
   frame labels.
6. getorf-modes overload (`getorf_FindORFs` port, single pass over every codon position,
   frame = pos mod 3): STOP→STOP modes (0, 2, 5, 6) start with every frame open at offset
   0/1/2; a STOP, or the last complete codon of a frame, writes the current region (the last
   codon is appended first if it is not a STOP) and reopens the frame after it. START→STOP
   modes (1, 3, 4) behave as step 4. Each written region passes the size filter, then is
   output as protein / nucleotides (`getorf_AppORF`) or as the flank window around the START
   (4), the STOP before the region (5; the leading region has none) or the terminating STOP
   (6; open regions have none) — `getorf_WriteORF`.
7. `circular` (getorf `-circular` [6]): scan the sequence tripled (and its reverse complement
   tripled) with every frame closed at the start; a START at position ≥ L is ignored; the
   first STOP / end codon at ≥ L closes the frame for good (`LASTORF`); a region starting at
   ≥ L is shifted back by L (`while(s > len)`), flank windows are cut from the tripled
   sequence. In modes 0/2 a frame with no STOP at all is written after the scan as one open
   record from offset f to 3L − 1 (mode 0: translated with EMBOSS's dangling-base rule —
   2 bases b1 b2 as b1 b2 N, 1 base as X [8]).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Standard genetic code (NCBI table 1): start codons {TTG, CTG, ATG}; stop codons
{TAA, TAG, TGA} [3]. Codon length = 3; three reading frames per strand [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| TranslateSixFrames | O(n) | O(n) | n = sequence length; one revcomp + six linear passes |
| FindOrfs | O(n) | O(n) | linear scan per frame on each searched strand |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [Translator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs)

- `Translator.TranslateSixFrames(DnaSequence, GeneticCode?)`: six-frame translation.
- `Translator.FindOrfs(DnaSequence, GeneticCode?, int minLength, bool searchBothStrands)`: ORF finding (getorf `-find 1`).
- `Translator.FindOrfs(DnaSequence, OrfFindMode, GeneticCode?, int minSize, int maxSize, int flanking, bool searchBothStrands, bool methionine, bool circular = false)`: all getorf `-find` modes 0–6, linear or `-circular` → `OrfFeature`.
- `Translator.TranslateSixFrames(string, GeneticCode?)`, `Translator.FindOrfs(string, …)` (both signatures): IUPAC DNA/RNA strings; reverse complement via the canonical `DnaSequence.GetReverseComplementString`.
- `OrfResult` (record struct): ORF coordinates, frame, protein, derived lengths.

### 5.2 Current Behavior

Translation uses a direct O(n) codon scan; the repository suffix tree is not applicable
(this is sequential triplet decoding, not substring search). Frames are computed in a
single offset loop with one shared reverse complement. `TranslateSixFrames` renders
internal stop codons as `*` (it does not terminate early), unlike `Translate(...,
toFirstStop: true)`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Six frames keyed +1..+3 / −1..−3, forward offsets 0/1/2 and reverse-complement offsets 0/1/2 (Biopython six_frame_translations) [2].
- Trailing partial codon ignored (Biopython `fragment_length` truncation) [2].
- Standard code start/stop codons {TTG,CTG,ATG} / {TAA,TAG,TGA} (NCBI table 1) [3].
- ORF = START→STOP region (EMBOSS getorf `-find 1`, `getorf_FindORFs`) [4][6]; initiator → `M` (`-methionine` default Y) [6]; open ORF ends at last complete codon (`WriteORF(start, pos+2)`) [6].
- End coordinate includes the STOP codon (INSDC feature table CDS: "location includes stop codon") [7]; getorf's printed end excludes it (`WriteORF(start, pos-1)`) — a coordinate-convention difference of exactly 3.
- getorf `-find` 0–6 with `-minsize`/`-maxsize` (nucleotides; ÷3 for modes 0/1/4), `-flanking`, `-reverse`, `-methionine` (getorf.c / getorf.acd) [6]: 0 mismatches against the EMBOSS 6.6.0 getorf binary over 2 500 random sequences × 290 parameter sets (3.47 M reference records), apart from divergence #5.
- IUPAC strings: six frames = Biopython `Seq(s)[f:].translate(t)` / `reverse_complement()` [2]: 0 mismatching frames over 20 000 random IUPAC strings × all 27 tables. ORFs on IUPAC strings vs the EMBOSS 6.6.0 getorf binary: 0 unexplained mismatches over 1 000 sequences × 93 runs, linear (384 927 reference records) and circular (344 889) — differences limited to #7.
- getorf `-circular` [6]: 0 mismatches against the binary over 1 000 random A/C/G/T sequences × 93 runs (447 526 reference records), apart from #8.
- Tables 27/28/31: every stop codon is dual-coding (Biopython `CodonTable`: `stop_codons ⊂ forward_table`), so `FindOrfs` throws `ArgumentException` — as Biopython refuses `to_stop=True` for these tables [2].

**Intentionally simplified:**

- Legacy-overload `minLength` is measured in amino acids; **consequence:** users converting from EMBOSS getorf `-minsize` (nucleotides) must divide by 3 — or use the getorf-modes overload, whose `minSize` is getorf's nucleotide value.

- EMBOSS transeq frame conventions (`SixFrameConvention.EmbossTranseq` / `EmbossTranseqAlternative`; transeq.c, transeq.acd `-alternative`, ajtranslate.c `ajTrnSeqFrameC`/`ajTrnSeqRevC`/`ajTrnSeqAltRevC`/`ajTrnSeqDangleS`) [1][8]: default frame −k = reverse complement of the frame +k codons ("frame -1 is the reverse-complement of the set of codons used in frame 1"), the k − 1 leading forward bases forming its dangling end; `-alternative Y` frame −k starts at the k-th base of the reverse complement (= Biopython phasing). Dangling bases are translated in every frame (two bases `b1 b2` as codon `b1 b2 N`, one base as X); stops print `*` (transeq `-clean` / `-trim` are post-processing: replace `*` by X / strip trailing X and `*`). 0 mismatching frames against the EMBOSS 6.6.0 transeq binary over 2 000 random A/C/G/T(U) sequences (1–150 nt, lower case / U) × tables 1/2/11 × both conventions (72 000 frames); on 600 IUPAC sequences (21 600 frames) every residue equals the per-codon table of its tool, so the frames differ only at codons where transeq's ambiguity table differs (#9).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default reverse-frame numbering = Biopython independent-offset (not EMBOSS phase-locked) | Assumption | −1/−2/−3 labels differ from EMBOSS default unless `SixFrameConvention.EmbossTranseq` is passed | accepted | Documented alternative in EMBOSS transeq [1]; see INV-03 |
| 2 | `minLength` in amino acids, not nucleotides | Deviation | Different parameter unit than getorf | accepted | §3.1 |
| 3 | End coordinate includes the STOP (INSDC) | Convention | getorf's printed end = EndPosition − 3 | accepted | [7] |
| 4 | Reverse-strand ORF coordinates in reverse-complement coordinates | Convention | getorf prints forward-strand coordinates (`len−s+1`) | accepted | §6.2 |
| 5 | A START in the very last complete codon of a frame (no STOP) is reported as a 1-aa open ORF when `minLength ≤ 1` (getorf-modes overload, modes 1/3/4: when it passes `minSize` and, mode 4, the window fits) | Divergence | getorf never writes it (branch-order artefact of `getorf_FindORFs`: the START branch opens the ORF and the loop ends) | accepted | unobservable for `minLength ≥ 2` / `minSize ≥ 6` (mode 3: ≥ 4) (getorf default minsize 30 nt = 10 aa) |
| 6 | Nucleotide output of the string overloads is upper-case DNA | Convention | getorf keeps the input case (and writes U as T) | accepted | §3.1 |
| 7 | Ambiguous codons in ORF finding resolved per codon by the all-expansions rule (STOP iff all expansions stop; residue = Biopython) | Divergence | EMBOSS's precomputed ambiguous table (`trnComplete`) gives X for most ambiguous codons and is inconsistent for a few: table 1/11 `TRR` → STOP (TGG = W), `TRA` (all stops) not a STOP, `YTN` → L, `MGN` → R, `ATS` → I, `WSN` → S; table 2 `WRR` → STOP. ORFs differ from getorf only around such codons; residues differ where getorf prints X | accepted | [8]; Evidence |
| 9 | EMBOSS conventions keep this library's IUPAC codon resolution (Biopython tables) | Divergence | transeq's ambiguity table never emits B/Z/J (RAY → X, not B) and resolves a few codons differently (table 1: `ATS` → I, `ATW` → X, `TRA` → X, `TRR` → `*`, `YTN` → L, `MGN` → R; 28 of 3 375 IUPAC codons in table 1, 15 in table 2, 28 in table 11); residues differ only there | accepted | [8]; Evidence; as #7 |
| 8 | Circular full-length record (STOP-free frame, modes 0/2): mode 2 reports the frame's nucleotides and filters on nucleotide length; frame 3 uses `≤ maxSize` | Divergence | getorf.c prints the protein in mode 2 (comparing its residue count with nucleotide bounds) and tests frame 3 with `>= maxsize` (typo), so it normally omits that record | accepted | [6] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `dna` | `ArgumentNullException` | contract |
| Empty sequence | six empty frames; no ORFs | no complete codons [2] |
| Length not multiple of 3 | trailing 1–2 nt ignored | Biopython truncation [2] |
| No START codon | no ORF | getorf START→STOP [4] |
| ORF runs to sequence end without STOP | open ORF emitted if ≥ minLength; ends at last complete codon | getorf `WriteORF(start, pos+2)` [6] |
| Alternative start codon (TTG, CTG, GTG…) | initiator reported as `M` | getorf `-methionine` [6]; Biopython `cds=True` [2] |
| Genetic code 27, 28 or 31 in `FindOrfs` | `ArgumentException` | no unambiguous stop codon; Biopython `to_stop` refusal [2] |

### 6.2 Limitations

Reverse-frame labels follow the Biopython convention by default; pass `SixFrameConvention.EmbossTranseq`
to reproduce EMBOSS transeq's phase-locked default (or `EmbossTranseqAlternative` for `-alternative Y`).
transeq's `-regions`, `-trim`, `-clean` and `-methionine` (regions only) are not exposed. ORF finding provides all seven getorf `-find` modes (legacy overload:
`-find 1` only), on linear or circular (`-circular`) sequences, A/C/G/T or IUPAC strings. ORF coordinates for
reverse-strand hits are expressed in the reverse complement's coordinate frame
(`OrfFeature.ToGetOrfRange` gives getorf's forward-strand header).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var dna = new DnaSequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
var frames = Translator.TranslateSixFrames(dna);
// frames[1]  == "MAIVMGR*KGAR*"
// frames[-1] == "LSGTLSAAHYNGH"

var orfs = Translator.FindOrfs(new DnaSequence("GGGATGAAACCCTAAGGG"),
                               minLength: 1, searchBothStrands: false);
// one ORF: Start=3, End=14, Frame=1, Protein="MKP"

// getorf -find 0 -minsize 9 (EMBOSS 6.6.0): six STOP→STOP regions, first one
var regions = Translator.FindOrfs(new DnaSequence("GGGATGAAACCCTAAGGG"),
                                  OrfFindMode.TranslationStopToStop, minSize: 9).ToList();
// regions[0]: Start=0, End=14, Frame=1, IsTerminated=true, Sequence="GMKP"; ToGetOrfRange(18) == (1, 12)
// FlankingInitialStop, flanking: 2 → [3 - 9] "GATGAAA" and [16 - 10] (REVERSE SENSE) "CTTAGGG"

// IUPAC genome (string overload): N codon read through as X, as getorf prints it
var amb = Translator.FindOrfs("ATGNNNCCCTAA", OrfFindMode.TranslationStartToStop, minSize: 3,
                              searchBothStrands: false).Single();   // [1 - 9] "MXP"

// getorf -circular: the ORF crosses the origin (EndPosition 20 ≥ L = 15)
var wrapped = Translator.FindOrfs("CCCTAAGGGATGAAA", OrfFindMode.TranslationStartToStop, minSize: 3,
                                  searchBothStrands: false, circular: true).Single();
// Start=9, End=20, Sequence="MKP"; ToGetOrfRange(15) == (10, 18) "(ORF crosses the breakpoint)"
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [Translator_SixFrames_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/Translator_SixFrames_Tests.cs) — covers `INV-01`..`INV-06`
- Tests: [Translator_GetOrfModes_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/Translator_GetOrfModes_Tests.cs) — getorf `-find` 0–6 against the EMBOSS 6.6.0 binary oracle `TestData/Getorf/getorf_find_modes.tsv`
- Tests: [Translator_IupacCircular_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/Translator_IupacCircular_Tests.cs) — IUPAC string overloads (Biopython six frames; getorf ambiguity rules) and `-circular`, oracle `TestData/Getorf/getorf_iupac_circular.tsv`
- Tests: [Translator_SixFrameConvention_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/Translator_SixFrameConvention_Tests.cs) — EMBOSS transeq conventions against the transeq 6.6.0 binary oracle `TestData/Transeq/transeq_six_frames.tsv`
- Evidence: [TRANS-SIXFRAME-001-Evidence.md](../../../docs/Evidence/TRANS-SIXFRAME-001-Evidence.md)

## 8. References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: transeq application documentation. EMBOSS. https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
2. Cock PJA et al. 2009. Biopython. Bioinformatics 25(11):1422–1423. Source Bio/SeqUtils/__init__.py (six_frame_translations). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
3. NCBI. The Genetic Codes — Standard Code (transl_table=1). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
4. Rice P, Longden I, Bleasby A. 2000. EMBOSS: getorf application documentation. EMBOSS. https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html
5. Wikipedia contributors. Reading frame (cites Lodish 2007; Pierce 2012). https://en.wikipedia.org/wiki/Reading_frame
6. EMBOSS getorf source: `emboss/getorf.c` (`getorf_FindORFs`, `getorf_WriteORF`, `getorf_AppORF`) and `emboss/acd/getorf.acd` (`-find` 0–6, `-minsize` 30, `-maxsize` 1000000, `-flanking` 100, `-methionine` Y, `-reverse` Y). https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/getorf.c
7. INSDC. The DDBJ/ENA/GenBank Feature Table Definition — CDS: "sequence of nucleotides that corresponds with the sequence of amino acids in a protein (location includes stop codon)". https://www.insdc.org/submitting-standards/feature-table/
8. EMBOSS translation library: `ajax/core/ajtranslate.c` (`ajTrnCodonstrTypeC`, `ajTrnReadFile`, `trnComplete`, `ajTrnSeqOrig`, `ajTrnSeqDangleC`). https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajtranslate.c; transeq source `emboss/transeq.c`, `emboss/acd/transeq.acd` (`-alternative`, `-trim`, `-clean`)
