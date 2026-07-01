using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Annotation area — open reading frame (ORF) detection
/// (ANNOT-ORF-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, malformed and boundary inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// out-of-bounds indexing (IndexOutOfRangeException / ArgumentOutOfRangeException
/// leaking from internal Substring), no DivideByZero, no NullReferenceException,
/// and no nonsense output (an ORF whose span is not codon-aligned, or that runs off
/// the end of the sequence, or that is reported when the biological contract says it
/// must not be). Every input must resolve to EITHER a well-defined, theory-correct
/// result, OR a documented, intentional validation outcome. A raw runtime exception,
/// a hang, or a blow-up of bogus ORFs on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-ORF-001 — ORF finding
/// Checklist: docs/checklists/03_FUZZING.md, row 28.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: a sequence with NO start codon, an all-start-no-stop
///          sequence, and minLength = 0.
///   • MC = Malformed Content — non-DNA characters (N, digits, IUPAC ambiguity
///          codes, garbage) embedded in the scanned sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The ORF-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// An ORF is a contiguous stretch of codons beginning with a start codon and ending
/// with an in-frame stop codon, with no intervening in-frame stop. The canonical
/// Annotation-area detector is
///   GenomeAnnotator.FindOrfs(string dnaSequence,
///                            int minLength = 100,
///                            bool searchBothStrands = true,
///                            bool requireStartCodon = true)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs lines 91–120).
/// This is the surface the checklist row probes (it is the canonical annotation ORF
/// finder named in ORF_Detection.md §1, §5.1). NOTE there is a distinct Analysis-area
/// surface GenomicAnalyzer.FindOpenReadingFrames (GENOMIC-ORF-001, row 177) with
/// different semantics (ATG-only start, minLength compared to nucleotide span); it is
/// NOT the unit under fuzz here.
///
/// Fixed codon sets (ORF_Detection.md §4.2; GenomeAnnotator.cs StartCodons/StopCodons):
///   • Start codons: ATG, GTG, TTG (common prokaryotic starts).
///   • Stop codons : TAA, TAG, TGA (standard).
/// The scan walks each of the three forward frames (plus three reverse-complement
/// frames when searchBothStrands = true) in codon-sized steps; it accumulates pending
/// start positions and, on reaching an in-frame stop, emits every pending start whose
/// stop-delimited amino-acid length L_aa = (t − s)/3 is ≥ minLength
/// (GenomeAnnotator.cs lines 139–213; ORF_Detection.md §2.2, §4.1).
///
/// Documented parameter / validation contract (ORF_Detection.md §3.3, §6.1):
///   • dnaSequence null or empty → NO ORFs (the method `yield break`s on
///     string.IsNullOrEmpty; GenomeAnnotator.cs lines 97–98). NOT an exception — this
///     is the documented contract for the canonical surface, so the fuzz tests pin the
///     empty/null result, NOT a throw.
///   • The canonical finder performs NO numeric validation of minLength and NO
///     alphabet pre-validation. minLength is simply compared with `>=` against the
///     amino-acid length, and each codon is upper-cased then tested for exact
///     membership in the fixed A/C/G/T codon sets. The fuzz tests below therefore
///     assert the THEORY-CORRECT consequences of that contract, not an exception:
///       – No start codon with requireStartCodon = true → no pending start is ever
///         recorded → empty result (ORF_Detection.md §6.1 row 3).
///       – ALL-START-NO-STOP (e.g. "ATGATGATG…" with no in-frame stop) with the
///         DEFAULT requireStartCodon = true → NOT reported: a terminating stop is
///         REQUIRED before emission, so a start that never reaches a stop yields
///         nothing (ORF_Detection.md §6.1 row 4, INV-02). This is THE key biological
///         boundary for this unit — a partial/unterminated ORF is suppressed by
///         default. Only when requireStartCodon = false does the finder additionally
///         emit trailing frame segments that run to the (codon-aligned) sequence end
///         (GenomeAnnotator.cs lines 183–212). The fuzz tests pin BOTH halves of that
///         contract so the default-suppression cannot silently drift.
///       – minLength = 0 is degenerate: `orfAaLength >= 0` is ALWAYS true, so a
///         zero-amino-acid ORF (a start codon immediately followed by an in-frame
///         stop, e.g. "ATGTAA", with L_aa = 0) IS reported. The canonical contract
///         does not reject minLength = 0 (no validation gate, ORF_Detection.md §3.3);
///         the disciplined outcome here is that the degenerate threshold produces a
///         WELL-FORMED trivial ORF (codon-aligned, in-bounds, correct start/stop),
///         not a crash and not a blow-up of garbage. The fuzz test pins that
///         well-formedness rather than asserting an exception the contract never
///         promises.
///       – Non-DNA characters (N, digits, IUPAC ambiguity codes, punctuation) are
///         NOT pre-validated and NOT rejected: a codon containing such a character
///         simply fails exact start/stop membership and is treated as an ordinary
///         non-boundary codon (ORF_Detection.md §3.3, §6.1 rows 6–7). The fuzz tests
///         pin that malformed content never crashes (no out-of-range Substring on the
///         3-char codon window), never invents a spurious start/stop from a non-ACGT
///         triplet, and that an ORF surrounding internal garbage can still close at a
///         genuine downstream stop.
///
/// Documented invariants pinned on every positive result (ORF_Detection.md §2.4):
///   INV-01 every reported ORF begins with ATG/GTG/TTG (requireStartCodon = true);
///   INV-02 every reported ORF ends with TAA/TAG/TGA (requireStartCodon = true);
///   INV-03 every reported ORF span (End − Start) is divisible by 3;
///   INV-04 coordinates satisfy 0 ≤ Start < End ≤ dnaSequence.Length.
///
/// FindOrfs is a LAZY iterator (`yield`); even the null/empty short-circuit sits
/// inside the iterator body, so it only fires on enumeration. Every test therefore
/// forces enumeration (`.ToList()`) so the documented behavior actually surfaces and
/// any hang would manifest as a non-terminating materialization. Deterministic only:
/// any random fuzz input is generated from a locally fixed `new Random(seed)`.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-GENE-001 — gene prediction (prokaryotic ORF-first heuristic)
/// Checklist: docs/checklists/03_FUZZING.md, row 29.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: an empty sequence, a non-coding-only sequence (no ORF,
///          no RBS), and overlapping genes (overlap bookkeeping must stay valid).
///   • MC = Malformed Content — a sequence that contains NO Shine-Dalgarno (RBS)
///          motif anywhere, so the RBS helper must invent nothing.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The gene-prediction contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Gene prediction here is a prokaryote-oriented, ORF-first heuristic split across two
/// helpers (Gene_Prediction.md §1, §5.1):
///   GenomeAnnotator.PredictGenes(string dnaSequence,
///                                int minOrfLength = 100,
///                                string prefix = "gene")
///     — converts every QUALIFYING ORF (translated length ≥ minOrfLength, found with
///       searchBothStrands: true, requireStartCodon: true) into a CDS GeneAnnotation,
///       ordered by genomic Start, with sequential IDs (Gene_Prediction.md §4.1, §5.2;
///       GenomeAnnotator.cs lines 389–421).
///   GenomeAnnotator.FindRibosomeBindingSites(string dnaSequence,
///                                            int upstreamWindow = 20,
///                                            int minDistance = 4,
///                                            int maxDistance = 15)
///     — independently scans the FORWARD-strand upstream window of every ORF (found with
///       an internal minLength = 30 aa) for Shine-Dalgarno-like motifs at an aligned
///       spacing within [minDistance, maxDistance] (GenomeAnnotator.cs lines 267–384).
///
/// CRITICAL: the two helpers are INDEPENDENT. PredictGenes does NOT require, call, or
/// consider an RBS — it emits genes purely from ORF structure (Gene_Prediction.md §5.2,
/// §5.4 deviation 1). So "no RBS in seq" suppresses RBS HITS, not gene predictions. The
/// fuzz tests below pin that exact split rather than a coupled model the code never
/// implements. The SD-like motif library is AGGAGG, GGAGG, AGGAG, GAGG, AGGA; the RBS
/// score is motif.Length / 6.0 (Gene_Prediction.md §4.2).
///
/// Documented parameter / validation contract (Gene_Prediction.md §3.3, §6.1):
///   • dnaSequence null or empty → BOTH helpers return an empty sequence (PredictGenes
///     via the underlying FindOrfs short-circuit; FindRibosomeBindingSites via its own
///     IsNullOrEmpty `yield break`). NOT an exception — pinned as the empty result.
///   • No valid ORF (non-coding-only input) → no genes AND no RBS hits, because RBS
///     scanning is ORF-driven (Gene_Prediction.md §6.1 row 2).
///   • Overlapping / nested ORFs → every qualifying ORF is emitted as its OWN CDS;
///     there is NO best-model selection or overlap suppression (Gene_Prediction.md §5.2
///     bullet 2, §5.3 "Intentionally simplified", §6.2). The disciplined outcome is that
///     overlap bookkeeping never crashes and every emitted gene keeps VALID coordinates.
///   • Non-DNA / no-motif content is not pre-validated; a region lacking any SD-like
///     motif simply produces no RBS hit (Gene_Prediction.md §3.3).
///
/// Documented invariants pinned on every predicted gene (Gene_Prediction.md §2.4):
///   GENE-INV-01 every gene derives from an ORF starting ATG/GTG/TTG and ending
///               TAA/TAG/TGA (PredictGenes delegates with requireStartCodon: true);
///   GENE-INV-02 every gene has valid bounds 0 ≤ Start < End ≤ length and Frame ∈ {1,2,3},
///               Strand ∈ {'+','-'}, Type == "CDS";
///   GENE-INV-03 every reported RBS hit lies at aligned distance ∈ [minDistance,
///               maxDistance] from its ORF start;
///   GENE-INV-04 every RBS score ∈ [4/6, 1.0] for the supported motif library.
/// PredictGenes / FindRibosomeBindingSites are LAZY iterators (`yield`), so every test
/// forces enumeration (`.ToList()`). Deterministic only: random fuzz input comes from a
/// locally fixed `new Random(seed)`.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-PROM-001 — promoter (−35 / −10 box) motif detection
/// Checklist: docs/checklists/03_FUZZING.md, row 30.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: NO −10 box present, the (inapplicable) threshold = 0,
///          an empty sequence, and a sequence SHORTER than the motif window.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The promoter-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Promoter detection here is a bacterial MOTIF-SEARCH helper, NOT a full TSS predictor
/// and NOT a −35/−10 pairing model (Promoter_Detection.md §1, §5.2):
///   GenomeAnnotator.FindPromoterMotifs(string dnaSequence)
///     → IEnumerable&lt;(int position, string type, string sequence, double score)&gt;
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs lines 561–615).
/// It uppercases the input, then scans it INDEPENDENTLY for each member of two small fixed
/// motif libraries, emitting EVERY exact substring hit with a hard-coded, literature-derived
/// score (Promoter_Detection.md §4.1, §4.2):
///   • −35 box variants: TTGACA(1.000), TTGAC(0.855), TGACA(0.815), TTGA(0.710);
///   • −10 box (Pribnow) variants: TATAAT(1.000), TATAA(0.801), ATAAT(0.813), TATA(0.665).
/// The scan loop is `for (int i = 0; i &lt;= seq.Length - motif.Length; i++)`, so when the
/// sequence is shorter than a motif the upper bound is NEGATIVE and the loop body never
/// runs — there is NO out-of-range Substring (Promoter_Detection.md §6.1). This is THE key
/// boundary for this unit (checklist target "seq shorter than motif").
///
/// CRITICAL — there is NO threshold parameter. FindPromoterMotifs takes ONLY a sequence;
/// it has no score gate, no minimum-score cutoff, no −35/−10 spacing check, and no pairing
/// (Promoter_Detection.md §3.1, §5.2, §6.2). The checklist's "threshold = 0" target is
/// therefore DEGENERATE/INAPPLICABLE for this surface: there is no threshold to set to 0,
/// so the theory-correct contract is that the method ALWAYS emits every exact motif hit at
/// its full hard-coded score, and every emitted score is strictly POSITIVE (∈ [0.665, 1.0],
/// the variant library's range) — i.e. a "threshold of 0" admits exactly the same hits the
/// method already emits unconditionally. The fuzz test pins that no score-gating exists and
/// no hit is suppressed or zero-scored, rather than asserting a parameter the contract never
/// exposes (never weaken a test to a no-op: it asserts a real, falsifiable consequence).
///
/// Documented parameter / validation contract (Promoter_Detection.md §3.3, §6.1):
///   • Empty sequence → NO motifs: the scan loops have no valid start position
///     (`0 <= 0 - motifLen` is false). NOT an exception — pinned as the empty result.
///   • Sequence shorter than a motif → NO motifs and NO IndexOutOfRange: the negative
///     upper bound suppresses the loop body for that motif (KEY boundary).
///   • No −10 box present → NO −10 hits: matching is EXACT against the fixed library, so a
///     sequence lacking TATAAT/TATAA/ATAAT/TATA yields no −10 motif (it may still legitimately
///     yield −35 hits; the two libraries are scanned independently).
///   • Null input → the method dereferences `dnaSequence` immediately via ToUpperInvariant()
///     with NO null guard, so null throws (Promoter_Detection.md §3.3, §6.1 "Null input — not
///     explicitly handled"). The fuzz test pins the DOCUMENTED throw rather than weakening it.
///   • Non-ACGT / lowercase content is not rejected: input is uppercased, then only exact
///     library substrings match, so non-library characters simply never match.
///
/// Documented invariants pinned on every positive result (Promoter_Detection.md §2.4):
///   PROM-INV-01 every reported score ∈ [0, 1] (a normalized fraction of full consensus);
///   PROM-INV-02 a full TTGACA / TATAAT match scores exactly 1.0;
///   PROM-INV-03 every reported motif belongs to the fixed −35 / −10 variant library;
///   PROM-INV-04 every reported position is a valid 0-based index whose motif window stays
///               in bounds: 0 ≤ position ≤ length − motif.Length.
/// FindPromoterMotifs is a LAZY iterator (`yield`), so every test forces enumeration
/// (`.ToList()`); any hang would manifest as a non-terminating materialization. Deterministic
/// only: random fuzz input comes from a locally fixed `new Random(seed)`.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-GFF-001 — GFF3 I/O (lightweight Annotation-area parser + serializer)
/// Checklist: docs/checklists/03_FUZZING.md, row 31.
/// Fuzz strategies exercised for THIS unit:
///   • TF  = Truncated Fields — a data line with FEWER than 9 tab-separated columns
///           (missing fields), plus a present-but-EMPTY strand column.
///   • MC  = Malformed Content — non-numeric start/end/score/phase, an out-of-domain
///           strand symbol, and garbage data lines mixed with valid ones.
///   • INJ = Injection — TAB injection: an extra literal tab inside a field (read side),
///           and attribute payloads carrying tab / newline / ';' / '=' / '&' / ',' that
///           would break GFF3 column structure if not escaped (write side).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Which surface is fuzzed, and why ANNOT (not FileIO PARSE-GFF-001 / row 68)
/// ───────────────────────────────────────────────────────────────────────────
/// ANNOT-GFF-001 is the Annotation-area GFF3 surface documented in GFF3_IO.md §5.1: the
/// LIGHTWEIGHT helpers on GenomeAnnotator —
///   • GenomeAnnotator.ParseGff3(IEnumerable&lt;string&gt;) → IEnumerable&lt;GenomicFeature&gt;
///     (reader; GenomeAnnotator.cs lines 426–481), and
///   • GenomeAnnotator.ToGff3(IEnumerable&lt;GeneAnnotation&gt;, string seqId) → IEnumerable&lt;string&gt;
///     plus its column-9 EncodeGff3Value helper (serializer; GenomeAnnotator.cs lines 488–554).
/// Both reading AND writing live in the Annotation area, so per the checklist this unit fuzzes
/// the Annotation-area helper surface. The SEPARATE, fuller FileIO parser
/// Seqeron.Genomics.IO.GffParser (PARSE-GFF-001, row 68) is deliberately NOT touched here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GFF3 I/O contract under test (GFF3_IO.md §2.2, §3.3, §4.2, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// A GFF3 data line is NINE TAB-separated columns: seqid, source, type, start, end, score,
/// strand (+/-/./?), phase, attributes. The lightweight reader:
///   • skips blank and '#'-prefixed lines, and skips any line with FEWER than 9 tab fields
///     (TF/missing-fields → SKIP, never an IndexOutOfRange; §4.2, §6.1);
///   • parses '.' as the null sentinel for score (col 6) and phase (col 8);
///   • fills a missing ID attribute with a running `feature_n` id;
///   • percent-DECODES attribute values with Uri.UnescapeDataString.
///
/// CALLER-SAFETY HARDENING applied for this unit (source fix — see report): the lightweight
/// reader previously used int.Parse / double.Parse on the numeric columns and `parts[6][0]`
/// on the strand column, so a MALFORMED-BUT-9-COLUMN data line crashed with an uncaught
/// FormatException (non-numeric start/end/score/phase) or IndexOutOfRangeException (a
/// present-but-EMPTY strand field) — exactly the undisciplined fuzz failure modes
/// (03_FUZZING.md §Description; ADVANCED_TESTING_CHECKLIST.md §8). These are now TOLERANT:
/// a malformed numeric column makes the reader SKIP that data line (mirroring the existing
/// &lt;9-field skip and the sibling GffParser.ParseLine discipline), and an empty strand column
/// defaults to the valid '.' sentinel — no throw, no garbage record. The fuzz tests below pin
/// that disciplined skip/normalize behavior, NOT a crash.
///
/// Strand DOMAIN: the spec strand alphabet is {+, -, ., ?}. The lightweight reader does NOT
/// reject an out-of-domain strand symbol; it stores `parts[6][0]` verbatim (GFF3_IO.md §3.2
/// "Parsed column 7 strand symbol"). So an invalid strand 'Z' is parsed as the char 'Z', not
/// rejected — the disciplined outcome is a deterministic verbatim parse, never a crash. The
/// fuzz test pins that documented verbatim contract rather than asserting a rejection the
/// helper never performs.
///
/// TAB INJECTION (read side): an EXTRA tab inside a field yields MORE than 9 columns; because
/// the eligibility gate is `parts.Length < 9` (not `!= 9`), the line is NOT skipped — the
/// reader consumes columns 0..8 and the injected-tab tail beyond column 9 is dropped. This is
/// deterministic and crash-free; the fuzz test pins that an injected tab cannot smuggle a tail
/// fragment into the parsed attributes nor split the record into wrong columns.
///
/// SERIALIZER (write side, the INJ high-value target): ToGff3 emits `##gff-version 3` then one
/// 9-column row per GeneAnnotation, writing col 4 as Start + 1 (0-based internal → 1-based GFF;
/// INV-01), source '.', phase '0' for CDS else '.', and percent-ENCODING every column-9 reserved
/// character via EncodeGff3Value: tab→%09, newline→%0A, CR→%0D, '%'→%25, ';'→%3B, '='→%3D,
/// '&'→%26, ','→%2C, and control chars → %XX (GFF3_IO.md §4.2; GenomeAnnotator.cs lines 530–553).
/// A well-behaved serializer must NEVER emit a structurally broken line when an attribute value
/// contains a tab/newline/';'/'='/'&'/',' — it must ESCAPE them. The fuzz tests pin that even an
/// attribute value stuffed with those exact GFF-breaking characters serializes to a SINGLE
/// 9-column row, and that the value round-trips byte-exact back through ParseGff3.
///
/// GFF3 invariants pinned on positive results (GFF3_IO.md §2.4):
///   GFF-INV-01 external GFF3 coordinates are 1-based inclusive (ToGff3 writes Start + 1);
///   GFF-INV-02 attribute keys are case-sensitive;
///   GFF-INV-03 phase is '0' for CDS, '.' otherwise.
/// ParseGff3 / ToGff3 are LAZY iterators (`yield`), so every test forces enumeration (`.ToList()`)
/// to surface the documented behavior; any hang would manifest as a non-terminating
/// materialization. Deterministic only: no randomness is required for this text-shaped unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-CODING-001 — coding-potential scoring (CPAT hexamer usage-bias)
/// Checklist: docs/checklists/03_FUZZING.md, row 216.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: a RANDOM sequence (no coding bias → classified non-coding
///          under non-coding-biased tables), a PERFECT ORF (strong coding bias →
///          classified coding under coding-biased tables), and the EMPTY/too-short
///          sequence (defined degenerate score 0, NO divide-by-zero).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The coding-potential contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// CPAT's hexamer usage-bias score is the mean per-hexamer log-likelihood ratio of
/// in-frame hexamer frequencies between a coding and a non-coding training set
/// (Coding_Potential_Calculation.md §1, §2.2). The surface is
///   GenomeAnnotator.CalculateCodingPotential(
///       string sequence,
///       IReadOnlyDictionary&lt;string,double&gt; codingHexamerFrequencies,
///       IReadOnlyDictionary&lt;string,double&gt; noncodingHexamerFrequencies,
///       int wordSize = 6, int stepSize = 3) → double
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs lines 693–743).
/// It up-cases the sequence and walks frame-0 in-frame hexamers (offset 0, step
/// stepSize, full-length words only). For each hexamer present in BOTH tables it adds
/// (Coding_Potential_Calculation.md §2.2, verified element-by-element against the source
/// and the canonical CPAT/lncScore FrameKmer.kmer_ratio branches, NOT echoed off the code):
///   • F&gt;0 and F'&gt;0  → ln(F/F')   (NATURAL log);
///   • F&gt;0 and F'==0 → +1;
///   • F==0 and F'==0 → skipped (`continue`), NOT counted;
///   • F==0 and F'&gt;0  → −1;
///   • a hexamer missing from EITHER table → skipped, NOT counted.
/// Score = (Σ contributions) / (count of scored hexamers); 0 when nothing is scorable
/// (Coding_Potential_Calculation.md §2.2, §5.4 deviation 1 — the reference's −1 sentinel
/// is intentionally replaced by 0 here; both are non-scores).
///
/// Documented parameter / validation contract (Coding_Potential_Calculation.md §3.3, §6.1):
///   • null sequence or either table → ArgumentNullException; wordSize ≤ 0 / stepSize ≤ 0
///     → ArgumentOutOfRangeException. These are the DOCUMENTED throws and are pinned as such.
///   • sequence.Length &lt; wordSize (INCLUDING empty) → returns 0, NOT a throw and NOT a
///     DivideByZero: the short-circuit fires before any hexamer is scored (INV-06). This is
///     THE empty/degenerate boundary the checklist row targets.
///   • Non-ACGT content is NOT rejected: a hexamer containing N/garbage simply is not a table
///     key, so it is skipped — no exception (§3.3, §6.1).
///
/// Documented invariants pinned on positive results (Coding_Potential_Calculation.md §2.4):
///   INV-01 score = (Σ per-hexamer contributions) / (count of scored hexamers);
///   INV-02 per-hexamer contribution = ln(F/F') (natural log) when both &gt; 0;
///   INV-03 coding-only hexamer ⇒ +1, non-coding-only ⇒ −1;
///   INV-04 coding-biased tables ⇒ score &gt; 0 (classified coding); non-coding-biased ⇒
///          score &lt; 0 (classified non-coding) — the coding/non-coding decision threshold is 0;
///   INV-05 only in-frame full-length hexamers (offset 0, step stepSize, length = wordSize)
///          are scored;
///   INV-06 sequence.Length &lt; wordSize ⇒ score = 0.
/// The fuzz bar for this unit: no crash / hang / NaN / Infinity on degenerate or random
/// input, the empty/too-short boundary returns a defined 0 with no divide-by-zero, the random
/// sequence is classified NON-coding (score &lt; 0 against non-coding-biased tables) while a
/// perfect coding ORF is classified CODING (score &gt; 0 against coding-biased tables), and the
/// exact worked-example score 0.34657359027997264 (Coding_Potential_Calculation.md §7.1) is
/// reproduced from first principles. Deterministic only: random fuzz input and random
/// frequency tables come from a locally fixed `new Random(seed)` — NEVER a shared static Rng.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-CODONUSAGE-001 — Relative Synonymous Codon Usage (RSCU)
/// Checklist: docs/checklists/03_FUZZING.md, row 217.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: the EMPTY input (empty collection / all-empty strings),
///          a SINGLE codon, and a sequence whose LENGTH is NOT a multiple of 3
///          (a partial trailing codon). — docs/checklists/03_FUZZING.md §Description
///          ("BE = Boundary Exploitation: 0, -1, MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The RSCU contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// RSCU measures within-synonymous-family codon bias, normalised away from amino-acid
/// composition (Relative_Synonymous_Codon_Usage.md §1, §2.2). The surface is
///   GenomeAnnotator.GetCodonUsage(IEnumerable&lt;string&gt; codingSequences)              [Standard code]
///   GenomeAnnotator.GetCodonUsage(IEnumerable&lt;string&gt; codingSequences, GeneticCode)  [supplied code]
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs lines 1044–1114).
/// NOTE this is the ANNOTATION-module RSCU (ANNOT-CODONUSAGE-001), DISTINCT from the
/// Codon-area CODON-USAGE-001 (row 61). It is ALSO distinct from the pre-existing raw
/// integer-count GetCodonUsage(string) overload (line 1000), which is NOT RSCU and is not
/// the unit under fuzz. The string-collection overloads returning IReadOnlyDictionary&lt;string,double&gt;
/// are the RSCU surface named in the doc §5.1.
///
/// For amino acid i with n_i synonymous codons and observed counts x_(i,j), the doc §2.2
/// formula (Sharp &amp; Li 1986, verbatim from LIRMM/CodonU, NOT echoed off the code) is
///   RSCU_(i,j) = n_i · x_(i,j) / Σ_j x_(i,j)
/// computed over counts pooled across all input sequences. The result maps EVERY sense
/// codon (uppercase DNA) to its RSCU value; an entirely unobserved family yields 0.0 for
/// each member (no division by zero; the CAI 0.5 pseudocount is deliberately NOT applied —
/// doc §3.3, §5.3, §5.4 deviation 1).
///
/// Documented framing / alphabet contract (Relative_Synonymous_Codon_Usage.md §3.1, §3.3, §4.1):
///   • Sequences are read IN FRAME from index 0 in non-overlapping triplets; a PARTIAL
///     trailing codon (length not a multiple of 3) is IGNORED — `i <= len - 3` step 3.
///   • Only A/C/G/T codons are counted (case-insensitive; input uppercased first); a codon
///     containing N/IUPAC/garbage is simply not counted.
///   • Null/empty individual sequences are SKIPPED; the genetic-code table is RNA-keyed (U)
///     and reconciled to DNA (T) internally; STOP codons are excluded (sense codons only).
///
/// Documented parameter / validation contract (doc §3.3, §6.1):
///   • Null `codingSequences` or null `code` → ArgumentNullException (DOCUMENTED throws,
///     pinned as such — never weakened to a tolerant no-op).
///   • Empty collection / all-empty strings → NOTHING observed: every family totals 0, so
///     EVERY sense codon maps to RSCU 0.0 (NO division by zero). This is THE empty boundary
///     the checklist row targets. The output still enumerates all 61 sense codons of the
///     Standard code (64 − 3 stops), so it is a defined, degenerate result, not an exception.
///
/// Documented invariants pinned on positive results (Relative_Synonymous_Codon_Usage.md §2.4):
///   INV-01 for every OBSERVED synonymous family, Σ RSCU over its codons = n_i;
///   INV-02 a single-codon amino acid (Met=ATG, Trp=TGG) always has RSCU = 1.0 (n_i = 1);
///   INV-03 every RSCU value lies in [0, n_i] (so within [0, 6] for the Standard code,
///          whose largest family — Leu/Ser/Arg — has n_i = 6);
///   INV-04 stop codons never appear in the output (sense codons only);
///   INV-05 uniform usage within a family ⇒ every member RSCU = 1.0.
/// The fuzz bar for this unit: no crash / hang / NaN / Infinity / corruption on degenerate
/// or random input; the EMPTY boundary returns a defined all-zero 61-codon map with NO
/// divide-by-zero; a SINGLE codon yields the exact n_i·x/Σ value for its family member; a
/// length NOT a multiple of 3 drops only the partial trailing codon and never crashes; and
/// the doc §7.1 Leucine worked example (CTT=3.0, CTG=1.5, TTA=1.5, others 0; Σ family = 6)
/// is reproduced from first principles. The dictionary values are pinned to be finite and
/// non-negative on a fixed-seed random sweep. Deterministic only: random fuzz input comes
/// from a locally fixed `new Random(seed)` — NEVER a shared static Rng.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-REPEAT-001 — repetitive-element detection (tandem + inverted repeats)
/// Checklist: docs/checklists/03_FUZZING.md, row 218.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row "no repeat, full repeat, minLen edge": a sequence with NO
///          repeat (empty result), a FULL repeat (the WHOLE sequence is one tandem
///          array), and the minLength EDGE (an array whose span is exactly the
///          threshold is reported; one base below the threshold is dropped). Plus the
///          empty/null sequence and the out-of-range minCopies/minRepeatLength guards.
/// — docs/checklists/03_FUZZING.md §Description ("BE = Boundary Exploitation: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The repetitive-element contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// This is the ANNOTATION-module repeat finder (Repetitive_Element_Detection.md,
/// Test Unit ID ANNOT-REPEAT-001), DISTINCT from the Repeats-area REP-* rows 13–17.
/// The surface is
///   GenomeAnnotator.FindRepetitiveElements(string dnaSequence,
///                                          int minRepeatLength = 10,
///                                          int minCopies = 2)
///     → IEnumerable&lt;(int start, int end, string type, string sequence)&gt;
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs lines 774–914).
/// It reports two structural classes (Repetitive_Element_Detection.md §1, §2.2):
///   • "tandem_repeat"  — a PRIMITIVE motif u repeated head-to-tail in ≥ minCopies
///                        directly adjacent copies; the array spans c·|u| bases and the
///                        reported `sequence` is the array slice dnaSequence[start..end]
///                        (0-based, end-EXCLUSIVE); non-primitive units (AA = A×2) are
///                        collapsed to the primitive period (INV-01, INV-02; §5.2).
///   • "inverted_repeat" — a left arm W followed by a gap G (|G| ≥ 0) then the right arm
///                        W̄ᴿ = revcomp(W); gap 0 ⇒ reverse-complement palindrome. The
///                        right arm is accepted only when it ordinal-equals revcomp(left)
///                        (INV-03; §2.2, §5.2).
///
/// Documented parameter / validation contract (Repetitive_Element_Detection.md §3.1, §3.3, §6.1):
///   • dnaSequence null → ArgumentNullException (DOCUMENTED throw, pinned as such).
///   • dnaSequence empty → EMPTY result, NOT a throw (the core `yield break`s on length 0).
///   • minCopies &lt; 2 → ArgumentOutOfRangeException (a tandem repeat needs "two or more"
///     copies; Wikipedia: Tandem repeat). minRepeatLength &lt; 1 → ArgumentOutOfRangeException.
///   • minRepeatLength is the MINIMUM TOTAL SPAN (bp) of a reported element: a tandem array
///     is dropped when end − start &lt; minRepeatLength (GenomeAnnotator.cs line 840). It is
///     ALSO reused as the inverted-repeat minimum ARM length (`minArmLength`), so an inverted
///     pair needs i + 2·minRepeatLength ≤ length to exist (line 887). The minLength-edge fuzz
///     tests therefore isolate the TANDEM contribution (filter type == "tandem_repeat") where
///     the span filter is the quantity under test.
///   • Non-ACGT / lowercase content is NOT pre-validated: the input is upper-cased, then exact
///     matching simply never fires for non-library/non-matching triplets (§3.3, §6.1).
///
/// Documented invariants pinned on positive results (Repetitive_Element_Detection.md §2.4):
///   INV-01 a tandem repeat's `sequence` equals dnaSequence[start..end] (0-based, exclusive
///          end) and its length is an integer multiple of the primitive unit length;
///   INV-02 a tandem repeat has ≥ minCopies (≥ 2) directly adjacent identical primitive copies;
///   INV-03 for every inverted repeat the right arm equals revcomp(left arm).
/// The fuzz bar for this unit: no crash / hang / corruption on degenerate or random input;
/// the NO-repeat boundary returns an empty result; the FULL-repeat boundary reports the WHOLE
/// sequence as one tandem array with the correct primitive unit; the minLength EDGE reports an
/// array of span exactly = minRepeatLength and DROPS the same array one base above the threshold;
/// the doc §7.1 worked example FindRepetitiveElements("ATTCGATTCGATTCG", 5, 2) →
/// (0, 15, "tandem_repeat", "ATTCGATTCGATTCG") (unit "ATTCG" ×3) is reproduced from first
/// principles; and a fixed-seed random sweep keeps every reported element structurally valid
/// (INV-01..INV-03, in-bounds, codon of the right type). Deterministic only: random fuzz input
/// comes from a locally fixed `new Random(seed)` — NEVER a shared static Rng.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AnnotationFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-ORF-001 — ORF finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-ORF-001 — ORF finding

    #region BE — Boundary: no start codon (no ATG/GTG/TTG) in the sequence

    /// <summary>
    /// BE: a sequence containing NO start codon at all. With the default
    /// requireStartCodon = true, no pending ORF start is ever recorded, so even a
    /// sequence riddled with stop codons emits NOTHING — no crash, no spurious ORF
    /// (ORF_Detection.md §6.1 row 3, INV-01). We use "CCC...TAA..." built entirely from
    /// codons that are neither starts nor (for the body) stops, plus a stop, and scan
    /// both strands so a start cannot sneak in on the reverse complement either. The
    /// body codon "CCC" reverse-complements to "GGG", still not a start; the stop "TAA"
    /// reverse-complements within frame to "TTA"/"TAA" variants — none are starts. The
    /// disciplined outcome is the empty result.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_NoStartCodon_ReturnsEmptyAndDoesNotThrow()
    {
        // No ATG/GTG/TTG anywhere; only C/A bodies and a TAA stop.
        const string seq = "CCCCCCCCCCCCCCCCCCTAA";

        var act = () => GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true).ToList();
        act.Should().NotThrow("a sequence with no start codon never records a pending ORF; it simply yields nothing");

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true).ToList();
        orfs.Should().BeEmpty(
            "with requireStartCodon = true and no ATG/GTG/TTG present, no ORF can be emitted on either strand");
    }

    #endregion

    #region BE — Boundary: all start, no stop (the key biological boundary)

    /// <summary>
    /// BE: "ATGATGATG…" — a start codon tiled repeatedly with NO in-frame stop. This is
    /// THE key biological boundary for ORF finding: a start that never reaches a stop is
    /// an UNTERMINATED ORF. Under the DEFAULT requireStartCodon = true the canonical
    /// finder REQUIRES a terminating stop before emission, so the unterminated frame is
    /// suppressed — empty result, no crash, no off-the-end Substring
    /// (ORF_Detection.md §6.1 row 4, INV-02; GenomeAnnotator.cs lines 158–187 emit only
    /// on a stop). We restrict to the forward strand so a reverse-complement stop cannot
    /// terminate it: revcomp("ATG") = "CAT", giving "CATCATCAT…", which has no stop
    /// either, so even both-strand search would stay empty — but pinning forward-only
    /// keeps the assertion about the all-start frame unambiguous.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_AllStartNoStop_RequireStartCodon_ReturnsEmpty()
    {
        // 30 nt of pure 'ATG' — in frame +1 it is ten start codons and never a stop.
        string allStart = string.Concat(Enumerable.Repeat("ATG", 10));

        var act = () => GenomeAnnotator.FindOrfs(allStart, minLength: 1, searchBothStrands: false).ToList();
        act.Should().NotThrow("an unterminated start-only frame is scanned to its end without indexing past it");

        var orfs = GenomeAnnotator.FindOrfs(allStart, minLength: 1, searchBothStrands: false).ToList();
        orfs.Should().BeEmpty(
            "INV-02 / §6.1 row 4: with requireStartCodon = true a stop is REQUIRED, so a start that never reaches a stop is not reported");
    }

    /// <summary>
    /// BE (other half of the contract): the SAME all-start-no-stop input WITH
    /// requireStartCodon = false. Here the canonical finder additionally emits trailing
    /// frame segments that run to the codon-aligned sequence end even without a stop
    /// (GenomeAnnotator.cs lines 183–212). Pinning this half proves the default
    /// (requireStartCodon = true) suppression above is a deliberate switch, not an
    /// accident, and that the partial-ORF path is itself disciplined: the reported span
    /// is codon-aligned (INV-03) and never runs past the sequence end (INV-04). For
    /// "ATGATGATG" (9 nt, frame +1) the trailing segment from the first pending start
    /// runs to the aligned end, length divisible by 3.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_AllStartNoStop_NoRequireStartCodon_EmitsCodonAlignedInBoundsPartialOrf()
    {
        string allStart = string.Concat(Enumerable.Repeat("ATG", 3)); // "ATGATGATG", 9 nt

        var orfs = GenomeAnnotator.FindOrfs(allStart, minLength: 0, searchBothStrands: false, requireStartCodon: false).ToList();

        orfs.Should().NotBeEmpty(
            "with requireStartCodon = false the canonical finder emits the trailing unterminated frame segment");
        orfs.Should().OnlyContain(o =>
            (o.End - o.Start) % 3 == 0 &&                       // INV-03: codon-aligned span
            o.Start >= 0 && o.Start < o.End && o.End <= allStart.Length, // INV-04: in bounds
            "even the partial-ORF path stays codon-aligned and within the sequence; it never runs off the end");
    }

    #endregion

    #region BE — Boundary: minLength = 0 (degenerate threshold)

    /// <summary>
    /// BE: minLength = 0 is the degenerate length floor. `orfAaLength >= 0` is ALWAYS
    /// true, so a zero-amino-acid ORF — a start codon immediately followed by an in-frame
    /// stop, e.g. "ATGTAA" with L_aa = (3 − 0)/3 − 1 = 0 — IS reported. The canonical
    /// contract has NO validation gate on minLength (ORF_Detection.md §3.3), so the
    /// disciplined outcome is NOT an exception but a WELL-FORMED trivial ORF: it must
    /// begin with a start, end with a stop, span exactly 6 nt (start + stop), be
    /// codon-aligned, and stay in bounds. We pin that the degenerate threshold yields a
    /// clean trivial result, never a crash and never a blow-up of garbage.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_MinLengthZero_EmitsWellFormedTrivialOrf()
    {
        const string seq = "ATGTAA"; // start immediately followed by stop → 0-aa ORF

        var act = () => GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: false).ToList();
        act.Should().NotThrow("minLength = 0 is not rejected by the canonical finder; it admits the 0-aa ORF without crashing");

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: false).ToList();

        var trivial = orfs.Should().ContainSingle(
            "exactly one frame (+1) contains the immediate start→stop pair").Subject;
        trivial.Sequence.Should().Be("ATGTAA", "the 0-aa ORF spans the start codon plus the terminal stop");
        trivial.Start.Should().Be(0);
        trivial.End.Should().Be(6);
        ((trivial.End - trivial.Start) % 3).Should().Be(0, "INV-03: the span is codon-aligned");
        trivial.End.Should().BeLessThanOrEqualTo(seq.Length, "INV-04: the ORF stays within the sequence bounds");
    }

    /// <summary>
    /// BE: minLength = 0 must NOT turn into a non-terminating or out-of-range scan on a
    /// longer mixed sequence — the threshold relaxation must not derail the codon walk.
    /// On a fixed-seed random 600-nt sequence the finder must complete promptly and every
    /// reported ORF (however short) must still satisfy all four invariants. This is the
    /// "degenerate parameter does not corrupt the scan" guard for minLength = 0.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindOrfs_MinLengthZero_RandomSequence_AllResultsWellFormed()
    {
        string seq = RandomDna(600, seed: 28_001);

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: true).ToList();

        orfs.Should().OnlyContain(o =>
            (o.End - o.Start) % 3 == 0 &&                           // INV-03
            o.Start >= 0 && o.Start < o.End && o.End <= seq.Length, // INV-04
            "even at the degenerate minLength = 0 every ORF on random input is codon-aligned and in bounds");
    }

    #endregion

    #region MC — Malformed content: non-DNA characters

    /// <summary>
    /// MC: non-DNA characters (N, IUPAC ambiguity codes, digits, punctuation) embedded
    /// in the sequence must NOT crash. The canonical finder does not pre-validate the
    /// alphabet (ORF_Detection.md §3.3): each 3-char codon window is upper-cased and
    /// tested for exact membership in the fixed A/C/G/T codon sets, so a triplet
    /// containing a non-ACGT character simply fails start/stop matching and is treated as
    /// an ordinary non-boundary codon — no out-of-range Substring, no spurious start or
    /// stop invented from garbage. Here the whole sequence is non-DNA filler with no
    /// valid start, so the result is empty on both strands.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_NonDnaCharacters_DoesNotCrashAndEmitsNoSpuriousOrf()
    {
        // Length divisible by 3 (12) but full of non-ACGT garbage and ambiguity codes;
        // no triplet here is a start codon, so nothing can be emitted.
        const string garbage = "NNNRYK123!?WSM";

        var act = () => GenomeAnnotator.FindOrfs(garbage, minLength: 0, searchBothStrands: true).ToList();
        act.Should().NotThrow("non-DNA codons fail set membership but never index past the 3-char window or throw");

        var orfs = GenomeAnnotator.FindOrfs(garbage, minLength: 0, searchBothStrands: true).ToList();
        orfs.Should().BeEmpty(
            "no non-ACGT triplet is a start codon, so malformed filler invents no spurious ORF");
    }

    /// <summary>
    /// MC: an ORF whose INTERNAL coding codons contain non-DNA characters must still
    /// close at a genuine downstream stop. Only start/stop boundaries are matched against
    /// the fixed codon sets, so an internal 'NNN' codon is just a non-boundary codon and
    /// the ORF continues until a recognized stop (ORF_Detection.md §6.1 row 7). We embed
    /// "ATG" + "NNN" + "TAA" in frame +1: the ORF must open at the ATG, ride through the
    /// ambiguous internal codon, and close at the TAA — codon-aligned and in bounds.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_NonDnaInternalCodon_OrfStillClosesAtGenuineStop()
    {
        const string seq = "ATGNNNTAA"; // start, ambiguous internal codon, stop — frame +1

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 0, searchBothStrands: false, requireStartCodon: true).ToList();

        var orf = orfs.Should().ContainSingle(o => o.Start == 0).Subject;
        orf.Sequence.Should().Be("ATGNNNTAA",
            "the ORF opens at ATG, rides through the ambiguous internal NNN codon, and closes at the genuine TAA stop");
        orf.End.Should().Be(9);
        ((orf.End - orf.Start) % 3).Should().Be(0, "INV-03: the span is codon-aligned despite the internal garbage");
        orf.End.Should().BeLessThanOrEqualTo(seq.Length, "INV-04: the ORF stays within the sequence bounds");
    }

    #endregion

    #region Positive sanity — a clean ORF is found with correct boundaries

    /// <summary>
    /// Positive sanity: alongside the degenerate/malformed probes, a textbook ORF —
    /// "ATG" + a run of coding codons + "TAA" — must be detected with the CORRECT start,
    /// stop, frame, codon-aligned length and translated protein, so the boundary
    /// hardening never silently breaks the core function. We build
    /// ATG (GCT)×6 TAA in frame +1: start at 0, stop closing the 8-codon ORF, span 24 nt
    /// (divisible by 3), 6 internal amino acids (L_aa = 6 ≥ minLength). Pinned per
    /// INV-01..INV-04 (ORF_Detection.md §2.4): starts ATG, ends TAA, span %3 == 0, in
    /// bounds; ProteinSequence includes the terminal '*' stop symbol (ORF_Detection.md
    /// §3.2, §5.2).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindOrfs_CleanOrf_DetectedWithCorrectStartStopFrameAndLength()
    {
        // ATG + (GCT)×6 [Ala] + TAA  → 1 + 6 + 1 = 8 codons, 24 nt, L_aa = 6.
        string coding = string.Concat(Enumerable.Repeat("GCT", 6));
        string seq = "ATG" + coding + "TAA";

        var orf = GenomeAnnotator.FindOrfs(seq, minLength: 6, searchBothStrands: false, requireStartCodon: true)
            .Should().ContainSingle("frame +1 holds exactly one clean ATG…TAA ORF at minLength = 6").Subject;

        orf.Start.Should().Be(0, "the ORF opens at the leading ATG");
        orf.End.Should().Be(seq.Length, "the ORF closes at the terminal TAA, spanning the whole sequence");
        orf.Frame.Should().Be(1, "an ORF starting at index 0 is in reading frame +1");
        orf.IsReverseComplement.Should().BeFalse("this ORF is on the forward strand");
        orf.Sequence.Should().StartWith("ATG", "INV-01: the ORF begins with a start codon");
        orf.Sequence.Should().EndWith("TAA", "INV-02: the ORF ends with a stop codon");
        ((orf.End - orf.Start) % 3).Should().Be(0, "INV-03: the ORF span is divisible by 3");
        (orf.Start >= 0 && orf.Start < orf.End && orf.End <= seq.Length).Should().BeTrue(
            "INV-04: 0 ≤ Start < End ≤ sequence length");
        orf.ProteinSequence.Should().EndWith("*",
            "the canonical translation includes the terminal '*' stop symbol for a stop-terminated ORF");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed results — every reported ORF begins with a start codon,
    /// ends with a stop codon, is codon-aligned and in bounds — so the degenerate-boundary
    /// and malformed-content probes do not corrupt the scan on ordinary input. Pinned per
    /// INV-01..INV-04 (ORF_Detection.md §2.4) with the default requireStartCodon = true.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindOrfs_RandomSequence_ProducesOnlyWellFormedOrfs()
    {
        var starts = new[] { "ATG", "GTG", "TTG" };
        var stops = new[] { "TAA", "TAG", "TGA" };
        string seq = RandomDna(1500, seed: 28_011);

        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true, requireStartCodon: true).ToList();

        orfs.Should().OnlyContain(o =>
            o.Sequence.Length >= 6 &&                                       // start + stop minimum
            starts.Contains(o.Sequence.Substring(0, 3)) &&                  // INV-01
            stops.Contains(o.Sequence.Substring(o.Sequence.Length - 3, 3)) && // INV-02
            (o.End - o.Start) % 3 == 0 &&                                   // INV-03
            o.Start >= 0 && o.Start < o.End && o.End <= seq.Length,         // INV-04
            "every ORF on random input begins with a start, ends with a stop, is codon-aligned, and stays in bounds");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-GENE-001 — gene prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-GENE-001 — gene prediction

    #region Helpers — gene-prediction constructions

    /// <summary>
    /// Builds a forward-strand coding cassette: a leading 5' pad, an in-frame ATG start,
    /// <paramref name="codingCodons"/> repeats of "GCT" (Ala), then a TAA stop. The pad is
    /// "AAA"-tiled so it introduces no start/stop in any frame and keeps the cassette
    /// frame-+1 aligned. Translated length = codingCodons aa (the terminal '*' is trimmed
    /// by PredictGenes for protein_length).
    /// </summary>
    private static string CodingCassette(int codingCodons, int padCodons)
    {
        string pad = string.Concat(Enumerable.Repeat("AAA", padCodons));
        string coding = string.Concat(Enumerable.Repeat("GCT", codingCodons));
        return pad + "ATG" + coding + "TAA";
    }

    #endregion

    #region BE — Boundary: empty sequence yields no genes and no RBS hits

    /// <summary>
    /// BE: an empty sequence must yield NO genes and NO RBS hits, never a crash
    /// (Gene_Prediction.md §3.3, §6.1 row 1). PredictGenes short-circuits through the
    /// underlying FindOrfs `yield break`; FindRibosomeBindingSites short-circuits on its
    /// own IsNullOrEmpty guard. Both are lazy iterators, so we force enumeration to make
    /// the documented empty result actually surface.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_EmptySequence_ReturnsNoGenesAndNoRbs()
    {
        var predict = () => GenomeAnnotator.PredictGenes(string.Empty, minOrfLength: 1).ToList();
        predict.Should().NotThrow("an empty sequence short-circuits through FindOrfs; it simply yields no genes");
        GenomeAnnotator.PredictGenes(string.Empty, minOrfLength: 1).ToList()
            .Should().BeEmpty("no ORF can exist in an empty sequence, so no gene is predicted");

        var rbs = () => GenomeAnnotator.FindRibosomeBindingSites(string.Empty).ToList();
        rbs.Should().NotThrow("the RBS helper `yield break`s on empty input before any indexing");
        GenomeAnnotator.FindRibosomeBindingSites(string.Empty).ToList()
            .Should().BeEmpty("RBS scanning is ORF-driven; no ORF means no upstream motif scan");
    }

    #endregion

    #region BE — Boundary: non-coding-only sequence (no ORF, no RBS) yields nothing

    /// <summary>
    /// BE: a non-coding-only sequence — no start codon in any frame, hence no ORF — must
    /// yield NO genes and NO RBS hits (Gene_Prediction.md §6.1 row 2). We tile "CCC"
    /// (revcomp "GGG"), which is neither a start nor a stop on either strand, so the ORF
    /// scan that BOTH helpers drive finds nothing. The disciplined outcome is the empty
    /// result on both surfaces, not a crash and not a spurious gene.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_NonCodingOnly_ReturnsNoGenesAndNoRbs()
    {
        // 60 nt of pure 'CCC' — no ATG/GTG/TTG and no TAA/TAG/TGA on either strand.
        string nonCoding = string.Concat(Enumerable.Repeat("CCC", 20));

        var predict = () => GenomeAnnotator.PredictGenes(nonCoding, minOrfLength: 1).ToList();
        predict.Should().NotThrow("a non-coding sequence records no pending ORF start; the scan terminates cleanly");
        GenomeAnnotator.PredictGenes(nonCoding, minOrfLength: 1).ToList()
            .Should().BeEmpty("with no start codon on either strand there is no ORF and therefore no gene");

        GenomeAnnotator.FindRibosomeBindingSites(nonCoding).ToList()
            .Should().BeEmpty("RBS scanning is ORF-driven, so a sequence with no ORF yields no RBS hits");
    }

    #endregion

    #region MC — Malformed content: no RBS motif anywhere upstream of a real gene

    /// <summary>
    /// MC: a sequence that DOES contain a genuine, qualifying ORF but has NO Shine-Dalgarno
    /// motif (AGGAGG / GGAGG / AGGAG / GAGG / AGGA) anywhere upstream must still predict the
    /// gene, yet report NO RBS hit. This pins the documented INDEPENDENCE of the two helpers
    /// (Gene_Prediction.md §5.2, §5.4 deviation 1): "no RBS in seq" suppresses RBS hits, NOT
    /// gene calls. The ORF is ≥ 30 aa so the RBS helper's internal FindOrfs(minLength: 30)
    /// discovers it and actually scans its upstream window; the pad is "AAA"-tiled so the
    /// upstream region contains no SD-like motif. The gene must appear with valid bounds; the
    /// RBS list must be empty.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_NoRbsMotifUpstream_PredictsGeneButReportsNoRbs()
    {
        // 35 coding codons (≥ 30 aa) with a 4-codon 'AAA' pad — no SD-like motif upstream.
        string seq = CodingCassette(codingCodons: 35, padCodons: 4);

        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();
        genes.Should().ContainSingle("the cassette holds exactly one qualifying forward-strand ORF")
            .Which.Strand.Should().Be('+', "the gene is on the forward strand");

        var rbs = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        rbs.Should().BeEmpty(
            "the upstream window holds no Shine-Dalgarno motif, so the RBS helper invents nothing — " +
            "yet the gene is still predicted because PredictGenes ignores RBS evidence entirely");
    }

    #endregion

    #region BE — Boundary: overlapping genes are all reported with valid coordinates

    /// <summary>
    /// BE: overlapping / nested ORFs must each be emitted as their OWN CDS — there is NO
    /// best-model selection or overlap suppression (Gene_Prediction.md §5.2 bullet 2,
    /// §5.3 "Intentionally simplified", §6.2). The disciplined boundary here is the overlap
    /// BOOKKEEPING: predicting multiple coincident genes must not crash, and EVERY emitted
    /// gene must keep valid coordinates (GENE-INV-02: 0 ≤ Start &lt; End ≤ length,
    /// Frame ∈ {1,2,3}, Strand ∈ {'+','-'}, Type == "CDS"), with IDs ordered by Start.
    ///
    /// We build a single coding cassette that is a palindromic-by-construction overlap
    /// generator: a forward ORF AND its reverse complement both score as ORFs over the SAME
    /// genomic span. Concretely "ATG...TAA" on the forward strand revcomps to "TTA...CAT";
    /// to force a genuine overlapping pair we instead concatenate two forward ORFs whose
    /// spans abut/overlap by sharing coordinates across frames. Simpler and fully
    /// deterministic: a cassette with an inner in-frame ATG produces a nested ORF (a second
    /// start before the shared stop), so PredictGenes emits BOTH the outer and the inner
    /// gene over overlapping coordinates.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_OverlappingGenes_AllReportedWithValidCoordinates()
    {
        // Outer ORF: ATG (GCT)x40 ATG (GCT)x5 TAA  → an INNER in-frame ATG creates a nested
        // ORF that shares the same TAA stop, so two overlapping genes are emitted in frame +1.
        string outerHead = "ATG" + string.Concat(Enumerable.Repeat("GCT", 40));
        string innerHead = "ATG" + string.Concat(Enumerable.Repeat("GCT", 5));
        string seq = outerHead + innerHead + "TAA";

        var predict = () => GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();
        predict.Should().NotThrow("emitting multiple overlapping ORFs must not crash the overlap bookkeeping");

        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();

        genes.Count.Should().BeGreaterThanOrEqualTo(2,
            "both the outer ORF and the nested inner ORF (sharing the terminal stop) are emitted; " +
            "there is no overlap suppression");

        genes.Should().OnlyContain(g =>
            g.Start >= 0 && g.Start < g.End && g.End <= seq.Length &&  // GENE-INV-02 bounds
            (g.Strand == '+' || g.Strand == '-') &&                    // GENE-INV-02 strand
            g.Type == "CDS",                                           // GENE-INV-02 type
            "every overlapping gene keeps valid coordinates, a ± strand, and Type=CDS");

        // GENE-INV-02 frame ∈ {1,2,3} — carried in the 'frame' attribute, not a record field.
        genes.Should().OnlyContain(g =>
            g.Attributes.ContainsKey("frame") &&
            (g.Attributes["frame"] == "1" || g.Attributes["frame"] == "2" || g.Attributes["frame"] == "3"),
            "every overlapping gene tags a reading frame in {1,2,3}");

        // Genes are ordered by genomic Start, and two distinct overlapping starts exist.
        genes.Select(g => g.Start).Should().BeInAscendingOrder("PredictGenes orders genes by genomic Start");
        genes.Select(g => g.Start).Distinct().Count().Should().BeGreaterThanOrEqualTo(2,
            "the outer and inner ORFs open at distinct positions yet overlap, confirming no overlap was suppressed");
    }

    #endregion

    #region Positive sanity — RBS upstream of a downstream ORF yields a gene + RBS hit

    /// <summary>
    /// Positive sanity: a textbook prokaryotic gene — a Shine-Dalgarno motif (AGGAGG) an
    /// aligned short distance upstream of an in-frame ATG that opens a ≥ 30 aa ORF — must
    /// produce BOTH a predicted gene with correct coordinates AND a matching RBS hit, so the
    /// boundary/malformed probes never silently break the core function.
    ///
    /// Layout (frame +1): [pad AAA] AGGAGG [8-nt 'AAAAAAAA' spacer] ATG (GCT)x35 TAA.
    /// The SD motif 3' end sits 8 nt upstream of ATG → aligned distance 8 ∈ [4,15]
    /// (GENE-INV-03), inside the default 20-nt upstream window. The ORF is 35 aa ≥ the RBS
    /// helper's internal 30-aa floor, so it is scanned. Pinned: the gene starts at the ATG
    /// with End at the TAA, valid bounds, '+' strand, CDS; the RBS hit reports AGGAGG with
    /// score 1.0 (= 6/6, GENE-INV-04) at the SD position.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void PredictGenes_RbsUpstreamOfOrf_YieldsGeneWithCorrectCoordsAndRbsHit()
    {
        const string pad = "AAA";               // 3 nt 5' pad (no start/stop, keeps frame)
        const string sd = "AGGAGG";             // full Shine-Dalgarno consensus, score 1.0
        const string spacer = "AAAAAAAA";       // 8 nt → aligned distance 8 ∈ [4,15]
        string coding = "ATG" + string.Concat(Enumerable.Repeat("GCT", 35)) + "TAA"; // 35 aa ORF
        string seq = pad + sd + spacer + coding;

        int expectedStart = pad.Length + sd.Length + spacer.Length; // start of ATG
        int expectedSdPos = pad.Length;                             // start of AGGAGG

        // ── Gene prediction ──────────────────────────────────────────────────────────
        var gene = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30)
            .Should().ContainSingle("exactly one qualifying forward-strand ORF is present").Subject;

        gene.Start.Should().Be(expectedStart, "the gene opens at the in-frame ATG downstream of the SD motif");
        gene.End.Should().Be(seq.Length, "the gene closes at the terminal TAA");
        gene.Strand.Should().Be('+', "the gene is on the forward strand");
        gene.Type.Should().Be("CDS", "PredictGenes emits CDS annotations");
        gene.Attributes["frame"].Should().BeOneOf("1", "2", "3", "GENE-INV-02: frame is in {1,2,3}");
        (gene.Start >= 0 && gene.Start < gene.End && gene.End <= seq.Length).Should().BeTrue(
            "GENE-INV-02: 0 ≤ Start < End ≤ sequence length");

        // ── RBS detection ────────────────────────────────────────────────────────────
        var rbs = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        rbs.Should().NotBeEmpty("a full AGGAGG motif sits at an aligned distance of 8 nt upstream of the start");

        var hit = rbs.Should().Contain(h => h.sequence == "AGGAGG",
            "the full Shine-Dalgarno consensus is the highest-scoring motif present").Subject;
        hit.position.Should().Be(expectedSdPos, "the RBS hit reports the genomic 5' position of the SD motif");
        hit.score.Should().BeApproximately(1.0, 1e-9, "GENE-INV-04: AGGAGG (6 nt) scores 6/6 = 1.0");

        int alignedDistance = gene.Start - hit.position - hit.sequence.Length;
        alignedDistance.Should().BeInRange(4, 15, "GENE-INV-03: aligned spacing lies within [minDistance, maxDistance]");
    }

    /// <summary>
    /// Positive sanity / robustness: on a fixed-seed random 1200-nt sequence both helpers
    /// must complete promptly and every result must be well-formed — every predicted gene
    /// satisfies GENE-INV-02 (valid bounds, frame ∈ {1,2,3}, ± strand, CDS) and every RBS
    /// hit satisfies GENE-INV-03/04 (aligned spacing in [4,15], score ∈ [4/6, 1]). This is
    /// the "degenerate-boundary probes do not corrupt the scan on ordinary input" guard.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void PredictGenes_RandomSequence_ProducesOnlyWellFormedGenesAndRbsHits()
    {
        string seq = RandomDna(1200, seed: 29_001);

        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();
        genes.Should().OnlyContain(g =>
            g.Start >= 0 && g.Start < g.End && g.End <= seq.Length &&
            (g.Strand == '+' || g.Strand == '-') &&
            g.Type == "CDS" &&
            g.Attributes.ContainsKey("frame") &&
            (g.Attributes["frame"] == "1" || g.Attributes["frame"] == "2" || g.Attributes["frame"] == "3"),
            "every predicted gene on random input is in bounds, codon-frame-tagged, ± stranded, and CDS");

        var rbs = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        rbs.Should().OnlyContain(h =>
            h.score >= 4.0 / 6.0 && h.score <= 1.0 + 1e-9 &&   // GENE-INV-04
            h.position >= 0 && h.position < seq.Length,        // in-bounds motif position
            "every RBS hit on random input carries a length-normalized score in [4/6, 1] at an in-bounds position");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-PROM-001 — promoter (−35 / −10 box) motif detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-PROM-001 — promoter prediction

    #region Helpers — promoter motif libraries

    /// <summary>The fixed −35 / −10 variant library scanned by FindPromoterMotifs (Promoter_Detection.md §4.2).</summary>
    private static readonly string[] PromoterMotifLibrary =
        { "TTGACA", "TTGAC", "TGACA", "TTGA", "TATAAT", "TATAA", "ATAAT", "TATA" };

    #endregion

    #region BE — Boundary: no -10 box present in the sequence

    /// <summary>
    /// BE: a sequence containing NO −10 box (none of TATAAT / TATAA / ATAAT / TATA) must
    /// emit NO −10 hits — matching is EXACT against the fixed library, so the absence of any
    /// −10 variant invents nothing (Promoter_Detection.md §3.3, §6.1 "No supported motifs",
    /// PROM-INV-03). The two motif libraries are scanned INDEPENDENTLY, so this is asserted
    /// specifically for the −10 class; the sequence here is built from G/C runs that contain
    /// no −10 variant (and, as it happens, no −35 variant either), so the whole result is
    /// empty — no crash, no spurious motif.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_NoMinusTenBox_EmitsNoMinusTenHits()
    {
        // Pure G/C run: contains no TATAAT/TATAA/ATAAT/TATA and no TTGACA-family variant.
        const string seq = "GGGGCCCCGGGGCCCCGGGGCCCC";

        var act = () => GenomeAnnotator.FindPromoterMotifs(seq).ToList();
        act.Should().NotThrow("an exact-match scan over a sequence with no library motif simply yields nothing");

        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        motifs.Where(m => m.type == "-10 box").Should().BeEmpty(
            "no −10 variant occurs in a G/C-only sequence, so the exact-match scan reports no Pribnow box");
        motifs.Should().BeEmpty(
            "neither library motif occurs, so the helper invents no spurious promoter hit on either class");
    }

    /// <summary>
    /// BE: the independence of the two libraries — a sequence carrying a genuine −35 box but
    /// NO −10 box must still report the −35 hits while reporting ZERO −10 hits. This pins that
    /// "no −10 box" suppresses ONLY the −10 class, not the whole scan (Promoter_Detection.md
    /// §5.2 "scanned independently"). We embed TTGACA in a G/C carrier that contains no −10
    /// variant; the −35 family (TTGACA, TTGAC, TGACA, TTGA) must appear and the −10 list must
    /// be empty.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_Minus35ButNoMinusTen_ReportsMinus35OnlyNoMinusTen()
    {
        // TTGACA flanked by G/C only — a real −35 box, but no −10 variant anywhere.
        const string seq = "GGGGGTTGACAGCGCGC";

        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        motifs.Where(m => m.type == "-35 box").Should().NotBeEmpty(
            "the genuine TTGACA box yields −35 hits even though no −10 box is present");
        motifs.Where(m => m.type == "-10 box").Should().BeEmpty(
            "the libraries are scanned independently, so the absence of a −10 box suppresses only the −10 class");
    }

    #endregion

    #region BE — Boundary: "threshold = 0" is inapplicable — no score gate exists

    /// <summary>
    /// BE: the checklist target "threshold = 0" is DEGENERATE for this surface —
    /// FindPromoterMotifs exposes NO threshold parameter and applies NO score gate
    /// (Promoter_Detection.md §3.1, §5.2). The theory-correct consequence is that the method
    /// ALWAYS emits every exact motif hit at its full hard-coded score, and EVERY emitted
    /// score is strictly positive (∈ [0.665, 1.0]); a hypothetical "threshold of 0" would
    /// admit exactly that same set (no hit has score &lt; 0). We pin the falsifiable contract:
    /// on a sequence carrying both full consensus boxes, every hit's score is &gt; 0 and ≤ 1,
    /// no hit is suppressed or zero-scored, so there is no degenerate-threshold blow-up or
    /// gate to mis-set.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_NoThresholdGate_AllHitsScorePositiveAndUnsuppressed()
    {
        // Both full consensus boxes present, independently — exercises every score path.
        const string seq = "GGGTTGACAGGGGGGGGGGGGGGGGGTATAATGGG";

        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        motifs.Should().NotBeEmpty("both full consensus boxes are present, so motif hits are emitted");
        motifs.Should().OnlyContain(m => m.score > 0.0 && m.score <= 1.0 + 1e-9,
            "PROM-INV-01 / no score gate: every emitted hit carries a strictly positive score in (0, 1], " +
            "so a degenerate 'threshold = 0' admits exactly the hits the gate-less scan already emits");

        // Both full consensus hits are present and unsuppressed at score 1.0 (PROM-INV-02).
        motifs.Should().Contain(m => m.sequence == "TTGACA" && Math.Abs(m.score - 1.0) < 1e-9,
            "PROM-INV-02: the full −35 consensus is emitted at score 1.0, never gated away");
        motifs.Should().Contain(m => m.sequence == "TATAAT" && Math.Abs(m.score - 1.0) < 1e-9,
            "PROM-INV-02: the full −10 consensus is emitted at score 1.0, never gated away");
    }

    #endregion

    #region BE — Boundary: empty sequence yields no motifs

    /// <summary>
    /// BE: an empty sequence must yield NO motifs and never crash — for every motif the scan
    /// bound `0 &lt;= 0 - motif.Length` is false, so no loop body executes and no Substring is
    /// taken (Promoter_Detection.md §6.1 "Empty sequence"). FindPromoterMotifs is a lazy
    /// iterator, so we force enumeration to make the documented empty result actually surface.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_EmptySequence_ReturnsNoMotifs()
    {
        var act = () => GenomeAnnotator.FindPromoterMotifs(string.Empty).ToList();
        act.Should().NotThrow("an empty sequence gives every scan loop a negative upper bound; no Substring is taken");

        GenomeAnnotator.FindPromoterMotifs(string.Empty).ToList()
            .Should().BeEmpty("no motif can occur in an empty sequence, so the scan yields nothing");
    }

    /// <summary>
    /// BE (documented null contract): null input is NOT guarded — FindPromoterMotifs
    /// dereferences `dnaSequence` immediately via ToUpperInvariant() (Promoter_Detection.md
    /// §3.3, §6.1 "Null input — not explicitly handled before uppercasing"). The disciplined
    /// fuzz outcome here is to PIN the documented throw, not to weaken the test to accept a
    /// silent empty result the contract does not promise. Because the throw fires inside the
    /// lazy iterator body, enumeration is forced to surface it.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_NullSequence_ThrowsPerDocumentedContract()
    {
        var act = () => GenomeAnnotator.FindPromoterMotifs(null!).ToList();
        act.Should().Throw<NullReferenceException>(
            "the method dereferences the sequence via ToUpperInvariant() with no null guard (Promoter_Detection.md §6.1)");
    }

    #endregion

    #region BE — Boundary: sequence shorter than the motif (KEY: no IndexOutOfRange)

    /// <summary>
    /// BE (KEY boundary): a sequence SHORTER than every library motif must yield NO motifs and
    /// — critically — must NOT raise IndexOutOfRange. The scan is
    /// `for (int i = 0; i &lt;= seq.Length - motif.Length; i++)`; when seq.Length &lt; motif.Length
    /// the upper bound `seq.Length - motif.Length` is NEGATIVE, so the loop body never runs and
    /// `seq.Substring(i, motif.Length)` is never reached past the end (Promoter_Detection.md
    /// §6.1, PROM-INV-04). The shortest library motif is 4 nt (TATA / TTGA), so a 3-nt sequence
    /// is shorter than EVERY motif — the hardest short-seq case.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_SequenceShorterThanMotif_ReturnsEmptyNoIndexOutOfRange()
    {
        // 3 nt: shorter than the shortest library motif (4 nt), so no motif window fits.
        const string seq = "TAT";

        var act = () => GenomeAnnotator.FindPromoterMotifs(seq).ToList();
        act.Should().NotThrow<IndexOutOfRangeException>(
            "a sequence shorter than every motif gives a negative scan bound; the Substring window is never taken");
        act.Should().NotThrow<ArgumentOutOfRangeException>(
            "no Substring(i, len) ever runs past the end of an undersized sequence");

        GenomeAnnotator.FindPromoterMotifs(seq).ToList()
            .Should().BeEmpty("no library motif (min length 4) can fit in a 3-nt sequence");
    }

    /// <summary>
    /// BE: the boundary right AT the shortest motif length — a sequence whose length exactly
    /// equals the 4-nt minimum motif. Here the scan bound is `i &lt;= 0`, so the single window at
    /// i = 0 IS evaluated: a 4-nt sequence that exactly IS a 4-nt library motif must match (one
    /// in-bounds hit), and a 4-nt non-motif must yield nothing — both without running off the
    /// end. This pins that the negative-bound guard does not also suppress the legitimate
    /// length-exact window (PROM-INV-04).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_SequenceEqualsShortestMotifLength_MatchesInBounds()
    {
        // Exactly the 4-nt −10 variant "TATA": the single i = 0 window must match in bounds.
        var hits = GenomeAnnotator.FindPromoterMotifs("TATA").ToList();
        hits.Should().ContainSingle("a 4-nt sequence that exactly is the TATA variant yields one in-bounds hit");
        var hit = hits.Single();
        hit.sequence.Should().Be("TATA");
        hit.type.Should().Be("-10 box");
        hit.position.Should().Be(0, "PROM-INV-04: the only fitting window starts at index 0");

        // A 4-nt non-motif must yield nothing — still no out-of-range access.
        GenomeAnnotator.FindPromoterMotifs("GGGG").ToList()
            .Should().BeEmpty("a 4-nt sequence matching no library motif yields no hit and no out-of-range Substring");
    }

    #endregion

    #region Positive sanity — consensus boxes detected with correct coords and score

    /// <summary>
    /// Positive sanity: a textbook bacterial promoter layout — a full −35 consensus (TTGACA)
    /// upstream of a full −10 consensus (TATAAT) at the classical ~17 bp spacing — must be
    /// detected with the CORRECT 0-based positions and full score 1.0 for both boxes, so the
    /// boundary hardening never silently breaks the core function. Layout (G filler):
    ///   [5 G] TTGACA [17 G spacer] TATAAT [3 G]
    /// The −35 box opens at index 5; the −10 box opens at 5 + 6 + 17 = 28. Pinned per
    /// PROM-INV-01..04 (Promoter_Detection.md §2.4): both full hits score 1.0, every hit's
    /// score ∈ [0,1] and belongs to the fixed library, and every position is in bounds. The
    /// spacing here is biologically faithful but NOT validated by the helper (it pairs nothing);
    /// it merely places two independent, individually-correct hits.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPromoterMotifs_ConsensusBoxesAtCanonicalSpacing_DetectedWithCorrectCoordsAndScore()
    {
        const string lead = "GGGGG";                       // 5 nt 5' filler (no motif)
        const string minus35 = "TTGACA";                   // full −35 consensus
        const string spacer = "GGGGGGGGGGGGGGGGG";          // 17 nt classical spacing
        const string minus10 = "TATAAT";                   // full −10 consensus (Pribnow)
        const string tail = "GGG";
        string seq = lead + minus35 + spacer + minus10 + tail;

        int expected35Pos = lead.Length;                              // 5
        int expected10Pos = lead.Length + minus35.Length + spacer.Length; // 28

        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        // ── −35 box ──────────────────────────────────────────────────────────────────
        var full35 = motifs.Should().ContainSingle(m => m.sequence == "TTGACA",
            "the full −35 consensus occurs exactly once").Subject;
        full35.type.Should().Be("-35 box");
        full35.position.Should().Be(expected35Pos, "the −35 box opens at the 0-based index after the 5-nt lead");
        full35.score.Should().BeApproximately(1.0, 1e-9, "PROM-INV-02: a full TTGACA match scores 1.0");

        // ── −10 box ──────────────────────────────────────────────────────────────────
        var full10 = motifs.Should().ContainSingle(m => m.sequence == "TATAAT",
            "the full −10 consensus occurs exactly once").Subject;
        full10.type.Should().Be("-10 box");
        full10.position.Should().Be(expected10Pos, "the −10 box opens 17 nt downstream of the −35 box end");
        full10.score.Should().BeApproximately(1.0, 1e-9, "PROM-INV-02: a full TATAAT match scores 1.0");

        // ── Whole-result invariants ──────────────────────────────────────────────────
        motifs.Should().OnlyContain(m =>
            m.score >= 0.0 && m.score <= 1.0 + 1e-9 &&                       // PROM-INV-01
            (m.type == "-35 box" || m.type == "-10 box") &&                  // type domain
            PromoterMotifLibrary.Contains(m.sequence) &&                     // PROM-INV-03
            m.position >= 0 && m.position <= seq.Length - m.sequence.Length, // PROM-INV-04
            "every hit scores in [0,1], names a known box type, is from the fixed library, and sits in bounds");
    }

    /// <summary>
    /// Positive sanity / robustness: on a fixed-seed random 1500-nt sequence the scan must
    /// complete promptly and produce ONLY well-formed hits — every reported motif belongs to
    /// the fixed library, scores in [0,1], names a valid box type, and sits at an in-bounds
    /// position whose motif window stays within the sequence (PROM-INV-01/03/04). This is the
    /// "degenerate-boundary probes do not corrupt the scan on ordinary input" guard. We also
    /// cross-check that the substring at each reported position actually equals the reported
    /// motif, so no position is ever mis-reported.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindPromoterMotifs_RandomSequence_ProducesOnlyWellFormedHits()
    {
        string seq = RandomDna(1500, seed: 30_001);

        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        motifs.Should().OnlyContain(m =>
            (m.type == "-35 box" || m.type == "-10 box") &&                  // type domain
            PromoterMotifLibrary.Contains(m.sequence) &&                     // PROM-INV-03
            m.score >= 0.0 && m.score <= 1.0 + 1e-9 &&                       // PROM-INV-01
            m.position >= 0 && m.position <= seq.Length - m.sequence.Length && // PROM-INV-04
            seq.Substring(m.position, m.sequence.Length) == m.sequence,      // position fidelity
            "every hit on random input is a known-library motif, scored in [0,1], in bounds, " +
            "and the reported position genuinely contains the reported motif");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-GFF-001 — GFF3 I/O (parser + serializer) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-GFF-001 — GFF annotation

    #region Helpers — GFF3 line construction

    /// <summary>
    /// Builds one tab-separated GFF3 data line from explicit column values, so a single
    /// injected/malformed column can be probed without ambiguity. Columns are joined with
    /// the literal tab GFF3 requires.
    /// </summary>
    private static string GffLine(
        string seqid, string source, string type, string start, string end,
        string score, string strand, string phase, string attributes) =>
        string.Join('\t', seqid, source, type, start, end, score, strand, phase, attributes);

    /// <summary>A well-formed 9-column GFF3 data line (valid in every column).</summary>
    private static string ValidGffLine() =>
        GffLine("seq1", "prog", "gene", "10", "99", ".", "+", ".", "ID=g1;Name=alpha");

    #endregion

    #region TF — Truncated / missing fields: a sub-9-column line is skipped, not crashed

    /// <summary>
    /// TF (missing fields): a truncated data line with FEWER than 9 tab-separated columns must
    /// be SKIPPED, never crash with IndexOutOfRange (GFF3_IO.md §4.2 "Fewer than 9 fields →
    /// skipped", §6.1). The eligibility gate `parts.Length &lt; 9` drops the line before any
    /// column is read, so no `parts[3]`/`parts[6]` access ever runs past the truncated array.
    /// We feed a 5-column fragment AMONG a valid line so the disciplined outcome is provable:
    /// the truncated line contributes nothing while the valid line still parses.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ParseGff3_TruncatedLine_IsSkippedNotIndexOutOfRange()
    {
        string truncated = "seq1\tprog\tgene\t10\t99"; // only 5 of 9 columns
        var lines = new[] { truncated, ValidGffLine() };

        var act = () => GenomeAnnotator.ParseGff3(lines).ToList();
        act.Should().NotThrow<IndexOutOfRangeException>(
            "a line with fewer than 9 fields is dropped by the `< 9` gate before any column index is read");

        var feats = GenomeAnnotator.ParseGff3(lines).ToList();
        feats.Should().ContainSingle("only the well-formed 9-column line survives; the 5-column fragment is skipped")
            .Which.FeatureId.Should().Be("g1", "the surviving feature is the valid line, identified by its ID attribute");
    }

    /// <summary>
    /// TF (present-but-empty strand column): a structurally-9-column line whose STRAND column is
    /// the empty string was a caller-safety crash — `parts[6][0]` raised IndexOutOfRange on the
    /// empty field. After the caller-safety hardening (mirroring GffParser.ParseLine) an empty
    /// strand column defaults to the valid '.' no-strand sentinel: the line parses, no throw.
    /// We pin BOTH the non-throw and the normalized '.' strand so the hardening cannot silently
    /// drift back to an indexing crash.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ParseGff3_EmptyStrandColumn_DefaultsToDotNotIndexOutOfRange()
    {
        string emptyStrand = GffLine("seq1", "prog", "gene", "10", "99", ".", "", ".", "ID=g1");

        var act = () => GenomeAnnotator.ParseGff3(new[] { emptyStrand }).ToList();
        act.Should().NotThrow<IndexOutOfRangeException>(
            "an empty (present-but-blank) strand column must not index past the empty field");

        var feat = GenomeAnnotator.ParseGff3(new[] { emptyStrand }).ToList()
            .Should().ContainSingle("a 9-column line with a blank strand still parses").Subject;
        feat.Strand.Should().Be('.', "an empty strand column defaults to the GFF3 '.' no-strand sentinel");
    }

    #endregion

    #region MC — Malformed content: non-numeric numeric columns are skipped, not thrown

    /// <summary>
    /// MC: a 9-column line whose START column is non-numeric ("ABC") was an uncaught
    /// FormatException (int.Parse on garbage). After hardening the malformed line is SKIPPED —
    /// mirroring the existing &lt;9-field skip discipline — so a single corrupt data line cannot
    /// abort the whole parse (GFF3_IO.md §6.1; 03_FUZZING.md §Description MC). We interleave the
    /// garbage line with a valid one and assert the valid line still parses while the garbage
    /// line contributes nothing and does NOT consume a `feature_n` id slot.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ParseGff3_NonNumericStart_LineIsSkippedNotFormatException()
    {
        string badStart = GffLine("seq1", "prog", "gene", "ABC", "99", ".", "+", ".", "ID=g1");
        // No ID on the valid line, so its auto feature id reveals whether the skipped line bumped the counter.
        string valid = GffLine("seq1", "prog", "gene", "10", "99", ".", "+", ".", "Name=alpha");

        var act = () => GenomeAnnotator.ParseGff3(new[] { badStart, valid }).ToList();
        act.Should().NotThrow<FormatException>("a non-numeric start column must skip the line, not throw");

        var feats = GenomeAnnotator.ParseGff3(new[] { badStart, valid }).ToList();
        feats.Should().ContainSingle("the malformed-start line is skipped; only the valid line survives");
        feats[0].FeatureId.Should().Be("feature_1",
            "the skipped malformed line did NOT advance the feature counter, so the first surviving feature is feature_1");
    }

    /// <summary>
    /// MC: non-numeric SCORE and PHASE columns are likewise skipped, never thrown. Score and
    /// phase use '.' as the documented null sentinel (GFF3_IO.md §4.2); any OTHER non-numeric
    /// token is malformed content. We assert each malformed line is dropped without a
    /// FormatException, and that a genuine '.' sentinel still parses to a null Score/Phase.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ParseGff3_NonNumericScoreOrPhase_LinesSkippedAndDotSentinelStillNull()
    {
        string badScore = GffLine("seq1", "prog", "gene", "10", "99", "XYZ", "+", ".", "ID=bad1");
        string badPhase = GffLine("seq1", "prog", "CDS", "10", "99", ".", "+", "ZZ", "ID=bad2");
        string sentinel = GffLine("seq1", "prog", "gene", "10", "99", ".", "+", ".", "ID=ok");

        var act = () => GenomeAnnotator.ParseGff3(new[] { badScore, badPhase, sentinel }).ToList();
        act.Should().NotThrow<FormatException>("malformed score/phase tokens skip the line; they never throw");

        var feats = GenomeAnnotator.ParseGff3(new[] { badScore, badPhase, sentinel }).ToList();
        var ok = feats.Should().ContainSingle("only the '.'-sentinel line is well-formed; both garbage lines are skipped").Subject;
        ok.FeatureId.Should().Be("ok");
        ok.Score.Should().BeNull("a '.' score column parses to the null sentinel, not 0");
        ok.Phase.Should().BeNull("a '.' phase column parses to the null sentinel, not 0");
    }

    /// <summary>
    /// MC (invalid strand symbol): an out-of-domain strand 'Z' (not in {+,-,.,?}) is NOT rejected
    /// by the lightweight reader — column 7 is parsed VERBATIM (GFF3_IO.md §3.2). The disciplined
    /// fuzz outcome is a deterministic verbatim parse, never a crash and never a silent drop. We
    /// pin that the line parses and the strand is exactly the supplied char, so the documented
    /// "no strand validation" contract cannot silently change to a rejection without this test
    /// failing.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ParseGff3_InvalidStrandSymbol_ParsedVerbatimNotRejected()
    {
        string invalidStrand = GffLine("seq1", "prog", "gene", "10", "99", ".", "Z", ".", "ID=g1");

        var feat = GenomeAnnotator.ParseGff3(new[] { invalidStrand }).ToList()
            .Should().ContainSingle("the lightweight reader applies no strand-domain validation").Subject;
        feat.Strand.Should().Be('Z',
            "column 7 is parsed verbatim, so an out-of-domain strand symbol is preserved exactly, not rejected or normalized");
    }

    #endregion

    #region INJ — Tab injection on the read side: extra tab cannot smuggle a tail fragment

    /// <summary>
    /// INJ (tab injection, read side): an EXTRA literal tab inside the attributes field produces
    /// MORE than 9 columns. Because the eligibility gate is `parts.Length &lt; 9` (not `!= 9`), the
    /// line is NOT skipped — the reader consumes columns 0..8 and the injected-tab tail beyond
    /// column 9 is DROPPED. The disciplined contract is that this is deterministic and crash-free,
    /// and critically that the smuggled tail fragment can NOT leak into the parsed attribute
    /// dictionary nor shift columns. We inject a tab after the real attributes followed by a
    /// fake "ID=evil" tail; the parsed attributes must contain ONLY the legitimate column-9
    /// attribute and must NOT pick up the post-tab "evil" id.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ParseGff3_TabInjectionInAttributes_DropsTailAndDoesNotCorruptColumns()
    {
        // 10 columns: a tab injected after the genuine attributes, with a fake ID tail.
        string injected = "seq1\tprog\tgene\t10\t99\t.\t+\t.\tID=real\tID=evil";

        var act = () => GenomeAnnotator.ParseGff3(new[] { injected }).ToList();
        act.Should().NotThrow("a >9-column line is parsed deterministically from columns 0..8; the tail is ignored");

        var feat = GenomeAnnotator.ParseGff3(new[] { injected }).ToList()
            .Should().ContainSingle("the injected tab does not split the record into multiple features").Subject;
        feat.FeatureId.Should().Be("real",
            "the genuine column-9 ID is used; the post-tab 'ID=evil' tail lies beyond column 9 and is dropped");
        feat.Attributes.Should().ContainKey("ID")
            .WhoseValue.Should().Be("real", "no smuggled attribute from the injected tail reaches the dictionary");
    }

    #endregion

    #region INJ — Tab/special-char injection on the WRITE side: serializer must escape, not break

    /// <summary>
    /// INJ (serializer, the high-value target): an attribute value stuffed with the exact
    /// characters that would break GFF3 column structure — TAB, newline, ';', '=', '&', ',', '%'
    /// — must serialize to a SINGLE 9-column row, NEVER a structurally broken line. ToGff3 percent-
    /// encodes every column-9 reserved character via EncodeGff3Value (GFF3_IO.md §4.2): tab→%09,
    /// newline→%0A, ';'→%3B, '='→%3D, '&'→%26, ','→%2C, '%'→%25. We pin that the emitted data row
    /// splits into exactly 9 tab columns (so the injected tab did NOT create a 10th column) and
    /// that the raw column-9 text carries no literal tab/newline/';'/'=' from the value.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ToGff3_AttributeValueWithGffBreakingChars_EscapedToSingleNineColumnRow()
    {
        var ann = new GenomeAnnotator.GeneAnnotation(
            GeneId: "g\t1",                 // tab injected into the ID itself
            Start: 0, End: 10, Strand: '+', Type: "CDS",
            Product: "p;a=b\nx",            // ';', '=', newline in the product
            Attributes: new Dictionary<string, string> { ["note"] = "has\ttab,comma&amp%pct" });

        var lines = GenomeAnnotator.ToGff3(new[] { ann }, "seqA").ToList();

        lines[0].Should().Be("##gff-version 3", "GFF3 output opens with the version directive");
        lines.Should().HaveCount(2, "one header line plus one feature row");

        string row = lines[1];
        row.Split('\t').Should().HaveCount(9,
            "the injected tab is percent-encoded (%09), so the row stays exactly 9 tab-separated columns and is not structurally broken");

        string col9 = row.Split('\t')[8];
        col9.Should().NotContain("\n", "a newline in an attribute value is escaped to %0A, never emitted raw");
        col9.Should().Contain("%09", "the injected tab is escaped to %09 inside column 9");
        col9.Should().Contain("%3B", "the ';' in the product value is escaped to %3B so it cannot start a spurious attribute");
        col9.Should().Contain("%26", "the '&' is escaped to %26");
        col9.Should().Contain("%2C", "the ',' is escaped to %2C");
        col9.Should().Contain("%25", "the literal '%' is escaped to %25 so encoding is unambiguous");
    }

    #endregion

    #region Positive sanity — a valid annotation serializes to 9 columns and round-trips

    /// <summary>
    /// Positive sanity: a clean GeneAnnotation must serialize to a well-formed 9-column GFF3 row
    /// AND round-trip back through ParseGff3 with faithful coordinates and attributes, so the
    /// malformed/injection hardening never silently breaks the core function. Pinned per the GFF3
    /// invariants (GFF3_IO.md §2.4): GFF-INV-01 the external start is 1-based (Start + 1);
    /// GFF-INV-03 phase is '0' for a CDS; and an attribute value carrying a tab round-trips
    /// byte-exact (escape on write, percent-decode on read).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GffRoundTrip_ValidAnnotation_NineColumnsAndFaithfulRoundTrip()
    {
        var ann = new GenomeAnnotator.GeneAnnotation(
            GeneId: "g1", Start: 0, End: 30, Strand: '+', Type: "CDS", Product: "alpha kinase",
            Attributes: new Dictionary<string, string> { ["note"] = "embedded\ttab" });

        var lines = GenomeAnnotator.ToGff3(new[] { ann }, "chr1").ToList();
        string row = lines[1];
        var cols = row.Split('\t');

        cols.Should().HaveCount(9, "a valid annotation serializes to exactly 9 GFF3 columns");
        cols[0].Should().Be("chr1", "column 1 is the supplied seqId");
        cols[2].Should().Be("CDS", "column 3 is the annotation type");
        cols[3].Should().Be("1", "GFF-INV-01: column 4 start is 1-based (internal 0-based Start 0 → '1')");
        cols[4].Should().Be("30", "column 5 end is written unchanged");
        cols[6].Should().Be("+", "column 7 carries the annotation strand");
        cols[7].Should().Be("0", "GFF-INV-03: phase is '0' for a CDS feature");

        // ── Round-trip through the lightweight reader ────────────────────────────────────
        var feat = GenomeAnnotator.ParseGff3(lines).ToList()
            .Should().ContainSingle("the single serialized feature parses back to one record").Subject;
        feat.Type.Should().Be("CDS");
        feat.Start.Should().Be(ann.Start + 1, "the reader preserves the 1-based start written by the exporter (Start + 1)");
        feat.End.Should().Be(ann.End, "the end coordinate round-trips unchanged");
        feat.Strand.Should().Be('+');
        feat.FeatureId.Should().Be("g1", "the ID attribute round-trips as the feature id");
        feat.Attributes.Should().ContainKey("note")
            .WhoseValue.Should().Be("embedded\ttab",
                "the tab-bearing attribute value round-trips byte-exact: escaped to %09 on write, percent-decoded on read");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-CODING-001 — coding-potential (CPAT hexamer) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-CODING-001 — coding-potential scoring

    #region Helpers — coding-potential constructions

    /// <summary>
    /// Builds a coding/non-coding hexamer-frequency table pair that is biased TOWARD coding:
    /// every in-frame hexamer of <paramref name="codingSeq"/> gets a high coding count and a
    /// low non-coding count, so any sequence built from those hexamers scores ln(high/low) &gt; 0.
    /// All frequencies are strictly positive (both-&gt;0 branch, INV-02), so the per-hexamer
    /// contribution is the natural log ln(codingValue/noncodingValue) — derived from the spec
    /// (Coding_Potential_Calculation.md §2.2), not from the code's own arrays.
    /// </summary>
    private static (Dictionary<string, double> coding, Dictionary<string, double> noncoding)
        CodingBiasedTables(string codingSeq, double codingValue, double noncodingValue)
    {
        var coding = new Dictionary<string, double>();
        var noncoding = new Dictionary<string, double>();
        for (int i = 0; i + 6 <= codingSeq.Length; i += 3)
        {
            string h = codingSeq.Substring(i, 6).ToUpperInvariant();
            coding[h] = codingValue;
            noncoding[h] = noncodingValue;
        }
        return (coding, noncoding);
    }

    #endregion

    #region BE — Boundary: empty / too-short sequence returns a defined 0 (no divide-by-zero)

    /// <summary>
    /// BE: the EMPTY sequence is the headline boundary for this unit. Per INV-06 /
    /// Coding_Potential_Calculation.md §3.3, §6.1, `sequence.Length &lt; wordSize` (empty is the
    /// extreme) short-circuits to a DEFINED score of 0 BEFORE any hexamer is scored — so there
    /// is no division by the zero scored-hexamer count, no NaN, no Infinity and no throw. The
    /// tables are non-empty and well-formed; the degeneracy is purely the input length. We pin
    /// the exact 0.0 (not merely "not NaN") because the contract promises a defined non-score.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateCodingPotential_EmptySequence_ReturnsDefinedZeroNoDivideByZero()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2 };

        double score = 0;
        var act = () => score = GenomeAnnotator.CalculateCodingPotential(string.Empty, coding, noncoding);
        act.Should().NotThrow("empty input short-circuits via the length guard before any scoring (INV-06)");

        score.Should().Be(0.0, "an empty sequence has no scorable hexamer, so the defined non-score is 0");
        double.IsNaN(score).Should().BeFalse("the length guard prevents any 0/0 division");
        double.IsInfinity(score).Should().BeFalse("no Infinity can arise from the guarded empty case");
    }

    /// <summary>
    /// BE: a sequence SHORTER than one hexamer window (length 5 &lt; wordSize 6) is the other half
    /// of the too-short boundary (INV-06). It must also resolve to a defined 0 with no division
    /// by zero — the loop `i + wordSize &lt;= length` never executes, so no hexamer is scored and
    /// the count-0 path returns 0 (Coding_Potential_Calculation.md §6.1 row 1). Pinned exactly.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateCodingPotential_SequenceShorterThanWord_ReturnsDefinedZero()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAA", coding, noncoding); // length 5 < 6

        score.Should().Be(0.0, "no full-length hexamer fits, so the scored-hexamer count is 0 → defined non-score 0");
        double.IsNaN(score).Should().BeFalse("no 0/0 division occurs for the too-short input");
    }

    #endregion

    #region BE — Boundary: random sequence is classified NON-coding (score < 0)

    /// <summary>
    /// BE (checklist "random seq"): a fixed-seed RANDOM DNA sequence carries no coding bias, so
    /// scored against tables that favour NON-coding (every in-frame hexamer present in both
    /// tables with a higher non-coding than coding value) it must land on the NON-coding side of
    /// the decision threshold — score &lt; 0 (INV-04, Coding_Potential_Calculation.md §2.4). We
    /// build coding-vs-noncoding tables from the SAME random sequence's hexamers but invert the
    /// bias (coding value 1, non-coding value 4), so every scored hexamer contributes
    /// ln(1/4) = −1.3862943611 &lt; 0 and the mean is exactly that negative value. The fuzz bar:
    /// no crash / NaN / Infinity, and the documented sign convention (non-coding ⇒ score &lt; 0)
    /// holds for random input.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void CalculateCodingPotential_RandomSequence_NoncodingBiasedTables_ClassifiedNonCoding()
    {
        string seq = RandomDna(600, seed: 216_001);
        var (coding, noncoding) = CodingBiasedTables(seq, codingValue: 1.0, noncodingValue: 4.0);

        double score = 0;
        var act = () => score = GenomeAnnotator.CalculateCodingPotential(seq, coding, noncoding);
        act.Should().NotThrow("a random in-frame sequence scored against well-formed tables never crashes");

        double.IsNaN(score).Should().BeFalse("a finite log-ratio mean can never be NaN here");
        double.IsInfinity(score).Should().BeFalse("all table values are strictly positive, so no log of 0 arises");
        score.Should().BeLessThan(0.0,
            "INV-04: non-coding-biased tables put a random sequence on the non-coding side of the threshold (0)");
        score.Should().BeApproximately(Math.Log(1.0 / 4.0), 1e-9,
            "every scored hexamer contributes ln(1/4); their mean is exactly ln(1/4) (INV-01, INV-02)");
    }

    /// <summary>
    /// BE / robustness sweep: across many fixed-seed random sequences AND randomly drawn table
    /// values, the score must stay FINITE (no NaN / Infinity), respect INV-06 for too-short
    /// inputs, and obey the sign convention dictated by which table is biased high. This is the
    /// "degenerate / random parameters do not derail the scan" boundary sweep. For each trial we
    /// pick a random bias direction: when coding &gt; noncoding the per-hexamer ln-ratio is &gt; 0 so
    /// the mean (when any hexamer scores) is &gt; 0; when noncoding &gt; coding it is &lt; 0. All values
    /// are strictly positive so no log-of-zero / Infinity can appear.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CalculateCodingPotential_RandomSweep_StaysFiniteAndRespectsSignConvention()
    {
        var rng = new Random(216_777);
        for (int trial = 0; trial < 40; trial++)
        {
            int len = rng.Next(0, 300);            // includes the too-short (0..5) boundary
            string seq = RandomDna(len, seed: 216_900 + trial);

            bool codingBias = rng.Next(2) == 0;
            double hi = 2.0 + rng.NextDouble() * 8.0;   // strictly positive
            double lo = 0.25 + rng.NextDouble() * 1.5;  // strictly positive, < hi
            double codingValue = codingBias ? hi : lo;
            double noncodingValue = codingBias ? lo : hi;

            var (coding, noncoding) = CodingBiasedTables(seq, codingValue, noncodingValue);

            double score = 0;
            var act = () => score = GenomeAnnotator.CalculateCodingPotential(seq, coding, noncoding);
            act.Should().NotThrow($"trial {trial}: well-formed tables + arbitrary length never crash");

            double.IsNaN(score).Should().BeFalse($"trial {trial}: score must be finite, never NaN");
            double.IsInfinity(score).Should().BeFalse($"trial {trial}: strictly positive table values forbid Infinity");

            // Count scorable in-frame hexamers from first principles (spec §2.2 INV-05/06).
            int scorable = 0;
            for (int i = 0; i + 6 <= seq.Length; i += 3) scorable++;

            if (scorable == 0)
            {
                score.Should().Be(0.0, $"trial {trial}: too-short / no scorable hexamer ⇒ defined 0 (INV-06)");
            }
            else if (codingBias)
            {
                score.Should().BeGreaterThan(0.0, $"trial {trial}: coding-biased tables ⇒ score > 0 (INV-04)");
            }
            else
            {
                score.Should().BeLessThan(0.0, $"trial {trial}: non-coding-biased tables ⇒ score < 0 (INV-04)");
            }
        }
    }

    #endregion

    #region BE — Boundary: perfect ORF is classified CODING (score > 0)

    /// <summary>
    /// BE (checklist "perfect ORF"): a clean, in-frame coding ORF (ATG + Ala codons + TAA),
    /// scored against tables trained to favour exactly its hexamers (high coding value, low
    /// non-coding value), must land on the CODING side of the decision threshold — score &gt; 0
    /// (INV-04). Every in-frame hexamer of the ORF is in both tables with coding 8 vs non-coding
    /// 2, so each contributes ln(8/2) = ln 4 = 1.3862943611 and the mean is exactly that
    /// positive value. We pin both the sign (classified coding) and the exact mean derived from
    /// the spec (Coding_Potential_Calculation.md §2.2), proving the perfect-ORF boundary is the
    /// mirror of the random/non-coding boundary above.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateCodingPotential_PerfectOrf_CodingBiasedTables_ClassifiedCoding()
    {
        // Perfect in-frame ORF: ATG + (GCT)x12 [Ala] + TAA — frame 0 throughout.
        string orf = "ATG" + string.Concat(Enumerable.Repeat("GCT", 12)) + "TAA";
        var (coding, noncoding) = CodingBiasedTables(orf, codingValue: 8.0, noncodingValue: 2.0);

        double score = GenomeAnnotator.CalculateCodingPotential(orf, coding, noncoding);

        double.IsNaN(score).Should().BeFalse("a coding ORF scored against well-formed tables yields a finite score");
        score.Should().BeGreaterThan(0.0,
            "INV-04: coding-biased tables place a perfect ORF on the coding side of the threshold (0)");
        score.Should().BeApproximately(Math.Log(8.0 / 2.0), 1e-9,
            "every in-frame hexamer contributes ln(8/2) = ln 4; their mean is exactly ln 4 (INV-01, INV-02)");
    }

    /// <summary>
    /// BE: the +1 (coding-only) and −1 (non-coding-only) sentinel branches (INV-03) are the
    /// degenerate-table boundary. A hexamer present ONLY in the coding table contributes +1; one
    /// present ONLY in the non-coding table contributes −1. We build a 2-hexamer sequence where
    /// the first hexamer is coding-only and the second non-coding-only, so the mean is
    /// (+1 + −1)/2 = 0 exactly — pinning that the sentinel branches are reached and combined per
    /// INV-01 without crashing or producing NaN/Infinity (no log of zero is taken on these paths).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateCodingPotential_TablePresenceSentinels_ProducePlusOneAndMinusOne()
    {
        const string seq = "ATGAAACCC"; // in-frame hexamers ATGAAA (i=0), AAACCC (i=3)
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 5.0 };   // coding-only ⇒ +1
        var noncoding = new Dictionary<string, double> { ["AAACCC"] = 5.0 }; // non-coding-only ⇒ −1

        double score = GenomeAnnotator.CalculateCodingPotential(seq, coding, noncoding);

        double.IsNaN(score).Should().BeFalse("the ±1 sentinel branches take no logarithm, so no NaN arises");
        score.Should().BeApproximately(0.0, 1e-12,
            "INV-03: coding-only ⇒ +1 and non-coding-only ⇒ −1; INV-01 mean = (+1 + −1)/2 = 0");
    }

    #endregion

    #region Positive sanity — exact worked example reproduced from first principles

    /// <summary>
    /// Positive sanity: the documented worked example (Coding_Potential_Calculation.md §7.1)
    /// scored from FIRST PRINCIPLES, independently of the implementation. Sequence "ATGAAACCC"
    /// (length 9), wordSize 6, stepSize 3 → in-frame hexamers ATGAAA (i=0) and AAACCC (i=3);
    /// i=6 yields "CCC" (length 3 &lt; 6) and is skipped (INV-05). With coding {ATGAAA:8, AAACCC:2}
    /// and non-coding {ATGAAA:2, AAACCC:4}:
    ///   ATGAAA: ln(8/2) = ln 4    =  1.3862943611198906
    ///   AAACCC: ln(2/4) = ln 0.5  = −0.6931471805599453
    ///   sum = 0.6931471805599453, count = 2 ⇒ 0.34657359027997264.
    /// The expected value is computed here from Math.Log so the assertion encodes the SPEC's
    /// arithmetic, not a number copied from the code. This guards the whole unit: if the scan,
    /// the log base, the denominator or the skip rule drifts, this exact value breaks.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateCodingPotential_WorkedExample_MatchesExactSpecScore()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8, ["AAACCC"] = 2 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2, ["AAACCC"] = 4 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        double expected = (Math.Log(8.0 / 2.0) + Math.Log(2.0 / 4.0)) / 2.0; // = 0.34657359027997264
        score.Should().BeApproximately(expected, 1e-12,
            "the score is the mean of ln(8/2) and ln(2/4) over the two in-frame hexamers (INV-01, INV-02; §7.1)");
        score.Should().BeApproximately(0.34657359027997264, 1e-12,
            "the documented worked-example score is reproduced exactly");
    }

    #endregion

    #region Validation — documented throws are pinned, not weakened

    /// <summary>
    /// Validation: null sequence or either table throws ArgumentNullException, and
    /// wordSize ≤ 0 / stepSize ≤ 0 throws ArgumentOutOfRangeException
    /// (Coding_Potential_Calculation.md §3.3, §6.1). These are the DOCUMENTED guards; the fuzz
    /// test pins the exact throw type rather than weakening it to a tolerant no-op.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void CalculateCodingPotential_InvalidArguments_ThrowDocumentedExceptions()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2 };

        ((Action)(() => GenomeAnnotator.CalculateCodingPotential(null!, coding, noncoding)))
            .Should().Throw<ArgumentNullException>("a null sequence is a documented ArgumentNullException");
        ((Action)(() => GenomeAnnotator.CalculateCodingPotential("ATGAAA", null!, noncoding)))
            .Should().Throw<ArgumentNullException>("a null coding table is a documented ArgumentNullException");
        ((Action)(() => GenomeAnnotator.CalculateCodingPotential("ATGAAA", coding, null!)))
            .Should().Throw<ArgumentNullException>("a null non-coding table is a documented ArgumentNullException");

        ((Action)(() => GenomeAnnotator.CalculateCodingPotential("ATGAAA", coding, noncoding, wordSize: 0)))
            .Should().Throw<ArgumentOutOfRangeException>("wordSize ≤ 0 is a documented ArgumentOutOfRangeException");
        ((Action)(() => GenomeAnnotator.CalculateCodingPotential("ATGAAA", coding, noncoding, stepSize: 0)))
            .Should().Throw<ArgumentOutOfRangeException>("stepSize ≤ 0 is a documented ArgumentOutOfRangeException");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-CODONUSAGE-001 — Relative Synonymous Codon Usage : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-CODONUSAGE-001 — Relative Synonymous Codon Usage (RSCU)

    #region Helpers — RSCU constants (NCBI translation table 1)

    /// <summary>
    /// The six synonymous Leucine codons of the Standard genetic code (NCBI table 1),
    /// n_i = 6. Derived from the genetic code, NOT from the implementation's arrays.
    /// </summary>
    private static readonly string[] LeucineFamily = { "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" };

    /// <summary>The Standard code's three stop codons (TAA, TAG, TGA) — excluded from RSCU output.</summary>
    private static readonly string[] StopCodons = { "TAA", "TAG", "TGA" };

    #endregion

    #region Positive sanity — the doc §7.1 Leucine worked example reproduced from first principles

    /// <summary>
    /// Positive sanity: the doc's hand-checkable worked example
    /// (Relative_Synonymous_Codon_Usage.md §7.1). CDS "CTTCTTCTGTTA" → in-frame codons
    /// CTT, CTT, CTG, TTA (all Leucine). Family counts CTT=2, CTG=1, TTA=1; Σ=4; n_i=6.
    /// RSCU = n_i·x/Σ ⇒ CTT = 6·2/4 = 3.0, CTG = 6·1/4 = 1.5, TTA = 6·1/4 = 1.5; the three
    /// UNobserved Leu codons (TTG, CTC, CTA) = 0.0. These expected values are derived
    /// independently from the §2.2 formula and the NCBI Leu family — NOT echoed off the
    /// code's arrays. INV-01 (Σ over the family = n_i = 6) is pinned alongside, so a wrong
    /// implementation that mis-normalised could not pass.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_RscuLeucineWorkedExample_ReproducesExactDocValues()
    {
        IReadOnlyDictionary<string, double> rscu = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTTCTGTTA" });

        rscu["CTT"].Should().BeApproximately(3.0, 1e-10, "§7.1: CTT = 6·2/4 = 3.0 (preferred Leu codon)");
        rscu["CTG"].Should().BeApproximately(1.5, 1e-10, "§7.1: CTG = 6·1/4 = 1.5");
        rscu["TTA"].Should().BeApproximately(1.5, 1e-10, "§7.1: TTA = 6·1/4 = 1.5");
        rscu["TTG"].Should().BeApproximately(0.0, 1e-10, "§7.1: TTG unobserved ⇒ 6·0/4 = 0.0");
        rscu["CTC"].Should().BeApproximately(0.0, 1e-10, "§7.1: CTC unobserved ⇒ 0.0");
        rscu["CTA"].Should().BeApproximately(0.0, 1e-10, "§7.1: CTA unobserved ⇒ 0.0");

        double leuSum = LeucineFamily.Sum(c => rscu[c]);
        leuSum.Should().BeApproximately(6.0, 1e-10,
            "INV-01: Σ RSCU over the observed 6-codon Leu family equals n_i = 6");
    }

    #endregion

    #region BE — Boundary: empty input (defined all-zero result, NO divide-by-zero)

    /// <summary>
    /// BE: the EMPTY boundary called out in the checklist row. An empty collection AND a
    /// collection of all-empty strings observe NOTHING, so every synonymous family totals
    /// Σ = 0. The doc's §3.3 / §5.4 deviation 1 contract is that the family-total-0 branch
    /// returns RSCU 0.0 for each member — a DEFINED degenerate result, with NO division by
    /// zero (n_i·0/0 is never evaluated). The output still enumerates all 61 sense codons of
    /// the Standard code (64 − 3 stops, INV-04), every value finite and exactly 0.0. We pin
    /// no throw, the 61-codon cardinality, and that every value is 0.0 (not NaN/Infinity).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_EmptyInput_ReturnsDefinedAllZeroMap_NoDivideByZero()
    {
        var fromEmptyCollection = (() => GenomeAnnotator.GetCodonUsage(Array.Empty<string>()));
        var fromAllEmptyStrings = (() => GenomeAnnotator.GetCodonUsage(new[] { "", "", "" }));
        fromEmptyCollection.Should().NotThrow("nothing observed ⇒ every family total is 0; the 0/0 branch returns 0.0, never divides");
        fromAllEmptyStrings.Should().NotThrow("all-empty strings are skipped; the result is the same defined all-zero map");

        IReadOnlyDictionary<string, double> emptyCollection = GenomeAnnotator.GetCodonUsage(Array.Empty<string>());
        IReadOnlyDictionary<string, double> allEmptyStrings = GenomeAnnotator.GetCodonUsage(new[] { "", "", "" });

        emptyCollection.Count.Should().Be(61, "the Standard code has 64 codons − 3 stops = 61 sense codons (INV-04)");
        emptyCollection.Values.Should().OnlyContain(v => v == 0.0,
            "with nothing observed every synonymous family totals 0 ⇒ RSCU 0.0 for every codon, no NaN/Infinity");

        allEmptyStrings.Count.Should().Be(61, "all-empty strings observe nothing yet still enumerate the 61 sense codons");
        allEmptyStrings.Values.Should().OnlyContain(v => v == 0.0,
            "all-empty input yields the same defined all-zero RSCU map");

        StopCodons.Should().OnlyContain(stop => !emptyCollection.ContainsKey(stop),
            "INV-04: stop codons (TAA/TAG/TGA) never appear in the RSCU output");
    }

    #endregion

    #region BE — Boundary: a single codon

    /// <summary>
    /// BE: a SINGLE codon. A lone codon is the only member of its family observed, so its
    /// family total Σ equals its own count and RSCU = n_i·x/Σ = n_i (the maximum allowed by
    /// INV-03, [0, n_i]). For a single Leucine codon "CTT" (n_i = 6, x = 1, Σ = 1) the doc
    /// formula gives RSCU = 6·1/1 = 6.0; the other five Leu codons are unobserved ⇒ 0.0, and
    /// Σ over the family = 6.0 = n_i (INV-01). The expected 6.0 is derived from the §2.2
    /// formula and the NCBI Leu family size, NOT echoed off the code.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_SingleLeucineCodon_RscuEqualsFamilySize()
    {
        IReadOnlyDictionary<string, double> rscu = GenomeAnnotator.GetCodonUsage(new[] { "CTT" });

        rscu["CTT"].Should().BeApproximately(6.0, 1e-10,
            "a single observed Leu codon: RSCU = n_i·x/Σ = 6·1/1 = 6.0 (the INV-03 upper bound n_i)");
        rscu["TTA"].Should().BeApproximately(0.0, 1e-10, "the other Leu codons are unobserved ⇒ 0.0");
        rscu["TTG"].Should().BeApproximately(0.0, 1e-10, "the other Leu codons are unobserved ⇒ 0.0");

        double leuSum = LeucineFamily.Sum(c => rscu[c]);
        leuSum.Should().BeApproximately(6.0, 1e-10, "INV-01: Σ RSCU over the Leu family equals n_i = 6");

        rscu.Values.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v) && v >= 0.0 && v <= 6.0,
            "INV-03: every RSCU lies in [0, n_i] (≤ 6 for the Standard code); all values finite");
    }

    /// <summary>
    /// BE: a single codon of a SINGLE-codon amino acid. Methionine (ATG) and Tryptophan
    /// (TGG) are the only n_i = 1 families in the Standard code, so by INV-02 their RSCU is
    /// ALWAYS 1.0 (n_i·x/Σ = 1·x/x = 1). A lone "ATG" must therefore yield RSCU 1.0 — a
    /// degenerate single-codon input that is fully defined, not a divide-by-zero.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_SingleMethionineCodon_RscuIsOne()
    {
        IReadOnlyDictionary<string, double> rscu = GenomeAnnotator.GetCodonUsage(new[] { "ATG" });

        rscu["ATG"].Should().BeApproximately(1.0, 1e-10,
            "INV-02: Met (n_i = 1) always has RSCU = 1·1/1 = 1.0");
    }

    #endregion

    #region BE — Boundary: length not a multiple of 3 (partial trailing codon ignored)

    /// <summary>
    /// BE: a sequence whose LENGTH is NOT a multiple of 3. The doc §3.3/§4.1 contract reads
    /// in-frame triplets from index 0 in steps of 3 and IGNORES a partial trailing codon
    /// (the loop bound is `i <= len - 3`). "CTTCTTCTGTTAC" is the §7.1 worked example
    /// (length 12) plus one trailing 'C' (length 13, 13 % 3 = 1). The trailing partial codon
    /// must be silently dropped, leaving exactly the §7.1 Leu counts ⇒ identical RSCU
    /// (CTT=3.0, CTG=1.5, TTA=1.5) and no crash / no out-of-range Substring on the 3-char
    /// window. This pins that a non-%3 length neither corrupts the count nor throws.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_LengthNotMultipleOfThree_IgnoresPartialTrailingCodon()
    {
        const string withTrailingPartial = "CTTCTTCTGTTAC"; // §7.1 example (12 nt) + trailing 'C' = 13 nt

        var act = (() => GenomeAnnotator.GetCodonUsage(new[] { withTrailingPartial }));
        act.Should().NotThrow("the partial trailing codon is dropped by the `i <= len - 3` framing; no out-of-range Substring");

        IReadOnlyDictionary<string, double> rscu = GenomeAnnotator.GetCodonUsage(new[] { withTrailingPartial });

        rscu["CTT"].Should().BeApproximately(3.0, 1e-10, "trailing 'C' dropped ⇒ same §7.1 Leu counts ⇒ CTT = 3.0");
        rscu["CTG"].Should().BeApproximately(1.5, 1e-10, "trailing 'C' dropped ⇒ CTG = 1.5");
        rscu["TTA"].Should().BeApproximately(1.5, 1e-10, "trailing 'C' dropped ⇒ TTA = 1.5");

        LeucineFamily.Sum(c => rscu[c]).Should().BeApproximately(6.0, 1e-10,
            "INV-01 still holds: the ignored partial codon does not perturb the family sum");
    }

    /// <summary>
    /// BE: the most degenerate non-%3 lengths — a 1-nt and a 2-nt sequence have NO complete
    /// in-frame codon at all (`0 &gt; len - 3`), so nothing is observed and the result is the
    /// defined all-zero 61-codon map, with no crash and no divide-by-zero. This is the
    /// length-not-%3 boundary at its smallest: a sub-codon input.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_SubCodonLength_ObservesNothing_DefinedAllZeroMap()
    {
        foreach (string subCodon in new[] { "A", "CT" })
        {
            var act = (() => GenomeAnnotator.GetCodonUsage(new[] { subCodon }));
            act.Should().NotThrow($"\"{subCodon}\" has no complete codon; the loop body never runs, no crash");

            IReadOnlyDictionary<string, double> rscu = GenomeAnnotator.GetCodonUsage(new[] { subCodon });
            rscu.Count.Should().Be(61, "the output still enumerates all 61 sense codons even with nothing observed");
            rscu.Values.Should().OnlyContain(v => v == 0.0,
                $"\"{subCodon}\" observes no codon ⇒ every family totals 0 ⇒ RSCU 0.0, no divide-by-zero");
        }
    }

    #endregion

    #region BE — Randomized boundary sweep (no crash / NaN / Infinity / out-of-range)

    /// <summary>
    /// BE / randomized sweep: across many fixed-seed random sequences — including ones whose
    /// length is deliberately NOT a multiple of 3 and short single-/sub-codon inputs — the RSCU
    /// computation must never crash, hang, or emit a NaN / Infinity / out-of-range value. For
    /// EVERY input the output enumerates the 61 sense codons (INV-04, no stops), every value is
    /// finite and within [0, 6] (INV-03, the Standard code's largest family n_i = 6), and every
    /// OBSERVED synonymous family sums to its n_i within tolerance (INV-01). The per-family sum
    /// check is computed independently from the NCBI Leu family here as a falsifiable witness
    /// that the normalisation is per-family (per the §2.2 formula), not per-total. Deterministic:
    /// the input lengths and bases come from a locally fixed `new Random(seed)`.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GetCodonUsage_RandomBoundarySweep_NeverCrashesAndStaysWithinInvariants()
    {
        var rng = new Random(217_001);

        for (int iter = 0; iter < 200; iter++)
        {
            // Lengths span 0..~150 including many non-%3 and sub-codon cases.
            int length = rng.Next(0, 151);
            string seq = RandomDna(length, seed: rng.Next());

            IReadOnlyDictionary<string, double>? rscu = null;
            var act = (() => rscu = GenomeAnnotator.GetCodonUsage(new[] { seq }));
            act.Should().NotThrow($"RSCU must not crash on random input (iter {iter}, length {length})");

            rscu!.Count.Should().Be(61, "every result enumerates exactly the 61 sense codons of the Standard code");
            rscu.Values.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v) && v >= 0.0 && v <= 6.0,
                $"INV-03: every RSCU is finite and in [0, n_i ≤ 6] (iter {iter})");
            StopCodons.Should().OnlyContain(stop => !rscu.ContainsKey(stop),
                "INV-04: no stop codon ever appears in the RSCU output");

            // INV-01 witness on the Leu family: if any Leu codon is observed (sum > 0), the
            // whole family must sum to n_i = 6; an entirely unobserved family sums to 0.0.
            double leuSum = LeucineFamily.Sum(c => rscu[c]);
            bool leuObserved = leuSum > 0.0;
            if (leuObserved)
                leuSum.Should().BeApproximately(6.0, 1e-9,
                    $"INV-01: an observed Leu family sums to n_i = 6 (per-family normalisation) (iter {iter})");
            else
                leuSum.Should().Be(0.0, $"an unobserved Leu family sums to exactly 0.0 (iter {iter})");
        }
    }

    #endregion

    #region Validation — documented null throws (pinned, not weakened)

    /// <summary>
    /// Validation: a null `codingSequences` collection (both overloads) and a null
    /// `GeneticCode` throw ArgumentNullException (Relative_Synonymous_Codon_Usage.md §3.3,
    /// §6.1). These are DOCUMENTED guards; the fuzz test pins the exact throw type rather
    /// than weakening it to a tolerant no-op.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GetCodonUsage_NullArguments_ThrowDocumentedArgumentNullException()
    {
        ((Action)(() => GenomeAnnotator.GetCodonUsage((IEnumerable<string>)null!)))
            .Should().Throw<ArgumentNullException>("a null coding-sequence collection is a documented ArgumentNullException");
        ((Action)(() => GenomeAnnotator.GetCodonUsage((IEnumerable<string>)null!, GeneticCode.Standard)))
            .Should().Throw<ArgumentNullException>("a null collection (genetic-code overload) is a documented ArgumentNullException");
        ((Action)(() => GenomeAnnotator.GetCodonUsage(new[] { "ATG" }, null!)))
            .Should().Throw<ArgumentNullException>("a null GeneticCode is a documented ArgumentNullException");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  ANNOT-REPEAT-001 — repetitive-element detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region ANNOT-REPEAT-001 — repetitive-element detection

    #region Positive sanity — doc §7.1 worked example (hand-checkable)

    /// <summary>
    /// Positive sanity: the Repetitive_Element_Detection.md §7.1 worked example,
    /// computed INDEPENDENTLY from the model (NOT echoed off the code). The motif
    /// "ATTCG" tiled three times head-to-tail gives "ATTCGATTCGATTCG" (15 bp). With
    /// minRepeatLength = 5 and minCopies = 2 the canonical full-span primitive array is
    /// (start = 0, end = 15, type = "tandem_repeat", sequence = "ATTCGATTCGATTCG") — the
    /// unit "ATTCG" (|u| = 5) repeated 3× over [0, 15). INV-01: sequence == dnaSequence
    /// [0..15] and 15 is a multiple of 5; INV-02: 3 ≥ minCopies adjacent copies. NOTE the
    /// detector ALSO legitimately reports the SHIFTED sub-arrays (e.g. "TTCGA" × 2 at
    /// start 1) — the §7.1 snippet shows the canonical element, NOT a uniqueness claim,
    /// matching the committed unit test which pins the (0, 15) element by membership. So
    /// we assert the canonical element is PRESENT with exact fields, and that EVERY
    /// reported tandem array is itself structurally valid (slice == source, span a
    /// multiple of its primitive period, ≥ minCopies copies). We isolate the tandem class
    /// because minRepeatLength doubles as the inverted-arm length, irrelevant here.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindRepetitiveElements_WorkedExample_ReportsCanonicalFullSpanTandemUnitAttcgTimesThree()
    {
        const string seq = "ATTCGATTCGATTCG"; // unit "ATTCG" × 3, 15 bp

        var tandem = GenomeAnnotator
            .FindRepetitiveElements(seq, minRepeatLength: 5, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        var tr = tandem.Should().ContainSingle(e => e.start == 0 && e.end == 15,
            "the §7.1 canonical element is the full-span \"ATTCG\" × 3 array over [0, 15)").Subject;
        tr.type.Should().Be("tandem_repeat");
        tr.sequence.Should().Be("ATTCGATTCGATTCG",
            "INV-01: the reported sequence is dnaSequence[0..15], the array slice");
        ((tr.end - tr.start) % 5).Should().Be(0, "INV-01: the span is an integer multiple of the unit length |u| = 5");
        ((tr.end - tr.start) / 5).Should().BeGreaterThanOrEqualTo(2, "INV-02: at least minCopies adjacent copies");

        // Every reported tandem array (canonical + shifted sub-arrays) is structurally valid.
        tandem.Should().OnlyContain(e =>
            e.start >= 0 && e.start < e.end && e.end <= seq.Length &&
            e.sequence == seq.Substring(e.start, e.end - e.start),
            "every reported tandem array stays in bounds and its slice equals the source (INV-01)");
    }

    #endregion

    #region BE — Boundary: no repeat present (empty result)

    /// <summary>
    /// BE ("no repeat"): a sequence with NO tandem array of ≥ 2 adjacent copies and NO
    /// inverted (revcomp-palindrome) arm must yield an EMPTY result — no crash, no spurious
    /// element (Repetitive_Element_Detection.md §6.1, INV-01..INV-03). "ACGT" has every base
    /// distinct (no head-to-tail repeat) and is its own revcomp only as a whole, but with
    /// minRepeatLength = 4 the inverted scan needs i + 2·4 ≤ 4 (impossible) and the tandem
    /// scan needs span ≥ 4 with ≥ 2 copies (impossible for 4 distinct bases). The disciplined
    /// outcome is the empty list.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindRepetitiveElements_NoRepeatPresent_ReturnsEmpty()
    {
        const string seq = "ACGT"; // all distinct, no adjacent repeat, no in-range inverted arm

        var act = () => GenomeAnnotator.FindRepetitiveElements(seq, minRepeatLength: 4, minCopies: 2).ToList();
        act.Should().NotThrow("a non-repetitive sequence is scanned cleanly to the end");

        GenomeAnnotator.FindRepetitiveElements(seq, minRepeatLength: 4, minCopies: 2).ToList()
            .Should().BeEmpty("no head-to-tail repeat and no in-range inverted arm exist in \"ACGT\"");
    }

    #endregion

    #region BE — Boundary: full repeat (whole sequence is one tandem array)

    /// <summary>
    /// BE ("full repeat"): the WHOLE sequence is a single tandem array. "GATCGATCGATC" is the
    /// primitive unit "GATC" tiled three times (12 bp). With minRepeatLength = 4, minCopies = 2
    /// the array spans 3·4 = 12 = the whole sequence, reported once with start = 0, end = length
    /// (INV-01, INV-02). The primitive-collapse rule means a homopolymeric/periodic whole sequence
    /// is reported by its shortest period, never double-counted. We pin that the full-span array
    /// is reported with the correct primitive unit and that NO tandem array extends past the
    /// sequence end.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindRepetitiveElements_FullRepeat_ReportsWholeSequenceAsSingleTandemArray()
    {
        const string seq = "GATCGATCGATC"; // unit "GATC" × 3, 12 bp — the entire sequence

        var tandem = GenomeAnnotator
            .FindRepetitiveElements(seq, minRepeatLength: 4, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        var full = tandem.Should().ContainSingle(o => o.start == 0 && o.end == seq.Length,
            "the whole sequence is one primitive tandem array spanning [0, length)").Subject;

        full.sequence.Should().Be(seq, "INV-01: the array slice is the entire dnaSequence[0..length]");
        // The primitive unit is the shortest period; for GATC×3 that is "GATC" (4 bp).
        (full.end - full.start).Should().Be(12, "the span is the full 12 bp");
        (12 % 4).Should().Be(0, "INV-01: 12 is a multiple of the primitive unit length 4");

        tandem.Should().OnlyContain(o => o.start >= 0 && o.start < o.end && o.end <= seq.Length,
            "every reported tandem array stays within the sequence bounds");
    }

    #endregion

    #region BE — Boundary: minLength edge (span exactly threshold vs one below)

    /// <summary>
    /// BE ("minLen edge"): minRepeatLength is the minimum TOTAL SPAN (bp) of a reported element
    /// (Repetitive_Element_Detection.md §3.1; the array is dropped when end − start &lt;
    /// minRepeatLength, GenomeAnnotator.cs line 840). The dinucleotide array "ATATAT" is the
    /// primitive unit "AT" × 3, spanning EXACTLY 6 bp. We pin the edge from both sides:
    ///   • minRepeatLength = 6 (span == threshold) → the array IS reported (6 ≥ 6);
    ///   • minRepeatLength = 7 (span one below threshold) → the SAME array is DROPPED (6 &lt; 7).
    /// This is the falsifiable boundary: an off-by-one in the span filter would either drop the
    /// exactly-threshold array or keep the one-below array. We isolate the tandem class because
    /// minRepeatLength doubles as the inverted-arm length, which here (needs i + 2·6 ≤ 6) admits
    /// no inverted repeat anyway, so the assertion targets the span filter cleanly.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindRepetitiveElements_MinLengthEdge_ReportsAtThresholdAndDropsOneBelow()
    {
        const string seq = "ATATAT"; // unit "AT" × 3, span exactly 6 bp

        var atThreshold = GenomeAnnotator
            .FindRepetitiveElements(seq, minRepeatLength: 6, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        var tr = atThreshold.Should().ContainSingle(
            "a tandem array whose span (6) is EXACTLY minRepeatLength (6) is reported").Subject;
        tr.start.Should().Be(0);
        tr.end.Should().Be(6, "the array spans the whole 6 bp");
        tr.sequence.Should().Be("ATATAT", "INV-01: the reported slice is dnaSequence[0..6]");

        var oneBelow = GenomeAnnotator
            .FindRepetitiveElements(seq, minRepeatLength: 7, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        oneBelow.Should().BeEmpty(
            "the same 6 bp array is DROPPED when its span (6) is one base below the threshold (7)");
    }

    #endregion

    #region BE — Boundary: empty / null sequence and out-of-range guards

    /// <summary>
    /// BE (empty / null + guards): the documented validation contract
    /// (Repetitive_Element_Detection.md §3.3, §6.1). A null sequence throws
    /// ArgumentNullException; an EMPTY sequence yields an EMPTY result (NOT a throw — the core
    /// `yield break`s on length 0); minCopies &lt; 2 and minRepeatLength &lt; 1 each throw
    /// ArgumentOutOfRangeException. The throws are eager (argument validation runs before the
    /// iterator body), but we force enumeration on the empty case to surface the lazy short-circuit.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindRepetitiveElements_EmptyAndInvalidArguments_HonourDocumentedContract()
    {
        ((Action)(() => GenomeAnnotator.FindRepetitiveElements(null!).ToList()))
            .Should().Throw<ArgumentNullException>("a null dnaSequence is a documented ArgumentNullException");

        GenomeAnnotator.FindRepetitiveElements(string.Empty).ToList()
            .Should().BeEmpty("an empty sequence has no repetitive elements; the core yield-breaks, it does NOT throw");

        ((Action)(() => GenomeAnnotator.FindRepetitiveElements("ACGTACGT", minRepeatLength: 4, minCopies: 1).ToList()))
            .Should().Throw<ArgumentOutOfRangeException>("minCopies < 2 is invalid: a tandem repeat needs ≥ 2 adjacent copies");

        ((Action)(() => GenomeAnnotator.FindRepetitiveElements("ACGTACGT", minRepeatLength: 0, minCopies: 2).ToList()))
            .Should().Throw<ArgumentOutOfRangeException>("minRepeatLength < 1 is invalid");
    }

    #endregion

    #region MC / RB — Malformed content + randomized boundary sweep

    /// <summary>
    /// MC: non-DNA characters (N, IUPAC ambiguity codes, digits, punctuation) embedded in the
    /// scanned sequence are NOT pre-validated and must NOT crash (Repetitive_Element_Detection.md
    /// §3.3). The detector upper-cases the input then matches exactly, so a garbage triplet simply
    /// never participates in a head-to-tail unit match or a revcomp arm match — no out-of-range
    /// Substring, no spurious element. We pin crash-freedom and that the result is structurally
    /// valid (every element in bounds and of a known type).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindRepetitiveElements_NonDnaContent_DoesNotCrashAndStaysWellFormed()
    {
        const string garbage = "NNNRYK123!?WSMxyz"; // no ACGT repeat, no revcomp arm

        List<(int start, int end, string type, string sequence)> elements = null!;
        var act = () => elements = GenomeAnnotator.FindRepetitiveElements(garbage, minRepeatLength: 2, minCopies: 2).ToList();
        act.Should().NotThrow("non-DNA content fails exact matching but never indexes past a window or throws");

        elements.Should().OnlyContain(e =>
            (e.type == "tandem_repeat" || e.type == "inverted_repeat") &&
            e.start >= 0 && e.start < e.end && e.end <= garbage.Length,
            "any reported element on malformed content is still a known type and stays in bounds");
    }

    /// <summary>
    /// RB / boundary sweep: a fixed-seed random sweep over varied lengths and parameters must
    /// complete promptly and keep EVERY reported element structurally valid — no crash, no hang,
    /// no corruption (Repetitive_Element_Detection.md §2.4 INV-01..INV-03). For each element:
    /// it is in bounds with start &lt; end; a tandem array's slice equals dnaSequence[start..end]
    /// and its span is an integer multiple of its primitive unit length with ≥ minCopies copies
    /// (INV-01, INV-02); an inverted repeat's first arm reverse-complements to its last arm
    /// (INV-03), checked on the literal arm slices (the composite `sequence` embeds the gap in
    /// brackets, so we read the arms back off the source sequence). Deterministic only:
    /// `new Random(seed)` is fixed locally.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindRepetitiveElements_RandomSweep_AllElementsStructurallyValid()
    {
        var rng = new Random(218_001);

        for (int iter = 0; iter < 120; iter++)
        {
            int length = rng.Next(0, 121);
            string seq = RandomDna(length, seed: rng.Next());
            int minRepeatLength = rng.Next(1, 8);
            int minCopies = rng.Next(2, 5);

            List<(int start, int end, string type, string sequence)> elements = null!;
            var act = () => elements = GenomeAnnotator
                .FindRepetitiveElements(seq, minRepeatLength, minCopies).ToList();
            act.Should().NotThrow(
                $"random input must never crash (iter {iter}, len {length}, minLen {minRepeatLength}, minCopies {minCopies})");

            foreach (var e in elements)
            {
                // In-bounds, exclusive-end, span ≥ minRepeatLength.
                (e.start >= 0 && e.start < e.end && e.end <= seq.Length).Should().BeTrue(
                    $"every element is in bounds (iter {iter})");
                (e.end - e.start).Should().BeGreaterThanOrEqualTo(minRepeatLength,
                    $"every reported element spans ≥ minRepeatLength (iter {iter})");

                if (e.type == "tandem_repeat")
                {
                    string slice = seq.ToUpperInvariant().Substring(e.start, e.end - e.start);
                    e.sequence.Should().Be(slice, $"INV-01: tandem slice == dnaSequence[start..end] (iter {iter})");

                    // INV-01/INV-02: the slice is a head-to-tail repeat of its primitive unit
                    // with ≥ minCopies copies. Derive the unit independently as the shortest
                    // period p of the slice such that slice == unit repeated (len % p == 0).
                    int span = slice.Length;
                    int period = SmallestTandemPeriod(slice);
                    period.Should().BePositive($"a tandem array has a head-to-tail period (iter {iter})");
                    (span % period).Should().Be(0, $"INV-01: span is a multiple of the unit length (iter {iter})");
                    (span / period).Should().BeGreaterThanOrEqualTo(minCopies,
                        $"INV-02: at least minCopies adjacent copies (iter {iter})");
                }
                else
                {
                    e.type.Should().Be("inverted_repeat", $"the only other class is inverted_repeat (iter {iter})");
                    // INV-03: the left arm's reverse complement equals the right arm. Arms are
                    // equal-length; the span minus the gap is 2·armLen, so each arm length is
                    // recovered by stripping the bracketed gap. Read arms off the source slice.
                    string upper = seq.ToUpperInvariant();
                    string composite = e.sequence;
                    int gapStart = composite.IndexOf('[');
                    int armLen;
                    if (gapStart < 0)
                    {
                        armLen = composite.Length / 2; // palindrome (gap 0): arm1arm2
                    }
                    else
                    {
                        armLen = gapStart; // arm1 precedes the '[' gap marker
                    }
                    armLen.Should().BePositive($"an inverted repeat has positive-length arms (iter {iter})");
                    string leftArm = upper.Substring(e.start, armLen);
                    string rightArm = upper.Substring(e.end - armLen, armLen);
                    DnaSequence.GetReverseComplementString(leftArm).Should().Be(rightArm,
                        $"INV-03: the right arm equals revcomp(left arm) (iter {iter})");
                }
            }
        }
    }

    /// <summary>
    /// Independently derives the shortest head-to-tail period p of <paramref name="s"/> such
    /// that s is exactly (s[0..p]) repeated len/p times. Used by the random sweep to verify
    /// INV-01/INV-02 from the slice alone, NOT from any value the code returned.
    /// </summary>
    private static int SmallestTandemPeriod(string s)
    {
        int n = s.Length;
        for (int p = 1; p <= n; p++)
        {
            if (n % p != 0) continue;
            bool ok = true;
            for (int i = p; i < n && ok; i++)
                if (s[i] != s[i - p]) ok = false;
            if (ok) return p;
        }
        return n;
    }

    #endregion

    #endregion
}
