using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the FileIO area — FASTA parsing (PARSE-FASTA-001) and FASTQ
/// parsing (PARSE-FASTQ-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, truncated, random and injected input to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no out-of-bounds indexing leaking from an internal Substring, no
/// NullReferenceException, no OutOfMemoryException or quadratic blow-up on a huge
/// header, and no silent emission of a corrupt record. Every input must resolve to
/// EITHER a well-defined, theory-correct parse, OR a DOCUMENTED, intentional
/// validation exception. A raw, unhandled runtime exception, a hang, or a silently
/// corrupt record on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing": parsers of file formats are
///   THE prime fuzz target because they sit on the untrusted-input boundary.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PARSE-FASTA-001 — FASTA parsing (FileIO)
/// Checklist: docs/checklists/03_FUZZING.md, row 64. This is the first FileIO-area
/// fuzz unit and the single highest-value fuzz target in the checklist.
/// Fuzz strategies exercised for THIS unit:
///   • RB  = Random Bytes — binary garbage / random characters fed as FASTA text.
///   • TF  = Truncated Fields — a record with no trailing newline; a header with no
///           sequence; a defline-less sequence block (the leading record header is
///           "truncated" away).
///   • MC  = Malformed Content — missing '>' marker, empty sequence body.
///   • INJ = Injection — embedded NUL bytes and other control / unicode characters
///           inside the defline and the sequence body.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The FASTA-parsing contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// FASTA is a line-oriented format: each record is a defline beginning with '>'
/// (its first whitespace-delimited token is the identifier, the remainder is the
/// optional description) followed by one or more sequence lines.[Lipman &amp; Pearson
/// 1985; NCBI BLAST FASTA guidance]
///   — docs/algorithms/FileIO/FASTA_Parsing.md §1, §2.1, §2.2.
///
/// API entry points (FastaParser.cs):
///   • IEnumerable&lt;FastaEntry&gt; FastaParser.Parse(string fastaContent)
///   • IEnumerable&lt;FastaEntry&gt; FastaParser.ParseFile(string filePath)
/// Both delegate to the SAME line-oriented state machine `ParseReader`, so the
/// string surface fully exercises the parsing contract; one file-path test below
/// additionally pins the `ParseFile` StreamReader path (a common real-world entry).
///
/// CRITICAL ENUMERATION NOTE: `Parse` is a `yield`-based iterator — NOTHING runs,
/// and no exception can surface, until the caller enumerates. Every test therefore
/// materializes the result with `.ToList()` so that any work (and any documented
/// throw) actually happens inside the assertion.
///
/// Documented, repository-specific behavior the tests pin (FASTA_Parsing.md §3.3,
/// §5.2, §5.3, §6.1; FastaParser.cs; DnaSequence.cs):
///   • DEFLINE REQUIRED (INV-01): an entry is emitted only after a '>' defline has
///     been seen AND at least one sequence character collected. A sequence block
///     with NO leading '>' therefore produces ZERO entries — the parser drops it
///     silently rather than crashing (header stays null). The fuzz test pins the
///     empty result, NOT a throw. (FASTA_Parsing.md §6.1 row 1, §2.2.)
///   • HEADER-ONLY DROPPED: a defline with no following sequence line is NOT
///     yielded, because emission is gated on `sequenceBuilder.Length > 0`
///     (FastaParser.cs lines 117 / 135). So "empty sequence" → record skipped, not
///     an empty-sequence record and not a crash. (FASTA_Parsing.md §5.2, §6.1.)
///   • DNA-ONLY MATERIALIZATION: each entry's payload is built with
///     `new DnaSequence(sequence)`, which upper-cases then REJECTS any character
///     outside A/C/G/T with an `ArgumentException` (DnaSequence.cs lines 112–124).
///     Therefore binary garbage / NUL bytes / non-DNA symbols that land in a
///     SEQUENCE body surface as a DOCUMENTED `ArgumentException` during enumeration
///     — never an unhandled IndexOutOfRange / NullReference / KeyNotFound. This is
///     the documented "rejected with ArgumentException" contract of §3.3, and the
///     fuzz tests assert exactly that exception type, never a weaker "NotThrow".
///   • WHITESPACE STRIPPED: inside sequence lines every `char.IsWhiteSpace` char is
///     dropped before materialization, so blank lines contribute nothing and a
///     record split across many lines concatenates cleanly (INV-02).
///   • HEADERS ARE NOT ALPHABET-VALIDATED: NUL bytes / control / unicode chars in a
///     DEFLINE are preserved verbatim (only the leading '>' is stripped and the
///     remainder trimmed and split on the first space/tab). A 1 MB defline is just a
///     long string token — O(n) Substring + Trim, no quadratic blow-up, no OOM.
///   • NO-TRAILING-NEWLINE: `TextReader.ReadLine` returns the final unterminated
///     line, so the last record is still parsed — the key real-world robustness case.
///   • NULL / EMPTY / WHITESPACE-ONLY input → no entries (`yield break` on
///     `string.IsNullOrWhiteSpace`, FastaParser.cs line 19).
///
/// Determinism note: random-byte payloads use a LOCAL fixed-seed `new Random(seed)`
/// (no shared static Rng), so every run is byte-for-byte reproducible. The
/// huge-header and binary-garbage tests carry `[CancelAfter]` to convert any
/// pathological hang into a deterministic failure rather than a stuck run.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: PARSE-FASTQ-001 — FASTQ parsing (FileIO)
/// ═══════════════════════════════════════════════════════════════════════════
/// Checklist: docs/checklists/03_FUZZING.md, row 65. The second FileIO-area fuzz
/// unit. Fuzz strategies exercised for THIS unit:
///   • RB  = Random Bytes — fixed-seed random characters fed into the sequence
///           AND quality fields (and a whole random-byte blob with no '@').
///   • TF  = Truncated Fields — a record that ends after only 1 / 2 / 3 of the 4
///           canonical lines (file truncated mid-record). THE key boundary case:
///           the classic IndexOutOfRange-on-a-truncated-4-line-record trap.
///   • MC  = Malformed Content — header line missing its leading '@'; separator
///           line missing its leading '+'; quality length ≠ sequence length.
///   • INJ = Injection — NUL bytes / control chars inside the quality string.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The FASTQ-parsing contract under test (the DOCUMENTED, repository-specific one)
/// ───────────────────────────────────────────────────────────────────────────
/// A canonical FASTQ record is FOUR logical lines: an '@'-header, a sequence, a
/// '+' separator, and a quality string whose length equals the sequence length;
/// quality is Phred-encoded printable ASCII (33..126).[Cock et al. 2009; NCBI SRA]
///   — docs/algorithms/FileIO/FASTQ_Parsing.md §2.1, §2.2, INV-01, INV-02.
///
/// API entry points (FastqParser.cs):
///   • IEnumerable&lt;FastqRecord&gt; FastqParser.Parse(string content, QualityEncoding)
///   • IEnumerable&lt;FastqRecord&gt; FastqParser.ParseFile(string path, QualityEncoding)
///   • IEnumerable&lt;FastqRecord&gt; FastqParser.Parse(TextReader, QualityEncoding)
/// All three delegate to the SAME `Parse(TextReader,...)` line state machine, so the
/// string surface fully exercises the contract; one file-path test additionally pins
/// the `ParseFile` StreamReader path. Like FASTA, `Parse` is a `yield`-based iterator,
/// so every test materializes with `.ToList()` to force the work to actually run.
///
/// CRITICAL CONTRACT DIVERGENCE FROM FASTA — TOLERANT, NOT STRICT.
/// FASTQ_Parsing.md §1, §5.2, §5.3, §6.2 state IN WRITING that this parser favors
/// "tolerant record assembly over strict format validation" and that "malformed
/// records can be skipped or partially assembled rather than rejected with an
/// error" — there is explicitly NO strict FASTQ validator. The DnaSequence alphabet
/// gate that made FASTA throw does NOT apply here: a `FastqRecord` stores the
/// sequence as a raw `string` (FastqParser.cs line 19–24), never materializing a
/// DnaSequence, so non-DNA / garbage sequence content does NOT throw. THEREFORE the
/// theory-correct contract these fuzz tests pin is the documented one:
///   • The BUSINESS GUARANTEE — the parser must NEVER crash, hang, throw an
///     undocumented runtime exception (IndexOutOfRange on a truncated 4-line
///     record!, NullReference, OOM), or loop forever on malformed/truncated/random/
///     injected input. EVERY fuzz test below asserts `NotThrow` + a pinned, exact
///     structural outcome, so the documented tolerance can never silently drift into
///     a crash, and a partially-assembled record's shape is nailed down precisely.
/// Documented behaviors pinned per target (FASTQ_Parsing.md §3.3, §4.1, §5.2):
///   • MISSING '@' HEADER (MC): a line not starting with '@' is SKIPPED
///     (FastqParser.cs line 96–97). A whole "record" whose header lacks '@' is
///     dropped — the parser does not crash and does not invent a record.
///   • MISSING '+' SEPARATOR (MC): the sequence loop accumulates every line until
///     one starts with '+' (line 104). With NO '+', the sequence absorbs the lines
///     that should have been the separator AND the quality, and the quality string
///     comes out EMPTY — a partially-assembled, NOT-rejected record (the documented
///     §5.2 tolerance). We pin that exact shape; the point is it does not crash.
///   • QUALITY-LENGTH ≠ SEQUENCE-LENGTH (MC): the quality loop reads lines only
///     until `qualityBuilder.Length >= sequence.Length` (line 115). A SHORT quality
///     run at EOF yields a record whose QualityString is shorter than the Sequence
///     (qual-len &lt; seq-len, NOT rejected); a quality line LONGER than the sequence
///     is appended whole (Trim never truncates), yielding qual-len &gt; seq-len. The
///     documented contract is tolerant assembly, so we pin BOTH mismatched shapes as
///     no-throw rather than asserting a rejection the parser explicitly does not do.
///   • INVALID QUALITY CHARS (INJ): `DecodeQualityScores` is `Math.Max(0, c-offset)`
///     for every char (line 178–181) — any byte, incl. NUL / control / chars outside
///     33..126, decodes (clamped at 0) without validation or crash. Pinned: a record
///     with control-char quality parses, QualityString preserved verbatim, scores
///     all ≥ 0 and never throwing.
///   • TRUNCATED RECORD (TF): every read loop guards `ReadLine() != null`, and
///     `DecodeQualityScores("")` returns an empty array (line 172–173), so a record
///     truncated to 1 / 2 / 3 lines yields a record with empty/short sequence and/or
///     empty quality and NEVER an IndexOutOfRange. This is the highest-value TF case.
///   • NULL / EMPTY input → no records (`yield break` on `IsNullOrEmpty`, line 73).
/// Determinism note: the random-byte FASTQ tests use a LOCAL fixed-seed
/// `new Random(seed)` and carry `[CancelAfter]`, exactly as the FASTA tests do.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
// I/O-bound fixture: the file-path tests share a per-fixture `_tempDir` set up in
// [SetUp] and removed in [TearDown]. Under the assembly-wide
// `Parallelizable(ParallelScope.Children)` policy a single fixture instance is shared
// across parallel test methods, so a concurrent [TearDown] can delete the temp dir
// another test is mid-write into (DirectoryNotFoundException). Per the documented
// house rule in Parallelization.cs ("I/O-bound fixtures should be marked
// [NonParallelizable]") this fixture runs its tests serially.
[NonParallelizable]
public class FileIoFuzzTests
{
    #region Helpers — temp-dir management for the ParseFile (file-path) surface

    private string _tempDir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(),
            "SeqeronFastaFuzz_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_tempDir);
    }

    [TearDown]
    public void TearDown()
    {
        try
        {
            if (Directory.Exists(_tempDir))
                Directory.Delete(_tempDir, recursive: true);
        }
        catch (IOException)
        {
            // Best-effort cleanup; a leaked temp file must not fail the suite.
        }
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp file and returns its path, so the file-path tests control byte
    /// layout exactly (in particular: whether a trailing newline exists).
    /// </summary>
    private string WriteTempFasta(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".fasta");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp `.fastq` file and returns its path, so the FASTQ file-path tests
    /// control byte layout exactly (in particular: truncation mid-record).
    /// </summary>
    private string WriteTempFastq(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".fastq");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    #endregion

    #region PARSE-FASTA-001 — FASTA parsing

    #region Positive sanity — a valid multi-record FASTA parses to the correct ids/sequences

    /// <summary>
    /// Positive control: a well-formed two-record FASTA (the worked example from
    /// FASTA_Parsing.md §7.1) must parse to exactly two entries with the correct
    /// id / description / sequence. If this fails, the fuzz suite is meaningless —
    /// it proves the parser is actually wired up and the happy path is intact.
    /// </summary>
    [Test]
    public void Parse_ValidMultiRecord_ParsesToCorrectEntries()
    {
        const string fasta =
            ">seq1 First sequence\n" +
            "ATGCATGC\n" +
            ">seq2 Second sequence\n" +
            "GGCCTTAA\n";

        List<FastaEntry> entries = FastaParser.Parse(fasta).ToList();

        entries.Should().HaveCount(2, "both well-formed records must be emitted");

        entries[0].Id.Should().Be("seq1");
        entries[0].Description.Should().Be("First sequence");
        entries[0].Sequence.Sequence.Should().Be("ATGCATGC");

        entries[1].Id.Should().Be("seq2");
        entries[1].Description.Should().Be("Second sequence");
        entries[1].Sequence.Sequence.Should().Be("GGCCTTAA");
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: a multi-line, multi-record FASTA
    /// written to disk (sequence wrapped across lines, exactly as real files do)
    /// parses via `ParseFile` to the same correct entries. This pins the StreamReader
    /// path and the multi-line concatenation invariant (INV-02).
    /// </summary>
    [Test]
    public void ParseFile_ValidMultiLineMultiRecord_ParsesToCorrectEntries()
    {
        const string fasta =
            ">chrA alpha contig\n" +
            "ATGC\n" +
            "ATGC\n" +
            ">chrB\n" +
            "GGGGCCCC\n";
        string path = WriteTempFasta(fasta);

        List<FastaEntry> entries = FastaParser.ParseFile(path).ToList();

        entries.Should().HaveCount(2);
        entries[0].Id.Should().Be("chrA");
        entries[0].Description.Should().Be("alpha contig");
        entries[0].Sequence.Sequence.Should().Be("ATGCATGC",
            "wrapped sequence lines must concatenate (INV-02)");
        entries[1].Id.Should().Be("chrB");
        entries[1].Description.Should().BeNull("a description-less defline yields a null Description");
        entries[1].Sequence.Sequence.Should().Be("GGGGCCCC");
    }

    #endregion

    #region MC — Missing '>': a defline-less sequence block yields NO entries (no crash)

    /// <summary>
    /// MC / TF: input that is a bare sequence block with NO leading '>' defline.
    /// Per INV-01 (FASTA_Parsing.md §2.2, §6.1 row 1) an entry is emitted only after
    /// a defline has been seen, so `header` never leaves null and the buffered
    /// characters are dropped. The documented, theory-correct outcome is ZERO
    /// entries — a clean parse, NOT an exception and certainly NOT a crash. We pin
    /// the empty result exactly so this silent-drop contract cannot drift.
    /// </summary>
    [Test]
    public void Parse_MissingHeaderMarker_YieldsNoEntries()
    {
        const string headerless = "ATGCATGCATGC\nGGCCTTAAGGCC\n";

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(headerless).ToList();

        act.Should().NotThrow("a defline-less block is dropped, not a crash");
        entries.Should().BeEmpty(
            "no entry can be emitted without a preceding '>' defline (INV-01)");
    }

    #endregion

    #region TF — Empty sequence: a header with no sequence line is dropped (no crash)

    /// <summary>
    /// TF / MC: a defline immediately followed by another defline (and a final
    /// defline with nothing after it) — i.e. headers with EMPTY sequence bodies.
    /// Emission is gated on `sequenceBuilder.Length > 0` (FastaParser.cs line 117 /
    /// 135), so header-only records are DROPPED (FASTA_Parsing.md §5.2, §6.1). The
    /// only record with sequence content survives; no empty-sequence record is
    /// emitted and nothing crashes. We verify the surviving record AND the dropping.
    /// </summary>
    [Test]
    public void Parse_HeaderWithEmptySequence_DropsEmptyRecords()
    {
        const string fasta =
            ">empty1\n" +     // header, then immediately another header → empty seq
            ">withSeq data\n" +
            "ACGT\n" +
            ">empty2\n";      // trailing header, no sequence at all

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(fasta).ToList();

        act.Should().NotThrow("empty-sequence records are dropped, not a crash");
        entries.Should().ContainSingle("only the one record with sequence content survives")
            .Which.Id.Should().Be("withSeq");
        entries[0].Sequence.Sequence.Should().Be("ACGT");
    }

    #endregion

    #region RB — Binary garbage / random bytes: documented ArgumentException, never a crash

    /// <summary>
    /// RB: a valid defline followed by a sequence body of fixed-seed RANDOM bytes
    /// (the full 0x00–0xFF range as chars). Because the body is materialized via
    /// `new DnaSequence(...)`, any non-A/C/G/T character is REJECTED with a
    /// DOCUMENTED `ArgumentException` (DnaSequence.cs lines 112–124; FASTA_Parsing.md
    /// §3.3). That is the contract: random garbage is rejected deterministically, it
    /// is NEVER an unhandled IndexOutOfRange / NullReference / OutOfMemory and never a
    /// silently-corrupt record. We assert the exact exception type. `[CancelAfter]`
    /// guards against any pathological hang while scanning the garbage.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Parse_BinaryGarbageSequenceBody_ThrowsDocumentedArgumentException(
        CancellationToken token)
    {
        var rng = new Random(0xFA57A);             // local fixed seed — fully reproducible
        var sb = new StringBuilder(">garbage record\n");
        for (int i = 0; i < 4096; i++)
            sb.Append((char)rng.Next(0, 256));     // arbitrary bytes incl. control chars
        sb.Append('\n');
        string fasta = sb.ToString();

        // Enumerate eagerly so the DnaSequence materialization (and its throw) runs.
        var act = () => FastaParser.Parse(fasta).ToList();

        act.Should().Throw<ArgumentException>(
            "non-A/C/G/T sequence content is rejected by DnaSequence, the documented contract")
           .And.Message.Should().Contain("nucleotide");
        token.IsCancellationRequested.Should().BeFalse("parsing garbage must not hang");
    }

    /// <summary>
    /// RB / robustness companion: pure random bytes that do NOT begin with '>' must
    /// be handled deterministically. The first random byte is almost never '>' (and
    /// we force it to be a non-'>' byte), so with no defline the whole block is a
    /// defline-less sequence buffer that is never emitted (INV-01) → ZERO entries,
    /// no crash, no hang. This pins the "random bytes that aren't even FASTA" case.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Parse_PureRandomBytesNoDefline_YieldsNoEntriesAndDoesNotCrash(
        CancellationToken token)
    {
        var rng = new Random(0xB10C);
        var sb = new StringBuilder();
        for (int i = 0; i < 8192; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '>') c = 'x';   // guarantee no defline is ever introduced
            if (c == '\n' || c == '\r') c = 'y'; // keep it a single buffered "line"
            sb.Append(c);
        }

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(sb.ToString()).ToList();

        act.Should().NotThrow("defline-less random bytes are buffered and dropped, not a crash");
        entries.Should().BeEmpty("with no '>' line no entry can ever be emitted (INV-01)");
        token.IsCancellationRequested.Should().BeFalse("parsing random bytes must not hang");
    }

    #endregion

    #region INJ — Embedded NUL bytes: handled deterministically

    /// <summary>
    /// INJ: NUL bytes (\0) injected into BOTH the defline and the sequence body.
    /// • In the DEFLINE, '\0' is NOT alphabet-validated — it is preserved verbatim
    ///   as part of the header token (only '>' is stripped, the rest trimmed/split).
    ///   So the parser must not crash on a NUL in the header.
    /// • In the SEQUENCE body, '\0' is not whitespace, so it reaches DnaSequence and
    ///   is REJECTED with the documented `ArgumentException` (it is not A/C/G/T).
    /// We pin BOTH halves: a NUL-bearing sequence throws the documented exception,
    /// while a NUL purely in the header (with a clean sequence) parses cleanly and
    /// the NUL survives in the id. Deterministic either way — never an unhandled crash.
    /// </summary>
    [Test]
    public void Parse_NullBytesInSequence_ThrowsDocumentedArgumentException()
    {
        const string fastaWithNulInSeq = ">rec\nAC\0GT\n";

        var act = () => FastaParser.Parse(fastaWithNulInSeq).ToList();

        act.Should().Throw<ArgumentException>(
            "a NUL byte in the sequence body is a non-DNA char rejected by DnaSequence");
    }

    /// <summary>
    /// INJ companion: a NUL byte injected only into the DEFLINE (header), with a
    /// clean A/C/G/T sequence, must parse deterministically — the header is not
    /// alphabet-validated, so the record is emitted and the NUL survives verbatim in
    /// the parsed id token. No crash, no silent loss of the record.
    /// </summary>
    [Test]
    public void Parse_NullByteInHeaderOnly_ParsesAndPreservesNul()
    {
        const string fasta = ">id\0X clean\nACGT\n";

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(fasta).ToList();

        act.Should().NotThrow("a NUL in the header is preserved, not validated, so it must not crash");
        entries.Should().ContainSingle();
        entries[0].Id.Should().Be("id\0X", "the first whitespace-delimited token, NUL included");
        entries[0].Sequence.Sequence.Should().Be("ACGT");
    }

    #endregion

    #region OVF — Huge header (1 MB defline): completes in linear time, no OOM/quadratic blow-up

    /// <summary>
    /// OVF / TF: a single defline whose header is 1 MB long (a header-only record,
    /// no sequence). The parser strips the leading '>', trims, and splits on the
    /// first space/tab — all O(n) over the header length, no quadratic re-scan, no
    /// unbounded allocation. Because the record has no sequence body it is dropped
    /// (INV / header-only rule) → ZERO entries. The point of THIS test is that the
    /// huge header is processed without OutOfMemory or a quadratic stall, guarded by
    /// `[CancelAfter]`. We add a sequence after the giant header to ALSO prove the
    /// parser recovers and still emits the following real record.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Parse_HugeHeader_CompletesWithoutBlowUpAndRecovers(CancellationToken token)
    {
        const int headerLen = 1_000_000;                 // 1 MB defline
        string hugeHeader = new string('H', headerLen);
        string fasta = ">" + hugeHeader + "\n" +         // header-only → dropped
                       ">realId real desc\n" +
                       "ACGTACGT\n";

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(fasta).ToList();

        act.Should().NotThrow("a 1 MB header must be processed in linear time, no OOM/quadratic blow-up");
        entries.Should().ContainSingle("the header-only giant record is dropped; the real record survives")
            .Which.Id.Should().Be("realId");
        entries[0].Sequence.Sequence.Should().Be("ACGTACGT");
        token.IsCancellationRequested.Should().BeFalse("a huge header must not stall the parser");
    }

    #endregion

    #region TF — No trailing newline: the last record is still parsed (key real-world case)

    /// <summary>
    /// TF: a FASTA whose final record has NO trailing newline — extremely common in
    /// real files and a classic off-by-one trap. `TextReader.ReadLine` returns the
    /// last unterminated line, and the end-of-input flush emits the final buffered
    /// record (FastaParser.cs lines 135–138). So the last record MUST still parse.
    /// We assert it over BOTH surfaces: the in-memory `Parse(string)` and the
    /// `ParseFile` StreamReader path (with a temp file written WITHOUT a trailing
    /// newline), since the file path is where this case bites in practice.
    /// </summary>
    [Test]
    public void Parse_NoTrailingNewline_StillParsesLastRecord()
    {
        const string noNewline = ">a first\nACGT\n>b second\nTTTT"; // no final '\n'

        List<FastaEntry> entries = FastaParser.Parse(noNewline).ToList();

        entries.Should().HaveCount(2, "the unterminated final record must still be emitted");
        entries[1].Id.Should().Be("b");
        entries[1].Description.Should().Be("second");
        entries[1].Sequence.Sequence.Should().Be("TTTT");
    }

    /// <summary>
    /// TF over the FILE-PATH surface: the same no-trailing-newline case written to
    /// disk and read back via `ParseFile`, pinning that the StreamReader flush emits
    /// the final unterminated record exactly as the string surface does.
    /// </summary>
    [Test]
    public void ParseFile_NoTrailingNewline_StillParsesLastRecord()
    {
        const string noNewline = ">x\nAAAA\n>y last\nCCCC"; // file ends mid-record
        string path = WriteTempFasta(noNewline);

        List<FastaEntry> entries = FastaParser.ParseFile(path).ToList();

        entries.Should().HaveCount(2);
        entries[1].Id.Should().Be("y");
        entries[1].Description.Should().Be("last");
        entries[1].Sequence.Sequence.Should().Be("CCCC",
            "ParseFile must flush the final unterminated record");
    }

    #endregion

    #endregion

    #region PARSE-FASTQ-001 — FASTQ parsing

    #region Positive sanity — a valid multi-record FASTQ parses to correct id/seq/quality

    /// <summary>
    /// Positive control: a well-formed two-record FASTQ must parse to exactly two
    /// records with the correct id / description / sequence / quality string, and the
    /// per-record quality-string length must equal the sequence length (INV-01). If
    /// this fails the whole FASTQ fuzz suite is meaningless — it proves the happy path
    /// is wired up before we throw garbage at it.
    /// </summary>
    [Test]
    public void ParseFastq_ValidMultiRecord_ParsesToCorrectRecords()
    {
        const string fastq =
            "@SEQ1 first read\n" +
            "ACGTACGT\n" +
            "+\n" +
            "IIIIIIII\n" +
            "@SEQ2\n" +
            "GGCCTTAA\n" +
            "+\n" +
            "!!!!!!!!\n";

        List<FastqParser.FastqRecord> records = FastqParser.Parse(fastq).ToList();

        records.Should().HaveCount(2, "both well-formed records must be emitted");

        records[0].Id.Should().Be("SEQ1");
        records[0].Description.Should().Be("first read");
        records[0].Sequence.Should().Be("ACGTACGT");
        records[0].QualityString.Should().Be("IIIIIIII");
        records[0].QualityScores.Should().HaveCount(records[0].Sequence.Length,
            "valid FASTQ has one quality score per base (INV-01)");
        records[0].QualityScores.Should().AllSatisfy(s => s.Should().BeGreaterThanOrEqualTo(0));

        records[1].Id.Should().Be("SEQ2");
        records[1].Description.Should().Be("", "a description-less header yields an empty Description");
        records[1].Sequence.Should().Be("GGCCTTAA");
        records[1].QualityString.Should().Be("!!!!!!!!");
        records[1].QualityScores.Should().HaveCount(8);
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: a valid two-record FASTQ written
    /// to disk parses via `ParseFile` to the same correct records, pinning the
    /// StreamReader path (the common real-world entry point).
    /// </summary>
    [Test]
    public void ParseFileFastq_ValidMultiRecord_ParsesToCorrectRecords()
    {
        const string fastq =
            "@readA alpha\n" +
            "ACGT\n" +
            "+\n" +
            "IIII\n" +
            "@readB\n" +
            "TTTT\n" +
            "+\n" +
            "####\n";
        string path = WriteTempFastq(fastq);

        List<FastqParser.FastqRecord> records = FastqParser.ParseFile(path).ToList();

        records.Should().HaveCount(2);
        records[0].Id.Should().Be("readA");
        records[0].Description.Should().Be("alpha");
        records[0].Sequence.Should().Be("ACGT");
        records[0].QualityString.Should().Be("IIII");
        records[1].Id.Should().Be("readB");
        records[1].Sequence.Should().Be("TTTT");
        records[1].QualityString.Should().Be("####");
    }

    #endregion

    #region MC — Missing '@' header: the malformed block is skipped, never a crash

    /// <summary>
    /// MC: the first "record" has a header line that does NOT begin with '@' (the '@'
    /// is missing). Per FastqParser.cs line 96–97 every line that does not start with
    /// '@' is SKIPPED, so the entire malformed block (header-less + its seq/+/qual
    /// lines, none of which start with '@') is dropped and NO record is invented for
    /// it. The following well-formed record is still parsed. The documented tolerant
    /// contract is "skip, do not crash" — we pin the surviving record exactly and
    /// assert no throw, so the skip cannot silently become a crash.
    /// </summary>
    [Test]
    public void ParseFastq_MissingAtHeaderMarker_SkipsMalformedBlockNoCrash()
    {
        const string fastq =
            "SEQ_no_at\n" +   // header missing '@' → skipped
            "ACGT\n" +        // these lines also lack '@' → skipped
            "+\n" +
            "IIII\n" +
            "@good\n" +       // a clean record after the garbage
            "GGCC\n" +
            "+\n" +
            "JJJJ\n";

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("a header without '@' is skipped, not a crash (tolerant parser)");
        records.Should().ContainSingle("only the well-formed '@good' record can be emitted")
            .Which.Id.Should().Be("good");
        records[0].Sequence.Should().Be("GGCC");
        records[0].QualityString.Should().Be("JJJJ");
    }

    #endregion

    #region MC — Missing '+' separator: sequence absorbs the rest, quality empty, no crash

    /// <summary>
    /// MC: a record with NO '+' separator line. The sequence-accumulation loop
    /// (FastqParser.cs line 104) reads lines until one begins with '+'; with no such
    /// line it absorbs BOTH the intended sequence AND the intended quality line, then
    /// hits EOF, leaving the quality string EMPTY. This is the documented §5.2
    /// "partially assembled rather than rejected" behavior — a tolerant parse, NOT a
    /// crash and NOT an exception. We pin the exact partially-assembled shape so the
    /// documented tolerance is nailed down and can never drift into an
    /// IndexOutOfRange/NullReference.
    /// </summary>
    [Test]
    public void ParseFastq_MissingPlusSeparator_AbsorbsIntoSequenceNoCrash()
    {
        const string fastq =
            "@rec\n" +
            "ACGT\n" +   // intended sequence
            "IIII\n";    // intended quality — but with no '+' it is absorbed too

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("a missing '+' separator is tolerated (partial assembly), not a crash");
        records.Should().ContainSingle("the record is partially assembled, not rejected");
        records[0].Id.Should().Be("rec");
        records[0].Sequence.Should().Be("ACGTIIII",
            "with no '+' the sequence loop absorbs the would-be quality line too (§5.2)");
        records[0].QualityString.Should().BeEmpty(
            "nothing is left for the quality loop after the sequence loop hits EOF");
        records[0].QualityScores.Should().BeEmpty("an empty quality string decodes to no scores");
    }

    #endregion

    #region MC — Quality length ≠ sequence length: tolerant partial assembly, never a crash

    /// <summary>
    /// MC (qual-len &lt; seq-len): the quality run is SHORTER than the sequence and the
    /// input ends before the quality loop fills up. The loop (line 115) stops at EOF,
    /// so the record carries a QualityString shorter than its Sequence — the parser
    /// does NOT reject the corruption (documented tolerant behavior, §5.2/§6.2), it
    /// emits the mismatched record. The business guarantee under fuzzing is that this
    /// must not crash and must not produce more scores than quality chars. We pin the
    /// exact short-quality shape so the documented tolerance is explicit and stable.
    /// </summary>
    [Test]
    public void ParseFastq_QualityShorterThanSequence_TolerantNoCrash()
    {
        const string fastq =
            "@rec\n" +
            "ACGTACGT\n" +  // 8-base sequence
            "+\n" +
            "II\n";          // only 2 quality chars, then EOF

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("a too-short quality run is tolerated at EOF, not a crash");
        records.Should().ContainSingle();
        records[0].Sequence.Should().Be("ACGTACGT");
        records[0].QualityString.Should().Be("II",
            "the quality loop stops at EOF with quality shorter than the sequence (no rejection)");
        records[0].QualityString.Length.Should().BeLessThan(records[0].Sequence.Length,
            "this is the documented qual-len < seq-len corruption the tolerant parser accepts");
        records[0].QualityScores.Should().HaveCount(records[0].QualityString.Length,
            "scores are decoded per quality char, never more than the quality length");
    }

    /// <summary>
    /// MC (qual-len &gt; seq-len): the quality LINE is LONGER than the sequence. The
    /// quality loop appends the whole trimmed line in one shot (Trim never truncates),
    /// so `qualityBuilder.Length` overshoots the sequence length and the QualityString
    /// comes out LONGER than the Sequence. Again tolerant assembly, not rejection — we
    /// pin the over-long shape and assert no crash. The following clean record proves
    /// the parser recovers after the over-long line is fully consumed.
    /// </summary>
    [Test]
    public void ParseFastq_QualityLongerThanSequence_TolerantNoCrash()
    {
        const string fastq =
            "@rec1\n" +
            "ACGT\n" +          // 4-base sequence
            "+\n" +
            "IIIIIIIIIIII\n" +  // 12 quality chars (over-long, appended whole)
            "@rec2\n" +
            "TTTT\n" +
            "+\n" +
            "JJJJ\n";

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("an over-long quality line is tolerated, not a crash");
        records.Should().HaveCount(2, "the parser recovers and still emits the following record");
        records[0].Sequence.Should().Be("ACGT");
        records[0].QualityString.Length.Should().BeGreaterThan(records[0].Sequence.Length,
            "the whole over-long quality line is appended (Trim does not truncate)");
        records[1].Id.Should().Be("rec2", "the parser recovers after the over-long quality line");
        records[1].QualityString.Should().Be("JJJJ");
    }

    #endregion

    #region INJ — Invalid / control quality chars: decoded without validation, never a crash

    /// <summary>
    /// INJ: quality characters OUTSIDE the canonical Phred printable range 33..126 —
    /// here NUL (\0) and other control bytes injected into the quality string (with a
    /// matching-length sequence). `DecodeQualityScores` is `Math.Max(0, c - offset)`
    /// for every char (line 178–181): it never validates the alphabet, so control /
    /// out-of-range chars decode (clamped at ≥ 0) without throwing. The documented
    /// contract is "decode, do not reject"; the fuzz guarantee is no crash. We pin
    /// that the record parses, the QualityString is preserved verbatim (NUL included),
    /// and every decoded score is ≥ 0.
    /// </summary>
    [Test]
    public void ParseFastq_InvalidControlCharQuality_DecodesWithoutCrash()
    {
        // 4-base sequence; 4 quality chars, all control bytes below ASCII 33.
        const string fastq = "@rec\nACGT\n+\n\0\n";

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("out-of-range/control quality chars are decoded, not validated (no crash)");
        records.Should().ContainSingle();
        records[0].Sequence.Should().Be("ACGT");
        records[0].QualityString.Should().Be("\0",
            "control-char quality is preserved verbatim, not sanitized");
        records[0].QualityScores.Should().HaveCount(4);
        records[0].QualityScores.Should().AllSatisfy(s => s.Should().BeGreaterThanOrEqualTo(0),
            "DecodeQualityScores clamps every value to ≥ 0, even for sub-offset bytes");
    }

    #endregion

    #region TF — Truncated record (1/2/3 of 4 lines): never IndexOutOfRange, defined result

    /// <summary>
    /// TF — THE key boundary case: a file truncated after only 1, 2, or 3 of the four
    /// canonical FASTQ lines. The classic trap is an IndexOutOfRange when code assumes
    /// four lines are always present. Every read loop here guards `ReadLine() != null`
    /// and `DecodeQualityScores("")` returns an empty array (line 172–173), so each
    /// truncation yields a well-DEFINED record (empty/short sequence and/or empty
    /// quality) and NEVER an unhandled exception. We assert no-throw plus the exact
    /// shape for each truncation point.
    /// </summary>
    [Test]
    public void ParseFastq_TruncatedRecords_NeverIndexOutOfRange()
    {
        // 1 line only: header, then EOF.
        const string oneLine = "@only_header\n";
        // 2 lines: header + sequence, then EOF (no '+', no quality).
        const string twoLines = "@hdr\nACGT\n";
        // 3 lines: header + sequence + '+' separator, then EOF (no quality).
        const string threeLines = "@hdr\nACGT\n+\n";

        List<FastqParser.FastqRecord> r1 = new(), r2 = new(), r3 = new();
        var act1 = () => r1 = FastqParser.Parse(oneLine).ToList();
        var act2 = () => r2 = FastqParser.Parse(twoLines).ToList();
        var act3 = () => r3 = FastqParser.Parse(threeLines).ToList();

        act1.Should().NotThrow("a header-only truncation must not IndexOutOfRange");
        act2.Should().NotThrow("a header+sequence truncation must not IndexOutOfRange");
        act3.Should().NotThrow("a header+sequence+'+' truncation must not IndexOutOfRange");

        // 1 line: header consumed; sequence loop hits EOF → empty seq; quality empty.
        r1.Should().ContainSingle();
        r1[0].Id.Should().Be("only_header");
        r1[0].Sequence.Should().BeEmpty();
        r1[0].QualityString.Should().BeEmpty();
        r1[0].QualityScores.Should().BeEmpty();

        // 2 lines: no '+', so the sequence absorbs "ACGT"; quality empty at EOF.
        r2.Should().ContainSingle();
        r2[0].Sequence.Should().Be("ACGT");
        r2[0].QualityString.Should().BeEmpty();

        // 3 lines: '+' stops the sequence loop at "ACGT"; quality loop hits EOF → empty.
        r3.Should().ContainSingle();
        r3[0].Sequence.Should().Be("ACGT");
        r3[0].QualityString.Should().BeEmpty(
            "the quality loop hits EOF immediately, leaving an empty (not crashed) quality");
        r3[0].QualityScores.Should().BeEmpty();
    }

    /// <summary>
    /// TF over the FILE-PATH surface: a multi-record FASTQ whose FINAL record is
    /// truncated mid-record (only 3 of the 4 lines, no trailing newline) is written to
    /// disk and read via `ParseFile`. The StreamReader flush must emit the complete
    /// first record AND the truncated final record without an IndexOutOfRange — the
    /// real-world shape where this bug actually bites.
    /// </summary>
    [Test]
    public void ParseFileFastq_FinalRecordTruncated_StillParsesNoCrash()
    {
        const string fastq =
            "@full\n" +
            "ACGT\n" +
            "+\n" +
            "IIII\n" +
            "@truncated\n" +
            "GGCC\n" +
            "+";              // file ends here: no quality line, no trailing newline
        string path = WriteTempFastq(fastq);

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.ParseFile(path).ToList();

        act.Should().NotThrow("a truncated final record on disk must not crash ParseFile");
        records.Should().HaveCount(2, "both the full and the truncated records are emitted");
        records[0].Id.Should().Be("full");
        records[0].QualityString.Should().Be("IIII");
        records[1].Id.Should().Be("truncated");
        records[1].Sequence.Should().Be("GGCC");
        records[1].QualityString.Should().BeEmpty("the truncated record has no quality at EOF");
    }

    #endregion

    #region RB — Random bytes: handled deterministically, no crash, no hang

    /// <summary>
    /// RB: a valid '@' header followed by fixed-seed RANDOM bytes (full 0x00–0xFF
    /// range) in the SEQUENCE/quality body. Unlike FASTA, the FASTQ sequence is stored
    /// as a raw string and never alphabet-validated, so random garbage does NOT throw —
    /// the documented tolerant contract. The fuzz guarantee is that the parser consumes
    /// the garbage deterministically: no crash, no hang, and any decoded quality scores
    /// are ≥ 0. `[CancelAfter]` converts any pathological hang into a deterministic
    /// failure; the LOCAL fixed seed makes the run byte-for-byte reproducible.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseFastq_RandomByteBody_DoesNotCrashOrHang(CancellationToken token)
    {
        var rng = new Random(0xFA5701);            // local fixed seed — fully reproducible
        var sb = new StringBuilder("@garbage record\n");
        for (int i = 0; i < 4096; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '\n' || c == '\r') c = 'N';   // keep it inside one logical body line
            sb.Append(c);
        }
        sb.Append('\n');
        string fastq = sb.ToString();

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow(
            "random bytes in a FASTQ body are tolerated (raw string, no alphabet gate), not a crash");
        records.Should().ContainSingle("one '@' header → one (tolerantly assembled) record");
        records[0].QualityScores.Should().AllSatisfy(s => s.Should().BeGreaterThanOrEqualTo(0),
            "every decoded score is clamped to ≥ 0 regardless of the byte");
        token.IsCancellationRequested.Should().BeFalse("parsing garbage must not hang");
    }

    /// <summary>
    /// RB companion: a whole blob of fixed-seed random bytes with NO '@' header line.
    /// Every line fails the `StartsWith('@')` gate (line 96–97) and is skipped, so the
    /// parser yields ZERO records — no crash, no hang. We force out '@', '\n' and '\r'
    /// so the blob is genuinely one header-less garbage line. `[CancelAfter]` guards
    /// against any stall while scanning the garbage.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseFastq_PureRandomBytesNoHeader_YieldsNoRecordsNoCrash(CancellationToken token)
    {
        var rng = new Random(0xB10C);
        var sb = new StringBuilder();
        for (int i = 0; i < 8192; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '@') c = 'x';                  // guarantee no header is ever introduced
            if (c == '\n' || c == '\r') c = 'y';    // keep it a single header-less line
            sb.Append(c);
        }

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(sb.ToString()).ToList();

        act.Should().NotThrow("header-less random bytes are skipped line by line, not a crash");
        records.Should().BeEmpty("with no '@' line no record can ever be emitted");
        token.IsCancellationRequested.Should().BeFalse("parsing random bytes must not hang");
    }

    #endregion

    #endregion
}
