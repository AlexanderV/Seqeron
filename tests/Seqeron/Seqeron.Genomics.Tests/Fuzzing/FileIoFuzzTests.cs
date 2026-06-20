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
/// Fuzz tests for the FileIO area — FASTA parsing (PARSE-FASTA-001).
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
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
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
}
