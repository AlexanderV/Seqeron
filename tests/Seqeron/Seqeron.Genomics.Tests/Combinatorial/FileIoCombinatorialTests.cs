namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the FileIO area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("FileIO")]
public class FileIoCombinatorialTests
{
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static string Wrap(string seq, int width) =>
        width <= 0
            ? seq
            : string.Join("\n", Enumerable.Range(0, (seq.Length + width - 1) / width)
                .Select(i => seq.Substring(i * width, Math.Min(width, seq.Length - i * width))));

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PARSE-FASTA-001 — FASTA parsing (FileIO)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 64.
    // Spec: tests/TestSpecs/PARSE-FASTA-001.md (canonical FastaParser.Parse).
    // Dimensions: multiRecord(2) × lineWrap(3) × headerStyle(3). Grid 2×3×3 = 18.
    //
    // Model (Pearson FASTA): each record is a '>' header line (id = first token, description =
    // the rest) followed by sequence lines; the residues are the concatenation of all sequence
    // lines with whitespace removed — so the parsed sequence is INVARIANT to line wrapping.
    //
    // The combinatorial point: record count, line-wrap width and header style interact, yet the
    // parser recovers every record's id, description and (wrap-independent) sequence in order.
    // ═══════════════════════════════════════════════════════════════════════

    public enum HeaderStyle { IdOnly, IdDescription, PipedAccession }

    private static (string Id, string? Desc, string Header) Header(int r, HeaderStyle style) => style switch
    {
        HeaderStyle.IdOnly => ($"seq{r}", null, $"seq{r}"),
        HeaderStyle.IdDescription => ($"seq{r}", $"description number {r}", $"seq{r} description number {r}"),
        // NCBI piped style has no space ⇒ the whole token is the id (parser splits on whitespace).
        _ => ($"sp|ACC{r}|seq{r}", null, $"sp|ACC{r}|seq{r}"),
    };

    [Test, Combinatorial]
    public void ParseFasta_RecoversRecords_AcrossWrapAndHeaderStyle(
        [Values(1, 3)] int nRecords,
        [Values(10, 25, 0)] int lineWrap, // 0 = single line
        [Values(HeaderStyle.IdOnly, HeaderStyle.IdDescription, HeaderStyle.PipedAccession)] HeaderStyle style)
    {
        var records = Enumerable.Range(0, nRecords)
            .Select(r => (Hdr: Header(r, style), Seq: DiverseDna(40 + r * 7, (uint)(0x100 + r))))
            .ToList();

        string content = string.Concat(records.Select(rec => $">{rec.Hdr.Header}\n{Wrap(rec.Seq, lineWrap)}\n"));

        var parsed = FastaParser.Parse(content).ToList();

        parsed.Should().HaveCount(nRecords, "every record is parsed, in order");
        for (int r = 0; r < nRecords; r++)
        {
            parsed[r].Id.Should().Be(records[r].Hdr.Id, "id is the first header token");
            parsed[r].Description.Should().Be(records[r].Hdr.Desc, "description is the header remainder");
            parsed[r].Sequence.Sequence.Should().Be(records[r].Seq, "line wrapping does not change the residues");
        }
    }

    /// <summary>
    /// Interaction witness: the same record wrapped at different widths (and unwrapped) parses to
    /// the identical sequence — parsing is line-wrap invariant.
    /// </summary>
    [Test]
    public void ParseFasta_WrapInvariance()
    {
        string seq = DiverseDna(120, 7u);
        string Parse(int w) => FastaParser.Parse($">s\n{Wrap(seq, w)}\n").Single().Sequence.Sequence;

        Parse(10).Should().Be(seq);
        Parse(60).Should().Be(seq);
        Parse(0).Should().Be(seq);
    }

    /// <summary>
    /// Interaction witness: a write→parse round-trip preserves ids and sequences for a
    /// multi-record file at any output line width.
    /// </summary>
    [Test]
    public void ParseFasta_RoundTripsThroughToFasta()
    {
        var entries = Enumerable.Range(0, 3)
            .Select(r => new FastaEntry($"id{r}", $"desc {r}", new DnaSequence(DiverseDna(50 + r * 11, (uint)(900 + r))))).ToList();

        string text = FastaParser.ToFasta(entries, lineWidth: 30);
        var parsed = FastaParser.Parse(text).ToList();

        parsed.Select(p => p.Id).Should().Equal(entries.Select(e => e.Id));
        parsed.Select(p => p.Sequence.Sequence).Should().Equal(entries.Select(e => e.Sequence.Sequence));
    }
}
