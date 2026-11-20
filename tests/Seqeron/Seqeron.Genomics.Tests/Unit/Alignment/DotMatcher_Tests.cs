namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A2-28d, B09 F152): <see cref="DotMatcher.Compute"/> = EMBOSS dotmatcher 6.6.0. Oracle (TestData/Alignment/dotmatcher_oracle.tsv, scripts/alignment/generate_dotmatcher_oracle.py): 600 related / unrelated
/// DNA (EDNAFULL, EDNASIMPLE) and protein (EBLOSUM62, EPAM250, EBLOSUM30) pairs, windows 3-20, thresholds around the identity score; the <c>Line</c> segments of <c>-graph data</c> in drawing order (x = the second sequence, y = the first: measured).
/// </summary>
[TestFixture]
public class DotMatcher_Tests
{
    private sealed record Row(string Kind, string Matrix, int Window, int Threshold, string Seq1, string Seq2, string Segments);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.dotmatcher_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], int.Parse(f[2]), int.Parse(f[3]), f[4], f[5], f[6]));
        }
        return rows;
    }

    private static string Format(IReadOnlyList<DotMatcherSegment> segments) =>
        segments.Count == 0 ? "-" : string.Join(" ", segments.Select(s => $"{s.Start2}:{s.Start1}:{s.End2}:{s.End1}"));   // the plot: x = the second sequence, y = the first

    [Test]
    public void Compute_EqualsDotmatcher_AllRows()
    {
        var rows = Oracle.Value;
        Assert.That(rows, Has.Count.EqualTo(600));
        var failures = new List<string>();
        int withSegments = 0, protein = 0, segmentsTotal = 0, simple = 0;
        foreach (var r in rows)
        {
            var actual = Format(DotMatcher.Compute(r.Seq1, r.Seq2, r.Window, r.Threshold, SubstitutionMatrix.Get(r.Matrix)));
            if (actual != r.Segments)
                failures.Add($"{r.Matrix} w={r.Window} t={r.Threshold}\n{r.Seq1}\n{r.Seq2}\n expected {r.Segments}\n actual   {actual}");
            if (r.Segments != "-") { withSegments++; segmentsTotal += r.Segments.Split(' ').Length; }
            if (r.Kind == "prot") protein++;
            if (r.Matrix == "EDNASIMPLE") simple++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        Assert.That(withSegments, Is.GreaterThan(300));
        Assert.That(segmentsTotal, Is.GreaterThan(2000));
        Assert.That(protein, Is.GreaterThan(200));
        Assert.That(simple, Is.GreaterThan(50));
    }

    [Test]
    public void DefaultMatrix_IsChosenByTheSequenceType_AndShortSequencesGiveNoSegments()
    {
        // identical 12-mers: one segment over the whole diagonal (the same one from both loops)
        var dna = DotMatcher.Compute("ACGTACGTACGT", "ACGTACGTACGT", 6, 20);
        Assert.That(dna, Does.Contain(new DotMatcherSegment(1, 1, 12, 12)));
        var protein = DotMatcher.Compute("MKVLAAGIVGLL", "MKVLAAGIVGLL", 6, 20);
        Assert.That(protein, Does.Contain(new DotMatcherSegment(1, 1, 12, 12)));
        Assert.That(DotMatcher.Compute("ACGT", "ACGTACGT", 10, 5), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => DotMatcher.Compute("ACGT", "ACGT", 0));
        Assert.Throws<ArgumentNullException>(() => DotMatcher.Compute(null!, "ACGT"));
    }
}
