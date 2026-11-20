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

    // ------------------------------------------------------------------ A4-9: the edge cases the XML docs state

    private sealed record EdgeRow(string Case, string Matrix, int Window, int Threshold, string Seq1, string Seq2, string Segments);

    private static readonly Lazy<List<EdgeRow>> EdgeOracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.dotmatcher_edge_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<EdgeRow>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new EdgeRow(f[0], f[1], int.Parse(f[2]), int.Parse(f[3]), f[4], f[5], f[6]));
        }

        return rows;
    });

    [Test]
    public void Compute_EqualsDotmatcher_OnTheEdgeCaseRows()
    {
        // 720 rows of the installed 6.6.0 binary: lower / mixed case, IUPAC ambiguity codes (and the letters J O U that are not in a matrix: score 0), '*', gap characters - . ~ (removed by the
        // reader: positions count the cleaned sequences), the default matrix by the type of the first sequence (also when the second is of the other type) and a window as long as / up to 8 residues
        // longer than a sequence (the program scores the zero terminator beyond the end as 0)
        var rows = EdgeOracle.Value;
        Assert.That(rows, Has.Count.EqualTo(720));
        var failures = new List<string>();
        var perCase = new Dictionary<string, int>();
        int withSegments = 0, refused = 0, overhanging = 0;
        foreach (var r in rows)
        {
            perCase[r.Case] = perCase.GetValueOrDefault(r.Case) + 1;
            string actual;
            try
            {
                var matrix = r.Matrix == "-" ? null : SubstitutionMatrix.Get(r.Matrix);
                actual = Format(DotMatcher.Compute(r.Seq1, r.Seq2, r.Window, r.Threshold, matrix));
            }
            catch (ArgumentException ex) when (ex is not ArgumentOutOfRangeException)
            {
                actual = "ERR:notnucleic";
            }

            if (actual != r.Segments)
                failures.Add($"{r.Case} {r.Matrix} w={r.Window} t={r.Threshold}\n{r.Seq1}\n{r.Seq2}\n expected {r.Segments}\n actual   {actual}");
            if (r.Segments == "ERR:notnucleic") refused++;
            else if (r.Segments != "-") withSegments++;
            if (r.Case == "overhang" && r.Window >= Math.Min(Clean(r.Seq1).Length, Clean(r.Seq2).Length)) overhanging++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}: " + string.Join("\n", failures.Take(3)));
        foreach (var c in new[] { "lower", "ambig", "star", "gap", "default", "overhang" })
            Assert.That(perCase[c], Is.GreaterThanOrEqualTo(100), c);
        Assert.That(withSegments, Is.GreaterThan(350));
        Assert.That(refused, Is.GreaterThan(30));
        Assert.That(overhanging, Is.GreaterThanOrEqualTo(100));
    }

    private static string Clean(string s) => EmbossSequenceReader.Clean(s);

    [Test]
    public void WindowEqualToASequence_ScoresTheDiagonalsThatFit_AndALongerWindowReadsZerosBeyondTheEnd()
    {
        // measured: a = 12 nt, b = 10 nt, threshold 5, EDNAFULL: w = 10 -> one segment (1,1,10,10); w = 11 -> (1,1,11,11), the window counted whole with the terminator scoring 0; w = 12 -> none
        string a = "ACGTACGTACGT", b = "ACGTACGTAC";
        Assert.That(Format(DotMatcher.Compute(a, b, 10, 5)), Is.EqualTo("1:1:10:10"));
        Assert.That(Format(DotMatcher.Compute(a, b, 11, 5)), Is.EqualTo("1:1:11:11"));
        Assert.That(DotMatcher.Compute(a, b, 12, 5), Is.Empty);
        // a nucleic first sequence needs a nucleic second one (ACD)
        Assert.Throws<ArgumentException>(() => DotMatcher.Compute("ACGTACGTACGT", "MKVLAAGIVGLL", 3, 5));
        Assert.DoesNotThrow(() => DotMatcher.Compute("MKVLAAGIVGLL", "ACGTACGTACGT", 3, 5));
    }
}
