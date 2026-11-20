// GENOMIC-SIMILARITY-001 (B09 F32-F34) — GenomicAnalyzer.CalculateSimilarity is a facade over the canonical
// KmerAnalyzer.JaccardSimilarity (DUP_MAP §16) x 100. Evidence: docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-SIMILARITY-001.md (section 8)
// Oracle: scripts/analysis/generate_similarity_oracle.py (TestData/Analysis/similarity_oracle.tsv, 2086 rows):
//   literal   = Python set == scipy.spatial.distance.jaccard on presence vectors (scipy 1.18.1);
//   canonical = Python set over min(w, revcomp(w)) == sourmash 4.9.4 MinHash(n=0, scaled=1).jaccard (exact).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_CalculateSimilarity_Oracle_Tests
{
    private sealed record Row(string Seq1, string Seq2, int K, double Literal, double Canonical);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracle() => OracleRows.Value;

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.similarity_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0] == "-" ? "" : f[0], f[1] == "-" ? "" : f[1], int.Parse(f[2]),
                double.Parse(f[3], System.Globalization.CultureInfo.InvariantCulture),
                double.Parse(f[4], System.Globalization.CultureInfo.InvariantCulture)));
        }
        return rows;
    }

    [Test]
    public void CalculateSimilarity_Literal_MatchesPythonSetAndScipyJaccard()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(2086));
        foreach (var r in rows)
        {
            double actual = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(r.Seq1), new DnaSequence(r.Seq2), r.K);
            Assert.That(actual, Is.EqualTo(r.Literal).Within(1e-9), $"{r.Seq1} / {r.Seq2} k={r.K}");
        }
    }

    [Test]
    public void CalculateSimilarity_Canonical_MatchesPythonSetAndSourmashExactJaccard()
    {
        int differing = 0;
        foreach (var r in ReadOracle())
        {
            double actual = GenomicAnalyzer.CalculateSimilarity(
                new DnaSequence(r.Seq1), new DnaSequence(r.Seq2), r.K, canonical: true);
            Assert.That(actual, Is.EqualTo(r.Canonical).Within(1e-9), $"{r.Seq1} / {r.Seq2} k={r.K}");
            if (Math.Abs(r.Canonical - r.Literal) > 1e-9) differing++;
        }
        Assert.That(differing, Is.GreaterThan(500), "the oracle must discriminate canonical from literal k-mers");
    }

    [Test]
    public void CalculateSimilarity_EqualsCanonicalKmerAnalyzerJaccardTimes100()
    {
        // DUP_MAP §16: the facade is KmerAnalyzer.JaccardSimilarity * 100 (literal and canonical modes).
        foreach (var r in ReadOracle().Take(600))
        {
            var a = new DnaSequence(r.Seq1);
            var b = new DnaSequence(r.Seq2);
            Assert.That(GenomicAnalyzer.CalculateSimilarity(a, b, r.K),
                Is.EqualTo(KmerAnalyzer.JaccardSimilarity(a.Sequence, b.Sequence, r.K) * 100.0));
            Assert.That(GenomicAnalyzer.CalculateSimilarity(a, b, r.K, canonical: true),
                Is.EqualTo(KmerAnalyzer.JaccardSimilarity(a.Sequence, b.Sequence, r.K,
                    new KmerCountingOptions(Canonical: true)) * 100.0));
        }
    }

    // AACCCGGG vs its reverse complement CCCGGGTT: literal 3-mer sets share 4 of 8 (50 %); canonical sets are equal (100 %).
    [Test]
    public void CalculateSimilarity_ReverseComplement_LiteralLowCanonical100()
    {
        var s = new DnaSequence("AACCCGGG");
        var rc = s.ReverseComplement();
        Assert.That(GenomicAnalyzer.CalculateSimilarity(s, rc, 3, canonical: true), Is.EqualTo(100.0).Within(1e-12));
        Assert.That(GenomicAnalyzer.CalculateSimilarity(s, rc, 3), Is.EqualTo(50.0).Within(1e-12));
    }

    [Test]
    public void CalculateSimilarity_CaseInsensitive_DnaSequenceFoldsCase()
    {
        Assert.That(GenomicAnalyzer.CalculateSimilarity(new DnaSequence("acgtacgt"), new DnaSequence("ACGTACGA"), 3),
            Is.EqualTo(80.0).Within(1e-10));
    }

    [Test]
    public void CalculateSimilarity_KLargerThanBothSequences_ReturnsZeroInBothModes()
    {
        var a = new DnaSequence("ACGT");
        var b = new DnaSequence("ACGT");
        Assert.That(GenomicAnalyzer.CalculateSimilarity(a, b, 9), Is.EqualTo(0.0));
        Assert.That(GenomicAnalyzer.CalculateSimilarity(a, b, int.MaxValue, canonical: true), Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateSimilarity_CanonicalOverload_ValidatesArguments()
    {
        var s = new DnaSequence("ACGT");
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.CalculateSimilarity(null!, s, 3, true));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.CalculateSimilarity(s, null!, 3, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.CalculateSimilarity(s, s, 0, true));
    }

    [Test]
    public void CalculateSimilarity_LargeSequences_Linear()
    {
        var rnd = new Random(7);
        string Make() => new string(Enumerable.Range(0, 300_000).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        var a = new DnaSequence(Make());
        var sw = System.Diagnostics.Stopwatch.StartNew();
        double self = GenomicAnalyzer.CalculateSimilarity(a, a, 21, canonical: true);
        sw.Stop();
        Assert.That(self, Is.EqualTo(100.0));
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(30));
    }
}
