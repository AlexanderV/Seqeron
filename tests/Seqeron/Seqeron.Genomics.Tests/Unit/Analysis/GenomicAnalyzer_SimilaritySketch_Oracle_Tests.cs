// GENOMIC-SIMILARITY-001 (B09 A2-31 / A2-35) — sketch estimators, containment, Jaccard distance and ANI of GenomicAnalyzer.CompareSimilarity.
// Oracle: scripts/analysis/generate_similarity_sketch_oracle.py, 4290 rows — exact canonical k-mer sets (Python sets, == sourmash MinHash(scaled=1)),
// sourmash 4.9.4 FracMinHash (scaled 2..10: jaccard / contained_by / max_containment) and bottom-s MinHash (num, k >= 17 = the mash dist estimator),
// sourmash distance_utils containment_to_distance / jaccard_to_distance point estimates. Evidence: docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_SimilaritySketch_Oracle_Tests
{
    private sealed record Row(string A, string B, int K, string Kind, int Param, double?[] V);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracle() => OracleRows.Value;

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.similarity_sketch_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], int.Parse(f[2]), f[3], int.Parse(f[4]),
                f.Skip(5).Select(x => x == "." ? (double?)null : double.Parse(x, System.Globalization.CultureInfo.InvariantCulture)).ToArray()));
        }
        return rows;
    }

    private const double Tol = 1e-12;

    [Test]
    public void Exact_MatchesPythonSetsAndSourmashScaled1()
    {
        var rows = ReadOracle().Where(r => r.Kind == "exact").ToList();
        Assert.That(rows, Has.Count.EqualTo(990));
        foreach (var r in rows)
        {
            var e = GenomicAnalyzer.CompareSimilarity(r.A, r.B, r.K, canonical: true);
            Assert.That(e.Estimator, Is.EqualTo(SimilarityEstimator.Exact));
            string id = $"{r.A} {r.B} k={r.K}";
            Assert.That(e.Jaccard, Is.EqualTo(r.V[0]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAInB!.Value, Is.EqualTo(r.V[1]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentBInA!.Value, Is.EqualTo(r.V[2]!.Value).Within(Tol), id);
            Assert.That(e.MaxContainment!.Value, Is.EqualTo(r.V[3]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAniAInB!.Value, Is.EqualTo(r.V[4]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAniBInA!.Value, Is.EqualTo(r.V[5]!.Value).Within(Tol), id);
            Assert.That(e.JaccardAni, Is.EqualTo(r.V[6]!.Value).Within(Tol), id);
            Assert.That(e.JaccardDistance, Is.EqualTo(1.0 - r.V[0]!.Value).Within(Tol), id);
        }
    }

    [Test]
    public void FracMinHash_MatchesSourmashScaled()
    {
        var rows = ReadOracle().Where(r => r.Kind == "scaled").ToList();
        Assert.That(rows, Has.Count.EqualTo(1980));
        foreach (var r in rows)
        {
            var e = GenomicAnalyzer.CompareSimilarity(r.A, r.B, r.K, canonical: true, scaled: r.Param);
            string id = $"{r.A} {r.B} k={r.K} scaled={r.Param}";
            Assert.That(e.Estimator, Is.EqualTo(SimilarityEstimator.FracMinHash));
            Assert.That(e.Jaccard, Is.EqualTo(r.V[0]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAInB!.Value, Is.EqualTo(r.V[1]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentBInA!.Value, Is.EqualTo(r.V[2]!.Value).Within(Tol), id);
            Assert.That(e.MaxContainment!.Value, Is.EqualTo(r.V[3]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAniAInB!.Value, Is.EqualTo(r.V[4]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAniBInA!.Value, Is.EqualTo(r.V[5]!.Value).Within(Tol), id);
            Assert.That(e.JaccardAni, Is.EqualTo(r.V[6]!.Value).Within(Tol), id);
        }
    }

    [Test]
    public void BottomSketch_MatchesSourmashNumSketch()
    {
        var rows = ReadOracle().Where(r => r.Kind == "num").ToList();
        Assert.That(rows, Has.Count.EqualTo(1320));
        int informative = 0;
        foreach (var r in rows)
        {
            var e = GenomicAnalyzer.CompareSimilarity(r.A, r.B, r.K, canonical: true, sketchSize: r.Param);
            string id = $"{r.A} {r.B} k={r.K} s={r.Param}";
            Assert.That(e.Estimator, Is.EqualTo(SimilarityEstimator.MashBottomSketch));
            Assert.That(e.Jaccard, Is.EqualTo(r.V[0]!.Value).Within(Tol), id);
            Assert.That(e.ContainmentAInB, Is.Null);
            Assert.That(GenomicAnalyzer.CalculateSimilarity(r.A, r.B, r.K, true, r.Param), Is.EqualTo(r.V[0]!.Value * 100.0).Within(1e-9), id);
            if (r.V[0]!.Value > 0 && r.V[0]!.Value < 1) informative++;
        }

        Assert.That(informative, Is.GreaterThan(300), "the corpus must contain partial-overlap estimates");
    }

    [Test]
    public void CalculateSimilarity_Sketch_IsThePercentOfTheEstimate()
    {
        var a = new DnaSequence("ACGTACGTTAGCTAGCTAGGATCGATCGATTTAGCGCGATATCGGCTAGCTAGCTAGGGATCTAGCTAGCTTAGCGA");
        var b = new DnaSequence("ACGTACGTTAGCTAGCTAGGATCGATCGATTTAGCGCGATATCGGCTAGGGGATCTAGCTAGCTTAGCGATTTTTT");
        double pct = GenomicAnalyzer.CalculateSimilarity(a, b, 21, true, 50);
        Assert.That(pct, Is.EqualTo(GenomicAnalyzer.CompareSimilarity(a, b, 21, true, sketchSize: 50).Jaccard * 100.0).Within(1e-9));
        // a sketch at least as large as the union is the exact index
        Assert.That(GenomicAnalyzer.CalculateSimilarity(a, b, 21, true, 100000), Is.EqualTo(GenomicAnalyzer.CalculateSimilarity(a, b, 21, true)).Within(1e-9));
    }

    [Test]
    public void JaccardDistance_IsOneMinusTheIndex_AndStringOverloadsSkipAmbiguousWindows()
    {
        var rnd = new Random(5);
        for (int iter = 0; iter < 200; iter++)
        {
            string x = new string(Enumerable.Range(0, rnd.Next(1, 60)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string y = new string(Enumerable.Range(0, rnd.Next(1, 60)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            int k = rnd.Next(1, 7);
            bool canonical = rnd.Next(2) == 0;
            var dx = new DnaSequence(x);
            var dy = new DnaSequence(y);
            Assert.That(GenomicAnalyzer.CalculateJaccardDistance(dx, dy, k, canonical), Is.EqualTo(1.0 - GenomicAnalyzer.CalculateSimilarity(dx, dy, k, canonical) / 100.0).Within(1e-12));
            Assert.That(GenomicAnalyzer.CalculateSimilarity(x.ToLowerInvariant(), y, k, canonical), Is.EqualTo(GenomicAnalyzer.CalculateSimilarity(dx, dy, k, canonical)).Within(1e-12));
            Assert.That(GenomicAnalyzer.CalculateJaccardDistance(x, y, k, canonical), Is.EqualTo(GenomicAnalyzer.CalculateJaccardDistance(dx, dy, k, canonical)).Within(1e-12));
        }

        // N windows are skipped: ACNGT has no 3-mer window without N except none; AC-GT-AC-GT: k=2 windows AC, GT only
        Assert.That(GenomicAnalyzer.CalculateSimilarity("ACNGT", "ACGT", 2), Is.EqualTo(100.0 * 2 / 3).Within(1e-12)); // {AC,GT} vs {AC,CG,GT}
    }

    [Test]
    public void CompareSimilarity_Boundaries_AndValidation()
    {
        var same = GenomicAnalyzer.CompareSimilarity("ACGTACGGTTCA", "ACGTACGGTTCA", 4, canonical: true);
        Assert.That((same.Jaccard, same.JaccardDistance, same.MashDistance, same.JaccardAni, same.ContainmentAniAInB), Is.EqualTo((1.0, 0.0, 0.0, 1.0, 1.0)));
        var none = GenomicAnalyzer.CompareSimilarity("AAAAAAAA", "CCCCCCCC", 4);
        Assert.That((none.Jaccard, none.MashDistance, none.JaccardAni, none.ContainmentAniAInB), Is.EqualTo((0.0, 1.0, 0.0, 0.0)));
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.CompareSimilarity("ACGT", "ACGT", 2, false, 10, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.CompareSimilarity("ACGT", "ACGT", 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.CompareSimilarity("ACGT", "ACGT", 4, false, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.CompareSimilarity("ACGT", "ACGT", 4, false, null, 0));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.CompareSimilarity((string)null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.CalculateJaccardDistance((DnaSequence)null!, new DnaSequence("ACGT")));
    }
}
