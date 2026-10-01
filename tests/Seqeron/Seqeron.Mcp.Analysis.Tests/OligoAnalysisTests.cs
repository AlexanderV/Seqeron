using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// oligo_analysis — delegation to the RSAT overload
/// <c>MotifFinder.DiscoverMotifs(seq, k, minCount, OligoBackgroundModel, OligoStrandMode, countOverlapping)</c>.
/// Expected values: RSAT oligo-analysis outputs locked in MotifFinder_OligoAnalysis_Tests (t2.fa).
/// </summary>
[TestFixture]
public class OligoAnalysisTests
{
    private const string T2 = "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC";
    private const double Rel = 1e-10;

    [Test]
    public void OligoAnalysis_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.OligoAnalysis(T2, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("", 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis("ACGU", 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.OligoAnalysis(T2, 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis(T2, 4, background: "uniform"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis(T2, 4, background: "bernoulli"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis(T2, 4, background: "markov_table"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.OligoAnalysis(T2, 4, background: "markov", markovOrder: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.OligoAnalysis(T2, 4, background: "markov", markovOrder: 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.OligoAnalysis(T2, 4, strands: "minus"));
    }

    [Test]
    public void OligoAnalysis_Equiprobable_SingleStrand_EqualsRsat()
    {
        // oligo-analysis -i t2.fa -l 4 -1str -bg equi -lth occ 2 -return occ,proba
        var r = AnalysisTools.OligoAnalysis(T2, 4, 2, "equiprobable");
        var atgc = r.Motifs.Single(m => m.Sequence == "ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(r.TotalOccurrences, Is.EqualTo(60));
            Assert.That(r.TestedPatterns, Is.EqualTo(10));
            Assert.That(r.PossibleOligos, Is.EqualTo(256));
            Assert.That(r.Strands, Is.EqualTo("single"));
            Assert.That(atgc.Count, Is.EqualTo(6));
            Assert.That(atgc.ReverseComplement, Is.Null);
            Assert.That(atgc.ExpectedFrequency, Is.EqualTo(0.00390625).Within(Rel).Percent);
            Assert.That(atgc.OccurrenceProbability, Is.EqualTo(1.4844942103072834e-07).Within(Rel).Percent);
            Assert.That(atgc.OccurrenceEValue, Is.EqualTo(1.4844942103072835e-06).Within(Rel).Percent);
            Assert.That(atgc.OccurrenceSignificance, Is.EqualTo(5.8284214918614703).Within(1e-11));
            Assert.That(atgc.Ratio, Is.EqualTo(6 / (60 / 256.0)).Within(Rel).Percent);
        });
    }

    [Test]
    public void OligoAnalysis_BothStrandsInput_AndMarkov1_EqualRsat()
    {
        // -l 4 -2str (input Bernoulli) -lth occ 2: atgc|gcat occ 9, occ_P 1.0760647330191343e-09, occ_sig 7.7920703428970466
        var both = AnalysisTools.OligoAnalysis(T2, 4, 2, "input", strands: "both");
        var pair = both.Motifs.Single(m => m.Sequence == "ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(both.PossibleOligos, Is.EqualTo(136));
            Assert.That(both.TestedPatterns, Is.EqualTo(15));
            Assert.That(pair.ReverseComplement, Is.EqualTo("GCAT"));
            Assert.That(pair.Count, Is.EqualTo(9));
            Assert.That(pair.OccurrenceProbability, Is.EqualTo(1.0760647330191343e-09).Within(Rel).Percent);
            Assert.That(pair.OccurrenceSignificance, Is.EqualTo(7.7920703428970466).Within(1e-11));
        });

        // -l 4 -1str -markov 1 (no threshold): atgc exp_freq 0.032915191520534237
        var m1 = AnalysisTools.OligoAnalysis(T2, 4, 1, "markov", markovOrder: 1);
        Assert.That(m1.TestedPatterns, Is.EqualTo(40));
        Assert.That(m1.Motifs.Single(m => m.Sequence == "ATGC").ExpectedFrequency,
            Is.EqualTo(0.032915191520534237).Within(Rel).Percent);
    }

    [Test]
    public void OligoAnalysis_MarkovTable_EqualsLibrary()
    {
        var table = new Dictionary<string, double>();
        int i = 1;
        foreach (char a in "ACGT")
            foreach (char b in "ACGT")
                table[$"{a}{b}"] = i++;

        var r = AnalysisTools.OligoAnalysis(T2, 4, 5, "markov_table", oligoFrequencies: table);
        var lib = global::Seqeron.Genomics.Analysis.MotifFinder.DiscoverMotifs(
            new global::Seqeron.Genomics.Core.DnaSequence(T2), 4, 5,
            global::Seqeron.Genomics.Analysis.OligoBackgroundModel.MarkovFromOligoFrequencies(table));
        // RSAT -bgfile (psi 0.01) -l 4 -1str -lth occ 5, unrounded: atgc occ_sig 8.436534013635193
        Assert.That(r.Motifs.Single(m => m.Sequence == "ATGC").OccurrenceSignificance, Is.EqualTo(8.436534013635193).Within(1e-11));
        Assert.That(r.Motifs.Select(m => (m.Sequence, m.Count, m.ExpectedFrequency, m.OccurrenceProbability, m.OccurrenceSignificance)),
            Is.EqualTo(lib.Motifs.Select(m => (m.Sequence, m.Count, m.ExpectedFrequency, m.OccurrenceProbability, m.OccurrenceSignificance))));
        Assert.That(r.TestedPatterns, Is.EqualTo(lib.TestedPatterns));
    }

    [Test]
    public void OligoAnalysis_PossibleOligosBeyondDoubleRange_MapsToNull()
    {
        // k = 600: 4^600 exceeds the double range (library returns +Infinity, not representable in JSON).
        var rng = new Random(5);
        string x = new(Enumerable.Range(0, 600).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        var r = AnalysisTools.OligoAnalysis(x + x, 600, 1, "equiprobable");
        Assert.That(r.PossibleOligos, Is.Null);
        Assert.That(r.Motifs.Single(m => m.Sequence == x).Count, Is.EqualTo(2));
    }
}
