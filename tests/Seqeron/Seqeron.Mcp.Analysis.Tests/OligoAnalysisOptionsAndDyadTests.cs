using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// oligo_analysis RSAT options (zscore, expectedFrequencyPseudo, degenerate, lexicon, calibration, extraSequences) →
/// <c>MotifFinder.AnalyzeOligos</c>; dyad_analysis → <c>MotifFinder.AnalyzeDyads</c>.
/// Expected values: RSAT oligo-analysis / dyad-analysis runs locked in MotifFinder_OligoAnalysisOptions_Tests and
/// MotifFinder_DyadAnalysis_Tests (t2.fa).
/// </summary>
[TestFixture]
public class OligoAnalysisOptionsAndDyadTests
{
    private const string T2 = "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC";
    private const double Rel = 1e-10;

    [Test]
    public void OligoAnalysis_ZScore_EqualsRsat()
    {
        // oligo-analysis -i t2.fa -l 4 -1str -bg equi -return occ,proba,zscore -lth occ 2
        var r = AnalysisTools.OligoAnalysis(T2, 4, 2, "equiprobable", zscore: true);
        var tgct = r.Motifs.Single(m => m.Sequence == "TGCT");
        Assert.Multiple(() =>
        {
            Assert.That(r.TestedPatterns, Is.EqualTo(10));
            Assert.That(r.SequenceCount, Is.EqualTo(1));
            Assert.That(tgct.OverlapCoefficient, Is.EqualTo(1.015625));
            Assert.That(tgct.ExpectedVariance, Is.EqualTo(0.23345947265625).Within(Rel).Percent);
            Assert.That(tgct.ZScore, Is.EqualTo(3.6542034172140836).Within(Rel).Percent);
            Assert.That(tgct.FittedDistribution, Is.EqualTo("Binomial"));
            Assert.That(r.Motifs.Single(m => m.Sequence == "ATGC").OccurrenceProbability,
                Is.EqualTo(1.4844942103072834e-07).Within(Rel).Percent);
        });
    }

    [Test]
    public void OligoAnalysis_PseudoLexiconDegenerateCalibration_EqualRsat()
    {
        // -l 4 -1str -markov 1 -pseudo 0.1 -lth occ 2: atgc exp_freq 0.030014297368480814, occ_P 0.0091641466020121916
        var pseudo = AnalysisTools.OligoAnalysis(T2, 4, 2, "markov", markovOrder: 1, expectedFrequencyPseudo: 0.1);
        // -l 4 -1str -lexicon -lth occ 2: atgc exp_freq 0.025, segments a|tgc, occ_P 0.003853224500294343
        var lex = AnalysisTools.OligoAnalysis(T2, 4, 2, "lexicon");
        // Fixed-oracle -onedeg (RSAT 1.169 returns nothing, see library tests): -l 3 -lth occ 6 → 71 tested, NPO 528
        var deg = AnalysisTools.OligoAnalysis(T2, 3, 6, degenerate: "onedeg");
        // -calibN cal3.tab -l 3 -1str: tgc Poisson(1.5) exact 0.0044559807752478468
        const string cal = "atg\t0.9\t1.2\t1.44\ntgc\t1.5\t1.1\t1.21\ncat\t2.0\t1.5\t2.25\ngca\t0.8\t0.8\t0.64\n";
        var calN = AnalysisTools.OligoAnalysis(T2, 3, 1, calibrationTable: cal);
        Assert.Multiple(() =>
        {
            Assert.That(pseudo.Motifs.Single(m => m.Sequence == "ATGC").ExpectedFrequency, Is.EqualTo(0.030014297368480814).Within(Rel).Percent);
            Assert.That(pseudo.Motifs.Single(m => m.Sequence == "ATGC").OccurrenceProbability, Is.EqualTo(0.0091641466020121916).Within(Rel).Percent);
            Assert.That(lex.Motifs.Single(m => m.Sequence == "ATGC").ExpectedFrequency, Is.EqualTo(0.025).Within(Rel).Percent);
            Assert.That(lex.Motifs.Single(m => m.Sequence == "ATGC").LexiconSegmentation, Is.EqualTo("A|TGC"));
            Assert.That(lex.Motifs.Single(m => m.Sequence == "ATGC").OccurrenceProbability, Is.EqualTo(0.003853224500294343).Within(Rel).Percent);
            Assert.That(deg.PossibleOligos, Is.EqualTo(528));
            Assert.That(deg.TestedPatterns, Is.EqualTo(71));
            Assert.That(deg.Degenerate, Is.EqualTo("onedeg"));
            Assert.That(deg.Motifs.Single(m => m.Sequence == "AYG").OccurrenceProbability, Is.EqualTo(0.00054574388255023234).Within(Rel).Percent);
            Assert.That(calN.TestedPatterns, Is.EqualTo(4));
            Assert.That(calN.Motifs.Single(m => m.Sequence == "TGC").FittedDistribution, Is.EqualTo("Poisson"));
            Assert.That(calN.Motifs.Single(m => m.Sequence == "TGC").OccurrenceProbability, Is.EqualTo(0.0044559807752478468).Within(Rel).Percent);
            Assert.That(calN.Motifs.Single(m => m.Sequence == "ATG").FittedDistribution, Is.EqualTo("NegativeBinomial"));
            Assert.That(calN.Motifs.Single(m => m.Sequence == "ATG").ExpectedVariance, Is.EqualTo(1.44));
        });
    }

    [Test]
    public void OligoAnalysis_ExtraSequences_AndValidation()
    {
        var r = AnalysisTools.OligoAnalysis("ATGCATGCATGCAAATTTGGGCC", 3, 1, strands: "both", countOverlapping: false,
            extraSequences: new[] { "CATGCTTAGCGGATCCATGCATGCTTTAAACG" },
            calibrationTable: "atg\t0.9\t1.2\t1.44\ngca\t0.8\t0.8\t0.64\n", calibrationMode: "sequence");
        var atg = r.Motifs.Single(m => m.Sequence == "ATG");
        Assert.Multiple(() =>
        {
            // RSAT (exact negbin) -calib1, 2 sequences, -2str -noov: atg|cat occ 6, overlaps 5, occ_P 0.03602153062820438
            Assert.That(r.SequenceCount, Is.EqualTo(2));
            Assert.That(atg.Count, Is.EqualTo(6));
            Assert.That(atg.Overlaps, Is.EqualTo(5));
            Assert.That(atg.SequenceIndices!.Distinct().Order(), Is.EqualTo(new[] { 0, 1 }));
            Assert.That(atg.OccurrenceProbability, Is.EqualTo(0.03602153062820438).Within(Rel).Percent);
            Assert.That(() => AnalysisTools.OligoAnalysis(T2, 3, degenerate: "two"), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.OligoAnalysis(T2, 3, expectedFrequencyPseudo: 2), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => AnalysisTools.OligoAnalysis(T2, 3, calibrationTable: "atg 1"), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.OligoAnalysis(T2, 3, calibrationTable: "atg\t1\t1\t1\n", calibrationMode: "x"), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.OligoAnalysis(T2, 3, extraSequences: new[] { "ACGU" }), Throws.ArgumentException);
        });
    }

    [Test]
    public void DyadAnalysis_EqualsRsat()
    {
        // dyad-analysis -i t2.fa -l 3 -sp 0-2 -1str -lth occ 2 (-noov): 13 tested, atgn{0}cat occ_P 0.079911679827561213
        var r = AnalysisTools.DyadAnalysis(new[] { T2 }, 3, 0, 2, strands: "single", minCount: 2);
        var d = r.Dyads.Single(x => x.Pattern == "ATGn{0}CAT");
        // RSAT defaults (-2str -noov): atgn{0}cat exp_occ 0.40311744154797091, per spacing (58,50,8)
        var both = AnalysisTools.DyadAnalysis(new[] { T2 }, 3, 0, 2, minCount: 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.TestedPatterns, Is.EqualTo(13));
            Assert.That(r.PossibleDyads, Is.EqualTo(12288));
            Assert.That(d.Occurrences, Is.EqualTo(2));
            Assert.That(d.Overlaps, Is.EqualTo(1));
            Assert.That(d.OccurrenceProbability, Is.EqualTo(0.079911679827561213).Within(Rel).Percent);
            Assert.That(d.ZScore, Is.EqualTo(2.41).Within(0.005));
            Assert.That(both.Strands, Is.EqualTo("both"));
            Assert.That(both.Spacings[0], Is.EqualTo(new DyadSpacingItem(0, 58, 50, 8)));
            Assert.That(both.Dyads.Single(x => x.Pattern == "ATGn{0}CAT").ExpectedOccurrences, Is.EqualTo(0.40311744154797091).Within(Rel).Percent);
            Assert.That(both.Dyads.Single(x => x.Pattern == "ATGn{0}CAT").IsReversePalindrome, Is.True);
        });
    }

    [Test]
    public void DyadAnalysis_TypesBackgroundsAndValidation()
    {
        // dyad-analysis -i t2.fa -l 2 -sp 0-6 -type dr -1str -ovlp -lth occ 2: 6 tested, atn{2}at occ_P 0.015985324783370104
        var dr = AnalysisTools.DyadAnalysis(new[] { T2 }, 2, 0, 6, "dr", "single", true, 2);
        var table = AnalysisTools.DyadAnalysis(new[] { T2 }, 3, 0, 0, strands: "single", minCount: 2,
            background: "dyad_table", dyadFrequencies: new Dictionary<string, double> { ["atgn{0}cat"] = 0.01 });
        Assert.Multiple(() =>
        {
            Assert.That(dr.TestedPatterns, Is.EqualTo(6));
            Assert.That(dr.Dyads.Single(x => x.Pattern == "ATn{2}AT").OccurrenceProbability, Is.EqualTo(0.015985324783370104).Within(Rel).Percent);
            Assert.That(table.TestedPatterns, Is.EqualTo(1));
            Assert.That(table.Dyads.Single(x => x.Pattern == "ATGn{0}CAT").ExpectedFrequency, Is.EqualTo(0.01));
            Assert.That(table.Dyads.Where(x => x.Pattern != "ATGn{0}CAT").All(x => x.OccurrenceProbability is null), Is.True);
            Assert.That(() => AnalysisTools.DyadAnalysis(Array.Empty<string>()), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.DyadAnalysis(new[] { T2 }, monadLength: 0), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => AnalysisTools.DyadAnalysis(new[] { T2 }, minSpacing: 3, maxSpacing: 2), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => AnalysisTools.DyadAnalysis(new[] { T2 }, dyadType: "xx"), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.DyadAnalysis(new[] { T2 }, background: "monad_table"), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.DyadAnalysis(new[] { T2 }, strands: "minus"), Throws.ArgumentException);
            Assert.That(() => AnalysisTools.DyadAnalysis(new[] { "ACGU" }), Throws.ArgumentException);
        });
    }
}
