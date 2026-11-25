using System.Globalization;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// ANNOT-GENE-001 (B11 A2): genome-trained Prodigal RBS / start model (Hyatt et al. 2010, BMC Bioinformatics 11:119).
/// Oracle: the Prodigal V2.6.3 binary run on E. coli K-12 MG1655 (U00096.2) nt 1..150000 — its binary training file
/// (<c>-t</c>: rbs_wt[28], type_wt[3]) and its <c>-f gff</c> gene set (see docs/Evidence/ANNOT-GENE-001-Evidence.md).
/// </summary>
[TestFixture]
[Category("Annotation")]
[Category("GenePrediction")]
public class GenomeAnnotator_ProdigalStart_Tests
{
    private static readonly Lazy<string> Genome = new(() => Resource("ecoli_U00096.2_1-150000.txt").Trim());
    private static readonly Lazy<RbsModel> Model = new(() => GenomeAnnotator.TrainRbsModel(Genome.Value));

    private static string Resource(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Prodigal." + file;
        using var stream = typeof(GenomeAnnotator_ProdigalStart_Tests).Assembly.GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        return reader.ReadToEnd();
    }

    private static IEnumerable<string[]> OracleRows(string file) =>
        Resource(file).Split('\n').Where(l => l.Length > 0 && l[0] != '#').Select(l => l.TrimEnd('\r').Split('\t'));

    private static List<(int begin, int end, char strand)> OracleGenes() =>
        OracleRows("prodigal_slice_genes.tsv").Select(f => (int.Parse(f[0]), int.Parse(f[1]), f[2][0])).ToList();

    private static (int stop, char strand) StopKey(int begin, int end, char strand) => (strand == '+' ? end : begin, strand);

    #region Training = Prodigal training file

    [Test]
    public void TrainRbsModel_Weights_EqualProdigalTrainingFile()
    {
        var rows = OracleRows("prodigal_slice_model.tsv").ToList();
        var rbs = rows.Where(r => r[0] == "rbs").Select(r => double.Parse(r[2], CultureInfo.InvariantCulture)).ToArray();
        var type = rows.Where(r => r[0] == "type").Select(r => double.Parse(r[2], CultureInfo.InvariantCulture)).ToArray();

        Assert.That(Model.Value.BinWeights, Has.Count.EqualTo(28));
        for (int i = 0; i < 28; i++)
            Assert.That(Model.Value.BinWeights[i], Is.EqualTo(rbs[i]).Within(1e-9), $"rbs_wt[{i}]");
        for (int i = 0; i < 3; i++)
            Assert.That(Model.Value.StartTypeWeights[i], Is.EqualTo(type[i]).Within(1e-9), $"type_wt[{i}]");
        Assert.That(Model.Value.UsesShineDalgarno, Is.True);
    }

    [Test]
    public void TrainRbsModel_BinFrequencies_AreDistributions()
    {
        Assert.That(Model.Value.BinFrequencyInGenes.Sum(), Is.EqualTo(1.0).Within(1e-9));
        Assert.That(Model.Value.BinFrequencyBackground.Sum(), Is.EqualTo(1.0).Within(1e-9));
        Assert.That(Model.Value.TypeFrequencyInGenes.Sum(), Is.EqualTo(1.0).Within(1e-9));
        // Strong SD (AGGAGG, 5-10 nt spacer = bin 27) is enriched in real starts over the background.
        Assert.That(Model.Value.BinWeights[27], Is.GreaterThan(3.0));
    }

    [Test]
    public void TrainRbsModel_ShortSequence_Throws() =>
        Assert.Throws<ArgumentException>(() => GenomeAnnotator.TrainRbsModel(new string('A', 19999)));

    #endregion

    #region Prodigal gene model (dynamic programming) = Prodigal gff

    [Test]
    public void PredictGenes_ProdigalDynamicProgramming_GeneSetIdenticalToProdigal()
    {
        var expected = OracleGenes();
        var genes = GenomeAnnotator.PredictGenes(Genome.Value, 30, "g", StartSelection.ProdigalDynamicProgramming, Model.Value)
            .Select(g => (begin: g.Start + 1, end: g.End, strand: g.Strand)).ToList();

        Assert.That(expected, Has.Count.EqualTo(129));
        Assert.That(genes, Is.EqualTo(expected));
    }

    #endregion

    #region Start selection among nested starts

    [Test]
    public void PredictGenes_ProdigalScore_StartsAgreeWithProdigal_BetterThanLongestOrf()
    {
        var oracle = OracleGenes().ToDictionary(g => StopKey(g.begin, g.end, g.strand));

        int Agree(StartSelection mode, out int shared)
        {
            var mine = GenomeAnnotator.PredictGenes(Genome.Value, 30, "g", mode, Model.Value)
                .Select(g => (begin: g.Start + 1, end: g.End, strand: g.Strand))
                .ToDictionary(g => StopKey(g.begin, g.end, g.strand));
            var common = oracle.Keys.Where(mine.ContainsKey).ToList();
            shared = common.Count;
            return common.Count(k => mine[k] == oracle[k]);
        }

        int scoreAgree = Agree(StartSelection.ProdigalScore, out int shared);
        int longestAgree = Agree(StartSelection.LongestOrf, out _);

        Assert.That(shared, Is.GreaterThanOrEqualTo(125));
        Assert.That(scoreAgree, Is.EqualTo(126), "Prodigal-score start selection: 126 of 128 stop-shared genes");
        Assert.That(longestAgree, Is.EqualTo(81), "longest-ORF rule: 81 of 128");
    }

    [Test]
    public void PredictGenes_DefaultStartSelection_IsLongestOrf()
    {
        string seq = "AAAA" + "ATG" + "GCT" + "ATG" + string.Concat(Enumerable.Repeat("GCT", 40)) + "TAA" + "CCCC";
        var def = GenomeAnnotator.PredictGenes(seq, 30).ToList();
        var explicitLongest = GenomeAnnotator.PredictGenes(seq, 30, "gene", StartSelection.LongestOrf).ToList();

        Assert.That(def, Has.Count.EqualTo(1));
        Assert.That(def[0].Start, Is.EqualTo(4));
        Assert.That((explicitLongest[0].Start, explicitLongest[0].End), Is.EqualTo((def[0].Start, def[0].End)));
    }

    [Test]
    public void PredictGenes_ProdigalScore_ShortGenome_Throws() =>
        Assert.Throws<ArgumentException>(() => GenomeAnnotator.PredictGenes(
            "ATG" + string.Concat(Enumerable.Repeat("GCT", 40)) + "TAA", 30, "g", StartSelection.ProdigalScore).ToList());

    #endregion

    #region RBS bins (formula-derived, cross-checked with an independent Python port)

    private static string Spaced(string motif, int spacer) => new string('C', 25) + motif + new string('C', spacer) + "ATG";

    [TestCase("AGGAGG", 7, 27, 0)]      // 6-mer score 14, spacer class 0 (5-10 nt)
    [TestCase("AGGAGG", 4, 26, 0)]      // 6-mer, rdis < 5 -> class 1
    [TestCase("AGGAGG", 11, 25, 0)]     // 6-mer, rdis 11-12 -> class 2
    [TestCase("AGGAGG", 13, 10, 0)]     // rdis >= 13 -> class 3
    [TestCase("GGAGG", 7, 24, 0)]       // 5-mer score 12, class 0
    [TestCase("TGGAGT", 7, 16, 0)]      // GGAG: 4-mer score 9, class 0 (Prodigal bins 11/12/16)
    [TestCase("AGTAGG", 7, 13, 19)]     // exact AGG (bin 13) vs 6-mer with one internal mismatch (bin 19)
    [TestCase("TGGAT", 5, 13, 0)]       // GGA: 3-mer score 6, class 0
    public void ClassifyRbsStart_Bins_MatchFormulaTable(string motif, int spacer, int exactBin, int mismatchBin)
    {
        string seq = Spaced(motif, spacer);
        var c = GenomeAnnotator.ClassifyRbsStart(seq, seq.Length - 3);

        Assert.That(c.ExactBin, Is.EqualTo(exactBin));
        Assert.That(c.MismatchBin, Is.EqualTo(mismatchBin));
        Assert.That(c.Bin, Is.EqualTo(Math.Max(exactBin, mismatchBin)), "all-zero weights: higher bin wins the tie");
    }

    [Test]
    public void ClassifyRbsStart_NoMotif_Bin0() =>
        Assert.That(GenomeAnnotator.ClassifyRbsStart(new string('C', 40) + "ATG", 40).Bin, Is.EqualTo(0));

    [Test]
    public void ClassifyRbsStart_ReverseStrand_MirrorsForward()
    {
        string fwd = Spaced("AGGAGG", 7);
        string rev = DnaSequence.GetReverseComplementString(fwd);
        var f = GenomeAnnotator.ClassifyRbsStart(fwd, fwd.Length - 3, '+');
        // Reverse strand: the start codon's first base (on that strand) is the highest forward coordinate of its CAT.
        var r = GenomeAnnotator.ClassifyRbsStart(rev, 2, '-');

        Assert.That(r, Is.EqualTo(f));
    }

    [Test]
    public void ClassifyRbsStart_UsesModelWeights()
    {
        var c = GenomeAnnotator.ClassifyRbsStart(Spaced("AGGAGG", 7), 25 + 6 + 7, '+', Model.Value);

        Assert.That(c.Bin, Is.EqualTo(27));
        Assert.That(c.Weight, Is.EqualTo(Model.Value.BinWeights[27]));
    }

    #endregion

    #region FindRibosomeBindingSites score modes

    [Test]
    public void FindRibosomeBindingSites_LogOdds_ScoreIsBinWeight()
    {
        string seq = Spaced("AGGAGG", 7) + string.Concat(Enumerable.Repeat("GCT", 40)) + "TAA";
        var legacy = GenomeAnnotator.FindRibosomeBindingSites(seq).ToList();
        var logOdds = GenomeAnnotator.FindRibosomeBindingSites(seq, scoreMode: RbsScoreMode.LogOdds, model: Model.Value).ToList();

        Assert.That(legacy.Single().score, Is.EqualTo(1.0), "legacy default: motif length / 6");
        Assert.That(logOdds.Single().score, Is.EqualTo(Model.Value.BinWeights[27]));
        Assert.That(logOdds.Single().position, Is.EqualTo(legacy.Single().position));
    }

    [Test]
    public void FindRibosomeBindingSites_LogOdds_SpacerChangesScore()
    {
        string tail = string.Concat(Enumerable.Repeat("GCT", 40)) + "TAA";
        double At(int spacer) => GenomeAnnotator.FindRibosomeBindingSites(Spaced("AGGAGG", spacer) + tail,
            scoreMode: RbsScoreMode.LogOdds, model: Model.Value).Single().score;

        Assert.That(At(7), Is.EqualTo(Model.Value.BinWeights[27]));
        Assert.That(At(11), Is.EqualTo(Model.Value.BinWeights[25]));
        Assert.That(At(13), Is.EqualTo(Model.Value.BinWeights[10]));
        Assert.That(At(7), Is.Not.EqualTo(At(13)));
    }

    [Test]
    public void FindRibosomeBindingSitesBothStrands_LogOdds_ReverseHitsScored()
    {
        string fwd = Spaced("AGGAGG", 7) + string.Concat(Enumerable.Repeat("GCT", 40)) + "TAA";
        string rev = DnaSequence.GetReverseComplementString(fwd);
        var hits = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(rev, scoreMode: RbsScoreMode.LogOdds, model: Model.Value).ToList();

        Assert.That(hits.Single().strand, Is.EqualTo('-'));
        Assert.That(hits.Single().score, Is.EqualTo(Model.Value.BinWeights[27]));
    }

    #endregion

    #region Training from a supplied gene set

    [Test]
    public void TrainRbsModel_FromGenes_ProducesEnrichedSdBins()
    {
        var genes = OracleGenes().Select((g, i) => new GenomeAnnotator.GeneAnnotation(
            $"g{i}", g.begin - 1, g.end, g.strand, "CDS", "", new Dictionary<string, string>())).ToList();
        var model = GenomeAnnotator.TrainRbsModel(Genome.Value, genes);

        Assert.That(model.BinWeights, Has.Count.EqualTo(28));
        Assert.That(model.BinWeights[27], Is.GreaterThan(2.0), "AGGAGG at 5-10 nt is enriched among the supplied starts");
        Assert.That(model.BinWeights[0], Is.LessThan(0.0), "no motif is depleted among the supplied starts");
        Assert.That(model.StartTypeWeights[0], Is.GreaterThan(model.StartTypeWeights[2]), "ATG favoured over TTG");
    }

    #endregion
}
