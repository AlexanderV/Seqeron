// CODON-RARE-001 — MCP: Sherlocc threshold comparison (inclusive_threshold) and %MinMax usage basis
// (FIN-B02 WP7, A1-1 / A1-3). Values are the reference locks of
// tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_SherloccThreshold_PerThousandPresets_Tests.cs
// (Sherlocc sherlocc.pl `$avrg <= $avrg_lim`; Clark-lab calculateMinMax over Kazusa 4932 / 9606 per thousand).
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class SherloccThresholdPerThousandToolsTests
{
    private static CodonUsageTableInput Preset(string id) => new(Preset: id);
    private const string SixAtThresholdPlusCommon = "UCUUCCUCGAGUGUAGGGCGC";
    private const string Gene20 = "AUGCUGAGAAAAGAAGCUCUGAGGUUUCCCGGAAUAUCUCGAUGGAAGACCCUAGUUCAG";

    [Test]
    public void FindRareCodonClusters_DefaultInclusive_StrictOnRequest()
    {
        var inclusive = MolToolsTools.find_rare_codon_clusters(SixAtThresholdPlusCommon, Preset("EColiK12")).Clusters;
        Assert.That(inclusive, Is.EqualTo(new[] { new RareCodonClusterDto(0, 6, 6) }));

        var strict = MolToolsTools.find_rare_codon_clusters(SixAtThresholdPlusCommon, Preset("EColiK12"), inclusive_threshold: false).Clusters;
        Assert.That(strict, Is.Empty);
    }

    [Test]
    public void FindRareCodons_DefaultStrict_InclusiveOnRequest()
    {
        var strict = MolToolsTools.find_rare_codons("UCUAGA", Preset("EColiK12")).RareCodons;
        Assert.That(strict.Select(r => r.Position), Is.EqualTo(new[] { 3 }));

        var inclusive = MolToolsTools.find_rare_codons("UCUAGA", Preset("EColiK12"), inclusive_threshold: true).RareCodons;
        Assert.That(inclusive.Select(r => (r.Position, r.Codon, r.Frequency)),
            Is.EqualTo(new[] { (0, "UCU", 0.15), (3, "AGA", 0.04) }));
    }

    [Test]
    public void CodonMinMaxProfile_PresetsUsePerThousandByDefault()
    {
        Assert.That(MolToolsTools.codon_minmax_profile("CUGAGA", Preset("Yeast"), 2).Windows[0].PercentMinMax,
            Is.EqualTo(33.99209486166008).Within(1e-9));
        Assert.That(MolToolsTools.codon_minmax_profile(Gene20, Preset("Human")).Windows.Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 28.796581754809075, 23.264123403923747, 16.195449171123403 }).Within(1e-9));

        Assert.That(MolToolsTools.codon_minmax_profile("CUGAGA", Preset("Yeast"), 2,
                CodonOptimizer.MinMaxUsageBasis.CodonFrequencies).Windows[0].PercentMinMax,
            Is.EqualTo(58.620689655172406).Within(1e-9));
        // SRC RO1 (2026-10-14): EColiK12 now carries the Kazusa 316407 per-thousand column → CHARMING
        // value 58.1 (port 58.10113519091846); was 36.470588235294116, still the CodonFrequencies value.
        Assert.That(MolToolsTools.codon_minmax_profile("CUGAGA", Preset("EColiK12"), 2).Windows[0].PercentMinMax,
            Is.EqualTo(58.10113519091846).Within(1e-9));
        Assert.That(MolToolsTools.codon_minmax_profile("CUGAGA", Preset("EColiK12"), 2,
                CodonOptimizer.MinMaxUsageBasis.CodonFrequencies).Windows[0].PercentMinMax,
            Is.EqualTo(36.470588235294116).Within(1e-10));
    }

    [Test]
    public void CodonMinMaxProfile_UnknownBasis_Rejected()
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            MolToolsTools.codon_minmax_profile("CUGAGA", Preset("Yeast"), 2, (CodonOptimizer.MinMaxUsageBasis)9));
        Assert.That(ex!.ParamName, Is.EqualTo("usage_basis"));
    }
}
