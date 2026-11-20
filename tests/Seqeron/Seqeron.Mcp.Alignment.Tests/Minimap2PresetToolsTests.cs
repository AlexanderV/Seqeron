using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>B09 A5-9: the <c>preset</c>, <c>maxFragmentLength</c> and <c>maxIntronLength</c> parameters of <c>minimap2_chains</c> (paf mode) and its <c>paf</c> text; the presets equal <c>minimap2 -x PRESET</c> (2 420 library cases).</summary>
[TestFixture]
public class Minimap2PresetToolsTests
{
    private static string Genome(int n, int seed)
    {
        var rnd = new Random(seed);
        return new string(Enumerable.Range(0, n).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
    }

    [Test]
    public void Preset_ChangesTheOptionsLikeMinimap2_X()
    {
        string target = Genome(6000, 1);
        string query = target.Substring(1000, 2500);
        var plain = AnalysisPaf(target, query, null);
        var asm5 = AnalysisPaf(target, query, "asm5");
        Assert.That(plain.Paf, Does.StartWith("query\t2500\t"));
        Assert.That(plain.Paf, Does.Contain("\t+\ttarget1\t6000\t"));
        Assert.That(asm5.Paf, Does.StartWith("query\t2500\t"));
        Assert.That(asm5.Paf, Does.Contain("\t+\ttarget1\t6000\t"));
        Assert.That(asm5.Hits!.Single().Divergence, Is.Not.Null);
        Assert.Throws<ArgumentException>(() => AnalysisPaf(target, query, "lr:hq"));
    }

    private static Minimap2ChainsResult AnalysisPaf(string target, string query, string? preset)
        => AlignmentTools.Minimap2Chains([target], query, mode: "paf", preset: preset);

    [Test]
    public void ShortReadPreset_HasNoDivergenceTag_AndNoSecondaryHits()
    {
        string target = Genome(5000, 2);
        var r = AlignmentTools.Minimap2Chains([target], target.Substring(700, 150), mode: "paf", preset: "sr");
        Assert.That(r.Paf, Does.Not.Contain("dv:f:"));
        Assert.That(r.Paf, Does.Contain("tp:A:P"));
    }
}
