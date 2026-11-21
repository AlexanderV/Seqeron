using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>B09 A5-9: the <c>preset</c>, <c>maxFragmentLength</c> and <c>maxIntronLength</c> parameters of <c>minimap2_chains</c> (paf mode) and its <c>paf</c> text; the presets equal <c>minimap2 -x PRESET</c> (2 420 library cases).</summary>
[TestFixture]
public class Minimap2PresetToolsTests
{
    // a 64-bit LCG (x <- x * 6364136223846793005 + 1442695040888963407, base = x >> 62), the same as the Python script that produced the expected lines
    private static string Genome(int n, int seed)
    {
        ulong x = (ulong)seed;
        var chars = new char[n];
        for (int i = 0; i < n; i++)
        {
            x = unchecked(x * 6364136223846793005UL + 1442695040888963407UL);
            chars[i] = "ACGT"[(int)(x >> 62)];
        }

        return new string(chars);
    }

    private static readonly string Target = Genome(6000, 1);

    /// <summary>target[1000, 3500) with every 97th base changed (A becomes C, any other base A).</summary>
    private static readonly string Query = string.Concat(Target.Substring(1000, 2500).Select((c, i) => i % 97 == 0 ? (c != 'A' ? 'A' : 'C') : c));

    private static string Paf(string? preset, int? k = null, int? w = null, int? minChainScore = null, int? minAnchors = null, bool? rmq = null, bool? outputSecondary = null, bool homopolymerCompressed = false)
        => AlignmentTools.Minimap2Chains([Target], Query, k: k, w: w, minChainScore: minChainScore, minAnchors: minAnchors, mode: "paf", rmq: rmq, outputSecondary: outputSecondary, homopolymerCompressed: homopolymerCompressed, preset: preset).Paf!;

    [Test]
    public void Preset_ChangesTheOptionsLikeMinimap2_X()
    {
        // minimap2 2.26 `minimap2 target1.fa query.fa` and `-x asm5` on the same pair (whole PAF lines, installed program)
        Assert.That(Paf(null), Is.EqualTo("query\t2500\t9\t2494\t+\ttarget1\t6000\t1009\t3494\t2302\t2485\t60\ttp:A:P\tcm:i:386\ts1:i:2302\ts2:i:0\tdv:f:0.0002\trl:i:0\n"));
        Assert.That(Paf("asm5"), Is.EqualTo("query\t2500\t9\t2486\t+\ttarget1\t6000\t1009\t3486\t2069\t2477\t60\ttp:A:P\tcm:i:188\ts1:i:2069\ts2:i:0\tdv:f:0.0003\trl:i:0\n"));
        Assert.That(Paf("map-pb"), Is.EqualTo("query\t2500\t9\t2500\t+\ttarget1\t6000\t1009\t3500\t2283\t2496\t60\ttp:A:P\tcm:i:262\ts1:i:2278\ts2:i:0\tdv:f:0.0003\trl:i:0\n"));
        Assert.That(AlignmentTools.Minimap2Chains([Target], Query, mode: "paf", preset: "asm5").Hits!.Single().Divergence, Is.Not.Null);
        Assert.Throws<ArgumentException>(() => Paf("lr:hq"));
    }

    [Test]
    public void ExplicitValue_EqualToTheDefault_OverridesThePreset()
    {
        // minimap2 -x asm5 -k 15 is k = 15 (the preset's 19 is replaced by a value that equals the DEFAULT): a different chain than -x asm5
        Assert.That(Paf("asm5", k: 15), Is.EqualTo("query\t2500\t9\t2494\t+\ttarget1\t6000\t1009\t3494\t2051\t2485\t60\ttp:A:P\tcm:i:213\ts1:i:2051\ts2:i:0\tdv:f:0.0003\trl:i:0\n"));
        Assert.That(Paf("asm5", k: 15), Is.Not.EqualTo(Paf("asm5")));
        // -x asm5 -w 10 -m 40 -n 3 (all three equal the defaults, the preset has w = 19)
        Assert.That(Paf("asm5", w: 10, minChainScore: 40, minAnchors: 3), Is.EqualTo("query\t2500\t9\t2496\t+\ttarget1\t6000\t1009\t3496\t2314\t2487\t60\ttp:A:P\tcm:i:348\ts1:i:2314\ts2:i:0\tdv:f:0.0002\trl:i:0\n"));
        // an omitted option keeps the preset's value, also when the option is passed as null
        Assert.That(Paf("asm5", k: null, w: null), Is.EqualTo(Paf("asm5")));
        // -m 40 after -x ava-ont (preset: 100) is a primary-less all-vs-all hit with the longer chain list (tp:A:S, no s2)
        Assert.That(AlignmentTools.Minimap2Chains([Target], Query, minChainScore: 40, mode: "paf", preset: "ava-ont").Paf,
            Is.EqualTo("query\t2500\t5\t2499\t+\ttarget1\t6000\t1005\t3499\t2400\t2494\t0\ttp:A:S\tcm:i:712\ts1:i:2400\tdv:f:0.0001\trl:i:0\n"));
        // the flags that minimap2 can switch off (--rmq=no, --secondary=yes) are explicit booleans; -H only switches on (false = not given: map-pb keeps it)
        Assert.That(Paf("asm5", rmq: false), Is.EqualTo(Paf("asm5")));                      // equal on this pair; the 440 override cases of the library oracle (MinimapPresets.Apply) tell them apart
        Assert.That(Paf("map-pb", homopolymerCompressed: false), Is.EqualTo(Paf("map-pb")));
        Assert.That(Paf("sr", outputSecondary: true), Does.Contain("tp:A:P"));
    }

    [Test]
    public void PresetOverrides_AreTheLibrarysApply_NotAThirdCopyOfTheMerge()
    {
        // the tool and the oracle-checked library function build the same MinimapMapOptions: every given parameter is an override, the preset supplies the rest
        var viaTool = AlignmentTools.Minimap2Chains([Target], Query, k: 15, rmq: false, outputSecondary: true, bestN: 3, mode: "paf", preset: "asm5").Paf;
        var options = global::Seqeron.Genomics.Alignment.MinimapPresets.Apply("asm5", new global::Seqeron.Genomics.Alignment.MinimapOverrides { KmerSize = 15, UseRmq = false, OutputSecondary = true, BestN = 3 });
        Assert.That((options.Seed.KmerSize, options.Seed.WindowSize, options.UseRmq, options.OutputSecondary, options.BestN, options.Bandwidth), Is.EqualTo((15, 19, false, true, 3, 1000)));
        var mapped = global::Seqeron.Genomics.Alignment.Minimap2Chaining.Map([Target], Query, options, null, null);
        Assert.That(viaTool, Is.EqualTo(global::Seqeron.Genomics.Alignment.Minimap2Chaining.FormatPaf(mapped, "query", Query.Length, ["target1"], [Target.Length])));
    }

    [Test]
    public void ShortReadPreset_HasNoDivergenceTag_AndNoSecondaryHits()
    {
        string target = Genome(5000, 2);
        var r = AlignmentTools.Minimap2Chains([target], target.Substring(700, 150), mode: "paf", preset: "sr");
        Assert.That(r.Paf, Does.Not.Contain("dv:f:"));
        Assert.That(r.Paf, Does.Contain("tp:A:P"));
        Assert.That(r.Hits!.All(h => h.Primary), Is.True, "-x sr prints no secondary hit");
    }
}
