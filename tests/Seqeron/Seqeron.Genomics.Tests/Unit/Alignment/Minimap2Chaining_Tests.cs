using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A3-9c, B09 F151): <see cref="Minimap2Chaining.Chain"/> = <c>mg_lchain_dp</c> + <c>mg_chain_backtrack</c> + <c>compact_a</c> of minimap2 2.26.
/// Oracle (TestData/Alignment/minimap2_chain_oracle.jsonl, scripts/alignment/generate_minimap2_chain_oracle.py): 450 related random pairs (errors 1-15 %, inversions, tandem repeats, two reference sequences) mapped by the installed
/// minimap2 2.26 with random <c>-k -w -g -r b,b -m -n</c>; <c>--print-seeds</c> gives minimap2's own sorted anchor array, <c>--print-chains</c> the chains it builds from it, the plain PAF the chain scores (<c>s1</c>) and anchor counts (<c>cm</c>).
/// </summary>
[TestFixture]
public class Minimap2Chaining_Tests
{
    private sealed record Case(int K, int W, int G, int Bw, int M, int N, string[] Targets, string Query, MinimapAnchor[] Seeds, int[][] Chains, (int Cm, int S1)[] Paf);

    private static readonly Lazy<List<Case>> Oracle = new(ReadOracle);

    private static List<Case> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.minimap2_chain_oracle.jsonl";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var cases = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            using var doc = JsonDocument.Parse(line);
            var r = doc.RootElement;
            var seeds = r.GetProperty("seeds").EnumerateArray()
                .Select(s => new MinimapAnchor(s[0].GetInt32(), s[1].GetInt32() == 1, s[2].GetInt32(), s[3].GetInt32(), s[4].GetInt32())).ToArray();
            var chains = r.GetProperty("chains").EnumerateArray().Select(c => c.EnumerateArray().Select(e => e.GetInt32()).ToArray()).ToArray();
            var paf = r.GetProperty("paf").EnumerateArray().Select(c => (c[0].GetInt32(), c[1].GetInt32())).ToArray();
            cases.Add(new Case(r.GetProperty("k").GetInt32(), r.GetProperty("w").GetInt32(), r.GetProperty("g").GetInt32(), r.GetProperty("bw").GetInt32(),
                r.GetProperty("m").GetInt32(), r.GetProperty("n").GetInt32(),
                r.GetProperty("targets").EnumerateArray().Select(e => e.GetString()!).ToArray(), r.GetProperty("query").GetString()!, seeds, chains, paf));
        }
        return cases;
    }

    private static MinimapChainOptions OptionsOf(Case c) => new()
    {
        KmerSize = c.K, MaxTargetGap = c.G, MaxQueryGap = c.G, Bandwidth = c.Bw, MinCount = c.N, MinScore = c.M,
    };

    [Test]
    public void Chain_EqualsMinimap2_AllCases()
    {
        var cases = Oracle.Value;
        Assert.That(cases, Has.Count.EqualTo(450));
        var failures = new List<string>();
        int withChains = 0, multi = 0, reverse = 0, twoTargets = 0, scoresChecked = 0, many = 0;
        foreach (var c in cases)
        {
            var chains = Minimap2Chaining.Chain(c.Seeds, OptionsOf(c));
            var actual = chains.Select(x => x.Anchors.ToArray()).OrderBy(x => x[0]).ToArray();
            bool equal = actual.Length == c.Chains.Length && actual.Zip(c.Chains).All(p => p.First.SequenceEqual(p.Second));
            if (!equal)
                failures.Add($"k={c.K} w={c.W} g={c.G} bw={c.Bw} m={c.M} n={c.N} seeds={c.Seeds.Length}\n expected {string.Join(" | ", c.Chains.Select(x => string.Join(",", x)))}\n actual   {string.Join(" | ", actual.Select(x => string.Join(",", x)))}");
            // every PAF line (cm anchors, chain score s1) is one of the chains
            var mine = chains.Select(x => (x.Anchors.Count, x.Score)).ToList();
            foreach (var line in c.Paf)
            {
                if (!mine.Remove(line))
                    failures.Add($"PAF line cm={line.Cm} s1={line.S1} is not a chain score of the port (k={c.K} bw={c.Bw})");
                scoresChecked++;
            }

            if (chains.Count > 0) withChains++;
            if (chains.Count > 1) multi++;
            if (c.Seeds.Any(s => s.Reverse)) reverse++;
            if (c.Seeds.Any(s => s.TargetId > 0)) twoTargets++;
            if (c.Seeds.Length > 100) many++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        Assert.That(withChains, Is.GreaterThan(350));
        Assert.That(multi, Is.GreaterThan(100));
        Assert.That(reverse, Is.GreaterThan(100));
        Assert.That(twoTargets, Is.GreaterThan(50));
        Assert.That(scoresChecked, Is.GreaterThan(350));
        Assert.That(many, Is.GreaterThan(100));
    }

    [Test]
    public void CollectAnchors_EqualsTheSeedArrayOfMinimap2_AllCases()
    {
        // minimizer sketch (mm_sketch), hash index, mm_collect_matches / collect_seed_hits and the sort: minimap2 --print-seeds on the same sequences
        int checkedSeeds = 0;
        foreach (var c in Oracle.Value)
        {
            var actual = Minimap2Chaining.CollectAnchors(c.Targets, c.Query, c.K, c.W);
            Assert.That(actual, Is.EqualTo(c.Seeds), $"k={c.K} w={c.W} targets={c.Targets.Length} query={c.Query.Length}");
            checkedSeeds += actual.Count;
        }

        Assert.That(checkedSeeds, Is.GreaterThan(30_000));
    }

    [Test]
    public void MapChains_EqualsTheChainsOfMinimap2_FromTheSequences()
    {
        int chains = 0, reverse = 0;
        foreach (var c in Oracle.Value)
        {
            var hits = Minimap2Chaining.MapChains(c.Targets, c.Query, c.K, c.W, OptionsOf(c));
            var actual = hits.Select(h => h.Anchors.Select(a => (a.TargetId, a.Reverse, a.TargetEnd, a.QueryEnd, a.Span)).ToArray()).OrderBy(x => x[0].TargetEnd).ToArray();
            var expected = c.Chains.Select(ch => ch.Select(i => (c.Seeds[i].TargetId, c.Seeds[i].Reverse, c.Seeds[i].TargetEnd, c.Seeds[i].QueryEnd, c.Seeds[i].Span)).ToArray()).ToArray();
            Assert.That(actual.Length, Is.EqualTo(expected.Length), $"k={c.K} w={c.W}");
            var byFirst = expected.OrderBy(x => x[0].TargetEnd).ToArray();
            for (int i = 0; i < actual.Length; i++)
                Assert.That(actual[i], Is.EqualTo(byFirst[i]), $"k={c.K} w={c.W} chain {i}");
            chains += hits.Count;
            reverse += hits.Count(h => h.Reverse);
        }

        Assert.That(chains, Is.GreaterThan(800));
        Assert.That(reverse, Is.GreaterThan(100));
    }

    [Test]
    public void CollectAnchors_FrequentMinimizers_AreNotPorted_AndValidation()
    {
        // mid_occ = max(10, the 0.9998 quantile of the occurrence counts + 1): on a 60 kb reference with 3 000 poly-A the poly-A minimizer is the one frequent key (3 000 > 10 occurrences)
        var rng = new Random(5);
        string reference = new string(Enumerable.Range(0, 60_000).Select(_ => "ACGT"[rng.Next(4)]).ToArray()) + new string('A', 3000);
        Assert.Throws<NotSupportedException>(() => Minimap2Chaining.CollectAnchors([reference], new string('A', 24) + reference[..200]));
        Assert.DoesNotThrow(() => Minimap2Chaining.CollectAnchors([reference], reference.Substring(1000, 300)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Minimap2Chaining.CollectAnchors(["ACGT"], "ACGT", 29, 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => Minimap2Chaining.CollectAnchors(["ACGT"], "ACGT", 15, 256));
        Assert.That(Minimap2Chaining.CollectAnchors(["ACGT"], "ACGT"), Is.Empty);
    }

    [Test]
    public void SmallChain_ScoreIsTheSumOfAnchorSpansMinusGapPenalties()
    {
        // five colinear 15-mers every 20 bases on the same diagonal: score = 15 + 4 * min(20, 15) = 75, one chain of five
        var anchors = Enumerable.Range(0, 5).Select(i => new MinimapAnchor(0, false, 100 + 20 * i, 200 + 20 * i, 15)).ToList();
        var chains = Minimap2Chaining.Chain(anchors);
        Assert.That(chains, Has.Count.EqualTo(1));
        Assert.That((chains[0].Score, chains[0].Anchors), Is.EqualTo((75, new[] { 0, 1, 2, 3, 4 })));
        // below -n anchors or -m score: no chain
        Assert.That(Minimap2Chaining.Chain(anchors.Take(2).ToList()), Is.Empty);
        Assert.That(Minimap2Chaining.Chain(anchors, new MinimapChainOptions { MinScore = 76 }), Is.Empty);
    }

    [Test]
    public void Validation_AndSorting()
    {
        Assert.Throws<ArgumentNullException>(() => Minimap2Chaining.Chain(null!));
        Assert.That(Minimap2Chaining.Chain([]), Is.Empty);
        var unsorted = new[] { new MinimapAnchor(0, false, 50, 60, 15), new MinimapAnchor(0, false, 10, 20, 15) };
        Assert.Throws<ArgumentException>(() => Minimap2Chaining.Chain(unsorted));
        Assert.That(Minimap2Chaining.SortAnchors(unsorted).Select(a => a.TargetEnd), Is.EqualTo(new[] { 10, 50 }));
    }
}
