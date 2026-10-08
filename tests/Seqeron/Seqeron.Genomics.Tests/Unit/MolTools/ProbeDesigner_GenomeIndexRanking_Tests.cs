// PROBE-DESIGN-001 — DesignProbes(target, ISuffixTree genomeIndex, …) (audit round 4, A4-1, B07.md F59): the overload
// checked only the top maxProbes × 5 raw-score candidates for specificity and, with requireUnique = false, scaled Score
// by the specificity without re-ranking. Contract (XML / Hybridization_Probe_Design.md §3.1): "maximum probes to return",
// "only unique probes" — i.e. the best maxProbes of ALL candidates after the specificity step, in the documented order.
// Oracles (independent of the overload under test): the exhaustive non-index DesignProbes (maxProbes = int.MaxValue =
// every score-positive candidate in the documented order) + an occurrence count by naive string scanning (overlapping
// matches), cross-checked with Python str.find on the same sequences.
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class ProbeDesigner_GenomeIndexRanking_Tests
{
    // Auditor's repro: X (400 nt) and Y (80 nt) drawn from System.Random(7) as uniform ACGT, in that order;
    // genome = X + T10 + X + T10 + Y, target = X + Y. Every 50–60-mer inside X occurs twice, every one inside Y once,
    // and the top raw-score candidates all lie in X — so the old maxProbes × 5 shortlist contained no unique probe.
    private static readonly (string Target, string Genome) Fixture = BuildFixture();
    private static readonly global::SuffixTree.ISuffixTree Index = global::SuffixTree.SuffixTree.Build(Fixture.Genome);

    private static (string, string) BuildFixture()
    {
        var rng = new Random(7);
        string Draw(int n) => new(Enumerable.Range(0, n).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        string x = Draw(400), y = Draw(80);
        string t10 = new('T', 10);
        return (x + y, x + t10 + x + t10 + y);
    }

    // Microarray preset with the cheap (Primer3 alignment / inverted-repeat) self-structure screens instead of the
    // per-candidate ntthal runs: the screens are orthogonal to the specificity walk under test, and the exhaustive
    // oracle has to evaluate every one of the ~4 000 candidate windows.
    private static ProbeDesigner.ProbeParameters FastScreen => ProbeDesigner.Defaults.Microarray with
    {
        StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic,
    };

    // Exhaustive candidate list per parameter set (one full non-index design), reused by every oracle below
    // (this assembly runs fixture tests in parallel, so the memo is locked).
    private static readonly Dictionary<ProbeDesigner.ProbeParameters, List<ProbeDesigner.Probe>> AllCandidates = new();

    private static List<ProbeDesigner.Probe> Candidates(ProbeDesigner.ProbeParameters param)
    {
        lock (AllCandidates)
        {
            if (!AllCandidates.TryGetValue(param, out var all))
                AllCandidates[param] = all = ProbeDesigner.DesignProbes(Fixture.Target, param, int.MaxValue).ToList();
            return all;
        }
    }

    // Overlapping exact occurrences by naive scanning (independent of the suffix tree).
    private static int NaiveCount(string probe)
    {
        int count = 0;
        for (int i = Fixture.Genome.IndexOf(probe, StringComparison.Ordinal); i >= 0;
             i = Fixture.Genome.IndexOf(probe, i + 1, StringComparison.Ordinal))
            count++;
        return count;
    }

    private static double NaiveSpecificity(string probe)
    {
        int hits = NaiveCount(probe);
        return hits == 0 ? 0 : 1.0 / hits;
    }

    // Brute force of the contract: every candidate, the specificity of each, then the documented order —
    // requireUnique: the unique ones in the base ranking order; otherwise Score × specificity, re-ranked stably
    // (score descending, ties in the base order = (length, start) ascending for the additive ranking; the Primer3
    // penalty is not touched by the specificity, so that ranking keeps its order).
    private static List<ProbeDesigner.Probe> Oracle(
        ProbeDesigner.ProbeParameters param, int maxProbes, bool requireUnique)
    {
        var all = Candidates(param);
        if (requireUnique)
            return all.Where(p => NaiveSpecificity(p.Sequence) == 1.0).Take(maxProbes).ToList();
        var scaled = all.Select(p => p with { Score = p.Score * NaiveSpecificity(p.Sequence) });
        return (param.Ranking == ProbeDesigner.ProbeRanking.AdditiveScore
                ? scaled.OrderByDescending(p => p.Score).ThenBy(p => p.Sequence.Length).ThenBy(p => p.Start)
                : scaled)
            .Take(maxProbes).ToList();
    }

    private static void AssertSameProbes(
        IReadOnlyList<ProbeDesigner.Probe> actual, IReadOnlyList<ProbeDesigner.Probe> expected)
    {
        Assert.That(actual.Select(p => (p.Start, p.Sequence)), Is.EqualTo(expected.Select(p => (p.Start, p.Sequence))));
        for (int i = 0; i < expected.Count; i++)
            Assert.That(actual[i].Score, Is.EqualTo(expected[i].Score).Within(1e-12), $"Score of probe {i}");
    }

    [Test]
    public void Repro_SequencesMatchTheIndependentOccurrenceCounts()
    {
        // Python cross-check on the same two strings (str.find loop): the 50-mer at target 402 occurs once, the one at
        // 100 twice, and a window crossing the X|Y junction not at all.
        Assert.Multiple(() =>
        {
            Assert.That(Fixture.Target, Has.Length.EqualTo(480));
            Assert.That(Fixture.Genome, Has.Length.EqualTo(900));
            Assert.That(NaiveCount(Fixture.Target.Substring(402, 50)), Is.EqualTo(1));
            Assert.That(NaiveCount(Fixture.Target.Substring(100, 50)), Is.EqualTo(2));
            Assert.That(NaiveCount(Fixture.Target.Substring(380, 50)), Is.EqualTo(0));
            Assert.That(ProbeDesigner.CheckSpecificity(Fixture.Target.Substring(402, 50), Index), Is.EqualTo(1.0));
            Assert.That(ProbeDesigner.CheckSpecificity(Fixture.Target.Substring(100, 50), Index), Is.EqualTo(0.5));
        });
    }

    [Test]
    public void RequireUnique_MicroarrayPreset_FindsTheUniqueProbesBeyondTheOldShortlist()
    {
        // Auditor's repro with the preset exactly as documented (thermodynamic screen): maxProbes = 2 used to return 0
        // probes because the top 10 raw-score candidates all lie in the duplicated X; the unique probes are the
        // 50-mers at target 402 and 417, additive score 0.85.
        var probes = ProbeDesigner.DesignProbes(
            Fixture.Target, Index, ProbeDesigner.Defaults.Microarray, maxProbes: 2, requireUnique: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(probes.Select(p => p.Start), Is.EqualTo(new[] { 402, 417 }));
            Assert.That(probes.Select(p => p.Sequence.Length), Is.All.EqualTo(50));
            Assert.That(probes.Select(p => p.Score), Is.All.EqualTo(0.85).Within(1e-12));
            foreach (var p in probes)
            {
                Assert.That(NaiveCount(p.Sequence), Is.EqualTo(1), $"probe at {p.Start} unique (naive scan)");
                Assert.That(ProbeDesigner.CheckSpecificity(p.Sequence, Index), Is.EqualTo(1.0));
            }
        });
    }

    [TestCase(1)]
    [TestCase(3)]
    [TestCase(100_000)]
    public void RequireUnique_EqualsBruteForceOverAllCandidates(int maxProbes)
    {
        var probes = ProbeDesigner.DesignProbes(Fixture.Target, Index, FastScreen, maxProbes, requireUnique: true)
            .ToList();

        Assert.Multiple(() =>
        {
            AssertSameProbes(probes, Oracle(FastScreen, maxProbes, requireUnique: true));
            Assert.That(probes, Is.Not.Empty);
            foreach (var p in probes)
                Assert.That(NaiveCount(p.Sequence), Is.EqualTo(1));
        });
    }

    [TestCase(2)]
    [TestCase(100_000)]
    public void NotRequireUnique_ReRanksByTheSpecificityScaledScore(int maxProbes)
    {
        // Windows inside X are scaled × 0.5, junction windows × 0 and Y windows × 1, so the Y probes lead: the output
        // must be the brute-force ranking of the scaled scores (descending, ties by length then start), not the raw one.
        var probes = ProbeDesigner.DesignProbes(Fixture.Target, Index, FastScreen, maxProbes, requireUnique: false)
            .ToList();

        Assert.Multiple(() =>
        {
            AssertSameProbes(probes, Oracle(FastScreen, maxProbes, requireUnique: false));
            for (int i = 1; i < probes.Count; i++)
                Assert.That(probes[i - 1].Score, Is.GreaterThanOrEqualTo(probes[i].Score), $"score-descending at {i}");
            Assert.That(probes[0].Start, Is.GreaterThanOrEqualTo(400), "a unique Y probe outranks the halved X probes");
        });
    }

    [Test]
    public void NotRequireUnique_TiesKeepTheDocumentedLengthThenStartOrder()
    {
        // All X windows of equal raw score are halved alike: among equal scaled scores the order stays
        // (length, start) ascending — the stable re-rank of the documented tie order.
        var probes = ProbeDesigner.DesignProbes(Fixture.Target, Index, FastScreen, 100_000, requireUnique: false)
            .ToList();

        bool sawTie = false;
        for (int i = 1; i < probes.Count; i++)
        {
            if (probes[i - 1].Score != probes[i].Score)
                continue;
            sawTie = true;
            var a = (probes[i - 1].Sequence.Length, probes[i - 1].Start);
            var b = (probes[i].Sequence.Length, probes[i].Start);
            Assert.That(a.CompareTo(b), Is.LessThan(0), $"tie order at {i}");
        }
        Assert.That(sawTie, Is.True, "the fixture has scaled-score ties");
    }

    [TestCase(true)]
    [TestCase(false)]
    public void Primer3Ranking_WalksAllCandidatesInPenaltyOrder(bool requireUnique)
    {
        var param = FastScreen with
        {
            Ranking = ProbeDesigner.ProbeRanking.Primer3Penalty, OptTm = 86, OptLength = 55,
        };

        var probes = ProbeDesigner.DesignProbes(Fixture.Target, Index, param, maxProbes: 3, requireUnique).ToList();

        Assert.Multiple(() =>
        {
            AssertSameProbes(probes, Oracle(param, 3, requireUnique));
            Assert.That(probes, Has.Count.EqualTo(3));
            if (requireUnique)
                Assert.That(probes.Select(p => p.Start), Is.All.InRange(400, 430),
                    "only the Y windows are unique, and F51's penalty order is kept");
            foreach (var p in probes)
                Assert.That(p.Primer3Penalty, Is.Not.Null);
        });
    }

    [Test]
    public void NoUniqueCandidate_ReturnsEmpty_AndNonPositiveMaxProbesReturnsEmpty()
    {
        string x = Fixture.Target[..400];
        var duplicated = global::SuffixTree.SuffixTree.Build(x + new string('T', 10) + x);
        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.DesignProbes(x, duplicated, FastScreen, 5, requireUnique: true), Is.Empty);
            Assert.That(ProbeDesigner.DesignProbes(x, duplicated, FastScreen, 0, requireUnique: false), Is.Empty);
            Assert.That(ProbeDesigner.DesignProbes(x, duplicated, FastScreen, -1, requireUnique: false), Is.Empty);
        });
    }
}
