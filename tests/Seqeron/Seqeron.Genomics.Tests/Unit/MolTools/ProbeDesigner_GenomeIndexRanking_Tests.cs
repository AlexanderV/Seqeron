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

    // --- bothStrands (audit round 4, A4-2, B07.md F60) ---
    // Same X / Y; the index holds X + T10 + Y + T10 + rc(X): every X window occurs once on the indexed strand and once
    // as its reverse complement, every Y window once. Python cross-check (str.find over the dumped strings, probe and
    // Biopython reverse complement, union of positions): of the 4686 candidate windows 3806 have 1 site on the indexed
    // strand and 2 on both strands, 297 (all inside Y) 1 / 1, 561 (crossing X|Y) 0 / 0 and 22 (crossing X|Y with only
    // Y's first 1–2 nt "AA", whose reverse complement "TT" lies in the T10 before rc(X)) 0 / 1.
    private static readonly string StrandGenome =
        Fixture.Target[..400] + new string('T', 10) + Fixture.Target[400..] + new string('T', 10)
        + DnaSequence.GetReverseComplementString(Fixture.Target[..400]);
    private static readonly global::SuffixTree.ISuffixTree StrandIndex = global::SuffixTree.SuffixTree.Build(StrandGenome);

    // Binding sites on both strands by naive scanning (independent of the suffix tree): the union of the overlapping
    // occurrence positions of the probe and of its reverse complement.
    private static int NaiveBothStrandSites(string probe)
    {
        var positions = new HashSet<int>();
        foreach (string pattern in new[] { probe, DnaSequence.GetReverseComplementString(probe) })
            for (int i = StrandGenome.IndexOf(pattern, StringComparison.Ordinal); i >= 0;
                 i = StrandGenome.IndexOf(pattern, i + 1, StringComparison.Ordinal))
                positions.Add(i);
        return positions.Count;
    }

    [Test]
    public void BothStrands_CandidateSiteCounts_MatchThePythonCrossCheck()
    {
        // Every candidate window of lengths 50–60 (the Microarray preset range): (indexed-strand sites, both-strand
        // sites) tallied; Python str.find on the dumped strings gives 3806 × (1, 2), 297 × (1, 1), 561 × (0, 0), 22 × (0, 1).
        var tally = new Dictionary<(int, int), int>();
        string target = Fixture.Target;
        for (int len = 50; len <= 60; len++)
            for (int start = 0; start + len <= target.Length; start++)
            {
                string w = target.Substring(start, len);
                double s1 = ProbeDesigner.CheckSpecificity(w, StrandIndex);
                double s2 = ProbeDesigner.CheckSpecificity(w, StrandIndex, bothStrands: true);
                var key = (s1 == 0 ? 0 : (int)Math.Round(1 / s1), s2 == 0 ? 0 : (int)Math.Round(1 / s2));
                tally[key] = tally.GetValueOrDefault(key) + 1;
            }

        Assert.That(tally, Is.EquivalentTo(new Dictionary<(int, int), int>
        {
            [(1, 2)] = 3806, [(1, 1)] = 297, [(0, 0)] = 561, [(0, 1)] = 22,
        }));
    }
}
