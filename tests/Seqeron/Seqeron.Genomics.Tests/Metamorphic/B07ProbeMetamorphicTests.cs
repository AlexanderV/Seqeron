using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.ProbeDesigner;

namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic tests for the genome-index <c>DesignProbes</c> overload after review batch B07
/// (docs/Validation/review-2026-09/B07.md F59 — every candidate walked lazily and the specificity-scaled scores
/// re-ranked; F60 — both-strand specificity). The source relation is the base design without an index: the
/// index overload must equal that full ranking filtered (requireUnique) or rescored and stably re-ranked
/// (requireUnique = false), truncated to maxProbes, for both rankings. Genomes are built from the target's own pieces
/// (repeats on either strand) plus random DNA, so specificities 0, 1 and 1/n all occur. Short targets and the cheap
/// fallback screen keep each case fast (the screen does not enter the relations).
///
/// Test Units: PROBE-DESIGN-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("MolTools")]
public class B07ProbeMetamorphicTests
{
    private static Gen<string> AcgtGen(int minLen, int maxLen) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
        select new string(chars);

    private static string RevComp(string s)
    {
        var r = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(r);
    }

    /// <summary>Independent occurrence count (overlapping) of a pattern in a text.</summary>
    private static int Count(string text, string pattern)
    {
        int c = 0;
        for (int i = text.IndexOf(pattern, StringComparison.Ordinal); i >= 0; i = text.IndexOf(pattern, i + 1, StringComparison.Ordinal))
            c++;
        return c;
    }

    private static (int, string) Id(Probe p) => (p.Start, p.Sequence);

    private sealed record Scenario(string Target, string Genome, ProbeParameters Param, int MaxProbes)
    {
        public override string ToString() => $"target {Target} genome({Genome.Length}) {Param.MinLength}-{Param.MaxLength} {Param.Ranking} k={MaxProbes}";
    }

    /// <summary>Target 40–80 nt; genome = random + target (optionally) + 0–3 copies of target pieces, each on either strand.</summary>
    private static Gen<Scenario> ScenarioGen() =>
        from target in AcgtGen(40, 80)
        from includeTarget in Gen.Elements(true, true, true, false)
        from nPieces in Gen.Choose(0, 3)
        from pieces in (
            from start in Gen.Choose(0, 1000)
            from len in Gen.Choose(14, 40)
            from rc in Gen.Elements(false, true)
            select (start, len, rc)).ArrayOf(nPieces)
        from fills in AcgtGen(0, 30).ArrayOf(nPieces + 2)
        from minLen in Gen.Choose(14, 18)
        from extra in Gen.Choose(0, 3)
        from minTm in Gen.Choose(30, 55)
        from span in Gen.Choose(0, 20)
        from primer3 in Gen.Elements(false, false, true)
        from k in Gen.Choose(1, 15)
        let copies = pieces.Select(p =>
        {
            int s = p.start % target.Length;
            string piece = target.Substring(s, Math.Min(p.len, target.Length - s));
            return p.rc ? RevComp(piece) : piece;
        }).ToArray()
        let genome = fills[0] + (includeTarget ? target : "") + string.Concat(copies.Select((c, i) => fills[i + 1] + c)) + fills[^1]
        let param = new ProbeParameters(minLen, minLen + extra, minTm, minTm + span, 0.3, 0.7, 5, false, 0.5)
        {
            StructureScreen = ProbeStructureScreen.Heuristic,
            Ranking = primer3 ? ProbeRanking.Primer3Penalty : ProbeRanking.AdditiveScore,
            OptTm = minTm + span / 2.0,
            OptLength = minLen,
        }
        select new Scenario(target, genome.Length == 0 ? "A" : genome, param, k);

    /// <summary>
    /// F59 (requireUnique): the index overload = the full base ranking (no index) filtered to the probes with
    /// <c>CheckSpecificity == 1</c> (independently: exactly one occurrence), truncated to maxProbes, records unchanged.
    /// F60: with bothStrands the filter counts the probe and its reverse complement (a reverse palindrome once); every
    /// such probe has <c>CheckSpecificity(…, true) == 1</c> and occurs at most once on the indexed strand (it may occur
    /// there only as its reverse complement).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property GenomeIndex_RequireUnique_EqualsFilteredBaseRanking_Truncated()
    {
        return Prop.ForAll(ScenarioGen().ToArbitrary(), sc =>
        {
            var index = global::SuffixTree.SuffixTree.Build(sc.Genome);
            var full = DesignProbes(sc.Target, sc.Param, 100_000).ToList();
            foreach (bool both in new[] { false, true })
            {
                var got = DesignProbes(sc.Target, index, sc.Param, sc.MaxProbes, requireUnique: true, bothStrands: both).ToList();
                var expected = full.Where(p =>
                {
                    string rc = RevComp(p.Sequence);
                    int n = Count(sc.Genome, p.Sequence) + (both && rc != p.Sequence ? Count(sc.Genome, rc) : 0);
                    return n == 1;
                }).Take(sc.MaxProbes).ToList();
                if (!got.Select(Id).SequenceEqual(expected.Select(Id))
                    || got.Zip(expected).Any(t => t.First.Score != t.Second.Score || t.First.Primer3Penalty != t.Second.Primer3Penalty))
                    return false.Label($"both={both}: got [{string.Join(",", got.Select(Id))}] expected [{string.Join(",", expected.Select(Id))}] {sc}");
                if (got.Any(p => CheckSpecificity(p.Sequence, index, both) != 1.0))
                    return false.Label($"both={both}: a returned probe is not unique {sc}");
                if (both && got.Any(p => CheckSpecificity(p.Sequence, index, false) is not (0.0 or 1.0)))
                    return false.Label($"both-strand-unique probe occurs more than once on the indexed strand {sc}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// F59 (requireUnique = false): the index overload = the full base ranking with <c>Score × CheckSpecificity</c>,
    /// re-ranked score-descending with the documented (length, start) tie order (additive ranking) or kept in the
    /// Primer3 penalty order (the specificity does not enter the penalty), truncated to maxProbes; the output is a prefix
    /// of the uncapped output; with the additive ranking it is score-descending. F60: both-strand scores never exceed
    /// the single-strand scores of a probe that occurs on the indexed strand (a reverse-complement occurrence only adds
    /// to its count).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property GenomeIndex_Lenient_EqualsRescoredStableReranking_PrefixClosed()
    {
        return Prop.ForAll(ScenarioGen().ToArbitrary(), sc =>
        {
            var index = global::SuffixTree.SuffixTree.Build(sc.Genome);
            var full = DesignProbes(sc.Target, sc.Param, 100_000).ToList();
            var byBoth = new Dictionary<bool, Dictionary<(int, string), double>>();
            foreach (bool both in new[] { false, true })
            {
                var rescored = full.Select(p =>
                {
                    string rc = RevComp(p.Sequence);
                    int n = Count(sc.Genome, p.Sequence) + (both && rc != p.Sequence ? Count(sc.Genome, rc) : 0);
                    return p with { Score = p.Score * (n == 0 ? 0 : 1.0 / n) };
                }).ToList();
                var expectedAll = sc.Param.Ranking == ProbeRanking.Primer3Penalty
                    ? rescored
                    : rescored.OrderByDescending(p => p.Score).ThenBy(p => p.Sequence.Length).ThenBy(p => p.Start).ToList();
                var all = DesignProbes(sc.Target, index, sc.Param, 100_000, requireUnique: false, bothStrands: both).ToList();
                var got = DesignProbes(sc.Target, index, sc.Param, sc.MaxProbes, requireUnique: false, bothStrands: both).ToList();
                if (!all.Select(Id).SequenceEqual(expectedAll.Select(Id)) || all.Zip(expectedAll).Any(t => t.First.Score != t.Second.Score))
                    return false.Label($"both={both}: full re-ranking differs {sc}\n got [{string.Join(",", all.Select(p => $"{p.Start}/{p.Sequence.Length}:{p.Score}"))}]\n exp [{string.Join(",", expectedAll.Select(p => $"{p.Start}/{p.Sequence.Length}:{p.Score}"))}]");
                if (!got.Select(Id).SequenceEqual(all.Take(sc.MaxProbes).Select(Id)))
                    return false.Label($"both={both}: top-{sc.MaxProbes} is not a prefix {sc}");
                if (sc.Param.Ranking == ProbeRanking.AdditiveScore && all.Zip(all.Skip(1)).Any(t => t.First.Score < t.Second.Score))
                    return false.Label($"both={both}: not score-descending {sc}");
                byBoth[both] = all.ToDictionary(Id, p => p.Score);
            }
            return byBoth[true].All(kv => byBoth[false][kv.Key] == 0 || kv.Value <= byBoth[false][kv.Key])
                .Label($"a both-strand score exceeds the single-strand score {sc}");
        });
    }
}
