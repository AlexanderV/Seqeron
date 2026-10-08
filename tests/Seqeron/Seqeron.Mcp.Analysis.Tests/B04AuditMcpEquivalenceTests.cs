using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Heavy-tier wrapper equivalence for the MCP tools added or extended by the B04 completeness audit (B04 F49 and the
/// additive options of F34–F65): on fixed-seed random inputs every wrapper returns exactly what its library call returns,
/// and the cross-tool identities hold — lempel_ziv_complexity.normalized = compression_ratio, the masked positions of
/// mask_low_complexity = find_low_complexity_intervals, and mask_approximate_tandem_repeats masks exactly the union of the
/// find_approximate_tandem_repeats spans. Expected values are the library's (themselves reference-locked in
/// Seqeron.Genomics.Tests), never re-derived from the wrapper.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class B04AuditMcpEquivalenceTests
{
    private const int Seed = 4004;

    private static string Segmented(Random rng, int length, string extra)
    {
        var sb = new System.Text.StringBuilder();
        while (sb.Length < length)
        {
            int kind = rng.Next(extra.Length > 0 ? 3 : 2);
            if (kind == 0)
                for (int i = rng.Next(2, 30); i > 0; i--) sb.Append("ACGT"[rng.Next(4)]);
            else if (kind == 2)
                for (int i = rng.Next(1, 8); i > 0; i--) sb.Append(extra[rng.Next(extra.Length)]);
            else
            {
                string unit = new(Enumerable.Range(0, rng.Next(1, 7)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
                for (int c = rng.Next(3, 12); c > 0; c--)
                    foreach (char ch in unit) sb.Append(rng.Next(20) == 0 ? "ACGT"[rng.Next(4)] : ch);
            }
        }
        return sb.ToString(0, length);
    }

    private static IEnumerable<int> Positions(IEnumerable<(int Start, int End)> intervals) =>
        intervals.SelectMany(iv => Enumerable.Range(iv.Start, iv.End - iv.Start));

    [Test]
    public void DustLongdustAndLzTools_EqualLibraryCalls()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 150; trial++)
        {
            string s = Segmented(rng, rng.Next(1, 300), "NNRY");
            int w = rng.Next(8, 65), level = rng.Next(2, 65), linker = rng.Next(1, 33);
            double t = level / 10.0;
            foreach (string engine in new[] { "sdust", "dustmasker" })
            {
                var dustEngine = SequenceComplexity.ParseDustEngine(engine);
                var lib = SequenceComplexity.FindLowComplexityIntervals(s, w, t, linker, dustEngine);
                var tool = AnalysisTools.FindLowComplexityIntervals(s, w, t, linker, engine).Items;
                Assert.That(tool.Select(x => (x.Start, x.End)), Is.EqualTo(lib), s);
                Assert.That(tool.All(x => x.Length == x.End - x.Start), Is.True);

                bool soft = trial % 2 == 0;
                string masked = AnalysisTools.MaskLowComplexity(s, w, t, 'X', linker, soft, engine).Masked;
                Assert.That(masked, Is.EqualTo(SequenceComplexity.MaskLowComplexity(s, w, t, 'X', linker, soft, dustEngine)), s);
                var maskedPositions = Enumerable.Range(0, s.Length).Where(i => soft ? char.IsLower(masked[i]) : masked[i] == 'X');
                Assert.That(maskedPositions, Is.EqualTo(Positions(lib)), s);
            }

            int k = rng.Next(1, 9), lw = rng.Next(10, 400);
            Assert.That(AnalysisTools.LongdustScore(s, k).Score, Is.EqualTo(SequenceComplexity.CalculateLongdustScore(s, k)));
            Assert.That(AnalysisTools.FindLongdustRegions(s, k, lw).Items.Select(x => (x.Start, x.End)),
                Is.EqualTo(SequenceComplexity.FindLongdustRegions(s, k, lw)), s);

            var lz = AnalysisTools.LempelZivComplexity(s);
            Assert.That(lz.Complexity, Is.EqualTo(SequenceComplexity.CalculateLempelZivComplexity(s)));
            Assert.That(lz.Normalized, Is.EqualTo(AnalysisTools.CompressionRatio(s).Ratio), s);
        }
    }

    [Test]
    public void LowComplexityRegionTools_EqualLibraryCalls()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 150; trial++)
        {
            string s = Segmented(rng, rng.Next(1, 300), "N");
            int w = rng.Next(6, 64), k = rng.Next(1, 6);
            double t = rng.Next(1, 20) / 10.0, cutoff = rng.Next(0, 101) / 100.0;
            Assert.That(AnalysisTools.FindLowComplexityRegions(s, w, t).Items.Select(x => (x.Start, x.End, x.Sequence)),
                Is.EqualTo(SequenceComplexity.FindLowComplexityRegions(s, w, t).Select(x => (x.Start, x.End, x.Sequence))), s);
            Assert.That(AnalysisTools.FindLowComplexityRegions(s, w, cutoff, "bbduk", k).Items.Select(x => (x.Start, x.End, x.MinEntropy)),
                Is.EqualTo(SequenceComplexity.FindLowEntropyRegionsBbduk(s, cutoff, w, k).Select(x => (x.Start, x.End, x.MinEntropy))), s);
        }
    }

    [Test]
    public void TandemRepeatTools_EqualLibraryCalls_AndMaskIsUnionOfSpans()
    {
        var rng = new Random(Seed + 2);
        for (int trial = 0; trial < 60; trial++)
        {
            string s = Segmented(rng, rng.Next(10, 400), trial % 2 == 0 ? "N" : "");
            int pm = trial % 3 == 0 ? 75 : 80, flank = trial % 4 == 0 ? 50 : 0;
            var p = TandemRepeatsFinderParameters.Recommended with { MatchProbability = pm, FlankLength = flank };
            var lib = RepeatFinder.FindApproximateTandemRepeats(s, p, out int count);
            var tool = AnalysisTools.FindApproximateTandemRepeats(s, matchProbability: pm, flankLength: flank).Items;
            Assert.That(tool.Select(x => (x.Start, x.SpanLength, x.Period, x.Consensus, x.AlignmentScore, x.EntropyTrf, x.AlignedSequence,
                    x.LeftFlank, x.CopyMatches, x.OutputIndex, x.OutputCount, x.DetectionPosition)),
                Is.EqualTo(lib.Select(x => (x.Start, x.SpanLength, x.Period, x.Consensus, x.AlignmentScore, x.EntropyTrf, x.AlignedSequence,
                    x.LeftFlank, x.CopyMatches, x.OutputIndex, x.OutputCount, x.DetectionPosition))), s);
            Assert.That(lib.All(r => r.OutputCount == count), Is.True);

            string masked = AnalysisTools.MaskApproximateTandemRepeats(s, matchProbability: pm).Masked;
            var spans = RepeatFinder.FindApproximateTandemRepeats(s, p with { FlankLength = 0 })
                .SelectMany(r => Enumerable.Range(r.Start, r.SpanLength)).ToHashSet();
            Assert.That(Enumerable.Range(0, s.Length).Where(i => masked[i] == 'N' && s[i] != 'N'),
                Is.EqualTo(Enumerable.Range(0, s.Length).Where(i => spans.Contains(i) && s[i] != 'N')), s);
            Assert.That(Enumerable.Range(0, s.Length).All(i => spans.Contains(i) || masked[i] == s[i]), Is.True, s);

            var dat = AnalysisTools.FindApproximateTandemRepeats(s, matchProbability: pm, flankLength: flank, format: "dat", sequenceName: "q").Formatted;
            Assert.That(dat, Is.EqualTo(RepeatFinder.FormatTrfDatFileHeader() + RepeatFinder.FormatTrfDatLines(s, lib, "q", p)), s);
            var html = AnalysisTools.FindApproximateTandemRepeats(s, matchProbability: pm, flankLength: flank, format: "html", sequenceName: "q");
            Assert.That(html.AlignmentPages!.Select(x => (x.FileName, x.Html)),
                Is.EqualTo(RepeatFinder.FormatTrfAlignmentPages(s, lib, "q", p, "q", count).Select(x => (x.FileName, x.Html))), s);
        }
    }

    [Test]
    public void DegenerateAndMicrosatelliteTools_EqualLibraryCalls()
    {
        var rng = new Random(Seed + 3);
        string[] definitions = ["1-10 2-6 3-5 4-5 5-5 6-5", "1-3 2-2 3-2 4-2 5-2 6-2", "2-3 7-2 9-2 10-2"];
        for (int trial = 0; trial < 120; trial++)
        {
            string s = Segmented(rng, rng.Next(10, 120), "N");
            int k = rng.Next(1, 3), minLength = k + rng.Next(2, 9);
            string distance = trial % 2 == 0 ? "edit" : "hamming";
            bool pal = trial % 4 >= 2;
            string reporting = trial % 3 == 0 ? "bestPerSeed" : "allMaximal";
            bool compat = trial % 5 == 0;
            var lib = RepeatFinder.FindDegenerateRepeats(s, minLength, k,
                distance == "edit" ? ApproximateRepeatDistance.Edit : ApproximateRepeatDistance.Hamming, pal, int.MaxValue, int.MinValue,
                reporting == "bestPerSeed" ? DegenerateRepeatReporting.BestPerSeed : DegenerateRepeatReporting.AllMaximal, compat);
            var tool = AnalysisTools.FindDegenerateRepeats(s, minLength, k, distance, pal, int.MaxValue, int.MinValue, reporting, compat).Items;
            Assert.That(tool.Select(x => (x.FirstPosition, x.FirstLength, x.SecondPosition, x.SecondLength, x.Distance, x.FirstCopy, x.SecondCopy)),
                Is.EqualTo(lib.Select(x => (x.FirstPosition, x.FirstLength, x.SecondPosition, x.SecondLength, x.Distance, x.FirstCopy, x.SecondCopy))), s);

            string def = definitions[trial % definitions.Length];
            var map = RepeatFinder.ParseMisaDefinition(def);
            foreach (bool misaScan in new[] { false, true })
            {
                var mode = misaScan ? MicrosatelliteScanMode.MisaRegex : MicrosatelliteScanMode.MaximalRuns;
                Assert.That(AnalysisTools.FindMicrosatellites(s, misaScan: misaScan, misaDefinition: def).Items.Select(x => (x.Position, x.RepeatUnit, x.RepeatCount)),
                    Is.EqualTo(RepeatFinder.FindMicrosatellites(s, map, mode).Select(x => (x.Position, x.RepeatUnit, x.RepeatCount))), s);
                var summary = AnalysisTools.TandemRepeatSummary(s, misaScan: misaScan, misaDefinition: def);
                var libSummary = RepeatFinder.GetTandemRepeatSummary(s, map, mode);
                Assert.That((summary.TotalRepeats, summary.TotalRepeatBases, summary.PercentageOfSequence),
                    Is.EqualTo((libSummary.TotalRepeats, libSummary.TotalRepeatBases, libSummary.PercentageOfSequence)), s);
                Assert.That(summary.CountsByUnitLength, Is.EquivalentTo(libSummary.CountsByUnitLength), s);
            }
        }
    }
}
