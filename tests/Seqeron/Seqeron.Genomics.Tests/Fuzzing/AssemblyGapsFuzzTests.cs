using System.Text;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for ASSEMBLY-GAPS-001 (GenomeAssemblyAnalyzer FindGaps / AnalyzeGapDistribution /
/// AnalyzeScaffolds / ExtractContigs). Random N-run-rich sequences are checked against an independent
/// regex oracle of QUAST <c>split_by_ns</c> (cut at maximal runs of ≥ minGapLength N/n) and of
/// <c>seqtk gap</c> (maximal N-runs), plus the invariants ContigLength + GapLength = TotalLength,
/// contigs/gaps tile the scaffold, and ExtractContigs = AnalyzeScaffolds contigs filtered by length.
/// The Python transcription (tools/review/oracles/B10/gaps_oracle.py) agreed 3000/3000 (Evidence).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyGapsFuzzTests
{
    [Test]
    public void GapsScaffoldsContigs_RandomSequences_MatchRegexOracle()
    {
        var rng = new Random(2026);
        int[] runs = { 1, 1, 2, 3, 5, 9, 10, 11, 15, 25 };
        for (int trial = 0; trial < 2000; trial++)
        {
            int length = rng.Next(0, 61);
            var sb = new StringBuilder();
            while (sb.Length < length)
            {
                if (rng.NextDouble() < 0.3)
                    for (int r = runs[rng.Next(runs.Length)]; r > 0; r--) sb.Append(rng.Next(2) == 0 ? 'N' : 'n');
                else
                    sb.Append("ACGTacgtRY"[rng.Next(10)]);
            }
            string seq = sb.ToString();
            int minGap = new[] { 1, 2, 3, 5, 10, 11, 20 }[rng.Next(7)];
            int minContig = new[] { 0, 1, 2, 3, 5, 8, 20 }[rng.Next(7)];
            var input = new[] { ("s", seq) };
            // A3-2: QUAST correct_seq reads every non-ACGT symbol (R, Y here) as N; false = only N/n.
            bool asGap = rng.Next(2) == 0;
            string gapClass = asGap ? "[^ACGTacgt]" : "[Nn]";
            string ctx = $"seq={seq} minGap={minGap} minContig={minContig} nonAcgtAsGap={asGap}";

            var expectedGaps = Regex.Matches(seq, gapClass + "+").Where(m => m.Length >= minGap)
                .Select(m => (m.Index, m.Index + m.Length - 1, m.Length)).ToList();
            var gaps = GenomeAssemblyAnalyzer.FindGaps(input, minGap, asGap).Select(g => (g.Start, g.End, g.Length)).ToList();
            Assert.That(gaps, Is.EqualTo(expectedGaps), ctx);

            var cuts = Regex.Matches(seq, $"{gapClass}{{{minGap},}}").ToList();
            var expectedContigs = new List<(int, int)>();
            int pos = 0;
            foreach (var m in cuts)
            {
                if (m.Index > pos) expectedContigs.Add((pos, m.Index - 1));
                pos = m.Index + m.Length;
            }
            if (seq.Length > pos) expectedContigs.Add((pos, seq.Length - 1));

            var st = GenomeAssemblyAnalyzer.AnalyzeScaffolds(input, minGap, asGap).Single();
            Assert.That(st.Contigs.Select(c => (c.Start, c.End)), Is.EqualTo(expectedContigs), ctx);
            Assert.That(st.Gaps.Select(g => (g.Start, g.End, g.Length)), Is.EqualTo(expectedGaps), ctx);
            Assert.That(st.ContigLength + st.GapLength, Is.EqualTo(seq.Length), ctx);

            var extracted = GenomeAssemblyAnalyzer.ExtractContigs(input, minContig, minGap, asGap).Select(c => c.Sequence);
            Assert.That(extracted, Is.EqualTo(expectedContigs.Select(c => seq[c.Item1..(c.Item2 + 1)])
                .Where(c => c.Length >= minContig)), ctx);

            var d = GenomeAssemblyAnalyzer.AnalyzeGapDistribution(gaps.Select(g => new GenomeAssemblyAnalyzer.GapInfo("s", g.Start, g.End, g.Length, "x")));
            if (gaps.Count > 0)
            {
                var sorted = gaps.Select(g => (double)g.Length).OrderBy(x => x).ToArray();
                double median = sorted.Length % 2 == 1 ? sorted[sorted.Length / 2]
                    : (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]) / 2;
                Assert.That(d.MedianLength, Is.EqualTo(median), ctx);
            }
        }
    }
}
