using System.IO.Compression;
using System.Reflection;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// B09 WP67 (F307, R14-10b): <see cref="PairwiseAligner"/> with gap-score functions of (position, length) (Biopython's Waterman-Smith-Beyer algorithm: <c>insertion_score</c> / <c>deletion_score</c> / <c>gap_score</c>, strands + and -) and with gap scores
/// above zero (or an open score above its extend score: Gotoh, Needleman-Wunsch, Smith-Waterman) equals Bio.Align.PairwiseAligner 1.88 on the rows of TestData/Alignment/pairwise_wsb_oracle.jsonl.gz (scripts/alignment/generate_pairwise_wsb_oracle.py): the
/// name of the algorithm, <c>score()</c>, the score of <c>align()</c> and the SET of alignments (as coordinates; local alignments with gap scores above zero include the ones the reference lists for end points it did not clear, e.g. a lower-scoring alignment ending at an earlier cell), and the RuntimeError the reference raises on one row.
/// </summary>
[TestFixture]
public class PairwiseAligner_GapFunctions_Oracle_Tests
{
    private static Func<int, int, double> Function(JsonElement e) => GapScoreFunctions.Create(e.GetProperty("kind").GetString()!, [.. e.GetProperty("p").EnumerateArray().Select(x => x.GetDouble())]);

    [Test]
    public void Aligner_EqualsBiopython_AllRows()
    {
        const string resource = "Seqeron.Genomics.Tests.TestData.Alignment.pairwise_wsb_oracle.jsonl.gz";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, sets = 0, wsb = 0, positive = 0, errors = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var aligner = new PairwiseAligner
            {
                Mode = j.GetProperty("mode").GetString() == "local" ? PairwiseAlignerMode.Local : PairwiseAlignerMode.Global,
                MatchScore = j.GetProperty("match").GetDouble(),
                MismatchScore = j.GetProperty("mismatch").GetDouble(),
            };
            var cfg = j.GetProperty("cfg");
            bool function = false, positiveGap = false;
            if (cfg.TryGetProperty("gap", out var g)) { aligner.SetGapFunction("gap_score", Function(g)); function = true; positiveGap |= g.GetProperty("kind").GetString() == "bonus"; }
            if (cfg.TryGetProperty("ins", out var ins)) { aligner.SetGapFunction("insertion_score", Function(ins)); function = true; positiveGap |= ins.GetProperty("kind").GetString() == "bonus"; }
            if (cfg.TryGetProperty("del", out var del)) { aligner.SetGapFunction("deletion_score", Function(del)); function = true; positiveGap |= del.GetProperty("kind").GetString() == "bonus"; }
            if (cfg.TryGetProperty("scores", out var sc))
            {
                foreach (var kv in sc.EnumerateObject())
                {
                    aligner.SetGapScore(kv.Name, kv.Value.GetDouble());
                    positiveGap |= kv.Value.GetDouble() > 0;
                }
            }

            if (function)
                wsb++;
            string a = j.GetProperty("a").GetString()!, b = j.GetProperty("b").GetString()!;
            char strand = j.GetProperty("strand").GetString()![0];
            string id = $"row {j.GetProperty("id")} {a} {b} {j.GetProperty("mode")} {strand} {cfg}";
            try
            {
                if (aligner.Algorithm != j.GetProperty("algorithm").GetString())
                {
                    failures.Add($"{id}: algorithm {aligner.Algorithm} expected {j.GetProperty("algorithm").GetString()}");
                    continue;
                }

                double expected = j.GetProperty("score").GetDouble();
                double score = aligner.Score(a, b, strand);
                if (Math.Abs(score - expected) > 1e-9)
                {
                    failures.Add($"{id}: score {score} expected {expected}");
                    continue;
                }

                if (j.TryGetProperty("error", out var reference))
                {
                    errors++;
                    string expectedMessage = reference.GetString()!["RuntimeError: ".Length..];
                    try
                    {
                        aligner.Align(a, b, strand);
                        failures.Add($"{id}: no error, expected {expectedMessage}");
                    }
                    catch (InvalidOperationException e) when (e.Message == expectedMessage)
                    {
                    }

                    continue;
                }

                var alignments = aligner.Align(a, b, strand);
                double aligned = j.GetProperty("aligned_score").GetDouble();
                if (Math.Abs(alignments.Score - aligned) > 1e-9)
                {
                    failures.Add($"{id}: aligned score {alignments.Score} expected {aligned}");
                    continue;
                }

                bool localPositive = positiveGap && aligner.Mode == PairwiseAlignerMode.Local;
                int total = j.GetProperty("total").GetInt32();
                if (total > 400)
                    continue;
                var got = alignments.Select(x => string.Join("|", x.Coordinates.Select(r => string.Join(",", r)))).OrderBy(x => x, StringComparer.Ordinal).ToList();
                var want = j.GetProperty("alignments").EnumerateArray().Select(x => string.Join("|", x.EnumerateArray().Select(r => string.Join(",", r.EnumerateArray().Select(v => v.GetInt64()))))).OrderBy(x => x, StringComparer.Ordinal).ToList();
                if (!got.SequenceEqual(want))
                    failures.Add($"{id}: {got.Count} alignments expected {want.Count}: got [{string.Join(" ; ", got.Take(4))}] expected [{string.Join(" ; ", want.Take(4))}]");
                else
                {
                    sets++;
                    if (localPositive)
                        positive++;
                }
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(1500));
        Assert.That(sets, Is.GreaterThan(1000), "rows whose set of alignments was compared");
        Assert.That(wsb, Is.GreaterThan(600), "rows with a gap-score function");
        Assert.That(positive, Is.GreaterThan(250), "local rows with gap scores above zero whose set of alignments equals the reference");
        Assert.That(errors, Is.EqualTo(1), "rows on which the reference raises RuntimeError");
    }
}
