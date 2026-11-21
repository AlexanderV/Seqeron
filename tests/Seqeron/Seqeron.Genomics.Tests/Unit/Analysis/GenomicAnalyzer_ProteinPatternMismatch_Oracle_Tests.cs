// GENOMIC-MOTIFS-001 (B09 A5-8c2, F195) — EMBOSS fuzzpro / fuzztran with -pmismatch n: GenomicAnalyzer.FindProteinPatternHits / FindTranslatedPatternHits(maxMismatches).
// Oracle: scripts/analysis/generate_fuzzpro_mismatch_oracle.py (1 000 random syntax-tree patterns, -pmismatch 1 - 3, fuzzpro 750 + fuzztran 250 rows, 6.6.0).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_ProteinPatternMismatch_Oracle_Tests
{
    private static List<string[]> Rows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.fuzzpro_mismatch_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0) rows.Add(line.Split('\t'));
        return rows;
    }

    [Test]
    public void FindProteinPatternHits_WithMismatches_EqualsFuzzpro_AllRows()
    {
        var rows = Rows().Where(r => r[0] == "pro").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(700));
        var failures = new List<string>();
        int hitRows = 0, withMismatch = 0, rejected = 0, ranges = 0;
        foreach (var f in rows)
        {
            string label = $"{f[1]} / {f[2]} mm {f[3]}";
            int mm = int.Parse(f[3]);
            if (f[4] == "ERR")
            {
                try
                {
                    GenomicAnalyzer.FindProteinPatternHits(f[1], f[2], default, mm);
                    failures.Add("accepted: " + label);
                }
                catch (ArgumentException)
                {
                    rejected++;
                }

                continue;
            }

            string[] hits;
            try
            {
                hits = [.. GenomicAnalyzer.FindProteinPatternHits(f[1], f[2], default, mm).Select(h => $"{h.Position}:{h.Length}:{h.Mismatches}").OrderBy(x => x, StringComparer.Ordinal)];
            }
            catch (ArgumentException e)
            {
                failures.Add($"rejected ({e.Message}): {label}");
                continue;
            }

            var expected = f[4] == "-" ? [] : f[4].Split(' ').OrderBy(x => x, StringComparer.Ordinal).ToArray();
            if (!hits.SequenceEqual(expected))
                failures.Add($"{label}: {hits.Length} vs {expected.Length}; first difference at {Enumerable.Range(0, Math.Min(hits.Length, expected.Length)).FirstOrDefault(i => hits[i] != expected[i], -1)}: ours {string.Join(' ', hits.Skip(Math.Max(0, Enumerable.Range(0, Math.Min(hits.Length, expected.Length)).FirstOrDefault(i => hits[i] != expected[i], 0) - 1)).Take(4))} expected {string.Join(' ', expected.Skip(Math.Max(0, Enumerable.Range(0, Math.Min(hits.Length, expected.Length)).FirstOrDefault(i => hits[i] != expected[i], 0) - 1)).Take(4))}");
            if (expected.Length > 0) hitRows++;
            if (expected.Any(e => !e.EndsWith(":0"))) withMismatch++;
            if (f[2].Contains(',')) ranges++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(hitRows, Is.GreaterThan(300));
        Assert.That(withMismatch, Is.GreaterThan(250));
        Assert.That(rejected, Is.GreaterThan(30));
        Assert.That(ranges, Is.GreaterThan(200));
    }

    [Test]
    public void FindTranslatedPatternHits_WithMismatches_EqualsFuzztran_AllRows()
    {
        var rows = Rows().Where(r => r[0] == "tran").ToList();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(200));
        var failures = new List<string>();
        int hitRows = 0, withMismatch = 0;
        foreach (var f in rows)
        {
            string label = $"{f[1]} / {f[2]} mm {f[3]}";
            int mm = int.Parse(f[3]);
            if (f[4] == "ERR")
            {
                try
                {
                    GenomicAnalyzer.FindTranslatedPatternHits(f[1], f[2], "6", 0, default, mm);
                    failures.Add("accepted: " + label);
                }
                catch (ArgumentException)
                {
                }

                continue;
            }

            string[] hits;
            try
            {
                hits = [.. GenomicAnalyzer.FindTranslatedPatternHits(f[1], f[2], "6", 0, default, mm).Select(h => $"{h.ProteinStart}:{h.ProteinEnd}:{h.Mismatches}@{h.Frame}").OrderBy(x => x, StringComparer.Ordinal)];
            }
            catch (ArgumentException e)
            {
                failures.Add($"rejected ({e.Message}): {label}");
                continue;
            }

            var expected = f[4] == "-" ? [] : f[4].Split(' ').OrderBy(x => x, StringComparer.Ordinal).ToArray();
            if (!hits.SequenceEqual(expected))
                failures.Add($"{label}: {hits.Length} vs {expected.Length}; first difference at {Enumerable.Range(0, Math.Min(hits.Length, expected.Length)).FirstOrDefault(i => hits[i] != expected[i], -1)}: ours {string.Join(' ', hits.Skip(Math.Max(0, Enumerable.Range(0, Math.Min(hits.Length, expected.Length)).FirstOrDefault(i => hits[i] != expected[i], 0) - 1)).Take(4))} expected {string.Join(' ', expected.Skip(Math.Max(0, Enumerable.Range(0, Math.Min(hits.Length, expected.Length)).FirstOrDefault(i => hits[i] != expected[i], 0) - 1)).Take(4))}");
            if (expected.Length > 0) hitRows++;
            if (expected.Any(e => !e.Contains(":0@"))) withMismatch++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(hitRows, Is.GreaterThan(80));
        Assert.That(withMismatch, Is.GreaterThan(60));
    }
}
