using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c2b (B09 F195): fuzznuc with a range pattern and <c>-pmismatch</c> lists a window once per successful unrolling of the pattern (<c>patOUBrute</c> pushes a hit per unrolling): TestData/Alignment/fuzznuc_dup_oracle.jsonl.gz
/// (scripts/alignment/generate_fuzznuc_dup_oracle.py, 300 patterns, both strands, excel report byte for byte).
/// </summary>
[TestFixture]
public class EmbossFuzznucRangeMismatch_Oracle_Tests
{
    [Test]
    public void FuzznucReport_RangePatternWithMismatches_EqualsTheProgram()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.fuzznuc_dup_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, duplicates = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string expected = j.GetProperty("report").GetString()!;
            string pattern = j.GetProperty("pattern").GetString()!;
            int mm = j.GetProperty("mismatch").GetInt32();
            bool comp = j.GetProperty("complement").GetBoolean();
            if (expected.Split('\n', StringSplitOptions.RemoveEmptyEntries).Skip(1).GroupBy(l => l).Any(g => g.Count() > 1))
                duplicates++;
            string got;
            try
            {
                got = EmbossReportPrograms.Fuzznuc(new EmbossReportRequest(j.GetProperty("seq").GetString()!, ReportQualifiers: "-rformat excel", RunDate: new DateTime(2020, 1, 1, 0, 0, 0, DateTimeKind.Utc)), pattern, mm, comp).Report;
            }
            catch (ArgumentException e)
            {
                failures.Add($"row {rows - 1} {pattern} mm {mm}: {e.Message}");
                continue;
            }

            if (got != expected)
                failures.Add($"row {rows - 1} {pattern} mm {mm} comp {comp}: report differs ({got.Split('\n').Length} vs {expected.Split('\n').Length} lines)");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(duplicates, Is.GreaterThan(80), "rows whose report lists a window more than once");
    }
}
