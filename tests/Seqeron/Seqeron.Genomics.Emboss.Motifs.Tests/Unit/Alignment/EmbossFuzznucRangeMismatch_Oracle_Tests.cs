// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.FeatureTools;
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
                got = EmbossMotifReportPrograms.Fuzznuc(new EmbossReportRequest(j.GetProperty("seq").GetString()!, ReportQualifiers: "-rformat excel", RunDate: new DateTime(2020, 1, 1, 0, 0, 0, DateTimeKind.Utc)), pattern, mm, comp).Report;
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

    /// <summary>
    /// A5-8c3b (F217): patterns whose elements are ALL optional ranges (<c>C(0,2)C(0,1)</c>, <c>X(0,2)[AC](0,1)</c>) of fuzznuc, fuzzpro and fuzztran, with and without <c>-pmismatch</c>, anchors <c>&lt;</c> / <c>&gt;</c>, both strands: the programs also list
    /// the EMPTY hits (<c>end = start - 1</c>, score 0) at every start of the sequence, and a pattern anchored at the end the empty match at the end when no match began at the last residue; the translation of an empty fuzztran hit is the
    /// <c>ajStrAssignSubS</c> clipping of its positions (TestData/Alignment/fuzznuc_optional_oracle.jsonl.gz, scripts/alignment/generate_fuzznuc_optional_oracle.py: 400 reports, whole excel files).
    /// </summary>
    [Test]
    public void AllOptionalRangePatterns_ListTheEmptyHits_EqualsTheProgram()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.fuzznuc_optional_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, empty = 0, anchoredEnd = 0, withMismatch = 0;
        var perProgram = new Dictionary<string, int>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string program = j.GetProperty("program").GetString()!;
            perProgram[program] = perProgram.GetValueOrDefault(program) + 1;
            string expected = j.GetProperty("report").GetString()!;
            string pattern = j.GetProperty("pattern").GetString()!;
            int mm = j.GetProperty("mismatch").GetInt32();
            var request = new EmbossReportRequest(j.GetProperty("seq").GetString()!, ReportQualifiers: "-rformat excel", RunDate: new DateTime(2020, 1, 1, 0, 0, 0, DateTimeKind.Utc));
            string got = program switch
            {
                "fuzznuc" => EmbossMotifReportPrograms.Fuzznuc(request, pattern, mm, j.GetProperty("complement").GetBoolean()).Report,
                "fuzzpro" => EmbossMotifReportPrograms.Fuzzpro(request, pattern, null, mm).Report,
                _ => EmbossMotifReportPrograms.Fuzztran(request, pattern, j.GetProperty("frame").GetString()!, 0, null, mm).Report,
            };
            if (got != expected)
                failures.Add($"row {rows - 1} {program} {pattern} mm {mm}: report differs ({got.Split('\n').Length} vs {expected.Split('\n').Length} lines)");
            if (expected.Split('\n').Any(l => l.Split('\t').Length > 3 && l.Split('\t')[2] != "End" && int.TryParse(l.Split('\t')[1], out int a) && int.TryParse(l.Split('\t')[2], out int b) && b == a - 1))
                empty++;
            if (pattern.EndsWith('>'))
                anchoredEnd++;
            if (mm > 0)
                withMismatch++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(400));
        Assert.That(perProgram["fuzznuc"], Is.GreaterThanOrEqualTo(180));
        Assert.That(perProgram["fuzzpro"], Is.GreaterThanOrEqualTo(90));
        Assert.That(perProgram["fuzztran"], Is.GreaterThanOrEqualTo(90));
        Assert.That(empty, Is.GreaterThan(250), "reports with an empty hit");
        Assert.That(anchoredEnd, Is.GreaterThan(30));
        Assert.That(withMismatch, Is.GreaterThan(100));
    }
}
