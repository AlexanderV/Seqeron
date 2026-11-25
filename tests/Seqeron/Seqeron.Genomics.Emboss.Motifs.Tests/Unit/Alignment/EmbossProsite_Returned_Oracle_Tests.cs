// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.SequenceFormats;
using System.IO.Compression;
using System.Reflection;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// B09 WP67 (F306): <c>patmatmotifs</c> of EMBOSS 6.6.0 on the REAL PROSITE release 2026_04 (the 85 runs of the user-run Docker package <c>handoff/WP62-patmatmotifs-prosite</c>: <c>prosextract</c> of the SIB files, <c>patmatmotifs -auto ARGS</c> on 6 input
/// files x -full x -prune, 22 report formats, the report qualifiers, ranges, damaged / empty / missing input; every run also under valgrind: no run reads invalid memory). The data directory is the one of the runs (<c>/usr/share/EMBOSS/data/PROSITE</c>:
/// the files <c>prosite.lines</c> of the entries the reports name and the documentation files of the <c>-full</c> runs, TestData/Alignment/emboss_prosite_release_data.json.gz); the port equals the program on the standard output, the error stream, the
/// exit status and every file written, on every run.
/// </summary>
[TestFixture]
public class EmbossProsite_Returned_Oracle_Tests
{
    private static JsonDocument Load(string name)
    {
        string resource = $"Seqeron.Genomics.Tests.TestData.Alignment.{name}";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        return JsonDocument.Parse(new StreamReader(gz).ReadToEnd());
    }

    [Test]
    public void PatmatmotifsOnTheRealPrositeRelease_EqualsTheProgram_AllRuns()
    {
        var data = new Dictionary<string, string>();
        using (var d = Load("emboss_prosite_release_data.json.gz"))
            foreach (var p in d.RootElement.EnumerateObject())
                data[p.Name] = p.Value.GetString()!;
        string resource = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_prosite_returned_oracle.jsonl.gz";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, hits = 0, died = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            var files = new Dictionary<string, string>(data);
            foreach (var p in j.GetProperty("files").EnumerateObject())
                files[p.Name] = p.Value.GetString()!;
            int rc = j.GetProperty("rc").GetInt32();
            // A row without an epoch was recorded on a real day (the GFF header prints "#!Date yyyy-MM-dd"): the clock of the port is set to that day, not to today (the test was date dependent: audit WP79).
            DateTime? recorded = null;
            foreach (var of in j.GetProperty("outfiles").EnumerateObject())
            {
                var dm = System.Text.RegularExpressions.Regex.Match(of.Value.GetString() ?? "", @"#!Date (\d{4})-(\d{2})-(\d{2})");
                if (dm.Success)
                    recorded = new DateTime(int.Parse(dm.Groups[1].Value), int.Parse(dm.Groups[2].Value), int.Parse(dm.Groups[3].Value), 12, 0, 0, DateTimeKind.Utc);
            }

            var clock = new EmbossRunEnvironment(j.TryGetProperty("epoch", out var epoch) ? DateTimeOffset.FromUnixTimeSeconds(epoch.GetInt64()).UtcDateTime : recorded);
            string id = $"row {j.GetProperty("id")} [{string.Join(' ', args)}]";
            var r = EmbossMotifProgramRunner.Run("patmatmotifs", args, files, clock);
            var expectedFiles = new SortedDictionary<string, string>(StringComparer.Ordinal);
            foreach (var p in j.GetProperty("outfiles").EnumerateObject())
                expectedFiles[p.Name] = p.Value.GetString()!;
            string why = "";
            if (r.ExitCode != rc)
                why = $"exit status {r.ExitCode} expected {rc}";
            else if (r.Stdout != j.GetProperty("out").GetString())
                why = "standard output differs";
            else if (r.Stderr != j.GetProperty("err").GetString())
                why = $"error stream '{r.Stderr.Replace("\n", "|")}' expected '{j.GetProperty("err").GetString()!.Replace("\n", "|")}'";
            else if (!expectedFiles.SequenceEqual(new SortedDictionary<string, string>(r.Files.ToDictionary(k => k.Key, k => k.Value), StringComparer.Ordinal)))
                why = $"files {string.Join(",", r.Files.Keys)} expected {string.Join(",", expectedFiles.Keys)} (or their text differs)";
            if (why.Length > 0)
                failures.Add($"{id}: {why}");
            else
            {
                if ((r.Stdout + string.Concat(r.Files.Values)).Contains("HitCount: ", StringComparison.Ordinal) && !(r.Stdout + string.Concat(r.Files.Values)).Contains("HitCount: 0", StringComparison.Ordinal))
                    hits++;
                if (rc != 0)
                    died++;
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.EqualTo(85));
        Assert.That(hits, Is.GreaterThan(20), "runs with hits");
        Assert.That(died, Is.EqualTo(6), "the runs that end with an error of the library");
    }
}
