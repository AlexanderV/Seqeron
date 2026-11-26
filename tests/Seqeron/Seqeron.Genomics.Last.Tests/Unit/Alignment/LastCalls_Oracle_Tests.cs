// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text.Json;
using Seqeron.Genomics.Last;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-last, B09 F376): every call of <c>lastdb</c> / <c>lastal</c> that the 204 runs of the user's Docker package made (LAST 1256 built from the sources, behind logging wrappers; <c>last_calls_oracle.tsv.gz</c>):
/// the port gives the same standard output (the raw text, byte for byte: the header with the derived parameters and the scoring matrix with its ambiguity scores, every alignment line), the same standard error and the same status.
/// </summary>
[TestFixture]
public class LastCalls_Oracle_Tests
{
    internal sealed record Call(string Set, int Index, string Program, string[] Arguments, Dictionary<string, string> Files, string StandardOutput, string StandardError, int ExitCode);

    internal static IReadOnlyList<Call> Load()
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.last_calls_oracle.tsv")!;
        using var reader = new StreamReader(stream);
        var list = new List<Call>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            string[] f = line.Split('\t');
            list.Add(new Call(f[0], int.Parse(f[1]), f[2], JsonSerializer.Deserialize<string[]>(OracleEscapes.Unescape(f[3]))!,
                JsonSerializer.Deserialize<Dictionary<string, string>>(OracleEscapes.Unescape(f[4]))!, OracleEscapes.Unescape(f[5]), OracleEscapes.Unescape(f[6]), int.Parse(f[7])));
        }

        return list;
    }

    private static readonly Lazy<IReadOnlyList<Call>> Calls = new(Load);

    private static LastProcessResult RunCall(Call c)
    {
        var dir = new Dictionary<string, string>(c.Files);
        var last = new LastPrograms();
        if (c.Program == "lastal")
        {
            // the database of the call was made by the preceding lastdb call of the run: <name> <name>
            string db = c.Arguments[^2];
            var made = last.Run("/opt/last/bin/lastdb", [db, db], dir);
            Assert.That(made.ExitCode, Is.Zero, c.Set + " lastdb " + db);
        }

        return last.Run("/opt/last/bin/" + c.Program, c.Arguments, dir);
    }

    private static string FirstDifference(Call c, LastProcessResult r)
    {
        string[] a = c.StandardOutput.Split('\n'), b = r.StandardOutput.Split('\n');
        for (int i = 0; i < Math.Max(a.Length, b.Length); i++)
        {
            string x = i < a.Length ? a[i] : "<none>", y = i < b.Length ? b[i] : "<none>";
            if (x != y)
                return $"{c.Set}#{c.Index} line {i}\n expected: {x}\n actual:   {y}\n rc {r.ExitCode}/{c.ExitCode} stderr [{r.StandardError}]";
        }

        return $"{c.Set}#{c.Index} rc {r.ExitCode}/{c.ExitCode} stderr [{r.StandardError}] expected [{c.StandardError}]";
    }

    [Test]
    public void Oracle_IsComplete()
    {
        Assert.That(Calls.Value.Count(c => c.Program == "lastal"), Is.GreaterThan(700));
        Assert.That(Calls.Value.Count(c => c.Program == "lastdb"), Is.GreaterThan(800));
        Assert.That(Calls.Value.Select(c => c.Set).Distinct().Count(), Is.GreaterThan(100));
    }

    [Test]
    public void Lastdb_EveryCall_WritesNothingAndSucceeds()
    {
        var bad = new List<string>();
        foreach (var c in Calls.Value.Where(c => c.Program == "lastdb"))
        {
            var r = new LastPrograms().Run("/opt/last/bin/lastdb", c.Arguments, new Dictionary<string, string>(c.Files));
            if (r.ExitCode != c.ExitCode || r.StandardOutput != c.StandardOutput || r.StandardError != c.StandardError)
                bad.Add($"{c.Set}#{c.Index}");
        }

        Assert.That(bad, Is.Empty);
    }

    [Test]
    public void Lastal_EveryCall_EqualsTheRawOutput()
    {
        var bad = new List<string>();
        string first = "";
        foreach (var c in Calls.Value.Where(c => c.Program == "lastal"))
        {
            var r = RunCall(c);
            if (r.ExitCode != c.ExitCode || r.StandardOutput != c.StandardOutput || r.StandardError != c.StandardError)
            {
                bad.Add($"{c.Set}#{c.Index}");
                if (first.Length == 0)
                    first = FirstDifference(c, r);
            }
        }

        Assert.That(bad, Is.Empty, bad.Count + " of " + Calls.Value.Count(c => c.Program == "lastal") + " calls differ; first: " + first);
    }
}
