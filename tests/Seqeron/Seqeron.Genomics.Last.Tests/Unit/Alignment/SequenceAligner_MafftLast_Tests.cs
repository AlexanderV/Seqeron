// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Last;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-last, B09 WP85, F376): <c>mafft --lastpair</c> / <c>--lastmultipair</c> of the script <c>mafft</c> (MAFFT 7.505-1) with LAST 1256 behind it (the user's Docker package, run 2, 204 sets,
/// <c>scripts/alignment/convert_last_returned.py</c>, <c>mafft_last_oracle.tsv.gz</c>): return code, standard output, files and standard error of every run equal the port with its own LAST (<see cref="LastPrograms"/>), except the
/// runs in which <c>pairlocalalign -R</c> overflows its heap (valgrind proof in the package): the port refuses them with <see cref="NotSupportedException"/> instead of guessing the place of the crash.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftLast_Tests
{
    internal sealed record Case(string Id, string Group, string[] Args, string Input, Dictionary<string, string> Files, int Rc, string Out, string ErrSha, string Err, Dictionary<string, string> OutFiles);

    private static Dictionary<string, string> JsonMap(string json)
    {
        using var doc = JsonDocument.Parse(json);
        return doc.RootElement.EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
    }

    internal static List<Case> Load()
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.mafft_last_oracle.tsv") ?? throw new InvalidOperationException("oracle resource missing");
        using var reader = new StreamReader(stream, Encoding.UTF8);
        var list = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Case(f[0], f[1], JsonSerializer.Deserialize<string[]>(OracleEscapes.Unescape(f[2]))!, OracleEscapes.Unescape(f[3]), JsonMap(OracleEscapes.Unescape(f[4])), int.Parse(f[5], CultureInfo.InvariantCulture),
                OracleEscapes.Unescape(f[6]), f[7], OracleEscapes.Unescape(f[8]), JsonMap(OracleEscapes.Unescape(f[9]))));
        }

        return list;
    }

    private static readonly Lazy<List<Case>> Oracle = new(Load);

    private static readonly MafftInstallation WithLast = MafftInstallation.Debian with { Last = new LastPrograms() };

    private const string InputName = "in.fa";

    private static MafftProgramResult Run(Case c, MafftInstallation installation)
    {
        var files = new Dictionary<string, string>(c.Files) { [InputName] = c.Input };
        return MafftProgram.Run([.. c.Args, InputName], installation, files);
    }

    /// <summary>The standard error as the converter normalised it for the groups with threads: the number of the thread that prints a progress line is a race of the program.</summary>
    private static string NormalizeThreads(string t) => Regex.Replace(Regex.Replace(t, @"by thread +\d+", "by thread T"), @"\((thread|threads) +\d+", "($1 T");

    public static IEnumerable<TestCaseData> Groups() => new[] { "last-plain", "last-plain-thread", "last-add", "last-protein", "last-stop" }.Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        var by = Oracle.Value.GroupBy(r => r.Group).ToDictionary(g => g.Key, g => g.Count());
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value, Has.Count.EqualTo(204));
            Assert.That(by.Keys, Is.EquivalentTo(new[] { "last-plain", "last-plain-thread", "last-crash", "last-add", "last-protein", "last-stop" }));
            Assert.That(by["last-crash"], Is.EqualTo(100), "the runs in which pairlocalalign -R overflows its heap");
            Assert.That(by["last-plain"] + by["last-plain-thread"], Is.EqualTo(35), "the runs that reach the end");
            Assert.That(Oracle.Value.Where(r => r.Group == "last-crash").Select(r => r.Rc), Is.All.EqualTo(1));
        });
    }

    [TestCaseSource(nameof(Groups))]
    public void ReturnCodeOutputAndStandardError_EqualMafft(string group)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var c in Oracle.Value.Where(x => x.Group == group))
        {
            n++;
            string what = c.Id + " (" + string.Join(' ', c.Args) + ")";
            MafftProgramResult r;
            try
            {
                r = Run(c, WithLast);
            }
            catch (Exception e)
            {
                failures.Add(what + ": " + e.GetType().Name + " " + e.Message);
                continue;
            }

            string err = group.EndsWith("-thread", StringComparison.Ordinal) ? NormalizeThreads(r.Stderr) : r.Stderr;
            string sha = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(err))).ToLowerInvariant();
            if (r.ExitCode != c.Rc)
                failures.Add(what + ": return code " + r.ExitCode + " != " + c.Rc);
            else if (r.Stdout != c.Out)
                failures.Add(what + ": standard output differs");
            else if (!r.Files.OrderBy(k => k.Key, StringComparer.Ordinal).SequenceEqual(c.OutFiles.OrderBy(k => k.Key, StringComparer.Ordinal)))
                failures.Add(what + ": files differ");
            else if (sha != c.ErrSha || (c.Err.Length > 0 && c.Err != err))
                failures.Add(what + ": standard error differs: " + err.Replace("\n", "|").Replace("\r", "~"));
        }

        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    /// <summary>
    /// The runs in which the real program crashes: <c>putlocalhom_last</c> writes behind the row i of the table of local homologies for a pair (i &gt;= 1, j &gt;= n - i) that has an alignment. The place of the crash (and whether the
    /// message says "Segmentation fault" or "double free or corruption") depends on the heap of the C library; the port stops with the refusal at the first such pair, and every one of the oracle's runs that crashes has such a pair.
    /// </summary>
    [Test]
    public void HeapOverflowRuns_AreRefusedWithTheProof_AndEveryOtherRunHasNoSuchPair()
    {
        var crash = Oracle.Value.Where(c => c.Group == "last-crash").ToList();
        foreach (var c in crash)
        {
            var e = Assert.Throws<NotSupportedException>(() => Run(c, WithLast), c.Id)!;
            Assert.That(e.Message, Does.Contain("putlocalhom_last").And.Contain("valgrind").And.Contain("Invalid write of size 4"), c.Id);
        }

        Assert.That(crash, Is.Not.Empty);
    }

    /// <summary>
    /// The proof of the refusal: valgrind on the installed <c>pairlocalalign</c> with the recorded output of LAST (<c>last_valgrind_oracle.tsv</c>, <c>scripts/alignment/generate_last_valgrind_oracle.py</c>): every run that crashes has invalid writes, no run that
    /// ends has one - and the port refuses exactly the runs with invalid writes.
    /// </summary>
    [Test]
    public void ThePortRefusesExactlyTheRunsInWhichValgrindSeesInvalidWrites()
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.last_valgrind_oracle.tsv")!;
        using var reader = new StreamReader(stream);
        var vg = reader.ReadToEnd().Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.Split('\t')).ToDictionary(f => f[0], f => (Invalid: int.Parse(f[1], CultureInfo.InvariantCulture), Crashed: f[2] == "1"));
        Assert.That(vg, Has.Count.EqualTo(204));
        var disagreements = new List<string>();
        foreach (var c in Oracle.Value)
        {
            var (invalid, crashed) = vg[c.Id];
            if (crashed != (invalid > 0))
                disagreements.Add(c.Id + " valgrind " + invalid + " crashed " + crashed);
            bool refused;
            try
            {
                Run(c, WithLast);
                refused = false;
            }
            catch (NotSupportedException)
            {
                refused = true;
            }

            if (refused != invalid > 0)
                disagreements.Add(c.Id + " refused " + refused + " invalid writes " + invalid);
        }

        Assert.That(disagreements, Is.Empty, string.Join("\n", disagreements.Take(5)));
        Assert.That(vg.Values.Count(v => v.Invalid > 0), Is.EqualTo(100));
    }

    [Test]
    public void WithoutLast_TheScriptStopsWithItsInstallText()
    {
        var r = MafftProgram.Run(["--lastpair", InputName], new Dictionary<string, string> { [InputName] = ">a\nACGTACGTACGTACGTAAAC\n>b\nACGTACGTACGTACGTAAAG\n" });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stderr, Does.Contain("== Install LAST ="));
        });
    }
}
