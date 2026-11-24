using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c, B09 F324): <see cref="MafftProgram"/> <c>--fastaparttree</c> against 200 sets of the user's Docker run <c>handoff/A2-19c-fasta34-battery</c> (MAFFT 7.505-1 + fasta34 34.26.5; mafft_fastaparttree_full_oracle.tsv.gz,
/// scripts/alignment/convert_fasta34_battery_returned.py): <c>mafft &lt;options&gt; in.fa</c> WITHOUT <c>--quiet</c>, 25 - 260 protein / DNA sequences (the 60 sets of the earlier run + 140), the options <c>--partsize</c>, <c>--groupsize</c>, <c>--retree 1</c>,
/// <c>--reorder</c>, <c>--treeout</c>: the return code, the aligned rows, the sequence names (their order), the tree file and the WHOLE standard error (progress counters, the matrix log, the strategy text) are equal; the one line of the standard error
/// that carries the process number (<c>rm -f /tmp/data-&lt;pid&gt;* ...</c>) is compared with the number removed (the port writes <c>&lt;pid&gt;</c>).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftFastaPartTreeFull_Tests
{
    internal sealed record Set(string Id, string Options, string[] Seqs, int Rc, string[] Rows, string[] Names, string Tree, string Stderr);

    internal static readonly Lazy<List<Set>> Oracle = new(() =>
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream("Seqeron.Genomics.Tests.TestData.Alignment.mafft_fastaparttree_full_oracle.tsv.gz") ?? throw new InvalidOperationException("oracle");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var list = new List<Set>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            var extra = JsonDocument.Parse(f[5]).RootElement;
            list.Add(new Set(f[0], f[1], f[2].Split('|'), int.Parse(f[3], CultureInfo.InvariantCulture), f[4].Length == 0 ? [] : f[4].Split('|'),
                extra.GetProperty("names").EnumerateArray().Select(e => e.GetString()!).Where(n => n.Length > 0).ToArray(), extra.GetProperty("tree").GetString()!, extra.GetProperty("stderr").GetString()!));
        }

        return list;
    });

    private static (string[] Names, string[] Rows) ParseFasta(string text)
    {
        var names = new List<string>();
        var rows = new List<string>();
        foreach (string l in text.Split('\n'))
        {
            if (l.StartsWith('>'))
            {
                names.Add(l[1..]);
                rows.Add("");
            }
            else if (l.Trim().Length > 0)
            {
                rows[^1] += l.Trim();
            }
        }

        return ([.. names], [.. rows]);
    }

    internal static MafftProgramResult Run(Set s)
    {
        string fasta = string.Concat(s.Seqs.Select((q, i) => $">s{i}\n{q}\n"));
        return MafftProgram.Run([.. s.Options.Split(' ', StringSplitOptions.RemoveEmptyEntries), "in.fa"], new Dictionary<string, string> { ["in.fa"] = fasta });
    }

    private static string FirstDifference(string ours, string theirs)
    {
        var a = ours.Split('\n');
        var b = theirs.Split('\n');
        for (int i = 0; i < Math.Max(a.Length, b.Length); i++)
        {
            string x = i < a.Length ? a[i] : "<end>";
            string y = i < b.Length ? b[i] : "<end>";
            if (x != y)
                return $"line {i + 1}: port '{x}', program '{y}'";
        }

        return "none";
    }

    [Test]
    public void Oracle_Shape()
    {
        var sets = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(sets, Has.Count.EqualTo(200));
            Assert.That(sets.Select(s => s.Options).Distinct().Count(), Is.EqualTo(7));
            Assert.That(sets.Count(s => s.Rc == 1), Is.EqualTo(1), "the set that the program refuses (\"Seq 31, too short.\")");
            Assert.That(sets.Max(s => s.Seqs.Length), Is.GreaterThanOrEqualTo(200));
            Assert.That(sets.Count(s => s.Tree.Length > 0), Is.EqualTo(20), "--treeout sets");
            Assert.That(sets.Count(s => s.Seqs[0].All(c => "ACGTacgtNn".Contains(c))), Is.GreaterThan(50));
            Assert.That(sets.All(s => !s.Stderr.Contains("/tmp/data-", StringComparison.Ordinal)), "the process-number line is removed");
        });
    }

    [Test]
    public void EverySet_ReturnCodeRowsNamesTreeAndStandardError_EqualTheProgram()
    {
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        Parallel.ForEach(Oracle.Value, s =>
        {
            MafftProgramResult res;
            try
            {
                res = Run(s);
            }
            catch (Exception e)
            {
                failures.Add($"{s.Id} ({s.Options}): {e.GetType().Name} {e.Message}");
                return;
            }

            string stderr = Regex.Replace(res.Stderr, @"^rm -f /tmp/data-<pid>\*  /tmp/query-<pid> /tmp/fasta-<pid>\n", "", RegexOptions.Multiline);
            var (names, rows) = ParseFasta(res.Stdout);
            string tree = res.Files.TryGetValue("in.fa.tree", out var t) ? t : "";
            if (res.ExitCode != s.Rc)
                failures.Add($"{s.Id} ({s.Options}): exit status {res.ExitCode}, program {s.Rc}");
            else if (!rows.SequenceEqual(s.Rows))
                failures.Add($"{s.Id} ({s.Options}): rows differ");
            else if (!names.SequenceEqual(s.Names))
                failures.Add($"{s.Id} ({s.Options}): names differ");
            else if (tree != s.Tree)
                failures.Add($"{s.Id} ({s.Options}): tree file differs");
            else if (stderr != s.Stderr)
                failures.Add($"{s.Id} ({s.Options}): standard error differs: {FirstDifference(stderr, s.Stderr)}");
        });

        Assert.That(failures, Is.Empty, $"{failures.Count} of 200 sets differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
    }
}
