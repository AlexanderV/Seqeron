using System.Globalization;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c stage 4, B09 F312): <see cref="MafftOptions.PartTree"/> = <c>mafft --parttree</c> / <c>--dpparttree</c> of MAFFT 7.505 (<c>splittbfast</c>: the sequences are classified by their distance to
/// <c>--partsize</c> representatives, followed recursively along the UPGMA tree of the representatives and aligned group by group; <c>--retree</c> passes, <c>--groupsize</c>, <c>--reorder</c>, <c>--treeout</c>) against the installed
/// MAFFT 7.505 (Debian 7.505-1; scripts/alignment/generate_mafft_parttree_oracle.py, mafft_parttree_oracle.tsv.gz). Every committed run must give the program's rows, the order and the names of the output (<c>Group-n</c> labels
/// of <c>--groupsize</c>), the guide tree file and the program's standard error (SHA-256 of the whole text, the full text for every 25th run).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftPartTree_Tests
{
    internal sealed record Row(string Id, string Strategy, string Options, string[] Seqs, int Rc, string[] Rows, string[] Names, string Sha, string Stderr, string Tree);

    internal static List<Row> LoadOracle(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + file;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], f[3].Split('|'), int.Parse(f[4], CultureInfo.InvariantCulture), f[5].Length == 0 ? [] : f[5].Split('|'), f[6].Length == 0 ? [] : f[6].Split('|'), f[7], Unescape(f[8]), Unescape(f[9])));
        }
        return list;
    }

    internal static string Unescape(string s)
    {
        var sb = new StringBuilder(s.Length);
        int i = 0;
        while (i < s.Length)
        {
            if (s[i] == '\\' && i + 1 < s.Length)
            {
                sb.Append(s[i + 1] switch { 'n' => '\n', 'r' => '\r', 'b' => '\b', 't' => '\t', _ => s[i + 1] });
                i += 2;
            }
            else
            {
                sb.Append(s[i]);
                i++;
            }
        }
        return sb.ToString();
    }

    internal static readonly Lazy<List<Row>> Oracle = new(() => LoadOracle("mafft_parttree_oracle.tsv.gz"));

    public static IEnumerable<TestCaseData> Strategies() => Oracle.Value.Select(r => r.Strategy).Distinct().Select(s => new TestCaseData(s).SetArgDisplayNames(s));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(2500));
            foreach (var g in rows.GroupBy(r => r.Strategy))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(rows.Count(r => r.Options.Contains("--dpparttree")), Is.GreaterThanOrEqualTo(1200));
            Assert.That(rows.Count(r => r.Options.Contains("--groupsize")), Is.GreaterThanOrEqualTo(400));
            Assert.That(rows.Count(r => r.Seqs.Length >= 100), Is.GreaterThanOrEqualTo(50), "sets above the size at which the classification collapses sequences (partsize 50)");
            Assert.That(rows.Count(r => r.Tree.Length > 0), Is.GreaterThanOrEqualTo(400), "guide tree files");
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(60), "full texts of the standard error");
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsNamesTreeAndStandardError_EqualMafft(string strategy) => CompareStrategy(Oracle.Value, strategy);

    internal static void CompareStrategy(List<Row> oracle, string strategy)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var r in oracle.Where(x => x.Strategy == strategy))
        {
            n++;
            var o = SequenceAligner_MafftProgressive_Tests.ParseOptions(r.Options);
            if (r.Rc != 0)
            {
                if (!SequenceAligner_MafftProgressive_Tests.Throws<ArgumentException>(() => SequenceAligner.RunMafft(r.Seqs, o)))
                    failures.Add(r.Id + ": the program refuses (rc " + r.Rc + "), the port does not");
                continue;
            }
            MafftRunResult res;
            try
            {
                res = SequenceAligner.RunMafft(r.Seqs, o);
            }
            catch (Exception e)
            {
                failures.Add(r.Id + ": " + e.GetType().Name + " " + e.Message);
                continue;
            }
            if (!res.Alignment.AlignedSequences.SequenceEqual(r.Rows))
            {
                failures.Add(r.Id + " (" + r.Options + "): rows differ");
                continue;
            }
            var names = res.InputIndex.Select((idx, k) => (res.GroupLabels is { } l ? l[k] + " " : "") + "s" + idx.ToString(CultureInfo.InvariantCulture)).ToArray();
            if (!names.SequenceEqual(r.Names))
                failures.Add(r.Id + " (" + r.Options + "): names / order differ");
            if ((res.TreeFile ?? "") != r.Tree)
                failures.Add(r.Id + " (" + r.Options + "): guide tree differs");
            string sha = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(res.StandardError))).ToLowerInvariant();
            if (sha != r.Sha)
                failures.Add(r.Id + " (" + r.Options + "): standard error differs");
            else if (r.Stderr.Length > 0 && r.Stderr != res.StandardError)
                failures.Add(r.Id + ": standard error text differs although the hash is equal");
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void Contracts()
    {
        string[] seqs = ["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQA", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLQA"];
        var o = new MafftOptions { AnySymbol = false, PartTree = MafftPartTree.KmerDistance };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.RunMafft(seqs, o).StandardError, Does.Contain("Strategy:\n FFT-NS-PartTree-2 (Not tested.)\n ?\n"));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { PartTree = MafftPartTree.DynamicProgramming }).StandardError, Does.Contain("FFT-NS-DPPartTree-1"));
            Assert.That(SequenceAligner.RunMafft(seqs, o with { Retree = 3, Fft = false }).StandardError, Does.Contain("NW-NS-PartTree-3"));
            Assert.Throws<ArgumentException>(() => SequenceAligner.RunMafft(seqs, o with { MaxIterate = 1 }), "dvtditr cannot read the missing distance matrix");
            Assert.Throws<ArgumentException>(() => SequenceAligner.RunMafft(seqs, o with { DistOut = true }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.RunMafft([.. seqs, "ACDE"], o), "sequences of fewer than 6 letters");
            Assert.That(SequenceAligner.RunMafft(seqs, o with { TreeOut = true }).TreeFile, Is.Not.Null);
            Assert.That(SequenceAligner.RunMafft(seqs, o).TreeFile, Is.Null);
        });
    }
}
