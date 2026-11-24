using System.Globalization;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c stage 4, B09 F313): <c>mafft --retree 0</c> (the guide tree only: the sequences come back not aligned), <c>--treeout</c>, <c>--nodeout</c> and <c>--distout</c> of MAFFT 7.505 against the installed program
/// (scripts/alignment/generate_mafft_treeout_oracle.py, mafft_treeout_oracle.tsv.gz): rows, names, standard error, the guide tree file <c>input.tree</c> and the distance file <c>input.hat2</c>, byte for byte.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftTreeOut_Tests
{
    internal sealed record Row(string Id, string Strategy, string Options, string[] Seqs, int Rc, string[] Rows, string[] Names, string Sha, string Stderr, string Tree, string Hat2);

    internal static List<Row> Load()
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream("Seqeron.Genomics.Tests.TestData.Alignment.mafft_treeout_oracle.tsv.gz") ?? throw new InvalidOperationException("oracle");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            var u = SequenceAligner_MafftPartTree_Tests.Unescape;
            list.Add(new Row(f[0], f[1], f[2], f[3].Split('|'), int.Parse(f[4], CultureInfo.InvariantCulture), f[5].Length == 0 ? [] : f[5].Split('|'), f[6].Length == 0 ? [] : f[6].Split('|'), f[7], u(f[8]), u(f[9]), u(f[10])));
        }
        return list;
    }

    internal static readonly Lazy<List<Row>> Oracle = new(Load);

    public static IEnumerable<TestCaseData> Strategies() => Oracle.Value.Select(r => r.Strategy).Distinct().Select(s => new TestCaseData(s).SetArgDisplayNames(s));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(1000));
            foreach (var g in rows.GroupBy(r => r.Strategy))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(rows.Count(r => r.Tree.Length > 0), Is.GreaterThanOrEqualTo(500));
            Assert.That(rows.Count(r => r.Hat2.Length > 0), Is.GreaterThanOrEqualTo(250));
            Assert.That(rows.Count(r => r.Tree.Contains("Density:")), Is.GreaterThanOrEqualTo(100), "--nodeout");
            Assert.That(rows.Count(r => r.Options.Contains("--retree 0")), Is.GreaterThanOrEqualTo(400));
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(30), "full texts of the standard error");
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsTreeDistancesAndStandardError_EqualMafft(string strategy)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var r in Oracle.Value.Where(x => x.Strategy == strategy))
        {
            n++;
            var o = SequenceAligner_MafftProgressive_Tests.ParseOptions(r.Options) with { Names = r.Seqs.Select((_, i) => "s" + i.ToString(CultureInfo.InvariantCulture)).ToArray() };
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
                failures.Add(r.Id + " (" + r.Options + "): rows differ");
            else if (!res.InputIndex.Select(i => "s" + i.ToString(CultureInfo.InvariantCulture)).SequenceEqual(r.Names))
                failures.Add(r.Id + ": order differs");
            else if ((res.TreeFile ?? "") != r.Tree)
                failures.Add(r.Id + " (" + r.Options + "): tree file differs");
            else if ((res.DistanceFile ?? "") != r.Hat2)
                failures.Add(r.Id + " (" + r.Options + "): distance file differs");
            else
            {
                string sha = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(res.StandardError))).ToLowerInvariant();
                if (sha != r.Sha)
                    failures.Add(r.Id + " (" + r.Options + "): standard error differs");
                else if (r.Stderr.Length > 0 && r.Stderr != res.StandardError)
                    failures.Add(r.Id + ": standard error text differs although the hash is equal");
            }
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void Contracts()
    {
        string[] seqs = ["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQA", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLQA"];
        var o = new MafftOptions { AnySymbol = false, Names = ["a", "b c", "d", "e"] };
        var res = SequenceAligner.RunMafft(seqs, o with { Retree = 0, TreeOut = true, DistOut = true });
        Assert.Multiple(() =>
        {
            Assert.That(res.Alignment.AlignedSequences, Is.EqualTo(seqs), "--retree 0 does not align");
            Assert.That(res.TreeFile, Does.Contain("1_a").And.Contain("2_b_c").And.EndWith(";\n"));
            Assert.That(res.DistanceFile, Does.StartWith("    1\n    4\n").And.Contain("   2. =b c\n"));
            Assert.That(res.StandardError, Does.Contain("FFT-NS-0"));
            Assert.That(SequenceAligner.RunMafft(seqs, o).TreeFile, Is.Null);
            Assert.That(SequenceAligner.RunMafft(seqs, o with { NodeOut = true }).TreeFile, Does.Contain("Density:").And.Contain("Node info:"));
            Assert.Catch<ArgumentException>(() => SequenceAligner.RunMafft(seqs, o with { NodeOut = true, MaxIterate = 2 }));
            Assert.Catch<ArgumentException>(() => SequenceAligner.RunMafft(seqs, o with { Names = ["a"] }));
        });
    }
}
