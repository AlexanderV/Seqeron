using System.Globalization;
using System.IO.Compression;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c stage 4, B09 F312): <see cref="MafftPartTree.Fasta"/> = <c>mafft --fastaparttree</c> of MAFFT 7.505 (<c>splittbfast -S</c> with the optimised scores of <c>fasta34</c>, FASTA 3.4: ported). The oracle is the
/// user's returned Docker run 3 (docs/Validation/review-2026-09/handoff/A2-19c-fastaparttree/returned/: MAFFT 7.505-1 + fasta34 built from fasta-34.26.5; deterministic: the repeat run is byte-identical; scripts/alignment/import_fastaparttree_run.py,
/// mafft_fastaparttree_oracle.tsv.gz: 60 sets of 25 - 60 protein / DNA sequences, the options --partsize, --groupsize, --retree): the rows of every set must be equal.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftFastaPartTree_Tests
{
    internal sealed record Row(string Id, string Options, string[] Seqs, int Rc, string[] Rows);

    internal static List<Row> Load()
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream("Seqeron.Genomics.Tests.TestData.Alignment.mafft_fastaparttree_oracle.tsv.gz") ?? throw new InvalidOperationException("oracle");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2].Split('|'), int.Parse(f[3], CultureInfo.InvariantCulture), f[4].Length == 0 ? [] : f[4].Split('|')));
        }
        return list;
    }

    internal static readonly Lazy<List<Row>> Oracle = new(Load);

    [Test]
    public void Oracle_Shape()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(60));
            Assert.That(rows.Select(r => r.Options.Replace(" --quiet", "")).Distinct().Count(), Is.EqualTo(4));
            Assert.That(rows.Count(r => r.Seqs[0].All(c => "ACGT".Contains(c))), Is.GreaterThanOrEqualTo(25), "DNA sets (the score of the nucleotide alphabet)");
            Assert.That(rows.Count(r => !r.Seqs[0].All(c => "ACGT".Contains(c))), Is.GreaterThanOrEqualTo(25), "protein sets");
            Assert.That(rows.All(r => r.Rc == 0));
        });
    }

    [Test]
    public void Rows_EqualTheProgramRun()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var o = SequenceAligner_MafftProgressive_Tests.ParseOptions(r.Options);
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
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void Contracts()
    {
        string[] seqs = ["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQA", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLQA"];
        var o = new MafftOptions { AnySymbol = false, PartTree = MafftPartTree.Fasta };
        var res = SequenceAligner.RunMafft(seqs, o);
        Assert.Multiple(() =>
        {
            Assert.That(res.StandardError, Does.Contain("Strategy:\n FFT-NS-FastaPartTree-1 (Not tested.)"));
            Assert.That(res.StandardError, Does.Contain("fastapath=fasta34\n").And.Contain("done.\nA: "));
            Assert.That(res.Alignment.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(seqs));
        });
    }
}
