using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19, B09 F119): <see cref="SequenceAligner.MultipleAlignMafft"/> (stage 1: two sequences) against the installed MAFFT 7.505 (Debian 7.505-1;
/// scripts/alignment/generate_mafft_pair_oracle.py, mafft_pair_oracle.tsv). Every committed run must give the program's rows (original letters restored, gaps '-').
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftPair_Tests
{
    internal sealed record Row(string Id, string Options, string[] Seqs, int Rc, string[] Rows);

    internal static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mafft_pair_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#' || line.StartsWith("id\t", StringComparison.Ordinal))
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2].Split('|'), int.Parse(f[3]), f[4].Length == 0 ? [] : f[4].Split('|')));
        }
        return list;
    });

    internal static MafftOptions ParseOptions(string options)
    {
        if (options == "plain")
            return new MafftOptions { AnySymbol = false };
        var o = new MafftOptions();
        if (options == "-")
            return o;
        var t = options.Split(' ');
        int i = 0;
        while (i < t.Length)
        {
            switch (t[i++])
            {
                case "--op": o = o with { GapOpen = double.Parse(t[i++], CultureInfo.InvariantCulture) }; break;
                case "--ep": o = o with { Offset = double.Parse(t[i++], CultureInfo.InvariantCulture) }; break;
                case "--bl": o = o with { Blosum = int.Parse(t[i++], CultureInfo.InvariantCulture) }; break;
                case "--nuc": o = o with { SequenceType = MafftSequenceType.Nucleotide }; break;
                case "--amino": o = o with { SequenceType = MafftSequenceType.Amino }; break;
                case "--retree": o = o with { Retree = int.Parse(t[i++], CultureInfo.InvariantCulture) }; break;
                case "--nofft": o = o with { Fft = false }; break;
                case "--fft": o = o with { Fft = true }; break;
                case "--reorder": o = o with { Reorder = true }; break;
                default: throw new InvalidOperationException(t[i - 1]);
            }
        }
        return o;
    }

    public static IEnumerable<TestCaseData> Sets() => Oracle.Value.Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    [Test]
    public void Oracle_HasEveryKindOfRun()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(600));
            Assert.That(rows.Select(r => r.Options).Distinct().Count(), Is.GreaterThan(30));
            Assert.That(rows.Count(r => r.Options == "plain"), Is.GreaterThanOrEqualTo(50));
        });
    }

    [TestCaseSource(nameof(Sets))]
    public void Rows_EqualMafft(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var o = ParseOptions(r.Options);
        if (r.Rc != 0)
        {
            Assert.Catch<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(r.Seqs, o));
            return;
        }
        var res = SequenceAligner.MultipleAlignMafft(r.Seqs, o);
        Assert.That(res.AlignedSequences, Is.EqualTo(r.Rows), r.Id + " " + r.Options + " " + string.Join(" / ", r.Seqs));
    }
}
