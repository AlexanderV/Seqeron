using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-20, B09 F118): <see cref="SequenceAligner.MultipleAlignKalign"/> against the installed Kalign 3.3.5 (Debian 1:3.3.5-1, AVX2 build).
/// Every committed run (scripts/alignment/generate_kalign_oracle.py, kalign_oracle.tsv) must give the program's rows - or, where the program refuses (alphabet / type mismatch,
/// fewer than two sequences with residues), an <see cref="ArgumentException"/>.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_KalignOracle_Tests
{
    internal sealed record Row(string Id, string Kind, string Options, int Threads, string[] Names, string[] Seqs, int Rc, string[] Rows);

    internal static readonly Lazy<List<Row>> Oracle = new(Load);

    private static List<Row> Load()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.kalign_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#' || line.StartsWith("id\t", StringComparison.Ordinal))
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], int.Parse(f[3]), f[4].Split('|'), f[5].Split('|'), int.Parse(f[6]), f[7].Length == 0 ? [] : f[7].Split('|')));
        }
        return list;
    }

    internal static KalignOptions ParseOptions(string options, int threads)
    {
        var t = options == "-" ? [] : options.Split(' ');
        var o = new KalignOptions { Threads = threads };
        for (int i = 0; i < t.Length; i += 2)
        {
            switch (t[i])
            {
                case "--type":
                    // kalign's set_aln_type: strstr(in, "rna"), "dna", "internal", "protein", "divergent" in this order ("internal" contains "rna")
                    var type = t[i + 1];
                    o = o with { Type = type.Contains("rna") ? KalignType.Rna : type.Contains("dna") ? KalignType.Dna : type.Contains("internal") ? KalignType.DnaInternal : type.Contains("protein") ? KalignType.Protein : KalignType.ProteinDivergent };
                    break;
                case "--gpo": o = o with { GapOpen = float.Parse(t[i + 1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "--gpe": o = o with { GapExtend = float.Parse(t[i + 1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                case "--tgpe": o = o with { TerminalGapExtend = float.Parse(t[i + 1], System.Globalization.CultureInfo.InvariantCulture) }; break;
                default: throw new InvalidOperationException(t[i]);
            }
        }
        return o;
    }

    public static IEnumerable<TestCaseData> Sets() => Oracle.Value.Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    [Test]
    public void Oracle_Has430Runs_Of5RefusalsAndEveryOption()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(300));
            Assert.That(rows.Count(r => r.Rc != 0), Is.EqualTo(5));
            Assert.That(rows.Select(r => r.Kind).Distinct(), Is.EquivalentTo(new[] { "dna", "rna", "prot" }));
            Assert.That(rows.Select(r => r.Threads).Distinct(), Is.EquivalentTo(new[] { 1, 2, 4 }));
        });
    }

    [TestCaseSource(nameof(Sets))]
    public void Rows_EqualKalign(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var o = ParseOptions(r.Options, r.Threads);
        if (r.Rc != 0)
        {
            Assert.Throws<ArgumentException>(() => SequenceAligner.KalignRun(r.Seqs, o));
            return;
        }
        var res = SequenceAligner.KalignRun(r.Seqs, o);
        Assert.That(res.Alignment.AlignedSequences, Is.EqualTo(r.Rows), r.Id + " " + r.Options + " -n " + r.Threads);
    }
}
