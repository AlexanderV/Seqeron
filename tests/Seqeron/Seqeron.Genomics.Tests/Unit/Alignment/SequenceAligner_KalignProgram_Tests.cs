using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-20, B09 F118): <see cref="KalignProgram"/> against whole command lines of the installed Kalign 3.3.5
/// (scripts/alignment/generate_kalign_command_oracle.py, kalign_command_oracle.tsv: 83 runs - every type word, gap options incl. C <c>atof</c> prefixes, thread counts, the <c>getopt_long_only</c>
/// forms, abbreviations, two input files, a missing file, type / alphabet mismatch, one sequence only (the program segfaults), 60-column output).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_KalignProgram_Tests
{
    private sealed record Row(string Id, string[] Args, string In, string In2, int Rc, string Out);

    private static string Unesc(string s) => System.Text.RegularExpressions.Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", var c => c });

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.kalign_command_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1].Split('|'), Unesc(f[2]), Unesc(f[3]), int.Parse(f[4]), Unesc(f[5])));
        }
        return list;
    });

    public static IEnumerable<TestCaseData> Runs() => Oracle.Value.Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    [Test]
    public void Oracle_Has83Runs_With27Refusals()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value, Has.Count.EqualTo(83));
            Assert.That(Oracle.Value.Count(r => r.Rc != 0), Is.EqualTo(27));
        });
    }

    [TestCaseSource(nameof(Runs))]
    public void Run_EqualsTheProgram(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var files = new Dictionary<string, string> { ["in"] = r.In };
        if (r.In2.Length > 0)
            files["in2"] = r.In2;
        var res = SequenceAligner_KalignProgramHelper.Run(r.Args, files);
        if (r.Rc != 0)
        {
            Assert.That(res.ExitCode, Is.Not.Zero, id + " " + string.Join(' ', r.Args));
            return;
        }
        string actual = res.Files.TryGetValue("out", out var text) ? text : res.Stdout;
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.Zero, id + " " + string.Join(' ', r.Args) + " " + string.Join("; ", res.Messages));
            Assert.That(actual, Is.EqualTo(r.Out), id + " " + string.Join(' ', r.Args));
        });
    }

    [Test]
    public void ParseType_InternalMeansRna_BecauseStrstrFindsRnaFirst()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KalignProgram.ParseType("internal"), Is.EqualTo(KalignType.Rna));
            Assert.That(KalignProgram.ParseType("dna"), Is.EqualTo(KalignType.Dna));
            Assert.That(KalignProgram.ParseType("xrnax"), Is.EqualTo(KalignType.Rna));
            Assert.That(KalignProgram.ParseType("protein"), Is.EqualTo(KalignType.Protein));
            Assert.That(KalignProgram.ParseType("divergent"), Is.EqualTo(KalignType.ProteinDivergent));
            Assert.Throws<ArgumentException>(() => KalignProgram.ParseType("RNA"));
        });
    }

    [Test]
    public void Parse_GetoptLongOnly_Forms()
    {
        var c = KalignProgram.Parse(["-i", "a.fa", "b.fa", "-type", "dna", "--gpo=3.5", "-gpe", "2", "--tgpe", "0.5x", "-n4", "-q", "-f", "fasta", "-o", "o.fa"]);
        Assert.Multiple(() =>
        {
            Assert.That(c.Options, Is.EqualTo(new KalignOptions { Type = KalignType.Dna, GapOpen = 3.5f, GapExtend = 2f, TerminalGapExtend = 0.5f, Threads = 4 }));
            Assert.That(c.InputFiles, Is.EqualTo(new[] { "a.fa", "b.fa" }));
            Assert.That(c.OutFile, Is.EqualTo("o.fa"));
            Assert.That(c.Quiet, Is.True);
            Assert.Throws<ArgumentException>(() => KalignProgram.Parse(["--gp", "3"]), "--gp is a prefix of --gpo and --gpe");
            Assert.Throws<ArgumentException>(() => KalignProgram.Parse(["-n", "0"]));
            Assert.Throws<ArgumentException>(() => KalignProgram.Parse(["--format", "xyz"]));
        });
    }

    [Test]
    public void MsfAndClustalOutput_AreAcceptedByTheParser_ButNotWritten()
    {
        var files = new Dictionary<string, string> { ["in"] = ">a\nACGTACGT\n>b\nACGTCGT\n" };
        Assert.Multiple(() =>
        {
            Assert.That(KalignProgram.Parse(["-i", "in", "--format", "msf"]).Format, Is.EqualTo("msf"));
            Assert.Throws<NotSupportedException>(() => KalignProgram.Run(["-q", "-i", "in", "--format", "clu"], files));
        });
    }
}

internal static class SequenceAligner_KalignProgramHelper
{
    public static KalignProgramResult Run(string[] args, Dictionary<string, string> files) => KalignProgram.Run(args, files);
}
