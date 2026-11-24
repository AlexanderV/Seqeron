using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c stage 4, B09 WP70, F315 ff.): the script <c>mafft</c> of MAFFT 7.505 as a command line (<see cref="MafftProgram"/>): the options as the script parses them, the messages and the exit status of every check,
/// the output formats, the files <c>input.tree</c> / <c>input.hat2</c>, against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_command_oracle.py</c>, <c>mafft_command_oracle.tsv.gz</c>): return code,
/// standard output and the files byte for byte, the standard error by SHA-256 (the full text for every 25th case and for the short ones).
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftProgram_Tests
{
    internal sealed record Case(string Id, string Group, string[] Args, string Input, Dictionary<string, string> Files, int Rc, string Out, string ErrSha, string Err, Dictionary<string, string> OutFiles);

    private static Dictionary<string, string> JsonMap(string json)
    {
        using var doc = System.Text.Json.JsonDocument.Parse(json);
        return doc.RootElement.EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
    }

    internal static List<Case> Load() => Load("mafft_command_oracle.tsv");

    internal static List<Case> Load(string file)
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment." + file) ?? throw new InvalidOperationException("oracle resource missing: " + file);
        using var reader = new StreamReader(stream, Encoding.UTF8);
        var list = new List<Case>();
        string? line;
        var u = SequenceAligner_MafftPartTree_Tests.Unescape;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            string[] args = System.Text.Json.JsonSerializer.Deserialize<string[]>(u(f[2]))!;
            list.Add(new Case(f[0], f[1], args, u(f[3]), JsonMap(u(f[4])), int.Parse(f[5], CultureInfo.InvariantCulture), u(f[6]), f[7], u(f[8]), JsonMap(u(f[9]))));
        }
        return list;
    }

    internal static readonly Lazy<List<Case>> Oracle = new(Load);

    public static IEnumerable<TestCaseData> Groups() => Oracle.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    internal const string InputName = "in.fa";

    internal static MafftProgramResult RunCase(Case c)
    {
        var files = new Dictionary<string, string>(c.Files) { [InputName] = c.Input };
        return MafftProgram.Run([.. c.Args, InputName], files);
    }

    [Test]
    public void Oracle_HasEveryGroup()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(2900));
            foreach (var g in rows.Where(r => r.Group != "contract").GroupBy(r => r.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(rows.Count(r => r.Group == "contract"), Is.GreaterThanOrEqualTo(80));
            Assert.That(rows.Count(r => r.OutFiles.Keys.Any(k => k.EndsWith(".tree", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(500), "guide tree files");
            Assert.That(rows.Count(r => r.OutFiles.Keys.Any(k => k.EndsWith(".hat2", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(300), "distance files");
            Assert.That(rows.Count(r => r.Args.Contains("--clustalout")), Is.GreaterThanOrEqualTo(200));
            Assert.That(rows.Count(r => r.Args.Contains("--phylipout")), Is.GreaterThanOrEqualTo(80));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(60), "the refusals of the script");
            Assert.That(rows.Count(r => r.Err.Length > 0), Is.GreaterThanOrEqualTo(150), "full texts of the standard error");
        });
    }

    /// <summary>Runs one case of an oracle and compares return code, standard output, the files and the standard error (by SHA-256; the text where the oracle holds it) with the installed program: null when equal.</summary>
    internal static string? Verify(Case c)
    {
        MafftProgramResult r;
        try
        {
            r = RunCase(c);
        }
        catch (Exception e)
        {
            return c.Id + " (" + string.Join(' ', c.Args) + "): " + e.GetType().Name + " " + e.Message;
        }
        string what = c.Id + " (" + string.Join(' ', c.Args) + ")";
        if (r.ExitCode != c.Rc)
            return what + ": return code " + r.ExitCode + " != " + c.Rc;
        if (r.Stdout != c.Out)
            return what + ": standard output differs";
        if (!r.Files.OrderBy(k => k.Key, StringComparer.Ordinal).SequenceEqual(c.OutFiles.OrderBy(k => k.Key, StringComparer.Ordinal)))
            return what + ": files differ (" + string.Join(',', r.Files.Keys) + " vs " + string.Join(',', c.OutFiles.Keys) + ")";
        if (c.ErrSha != "-")
        {
            string sha = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(r.Stderr))).ToLowerInvariant();
            if (sha != c.ErrSha)
                return what + ": standard error differs";
            if (c.Err.Length > 0 && c.Err != r.Stderr)
                return what + ": standard error text differs although the hash is equal";
        }
        return null;
    }

    [TestCaseSource(nameof(Groups))]
    public void ReturnCodeOutputFilesAndStandardError_EqualMafft(string group)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var c in Oracle.Value.Where(x => x.Group == group))
        {
            n++;
            if (Verify(c) is { } failure)
                failures.Add(failure);
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void SplitCommandLine_FollowsTheShell()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MafftProgram.SplitCommandLine("--op 1.5 --anysymbol  in.fa"), Is.EqualTo(new[] { "--op", "1.5", "--anysymbol", "in.fa" }));
            Assert.That(MafftProgram.SplitCommandLine("--out 'my file.txt' \"a b\" c\\ d"), Is.EqualTo(new[] { "--out", "my file.txt", "a b", "c d" }));
            Assert.That(MafftProgram.SplitCommandLine(""), Is.Empty);
            Assert.That(MafftProgram.SplitCommandLine("''"), Is.EqualTo(new[] { "" }));
            Assert.Catch<ArgumentException>(() => MafftProgram.SplitCommandLine("--op '1"));
            Assert.Catch<ArgumentException>(() => MafftProgram.SplitCommandLine("--op \"1"));
        });
    }

    [Test]
    public void Contracts()
    {
        string fasta = ">a\nMKVLAAGIVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n";
        var files = new Dictionary<string, string> { [InputName] = fasta };
        var plain = MafftProgram.Run(["--quiet", InputName], files);
        Assert.Multiple(() =>
        {
            Assert.That(plain.ExitCode, Is.EqualTo(0));
            Assert.That(plain.Stdout, Is.EqualTo(">a\nMKVLAAGIVGLLLA\n>b\nMKVL-AGIVGLLA-\n>c\nMKVIAAGLVGLLLA\n"));
            Assert.That(plain.Stderr, Is.Empty);
            Assert.That(MafftProgram.Run(["--version"]).Stderr, Is.EqualTo("v7.505 (2022/Apr/10)\n"));
            Assert.That(MafftProgram.Run(["--foo", InputName], files).ExitCode, Is.EqualTo(1));
            Assert.That(MafftProgram.Run(["--foo", InputName], files).Stderr, Does.StartWith("Unknown option:  --foo\n---"));
            Assert.That(MafftProgram.Run(["missing.fa"], files).Stderr, Does.StartWith("/usr/bin/mafft: Cannot open missing.fa.\n\n"));
            Assert.That(MafftProgram.Run(["--quiet", "--out", "result.txt", InputName], files).Files["result.txt"], Is.EqualTo(plain.Stdout));
            Assert.That(MafftProgram.Run(["--quiet", "--treeout", InputName], files).Files, Does.ContainKey("in.fa.tree"));
            Assert.That(MafftProgram.Run(["--quiet", "-"], files, stdin: fasta).Stdout, Is.EqualTo(plain.Stdout), "the input file - is standard input");
            Assert.That(MafftProgram.Run(["--quiet", "--progress", "/tmp/p.txt", InputName], files).Files["/tmp/p.txt"], Is.Empty.Or.Not.Null, "the progress file is a file of the run");
            Assert.Throws<ArgumentNullException>(() => MafftProgram.Run(null!));
            Assert.Throws<ArgumentNullException>(() => MafftProgram.Run(["a", null!]));
            Assert.Throws<NotSupportedException>(() => MafftProgram.Run([]), "no arguments: the interactive wizard");
        });
    }
}
