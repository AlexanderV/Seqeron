using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-20b, B09 F222): <see cref="KalignProgram"/> output formats <c>msf</c> / <c>clu</c> (GCG checksums, date line, blocks), MSF / Clustal inputs with the program's reader quirks, several
/// inputs, and the non-quiet standard output (banner, log lines, warnings, <c>-n 1</c> alignment lines) against the installed Kalign 3.3.5 run under the time shim
/// (scripts/alignment/generate_kalign_formats_oracle.py, kalign_formats_oracle.tsv).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_KalignFormatsOracle_Tests
{
    private sealed record Row(string Id, string[] Args, string In, string In2, string In3, long Epoch, int Rc, string Stdout, string Out);

    private static string Unesc(string s) => Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", var c => c });

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.kalign_formats_oracle.tsv";
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
            list.Add(new Row(f[0], f[1].Length == 0 ? [] : f[1].Split('|'), Unesc(f[2]), Unesc(f[3]), Unesc(f[4]), long.Parse(f[5]), int.Parse(f[6]), Unesc(f[7]), Unesc(f[8])));
        }
        return list;
    });

    private static string Normalize(string stdout) => Regex.Replace(stdout, @"CPU Time: .*", "CPU Time: ~");

    [Test]
    public void Oracle_HasTheRecordedRuns()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value, Has.Count.GreaterThanOrEqualTo(250));
            Assert.That(Oracle.Value.Count(r => r.Args.Contains("msf") || r.Args.Contains("clu")), Is.GreaterThan(30));
            Assert.That(Oracle.Value.Count(r => r.Rc != 0), Is.GreaterThan(20));
        });
    }

    [Test]
    public void Run_EveryCommandLine_EqualsTheProgram()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var files = new Dictionary<string, string> { ["in"] = r.In };
            if (r.In2.Length > 0)
                files["in2"] = r.In2;
            if (r.In3.Length > 0)
                files["in3"] = r.In3;
            var runtime = new KalignRuntime { Now = () => DateTime.UnixEpoch.AddSeconds(r.Epoch), CpuSeconds = () => 0 };
            var res = KalignProgram.Run(r.Args, files, "", runtime);
            string outFile = res.Files.Values.FirstOrDefault() ?? "";
            if ((r.Rc != 0) != (res.ExitCode != 0))
                failures.Add($"{r.Id} exit {res.ExitCode} (program {r.Rc}) {string.Join(' ', r.Args)} {string.Join("; ", res.Messages)}");
            else if (r.Rc is >= 0 and < 128 && Normalize(res.Stdout) != r.Stdout) // after a crash of the program (signal) its buffered text is not meaningful
                failures.Add($"{r.Id} stdout {string.Join(' ', r.Args)}\n--- expected\n{r.Stdout}\n--- actual\n{Normalize(res.Stdout)}");
            else if (r.Rc == 0 && outFile != r.Out)
                failures.Add($"{r.Id} file {string.Join(' ', r.Args)}\n--- expected\n{r.Out}\n--- actual\n{outFile}");
        }
        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(5)) + $"\n({failures.Count} of {Oracle.Value.Count})");
    }
}
