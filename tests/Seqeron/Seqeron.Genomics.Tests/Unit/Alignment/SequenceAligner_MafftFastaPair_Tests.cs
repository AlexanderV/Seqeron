using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c, B09 F346 - F349): <c>mafft --fastapair</c> / <c>--fastswpair</c> (<c>dndfast7</c> and the alignment half of FASTA 3.4) against the user's Docker run <c>handoff/A2-19c-fastapair</c> (MAFFT 7.505-1 + fasta34 34.26.5 sha256 7491909b...;
/// scripts/alignment/convert_fastapair_returned.py; fasta34_align_oracle.tsv.gz, mafft_fastapair_oracle.tsv.gz): 180 sets (2 - 40 protein / DNA sequences, duplicates, unknown letters, low complexity; the options <c>--retree 1</c>, <c>--maxiterate 2 / 3</c>, <c>--treeout --distout</c>,
/// <c>--reorder</c>, <c>--anysymbol</c>, <c>--op --ep</c>, <c>--nofft</c>).
/// (1) the RAW output of every call of <c>fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000 query library ktup</c> that <c>dndfast7</c> makes (886 calls: the list of the best scores, the alignment blocks with their <c>fa_</c> / <c>sw_</c> / <c>bs_</c> fields and the aligned rows with their context, the histogram,
/// the Kolmogorov-Smirnov statistic) equals the port byte for byte (the paths, the clock lines and the two numbers that the program prints from uninitialised memory are normalised);
/// (2) <c>dndfast7</c> run directly: <c>hat2</c>, <c>hat3</c> (the address of the next block is normalised), the standard error (every set of the battery, including the 62 DNA sets that end with a segmentation fault of the program and the sets that stop with "count (fasta34 -z 3) = n") and the exit status are equal;
/// (3) the complete <c>mafft</c> run: exit status, rows, names, tree and hat2 files and the whole standard error are equal.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftFastaPair_Tests
{
    internal sealed record Call(string Id, string Args, string Query, string Library, string Out);

    internal sealed record Set(string Id, string Options, string[] Seqs, int Rc, string[] Rows, string[] Names, string Tree, string Hat2, string Stderr, int DirectRc, string DirectStderr, string DirectHat2, string DirectHat3);

    private static string Unescape(string s) => Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "n" => "\n", "t" => "\t", "r" => "\r", var other => other }, RegexOptions.Singleline);

    private static IEnumerable<string[]> Lines(string resource)
    {
        using var stream = typeof(SequenceAligner_MafftFastaPair_Tests).Assembly.GetManifestResourceStream("Seqeron.Genomics.Tests.TestData.Alignment." + resource) ?? throw new InvalidOperationException(resource);
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length > 0)
                yield return line.Split('\t');
        }
    }

    internal static readonly Lazy<List<Call>> Calls = new(() => Lines("fasta34_align_oracle.tsv.gz").Select(f => new Call(f[0], f[1], Unescape(f[2]), Unescape(f[3]), Unescape(f[4]))).ToList());

    internal static readonly Lazy<List<Set>> Sets = new(() => Lines("mafft_fastapair_oracle.tsv.gz").Select(f =>
    {
        var extra = JsonDocument.Parse(f[5]).RootElement;
        var direct = extra.GetProperty("direct");
        return new Set(f[0], f[1], f[2].Split('|'), int.Parse(f[3], CultureInfo.InvariantCulture), f[4].Length == 0 ? [] : f[4].Split('|'),
            extra.GetProperty("names").EnumerateArray().Select(e => e.GetString()!).Where(n => n.Length > 0).ToArray(), extra.GetProperty("tree").GetString()!, extra.GetProperty("hat2").GetString()!, extra.GetProperty("stderr").GetString()!,
            direct.GetProperty("rc").GetInt32(), direct.GetProperty("stderr").GetString()!, direct.GetProperty("hat2").GetString()!, direct.GetProperty("hat3").GetString()!);
    }).ToList());

    private static (List<string> Titles, List<string> Seqs) ParseFasta(string text)
    {
        var titles = new List<string>();
        var seqs = new List<string>();
        foreach (string l in text.Split('\n'))
        {
            if (l.StartsWith('>'))
            {
                titles.Add(l[1..]);
                seqs.Add("");
            }
            else if (l.Length > 0)
            {
                seqs[^1] += l;
            }
        }

        return (titles, seqs);
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

    private static string NormaliseSegfault(string stderr) => Regex.Replace(stderr, @"(/usr/bin/mafft: line \d+:)\s+(?:\d+|<pid>) (Segmentation fault)", "$1 P $2");

    private static string NormaliseHat3(string hat3) => Regex.Replace(hat3, @" (\(nil\)|0x[0-9a-f]+)$", m => m.Groups[1].Value == "(nil)" ? " NIL" : " PTR", RegexOptions.Multiline);

    [Test]
    public void Oracle_Shape()
    {
        var sets = Sets.Value;
        var calls = Calls.Value;
        Assert.Multiple(() =>
        {
            Assert.That(sets, Has.Count.EqualTo(180));
            Assert.That(calls, Has.Count.EqualTo(886));
            Assert.That(sets.Select(s => s.Options).Distinct().Count(), Is.EqualTo(12));
            Assert.That(sets.Count(s => s.Rc == 0), Is.EqualTo(77), "the sets that finish (one of them, 161, only because the program stores an unassigned number)");
            Assert.That(sets.Count(s => s.DirectRc == -11), Is.EqualTo(62), "DNA sets: dndfast7 stops with a segmentation fault (a null pointer from cutal, then putlocalhom2)");
            Assert.That(sets.Count(s => s.DirectStderr.Contains("count (fasta34 -z 3) =", StringComparison.Ordinal)), Is.EqualTo(41));
            Assert.That(calls.Count(c => c.Args.Contains("-A", StringComparison.Ordinal)), Is.EqualTo(276), "calls with -A (--fastswpair)");
            Assert.That(calls.Count(c => c.Out.Contains("       opt      E()", StringComparison.Ordinal)), Is.EqualTo(162), "calls with a histogram (more than 20 library sequences)");
            Assert.That(calls.Count(c => c.Out.Contains("; mp_KS: -0.0000 (N=UNINIT) at UNINIT", StringComparison.Ordinal)), Is.EqualTo(724), "calls whose statistic comes from uninitialised memory in the program");
        });
    }

    [Test]
    public void EveryRawFasta34Output_EqualsThePort_ByteForByte()
    {
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        Parallel.ForEach(Calls.Value, c =>
        {
            var a = c.Args.Split(' ');
            var (qt, qs) = ParseFasta(c.Query);
            var (lt, ls) = ParseFasta(c.Library);
            var run = SequenceAligner.MafftFasta34.Run(qt[0], qs[0], lt, ls, a.Contains("-n"), a.Contains("-A"), int.Parse(a[^1], CultureInfo.InvariantCulture));
            if (run.ExitCode != 0 || run.Text != c.Out)
                failures.Add($"{c.Id} ({c.Args}): exit {run.ExitCode}, {FirstDifference(run.Text, c.Out)}");
        });

        Assert.That(failures, Is.Empty, $"{failures.Count} of 886 calls differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
    }

    [Test]
    public void Dndfast7_Hat2Hat3StandardErrorAndExitStatus_EqualTheProgram_OnEverySet()
    {
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        int refused = 0;
        Parallel.ForEach(Sets.Value, s =>
        {
            var names = Enumerable.Range(0, s.Seqs.Length).Select(i => "=s" + i).ToList();
            SequenceAligner.MafftRun.Dndfast7Result r;
            try
            {
                r = SequenceAligner.MafftRun.Dndfast7(names, s.Seqs, s.Options.Contains("--fastswpair", StringComparison.Ordinal));
            }
            catch (NotSupportedException e)
            {
                // the one set whose program reads an unassigned variable (see Dndfast7_TheProgramReadsAnUnassignedScore_IsRefused)
                if (s.Id == "161" && e.Message.Contains("unassigned variable", StringComparison.Ordinal))
                    Interlocked.Increment(ref refused);
                else
                    failures.Add($"{s.Id} ({s.Options}): {e.GetType().Name} {e.Message}");
                return;
            }

            int expected = s.DirectRc == -11 ? 139 : s.DirectRc;
            if (r.ExitCode != expected)
                failures.Add($"{s.Id} ({s.Options}): exit status {r.ExitCode}, program {s.DirectRc}");
            else if (r.Err != s.DirectStderr)
                failures.Add($"{s.Id} ({s.Options}): standard error differs: {FirstDifference(r.Err, s.DirectStderr)}");
            else if (r.ExitCode == 0 && r.Hat2 != s.DirectHat2)
                failures.Add($"{s.Id} ({s.Options}): hat2 differs: {FirstDifference(r.Hat2, s.DirectHat2)}");
            else if (r.ExitCode == 0 && NormaliseHat3(r.Hat3) != s.DirectHat3)
                failures.Add($"{s.Id} ({s.Options}): hat3 differs: {FirstDifference(NormaliseHat3(r.Hat3), s.DirectHat3)}");
        });

        Assert.Multiple(() =>
        {
            Assert.That(failures, Is.Empty, $"{failures.Count} of 180 sets differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
            Assert.That(refused, Is.EqualTo(1));
        });
    }

    [Test]
    public void Dndfast7_TheProgramReadsAnUnassignedScore_IsRefused()
    {
        // set 161: the 6th sequence (40 times A) is a DNA query for fasta34 (more than 85 percent ACGTUN in more than 20 residues) in a protein run; the list of the best scores has the lines of DNA ("(  40) [f]  200 ..."),
        // sscanf( pt, "%d ..." ) fails on every one of them and the variable opt of ReadFasta34m10 was never assigned: the program stores 32767 (the upper half of a stack address under the emulator of the oracle run) and prints
        // distances of -25 ... -2728 to hat2. The value depends on the stack of the process, so the port refuses it instead of inventing a number.
        var s = Sets.Value.Single(x => x.Id == "161");
        var names = Enumerable.Range(0, s.Seqs.Length).Select(i => "=s" + i).ToList();
        var e = Assert.Throws<NotSupportedException>(() => SequenceAligner.MafftRun.Dndfast7(names, s.Seqs, true));
        Assert.That(e!.Message, Does.Contain("unassigned variable"));
        Assert.That(s.DirectHat2, Does.Contain("-2728.583"), "the program's distances of that run");
        var guess = s.Seqs[5];
        Assert.That(guess.All(c => c == 'A') && guess.Length == 40, Is.True);
    }

    [Test]
    public void EveryMafftRun_ReturnCodeRowsNamesTreeHat2AndStandardError_EqualTheProgram()
    {
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        int refused = 0;
        Parallel.ForEach(Sets.Value, s =>
        {
            MafftProgramResult res;
            try
            {
                string fasta = string.Concat(s.Seqs.Select((q, i) => $">s{i}\n{q}\n"));
                res = MafftProgram.Run([.. s.Options.Split(' ', StringSplitOptions.RemoveEmptyEntries), "in.fa"], new Dictionary<string, string> { ["in.fa"] = fasta });
            }
            catch (NotSupportedException e)
            {
                if (s.Id == "161" && e.Message.Contains("unassigned variable", StringComparison.Ordinal))
                    Interlocked.Increment(ref refused);
                else
                    failures.Add($"{s.Id} ({s.Options}): {e.GetType().Name} {e.Message}");
                return;
            }

            var (titles, seqs) = ParseFasta(res.Stdout);
            string tree = res.Files.TryGetValue("in.fa.tree", out var t) ? t : "";
            string hat2 = res.Files.TryGetValue("in.fa.hat2", out var h) ? h : "";
            string stderr = NormaliseSegfault(res.Stderr);
            if (res.ExitCode != s.Rc)
                failures.Add($"{s.Id} ({s.Options}): exit status {res.ExitCode}, program {s.Rc}");
            else if (!seqs.SequenceEqual(s.Rows))
                failures.Add($"{s.Id} ({s.Options}): rows differ");
            else if (!titles.SequenceEqual(s.Names))
                failures.Add($"{s.Id} ({s.Options}): names differ");
            else if (tree != s.Tree)
                failures.Add($"{s.Id} ({s.Options}): tree file differs");
            else if (hat2 != s.Hat2)
                failures.Add($"{s.Id} ({s.Options}): hat2 file differs: {FirstDifference(hat2, s.Hat2)}");
            else if (stderr != NormaliseSegfault(s.Stderr))
                failures.Add($"{s.Id} ({s.Options}): standard error differs: {FirstDifference(stderr, NormaliseSegfault(s.Stderr))}");
        });

        Assert.Multiple(() =>
        {
            Assert.That(failures, Is.Empty, $"{failures.Count} of 180 sets differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
            Assert.That(refused, Is.EqualTo(1));
        });
    }
}
