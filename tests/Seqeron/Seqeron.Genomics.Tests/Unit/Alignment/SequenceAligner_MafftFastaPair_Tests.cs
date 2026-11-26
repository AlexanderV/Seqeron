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

    internal sealed record Set(string Id, string Options, string[] Seqs, int Rc, string[] Rows, string[] Names, string Tree, string Hat2, string Stderr, int DirectRc, string DirectStderr, string DirectHat2, string DirectHat3, bool Full, Dictionary<string, string> Files, string? Stdout, Dictionary<string, string> OutFiles);

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

    // run "" = handoff/A2-19c-fastapair (180 sets), run "2" = handoff/A2-19c-fastapair2 (all 288 sets: run 2 supplies the raw calls and the direct runs of dndfast7 of sets 1 - 253, run 3 the complete runs of mafft of 1 - 253, run 4 everything of sets 254 - 288)
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<string, List<Call>> CallCache = new();
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<string, List<Set>> SetCache = new();

    internal static List<Call> CallsOf(string run) => CallCache.GetOrAdd(run, r => Lines($"fasta34_align{r}_oracle.tsv.gz").Select(f => new Call(f[0], f[1], Unescape(f[2]), Unescape(f[3]), Unescape(f[4]))).ToList());

    internal static List<Set> SetsOf(string run) => SetCache.GetOrAdd(run, r => Lines($"mafft_fastapair{r}_oracle.tsv.gz").Select(f =>
    {
        var extra = JsonDocument.Parse(f[5]).RootElement;
        var direct = extra.GetProperty("direct");
        bool full = !extra.TryGetProperty("full", out var fl) || fl.GetBoolean();
        var files = new Dictionary<string, string>();
        if (extra.TryGetProperty("files", out var fs))
        {
            foreach (var p in fs.EnumerateObject())
                files[p.Name] = p.Value.GetString()!;
        }

        var outFiles = new Dictionary<string, string>();
        if (extra.TryGetProperty("outfiles", out var of))
        {
            foreach (var p in of.EnumerateObject())
                outFiles[p.Name] = p.Value.GetString()!;
        }

        string? stdout = extra.TryGetProperty("stdout", out var so) ? so.GetString() : null;      // run 4 of the package: the raw standard output of the program (CLUSTAL, PHYLIP ...); runs 1 - 3 have the rows only
        return new Set(f[0], f[1], f[2].Split('|'), full ? int.Parse(f[3], CultureInfo.InvariantCulture) : 0, f[4].Length == 0 ? [] : f[4].Split('|'),
            extra.GetProperty("names").EnumerateArray().Select(e => e.GetString()!).Where(n => n.Length > 0).ToArray(), extra.GetProperty("tree").GetString()!, extra.GetProperty("hat2").GetString()!, extra.GetProperty("stderr").GetString()!,
            direct.GetProperty("rc").GetInt32(), direct.GetProperty("stderr").GetString()!, direct.GetProperty("hat2").GetString()!, direct.GetProperty("hat3").GetString()!, full, files, stdout, outFiles);
    }).ToList());

    internal static List<Call> Calls => CallsOf("");

    internal static List<Set> Sets => SetsOf("");

    /// <summary>The fasta text of a set: the file <c>in.fa</c> that the generator named, or the sequences as <c>s0, s1, ...</c>.</summary>
    private static string FastaOf(Set s) => s.Files.TryGetValue("in.fa", out var f) ? f : string.Concat(s.Seqs.Select((q, i) => $">s{i}\n{q}\n"));

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
            else if (l.Length > 0 && seqs.Count > 0)      // lines before the first header (CLUSTAL, PHYLIP ...) are not rows
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

    /// <summary>
    /// With <c>--thread n</c> (n &gt; 1) the steps of the progressive alignment (<c>treebasethread</c> of <c>tbfast</c>) run on n threads: the program prints "STEP k /m (thread t) " when a thread has finished a step and "Reallocating (by thread t) ..done. *alloclen = a" when the thread that took the step has to grow the buffers.
    /// The steps are taken in their order, so the sizes of the buffers (a) and the list of the reallocations are fixed by the guide tree, but WHICH thread takes a step and in which order the threads print is up to the scheduler: ten runs of the installed <c>mafft --retree 2 --thread 2</c> on one input gave ten different
    /// orders and thread numbers, and sets 269 and 285 of the oracle show two different assignments of the five steps of the same tree shape. Only that is removed from both sides: the thread numbers and the order of the "STEP" lines (their set of numbers k /m stays), the order of the "Reallocating" lines relative to the
    /// "STEP" lines (their list of sizes stays); every other character of the standard error, and the rows, are compared byte for byte. With one thread (or none) nothing is normalised.
    /// </summary>
    private static string NormaliseThreadSchedule(string stderr)
    {
        var steps = new List<(int K, string Total)>();
        var reallocs = new List<string>();
        string rest = Regex.Replace(stderr, @"\rSTEP\s+(\d+) /(\d+) \(thread\s+\d+\) ", m =>
        {
            steps.Add((int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture), m.Groups[2].Value));
            return "";
        });
        rest = Regex.Replace(rest, @"\nReallocating \(by thread \d+\) \.\.done\. \*alloclen = (\d+)\n", m =>
        {
            reallocs.Add(m.Groups[1].Value);
            return "";
        });
        return rest + "\n#STEPS " + string.Join(",", steps.OrderBy(x => x.K).Select(x => $"{x.K}/{x.Total}")) + "\n#REALLOC " + string.Join(",", reallocs);
    }

    private static bool RunsOnSeveralThreads(string options)
    {
        var parts = options.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        int i = Array.IndexOf(parts, "--thread");
        return i >= 0 && i + 1 < parts.Length && int.TryParse(parts[i + 1], CultureInfo.InvariantCulture, out int n) && n > 1;
    }

    private static string NormaliseHat3(string hat3) => Regex.Replace(hat3, @" (\(nil\)|0x[0-9a-f]+)$", m => m.Groups[1].Value == "(nil)" ? " NIL" : " PTR", RegexOptions.Multiline);

    [Test]
    public void Oracle_Shape()
    {
        var sets = Sets;
        var calls = Calls;
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
    public void Oracle2_Shape()
    {
        var sets = SetsOf("2");
        var calls = CallsOf("2");
        Assert.Multiple(() =>
        {
            Assert.That(sets, Has.Count.EqualTo(288), "all sets of handoff/A2-19c-fastapair2: 1 - 253 from runs 2 + 3, 254 - 288 from run 4");
            Assert.That(calls, Has.Count.EqualTo(695), "raw fasta34 calls: 551 of run 2 + 144 of run 4 (the sets 254 - 288 whose dndfast7 reached fasta34)");
            Assert.That(sets.Count(s => s.Full), Is.EqualTo(288), "every set has the complete mafft run (results row, standard error); no set of run 4 recorded an error (results.tsv column 10 empty in all 35 rows), none ended with 124 / 127");
            Assert.That(sets.Count(s => s.Full && s.Rc == 0), Is.EqualTo(69), "complete runs that finish with exit status 0 (47 of run 3 + 22 of run 4)");
            Assert.That(sets.Count(s => s.Full && s.Rc == 1), Is.EqualTo(219), "complete runs that stop with exit status 1 (the program's segmentation fault / count message is turned into exit 1 by the script)");
            Assert.That(sets.Count(s => s.Full && s.Rows.Length > 0), Is.EqualTo(65), "69 finished runs minus the 4 CLUSTAL / PHYLIP runs (no '>' rows)");
            Assert.That(sets.Count(s => s.Tree.Length > 0 || s.Hat2.Length > 0), Is.EqualTo(0), "no complete run left a tree or hat2 file (no --treeout / --distout set exists among 254 - 288; the DNA sets of 1 - 253 stop in dndfast7, mafft removes hat2): the comparison of both files is vacuous");
            Assert.That(sets.Count(s => s.Files.Count > 0), Is.EqualTo(27), "sets whose options name further files");
            Assert.That(sets.Select(s => s.Options).Distinct().Count(), Is.EqualTo(41));
            Assert.That(sets.Count(s => s.Stdout is not null), Is.EqualTo(35), "run 4 stores the raw standard output of its 35 sets");
            Assert.That(sets.Count(s => s.Stdout is { Length: > 0 }), Is.EqualTo(22));
            Assert.That(sets.Count(s => s.Stdout is not null && s.Stdout.StartsWith("CLUSTAL format alignment by MAFFT", StringComparison.Ordinal)), Is.EqualTo(2), "--clustalout");
            Assert.That(sets.Count(s => s.Stdout is not null && s.Stdout.StartsWith(" 3 ", StringComparison.Ordinal)), Is.EqualTo(2), "--phylipout (sequential PHYLIP: ' 3 150')");
            Assert.That(sets.Count(s => s.OutFiles.Count > 0), Is.EqualTo(1), "set 261: new.fa.map of --mapout");
            Assert.That(sets.Count(s => s.DirectRc == -11), Is.EqualTo(210), "dndfast7 stops with a segmentation fault (a null pointer from cutal, then putlocalhom2)");
            Assert.That(sets.Count(s => s.DirectRc == 1 && s.DirectStderr.Contains("count (fasta34 -z 3) =", StringComparison.Ordinal)), Is.EqualTo(6));
            Assert.That(sets.Count(s => s.DirectRc == 0), Is.EqualTo(72));
            Assert.That(calls.Count(c => c.Args.Split(' ').Contains("-n")), Is.EqualTo(219), "nucleotide calls of fasta34");
            Assert.That(calls.Count(c => c.Args.Split(' ').Contains("-A")), Is.EqualTo(165), "calls with -A (--fastswpair)");
        });
    }

    [TestCase("")]
    [TestCase("2")]
    public void EveryRawFasta34Output_EqualsThePort_ByteForByte(string run)
    {
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        var allCalls = CallsOf(run);
        Parallel.ForEach(allCalls, c =>
        {
            var a = c.Args.Split(' ');
            var (qt, qs) = ParseFasta(c.Query);
            var (lt, ls) = ParseFasta(c.Library);
            var run = SequenceAligner.MafftFasta34.Run(qt[0], qs[0], lt, ls, a.Contains("-n"), a.Contains("-A"), int.Parse(a[^1], CultureInfo.InvariantCulture));
            if (run.ExitCode != 0 || run.Text != c.Out)
                failures.Add($"{c.Id} ({c.Args}): exit {run.ExitCode}, {FirstDifference(run.Text, c.Out)}");
        });

        Assert.That(failures, Is.Empty, $"{failures.Count} of {allCalls.Count} calls differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
    }

    [TestCase("", 1)]
    [TestCase("2", 0)]
    public void Dndfast7_Hat2Hat3StandardErrorAndExitStatus_EqualTheProgram_OnEverySet(string run, int expectedRefused)
    {
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        int refused = 0;
        var allSets = SetsOf(run);
        Parallel.ForEach(allSets, s =>
        {
            var (titles, seqs) = ParseFasta(FastaOf(s));
            var names = titles.Select(t => "=" + t).ToList();
            SequenceAligner.MafftRun.Dndfast7Result r;
            try
            {
                r = SequenceAligner.MafftRun.Dndfast7(names, seqs.ToArray(), s.Options.Contains("--fastswpair", StringComparison.Ordinal));
            }
            catch (NotSupportedException e)
            {
                // the one set whose program reads an unassigned variable (see Dndfast7_TheProgramReadsAnUnassignedScore_IsRefused)
                if (e.Message.Contains("unassigned variable", StringComparison.Ordinal))
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
            Assert.That(failures, Is.Empty, $"{failures.Count} of {allSets.Count} sets differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
            Assert.That(refused, Is.EqualTo(expectedRefused));
        });
    }

    [Test]
    public void Dndfast7_TheProgramReadsAnUnassignedScore_IsRefused()
    {
        // set 161: the 6th sequence (40 times A) is a DNA query for fasta34 (more than 85 percent ACGTUN in more than 20 residues) in a protein run; the list of the best scores has the lines of DNA ("(  40) [f]  200 ..."),
        // sscanf( pt, "%d ..." ) fails on every one of them and the variable opt of ReadFasta34m10 was never assigned: the program stores 32767 (the upper half of a stack address under the emulator of the oracle run) and prints
        // distances of -25 ... -2728 to hat2. The value depends on the stack of the process, so the port refuses it instead of inventing a number.
        var s = Sets.Single(x => x.Id == "161");
        var names = Enumerable.Range(0, s.Seqs.Length).Select(i => "=s" + i).ToList();
        var e = Assert.Throws<NotSupportedException>(() => SequenceAligner.MafftRun.Dndfast7(names, s.Seqs, true));
        Assert.That(e!.Message, Does.Contain("unassigned variable"));
        Assert.That(s.DirectHat2, Does.Contain("-2728.583"), "the program's distances of that run");
        var guess = s.Seqs[5];
        Assert.That(guess.All(c => c == 'A') && guess.Length == 40, Is.True);
    }

    [TestCase("", 1)]
    [TestCase("2", 0)]
    public void EveryMafftRun_ReturnCodeRowsNamesTreeHat2AndStandardError_EqualTheProgram(string run, int expectedRefused)
    {
        // run "2": all 288 sets have the complete run (status, rows, names, tree, hat2 file, standard error; raw standard output and output files where run 4 stored them)
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        int refused = 0;
        var allSets = SetsOf(run);
        Parallel.ForEach(allSets, s =>
        {
            MafftProgramResult res;
            try
            {
                var files = new Dictionary<string, string>(s.Files) { ["in.fa"] = FastaOf(s) };
                res = MafftProgram.Run([.. s.Options.Split(' ', StringSplitOptions.RemoveEmptyEntries), "in.fa"], files);
            }
            catch (NotSupportedException e)
            {
                if (e.Message.Contains("unassigned variable", StringComparison.Ordinal))
                    Interlocked.Increment(ref refused);
                else
                    failures.Add($"{s.Id} ({s.Options}): {e.GetType().Name} {e.Message}");
                return;
            }

            var (titles, seqs) = ParseFasta(res.Stdout);
            string tree = res.Files.TryGetValue("in.fa.tree", out var t) ? t : "";
            string hat2 = res.Files.TryGetValue("in.fa.hat2", out var h) ? h : "";
            bool threaded = RunsOnSeveralThreads(s.Options);
            string stderr = NormaliseSegfault(res.Stderr);
            string expectedStderr = NormaliseSegfault(s.Stderr);
            if (threaded)
            {
                stderr = NormaliseThreadSchedule(stderr);
                expectedStderr = NormaliseThreadSchedule(expectedStderr);
            }

            if (s.Full && res.ExitCode != s.Rc)
                failures.Add($"{s.Id} ({s.Options}): exit status {res.ExitCode}, program {s.Rc}");
            else if (s.Full && !seqs.SequenceEqual(s.Rows))
                failures.Add($"{s.Id} ({s.Options}): rows differ");
            else if (s.Full && !titles.SequenceEqual(s.Names))
                failures.Add($"{s.Id} ({s.Options}): names differ");
            else if (s.Full && tree != s.Tree)
                failures.Add($"{s.Id} ({s.Options}): tree file differs");
            else if (s.Full && hat2 != s.Hat2)
                failures.Add($"{s.Id} ({s.Options}): hat2 file differs: {FirstDifference(hat2, s.Hat2)}");
            else if (s.Full && s.Stdout is not null && res.Stdout != s.Stdout)
                failures.Add($"{s.Id} ({s.Options}): standard output differs: {FirstDifference(res.Stdout, s.Stdout)}");
            else if (s.Full && s.OutFiles.FirstOrDefault(kv => !res.Files.TryGetValue(kv.Key, out var content) || content != kv.Value) is { Key: not null } bad)
                failures.Add($"{s.Id} ({s.Options}): output file {bad.Key} differs or is missing");
            else if (stderr != expectedStderr)
                failures.Add($"{s.Id} ({s.Options}): standard error differs: {FirstDifference(stderr, expectedStderr)}");
        });

        Assert.Multiple(() =>
        {
            Assert.That(failures, Is.Empty, $"{failures.Count} of {allSets.Count} sets differ; first: {failures.OrderBy(f => f, StringComparer.Ordinal).FirstOrDefault()}");
            Assert.That(refused, Is.EqualTo(expectedRefused));
        });
    }
}
