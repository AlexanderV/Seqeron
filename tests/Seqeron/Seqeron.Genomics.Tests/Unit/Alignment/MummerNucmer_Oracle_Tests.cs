using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9b (B09 F199): the whole <c>nucmer</c> 3.23 pipeline - prenuc | mummer | mgaps | postnuc - against the installed program: the text of <c>out.delta</c> (or <c>out.cluster</c> with --nodelta) of
/// 434 runs, 20+ with the coords file (TestData/Alignment/nucmer_oracle.jsonl.gz, scripts/alignment/generate_nucmer_oracle.py): multi-record references and queries, both strands, repeats, N runs, IUPAC and non-letter characters,
/// every option family (--mum / --mumreference / --maxmatch, -b -c -g -l -D -d, --banded --noextend --nooptimize --nosimplify -f -r --nodelta) and the inputs on which the programs abort.
/// </summary>
[TestFixture]
public class MummerNucmer_Oracle_Tests
{
    internal static NucmerOptions Parse(IReadOnlyList<string> a)
    {
        var o = new NucmerOptions();
        var algorithms = new HashSet<string>();
        var ci = CultureInfo.InvariantCulture;
        int i = 0;
        while (i < a.Count)
        {
            string opt = a[i++];
            switch (opt)
            {
                case "--mum": o = o with { MatchMode = CommonRegionMode.MaximalUnique }; algorithms.Add("mum"); break;
                case "--maxmatch": o = o with { MatchMode = CommonRegionMode.MaximalExactMatches }; algorithms.Add("maxmatch"); break;
                case "--mumreference":
                case "--mumcand": o = o with { MatchMode = CommonRegionMode.MaximalUniqueReference }; algorithms.Add("mumreference"); break;
                case "-b": o = o with { BreakLength = int.Parse(a[i++], ci) }; break;
                case "-c": o = o with { MinCluster = long.Parse(a[i++], ci) }; break;
                case "-g": o = o with { MaxGap = long.Parse(a[i++], ci) }; break;
                case "-l": o = o with { MinMatch = int.Parse(a[i++], ci) }; break;
                case "-D": o = o with { DiagonalDifference = int.Parse(a[i++], ci) }; break;
                case "-d": o = o with { DiagonalFactor = double.Parse(a[i++], ci) }; break;
                case "--noextend": o = o with { Extend = false }; break;
                case "--nooptimize": o = o with { Optimize = false }; break;
                case "--nosimplify": o = o with { Simplify = false }; break;
                case "--banded": o = o with { Banded = true }; break;
                case "-f": o = o with { ForwardOnly = true }; break;
                case "-r": o = o with { ReverseOnly = true }; break;
                case "--nodelta": o = o with { Delta = false }; break;
                case "-o" or "--coords": o = o with { Coords = true }; break;
                default: throw new ArgumentException("option " + opt);
            }
        }

        if (algorithms.Count > 1)
            throw new ArgumentException("ERROR: Multiple matching algorithms selected");
        return o;
    }

    [Test]
    public void Nucmer_FilesEqualTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.nucmer_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int coordsRows = 0, rows = 0, errors = 0, reverseAlignments = 0, multi = 0, nodelta = 0, banded = 0, noextend = 0, nooptimize = 0, nosimplify = 0, maxmatch = 0, mum = 0, multiRef = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string reference = j.GetProperty("ref").GetString()!, query = j.GetProperty("qry").GetString()!;
            string? expected = j.GetProperty("out").ValueKind == JsonValueKind.Null ? null : j.GetProperty("out").GetString();
            string? got;
            try
            {
                var r = GenomicAnalyzer.Nucmer(reference, query, "/tmp/b09nuc/ref.fa", "/tmp/b09nuc/qry.fa", Parse(args));
                got = r.Delta ?? r.Cluster;
                if (j.TryGetProperty("coords", out var ec))
                {
                    coordsRows++;
                    if (r.Coords != ec.GetString())
                        failures.Add($"row {rows} [{string.Join(' ', args)}]: coords file differs");
                }
            }
            catch (Exception e) when (e is FormatException or ArgumentException)
            {
                got = null;
            }

            if (expected is null)
                errors++;
            else
            {
                if (args.Contains("--nodelta")) nodelta++;
                if (args.Contains("--banded")) banded++;
                if (args.Contains("--noextend")) noextend++;
                if (args.Contains("--nooptimize")) nooptimize++;
                if (args.Contains("--nosimplify")) nosimplify++;
                if (args.Contains("--maxmatch")) maxmatch++;
                if (args.Contains("--mum")) mum++;
                if (reference.Split('>').Length > 2) multiRef++;
                var lines = expected.Split('\n');
                int alignments = lines.Count(l => l.Split(' ').Length == 7);
                int headers = lines.Count(l => l.StartsWith('>'));
                if (alignments > headers) multi++;
                reverseAlignments += lines.Count(l => { var p = l.Split(' '); return p.Length == 7 && long.Parse(p[2], CultureInfo.InvariantCulture) > long.Parse(p[3], CultureInfo.InvariantCulture); });
            }

            if (got != expected)
                failures.Add($"row {rows} [{string.Join(' ', args)}]: {Diff(got ?? "<error>", expected ?? "<error>")}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(420));
        Assert.That(errors, Is.GreaterThanOrEqualTo(5));
        Assert.That(coordsRows, Is.GreaterThan(20));
        Assert.That(reverseAlignments, Is.GreaterThan(500));
        Assert.That(multi, Is.GreaterThan(150));
        Assert.That(multiRef, Is.GreaterThan(100));
        Assert.That(new[] { nodelta, banded, noextend, nooptimize, nosimplify, maxmatch, mum }, Has.All.GreaterThan(15));
    }

    private static string Diff(string a, string b)
    {
        int i = 0;
        while (i < Math.Min(a.Length, b.Length) && a[i] == b[i])
            i++;
        string Ctx(string s) => i >= s.Length ? "<end>" : s.Substring(Math.Max(0, i - 30), Math.Min(s.Length - Math.Max(0, i - 30), 80)).Replace("\n", "\\n", StringComparison.Ordinal);
        return $"{i}: ours '{Ctx(a)}' expected '{Ctx(b)}'";
    }

    [Test]
    public void LibstdcxxSort_SmallRangesAreStableAndLargeOnesDoNotChangeTheMultiset()
    {
        var small = Enumerable.Range(0, 16).Select(i => (Key: i % 3, Id: i)).ToList();
        LibstdcxxSort_Sort(small);
        Assert.That(small.Select(x => x.Id), Is.EqualTo(small.OrderBy(x => x.Key).ThenBy(x => x.Id).Select(x => x.Id)), "up to 16 elements the insertion sort is stable");
        var rnd = new Random(7);
        var big = Enumerable.Range(0, 5000).Select(i => (Key: rnd.Next(50), Id: i)).ToList();
        LibstdcxxSort_Sort(big);
        Assert.That(big.Select(x => x.Key), Is.Ordered);
        Assert.That(big.Select(x => x.Id).Order(), Is.EqualTo(Enumerable.Range(0, 5000)));
    }

    private static void LibstdcxxSort_Sort(List<(int Key, int Id)> list)
        => LibstdcxxSort.Sort(list, (a, b) => a.Key < b.Key);

    [Test]
    public void SwAlign_MatchScoreIsThreeForEqualBasesAndMinusSevenOtherwise()
    {
        Assert.That(NucmerSwAlign.MatchScore('a', 'A'), Is.EqualTo(3));
        Assert.That(NucmerSwAlign.MatchScore('c', 'g'), Is.EqualTo(-7));
        Assert.That(NucmerSwAlign.MatchScore('n', 'n'), Is.EqualTo(-7), "N never matches");
        Assert.That(NucmerSwAlign.MatchScore('*', '*'), Is.EqualTo(-7));
    }

    [Test]
    public void Postnuc_ReadsTheMgapsTextOfTheProgram()
    {
        // two matches with a one-base insertion in the query between them (the program's output for exactly this input: `postnuc r.fa q.fa out < mgaps`)
        string s = "ACGTTGCAAGCTTCGATCGGATCCATGCAGTCAAGTCCATGGATCAGTCCAGTCAATGCCATGAGGCTAACGTTAGCCATGCAAGTCCATTGACCA";
        string r = ">r\n" + s + "\n", q = ">q\n" + s[..45] + "T" + s[45..] + "\n";
        string mgaps = "> q\n       1        1     45    none      -      -\n      46       47     48    none      0      1\n";
        var res = NucmerPostnuc.Run("r.fa", "q.fa", r, q, mgaps);
        Assert.That(res.Delta, Is.EqualTo("r.fa q.fa\nNUCMER\n>r q 96 97\n1 96 1 97 1 1 0\n-46\n0\n"));
    }
}
