using System.Globalization;
using System.Text.Json;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>Runs one recorded MUMmer program invocation (TestData/Alignment/mummer_show_oracle.jsonl.gz, scripts/alignment/generate_mummer_show_oracle.py) through the ports: the command-line switches are mapped to the options exactly.</summary>
public static class MummerShowOracleRunner
{
    private static readonly CultureInfo Ci = CultureInfo.InvariantCulture;

    public const uint FakeTime = 1700000000;

    public static readonly DateTime FakeNow = new(2023, 11, 14, 22, 13, 20, DateTimeKind.Utc);

    /// <summary>Returns the program's standard output, or null when the program exits with an error.</summary>
    public static string? Run(string prog, IReadOnlyList<string> a, string? stdin, string delta, string refFasta, string qryFasta)
        => Run(prog, a, stdin, delta, refFasta, qryFasta, out _);

    /// <summary>As <see cref="Run(string, IReadOnlyList{string}, string?, string, string, string)"/>; <paramref name="files"/> receives the files of show-tiling (-t, -u, -p) by option.</summary>
    public static string? Run(string prog, IReadOnlyList<string> a, string? stdin, string delta, string refFasta, string qryFasta, out Dictionary<string, string?> files)
    {
        files = [];
        try
        {
            return Dispatch(prog, a, stdin, delta, refFasta, qryFasta, files);
        }
        catch (Exception e) when (e is FormatException or ArgumentException or InvalidOperationException)
        {
            return null;
        }
    }

    private static string Dispatch(string prog, IReadOnlyList<string> a, string? stdin, string delta, string refFasta, string qryFasta, Dictionary<string, string?> files)
    {
        switch (prog)
        {
            case "delta-filter":
            {
                var o = new MummerDeltaFilterOptions();
                var r = new Args(a);
                while (r.Next(out string opt))
                {
                    switch (opt)
                    {
                        case "-q": o = o with { QueryLis = true }; break;
                        case "-r": o = o with { ReferenceLis = true }; break;
                        case "-g": o = o with { GlobalLis = true }; break;
                        case "-1": o = o with { OneToOne = true }; break;
                        case "-m": o = o with { ManyToMany = true }; break;
                        case "-i": o = o with { MinIdentity = float.Parse(r.Value(), Ci) }; break;
                        case "-l": o = o with { MinLength = long.Parse(r.Value(), Ci) }; break;
                        case "-u": o = o with { MinUnique = float.Parse(r.Value(), Ci) }; break;
                        case "-o": o = o with { MaxOverlap = float.Parse(r.Value(), Ci) }; break;
                        case "-e": o = o with { Epsilon = float.Parse(r.Value(), Ci) }; break;
                        default: throw new ArgumentException(opt);
                    }
                }

                return MummerDeltaFilter.Run(delta, o);
            }

            case "show-coords":
            {
                var o = new MummerShowCoordsOptions { RandomSeed = FakeTime, Now = FakeNow };
                var r = new Args(a);
                while (r.Next(out string opt))
                {
                    switch (opt)
                    {
                        case "-b": o = o with { Brief = true }; break;
                        case "-B": o = o with { Btab = true }; break;
                        case "-c": o = o with { ShowCoverage = true }; break;
                        case "-d": o = o with { ShowDirection = true }; break;
                        case "-H": o = o with { PrintHeader = false }; break;
                        case "-I": o = o with { MinIdentity = float.Parse(r.Value(), Ci) }; break;
                        case "-k": o = o with { Knockout = true }; break;
                        case "-l": o = o with { ShowSeqLengths = true }; break;
                        case "-L": o = o with { MinLength = long.Parse(r.Value(), Ci) }; break;
                        case "-o": o = o with { AnnotateOverlaps = true }; break;
                        case "-q": o = o with { SortByQuery = true }; break;
                        case "-r": o = o with { SortByReference = true }; break;
                        case "-g": o = o with { Las = true }; break;
                        case "-G": o = o with { WeightedLas = true }; break;
                        case "-T": o = o with { Tabular = true }; break;
                        case "-w": o = o with { ShowWarnings = true }; break;
                        default: throw new ArgumentException(opt);
                    }
                }

                return MummerShowCoords.Run(delta, o).Text;
            }

            case "show-aligns":
            {
                var o = new MummerShowAlignsOptions();
                var r = new Args(a.Take(a.Count - 2).ToList());
                while (r.Next(out string opt))
                {
                    switch (opt)
                    {
                        case "-q": o = o with { SortByQuery = true }; break;
                        case "-r": o = o with { SortByReference = true }; break;
                        case "-w": o = o with { ScreenWidth = int.Parse(r.Value(), Ci) }; break;
                        case "-x": o = o with { MatrixType = int.Parse(r.Value(), Ci) }; break;
                        default: throw new ArgumentException(opt);
                    }
                }

                return MummerShowAligns.Run(delta, refFasta, qryFasta, a[^2], a[^1], o).Text;
            }

            case "show-snps":
            {
                var o = new MummerShowSnpsOptions();
                var r = new Args(a);
                while (r.Next(out string opt))
                {
                    switch (opt)
                    {
                        case "-C": o = o with { ShowConflict = false }; break;
                        case "-H": o = o with { PrintHeader = false }; break;
                        case "-I": o = o with { ShowIndels = false }; break;
                        case "-l": o = o with { ShowLength = true }; break;
                        case "-q": o = o with { SortQuery = true }; break;
                        case "-r": o = o with { SortReference = true }; break;
                        case "-T": o = o with { Tabular = true }; break;
                        case "-x": o = o with { Context = int.Parse(r.Value(), Ci) }; break;
                        case "-S": o = o with { SelectedAlignments = (stdin ?? string.Empty).Split('\n').ToList() }; break;
                        default: throw new ArgumentException(opt);
                    }
                }

                return MummerShowSnps.Run(delta, refFasta, qryFasta, o).Text;
            }

            case "show-diff":
            {
                var o = new MummerShowDiffOptions();
                var r = new Args(a);
                while (r.Next(out string opt))
                {
                    switch (opt)
                    {
                        case "-f": o = o with { Amos = true }; break;
                        case "-H": o = o with { PrintHeader = false }; break;
                        case "-q": o = o with { QueryDiff = true }; break;
                        case "-r": o = o with { ReferenceDiff = true }; break;
                        default: throw new ArgumentException(opt);
                    }
                }

                return MummerShowDiff.Run(delta, o).Text;
            }

            case "dnadiff":
            {
                var res = GenomicAnalyzer.Dnadiff(refFasta, qryFasta, "/tmp/b09nuc/ref.fa", "/tmp/b09nuc/qry.fa", a.Contains("-d") ? delta : null);
                files["report"] = res.Report;
                files["delta"] = a.Contains("-d") ? null : res.Delta;
                files["1delta"] = res.Delta1;
                files["mdelta"] = res.DeltaM;
                files["1coords"] = res.Coords1;
                files["mcoords"] = res.CoordsM;
                files["snps"] = res.Snps;
                files["rdiff"] = res.ReferenceDiff;
                files["qdiff"] = res.QueryDiff;
                files["unref"] = res.UnalignedReference;
                files["unqry"] = res.UnalignedQuery;
                return string.Empty;
            }

            case "show-tiling":
            {
                var o = new MummerShowTilingOptions { RandomSeed = FakeTime, Now = FakeNow, CommandLine = ["show-tiling", .. a, "out.delta"] };
                var r = new Args(a);
                while (r.Next(out string opt))
                {
                    switch (opt)
                    {
                        case "-a": o = o with { PrintAlignments = true }; break;
                        case "-c": o = o with { Circular = true }; break;
                        case "-g": o = o with { MaxGap = int.Parse(r.Value(), Ci) }; break;
                        case "-i": o = o with { MinIdentity = float.Parse(r.Value(), Ci) }; break;
                        case "-l": o = o with { MinContigLength = int.Parse(r.Value(), Ci) }; break;
                        case "-R": o = o with { RandomRepeats = true }; break;
                        case "-v": o = o with { MinCoverage = float.Parse(r.Value(), Ci) }; break;
                        case "-V": o = o with { MinCoverageDifference = float.Parse(r.Value(), Ci) }; break;
                        case "-x": o = o with { Xml = true }; break;
                        case "-t": r.Value(); o = o with { ContigList = true }; break;
                        case "-u": r.Value(); o = o with { Unusable = true }; break;
                        case "-p": r.Value(); o = o with { PseudoMolecule = true }; break;
                        default: throw new ArgumentException(opt);
                    }
                }

                var res = MummerShowTiling.Run(delta, qryFasta, o);
                if (o.ContigList) files["-t"] = res.ContigList;
                if (o.Unusable) files["-u"] = res.Unusable;
                if (o.PseudoMolecule) files["-p"] = res.PseudoMolecule;
                return res.Text;
            }

            default:
                throw new ArgumentException(prog);
        }
    }

    private sealed class Args(IReadOnlyList<string> words)
    {
        private int _i;

        public bool Next(out string opt)
        {
            if (_i < words.Count)
            {
                opt = words[_i++];
                return true;
            }

            opt = string.Empty;
            return false;
        }

        public string Value() => words[_i++];
    }

    /// <summary>Compares a file of the oracle with ours; the lines of the hash-ordered blocks of dnadiff (report, .unref, .unqry) are compared as sets of lines per blank-line separated block.</summary>
    public static bool FileEqual(string name, string? ours, string? expected, bool orderInsensitive)
    {
        if (ours == expected)
            return true;
        if (!orderInsensitive || ours is null || expected is null)
            return false;
        static string Canon(string t) => string.Join("\n\n", t.Split("\n\n").Select(b => string.Join("\n", b.Split('\n').Order(StringComparer.Ordinal))));
        return Canon(ours) == Canon(expected);
    }

    public static string Diff(string a, string b)
    {
        int i = 0;
        while (i < Math.Min(a.Length, b.Length) && a[i] == b[i])
            i++;
        string Ctx(string s) => i >= s.Length ? "<end>" : s.Substring(Math.Max(0, i - 40), Math.Min(s.Length - Math.Max(0, i - 40), 100)).Replace("\n", "\\n", StringComparison.Ordinal);
        return $"{i}: ours '{Ctx(a)}' expected '{Ctx(b)}'";
    }

    public static JsonElement Parse(string line) => JsonDocument.Parse(line).RootElement;
}
