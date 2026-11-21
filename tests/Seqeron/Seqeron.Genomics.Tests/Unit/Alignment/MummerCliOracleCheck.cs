using System.Globalization;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of the <c>mummer</c> command (no test framework: also run from a scratch console).</summary>
public static class MummerCliOracleCheck
{
    internal static MummerOptions ParseArgs(IReadOnlyList<string> args)
    {
        var o = new MummerOptions { ReferenceName = "reference.fa" };
        using var e = args.GetEnumerator();
        while (e.MoveNext())
        {
            switch (e.Current)
            {
                case "-mum": o = o with { Mode = MummerMatchMode.Mum }; break;
                case "-mumcand" or "-mumreference": o = o with { Mode = MummerMatchMode.MumReference }; break;
                case "-maxmatch": o = o with { Mode = MummerMatchMode.MaxMatch }; break;
                case "-n": o = o with { MatchNucleotidesOnly = true }; break;
                case "-l": o = o with { MinLength = long.Parse(Next(e), CultureInfo.InvariantCulture) }; break;
                case "-b": o = o with { Both = true }; break;
                case "-r": o = o with { OnlyReverse = true }; break;
                case "-s": o = o with { ShowString = true }; break;
                case "-c": o = o with { ShowReversePositions = true }; break;
                case "-F": o = o with { FourColumn = true }; break;
                case "-L": o = o with { ShowSequenceLengths = true }; break;
                default: throw new ArgumentException(e.Current);
            }
        }

        return o;

        static string Next(IEnumerator<string> e) => e.MoveNext() ? e.Current : throw new ArgumentException("missing value");
    }

    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll()
    {
        var failures = new List<string>();
        int rows = 0, errors = 0, multiRef = 0, strings = 0, reverse = 0, nonEmpty = 0, multiQuery = 0;
        var modes = new Dictionary<string, int>();
        foreach (var j in MummerOracleData.Rows("mummer_cli_oracle.jsonl.gz"))
        {
            rows++;
            var args = j.StrList("args");
            var qrys = j.StrList("qrys");
            var o = ParseArgs(args) with { QueryNames = Enumerable.Range(0, qrys.Count).Select(k => $"query{k}.fa").ToList() };
            string label = $"row {rows} [{string.Join(' ', args)}]";
            string? got = null, messages = null;
            try
            {
                var r = MummerProgram.Run(j.Str("ref"), qrys, o);
                got = r.Text;
                messages = r.Messages;
            }
            catch (Exception e) when (e is ArgumentException or FormatException)
            {
            }

            if (j.GetProperty("rc").GetInt32() != 0)
            {
                errors++;
                if (got is not null)
                    failures.Add(label + ": the program fails, the port returns text");
                continue;
            }

            string expected = j.Str("out");
            if (got != expected)
            {
                failures.Add($"{label}: {(got is null ? "<error>" : MummerOracleData.Diff(got, expected))}");
                continue;
            }

            if (messages!.TrimEnd('\n') != j.Str("err").TrimEnd('\n'))
                failures.Add($"{label}: stderr {MummerOracleData.Diff(messages, j.Str("err"))}");
            modes[o.Mode.ToString()] = modes.GetValueOrDefault(o.Mode.ToString()) + 1;
            if (expected.Split('\n').Length > qrys.Count * 2 + 2) nonEmpty++;
            if (j.Str("ref").Split('>').Length > 2) multiRef++;
            if (o.ShowString) strings++;
            if (o.Both || o.OnlyReverse) reverse++;
            if (qrys.Count > 1) multiQuery++;
        }

        var counts = new Dictionary<string, int> { ["rows"] = rows, ["errors"] = errors, ["nonEmpty"] = nonEmpty, ["multiRef"] = multiRef, ["strings"] = strings, ["reverse"] = reverse, ["multiQuery"] = multiQuery, ["Mum"] = modes.GetValueOrDefault("Mum"), ["MumReference"] = modes.GetValueOrDefault("MumReference"), ["MaxMatch"] = modes.GetValueOrDefault("MaxMatch") };
        return (failures, counts);
    }

}
