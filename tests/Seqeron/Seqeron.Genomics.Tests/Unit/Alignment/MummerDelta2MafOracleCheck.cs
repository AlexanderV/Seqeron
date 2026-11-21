using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The oracle comparison of <c>delta2maf</c> (no test framework: also run from a scratch console).</summary>
public static class MummerDelta2MafOracleCheck
{
    public static (List<string> Failures, Dictionary<string, int> Counts) CheckAll(string resource = "delta2maf_oracle.jsonl.gz")
    {
        var failures = new List<string>();
        int rows = 0, errors = 0, equal = 0, withIds = 0, sorted = 0, widths = 0, multiRecordMaf = 0;
        foreach (var j in MummerOracleData.Rows(resource))
        {
            rows++;
            var args = j.StrList("args");
            var o = new MummerDelta2MafOptions();
            var positional = new List<string>();
            using var e = args.GetEnumerator();
            while (e.MoveNext())
            {
                switch (e.Current)
                {
                    case "-q": o = o with { SortByQuery = true }; break;
                    case "-r": o = o with { SortByReference = true }; break;
                    case "-w": e.MoveNext(); o = o with { ScreenWidth = int.Parse(e.Current, System.Globalization.CultureInfo.InvariantCulture) }; widths++; break;
                    case "-x": e.MoveNext(); o = o with { MatrixType = int.Parse(e.Current, System.Globalization.CultureInfo.InvariantCulture) }; break;
                    default: positional.Add(e.Current); break;
                }
            }

            string refId = positional.Count > 1 ? positional[1] : string.Empty, qryId = positional.Count > 2 ? positional[2] : string.Empty;
            string label = $"row {rows} [{string.Join(' ', args)}]";
            string? got = null;
            try
            {
                got = MummerDelta2Maf.Run(j.Str("delta"), j.Str("ref"), j.Str("qry"), refId, qryId, o).Text;
            }
            catch (FormatException)
            {
            }

            if (j.GetProperty("rc").GetInt32() != 0)
            {
                errors++;
                if (got is not null)
                    failures.Add(label + ": the program fails, the port returns text");
                continue;
            }

            if (got != j.Str("out"))
            {
                failures.Add($"{label}: {MummerOracleData.Diff(got ?? "<error>", j.Str("out"))}");
                continue;
            }

            equal++;
            if (positional.Count > 1) withIds++;
            if (o.SortByQuery || o.SortByReference) sorted++;
            if (got.Split("a score=").Length > 3) multiRecordMaf++;
        }

        return (failures, new Dictionary<string, int> { ["rows"] = rows, ["errors"] = errors, ["equal"] = equal, ["withIds"] = withIds, ["sorted"] = sorted, ["widths"] = widths, ["manyAlignments"] = multiRecordMaf });
    }
}
