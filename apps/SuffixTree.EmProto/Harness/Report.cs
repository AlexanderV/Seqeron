using System.Globalization;

namespace SuffixTree.EmProto.Harness;

/// <summary>
/// Raw-row writer of the experiments: one tab-separated line per metric, <c>R&lt;TAB&gt;experiment&lt;TAB&gt;case&lt;TAB&gt;[run&lt;TAB&gt;]metric&lt;TAB&gt;value&lt;TAB&gt;unit&lt;TAB&gt;notes</c> (the run column only with <c>--run</c>),
/// to stdout and, when given, appended to a file (the TSVs in <c>docs/refactoring/suffix-tree-upgrade/analysis/em</c> are such files, concatenated by the runner scripts).
/// </summary>
internal sealed class Report(string? path, string? run = null)
{
    private readonly object _gate = new();

    internal void Row(string experiment, string caseName, string metric, double value, string unit, string notes = "")
        => Row(experiment, caseName, metric, value.ToString("0.######", CultureInfo.InvariantCulture), unit, notes);

    internal void Row(string experiment, string caseName, string metric, long value, string unit, string notes = "")
        => Row(experiment, caseName, metric, value.ToString(CultureInfo.InvariantCulture), unit, notes);

    internal void Row(string experiment, string caseName, string metric, string value, string unit, string notes = "")
    {
        string line = run is null ? string.Join('\t', "R", experiment, caseName, metric, value, unit, notes) : string.Join('\t', "R", experiment, caseName, run, metric, value, unit, notes);
        lock (_gate)
        {
            Console.WriteLine(line);
            if (path is not null) File.AppendAllText(path, line + "\n");
        }
    }
}
