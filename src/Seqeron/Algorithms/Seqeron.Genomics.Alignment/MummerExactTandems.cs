using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>Output of <see cref="MummerExactTandems"/>.</summary>
/// <param name="Text">What the script writes to standard output.</param>
/// <param name="Log">What the script and its programs write to standard error.</param>
/// <param name="ExitCode">The exit status of the script (0; the installed script stops with the status 2 of <c>awk</c>).</param>
/// <param name="LeftBehind">The files the script leaves in the current directory (the installed script keeps <c>&lt;pid&gt;.tmp.matches</c>; the name is given with the process id replaced by <c>$$</c>).</param>
public sealed record MummerExactTandemsResult(string Text, string Log, int ExitCode, IReadOnlyList<string> LeftBehind);

/// <summary>
/// <c>exact-tandems</c> of MUMmer 3.23: the exact tandem repeats of a sequence with a duplicated unit of at least a given length. The script runs
/// <c>repeat-match -t -n N file | tail +3 &gt; $$.tmp.matches</c>, sorts the matches (<c>sort -k1n -k2n</c>) and feeds them to <c>tandem-repeat.awk</c>, which writes one line per repeat region (start, extent, unit length, copies).
/// The script as installed by the Debian package is broken: its last command is split over two lines (<c>awk -f $</c> and, on the next line, <c>$scriptdir/tandem-repeat.awk</c>), so <c>awk</c> looks for a file named <c>$</c>,
/// fails with status 2, and the script - run with <c>sh -e</c> - stops there without writing a table and without removing its temporary file. <see cref="Run"/> is the pipeline the script was meant to run,
/// <see cref="RunInstalled"/> the script as it behaves.
/// </summary>
public static class MummerExactTandems
{
    /// <summary>The message <c>awk</c> (gawk) of the installed script prints for the file named <c>$</c>.</summary>
    public const string InstalledAwkError = "awk: fatal: cannot open source file `$' for reading: No such file or directory\n";

    /// <summary>The intended pipeline: <c>repeat-match -t -n N file | tail +3 | sort -k1n -k2n | awk -f tandem-repeat.awk</c>.</summary>
    /// <param name="genomeFasta">The FASTA text (only the first record is used).</param>
    /// <param name="minMatchLength">The minimum length of the duplicated unit (<c>-n</c>; a value below 1 makes <c>repeat-match</c> fail and the table empty).</param>
    /// <exception cref="FormatException">The program reads uninitialised memory or aborts (see <see cref="MummerRepeatMatch"/>).</exception>
    public static MummerExactTandemsResult Run(string genomeFasta, int minMatchLength)
    {
        var (matches, log) = Matches(genomeFasta, minMatchLength);
        return new MummerExactTandemsResult(Table(matches), log, 0, []);
    }

    /// <summary>The script as installed: <c>Finding matches</c>, <c>Tandem repeats</c>, then the failure of <c>awk -f $</c>.</summary>
    /// <param name="genomeFasta">The FASTA text.</param>
    /// <param name="minMatchLength">The minimum length of the duplicated unit.</param>
    /// <exception cref="FormatException">The program reads uninitialised memory or aborts (see <see cref="MummerRepeatMatch"/>).</exception>
    public static MummerExactTandemsResult RunInstalled(string genomeFasta, int minMatchLength)
    {
        var (_, log) = Matches(genomeFasta, minMatchLength);
        return new MummerExactTandemsResult("Finding matches\nTandem repeats\n", log + InstalledAwkError, 2, ["$$.tmp.matches"]);
    }

    private static (List<string> Lines, string Log) Matches(string genomeFasta, int minMatchLength)
    {
        ArgumentNullException.ThrowIfNull(genomeFasta);
        if (minMatchLength < 1)
            return ([], "ERROR:  Illegal min match length \"" + minMatchLength.ToString(CultureInfo.InvariantCulture) + "\"\n");
        var r = MummerRepeatMatch.Run(genomeFasta, new MummerRepeatMatchOptions { MinMatchLength = minMatchLength, TandemOnly = true });
        var lines = r.Text.Split('\n', StringSplitOptions.None).ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);                         // tail +3: from the third line on
        return (lines.Skip(2).ToList(), r.Messages);
    }

    private static string Table(List<string> matches)
    {
        var ci = CultureInfo.InvariantCulture;
        var rows = matches.Select(l => (Line: l, F: l.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries))).Where(t => t.F.Length >= 3)
            .Select(t => (t.Line, Start: long.Parse(t.F[0], ci), Second: long.Parse(t.F[1], ci), Length: long.Parse(t.F[2], ci))).ToList();
        // sort -k1n -k2n: the numeric keys, then the whole line (last-resort comparison)
        rows.Sort((a, b) =>
        {
            int c = a.Start.CompareTo(b.Start);
            if (c == 0)
                c = a.Second.CompareTo(b.Second);
            return c != 0 ? c : string.CompareOrdinal(a.Line, b.Line);
        });
        var sb = new StringBuilder();
        sb.Append("Start".PadLeft(8)).Append(' ').Append("Extent".PadLeft(8)).Append(' ').Append("UnitLen".PadLeft(8)).Append(' ').Append("Copies".PadLeft(10)).Append('\n');
        bool havePrev = false;
        long prev = 0;
        foreach (var (_, start, second, length) in rows)
        {
            if (start + length < second)
                continue;
            if (havePrev && start == prev)
                continue;
            long extent = second + length - start, unit = second - start;
            sb.Append(start.ToString(ci).PadLeft(8)).Append(' ').Append(extent.ToString(ci).PadLeft(8)).Append(' ').Append(unit.ToString(ci).PadLeft(8)).Append(' ')
                .Append(CFormat.Fixed((double)extent / unit, 1).PadLeft(10)).Append('\n');
            prev = start;
            havePrev = true;
        }

        return sb.ToString();
    }
}
