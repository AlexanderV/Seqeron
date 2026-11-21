using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The clock of a Kalign run: the time of the day printed in front of every log line (<c>get_time</c>: <c>[%F %H:%M:%S] </c> of <c>localtime_r(time())</c>) and in the MSF header (<c>%B %d, %Y %H:%M</c>), and the
/// CPU seconds of the stopwatch lines (<c>CPU Time: 0.00u 00:00:00.00 Elapsed: 00:00:00.00</c>; the installed build has no <c>times()</c>, so <c>user</c> is <c>clock()</c> and elapsed is whole seconds of <c>time()</c>).
/// The defaults read the machine clock; a test passes fixed values.
/// </summary>
public sealed class KalignRuntime
{
    /// <summary>The local time of the day.</summary>
    public Func<DateTime> Now { get; init; } = () => DateTime.Now;

    /// <summary>The CPU time of the process in seconds (<c>clock() / CLOCKS_PER_SEC</c>).</summary>
    public Func<double> CpuSeconds { get; init; } = () =>
    {
        var u = Environment.CpuUsage;
        return (u.UserTime + u.PrivilegedTime).TotalSeconds;
    };
}

/// <summary>The messages of <c>tldevel.c</c> (<c>log_message</c>, <c>warning</c>: both to standard output) and the stopwatch of <c>esl_stopwatch.c</c> of Kalign 3.3.5.</summary>
internal sealed class KalignLog(KalignRuntime runtime, StringBuilder stdout, bool quiet)
{
    private double _cpu0;
    private long _t0;
    private double _user;
    private long _elapsedSeconds;

    internal bool Quiet { get; } = quiet;

    internal KalignRuntime Runtime { get; } = runtime;

    private string Stamp() => "[" + Runtime.Now().ToString("yyyy-MM-dd HH:mm:ss", CultureInfo.InvariantCulture) + "] ";

    /// <summary><c>LOG_MSG</c>.</summary>
    internal void Log(string message) => stdout.Append(Stamp().PadLeft(22)).Append(": ").Append("LOG ".PadLeft(8)).Append(": ").Append(message).Append('\n');

    /// <summary><c>WARNING_MSG</c> in <paramref name="file"/> at <paramref name="line"/> (the program prints the path of its build tree).</summary>
    internal void Warning(string message, string file, int line)
        => stdout.Append(Stamp().PadLeft(22)).Append(": ").Append("WARNING ".PadLeft(8)).Append(": ").Append(message).Append(" (./lib/src/").Append(file).Append(" line ").Append(line.ToString(CultureInfo.InvariantCulture)).Append(")\n");

    private static long Seconds(DateTime t) => (long)Math.Floor((t - DateTime.UnixEpoch).TotalSeconds);

    /// <summary><c>START_TIMER</c>.</summary>
    internal void Start()
    {
        _t0 = Seconds(Runtime.Now());
        _cpu0 = Runtime.CpuSeconds();
    }

    /// <summary><c>STOP_TIMER</c>.</summary>
    internal void Stop()
    {
        _elapsedSeconds = Seconds(Runtime.Now()) - _t0;
        _user = Runtime.CpuSeconds() - _cpu0;
    }

    private static string TimeString(double sec)
    {
        int h = (int)(sec / 3600.0);
        int m = (int)(sec / 60.0) - h * 60;
        int s = (int)sec - h * 3600 - m * 60;
        int hs = (int)(sec * 100.0) - h * 360000 - m * 6000 - s * 100;
        return string.Create(CultureInfo.InvariantCulture, $"{h:00}:{m:00}:{s:00}.{hs:00}");
    }

    /// <summary><c>GET_TIMING</c> (<c>tl_stopwatch_Display</c>).</summary>
    internal void Timing()
        => Log(string.Create(CultureInfo.InvariantCulture, $"CPU Time: {_user:F2}u {TimeString(_user)} Elapsed: {TimeString(_elapsedSeconds)}"));
}

/// <summary>One sequence read by Kalign's alignment readers: its name and the characters after it as the program counts them (letters and punctuation are the sequence, all of them enter the letter frequencies).</summary>
internal sealed record KalignInputSequence(string Name, string Raw);

/// <summary>
/// The alignment readers and writers of Kalign 3.3.5 <c>lib/src/msa_io.c</c>: the line buffer (a line ends at its first control character), the format detection of the first 100 lines, the FASTA, MSF and
/// Clustal readers and the FASTA, MSF (with the GCG checksums and the date line) and Clustal writers.
/// </summary>
internal static class KalignMsaIo
{
    internal const int NameLen = 256;

    internal enum Format { Fasta, Msf, Clustal, Fail }

    private static bool IsCntrl(char c) => c < 32 || c == 127;

    private static bool IsSpace(char c) => c is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';

    private static bool IsAlpha(char c) => c is >= 'a' and <= 'z' or >= 'A' and <= 'Z';

    /// <summary><c>read_file_stdin</c>: the lines of the text (as <c>getline</c> splits them), each cut at its first control character.</summary>
    internal static List<string> ReadLines(string text)
    {
        var lines = new List<string>();
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            string raw = nl < 0 ? text[pos..] : text[pos..nl];
            pos = nl < 0 ? text.Length : nl + 1;
            int cut = raw.Length;
            for (int i = 0; i < raw.Length; i++)
            {
                if (IsCntrl(raw[i]))
                {
                    cut = i;
                    break;
                }
            }
            lines.Add(raw[..cut]);
        }
        return lines;
    }

    /// <summary>The "no input" test of <c>kalign_read_input</c>: the first line minus one character is zero, or there is no line.</summary>
    internal static bool IsEmptyInput(List<string> lines) => lines.Count == 0 || lines[0].Length - 1 == 0;

    /// <summary><c>detect_alignment_format</c>.</summary>
    internal static Format Detect(List<string> lines)
    {
        bool fasta = false, msf = false, clu = false;
        for (int i = 0; i < Math.Min(lines.Count, 100); i++)
        {
            string line = lines[i];
            if (line.Length > 0 && line[0] == '>')
                fasta = true;
            if (line.Contains("multiple sequence alignment", StringComparison.Ordinal) || line.Contains("CLUSTAL W", StringComparison.Ordinal) || line.Contains("CLUSTAL O", StringComparison.Ordinal))
                clu = true;
            if (line.Contains("!!AA_MULTIPLE_ALIGNMENT", StringComparison.Ordinal) || line.Contains("!!NA_MULTIPLE_ALIGNMENT", StringComparison.Ordinal) || line.Contains("MSF:", StringComparison.Ordinal))
                msf = true;
        }
        if (clu)
            return Format.Clustal;
        if (msf)
            return Format.Msf;
        return fasta ? Format.Fasta : Format.Fail;
    }

    /// <summary><c>read_fasta</c>.</summary>
    internal static List<KalignInputSequence> ReadFasta(List<string> lines)
    {
        var names = new List<string>();
        var raws = new List<StringBuilder>();
        foreach (string line in lines)
        {
            if (line.Length > 0 && line[0] == '>')
            {
                names.Add(line[1..]);
                raws.Add(new StringBuilder());
                continue;
            }
            if (line.Length == 0)
                continue;
            if (raws.Count == 0)
            {
                if (line.Any(IsAlpha))
                    throw new ArgumentException("Encountered a sequence before encountering it's name");
                throw new ArgumentException("A line of punctuation before the first sequence name (the program dereferences a null sequence).");
            }
            raws[^1].Append(line);
        }
        return [.. names.Select((n, i) => new KalignInputSequence(n, raws[i].ToString()))];
    }

    /// <summary><c>read_clu</c>: the first line is skipped; a line that does not start with a blank is the next sequence of the block (name = up to 255 characters up to the first white space), an empty line ends the block.</summary>
    internal static List<KalignInputSequence> ReadClustal(List<string> lines)
    {
        var names = new List<string>();
        var raws = new List<StringBuilder>();
        int active = 0;
        int numseq = 0;
        for (int nl = 1; nl < lines.Count; nl++)
        {
            string line = lines[nl];
            if (line.Length == 0)
            {
                active = 0;
                continue;
            }
            if (IsSpace(line[0]))
                continue;
            while (names.Count <= active)
            {
                names.Add("");
                raws.Add(new StringBuilder());
            }
            int j = 0;
            var name = new StringBuilder();
            for (int i = 0; i < line.Length; i++)
            {
                if (i == NameLen - 1 || IsSpace(line[i]))
                {
                    j = i;
                    break;
                }
                name.Append(line[i]);
            }
            // a line without any white space leaves j = 0: the name is empty and the whole line is read as residues (as in the program)
            names[active] = j == 0 ? "" : name.ToString();
            raws[active].Append(line, j, line.Length - j);
            active++;
            numseq = Math.Max(numseq, active);
        }
        return [.. Enumerable.Range(0, numseq).Select(i => new KalignInputSequence(names[i], raws[i].ToString()))];
    }

    /// <summary><c>read_msf</c>: the names of the <c>Name: ... Len:</c> lines before the <c>//</c> line, then the blocks (the first <c>strlen(name)</c> characters of a line are skipped).</summary>
    internal static List<KalignInputSequence> ReadMsf(List<string> lines)
    {
        var names = new List<string>();
        int numseq = 0;
        int li = 0;
        for (int nl = 0; nl < lines.Count; nl++)
        {
            string line = lines[nl];
            li++;
            if (line.Contains("//", StringComparison.Ordinal))
                break;
            int p = line.IndexOf("Name:", StringComparison.Ordinal);
            if (p >= 0 && line.Contains("Len:", StringComparison.Ordinal))
            {
                p += 5;
                while (p < line.Length && IsSpace(line[p]))
                    p++;
                var name = new StringBuilder();
                for (int i = 0; i < line.Length; i++)
                {
                    if (i == NameLen - 1 || p + i >= line.Length || IsSpace(line[p + i]))
                        break;
                    name.Append(line[p + i]);
                }
                while (names.Count <= numseq)
                    names.Add("");
                names[numseq] = name.ToString();
                numseq++;
            }
        }
        var raws = new List<StringBuilder>();
        int active = 0;
        for (int nl = li; nl < lines.Count; nl++)
        {
            string line = lines[nl];
            if (line.Length == 0)
            {
                active = 0;
                continue;
            }
            if (IsSpace(line[0]))
                continue;
            while (raws.Count <= active)
                raws.Add(new StringBuilder());
            int j = active < names.Count ? Math.Min(names[active].Length, NameLen) : 0;
            if (line.Length > j)
                raws[active].Append(line, j, line.Length - j);
            active++;
        }
        return [.. Enumerable.Range(0, numseq).Select(i => new KalignInputSequence(names[i], i < raws.Count ? raws[i].ToString() : ""))];
    }

    // ------------------------------------------------------------------ writers

    /// <summary><c>write_msa_fasta</c>: 60 characters per line.</summary>
    internal static string WriteFasta(IReadOnlyList<string> names, IReadOnlyList<string> rows)
    {
        var sb = new StringBuilder();
        for (int i = 0; i < rows.Count; i++)
        {
            sb.Append('>').Append(names[i]).Append('\n');
            string row = rows[i];
            for (int p = 0; p < row.Length; p += 60)
                sb.Append(row, p, Math.Min(60, row.Length - p)).Append('\n');
        }
        return sb.ToString();
    }

    private static string Cut(string name) => name.Length > NameLen ? name[..NameLen] : name;

    private static void Blocks(StringBuilder sb, IReadOnlyList<string> names, IReadOnlyList<string> rows, int alnLen, int maxName)
    {
        // the lines are collected (block, sequence) and sorted: block 0 of every sequence, an empty line ("\n" + the line feed), block 1 ...
        int blocks = Math.Max(1, (alnLen + 59) / 60);
        for (int b = 0; b < blocks; b++)
        {
            for (int i = 0; i < rows.Count; i++)
            {
                sb.Append(Cut(names[i]).PadRight(maxName + 5));
                int f = b * 60;
                sb.Append(rows[i], f, Math.Min(60, Math.Max(0, alnLen - f))).Append('\n');
            }
            sb.Append("\n\n");
        }
    }

    /// <summary><c>write_msa_clu</c>.</summary>
    internal static string WriteClustal(IReadOnlyList<string> names, IReadOnlyList<string> rows)
    {
        var sb = new StringBuilder();
        int maxName = names.Max(n => Cut(n).Length);
        sb.Append("Kalign (3.3.5) multiple sequence alignment\n\n");
        Blocks(sb, names, rows, rows[0].Length, maxName);
        return sb.ToString();
    }

    /// <summary><c>GCGchecksum</c> of the first <paramref name="len"/> characters.</summary>
    internal static int GcgChecksum(string seq, int len)
    {
        int chk = 0;
        for (int i = 0; i < len; i++)
            chk = (chk + (i % 57 + 1) * char.ToUpperInvariant(i < seq.Length ? seq[i] : '\0')) % 10000;
        return chk;
    }

    /// <summary>
    /// <c>write_msa_msf</c>. Reproduced as written: the header and the <c>Len:</c> fields print the number of RESIDUES of the first sequence (<c>sequences[0]->len</c>, never the aligned length), the checksums
    /// run over that many leading characters of the ALIGNED rows, and the header always says <c>!!NA_MULTIPLE_ALIGNMENT</c> / <c>Type: N</c> (the alphabet size compared with <c>ALPHA_defPROTEIN</c> is 5 or 23).
    /// </summary>
    internal static string WriteMsf(IReadOnlyList<string> names, IReadOnlyList<string> rows, IReadOnlyList<int> residueCounts, string? outFile, DateTime now)
    {
        var sb = new StringBuilder();
        int maxName = names.Max(n => Cut(n).Length);
        int len0 = residueCounts[0];
        int mult = 0;
        for (int i = 0; i < rows.Count; i++)
            mult = (mult + GcgChecksum(rows[i], residueCounts[i])) % 10000;
        string baseName = "stdout";
        if (outFile != null)
        {
            baseName = outFile[(outFile.LastIndexOf('/') + 1)..];
            if (baseName.Length == 0)
                throw new ArgumentException($"No filename found in: {outFile}");
        }
        string date = now.ToString("MMMM dd, yyyy HH:mm", CultureInfo.InvariantCulture);
        sb.Append("!!NA_MULTIPLE_ALIGNMENT 1.0\n\n");
        sb.Append(' ').Append(baseName).Append("  MSF: ").Append(len0.ToString(CultureInfo.InvariantCulture)).Append("  Type: N  ").Append(date).Append("  Check: ").Append(mult.ToString(CultureInfo.InvariantCulture)).Append("  ..\n\n");
        for (int i = 0; i < rows.Count; i++)
        {
            string n = Cut(names[i]);
            sb.Append(" Name: ").Append(n.PadRight(maxName)).Append("  Len:  ").Append(len0.ToString(CultureInfo.InvariantCulture).PadLeft(5))
              .Append("  Check: ").Append(GcgChecksum(rows[i], residueCounts[i]).ToString(CultureInfo.InvariantCulture).PadLeft(4)).Append("  Weight: 1.00\n");
        }
        sb.Append("\n//\n\n");
        Blocks(sb, names, rows, rows[0].Length, maxName);
        return sb.ToString();
    }
}
