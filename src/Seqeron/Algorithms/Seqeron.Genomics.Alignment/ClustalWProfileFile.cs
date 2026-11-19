using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>A profile file as ClustalW reads it: the sequences and the structure information of the file (<see cref="ClustalWProfileStructure"/>), if any.</summary>
/// <param name="Sequences">The sequences (names and residues; '-' gaps).</param>
/// <param name="Structure">The first structure line or gap-penalty mask of the file; null = none (or refused, see <paramref name="StructureProblem"/>).</param>
/// <param name="StructureProblem">"bad secondary structure format" / "bad gap penalty mask format": the program printed this and ignored the structure.</param>
public sealed record ClustalWProfileFile(IReadOnlyList<ClustalWSequence> Sequences, ClustalWProfileStructure? Structure, string? StructureProblem);

/// <summary>
/// The input files of the profile commands in the two formats that can carry a secondary-structure gap-penalty mask: CLUSTAL (<c>ClustalFileParser::countSeqs / getSeq / getSecStructure</c>: structure lines
/// <c>!SS_name</c> / <c>!GM_name</c> before the sequence lines of every block) and GDE (<c>GDEFileParser::getSeq / getSecStructure</c>: <c>"SS_name</c> / <c>"GM_name</c> entries). The line scanning is the program's,
/// including how it skips blank, conservation and comment lines.
/// </summary>
public static class ClustalWProfileFileReader
{
    private const int MaxNames = 150;

    /// <summary>Reads a CLUSTAL-format file (<c>.aln</c>) as a profile: the sequences (the first block gives their number) and the structure line / mask.</summary>
    /// <param name="text">The file.</param>
    public static ClustalWProfileFile ParseClustal(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = SplitLines(text);
        int n = CountClustalSequences(lines);
        var sequences = new List<ClustalWSequence>();
        for (int s = 1; s <= n; s++)
            sequences.Add(ReadClustalSequence(lines, s));
        var (structure, problem) = ReadClustalStructure(lines);
        return new ClustalWProfileFile(sequences, structure, problem);
    }

    /// <summary>Reads a GDE-format file (<c>.gde</c>) as a profile: the sequences of the type (<c>#</c> for DNA, <c>%</c> for protein) and the <c>"SS_</c> / <c>"GM_</c> entry.</summary>
    /// <param name="text">The file.</param>
    /// <param name="dna">True for DNA (<c>#name</c> entries), false for protein (<c>%name</c>).</param>
    public static ClustalWProfileFile ParseGde(string text, bool dna)
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = SplitLines(text);
        char marker = dna ? '#' : '%';
        var sequences = new List<ClustalWSequence>();
        int at = 0;
        while (at < lines.Count)
        {
            if (lines[at].Length == 0 || lines[at][0] != marker)
            {
                at++;
                continue;
            }

            string name = GdeName(lines[at], 1);
            var residues = new System.Text.StringBuilder();
            int j = at + 1;
            for (; j < lines.Count; j++)
            {
                string line = lines[j];
                if (line.Length > 0 && (line[0] is '%' or '#' or '"'))
                    break;
                AppendResidues(residues, line);
            }

            sequences.Add(new ClustalWSequence(name, residues.ToString()));
            at = j;
        }

        return new ClustalWProfileFile(sequences, ReadGdeStructure(lines), null);
    }

    // std::getline semantics: lines end at '\n'; a final line without '\n' is a line; nothing after the last '\n'
    private static List<string> SplitLines(string text)
    {
        var lines = text.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        return lines;
    }

    // Utility::blankLine: only digits, white space, '*', ':' and '.'
    private static bool BlankLine(string line)
    {
        foreach (char c in line)
        {
            if (c == '\n')
                break;
            if (!(char.IsAsciiDigit(c) || char.IsWhiteSpace(c) || c is '*' or ':' or '.'))
                return false;
        }

        return true;
    }

    private static bool ClustalBlankLine(string line) => (line.Length > 0 && line[0] == '!') || BlankLine(line);

    private static int CountClustalSequences(List<string> lines)
    {
        int i = 0;
        while (i < lines.Count && BlankLine(lines[i]))
            i++;
        if (i >= lines.Count)
            return 0;
        i++;   // the line that ended the first loop (the title) is consumed
        while (i < lines.Count && ClustalBlankLine(lines[i]))
            i++;
        if (i >= lines.Count)
            return 0;
        int count = 1;
        for (i++; i < lines.Count; i++)
        {
            if (ClustalBlankLine(lines[i]))
                return count;
            count++;
        }

        return 0;
    }

    private static string BlankToUnderscore(string s) => s.Replace(' ', '_');

    private static void AppendResidues(System.Text.StringBuilder sb, string token)
    {
        foreach (char c in token)
        {
            if (char.IsWhiteSpace(c))
                break;
            char u = char.ToUpperInvariant(c);
            if (ClustalWUserMatrix.ResidueCodes.Contains(u, StringComparison.Ordinal))
                sb.Append(u);
        }
    }

    private static string[] Tokens(string line) => line.Split([' ', '\t', '\r', '\v', '\f'], StringSplitOptions.RemoveEmptyEntries);

    private static ClustalWSequence ReadClustalSequence(List<string> lines, int seqNum)
    {
        string name = "";
        var residues = new System.Text.StringBuilder();
        int idx = 1;   // the title line is ignored
        while (idx < lines.Count)
        {
            string line = lines[idx++];
            if (ClustalBlankLine(line))
                continue;
            for (int i = 1; i < seqNum && idx < lines.Count; i++)
                line = lines[idx++];
            var tokens = Tokens(line);
            if (tokens.Length > 0)
            {
                name = tokens[0].Length > MaxNames ? tokens[0][..MaxNames] : tokens[0];
                if (tokens.Length > 1)
                    AppendResidues(residues, tokens[1]);
            }

            // the rest of the block
            while (true)
            {
                if (idx >= lines.Count)
                    return new ClustalWSequence(BlankToUnderscore(name), residues.ToString());
                if (ClustalBlankLine(lines[idx++]))
                    break;
            }
        }

        return new ClustalWSequence(BlankToUnderscore(name), residues.ToString());
    }

    // ClustalFileParser::getSecStructure with an unlimited mask length
    private static (ClustalWProfileStructure? Structure, string? Problem) ReadClustalStructure(List<string> lines)
    {
        int idx = 0;
        string line = "";
        bool Next()
        {
            if (idx >= lines.Count)
                return false;
            line = lines[idx++];
            return true;
        }

        var secondary = new System.Text.StringBuilder();
        var mask = new System.Text.StringBuilder();
        int structPenalties = 0;   // 0 none, 1 SECST, 2 GMASK
        string name = "";
        int structIndex = 0;

        void Take(string entryLine, System.Text.StringBuilder target)
        {
            var tokens = Tokens(entryLine.Length > 4 ? entryLine[4..] : "");
            if (tokens.Length > 1)
            {
                foreach (char c in tokens[1])
                {
                    if (!char.IsWhiteSpace(c))
                        target.Append(c);
                }
            }
        }

        string EntryName(string entryLine)
        {
            var tokens = Tokens(entryLine.Length > 4 ? entryLine[4..] : "");
            string raw = tokens.Length > 0 ? tokens[0] : "";
            return BlankToUnderscore(raw.Length > MaxNames ? raw[..MaxNames] : raw);
        }

        // the title line, then the line behind it
        if (!Next())
            return (null, null);
        if (!Next())
            return (null, null);
        do
        {
            if (!Next())
                return (null, null);
        }
        while (BlankLine(line));

        int ix = -1;
        while (line.Length > 0 && line[0] == '!')
        {
            if (line.StartsWith("!SS", StringComparison.Ordinal))
            {
                ix++;
                name = EntryName(line);
                structPenalties = 1;
                structIndex = ix;
                Take(line, secondary);
            }
            else if (line.StartsWith("!GM", StringComparison.Ordinal))
            {
                ix++;
                name = EntryName(line);
                structPenalties = 2;
                structIndex = ix;
                Take(line, mask);
            }

            if (structPenalties != 0)
                break;
            if (!Next())
                return (null, null);
        }

        if (structPenalties == 0)
            return (null, null);
        while (line.Length > 0 && line[0] == '!')
        {
            if (!Next())
                return (Build(), null);
        }

        while (true)
        {
            if (line.Length == 0 || char.IsWhiteSpace(line[0]))
                break;
            if (!Next())
                return (Build(), null);
        }

        while (true)
        {
            while (BlankLine(line))
            {
                if (!Next())
                    return (Build(), null);
            }

            for (ix = 0; ix < structIndex; ix++)
            {
                if (line.Length == 0 || line[0] != '!')
                    return (null, structPenalties == 1 ? "bad secondary structure format" : "bad gap penalty mask format");
                if (!Next())
                    return (Build(), null);
            }

            if (structPenalties == 1)
            {
                if (!line.StartsWith("!SS", StringComparison.Ordinal))
                    return (null, "bad secondary structure format");
                Take(line, secondary);
            }
            else
            {
                if (!line.StartsWith("!GM", StringComparison.Ordinal))
                    return (null, "bad gap penalty mask format");
                Take(line, mask);
            }

            while (line.Length > 0 && line[0] == '!')
            {
                if (!Next())
                    return (Build(), null);
            }

            while (true)
            {
                if (line.Length == 0 || char.IsWhiteSpace(line[0]))
                    break;
                if (!Next())
                    return (Build(), null);
            }
        }

        ClustalWProfileStructure Build() => new(name, structPenalties == 1 ? secondary.ToString() : null, structPenalties == 2 ? mask.ToString() : null);
    }

    // GDEFileParser::getSeq / getSecStructure: the name is the text after the marker up to '(' or the end of the line, trailing blanks removed
    private static string GdeName(string line, int start)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = start; i < line.Length && i < start + MaxNames; i++)
        {
            if (line[i] is '(' or '\r')
                break;
            sb.Append(line[i]);
        }

        return BlankToUnderscore(sb.ToString().TrimEnd());
    }

    private static ClustalWProfileStructure? ReadGdeStructure(List<string> lines)
    {
        for (int i = 0; i < lines.Count; i++)
        {
            string line = lines[i];
            if (line.Length == 0 || line[0] != '"')
                continue;
            bool isSs = line.StartsWith("\"SS_", StringComparison.Ordinal);
            if (!isSs && !line.StartsWith("\"GM_", StringComparison.Ordinal))
                continue;
            string name = GdeName(line, 4);
            var text = new System.Text.StringBuilder();
            for (int j = i + 1; j < lines.Count; j++)
            {
                string l = lines[j];
                if (l.Length > 0 && (l[0] is '%' or '#' or '"'))
                    break;
                text.Append(l.TrimEnd('\r'));
            }

            return new ClustalWProfileStructure(name, isSs ? text.ToString() : null, isSs ? null : text.ToString());
        }

        return null;
    }
}
