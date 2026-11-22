#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c / ajnexus.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

// seqReadNexus of ajseqread.c with the parser of ajnexus.c (ajNexusParse, nexusBlockSave, nexusCommand, the block parsers and nexusSetSequences): the blocks of a NEXUS file are read into one object, the sequences are the rows of the
// matrix of the characters (data) block under the taxa labels.
internal sealed partial class EmbossSeqin
{
    private const RegexOptions CaseC = RegexOptions.CultureInvariant | RegexOptions.IgnoreCase;
    private static readonly Regex NxCommExp = new(@"[\[]([^&][^\[\]]+)[\]]$", RegexOptions.CultureInvariant);
    private static readonly Regex NxEndCommExp = new(@"^([^\]]*)[\]]", RegexOptions.CultureInvariant);
    private static readonly Regex NxBegCommExp = new(@"[\[]([^\[\]]*)$", RegexOptions.CultureInvariant);
    private static readonly Regex NxBegExp = new(@"^\s*[Bb][Ee][Gg][Ii][Nn]\s+(\S+)\s*;\s*$", RegexOptions.CultureInvariant);
    private static readonly Regex NxEndExp = new(@"^\s*[Ee][Nn][Dd]\S*;\s*$", RegexOptions.CultureInvariant);
    private static readonly Regex NxBegCommandExp = new(@"^\s*(\S*)\s*", RegexOptions.CultureInvariant);
    private static readonly Regex NxEndCommandExp = new(@";\s*$", RegexOptions.CultureInvariant);
    private static readonly Regex NxNewlineExp = new(@"[^\n]+", RegexOptions.CultureInvariant);
    private static readonly Regex NxTaxDimExp = new(@"\bntax\s*=\s*(\d+)", CaseC);
    private static readonly Regex NxNcharExp = new(@"\bnchar\s*=\s*(\d+)", CaseC);
    private static readonly Regex NxDatatypeExp = new(@"\bdatatype\s*=\s*(\S+)", CaseC);
    private static readonly Regex NxGapExp = new(@"\bgap\s*=\s*(\S+)", CaseC);
    private static readonly Regex NxInterleaveExp = new(@"\binterleave\b", CaseC);
    private static readonly Regex NxWordExp = new(@"\S+", CaseC);

    private sealed class NexusObject
    {
        public long Ntax;
        public bool HasTaxa;
        public long TaxaNtax;
        public List<string>? TaxLabels;
        public bool HasCharacters;
        public long CharNtax, Nchar;
        public string? DataType = "standard";
        public char Gap;
        public bool Interleave;
        public List<string>? Matrix;
        public string?[]? Sequences;
    }

    // nexusGetArray: every match is added to the array; a second call on an array that exists reads memory past the end of the structure in the library (ajListstrToarrayAppend counts "array[j]" instead of "(*array)[j]")
    private static int NexusGetArray(string src, Regex exp, ref List<string>? dest)
    {
        if (dest is not null)
            throw new InvalidOperationException("The reference program appends to an array it counts with a wrong pointer (ajListstrToarrayAppend reads memory past the end of the structure).");
        var list = new List<string>();
        string tmp = src;
        for (var m = exp.Match(tmp); m.Success; m = exp.Match(tmp))
        {
            list.Add(m.Value);
            tmp = tmp[(m.Index + m.Length)..];
        }

        dest = list.Count > 0 ? list : null;
        return list.Count;
    }

    // nexusGetInt (ajStrToUint: the value is 0 when the text does not scan)
    private static bool NexusGetInt(string src, Regex exp, ref long dest)
    {
        var m = exp.Match(src);
        if (!m.Success)
            return false;
        bool ok = CStr.ToUint(m.Groups[1].Value, out uint v);
        dest = v;
        return ok;
    }

    // nexusCommand: the next command of the block records (the list is shared by all the blocks of the file)
    private static bool NexusCommand(Queue<string> list, ref string command, ref string cmdstr)
    {
        bool incommand = false;
        while (list.Count > 0)
        {
            string rdline = list.Dequeue();
            if (!incommand)
            {
                var m = NxBegCommandExp.Match(rdline);
                if (m.Success)
                {
                    command = m.Groups[1].Value;
                    cmdstr = rdline[(m.Index + m.Length)..];
                    var e = NxEndCommandExp.Match(cmdstr);
                    if (e.Success)
                    {
                        cmdstr = cmdstr[..e.Index];
                        return true;
                    }

                    incommand = true;
                }
            }
            else
            {
                var e = NxEndCommandExp.Match(rdline);
                var sb = new StringBuilder(cmdstr);
                if (e.Success)
                {
                    cmdstr = sb.Append(rdline, 0, e.Index).ToString();
                    return true;
                }

                if (sb.Length > 0)
                    sb.Append('\n');
                cmdstr = sb.Append(rdline).ToString();
            }
        }

        return false;
    }

    private static bool NxIs(string command, string name) => string.Equals(command, name, StringComparison.OrdinalIgnoreCase);

    // nexusParseTaxa
    private static bool NexusParseTaxa(NexusObject nx, Queue<string> list)
    {
        string command = "", cmdstr = "";
        nx.HasTaxa = true;
        while (NexusCommand(list, ref command, ref cmdstr))
        {
            if (NxIs(command, "dimensions"))
            {
                if (!NxTaxDimExp.IsMatch(cmdstr))
                    break;
                if (NexusGetInt(cmdstr, NxTaxDimExp, ref nx.TaxaNtax) && nx.Ntax == 0)
                    nx.Ntax = nx.TaxaNtax;
            }
            else if (NxIs(command, "taxlabels"))
            {
                if (nx.TaxaNtax == 0)
                    return false;
                int i = NexusGetArray(cmdstr, NxWordExp, ref nx.TaxLabels);
                if (i != nx.TaxaNtax)
                    return false;
            }
        }

        return nx.TaxaNtax != 0;
    }

    // nexusParseCharacters
    private static bool NexusParseCharacters(NexusObject nx, Queue<string> list)
    {
        string command = "", cmdstr = "";
        nx.HasCharacters = true;
        while (NexusCommand(list, ref command, ref cmdstr))
        {
            if (NxIs(command, "dimensions"))
            {
                if (NexusGetInt(cmdstr, NxTaxDimExp, ref nx.CharNtax))
                {
                    if (nx.Ntax == 0)
                    {
                        nx.HasTaxa = true;
                        nx.Ntax = nx.TaxaNtax = nx.CharNtax;
                    }
                }

                if (!NexusGetInt(cmdstr, NxNcharExp, ref nx.Nchar))
                    return false;
            }
            else if (NxIs(command, "format"))
            {
                var dt = NxDatatypeExp.Match(cmdstr);
                if (dt.Success)
                    nx.DataType = dt.Groups[1].Value;
                var gap = NxGapExp.Match(cmdstr);
                if (gap.Success && gap.Groups[1].Length > 0)
                    nx.Gap = gap.Groups[1].Value[0];
                else if (gap.Success)
                    nx.Gap = '\0';
                if (NxInterleaveExp.IsMatch(cmdstr))
                    nx.Interleave = true;
            }
            else if (NxIs(command, "taxlabels"))
            {
                nx.HasTaxa = true;
                int i = NexusGetArray(cmdstr, NxWordExp, ref nx.TaxLabels);
                if (i != nx.TaxaNtax)
                    return false;
            }
            else if (NxIs(command, "matrix"))
            {
                NexusGetArray(cmdstr, NxNewlineExp, ref nx.Matrix);
            }
            else if (NxIs(command, "charstatelabels") || NxIs(command, "charlabels") || NxIs(command, "statelabels"))
            {
                // the labels are kept by the library and not used by the reader (their array is not read again)
            }
        }

        return nx.Nchar != 0;
    }

    // nexusParseUnaligned: the number of taxa of the block is kept
    private static bool NexusParseUnaligned(NexusObject nx, Queue<string> list)
    {
        string command = "", cmdstr = "";
        long unNtax = 0;
        while (NexusCommand(list, ref command, ref cmdstr))
        {
            if (NxIs(command, "dimensions") && NexusGetInt(cmdstr, NxTaxDimExp, ref unNtax) && nx.Ntax == 0)
            {
                nx.HasTaxa = true;
                nx.Ntax = nx.TaxaNtax = unNtax;
            }
        }

        return true;
    }

    // nexusParseDistances
    private static bool NexusParseDistances(NexusObject nx, Queue<string> list)
    {
        string command = "", cmdstr = "";
        long dNtax = 0;
        while (NexusCommand(list, ref command, ref cmdstr))
        {
            if (NxIs(command, "dimensions"))
            {
                if (NexusGetInt(cmdstr, NxTaxDimExp, ref dNtax) && nx.Ntax == 0)
                {
                    nx.HasTaxa = true;
                    nx.Ntax = nx.TaxaNtax = dNtax;
                }
            }
            else if (NxIs(command, "taxlabels"))
            {
                nx.HasTaxa = true;
                int i = NexusGetArray(cmdstr, NxWordExp, ref nx.TaxLabels);
                if (i != nx.TaxaNtax)
                    return false;
            }
            else if (NxIs(command, "matrix"))
            {
                if (!nx.HasCharacters)
                    throw new InvalidOperationException("The reference program dereferences a null pointer (a distances matrix with no characters block).");
                NexusGetArray(cmdstr, NxNewlineExp, ref nx.Matrix);
            }
        }

        return true;
    }

    // ajStrQuoteStripAll
    private static string QuoteStripAll(string s)
    {
        if (s.Length > 0 && s[^1] == '"')
        {
            s = s[..^1];
            if (s.Length > 0 && s[0] == '"')
                s = s[1..];
            s = s.Replace("\"\"", "\"").Replace("\\\"", "\"");
        }
        else if (s.Length > 0 && s[^1] == '\'')
        {
            s = s[..^1];
            if (s.Length > 0 && s[0] == '\'')
                s = s[1..];
            s = s.Replace("''", "'").Replace("\\'", "'");
        }

        return s;
    }

    // nexusSetSequences
    private bool NexusSetSequences(NexusObject nx)
    {
        if (!nx.HasCharacters)
            return false;
        if (nx.Sequences is not null)
            return true;
        if (nx.Nchar == 0 || nx.Ntax == 0 || !nx.HasTaxa || nx.Matrix is null)
            return false;
        if (string.Equals(nx.DataType, "continuous", StringComparison.OrdinalIgnoreCase) || string.Equals(nx.DataType, "standard", StringComparison.OrdinalIgnoreCase))
            return false;
        char gapch = nx.Gap != '\0' ? nx.Gap : '.';
        var seqtab = new Dictionary<string, StringBuilder>();
        // ajTablePut of a label that is in the table already deletes the old key, which is the string of the label list: the library then reads and deletes freed memory (valgrind: Invalid read / write of size 4 in ajStrDel
        // and ajStrGetPtr, the output depends on how the heap is reused), so a repeated label is refused
        void Put(string label)
        {
            if (seqtab.ContainsKey(label))
                throw new InvalidOperationException("The reference program reads and deletes a string it freed (ajTablePut of a taxon label that is repeated: use after free, valgrind: Invalid read).");
        }

        bool havetaxa = false;
        var filled = new List<string>();
        if (nx.TaxLabels is not null)
        {
            havetaxa = true;
            foreach (string label in nx.TaxLabels)
            {
                Put(label);
                filled.Add(label);
                seqtab[label] = new StringBuilder();
            }
        }

        string? taxlabel = null;
        StringBuilder? seqstr = null;
        long itax = 0;
        for (int i = 0; i < nx.Matrix.Count; i++)
        {
            string rdline = RemoveWhiteExcess(nx.Matrix[i]);
            if (taxlabel is null || nx.Interleave || seqstr!.Length >= nx.Nchar)
            {
                var w = NxWordExp.Match(rdline);
                if (!w.Success)
                    continue;
                taxlabel = QuoteStripAll(w.Value);
                rdline = rdline[(w.Index + w.Length)..];
                if (!havetaxa)
                {
                    Put(taxlabel);
                    filled.Add(taxlabel);
                    seqtab[taxlabel] = new StringBuilder();
                    itax++;
                    if (itax >= nx.Ntax)
                        havetaxa = true;
                }

                if (!seqtab.TryGetValue(taxlabel, out seqstr))
                {
                    _context.Error($"Unknown taxon '{taxlabel}' in nexus data");
                    return false;
                }
            }

            for (var m = NxWordExp.Match(rdline); m.Success; m = NxWordExp.Match(rdline))
            {
                seqstr!.Append(m.Value);
                rdline = rdline[(m.Index + m.Length)..];
            }

            if (gapch != '-')
            {
                string replaced = seqstr!.ToString().Replace(gapch, '-');
                seqstr.Clear().Append(replaced);
            }
        }

        nx.TaxLabels = filled;
        var seqs = new string?[nx.Ntax + 1];
        for (int i = 0; i < filled.Count && i < seqs.Length; i++)
        {
            var s = seqtab[filled[i]];
            seqs[i] = s.ToString();
            if (s.Length != nx.Nchar)
            {
                _context.Error($"Nexus sequence length for '{filled[i]}' is {s.Length}, expected {nx.Nchar}\n");
                return false;
            }
        }

        nx.Sequences = seqs;
        return true;
    }

    // ajNexusParse: null when the text is not a NEXUS file, a block is not valid or a block is not closed
    private NexusObject? NexusParse()
    {
        string? first = _buff.Next();
        if (first is null || !PrefixC(first, "#NEXUS"))
        {
            _buff.Reset();
            return null;
        }

        string? line = _buff.Next();
        if (line is null)
            return null;
        bool block = false;
        int commentlevel = 0;
        var blocklist = new Queue<string>();
        var nx = new NexusObject();
        string blockname = "";
        while (line is not null)
        {
            string rdline = line.TrimEnd(' ', '\t', '\n', '\r', '\f', '\v');
            while (commentlevel > 0)
            {
                var m = NxEndCommExp.Match(rdline);
                if (!m.Success)
                    break;
                rdline = rdline[(m.Index + m.Length)..];
                --commentlevel;
            }

            for (var m = NxCommExp.Match(rdline); m.Success; m = NxCommExp.Match(rdline))
                rdline = rdline[..m.Index] + rdline[(m.Index + m.Length)..];
            for (var m = NxBegCommExp.Match(rdline); m.Success; m = NxBegCommExp.Match(rdline))
            {
                rdline = rdline[..m.Index];
                ++commentlevel;
            }

            Match b;
            if (commentlevel > 0)
            {
                // inside a comment: the line is skipped
            }
            else if (!block && (b = NxBegExp.Match(rdline)).Success)
            {
                block = true;
                blockname = b.Groups[1].Value;
            }
            else if (block && NxEndExp.IsMatch(rdline))
            {
                if (!NexusBlockSave(nx, blocklist, blockname))
                    return null;
                blockname = "";
                block = false;
            }
            else if (block)
            {
                rdline = RemoveWhiteExcess(rdline);
                if (rdline.Length > 0)
                    blocklist.Enqueue(rdline);
            }

            line = _buff.Next();
        }

        return block ? null : nx;
    }

    private static bool NexusBlockSave(NexusObject nx, Queue<string> list, string blockname)
    {
        if (NxIs(blockname, "taxa"))
            return NexusParseTaxa(nx, list);
        if (NxIs(blockname, "characters") || NxIs(blockname, "data"))
            return NexusParseCharacters(nx, list);
        if (NxIs(blockname, "unaligned"))
            return NexusParseUnaligned(nx, list);
        if (NxIs(blockname, "distances"))
            return NexusParseDistances(nx, list);
        if (NxIs(blockname, "sets") || NxIs(blockname, "assumptions") || NxIs(blockname, "codons") || NxIs(blockname, "trees") || NxIs(blockname, "notes"))
        {
            string command = "", cmdstr = "";
            while (NexusCommand(list, ref command, ref cmdstr))
            {
                // the commands of these blocks are read and not kept by the reader
            }

            return true;
        }

        return true;
    }

    private sealed class NexusState
    {
        public NexusObject Nexus = null!;
        public string?[] Names = [];
        public int Count, Nseq;
    }

    private NexusState? _nexusData;

    // seqReadNexus
    private bool ReadNexus(EmbossSeq t)
    {
        if (_nexusData is null)
        {
            _buff.SetBuffered();
            bool ok = ReadLine();
            if (!ok)
                return false;
            if (!PrefixC(_line, "#NEXUS"))
            {
                StoreReset();
                return false;
            }

            ok = ReadLine();
            while (ok && !PrefixC(_line, "#NEXUS"))
                ok = ReadLine();
            StoreReset();
            var nexus = NexusParse();
            if (nexus is null)
            {
                StoreReset();
                return false;
            }

            var st = new NexusState { Nexus = nexus, Count = 0, Nseq = nexus.HasTaxa ? (int)Math.Min(nexus.TaxaNtax, 1 << 24) : 0 };
            st.Names = new string?[st.Nseq];
            if (nexus.TaxLabels is { } labels)
            {
                for (int j = 0; j < st.Nseq; j++)
                    st.Names[j] = j < labels.Count ? labels[j] : throw new InvalidOperationException("The reference program reads past the end of the taxa labels.");
            }

            _nexusData = st;
        }

        var d = _nexusData;
        int i = d.Count;
        if (!NexusSetSequences(d.Nexus) || d.Nexus.Sequences is null)
        {
            _nexusData = null;
            return false;
        }

        var seqs = d.Nexus.Sequences;
        t.Seq.Clear();
        if (i < seqs.Length && seqs[i] is { } text)
            t.Seq.Append(text);
        if (d.Names.Length == 0 || d.Names[0] is null)
        {
            var names = d.Nexus.TaxLabels;
            for (int j = 0; j < d.Nseq; j++)
            {
                if (names is null || j >= names.Count)
                    throw new InvalidOperationException("The reference program reads past the end of the taxa labels.");
                d.Names[j] = names[j];
            }
        }

        if (i < d.Names.Length && d.Names[i] is { } name)
            t.Name = name;
        d.Count++;
        if (d.Count >= d.Nseq)
        {
            _buff.Clear(0);
            _nexusData = null;
        }

        return true;
    }
}
