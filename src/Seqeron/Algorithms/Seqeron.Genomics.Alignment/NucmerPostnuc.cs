using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="NucmerPostnuc.Run"/> = the switches of <c>postnuc</c> that <c>nucmer</c> passes (<c>-b</c>, <c>-B</c>, <c>-d</c>, <c>-e</c>, <c>-s</c>, <c>-t</c>).</summary>
public sealed record NucmerPostnucOptions
{
    /// <summary>Break (give-up) length of the extension (<c>-b</c>, <c>nucmer -b|--breaklen</c>, default 200). Values outside 1 - 10 000 are ignored with a warning, as the program does.</summary>
    public int BreakLength { get; init; } = 200;

    /// <summary>Absolute diagonal banding of the extension (<c>-B</c>, <c>nucmer --banded</c> passes the diagonal difference; 0 = none). A negative value is ignored with a warning.</summary>
    public int Banding { get; init; }

    /// <summary>Create the delta file (default); false = <c>-d</c> (<c>nucmer --nodelta</c>): output only the clusters.</summary>
    public bool Delta { get; init; } = true;

    /// <summary>Extend the alignments outward from the clusters (default); false = <c>-e</c> (<c>nucmer --noextend</c>).</summary>
    public bool Extend { get; init; } = true;

    /// <summary>Force the alignment to the ends of the sequences when within the break length (<c>-t</c>, <c>nucmer --nooptimize</c>).</summary>
    public bool ToSequenceEnd { get; init; }

    /// <summary>Do not remove shadowed alignments (<c>-s</c>, <c>nucmer --nosimplify</c>).</summary>
    public bool KeepShadowed { get; init; }
}

/// <summary>Output of <see cref="NucmerPostnuc.Run"/>.</summary>
/// <param name="Delta">The text of <c>&lt;pfx&gt;.delta</c> (null with <c>-d</c>).</param>
/// <param name="Cluster">The text of <c>&lt;pfx&gt;.cluster</c> (null unless <c>-d</c>).</param>
/// <param name="Messages">What the program writes to standard error (warnings), in order.</param>
public sealed record NucmerPostnucResult(string? Delta, string? Cluster, string Messages);

/// <summary>
/// The <c>postnuc</c> program of MUMmer 3.23 (<c>src/tigr/postnuc.cc</c>, Phillippy, TIGR): reads the clusters <c>mgaps</c> found on the concatenated reference (the <c>prenuc</c> file: the sequences joined with an
/// <c>x</c> after each), maps their coordinates back to the individual reference records, extends every cluster by the banded Smith-Waterman engine of <see cref="NucmerSwAlign"/> (the gaps between the matches of
/// a cluster, then outward off both ends until the score stops improving within the break length, fusing it with a neighbouring cluster or alignment when they meet, dropping clusters that are shadowed by an
/// earlier alignment) and writes the <c>.delta</c> file. Clusters of one (reference record, query record) pair are sorted by their first start with the <c>std::sort</c> of libstdc++
/// (<see cref="LibstdcxxSort"/>: not stable, equal starts keep the order that algorithm gives them).
/// </summary>
/// <remarks>
/// The reference and query are FASTA texts read as the program's <c>Read_String</c> does (header = first word, residues lower-cased, anything that is not a letter or <c>*</c> becomes <c>x</c>, white space skipped).
/// The mgaps text is the one <c>mgaps</c> prints (<c>&gt; id</c> / <c>&gt; id Reverse</c> header lines, <c>#</c> between clusters, then <c>start1 start2 length ...</c> lines).
/// </remarks>
public static class NucmerPostnuc
{
    private const sbyte Forward = 1;
    private const sbyte Reverse = -1;

    private sealed class FastaRecord
    {
        public string Id = string.Empty;
        public long Len;
        public byte[] Seq = [0, 0];   // 1-based, Seq[0] = Seq[Len + 1] = 0
    }

    private struct Match
    {
        public long SA, SB, Len;
    }

    private sealed class Cluster
    {
        public bool WasFused;
        public sbyte DirB;
        public List<Match> Matches = [];
    }

    private sealed class Synteny
    {
        public FastaRecord AfP = null!;
        public string BId = string.Empty;
        public long BLen = -1;
        public List<Cluster> Clusters = [];
    }

    private sealed class Alignment
    {
        public sbyte DirB;
        public long SA, SB, EA, EB;
        public List<long> Delta = [];
        public long DeltaApos;
        public long Errors, SimErrors, NonAlphas;
    }

    private enum LineType { None, Header, MatchLine }

    private sealed class State
    {
        public NucmerPostnucOptions O = null!;
        public int BreakLen = 200;
        public int Banding;
        public StringBuilder Err = new();
        public StringBuilder? Delta;
        public StringBuilder? ClusterFile;
    }

    private static bool IsSpace(int c) => c == ' ' || (c >= '\t' && c <= '\r');

    private static bool IsAlpha(int c) => (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');

    /// <summary><c>Read_String (fp, T, Size, Name, FALSE)</c>: a sequential reader of FASTA records.</summary>
    private sealed class FastaReader(string text)
    {
        private int _pos;

        public bool Read(out string name, out byte[] seq, out long len, StringBuilder err)
        {
            name = string.Empty;
            seq = [0, 0];
            len = 0;
            int ch;
            while (_pos < text.Length && text[_pos] != '>')
                _pos++;
            if (_pos >= text.Length)
                return false;
            _pos++; // the '>'

            // fgets (Line, MAX_LINE = 1024): at most 1023 characters up to and including the newline
            int start = _pos;
            int nl = text.IndexOf('\n', _pos);
            if (nl < 0 || nl - start + 1 > 1023)
                throw new FormatException("Assertion `Len > 0 && Line [Len - 1] == '\\n'' failed: the header line of a FASTA record is missing its newline or is longer than 1023 characters.");
            string line = text.Substring(start, nl - start + 1);
            _pos = nl + 1;
            // strtok (Line, " \t\n")
            int i = 0;
            while (i < line.Length && (line[i] == ' ' || line[i] == '\t' || line[i] == '\n'))
                i++;
            int j = i;
            while (j < line.Length && line[j] != ' ' && line[j] != '\t' && line[j] != '\n')
                j++;
            name = line.Substring(i, j - i);

            var bytes = new List<byte>(1024) { 0 };
            while (_pos < text.Length)
            {
                ch = text[_pos];
                if (ch == '>')
                    break;
                _pos++;
                if (IsSpace(ch))
                    continue;
                ch = ch >= 'A' && ch <= 'Z' ? ch + 32 : ch;
                if (!IsAlpha(ch) && ch != '*')
                {
                    err.Append("Unexpected character `").Append((char)ch).Append("' in string ").Append(name).Append('\n');
                    ch = 'x';
                }

                bytes.Add((byte)ch);
            }

            bytes.Add(0);
            seq = [.. bytes];
            len = bytes.Count - 2;
            return true;
        }
    }

    private static long RevC(long coord, long len) => len - coord + 1;

    /// <summary>
    /// The FASTA records as the MUMmer programs read them (<c>Read_String</c>): the first word of the header, the residues lower-cased, white space skipped, every character that is not a letter or <c>*</c> replaced by <c>x</c>.
    /// </summary>
    /// <param name="fasta">FASTA text.</param>
    /// <param name="messages">What the reader writes to standard error.</param>
    /// <exception cref="FormatException">A header line is longer than 1023 characters or has no newline (the program's assertion).</exception>
    public static IReadOnlyList<(string Id, string Sequence)> ReadFasta(string fasta, out string messages)
    {
        ArgumentNullException.ThrowIfNull(fasta);
        var err = new StringBuilder();
        var reader = new FastaReader(fasta);
        var list = new List<(string, string)>();
        while (reader.Read(out string id, out byte[] seq, out long len, err))
            list.Add((id, Encoding.Latin1.GetString(seq, 1, (int)len)));
        messages = err.ToString();
        return list;
    }

    /// <summary>
    /// <c>postnuc [-b int] [-B int] [-d] [-e] [-s] [-t] &lt;reference&gt; &lt;query&gt; &lt;pfx&gt; &lt; mgaps</c>: the delta text (or, with <c>-d</c>, the cluster text).
    /// </summary>
    /// <param name="referenceName">Reference file name as printed in the header of the output (nucmer passes an absolute path).</param>
    /// <param name="queryName">Query file name as printed in the header of the output.</param>
    /// <param name="referenceFasta">Contents of the reference FASTA file.</param>
    /// <param name="queryFasta">Contents of the query FASTA file.</param>
    /// <param name="mgaps">The output of <c>mgaps</c> (clusters in the coordinates of the <c>prenuc</c> concatenation).</param>
    /// <param name="options">The switches (null = defaults).</param>
    /// <exception cref="FormatException">The program aborts: the input cannot be parsed, no reference sequence, the query records are not in the order of the clusters ...</exception>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static NucmerPostnucResult Run(string referenceName, string queryName, string referenceFasta, string queryFasta, string mgaps, NucmerPostnucOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(referenceName);
        ArgumentNullException.ThrowIfNull(queryName);
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(mgaps);
        var st = new State { O = options ?? new NucmerPostnucOptions() };
        return RunCore(st, referenceName, queryName, referenceFasta, queryFasta, mgaps);
    }

    private static NucmerPostnucResult RunCore(State st, string refName, string qryName, string refFasta, string qryFasta, string mgaps)
    {
        var o = st.O;
        // setBreakLen / setBanding
        if (o.BreakLength < 1 || o.BreakLength > NucmerSwAlign.MaxAlignmentLength)
            st.Err.Append("WARNING: Invalid break length ").Append(o.BreakLength.ToString(CultureInfo.InvariantCulture)).Append(", ignoring\n");
        else
            st.BreakLen = o.BreakLength;
        if (o.Banding < 0)
            st.Err.Append("WARNING: Invalid banding ").Append(o.Banding.ToString(CultureInfo.InvariantCulture)).Append(", ignoring\n");
        else
            st.Banding = o.Banding;

        if (o.Delta)
            st.Delta = new StringBuilder();
        else
            st.ClusterFile = new StringBuilder();
        (o.Delta ? st.Delta! : st.ClusterFile!).Append(refName).Append(' ').Append(qryName).Append("\nNUCMER\n");

        // the reference records
        var af = new List<FastaRecord>();
        var rr = new FastaReader(refFasta);
        while (rr.Read(out string id, out byte[] seq, out long len, st.Err))
            af.Add(new FastaRecord { Id = id, Len = len, Seq = seq });
        if (af.Count <= 0)
            ParseAbort(refName);

        var qr = new FastaReader(qryFasta);
        var syntenys = new List<Synteny>();
        int currSp = -1;
        var prevLine = LineType.None;
        string idA = string.Empty, idB = string.Empty, currIdB = string.Empty;
        sbyte dirB = Forward;

        foreach (string rawLine in SplitLines(mgaps))
        {
            string line = rawLine;
            if (line.Length > 0 && line[0] == '>')
            {
                string? tok = FirstWord(line, 1);
                if (tok == null)
                    ParseAbort("stdin");
                currIdB = tok!;
                dirB = line.Contains(" Reverse", StringComparison.Ordinal) ? Reverse : Forward;
                prevLine = LineType.Header;
            }
            else if (line.Length > 0 && line[0] == '#')
            {
                prevLine = LineType.Header;
            }
            else
            {
                if (!ScanThree(line, out long sA, out long sB, out long len))
                    ParseAbort("stdin");

                int seqi = 0;
                for (; seqi < af.Count && sA > af[seqi].Len; seqi++)
                    sA -= af[seqi].Len + 1;
                if (seqi >= af.Count)
                    throw new FormatException("ERROR: A MUM was found with a start coordinate greater than\n       the sequence length, a serious error has occured.\n       Please file a bug report\n");

                if (sA + len - 1 > af[seqi].Len || sA <= 0)
                {
                    st.Err.Append("WARNING: A MUM was found extending beyond the boundry of:\n         Reference sequence '>").Append(af[seqi].Id)
                        .Append("'\n\nPlease check that the '-n' option is activated on 'mummer2'\nand try again, or file a bug report\nAttempting to continue.\n");
                    continue;
                }

                if (idA != af[seqi].Id || idB != currIdB)
                {
                    bool found = false;
                    if (idB == currIdB)
                    {
                        for (int sp = syntenys.Count - 1; sp >= 0; sp--)
                        {
                            if (syntenys[sp].AfP.Id == af[seqi].Id)
                            {
                                if (syntenys[sp].AfP.Len != af[seqi].Len)
                                    throw new FormatException("ERROR: The reference file may contain sequences with non-unique\n       header Ids, please check your input files and try again\n");
                                currSp = sp;
                                found = true;
                                break;
                            }
                        }
                    }
                    else
                    {
                        ProcessSyntenys(st, syntenys, qr);
                    }

                    idA = af[seqi].Id;
                    idB = currIdB;

                    if (!found)
                    {
                        syntenys.Add(new Synteny { AfP = af[seqi], BId = idB, BLen = -1 });
                        currSp = syntenys.Count - 1;
                    }

                    if (syntenys.Count > 0 && syntenys[currSp].Clusters.Count > 0 && syntenys[currSp].Clusters[^1].Matches.Count == 0)
                        syntenys[currSp].Clusters.RemoveAt(syntenys[currSp].Clusters.Count - 1);
                    syntenys[currSp].Clusters.Add(new Cluster { WasFused = false, DirB = dirB });
                }
                else if (prevLine == LineType.Header)
                {
                    if (syntenys.Count > 0 && syntenys[currSp].Clusters.Count > 0 && syntenys[currSp].Clusters[^1].Matches.Count == 0)
                        syntenys[currSp].Clusters.RemoveAt(syntenys[currSp].Clusters.Count - 1);
                    syntenys[currSp].Clusters.Add(new Cluster { WasFused = false, DirB = dirB });
                }

                if (len > 1)
                    syntenys[currSp].Clusters[^1].Matches.Add(new Match { SA = sA, SB = sB, Len = len });

                prevLine = LineType.MatchLine;
            }
        }

        if (syntenys.Count > 0 && syntenys[currSp].Clusters.Count > 0 && syntenys[currSp].Clusters[^1].Matches.Count == 0)
            syntenys[currSp].Clusters.RemoveAt(syntenys[currSp].Clusters.Count - 1);
        ProcessSyntenys(st, syntenys, qr);

        return new NucmerPostnucResult(st.Delta?.ToString(), st.ClusterFile?.ToString(), st.Err.ToString());
    }

    private static IEnumerable<string> SplitLines(string text)
    {
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            if (nl < 0)
            {
                yield return text.Substring(pos);
                yield break;
            }

            yield return text.Substring(pos, nl - pos);
            pos = nl + 1;
        }
    }

    private static string? FirstWord(string line, int from)
    {
        int i = from;
        while (i < line.Length && IsSpace(line[i]))
            i++;
        int j = i;
        while (j < line.Length && !IsSpace(line[j]))
            j++;
        return j > i ? line.Substring(i, j - i) : null;
    }

    /// <summary><c>sscanf (Line, "%ld %ld %ld", ...) == 3</c>.</summary>
    private static bool ScanThree(string line, out long a, out long b, out long c)
    {
        a = b = c = 0;
        int pos = 0;
        return ScanLong(line, ref pos, out a) && ScanLong(line, ref pos, out b) && ScanLong(line, ref pos, out c);
    }

    private static bool ScanLong(string s, ref int pos, out long value)
    {
        value = 0;
        while (pos < s.Length && IsSpace(s[pos]))
            pos++;
        int start = pos;
        if (pos < s.Length && (s[pos] == '+' || s[pos] == '-'))
            pos++;
        int digits = pos;
        while (pos < s.Length && s[pos] >= '0' && s[pos] <= '9')
            pos++;
        if (pos == digits)
        {
            pos = start;
            return false;
        }

        return long.TryParse(s.AsSpan(start, pos - start), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out value);
    }

    private static void ParseAbort(string s)
        => throw new FormatException("ERROR: Could not parse input from '" + s + "'. \nPlease check the filename and format, or file a bug report\n");

    private static void ProcessSyntenys(State st, List<Synteny> syntenys, FastaReader qr)
    {
        string bId = string.Empty;
        byte[] bSeq = [0, 0];
        long bLen = 0;
        for (int i = 0; i < syntenys.Count; i++)
        {
            var cur = syntenys[i];
            if (cur.Clusters.Count == 0)
                continue;

            if (i == 0 || cur.BId != syntenys[i - 1].BId)
            {
                while (qr.Read(out bId, out bSeq, out bLen, st.Err))
                    if (cur.BId == bId)
                        break;
                if (cur.BId != bId)
                    ParseAbort("Query File");
            }

            cur.BLen = bLen;
            ExtendClusters(st, cur.Clusters, cur.AfP, bId, bSeq, bLen);
        }

        FlushSyntenys(st, syntenys);
        syntenys.Clear();
    }

    private static void FlushSyntenys(State st, List<Synteny> syntenys)
    {
        var sb = st.ClusterFile;
        if (sb == null)
            return;
        foreach (var sp in syntenys)
        {
            sb.Append('>').Append(sp.AfP.Id).Append(' ').Append(sp.BId).Append(' ').Append(sp.AfP.Len.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(sp.BLen.ToString(CultureInfo.InvariantCulture)).Append('\n');
            foreach (var cp in sp.Clusters)
            {
                sb.Append(Forward.ToString(CultureInfo.InvariantCulture).PadLeft(2)).Append(' ').Append(cp.DirB.ToString(CultureInfo.InvariantCulture).PadLeft(2)).Append('\n');
                for (int k = 0; k < cp.Matches.Count; k++)
                {
                    var mp = cp.Matches[k];
                    long sb2 = cp.DirB == Forward ? mp.SB : RevC(mp.SB, sp.BLen);
                    sb.Append(mp.SA.ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append(' ').Append(sb2.ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append(' ').Append(mp.Len.ToString(CultureInfo.InvariantCulture).PadLeft(6));
                    if (k > 0)
                    {
                        var pm = cp.Matches[k - 1];
                        sb.Append((mp.SA - pm.SA - pm.Len).ToString(CultureInfo.InvariantCulture).PadLeft(6)).Append(' ').Append((mp.SB - pm.SB - pm.Len).ToString(CultureInfo.InvariantCulture).PadLeft(6)).Append('\n');
                    }
                    else
                    {
                        sb.Append("-".PadLeft(6)).Append(' ').Append("-".PadLeft(6)).Append('\n');
                    }
                }
            }
        }
    }

    private static void ExtendClusters(State st, List<Cluster> clusters, FastaRecord af, string bId, byte[] bSeqIn, long bLen)
    {
        var o = st.O;
        LibstdcxxSort.Sort(clusters, (pa, pb) => pa.Matches[0].SA < pb.Matches[0].SA);

        if (!o.Delta)
            return;

        bool targetReached = false;
        byte[] a = af.Seq;
        byte[]? bRev = null;
        byte[] b;
        long targetA, targetB;
        var bRec = new FastaRecord { Id = bId, Len = bLen, Seq = bSeqIn };

        int count = clusters.Count;
        int prevCp = 0, currCp = 0, targetCp = count;
        var alignments = new List<Alignment>();
        int currAp = 0;

        while (currCp < count)
        {
            var cc = clusters[currCp];
            if (o.Extend && !targetReached && (cc.WasFused || (!o.KeepShadowed && IsShadowedCluster(cc, alignments, currAp))))
            {
                cc.WasFused = true;
                currCp = ++prevCp;
                continue;
            }

            if (cc.DirB == Forward)
            {
                b = bSeqIn;
            }
            else
            {
                bRev ??= MummerSequence.ReverseComplement(bRec.Seq);
                b = bRev;
            }

            for (int mp = 0; mp < cc.Matches.Count; mp++)
            {
                var m = cc.Matches[mp];
                if (targetReached)
                {
                    var cur = alignments[currAp];
                    if (cur.EA != m.SA || cur.EB != m.SB)
                    {
                        if (mp >= cc.Matches.Count - 1)
                            throw new FormatException("ERROR: Target match does not exist, please\n       file a bug report\n");
                        continue;
                    }

                    cur.EA += m.Len - 1;
                    cur.EB += m.Len - 1;
                }
                else
                {
                    alignments.Add(new Alignment { SA = m.SA, SB = m.SB, EA = m.SA + m.Len - 1, EB = m.SB + m.Len - 1, DirB = cc.DirB });
                    currAp = alignments.Count - 1;

                    if (o.Extend || mp != 0)
                    {
                        int targetAp = GetReverseTargetAlignment(st, alignments, currAp);
                        if (ExtendBackward(st, alignments, currAp, targetAp, a, b))
                            currAp = targetAp;
                    }
                }

                uint mo = NucmerSwAlign.ForwardAlign;

                if (mp < cc.Matches.Count - 1)
                {
                    targetA = cc.Matches[mp + 1].SA;
                    targetB = cc.Matches[mp + 1].SB;
                    targetReached = ExtendForward(st, alignments[currAp], a, targetA, b, targetB, mo);
                }
                else if (o.Extend)
                {
                    targetA = af.Len;
                    targetB = bLen;

                    targetCp = GetForwardTargetCluster(st, clusters, currCp, ref targetA, ref targetB);
                    if (targetCp == count)
                    {
                        mo |= NucmerSwAlign.OptimalBit;
                        if (o.ToSequenceEnd)
                            mo |= NucmerSwAlign.SeqEndBit;
                    }

                    targetReached = ExtendForward(st, alignments[currAp], a, targetA, b, targetB, mo);
                }
            }

            if (targetCp == count)
                targetReached = false;

            cc.WasFused = true;

            if (!targetReached)
                currCp = ++prevCp;
            else
                currCp = targetCp;
        }

        FlushAlignments(st, alignments, af, bRec, bRev);
    }

    private static bool ExtendBackward(State st, List<Alignment> alignments, int currAp, int targetAp, byte[] a, byte[] b)
    {
        bool targetReached;
        bool overflow = false, doubleFlag = false;
        var cur = alignments[currAp];
        uint mo = NucmerSwAlign.BackwardSearch;
        long targetA, targetB;
        bool noTarget = targetAp == alignments.Count;

        if (!noTarget)
        {
            targetA = alignments[targetAp].EA;
            targetB = alignments[targetAp].EB;
        }
        else
        {
            targetA = 1;
            targetB = 1;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (cur.SA - targetA + 1 > NucmerSwAlign.MaxAlignmentLength)
        {
            targetA = cur.SA - NucmerSwAlign.MaxAlignmentLength + 1;
            overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (cur.SB - targetB + 1 > NucmerSwAlign.MaxAlignmentLength)
        {
            targetB = cur.SB - NucmerSwAlign.MaxAlignmentLength + 1;
            if (overflow)
                doubleFlag = true;
            else
                overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (st.O.ToSequenceEnd && !doubleFlag)
            mo |= NucmerSwAlign.SeqEndBit;

        targetReached = NucmerSwAlign.AlignEngine(a, cur.SA, ref targetA, b, cur.SB, ref targetB, [], mo, st.BreakLen, st.Banding);

        if (overflow || noTarget)
            targetReached = false;

        if (targetReached)
        {
            ExtendForward(st, alignments[targetAp], a, cur.SA, b, cur.SB, NucmerSwAlign.ForcedForwardAlign);
            alignments[targetAp].EA = cur.EA;
            alignments[targetAp].EB = cur.EB;
            alignments.RemoveAt(alignments.Count - 1);
        }
        else
        {
            long endA = cur.SA, endB = cur.SB;
            NucmerSwAlign.AlignEngine(a, targetA, ref endA, b, targetB, ref endB, cur.Delta, NucmerSwAlign.ForcedForwardAlign, st.BreakLen, st.Banding);
            cur.SA = targetA;
            cur.SB = targetB;
            foreach (long d in cur.Delta)
                cur.DeltaApos += d > 0 ? d : Math.Abs(d) - 1;
        }

        return targetReached;
    }

    private static bool ExtendForward(State st, Alignment cur, byte[] a, long targetA, byte[] b, long targetB, uint mo)
    {
        bool overflow = false, doubleFlag = false;
        int di = cur.Delta.Count;

        long valA = targetA - cur.EA + 1;
        long valB = targetB - cur.EB + 1;

        if (valA > NucmerSwAlign.MaxAlignmentLength)
        {
            targetA = cur.EA + NucmerSwAlign.MaxAlignmentLength - 1;
            overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (valB > NucmerSwAlign.MaxAlignmentLength)
        {
            targetB = cur.EB + NucmerSwAlign.MaxAlignmentLength - 1;
            if (overflow)
                doubleFlag = true;
            else
                overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (doubleFlag)
            mo &= ~NucmerSwAlign.SeqEndBit;

        bool targetReached = NucmerSwAlign.AlignEngine(a, cur.EA, ref targetA, b, cur.EB, ref targetB, cur.Delta, mo, st.BreakLen, st.Banding);

        if (targetReached && overflow)
            targetReached = false;

        if (di < cur.Delta.Count)
        {
            valA = (cur.EA - cur.SA + 1) - cur.DeltaApos - 1;
            cur.Delta[di] += cur.Delta[di] > 0 ? valA : -valA;
            if (cur.Delta[di] == 0 || valA < 0)
                throw new FormatException("ERROR: failed to merge alignments at position " + cur.EA.ToString(CultureInfo.InvariantCulture) + "\n       Please file a bug report\n");

            for (int k = di; k < cur.Delta.Count; k++)
                cur.DeltaApos += cur.Delta[k] > 0 ? cur.Delta[k] : Math.Abs(cur.Delta[k]) - 1;
        }

        cur.EA = targetA;
        cur.EB = targetB;
        return targetReached;
    }

    /// <summary>Returns the index of the cluster that is most likely to join with the current one (clusters.Count when none).</summary>
    private static int GetForwardTargetCluster(State st, List<Cluster> clusters, int currCp, ref long targetA, ref long targetB)
    {
        var cur = clusters[currCp];
        var last = cur.Matches[^1];
        long sA = last.SA + last.Len - 1;
        long sB = last.SB + last.Len - 1;

        long dist = targetA - sA < targetB - sB ? targetA - sA : targetB - sB;

        int cp = clusters.Count;
        for (int cip = currCp + 1; cip < clusters.Count; cip++)
        {
            var c = clusters[cip];
            if (cur.DirB == c.DirB)
            {
                long eA = c.Matches[0].SA;
                long eB = c.Matches[0].SB;

                if ((eA < sA || eB < sB) && c.Matches[^1].SA >= sA && c.Matches[^1].SB >= sB)
                {
                    for (int mip = 0; mip < c.Matches.Count && (eA < sA || eB < sB); mip++)
                    {
                        eA = c.Matches[mip].SA;
                        eB = c.Matches[mip].SB;
                    }
                }

                if (eA >= sA && eB >= sB)
                {
                    long greater, lesser;
                    if (eA - sA > eB - sB)
                    {
                        greater = eA - sA;
                        lesser = eB - sB;
                    }
                    else
                    {
                        lesser = eA - sA;
                        greater = eB - sB;
                    }

                    if (greater < st.BreakLen || lesser * 3 + (greater - lesser) * -7 >= 0)
                    {
                        cp = cip;
                        targetA = eA;
                        targetB = eB;
                        break;
                    }

                    if ((greater << 1) - lesser < dist)
                    {
                        cp = cip;
                        targetA = eA;
                        targetB = eB;
                        dist = (greater << 1) - lesser;
                    }
                }
            }
        }

        return cp;
    }

    /// <summary>Returns the index of the earlier alignment that is most likely to join with the current one (alignments.Count when none).</summary>
    private static int GetReverseTargetAlignment(State st, List<Alignment> alignments, int currAp)
    {
        var cur = alignments[currAp];
        long sA = cur.SA, sB = cur.SB;
        long dist = sA < sB ? sA : sB;

        int ap = alignments.Count;
        for (int aip = currAp - 1; aip >= 0; aip--)
        {
            var al = alignments[aip];
            if (cur.DirB == al.DirB)
            {
                long eA = al.EA, eB = al.EB;
                if (eA <= sA && eB <= sB)
                {
                    long greater, lesser;
                    if (sA - eA > sB - eB)
                    {
                        greater = sA - eA;
                        lesser = sB - eB;
                    }
                    else
                    {
                        lesser = sA - eA;
                        greater = sB - eB;
                    }

                    if (greater < st.BreakLen || lesser * 3 + (greater - lesser) * -7 >= 0)
                    {
                        ap = aip;
                        break;
                    }

                    if ((greater << 1) - lesser < dist)
                    {
                        ap = aip;
                        dist = (greater << 1) - lesser;
                    }
                }
            }
        }

        return ap;
    }

    private static bool IsShadowedCluster(Cluster cur, List<Alignment> alignments, int ap)
    {
        long sA = cur.Matches[0].SA;
        long eA = cur.Matches[^1].SA + cur.Matches[^1].Len - 1;
        long sB = cur.Matches[0].SB;
        long eB = cur.Matches[^1].SB + cur.Matches[^1].Len - 1;

        if (alignments.Count > 0)
        {
            int aip;
            for (aip = ap; aip >= 0; aip--)
            {
                var al = alignments[aip];
                if (al.DirB == cur.DirB && al.EA >= eA && al.EB >= eB && al.SA <= sA && al.SB <= sB)
                    break;
            }

            if (aip >= 0)
                return true;
        }

        return false;
    }

    private static void FlushAlignments(State st, List<Alignment> alignments, FastaRecord af, FastaRecord bf, byte[]? bRev)
    {
        var sb = st.Delta!;
        sb.Append('>').Append(af.Id).Append(' ').Append(bf.Id).Append(' ').Append(af.Len.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(bf.Len.ToString(CultureInfo.InvariantCulture)).Append('\n');

        ParseDelta(alignments, af, bf, bRev);

        foreach (var ap in alignments)
        {
            if (ap.DirB == Forward)
                sb.Append(Line7(ap.SA, ap.EA, ap.SB, ap.EB, ap));
            else
                sb.Append(Line7(ap.SA, ap.EA, RevC(ap.SB, bf.Len), RevC(ap.EB, bf.Len), ap));
            foreach (long d in ap.Delta)
                sb.Append(d.ToString(CultureInfo.InvariantCulture)).Append('\n');
            sb.Append("0\n");
        }
    }

    private static string Line7(long a, long b, long c, long d, Alignment ap)
    {
        var ci = CultureInfo.InvariantCulture;
        return string.Concat(a.ToString(ci), " ", b.ToString(ci), " ", c.ToString(ci), " ", d.ToString(ci), " ",
            ap.Errors.ToString(ci), " ", ap.SimErrors.ToString(ci), " ", ap.NonAlphas.ToString(ci), "\n");
    }

    private static void ParseDelta(List<Alignment> alignments, FastaRecord af, FastaRecord bf, byte[]? bRevIn)
    {
        byte[]? brev = bRevIn;
        foreach (var ap in alignments)
        {
            byte[] A = af.Seq;
            byte[] B = bf.Seq;
            if (ap.DirB == Reverse)
            {
                brev ??= MummerSequence.ReverseComplement(bf.Seq);
                B = brev;
            }

            long apos = ap.SA, bpos = ap.SB;
            long errors = 0, simErrors = 0, nonAlphas = 0;
            long remain = ap.EA - ap.SA + 1;

            void Score(byte c1, byte c2)
            {
                int ch1 = c1, ch2 = c2;
                if (!IsAlpha(ch1))
                {
                    ch1 = 'O';
                    nonAlphas++;
                }

                if (!IsAlpha(ch2))
                {
                    ch2 = 'O';
                    nonAlphas++;
                }

                ch1 = ch1 >= 'a' && ch1 <= 'z' ? ch1 - 32 : ch1;
                ch2 = ch2 >= 'a' && ch2 <= 'z' ? ch2 - 32 : ch2;
                if (1 > NucmerSwAlign.MatchScore(ch1, ch2))
                    simErrors++;
                if (ch1 != ch2)
                    errors++;
            }

            foreach (long dp in ap.Delta)
            {
                long delta = dp;
                int sign = delta > 0 ? 1 : -1;
                delta = Math.Abs(delta);
                long i;
                for (i = 1; i < delta; i++)
                    Score(A[apos++], B[bpos++]);

                remain -= i - 1;
                errors++;
                simErrors++;

                if (sign == 1)
                {
                    apos++;
                    remain--;
                }
                else
                {
                    bpos++;
                }
            }

            for (long i = 0; i < remain; i++)
                Score(A[apos++], B[bpos++]);

            ap.Errors = errors;
            ap.SimErrors = simErrors;
            ap.NonAlphas = nonAlphas;
        }
    }
}
