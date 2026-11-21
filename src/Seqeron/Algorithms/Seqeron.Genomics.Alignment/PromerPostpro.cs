using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="PromerPostpro.Run"/> = the switches of <c>postpro</c> that <c>promer</c> passes (<c>-b</c>, <c>-d</c>, <c>-e</c>, <c>-t</c>, <c>-x</c>).</summary>
public sealed record PromerPostproOptions
{
    /// <summary>Break (give-up) length of the extension in amino acids (<c>-b</c>, <c>promer -b|--breaklen</c>, default 60). Values outside 1 - 10 000 are ignored with a warning.</summary>
    public int BreakLength { get; init; } = 60;

    /// <summary>Create the delta file (default); false = <c>-d</c> (<c>promer --nodelta</c>): output only the clusters.</summary>
    public bool Delta { get; init; } = true;

    /// <summary>Extend the alignments outward from the clusters (default); false = <c>-e</c> (<c>promer --noextend</c>).</summary>
    public bool Extend { get; init; } = true;

    /// <summary>Force the alignment to the ends of the sequences when within the break length (<c>-t</c>, <c>promer --nooptimize</c>).</summary>
    public bool ToSequenceEnd { get; init; }

    /// <summary>The alignment matrix (<c>-x</c>, <c>promer -x|--matrix</c>): 1 = BLOSUM 45, 2 = BLOSUM 62 (default), 3 = BLOSUM 80; 0 is refused with a warning (the nucleotide matrix) and 2 is used, any other value is ignored with a warning.</summary>
    public int Matrix { get; init; } = 2;
}

/// <summary>
/// The <c>postpro</c> program of MUMmer 3.23 (<c>src/tigr/postpro.cc</c>, Phillippy, TIGR): reads the clusters <c>mgaps</c> found on the amino-acid concatenation of the six-frame translations of the reference
/// (<see cref="PromerPrepro"/>), maps their coordinates back to the individual reference sequences and reading frames, extends every cluster in amino-acid space with the affine Smith-Waterman engine
/// of <see cref="NucmerSwAlign"/> (BLOSUM 45 / 62 / 80, break length in amino acids; the same routines as <see cref="NucmerPostnuc"/>, shared in <c>MummerExtension</c>) and writes the <c>.delta</c> file whose coordinates
/// reference the original DNA while the delta values count amino acids (one delta = three nucleotides).
/// </summary>
public static class PromerPostpro
{
    private sealed class FastaRecord
    {
        public string Id = string.Empty;
        public long Len;
        public byte[] Seq = [0, 0];   // 1-based
    }

    private struct Match
    {
        public long SA, SB, Len;
    }

    private sealed class Cluster
    {
        public bool WasFused;
        public int FrameA, FrameB;
        public List<Match> Matches = [];
    }

    private sealed class Synteny
    {
        public FastaRecord AfP = null!;
        public string BId = string.Empty;
        public long BLen = -1;
        public List<Cluster> Clusters = [];
    }

    private sealed class Alignment : MummerAlignment
    {
        public int FrameA, FrameB;
    }

    private enum LineType { None, Header, MatchLine }

    private sealed class State
    {
        public PromerPostproOptions O = null!;
        public int BreakLen = 60;
        public int Matrix = 2;
        public StringBuilder Err = new();
        public StringBuilder? Delta;
        public StringBuilder? ClusterFile;

        public MummerExtension.Settings Extension => new(O.ToSequenceEnd, BreakLen, 0, Matrix);
    }

    /// <summary><c>postpro [-b int] [-d] [-e] [-t] [-x type] &lt;reference&gt; &lt;query&gt; &lt;pfx&gt; &lt; mgaps</c>: the delta text (or, with <c>-d</c>, the cluster text).</summary>
    /// <param name="referenceName">Reference file name as printed in the first line of the output.</param>
    /// <param name="queryName">Query file name as printed in the first line of the output.</param>
    /// <param name="referenceFasta">Contents of the (DNA) reference FASTA file.</param>
    /// <param name="queryFasta">Contents of the (DNA) query FASTA file.</param>
    /// <param name="mgaps">The output of <c>mgaps</c> run over the output of <c>mummer</c> on the translated sequences (query headers <c>&gt; id.frame</c>).</param>
    /// <param name="options">The switches (null = defaults).</param>
    /// <exception cref="FormatException">The program aborts (no reference sequence, unparsable input, a failed assertion, the query records are not in the order of the clusters ...).</exception>
    public static NucmerPostnucResult Run(string referenceName, string queryName, string referenceFasta, string queryFasta, string mgaps, PromerPostproOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(referenceName);
        ArgumentNullException.ThrowIfNull(queryName);
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(mgaps);
        var o = options ?? new PromerPostproOptions();
        var st = new State { O = o };

        // setMatrixType (BLOSUM62); setBreakLen (60) are the defaults; promer passes -x before -b, so the matrix warnings come first
        if (o.Matrix < 0 || o.Matrix > 3)
            st.Err.Append("WARNING: Invalid matrix type ").Append(o.Matrix.ToString(CultureInfo.InvariantCulture)).Append(", ignoring\n");
        else
        {
            st.Matrix = o.Matrix;
            if (st.Matrix == 0)
            {
                st.Err.Append("WARNING: invalid matrix type 0, ignoring\n");
                st.Matrix = 2;
            }
        }

        if (o.BreakLength < 1 || o.BreakLength > NucmerSwAlign.MaxAlignmentLength)
            st.Err.Append("WARNING: Invalid break length ").Append(o.BreakLength.ToString(CultureInfo.InvariantCulture)).Append(", ignoring\n");
        else
            st.BreakLen = o.BreakLength;

        if (o.Delta)
            st.Delta = new StringBuilder();
        else
            st.ClusterFile = new StringBuilder();
        (o.Delta ? st.Delta! : st.ClusterFile!).Append(referenceName).Append(' ').Append(queryName).Append("\nPROMER\n");

        var af = new List<FastaRecord>();
        var rr = new MummerFastaReader(referenceFasta);
        while (rr.Read(out string id, out byte[] seq, out long len, st.Err))
            af.Add(new FastaRecord { Id = id, Len = len, Seq = seq });
        if (af.Count <= 0)
            ParseAbort(referenceName);

        var qr = new MummerFastaReader(queryFasta);
        var syntenys = new List<Synteny>();
        int currSp = -1;
        var prevLine = LineType.None;
        string idA = string.Empty, idB = string.Empty, currIdB = string.Empty;
        int frameA = 0, frameB = 0, currFrameB = -1;

        foreach (string line in SplitLines(mgaps))
        {
            if (line.Length > 0 && line[0] == '>')
            {
                string? tok = FirstWord(line, 1);
                if (tok == null)
                    ParseAbort("stdin");
                currIdB = tok!;
                // sscanf (CurrIdB + strlen - 1, "%d", &CurrFrameB): the last character as a number (unchanged when it is not a digit)
                if (char.IsAsciiDigit(currIdB[^1]))
                    currFrameB = currIdB[^1] - '0';
                currIdB = currIdB.Length >= 2 ? currIdB[..^2] : string.Empty;
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

                // map the coordinate in the concatenated frames back to its sequence: (len + 1) * 2 amino acids per sequence including the mask characters
                int seqi = 0;
                for (; seqi < af.Count && sA > RefLen(af[seqi].Len); seqi++)
                    sA -= RefLen(af[seqi].Len);
                if (seqi >= af.Count)
                    throw new FormatException("Assertion `Seqi < As' failed: a match lies beyond the reference sequences.");

                int i = 0;
                for (; sA > (af[seqi].Len - (i % 3)) / 3 + 1; i++)
                    sA -= (af[seqi].Len - (i % 3)) / 3 + 1;
                int currFrameA = i + 1;
                if (currFrameA <= 0 || currFrameA >= 7)
                    throw new FormatException("Assertion `CurrFrameA > 0 && CurrFrameA < 7' failed.");

                if (sA + len - 1 > (af[seqi].Len - ((currFrameA - 1) % 3)) / 3 + 1 || sA <= 0)
                {
                    st.Err.Append("\nWARNING: A MUM was found extending beyond the boundary of:\n         Reference sequence '>").Append(af[seqi].Id).Append("', frame ")
                        .Append(currFrameA.ToString(CultureInfo.InvariantCulture)).Append("\nPlease file a bug report\nAttempting to continue.\n");
                    continue;
                }

                if (sA + len - 1 > RefLen(af[seqi].Len) || sA <= 0)
                {
                    st.Err.Append("\nWARNING: A MUM was found extending beyond the boundary of:\n         Reference sequence '>").Append(af[seqi].Id).Append("'\nPlease file a bug report\nAttempting to continue.\n");
                    continue;
                }

                if (idA != af[seqi].Id || idB != currIdB || currFrameA != frameA || currFrameB != frameB)
                {
                    bool found = false;
                    if (idB == currIdB)
                    {
                        for (int sp = syntenys.Count - 1; sp >= 0; sp--)
                        {
                            if (syntenys[sp].AfP.Id == af[seqi].Id)
                            {
                                if (syntenys[sp].AfP.Len != af[seqi].Len)
                                    throw new FormatException("Assertion `Sp->AfP->len == Af[Seqi].len' failed: the reference file contains sequences with non-unique header ids.");
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
                    frameA = currFrameA;
                    frameB = currFrameB;

                    if (!found)
                    {
                        syntenys.Add(new Synteny { AfP = af[seqi], BId = idB, BLen = -1 });
                        currSp = syntenys.Count - 1;
                    }

                    PopEmptyCluster(syntenys, currSp);
                    syntenys[currSp].Clusters.Add(new Cluster { WasFused = false, FrameA = frameA, FrameB = frameB });
                }
                else if (prevLine == LineType.Header)
                {
                    PopEmptyCluster(syntenys, currSp);
                    syntenys[currSp].Clusters.Add(new Cluster { WasFused = false, FrameA = frameA, FrameB = frameB });
                }

                if (len > 1)
                    syntenys[currSp].Clusters[^1].Matches.Add(new Match { SA = sA, SB = sB, Len = len });

                prevLine = LineType.MatchLine;
            }
        }

        PopEmptyCluster(syntenys, currSp);
        ProcessSyntenys(st, syntenys, qr);
        return new NucmerPostnucResult(st.Delta?.ToString(), st.ClusterFile?.ToString(), st.Err.ToString());
    }

    private static void PopEmptyCluster(List<Synteny> syntenys, int currSp)
    {
        if (syntenys.Count > 0 && syntenys[currSp].Clusters.Count > 0 && syntenys[currSp].Clusters[^1].Matches.Count == 0)
            syntenys[currSp].Clusters.RemoveAt(syntenys[currSp].Clusters.Count - 1);
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
        while (i < line.Length && MummerFastaReader.IsSpace(line[i]))
            i++;
        int j = i;
        while (j < line.Length && !MummerFastaReader.IsSpace(line[j]))
            j++;
        return j > i ? line.Substring(i, j - i) : null;
    }

    private static bool ScanThree(string line, out long a, out long b, out long c)
    {
        a = b = c = 0;
        int pos = 0;
        return ScanLong(line, ref pos, out a) && ScanLong(line, ref pos, out b) && ScanLong(line, ref pos, out c);
    }

    private static bool ScanLong(string s, ref int pos, out long value)
    {
        value = 0;
        while (pos < s.Length && MummerFastaReader.IsSpace(s[pos]))
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

    /// <summary><c>refLen</c>: the length of the concatenated amino-acid frames of a sequence including the appended mask character of each frame.</summary>
    private static long RefLen(long ntLen) => (ntLen + 1) << 1;

    /// <summary><c>transC</c>: an amino-acid coordinate of a frame as a nucleotide coordinate of the forward strand.</summary>
    private static long TransC(long coord, int frame, long len)
        => frame > 3 ? len - (coord * 3 - (3 - (frame - 3))) + 1 : coord * 3 - (3 - frame);

    private static void ProcessSyntenys(State st, List<Synteny> syntenys, MummerFastaReader qr)
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
                {
                    st.Err.Append(cur.BId).Append('\n');
                    ParseAbort("Query File");
                }
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
        var ci = CultureInfo.InvariantCulture;
        foreach (var sp in syntenys)
        {
            sb.Append('>').Append(sp.AfP.Id).Append(' ').Append(sp.BId).Append(' ').Append(sp.AfP.Len.ToString(ci)).Append(' ').Append(sp.BLen.ToString(ci)).Append('\n');
            foreach (var cp in sp.Clusters)
            {
                sb.Append((cp.FrameA > 3 ? -(cp.FrameA - 3) : cp.FrameA).ToString(ci).PadLeft(2)).Append(' ').Append((cp.FrameB > 3 ? -(cp.FrameB - 3) : cp.FrameB).ToString(ci).PadLeft(2)).Append('\n');
                for (int k = 0; k < cp.Matches.Count; k++)
                {
                    var mp = cp.Matches[k];
                    sb.Append(TransC(mp.SA, cp.FrameA, sp.AfP.Len).ToString(ci).PadLeft(8)).Append(' ').Append(TransC(mp.SB, cp.FrameB, sp.BLen).ToString(ci).PadLeft(8)).Append(' ').Append((mp.Len * 3).ToString(ci).PadLeft(6));
                    if (k > 0)
                    {
                        var pm = cp.Matches[k - 1];
                        sb.Append(((mp.SA - pm.SA - pm.Len) * 3).ToString(ci).PadLeft(6)).Append(' ').Append(((mp.SB - pm.SB - pm.Len) * 3).ToString(ci).PadLeft(6)).Append('\n');
                    }
                    else
                    {
                        sb.Append("-".PadLeft(6)).Append(' ').Append("-".PadLeft(6)).Append('\n');
                    }
                }
            }
        }
    }

    private static byte[] Translation(State st, FastaRecord rec, int frame, out long length)
    {
        var (t, n) = MummerTranslate.Translate(rec.Seq, rec.Len, frame, st.Err);
        length = n;
        return t;
    }

    private static void ExtendClusters(State st, List<Cluster> clusters, FastaRecord af, string bId, byte[] bSeq, long bLen)
    {
        var o = st.O;
        LibstdcxxSort.Sort(clusters, (pa, pb) => pa.Matches[0].SA < pb.Matches[0].SA);

        if (!o.Delta)
            return;

        var bf = new FastaRecord { Id = bId, Len = bLen, Seq = bSeq };
        bool targetReached = false;
        var a = new byte[7][];
        var b = new byte[7][];
        var aLen = new long[7];
        var bLenArr = new long[7];
        for (int i = 0; i < 6; i++)
        {
            aLen[i + 1] = (af.Len - (i % 3)) / 3;
            bLenArr[i + 1] = (bf.Len - (i % 3)) / 3;
        }

        int count = clusters.Count;
        int prevCp = 0, currCp = 0, targetCp = count;
        var alignments = new List<Alignment>();
        int currAp = 0;

        while (currCp < count)
        {
            var cc = clusters[currCp];
            if (o.Extend && !targetReached && (cc.WasFused || IsShadowedCluster(cc, alignments, currAp)))
            {
                cc.WasFused = true;
                currCp = ++prevCp;
                continue;
            }

            a[cc.FrameA] ??= Translation(st, af, cc.FrameA, out aLen[cc.FrameA]);
            b[cc.FrameB] ??= Translation(st, bf, cc.FrameB, out bLenArr[cc.FrameB]);

            for (int mp = 0; mp < cc.Matches.Count; mp++)
            {
                var m = cc.Matches[mp];
                if (targetReached)
                {
                    var cur = alignments[currAp];
                    if (cur.EA != m.SA || cur.EB != m.SB)
                    {
                        if (mp >= cc.Matches.Count - 1)
                            throw new FormatException("Assertion `Mp < CurrCp->matches.end( ) - 1' failed: the target match does not exist.");
                        continue;
                    }

                    cur.EA += m.Len - 1;
                    cur.EB += m.Len - 1;
                }
                else
                {
                    alignments.Add(new Alignment { SA = m.SA, SB = m.SB, EA = m.SA + m.Len - 1, EB = m.SB + m.Len - 1, FrameA = cc.FrameA, FrameB = cc.FrameB });
                    currAp = alignments.Count - 1;

                    if (o.Extend || mp != 0)
                    {
                        int targetAp = GetReverseTargetAlignment(st, alignments, currAp);
                        if (MummerExtension.ExtendBackward(st.Extension, alignments, currAp, targetAp, a[cc.FrameA], b[cc.FrameB]))
                            currAp = targetAp;
                    }
                }

                uint mo = NucmerSwAlign.ForwardAlign;
                long targetA, targetB;

                if (mp < cc.Matches.Count - 1)
                {
                    targetA = cc.Matches[mp + 1].SA;
                    targetB = cc.Matches[mp + 1].SB;
                    targetReached = MummerExtension.ExtendForward(st.Extension, alignments[currAp], a[cc.FrameA], targetA, b[cc.FrameB], targetB, mo);
                }
                else if (o.Extend)
                {
                    targetA = aLen[cc.FrameA];
                    targetB = bLenArr[cc.FrameB];

                    targetCp = GetForwardTargetCluster(st, clusters, currCp, ref targetA, ref targetB);
                    if (targetCp == count)
                    {
                        mo |= NucmerSwAlign.OptimalBit;
                        if (o.ToSequenceEnd)
                            mo |= NucmerSwAlign.SeqEndBit;
                    }

                    targetReached = MummerExtension.ExtendForward(st.Extension, alignments[currAp], a[cc.FrameA], targetA, b[cc.FrameB], targetB, mo);
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

        FlushAlignments(st, alignments, af, bf);
    }

    private static void FlushAlignments(State st, List<Alignment> alignments, FastaRecord af, FastaRecord bf)
    {
        var sb = st.Delta!;
        var ci = CultureInfo.InvariantCulture;
        sb.Append('>').Append(af.Id).Append(' ').Append(bf.Id).Append(' ').Append(af.Len.ToString(ci)).Append(' ').Append(bf.Len.ToString(ci)).Append('\n');

        ParseDelta(st, alignments, af, bf);

        foreach (var ap in alignments)
        {
            sb.Append(TransC(ap.SA, ap.FrameA, af.Len).ToString(ci)).Append(' ')
                .Append((TransC(ap.EA, ap.FrameA, af.Len) + (ap.FrameA > 3 ? -2 : 2)).ToString(ci)).Append(' ')
                .Append(TransC(ap.SB, ap.FrameB, bf.Len).ToString(ci)).Append(' ')
                .Append((TransC(ap.EB, ap.FrameB, bf.Len) + (ap.FrameB > 3 ? -2 : 2)).ToString(ci)).Append(' ')
                .Append(ap.Errors.ToString(ci)).Append(' ').Append(ap.SimErrors.ToString(ci)).Append(' ').Append(ap.NonAlphas.ToString(ci)).Append('\n');
            foreach (long d in ap.Delta)
                sb.Append(d.ToString(ci)).Append('\n');
            sb.Append("0\n");
        }
    }

    private static void ParseDelta(State st, List<Alignment> alignments, FastaRecord af, FastaRecord bf)
    {
        var a = new byte[7][];
        var b = new byte[7][];
        foreach (var ap in alignments)
        {
            a[ap.FrameA] ??= Translation(st, af, ap.FrameA, out _);
            b[ap.FrameB] ??= Translation(st, bf, ap.FrameB, out _);
            byte[] aa = a[ap.FrameA], bb = b[ap.FrameB];

            long apos = ap.SA, bpos = ap.SB;
            long errors = 0, simErrors = 0, nonAlphas = 0;
            long remain = ap.EA - ap.SA + 1;

            void Score(byte c1, byte c2)
            {
                int ch1 = c1, ch2 = c2;
                if (!MummerFastaReader.IsAlpha(ch1))
                {
                    ch1 = 'O';
                    nonAlphas++;
                }

                if (!MummerFastaReader.IsAlpha(ch2))
                {
                    ch2 = 'O';
                    nonAlphas++;
                }

                ch1 = ch1 >= 'a' && ch1 <= 'z' ? ch1 - 32 : ch1;
                ch2 = ch2 >= 'a' && ch2 <= 'z' ? ch2 - 32 : ch2;
                if (1 > NucmerSwAlign.MatchScore(ch1, ch2, st.Matrix))
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
                    Score(aa[apos++], bb[bpos++]);

                remain -= i - 1;
                errors++;
                simErrors++;

                if (sign == 1)
                {
                    if (!MummerFastaReader.IsAlpha(aa[apos++]))
                        nonAlphas++;
                    remain--;
                }
                else
                {
                    if (!MummerFastaReader.IsAlpha(bb[bpos++]))
                        nonAlphas++;
                }
            }

            for (long i = 0; i < remain; i++)
                Score(aa[apos++], bb[bpos++]);

            ap.Errors = errors;
            ap.SimErrors = simErrors;
            ap.NonAlphas = nonAlphas;
        }
    }

    private static int GetForwardTargetCluster(State st, List<Cluster> clusters, int currCp, ref long targetA, ref long targetB)
    {
        var cur = clusters[currCp];
        var last = cur.Matches[^1];
        long sA = last.SA + last.Len - 1;
        long sB = last.SB + last.Len - 1;
        int good = MummerPromerTables.GoodScore[st.Matrix], cont = MummerPromerTables.ContGapScore[st.Matrix];

        long dist = targetA - sA < targetB - sB ? targetA - sA : targetB - sB;

        int cp = clusters.Count;
        for (int cip = currCp + 1; cip < clusters.Count; cip++)
        {
            var c = clusters[cip];
            if (cur.FrameA == c.FrameA && cur.FrameB == c.FrameB)
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

                    if (greater <= st.BreakLen || lesser * good + (greater - lesser) * cont > 0)
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

    private static int GetReverseTargetAlignment(State st, List<Alignment> alignments, int currAp)
    {
        var cur = alignments[currAp];
        long sA = cur.SA, sB = cur.SB;
        int good = MummerPromerTables.GoodScore[st.Matrix], cont = MummerPromerTables.ContGapScore[st.Matrix];
        long dist = sA < sB ? sA : sB;

        int ap = alignments.Count;
        for (int aip = currAp - 1; aip >= 0; aip--)
        {
            var al = alignments[aip];
            if (cur.FrameA == al.FrameA && cur.FrameB == al.FrameB)
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

                    if (greater <= st.BreakLen || lesser * good + (greater - lesser) * cont > 0)
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
                if (al.FrameA == cur.FrameA && al.FrameB == cur.FrameB && al.EA >= eA && al.EB >= eB && al.SA <= sA && al.SB <= sB)
                    break;
            }

            if (aip >= 0)
                return true;
        }

        return false;
    }
}
