#pragma warning disable S3776, S107, S1541, S138, S1854, S1121, S127, S3358, S1066, S4144, S2325, S3871
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: dndfast7 (--fastapair / --fastswpair): fasta34 for every sequence against the following ones, the alignment blocks as pair library (hat3) and the distances of the scores (hat2)

    internal sealed partial class MafftGlobals
    {
        /// <summary>The globals of <c>dndfast7</c>: the program never reads the penalties, the offset or the Kimura parameter, so the zero of the C globals stays (the model line shows "-0.00"); <c>alg</c> is the NUL character. <paramref name="dna"/> is the decision of <c>getnumlen</c>.</summary>
        internal static MafftGlobals ForDndfast7(bool dna)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = 0, PpenaltyEx = 0, Poffset = 0, PpenaltyDist = 0, PpenaltyOp = 0, PpenaltyExGenaf = 0, KimuraR = 0, PamN = 0, Alg = '\0', Outgap = 1 };
            g.PenaltyShiftFactor = 1000.0;
            if (dna)
            {
                g.Dorp = 'd';
                g.Scoremtx = -1;
                g.PamN = NotSpecified;
            }
            else
            {
                g.Dorp = 'p';
                g.Nblosum = 62;
                g.Scoremtx = 1;
            }

            return g;
        }
    }

    internal sealed partial class MafftRun
    {
        /// <summary>What <c>dndfast7</c> leaves: the files <c>hat2</c> and <c>hat3</c> (the address of the next block that <c>%p</c> prints is a synthetic <c>0x...</c> value: the program prints a heap address), the standard error and the exit status; <see cref="Segfault"/>: the program stopped with a segmentation fault (the standard error is that of the run up to it).</summary>
        internal sealed record Dndfast7Result(string Hat2, string Hat3, string Err, int ExitCode, bool Segfault);

        /// <summary>The state that survives from one call of <c>ReadFasta34m10</c> to the next: the static buffers <c>qal</c> / <c>tal</c> of the program (the part behind a shorter string keeps the letters of earlier alignments).</summary>
        private sealed class Fasta34Reader
        {
            internal readonly char[] Qal = new char[1 << 20];
            internal readonly char[] Tal = new char[1 << 20];
        }

        /// <summary>A text as the C library reads it: <c>fgets</c>, <c>fgetc</c>, <c>ungetc</c>, <c>feof</c>.</summary>
        private sealed class CStream
        {
            private readonly string _t;
            private int _pos;
            private int _pushed = -2;

            internal CStream(string t) => _t = t;

            internal bool Eof { get; private set; }

            /// <summary><c>fgets( b, size, fp )</c>: null (and the buffer stays) when nothing could be read.</summary>
            internal string? Fgets(int size)
            {
                if (_pushed >= 0 && size > 1)
                {
                    // the character of ungetc comes first
                    var sb0 = new StringBuilder();
                    sb0.Append((char)_pushed);
                    _pushed = -2;
                    if (sb0[0] != '\n')
                    {
                        string? rest = Fgets(size - 1);
                        if (rest != null)
                            sb0.Append(rest);
                    }

                    return sb0.ToString();
                }

                if (_pos >= _t.Length)
                {
                    Eof = true;
                    return null;
                }

                int take = 0;
                while (_pos + take < _t.Length && take < size - 1)
                {
                    take++;
                    if (_t[_pos + take - 1] == '\n')
                        break;
                }

                string line = _t.Substring(_pos, take);
                _pos += take;
                if (_pos >= _t.Length && !line.EndsWith('\n'))
                    Eof = true;
                return line;
            }

            internal int Fgetc()
            {
                if (_pushed >= 0)
                {
                    int c = _pushed;
                    _pushed = -2;
                    return c;
                }

                if (_pos >= _t.Length)
                {
                    Eof = true;
                    return -1;
                }

                return _t[_pos++];
            }

            internal void Ungetc(int c)
            {
                _pushed = c;
                Eof = false;
            }
        }

        private static bool IsAlpha(int c) => c is >= 'a' and <= 'z' or >= 'A' and <= 'Z';

        /// <summary>The <c>strstr( b, ":" ) + 1</c> then <c>atoi</c> of the block readers (a line without a colon is a null pointer: the program stops with a segmentation fault).</summary>
        private static int AtoiAfterColon(string b)
        {
            int c = b.IndexOf(':');
            if (c < 0)
                throw new Dndfast7CrashException();
            return AtoiAt(b, c + 1);
        }

        private sealed class Dndfast7CrashException : Exception
        {
        }

        /// <summary><c>cutal</c> (<c>io.c</c>): the position (in <paramref name="al"/>) of the residue <paramref name="start"/> (or -1: a null pointer) and the string is cut after the position of residue <paramref name="end"/>.</summary>
        private static int Cutal(char[] al, int alDisplayStart, int start, int end)
        {
            int pos = alDisplayStart;
            int pt = 0;
            int val = -1;
            do
            {
                if (start == pos)
                    val = pt;
                if (end == pos)
                    break;
                if (al[pt] != '-')
                    pos++;
            }
            while (al[pt++] != 0);
            al[pt + 1] = '\0';
            return val;
        }

        /// <summary><c>putlocalhom2</c> (<c>divpairscore</c> 0) over the C strings <c>al1[p1..]</c> and <c>al2[p2..]</c>: the blocks of the ungapped runs, each with the mean score per aligned column of all blocks.</summary>
        private void PutLocalHom2C(char[] al1, int p1, char[] al2, int p2, MafftLocalHom head, int off1, int off2)
        {
            if (p1 < 0 || p2 < 0)
                throw new Dndfast7CrashException();
            var tmppt = head;
            int pos1 = off1, pos2 = off2, start1 = 0, start2 = 0, end1, end2;
            int isumscore = 0, sumoverlap = 0;
            bool first = true;
            int st = 0, iscore = 0;
            int pt1 = p1, pt2 = p2;
            while (al1[pt1] != 0)
            {
                char c1 = al1[pt1], c2 = al2[pt2];
                if (st == 1 && (c1 == '-' || c2 == '-'))
                {
                    end1 = pos1 - 1;
                    end2 = pos2 - 1;
                    if (!first)
                    {
                        tmppt.Next = new MafftLocalHom();
                        tmppt = tmppt.Next;
                    }

                    first = false;
                    tmppt.Start1 = start1;
                    tmppt.Start2 = start2;
                    tmppt.End1 = end1;
                    tmppt.End2 = end2;
                    tmppt.Korh = 'h';
                    isumscore += iscore;
                    sumoverlap += end2 - start2 + 1;
                    iscore = 0;
                    st = 0;
                }
                else if (c1 != '-' && c2 != '-')
                {
                    if (st == 0)
                    {
                        start1 = pos1;
                        start2 = pos2;
                        st = 1;
                    }

                    iscore += NDisOf(c1, c2);
                }

                if (c1 != '-')
                    pos1++;
                if (c2 != '-')
                    pos2++;
                pt1++;
                pt2++;
            }

            // *(pt1 - 1) != '-' && *(pt2 - 1) != '-'  (for an empty string the program reads the byte before the buffer; here the letter before the start)
            char l1 = pt1 > 0 ? al1[pt1 - 1] : '\0', l2 = pt2 > 0 ? al2[pt2 - 1] : '\0';
            if (l1 != '-' && l2 != '-')
            {
                if (!first)
                {
                    tmppt.Next = new MafftLocalHom();
                    tmppt = tmppt.Next;
                }

                first = false;
                end1 = pos1 - 1;
                end2 = pos2 - 1;
                tmppt.Start1 = start1;
                tmppt.Start2 = start2;
                tmppt.End1 = end1;
                tmppt.End2 = end2;
                tmppt.Korh = 'h';
                isumscore += iscore;
                sumoverlap += end2 - start2 + 1;
            }

            for (var t = head; t != null; t = t.Next)
            {
                t.Overlapaa = sumoverlap;
                t.Opt = sumoverlap > 0 ? (double)isumscore * 5.8 / (600 * sumoverlap) : -1.0;
            }
        }

        private int NDisOf(char c1, char c2)
        {
            int a = c1 < G.AminoN.Length ? G.AminoN[c1] : -1;
            int b = c2 < G.AminoN.Length ? G.AminoN[c2] : -1;
            if (a < 0)
                throw new Dndfast7CrashException(); // n_dis[-1]: the program reads the header of the allocation as a pointer
            if (b < 0)
                return 0; // n_dis[a][-1]: the four bytes before a row are the upper half of the size field of its allocation (0 below 4 GB)
            return G.NDis[a][b];
        }

        /// <summary><c>ReadFasta34m10</c> / <c>ReadFasta34m10_nuc</c> (<c>io.c</c>): the optimised scores of the list of the best scores into <paramref name="dis"/> and the alignment blocks into the chains of <paramref name="localhomlist"/> (indexed by library sequence); returns the number of blocks (or -1 without a list).</summary>
        private int ReadFasta34m10(Fasta34Reader st, string text, double[] dis, MafftLocalHom[] localhomlist, bool nuc)
        {
            var fp = new CStream(text);
            const int B = 256;
            string b = "";
            var junban = new List<int>();
            int count = 0;
            int optValue = 0;
            bool optDefined = false;
            while (!fp.Eof)
            {
                string? line = fp.Fgets(B - 1);
                if (line != null)
                    b = line;
                if (b.StartsWith("+==========+", StringComparison.Ordinal))
                {
                    while (junban.Count <= count)
                        junban.Add(0);
                    junban[count] = AtoiAt(b, 12);
                    if (nuc && b.Contains('r'))
                        continue;
                    int pp = b.IndexOf(nuc ? ']' : ')');
                    if (pp < 0)
                        throw new Dndfast7CrashException();
                    var sc = new CScanner(b[(pp + 1)..]);
                    if (sc.Int(out int scanned))
                    {
                        optValue = scanned;
                        optDefined = true;
                    }
                    else if (!optDefined)
                    {
                        // sscanf( pt, "%d %lf %lf", &opt, &bits, &z ) fails on a list line of the other kind (fasta34 recognised a DNA query in a protein run: "(  40) [f]  200 ...") and the variable opt of the function was never assigned:
                        // the program stores what its stack holds (32767 = the upper half of a stack address under the emulator of the oracle run; another value elsewhere)
                        throw UndefinedBehaviour("the list of fasta34 has the format of the other kind of sequence (fasta34 recognised a DNA query in a protein run: lines with '[f]' / '[r]'): the program reads the score into an unassigned variable (the upper half of a stack address)");
                    }

                    dis[junban[count]] = optValue;
                    count++;
                }
                else if (b.StartsWith(">>+==========+", StringComparison.Ordinal))
                {
                    break;
                }
            }

            if (count == 0)
                return -1;
            count = 0;
            while (true)
            {
                if (!b.StartsWith(">>+==========+", StringComparison.Ordinal))
                {
                    string? line = fp.Fgets(B - 1);
                    if (line != null)
                        b = line;
                    if (fp.Eof)
                        break;
                    continue;
                }

                while (junban.Count <= count)
                    junban.Add(0);
                junban[count++] = AtoiAt(b, 14);
                string? l;
                while ((l = fp.Fgets(B - 1)) != null)
                {
                    b = l;
                    if (b.StartsWith("; fa_opt:", StringComparison.Ordinal) || b.StartsWith("; sw_s-w opt:", StringComparison.Ordinal))
                        break;
                }

                int opt = AtoiAfterColon(b);
                int overlapaa = ReadField(fp, ref b, "_overlap:", 9);
                int qstart = ReadField(fp, ref b, "_start:", 7) - 1;
                int qend = ReadField(fp, ref b, "_stop:", 6) - 1;
                int qalDisplayStart = ReadField(fp, ref b, "_display_start:", 15) - 1;
                ReadResidues(fp, st.Qal);
                int tstart = ReadField(fp, ref b, "_start:", 7) - 1;
                int tend = ReadField(fp, ref b, "_stop:", 6) - 1;
                int talDisplayStart = ReadField(fp, ref b, "_display_start:", 15) - 1;
                ReadResidues(fp, st.Tal);
                _ = overlapaa;
                _ = opt;
                int qal2 = Cutal(st.Qal, qalDisplayStart, qstart, qend);
                int tal2 = Cutal(st.Tal, talDisplayStart, tstart, tend);
                PutLocalHom2C(st.Qal, qal2, st.Tal, tal2, localhomlist[junban[count - 1]], qstart, tstart);
            }

            return count;
        }

        /// <summary><c>while( fgets( b, B-1, fp ) ) if( !strncmp( key, b+4, n ) ) break;</c> then <c>atoi( strstr( b, ":" ) + 1 )</c>.</summary>
        private static int ReadField(CStream fp, ref string b, string key, int n)
        {
            string? l;
            while ((l = fp.Fgets(255)) != null)
            {
                b = l;
                if (b.Length > 4 && string.CompareOrdinal(b, 4, key, 0, n) == 0)
                    break;
            }

            return AtoiAfterColon(b);
        }

        /// <summary>The loop <c>while( (c = fgetc( fp )) ) { if( c == '&gt;' ) { ungetc; break; } if( isalpha( c ) || c == '-' ) *pt++ = c; } *pt = 0;</c>.</summary>
        private static void ReadResidues(CStream fp, char[] buf)
        {
            int pt = 0;
            int c;
            while ((c = fp.Fgetc()) != 0)
            {
                if (c == -1)
                    throw UndefinedBehaviour("the alignment of fasta34 ends the file inside a sequence: the program reads the end of the file for ever");
                if (c == '>')
                {
                    fp.Ungetc(c);
                    break;
                }

                if (IsAlpha(c) || c == '-')
                    buf[pt++] = (char)c;
            }

            buf[pt] = '\0';
        }

        /// <summary>
        /// <c>dndfast7 [-A] &lt; infile</c> (<c>dndfast7.c</c>): for every sequence <c>i</c> the call <c>fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000 query library ktup</c> (library: the sequences <c>i ..</c>; ktup 1 for protein, 6 for DNA), the scores of
        /// the list of the best scores (<c>hat2</c>: 2 (1 - score / self score), the smaller of the two directions) and the alignment blocks (<c>hat3</c>).
        /// </summary>
        internal static Dndfast7Result Dndfast7(IReadOnlyList<string> names, IReadOnlyList<string> seqs, bool sw)
        {
            int njob = seqs.Count;
            bool dna = MafftEngine.GuessNucleotide(seqs, new MafftOptions());
            var g = MafftGlobals.ForDndfast7(dna);
            var run = new MafftRun(g);
            return run.Dndfast7Core(names, seqs, sw, njob);
        }

        private Dndfast7Result Dndfast7Core(IReadOnlyList<string> names, IReadOnlyList<string> seqs, bool sw, int njob)
        {
            var fastaFilter = seqs.Select(s => s.Replace("-", "", StringComparison.Ordinal)).ToArray();
            G.Constants(fastaFilter);
            string constantsLog = G.Log.ToString();
            var err = new StringBuilder(constantsLog);
            G.Log.Clear();
            bool dna = G.Dorp == 'd';
            int ktuple = dna ? 6 : 1;
            string[] name1 = new string[njob];
            for (int j = 0; j < njob; j++)
                name1[j] = "+==========+" + j.ToString(CultureInfo.InvariantCulture) + new string(' ', 22);
            var mtx = new double[njob][];
            var mtx2 = new double[njob][];
            for (int i = 0; i < njob; i++)
            {
                mtx[i] = new double[njob];
                mtx2[i] = new double[njob];
            }

            var localhomtable = new MafftLocalHom[njob][];
            for (int i = 0; i < njob; i++)
            {
                localhomtable[i] = new MafftLocalHom[njob];
                for (int j = 0; j < njob; j++)
                    localhomtable[i][j] = new MafftLocalHom();
            }

            var reader = new Fasta34Reader();
            for (int i = 0; i < njob; i++)
            {
                var libTitles = new List<string>();
                var library = new List<string>();
                for (int j = i; j < njob; j++)
                {
                    libTitles.Add(name1[j]);
                    library.Add(fastaFilter[j]);
                }

                var fa = MafftFasta34.Run(name1[i], fastaFilter[i], libTitles, library, dna, sw, ktuple);
                err.Append(fa.Err);
                if (fa.ExitCode != 0)
                {
                    err.Append("error in fasta\n");
                    return new Dndfast7Result("", "", err.ToString(), 1, false);
                }

                int res;
                try
                {
                    res = ReadFasta34m10(reader, fa.Text, mtx[i], localhomtable[i], dna);
                }
                catch (Dndfast7CrashException)
                {
                    return new Dndfast7Result("", "", err.ToString(), 139, true);
                }

                if (res < njob - i)
                {
                    err.Append("count (fasta34 -z 3) = ").Append(res.ToString(CultureInfo.InvariantCulture)).Append('\n');
                    return new Dndfast7Result("", "", err.ToString(), 1, false);
                }

                err.Append("query : ").Append(MafftGlobals.F(i + 1, 4, 0)).Append(" / ").Append(MafftGlobals.F(njob, 5, 0)).Append('\r');
            }

            for (int i = 0; i < njob; i++)
            {
                double max = mtx[i][i];
                if (max == 0.0)
                {
                    for (int j = 0; j < njob; j++)
                        mtx2[i][j] = 2.0;
                }
                else
                {
                    for (int j = 0; j < njob; j++)
                        mtx2[i][j] = (max - mtx[Math.Min(i, j)][Math.Max(i, j)]) / max * 2.0;
                }
            }

            for (int i = 0; i < njob - 1; i++)
            {
                for (int j = i + 1; j < njob; j++)
                    mtx2[i][j] = Math.Min(mtx2[i][j], mtx2[j][i]);
            }

            var half = new double[njob][];
            for (int i = 0; i < njob; i++)
            {
                half[i] = new double[njob - i];
                for (int j = i; j < njob; j++)
                    half[i][j - i] = mtx2[i][j];
            }

            string hat2 = WriteHat2(half, names);
            err.Append("##### writing hat3\n");
            var hat3 = new StringBuilder();
            long address = 0x55555bf86480;
            for (int i = 0; i < njob - 1; i++)
            {
                for (int j = i + 1; j < njob; j++)
                {
                    for (var t = localhomtable[i][j]; t != null; t = t.Next)
                    {
                        if (t.Opt == -1.0)
                            continue;
                        hat3.Append(i.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(j.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(t.Overlapaa.ToString(CultureInfo.InvariantCulture)).Append(' ')
                            .Append(MafftGlobals.F(t.Opt, 6, 3)).Append(' ').Append(t.Start1.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(t.End1.ToString(CultureInfo.InvariantCulture)).Append(' ')
                            .Append(t.Start2.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(t.End2.ToString(CultureInfo.InvariantCulture)).Append(' ');
                        if (t.Next == null)
                        {
                            hat3.Append("(nil)");
                        }
                        else
                        {
                            hat3.Append("0x").Append(address.ToString("x", CultureInfo.InvariantCulture));
                            address += 0x50;
                        }

                        hat3.Append('\n');
                    }
                }
            }

            err.Append("dndfast7 (").Append(G.KindName).Append(") Version 7.505\nalg=").Append(G.Alg).Append(", model=").Append(G.ModelName).Append(", amax=").Append(G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture))
                .Append('\n').Append(G.Nthread.ToString(CultureInfo.InvariantCulture)).Append(" thread(s)\n\n");
            return new Dndfast7Result(hat2, hat3.ToString(), err.ToString(), 0, false);
        }
    }

    #endregion
}
