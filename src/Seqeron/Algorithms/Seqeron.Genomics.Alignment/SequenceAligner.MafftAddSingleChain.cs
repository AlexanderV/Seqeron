#pragma warning disable S127, S3358, S1066, S3776, S107, S1854, IDE0059, S2325, S138, S1541
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: the programs that make the distance files and the pair library of addsingle (dndpre, pairlocalalign -Y / -L, disttbfast -T -y)

    internal sealed partial class MafftRun
    {
        /// <summary>The files and the standard error that a program of the chain leaves.</summary>
        internal sealed record PairlocalalignResult(string Hat2, string Hat3);

        /// <summary>
        /// <c>pairlocalalign</c> (the program, <c>pairlocalalignmain.c</c>, with <c>-Y</c> or <c>-L</c> and <c>-I nadd</c>): the local alignment of the pairs of sequences (<c>-Y</c>: every old sequence with every new one; <c>-L</c>: every pair, the pairs of two
        /// new sequences get the score 0), the distances to the file <c>hat2</c> and the ungapped blocks of the alignments to the file <c>hat3</c>. The messages go to <see cref="MafftGlobals.Log"/>.
        /// </summary>
        internal PairlocalalignResult PairlocalalignProgram(IReadOnlyList<string> seq0, IReadOnlyList<string> names)
        {
            int njob = seq0.Count;
            int nadd = G.Nadd;
            if (njob < 2)
                throw new MafftProgramFailureException("At least 2 sequences should be input!\nOnly " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
            if (G.Alg is not ('L' or 'Y'))
                throw new NotSupportedException("pairlocalalign of this algorithm is not available in the chain of addsingle.");
            bool part = G.Alg == 'Y';
            if (part && nadd == 0 && G.Nthread > 0)
                throw UndefinedBehaviour("pairlocalalign -Y with threads and no sequence added (--hybridpair --thread n without --add): the threads take the pair (njob - nadd - 1, njob) of the job list that does not exist (valgrind: 'Invalid read of size 1', address 0x0 in a thread; segmentation fault whose message names a process id)");
            var seq = seq0.ToArray();
            _njob = njob;
            _nlenmax = seq.Max(x => x.Length);
            int alloclen = _nlenmax * 2;
            G.Constants(seq);
            if (CheckAlphabet(seq) is { } unknownLetter)
                throw unknownLetter;
            var bseq = new string[njob];
            var dseq = new string[njob];
            var thereisx = new int[njob];
            for (int i = 0; i < njob; i++)
            {
                bseq[i] = GapPick0(seq[i]);
                if (G.Dorp == 'p' && G.Scoremtx == 1 && G.Nblosum > 0)
                {
                    int count = 0;
                    var sb = new StringBuilder(bseq[i].Length);
                    foreach (char c in bseq[i])
                    {
                        if (c is 'X' or 'x')
                            count++;
                        else
                            sb.Append(c);
                    }

                    dseq[i] = sb.ToString();
                    thereisx[i] = count;
                }
                else
                {
                    dseq[i] = bseq[i];
                }
            }

            int norg = njob - nadd;
            var localhom = new MafftLocalHom[part ? norg : njob - 1][];
            for (int i = 0; i < localhom.Length; i++)
            {
                int cols = part ? nadd : njob - i;
                localhom[i] = new MafftLocalHom[cols];
                for (int j = 0; j < cols; j++)
                    localhom[i][j] = new MafftLocalHom();
            }

            var distance = new double[part ? norg : njob][];
            for (int i = 0; i < distance.Length; i++)
                distance[i] = new double[part ? nadd : njob - i];
            var selfscore = new double[njob];
            Report("All-to-all alignment.\n");
            for (int i = 0; i < njob; i++)
            {
                double s = 0.0;
                foreach (char c in bseq[i])
                    s += G.AminoDis[c][c];
                selfscore[i] = s;
            }

            int ilim = part ? norg : njob - 1;
            for (int i = 0; i < ilim; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                if (G.Nthread > 0)
                {
                    // the threads: a progress line at the first column of the first row of every tenth row (the first thread does all the work)
                    if (!part || i == norg - 1)
                        PairlocalalignThreadProgress(i, njob, part, norg);
                }
                else
                {
                    Report(MafftGlobals.D5(i) + " / " + norg.ToString(CultureInfo.InvariantCulture) + "\r");
                }

                int jst = part ? norg : i + 1;
                for (int j = jst; j < njob; j++)
                {
                    if (bseq[i].Length == 0 || bseq[j].Length == 0)
                    {
                        if (part)
                            distance[i][j - norg] = 3.0;
                        else
                            distance[i][j - i] = 3.0;
                        continue;
                    }

                    var m1 = new[] { bseq[i] };
                    var m2 = new[] { bseq[j] };
                    int tx = thereisx[i] + thereisx[j];
                    double pscore;
                    int off1 = 0, off2 = 0;
                    bool newNew = !part && nadd > 0 && j >= norg && i >= norg;
                    if (newNew)
                    {
                        pscore = 0.0; // the pair of two new sequences: not aligned (the program then takes the two unaligned sequences as an alignment and reads past the end of the shorter one; addsingle never looks at these blocks)
                    }
                    else if (G.UseNaiveScore)
                    {
                        L__align11(G.NDisConsweightMulti, 0.0, m1, m2, alloclen, true, out off1, out off2);
                        pscore = NaivePairScore11(m1[0], m2[0], 0);
                    }
                    else if (G.StoreLocalhom != 0)
                    {
                        pscore = L__align11(G.NDisConsweightMulti, 0.0, m1, m2, alloclen, true, out off1, out off2);
                        if (tx != 0)
                            pscore = L__align11(G.NDisConsweightMulti, 0.0, new[] { dseq[i] }, new[] { dseq[j] }, alloclen, false, out _, out _);
                    }
                    else
                    {
                        pscore = L__align11(G.NDisConsweightMulti, 0.0, new[] { dseq[i] }, new[] { dseq[j] }, alloclen, false, out _, out _);
                    }

                    if (m1[0].Length != 0 && m2[0].Length != 0)
                    {
                        if (G.StoreLocalhom != 0 && !newNew)
                        {
                            if (part)
                                PutLocalHom2(m1[0], m2[0], localhom[i][j - norg], off1, off2);
                            else
                                PutLocalHom2(m1[0], m2[0], localhom[i][j - i], off1, off2);
                        }

                        pscore = Score2Dist(pscore, selfscore[i], selfscore[j]);
                    }
                    else
                    {
                        pscore = 2.0;
                    }

                    if (part)
                        distance[i][j - norg] = pscore;
                    else
                        distance[i][j - i] = pscore;
                }
            }

            string hat2 = part ? WriteHat2Part(njob, nadd, names, distance) : WriteHat2(distance, names);
            var hat3 = new StringBuilder();
            if (G.StoreLocalhom != 0)
            {
                Report("\n\n##### writing hat3\n");
                for (int i = 0; i < ilim; i++)
                {
                    int jst = part ? norg : i;
                    int jj = 0;
                    for (int j = jst; j < njob; j++, jj++)
                    {
                        for (var t = localhom[i][jj]; t != null; t = t.Next)
                        {
                            if (t.Opt == -1.0)
                                continue;
                            if (i < j)
                                hat3.Append(i).Append(' ').Append(j).Append(' ').Append(t.Overlapaa).Append(' ').Append(MafftGlobals.F(t.Opt, 7, 5)).Append(' ')
                                    .Append(t.Start1).Append(' ').Append(t.End1).Append(' ').Append(t.Start2).Append(' ').Append(t.End2).Append(" h\n");
                        }
                    }
                }
            }

            Report("pairlocalalign (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new PairlocalalignResult(hat2, hat3.ToString());
        }

        /// <summary>The line of the thread that takes the pair (i, i + 1) of every tenth row (<c>"% 5d / %d (by thread %3d) \r"</c>).</summary>
        private void PairlocalalignThreadProgress(int i, int njob, bool part, int norg)
        {
            int first = part ? norg : i + 1;
            if (first == i + 1 && i % 10 == 0)
                Report(MafftGlobals.D5(i) + " / " + norg.ToString(CultureInfo.InvariantCulture) + " (by thread " + 0.ToString(CultureInfo.InvariantCulture).PadLeft(3) + ") \r");
        }

        /// <summary><c>WriteHat2_part_pointer</c>: the distances of the old sequences to the new ones (the file <c>hat2</c> of <c>pairlocalalign -Y</c>).</summary>
        private static string WriteHat2Part(int njob, int nadd, IReadOnlyList<string> names, double[][] mtx)
        {
            var sb = new StringBuilder();
            int norg = njob - nadd;
            sb.Append(MafftGlobals.F(1, 5, 0)).Append('\n').Append(njob.ToString(CultureInfo.InvariantCulture).PadLeft(5)).Append('\n');
            sb.Append(' ').Append(Hat2Number(0.0)).Append('\n');
            for (int i = 0; i < njob; i++)
                sb.Append((i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4)).Append(". ").Append(names[i]).Append('\n');
            for (int i = 0; i < norg; i++)
            {
                for (int j = 0; j < nadd; j++)
                {
                    sb.Append(Hat2Number(mtx[i][j]));
                    if ((j + 1) % 12 == 0 || j == nadd - 1)
                        sb.Append('\n');
                }
            }

            return sb.ToString();
        }

        /// <summary><c>readhat2_doublehalf_part_pointer</c>: the distances of the old sequences to the new ones (fields of 6 characters).</summary>
        private static void ReadHat2PartText(string text, int njob, int nadd, double[][] nmtx)
        {
            int pos = 0;
            string Line()
            {
                int e = text.IndexOf('\n', pos);
                string l = e < 0 ? text[pos..] : text[pos..e];
                pos = e < 0 ? text.Length : e + 1;
                return l;
            }

            Line();
            int nseq0 = MafftGlobals.CAtoi(Line());
            if (nseq0 != njob)
                throw new MafftProgramFailureException(njob.ToString(CultureInfo.InvariantCulture) + " != " + nseq0.ToString(CultureInfo.InvariantCulture) + "\nhat2 is wrong.\n");
            Line();
            for (int i = 0; i < njob; i++)
                Line();
            int norg = njob - nadd;
            for (int i = 0; i < norg; i++)
            {
                for (int j = 0; j < nadd; j++)
                {
                    if (pos < text.Length && text[pos] == '\n')
                        pos++;
                    int take = Math.Min(6, Math.Max(0, text.Length - pos));
                    nmtx[i][j] = MafftGlobals.CAtof(text.Substring(pos, take));
                    pos += take;
                }
            }
        }

        // -------------------------------------------------------------------- hat3 as addsingle reads it

        /// <summary>The table of blocks that <c>addsingle</c> reads from the file <c>hat3</c>: <c>readlocalhomtable_one</c> (<c>-d</c>: the old sequences against the new ones, one column per new sequence) or <c>readlocalhomtable</c> (every pair, both directions).</summary>
        private static MafftLocalHom[][] ReadHat3ForAddsingle(string hat3Text, int norg, int nadd, bool one, bool[] kozo)
        {
            int njob = norg + nadd;
            int rows = one ? norg : njob;
            int cols = one ? nadd : njob;
            var table = new MafftLocalHom[rows][];
            var count = new int[rows][];
            var last = new MafftLocalHom[rows][];
            for (int i = 0; i < rows; i++)
            {
                table[i] = new MafftLocalHom[cols];
                count[i] = new int[cols];
                last[i] = new MafftLocalHom[cols];
                for (int j = 0; j < cols; j++)
                    table[i][j] = new MafftLocalHom();
            }

            foreach (string buff in SplitFgetsLines(hat3Text, 256))
            {
                var sc = new CScanner(buff);
                if (!sc.Int(out int i) || !sc.Int(out int j) || !sc.Int(out int overlapaa) || !sc.Double(out double opt) || !sc.Int(out int start1) || !sc.Int(out int end1) || !sc.Int(out int start2) || !sc.Int(out int end2) || !sc.Word(out string infor))
                    throw new NotSupportedException("A line of the file hat3 has fewer than nine fields: the program reads the variables of the line before.");
                if (one)
                {
                    if (infor.Length > 0 && infor[0] == 'k')
                        throw new MafftProgramFailureException("Not supported!\n");
                    j -= norg;
                    Append(table[i], count[i], last[i], j, start1, end1, start2, end2, opt, overlapaa, infor[0]);
                }
                else
                {
                    if (infor.Length > 0 && infor[0] == 'k')
                        kozo[i] = kozo[j] = true;
                    if (j <= i)
                        throw new MafftProgramFailureException("Check hat3.  The first sequence must be younger than the second one.\n");
                    Append(table[i], count[i], last[i], j, start1, end1, start2, end2, opt, overlapaa, infor[0]);
                    Append(table[j], count[j], last[j], i, start2, end2, start1, end1, opt, overlapaa, infor[0]);
                }
            }

            return table;

            static void Append(MafftLocalHom[] row, int[] cnt, MafftLocalHom[] lst, int j, int s1, int e1, int s2, int e2, double opt, int overlapaa, char korh)
            {
                MafftLocalHom node;
                if (cnt[j]++ > 0)
                {
                    node = new MafftLocalHom { Importance = 0.0 };
                    lst[j].Next = node;
                }
                else
                {
                    node = row[j];
                }

                lst[j] = node;
                node.Start1 = s1;
                node.Start2 = s2;
                node.End1 = e1;
                node.End2 = e2;
                node.Opt = (opt + 0.00) / 5.8 * 600;
                node.Overlapaa = overlapaa;
                node.Korh = korh;
            }
        }
    }

    #endregion
}
