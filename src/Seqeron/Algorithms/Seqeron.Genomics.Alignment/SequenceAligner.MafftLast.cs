using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: pairlocalalign -R (--lastpair): the alignments of LAST

    internal sealed partial class MafftRun
    {
        private sealed class LastAlignmentRecord
        {
            internal int Score;
            internal readonly List<(int Start1, int End1, int Start2, int End2)[]> Alignments = [];
        }

        /// <summary>
        /// The refusal for the heap overflow of <c>pairlocalalign -R</c> with three or more related sequences: <c>putlocalhom_last( ..., localhomtable[i]+j, lastresx[i]+j, 'h' )</c> indexes the row i of the table with the
        /// column j, but the row has only <c>njob - i</c> entries (<c>jst--</c> for every algorithm but <c>-Y</c> / <c>-r</c>): a pair (i &gt;= 1, j &gt;= njob - i) that has an alignment writes behind the block.
        /// </summary>
        private static NotSupportedException LastHeapOverflow(int i, int j, int njob) => UndefinedBehaviour(
            $"pairlocalalign -R writes the alignments of the pair ({i}, {j}) of {njob} sequences behind the row {i} of its table of local homologies (the row holds {njob - i} entries; putlocalhom_last indexes it with j, not j - i). "
            + "valgrind on the installed program: 'Invalid write of size 4 ... 16 bytes after a block of size ...'; the program then stops with a segmentation fault or 'double free or corruption' (a message with a process id; 100 of 100 runs of the oracle)");

        /// <summary>The text of the file <c>_scoringmatrixforlast</c>: the matrix of the sequence type in the letters of the alphabet.</summary>
        private string LastScoringMatrixText()
        {
            int size = G.Dorp == 'd' ? 4 : 20;
            var sb = new StringBuilder("      ");
            for (int j = 0; j < size; j++)
                sb.Append(' ').Append(G.Amino[j]).Append(' ');
            sb.Append('\n');
            for (int i = 0; i < size; i++)
            {
                sb.Append(G.Amino[i]).Append(' ');
                for (int j = 0; j < size; j++)
                    sb.Append(' ').Append(G.NDis[i][j].ToString(CultureInfo.InvariantCulture)).Append(' ');
                sb.Append('\n');
            }

            return sb.ToString();
        }

        /// <summary><c>block2reg</c>: the blocks of a line of the TAB format (<c>len,gap1:gap2,len ...</c>) as regions of the two sequences.</summary>
        private static (int Start1, int End1, int Start2, int End2)[] Block2Reg(string blocks, int start1, int start2)
        {
            var regs = new List<(int, int, int, int)>();
            int pos1 = start1, pos2 = start2;
            int p = 0;
            while (true)
            {
                string rest = blocks[p..];
                int tpt = rest.IndexOf(':');
                int npt = rest.IndexOf(',');
                if (tpt < 0 || (npt >= 0 && tpt > npt) || npt < 0)
                {
                    int len = MafftGlobals.CAtoi(rest);
                    int s1 = pos1, s2 = pos2;
                    pos1 += len - 1;
                    pos2 += len - 1;
                    regs.Add((s1, pos1, s2, pos2));
                }
                else
                {
                    string[] parts = rest[..npt].Split(':');
                    int glen1 = MafftGlobals.CAtoi(parts[0]);
                    int glen2 = MafftGlobals.CAtoi(parts[1]);
                    pos1 += glen1 + 1;
                    pos2 += glen2 + 1;
                }

                if (npt < 0)
                    break;
                p += npt + 1;
            }

            return regs.ToArray();
        }

        /// <summary><c>readlastresx_singleq</c>: the lines of the output of <c>lastal</c> for the query <paramref name="nameq"/> are added to the table of the pairs.</summary>
        private void ReadLastResultSingleQuery(string text, int nameq, LastAlignmentRecord[][] lastresx, bool subopt, int penalty)
        {
            foreach (string raw in text.Split('\n'))
            {
                if (raw.Length == 0 || raw[0] == '#')
                    continue;
                string[] f = raw.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
                int score = MafftGlobals.CAtoi(f[0]);
                int name1 = MafftGlobals.CAtoi(f[1]);
                int start1 = MafftGlobals.CAtoi(f[2]);
                int name2 = MafftGlobals.CAtoi(f[6]);
                int start2 = MafftGlobals.CAtoi(f[7]);
                if (G.Alg == 'R' && name2 <= name1)
                    continue;
                if (name2 != nameq)
                {
                    throw new MafftProgramFailureException("BUG!!!\n", 1);
                }

                var rec = lastresx[name1][name2];
                bool includeIntoScore = rec.Alignments.Count == 0;
                if (!includeIntoScore && !subopt)
                    continue;
                int tab = raw.LastIndexOf('\t');
                rec.Alignments.Add(Block2Reg(raw[(tab + 1)..], start1, start2));
                if (includeIntoScore)
                {
                    if (rec.Score != 0)
                        score += penalty;
                    rec.Score += score;
                }
            }
        }

        /// <summary><c>putlocalhom_last</c>: the regions of the alignments of a pair as blocks of the local homology (all the blocks of an alignment have the sum of the scores and of the lengths of its blocks).</summary>
        private void PutLocalHomLast(string s1, string s2, MafftLocalHom first, LastAlignmentRecord rec)
        {
            if (rec.Alignments.Count == 0)
                return;
            var tmppt = first;
            int nlocalhom = 0;
            foreach (var alignment in rec.Alignments)
            {
                int isumscore = 0, sumoverlap = 0;
                MafftLocalHom? first0 = null;
                for (int r = 0; r < alignment.Length; r++)
                {
                    if (nlocalhom++ > 0)
                    {
                        tmppt.Next = new MafftLocalHom { Start1 = 0, End1 = 0, Start2 = 0, End2 = 0, Opt = 0.0, Overlapaa = 0, Korh = '\0' };
                        tmppt = tmppt.Next;
                    }

                    var (a1, e1, a2, e2) = alignment[r];
                    tmppt.Start1 = a1;
                    tmppt.Start2 = a2;
                    tmppt.End1 = e1;
                    tmppt.End2 = e2;
                    tmppt.Korh = 'h';
                    if (r == 0)
                        first0 = tmppt;
                    int len = e1 - a1 + 1;
                    int iscore = 0;
                    for (int k = 0; k < len; k++)
                        iscore += G.NDis[G.AminoN[s1[a1 + k]]][G.AminoN[s2[a2 + k]]];
                    isumscore += iscore;
                    sumoverlap += e2 - a2 + 1;
                }

                for (var t = first0; t != null; t = t.Next)
                {
                    t.Overlapaa = sumoverlap;
                    t.Opt = (double)isumscore * 5.8 / (600 * sumoverlap);
                }
            }
        }

        /// <summary>
        /// <c>pairlocalalign -R</c> (<c>--lastpair</c>): every sequence is the query of one call of <c>lastal</c> against the database of the sequences before it (<c>calllast_fast</c>, <c>lastcallthread</c>), the first
        /// alignment of a pair gives its score (the distance of <c>hat2</c>) and the blocks of its alignment <c>hat3</c> (with the table index of the program: see <see cref="LastHeapOverflow"/>). The programs of LAST are run by <paramref name="last"/> in a working
        /// directory (the files <c>_dbd</c>, <c>_db&lt;k&gt;d</c>, <c>_q&lt;k&gt;</c>, <c>_scoringmatrixforlast</c>).
        /// </summary>
        internal PairlocalalignResult PairlocalalignLastProgram(IReadOnlyList<string> seq0, IReadOnlyList<string> names, ILastPrograms last, string lastDirectory, string whereIsPairAlign, int laste, int lastm, bool subopt)
        {
            int njob = seq0.Count;
            Report("laste = " + laste.ToString(CultureInfo.InvariantCulture) + "\nlastm = " + lastm.ToString(CultureInfo.InvariantCulture) + "\nwhereispairalign = " + whereIsPairAlign + "\n");
            if (njob < 2)
                throw new MafftProgramFailureException("At least 2 sequences should be input!\nOnly " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
            if (G.Alg != 'R')
                throw new NotSupportedException("pairlocalalign -r (--lastmultipair) is not available: the script stops in dndpre before.");
            if (G.Dorp == 'p')
            {
                throw new MafftProgramFailureException("Not yet supported\n", 1);
            }

            var seq = seq0.ToArray();
            _njob = njob;
            _nlenmax = seq.Max(x => x.Length);
            G.Constants(seq);
            if (CheckAlphabet(seq) is { } unknownLetter)
                throw unknownLetter;
            var bseq = new string[njob];
            for (int i = 0; i < njob; i++)
                bseq[i] = GapPick0(seq[i]);
            var lastresx = new LastAlignmentRecord[njob][];
            for (int i = 0; i < njob; i++)
            {
                lastresx[i] = new LastAlignmentRecord[njob];
                for (int j = 0; j < njob; j++)
                    lastresx[i][j] = new LastAlignmentRecord();
            }

            // the table of local homologies: the row i has njob - i entries (jst--)
            var localhom = new MafftLocalHom[njob - 1][];
            for (int i = 0; i < localhom.Length; i++)
            {
                localhom[i] = new MafftLocalHom[njob - i];
                for (int j = 0; j < localhom[i].Length; j++)
                    localhom[i][j] = new MafftLocalHom();
            }

            var distance = new double[njob][];
            for (int i = 0; i < njob; i++)
                distance[i] = new double[njob - i];
            Report("All-to-all alignment.\n");
            Report("Calling last (http://last.cbrc.jp/)\n");
            RunLastQueries(bseq, G.Nadd, last, lastDirectory, laste, lastm, subopt, lastresx);
            Report("done.\n");
            var selfscore = new double[njob];
            for (int i = 0; i < njob; i++)
            {
                double s = 0.0;
                foreach (char c in bseq[i])
                    s += G.AminoDis[c][c];
                selfscore[i] = s;
            }

            for (int i = 0; i < njob - 1; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                if (G.Nthread > 0)
                    PairlocalalignThreadProgress(i, njob, false, njob);
                else
                    Report(MafftGlobals.D5(i) + " / " + njob.ToString(CultureInfo.InvariantCulture) + "\r");
                for (int j = i + 1; j < njob; j++)
                {
                    if (bseq[i].Length == 0 || bseq[j].Length == 0)
                    {
                        distance[i][j - i] = 3.0;
                        continue;
                    }

                    int norg = njob - G.Nadd;
                    double pscore = G.Nadd > 0 && norg <= j && norg <= i ? 0.0 : lastresx[i][j].Score;
                    if (lastresx[i][j].Alignments.Count > 0 && j >= njob - i)
                        throw LastHeapOverflow(i, j, njob);
                    PutLocalHomLast(bseq[i], bseq[j], localhom[i][j < localhom[i].Length ? j : 0], lastresx[i][j]);
                    distance[i][j - i] = Score2Dist(pscore, selfscore[i], selfscore[j]);
                }
            }

            string hat2 = WriteHat2(distance, names);
            var hat3 = new StringBuilder();
            Report("\n\n##### writing hat3\n");
            for (int i = 0; i < njob - 1; i++)
            {
                int jj = 0;
                for (int j = i; j < njob; j++, jj++)
                {
                    for (var t = localhom[i][jj]; t != null; t = t.Next)
                    {
                        if (t.Opt == -1.0)
                            continue;
                        hat3.Append(i).Append(' ').Append(j).Append(' ').Append(t.Overlapaa).Append(' ').Append(MafftGlobals.F(t.Opt, 7, 5)).Append(' ')
                            .Append(t.Start1).Append(' ').Append(t.End1).Append(' ').Append(t.Start2).Append(' ').Append(t.End2).Append(" h\n");
                    }
                }
            }

            Report("pairlocalalign (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new PairlocalalignResult(hat2, hat3.ToString());
        }

        /// <summary><c>calllast_fast</c> / <c>lastcallthread</c>: the database of all the sequences, then for every query k the database of the k sequences before it, the query file and the call of <c>lastal</c>.</summary>
        private void RunLastQueries(string[] dseq, int nadd, ILastPrograms last, string lastDirectory, int laste, int lastm, bool subopt, LastAlignmentRecord[][] lastresx)
        {
            int njob = dseq.Length;
            var dir = new Dictionary<string, string>
            {
                ["_scoringmatrixforlast"] = LastScoringMatrixText(),
            };
            var db = new StringBuilder();
            for (int i = 0; i < njob - nadd; i++)
                db.Append('>').Append(i.ToString(CultureInfo.InvariantCulture)).Append('\n').Append(dseq[i]).Append('\n');
            dir["_dbd"] = db.ToString();
            var made = last.Run(lastDirectory + "/lastdb", G.Dorp == 'd' ? ["_dbd", "_dbd"] : ["-p", "_dbd", "_dbd"], dir);
            Report(made.StandardError);
            for (int k = 0; k < njob; k++)
            {
                if (G.Nthread > 0)
                    Report("\r" + k.ToString(CultureInfo.InvariantCulture) + " / " + njob.ToString(CultureInfo.InvariantCulture) + " (thread " + 0.ToString(CultureInfo.InvariantCulture) + ")                    \r");
                else
                    Report("\r" + k.ToString(CultureInfo.InvariantCulture) + " / " + njob.ToString(CultureInfo.InvariantCulture) + "                    \r");
                string kq = k.ToString(CultureInfo.InvariantCulture);
                string kd = kq;
                if (Math.Min(k, njob - nadd) == k)
                {
                    var dbk = new StringBuilder();
                    for (int i = 0; i < k; i++)
                        dbk.Append('>').Append(i.ToString(CultureInfo.InvariantCulture)).Append('\n').Append(dseq[i]).Append('\n');
                    dir["_db" + kd + "d"] = dbk.ToString();
                    var made1 = last.Run(lastDirectory + "/lastdb", G.Dorp == 'd' ? ["_db" + kd + "d", "_db" + kd + "d"] : ["-p", "_db" + kd + "d", "_db" + kd + "d"], dir);
                    Report(made1.StandardError);
                }
                else
                {
                    kd = ""; // the database of all the old sequences made before
                }

                dir["_q" + kq] = ">" + kq + "\n" + dseq[k] + "\n";
                int msize = Math.Max(10, k * lastm);
                var call = last.Run(lastDirectory + "/lastal", ["-m", msize.ToString(CultureInfo.InvariantCulture), "-e", laste.ToString(CultureInfo.InvariantCulture), "-f", "0", "-s", "1", "-p", "_scoringmatrixforlast",
                    "-a", (-G.Penalty).ToString(CultureInfo.InvariantCulture), "-b", (-G.PenaltyEx).ToString(CultureInfo.InvariantCulture), "_db" + kd + "d", "_q" + kq], dir);
                Report(call.StandardError);
                if (call.ExitCode != 0)
                    throw new MafftProgramFailureException("", 1);
                dir["_lastres" + kq] = call.StandardOutput;
                ReadLastResultSingleQuery(call.StandardOutput, k, lastresx, subopt, G.Penalty);
            }
        }
    }

    #endregion
}
