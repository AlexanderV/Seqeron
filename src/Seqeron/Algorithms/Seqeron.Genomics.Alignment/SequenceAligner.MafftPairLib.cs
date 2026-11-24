using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: tbfast with a pair library (--globalpair, --localpair, --genafpair without iteration)

    /// <summary><c>LocalHom</c> of <c>mltaln.h</c>: one ungapped block of a pairwise alignment (residue positions of the two sequences), its score (<c>opt</c>) and its importance; a chain of blocks per pair of sequences.</summary>
    internal sealed class MafftLocalHom
    {
        internal MafftLocalHom? Next;
        internal int Start1 = -1, End1 = -1, Start2 = -1, End2 = -1;
        internal double Opt = -1.0;
        internal int Overlapaa = -1;
        internal double Importance = -1.0, Rimportance;
        internal char Korh = 'h';
    }

    internal sealed partial class MafftGlobals
    {
        // pairlocalalign.c / tbfast.c
        internal int StoreDist = 1, StoreLocalhom = 1;
        internal int PpenaltyOp = NotSpecified, PpenaltyExGenaf = NotSpecified;
        internal bool UseNaiveScore;
        internal int Constraint;
        /// <summary><c>tbutree</c> of <c>tbfast</c>: 1 = the distances come from the pair stage / the files, 0 (<c>-J</c>) = from the alignment that is read.</summary>
        internal bool TbUtree = true;
        internal double FastaThreshold;
        internal int OutputHat23;

        /// <summary><c>arguments()</c> of <c>pairlocalalign.c</c> for the options the script <c>mafft</c> passes in <c>tbfast _ ... _</c>.</summary>
        internal static MafftGlobals ParsePairlocalalign(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified, Alg = 'A', Outgap = 1 };
            g.PenaltyShiftFactor = 1000.0;
            int k = 0;
            string Next() => k < argv.Count ? argv[k++] : throw new ArgumentException("A MAFFT option is missing its argument.");
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'O': g.PpenaltyOp = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'E': g.PpenaltyExGenaf = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q': g.PenaltyShiftFactor = CAtof(Next()); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; rest = true; break;
                        case 'l':
                            if (CAtof(Next()) < 0.00001)
                                g.StoreLocalhom = 0;
                            rest = true;
                            break;
                        case 'C': g.Nthread = CAtoi(Next()); rest = true; break;
                        case 'u': g.SpecificityConsideration = CAtof(Next()); rest = true; break;
                        case 'K': break;
                        case 'a': g.Fmodel = 1; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'L': g.Alg = 'L'; break;
                        case 'Y': g.Alg = 'Y'; break;
                        case 'I': g.Nadd = CAtoi(Next()); rest = true; break;
                        case 'N': g.Alg = 'N'; break;
                        case 'A': g.Alg = 'A'; break;
                        case 'Z': g.UseNaiveScore = true; break;
                        case ':': g.NWildcard = 1; break;
                        default: throw new NotSupportedException($"The MAFFT option -{c} of pairlocalalign is not available in this port.");
                    }
                }
            }
            return g;
        }

        /// <summary>
        /// <c>arguments()</c> of <c>tbfast.c</c> (the part after the pair-library arguments; the options that the script can produce for the progressive method with a pair library).
        /// Like the program, it reports <c>outputhat23</c> each time it reads <c>-+</c>.
        /// </summary>
        internal static MafftGlobals ParseTbfast(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified, Alg = 'A', Outgap = 1 };
            g.PenaltyShiftFactor = 1000.0;
            g.MinimumWeight = 0.00001;
            int k = 0;
            string Next() => k < argv.Count ? argv[k++] : throw new ArgumentException("A MAFFT option is missing its argument.");
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'V': g.PpenaltyDist = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q': g.PenaltyShiftFactor = CAtof(Next()); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; rest = true; break;
                        case 'l': g.FastaThreshold = CAtof(Next()); g.Constraint = 2; rest = true; break;
                        case 'C': g.NthreadPair = g.Nthread = CAtoi(Next()); rest = true; break;
                        case 's': g.SpecificityConsideration = CAtof(Next()); rest = true; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'y': g.DistOut = 1; break;
                        case 'S': g.SpScoreOut = 1; break;
                        case 't': g.TreeOut = 1; break;
                        case '^': g.TreeOut = 2; break;
                        case 'T': g.NoAlign = 1; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'L': g.LegacyGapCost = 1; break;
                        case 'O': g.Outgap = 0; break;
                        case 'X': g.TreeMethod = 'X'; g.SueffGlobal = CAtof(Next()); rest = true; break;
                        case 'n': g.OutNumber = 1; break;
                        case 'A': g.Alg = 'A'; break;
                        case 'M': g.Alg = 'M'; break;
                        case 'N': g.NeverMemSave = 1; break;
                        case 'B': break;
                        case 'K': g.AddProfile = false; break;
                        case 'H': g.Subalignment = true; g.SubalignmentOffset = CAtoi(Next()); rest = true; break;
                        case 'I': g.Nadd = CAtoi(Next()); rest = true; break;
                        case 'Y': g.KeepLength = 1; break;
                        case 'z': g.MapOut = 2; break;
                        case 'Z': g.MapOut = 1; break;
                        case 'p': g.Smoothing = 1; break;
                        case 'F': g.UseFft = true; break;
                        case 'G': g.ForceFft = true; g.UseFft = true; break;
                        case 'u': g.Tbrweight = 0; g.Weight = 0; break;
                        case 'U': g.TreeIn = 1; break;
                        case 'J': g.TbUtree = false; break;
                        case 'v': g.Tbrweight = 3; break;
                        case 'W': g.MinimumWeight = CAtof(Next()); rest = true; break;
                        case ':': g.NWildcard = 1; break;
                        case '@': g.Alg = 'd'; break;
                        case '+': g.OutputHat23 = CAtoi(Next()); g.Report("outputhat23=" + g.OutputHat23.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        default: throw new NotSupportedException($"The MAFFT option -{c} of tbfast is not available in this port.");
                    }
                }
            }
            if (g.FastaThreshold < 0.0001)
                g.Constraint = 0;
            return g;
        }
    }

    internal sealed partial class MafftRun
    {
        // -------------------------------------------------------------------- pairlocalalign.c (alg 'A': Galign11 on every pair)

        private static double Score2Dist(double pscore, double selfscore1, double selfscore2)
        {
            double bunbo = Math.Min(selfscore1, selfscore2);
            if (bunbo == 0.0)
                return 2.0;
            if (bunbo < pscore)
                return 0.0;
            return (1.0 - pscore / bunbo) * 2.0;
        }

        /// <summary>
        /// <c>pairlocalalign</c> as <c>tbfast</c> calls it (<c>pairlocalalign( njob, nlenmax, name, seq, iscore, localhomtable, ...)</c>): the pairwise alignment of every pair of sequences, the distance of every pair
        /// (<c>iscore[i][j - i]</c>) and the ungapped blocks of the alignments (<c>localhomtable[i][j - i]</c>). The messages go to <see cref="MafftGlobals.Log"/>.
        /// </summary>
        internal (double[][] Iscore, MafftLocalHom[][] LocalHom) PairLocalAlign(string[] seq, double[][]? expdist = null)
        {
            int njob = seq.Length;
            G.Constants(seq);
            var bseq = new string[njob];
            var dseq = new string[njob];
            var thereisx = new int[njob];
            if (CheckAlphabet(seq) is { } unknownLetter)
                throw unknownLetter;
            for (int i = 0; i < njob; i++)
            {
                bseq[i] = seq[i].Replace("-", "");
                if (G.Dorp == 'p' && G.Scoremtx == 1 && G.Nblosum > 0)
                {
                    int count = 0;
                    var sb = new System.Text.StringBuilder(bseq[i].Length);
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
            var localhom = new MafftLocalHom[njob][];
            for (int i = 0; i < njob; i++)
            {
                localhom[i] = new MafftLocalHom[njob - i];
                for (int j = 0; j < localhom[i].Length; j++)
                    localhom[i][j] = new MafftLocalHom();
            }
            var iscore = new double[njob][];
            for (int i = 0; i < njob; i++)
                iscore[i] = new double[njob - i];
            Report("All-to-all alignment.\n");
            var selfscore = new double[njob];
            for (int i = 0; i < njob; i++)
            {
                double s = 0.0;
                foreach (char c in bseq[i])
                    s += G.AminoDis[c][c];
                selfscore[i] = s;
            }
            if (G.Alg is not ('A' or 'L' or 'N'))
                throw new NotSupportedException("The pair library of this alignment method is not available in this port.");
            for (int i = 0; i < njob - 1; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                Report(MafftGlobals.D5(i) + " / " + njob + "\r");
                for (int j = i + 1; j < njob; j++)
                {
                    if (bseq[i].Length == 0 || bseq[j].Length == 0)
                    {
                        iscore[i][j - i] = 3.0;
                        continue;
                    }
                    var m1 = new[] { bseq[i] };
                    var m2 = new[] { bseq[j] };
                    int tx = thereisx[i] + thereisx[j];
                    double pscore;
                    int off1 = 0, off2 = 0;
                    if (G.Alg == 'A')
                    {
                        if (G.UseNaiveScore)
                        {
                            G__align11(G.NDisConsweightMulti, m1, m2, _nlenmax * 2, G.Outgap, G.Outgap);
                            pscore = NaivePairScore11(m1[0], m2[0], 0);
                        }
                        else if (G.StoreLocalhom != 0)
                        {
                            pscore = G__align11(G.NDisConsweightMulti, m1, m2, _nlenmax * 2, G.Outgap, G.Outgap);
                            if (tx != 0)
                                pscore = G__align11Noalign(G.NDisConsweightMulti, dseq[i], dseq[j]);
                            if (G.SpecificityConsideration > 0.0)
                            {
                                // --unalignlevel: a pair that is distant enough is aligned again with the scoring matrix shifted by the offset of its distance (the score stays that of the first alignment, the blocks come from the second)
                                double dist = expdist != null ? expdist[i][j] : Score2Dist(pscore, selfscore[i], selfscore[j]);
                                if (dist * 0.5 - G.SpecificityConsideration < 0.0)
                                {
                                    var dynamicmtx = MakeDynamicMtx(G.NDisConsweightMulti, 0.5 * dist);
                                    m1[0] = bseq[i];
                                    m2[0] = bseq[j];
                                    G__align11(dynamicmtx, m1, m2, _nlenmax * 2, G.Outgap, G.Outgap);
                                }
                            }
                        }
                        else
                        {
                            // no pair library is kept (--weighti 0): the score only, G__align11_noalign
                            pscore = G__align11Noalign(G.NDisConsweightMulti, dseq[i], dseq[j]);
                        }
                    }
                    else if (G.Alg == 'N')
                    {
                        if (G.UseNaiveScore)
                        {
                            GenL__align11(G.NDisConsweightMulti, m1, m2, _nlenmax * 2, out off1, out off2);
                            pscore = NaivePairScore11(m1[0], m2[0], 0);
                        }
                        else
                        {
                            pscore = GenL__align11(G.NDisConsweightMulti, m1, m2, _nlenmax * 2, out off1, out off2);
                            if (tx != 0)
                                pscore = GenL__align11(G.NDisConsweightMulti, new[] { dseq[i] }, new[] { dseq[j] }, _nlenmax * 2, out _, out _);
                        }
                    }
                    else if (G.UseNaiveScore)
                    {
                        L__align11(G.NDisConsweightMulti, 0.0, m1, m2, _nlenmax * 2, true, out off1, out off2);
                        pscore = NaivePairScore11(m1[0], m2[0], 0);
                    }
                    else if (G.StoreLocalhom != 0)
                    {
                        pscore = L__align11(G.NDisConsweightMulti, 0.0, m1, m2, _nlenmax * 2, true, out off1, out off2);
                        if (tx != 0)
                            pscore = L__align11(G.NDisConsweightMulti, 0.0, new[] { dseq[i] }, new[] { dseq[j] }, _nlenmax * 2, false, out _, out _);
                    }
                    else
                    {
                        // no pair library is kept (--weighti 0): the score only, L__align11_noalign (it keeps no trace-back matrix and does not touch commonIP)
                        pscore = L__align11(G.NDisConsweightMulti, 0.0, new[] { dseq[i] }, new[] { dseq[j] }, _nlenmax * 2, false, out _, out _);
                    }
                    if (m1[0].Length != 0 && m2[0].Length != 0)
                    {
                        if (G.StoreLocalhom != 0)
                            PutLocalHom2(m1[0], m2[0], localhom[i][j - i], off1, off2);
                        pscore = Score2Dist(pscore, selfscore[i], selfscore[j]);
                    }
                    else
                    {
                        pscore = 2.0;
                    }
                    iscore[i][j - i] = pscore;
                }
            }
            Report("tbfast-pair (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return (iscore, localhom);
        }

        /// <summary><c>G__align11_noalign</c> with the end gaps penalised (<c>outgap</c> 1): the score of the global alignment of two sequences without the alignment (the same dynamic program, the value of the last cell).</summary>
        private double G__align11Noalign(double[][] dynamicmtx, string s1, string s2)
        {
            if (s1.Length == 0 || s2.Length == 0)
            {
                // the program's noalign variant reads the terminating NUL: the score of an empty sequence is the value of the empty dynamic program
                return 0.0;
            }
            var a = new[] { s1 };
            var b = new[] { s2 };
            return G__align11(dynamicmtx, a, b, _nlenmax * 2, 1, 1, noalign: true);
        }

        /// <summary>
        /// <c>L__align11</c> / <c>L__align11_noalign</c> (<c>Lalign11.c</c>): the local alignment of two sequences (the best scoring cell of the matrix, cells below <c>-offset</c> restart); with <paramref name="withAlignment"/> the two
        /// sequences are replaced by the aligned segments (empty when nothing scores above the threshold) and <paramref name="off1"/> / <paramref name="off2"/> are the positions of their first residues. Returns the score.
        /// </summary>
        internal double L__align11(double[][] dynamicmtx, double scoreoffset, string[] seq1, string[] seq2, int alloclen, bool withAlignment, out int off1, out int off2)
        {
            off1 = off2 = 0;
            string s1 = seq1[0], s2 = seq2[0];
            int lgth1 = s1.Length, lgth2 = s2.Length;
            if (G.TryWarp)
                throw new NotSupportedException("--allowshift / --shiftpenalty (the warp of L__align11) is not available in this port.");
            var amino = new double[256][];
            for (int i = 0; i < 256; i++)
                amino[i] = new double[256];
            for (int i = 0; i < G.Nalphabets; i++)
                for (int j = 0; j < G.Nalphabets; j++)
                    amino[G.Amino[i]][G.Amino[j]] = dynamicmtx[i][j];
            double localthr = withAlignment ? -G.Offset + scoreoffset * 600 : -G.Offset;
            double fpenalty = G.Penalty, fpenaltyEx = G.PenaltyEx;
            int[][]? ijp = null;
            if (withAlignment)
            {
                RequireCells(((long)lgth1 + 1) * ((long)lgth2 + 1), MafftMaxPairCells, nameof(MafftMaxPairCells), "MAFFT", "keeps an int trace-back matrix ijp[len1 + 1][len2 + 1] (L__align11)",
                    string.Create(CultureInfo.InvariantCulture, $"{lgth1:N0} x {lgth2:N0} residues = {((long)lgth1 + 1) * ((long)lgth2 + 1):N0} cells"), "The program has no limit and is killed when the allocation (4 bytes per cell plus the row arrays) fails", "use Kalign for sequences of this length");
                _lCap.Fit(lgth1, lgth2);
                ijp = G.Common.IpFor(_lCap.Orlgth1, _lCap.Orlgth2);
            }
            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var m = new double[lgth2 + 2];
            var mp = new int[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            double[] currentw = w1, previousw = w2;
            // an empty sequence reads the terminating NUL of the C string (its row of the matrix is zero)
            char s1First = lgth1 > 0 ? s1[0] : '\0', s2First = lgth2 > 0 ? s2[0] : '\0';
            for (int j = 0; j < lgth1; j++)
                initverticalw[j] = amino[s2First][s1[j]];
            for (int j = 0; j < lgth2; j++)
                currentw[j] = amino[s1First][s2[j]];
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1];
                mp[j] = 0;
            }
            int localstop = lgth1 + lgth2 + 1;
            double maxwm = -999999999.9;
            int endali = 0, endalj = 0;
            for (int i = 1; i < lgth1 + 1; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                double[] aminoRow = amino[i < lgth1 ? s1[i] : '\0']; // the terminating NUL of the sequence scores 0
                for (int j = 0; j < lgth2; j++)
                    currentw[j] = aminoRow[s2[j]];
                currentw[0] = initverticalw[i];
                double mi = previousw[0];
                int mpi = 0;
                var ijpRow = ijp?[i];
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    double wm = previousw[j - 1];
                    int ptr = 0;
                    double g;
                    if ((g = mi + fpenalty) > wm)
                    {
                        wm = g;
                        ptr = -(j - mpi);
                    }
                    if (previousw[j - 1] > mi)
                    {
                        mi = previousw[j - 1];
                        mpi = j - 1;
                    }
                    mi += fpenaltyEx;
                    if ((g = m[j] + fpenalty) > wm)
                    {
                        wm = g;
                        ptr = +(i - mp[j]);
                    }
                    if (previousw[j - 1] > m[j])
                    {
                        m[j] = previousw[j - 1];
                        mp[j] = i - 1;
                    }
                    m[j] += fpenaltyEx;
                    if (maxwm < wm)
                    {
                        maxwm = wm;
                        endali = i;
                        endalj = j;
                    }
                    if (wm < localthr)
                    {
                        ptr = localstop;
                        wm = localthr;
                    }
                    if (ijpRow != null)
                        ijpRow[j] = ptr;
                    currentw[j] += wm;
                }
            }
            if (!withAlignment)
            {
                // L__align11_noalign
                if (maxwm < localthr)
                {
                    Report("maxwm <- 0.0 \n");
                    return 0.0;
                }
                return maxwm;
            }
            if (ijp![endali][endalj] == localstop)
            {
                seq1[0] = "";
                seq2[0] = "";
                Report("maxwm <- 0.0 \n");
                return 0.0;
            }
            // Ltracking
            for (int i = 0; i < lgth1 + 1; i++)
                ijp[i][0] = localstop;
            for (int j = 0; j < lgth2 + 1; j++)
                ijp[0][j] = localstop;
            var a = new char[lgth1 + lgth2 + 1];
            var b = new char[lgth1 + lgth2 + 1];
            int pa = lgth1 + lgth2, pb = lgth1 + lgth2;
            int iin = endali, jin = endalj;
            int limk = lgth1 + lgth2;
            int ifi = 0, jfi = 0;
            int k = 0;
            while (k <= limk)
            {
                int v = ijp[iin][jin];
                if (v < 0)
                {
                    ifi = iin - 1;
                    jfi = jin + v;
                }
                else if (v > 0)
                {
                    ifi = iin - v;
                    jfi = jin - 1;
                }
                else
                {
                    ifi = iin - 1;
                    jfi = jin - 1;
                }
                int l = iin - ifi;
                while (--l > 0)
                {
                    a[--pa] = s1[ifi + l];
                    b[--pb] = '-';
                    k++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    a[--pa] = '-';
                    b[--pb] = s2[jfi + l];
                    k++;
                }
                if (iin <= 0 || jin <= 0)
                    break;
                a[--pa] = s1[ifi];
                b[--pb] = s2[jfi];
                if (ijp[ifi][jfi] == localstop)
                    break;
                k++;
                iin = ifi;
                jin = jfi;
                k++;
            }
            off1 = ifi == -1 ? 0 : ifi;
            off2 = jfi == -1 ? 0 : jfi;
            if (alloclen < lgth1 + lgth2 - pa)
                throw new InvalidOperationException("LENGTH OVER!");
            seq1[0] = new string(a, pa, lgth1 + lgth2 - pa);
            seq2[0] = new string(b, pb, lgth1 + lgth2 - pb);
            return maxwm;
        }

        /// <summary>
        /// <c>genL__align11</c> (<c>genalign11.c</c>): the local alignment of two sequences with generalised affine gap costs (Altschul 1998): besides the gap of the ordinary affine model, a second, more costly opening
        /// (<c>penalty_OP</c>) lets the alignment jump from any cell to any earlier cell of the same row or column through the best scoring one (long unalignable regions). The same result convention as <see cref="L__align11"/>.
        /// </summary>
        internal double GenL__align11(double[][] dynamicmtx, string[] seq1, string[] seq2, int alloclen, out int off1, out int off2)
        {
            off1 = off2 = 0;
            string s1 = seq1[0], s2 = seq2[0];
            int lgth1 = s1.Length, lgth2 = s2.Length;
            var amino = new double[256][];
            for (int i = 0; i < 256; i++)
                amino[i] = new double[256];
            for (int i = 0; i < G.Nalphabets; i++)
                for (int j = 0; j < G.Nalphabets; j++)
                    amino[G.Amino[i]][G.Amino[j]] = dynamicmtx[i][j];
            double fpenalty = G.Penalty, fpenaltyOp = G.PenaltyOp, fpenaltyEx = G.PenaltyEx;
            double localthr = -(double)G.Offset;
            RequireCells(2 * ((long)lgth1 + 1) * ((long)lgth2 + 1), MafftMaxPairCells, nameof(MafftMaxPairCells), "MAFFT", "keeps two int trace-back matrices ijpi, ijpj[len1 + 1][len2 + 1] (genL__align11)",
                string.Create(CultureInfo.InvariantCulture, $"{lgth1:N0} x {lgth2:N0} residues = {((long)lgth1 + 1) * ((long)lgth2 + 1):N0} cells"), "The program has no limit and is killed when the allocation (8 bytes per cell plus the row arrays) fails", "use Kalign for sequences of this length");
            _genCap.Fit(lgth1, lgth2);
            var ijpi = G.Common.IpFor(_genCap.Orlgth1, _genCap.Orlgth2, withJp: true);
            var ijpj = G.Common.JpFor(lgth1, lgth2);
            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var m = new double[lgth2 + 2];
            var mp = new int[lgth2 + 2];
            var largeM = new double[lgth2 + 2];
            var largeMp = new int[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            double[] currentw = w1, previousw = w2;
            for (int j = 0; j < lgth1; j++)
                initverticalw[j] = amino[lgth2 > 0 ? s2[0] : '\0'][s1[j]];
            for (int j = 0; j < lgth2; j++)
                currentw[j] = amino[lgth1 > 0 ? s1[0] : '\0'][s2[j]];
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1];
                mp[j] = 0;
                largeM[j] = currentw[j - 1];
                largeMp[j] = 0;
            }
            int localstop = lgth1 + lgth2 + 1;
            double maxwm = -999999999.9;
            int endali = 0, endalj = 0;
            for (int i = 1; i < lgth1 + 1; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                double[] aminoRow = amino[i < lgth1 ? s1[i] : '\0'];
                for (int j = 0; j < lgth2; j++)
                    currentw[j] = aminoRow[s2[j]];
                currentw[0] = initverticalw[i];
                double mi = previousw[0];
                int mpi = 0;
                double largeMi = previousw[0];
                int largeMpi = 0;
                var rowI = ijpi[i];
                var rowJ = ijpj[i];
                double tbk = -999999.9;
                int tbki = 0, tbkj = 0;
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    double wm = previousw[j - 1];
                    rowI[j] = i - 1;
                    rowJ[j] = j - 1;
                    double g = mi + fpenalty;
                    if (g > wm)
                    {
                        wm = g;
                        rowJ[j] = mpi;
                    }
                    g = previousw[j - 1];
                    if (g > mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    mi += fpenaltyEx;
                    g = m[j] + fpenalty;
                    if (g > wm)
                    {
                        wm = g;
                        rowI[j] = mp[j];
                        rowJ[j] = j - 1;
                    }
                    g = previousw[j - 1];
                    if (g > m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    m[j] += fpenaltyEx;
                    g = tbk + fpenaltyOp;
                    if (g > wm)
                    {
                        wm = g;
                        rowI[j] = tbki;
                        rowJ[j] = tbkj;
                    }
                    if (largeMi > tbk)
                    {
                        tbk = largeMi;
                        tbki = i - 1;
                        tbkj = largeMpi;
                    }
                    if (largeM[j] > tbk)
                    {
                        tbk = largeM[j];
                        tbki = largeMp[j];
                        tbkj = j - 1;
                    }
                    if (previousw[j - 1] > largeM[j])
                    {
                        largeM[j] = previousw[j - 1];
                        largeMp[j] = i - 1;
                    }
                    if (previousw[j - 1] > largeMi)
                    {
                        largeMi = previousw[j - 1];
                        largeMpi = j - 1;
                    }
                    if (maxwm < wm)
                    {
                        maxwm = wm;
                        endali = i;
                        endalj = j;
                    }
                    if (wm < localthr)
                    {
                        rowI[j] = localstop;
                        wm = localthr;
                    }
                    currentw[j] += wm;
                }
            }
            if (ijpi[endali][endalj] == localstop)
            {
                seq1[0] = "";
                seq2[0] = "";
                return 0.0;
            }
            for (int i = 0; i < lgth1 + 1; i++)
            {
                ijpi[i][0] = localstop;
                ijpj[i][0] = localstop;
            }
            for (int j = 0; j < lgth2 + 1; j++)
            {
                ijpi[0][j] = localstop;
                ijpj[0][j] = localstop;
            }
            var a = new char[lgth1 + lgth2 + 1];
            var b = new char[lgth1 + lgth2 + 1];
            int pa = lgth1 + lgth2, pb = lgth1 + lgth2;
            int iin = endali, jin = endalj;
            int limk = lgth1 + lgth2;
            int ifi = 0, jfi = 0;
            int k = 0;
            while (k <= limk)
            {
                ifi = ijpi[iin][jin];
                jfi = ijpj[iin][jin];
                int l = iin - ifi;
                while (--l > 0)
                {
                    a[--pa] = s1[ifi + l];
                    b[--pb] = '-';
                    k++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    a[--pa] = '-';
                    b[--pb] = s2[jfi + l];
                    k++;
                }
                if (iin <= 0 || jin <= 0)
                    break;
                a[--pa] = s1[ifi];
                b[--pb] = s2[jfi];
                if (ijpi[ifi][jfi] == localstop || ijpj[ifi][jfi] == localstop)
                    break;
                k++;
                iin = ifi;
                jin = jfi;
                k++;
            }
            off1 = ifi == -1 ? 0 : ifi;
            off2 = jfi == -1 ? 0 : jfi;
            if (alloclen < lgth1 + lgth2 - pa)
                throw new InvalidOperationException("LENGTH OVER!");
            seq1[0] = new string(a, pa, lgth1 + lgth2 - pa);
            seq2[0] = new string(b, pb, lgth1 + lgth2 - pb);
            return maxwm;
        }

        /// <summary><c>putlocalhom2</c> (<c>io.c</c>, <c>divpairscore</c> 0): the ungapped runs of the aligned pair become blocks; every block gets the mean score per aligned column of all blocks (<c>opt</c>).</summary>
        internal void PutLocalHom2(string al1, string al2, MafftLocalHom head, int off1, int off2, char korh = 'h')
        {
            var tmppt = head;
            int pos1 = off1, pos2 = off2, start1 = 0, start2 = 0, end1, end2;
            int isumscore = 0, sumoverlap = 0;
            bool first = true;
            int st = 0, iscore = 0;
            int n = al1.Length;
            for (int p = 0; p < n; p++)
            {
                char c1 = al1[p], c2 = al2[p];
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
                    tmppt.Korh = korh;
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
                    iscore += G.NDis[G.AminoN[c1]][G.AminoN[c2]];
                }
                if (c1 != '-')
                    pos1++;
                if (c2 != '-')
                    pos2++;
            }
            // for an empty alignment the program reads the byte before its buffer (the high byte of the malloc size field, 0): not a gap
            if (n == 0 || (al1[n - 1] != '-' && al2[n - 1] != '-'))
            {
                if (!first)
                {
                    tmppt.Next = new MafftLocalHom();
                    tmppt = tmppt.Next;
                }
                end1 = pos1 - 1;
                end2 = pos2 - 1;
                tmppt.Start1 = start1;
                tmppt.Start2 = start2;
                tmppt.End1 = end1;
                tmppt.End2 = end2;
                tmppt.Korh = korh;
                isumscore += iscore;
                sumoverlap += end2 - start2 + 1;
            }
            for (tmppt = head; tmppt != null; tmppt = tmppt.Next)
            {
                tmppt.Overlapaa = sumoverlap;
                if (tmppt.Overlapaa > 0)
                    tmppt.Opt = (double)isumscore * 5.8 / (600 * sumoverlap);
                else
                    tmppt.Opt = -1.0;
            }
        }

        // -------------------------------------------------------------------- mltaln9.c: the importance of the blocks and the match matrix of a pair of groups

        /// <summary><c>calcimportance_half</c>: the importance of every block = the mean over its residues of the weighted number of sequences that have a block there, times the score of the block (the two sequences' values averaged).</summary>
        private static void CalcImportanceHalf(int nseq, double[] eff, string[] seq, MafftLocalHom[][] localhom, int alloclen)
        {
            var importance = new double[Math.Max(alloclen, 1)];
            var ieff = new double[nseq];
            double totaleff = 0.0;
            for (int i = 0; i < nseq; i++)
            {
                int nogaplen = 0;
                foreach (char c in seq[i])
                    if (c != '-')
                        nogaplen++;
                ieff[i] = nogaplen == 0 ? 0.0 : eff[i];
                totaleff += ieff[i];
            }
            for (int i = 0; i < nseq; i++)
                ieff[i] /= totaleff;
            for (int i = 0; i < nseq; i++)
            {
                Array.Clear(importance);
                for (int j = 0; j < nseq; j++)
                {
                    if (i == j)
                        continue;
                    if (i < j)
                    {
                        for (var t = localhom[i][j - i]; t != null; t = t.Next)
                        {
                            if (t.Opt == -1)
                                continue;
                            for (int pos = t.Start1; pos <= t.End1; pos++)
                                importance[pos] += ieff[j];
                        }
                    }
                    else
                    {
                        for (var t = localhom[j][i - j]; t != null; t = t.Next)
                        {
                            if (t.Opt == -1)
                                continue;
                            for (int pos = t.Start2; pos <= t.End2; pos++)
                                importance[pos] += ieff[j];
                        }
                    }
                }
                for (int j = 0; j < nseq; j++)
                {
                    if (i == j)
                        continue;
                    if (i < j)
                    {
                        if (localhom[i][j - i].Opt == -1.0)
                            continue;
                        for (var t = localhom[i][j - i]; t != null; t = t.Next)
                        {
                            if (t.Opt == -1.0)
                                continue;
                            double tmpdouble = 0.0;
                            int len = 0;
                            for (int pos = t.Start1; pos <= t.End1; pos++)
                            {
                                tmpdouble += importance[pos];
                                len++;
                            }
                            tmpdouble /= len;
                            t.Importance = tmpdouble * t.Opt;
                        }
                    }
                    else
                    {
                        if (localhom[j][i - j].Opt == -1.0)
                            continue;
                        for (var t = localhom[j][i - j]; t != null; t = t.Next)
                        {
                            if (t.Opt == -1.0)
                                continue;
                            double tmpdouble = 0.0;
                            int len = 0;
                            for (int pos = t.Start2; pos <= t.End2; pos++)
                            {
                                tmpdouble += importance[pos];
                                len++;
                            }
                            tmpdouble /= len;
                            t.Rimportance = tmpdouble * t.Opt;
                        }
                    }
                }
            }
            for (int i = 0; i < nseq - 1; i++)
            {
                for (int j = i + 1; j < nseq; j++)
                {
                    for (var t = localhom[i][j - i]; t != null; t = t.Next)
                    {
                        if (t.Opt == -1.0)
                            continue;
                        double imp = 0.5 * (t.Importance + t.Rimportance);
                        t.Importance = t.Rimportance = imp;
                    }
                }
            }
        }

        /// <summary><c>movereg</c> / <c>movereg_swap</c>: the columns of the two group rows that hold the residue positions of a block (the first sequence is the one of the first residue range unless <paramref name="swap"/>).</summary>
        private static void MoveReg(string seq1, string seq2, MafftLocalHom t, bool swap, out int start1, out int start2, out int end1, out int end2)
        {
            int s1 = swap ? t.Start2 : t.Start1, e1 = swap ? t.End2 : t.End1;
            int s2 = swap ? t.Start1 : t.Start2, e2 = swap ? t.End1 : t.End2;
            ScanRange(seq1, s1, e1, out start1, out end1);
            ScanRange(seq2, s2, e2, out start2, out end2);
        }

        private static void ScanRange(string seq, int startRes, int endRes, out int startCol, out int endCol)
        {
            int pt = 0;
            int tmpint = -1;
            while (pt < seq.Length)
            {
                if (seq[pt++] != '-')
                    tmpint++;
                if (tmpint == startRes)
                    break;
            }
            startCol = pt - 1;
            if (startRes == endRes)
            {
                endCol = startCol;
                return;
            }
            while (pt < seq.Length)
            {
                if (seq[pt++] != '-')
                    tmpint++;
                if (tmpint == endRes)
                    break;
            }
            endCol = pt - 1;
        }

        private double[][]? _impMtx;
        private int _impAllocLen;
        private double[][]? _partImpMtx;
        private int _partImpAllocLen;
        private double[][]? _dImpMtx;
        private int _dImpAllocLen;

        /// <summary><c>imp_match_init_strict</c> / <c>fillimp</c>: the matrix of the importance of every pair of columns of two groups, from the blocks of every pair of their sequences, weighted by the weights of the sequences and <c>fastathreshold</c>.</summary>
        private double[][] FillImp(int clus1, int clus2, int lgth1, int lgth2, string[] seq1, string[] seq2, double[] eff1, double[] eff2, MafftLocalHom?[][] localhom, int[] orinum1, int[] orinum2, bool part = false, double[]? eff1Kozo = null, double[]? eff2Kozo = null, MafftLocalHom?[][]? seedShrink = null, int[]? seedIn1 = null, int[]? seedIn2 = null, bool dmtx = false)
        {
            // impmtx of Salignmm.c (partSalignmm.c for the refinement, <paramref name="part"/>) is static: it is allocated (zeroed) when a group is longer than any before and only [0, lgth1) x [0, lgth2) is cleared by a call, so the row lgth1 read by the
            // dynamic program after the last residue (end gaps penalised) holds what earlier, larger calls left there
            if (dmtx)
            {
                // imp_match_init_strictD: the matrix of Dalignmm.c (its own static)
                if (_dImpMtx == null || _dImpAllocLen < lgth1 + 2 || _dImpAllocLen < lgth2 + 2)
                {
                    _dImpAllocLen = Math.Max(lgth1, lgth2) + 2;
                    _dImpMtx = new double[_dImpAllocLen][];
                    for (int i = 0; i < _dImpAllocLen; i++)
                        _dImpMtx[i] = new double[_dImpAllocLen];
                }
            }
            else if (part)
            {
                if (_partImpMtx == null || _partImpAllocLen < lgth1 + 2 || _partImpAllocLen < lgth2 + 2)
                {
                    _partImpAllocLen = Math.Max(lgth1, lgth2) + 2;
                    _partImpMtx = new double[_partImpAllocLen][];
                    for (int i = 0; i < _partImpAllocLen; i++)
                        _partImpMtx[i] = new double[_partImpAllocLen];
                }
            }
            else if (_impMtx == null || _impAllocLen < lgth1 + 2 || _impAllocLen < lgth2 + 2)
            {
                _impAllocLen = Math.Max(lgth1, lgth2) + 2;
                _impMtx = new double[_impAllocLen][];
                for (int i = 0; i < _impAllocLen; i++)
                    _impMtx[i] = new double[_impAllocLen];
            }
            double[][] impmtx;
            if (dmtx)
                impmtx = _dImpMtx!;
            else
                impmtx = part ? _partImpMtx! : _impMtx!;
            for (int i = 0; i < lgth1; i++)
                Array.Clear(impmtx[i], 0, lgth2);
            double effijx = 1.0 * G.FastaThreshold;
            if (seedShrink != null)
            {
                // fillimp_file (nodepair): the blocks of the seeds come first, weighted by the structure weights; a block of another kind is a defect of the file (the program stops)
                for (int i = 0; i < clus1; i++)
                {
                    if (seedIn1![i] == -1)
                        continue;
                    for (int j = 0; j < clus2; j++)
                    {
                        if (seedIn2![j] == -1)
                            continue;
                        double effijKozo = (eff1Kozo?[i] ?? 0.0) * (eff2Kozo?[j] ?? 0.0) * effijx;
                        bool swapSeed = orinum1[i] > orinum2[j];
                        for (var t = seedShrink[seedIn1[i]][seedIn2[j]]; t != null; t = t.Next)
                        {
                            MoveReg(seq1[i], seq2[j], t, swapSeed, out int start1, out int start2, out int end1, out int end2);
                            int k1 = start1, k2 = start2;
                            string s1 = seq1[i], s2 = seq2[j];
                            while (k1 < s1.Length && k2 < s2.Length)
                            {
                                char c1 = s1[k1], c2 = s2[k2];
                                if (c1 != '-' && c2 != '-')
                                {
                                    if (t.Korh != 'k')
                                        throw new MafftProgramFailureException("okashii\n");
                                    impmtx[k1][k2] += t.Importance * effijKozo;
                                    k1++;
                                    k2++;
                                }
                                else if (c1 != '-' && c2 == '-')
                                {
                                    k2++;
                                }
                                else if (c1 == '-' && c2 != '-')
                                {
                                    k1++;
                                }
                                else
                                {
                                    k1++;
                                    k2++;
                                }
                                if (k1 > end1 || k2 > end2)
                                    break;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < clus1; i++)
            {
                bool swap = false;
                for (int j = 0; j < clus2; j++)
                {
                    swap = orinum1[i] > orinum2[j];
                    double effij = eff1[i] * eff2[j] * effijx;
                    double effijKozo = (eff1Kozo?[i] ?? 0.0) * (eff2Kozo?[j] ?? 0.0) * effijx; // the blocks of kind 'k' (a structure) are weighted by the structure weights
                    for (var t = localhom[i][j]; t != null; t = t.Next)
                    {
                        MoveReg(seq1[i], seq2[j], t, swap, out int start1, out int start2, out int end1, out int end2);
                        int k1 = start1, k2 = start2;
                        string s1 = seq1[i], s2 = seq2[j];
                        while (k1 < s1.Length && k2 < s2.Length)
                        {
                            char c1 = s1[k1], c2 = s2[k2];
                            if (c1 != '-' && c2 != '-')
                            {
                                impmtx[k1][k2] += t.Importance * (t.Korh == 'k' ? effijKozo : effij);
                                k1++;
                                k2++;
                            }
                            else if (c1 != '-' && c2 == '-')
                            {
                                k2++;
                            }
                            else if (c1 == '-' && c2 != '-')
                            {
                                k1++;
                            }
                            else
                            {
                                k1++;
                                k2++;
                            }
                            if (k1 > end1 || k2 > end2)
                                break;
                        }
                    }
                }
            }
            return impmtx;
        }

        /// <summary><c>fastshrinklocalhom_half_seed</c>: the positions of the members of two groups among the seeds (the members below <paramref name="nseed"/>, -1 for the others; the entries behind the groups keep their old values, as the program's
        /// arrays do) and the blocks of every pair of seeds (null when the pair has none).</summary>
        private static MafftLocalHom?[][] ShrinkSeedHom(int[] mem1, int[] mem2, int nseed, int[] pos1, int[] pos2, MafftLocalHom[][] localhom)
        {
            int sk1 = 0;
            for (int k1 = 0; k1 < mem1.Length; k1++)
                pos1[k1] = mem1[k1] >= nseed ? -1 : sk1++;
            int sk2 = 0;
            for (int k2 = 0; k2 < mem2.Length; k2++)
                pos2[k2] = mem2[k2] >= nseed ? -1 : sk2++;
            var shrink = new MafftLocalHom?[sk1][];
            int r = 0;
            foreach (int a in mem1)
            {
                if (a >= nseed)
                    continue;
                shrink[r] = new MafftLocalHom?[sk2];
                int c = 0;
                foreach (int b in mem2)
                {
                    if (b >= nseed)
                        continue;
                    var head = a < b ? localhom[a][b - a] : localhom[b][a - b];
                    shrink[r][c++] = head.Opt == -1 ? null : head;
                }
                r++;
            }
            return shrink;
        }

        /// <summary><c>fastshrinklocalhom_half</c>: the blocks of every pair (member of group 1, member of group 2); null when the pair has none.</summary>
        private static MafftLocalHom?[][] ShrinkLocalHom(int[] mem1, int[] mem2, MafftLocalHom[][] localhom)
        {
            var shrink = new MafftLocalHom?[mem1.Length][];
            for (int k1 = 0; k1 < mem1.Length; k1++)
            {
                shrink[k1] = new MafftLocalHom?[mem2.Length];
                for (int k2 = 0; k2 < mem2.Length; k2++)
                {
                    int a = mem1[k1], b = mem2[k2];
                    var head = a < b ? localhom[a][b - a] : localhom[b][a - b];
                    shrink[k1][k2] = head.Opt == -1 ? null : head;
                }
            }
            return shrink;
        }

        // -------------------------------------------------------------------- tbfast.c

        /// <summary>The pair library as <c>dvtditr</c> reads it from the file <c>hat3</c> (<c>-+ n</c> with n &gt; 0 writes it, five decimals), null without <c>-+</c>.</summary>
        internal MafftLocalHom[][]? Hat3Table;

        /// <summary><c>kozoarivec</c>: the sequences that a block of kind 'k' (a structure) names; null when there is none.</summary>
        internal bool[]? Kozo;

        /// <summary>The distances of the pair library as <c>dvtditr</c> reads them from the file <c>hat2</c> (six-character fields with three decimals), [i][j] for i &lt; j.</summary>
        internal double[][]? Hat2Table;

        /// <summary><c>tbfast</c> writes the blocks as <c>"%d %d %d %7.5f %d %d %d %d %c"</c> with <c>opt / 600 * 5.8</c>; <c>readlocalhomtable2_half</c> reads <c>opt / 5.8 * 600</c> back into fresh blocks.</summary>
        private static MafftLocalHom[][] RoundTripHat3(MafftLocalHom[][] lh, int njob)
        {
            var res = new MafftLocalHom[njob][];
            for (int i = 0; i < njob; i++)
            {
                res[i] = new MafftLocalHom[njob - i];
                for (int j = 0; j < res[i].Length; j++)
                    res[i][j] = new MafftLocalHom();
            }
            for (int i = 0; i < njob - 1; i++)
            {
                for (int j = i; j < njob; j++)
                {
                    bool first = true;
                    var tail = res[i][j - i];
                    for (var t = lh[i][j - i]; t != null; t = t.Next)
                    {
                        if (t.Opt == -1.0)
                            continue;
                        if (j <= i)
                            continue; // the program writes the pairs of different sequences with the first younger than the second only
                        double opt = MafftGlobals.CAtof((t.Opt / 600 * 5.8).ToString("F5", CultureInfo.InvariantCulture));
                        MafftLocalHom node;
                        if (first)
                        {
                            node = tail;
                            first = false;
                        }
                        else
                        {
                            node = new MafftLocalHom { Importance = 0.0, Overlapaa = 0 };
                            tail.Next = node;
                            tail = node;
                        }
                        node.Start1 = t.Start1;
                        node.Start2 = t.Start2;
                        node.End1 = t.End1;
                        node.End2 = t.End2;
                        node.Opt = (opt + 0.00) / 5.8 * 600;
                        node.Overlapaa = t.Overlapaa;
                        node.Korh = t.Korh;
                    }
                }
            }
            return res;
        }

        /// <summary><c>WriteFloatHat2_pointer_halfmtx</c> (<c>%#6.3f</c> values, twelve to a line) and <c>readhat2_pointer</c> (<c>input_new</c>: an optional line end, then six characters): the distances as <c>dvtditr</c> gets them.</summary>
        private static double[][] RoundTripHat2(double[][] iscore, int njob)
        {
            var sb = new System.Text.StringBuilder();
            for (int i = 0; i < njob; i++)
            {
                for (int j = i + 1; j < njob; j++)
                {
                    sb.Append(iscore[i][j - i].ToString("F3", CultureInfo.InvariantCulture).PadLeft(6));
                    int ijsa = j - i;
                    if (ijsa % 12 == 0 || ijsa == njob - i - 1)
                        sb.Append('\n');
                }
            }
            string text = sb.ToString();
            int pos = 0;
            var res = new double[njob][];
            for (int i = 0; i < njob; i++)
                res[i] = new double[njob];
            for (int i = 0; i < njob - 1; i++)
            {
                for (int j = i + 1; j < njob; j++)
                {
                    if (pos < text.Length && text[pos] == '\n')
                        pos++;
                    int take = Math.Min(6, text.Length - pos);
                    res[i][j] = MafftGlobals.CAtof(text.Substring(pos, take));
                    pos += take;
                }
            }
            return res;
        }

        /// <summary>The distance matrix of <c>tbfast -J</c> from the alignment that is read (<c>naivepairscorefast</c> of every pair, 2 (1 - score / min(self scores)), at most 10); the rows must be equally long.</summary>
        private double[][] TbfastMsaDistances(string[] seq, int njob)
        {
            Report("Making a distance matrix from msa .. \n");
            for (int i = 1; i < njob; i++)
                if (_nlen[i] != _nlen[0])
                    throw new MafftProgramFailureException("Input pre-aligned seqences or make hat2.\n");
            var skip = new int[njob][];
            for (int i = 0; i < njob; i++)
                skip[i] = MakeSkipTable(seq[i]);
            var selfscore = new int[njob];
            for (int i = 0; i < njob; i++)
                selfscore[i] = (int)NaivePairScoreFast(seq[i], seq[i], skip[i], skip[i], G.PenaltyDist);
            var iscore = new double[njob][];
            for (int i = 0; i < njob; i++)
                iscore[i] = new double[njob - i];
            int ien = njob - 1;
            for (int i = 0; i < ien; i++)
            {
                if (i % 10 == 0)
                    Report("\r" + MafftGlobals.D5(i) + " / " + ien.ToString(CultureInfo.InvariantCulture));
                double ssi = selfscore[i];
                for (int j = i + 1; j < njob; j++)
                {
                    double ssj = selfscore[j];
                    double bunbo = Math.Min(ssi, ssj);
                    if (bunbo == 0.0)
                    {
                        iscore[i][j - i] = 2.0;
                    }
                    else
                    {
                        iscore[i][j - i] = (1.0 - NaivePairScoreFast(seq[i], seq[j], skip[i], skip[j], G.PenaltyDist) / bunbo) * 2.0;
                    }
                    if (iscore[i][j - i] > 10)
                        iscore[i][j - i] = 10.0;
                }
            }
            Report("\rdone.                                           \n");
            return iscore;
        }

        /// <summary>
        /// <c>tbfast _ (pair library arguments) _ -+ n ...</c> (<c>main</c> of <c>tbfast.c</c> with <c>pairlocalalign</c> called inside): the pairwise alignments of all pairs give the guide tree (UPGMA of their distances) and the
        /// constraints (the importance of the aligned blocks) of the progressive alignment of the groups (<c>A__align</c> with <c>constraint</c> 2).
        /// </summary>
        internal Result Tbfast(IReadOnlyList<string> seq0, MafftGlobals? pairG, IReadOnlyList<string>? names = null, NodePairData? nodePair = null, string? hat3Text = null, string? hat2Text = null)
        {
            _isTbfast = true;
            int njob = seq0.Count;
            _njob = njob;
            var seq = new string[njob];
            for (int i = 0; i < njob; i++)
                seq[i] = seq0[i];
            _nlenmax = 0;
            for (int i = 0; i < njob; i++)
                _nlenmax = Math.Max(_nlenmax, seq[i].Length);
            _nlen = new int[njob];
            for (int i = 0; i < njob; i++)
                _nlen[i] = seq[i].Length;
            _nogaplen = new int[njob];
            var fileNames = names ?? DefaultNames(njob);
            char treein = '\0';
            if (G.TreeIn != 0)
            {
                treein = CheckGuideTreeFile().Kind;
                if (treein == 'C')
                {
                    _compactTree = 2;
                    treein = '\0';
                    G.UseFft = false;
                }
                else if (treein == 'n')
                {
                    if (nodePair == null)
                        throw UndefinedBehaviour("the guide tree file says nodepair: the program reads the directory hat3dir that nodepair writes, and it does not exist (segmentation fault, process id in the message)");
                    _compactTree = 3;
                    treein = '\0';
                    G.UseFft = false;
                }
            }
            Report("treein = " + ((int)treein).ToString(CultureInfo.InvariantCulture) + "\ncompacttree = " + _compactTree.ToString(CultureInfo.InvariantCulture) + "\n");
            if (njob < 2)
            {
                // an empty file (the stopped disttbfast of --seed with --merge leaves an empty 'pre'): the program says so and exits with status 0
                Report("Warning: Only " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
                return new Result(Array.Empty<string>(), Array.Empty<int>());
            }
            Report(StackSizeText(njob));
            MafftSubalignments? merge = G.Subalignment ? ReadSubalignmentsTable(njob) : null;
            var topol0 = new int[njob][];
            var topol1 = new int[njob][];
            var len = new double[njob][];
            for (int i = 0; i < njob; i++)
                len[i] = new double[2];
            var dep = new Treedep[njob];
            double[][]? expdist = null;
            if (treein != '\0')
            {
                LoadTree(njob, topol0, topol1, len, dep, G.TreeOut != 0, names ?? DefaultNames(njob));
                Report("\ndone.\n\n");
                if (pairG != null && G.SpecificityConsideration > 0.0)
                {
                    // the distances of the pairs along the user tree (the pair stage needs them for --unalignlevel)
                    var userTree = new GuideTree(topol0, topol1, dep, len);
                    expdist = new double[njob][];
                    for (int i = 0; i < njob; i++)
                        expdist[i] = new double[njob];
                    for (int i = 0; i < njob - 1; i++)
                    {
                        foreach (int a in userTree.Order(i, 0))
                            foreach (int b in userTree.Order(i, 1))
                            {
                                expdist[a][b] += len[i][0] + len[i][1];
                                expdist[b][a] += len[i][0] + len[i][1];
                            }
                    }
                }
            }

            MafftLocalHom[][] localhom = null!;
            double[][] iscore = null!;
            var kozo = new bool[njob];
            int nkozo = 0;
            int nseed = 0;
            MafftLocalHom[][]? seedHom = null;
            if (pairG == null)
            {
                // callpairlocalalign == 0 (the later cycles of the script: tbfast -J): the blocks come from the file hat3 (--seed: the structures of the seed alignments), the distances from the alignment that is read
                if (G.Constraint == 2)
                {
                    Report("Loading 'hat3' ... ");
                    if (hat3Text == null)
                        throw new MafftProgramFailureException("Make hat3.\n");
                    localhom = ReadLocalHomTable2Half(hat3Text, njob, kozo);
                    Report("\ndone.\n");
                    nkozo = kozo.Count(x => x);
                    Kozo = nkozo > 0 ? kozo : null;
                }
                else if (G.DistOut != 0 || G.OutputHat23 != 0)
                {
                    // constraint == 0 (--weighti 0, --retree 0): the program writes the distances it holds before it loads the file (none yet: all zeros) to hat2
                    Report("\nwriting hat2 (1)\n");
                    var zero = new double[njob][];
                    for (int i = 0; i < njob; i++)
                        zero[i] = new double[njob - i];
                    DistFile = WriteHat2(zero, fileNames);
                }
            }
            else if (_compactTree == 3)
            {
                // the pair library is the blocks that nodepair wrote for every node, the guide tree is nodepair's
                Report("Checking 'hat3.seed' ... ");
                if (!string.IsNullOrEmpty(hat3Text) && (char.IsAsciiLetterOrDigit(hat3Text[0]) || hat3Text[0] == ' '))
                {
                    // the sequences of the seed alignments are the first ones of the input (names that start with "_see" after the first character)
                    nseed = 0;
                    while (nseed < njob && fileNames[nseed].Length >= 5 && string.CompareOrdinal(fileNames[nseed], 1, "_seed", 0, 4) == 0)
                        nseed++;
                    seedHom = ReadLocalHomTable2Half(hat3Text.Substring(1), nseed, kozo, null, afterFirstCharacter: true); // fgetc( prep ) has consumed the first character
                    nkozo = 0;
                    for (int i = 0; i < njob; i++)
                    {
                        Report("kozoarivec[" + i.ToString(CultureInfo.InvariantCulture) + "] = " + (kozo[i] ? "1" : "0") + "\n");
                        if (kozo[i])
                            nkozo++;
                    }
                    if (nkozo != nseed)
                        throw new MafftProgramFailureException("problem in input file?  nkozo != nseed\n");
                }
                if (G.Constraint != 2 && G.Constraint != 0)
                    throw new NotSupportedException("The progressive method of tbfast without a pair library is not available in this port.");
            }
            else
            {
                if (G.Constraint != 2 && G.Constraint != 0)
                    throw new NotSupportedException("The progressive method of tbfast without a pair library is not available in this port.");
                pairG.Dorp = G.Dorp;
                pairG.Common = G.Common; // one process: tbfast calls pairlocalalign itself, the trace-back matrices are shared
                var pairRun = new MafftRun(pairG) { _njob = njob, _nlenmax = _nlenmax };
                (double[][] Iscore, MafftLocalHom[][] LocalHom) pairResult;
                try
                {
                    pairResult = pairRun.PairLocalAlign(seq, expdist);
                }
                catch (MafftProgramFailureException)
                {
                    Report(pairRun.Log.ToString());
                    throw;
                }
                (iscore, localhom) = (pairResult.Iscore, pairResult.LocalHom);
                Report(pairRun.Log.ToString());
                Report("outputhat23=" + G.OutputHat23.ToString(CultureInfo.InvariantCulture) + "\n");
                if (G.Constraint == 2)
                {
                    for (int i = 0; i < njob; i++)
                    {
                        for (int j = 0; j < njob - i; j++)
                        {
                            for (var t = localhom[i][j]; t != null; t = t.Next)
                            {
                                if (t.Opt == -1.0)
                                    continue;
                                t.Opt = t.Opt / 5.8 * 600;
                            }
                        }
                    }
                    Report("Loading 'hat3.seed' ... ");
                    if (!string.IsNullOrEmpty(hat3Text))
                        ReadLocalHomTable2Half(hat3Text, njob, kozo, localhom); // --seed: the blocks of the seeds are added to the chains of the pair stage
                    Report("\ndone.\n");
                    nkozo = kozo.Count(x => x);
                    Kozo = nkozo > 0 ? kozo : null;
                    if (G.OutputHat23 != 0)
                    {
                        // the files hat3 and hat2 that dvtditr reads: the values go through their text
                        Report("Writing hat3 for iterative refinement\n");
                        Hat3Table = RoundTripHat3(localhom, njob);
                        Hat2Table = RoundTripHat2(iscore, njob);
                        DistFile = WriteHat2(iscore, fileNames);
                    }
                    else if (G.DistOut != 0)
                    {
                        DistFile = WriteHat2(iscore, fileNames);
                    }
                }
                else
                {
                    // constraint == 0 (--weighti 0, --retree 0): no blocks are kept; the pair stage only gives the distances
                    Report("blosum " + G.Nblosum.ToString(CultureInfo.InvariantCulture) + " / kimura 200\nscoremtx=" + G.Scoremtx.ToString(CultureInfo.InvariantCulture) + "\nfastathreshold=" + G.FastaThreshold.ToString("F6", CultureInfo.InvariantCulture) + "\n");
                    if (G.DistOut != 0 || G.OutputHat23 != 0)
                    {
                        Report("\nwriting hat2 (1)\n");
                        DistFile = WriteHat2(iscore, fileNames);
                        if (G.OutputHat23 != 0)
                            Hat2Table = RoundTripHat2(iscore, njob);
                    }
                    localhom = null!;
                }
            }
            G.Constants(seq);
            Report("Gap Penalty = " + MafftGlobals.F((double)G.Ppenalty / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.PpenaltyEx / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.Poffset / 1000, 5, 2, true) + "\n");
            if (G.DistOut != 0 && G.TreeOut == 0 && G.NoAlign != 0)
                return NoAlignResult(seq, njob, "tbfast");
            if (treein == '\0' && _compactTree != 3)
            {
                Report("tbutree = " + (G.TbUtree ? "1" : "0") + ", compacttree = " + _compactTree.ToString(CultureInfo.InvariantCulture) + "\n");
                if (!G.TbUtree)
                    iscore = TbfastMsaDistances(seq, njob);
                else if (pairG == null && hat2Text != null)
                {
                    // callpairlocalalign == 0 with the distances of the files (dndfast7: --fastapair / --fastswpair): readhat2_doublehalf_pointer, and --distout writes the file again
                    Report("Loading 'hat2' ... ");
                    var full = ReadHat2Text(hat2Text, njob);
                    iscore = new double[njob][];
                    for (int i = 0; i < njob; i++)
                    {
                        iscore[i] = new double[njob - i];
                        for (int j = i + 1; j < njob; j++)
                            iscore[i][j - i] = full[i][j];
                    }

                    Report("done.\n");
                    if (G.DistOut != 0)
                    {
                        Report("\nwriting hat2 (2)\n");
                        DistFile = WriteHat2(iscore, fileNames);
                    }
                }

                double[][]? iscoreKozo = null;
                if (nkozo > 0)
                {
                    // the distances of the sequences that have a structure
                    iscoreKozo = new double[nkozo][];
                    for (int i = 0; i < nkozo; i++)
                        iscoreKozo[i] = new double[nkozo - i];
                    int ik = 0;
                    for (int i = 0; i < njob - 1; i++)
                    {
                        int jk = ik + 1;
                        for (int j = i + 1; j < njob; j++)
                        {
                            if (kozo[i] && kozo[j])
                                iscoreKozo[ik][jk - ik] = iscore[i][j - i];
                            if (kozo[j])
                                jk++;
                        }
                        if (kozo[i])
                            ik++;
                    }
                }
                Report("Constructing a UPGMA tree ... ");
                if (merge != null)
                    ConstrainedSupg(njob, iscore, topol0, topol1, len, dep, merge, fileNames); // --merge: the UPGMA tree in which every subalignment is a group (the branch of the program that precedes compacttree == 3)
                else if (G.TreeOut != 0)
                    MuscleSupg(njob, iscore, topol0, topol1, len, dep, false, fileNames);
                else
                    MuscleSupg(njob, iscore, topol0, topol1, len, dep);
                if (nkozo > 0)
                {
                    // fixed_musclesupg_double_realloc_nobk_halfmtx( nkozo ... ): the tree of the structures is built (its progress is written) and not used
                    var t0 = new int[nkozo][];
                    var t1 = new int[nkozo][];
                    var lk = new double[nkozo][];
                    for (int i = 0; i < nkozo; i++)
                        lk[i] = new double[2];
                    MuscleSupg(nkozo, iscoreKozo!, t0, t1, lk, new Treedep[nkozo]);
                }
                Report("\ndone.\n\n");
            }
            else if (_compactTree == 3)
            {
                // treein_bin: the tree that nodepair wrote
                Report("tbutree = 1, compacttree = 3\n");
                double[][]? iscoreKozo3 = null;
                if (nkozo > 0)
                {
                    // the distances between the seeds: the dynamic program without the alignment on the rows of the seed alignments (gaps included); the program prints every pair
                    var selfscore = new int[njob];
                    for (int i = 0; i < njob; i++)
                        selfscore[i] = kozo[i] ? (int)NaivePairScore11(seq[i], seq[i], 0) : -1;
                    iscoreKozo3 = new double[nkozo][];
                    for (int i = 0; i < nkozo; i++)
                        iscoreKozo3[i] = new double[nkozo - i];
                    int ik = 0;
                    for (int i = 0; i < njob - 1; i++)
                    {
                        int jk = ik + 1;
                        for (int j = i + 1; j < njob; j++)
                        {
                            if (kozo[i] && kozo[j])
                            {
                                Report("seq0=" + seq[i] + "\n");
                                Report("seq1=" + seq[j] + "\n");
                                Report("selfscore0=" + selfscore[0].ToString(CultureInfo.InvariantCulture) + "\n");
                                Report("selfscore1=" + selfscore[1].ToString(CultureInfo.InvariantCulture) + "\n");
                                iscoreKozo3[ik][jk - ik] = Score2Dist(G__align11Noalign(G.NDisConsweightMulti, seq[i], seq[j]), selfscore[i], selfscore[j]);
                                Report("iscore_kozo[" + ik.ToString(CultureInfo.InvariantCulture) + "][" + jk.ToString(CultureInfo.InvariantCulture) + "]=" + MafftGlobals.F(iscoreKozo3[ik][jk - ik], 0, 6) + "\n");
                            }
                            if (kozo[j])
                                jk++;
                        }
                        if (kozo[i])
                            ik++;
                    }
                }
                for (int i = 0; i < njob - 1; i++)
                {
                    topol0[i] = nodePair!.Topol0[i];
                    topol1[i] = nodePair.Topol1[i];
                    len[i][0] = nodePair.Len[i][0];
                    len[i][1] = nodePair.Len[i][1];
                    dep[i] = nodePair.Dep[i];
                }
                if (nkozo > 0)
                {
                    var t0 = new int[nkozo][];
                    var t1 = new int[nkozo][];
                    var lk = new double[nkozo][];
                    for (int i = 0; i < nkozo; i++)
                        lk[i] = new double[2];
                    MuscleSupg(nkozo, iscoreKozo3!, t0, t1, lk, new Treedep[nkozo]);
                }
                Report("\ndone.\n\n");
            }
            var tree = new GuideTree(topol0, topol1, dep, len);
            int[] order = tree.Order(njob - 2, 2);
            if (G.TreeOut != 0 && G.NoAlign != 0)
                return NoAlignResult(seq, njob, "tbfast");
            double[] eff = new double[njob];
            if (G.Tbrweight != 0)
            {
                G.Weight = 3;
                CountEffSimple(njob, tree, eff);
                for (int i = njob - G.Nadd; i < njob; i++)
                    eff[i] /= 100.0;
            }
            else
            {
                for (int i = 0; i < njob; i++)
                    eff[i] = 1.0;
            }
            double[]? effKozo = null;
            if (nkozo > 0)
            {
                effKozo = new double[njob];
                for (int i = 0; i < njob; i++)
                {
                    if (!kozo[i])
                        effKozo[i] = 0.0;
                    else
                        effKozo[i] = G.Tbrweight != 0 ? eff[i] : 1.0;
                }
            }
            var bseq = new string[njob];
            string? originalGaps = null;
            string[]? addbk = null;
            if (G.Nadd > 0 && G.KeepLength != 0)
            {
                originalGaps = RecordOriginalGaps(njob - G.Nadd, seq);
                if (G.MapOut != 0)
                    addbk = seq.Skip(njob - G.Nadd).Select(x => x.Replace("-", "", StringComparison.Ordinal)).ToArray();
            }
            MafftAddState? addState = null;
            if (G.Nadd > 0)
            {
                addState = PrepareAdd(njob, G.Nadd, seq, bseq, tree, fileNames, tbfastMessages: true);
            }
            else if (merge != null)
            {
                addState = PrepareMerge(njob, seq, bseq, topol0, topol1, merge, fileNames, tbfastMessages: true);
            }
            else
            {
                for (int i = 0; i < njob; i++)
                    bseq[i] = seq[i].Replace("-", "");
            }
            int alloclen = _nlenmax * 2 + 1;
            Report("Progressive alignment ... \n");
            TreeBase(bseq, tree, eff, ref alloclen, localhom, tbfast: G.Constraint == 0, add: addState, nodePair: _compactTree == 3 ? nodePair : null, effarrKozo: effKozo, seedHom: seedHom, nseed: nseed);
            Report("\ndone.\n");
            var (deletedRuns, deleteMap, newNames, ndeleted) = KeepLength(njob, seq, bseq, originalGaps, addbk, fileNames);
            if (G.SpScoreOut != 0)
                throw UndefinedBehaviour("tbfast with -S (--scoreout with a pair library, --hybridpair, --seed ...: every run of tbfast): sumofpairsscore runs after freeconstants() and reads a table that the program has released (valgrind: 'Invalid read of size 8', address 0x308 = NULL + 97 * 8; the program ends with a segmentation fault whose message names a process id)");
            Report("tbfast (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + Math.Max(G.Nthread, 1).ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            if (ndeleted > 0)
            {
                Report("\nTo keep the alignment length, " + ndeleted.ToString(CultureInfo.InvariantCulture) + " letters were DELETED.\n");
                Report(G.MapOut != 0 ? "The deleted letters are shown in the (filename).map file.\n" : "To know the positions of deleted letters, rerun the same command with the --mapout option.\n");
            }
            return new Result(bseq, order) { Deleted = deletedRuns, DeleteMap = deleteMap, NewNames = newNames };
        }
    }

    #endregion
}
