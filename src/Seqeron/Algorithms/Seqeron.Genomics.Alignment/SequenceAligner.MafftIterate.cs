using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: dndpre and dvtditr (the tree-dependent iterative refinement, --maxiterate)

    internal sealed partial class MafftGlobals
    {
        // dndpre.c / dvtditr.c / tditeration.c
        internal int MaxDist = 1;
        internal int Niter = 100;
        internal int Cooling = 1;
        internal int UTree = 1;
        internal int Bunkatsu = 1;
        internal int Kobetsubunkatsu;
        internal int FftKeika;
        internal int ScoreCheck = 1;
        internal int DivWinSize = 20;
        internal int DivThreshold = 65;
        internal double MinimumWeight = 0.0005;
        internal double AutoSubalignment;
        internal double Cut;
        internal int RandomSeed;

        /// <summary><c>arguments()</c> of <c>dndpre.c</c>: the penalties are not given by the script (the program uses the defaults of the scoring model).</summary>
        internal static MafftGlobals ParseDndpre(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified };
            int k = 0;
            string Next() => k < argv.Count ? argv[k++] : throw new ArgumentException("A MAFFT option is missing its argument.");
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false;
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    switch (a[ci])
                    {
                        case 'Z': g.UseNaiveScore = true; break;
                        case 't': break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'P': g.Dorp = 'p'; break;
                        case ':': g.NWildcard = 1; break;
                        case 'K': break;
                        case 'I': g.Nadd = CAtoi(Next()); g.Report("nadd = " + g.Nadd.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; rest = true; break;
                        case 'M': g.MaxDist = CAtoi(Next()); rest = true; break;
                        case 'C': g.Nthread = CAtoi(Next()); rest = true; break;
                        default: break;
                    }
                }
            }
            return g;
        }

        /// <summary><c>arguments()</c> of <c>dvtditr.c</c> (the options the script can produce for the iterative refinement without a pair library).</summary>
        internal static MafftGlobals ParseDvtditr(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Nthread = 1, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified, Weight = 4, Alg = 'A', Kobetsubunkatsu = 1, FftKeika = 1 };
            g.UseFft = false;
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
                        case 'H': g.Subalignment = true; g.SubalignmentOffset = CAtoi(Next()); rest = true; break;
                        case 'K': g.Nadd = CAtoi(Next()); g.Report("nadd = " + g.Niter.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'I': g.Niter = CAtoi(Next()); g.Report("niter = " + g.Niter.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q':
                            g.PenaltyShiftFactor = CAtof(Next());
                            if (g.PenaltyShiftFactor < 100.0 && g.PenaltyShiftFactor != 2.0)
                            {
                                g.Report(g.PenaltyShiftFactor.ToString("F6", CultureInfo.InvariantCulture) + ", penalty_shift is fixed to penalty x 2 in the iterative refinement phase.\n");
                                g.PenaltyShiftFactor = 2.0;
                            }
                            rest = true;
                            break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); g.Report("poffset = " + g.Poffset.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); g.Report("kappa = " + g.KimuraR.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; g.Report("blosum " + g.Nblosum.ToString(CultureInfo.InvariantCulture) + " / kimura 200\n"); rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; g.Report("jtt/kimura " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; g.Report("tm " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'C': g.Nthread = CAtoi(Next()); g.Report("nthread = " + g.Nthread.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 't': g.RandomSeed = CAtoi(Next()); g.Report("randomseed = " + g.RandomSeed.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'p': Next(); rest = true; break;
                        case 's': g.SpecificityConsideration = CAtof(Next()); rest = true; break;
                        case 'S': g.SpScoreOut = 1; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'N': g.NeverMemSave = 1; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'l': g.Bunkatsu = 0; g.FastaThreshold = CAtof(Next()); g.Constraint = 2; rest = true; break;
                        case 'T': g.Kobetsubunkatsu = 0; break;
                        case 'B': g.Bunkatsu = 0; break;
                        case 'A': g.Alg = 'A'; break;
                        case 'M': g.Alg = 'M'; break;
                        case '@': g.Alg = 'd'; break;
                        case 'F': g.UseFft = true; break;
                        case 'u': g.Weight = 0; break;
                        case 'U': g.TreeIn = 1; break;
                        case 'J': g.UTree = 0; break;
                        case 'Z': g.ScoreCheck = 0; break;
                        case 'Y': g.ScoreCheck = 2; break;
                        case 'L': g.LegacyGapCost = 1; break;
                        case 'n': g.OutNumber = 1; break;
                        case 'X': g.TreeMethod = 'X'; g.SueffGlobal = CAtof(Next()); g.Report("sueff_global = " + g.SueffGlobal.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'E': g.AutoSubalignment = CAtof(Next()); g.Report("autosubalignment = " + g.AutoSubalignment.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'W': g.MinimumWeight = CAtof(Next()); g.Report("minimumweight = " + g.MinimumWeight.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'z': g.FftThreshold = CAtoi(Next()); rest = true; break;
                        case 'w': g.FftWinSize = CAtoi(Next()); rest = true; break;
                        case ':': g.NWildcard = 1; break;
                        case 'O': break;
                        default: throw new NotSupportedException($"The MAFFT option -{c} of dvtditr is not available in this port.");
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
        // -------------------------------------------------------------------- dndpre.c

        /// <summary>
        /// <c>dndpre -M 2</c> (<c>dndpre.c</c>): the distance of every pair of rows of the alignment <paramref name="rows"/> (<c>2 (1 - score / min(self scores))</c> with the gap penalty <c>penalty</c> of the default scoring,
        /// at most 9.9) written to the file <c>hat2</c>, which <c>dvtditr</c> reads: the values come back with three decimals. Returns the distances as <c>dvtditr</c> reads them (upper triangle, i &lt; j).
        /// </summary>
        internal double[][] Dndpre(string[] rows, IReadOnlyList<string>? names = null)
        {
            var g2 = G;
            int njobAll = rows.Length;
            if (g2.Nadd > 0 && njobAll < 2)
                throw new MafftProgramFailureException("At least 2 sequences should be input!\nOnly " + njobAll.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
            if (g2.Dorp == '\0')
                g2.Dorp = GuessNucleotideFromRows(rows) ? 'd' : 'p';
            int njob = njobAll - g2.Nadd; // the new sequences are ignored
            var seq = new string[njob];
            for (int i = 0; i < njob; i++)
                seq[i] = rows[i];
            for (int i = 1; i < njob; i++)
                if (seq[i].Length != seq[0].Length)
                    throw new MafftProgramFailureException("Not aligned!\n");
            if (g2.Dorp == 'd')
                for (int i = 0; i < njob; i++)
                    seq[i] = seq[i].ToLowerInvariant();
            g2.Constants(seq);
            if (new MafftRun(g2).CheckAlphabet(seq) is { } unknownLetter)
                throw unknownLetter;
            var helper = new MafftRun(g2);
            var selfscore = new double[njob];
            for (int i = 0; i < njob; i++)
                selfscore[i] = helper.NaivePairScore11(seq[i], seq[i], g2.Penalty);
            var skip = new int[njob][];
            for (int i = 0; i < njob; i++)
                skip[i] = MakeSkipTable(seq[i]);
            var mtx = new double[njob][];
            for (int i = 0; i < njob; i++)
                mtx[i] = new double[njob];
            int ilim = njob - 1;
            for (int i = 0; i < ilim; i++)
            {
                double ssi = selfscore[i];
                if (g2.Nthread > 0)
                {
                    // the threads take the pairs in turn (the first one does the rows, every other one that arrives after the last pair prints the progress of the last row once more)
                    if (njob == 1)
                        break;
                }
                else
                {
                    Log.Append((i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "/" + njob.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "\r");
                }
                for (int j = i + 1; j < njob; j++)
                {
                    double ssj = selfscore[j];
                    double bunbo = Math.Min(ssi, ssj);
                    double mtxv;
                    if (bunbo == 0.0)
                        mtxv = g2.MaxDist;
                    else
                        mtxv = g2.MaxDist * (1.0 - helper.NaivePairScoreFast(seq[i], seq[j], skip[i], skip[j], g2.UseNaiveScore ? 0 : g2.Penalty) / bunbo);
                    if (mtxv < 0.0)
                    {
                        Log.Append("WARNING: negative distance, mtxv = " + mtxv.ToString("F6", CultureInfo.InvariantCulture) + "\n");
                        mtxv = 0.0;
                    }
                    if (mtxv > 9.0)
                    {
                        Log.Append("WARNING: Distance " + i + "-" + j + " is strange, " + mtxv.ToString("F6", CultureInfo.InvariantCulture) + ".\n");
                        mtxv = 9.9;
                    }
                    mtx[i][j] = mtxv;
                }

                if (g2.Nthread > 0)
                    Log.Append(D4(i + 1) + "/" + D4(njob) + " (thread " + D4(0) + "), dndpre\r");
            }

            for (int t = 1; t < g2.Nthread && njob > 1; t++)
                Log.Append(D4(njob - 1) + "/" + D4(njob) + " (thread " + D4(t) + "), dndpre\r");
            Log.Append("dndpre (" + g2.KindName + ") Version 7.505\nalg=X, model=" + g2.ModelName + ", amax=" + g2.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + g2.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            // WriteHat2_pointer writes every value as "%#6.3f" (six characters), readhat2_pointer reads six characters back with atof
            var hat2 = new double[njob][];
            for (int i = 0; i < njob; i++)
                hat2[i] = new double[njob];
            for (int i = 0; i < njob - 1; i++)
                for (int j = i + 1; j < njob; j++)
                    hat2[i][j] = MafftGlobals.CAtof(mtx[i][j].ToString("F3", CultureInfo.InvariantCulture).PadLeft(6));
            // the file hat2 itself (WriteHat2_pointer): what --distout copies when the iterative refinement ran
            var half = new double[njob][];
            for (int i = 0; i < njob; i++)
            {
                half[i] = new double[njob - i];
                for (int j = i; j < njob; j++)
                    half[i][j - i] = mtx[i][j];
            }
            DistFile = WriteHat2(half, names?.Take(njob).ToArray() ?? DefaultNames(njob));
            return hat2;
        }

        // -------------------------------------------------------------------- dvtditr.c / tditeration.c

        /// <summary><c>intergroup_score</c>: the sum-of-pairs score between the rows of two groups (a run of gaps costs <c>penalty</c> once), weighted by <paramref name="eff1"/> x <paramref name="eff2"/>.</summary>
        private double IntergroupScore(string[] seq1, string[] seq2, double[] eff1, double[] eff2, int clus1, int clus2, int len)
        {
            int len2 = len - 2;
            double value = 0.0;
            for (int i = 0; i < clus1; i++)
            {
                for (int j = 0; j < clus2; j++)
                {
                    double efficient = eff1[i] * eff2[j];
                    string mseq1 = seq1[i];
                    string mseq2 = seq2[j];
                    double tmpscore = 0.0;
                    int k = 0;
                    while (k < len)
                    {
                        char ms1 = mseq1[k];
                        char ms2 = mseq2[k];
                        if (ms1 == '-' && ms2 == '-')
                        {
                            k++;
                            continue;
                        }
                        tmpscore += G.AminoDisConsweightMulti[ms1][ms2];
                        if (ms1 == '-')
                        {
                            tmpscore += G.Penalty;
                            tmpscore += G.AminoDisConsweightMulti[ms1][ms2];
                            k++;
                            while ((ms1 = CharAt(mseq1, k)) == '-')
                            {
                                tmpscore += G.AminoDisConsweightMulti[ms1][ms2];
                                k++;
                            }
                            k--;
                            if (k > len2)
                                break;
                        }
                        else if (ms2 == '-')
                        {
                            tmpscore += G.Penalty;
                            tmpscore += G.AminoDisConsweightMulti[ms1][ms2];
                            k++;
                            while ((ms2 = CharAt(mseq2, k)) == '-')
                            {
                                tmpscore += G.AminoDisConsweightMulti[ms1][ms2];
                                k++;
                            }
                            k--;
                            if (k > len2)
                                break;
                        }
                        k++;
                    }
                    value += tmpscore * efficient;
                }
            }
            return value;
        }

        private static char CharAt(string s, int i) => i < s.Length ? s[i] : '\0';

        /// <summary><c>searchAnchors</c> (<c>mltaln9.c</c>): the windows of the alignment whose average column score exceeds the threshold (the alignment is cut at their centres; <c>dvtditr</c> refines the pieces one by one).</summary>
        private List<int> SearchAnchors(string[] seq)
        {
            int nseq = seq.Length;
            int len = seq[0].Length;
            int divWinSize = G.DivWinSize;
            double threshold = (double)G.DivThreshold / 100.0 * 600.0 * divWinSize;
            var stra = new double[len];
            for (int i = 0; i < len; i++)
            {
                stra[i] = 0.0;
                int kcyc = nseq - 1;
                for (int k = 0; k < kcyc; k++)
                    for (int j = k + 1; j < nseq; j++)
                        stra[i] += G.NDis[G.AminoN[seq[k][i]]][G.AminoN[seq[j][i]]];
                stra[i] /= (double)nseq * (nseq - 1) / 2;
            }
            var centers = new List<int>();
            int status = 0, length = 0, start = 0;
            double cumscore = 0.0, score = 0.0;
            for (int j = 0; j < Math.Min(divWinSize, len); j++) // the program reads past the end for a shorter alignment; the sum is not used then (the loop below does not run)
                score += stra[j];
            int ii;
            for (ii = 1; ii < len - divWinSize; ii++)
            {
                score = score - stra[ii - 1] + stra[ii + divWinSize - 1];
                if (score > threshold)
                {
                    if (status == 0)
                    {
                        status = 1;
                        start = ii;
                        length = 0;
                        cumscore = 0.0;
                    }
                    length++;
                    cumscore += score;
                }
                if ((score <= threshold || length > 150) && status != 0)
                {
                    centers.Add((start + ii + divWinSize) / 2);
                    length = 0;
                    cumscore = 0.0;
                    status = 0;
                }
            }
            _ = cumscore;
            if (status != 0)
                centers.Add((start + ii + divWinSize) / 2);
            return centers;
        }

        /// <summary>
        /// <c>TreeDependentIteration</c> (<c>tditeration.c</c>, one thread, no pair library): for every branch of the guide tree (the two groups it separates, the sequences weighted by the branch weights of the tree when <c>weight == 4</c>)
        /// the two groups are re-aligned (<see cref="Falign"/> with the groups' own columns as anchors), and the new alignment is kept when its score is not lower than the old one; a pass over all branches is an iteration (forward
        /// and backward alternately) and the refinement ends when <c>2 n</c> branches in a row do not change or the score oscillates. <paramref name="aseq"/> is the alignment of the segment, the result is returned.
        /// </summary>
        private int _maxDistClass = -1;

        private string[] TreeDependentIteration(int locnjob, string[] segment, int[][][] topol, double[][] len, double[][] distmtx, int[][] skipthisbranch, int alloclen, MafftLocalHom[][]? localhom = null, bool[]? kozo = null)
        {
            bool anyKozo = kozo != null && kozo.Any(x => x);
            var effarrKozo = new double[locnjob];
            var aseq = (string[])segment.Clone();
            var bseq = (string[])segment.Clone();
            // --unalignlevel: one scoring matrix per distance class of the pair of groups (makescoringmatrices); the distance of a pair is the sum of the distances of its two sequences from the branch (<c>USEDISTONTREE</c>)
            VariousDist? various = null;
            double[][]? smalldistmtx = null;
            if (G.SpecificityConsideration != 0.0)
            {
                various = new VariousDist { MaxDistClass = _maxDistClass, Matrices = new double[_maxDistClass][][], Eff1s = new double[_maxDistClass][], Eff2s = new double[_maxDistClass][] };
                for (int c = 0; c < _maxDistClass; c++)
                {
                    double rep = 2.0 * c / G.NDistClass;
                    various.Matrices[c] = MakeDynamicMtx(G.NDisConsweightMulti, rep * 0.5);
                    various.Eff1s[c] = new double[locnjob];
                    various.Eff2s[c] = new double[locnjob];
                }
                smalldistmtx = new double[locnjob][];
                for (int i = 0; i < locnjob; i++)
                    smalldistmtx[i] = new double[locnjob];
            }
            _ = distmtx;
            var effarr = new double[locnjob];
            var distarr = new double[locnjob];
            var branchWeight = new double[locnjob][];
            for (int i = 0; i < locnjob; i++)
                branchWeight[i] = new double[2];
            if (localhom != null && G.UTree != 0)
            {
                // constraint (before treeCnv, which raises short branch lengths in place): the importance of the blocks of the pair library from the weights of the tree (counteff_simple of the full member lists)
                var effarrforlocalhom = new double[locnjob];
                CountEffSimpleFull(locnjob, topol, len, effarrforlocalhom);
                CalcImportanceHalf(locnjob, effarrforlocalhom, aseq, localhom, alloclen);
            }
            MafftUnrootedTree? stopol = null;
            if (G.UTree != 0 && locnjob > 2 && (G.Weight == 4 || G.Weight == 0))
                stopol = new MafftUnrootedTree(locnjob, topol, len, branchWeight);
            var history = new double[G.Niter][][];
            for (int i = 0; i < G.Niter; i++)
            {
                history[i] = new double[locnjob][];
                for (int j = 0; j < locnjob; j++)
                    history[i][j] = new double[2];
            }
            int converged = 0;
            double cut = G.Cut;
            if (G.Cooling != 0)
                cut *= 2.0;
            for (int iterate = 0; iterate < G.Niter; iterate++)
            {
                if (G.Cooling != 0)
                    cut *= 0.5;
                int lin, ldf;
                if (iterate % 2 == 0)
                {
                    lin = 0;
                    ldf = +1;
                }
                else
                {
                    lin = locnjob - 2;
                    ldf = -1;
                }
                for (int l = lin; l < locnjob - 1 && l >= 0; l += ldf)
                {
                    foreach (int k in l == locnjob - 2 ? new[] { 1 } : new[] { 0, 1 })
                    {
                        Report("STEP " + (iterate + 1).ToString("000", CultureInfo.InvariantCulture) + "-" + (l + 1).ToString("000", CultureInfo.InvariantCulture) + "-" + k + " ");
                        if (skipthisbranch[l][k] != 0)
                        {
                            Report(" skip.      \r");
                            continue;
                        }
                        if (locnjob == 2)
                        {
                            // distFromABranch for two sequences: the two branch lengths of the only node
                            distarr[0] = len[0][0];
                            distarr[1] = len[0][1];
                        }
                        else
                        {
                            stopol?.DistFromABranch(distarr, len, l, k);
                        }
                        var (memlist0, memlist1) = OneClusterAndTheOther(locnjob, topol, l, k);
                        int s1 = memlist0[0], s2 = memlist1[0];
                        if (smalldistmtx != null)
                        {
                            for (int i = 0; i < memlist0.Length; i++)
                                for (int j = 0; j < memlist1.Length; j++)
                                {
                                    int im = memlist0[i], jm = memlist1[j];
                                    smalldistmtx[i][j] = distarr[im] + distarr[jm]; // USEDISTONTREE (tddis.c): the distance of the two sequences from the branch
                                }
                        }
                        if (G.Weight == 0)
                        {
                            for (int i = 0; i < locnjob; i++)
                                effarr[i] = 1.0;
                            if (anyKozo)
                                for (int i = 0; i < locnjob; i++)
                                    effarrKozo[i] = kozo![i] ? 1.0 : 0.0;
                        }
                        else if (G.Weight == 4)
                        {
                            stopol!.WeightFromABranch(effarr, l, k);
                            if (anyKozo)
                                for (int i = 0; i < locnjob; i++)
                                    effarrKozo[i] = kozo![i] ? effarr[i] : 0.0; // a single weight
                        }
                        for (int i = 0; i < locnjob; i++)
                            aseq[i] = bseq[i];
                        int length = aseq[0].Length;
                        var mseq1 = new string[memlist0.Length];
                        var mseq2 = new string[memlist1.Length];
                        var effarr1 = new double[memlist0.Length];
                        var effarr2 = new double[memlist1.Length];
                        double[]? effarr1Kozo = null, effarr2Kozo = null;
                        if (anyKozo)
                        {
                            effarr1Kozo = new double[memlist0.Length];
                            effarr2Kozo = new double[memlist1.Length];
                            FastConjunctionKozo(memlist0, aseq, mseq1, effarr1, effarr, effarr1Kozo, effarrKozo);
                            FastConjunctionKozo(memlist1, aseq, mseq2, effarr2, effarr, effarr2Kozo, effarrKozo);
                        }
                        else
                        {
                            FastConjunction(memlist0, aseq, mseq1, effarr1, effarr, G.MinimumWeight);
                            FastConjunction(memlist1, aseq, mseq2, effarr2, effarr, G.MinimumWeight);
                        }
                        int clus1 = memlist0.Length, clus2 = memlist1.Length;
                        if (various != null)
                            ClassifyPairs(various, clus1, effarr1, clus2, effarr2, smalldistmtx!);
                        string[] scoreSeq1 = mseq1, scoreSeq2 = mseq2;
                        if (G.Alg == 'd')
                        {
                            // "Ichijiteki, koredeha scorecheck ga dekinai": the common gaps of both groups are removed in place, the score of the old alignment is taken over the old length from the buffers of the program (the NUL of
                            // the shortened row and what the rows held behind it before)
                            (scoreSeq1, mseq1) = CommonGapPickInPlace(mseq1);
                            (scoreSeq2, mseq2) = CommonGapPickInPlace(mseq2);
                        }
                        double mscore = IntergroupScore(scoreSeq1, scoreSeq2, effarr1, effarr2, clus1, clus2, length); // computed with -Z too ("Oscilation check no tame hitsuyou"): -Z only skips the score of the new alignment (tscore = mscore + 1)
                        double impmatchdouble = 0.0;
                        if (localhom != null && G.Alg == 'd')
                        {
                            // constraint == 2 with D__align: the importance matrix of the groups (imp_match_init_strictD, summed over the diagonal only when the score is checked), the dynamic program with the importance
                            if (G.ScoreCheck != 0)
                                throw new NotSupportedException("The iterative refinement with D__align and the score check is not reachable from the script (--sp passes -Z).");
                            var shrink = ShrinkLocalHom(memlist0, memlist1, localhom);
                            var imp = FillImp(clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, shrink, memlist0, memlist1, false, effarr1Kozo, effarr2Kozo, dmtx: true);
                            DAlign(G.NDisConsweightMulti, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, imp, ref impmatchdouble, 1, 1, various);
                        }
                        else if (localhom != null)
                        {
                            // msshrinklocalhom_fast_half, part_imp_match_init_strict, the importance on the diagonal of the groups' own alignment
                            var shrink = ShrinkLocalHom(memlist0, memlist1, localhom);
                            double oimpmatchdouble = 0.0;
                            PartFillImp(clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, shrink, memlist0, memlist1, effarr1Kozo, effarr2Kozo);
                            for (int i = length - 1; i >= 0; i--)
                                oimpmatchdouble += _partImpMtx![i][i];
                            if (G.ScoreCheck != 0)
                                mscore = oimpmatchdouble + mscore;
                            if (!G.UseFft)
                                throw new InvalidOperationException("Unexpected error.  Please contact katoh@ifrec.osaka-u.ac.jp");
                            FalignLocalHom(mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, out impmatchdouble, various);
                        }
                        else if (G.UseFft)
                            Falign(G.NDisConsweightMulti, mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen);
                        else
                            throw new InvalidOperationException("Unexpected error.  Please contact katoh@ifrec.osaka-u.ac.jp");
                        for (int i = 0; i < clus1; i++)
                            aseq[memlist0[i]] = mseq1[i];
                        for (int i = 0; i < clus2; i++)
                            aseq[memlist1[i]] = mseq2[i];
                        bool identity = aseq[s1] == bseq[s1] && aseq[s2] == bseq[s2];
                        length = mseq1[0].Length;
                        double tscore;
                        if (identity)
                        {
                            tscore = mscore;
                            Report(" identical.   ");
                            converged++;
                        }
                        else
                        {
                            if (G.ScoreCheck != 0)
                            {
                                tscore = IntergroupScore(mseq1, mseq2, effarr1, effarr2, clus1, clus2, length);
                                if (localhom != null)
                                    tscore = impmatchdouble + tscore;
                            }
                            else
                            {
                                tscore = mscore + 1.0;
                            }
                            if (double.IsNaN(mscore))
                                throw new MafftProgramFailureException("\n\nmscore became NaN\n");
                            if (double.IsNaN(tscore))
                                throw new MafftProgramFailureException("\n\ntscore became NaN\n");
                            if (tscore > mscore - cut / 100.0 * mscore)
                            {
                                for (int i = 0; i < locnjob; i++)
                                    bseq[i] = aseq[i];
                                Report(" accepted.");
                                converged = 0;
                            }
                            else
                            {
                                Report(" rejected.");
                                tscore = mscore;
                                converged++;
                            }
                        }
                        Report(G.Alg == 'd' ? "\n" : "\r");
                        history[iterate][l][k] = tscore;
                        if (converged >= locnjob * 2)
                        {
                            Report("\nConverged.\n\n");
                            return bseq;
                        }
                        if (iterate >= 1)
                        {
                            bool oscillating = false;
                            for (int ii = iterate - 2; ii >= 0; ii -= 2)
                            {
                                if (tscore == history[ii][l][k])
                                {
                                    oscillating = true;
                                    break;
                                }
                            }
                            if ((oscillating && G.Cooling == 0) || (oscillating && cut < 0.001 && G.Cooling != 0))
                            {
                                Report("\nOscillating.\n\n");
                                return bseq;
                            }
                        }
                    }
                }
            }
            return bseq;
        }

        /// <summary><c>classifypairs</c>: the class of every pair of the two groups from its distance (<c>(int)( dist / 2 * ndistclass )</c>, at most <c>maxdistclass - 1</c>); the weights of a sequence are kept for the classes in which it has a pair.</summary>
        private static void ClassifyPairs(VariousDist vd, int n1, double[] eff1, int n2, double[] eff2, double[][] smalldistmtx)
        {
            int nc = vd.MaxDistClass;
            for (int c = 0; c < nc; c++)
            {
                Array.Clear(vd.Eff1s[c], 0, n1);
                Array.Clear(vd.Eff2s[c], 0, n2);
            }
            vd.Which = new int[n1][];
            for (int i = 0; i < n1; i++)
            {
                vd.Which[i] = new int[n2];
                for (int j = 0; j < n2; j++)
                {
                    int c = (int)(smalldistmtx[i][j] / 2.0 * 10);
                    if (c >= nc)
                        c = nc - 1;
                    vd.Eff1s[c][i] = eff1[i];
                    vd.Eff2s[c][j] = eff2[j];
                    vd.Which[i][j] = c;
                }
            }
        }

        /// <summary><c>counteff_simple</c> (<c>mltaln9.c</c>): the weight of a sequence from the lengths of the branches between it and the root of the tree, with the complete member lists of the groups.</summary>
        private static void CountEffSimpleFull(int nseq, int[][][] topol, double[][] len, double[] node)
        {
            var rootnode = new double[nseq];
            var eff = new double[nseq];
            for (int i = 0; i < nseq; i++)
                eff[i] = 1.0;
            for (int i = 0; i < nseq - 1; i++)
            {
                foreach (int s1 in topol[i][0])
                {
                    rootnode[s1] += len[i][0] * eff[s1];
                    eff[s1] *= 0.5;
                }
                foreach (int s2 in topol[i][1])
                {
                    rootnode[s2] += len[i][1] * eff[s2];
                    eff[s2] *= 0.5;
                }
            }
            for (int i = 0; i < nseq; i++)
                rootnode[i] += 0.001;
            double total = 0.0;
            for (int i = 0; i < nseq; i++)
                total += rootnode[i];
            for (int i = 0; i < nseq; i++)
                node[i] = rootnode[i] / total;
        }

        private void PartFillImp(int clus1, int clus2, int lgth1, int lgth2, string[] seq1, string[] seq2, double[] eff1, double[] eff2, MafftLocalHom?[][] localhom, int[] orinum1, int[] orinum2, double[]? eff1Kozo = null, double[]? eff2Kozo = null)
            => FillImp(clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, orinum1, orinum2, true, eff1Kozo, eff2Kozo);

        /// <summary><c>OneClusterAndTheOther_fast</c>: the sequences of the group (step, branch) and all the others.</summary>
        private static (int[] Mem1, int[] Mem2) OneClusterAndTheOther(int locnjob, int[][][] topol, int step, int branch)
        {
            var pair = new bool[locnjob];
            var m1 = new List<int>();
            foreach (int r1 in topol[step][branch])
            {
                pair[r1] = true;
                m1.Add(r1);
            }
            var m2 = new List<int>();
            for (int i = 0; i < locnjob; i++)
                if (!pair[i])
                    m2.Add(i);
            return (m1.ToArray(), m2.ToArray());
        }

        /// <summary>
        /// <c>dvtditr</c> (<c>dvtditr.c</c>): the iterative refinement of the alignment <paramref name="rows"/> (the output of <c>disttbfast</c>): the guide tree is the UPGMA of the distances <paramref name="hat2"/> of <c>dndpre</c>, the alignment is cut
        /// into segments at the windows with a high column score (<see cref="SearchAnchors"/>; <c>--bunkatsunashi</c> = <c>-B</c> keeps one), and every segment is refined by <see cref="TreeDependentIteration"/>.
        /// </summary>
        internal Result Dvtditr(string[] rows, double[][] hat2, MafftLocalHom[][]? localhom = null, IReadOnlyList<string>? names = null, bool[]? kozo = null)
        {
            int njob = rows.Length;
            if (G.SpScoreOut != 0)
                throw UndefinedBehaviour("dvtditr with -S (--maxiterate n --scoreout, every algorithm): sumofpairsscore reads a block of 48 bytes that the program has freed (valgrind: invalid read; the program ends with a segmentation fault whose message names a process id)");
            MafftSubalignments? merge = G.Subalignment ? ReadSubalignmentsTable(njob) : null;
            for (int i = 1; i < njob; i++)
                if (rows[i].Length != rows[0].Length)
                    throw new ArgumentException("Input pre-aligned data");
            var seqG = (string[])rows.Clone();
            if (G.Dorp == '\0')
                G.Dorp = GuessNucleotideFromRows(rows) ? 'd' : 'p';
            if (G.Dorp == 'd')
                for (int i = 0; i < njob; i++)
                    seqG[i] = seqG[i].ToLowerInvariant();
            if (localhom != null)
                Report("Loading 'hat3' ... done.\n");
            if (G.SpecificityConsideration != 0.0)
            {
                // calcmaxdistclass: the classes of distance up to the one whose offset is 0
                int cmax = 0;
                for (; cmax < G.NDistClass; cmax++)
                {
                    double rep = 2.0 * cmax / G.NDistClass;
                    if (Math.Min(0.0, rep * 0.5 - G.SpecificityConsideration) == 0.0)
                        break;
                }
                Report("ndistclass = " + G.NDistClass.ToString(CultureInfo.InvariantCulture) + ", maxdistclass = " + (cmax + 1).ToString(CultureInfo.InvariantCulture) + "\n");
                _maxDistClass = cmax + 1;
            }
            G.Constants(seqG);
            if (CheckAlphabet(seqG) is { } unknownLetter)
                throw unknownLetter;
            if (njob == 2)
            {
                G.Weight = 0;
                G.Niter = 1;
            }
            CommonGapPick(seqG);
            var eff = new double[njob][];
            for (int i = 0; i < njob; i++)
            {
                eff[i] = new double[njob - i];
                for (int j = i; j < njob; j++)
                    eff[i][j - i] = j > i ? hat2[i][j] : 0.0;
            }
            var topol0 = new int[njob][];
            var topol1 = new int[njob][];
            var len = new double[njob][];
            for (int i = 0; i < njob; i++)
                len[i] = new double[2];
            if (G.TreeIn != 0)
                VeryFastSupgLoadTree(njob, eff, topol0, topol1, len, names ?? DefaultNames(njob));
            else if (merge != null)
                ConstrainedSupg(njob, eff, topol0, topol1, len, new Treedep[njob], merge, names ?? DefaultNames(njob));
            else if (G.TreeMethod is 'X' or 'E' or 'q')
                MuscleSupg(njob, eff, topol0, topol1, len, new Treedep[njob], fullMembership: true, treeNames: names ?? DefaultNames(njob));
            else
                throw new InvalidOperationException("Incorrect treemethod.");
            List<int[]>? subtable = merge?.Groups;
            if (G.AutoSubalignment > 0.0 && merge == null)
            {
                // --skipiterate x (dvtditr -E x): the groups of the tree whose branches (summed from a tip) pass x are subalignments that the refinement does not change inside
                subtable = new List<int[]>();
                if (GenerateSubalignmentsTable(njob, topol0, topol1, len, G.AutoSubalignment, subtable, out int maxmem))
                {
                    Report("################################################################################################ \n#\n# WARNING: Iterative refinment was not done because you gave a large --fix value ("
                        + MafftGlobals.F(G.AutoSubalignment, 6, 3) + ").\n#\n################################################################################################ \n");
                    Report("dvtditr (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                        + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
                    // the program leaves the file 'order' of the previous program and writes the alignment as it was (without the common gaps)
                    return new Result(seqG, Array.Empty<int>());
                }
                Report("nsubalignments = " + subtable.Count.ToString(CultureInfo.InvariantCulture) + ", maxmem = " + maxmem.ToString(CultureInfo.InvariantCulture) + "\n");
            }
            var order = new int[njob];
            Array.Copy(topol0[njob - 2], 0, order, 0, topol0[njob - 2].Length);
            Array.Copy(topol1[njob - 2], 0, order, topol0[njob - 2].Length, topol1[njob - 2].Length);
            Report("\n");
            var topol = new int[njob][][];
            for (int i = 0; i < njob - 1; i++)
                topol[i] = new[] { topol0[i], topol1[i] };
            var anchors = new List<int> { 0 };
            if (G.Bunkatsu == 0 || localhom != null)
            {
                anchors.Add(seqG[0].Length);
            }
            else
            {
                foreach (int c in SearchAnchors(seqG))
                    anchors.Add(c);
                anchors.Add(seqG[0].Length);
            }
            var skipthisbranch = new int[njob][];
            for (int i = 0; i < njob; i++)
                skipthisbranch[i] = new int[2];
            if (subtable != null)
                CheckSubalignments(njob, topol, subtable, skipthisbranch, merge != null);
            if (G.Nadd > 0)
            {
                // --add with --maxiterate (dvtditr -K n): the branches that cut the original alignment (the first njob - n sequences) are not refined
                var reftable = Enumerable.Range(0, njob - G.Nadd).ToArray();
                for (int j = 0; j < njob - 1; j++)
                {
                    for (int k = 0; k < 2; k++)
                        if (topol[j][k].Any(reftable.Contains) && !SameMember(topol[j][k], reftable))
                            skipthisbranch[j][k] = 1;
                }
            }
            int nseg = anchors.Count;
            var resG = new StringBuilder[njob];
            for (int i = 0; i < njob; i++)
                resG[i] = new StringBuilder();
            int offset = 0;
            int alloclen = seqG[0].Length * 9;
            for (int iseg = 0; iseg < nseg - 1; iseg++)
            {
                int tmplen = anchors[iseg + 1] - anchors[iseg];
                int pos = resG[0].Length;
                var seg = new string[njob];
                for (int j = 0; j < njob; j++)
                    seg[j] = seqG[j].Substring(offset, tmplen);
                offset += tmplen;
                Report("Segment " + (iseg + 1).ToString(CultureInfo.InvariantCulture).PadLeft(3) + "/" + (nseg - 1).ToString(CultureInfo.InvariantCulture).PadLeft(3) + " " + (pos + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "-" + (pos + 1 + tmplen).ToString(CultureInfo.InvariantCulture).PadLeft(4) + "\n");
                var bseq = TreeDependentIteration(njob, seg, topol, len, eff, skipthisbranch, alloclen, localhom, kozo);
                for (int i = 0; i < njob; i++)
                    resG[i].Append(bseq[i]);
            }
            Report("done\n");
            if (G.SpScoreOut != 0)
                ReportSumOfPairs(resG.Select(b => b.ToString()).ToArray());
            Report("dvtditr (" + G.KindName + ") Version 7.505\nalg=" + G.Alg + ", model=" + G.ModelName + ", amax=" + G.SpecificityConsideration.ToString("F1", CultureInfo.InvariantCulture) + "\n"
                + G.Nthread.ToString(CultureInfo.InvariantCulture) + " thread(s)\n\n");
            return new Result(resG.Select(b => b.ToString()).ToArray(), order);
        }

        /// <summary><c>generatesubalignmentstable</c> (mltaln9.c): true when the whole tree lies below the threshold (nothing is refined).</summary>
        private bool GenerateSubalignmentsTable(int nseq, int[][] topol0, int[][] topol1, double[][] len, double threshold, List<int[]> table, out int maxmem)
        {
            Report("\n\n\n");
            maxmem = 0;
            var distfromtip = new double[nseq];
            for (int i = 0; i < nseq - 1; i++)
            {
                int rep0 = topol0[i][0];
                double distfromtip0 = distfromtip[rep0];
                distfromtip[rep0] += len[i][0];
                int rep1 = topol1[i][0];
                double distfromtip1 = distfromtip[rep1];
                distfromtip[rep1] += len[i][1];
                if (topol0[i].Length > 1 && distfromtip0 <= threshold && threshold < distfromtip[rep0])
                {
                    table.Add((int[])topol0[i].Clone());
                    maxmem = Math.Max(maxmem, topol0[i].Length);
                }
                if (topol1[i].Length > 1 && distfromtip1 <= threshold && threshold < distfromtip[rep1])
                {
                    table.Add((int[])topol1[i].Clone());
                    maxmem = Math.Max(maxmem, topol1[i].Length);
                }
            }
            return distfromtip[0] <= threshold;
        }

        /// <summary>The check of dvtditr for every subalignment: the branches that lie inside it are skipped by the refinement, and it must be a group of the tree.</summary>
        private void CheckSubalignments(int njob, int[][][] topol, List<int[]> subtable, int[][] skipthisbranch, bool fromFile = false)
        {
            for (int i = 0; i < subtable.Count; i++)
            {
                Report("\nChecking subalignment " + (i + 1).ToString(CultureInfo.InvariantCulture) + ":\n");
                foreach (int m in subtable[i])
                    Report(" " + (m + 1).ToString(CultureInfo.InvariantCulture) + " ");
                Report("\n");
                for (int j = 0; j < njob - 1; j++)
                {
                    for (int k = 0; k < 2; k++)
                        if (IncludeMember(topol[j][k], subtable[i]) && !SameMember(topol[j][k], subtable[i]))
                            skipthisbranch[j][k] = 1;
                }
                bool found = false;
                for (int j = 0; j < njob - 1 && !found; j++)
                    found = SameMember(topol[j][0], subtable[i]) || SameMember(topol[j][1], subtable[i]);
                if (!found)
                {
                    if (!fromFile)
                        throw new InvalidOperationException("A subalignment made from the guide tree is a group of the tree.");
                    Report(TreeFile != null ? "" : "cp: cannot stat 'infile.tree': No such file or directory\n");
                    int off = G.SubalignmentOffset;
                    string note = off == 0 ? "" : "#\n# You specified seed alignment(s) consisting of " + off.ToString(CultureInfo.InvariantCulture) + " sequences.\n# In this case, the rule of numbering is:\n#   The aligned seed sequences are numbered as 1 .. "
                        + off.ToString(CultureInfo.InvariantCulture) + "\n#   The input sequences to be aligned are numbered as " + (off + 1).ToString(CultureInfo.InvariantCulture) + " .. " + (off + njob).ToString(CultureInfo.InvariantCulture) + "\n";
                    throw new MafftProgramFailureException("\n###############################################################################\n# ERROR!\n# Subalignment " + (i + 1).ToString(CultureInfo.InvariantCulture)
                        + " does not seem to form a monophyletic cluster\n# in the guide tree ('GuideTree' in this directory) internally computed.\n# If you really want to use this subalignment, pelase give a tree with --treein \n# http://mafft.cbrc.jp/alignment/software/treein.html\n"
                        + "# http://mafft.cbrc.jp/alignment/software/merge.html\n" + note + "############################################################################### \n\n");
                }
                Report(" -> OK\n");
            }
        }

        private static bool IncludeMember(int[] mem, int[] cand) => mem.All(cand.Contains);

        private static bool SameMember(int[] mem, int[] cand) => mem.Length == cand.Length && mem.All(cand.Contains);

        /// <summary><c>commongappick</c> of the rows of a group in place (<c>io.c</c>): returns the compacted rows and the rows as the program's buffers hold them afterwards (the shortened row, its NUL, then the old characters behind the new end).</summary>
        private static (string[] Buffers, string[] Rows) CommonGapPickInPlace(string[] seq)
        {
            int len = seq[0].Length;
            var keep = new List<int>();
            for (int i = 0; i < len; i++)
            {
                bool any = false;
                for (int j = 0; j < seq.Length; j++)
                {
                    if (seq[j][i] != '-')
                    {
                        any = true;
                        break;
                    }
                }
                if (any)
                    keep.Add(i);
            }
            var buffers = new string[seq.Length];
            var rows = new string[seq.Length];
            for (int j = 0; j < seq.Length; j++)
            {
                var sb = new StringBuilder(len);
                foreach (int i in keep)
                    sb.Append(seq[j][i]);
                rows[j] = sb.ToString();
                if (keep.Count < len)
                    sb.Append('\0').Append(seq[j], keep.Count + 1, len - keep.Count - 1);
                buffers[j] = sb.ToString();
            }
            return (buffers, rows);
        }

        private static void CommonGapPick(string[] seq)
        {
            int len = seq[0].Length;
            var keep = new List<int>();
            for (int i = 0; i < len; i++)
            {
                bool any = false;
                for (int j = 0; j < seq.Length; j++)
                {
                    if (seq[j][i] != '-')
                    {
                        any = true;
                        break;
                    }
                }
                if (any)
                    keep.Add(i);
            }
            for (int j = 0; j < seq.Length; j++)
            {
                var sb = new StringBuilder(keep.Count);
                foreach (int i in keep)
                    sb.Append(seq[j][i]);
                seq[j] = sb.ToString();
            }
        }

        private static bool GuessNucleotideFromRows(string[] rows)
        {
            int atgc = 0, total = 0;
            foreach (var s in rows)
            {
                if (total >= 1000000)
                    break;
                int n = 0, a = 0;
                foreach (char ch in s)
                {
                    if (!((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')))
                        continue;
                    n++;
                    char c = char.ToLowerInvariant(ch);
                    if (c is 'a' or 't' or 'g' or 'c' or 'u' or 'n')
                        a++;
                }
                atgc += a;
                total += n;
            }
            return (double)atgc / total > 0.75;
        }
    }

    #endregion
}
