using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: setcore (--core, --coreext, --corethr, --corewin)

    internal sealed partial class MafftGlobals
    {
        // setcore.c
        internal double CoreThreshold = 0.5;
        internal bool CoreExtend;

        /// <summary><c>arguments()</c> of <c>setcore.c</c> for the command line that the script writes: <c>-w corewin -i corethr [-c]</c>.</summary>
        internal static MafftGlobals ParseSetcore(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals { Weight = 3, Scoremtx = 0, Nblosum = 62, Ppenalty = NotSpecified, PpenaltyEx = NotSpecified, Poffset = NotSpecified, PpenaltyDist = NotSpecified, KimuraR = NotSpecified, PamN = NotSpecified };
            g.FftWinSize = NotSpecified;
            g.FftThreshold = NotSpecified;
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
                        case 'w': g.FftWinSize = CAtoi(Next()); rest = true; break;
                        case 'i': g.CoreThreshold = CAtof(Next()); g.Report("corethr = " + g.CoreThreshold.ToString("F6", CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'c': g.CoreExtend = true; break;
                        default: break;
                    }
                }
            }
            return g;
        }
    }

    internal sealed partial class MafftRun
    {
        /// <summary>
        /// <c>setcore</c> (<c>setcore.c</c>): reads the alignment <c>pre</c>, the distances of the file <c>hat2</c> (null: there is none, "Make hat2." and the program stops with an empty output), builds the UPGMA tree
        /// (<c>supg</c>), weights the sequences (<c>counteff_simple</c>), and writes the columns between the first and the last window whose average weighted gap content is below the threshold, with <c>corewin</c> gap columns
        /// on both sides (<c>--coreext</c>: with the next <c>corewin</c> residues of every sequence that lie outside the core). The result is the new <c>pre</c>: the names keep their leading '=', every sequence is on one line.
        /// </summary>
        internal (string[] Rows, string[] Names) SetCore(string[] rows, IReadOnlyList<string> names, string? hat2Text)
        {
            int njob = rows.Length;
            if (njob < 2)
            {
                Report("At least 2 sequences should be input!\nOnly " + njob.ToString(CultureInfo.InvariantCulture) + " sequence found.\n");
                return (Array.Empty<string>(), Array.Empty<string>());
            }
            int nlenmax = rows.Max(r => r.Length);
            if (G.FftWinSize < 0)
                throw new NotSupportedException("--corewin with a negative window: the program writes outside its arrays.");
            var seq = (string[])rows.Clone();
            if (G.Dorp == '\0')
                G.Dorp = GuessNucleotideFromRows(rows) ? 'd' : 'p';
            if (G.Dorp == 'd')
                for (int i = 0; i < njob; i++)
                    seq[i] = seq[i].ToLowerInvariant();
            G.Constants(seq);
            Report("Gap Penalty = " + MafftGlobals.F((double)G.Ppenalty / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.PpenaltyEx / 1000, 5, 2, true) + ", " + MafftGlobals.F((double)G.Poffset / 1000, 5, 2, true) + "\n");
            if (CheckAlphabet(seq, "Illeagal") != null)
                throw new NotSupportedException("--core: the alignment holds a letter that setcore does not know (it guesses the alphabet again from the letters; with --anysymbol and a forced --nuc / --amino the symbols of replaceu differ from the program's): registered in B09.md Leftovers (A2-19c-s5 --core).");
            Report("Loading 'hat2' ... ");
            if (hat2Text == null)
            {
                Report("Make hat2.\n");
                return (Array.Empty<string>(), Array.Empty<string>());
            }
            double[][] pscore;
            try
            {
                pscore = ReadHat2Text(hat2Text, njob);
            }
            catch (MafftProgramFailureException e)
            {
                Report(e.Text);
                return (Array.Empty<string>(), Array.Empty<string>());
            }
            Report("done.\n");
            Report("Constructing dendrogram ... ");
            var topol = new int[njob][][];
            var len = new double[njob][];
            for (int i = 0; i < njob; i++)
            {
                topol[i] = new int[2][];
                len[i] = new double[2];
            }
            Supg(njob, pscore, topol, len);
            Report("done.\n");
            var eff = new double[njob];
            CountEffSimpleFull(njob, topol, len, eff);
            int winsize = G.FftWinSize;
            var gapc = new double[Math.Max(nlenmax, winsize) + 1];
            for (int i = 0; i < nlenmax; i++)
            {
                gapc[i] = 0.0;
                for (int j = 0; j < njob; j++)
                    if (i < seq[j].Length && seq[j][i] == '-')
                        gapc[i] += eff[j];
            }
            double tmpavgap = 0.0;
            int corestart = -1, coreend = -1;
            for (int i = 0; i < winsize; i++)
                tmpavgap += gapc[i];
            for (int i = winsize; i < nlenmax; i++)
            {
                double avgap = tmpavgap / winsize;
                if (avgap < G.CoreThreshold)
                {
                    if (corestart == -1)
                        corestart = i - winsize;
                    coreend = i;
                }
                tmpavgap -= gapc[i - winsize];
                tmpavgap += gapc[i];
            }
            if (corestart == -1 || coreend == -1)
            {
                corestart = 0;
                coreend = nlenmax - 1;
            }
            var outRows = new string[njob];
            for (int i = 0; i < njob; i++)
            {
                string s = seq[i];
                var o = new char[winsize + (coreend - corestart + 1) + winsize];
                Array.Fill(o, '-');
                for (int j = corestart; j <= coreend; j++)
                    o[winsize + j - corestart] = j < s.Length ? s[j] : '\0';
                if (G.CoreExtend)
                {
                    int ot = winsize - 1;
                    int m = winsize;
                    int pt = corestart - 1;
                    while (m != 0 && --pt > 0)
                    {
                        if (s[pt] != '-')
                        {
                            o[ot--] = s[pt];
                            m--;
                        }
                    }
                    ot = winsize + coreend - corestart + 1;
                    pt = coreend;
                    m = winsize;
                    while (m != 0 && ++pt < s.Length)
                    {
                        if (s[pt] != '-')
                        {
                            o[ot++] = s[pt];
                            m--;
                        }
                    }
                }
                outRows[i] = new string(o);
            }
            return (outRows, names.ToArray());
        }

        /// <summary><c>readhat2_pointer</c> (the matrix of the file <c>hat2</c>: fixed fields of 6 characters, a line break after 12 values and after the last of a row).</summary>
        private static double[][] ReadHat2Text(string text, int njob)
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
            string l2 = Line();
            int nseq0 = MafftGlobals.CAtoi(l2.Length > 5 ? l2[..5] : l2);
            if (nseq0 != njob)
                throw new MafftProgramFailureException("hat2 is wrong.\n");
            Line();
            for (int i = 0; i < njob; i++)
                Line();
            var mtx = new double[njob][];
            for (int i = 0; i < njob; i++)
                mtx[i] = new double[njob];
            for (int i = 0; i < njob - 1; i++)
            {
                for (int j = i + 1; j < njob; j++)
                {
                    if (pos < text.Length && text[pos] == '\n')
                        pos++;
                    int take = Math.Min(6, Math.Max(0, text.Length - pos));
                    mtx[i][j] = MafftGlobals.CAtof(text.Substring(pos, take));
                    pos += take;
                }
            }
            return mtx;
        }

        /// <summary><c>supg</c> (<c>mltaln9.c</c>): the UPGMA of the program with the mixture <c>sueff_global</c> of the minimum and the average linkage; every group lists all its members.</summary>
        private void Supg(int nseq, double[][] oeff, int[][][] topol, double[][] len)
        {
            var eff = new double[nseq][];
            var pair = new bool[nseq][];
            var tmplen = new double[nseq];
            for (int i = 0; i < nseq; i++)
            {
                eff[i] = (double[])oeff[i].Clone();
                pair[i] = new bool[nseq];
                pair[i][i] = true;
            }
            double sueff = G.SueffGlobal;
            for (int k = 0; k < nseq - 1; k++)
            {
                double minscore = 9999.0;
                int im = -1, jm = -1;
                for (int i = 0; i < nseq - 1; i++)
                {
                    for (int j = i + 1; j < nseq; j++)
                    {
                        double t = eff[i][j];
                        if (t < minscore)
                        {
                            minscore = t;
                            im = i;
                            jm = j;
                        }
                    }
                }
                topol[k][0] = Enumerable.Range(0, nseq).Where(i => pair[im][i]).ToArray();
                topol[k][1] = Enumerable.Range(0, nseq).Where(i => pair[jm][i]).ToArray();
                len[k][0] = minscore / 2.0 - tmplen[im];
                len[k][1] = minscore / 2.0 - tmplen[jm];
                tmplen[im] = minscore / 2.0;
                for (int i = 0; i < nseq; i++)
                    pair[im][i] |= pair[jm][i];
                for (int i = 0; i < nseq; i++)
                    pair[jm][i] = false;
                for (int i = 0; i < nseq; i++)
                {
                    if (i == im || i == jm)
                        continue;
                    int miniim, maxiim, minijm, maxijm;
                    if (i < im)
                    {
                        miniim = i; maxiim = im; minijm = i; maxijm = jm;
                    }
                    else if (i < jm)
                    {
                        miniim = im; maxiim = i; minijm = i; maxijm = jm;
                    }
                    else
                    {
                        miniim = im; maxiim = i; minijm = jm; maxijm = i;
                    }
                    double eff0 = eff[miniim][maxiim];
                    double eff1 = eff[minijm][maxijm];
                    eff[miniim][maxiim] = Math.Min(eff0, eff1) * (1.0 - sueff) + (eff0 + eff1) * 0.5 * sueff;
                    eff[minijm][maxijm] = 9999.0;
                    eff[im][jm] = 9999.0;
                }
            }
        }
    }

    #endregion
}
