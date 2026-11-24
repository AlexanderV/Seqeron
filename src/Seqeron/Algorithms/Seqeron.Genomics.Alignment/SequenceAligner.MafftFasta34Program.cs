#pragma warning disable S3776, S107, S1541, S138, S1854, S1121, S127, S3358, S1066, S4144, S5443, S1125
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region FASTA 3.4 (fasta34): the whole program as dndfast7 calls it

    internal static partial class MafftFasta34
    {
        /// <summary>What a run of <c>fasta34</c> writes: standard output (with the process-dependent parts in the canonical form of the oracle: <c>/tmp/query-P</c>, <c>/tmp/data-P</c>, times <c>T</c>), standard error and the exit status.</summary>
        internal sealed record ProgramRun(string Text, string Err, int ExitCode);

        private const string QueryPath = "/tmp/query-P", LibraryPath = "/tmp/data-P";

        /// <summary><c>nrv[]</c> of <c>scaleswn.c</c>: the random offsets of <c>addhistz</c> (the array has 99 initialisers: the program lacks a comma after the 96th, "0.0219427302 -0.1774469919" is one value, and the last element is 0).</summary>
        private static readonly double[] Nrv =
        {
            0.3098900570, -0.0313400923, 0.1131975903, -0.2832547606, 0.0073672659,
            0.2914489107, 0.4209306311, -0.4630181404, 0.3326537896, 0.0050140359,
            -0.1117435426, -0.2835630301, 0.2302997065, -0.3102716394, 0.0819894916,
            -0.1676455701, -0.3782225018, -0.3204509938, -0.3594969187, -0.0308950398,
            0.2922813812, 0.1337170751, 0.4666577031, -0.2917784349, -0.2438179916,
            0.3002301394, 0.0231147123, 0.5687927366, -0.2318208709, -0.1476839273,
            -0.0385043851, -0.1213476523, 0.1486341995, 0.1027917167, 0.1409192644,
            -0.3280652579, 0.4232041455, 0.0775993309, 0.1159071787, 0.2769424442,
            0.3197284751, 0.1507346903, 0.0028580909, 0.4825103412, -0.0496843610,
            -0.2754357656, 0.6021881753, -0.0816123956, -0.0899148991, 0.4847183201,
            0.2151621865, -0.4542246220, 0.0690709102, 0.2461894193, 0.2126042295,
            -0.0767060668, 0.4819746149, 0.3323031326, 0.0177600676, 0.1143185210,
            0.2653977455, 0.0921872958, -0.1330986718, 0.0412287716, -0.1691604748,
            -0.0529679078, -0.0194157955, -0.6117493924, 0.1199067932, 0.0210243193,
            -0.5832259838, -0.1685528664, 0.0008591271, -0.1120347822, 0.0839125069,
            -0.2787486831, -0.1937017962, -0.1915733940, -0.7888453635, -0.3316745163,
            0.1180885226, -0.3347001067, -0.2477492636, -0.2445697600, 0.0001342482,
            -0.0015759812, -0.1516473992, -0.5202267615, 0.2136975210, 0.2500423188,
            -0.2402926401, -0.1094186280, -0.0618869933, -0.0815221188, 0.2623337275,
            0.0219427302 - 0.1774469919, 0.0828245026, -0.3271952808, -0.0632898028, 0.0,
        };

        private static Molecule WithKtup(Molecule m, int ktup)
        {
            var c = m.Clone();
            c.Ktup = ktup;
            c.OptWid = m.Nucleotide ? 16 : ktup == 1 ? 32 : 16;
            return c;
        }

        /// <summary><c>zs_to_Ec</c> (<c>scaleswn.c</c>).</summary>
        private static double ZsToEc(double zs, long entries)
        {
            if (entries < 5)
                return 0.0;
            double z = (zs - 50.0) / 10.0;
            if (z > ZsMax)
                return 1.0;
            double e = Math.Exp(-PiSqrt6 * z - EulerG);
            return entries * (e > .01 ? Math.Exp(-e) : 1.0 - e);
        }

        private static string StatInfo(Molecule m, int n0) =>
            "Altschul/Gish params: n0: " + I(n0) + " Lambda: " + F(m.Lambda, 3).PadLeft(5) + " K: " + F(m.K, 3).PadLeft(5) + " H: " + F(m.H, 3).PadLeft(5);

        /// <summary><c>prhist</c> (<c>compacc.c</c>) for <c>-z 3</c> without <c>-H</c>: the histogram of the Z-scores of the library (only for more than 20 library sequences), the size of the library, the parameters and the Kolmogorov-Smirnov statistic; returns the <c>hstring1</c> of <c>-m 10</c>.</summary>
        private static string PrHist(StringBuilder o, Molecule m, CoreResult core)
        {
            int nstats = core.Stats.Count;
            int n0 = core.N0;
            // inithistz( MAXHIST ) and addhistz for every library sequence
            const int minHist = 20, maxHistValue = 120, maxHist = 50;
            int histint = (int)((maxHistValue - minHist + 2) / (double)maxHist + 0.5);
            int maxh = (int)((maxHistValue - minHist + 2) / (double)histint + 0.5);
            var histA = new int[maxh];
            long histEntries = 0;
            int zCalls = 0;
            double a_n0 = n0;
            double a_n0f = Math.Log(m.K * a_n0) / m.H;
            foreach (var (score, n1) in core.Stats)
            {
                double zs = FindZa(m, a_n0, a_n0f, score, n1);
                double rv = Nrv[zCalls++ % 100];
                int zi = (int)(zs + 0.5 + rv);
                if (zi >= 0 && zi <= 120)
                    histEntries++;
                if (zi < minHist)
                    zi = minHist;
                if (zi > maxHistValue)
                    zi = maxHistValue;
                histA[(zi - minHist) / histint]++;
            }

            o.Append('\n');
            double maxDev = 0.0;
            int nChiSq = 0, cumHl = 0, maxI = 0;
            bool maxIset = false;
            if (nstats > 20)
            {
                int mh1 = maxh - 1;
                int mht = (3 * maxh - 3) / 4 - 1;
                if (mh1 > 0)
                {
                    int maxval = 0, maxvalt = 0;
                    for (int i = 0; i < maxh; i++)
                    {
                        if (histA[i] > maxval)
                            maxval = histA[i];
                        if (i >= mht && histA[i] > maxvalt)
                            maxvalt = histA[i];
                    }

                    nChiSq = 0;
                    cumHl = -histA[0];
                    int dotsiz = (maxval - 1) / 60 + 1;
                    int ddotsiz = (maxvalt - 1) / 50 + 1;
                    bool doinset = ddotsiz < dotsiz && dotsiz > 2;
                    o.Append("       opt      E()\n");
                    double prevE = ZsToEc(minHist - histint / 2, histEntries);
                    double curE = 0;
                    var hline = new char[160];
                    for (int i = 0; i <= mh1; i++)
                    {
                        char pch = i == mh1 ? '>' : ' ';
                        pch = i == 0 ? '<' : pch;
                        int hl = histA[i];
                        int hll = hl;
                        cumHl += hl;
                        double fInt = i * histint + minHist + histint / 2.0;
                        curE = ZsToEc(fInt, histEntries);
                        int el, ell, ev;
                        ev = el = ell = (int)(curE - prevE + 0.5);
                        if (hl > 0 && i > 5 && i < (90 - minHist) / histint)
                        {
                            double xTmp = Math.Abs(cumHl - curE);
                            if (xTmp > maxDev)
                            {
                                maxDev = xTmp;
                                maxI = i;
                                maxIset = true;
                            }

                            nChiSq++;
                        }

                        if ((el = (el + dotsiz - 1) / dotsiz) > 60)
                            el = 60;
                        if ((ell = (ell + ddotsiz - 1) / ddotsiz) > 40)
                            ell = 40;
                        o.Append(pch).Append(Iw(i < mh1 ? i * histint + minHist : mh1 * histint + minHist, 3)).Append(' ').Append(Iw(hl, 5)).Append(' ').Append(Iw(ev, 5)).Append(':');
                        if ((hl = (hl + dotsiz - 1) / dotsiz) > 60)
                            hl = 60;
                        if ((hll = (hll + ddotsiz - 1) / ddotsiz) > 40)
                            hll = 40;
                        for (int j = 0; j < hl; j++)
                            hline[j] = '=';
                        if (el <= hl)
                        {
                            if (el > 0)
                                hline[el - 1] = '*';
                            hline[hl] = '\0';
                        }
                        else
                        {
                            for (int j = hl; j < el; j++)
                                hline[j] = ' ';
                            hline[el - 1] = '*';
                            hline[hl = el] = '\0';
                        }

                        if (i == 1)
                        {
                            for (int j = hl; j < 10; j++)
                                hline[j] = ' ';
                            PutC(hline, 10, " one = represents " + I(dotsiz) + " library sequences");
                        }

                        if (doinset && i == mht - 2)
                        {
                            for (int j = hl; j < 10; j++)
                                hline[j] = ' ';
                            PutC(hline, 10, " inset = represents " + I(ddotsiz) + " library sequences");
                        }

                        if (i >= mht && doinset)
                        {
                            for (int j = hl; j < 10; j++)
                                hline[j] = ' ';
                            hline[10] = ':';
                            for (int j = 11; j < 11 + hll; j++)
                                hline[j] = '=';
                            hline[11 + hll] = '\0';
                            if (ell <= hll)
                            {
                                hline[10 + ell] = '*';
                            }
                            else
                            {
                                for (int j = 11 + hll; j < 10 + ell; j++)
                                    hline[j] = ' ';
                                hline[10 + ell] = '*';
                                hline[11 + ell] = '\0';
                            }
                        }

                        o.Append(CString(hline)).Append('\n');
                        prevE = curE;
                    }
                }
            }

            o.Append(Iw(core.DbLength, 7)).Append(" residues in ").Append(Iw(core.DbEntries, 5)).Append(" sequences\n");
            string statInfo = StatInfo(m, n0);
            o.Append(' ').Append(statInfo).Append('\n');
            if (cumHl > 0)
                o.Append(" Kolmogorov-Smirnov  statistic: ").Append(F(maxDev / (float)cumHl, 4).PadLeft(6)).Append(" (N=").Append(I(nChiSq)).Append(") at ").Append(Iw(maxIset ? maxI * histint + minHist : 0, 3)).Append('\n');
            if (cumHl <= 0)
                cumHl = -1;
            string ks = nstats > 20 && maxIset ? "(N=" + I(nChiSq) + ") at " + Iw(maxI * histint + minHist, 3) : "(N=UNINIT) at UNINIT";
            return "; mp_extrap: 60000 " + I(histEntries) + "\n; mp_stats: " + statInfo + "\n; mp_KS: " + F(maxDev / (float)cumHl, 4).PadLeft(6) + " " + ks + "\n";
        }

        private static void PutC(char[] buf, int pos, string s)
        {
            for (int i = 0; i < s.Length; i++)
                buf[pos + i] = s[i];
            buf[pos + s.Length] = '\0';
        }

        /// <summary>
        /// <c>fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000 query library ktup</c> (the command of <c>dndfast7</c>): the standard output as the oracle keeps it, the standard error and the exit status. The query file holds one record
        /// (<paramref name="queryTitle"/>, <paramref name="query"/>), the library file the records of <paramref name="libTitles"/> / <paramref name="library"/> in this order.
        /// </summary>
        internal static ProgramRun Run(string queryTitle, string query, IReadOnlyList<string> libTitles, IReadOnlyList<string> library, bool nFlag, bool sw, int ktup, int mshow = 500000, int ashow = 500000, double eCut = 500000)
        {
            var o = new StringBuilder();
            var err = new StringBuilder();
            string argv = (sw ? "-A " : "") + "-z3 -m10 " + (nFlag ? "-n " : "") + "-Q -b" + I(mshow) + " -E" + I((long)eCut) + " -d" + I(ashow) + " " + QueryPath + " " + LibraryPath + " " + I(ktup);
            o.Append("# fasta34 ").Append(argv).Append('\n');
            o.Append("FASTA searches a protein or DNA sequence data bank\n version 34.26.5 April 26, 2007\nPlease cite:\n W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448\n\n");
            o.Append("Query library ").Append(QueryPath).Append(" vs ").Append(LibraryPath).Append(" library\nsearching ").Append(LibraryPath).Append(" library\n\n");

            // the query: without -n the program reads the letters of both alphabets and decides on the composition (more than 85 percent of ACGTUN in more than 20 residues: DNA)
            bool dna = nFlag;
            int[] q0;
            if (nFlag)
            {
                q0 = Encode(query, Dna.Value.Ascii);
            }
            else
            {
                var raw = new List<char>();
                foreach (char ch in query)
                {
                    if ((char.IsAsciiLetter(ch) && "JOjo".IndexOf(ch) < 0) || ch == '*')
                        raw.Add(ch);
                }

                int nraw = raw.Count;
                int acgtun = raw.Count(ch => "ACGTUNacgtun".Contains(ch));
                if (nraw > 20 && acgtun / (float)nraw > 0.85)
                {
                    dna = true;
                    var codes = new List<int> { 0 };
                    var dm = Dna.Value;
                    for (int i = 0; i < raw.Count; i++)
                    {
                        char ch = raw[i];
                        int code = ch == '*' ? 24 : dm.Ascii[ch];
                        if (code > 34 && code != 49)
                        {
                            err.Append("*** Warning - unrecognized residue at ").Append(I(i)).Append(':').Append(ch).Append(" - ").Append(Iw(ch, 2)).Append('\n');
                            code = dm.Ascii['X'];
                        }

                        codes.Add(code);
                    }

                    codes.Add(0);
                    q0 = codes.ToArray();
                }
                else
                {
                    q0 = Encode(query, Protein.Value.Ascii);
                }
            }

            var m = WithKtup(dna ? Dna.Value : Protein.Value, ktup);
            int n0 = q0.Length - 2;
            if (n0 <= 0)
            {
                err.Append("\n***[fasta34] Query sequence length <= 0: ").Append(QueryPath).Append("***\n");
                return new ProgramRun(o.ToString(), err.ToString(), 1);
            }

            var lib = library.Select(x => Encode(x, m.Ascii)).ToList();
            var core = SearchCore(m, q0, lib, mshow, eCut);
            string sqnam = dna ? "nt" : "aa";
            string sqtype = dna ? "DNA" : "protein";
            string qtitle = queryTitle + " " + I(n0) + " " + sqnam;
            string qlabel = queryTitle.Split(' ')[0];
            o.Append(Iw(1, 3)).Append(">>>").Append(qtitle).Append(" - ").Append(I(n0)).Append(' ').Append(sqnam).Append("\n vs  ").Append(LibraryPath).Append(" library\n");
            string hstring1 = PrHist(o, m, core);
            int pamH = m.PamH, pamL = m.PamL;
            string pamfile = dna ? "+5/-4" : "BL50";
            var f0 = core.Frames[0];
            string gstring2 = "FASTA (3.5 Sept 2006) function [optimized, " + pamfile + " matrix (" + I(pamH) + ":" + I(pamL) + ")] ktup: " + I(ktup) + "\n join: " + I(f0.Cgap) + ", opt: " + I(f0.OptCut) + ", open/ext: " + I(m.GDel) + "/" + I(m.GGap) + ", width: " + Iw(m.OptWid, 3);
            string gstring3 = "; pg_name: FASTA\n; pg_ver: 3.5 Sept 2006\n; pg_matrix: " + pamfile + " (" + I(pamH) + ":" + I(pamL) + ")\n; pg_open-ext: " + I(m.GDel) + " " + I(m.GGap) + "\n; pg_ktup: " + I(ktup) + "\n; pg_optcut: " + I(f0.OptCut) + "\n; pg_cgap: " + I(f0.Cgap) + "\n";
            o.Append('\n').Append(gstring2).Append('\n');
            o.Append(" Scan time:  T\n");

            // showbest
            var listing = core.Listing;
            int zdbSize = core.DbEntries;
            const int llen = 60, rMargin = 10;
            int width = dna ? llen - (rMargin + 4) : llen - rMargin;
            if (!(core.Entries.Count == 0 || core.Entries[0].E < eCut))
            {
                o.Append("!! No library sequences with E() < ").Append(FormatG(eCut, 2)).Append('\n');
            }
            else
            {
                o.Append("The best scores are:").Append(new string(' ', llen - (rMargin + 12))).Append("opt bits E(").Append(I(zdbSize)).Append(")\n");
                foreach (var h in listing.Hits)
                {
                    string title = libTitles[h.Sequence];
                    if (title.Length > width)
                        title = title[..width];
                    o.Append(title.PadRight(width)).Append(" (").Append(Iw(h.Length, 4)).Append(')');
                    if (dna)
                        o.Append(h.Frame > 0 ? " [r]" : " [f]");
                    o.Append(' ').Append(Iw(h.Opt, 4)).Append(' ').Append(F(h.Bits, 1)).Append(' ').Append(FormatG(h.EValue, 2).PadLeft(7)).Append('\n');
                }

                if (core.Nshow > 0)
                {
                    o.Append("\n>>>").Append(qlabel).Append(", ").Append(I(n0)).Append(' ').Append(sqnam).Append(" vs ").Append(LibraryPath).Append(" library\n");
                    o.Append("; pg_name: fasta34\n; pg_ver: 34.26\n; pg_argv: fasta34 ").Append(argv).Append('\n');
                    o.Append(gstring3).Append(hstring1);
                    ShowAlign(o, err, core, m, dna, sw, libTitles, qlabel, sqtype, ashow);
                }
            }

            o.Append(">>><<<\n");
            o.Append("\n\n").Append(I(n0)).Append(" residues in 1 query   sequences\n").Append(I(core.DbLength)).Append(" residues in ").Append(I(core.DbEntries)).Append(" library sequences\n");
            o.Append(" Scomplib [34.26]\n start: T done: T\n Total Scan time:  T Total Display time:  T\n\nFunction used was FASTA [version 34.26.5 April 26, 2007]\n");
            return new ProgramRun(o.ToString(), err.ToString(), 0);
        }

        /// <summary><c>showalign</c> of <c>mshowalign.c</c> (<c>-m 10</c>): for every entry of the list that is shown, the scores and the alignment as two records.</summary>
        private static void ShowAlign(StringBuilder o, StringBuilder err, CoreResult core, Molecule m, bool dna, bool sw, IReadOnlyList<string> libTitles, string qlabel, string sqtype, int ashow)
        {
            int n0 = core.N0;
            string name0 = qlabel[..Math.Min(12, qlabel.Length)];
            string f_id0 = "fa";
            string f_id1 = dna && !sw ? "bs" : "sw";
            bool swFlag = !dna || sw;
            int istop = Math.Min(Math.Min(core.Entries.Count, ashow), core.Nshow);
            for (int ib = 0; ib < istop; ib++)
            {
                var bbp = core.Entries[ib];
                if (bbp.Initn <= 0)
                    break;
                int[] aa1 = core.Library[bbp.Seq];
                int n1 = bbp.N1;
                var f = core.Frames[bbp.Frame];
                var ares = new AlnRes();
                int swScore = DoWalign(m, f, aa1, n1, swFlag, ares, err);
                var aln = new Aln { Qlrev = bbp.Frame > 0 };
                string bline = libTitles[bbp.Seq].Replace('\n', ' ');
                o.Append(">>").Append(bline).Append('\n');
                o.Append("; ").Append(f_id0).Append("_frame: ").Append(bbp.Frame > 0 ? 'r' : 'f').Append('\n');
                o.Append("; ").Append(f_id0).Append("_initn: ").Append(Iw(bbp.Initn, 3)).Append('\n');
                o.Append("; ").Append(f_id0).Append("_init1: ").Append(Iw(bbp.Init1, 3)).Append('\n');
                o.Append("; ").Append(f_id0).Append("_opt: ").Append(Iw(bbp.Opt, 3)).Append('\n');
                o.Append("; ").Append(f_id0).Append("_z-score: ").Append(F(bbp.Z, 1).PadLeft(4)).Append('\n');
                o.Append("; ").Append(f_id0).Append("_bits: ").Append(F(ZsToBit(bbp.Z, n0, bbp.N1), 1).PadLeft(3)).Append('\n');
                o.Append("; ").Append(f_id0).Append("_expect: ").Append(FormatG(bbp.E, 2).PadLeft(6)).Append('\n');
                int maxc = ares.Nres + n0 + n1 + 4 * aln.Llen + 64;
                var seqc0 = new char[maxc];
                var seqc1 = new char[maxc];
                var seqca = new char[maxc];
                int[] aa0 = bbp.Frame > 0 ? core.Frames[1].Aa0 : core.Frames[0].Aa0;
                int nc = Calcons(m, aa0, n0, aa1, n1, out int lc, aln, ares, seqc0, seqc1, seqca);
                float percent = lc > 0 ? (float)((100.0 * (float)aln.Nident) / (float)lc) : -1.00f;
                float gpercent = lc > 0 ? (float)((100.0 * (float)aln.Nsim) / (float)lc) : -1.00f;
                CalCoord(n0, bbp.N1, aln);
                if (ares.Nres > 0)
                {
                    if (swFlag)
                        o.Append("; ").Append(f_id1).Append("_score: ").Append(I(swScore)).Append('\n');
                    o.Append("; ").Append(f_id1).Append("_ident: ").Append(F(percent / 100.0, 3).PadLeft(5)).Append('\n');
                    o.Append("; ").Append(f_id1).Append("_sim: ").Append(F(gpercent / 100.0, 3).PadLeft(5)).Append('\n');
                    o.Append("; ").Append(f_id1).Append("_overlap: ").Append(I(lc)).Append('\n');
                    string name1 = bline.Split(' ')[0];
                    Discons(o, seqc0, seqc1, nc, n0, bbp.N1, name0, name1, sqtype, aln);
                }
            }
        }
    }

    #endregion
}
