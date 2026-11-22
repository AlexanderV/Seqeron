using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: the arguments of the binaries and constants.c

    /// <summary>
    /// The state that MAFFT's binaries keep in global variables (<c>defs.c</c>, <c>arguments()</c> of <c>disttbfast.c</c>, <c>constants()</c> of <c>constants.c</c>): the parsed
    /// command line of <c>disttbfast</c> (the script <c>mafft</c> turns every option into such an argument list), the scoring matrices and the penalties scaled to integers.
    /// </summary>
    internal sealed partial class MafftGlobals
    {
        internal const int NotSpecified = 100009;

        /// <summary>The messages the programs write to standard error (<c>reporterr</c>).</summary>
        internal readonly System.Text.StringBuilder Log = new();

        /// <summary><c>modelname</c> of <c>constants()</c> (shown at the end of the run).</summary>
        internal string ModelName = "";

        internal void Report(string s) => Log.Append(s);

        internal static string F(double v, int width, int dec, bool plus = false)
        {
            string r = v.ToString("F" + dec.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
            if (plus && !r.StartsWith('-'))
                r = "+" + r;
            return r.PadLeft(width);
        }

        // ---- arguments() defaults (disttbfast.c)
        internal int Nthread = 1, NthreadPair = 1, NthreadTb = 1;
        internal int Nadd;
        internal bool AddProfile = true;
        internal int Nblosum = 62;
        internal int Fmodel;
        internal bool UseFft, ForceFft;
        internal int FftScore = 1;
        internal int Weight = 3;
        internal char Alg = 'A';
        internal int Outgap = 1;
        internal int Scoremtx = 1;
        internal int TreeIn, TreeOut, DistOut, NoAlign;
        internal int Tbrweight = 3;
        internal char TreeMethod = 'X';
        internal double SueffGlobal = 0.1;
        internal char Dorp;
        internal int PpenaltyDist = NotSpecified;
        internal int Ppenalty = -1530;
        internal int PpenaltyEx = NotSpecified;
        internal double PenaltyShiftFactor = 1000.0;
        internal int Poffset = -123;
        internal int KimuraR = NotSpecified;
        internal int PamN = NotSpecified;
        internal bool Tm;
        internal int FftWinSize = NotSpecified, FftThreshold = NotSpecified;
        internal int TupleSize = 6;
        internal int LegacyGapCost;
        internal double SpecificityConsideration;
        internal int NGuideTree = 1;
        internal int NWildcard;
        internal int NeverMemSave;
        internal bool OneIteration;
        internal int Disp;
        internal int KeepLength, MapOut, Smoothing;
        internal int SpScoreOut;
        internal int OutNumber;
        internal int Cycle = 1;

        // ---- constants()
        internal int Nalphabets = 26;
        internal int NscoredAlphabets = 20;
        internal int Penalty, PenaltyEx, PenaltyDist, PenaltyShift, Offset, OffsetFft;
        internal int[] AminoN = new int[256];
        internal int[] AminoGrp = new int[256];
        internal char[] Amino = new char[26];
        internal int[][] NDis = null!;
        internal double[][] NDisConsweightMulti = null!;
        internal int[][] NDisFft = null!;
        internal double[][] AminoDisConsweightMulti = null!;
        internal int[][] AminoDis = null!;
        internal double[] Polarity = new double[256];
        internal double[] Volume = new double[256];
        internal bool TryWarp;
        internal double ConsweightMulti = 1.0;

        /// <summary>The arguments of <c>disttbfast</c> as the script <c>mafft</c> passes them (<c>arguments()</c> of <c>disttbfast.c</c>, the options that the script can produce).</summary>
        internal static MafftGlobals Parse(IReadOnlyList<string> argv)
        {
            var g = new MafftGlobals();
            int k = 0;
            string Next() => k < argv.Count ? argv[k++] : throw new ArgumentException("A MAFFT option is missing its argument.");
            while (k < argv.Count && argv[k].Length > 0 && argv[k][0] == '-')
            {
                string a = argv[k++];
                bool rest = false; // an option with an argument ends the bundle of flags ("goto nextoption")
                for (int ci = 1; ci < a.Length && !rest; ci++)
                {
                    char c = a[ci];
                    switch (c)
                    {
                        case 'V': g.PpenaltyDist = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'f': g.Ppenalty = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'Q': g.PenaltyShiftFactor = CAtof(Next()); rest = true; break;
                        case 'g': g.PpenaltyEx = (int)(CAtof(Next()) * 1000 - 0.5); g.Report("ppenalty_ex = " + g.PpenaltyEx.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'h': g.Poffset = (int)(CAtof(Next()) * 1000 - 0.5); rest = true; break;
                        case 'k': g.KimuraR = CAtoi(Next()); g.Report("kappa = " + g.KimuraR.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'b': g.Nblosum = CAtoi(Next()); g.Scoremtx = 1; rest = true; break;
                        case 'j': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = false; g.Report("jtt/kimura " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'm': g.PamN = CAtoi(Next()); g.Scoremtx = 0; g.Tm = true; g.Report("tm " + g.PamN.ToString(CultureInfo.InvariantCulture) + "\n"); rest = true; break;
                        case 'C':
                            {
                                string v = Next();
                                g.NthreadPair = g.Nthread = CAtoi(v);
                                g.Report("nthread = " + g.Nthread.ToString(CultureInfo.InvariantCulture) + "\nnthreadpair = " + g.Nthread.ToString(CultureInfo.InvariantCulture) + "\n");
                                int dash = v.IndexOf('-');
                                g.NthreadTb = dash >= 0 ? CAtoi(v[(dash + 1)..]) : g.Nthread;
                                g.Report("nthreadtb = " + g.NthreadTb.ToString(CultureInfo.InvariantCulture) + "\n");
                                rest = true;
                                break;
                            }
                        case 's': g.SpecificityConsideration = CAtof(Next()); rest = true; break;
                        case 'a': g.Fmodel = 1; break;
                        case 'K': g.AddProfile = false; break;
                        case 'y': g.DistOut = 1; break;
                        case 't': g.TreeOut = 1; break;
                        case '^': g.TreeOut = 2; break;
                        case 'T': g.NoAlign = 1; break;
                        case 'r': g.OneIteration = true; break;
                        case 'D': g.Dorp = 'd'; break;
                        case 'P': g.Dorp = 'p'; break;
                        case 'L': g.LegacyGapCost = 1; break;
                        case 'e': g.FftScore = 0; break;
                        case 'x': CAtof(Next()); rest = true; break;
                        case 'q': CAtoi(Next()); rest = true; break;
                        case 'X': g.TreeMethod = 'X'; g.SueffGlobal = CAtof(Next()); rest = true; break;
                        case 'E': g.NGuideTree = CAtoi(Next()); rest = true; break;
                        case 'A': g.Alg = 'A'; break;
                        case '&': g.Alg = 'a'; break;
                        case '@': g.Alg = 'd'; break;
                        case 'N': g.NeverMemSave = 1; break;
                        case 'M': g.Alg = 'M'; break;
                        case 'S': g.SpScoreOut = 1; break;
                        case 'B': break;
                        case 'F': g.UseFft = true; break;
                        case 'G': g.UseFft = true; g.ForceFft = true; break;
                        case 'U': g.TreeIn = 1; break;
                        case 'u': g.Weight = 0; g.Tbrweight = 0; break;
                        case 'v': g.Tbrweight = 3; break;
                        case 'd': g.Disp = 1; break;
                        case 'O': g.Outgap = 0; break;
                        case 'w': g.FftWinSize = CAtoi(Next()); rest = true; break;
                        case 'W': g.TupleSize = CAtoi(Next()); rest = true; break;
                        case 'n': g.OutNumber = 1; break;
                        case ':': g.NWildcard = 1; break;
                        default:
                            throw new NotSupportedException($"The MAFFT option -{c} is not available in this port.");
                    }
                }
            }
            return g;
        }

        /// <summary><c>atof</c> of the C library: the longest numeric prefix (decimal, optional sign, fraction and exponent), 0 when there is none.</summary>
        internal static double CAtof(string s)
        {
            int i = 0, n = s.Length;
            while (i < n && char.IsWhiteSpace(s[i]))
                i++;
            int start = i;
            if (i < n && (s[i] == '+' || s[i] == '-'))
                i++;
            int digits = 0;
            while (i < n && s[i] >= '0' && s[i] <= '9') { i++; digits++; }
            if (i < n && s[i] == '.')
            {
                i++;
                while (i < n && s[i] >= '0' && s[i] <= '9') { i++; digits++; }
            }
            if (digits == 0)
                return 0.0;
            if (i < n && (s[i] == 'e' || s[i] == 'E'))
            {
                int j = i + 1;
                if (j < n && (s[j] == '+' || s[j] == '-'))
                    j++;
                int ed = 0;
                while (j < n && s[j] >= '0' && s[j] <= '9') { j++; ed++; }
                if (ed > 0)
                    i = j;
            }
            return double.Parse(s[start..i], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture);
        }

        /// <summary><c>myatoi</c> / <c>atoi</c>: an optional sign and digits, 0 when there are none.</summary>
        internal static int CAtoi(string s)
        {
            int i = 0, n = s.Length;
            while (i < n && char.IsWhiteSpace(s[i]))
                i++;
            bool neg = false;
            if (i < n && (s[i] == '+' || s[i] == '-'))
            {
                neg = s[i] == '-';
                i++;
            }
            long v = 0;
            while (i < n && s[i] >= '0' && s[i] <= '9')
            {
                v = v * 10 + (s[i] - '0');
                if (v > int.MaxValue)
                    v = int.MaxValue;
                i++;
            }
            return (int)(neg ? -v : v);
        }

        private static int Shishagonyuu(double v)
        {
            if (v > 0.0)
                return (int)(v + 0.5);
            if (v == 0.0)
                return 0;
            if (v < 0.0)
                return (int)(v - 0.5);
            return 0;
        }

        /// <summary><c>constants( nseq, seq )</c> of <c>constants.c</c>: the scoring matrices, penalties, the FFT tables of the alphabet in use. <paramref name="seq"/> are the sequences (used by <c>--fmodel</c> only).</summary>
        internal void Constants(IReadOnlyList<string> seq)
        {
            TryWarp = PenaltyShiftFactor < 10;
            if (Nblosum < 0)
                Dorp = 'p';
            for (int i = 0; i < 256; i++)
                AminoN[i] = -1;
            int[,] nDis;
            int charsize;
            if (Dorp == 'd')
            {
                charsize = 0x80;
                NscoredAlphabets = 10;
                Scoremtx = -1;
                if (Ppenalty == NotSpecified) Ppenalty = -1530;
                if (PenaltyDistIsUnset()) PpenaltyDist = Ppenalty;
                if (PpenaltyEx == NotSpecified) PpenaltyEx = 0;
                if (Poffset == NotSpecified) Poffset = -369;
                if (PamN == NotSpecified) PamN = 200;
                if (KimuraR == NotSpecified) KimuraR = 2;
                Penalty = (int)(3 * 600.0 / 1000.0 * Ppenalty + 0.5);
                PenaltyDist = (int)(3 * 600.0 / 1000.0 * PpenaltyDist + 0.5);
                PenaltyShift = (int)(PenaltyShiftFactor * Penalty);
                PenaltyEx = (int)(3 * 600.0 / 1000.0 * PpenaltyEx + 0.5);
                Offset = (int)(1 * 600.0 / 1000.0 * Poffset + 0.5);
                OffsetFft = (int)(1 * 600.0 / 1000.0 * (-0) + 0.5);
                string amino = MafftNucAlphabet;
                for (int i = 0; i < 26; i++)
                    Amino[i] = amino[i];
                for (int i = 0; i < 26; i++)
                    AminoN[Amino[i]] = i;
                var freq = new double[4];
                if (Fmodel == 1)
                {
                    CalcFreqNuc(seq, freq);
                    Report("a, freq[0] = " + freq[0].ToString("F6", CultureInfo.InvariantCulture) + "\ng, freq[1] = " + freq[1].ToString("F6", CultureInfo.InvariantCulture) + "\nc, freq[2] = " + freq[2].ToString("F6", CultureInfo.InvariantCulture) + "\nt, freq[3] = " + freq[3].ToString("F6", CultureInfo.InvariantCulture) + "\n");
                }
                else
                {
                    freq[0] = freq[1] = freq[2] = freq[3] = 0.25;
                }
                if (KimuraR == 9999)
                    throw new NotSupportedException("--kappa 9999 (the fixed nucleotide matrix of MAFFT) is not available in this port.");
                Report("generating a scoring matrix for nucleotide (dist=" + PamN.ToString(CultureInfo.InvariantCulture) + ") ... ");
                nDis = NucleotideMatrix(freq);
                Report("done\n");
                ModelName = "DNA" + PamN.ToString(CultureInfo.InvariantCulture) + " (" + KimuraR.ToString(CultureInfo.InvariantCulture) + "), " + F(-(double)Ppenalty * 0.001, 4, 2) + " (" + F(-(double)Ppenalty * 0.003, 4, 2) + "), "
                    + F(-(double)Poffset * 0.001, 4, 2) + " (" + F(-(double)Poffset * 0.003, 4, 2) + "), " + (TryWarp ? F(-(double)PenaltyShift / 1800, 4, 2) + " (" + F(-(double)PenaltyShift / 600, 4, 2) + ")" : "noshift");
                int[] grp = { 0, 1, 2, 3, 3, 0, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };
                for (int i = 0; i < 26; i++)
                    AminoGrp[Amino[i]] = grp[i];
                AmbiguousScore(nDis);
                if (NWildcard != 0)
                    NScore(nDis);
            }
            else if (Dorp == 'p' && Scoremtx == 1 && Nblosum >= 0)
            {
                charsize = 0x80;
                if (Ppenalty == NotSpecified) Ppenalty = -1530;
                if (PenaltyDistIsUnset()) PpenaltyDist = Ppenalty;
                if (PpenaltyEx == NotSpecified) PpenaltyEx = 0;
                if (Poffset == NotSpecified) Poffset = -123;
                if (PamN == NotSpecified) PamN = 0;
                if (KimuraR == NotSpecified) KimuraR = 1;
                SetProteinPenalties();
                string amino = MafftAminoAlphabet;
                for (int i = 0; i < 26; i++)
                    Amino[i] = amino[i];
                int[] grp = { 0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2, 6, 6, 6, 1 };
                for (int i = 0; i < 26; i++)
                    AminoGrp[Amino[i]] = grp[i];
                for (int i = 0; i < 26; i++)
                    AminoN[Amino[i]] = i;
                Report("rescale = 1\n");
                nDis = BlosumMatrix(seq);
                ModelName = "BLOSUM" + Nblosum.ToString(CultureInfo.InvariantCulture) + ", " + F(-(double)Ppenalty / 1000, 4, 2) + ", " + F(-(double)Poffset / 1000, 4, 2, true) + ", " + F(-(double)PpenaltyEx / 1000, 4, 2, true) + ", "
                    + (TryWarp ? F(-(double)PenaltyShift / 600, 4, 2) : "noshift");
            }
            else if (Dorp == 'p' && Scoremtx == 1)
            {
                throw new NotSupportedException("A user-defined amino-acid matrix (--aamatrix) / the text mode of MAFFT is not available in this port.");
            }
            else if (Dorp == 'p' && Scoremtx == 2)
            {
                throw new NotSupportedException("Not supported (MAFFT: Miyata-Yasunaga).");
            }
            else
            {
                charsize = 0x80;
                if (Ppenalty == NotSpecified) Ppenalty = -1530;
                if (PenaltyDistIsUnset()) PpenaltyDist = Ppenalty;
                if (PpenaltyEx == NotSpecified) PpenaltyEx = 0;
                if (Poffset == NotSpecified) Poffset = -123;
                if (PamN == NotSpecified) PamN = 200;
                if (KimuraR == NotSpecified) KimuraR = 1;
                if (PamN == 0)
                    throw new ArgumentException("pamN=0??");
                bool makeAverage0 = true;
                int pamn = PamN;
                if (pamn < 0)
                {
                    pamn = -pamn;
                    makeAverage0 = false;
                }
                SetProteinPenalties();
                string amino = MafftAminoAlphabet;
                for (int i = 0; i < 26; i++)
                    Amino[i] = amino[i];
                int[] grp = { 0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2, 6, 6, 6, 1 };
                for (int i = 0; i < 26; i++)
                    AminoGrp[Amino[i]] = grp[i];
                for (int i = 0; i < 26; i++)
                    AminoN[Amino[i]] = i;
                Report("generating " + pamn.ToString(CultureInfo.InvariantCulture) + "PAM " + (Tm ? "Transmembrane" : "JTT") + " scoring matrix for amino acids ... ");
                nDis = JttMatrix(seq, pamn, makeAverage0);
                Report("done.\n");
                ModelName = (Tm ? "Transmembrane" : "JTT") + " " + pamn.ToString(CultureInfo.InvariantCulture) + "PAM, " + F(-(double)Ppenalty / 1000, 4, 2) + ", " + F(-(double)Poffset / 1000, 4, 2) + ", "
                    + (TryWarp ? F(-(double)PenaltyShift / 600, 4, 2) : "noshift");
            }
            _ = charsize;
            NDis = new int[26][];
            for (int i = 0; i < 26; i++)
            {
                NDis[i] = new int[26];
                for (int j = 0; j < 26; j++)
                    NDis[i][j] = nDis[i, j];
            }
            AminoDis = new int[256][];
            AminoDisConsweightMulti = new double[256][];
            for (int i = 0; i < 256; i++)
            {
                AminoDis[i] = new int[256];
                AminoDisConsweightMulti[i] = new double[256];
            }
            NDisConsweightMulti = new double[26][];
            NDisFft = new int[26][];
            for (int i = 0; i < 26; i++)
            {
                NDisConsweightMulti[i] = new double[26];
                NDisFft[i] = new int[26];
            }
            for (int i = 0; i < Nalphabets; i++)
            {
                for (int j = 0; j < Nalphabets; j++)
                {
                    AminoDis[Amino[i]][Amino[j]] = NDis[i][j];
                    NDisConsweightMulti[i][j] = NDis[i][j] * ConsweightMulti;
                    AminoDisConsweightMulti[Amino[i]][Amino[j]] = NDis[i][j] * ConsweightMulti;
                }
            }
            if (Dorp == 'd')
            {
                for (int i = 0; i < 10; i++)
                    for (int j = 0; j < 10; j++)
                        NDisFft[i][j] = NDis[i][j] + Offset - OffsetFft;
            }
            else
            {
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        NDisFft[i][j] = NDis[i][j] + Offset - OffsetFft;
            }
            if (FftThreshold == NotSpecified)
                FftThreshold = 80;
            if (FftWinSize == NotSpecified)
                FftWinSize = Dorp == 'd' ? 100 : 20;
            if (FftScore != 0)
            {
                Normalise(MafftTables2.Polarity, Polarity);
                Normalise(MafftTables2.Volume, Volume);
            }
        }

        private bool PenaltyDistIsUnset() => PpenaltyDist == NotSpecified;

        private static void Normalise(double[] source, double[] target)
        {
            for (int i = 0; i < 20; i++)
                target[i] = source[i];
            double av = 0.0;
            for (int i = 0; i < 20; i++)
                av += target[i];
            av /= 20.0;
            double sd = 0.0;
            for (int i = 0; i < 20; i++)
                sd += (target[i] - av) * (target[i] - av);
            sd /= 20.0;
            sd = Math.Sqrt(sd);
            for (int i = 0; i < 20; i++)
                target[i] -= av;
            for (int i = 0; i < 20; i++)
                target[i] /= sd;
        }

        private void SetProteinPenalties()
        {
            Penalty = (int)(600.0 / 1000.0 * Ppenalty + 0.5);
            PenaltyDist = (int)(600.0 / 1000.0 * PpenaltyDist + 0.5);
            PenaltyShift = (int)(PenaltyShiftFactor * Penalty);
            PenaltyEx = (int)(600.0 / 1000.0 * PpenaltyEx + 0.5);
            Offset = (int)(600.0 / 1000.0 * Poffset + 0.5);
            OffsetFft = (int)(600.0 / 1000.0 * (-0) + 0.5);
        }

        private void CalcFreqNuc(IReadOnlyList<string> seq, double[] datafreq)
        {
            double total = 0.0;
            foreach (var s in seq)
            {
                foreach (char c in s)
                {
                    int aan = AminoN[c];
                    if (aan == 4)
                        aan = 3;
                    if (aan >= 0 && aan < 4)
                    {
                        datafreq[aan] += 1.0;
                        total += 1.0;
                    }
                }
            }
            total = 0.0;
            for (int i = 0; i < 4; i++)
                total += datafreq[i];
            for (int i = 0; i < 4; i++)
                datafreq[i] /= total;
            for (int i = 0; i < 4; i++)
                if (datafreq[i] < 0.0001)
                    datafreq[i] = 0.0001;
            total = 0.0;
            for (int i = 0; i < 4; i++)
                total += datafreq[i];
            for (int i = 0; i < 4; i++)
                datafreq[i] /= total;
        }

        private void CalcFreq(IReadOnlyList<string> seq, double[] datafreq)
        {
            for (int i = 0; i < NscoredAlphabets; i++)
                datafreq[i] = 0.0;
            foreach (var s in seq)
            {
                foreach (char c in s)
                {
                    int aan = AminoN[c];
                    if (aan >= 0 && aan < NscoredAlphabets && c != '-')
                        datafreq[aan] += 1.0;
                }
            }
            double total = 0.0;
            for (int i = 0; i < NscoredAlphabets; i++)
                total += datafreq[i];
            for (int i = 0; i < NscoredAlphabets; i++)
                datafreq[i] /= total;
            for (int i = 0; i < NscoredAlphabets; i++)
                if (datafreq[i] < 0.0001)
                    datafreq[i] = 0.0001;
            total = 0.0;
            for (int i = 0; i < NscoredAlphabets; i++)
                total += datafreq[i];
            for (int i = 0; i < NscoredAlphabets; i++)
                datafreq[i] /= total;
        }

        private static void MtxMltDouble(double[,] a, double[,] b, int n)
        {
            var row = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < n; k++)
                    row[k] = a[i, k];
                for (int j = 0; j < n; j++)
                {
                    double s = 0.0;
                    for (int k = 0; k < n; k++)
                        s += row[k] * b[k, j];
                    a[i, j] = s;
                }
            }
        }

        /// <summary>The Kimura two-parameter PAM matrix (<c>generatenuc1pam</c>, <c>MtxmltDouble</c> pamN times, log10 * 1000, centred, scaled to 600 / mean of the diagonal).</summary>
        private int[,] NucleotideMatrix(double[] freq)
        {
            var r = new double[4, 4];
            double kr = KimuraR;
            r[0, 0] = 0.0; r[0, 1] = kr; r[0, 2] = 1.0; r[0, 3] = 1.0;
            r[1, 0] = kr; r[1, 1] = 0.0; r[1, 2] = 1.0; r[1, 3] = 1.0;
            r[2, 0] = 1.0; r[2, 1] = 1.0; r[2, 2] = 0.0; r[2, 3] = kr;
            r[3, 0] = 1.0; r[3, 1] = 1.0; r[3, 2] = kr; r[3, 3] = 0.0;
            var mut = new double[4];
            double total = 0.0;
            for (int i = 0; i < 4; i++)
            {
                double tmp = 0.0;
                for (int j = 0; j < 4; j++)
                    tmp += r[i, j] * freq[j];
                mut[i] = tmp;
                total += tmp * freq[i];
            }
            var pam1 = new double[4, 4];
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pam1[i, j] = i != j ? 0.01 / total * r[i, j] * freq[j] : 1.0 - 0.01 / total * mut[i];
            var pamx = new double[11, 11];
            for (int i = 0; i < 4; i++)
                pamx[i, i] = 1.0;
            for (int x = 0; x < PamN; x++)
                MtxMltDouble(pamx, pam1, 4);
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pamx[i, j] /= freq[j];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    if (pamx[i, j] == 0.0)
                        pamx[i, j] = 0.00001;
                    pamx[i, j] = Math.Log10(pamx[i, j]) * 1000.0;
                }
            }
            double average = 0.0;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    average += pamx[i, j] * freq[i] * freq[j];
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pamx[i, j] -= average;
            average = 0.0;
            for (int i = 0; i < 4; i++)
                average += pamx[i, i] * 1.0 / 4.0;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pamx[i, j] *= 600.0 / average;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pamx[i, j] -= Offset;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    pamx[i, j] = Shishagonyuu(pamx[i, j]);
            for (int i = 0; i < 5; i++)
            {
                pamx[4, i] = pamx[3, i];
                pamx[i, 4] = pamx[i, 3];
            }
            for (int i = 5; i < 10; i++)
                for (int j = 5; j < 10; j++)
                    pamx[i, j] = pamx[i - 5, j - 5];
            var n = new int[26, 26];
            for (int i = 0; i < 10; i++)
                for (int j = 0; j < 10; j++)
                    n[i, j] = Shishagonyuu(pamx[i, j]);
            return n;
        }

        private void NScore(int[,] n)
        {
            int nn = AminoN['n'];
            for (int i = 0; i < 26; i++)
            {
                n[i, nn] = Shishagonyuu(0.25 * n[i, i]);
                n[nn, i] = n[i, nn];
            }
            n[nn, nn] = Shishagonyuu(0.25 * (n[0, 0] + n[1, 1] + n[2, 2] + n[3, 3]));
        }

        private void AmbiguousScore(int[,] n)
        {
            int Idx(char c) => AminoN[c];
            int a = Idx('a'), g = Idx('g'), c = Idx('c'), t = Idx('t');
            int r = Idx('r'), y = Idx('y'), k = Idx('k'), m = Idx('m'), s = Idx('s'), w = Idx('w'), b = Idx('b'), d = Idx('d'), h = Idx('h'), v = Idx('v');
            for (int i = 0; i < 26; i++)
            {
                n[i, r] = Shishagonyuu((double)1 / 2 * (n[a, i] + n[g, i]));
                n[i, y] = Shishagonyuu((double)1 / 2 * (n[c, i] + n[t, i]));
                n[i, k] = Shishagonyuu((double)1 / 2 * (n[g, i] + n[t, i]));
                n[i, m] = Shishagonyuu((double)1 / 2 * (n[a, i] + n[c, i]));
                n[i, s] = Shishagonyuu((double)1 / 2 * (n[g, i] + n[c, i]));
                n[i, w] = Shishagonyuu((double)1 / 2 * (n[a, i] + n[t, i]));
                n[i, b] = Shishagonyuu((double)1 / 3 * (n[c, i] + n[g, i] + n[t, i]));
                n[i, d] = Shishagonyuu((double)1 / 3 * (n[a, i] + n[g, i] + n[t, i]));
                n[i, h] = Shishagonyuu((double)1 / 3 * (n[a, i] + n[c, i] + n[t, i]));
                n[i, v] = Shishagonyuu((double)1 / 3 * (n[a, i] + n[c, i] + n[g, i]));
                n[r, i] = n[i, r];
                n[y, i] = n[i, y];
                n[k, i] = n[i, k];
                n[m, i] = n[i, m];
                n[s, i] = n[i, s];
                n[w, i] = n[i, w];
                n[b, i] = n[i, b];
                n[d, i] = n[i, d];
                n[h, i] = n[i, h];
                n[v, i] = n[i, v];
            }
            n[r, r] = Shishagonyuu((double)1 / 2 * (n[a, a] + n[g, g]));
            n[y, y] = Shishagonyuu((double)1 / 2 * (n[c, c] + n[t, t]));
            n[k, k] = Shishagonyuu((double)1 / 2 * (n[g, g] + n[t, t]));
            n[m, m] = Shishagonyuu((double)1 / 2 * (n[a, a] + n[c, c]));
            n[s, s] = Shishagonyuu((double)1 / 2 * (n[g, g] + n[c, c]));
            n[w, w] = Shishagonyuu((double)1 / 2 * (n[a, a] + n[t, t]));
            n[b, b] = Shishagonyuu((double)1 / 3 * (n[c, c] + n[g, g] + n[t, t]));
            n[d, d] = Shishagonyuu((double)1 / 3 * (n[a, a] + n[g, g] + n[t, t]));
            n[h, h] = Shishagonyuu((double)1 / 3 * (n[a, a] + n[c, c] + n[t, t]));
            n[v, v] = Shishagonyuu((double)1 / 3 * (n[a, a] + n[c, c] + n[g, g]));
        }

        /// <summary>The BLOSUM branch of <c>constants()</c> (<c>--bl</c>, default 62): the matrix scaled to 600 / (frequency-weighted mean of the diagonal), shifted by the offset, rounded to integers.</summary>
        private int[,] BlosumMatrix(IReadOnlyList<string> seq)
        {
            var tri = MafftTables.BlosumLowerTriangle[Nblosum];
            double[] freq = MafftTables.BlosumFrequencies;
            double[] freq1 = freq;
            if (Fmodel == 1)
            {
                var datafreq = new double[20];
                CalcFreq(seq, datafreq);
                freq1 = datafreq;
            }
            var d = new double[20, 20];
            int count = 0;
            for (int i = 0; i < 20; i++)
                for (int j = 0; j <= i; j++)
                    d[i, j] = d[j, i] = tri[count++];
            double average = 0.0;
            if (Fmodel != -1)
            {
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        average += d[i, j] * freq1[i] * freq1[j];
            }
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    d[i, j] -= average;
            average = 0.0;
            for (int i = 0; i < 20; i++)
                average += d[i, i] * freq1[i];
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    d[i, j] *= 600.0 / average;
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    d[i, j] -= Offset;
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    d[i, j] = Shishagonyuu(d[i, j]);
            var n = new int[26, 26];
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    n[i, j] = (int)d[i, j];
            return n;
        }

        /// <summary>The JTT / transmembrane branch of <c>constants()</c> (<c>--jtt</c>, <c>--tm</c>): PAM matrix of <paramref name="pamN"/> steps.</summary>
        private int[,] JttMatrix(IReadOnlyList<string> seq, int pamN, bool makeAverage0)
        {
            var rsr = new double[20, 20];
            double[] freq = Tm ? MafftTables2.TmFrequencies : MafftTables2.JttFrequencies;
            var r = MafftTables2.JttAcceptedMutations;
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    rsr[i, j] = r[i][j];
            for (int i = 0; i < 20; i++)
                rsr[i, i] = 0.0;
            if (Tm)
            {
                for (int i = 1; i < 20; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        rsr[j, i] /= 400.0 * freq[i] * freq[j];
                        rsr[i, j] = rsr[j, i];
                    }
                }
            }
            else
            {
                for (int i = 1; i < 20; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        rsr[i, j] /= 400.0 * freq[i] * freq[j];
                        rsr[j, i] = rsr[i, j];
                    }
                }
            }
            double[] freq1 = freq;
            if (Fmodel == 1)
            {
                var datafreq = new double[20];
                CalcFreq(seq, datafreq);
                freq1 = datafreq;
            }
            var mutab = new double[20];
            double tmp = 0.0;
            for (int i = 0; i < 20; i++)
            {
                mutab[i] = 0.0;
                for (int j = 0; j < 20; j++)
                    mutab[i] += rsr[i, j] * freq1[j];
                tmp += mutab[i] * freq1[i];
            }
            double delta = 0.01 / tmp;
            var pam1 = new double[20, 20];
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    pam1[i, j] = i != j ? delta * rsr[i, j] * freq1[j] : 1.0 - delta * mutab[i];
            var pamx = new double[20, 20];
            for (int i = 0; i < 20; i++)
                pamx[i, i] = 1.0;
            for (int x = 0; x < pamN; x++)
                MtxMltDouble(pamx, pam1, 20);
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    pamx[i, j] /= freq1[j];
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j < 20; j++)
                {
                    if (pamx[i, j] == 0.0)
                        pamx[i, j] = 0.00001;
                    pamx[i, j] = Math.Log10(pamx[i, j]) * 1000.0;
                }
            }
            double average = 0.0;
            if (Fmodel != -1)
            {
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        average += pamx[i, j] * freq1[i] * freq1[j];
            }
            if (makeAverage0)
            {
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        pamx[i, j] -= average;
            }
            average = 0.0;
            for (int i = 0; i < 20; i++)
                average += pamx[i, i] * freq1[i];
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    pamx[i, j] *= 600.0 / average;
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    pamx[i, j] -= Offset;
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    pamx[i, j] = Shishagonyuu(pamx[i, j]);
            var n = new int[26, 26];
            for (int i = 0; i < 20; i++)
                for (int j = 0; j < 20; j++)
                    n[i, j] = (int)pamx[i, j];
            return n;
        }
    }

    #endregion
}
