#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --text / --textmatrix (constants.c: the extended alphabet of 256 characters, blosum.c: extendedmtx, overridematrix, calcfreq_extended)

    internal sealed partial class MafftGlobals
    {
        /// <summary>The kind of data that the version line of the programs shows (<c>SHOWVERSION</c> of <c>mltaln.h</c>): nuc, text (<c>nblosum == -2</c>) or aa.</summary>
        internal string KindName
        {
            get
            {
                if (Dorp == 'd')
                    return "nuc";
                return Nblosum == -2 ? "text" : "aa";
            }
        }

        /// <summary>The text mode (<c>-b -2</c>, <c>nblosum == -2</c>): every byte except gap, blank and line ends is a letter of the 256-letter alphabet.</summary>
        internal bool TextMode => Dorp == 'p' && Scoremtx == 1 && Nblosum == -2;

        /// <summary><c>extendedmtx</c> (<c>blosum.c</c>) and the extended branch of <c>constants()</c>: the matrix is -1 for a mismatch, +1 for a match, then the entries of the file <c>_aamtx</c> (<c>overridematrix</c>) replace them; the average
        /// over the frequencies (the letters of the data with <c>-a</c>, otherwise 1/256) is removed, the matrix is scaled so that the average of the diagonal is 600, the offset is taken off and the values are rounded.</summary>
        private int[,] ExtendedMatrix(IReadOnlyList<string> seq, out double[] freq)
        {
            Nalphabets = 0x100;
            NscoredAlphabets = 0x100;
            Report("nalphabets = " + Nalphabets.ToString(CultureInfo.InvariantCulture) + "\n");
            int na = Nalphabets;
            if (Ppenalty == NotSpecified) Ppenalty = -1530;
            if (PenaltyDistIsUnset()) PpenaltyDist = Ppenalty;
            if (PpenaltyEx == NotSpecified) PpenaltyEx = 0;
            if (PpenaltyOp == NotSpecified) PpenaltyOp = -1530;
            if (PpenaltyExGenaf == NotSpecified) PpenaltyExGenaf = 0;
            if (Poffset == NotSpecified) Poffset = -123;
            if (PamN == NotSpecified) PamN = 0;
            if (KimuraR == NotSpecified) KimuraR = 1;
            SetProteinPenalties();
            // extendedmtx
            for (int i = 0; i < na; i++)
                Amino[i] = (char)i;
            for (int i = 0; i < na; i++)
                AminoGrp[Amino[i]] = i % 6;
            var datafreqExtended = new double[na];
            var distmp = new double[na, na];
            for (int i = 0; i < na; i++)
                for (int j = 0; j <= i; j++)
                    distmp[i, j] = distmp[j, i] = -1.0;
            for (int i = 0; i < na; i++)
                distmp[i, i] = 1.0;
            OverrideMatrix(distmp);
            ModelName = "Extended, " + F(-(double)Ppenalty / 1000, 4, 2) + ", " + F(-(double)Poffset / 1000, 4, 2, true) + ", " + F(-(double)PpenaltyEx / 1000, 4, 2, true) + ", "
                + (TryWarp ? F(-(double)PenaltyShift / 600, 4, 2) : "noshift");
            for (int i = 0; i < 0x100; i++)
                AminoN[i] = -1;
            for (int i = 0; i < na; i++)
                AminoN[Amino[i]] = i;
            var freq0 = new double[na];
            for (int i = 0; i < na; i++)
                freq0[i] = 1.0 / na;
            double[] freq1;
            if (Fmodel == 1)
            {
                CalcFreqExtended(seq, datafreqExtended);
                freq1 = datafreqExtended;
            }
            else
            {
                freq1 = freq0;
            }
            double average;
            if (Fmodel == -1)
            {
                average = 0.0;
            }
            else
            {
                average = 0.0;
                for (int i = 0; i < na; i++)
                    for (int j = 0; j < na; j++)
                        average += distmp[i, j] * freq1[i] * freq1[j];
            }
            for (int i = 0; i < na; i++)
                for (int j = 0; j < na; j++)
                    distmp[i, j] -= average;
            average = 0.0;
            for (int i = 0; i < na; i++)
                average += distmp[i, i] * freq1[i];
            if (average < 0.0)
                throw new MafftProgramFailureException("\nUnrealistic scoring matrix.  Give larger positive values to matches (A/A, B/B, etc).\n\n");
            for (int i = 0; i < na; i++)
                for (int j = 0; j < na; j++)
                    distmp[i, j] *= 600.0 / average;
            for (int i = 0; i < na; i++)
                for (int j = 0; j < na; j++)
                    distmp[i, j] -= Offset;
            for (int i = 0; i < na; i++)
                for (int j = 0; j < na; j++)
                    distmp[i, j] = Shishagonyuu(distmp[i, j]);
            var nDis = new int[na, na];
            for (int i = 0; i < na; i++)
                for (int j = 0; j < na; j++)
                    nDis[i, j] = CInt(distmp[i, j]);
            int gap = AminoN['-'];
            for (int i = 0; i < na; i++)
                for (int j = 0; j < na; j++)
                    nDis[i, gap] = nDis[gap, i] = 0;
            freq = freq1;
            return nDis;
        }

        /// <summary><c>calcfreq_extended</c>: the frequencies of the letters of the data (gaps not counted), normalised to sum 1 (NaN when there is no letter).</summary>
        private void CalcFreqExtended(IReadOnlyList<string> seq, double[] datafreq)
        {
            int n = NscoredAlphabets;
            double total = 0.0;
            foreach (var s in seq)
            {
                foreach (char c in s)
                {
                    int aan = c < 256 ? AminoN[c] : -1;
                    if (aan >= 0 && aan < n && c != '-')
                    {
                        datafreq[aan] += 1.0;
                        total += 1.0;
                    }
                }
            }
            total = 0.0;
            for (int i = 0; i < n; i++)
                total += datafreq[i];
            for (int i = 0; i < n; i++)
                datafreq[i] /= total;
            total = 0.0;
            for (int i = 0; i < n; i++)
                total += datafreq[i];
            for (int i = 0; i < n; i++)
                datafreq[i] /= total;
        }

        /// <summary>
        /// <c>overridematrix</c> (<c>blosum.c</c>): the lines of the file <c>_aamtx</c> ("0x41 0x42 3.5", a comment after '#') set the score of a pair of letters (both orders); a line that is not three numbers stops the program, a letter that
        /// cannot be a letter (NUL, blank, '-', '&lt;', '=', '&gt;', CR, LF, above 0xff) as well. The text file is read with <c>fgets( buf, 499 )</c>; the last line counts only when it ends with a line end (the script's <c>grep</c> makes sure it does).
        /// </summary>
        private void OverrideMatrix(double[,] matrix)
        {
            string? text = AaMatrixText;
            if (text == null)
                return;
            int pos = 0;
            while (pos < text.Length)
            {
                int nl = text.IndexOf('\n', pos);
                int end = nl < 0 ? text.Length : nl + 1;
                int take = Math.Min(end - pos, 498);
                string buf = text.Substring(pos, take);
                pos += take;
                if (pos >= text.Length && !buf.EndsWith('\n'))
                    break; // feof after the read: the program leaves the loop before it looks at the line
                int hash = buf.IndexOf('#');
                if (hash >= 0)
                    buf = buf[..hash];
                var sc = new MafftRun.CScanner(buf);
                int nread = 0;
                uint i1 = 0, i2 = 0;
                double v = 0.0;
                bool eof = false;
                if (sc.Hex(out i1))
                {
                    nread = 1;
                    if (sc.Hex(out i2))
                    {
                        nread = 2;
                        if (sc.Double(out v))
                            nread = 3;
                    }
                }
                else if (sc.AtEnd)
                {
                    eof = true;
                }
                if (eof)
                    continue;
                if (nread != 3)
                    throw new MafftProgramFailureException("Format error in this line?\n" + buf + "\n");
                if (CheckChar(i1))
                    throw new MafftProgramFailureException((char)(i1 & 0xFF) + "=0x" + i1.ToString("x", CultureInfo.InvariantCulture) + " cannot be used (1)\n");
                if (CheckChar(i2))
                    throw new MafftProgramFailureException((char)(i2 & 0xFF) + "=0x" + i2.ToString("x", CultureInfo.InvariantCulture) + " cannot be used (2)\n");
                Report("Score(" + (char)(i1 & 0xFF) + "=0x" + i1.ToString("x", CultureInfo.InvariantCulture) + "," + (char)(i2 & 0xFF) + "=0x" + i2.ToString("x", CultureInfo.InvariantCulture) + ")=" + F(v, 0, 6) + "\n");
                matrix[i1, i2] = v;
                matrix[i2, i1] = v;
            }
        }

        private static bool CheckChar(uint i) => i > 0xff || i <= 0x00 || i == 0x3E || i == 0x3D || i == 0x3C || i == 0x2D || i == 0x20 || i == 0x0d || i == 0x0a;
    }

    internal sealed partial class MafftRun
    {
        /// <summary><c>charfilter</c> + the loader of the text mode (<c>load1SeqWithoutName_realloc</c> with <c>nblosum == -2</c>): every byte stays except blank, CR and LF; '=', '&lt;' and '&gt;' stop the program.</summary>
        internal static string TextResidues(string raw)
        {
            var sb = new StringBuilder(raw.Length);
            foreach (char ch in raw)
            {
                if (ch >= 256)
                    throw new ArgumentException("MAFFT works on bytes; the sequence contains a character above U+00FF.", nameof(raw));
                if (ch is '=' or '<' or '>')
                    throw new MafftProgramFailureException("\nCharacters '= < >' can be used only in the title lines in the --anysymbol or --text mode.\n\n");
                if (ch is not (' ' or '\r' or '\n'))
                    sb.Append(ch);
            }
            return sb.ToString();
        }
    }

    #endregion
}
