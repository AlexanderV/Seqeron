#pragma warning disable S127, S3358, S1066, S3776, S107, S1854, IDE0059, S2325
using System.Globalization;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --codonpos / --codonscore (Falign.c, Galign11.c G__align11psg, io.c loadcodonscore)

    /// <summary>The codon score matrix and the position-specific gap penalties of <c>G__align11psg</c> (<see cref="Matrix"/> null: <c>--codonpos</c> without <c>--codonscore</c>, the gap penalties only); <see cref="Offset"/> is the position of the segment in the reference.</summary>
    internal sealed record MafftCodonDp(double[][]? Matrix, double[] GStart, double[] GEnd, int Offset);

    /// <summary><c>codon2id</c> / <c>nuc2id</c> of <c>io.c</c>: t = 0, c = 1, a = 2, g = 3; -1 for a codon with another letter (also past the end of the sequence).</summary>
    internal static class MafftCodon
    {
        internal static int Nuc(char c) => c switch { 't' => 0, 'c' => 1, 'a' => 2, 'g' => 3, _ => -1 };

        internal static int Id(string s, int pos)
        {
            int id0 = pos < s.Length ? Nuc(s[pos]) : -1;
            int id1 = pos + 1 < s.Length ? Nuc(s[pos + 1]) : -1;
            int id2 = pos + 2 < s.Length ? Nuc(s[pos + 2]) : -1;
            if (id0 < 0 || id1 < 0 || id2 < 0)
                return -1;
            return id0 * 16 + id1 * 4 + id2;
        }

        private static string Letter(int id) => id switch { 0 => "t", 1 => "c", 2 => "a", _ => "g" };

        internal static string Name(int id) => Letter(id / 16) + Letter(id % 16 / 4) + Letter(id % 4);
    }

    internal sealed partial class MafftGlobals
    {
        /// <summary>The files <c>_codonpos</c> and <c>_codonscore</c> (the script writes them without empty lines and with a line break for every carriage return).</summary>
        internal string? CodonPosText, CodonScoreText;
    }

    internal sealed partial class MafftRun
    {
        private bool _codonReady;
        private double[]? _codonGStart, _codonGEnd;
        private double[][]? _codonMatrix;

        /// <summary>The first call of <c>Falign</c> reads the codon positions of the reference (the first group, one sequence) and the codon scores.</summary>
        private void InitCodon(string seq1)
        {
            int len1 = seq1.Length;
            const double strong = 3.0;
            _codonReady = true;
            _codonGStart = null;
            _codonGEnd = null;
            if (G.CodonPos)
            {
                if (G.Dorp != 'd')
                    throw new MafftProgramFailureException("\n\nThe --codonpos and --codonscore options are only for DNA data.\n\n");
                var gstart = new double[len1 + 1];
                var gend = new double[len1 + 1];
                for (int i = 0; i < len1 + 1; i++)
                {
                    gstart[i] = 1.0 * 0.5;
                    gend[i] = 1.0 * 0.5;
                }

                if (G.CodonPosText == null)
                    throw new MafftProgramFailureException("Cannot open _codonpos file\n");
                int n = 0;
                foreach (string line in SplitFgetsLines(G.CodonPosText, 1001))
                {
                    if (n == len1)
                        throw new MafftProgramFailureException("\n\nNumber of lines in the codonposition file must be the same as the length of the first sequences (" + len1.ToString(CultureInfo.InvariantCulture) + ").\n\n");
                    char c0 = line.Length > 0 ? line[0] : '\0';
                    if (c0 == '#')
                        continue;
                    if (c0 == '0')
                    {
                        gstart[n] = 1.0 * 0.5;
                        gend[n] = 1.0 * 0.5;
                    }
                    else if (c0 == '1')
                    {
                        gstart[n] = strong * 0.5;
                        gend[n] = 1.0 * 0.5;
                    }
                    else if (c0 == '2')
                    {
                        gstart[n] = strong * 0.5;
                        gend[n] = strong * 0.5;
                    }
                    else if (c0 == '3')
                    {
                        gstart[n] = 1.0 * 0.5;
                        gend[n] = strong * 0.5;
                    }
                    else
                    {
                        throw new MafftProgramFailureException("In the codonposition file, 1st letter in a line must be either of:\n\t0 (noncoding)\n\t1 (1st position in codon)\n\t2 (2nd position in codon)\n\t3 (3rd position in codon)\nWhen mutliple difference frames are used, set 2.\n");
                    }

                    n++;
                }

                if (n < len1)
                    throw new MafftProgramFailureException("\n\nNumber of lines in the codonposition file (" + n.ToString(CultureInfo.InvariantCulture) + ") is less than the length of the first sequences (" + len1.ToString(CultureInfo.InvariantCulture) + ").\n\n");
                _codonGStart = gstart;
                _codonGEnd = gend;
            }

            _codonMatrix = null;
            if (G.CodonScore)
            {
                if (!G.CodonPos)
                    throw new MafftProgramFailureException("\n\n --codonpos is necessary for --codonscore\n\n");
                var mtx = LoadCodonScore(G.CodonScoreText ?? throw new MafftProgramFailureException("Cannot open _codonscore file\n"));
                var codonfreq = new double[64];
                double totalcount = 0.0;
                for (int i = 3; i < len1; i++)
                {
                    if (_codonGStart![i - 2] == strong * 0.5 && _codonGEnd![i - 2] == 1.0 * 0.5 &&
                        _codonGStart[i - 1] == strong * 0.5 && _codonGEnd[i - 1] == strong * 0.5 &&
                        _codonGStart[i - 0] == 1.0 * 0.5 && _codonGEnd[i - 0] == strong * 0.5)
                    {
                        int id = MafftCodon.Id(seq1, i - 2);
                        codonfreq[id] += 1.0; // an id of -1 (a letter that is not t, c, a, g) is written before the vector: the program does so too
                        totalcount += 1.0;
                    }
                }

                for (int i = 0; i < 64; i++)
                    codonfreq[i] /= totalcount;
                double codonscore0av = 0.0;
                for (int i = 0; i < 64; i++)
                    for (int j = 0; j < 64; j++)
                        codonscore0av += codonfreq[i] * codonfreq[j] * mtx[i][j];
                for (int i = 0; i < 64; i++)
                    for (int j = 0; j < 64; j++)
                        mtx[i][j] -= codonscore0av;
                double codonscore1av = 0.0;
                for (int i = 0; i < 64; i++)
                    codonscore1av += codonfreq[i] * mtx[i][i];
                for (int i = 0; i < 64; i++)
                    for (int j = 0; j < 64; j++)
                        mtx[i][j] /= codonscore1av;
                _codonMatrix = mtx;
            }
        }

        /// <summary><c>loadcodonscore</c>: one line per pair of codons ("ttt F ttc F 6.0", the pair in both directions), the lines that start with '#' are comments; every pair of the 64 codons must be given.</summary>
        private static double[][] LoadCodonScore(string text)
        {
            var mtx = new double[64][];
            for (int i = 0; i < 64; i++)
            {
                mtx[i] = new double[64];
                Array.Fill(mtx[i], -99999.0);
            }

            foreach (string buf in SplitFgetsLines(text, 1001))
            {
                if (buf.Length > 0 && buf[0] == '#')
                    continue;
                if (buf.Length == 0 || buf[^1] != '\n')
                    throw new MafftProgramFailureException(buf + ": too long in codonscore file.\n");
                var sc = new CScanner(buf);
                double score = 0.0;
                if (!sc.Word(out string c1) || !sc.Word(out _) || !sc.Word(out string c2) || !sc.Word(out _) || !sc.Double(out score))
                    throw new NotSupportedException("A line of the codon score file has fewer than five fields: the program reads the variables of the line before.");
                string codon1 = c1.Length > 3 ? c1[..3] : c1;
                string codon2 = c2.Length > 3 ? c2[..3] : c2;
                int id1 = MafftCodon.Id(codon1, 0);
                int id2 = MafftCodon.Id(codon2, 0);
                if (id1 < 0 || id2 < 0)
                    throw new MafftProgramFailureException("Cannot use codon pair " + codon1 + " - " + codon2 + ": Use small letter, a, c, g, t (instead of u)\n");
                mtx[id1][id2] = mtx[id2][id1] = score;
            }

            for (int i = 0; i < 64; i++)
            {
                for (int j = 0; j < 64; j++)
                {
                    if (mtx[i][j] == -99999.0)
                        throw new MafftProgramFailureException("\nCodon score for " + MafftCodon.Name(i) + "-" + MafftCodon.Name(j) + " (id" + i.ToString(CultureInfo.InvariantCulture) + "-id" + j.ToString(CultureInfo.InvariantCulture) + ") is not given.\n");
                }
            }

            return mtx;
        }
    }

    #endregion
}
