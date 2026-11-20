using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Options of EMBOSS <c>distmat</c> (<see cref="MultipleAlignmentStatistics.CalculateDistanceMatrix(IReadOnlyList{string}, MsaDistanceOptions)"/>), qualifier for qualifier:
/// <c>-nucmethod</c> / <c>-protmethod</c> = <see cref="Method"/>, <c>-gapweight</c>, <c>-ambiguous</c>, <c>-position</c>, <c>-calculatea</c>, <c>-parametera</c>.
/// </summary>
public sealed record MsaDistanceOptions
{
    /// <summary>The correction (0-5 as the distmat menus: 0 uncorrected, 1 Jukes-Cantor, 2 Kimura = 2-parameter for nucleotides / "Kimura Protein" for proteins, 3 Tamura, 4 Tajima-Nei, 5 Jin-Nei Gamma; 3-5 are nucleotide only).</summary>
    public MsaDistanceMethod Method { get; init; } = MsaDistanceMethod.Uncorrected;

    /// <summary>Residue type; <see cref="ConsensusResidueType.Auto"/> = EMBOSS' type of the first row.</summary>
    public ConsensusResidueType ResidueType { get; init; } = ConsensusResidueType.Auto;

    /// <summary><c>-gapweight</c>: weight of the gap columns in the uncorrected and Jukes-Cantor methods (default 0).</summary>
    public double GapWeight { get; init; }

    /// <summary>
    /// <c>-ambiguous</c>: count ambiguity codes in the uncorrected / Jukes-Cantor / Tajima-Nei match count. For proteins <c>distmat_checkambigProt</c> (identical residues 1, B vs D/N, Z vs E/Q, J vs I/L 0.5,
    /// X vs X 0.0025, X vs anything 0.05); for nucleotides EMBOSS' <c>ajBaseAlphacharCompare</c>, whose own defects are reproduced (it examines all but the last alternative of the first code, so
    /// an unambiguous base never matches, and it does not terminate when the first code is not an IUB letter - gap columns - which is refused with <see cref="NotSupportedException"/>).
    /// </summary>
    public bool Ambiguous { get; init; }

    /// <summary><c>-position</c>: codon positions analysed: 123 (all, default), 12, or for nucleotides 1, 2, 3; anything else is EMBOSS' fatal "Choose base positions 1, 2, 3, 12, or 123".</summary>
    public int Position { get; init; } = 123;

    /// <summary><c>-calculatea</c>: force the calculation of the Jin-Nei parameter a (default: use <see cref="ParameterA"/>).</summary>
    public bool CalculateA { get; init; }

    /// <summary><c>-parametera</c>: the Jin-Nei Gamma parameter a (default 1.0).</summary>
    public double ParameterA { get; init; } = 1.0;
}

public static partial class MultipleAlignmentStatistics
{
    // ---------------------------------------------------------------------------------------------------------------------------------
    // EMBOSS 6.6.0 emboss/distmat.c, ported function by function with the C arithmetic types kept: `float` variables are single precision,
    // `log` / `pow` / `sqrt` take and return double, an integer or `float` operand of a double expression is promoted. The output file prints the
    // matrix x 100 with two decimals; the values returned here are the unscaled float results.

    /// <summary>Is <paramref name="c"/> found by C's <c>strchr(set, c)</c>? A NUL character matches (strchr finds the terminator).</summary>
    private static bool StrChr(string set, char c) => c == '\0' || set.Contains(c);

    private static char At(string s, int k) => k < s.Length ? s[k] : '\0';

    /// <summary>
    /// Pairwise distance matrix of aligned rows exactly as EMBOSS <c>distmat</c> 6.6.0 computes it (every <c>-nucmethod</c> / <c>-protmethod</c>, <c>-gapweight</c>, <c>-ambiguous</c>, <c>-position</c>,
    /// <c>-calculatea</c>, <c>-parametera</c>): the fraction D of the printed percentage (the program prints 100 D with two decimals), single-precision arithmetic as in the C source, symmetric with a
    /// zero diagonal. Non-finite entries (nan / inf) are what the program prints.
    /// </summary>
    /// <exception cref="ArgumentException">Fewer than two rows, a position the program rejects, a nucleotide-only method for proteins.</exception>
    /// <exception cref="NotSupportedException">The reference's behaviour is undefined or it does not terminate: Tajima-Nei and the protein Kimura method with a position other than 123 read beyond the end of
    /// their buffers; <see cref="MsaDistanceOptions.Ambiguous"/> on nucleotides with a gap in the first sequence of a pair never ends.</exception>
    public static double[,] CalculateDistanceMatrix(IReadOnlyList<string> alignedRows, MsaDistanceOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        var rows = Normalize(alignedRows);
        if (rows.Length < 2)
            throw new ArgumentException($"Insufficient sequences ({rows.Length}) to create a matrix; EMBOSS distmat requires at least 2.", nameof(alignedRows));
        if (double.IsNaN(options.GapWeight) || double.IsInfinity(options.GapWeight))
            throw new ArgumentOutOfRangeException(nameof(options), "GapWeight must be finite.");
        if (!Enum.IsDefined(options.Method))
            throw new ArgumentOutOfRangeException(nameof(options), "Unknown distance method.");
        bool nuc = IsNucleotide(rows, options.ResidueType);
        if (!nuc && options.Method is MsaDistanceMethod.Tamura or MsaDistanceMethod.TajimaNei or MsaDistanceMethod.JinNei)
            throw new ArgumentException($"{options.Method} is a nucleotide method (distmat -nucmethod); proteins have 0 (uncorrected), 1 (Jukes-Cantor) and 2 (Kimura Protein).", nameof(options));

        // The EMBOSS reader types every row of an alignment on its own (Auto) and a nucleotide row has X turned into N (ajSeqSetNuc); with an explicit residue type every row is of that type.
        for (int i = 0; i < rows.Length; i++)
        {
            bool rowNuc = options.ResidueType == ConsensusResidueType.Auto ? MotifFinder.IsEmbossNucleotideSequence(rows[i]) : nuc;
            if (rowNuc)
                rows[i] = rows[i].Replace('X', 'N');
        }

        // main(): codon positions
        int incr = 1, posn = options.Position;
        if (posn >= 1 && posn <= 3 && nuc)
        {
            incr = 3;
            posn--;
        }
        else if (posn == 123)
            posn = 0;
        else if (posn == 23 || posn == 13 || posn != 12)
            throw new ArgumentException("Choose base positions 1, 2, 3, 12, or 123.", nameof(options));

        int nseqs = rows.Length, mlen = rows[0].Length;
        var seq = GetSeq(rows, mlen, incr, posn, out int len);
        // distmat_TajimaNei and distmat_KimuraProt are called with mlen (the alignment length), not len: with a position other than 123 they read mlen - len characters beyond the end of the
        // position-selected sequence buffers, i.e. heap memory (the program's result depends on the allocator's layout: about half of the oracle's runs equal the NUL-padded model, the others do not)
        if (len < mlen && (options.Method == MsaDistanceMethod.TajimaNei || (options.Method == MsaDistanceMethod.Kimura && !nuc)))
            throw new NotSupportedException($"EMBOSS distmat -{(nuc ? "nuc" : "prot")}method {(int)options.Method} with -position {options.Position} reads beyond the end of its sequence buffers (undefined behaviour of the reference); use position 123.");
        float gapwt = (float)options.GapWeight;
        var method = options.Method;
        float[,] matchF = null!, gapF = null!;
        if (method is MsaDistanceMethod.Uncorrected or MsaDistanceMethod.JukesCantor or MsaDistanceMethod.TajimaNei)
            (matchF, gapF) = CalcMatch(seq, len, nseqs, options.Ambiguous, nuc);

        float[,] dist = method switch
        {
            MsaDistanceMethod.Uncorrected => DistUncorrected(matchF, gapF, gapwt, len, nseqs),
            MsaDistanceMethod.JukesCantor => DistJukesCantor(matchF, gapF, gapwt, len, nseqs, nuc),
            MsaDistanceMethod.Kimura => nuc ? DistKimura(seq, len, nseqs) : DistKimuraProt(seq, mlen, nseqs),
            MsaDistanceMethod.Tamura => DistTamura(seq, len, nseqs),
            MsaDistanceMethod.TajimaNei => DistTajimaNei(seq, matchF, mlen, nseqs),
            _ => DistJinNei(seq, len, nseqs, options.CalculateA, (float)options.ParameterA),
        };
        var d = new double[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                d[i, j] = dist[i, j];
                d[j, i] = dist[i, j];
            }

        return d;
    }

    // distmat_getSeq: the codon positions that enter the distances (position 12 takes the first two bases of each codon; a missing character is the NUL that strncpy copies from the terminator)
    private static string[] GetSeq(string[] rows, int mlen, int incr, int posn, out int len)
    {
        bool onetwo = false;
        if (posn == 12)
        {
            onetwo = true;
            posn = 0;
            incr = 3;
        }

        len = 0;
        for (int j = posn; j < mlen; j += incr)
        {
            len++;
            if (onetwo)
                len++;
        }

        var result = new string[rows.Length];
        for (int i = 0; i < rows.Length; i++)
        {
            var sb = new System.Text.StringBuilder(len);
            for (int j = posn; j < mlen; j += incr)
            {
                sb.Append(At(rows[i], j));
                if (onetwo)
                    sb.Append(At(rows[i], j + 1));
            }

            result[i] = sb.ToString();
        }

        return result;
    }

    // distmat_calc_match
    private static (float[,] Match, float[,] Gap) CalcMatch(string[] seq, int len, int nseqs, bool ambig, bool nuc)
    {
        var match = new float[nseqs, nseqs];
        var gap = new float[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per sequence (the pair loops below are O(length) each)
            for (int j = i + 1; j < nseqs; j++)
                for (int k = 0; k < len; k++)
                {
                    char m1 = At(seq[i], k), m2 = At(seq[j], k);
                    if (ambig && nuc)
                        match[i, j] = match[i, j] + BaseAlphacharCompare(m1, m2);
                    else if (ambig && !nuc && m1 != '-' && m2 != '-')
                        match[i, j] = match[i, j] + CheckAmbigProt(m1, m2);
                    else if (m2 == m1 && m1 != '-')
                        match[i, j] = match[i, j] + 1;
                    if (m1 == '-' || m2 == '-')
                        gap[i, j] = gap[i, j] + 1;
                }
        }

        return (match, gap);
    }

    // ajBaseAlphacharCompare (ajax/core/ajbase.c) with the IUB table of Ebases.iub: n = (1/len1)(1/len2) when any of the first len1 = |list1| - 1 letters of code 1 is in list 2, with the program's
    // own length arithmetic (len = |list| - 1: an unambiguous base has len 0 and examines nothing, 1/0 = +inf in float)
    private static float BaseAlphacharCompare(char c, char c2)
    {
        string l1 = IubList(c), l2 = IubList(c2);
        int len1 = l1.Length - 1, len2 = l2.Length - 1;
        if (len1 < 0)
            throw new NotSupportedException($"EMBOSS distmat -ambiguous on nucleotides does not terminate on the character '{c}' (it loops over (unsigned)(-1) alternatives): not an IUB nucleotide code.");
        bool pmatch = false;
        for (int i = 0; i < len1; i++)
            if (l2.Contains(l1[i]))
                pmatch = true;
        return pmatch ? (1f / len1) * (1f / len2) : 0f;
    }

    private static string IubList(char c) => c switch
    {
        'A' => "A", 'B' => "CGT", 'C' => "C", 'D' => "AGT", 'G' => "G", 'H' => "ACT", 'K' => "GT", 'M' => "AC", 'N' => "ACGT", 'R' => "AG", 'S' => "CG",
        'T' => "T", 'U' => "T", 'V' => "ACG", 'W' => "AT", 'X' => "ACGT", 'Y' => "CT",
        _ => "",
    };

    // distmat_uncorrected
    private static float[,] DistUncorrected(float[,] match, float[,] gap, float gapwt, int len, int nseqs)
    {
        var d = new float[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                float m = match[i, j], g = gap[i, j];
                d[i, j] = 1 - (m / ((float)len - g + (g * gapwt)));
            }

        return d;
    }

    // distmat_JukesCantor
    private static float[,] DistJukesCantor(float[,] match, float[,] gap, float gapwt, int mlen, int nseqs, bool nuc)
    {
        float b = (float)(19.0 / 20.0);
        if (nuc)
            b = (float)(3.0 / 4.0);
        var d = new float[nseqs, nseqs];
        for (int j = 0; j < nseqs; j++)
            for (int i = j; i < nseqs; i++)
            {
                float m = match[j, i], g = gap[j, i];
                float dd = 1 - (m / ((float)mlen - g + (g * gapwt)));
                float val = (float)Math.Log(1.0 - (double)(dd / b));
                val *= -b;
                d[j, i] = val;
            }

        return d;
    }

    // distmat_checkSubs / distmat_checkRY
    private static void CheckSubs(char t1, char t2, ref int trans, ref int tranv)
    {
        if (StrChr("A", t1))
        {
            if (StrChr("G", t2)) ++trans;
            else if (StrChr("CTU", t2)) ++tranv;
        }
        else if (StrChr("G", t1))
        {
            if (StrChr("A", t2)) ++trans;
            else if (StrChr("CTU", t2)) ++tranv;
        }
        else if (StrChr("C", t1))
        {
            if (StrChr("TU", t2)) ++trans;
            else if (StrChr("AG", t2)) ++tranv;
        }
        else if (StrChr("T", t1))
        {
            if (StrChr("CU", t2)) ++trans;
            else if (StrChr("AG", t2)) ++tranv;
        }
        else if (StrChr("U", t1))
        {
            if (StrChr("TC", t2)) ++trans;
            else if (StrChr("AG", t2)) ++tranv;
        }
    }

    private static void CheckRY(char t1, char t2, ref int trans, ref int tranv)
    {
        if (StrChr("R", t1))
        {
            if (StrChr("AGR", t2)) ++trans;
            else if (StrChr("CTUY", t2)) ++tranv;
        }
        else if (StrChr("AG", t1))
        {
            if (StrChr("R", t2)) ++trans;
            else if (StrChr("Y", t2)) ++tranv;
        }
        else if (StrChr("Y", t1))
        {
            if (StrChr("CTUY", t2)) ++trans;
            else if (StrChr("AGR", t2)) ++tranv;
        }
        else if (StrChr("CTU", t1))
        {
            if (StrChr("Y", t2)) ++trans;
            else if (StrChr("R", t2)) ++tranv;
        }
    }

    // distmat_checkambigProt
    private static float CheckAmbigProt(char t1, char t2)
    {
        float n = 0f;
        if (!StrChr("X", t1) && t1 == t2)
            n = 1.0f;
        else if ((StrChr("B", t1) && StrChr("DN", t2)) || (StrChr("B", t2) && StrChr("DN", t1)))
            n = 0.5f;
        else if ((StrChr("Z", t1) && StrChr("EQ", t2)) || (StrChr("Z", t2) && StrChr("EQ", t1)))
            n = 0.5f;
        else if ((StrChr("J", t1) && StrChr("IL", t2)) || (StrChr("J", t2) && StrChr("IL", t1)))
            n = 0.5f;
        else if (StrChr("X", t1) && StrChr("X", t2))
            n = 0.0025f;
        else if (StrChr("X", t1) || StrChr("X", t2))
            n = 0.05f;
        return n;
    }

    // distmat_Kimura (2-parameter, nucleotides)
    private static float[,] DistKimura(string[] seq, int len, int nseqs)
    {
        var d = new float[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                int ptrans = 0, qtranv = 0, match = 0;
                for (int k = 0; k < len; k++)
                {
                    char t1 = char.ToUpperInvariant(At(seq[i], k)), t2 = char.ToUpperInvariant(At(seq[j], k));
                    if (!StrChr("-", t2) && !StrChr("-", t1))
                    {
                        int trans = ptrans, tranv = qtranv;
                        CheckSubs(t1, t2, ref trans, ref tranv);
                        CheckRY(t1, t2, ref trans, ref tranv);
                        ptrans = trans;
                        qtranv = tranv;
                        match++;
                    }
                }

                float p = (float)ptrans / (float)match;
                float q = (float)qtranv / (float)match;
                d[i, j] = (float)(-0.5 * Math.Log((1 - (2 * p) - q) * Math.Sqrt(1 - (2 * q))));
            }

        return d;
    }

    // distmat_KimuraProt
    private static float[,] DistKimuraProt(string[] seq, int mlen, int nseqs)
    {
        var d = new float[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                float match = 0;
                int scored = 0;
                for (int k = 0; k < mlen; k++)
                {
                    char m1 = At(seq[i], k), m2 = At(seq[j], k);
                    if (m1 != '-' && m2 != '-')
                    {
                        match = match + CheckAmbigProt(char.ToUpperInvariant(m1), char.ToUpperInvariant(m2));
                        scored++;
                    }
                }

                float dd = (float)(1.0 - (match / (float)scored));
                dd = (float)(-Math.Log(1 - dd - (0.2 * dd * dd)));
                d[i, j] = dd;
            }

        return d;
    }

    // distmat_Tamura
    private static float[,] DistTamura(string[] seq, int len, int nseqs)
    {
        var gc = new float[nseqs, nseqs];
        var gap = new float[nseqs, nseqs];
        var ptrans = new int[nseqs, nseqs];
        var qtranv = new int[nseqs, nseqs];
        var score = new int[nseqs, nseqs];
        var d = new float[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
            for (int j = 0; j < nseqs; j++)
                for (int k = 0; k < len; k++)
                {
                    float cgap = gap[i, j], geecee = gc[i, j];
                    char t1 = char.ToUpperInvariant(At(seq[j], k)), t2 = char.ToUpperInvariant(At(seq[i], k));
                    if (StrChr("-NXWRYMKBVDH", t1) || StrChr("-NXWRYMKBVDH", t2))
                        ++cgap;
                    else if (StrChr("GCS", t2))
                        ++geecee;
                    gc[i, j] = geecee;
                    gap[i, j] = cgap;
                }

        for (int i = 0; i < nseqs; i++)
            for (int j = 0; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                float cgap = gap[i, j];
                gc[i, j] = gc[i, j] / ((float)len - cgap);
            }

        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
                for (int k = 0; k < len; k++)
                {
                    char t1 = char.ToUpperInvariant(At(seq[i], k)), t2 = char.ToUpperInvariant(At(seq[j], k));
                    if (!StrChr("-", t2) && !StrChr("-", t1))
                    {
                        int trans = ptrans[i, j], tranv = qtranv[i, j];
                        CheckSubs(t1, t2, ref trans, ref tranv);
                        ptrans[i, j] = trans;
                        qtranv[i, j] = tranv;
                        score[i, j] = score[i, j] + 1;
                    }
                }

        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                float c = gc[j, i] + gc[i, j];
                c = c - (2 * gc[j, i] * gc[i, j]);
                float p = (float)ptrans[i, j] / (float)score[i, j];
                float q = (float)qtranv[i, j] / (float)score[i, j];
                if (p != 0f)
                    p = p / c;
                d[i, j] = (float)(-(c * Math.Log(1 - p - q)) - (0.5 * (1 - c) * Math.Log(1 - 2 * q)));
            }

        return d;
    }

    // distmat_TajimaNei: counts over mlen characters (the C code passes mlen, not len, so a position other than 123 reads NUL padding beyond the buffer)
    private static float[,] DistTajimaNei(string[] seq, float[,] match, int mlen, int nseqs)
    {
        var len = new int[nseqs, nseqs];
        var pfreq = new int[nseqs, nseqs, 6];
        var cbase = new int[nseqs, nseqs, 4];
        var d = new float[nseqs, nseqs];
        for (int i = 0; i < nseqs; i++)
            for (int j = 0; j < nseqs; j++)
                for (int l = 0; l < mlen; l++)
                {
                    char ti = char.ToUpperInvariant(At(seq[i], l));
                    if (!StrChr("-NXWMKBVDH", ti))
                    {
                        char tj = char.ToUpperInvariant(At(seq[j], l));
                        if (!StrChr("-NXWMKBVDH", tj))
                        {
                            len[i, j] = len[i, j] + 1;
                            if (StrChr("G", ti)) cbase[i, j, 3]++;
                            else if (StrChr("C", ti)) cbase[i, j, 2]++;
                            else if (StrChr("T", ti)) cbase[i, j, 1]++;
                            else if (StrChr("A", ti)) cbase[i, j, 0]++;
                        }

                        if (ti != tj)
                        {
                            if ((StrChr("A", ti) && StrChr("T", tj)) || (StrChr("T", ti) && StrChr("A", tj)))
                                pfreq[i, j, 0]++;
                            else if ((StrChr("A", ti) && StrChr("C", tj)) || (StrChr("C", ti) && StrChr("A", tj)))
                                pfreq[i, j, 1]++;
                            else if ((StrChr("A", ti) && StrChr("G", tj)) || (StrChr("G", ti) && StrChr("A", tj)))
                                pfreq[i, j, 2]++;
                            else if ((StrChr("T", ti) && StrChr("C", tj)) || (StrChr("C", ti) && StrChr("T", tj)))
                                pfreq[i, j, 3]++;
                            else if ((StrChr("T", ti) && StrChr("G", tj)) || (StrChr("G", ti) && StrChr("T", tj)))
                                pfreq[i, j, 4]++;
                            else if ((StrChr("C", ti) && StrChr("G", tj)) || (StrChr("G", ti) && StrChr("C", tj)))
                                pfreq[i, j, 5]++;
                        }
                    }
                }

        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                float slen = (float)len[i, j];
                float fij2 = 0f;
                for (int bs = 0; bs < 4; bs++)
                {
                    float fi = (float)cbase[i, j, bs], fj = (float)cbase[j, i, bs];
                    float fij = 0f;
                    if (fi != 0f && fj != 0f)
                        fij = (fi + fj) / ((float)2.0 * slen);
                    fij2 += fij * fij;
                }

                int pair = 0;
                float h = 0f;
                for (int bs = 0; bs < 3; bs++)
                    for (int bs1 = bs + 1; bs1 < 4; bs1++)
                    {
                        float x1 = (float)pfreq[i, j, pair];
                        float fij = x1 / slen;
                        float ci1 = (float)cbase[j, i, bs], cj1 = (float)cbase[i, j, bs];
                        float ci2 = (float)cbase[j, i, bs1], cj2 = (float)cbase[i, j, bs1];
                        if (fij != 0f)
                            h += ((float)0.5 * fij * fij) / ((ci1 + cj1) / ((float)2.0 * slen) * (ci2 + cj2) / ((float)2.0 * slen));
                        pair++;
                    }

                float m = match[i, j];
                float dd = (float)1.0 - m / slen;
                float b = (float)0.5 * (1 - fij2 + ((dd * dd) / h));
                d[i, j] = (float)(-b * Math.Log(1.0 - (dd / b)));
            }

        return d;
    }

    // distmat_JinNei (the variable slen keeps its value from the previous pair when a pair has no scored position, as in the C source)
    private static float[,] DistJinNei(string[] seq, int mlen, int nseqs, bool calcA, float varA)
    {
        var ptrans = new int[nseqs, nseqs];
        var qtranv = new int[nseqs, nseqs];
        var lenArr = new int[nseqs, nseqs];
        var avL = new float[nseqs, nseqs];
        var cval = new float[nseqs, nseqs];
        var d = new float[nseqs, nseqs];
        int slen = 0;
        float xlen;
        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                float av = avL[i, j];
                for (int k = 0; k < mlen; k++)
                {
                    char t1 = char.ToUpperInvariant(At(seq[i], k)), t2 = char.ToUpperInvariant(At(seq[j], k));
                    if (!StrChr("-", t1) && !StrChr("-", t2))
                    {
                        slen = lenArr[i, j] + 1;
                        lenArr[i, j] = slen;
                        int trans = 0, tranv = 0;
                        CheckSubs(t1, t2, ref trans, ref tranv);
                        CheckRY(t1, t2, ref trans, ref tranv);
                        av += ((float)trans + ((float)2.0 * (float)tranv));
                        trans += ptrans[i, j];
                        tranv += qtranv[i, j];
                        ptrans[i, j] = trans;
                        qtranv[i, j] = tranv;
                    }
                }

                xlen = (float)slen;
                avL[i, j] = av / xlen;
            }

        if (calcA)
        {
            for (int i = 0; i < nseqs; i++)
                for (int j = i + 1; j < nseqs; j++)
                {
                    AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                    xlen = (float)lenArr[i, j];
                    float av = avL[i, j];
                    float variance = 0f;
                    for (int k = 0; k < mlen; k++)
                    {
                        char t1 = char.ToUpperInvariant(At(seq[i], k)), t2 = char.ToUpperInvariant(At(seq[j], k));
                        if (!StrChr("-", t1) && !StrChr("-", t2))
                        {
                            int trans = 0, tranv = 0;
                            CheckSubs(t1, t2, ref trans, ref tranv);
                            CheckRY(t1, t2, ref trans, ref tranv);
                            variance += (av - (float)(trans + (2 * tranv))) * (av - (float)(trans + (2 * tranv)));
                        }
                    }

                    variance = variance / xlen;
                    cval[i, j] = (av * av) / variance;
                }
        }

        for (int i = 0; i < nseqs; i++)
            for (int j = i + 1; j < nseqs; j++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per pair of sequences
                xlen = (float)lenArr[i, j];
                float p = (float)ptrans[i, j] / xlen;
                float q = (float)qtranv[i, j] / xlen;
                float aPar = calcA ? cval[i, j] : varA;
                d[i, j] = (float)(0.5 * aPar * (Math.Pow(1.0 - (2 * p) - q, -1.0 / aPar) + (0.5 * Math.Pow(1.0 - (2 * q), -1.0 / aPar)) - 1.5));
            }

        return d;
    }
}
