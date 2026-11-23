#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905, S2681
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>Bio.Align.CodonAligner</c> of Biopython 1.88 (B09 WP67, F307, R14-10b; <c>_codonaligner.c</c>): aligns a nucleotide sequence to its protein by dynamic programming over the three translated frames of the nucleotide sequence, with scores for
/// the frame shifts. Every residue is aligned to one codon (three nucleotides); between two codons the nucleotide sequence may advance by 1, 2, 4 or 5 instead of 3 (a frame shift: -2, -1, +1, +2) at the score of the shift; leading and trailing
/// nucleotides are free. Residues and translated codons score <see cref="MatchScore"/> / <see cref="MismatchScore"/>, 0 against the <see cref="Wildcard"/> letter (X). The query coordinates of a frame shift go back (overlapping codons: -2 and -1)
/// or skip (+1 and +2) inside the alignment, as in Biopython.
/// </summary>
public sealed class CodonAligner
{
    private const int FrameshiftMinusTwo = 1, FrameshiftMinusOne = 2, FrameshiftNone = 4, FrameshiftPlusOne = 8, FrameshiftPlusTwo = 16;
    private readonly GeneticCode _code;

    /// <summary>Creates the aligner for a genetic code (NCBI table number, default 1: the standard code).</summary>
    /// <exception cref="ArgumentException">Unknown table.</exception>
    public CodonAligner(int tableNumber = 1)
    {
        _code = GeneticCode.GetByTableNumberOrThrow(tableNumber, nameof(tableNumber));
        TableNumber = tableNumber;
    }

    /// <summary>The NCBI number of the genetic code.</summary>
    public int TableNumber { get; }

    /// <summary><c>match_score</c> (default 1).</summary>
    public double MatchScore { get; set; } = 1.0;

    /// <summary><c>mismatch_score</c> (default 0).</summary>
    public double MismatchScore { get; set; }

    /// <summary><c>epsilon</c>: two scores closer than this are equal (default 1e-6).</summary>
    public double Epsilon { get; set; } = 1e-6;

    /// <summary><c>wildcard</c>: the letter that scores 0 against every letter (default <c>X</c>).</summary>
    public char Wildcard { get; set; } = 'X';

    /// <summary><c>frameshift_minus_two_score</c> (default -3): the next codon overlaps the previous one by two nucleotides.</summary>
    public double FrameshiftMinusTwoScore { get; set; } = -3.0;

    /// <summary><c>frameshift_minus_one_score</c> (default -3).</summary>
    public double FrameshiftMinusOneScore { get; set; } = -3.0;

    /// <summary><c>frameshift_plus_one_score</c> (default -3): one nucleotide is skipped between two codons.</summary>
    public double FrameshiftPlusOneScore { get; set; } = -3.0;

    /// <summary><c>frameshift_plus_two_score</c> (default -3).</summary>
    public double FrameshiftPlusTwoScore { get; set; } = -3.0;

    /// <summary><c>frameshift_minus_score</c>: the score of both negative shifts; reading it while they differ is an error, as in Biopython.</summary>
    public double FrameshiftMinusScore
    {
        get => Same(FrameshiftMinusTwoScore, FrameshiftMinusOneScore) ? FrameshiftMinusTwoScore : throw new InvalidOperationException("frame shift scores are different");
        set => FrameshiftMinusTwoScore = FrameshiftMinusOneScore = value;
    }

    /// <summary><c>frameshift_plus_score</c>.</summary>
    public double FrameshiftPlusScore
    {
        get => Same(FrameshiftPlusTwoScore, FrameshiftPlusOneScore) ? FrameshiftPlusTwoScore : throw new InvalidOperationException("frame shift scores are different");
        set => FrameshiftPlusTwoScore = FrameshiftPlusOneScore = value;
    }

    /// <summary><c>frameshift_two_score</c>: the score of the -2 and +2 shifts.</summary>
    public double FrameshiftTwoScore
    {
        get => Same(FrameshiftMinusTwoScore, FrameshiftPlusTwoScore) ? FrameshiftMinusTwoScore : throw new InvalidOperationException("frame shift scores are different");
        set => FrameshiftMinusTwoScore = FrameshiftPlusTwoScore = value;
    }

    /// <summary><c>frameshift_one_score</c>: the score of the -1 and +1 shifts.</summary>
    public double FrameshiftOneScore
    {
        get => Same(FrameshiftMinusOneScore, FrameshiftPlusOneScore) ? FrameshiftMinusOneScore : throw new InvalidOperationException("frame shift scores are different");
        set => FrameshiftMinusOneScore = FrameshiftPlusOneScore = value;
    }

    /// <summary><c>frameshift_score</c>: the score of all four shifts.</summary>
    public double FrameshiftScore
    {
        get => Same(FrameshiftMinusOneScore, FrameshiftMinusTwoScore) && Same(FrameshiftMinusOneScore, FrameshiftPlusOneScore) && Same(FrameshiftMinusOneScore, FrameshiftPlusTwoScore) ? FrameshiftMinusOneScore : throw new InvalidOperationException("frame shift scores are different");
        set => FrameshiftMinusTwoScore = FrameshiftMinusOneScore = FrameshiftPlusOneScore = FrameshiftPlusTwoScore = value;
    }

    private static bool Same(double a, double b) => a == b;

    // the three translated frames of the nucleotide sequence (seqB[0:], seqB[1:], seqB[2:], each cut to whole codons) as Bio.Seq.translate(str, CodonTable) does: DNA or RNA codons (not mixed) of the table by name,
    // stop codons as '*', every other codon of IUPAC letters as 'X' (no ambiguity resolution), anything else is an error
    private string[] Frames(string nucleotides)
    {
        var frames = new string[3];
        for (int f = 0; f < 3; f++)
        {
            int n = nucleotides.Length - f;
            int whole = n <= 0 ? 0 : 3 * (n / 3);
            var sb = new System.Text.StringBuilder(whole / 3);
            for (int k = f; k + 3 <= f + whole; k += 3)
                sb.Append(TranslateCodon(nucleotides.Substring(k, 3).ToUpperInvariant()));
            frames[f] = sb.ToString();
        }

        return frames;
    }

    private const string CodonLetters = "ACGTURYSWKMBDHVN";

    private char TranslateCodon(string codon)
    {
        if (!(codon.Contains('T') && codon.Contains('U')) && _code.CodonTable.TryGetValue(codon.Replace('T', 'U'), out char aa))
            return aa;
        foreach (char c in codon)
            if (!CodonLetters.Contains(c))
                throw new ArgumentException($"Codon '{codon}' is invalid");
        return 'X';
    }

    private double Compare(char a, char b) => a == Wildcard || b == Wildcard ? 0 : a == b ? MatchScore : MismatchScore;

    /// <summary><c>aligner.score(protein, nucleotides)</c>.</summary>
    public double Score(string protein, string nucleotides)
    {
        ArgumentNullException.ThrowIfNull(protein);
        ArgumentNullException.ThrowIfNull(nucleotides);
        var f = Frames(nucleotides);
        int nA = protein.Length, nB = ColumnCount(nucleotides.Length);
        var row = new double[nB + 1];
        for (int i = 1; i <= nA; i++)
        {
            char cA = protein[i - 1];
            for (int j = nB; j > 0; j--)
            {
                double score = -double.MaxValue;
                if (j >= 3)
                {
                    char cB = Codon(f, j);
                    double s = Compare(cA, cB);
                    score = Math.Max(score, row[j - 1] + s + FrameshiftMinusTwoScore);
                    score = Math.Max(score, row[j - 2] + s + FrameshiftMinusOneScore);
                    score = Math.Max(score, row[j - 3] + s);
                    if (j >= 4)
                        score = Math.Max(score, row[j - 4] + s + FrameshiftPlusOneScore);
                    if (j >= 5)
                        score = Math.Max(score, row[j - 5] + s + FrameshiftPlusTwoScore);
                }

                row[j] = score;
            }
        }

        return row.Max();
    }

    private static char Codon(string[] frames, int j)
    {
        int div = (j - 3) / 3, mod = (j - 3) % 3;
        return frames[mod][div];
    }

    /// <summary>
    /// <c>aligner.align(protein, nucleotides)</c>: every alignment of the optimal score as coordinates of the protein (target, residues) and the nucleotide sequence (query).
    /// </summary>
    public PairwiseAlignments Align(string protein, string nucleotides)
    {
        ArgumentNullException.ThrowIfNull(protein);
        ArgumentNullException.ThrowIfNull(nucleotides);
        var frames = Frames(nucleotides);
        int nA = protein.Length, nB = ColumnCount(nucleotides.Length);
        double eps = Epsilon;
        var trace = new int[nA + 1, nB + 1];
        var row = new double[nB + 1];
        for (int i = 1; i <= nA; i++)
        {
            char cA = protein[i - 1];
            for (int j = nB; j > 0; j--)
            {
                double score = -double.MaxValue;
                int t = 0;
                if (j >= 3)
                {
                    double s = Compare(cA, Codon(frames, j));
                    Select(row[j - 1] + s + FrameshiftMinusTwoScore, FrameshiftMinusTwo, ref score, ref t, eps);
                    Select(row[j - 2] + s + FrameshiftMinusOneScore, FrameshiftMinusOne, ref score, ref t, eps);
                    Select(row[j - 3] + s, FrameshiftNone, ref score, ref t, eps);
                    if (j >= 4)
                        Select(row[j - 4] + s + FrameshiftPlusOneScore, FrameshiftPlusOne, ref score, ref t, eps);
                    if (j >= 5)
                        Select(row[j - 5] + s + FrameshiftPlusTwoScore, FrameshiftPlusTwo, ref score, ref t, eps);
                }

                trace[i, j] = t;
                row[j] = score;
            }
        }

        double best = row.Max();
        // the trace of the last residue keeps the end points within epsilon of the best score; the others lose it
        var path = new int[nA + 1, nB + 1];
        for (int j = 0; j <= nB; j++)
        {
            if (row[j] < best - eps)
                trace[nA, j] = 0;
            else
                path[nA, j] = 0;
        }

        // the generator keeps its state in the path field of the cells; a cell the library never wrote holds whatever the heap held, so a read of one has no defined result
        var written = new bool[nA + 1, nB + 1];
        for (int j = 0; j <= nB; j++)
        {
            written[0, j] = true;
            if (row[j] >= best - eps)
                written[nA, j] = true;
        }

        // the library's length() counts the paths that start in the first row; its iterator walks the end points whatever that count is, so at least one step is taken
        int limit = (int)Math.Max(1, Math.Min(MaxAlignments, PathCount(trace, nA, nB)));
        var coordinates = new List<long[][]>();
        while (coordinates.Count < limit)
        {
            if (path[0, 0] == Done)
                break;
            int i = 0, j, p = 0, t;
            for (j = 0; j <= nB; j++)
            {
                p = path[0, j];
                if (p != 0)
                {
                    path[0, j] = 0;
                    while (true)
                    {
                        j += p;
                        t = trace[i + 1, j];
                        if (p == 1 && (t & FrameshiftMinusOne) != 0) { p = 2; break; }
                        if (p <= 2 && (t & FrameshiftNone) != 0) { p = 3; break; }
                        if (p <= 3 && (t & FrameshiftPlusOne) != 0) { p = 4; break; }
                        if (p <= 4 && (t & FrameshiftPlusTwo) != 0) { p = 5; break; }
                        i++;
                        p = path[i, j];
                        if (p == 0)
                            break;                                              // the end of the alignment without an alternative
                    }

                    if (p != 0)
                    {
                        j -= p;
                        path[i, j] = p;
                        written[i, j] = true;
                    }

                    break;
                }
            }

            if (p == 0)
            {
                if (i == 0)
                {
                    j = 0;
                    i = nA;
                }
                else
                {
                    j++;
                }

                for (; j <= nB; j++)
                    if (trace[nA, j] != 0)
                        break;
                if (j > nB)
                {
                    path[0, 0] = Done;
                    break;
                }
            }

            while (true)
            {
                t = trace[i, j];
                if ((t & FrameshiftMinusTwo) != 0) p = 1;
                else if ((t & FrameshiftMinusOne) != 0) p = 2;
                else if ((t & FrameshiftNone) != 0) p = 3;
                else if ((t & FrameshiftPlusOne) != 0) p = 4;
                else if ((t & FrameshiftPlusTwo) != 0) p = 5;
                else break;
                j -= p;
                i--;
                path[i, j] = p;
                written[i, j] = true;
            }

            coordinates.Add(CreatePath(path, written, j, nA, nB));
        }

        var target = new BioCoordSequence { Id = null, Length = nA, Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, protein)], nA) };
        var query = new BioCoordSequence { Id = null, Length = nB, Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, nucleotides)], nB) };
        var attributes = new System.Text.Json.Nodes.JsonObject();
        var alignments = coordinates.Select(c => new BioCoordAlignment { Sequences = [target, query], Coordinates = c, Attributes = attributes }).ToList();
        return new PairwiseAlignments(best, alignments);
    }

    private const int Done = 6;

    // PathGenerator_create_path: the vertices of the path that starts in the first row at the nucleotide position j. The library starts it in the first row even when the traceback ended in a later row (a free start after skipped residues);
    // it then follows cells the traceback never wrote, past the end of its arrays and of the vertex tuple it allocated, and has no defined result: refused.
    private static long[][] CreatePath(int[,] path, bool[,] written, int j0, int nA, int nB)
    {
        const string Undefined = "The reference library (Bio.Align.CodonAligner) reads memory it never initialised on this input: the optimal alignment skips leading residues, its traceback ends below the first row and the path is read from the first row; the result is undefined.";
        // the first pass counts the vertices the tuples are allocated for
        int n = 1, k0 = 0, l = j0, dir0 = 0;
        while (true)
        {
            if (k0 > nA || l > nB || l < 0 || !written[k0, l])
                throw new InvalidOperationException(Undefined);
            int p0 = path[k0, l];
            if (p0 == 0)
                break;
            if (p0 % 3 != 0)
            {
                n += 2;
                dir0 = 3;
            }
            else if (p0 != dir0)
            {
                n++;
                dir0 = p0;
            }

            l += p0;
            k0++;
        }

        var target = new List<long>(n);
        var query = new List<long>(n);
        int direction = 0;
        int i = 0, j = j0;
        while (true)
        {
            if (i > nA || j > nB || j < 0 || !written[i, j])
                throw new InvalidOperationException(Undefined);
            int p = path[i, j];
            if (p % 3 != 0)
            {
                target.Add(i);
                query.Add(j);
                j += p - 3;
                target.Add(i);
                query.Add(j);
                direction = 3;
            }
            else if (p != direction)
            {
                target.Add(i);
                query.Add(j);
                direction = p;
                if (p == 0 && target.Count <= n)
                    return [[.. target], [.. query]];
            }

            if (target.Count > n)
                throw new InvalidOperationException(Undefined);
            i++;
            j += 3;
        }
    }

    // PathGenerator_length: the number of paths from the first row to the end points (saturating at the maximum of the library, which raises OverflowError there)
    private static long PathCount(int[,] trace, int nA, int nB)
    {
        var counts1 = new long[nB + 1];
        var counts2 = new long[nB + 1];
        Array.Fill(counts2, 1L);
        for (int i = 1; i <= nA; i++)
        {
            Array.Copy(counts2, counts1, nB + 1);
            for (int j = 0; j <= nB; j++)
            {
                int t = trace[i, j];
                long count = 0;
                if ((t & FrameshiftMinusTwo) != 0) count = SaturatingAdd(count, counts1[j - 1]);
                if ((t & FrameshiftMinusOne) != 0) count = SaturatingAdd(count, counts1[j - 2]);
                if ((t & FrameshiftNone) != 0) count = SaturatingAdd(count, counts1[j - 3]);
                if ((t & FrameshiftPlusOne) != 0) count = SaturatingAdd(count, counts1[j - 4]);
                if ((t & FrameshiftPlusTwo) != 0) count = SaturatingAdd(count, counts1[j - 5]);
                counts2[j] = count;
            }
        }

        long total = 0;
        for (int j = 0; j <= nB; j++)
            total = SaturatingAdd(total, counts2[j]);
        return total;
    }

    private static long SaturatingAdd(long a, long b) => a > long.MaxValue - b ? long.MaxValue : a + b;

    // the number of nucleotide positions of the library: the three frames cut to whole codons give 3 * codons + 2 for a sequence shorter than three letters
    private static int ColumnCount(int n)
    {
        int l0 = n / 3, l1 = n >= 1 ? (n - 1) / 3 : 0, l2 = n >= 2 ? (n - 2) / 3 : 0;
        if (l1 == l0 && l2 == l0)
            return 3 * l0 + 2;
        if (l1 == l0 && l2 == l0 - 1)
            return 3 * l0 + 1;
        return 3 * l0;
    }

    private const int MaxAlignments = 1000;

    private static void Select(double temp, int flag, ref double score, ref int trace, double eps)
    {
        if (temp > score + eps)
        {
            score = temp;
            trace = flag;
        }
        else if (temp > score - eps)
        {
            trace |= flag;
        }
    }
}
