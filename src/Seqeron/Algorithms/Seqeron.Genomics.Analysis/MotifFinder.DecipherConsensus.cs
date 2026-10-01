using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>Sequence alphabet of the input to <see cref="MotifFinder.GenerateDecipherConsensus"/>.</summary>
public enum DecipherSequenceType
{
    /// <summary>DNA (Biostrings <c>DNAStringSet</c>, <c>DNA_ALPHABET</c>).</summary>
    Dna,

    /// <summary>RNA (Biostrings <c>RNAStringSet</c>, <c>RNA_ALPHABET</c>): <c>U</c> instead of <c>T</c>.</summary>
    Rna,

    /// <summary>Amino acids (Biostrings <c>AAStringSet</c>, <c>AA_ALPHABET</c>).</summary>
    AminoAcid,
}

/// <summary>
/// Bioconductor DECIPHER <c>ConsensusSequence</c> (Erik Wright): degeneracy-code consensus that
/// loses less than a <c>threshold</c> fraction of the sequence information at every position.
/// </summary>
public static partial class MotifFinder
{
    // Biostrings alphabets (DNA_ALPHABET, RNA_ALPHABET, AA_ALPHABET).
    private const string DecipherDnaAlphabet = "ACGTMRWSYKVHDBN-+.";
    private const string DecipherRnaAlphabet = "ACGUMRWSYKVHDBN-+.";
    private const string DecipherAaAlphabet = "ARNDCQEGHILKMFPSTWYVUOBJZX*-+.";

    // makeConsensusAA index order (AAs[0..22]): 20 canonical, U, O, '*'.
    private const string DecipherAaSymbols = "ARNDCQEGHILKMFPSTWYVUO*";

    // makeConsensus IUPAC tests in source order; bit 0 = A, 1 = C, 2 = G, 3 = T.
    private static readonly (char Symbol, int Mask)[] DecipherDnaTests =
    {
        ('A', 0b0001), ('C', 0b0010), ('G', 0b0100), ('T', 0b1000),
        ('Y', 0b1010), ('K', 0b1100), ('W', 0b1001), ('S', 0b0110), ('R', 0b0101), ('M', 0b0011),
        ('B', 0b1110), ('D', 0b1101), ('H', 0b1011), ('V', 0b0111),
    };

    /// <summary>
    /// Creates a consensus sequence with the rule of Bioconductor DECIPHER
    /// <c>ConsensusSequence(myXStringSet, threshold, ambiguity, noConsensusChar, minInformation,
    /// includeNonLetters, includeTerminalGaps)</c> — a line-by-line port of DECIPHER 3.9.4
    /// <c>R/ConsensusSequence.R</c> and <c>src/ConsensusSequence.c</c>
    /// (<c>alphabetFrequency</c>/<c>makeConsensus</c> and their <c>AA</c> forms).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Per column the characters are tallied as fractions of the counted characters (IUPAC codes
    /// split equally between their bases when <paramref name="ambiguity"/> is true, e.g. <c>R</c> =
    /// ½ A + ½ G, <c>X</c> = 1/20 of each canonical amino acid). With t = 1 − <paramref name="threshold"/>
    /// the first test that holds, in source order, chooses the symbol: a single residue that is
    /// strictly the most frequent with fraction ≥ t; then (DNA) the two- and three-base codes
    /// <c>Y K W S R M B D H V</c> whose bases are all strictly more frequent than the excluded
    /// bases and sum to ≥ t, then <c>N</c> if A+C+G+T ≥ t; (protein) <c>B</c>/<c>Z</c>/<c>J</c>
    /// when N/Q/I is the unique (or with D/E/L the tied) maximum and the pair sums to ≥ t, then
    /// <c>X</c> if all residues sum to ≥ t. The chosen fraction must also be ≥ the gap and mask
    /// fractions, otherwise <c>-</c> (or <c>+</c> when masks outnumber gaps) is emitted; a column
    /// failing every test emits <c>-</c>/<c>+</c> if the gap/mask fraction is ≥ t. When the
    /// fraction behind the emitted symbol is below <paramref name="minInformation"/> (or nothing
    /// was chosen), <paramref name="noConsensusChar"/> is emitted.
    /// </para>
    /// <para>
    /// Without <paramref name="includeTerminalGaps"/>, leading and trailing <c>-</c>/<c>.</c> of
    /// each sequence are not counted and a column with nothing counted is <c>-</c>; with it, such a
    /// column gives <paramref name="noConsensusChar"/>. <c>.</c> is a gap. Sequences may differ in
    /// length (the consensus has the length of the longest; shorter sequences do not count beyond
    /// their end), as in DECIPHER.
    /// </para>
    /// <para>
    /// <paramref name="includeNonLetters"/> is passed to the C argument <c>ignoreNonLetters</c>
    /// exactly as DECIPHER does: with the default <c>false</c> gaps (<c>-</c> <c>.</c>) and masks
    /// (<c>+</c>) are counted and compete with the residues; with <c>true</c> they are left out of
    /// the column (DECIPHER manual example: <c>c("A-+.A","AAAAA")</c> → <c>ANNNA</c> with
    /// <c>noConsensusChar="N"</c>, and <c>AAAAA</c> with <c>includeNonLetters=TRUE</c>). With
    /// <paramref name="ambiguity"/> = false the DNA/RNA non-letters are always counted and IUPAC
    /// codes are skipped; for proteins <c>B Z J X</c> are skipped.
    /// </para>
    /// <para>
    /// Input is case-insensitive (Biostrings upper-cases); the RNA consensus writes <c>U</c>.
    /// Cross-checked against DECIPHER's own C/R source compiled in R 4.3.3 with Biostrings 2.70
    /// (docs/Evidence/MOTIF-GENERATE-001-Evidence.md).
    /// </para>
    /// </remarks>
    /// <param name="sequences">Aligned sequences (DECIPHER also accepts unequal lengths).</param>
    /// <param name="sequenceType">DNA, RNA or amino-acid alphabet.</param>
    /// <param name="threshold">
    /// Fraction of sequence information that may be lost at a position, in [0, 1); DECIPHER
    /// default 0.05.
    /// </param>
    /// <param name="ambiguity">Split IUPAC degeneracy codes between their residues (default true).</param>
    /// <param name="noConsensusChar">
    /// Character for positions without consensus; must belong to the input alphabet (default <c>+</c>).
    /// </param>
    /// <param name="minInformation">
    /// Minimum fraction carried by the consensus character, in (0, 1]; default 1 − <paramref name="threshold"/>.
    /// </param>
    /// <param name="includeNonLetters">DECIPHER <c>includeNonLetters</c> (see remarks; default false).</param>
    /// <param name="includeTerminalGaps">Count leading/trailing gaps (default false).</param>
    /// <returns>The consensus (length of the longest sequence); "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// A null element, a character outside the alphabet, or <paramref name="noConsensusChar"/>
    /// outside the alphabet.
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="threshold"/> outside [0, 1), <paramref name="minInformation"/> outside (0, 1],
    /// NaN, or an undefined <paramref name="sequenceType"/>.
    /// </exception>
    public static string GenerateDecipherConsensus(
        IEnumerable<string> sequences,
        DecipherSequenceType sequenceType = DecipherSequenceType.Dna,
        double threshold = 0.05,
        bool ambiguity = true,
        char noConsensusChar = '+',
        double? minInformation = null,
        bool includeNonLetters = false,
        bool includeTerminalGaps = false)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (sequenceType is not (DecipherSequenceType.Dna or DecipherSequenceType.Rna or DecipherSequenceType.AminoAcid))
            throw new ArgumentOutOfRangeException(nameof(sequenceType));
        if (double.IsNaN(threshold) || threshold >= 1 || threshold < 0)
            throw new ArgumentOutOfRangeException(nameof(threshold), "Threshold must be in [0, 1).");
        double minInfo = minInformation ?? 1 - threshold;
        if (double.IsNaN(minInfo) || minInfo > 1 || minInfo <= 0)
            throw new ArgumentOutOfRangeException(nameof(minInformation), "minInformation must be in (0, 1].");

        string alphabet = sequenceType switch
        {
            DecipherSequenceType.Dna => DecipherDnaAlphabet,
            DecipherSequenceType.Rna => DecipherRnaAlphabet,
            _ => DecipherAaAlphabet,
        };
        if (alphabet.IndexOf(noConsensusChar) < 0)
            throw new ArgumentException(
                $"noConsensusChar '{noConsensusChar}' is not in the {sequenceType} alphabet.", nameof(noConsensusChar));

        // Biostrings encoding: upper-case, validate against the alphabet.
        var rows = new List<string>();
        int seqLength = 0;
        foreach (string s in sequences)
        {
            if (s is null)
                throw new ArgumentException("Sequences cannot contain null elements.", nameof(sequences));
            var upper = new char[s.Length];
            for (int k = 0; k < s.Length; k++)
            {
                char c = s[k] is >= 'a' and <= 'z' ? (char)(s[k] - 'a' + 'A') : s[k];
                if (alphabet.IndexOf(c) < 0)
                    throw new ArgumentException(
                        $"Invalid character '{s[k]}' at position {k} in sequence {rows.Count} for {sequenceType}.",
                        nameof(sequences));
                upper[k] = c;
            }
            rows.Add(new string(upper));
            seqLength = Math.Max(seqLength, s.Length);
        }

        bool ignore = includeNonLetters; // DECIPHER passes includeNonLetters as C ignoreNonLetters
        double t = 1 - threshold;
        string consensus = sequenceType == DecipherSequenceType.AminoAcid
            ? DecipherConsensusAa(rows, seqLength, ambiguity, ignore, includeTerminalGaps, t, minInfo)
            : DecipherConsensusDna(rows, seqLength, ambiguity, ignore, includeTerminalGaps, t, minInfo);

        if (sequenceType == DecipherSequenceType.Rna)
            consensus = consensus.Replace('T', 'U');
        return consensus.Replace('?', noConsensusChar);
    }

    private static bool IsDecipherGap(char c) => c is '-' or '.';

    /// <summary>frontTerminalGaps/endTerminalGaps: leading/trailing '-'/'.' counts.</summary>
    private static (int Start, int End) DecipherCountedRange(string row, bool includeTerminalGaps)
    {
        if (includeTerminalGaps)
            return (0, row.Length);
        int front = 0;
        while (front < row.Length && IsDecipherGap(row[front]))
            front++;
        int back = 0;
        while (back < row.Length && IsDecipherGap(row[row.Length - 1 - back]))
            back++;
        return (front, row.Length - back);
    }

    /// <summary>Biostrings DNA/RNA code of an IUPAC letter (A=1, C=2, G=4, T/U=8, OR for codes).</summary>
    private static int DecipherNucleotideMask(char c) => c switch
    {
        'A' => 1, 'C' => 2, 'M' => 3, 'G' => 4, 'R' => 5, 'S' => 6, 'V' => 7,
        'T' or 'U' => 8, 'W' => 9, 'Y' => 10, 'H' => 11, 'K' => 12, 'D' => 13, 'B' => 14, 'N' => 15,
        _ => 0,
    };

    private static string DecipherConsensusDna(
        List<string> rows, int seqLength, bool degeneracy, bool ignore, bool tGaps, double threshold, double minInfo)
    {
        // bits rows: 0..3 A C G T, 4 gap, 5 mask, 6 total.
        var bits = new double[7, seqLength];
        const double weight = 1;
        foreach (string row in rows)
        {
            var (start, end) = DecipherCountedRange(row, tGaps);
            for (int j = start; j < end; j++)
            {
                char p = row[j];
                if (degeneracy)
                {
                    bits[6, j] += weight;
                    if (p is '-' or '.' or '+')
                    {
                        if (ignore)
                            bits[6, j] -= weight;
                        else
                            bits[p == '+' ? 5 : 4, j] += weight;
                        continue;
                    }

                    int mask = DecipherNucleotideMask(p);
                    int n = System.Numerics.BitOperations.PopCount((uint)mask);
                    // .5*weight, (double)1/3*weight, .25*weight as in alphabetFrequency.
                    double share = n switch { 1 => weight, 2 => .5 * weight, 3 => (double)1 / 3 * weight, _ => .25 * weight };
                    for (int b = 0; b < 4; b++)
                    {
                        if ((mask & (1 << b)) != 0)
                            bits[b, j] += share;
                    }
                }
                else
                {
                    int row6 = p switch
                    {
                        'A' => 0, 'C' => 1, 'G' => 2, 'T' or 'U' => 3, '-' or '.' => 4, '+' => 5,
                        _ => -1, // degeneracy code: not counted
                    };
                    if (row6 >= 0)
                    {
                        bits[row6, j] += weight;
                        bits[6, j] += weight;
                    }
                }
            }
        }

        var seq = new StringBuilder(seqLength);
        var pct = new double[4];
        for (int j = 0; j < seqLength; j++)
        {
            double total = bits[6, j];
            if (!tGaps && total == 0)
            {
                seq.Append('-');
                continue;
            }

            for (int b = 0; b < 4; b++)
                pct[b] = bits[b, j] / total;
            double percentGap = bits[4, j] / total;
            double percentMask = bits[5, j] / total;
            double information = 0;
            char symbol = '?';
            bool chosen = false;

            foreach (var (code, mask) in DecipherDnaTests)
            {
                double sum = 0;
                bool dominant = true;
                for (int b = 0; b < 4; b++)
                {
                    if ((mask & (1 << b)) == 0)
                        continue;
                    sum += pct[b];
                    for (int o = 0; o < 4; o++)
                    {
                        if ((mask & (1 << o)) == 0 && !(pct[b] > pct[o]))
                            dominant = false;
                    }
                }
                if (sum >= threshold && dominant)
                {
                    (symbol, information) = DecipherPick(code, sum, percentGap, percentMask);
                    chosen = true;
                    break;
                }
            }

            if (!chosen)
            {
                double all = pct[0] + pct[1] + pct[2] + pct[3];
                if (all >= threshold)
                    (symbol, information) = DecipherPick('N', all, percentGap, percentMask);
                else if (percentGap >= threshold || percentMask >= threshold)
                    (symbol, information) = percentGap >= percentMask ? ('-', percentGap) : ('+', percentMask);
            }

            if (information < minInfo)
                symbol = '?';
            seq.Append(symbol);
        }

        return seq.ToString();
    }

    /// <summary>The emission block shared by every makeConsensus branch.</summary>
    private static (char Symbol, double Information) DecipherPick(
        char symbol, double fraction, double percentGap, double percentMask)
    {
        if (fraction >= percentGap && fraction >= percentMask)
            return (symbol, fraction);
        return percentGap >= percentMask ? ('-', percentGap) : ('+', percentMask);
    }

    private static string DecipherConsensusAa(
        List<string> rows, int seqLength, bool degeneracy, bool ignore, bool tGaps, double threshold, double minInfo)
    {
        // bits rows: 0..19 canonical, 20 U, 21 O, 22 '*', 23 gap, 24 mask, 25 total.
        var bits = new double[26, seqLength];
        const double weight = 1;
        foreach (string row in rows)
        {
            var (start, end) = DecipherCountedRange(row, tGaps);
            for (int j = start; j < end; j++)
            {
                char p = row[j];
                bits[25, j] += weight;
                int index = DecipherAaSymbols.IndexOf(p);
                if (index >= 0)
                {
                    bits[index, j] += weight;
                    continue;
                }

                switch (p)
                {
                    case 'B': DecipherSplitPair(bits, 2, 3, j, degeneracy, weight); break;
                    case 'Z': DecipherSplitPair(bits, 5, 6, j, degeneracy, weight); break;
                    case 'J': DecipherSplitPair(bits, 9, 10, j, degeneracy, weight); break;
                    case 'X':
                        if (degeneracy)
                        {
                            for (int i = 0; i < 20; i++)
                                bits[i, j] += weight / 20;
                        }
                        else
                        {
                            bits[25, j] -= weight;
                        }
                        break;
                    default: // '-', '.', '+'
                        if (ignore)
                            bits[25, j] -= weight;
                        else
                            bits[p == '+' ? 24 : 23, j] += weight;
                        break;
                }
            }
        }

        var seq = new StringBuilder(seqLength);
        var aas = new double[23];
        for (int j = 0; j < seqLength; j++)
        {
            double total = bits[25, j];
            if (!tGaps && total == 0)
            {
                seq.Append('-');
                continue;
            }

            double sumAll = 0;
            int m = 0;
            int tied = 0;
            for (int i = 0; i < 23; i++)
            {
                aas[i] = bits[i, j] / total;
                sumAll += aas[i];
                if (aas[i] > aas[m])
                {
                    tied = 1;
                    m = i;
                }
                else if (aas[i] == aas[m])
                {
                    tied++;
                }
            }

            double percentGap = bits[23, j] / total;
            double percentMask = bits[24, j] / total;
            double information = 0;
            char symbol = '?';

            if (tied == 1 && aas[m] >= threshold)
                (symbol, information) = DecipherPick(DecipherAaSymbols[m], aas[m], percentGap, percentMask);
            else if (DecipherAaPair(aas, m, tied, 2, 3, threshold))
                (symbol, information) = DecipherPick('B', aas[2] + aas[3], percentGap, percentMask);
            else if (DecipherAaPair(aas, m, tied, 5, 6, threshold))
                (symbol, information) = DecipherPick('Z', aas[5] + aas[6], percentGap, percentMask);
            else if (DecipherAaPair(aas, m, tied, 9, 10, threshold))
                (symbol, information) = DecipherPick('J', aas[9] + aas[10], percentGap, percentMask);
            else if (sumAll >= threshold)
                (symbol, information) = DecipherPick('X', sumAll, percentGap, percentMask);
            else if (percentGap >= threshold || percentMask >= threshold)
                (symbol, information) = percentGap >= percentMask ? ('-', percentGap) : ('+', percentMask);

            if (information < minInfo)
                symbol = '?';
            seq.Append(symbol);
        }

        return seq.ToString();
    }

    private static void DecipherSplitPair(double[,] bits, int a, int b, int j, bool degeneracy, double weight)
    {
        if (degeneracy)
        {
            bits[a, j] += 0.5 * weight;
            bits[b, j] += 0.5 * weight;
        }
        else
        {
            bits[25, j] -= weight;
        }
    }

    /// <summary>makeConsensusAA B/Z/J test: M == a and (unique max, or tied only with b) and a+b ≥ t.</summary>
    private static bool DecipherAaPair(double[] aas, int m, int tied, int a, int b, double threshold) =>
        m == a && ((tied == 2 && aas[a] == aas[b]) || tied == 1) && aas[a] + aas[b] >= threshold;
}
