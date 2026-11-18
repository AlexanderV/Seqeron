using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region ClustalW PROTEIN pairwise stages (-TYPE=PROTEIN: full DP with GONNET / BLOSUM / PAM / ID matrices, fast k-tuple stage)

    // ClustalW 2.1 pairwise stages for protein sequences (audit A1-15b part b): the SAME members as the DNA stages
    // (FullPairwiseAlign: forwardPass / reversePass / diff / tracePath, shared through RunClustalWPairwise; FastPairwiseAlign: pairAlign, shared through
    // ClustalWQuickFragmentScore) with the protein parameters of SubMatrix.cpp::getMatrix / getPairwiseMatrix, UserParameters.cpp and FullPairwiseAlign.cpp:
    //
    //   * matrix: the lower-triangular table of matrices.h in the packing order "ABCDEFGHIKLMNPQRSTVWXYZ" (gon250mt, blosum30mt, pam350mt, idmat =
    //     SubstitutionMatrix.ClustalWGonnet250 / ClustalWBlosum30 / ClustalWPam350 / ClustalWIdentity) scattered to the residue CODES of
    //     UserParameters::aminoAcidCodes ("ABCDEFGHIKLMNOPQRSTUVWXYZ-"; O and U have no table row and score 0 against everything) and multiplied by
    //     intScale (100; GONNET 10 because gon250mt is in tenth-bits);
    //   * matrixAvgScore = -(int)(S / ((float)(maxRes^2 - maxRes) / 2)) with maxRes = 22 (the 23 table residues minus one) and S the sum of all
    //     off-diagonal cells (i != j, j < i) over the 25 residue codes, in single precision;
    //   * gap open  = 2 * (int)((PWGAPOPEN + ln MIN(len1, len2)) * intScale)                         when matrixAvgScore <= 0 (identity matrix),
    //                 (int)(2 * matrixAvgScore * (PWGAPOPEN + ln MIN(len1, len2)) * gapOpenScale)    otherwise  (gapOpenScale = 1),
    //     gap extend = (int)(PWGAPEXT * intScale) in single precision; a gap of length k costs open + k * extend as in the DNA stage;
    //   * fast stage: residues are their codes 0..24 (k-tuple code = sum residue * 25^(j-1), 25^ktup tuples), defaults KTUPLE 1 / WINDOW 5 / TOPDIAGS 5 /
    //     PAIRGAP 3 (the command line sets PAIRGAP = KTUPLE + 3 for protein), KTUPLE <= 2.
    //
    // Oracle: the installed `clustalw -TYPE=PROTEIN -ALIGN` / `-QUICKTREE -ALIGN` (see generate_clustalw_protein_pairwise_oracle.py).

    /// <summary>ClustalW 2.1 PROTEIN default pairwise gap-opening penalty (<c>-PWGAPOPEN</c>).</summary>
    public const double ClustalWProteinPairwiseGapOpen = 10.0;

    /// <summary>ClustalW 2.1 PROTEIN default pairwise gap-extension penalty (<c>-PWGAPEXT</c>).</summary>
    public const double ClustalWProteinPairwiseGapExtend = 0.1;

    private const int ClustalWMaxAminoAcidCode = 24;   // UserParameters::maxAA = aminoAcidCodes.Length - 2

    private static readonly Lazy<(int[,] Matrix, int IntScale, int MatrixAvgScore)>[] ClustalWProteinMatrices =
    {
        new(() => BuildClustalWProteinMatrix(SubstitutionMatrix.ClustalWGonnet250, 10)),
        new(() => BuildClustalWProteinMatrix(SubstitutionMatrix.ClustalWBlosum30, 100)),
        new(() => BuildClustalWProteinMatrix(SubstitutionMatrix.ClustalWPam350, 100)),
        new(() => BuildClustalWProteinMatrix(SubstitutionMatrix.ClustalWIdentity, 100)),
    };

    /// <summary><c>SubMatrix::getMatrix</c> (pairwise, <c>negFlag</c> = true) for a table of <c>matrices.h</c>: the matrix over residue codes and <c>matrixAvgScore</c>.</summary>
    private static (int[,] Matrix, int IntScale, int MatrixAvgScore) BuildClustalWProteinMatrix(SubstitutionMatrix table, int intScale)
    {
        int n = ClustalWResidueCodes.Length;
        var matrix = new int[n, n];
        string order = SubstitutionMatrix.ClustalWAminoAcidOrder;
        int maxRes = 0;
        for (int i = 0; i < order.Length; i++)
        {
            int ti = ClustalWResidueCodes.IndexOf(order[i]);
            for (int j = 0; j <= i; j++)
            {
                int tj = ClustalWResidueCodes.IndexOf(order[j]);
                int k = table.GetScore(order[i], order[j]) * intScale;
                matrix[ti, tj] = k;
                matrix[tj, ti] = k;
                if (i == j)
                    maxRes++;
            }
        }

        maxRes--;
        int av3 = 0;
        for (int i = 0; i <= ClustalWMaxAminoAcidCode; i++)
            for (int j = 0; j < i; j++)
                av3 += matrix[i, j];

        av3 = (int)((float)av3 / ((float)(maxRes * maxRes - maxRes) / 2));
        return (matrix, intScale, -av3);
    }

    /// <summary>ClustalW residue codes of a protein sequence; a symbol ClustalW does not read (J, '*', non-letters) is rejected (the program silently drops it).</summary>
    private static int[] EncodeClustalWProtein(string sequence, string paramName)
    {
        var codes = new int[sequence.Length];
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            int code = c is >= 'A' and <= 'Z' && c != 'J' ? ClustalWResidueCodes.IndexOf(c) : -1;
            if (code < 0)
                throw new ArgumentException($"'{sequence[i]}' at position {i} is not a residue code known to ClustalW's protein stage (A-Z except J).", paramName);
            codes[i] = code;
        }

        return codes;
    }

    /// <summary>
    /// ClustalW 2.1 pairwise stage for PROTEIN (<c>clustalw -TYPE=PROTEIN -ALIGN</c>; default <c>-PWMATRIX=GONNET -PWGAPOPEN=10 -PWGAPEXT=0.1</c>) with the percent
    /// identity and guide-tree distance ClustalW derives from it: the same affine local (Smith–Waterman–Gotoh / Myers–Miller) port as the DNA overload
    /// with ClustalW's protein matrix (× intScale), its length-dependent gap-opening penalty and its single-precision gap extension (see the region comment).
    /// </summary>
    /// <remarks>Verified against the installed <c>clustalw</c> 2.1 (<c>Sequences (i:j) Aligned. Score:</c>) on the pairs of <c>clustalw_protein_pairwise_oracle.tsv</c>
    /// (four matrices, gap variants, B/Z/X, low-complexity and repeat pairs).</remarks>
    /// <param name="sequence1">First sequence (letters A-Z except J; case-insensitive; O and U score 0 as in ClustalW).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="options">Matrix and pairwise gap penalties (default: GONNET, 10, 0.1).</param>
    /// <exception cref="ArgumentNullException">A sequence or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException">A symbol is not a residue ClustalW reads, or a penalty is negative / NaN.</exception>
    public static ClustalWPairwiseResult ClustalWPairwiseAlign(string sequence1, string sequence2, ClustalWProteinPairwiseOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(options);
        ValidateClustalWProteinOptions(options);

        var codes1 = EncodeClustalWProtein(sequence1, nameof(sequence1));
        var codes2 = EncodeClustalWProtein(sequence2, nameof(sequence2));

        // FullPairwiseAlign: an empty sequence gives distance 1.0 (percent identity 0) without any alignment.
        if (codes1.Length == 0 || codes2.Length == 0)
            return new ClustalWPairwiseResult(0, 0, 0.0, 1.0, 0, "", "", -1, -1, -1, -1);

        var (matrix, intScale, matAvg) = ClustalWProteinMatrices[(int)options.Matrix].Value;

        float pwGapOpen = (float)options.GapOpen;
        float pwGapExtend = (float)options.GapExtend;
        double logMin = Math.Log((double)Math.Min(codes1.Length, codes2.Length));
        int gapOpenInt = matAvg <= 0
            ? 2 * (int)(((double)pwGapOpen + logMin) * intScale)
            : (int)(2 * matAvg * ((double)pwGapOpen + logMin) * 1.0);
        float extend = pwGapExtend * (float)intScale;
        int gapExtendInt = (int)extend;

        return RunClustalWPairwise(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), codes1, codes2, matrix, gapOpenInt, gapExtendInt);
    }

    /// <summary><see cref="ClustalWPairwiseAlign(string, string, ClustalWProteinPairwiseOptions)"/> for <see cref="ProteinSequence"/> inputs.</summary>
    public static ClustalWPairwiseResult ClustalWPairwiseAlign(ProteinSequence sequence1, ProteinSequence sequence2, ClustalWProteinPairwiseOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ClustalWPairwiseAlign(sequence1.Sequence, sequence2.Sequence, options);
    }

    private static void ValidateClustalWProteinOptions(ClustalWProteinPairwiseOptions o)
    {
        if (!Enum.IsDefined(o.Matrix))
            throw new ArgumentOutOfRangeException(nameof(o), o.Matrix, "Unknown ClustalW pairwise matrix.");
        if (o.GapOpen < 0 || o.GapExtend < 0 || double.IsNaN(o.GapOpen) || double.IsNaN(o.GapExtend))
            throw new ArgumentException("ClustalW pairwise gap penalties are non-negative numbers.", nameof(o));
        ArgumentNullException.ThrowIfNull(o.Quick);
    }

    /// <summary>
    /// ClustalW's PROTEIN fast pairwise stage (<c>-TYPE=PROTEIN -QUICKTREE</c>): the port of <c>FastPairwiseAlign::pairAlign</c> with residue codes as k-tuple letters.
    /// Same score / distance conversion as <see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/>; defaults KTUPLE 1, WINDOW 5, TOPDIAGS 5, PAIRGAP 3.
    /// </summary>
    /// <param name="sequence1">First sequence (letters A-Z except J).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="type">Residue type; <see cref="ClustalWSequenceType.Dna"/> is the existing DNA stage.</param>
    /// <param name="options">k-tuple parameters (default: the defaults of <paramref name="type"/>).</param>
    /// <param name="maxSequenceLength">Longest sequence of the whole input set (see the DNA overload); <c>null</c> = the longer of the two.</param>
    public static ClustalWQuickPairwiseResult ClustalWQuickPairwiseAlign(
        string sequence1, string sequence2, ClustalWSequenceType type, ClustalWQuickOptions? options = null, int? maxSequenceLength = null)
    {
        if (!Enum.IsDefined(type))
            throw new ArgumentOutOfRangeException(nameof(type), type, "Unknown ClustalW sequence type.");
        if (type == ClustalWSequenceType.Dna)
            return ClustalWQuickPairwiseAlign(sequence1, sequence2, options, maxSequenceLength);

        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ClustalWQuickPairwiseAlignCore(
            EncodeClustalWProtein(sequence1, nameof(sequence1)), EncodeClustalWProtein(sequence2, nameof(sequence2)),
            protein: true, options ?? new ClustalWProteinPairwiseOptions().Quick, maxSequenceLength);
    }

    /// <summary><see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWSequenceType, ClustalWQuickOptions?, int?)"/> for protein sequences.</summary>
    public static ClustalWQuickPairwiseResult ClustalWQuickPairwiseAlign(
        ProteinSequence sequence1, ProteinSequence sequence2, ClustalWQuickOptions? options = null, int? maxSequenceLength = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ClustalWQuickPairwiseAlign(sequence1.Sequence, sequence2.Sequence, ClustalWSequenceType.Protein, options, maxSequenceLength);
    }

    /// <summary>
    /// Single-precision guide-tree distance of the full stage for a protein pair (as <c>ClustalWSinglePrecisionDistance</c> for DNA): <c>mmScore = 100f·count;
    /// mmScore /= MIN(len1, len2); d = (100f − mmScore) / 100f</c> (an empty sequence: 1).
    /// </summary>
    private static double ClustalWProteinSinglePrecisionDistance(string a, string b, ClustalWProteinPairwiseOptions options)
    {
        if (a.Length == 0 || b.Length == 0)
            return 1.0;

        var r = ClustalWPairwiseAlign(a, b, options);
        float mm = 100.0f * r.Identities;
        mm /= Math.Min(a.Length, b.Length);
        return (100.0f - mm) / 100.0f;
    }

    #endregion
}
