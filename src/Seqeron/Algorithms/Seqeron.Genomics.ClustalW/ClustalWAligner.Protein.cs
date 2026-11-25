// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

public static partial class ClustalWAligner
{
    #region ClustalW PROTEIN pairwise stages (-TYPE=PROTEIN: full DP with GONNET / BLOSUM / PAM / ID matrices, fast k-tuple stage)

    // ClustalW 2.1 pairwise stages for protein sequences (audit A1-15b part b): the SAME members as the DNA stages
    // (FullPairwiseAlign: forwardPass / reversePass / diff / tracePath, shared through RunClustalWPairwise; FastPairwiseAlign: pairAlign, shared through
    // ClustalWQuickFragmentScore) with the protein parameters of SubMatrix.cpp::getMatrix / getPairwiseMatrix, UserParameters.cpp and FullPairwiseAlign.cpp:
    //
    //   * matrix: the lower-triangular table of matrices.h in the packing order "ABCDEFGHIKLMNPQRSTVWXYZ" (gon250mt, blosum30mt, pam350mt, idmat =
    //     ClustalWMatrices.ClustalWGonnet250 / ClustalWBlosum30 / ClustalWPam350 / ClustalWIdentity) scattered to the residue CODES of
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
        new(() => BuildClustalWProteinMatrix(ClustalWMatrices.ClustalWGonnet250, 10)),
        new(() => BuildClustalWProteinMatrix(ClustalWMatrices.ClustalWBlosum30, 100)),
        new(() => BuildClustalWProteinMatrix(ClustalWMatrices.ClustalWPam350, 100)),
        new(() => BuildClustalWProteinMatrix(ClustalWMatrices.ClustalWIdentity, 100)),
    };

    // C++ static_cast<int>(floating point) on x86-64 (cvttsd2si / cvttss2si): NaN and out-of-range values give INT_MIN
    private static int ClustalWCastToInt(double v) => double.IsNaN(v) || v >= 2147483648.0 || v < -2147483648.0 ? int.MinValue : (int)v;

    /// <summary><c>defaultAAXref</c>: for each of the 23 columns of the packing order of <c>matrices.h</c> its residue code (-1 beyond).</summary>
    private static readonly int[] ClustalWDefaultAAXref =
        ClustalWMatrices.ClustalWAminoAcidOrder.Select(c => ClustalWResidueCodes.IndexOf(c)).ToArray();

    /// <summary>The stored lower triangle of a table of <c>matrices.h</c> (x <paramref name="multiplier"/>) as the cells <c>getMatrix</c> consumes.</summary>
    private static IEnumerable<(int Row, int Column, short Value)> ClustalWTableCells(SubstitutionMatrix table, int multiplier)
    {
        string order = ClustalWMatrices.ClustalWAminoAcidOrder;
        for (int i = 0; i < order.Length; i++)
            for (int j = 0; j <= i; j++)
                yield return (i, j, (short)(table.GetScore(order[i], order[j]) * multiplier));
    }

    /// <summary><c>SubMatrix::getMatrix</c> up to the average score (the part shared by the pairwise and the profile matrices) for a table of <c>matrices.h</c>:
    /// the matrix over residue codes (cell = table x <paramref name="multiplier"/> x <paramref name="intScale"/>) and <c>matrixAvgScore</c>.</summary>
    private static (int[,] Matrix, int IntScale, int MatrixAvgScore) BuildClustalWProteinMatrix(SubstitutionMatrix table, int intScale, int multiplier = 1)
    {
        var (matrix, avg) = ClustalWGetMatrix(ClustalWTableCells(table, multiplier), ClustalWDefaultAAXref, intScale);
        return (matrix, intScale, avg);
    }

    /// <summary>
    /// <c>SubMatrix::getMatrix</c> before the shift: the stored cells (row / column = positions of the file header or of the packing order, in the order
    /// <c>getMatrix</c> consumes them) are scattered through <paramref name="cross"/> (<c>Xref</c>: position -&gt; residue code) into the symmetric
    /// <see cref="CwNumRes"/> x <see cref="CwNumRes"/> matrix over residue codes (cell x <paramref name="intScale"/>; unlisted residues score 0), and
    /// <c>matrixAvgScore = -(int)(S / ((float)(maxRes^2 - maxRes) / 2))</c> with maxRes = listed residues - 1 and S the sum of the off-diagonal cells (j &lt; i)
    /// over the 25 residue codes, in single precision.
    /// </summary>
    private static (int[,] Matrix, int MatrixAvgScore) ClustalWGetMatrix(IEnumerable<(int Row, int Column, short Value)> cells, IReadOnlyList<int> cross, int intScale)
    {
        var matrix = new int[CwNumRes, CwNumRes];
        int maxRes = 0;
        foreach (var (i, j, value) in cells)
        {
            int ti = cross[i], tj = cross[j];
            int k = value * intScale;
            if (ti == tj)
            {
                matrix[ti, ti] = k;
                maxRes++;
            }
            else
            {
                matrix[ti, tj] = k;
                matrix[tj, ti] = k;
            }
        }

        maxRes--;
        int av3 = 0;
        for (int i = 0; i <= ClustalWMaxAminoAcidCode; i++)
            for (int j = 0; j < i; j++)
                av3 += matrix[i, j];

        av3 = ClustalWCastToInt((float)av3 / ((float)(maxRes * maxRes - maxRes) / 2));   // fewer than three residues: 0/0 or x/0 -> INT_MIN on x86, as in the program
        return (matrix, unchecked(-av3));
    }

    /// <summary>
    /// <c>SubMatrix::getProfileAlignMatrix</c>, protein branch: the table of the matrix series that ClustalW uses for two profiles with mean percent identity
    /// <paramref name="pcid"/> and shorter average length <paramref name="minLen"/>, with the scale factor of the gap penalties (a <c>float</c> in ClustalW: <c>PrfScaleValues.scale</c>).
    /// BLOSUM: 80 / 62 (x2) / 45 / 45 (scale 0.5) / 45 (0.6) / 30 (0.6) at &gt; 80 / 60 / 40 / 30 / 20 / else; PAM: 20 / 60 / 120 / 350 at &gt; 80 / 60 / 40 / else;
    /// GONNET (intScale 10, scale 0.375): 80 (scale 0.1875) above 35 %, between 25 and 35 % 250 (profiles shorter than 100) or 120, else 350 or 160;
    /// ID: the identity matrix. <c>-NEGATIVE</c> selects BLOSUM40 / PAM120 / GONNET250 for every pair of profiles.
    /// </summary>
    private static (SubstitutionMatrix Table, int Multiplier, int IntScale, float Scale) ClustalWSelectProfileMatrix(
        ClustalWProfileMatrix matrix, bool negative, double pcid, int minLen, bool distanceTree = true)
    {
        // SubMatrix::getProfileAlignMatrix: `_negMatrix || userParameters->getDistanceTree() == false` selects the single matrix (a profile alignment run while the guide tree has no distances)
        bool single = negative || !distanceTree;
        switch (matrix)
        {
            case ClustalWProfileMatrix.Blosum:
                if (single)
                    return (ClustalWMatrices.ClustalWBlosum40, 1, 100, 0.75f);
                if (pcid > 80.0)
                    return (ClustalWMatrices.ClustalWBlosum80, 1, 100, 0.75f);
                if (pcid > 60.0)
                    return (SubstitutionMatrix.Blosum62, 2, 100, 0.75f);          // blosum62mt2 = 2 x BLOSUM62
                if (pcid > 40.0)
                    return (ClustalWMatrices.ClustalWBlosum45, 1, 100, 0.75f);
                if (pcid > 30.0)
                    return (ClustalWMatrices.ClustalWBlosum45, 1, 100, 0.5f);
                if (pcid > 20.0)
                    return (ClustalWMatrices.ClustalWBlosum45, 1, 100, 0.6f);
                return (ClustalWMatrices.ClustalWBlosum30, 1, 100, 0.6f);
            case ClustalWProfileMatrix.Pam:
                if (single)
                    return (ClustalWMatrices.ClustalWPam120, 1, 100, 0.75f);
                if (pcid > 80.0)
                    return (ClustalWMatrices.ClustalWPam20, 1, 100, 0.75f);
                if (pcid > 60.0)
                    return (ClustalWMatrices.ClustalWPam60, 1, 100, 0.75f);
                if (pcid > 40.0)
                    return (ClustalWMatrices.ClustalWPam120, 1, 100, 0.75f);
                return (ClustalWMatrices.ClustalWPam350, 1, 100, 0.75f);
            case ClustalWProfileMatrix.Gonnet:
                if (single)
                    return (ClustalWMatrices.ClustalWGonnet250, 1, 10, 0.375f);
                if (pcid > 35.0)
                    return (ClustalWMatrices.ClustalWGonnet80, 1, 10, 0.1875f);
                if (pcid > 25.0)
                    return (minLen < 100 ? ClustalWMatrices.ClustalWGonnet250 : ClustalWMatrices.ClustalWGonnet120, 1, 10, 0.375f);
                return (minLen < 100 ? ClustalWMatrices.ClustalWGonnet350 : ClustalWMatrices.ClustalWGonnet160, 1, 10, 0.375f);
            case ClustalWProfileMatrix.Identity:
                return (ClustalWMatrices.ClustalWIdentity, 1, 100, 0.75f);
            default:
                throw new ArgumentOutOfRangeException(nameof(matrix), matrix, "Unknown ClustalW profile matrix.");
        }
    }

    /// <summary>
    /// The matrix of the profile stage for two profiles of mean percent identity <paramref name="pcid"/> and shorter length <paramref name="minLen"/>:
    /// <c>getProfileAlignMatrix</c>' protein branch with the user's matrix file (<c>-MATRIX=file</c>) when given. A single user matrix has scale 0.75
    /// and intScale 100; a user series picks the first entry with <c>llimit &lt;= pcid &lt;= ulimit</c> (the first entry when none matches) with the scale
    /// <c>0.5 + (pcid - llimit) / ((ulimit - llimit) * 2.0)</c> (stored in a float: 0.5 at the lower, 1.0 at the upper limit).
    /// </summary>
    private static (IEnumerable<(int Row, int Column, short Value)> Cells, IReadOnlyList<int> Cross, int IntScale, float Scale) ClustalWChooseProfileMatrix(
        ClustalWProfileMatrix matrix, bool negative, ClustalWUserMatrix? userMatrix, IReadOnlyList<ClustalWMatrixSeriesEntry>? userSeries, double pcid, int minLen, bool distanceTree = true)
    {
        if (userSeries is { Count: > 0 })
        {
            int j = 0;
            for (int i = 0; i < userSeries.Count; i++)
            {
                if (pcid >= userSeries[i].LowerLimit && pcid <= userSeries[i].UpperLimit)
                {
                    j = i;
                    break;
                }
            }

            var entry = userSeries[j];
            float scale = (float)(0.5 + (pcid - entry.LowerLimit) / ((entry.UpperLimit - entry.LowerLimit) * 2.0));
            return (entry.Matrix.Cells(), entry.Matrix.Cross, 100, scale);
        }

        if (userMatrix is not null)
            return (userMatrix.Cells(), userMatrix.Cross, 100, 0.75f);

        var (table, multiplier, intScale, builtinScale) = ClustalWSelectProfileMatrix(matrix, negative, pcid, minLen, distanceTree);
        return (ClustalWTableCells(table, multiplier), ClustalWDefaultAAXref, intScale, builtinScale);
    }

    /// <summary>
    /// <c>SubMatrix::getMatrix</c> for the profile stage: the 32 x 32 matrix over residue codes, shifted to non-negative scores unless <paramref name="negative"/>
    /// (the lowest cell of the lower triangle without column A, as in the C++, is subtracted from every pair of listed residues), gap rows / columns 0; and
    /// <c>matrixAvgScore</c>.
    /// </summary>
    private static (int[,] Matrix, int MatrixAvgScore) BuildClustalWProfileMatrix(
        IEnumerable<(int Row, int Column, short Value)> cells, IReadOnlyList<int> cross, int intScale, bool negative)
    {
        var (matrix, matAvg) = ClustalWGetMatrix(cells, cross, intScale);

        int min = matrix[0, 0];
        for (int i = 0; i <= ClustalWMaxAminoAcidCode; i++)
            for (int j = 1; j <= i; j++)
                min = Math.Min(min, matrix[i, j]);

        if (!negative && min < 0)
        {
            for (int i = 0; i <= ClustalWMaxAminoAcidCode && i < cross.Count; i++)
            {
                int ti = cross[i];
                if (ti < 0)
                    continue;
                for (int j = 0; j <= ClustalWMaxAminoAcidCode && j < cross.Count; j++)
                {
                    int tj = cross[j];
                    if (tj >= 0)
                        matrix[ti, tj] -= min;
                }
            }
        }

        // gap rows / columns (grScore = ggScore = 0): the 32 x 32 array already holds zeros for the codes 25 .. 31
        return (matrix, matAvg);
    }

    private static (int[,] Matrix, int IntScale, int MatrixAvgScore) BuildClustalWUserPairwiseMatrix(ClustalWUserMatrix user)
    {
        var (matrix, avg) = ClustalWGetMatrix(user.Cells(), user.Cross, 100);
        return (matrix, 100, avg);
    }

    /// <summary>ClustalW residue codes of a protein sequence; a symbol ClustalW does not read (J, '*', non-letters) is rejected (the program silently drops it).</summary>
    private static int[] EncodeClustalWProtein(string sequence, string paramName, bool allowGaps = false)
    {
        var codes = new int[sequence.Length];
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            if (allowGaps && c == '-')
            {
                codes[i] = CwGapPos2;   // a gap already in the alignment (gapPos2)
                continue;
            }

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

        // -PWMATRIX=file: getPairwiseMatrix reads the user's matrix (intScale 100, gapOpenScale 1) with the same getMatrix; the pairwise matrix is never shifted
        var (matrix, intScale, matAvg) = options.UserMatrix is { } user
            ? BuildClustalWUserPairwiseMatrix(user)
            : ClustalWProteinMatrices[(int)options.Matrix].Value;

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
