// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.PairwiseAlignment;

using Seqeron.Genomics.Analysis;

public static partial class EmbossAlignmentStatistics
{

    /// <summary>
    /// Pairwise distance matrix of aligned rows as EMBOSS <c>distmat</c> (uncorrected and Jukes-Cantor; no ambiguity option, all codon positions - the full option set is <see cref="CalculateDistanceMatrix(IReadOnlyList{string}, MsaDistanceOptions)"/>): matches =
    /// columns where both rows carry the same non-gap symbol; g = columns where either row has a gap;
    /// <see cref="MsaDistanceMethod.Uncorrected"/> D = 1 − matches/(L − g + g·<paramref name="gapWeight"/>);
    /// <see cref="MsaDistanceMethod.JukesCantor"/> D' = −b·ln(1 − D/b). A pair with a zero denominator is NaN
    /// (distmat divides by zero). The matrix is symmetric with a zero diagonal.
    /// </summary>
    /// <param name="alignedRows">At least two aligned rows.</param>
    /// <param name="method">Distance model.</param>
    /// <param name="residueType">Selects b = 3/4 (nucleotide) or 19/20 (protein) for Jukes-Cantor; Auto = EMBOSS type of the first row.</param>
    /// <param name="gapWeight">distmat <c>-gapweight</c> (default 0).</param>
    public static double[,] CalculateDistanceMatrix(
        IReadOnlyList<string> alignedRows, MsaDistanceMethod method = MsaDistanceMethod.Uncorrected,
        ConsensusResidueType residueType = ConsensusResidueType.Auto, double gapWeight = 0)
    {
        return CalculateDistanceMatrix(alignedRows, new MsaDistanceOptions { Method = method, ResidueType = residueType, GapWeight = gapWeight });
    }
}
