// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Infrastructure;

public static partial class ClustalWMatrices
{
    /// <summary>
    /// ClustalW 2.1's GONNET250 matrix (<c>gon250mt</c> of <c>matrices.h</c>; the default PROTEIN pairwise matrix <c>-PWMATRIX=GONNET</c>): Gonnet, Cohen &amp;
    /// Benner 1992 (Science 256:1443) log-odds at 250 PAM in tenth-bits, i.e. every cell equals <c>round(10 × Biopython GONNET1992)</c> on the 20 amino acids.
    /// Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (ClustalW's packing order; B = Asx, Z = Glx, X = unknown).
    /// </summary>
    public static SubstitutionMatrix ClustalWGonnet250 => LazyClustalWGonnet250.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet250 = new(() => BuildTriangle("CLUSTALW-GONNET250", ClustalWGon250Triangle()));

    /// <summary>
    /// ClustalW 2.1's BLOSUM30 matrix (<c>blosum30mt</c>; the PROTEIN pairwise matrix <c>-PWMATRIX=BLOSUM</c>): Henikoff &amp; Henikoff 1992, equal on all 23 × 23 cells to
    /// the EMBOSS data file <c>EBLOSUM30</c>.
    /// </summary>
    public static SubstitutionMatrix ClustalWBlosum30 => LazyClustalWBlosum30.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWBlosum30 = new(() => BuildTriangle("CLUSTALW-BLOSUM30", ClustalWBlosum30Triangle()));

    /// <summary>
    /// ClustalW 2.1's PAM350 matrix (<c>pam350mt</c>; the PROTEIN pairwise matrix <c>-PWMATRIX=PAM</c>): the Dayhoff mutation-data log-odds at 350 PAM as distributed
    /// with ClustalW; equal on all 23 × 23 cells to the EMBOSS data file <c>EPAM350</c>).
    /// </summary>
    public static SubstitutionMatrix ClustalWPam350 => LazyClustalWPam350.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWPam350 = new(() => BuildTriangle("CLUSTALW-PAM350", ClustalWPam350Triangle()));

    /// <summary>ClustalW 2.1's identity matrix (<c>idmat</c>; <c>-PWMATRIX=ID</c>): 10 for every identical pair of residues (incl. B, X, Z), 0 otherwise.</summary>
    public static SubstitutionMatrix ClustalWIdentity => LazyClustalWIdentity.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWIdentity = new(() => BuildTriangle("CLUSTALW-ID", ClustalWIdTriangle()));

    /// <summary>ClustalW's amino-acid packing order of the lower-triangular tables of <c>matrices.h</c> (23 residues).</summary>
    public const string ClustalWAminoAcidOrder = "ABCDEFGHIKLMNPQRSTVWXYZ";

    private static SubstitutionMatrix BuildTriangle(string name, short[] triangle)
    {
        const string alphabet = ClustalWAminoAcidOrder;
        int n = alphabet.Length;
        var m = new int[n, n];
        int ix = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                m[i, j] = triangle[ix];
                m[j, i] = triangle[ix];
                ix++;
            }

        return new SubstitutionMatrix(name, alphabet, m);
    }


    private static short[] ClustalWGon250Triangle() => new short[]
    {
        24,
        0, 0,
        5, 0, 115,
        -3, 0, -32, 47,
        0, 0, -30, 27, 36,
        -23, 0, -8, -45, -39, 70,
        5, 0, -20, 1, -8, -52, 66,
        -8, 0, -13, 4, 4, -1, -14, 60,
        -8, 0, -11, -38, -27, 10, -45, -22, 40,
        -4, 0, -28, 5, 12, -33, -11, 6, -21, 32,
        -12, 0, -15, -40, -28, 20, -44, -19, 28, -21, 40,
        -7, 0, -9, -30, -20, 16, -35, -13, 25, -14, 28, 43,
        -3, 0, -18, 22, 9, -31, 4, 12, -28, 8, -30, -22, 38,
        3, 0, -31, -7, -5, -38, -16, -11, -26, -6, -23, -24, -9, 76,
        -2, 0, -24, 9, 17, -26, -10, 12, -19, 15, -16, -10, 7, -2, 27,
        -6, 0, -22, -3, 4, -32, -10, 6, -24, 27, -22, -17, 3, -9, 15, 47,
        11, 0, 1, 5, 2, -28, 4, -2, -18, 1, -21, -14, 9, 4, 2, -2, 22,
        6, 0, -5, 0, -1, -22, -11, -3, -6, 1, -13, -6, 5, 1, 0, -2, 15, 25,
        1, 0, 0, -29, -19, 1, -33, -20, 31, -17, 18, 16, -22, -18, -15, -20, -10, 0, 34,
        -36, 0, -10, -52, -43, 36, -40, -8, -18, -35, -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -22, 0, -5, -28, -27, 51, -40, 22, -7, -21, 0, -2, -14, -31, -17, -18, -19, -19, -11, 41, 0, 78,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    private static short[] ClustalWBlosum30Triangle() => new short[]
    {
        4,
        0, 5,
        -3, -2, 17,
        0, 5, -3, 9,
        0, 0, 1, 1, 6,
        -2, -3, -3, -5, -4, 10,
        0, 0, -4, -1, -2, -3, 8,
        -2, -2, -5, -2, 0, -3, -3, 14,
        0, -2, -2, -4, -3, 0, -1, -2, 6,
        0, 0, -3, 0, 2, -1, -1, -2, -2, 4,
        -1, -1, 0, -1, -1, 2, -2, -1, 2, -2, 4,
        1, -2, -2, -3, -1, -2, -2, 2, 1, 2, 2, 6,
        0, 4, -1, 1, -1, -1, 0, -1, 0, 0, -2, 0, 8,
        -1, -2, -3, -1, 1, -4, -1, 1, -3, 1, -3, -4, -3, 11,
        1, -1, -2, -1, 2, -3, -2, 0, -2, 0, -2, -1, -1, 0, 8,
        -1, -2, -2, -1, -1, -1, -2, -1, -3, 1, -2, 0, -2, -1, 3, 8,
        1, 0, -2, 0, 0, -1, 0, -1, -1, 0, -2, -2, 0, -1, -1, -1, 4,
        1, 0, -2, -1, -2, -2, -2, -2, 0, -1, 0, 0, 1, 0, 0, -3, 2, 5,
        1, -2, -2, -2, -3, 1, -3, -3, 4, -2, 1, 0, -2, -4, -3, -1, -1, 1, 5,
        -5, -5, -2, -4, -1, 1, 1, -5, -3, -2, -2, -3, -7, -3, -1, 0, -3, -5, -3, 20,
        0, -1, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, -2, -1,
        -4, -3, -6, -1, -2, 3, -3, 0, -1, -1, 3, -1, -4, -2, -1, 0, -2, -1, 1, 5, -1, 9,
        0, 0, 0, 0, 5, -4, -2, 0, -3, 1, -1, -1, -1, 0, 4, 0, -1, -1, -3, -1, 0, -2, 4,
    };

    private static short[] ClustalWPam350Triangle() => new short[]
    {
        2,
        1, 3,
        -2, -5, 18,
        1, 3, -6, 4,
        1, 3, -6, 4, 4,
        -4, -5, -5, -6, -6, 13,
        2, 1, -4, 1, 1, -6, 5,
        -1, 1, -4, 1, 1, -2, -2, 7,
        0, -2, -3, -2, -2, 2, -2, -2, 5,
        -1, 1, -6, 1, 0, -6, -1, 1, -2, 5,
        -2, -4, -7, -4, -4, 3, -4, -2, 4, -3, 8,
        -1, -2, -6, -3, -2, 1, -3, -2, 3, 0, 5, 6,
        0, 2, -4, 2, 2, -4, 1, 2, -2, 1, -3, -2, 2,
        1, 0, -3, 0, 0, -5, 0, 0, -2, -1, -3, -2, 0, 6,
        0, 2, -6, 2, 3, -5, -1, 3, -2, 1, -2, -1, 1, 1, 4,
        -1, 0, -4, -1, 0, -5, -2, 2, -2, 4, -3, 0, 1, 0, 2, 7,
        1, 1, 0, 1, 0, -4, 1, -1, -1, 0, -3, -2, 1, 1, 0, 0, 1,
        1, 0, -2, 0, 0, -3, 1, -1, 0, 0, -2, -1, 1, 1, 0, -1, 1, 2,
        0, -2, -2, -2, -2, -1, -1, -2, 4, -2, 3, 2, -2, -1, -2, -3, -1, 0, 5,
        -7, -6, -10, -8, -8, 1, -8, -3, -6, -4, -2, -5, -5, -7, -5, 4, -3, -6, -7, 27,
        0, 0, -3, -1, 0, -2, -1, 0, 0, -1, -1, 0, 0, 0, 0, -1, 0, 0, 0, -5, -1,
        -4, -4, 1, -5, -5, 11, -6, 0, 0, -5, 0, -2, -3, -6, -5, -5, -3, -3, -2, 1, -2, 14,
        0, 2, -6, 3, 3, -6, 0, 2, -2, 1, -3, -2, 2, 0, 3, 1, 0, 0, -2, -7, 0, -5, 3,
    };

    private static short[] ClustalWIdTriangle() => new short[]
    {
        10,
        0, 10,
        0, 0, 10,
        0, 0, 0, 10,
        0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
    };
}
