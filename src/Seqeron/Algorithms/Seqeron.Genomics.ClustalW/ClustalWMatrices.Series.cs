// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Infrastructure;

// ClustalW 2.1's matrix SERIES for the PROTEIN profile stage (audit A1-10c part b). SubMatrix::getProfileAlignMatrix chooses one table of
// matrices.h by the percent identity of the two profiles being aligned; the tables are the lower-triangular arrays blosum40mt / blosum45mt /
// blosum62mt2 / blosum80mt, pam20mt / pam60mt / pam120mt, gon40mt / gon80mt / gon120mt / gon160mt / gon350mt in ClustalW's packing order
// ("ABCDEFGHIKLMNPQRSTVWXYZ", ClustalWAminoAcidOrder), re-derived by scripts/alignment/generate_clustalw_protein_msa_oracle.py from the
// clustalw 2.1+lgpl-7 source (sources.debian.org). The BLOSUM / PAM tables equal the EMBOSS data files on all 529 cells (checked by the oracle
// script and by ClustalWSeriesMatrices_*Tests); blosum62mt2 is exactly 2 x EBLOSUM62 (= 2 x SubstitutionMatrix.Blosum62).
public static partial class ClustalWMatrices
{

    /// <summary>ClustalW 2.1 <c>blosum40mt</c> (BLOSUM40; profile-stage matrix series): EMBOSS data file EBLOSUM40 (all 23 x 23 cells). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWBlosum40 => LazyClustalWBlosum40.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWBlosum40 = new(() => BuildTriangle("CLUSTALW-BLOSUM40", ClustalWBlosum40Triangle()));

    /// <summary>ClustalW 2.1 <c>blosum45mt</c> (BLOSUM45; profile-stage matrix series): EMBOSS data file EBLOSUM45 (all 23 x 23 cells). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWBlosum45 => LazyClustalWBlosum45.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWBlosum45 = new(() => BuildTriangle("CLUSTALW-BLOSUM45", ClustalWBlosum45Triangle()));

    /// <summary>ClustalW 2.1 <c>blosum80mt</c> (BLOSUM80; profile-stage matrix series): EMBOSS data file EBLOSUM80 (all 23 x 23 cells). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWBlosum80 => LazyClustalWBlosum80.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWBlosum80 = new(() => BuildTriangle("CLUSTALW-BLOSUM80", ClustalWBlosum80Triangle()));

    /// <summary>ClustalW 2.1 <c>pam20mt</c> (PAM20; profile-stage matrix series): EMBOSS data file EPAM20 (all 23 x 23 cells). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWPam20 => LazyClustalWPam20.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWPam20 = new(() => BuildTriangle("CLUSTALW-PAM20", ClustalWPam20Triangle()));

    /// <summary>ClustalW 2.1 <c>pam60mt</c> (PAM60; profile-stage matrix series): EMBOSS data file EPAM60 (all 23 x 23 cells). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWPam60 => LazyClustalWPam60.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWPam60 = new(() => BuildTriangle("CLUSTALW-PAM60", ClustalWPam60Triangle()));

    /// <summary>ClustalW 2.1 <c>pam120mt</c> (PAM120; profile-stage matrix series): EMBOSS data file EPAM120 (all 23 x 23 cells). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWPam120 => LazyClustalWPam120.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWPam120 = new(() => BuildTriangle("CLUSTALW-PAM120", ClustalWPam120Triangle()));

    /// <summary>ClustalW 2.1 <c>gon40mt</c> (GONNET40; profile-stage matrix series): ClustalW 2.1 matrices.h only (Gonnet, Cohen &amp; Benner 1992 at 40 PAM, tenth-bits; no local reference carries it). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWGonnet40 => LazyClustalWGonnet40.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet40 = new(() => BuildTriangle("CLUSTALW-GONNET40", ClustalWGonnet40Triangle()));

    /// <summary>ClustalW 2.1 <c>gon80mt</c> (GONNET80; profile-stage matrix series): ClustalW 2.1 matrices.h only (Gonnet et al. 1992 at 80 PAM, tenth-bits). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWGonnet80 => LazyClustalWGonnet80.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet80 = new(() => BuildTriangle("CLUSTALW-GONNET80", ClustalWGonnet80Triangle()));

    /// <summary>ClustalW 2.1 <c>gon120mt</c> (GONNET120; profile-stage matrix series): ClustalW 2.1 matrices.h only (Gonnet et al. 1992 at 120 PAM, tenth-bits). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWGonnet120 => LazyClustalWGonnet120.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet120 = new(() => BuildTriangle("CLUSTALW-GONNET120", ClustalWGonnet120Triangle()));

    /// <summary>ClustalW 2.1 <c>gon160mt</c> (GONNET160; profile-stage matrix series): ClustalW 2.1 matrices.h only (Gonnet et al. 1992 at 160 PAM, tenth-bits). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWGonnet160 => LazyClustalWGonnet160.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet160 = new(() => BuildTriangle("CLUSTALW-GONNET160", ClustalWGonnet160Triangle()));

    /// <summary>ClustalW 2.1 <c>gon350mt</c> (GONNET350; profile-stage matrix series): ClustalW 2.1 matrices.h only (Gonnet et al. 1992 at 350 PAM, tenth-bits). Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (<see cref="ClustalWAminoAcidOrder"/>).</summary>
    public static SubstitutionMatrix ClustalWGonnet350 => LazyClustalWGonnet350.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet350 = new(() => BuildTriangle("CLUSTALW-GONNET350", ClustalWGonnet350Triangle()));

    private static short[] ClustalWBlosum40Triangle() => new short[]
    {
        5,
        -1, 5,
        -2, -2, 16,
        -1, 6, -2, 9,
        -1, 1, -2, 2, 7,
        -3, -3, -2, -4, -3, 9,
        1, -1, -3, -2, -3, -3, 8,
        -2, 0, -4, 0, 0, -2, -2, 13,
        -1, -3, -4, -4, -4, 1, -4, -3, 6,
        -1, 0, -3, 0, 1, -3, -2, -1, -3, 6,
        -2, -3, -2, -3, -2, 2, -4, -2, 2, -2, 6,
        -1, -3, -3, -3, -2, 0, -2, 1, 1, -1, 3, 7,
        -1, 4, -2, 2, -1, -3, 0, 1, -2, 0, -3, -2, 8,
        -2, -2, -5, -2, 0, -4, -1, -2, -2, -1, -4, -2, -2, 11,
        0, 0, -4, -1, 2, -4, -2, 0, -3, 1, -2, -1, 1, -2, 8,
        -2, -1, -3, -1, -1, -2, -3, 0, -3, 3, -2, -1, 0, -3, 2, 9,
        1, 0, -1, 0, 0, -2, 0, -1, -2, 0, -3, -2, 1, -1, 1, -1, 5,
        0, 0, -1, -1, -1, -1, -2, -2, -1, 0, -1, -1, 0, 0, -1, -2, 2, 6,
        0, -3, -2, -3, -3, 0, -4, -4, 4, -2, 2, 1, -3, -3, -3, -2, -1, 1, 5,
        -3, -4, -6, -5, -2, 1, -2, -5, -3, -2, -1, -2, -4, -4, -1, -2, -5, -4, -3, 19,
        0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -2, -1, -1, 0, 0, -1, -2, -1,
        -2, -3, -4, -3, -2, 4, -3, 2, 0, -1, 0, 1, -2, -3, -1, -1, -2, -1, -1, 3, -1, 9,
        -1, 2, -3, 1, 5, -4, -2, 0, -4, 1, -2, -2, 0, -1, 4, 0, 0, -1, -3, -2, -1, -2, 5,
    };

    private static short[] ClustalWBlosum45Triangle() => new short[]
    {
        5,
        -1, 4,
        -1, -2, 12,
        -2, 5, -3, 7,
        -1, 1, -3, 2, 6,
        -2, -3, -2, -4, -3, 8,
        0, -1, -3, -1, -2, -3, 7,
        -2, 0, -3, 0, 0, -2, -2, 10,
        -1, -3, -3, -4, -3, 0, -4, -3, 5,
        -1, 0, -3, 0, 1, -3, -2, -1, -3, 5,
        -1, -3, -2, -3, -2, 1, -3, -2, 2, -3, 5,
        -1, -2, -2, -3, -2, 0, -2, 0, 2, -1, 2, 6,
        -1, 4, -2, 2, 0, -2, 0, 1, -2, 0, -3, -2, 6,
        -1, -2, -4, -1, 0, -3, -2, -2, -2, -1, -3, -2, -2, 9,
        -1, 0, -3, 0, 2, -4, -2, 1, -2, 1, -2, 0, 0, -1, 6,
        -2, -1, -3, -1, 0, -2, -2, 0, -3, 3, -2, -1, 0, -2, 1, 7,
        1, 0, -1, 0, 0, -2, 0, -1, -2, -1, -3, -2, 1, -1, 0, -1, 4,
        0, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 2, 5,
        0, -3, -1, -3, -3, 0, -3, -3, 3, -2, 1, 1, -3, -3, -3, -2, -1, 0, 5,
        -2, -4, -5, -4, -3, 1, -2, -3, -2, -2, -2, -2, -4, -3, -2, -2, -4, -3, -3, 15,
        0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -2, -1,
        -2, -2, -3, -2, -2, 3, -3, 2, 0, -1, 0, 0, -2, -3, -1, -1, -2, -1, -1, 3, -1, 8,
        -1, 2, -3, 1, 4, -3, -2, 0, -3, 1, -2, -1, 0, -1, 4, 0, 0, -1, -3, -2, -1, -2, 4,
    };

    private static short[] ClustalWBlosum80Triangle() => new short[]
    {
        7,
        -3, 6,
        -1, -6, 13,
        -3, 6, -7, 10,
        -2, 1, -7, 2, 8,
        -4, -6, -4, -6, -6, 10,
        0, -2, -6, -3, -4, -6, 9,
        -3, -1, -7, -2, 0, -2, -4, 12,
        -3, -6, -2, -7, -6, -1, -7, -6, 7,
        -1, -1, -6, -2, 1, -5, -3, -1, -5, 8,
        -3, -7, -3, -7, -6, 0, -7, -5, 2, -4, 6,
        -2, -5, -3, -6, -4, 0, -5, -4, 2, -3, 3, 9,
        -3, 5, -5, 2, -1, -6, -1, 1, -6, 0, -6, -4, 9,
        -1, -4, -6, -3, -2, -6, -5, -4, -5, -2, -5, -4, -4, 12,
        -2, -1, -5, -1, 3, -5, -4, 1, -5, 2, -4, -1, 0, -3, 9,
        -3, -2, -6, -3, -1, -5, -4, 0, -5, 3, -4, -3, -1, -3, 1, 9,
        2, 0, -2, -1, -1, -4, -1, -2, -4, -1, -4, -3, 1, -2, -1, -2, 7,
        0, -1, -2, -2, -2, -4, -3, -3, -2, -1, -3, -1, 0, -3, -1, -2, 2, 8,
        -1, -6, -2, -6, -4, -2, -6, -5, 4, -4, 1, 1, -5, -4, -4, -4, -3, 0, 7,
        -5, -8, -5, -8, -6, 0, -6, -4, -5, -6, -4, -3, -7, -7, -4, -5, -6, -5, -5, 16,
        -1, -3, -4, -3, -2, -3, -3, -2, -2, -2, -2, -2, -2, -3, -2, -2, -1, -1, -2, -5, -2,
        -4, -5, -5, -6, -5, 4, -6, 3, -3, -4, -2, -3, -4, -6, -3, -4, -3, -3, -3, 3, -3, 11,
        -2, 0, -7, 1, 6, -6, -4, 0, -6, 1, -5, -3, -1, -2, 5, 0, -1, -2, -4, -5, -1, -4, 6,
    };

    private static short[] ClustalWPam20Triangle() => new short[]
    {
        6,
        -5, 6,
        -8, -14, 10,
        -4, 6, -16, 8,
        -3, 0, -16, 2, 8,
        -9, -12, -15, -17, -16, 9,
        -3, -4, -11, -4, -5, -10, 7,
        -8, -2, -8, -5, -6, -7, -10, 9,
        -6, -7, -7, -9, -6, -3, -13, -11, 9,
        -8, -3, -16, -6, -5, -16, -8, -8, -7, 7,
        -7, -10, -17, -15, -10, -4, -12, -7, -2, -9, 7,
        -6, -12, -16, -13, -8, -5, -10, -13, -2, -3, 0, 11,
        -5, 6, -13, 1, -3, -10, -4, -1, -6, -2, -8, -11, 8,
        -2, -8, -9, -9, -7, -11, -7, -5, -10, -8, -8, -9, -7, 8,
        -5, -4, -16, -4, 0, -15, -8, 0, -9, -4, -6, -5, -5, -4, 9,
        -8, -9, -9, -12, -11, -10, -11, -3, -6, -1, -10, -5, -7, -5, -2, 9,
        -1, -2, -4, -5, -5, -7, -3, -7, -8, -5, -9, -6, -1, -3, -6, -4, 7,
        -1, -4, -9, -6, -7, -10, -7, -8, -3, -4, -8, -5, -3, -5, -7, -8, 0, 7,
        -3, -9, -7, -9, -8, -9, -7, -7, 1, -10, -3, -2, -9, -7, -8, -9, -8, -4, 7,
        -16, -11, -18, -17, -19, -6, -17, -8, -16, -14, -7, -15, -9, -16, -15, -3, -6, -15, -18, 13,
        -4, -6, -11, -7, -6, -9, -6, -6, -6, -6, -7, -6, -4, -6, -6, -7, -4, -5, -6, -13, -6,
        -9, -7, -5, -13, -9, 1, -16, -4, -7, -10, -8, -13, -5, -16, -14, -11, -8, -7, -8, -6, -9, 10,
        -4, -1, -16, 0, 6, -16, -6, -2, -7, -5, -8, -6, -4, -5, 7, -5, -6, -7, -8, -17, -6, -11, 6,
    };

    private static short[] ClustalWPam60Triangle() => new short[]
    {
        5,
        -2, 5,
        -5, -9, 9,
        -2, 5, -10, 7,
        -1, 2, -10, 3, 7,
        -6, -8, -9, -11, -10, 8,
        0, -2, -7, -2, -2, -7, 6,
        -5, 0, -6, -2, -3, -4, -6, 8,
        -3, -4, -4, -5, -4, -1, -7, -6, 7,
        -5, -1, -10, -2, -3, -10, -5, -4, -4, 6,
        -4, -7, -11, -9, -7, -1, -8, -4, 0, -6, 6,
        -3, -6, -10, -7, -5, -2, -6, -7, 1, 0, 2, 10,
        -2, 5, -7, 2, 0, -6, -1, 1, -4, 0, -5, -6, 6,
        0, -4, -6, -5, -3, -7, -4, -2, -6, -4, -5, -6, -4, 7,
        -3, -1, -10, -1, 2, -9, -5, 2, -5, -1, -3, -2, -2, -1, 7,
        -5, -5, -6, -6, -6, -7, -7, 0, -4, 2, -6, -2, -3, -2, 0, 8,
        1, 0, -1, -2, -2, -5, 0, -4, -4, -2, -6, -4, 1, 0, -3, -2, 5,
        1, -2, -5, -3, -4, -6, -3, -5, -1, -2, -5, -2, -1, -2, -4, -4, 1, 6,
        -1, -5, -4, -6, -4, -5, -4, -5, 3, -6, -1, 0, -5, -4, -5, -5, -4, -1, 6,
        -10, -8, -12, -11, -12, -3, -11, -5, -10, -8, -4, -9, -6, -10, -9, 0, -4, -9, -11, 13,
        -2, -3, -6, -3, -3, -5, -3, -3, -3, -3, -4, -3, -2, -3, -3, -4, -2, -2, -3, -8, -3,
        -6, -5, -2, -8, -7, 3, -10, -2, -4, -7, -5, -7, -3, -10, -8, -8, -5, -5, -5, -3, -5, 9,
        -2, 1, -10, 2, 5, -10, -3, 0, -4, -2, -5, -4, -1, -2, 6, -2, -3, -4, -5, -11, -3, -7, 5,
    };

    private static short[] ClustalWPam120Triangle() => new short[]
    {
        3,
        0, 4,
        -3, -6, 9,
        0, 4, -7, 5,
        0, 3, -7, 3, 5,
        -4, -5, -6, -7, -7, 8,
        1, 0, -4, 0, -1, -5, 5,
        -3, 1, -4, 0, -1, -3, -4, 7,
        -1, -3, -3, -3, -3, 0, -4, -4, 6,
        -2, 0, -7, -1, -1, -7, -3, -2, -3, 5,
        -3, -4, -7, -5, -4, 0, -5, -3, 1, -4, 5,
        -2, -4, -6, -4, -3, -1, -4, -4, 1, 0, 3, 8,
        -1, 3, -5, 2, 1, -4, 0, 2, -2, 1, -4, -3, 4,
        1, -2, -4, -3, -2, -5, -2, -1, -3, -2, -3, -3, -2, 6,
        -1, 0, -7, 1, 2, -6, -3, 3, -3, 0, -2, -1, 0, 0, 6,
        -3, -2, -4, -3, -3, -5, -4, 1, -2, 2, -4, -1, -1, -1, 1, 6,
        1, 0, 0, 0, -1, -3, 1, -2, -2, -1, -4, -2, 1, 1, -2, -1, 3,
        1, 0, -3, -1, -2, -4, -1, -3, 0, -1, -3, -1, 0, -1, -2, -2, 2, 4,
        0, -3, -3, -3, -3, -3, -2, -3, 3, -4, 1, 1, -3, -2, -3, -3, -2, 0, 5,
        -7, -6, -8, -8, -8, -1, -8, -3, -6, -5, -3, -6, -4, -7, -6, 1, -2, -6, -8, 12,
        -1, -1, -4, -2, -1, -3, -2, -2, -1, -2, -2, -2, -1, -2, -1, -2, -1, -1, -1, -5, -2,
        -4, -3, -1, -5, -5, 4, -6, -1, -2, -5, -2, -4, -2, -6, -5, -5, -3, -3, -3, -2, -3, 8,
        -1, 2, -7, 3, 4, -6, -2, 1, -3, -1, -3, -2, 0, -1, 4, -1, -1, -2, -3, -7, -1, -5, 4,
    };

    private static short[] ClustalWGonnet40Triangle() => new short[]
    {
        92,
        0, 0,
        -31, 0, 163,
        -56, 0, -135, 111,
        -37, 0, -140, 16, 105,
        -92, 0, -64, -152, -143, 126,
        -32, 0, -91, -51, -76, -152, 105,
        -65, 0, -67, -41, -40, -50, -81, 145,
        -76, 0, -87, -150, -106, -39, -158, -94, 104,
        -54, 0, -132, -47, -13, -127, -79, -34, -86, 103,
        -68, 0, -85, -155, -108, -13, -141, -85, 5, -85, 89,
        -45, 0, -63, -130, -80, -16, -114, -60, 10, -57, 16, 140,
        -62, 0, -83, 6, -38, -104, -40, -7, -99, -20, -112, -91, 115,
        -37, 0, -137, -69, -60, -128, -87, -71, -108, -62, -83, -119, -78, 124,
        -43, 0, -113, -32, 10, -100, -71, 0, -91, 2, -60, -35, -25, -46, 118,
        -61, 0, -86, -77, -50, -130, -69, -31, -103, 19, -84, -81, -47, -73, -6, 112,
        0, 0, -35, -36, -41, -111, -37, -48, -95, -43, -95, -64, -11, -35, -35, -51, 99,
        -25, 0, -59, -47, -52, -90, -85, -46, -51, -34, -78, -44, -27, -42, -39, -52, 13, 100,
        -22, 0, -43, -133, -74, -58, -122, -98, 28, -82, -18, -22, -103, -86, -79, -88, -74, -25, 97,
        -120, 0, -68, -171, -131, -6, -108, -70, -93, -127, -71, -72, -119, -149, -87, -63, -98, -120, -115, 181,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -95, 0, -56, -98, -107, 31, -129, 5, -76, -88, -64, -66, -62, -106, -81, -75, -69, -87, -73, 1, 0, 135,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    private static short[] ClustalWGonnet80Triangle() => new short[]
    {
        75,
        0, 0,
        -10, 0, 154,
        -31, 0, -93, 96,
        -17, 0, -94, 31, 88,
        -64, 0, -39, -111, -102, 114,
        -11, 0, -61, -26, -47, -115, 97,
        -39, 0, -43, -17, -17, -26, -53, 127,
        -43, 0, -54, -106, -73, -15, -114, -64, 86,
        -30, 0, -88, -21, 4, -89, -50, -12, -59, 85,
        -43, 0, -55, -109, -75, 7, -104, -57, 22, -58, 77,
        -26, 0, -39, -88, -53, 3, -83, -38, 25, -37, 31, 117,
        -34, 0, -55, 21, -13, -75, -18, 9, -71, -2, -79, -62, 97,
        -16, 0, -93, -42, -35, -93, -58, -45, -75, -37, -58, -78, -48, 114,
        -22, 0, -76, -9, 23, -70, -44, 14, -60, 17, -39, -19, -6, -24, 95,
        -36, 0, -60, -44, -23, -90, -43, -10, -71, 33, -58, -53, -22, -45, 11, 97,
        14, 0, -15, -14, -19, -77, -16, -25, -62, -20, -64, -41, 5, -14, -15, -27, 78,
        -5, 0, -34, -24, -27, -62, -52, -24, -28, -15, -49, -25, -7, -20, -18, -27, 25, 81,
        -6, 0, -21, -89, -51, -31, -86, -65, 41, -54, 3, 1, -69, -57, -51, -60, -43, -9, 80,
        -87, 0, -43, -124, -98, 16, -81, -43, -63, -89, -44, -45, -86, -112, -62, -41, -72, -87, -80, 173,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -65, 0, -32, -69, -74, 49, -94, 21, -47, -60, -35, -37, -39, -76, -53, -50, -46, -58, -47, 23, 0, 123,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    private static short[] ClustalWGonnet120Triangle() => new short[]
    {
        59,
        0, 0,
        -1, 0, 144,
        -18, 0, -69, 82,
        -9, 0, -68, 35, 72,
        -48, 0, -26, -87, -78, 102,
        -3, 0, -45, -14, -31, -92, 90,
        -26, 0, -31, -7, -6, -14, -37, 110,
        -27, 0, -36, -80, -55, -3, -87, -48, 72,
        -19, 0, -64, -8, 11, -67, -34, -2, -44, 69,
        -30, 0, -39, -82, -57, 15, -82, -42, 28, -44, 66,
        -17, 0, -26, -64, -40, 11, -65, -28, 29, -27, 34, 95,
        -20, 0, -41, 26, -1, -58, -7, 14, -55, 5, -61, -46, 80,
        -6, 0, -68, -28, -22, -72, -41, -31, -56, -24, -44, -56, -32, 105,
        -12, 0, -56, 1, 25, -53, -30, 17, -43, 20, -30, -14, 1, -14, 74,
        -23, 0, -45, -27, -10, -68, -30, -1, -53, 36, -44, -38, -10, -30, 16, 83,
        16, 0, -7, -5, -9, -58, -6, -14, -44, -10, -47, -29, 10, -5, -7, -15, 60,
        2, 0, -21, -13, -15, -47, -35, -14, -17, -6, -34, -16, 0, -10, -9, -16, 26, 64,
        0, 0, -11, -65, -38, -17, -65, -47, 42, -39, 13, 10, -50, -42, -36, -44, -28, -3, 65,
        -68, 0, -29, -96, -78, 27, -66, -28, -46, -68, -29, -31, -68, -89, -49, -30, -57, -67, -59, 166,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -48, 0, -20, -53, -56, 55, -74, 26, -31, -44, -20, -22, -28, -59, -38, -37, -35, -42, -33, 33, 0, 111,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    private static short[] ClustalWGonnet160Triangle() => new short[]
    {
        46,
        0, 0,
        3, 0, 135,
        -11, 0, -53, 70,
        -4, 0, -52, 34, 59,
        -38, 0, -18, -70, -62, 91,
        2, 0, -34, -7, -21, -76, 82,
        -18, 0, -23, -1, -1, -7, -27, 93,
        -18, 0, -25, -62, -43, 3, -70, -37, 59,
        -12, 0, -48, -1, 13, -53, -24, 2, -35, 55,
        -22, 0, -29, -65, -45, 19, -67, -32, 30, -34, 57,
        -12, 0, -19, -50, -31, 14, -52, -21, 29, -21, 34, 76,
        -12, 0, -31, 26, 5, -47, -2, 15, -44, 8, -48, -36, 65,
        -1, 0, -52, -19, -14, -58, -30, -22, -43, -16, -35, -42, -22, 96,
        -7, 0, -42, 6, 23, -41, -21, 17, -32, 20, -24, -12, 5, -8, 56,
        -16, 0, -35, -16, -3, -53, -21, 3, -41, 35, -35, -29, -4, -21, 17, 71,
        16, 0, -2, 0, -3, -45, -1, -8, -33, -4, -36, -23, 11, 0, -2, -9, 44,
        5, 0, -14, -6, -8, -36, -24, -8, -12, -2, -24, -11, 3, -4, -4, -9, 23, 50,
        1, 0, -6, -49, -30, -8, -52, -35, 40, -30, 17, 14, -38, -32, -27, -34, -20, 0, 53,
        -55, 0, -21, -78, -64, 32, -55, -19, -34, -54, -20, -22, -55, -74, -40, -24, -47, -54, -45, 158,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -37, 0, -13, -42, -44, 56, -60, 27, -20, -35, -11, -13, -22, -48, -29, -29, -28, -32, -24, 38, 0, 100,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    private static short[] ClustalWGonnet350Triangle() => new short[]
    {
        10,
        0, 0,
        4, 0, 93,
        0, 0, -19, 29,
        1, 0, -17, 19, 20,
        -14, 0, -3, -30, -25, 51,
        5, 0, -12, 4, -2, -35, 51,
        -4, 0, -8, 5, 4, 1, -7, 33,
        -4, 0, -5, -24, -17, 11, -29, -13, 27,
        -1, 0, -16, 6, 9, -21, -4, 5, -13, 18,
        -7, 0, -7, -25, -18, 18, -30, -11, 22, -14, 28,
        -4, 0, -4, -19, -13, 14, -23, -8, 19, -9, 21, 23,
        0, 0, -11, 15, 9, -20, 5, 8, -18, 7, -19, -14, 20,
        3, 0, -18, -2, 0, -25, -7, -5, -16, -2, -15, -14, -3, 56,
        0, 0, -14, 8, 11, -16, -4, 7, -11, 10, -11, -7, 6, 0, 12,
        -2, 0, -13, 2, 6, -20, -4, 6, -15, 18, -14, -11, 4, -4, 10, 28,
        6, 0, 1, 5, 3, -18, 5, 0, -11, 2, -13, -9, 6, 4, 2, 1, 10,
        4, 0, -2, 2, 1, -13, -5, -1, -3, 2, -7, -4, 4, 2, 1, 0, 8, 11,
        0, 0, 2, -18, -12, 5, -21, -11, 22, -10, 16, 14, -13, -11, -9, -12, -6, 0, 21,
        -24, 0, -4, -35, -29, 35, -30, -3, -9, -23, -1, -3, -24, -34, -19, -12, -22, -23, -14, 124,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -14, 0, -1, -18, -17, 42, -27, 15, -1, -14, 5, 2, -10, -20, -11, -12, -12, -12, -4, 39, 0, 57,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
}
