// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.TCoffee;

using Seqeron.Genomics.ClustalW;
using static Seqeron.Genomics.ClustalW.ClustalWAligner;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

/// <summary>
/// The triangular score vectors <c>blosum62mt</c> / <c>blosum50mt</c> of T-Coffee 13.41 (<c>matrices.h</c>), over the 23 letters
/// "abcdefghiklmnpqrstvwxyz" (row a holds a + 1 cells). Transcribed by a script from the Debian source package, 276 values each.
/// </summary>
internal static class TcMatrixTables
{
    public static readonly int[] Blosum62Mt =
    {
        4, -2, 4, 0, -3, 9, -2, 4, -3, 6, -1, 1, -4, 2, 5, -2, -3, -2, -3, -3, 6, 0, -1,
        -3, -1, -2, -3, 6, -2, 0, -3, -1, 0, -1, -2, 8, -1, -3, -1, -3, -3, 0, -4, -3, 4, -1,
        0, -3, -1, 1, -3, -2, -1, -3, 5, -1, -4, -1, -4, -3, 0, -4, -3, 2, -2, 4, -1, -3, -1,
        -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, 3, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -1,
        -2, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, 0, -3, 0, 2, -3, -2, 0, -3, 1,
        -2, 0, 0, -1, 5, -1, -1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, 1, 0,
        -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 0, -1, -1, -1, -1, -2, -2, -2,
        -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -3, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3,
        -2, -2, -3, -2, 0, 4, -3, -4, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3,
        -2, -3, 11, 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1, 0, 0, -1, -2,
        -1, -2, -3, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, -1, 7,
        -1, 1, -3, 1, 4, -3, -2, 0, -3, 1, -3, -1, 0, -1, 3, 0, 0, -1, -2, -3, -1, -2, 4,
    };

    public static readonly int[] Blosum50Mt =
    {
        5, -2, 5, -1, -3, 13, -2, 5, -4, 8, -1, 1, -3, 2, 6, -3, -4, -2, -5, -3, 8, 0, -1,
        -3, -1, -3, -4, 8, -2, 0, -3, -1, 0, -1, -2, 10, -1, -4, -2, -4, -4, 0, -4, -4, 5, -1,
        0, -3, -1, 1, -4, -2, 0, -3, 6, -2, -4, -2, -4, -3, 1, -4, -3, 2, -3, 5, -1, -3, -2,
        -4, -2, 0, -3, -1, 2, -2, 3, 7, -1, 4, -2, 2, 0, -4, 0, 1, -3, 0, -4, -2, 7, -1,
        -2, -4, -1, -1, -4, -2, -2, -3, -1, -4, -3, -2, 10, -1, 0, -3, 0, 2, -4, -2, 1, -3, 2,
        -2, 0, 0, -1, 7, -2, -1, -4, -2, 0, -3, -3, 0, -4, 3, -3, -2, -1, -3, 1, 7, 1, 0,
        -1, 0, -1, -3, 0, -1, -3, 0, -3, -2, 1, -1, 0, -1, 5, 0, 0, -1, -1, -1, -2, -2, -2,
        -1, -1, -1, -1, 0, -1, -1, -1, 2, 5, 0, -4, -1, -4, -3, -1, -4, -4, 4, -3, 1, 1, -3,
        -3, -3, -3, -2, 0, 5, -3, -5, -5, -5, -3, 1, -3, -3, -3, -3, -2, -1, -4, -4, -1, -3, -4,
        -3, -3, 15, -1, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, -2, -1, -1, -1, 0, -1, -3,
        -1, -2, -3, -3, -3, -2, 4, -3, 2, -1, -2, -1, 0, -2, -3, -1, -1, -2, -2, -1, 2, -1, 8,
        -1, 2, -3, 1, 5, -4, -2, 0, -3, 1, -3, -1, 0, -1, 4, 0, 0, -1, -3, -2, -1, -2, 5,
    };
}
