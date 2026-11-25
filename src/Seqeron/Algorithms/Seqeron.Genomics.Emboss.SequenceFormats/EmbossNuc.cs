// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>The nucleotide character functions of the AJAX library (<c>ajbase.c</c>, <c>ajseq.c</c>): the complement of a base character and <c>ajSeqstrReverse</c> (B09 WP61, F280; moved from <c>EmbossSeqRange</c>).</summary>
public static class EmbossNuc
{
    // ajBaseAlphacharComp of every character ("ACGTURYWSMKBDHVNXacgturywsmkbdhvnx" -> "TGCAAYRWSKMVHDBNXtgcaayrwskmvhdbnx", every other character is kept)
    private const string Forward = "ACGTURYWSMKBDHVNXacgturywsmkbdhvnx", Complement = "TGCAAYRWSKMVHDBNXtgcaayrwskmvhdbnx";

    /// <summary><c>ajBaseAlphacharComp</c>.</summary>
    public static char Comp(char c)
    {
        int k = Forward.IndexOf(c, StringComparison.Ordinal);
        return k >= 0 ? Complement[k] : c;
    }

    /// <summary><c>ajSeqstrReverse</c>: the reverse complement.</summary>
    public static string ReverseComplement(string s) => string.Create(s.Length, s, static (dest, src) =>
    {
        for (int i = 0; i < src.Length; i++)
            dest[i] = Comp(src[src.Length - 1 - i]);
    });
}
