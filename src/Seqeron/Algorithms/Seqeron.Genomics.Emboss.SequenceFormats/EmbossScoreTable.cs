// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

/// <summary>
/// The EMBOSS <c>AjPMatrixf</c> + <c>AjPSeqCvt</c> pair: the single-precision score of two residues through the matrix's own symbol codes. A residue that is not in the matrix (the program's conversion table
/// maps it to a symbol that scores 0) has the extra code <see cref="UnknownCode"/>; both cases of a letter have the same code. One table serves needle, water, supermatcher, needleall and the similarity count.
/// </summary>
public sealed class EmbossScoreTable
{
    private readonly string _alphabet;
    private readonly int _nsym;
    private readonly float[] _sub;

    public EmbossScoreTable(SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        _alphabet = matrix.Alphabet;
        _nsym = _alphabet.Length + 1;
        _sub = new float[_nsym * _nsym];
        for (int i = 0; i < _alphabet.Length; i++)
            for (int j = 0; j < _alphabet.Length; j++)
                _sub[(i * _nsym) + j] = matrix.TryGetScore(_alphabet[i], _alphabet[j], out int sc) ? sc : 0;
    }

    /// <summary>The code of a residue that the matrix does not contain.</summary>
    public int UnknownCode => _alphabet.Length;

    /// <summary>The matrix symbol code of a residue (<c>ajSeqcvtGetCodeK</c>).</summary>
    public int Code(char c)
    {
        int k = _alphabet.IndexOf(char.ToUpperInvariant(c), StringComparison.Ordinal);
        return k < 0 ? _alphabet.Length : k;
    }

    /// <summary>The codes of every residue of a sequence.</summary>
    public int[] Codes(string s)
    {
        var r = new int[s.Length];
        for (int i = 0; i < r.Length; i++)
            r[i] = Code(s[i]);
        return r;
    }

    /// <summary>The score of two symbol codes.</summary>
    public float Score(int x, int y) => _sub[(x * _nsym) + y];

    /// <summary>The score of two residues.</summary>
    public float Score(char a, char b) => _sub[(Code(a) * _nsym) + Code(b)];
}
