// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

/// <summary>The null checks of the EMBOSS input sequences taken by the programs (moved out of the MIT <c>SequenceAligner</c>, B09 WP80).</summary>
public static class EmbossInputContract
{
    /// <summary>
    /// A set of EMBOSS input sequences (<c>NeedleAll</c> / <c>Supermatcher</c> / ...): a null list, a null member, or a member without <c>Name</c> / <c>Residues</c> is an <see cref="ArgumentNullException"/> naming the parameter and the index
    /// (audit A6-4: the first dereference was a <see cref="NullReferenceException"/> in <c>EmbossSetInput.Any</c>).
    /// </summary>
    public static void RequireSet(IReadOnlyList<EmbossSequence>? set, string paramName)
    {
        SequenceAligner.RequireNoNullElements(set, paramName, "Sequence");
        for (int i = 0; i < set!.Count; i++)
            RequireSequence(set[i], paramName, i);
    }

    /// <summary>One EMBOSS input sequence: a null <c>Name</c> or <c>Residues</c> is an <see cref="ArgumentNullException"/> naming the parameter (and the index when it is a member of a set).</summary>
    public static void RequireSequence(EmbossSequence sequence, string paramName, int index = -1)
    {
        string which = index < 0 ? "The sequence" : $"Sequence {index}";
        if (sequence.Name is null)
            throw new ArgumentNullException(paramName, $"{which} has no Name.");
        if (sequence.Residues is null)
            throw new ArgumentNullException(paramName, $"{which} has no Residues.");
    }
}
