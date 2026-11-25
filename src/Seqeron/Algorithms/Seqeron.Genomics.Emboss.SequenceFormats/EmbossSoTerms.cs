// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>The Sequence Ontology of EMBOSS 6.6.0 (<c>data/OBO/so.obo</c>): the name of a term by its id, as <c>ajOboinQryS</c> finds it for the feature writers.</summary>
public static partial class EmbossSoTerms
{
    private static Dictionary<string, (string Name, bool Obsolete)>? _terms;

    /// <summary>The name of the term, or null when the id is no term of the file or the term is obsolete (the writers skip those).</summary>
    public static string? Name(string id)
    {
        _terms ??= Data.GroupBy(d => d.Id).ToDictionary(g => g.Key, g => (g.First().Name, g.First().Obsolete), StringComparer.Ordinal);
        return _terms.TryGetValue(id, out var t) && !t.Obsolete ? t.Name : null;
    }
}
