// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using System.Text;
using Seqeron.Genomics.Infrastructure;

/// <summary>The output of <c>seqret</c>: the sequence file, the feature file and the messages the library printed.</summary>
/// <param name="SequenceText">The text of the sequence output file.</param>
/// <param name="FeatureText">The text of the feature output file (empty when no feature was written).</param>
/// <param name="Sequences">The number of sequences written.</param>
/// <param name="Messages">What the library printed on its error stream, in order (the ones before the program died, too).</param>
public sealed record EmbossSeqretOutput(string SequenceText, string FeatureText, int Sequences, IReadOnlyList<string> Messages)
{
    /// <summary>The files of <c>-ossingle</c> by path (<c>directory/name.extension</c>); empty when the output is one file.</summary>
    public IReadOnlyDictionary<string, string> SingleFiles { get; init; } = new Dictionary<string, string>();

    /// <summary>The name of the feature output file.</summary>
    public string FeatureFileName { get; init; } = "";

    /// <summary>The path of the sequence output file (<c>-osdirectory</c> + <c>-outseq</c>).</summary>
    public string SequenceFileName { get; init; } = "";

    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}
