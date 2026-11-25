// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>A run of <c>entret</c> (A5-8c3h, B09 F257).</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA): a file name, a query on a file or a database, a list file, with the prefixes and ranges of the USA forms.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers.</param>
/// <param name="FirstOnly"><c>-firstonly</c>: read one entry and stop.</param>
public sealed record EmbossEntretRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, bool FirstOnly = false)
{
    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }
}

/// <summary>The output of <c>entret</c>.</summary>
/// <param name="Text">The text of the output file: the entries as the readers stored them.</param>
/// <param name="Entries">The number of entries written.</param>
/// <param name="Messages">What the library printed on its error stream.</param>
public sealed record EmbossEntretOutput(string Text, int Entries, IReadOnlyList<string> Messages)
{
    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}

/// <summary>
/// <c>entret</c> of EMBOSS 6.6.0 (<c>emboss/entret.c</c>): the sequences of the input are read with the text of the entries stored (the ACD <c>seqall</c> has <c>entry: Y</c>) and the text of every entry is written to the output file
/// (<c>ajSeqGetEntryS</c>: the lines the sequence reader read, <c>ajTextinStoreReadline</c>).
/// </summary>
public static class EmbossEntretProgram
{
    /// <summary>Runs <c>entret</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program does not stop or dies by a fault.</exception>
    public static EmbossEntretOutput Entret(EmbossEntretRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "entret", StoreText = true };
        var seqall = EmbossSeqall.Open("entret", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, false, context, "any", false, r.Definitions, true);
        if (seqall is null)
        {
            var died = EmbossSeqretPrograms.Died(context);
            return new EmbossEntretOutput("", 0, died.Messages) { Died = true };
        }

        var sb = new StringBuilder();
        int entries = 0;
        while (EmbossSeqretPrograms.TryNext(seqall, context, out var seq))
        {
            if (seq.TextPtr is { } text)
                sb.Append(text);
            entries++;
            if (r.FirstOnly)
                break;
        }

        return new EmbossEntretOutput(sb.ToString(), entries, [.. context.Messages]) { Died = context.Fatal };
    }
}
