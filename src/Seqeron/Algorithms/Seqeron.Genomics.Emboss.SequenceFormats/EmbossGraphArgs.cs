// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary><c>acdSetGraphxy</c> for the ACD item <c>xygraph</c> of a program (<c>-graph</c> and its associated qualifiers <c>-gdesc -gtitle -gsubtitle -gxtitle -gytitle -goutfile -gdirectory</c>; B09 WP63, F282).</summary>
public static class EmbossGraphArgs
{
    /// <summary>
    /// The graph of the item <paramref name="item"/>: the value is the device (default <c>x11</c>: the variable <c>EMBOSS_GRAPHICS</c> is not defined), then the title, description, subtitle, labels, the output file and its directory as <c>acdSetGraphxy</c> sets them.
    /// </summary>
    /// <returns>The graph; null when the program dies on the value (the messages are in <paramref name="context"/>).</returns>
    public static EmbossGraphXy? Open(EmbossAcdArgs a, EmbossFeatureReadContext context, string program, string item, string defaultDescription, int multi, string firstName, IReadOnlyDictionary<string, string> files)
    {
        string value = a.Raw(item) ?? "x11";
        if (value.Length == 0)
        {
            context.Error("Graph is required");
            return null;
        }

        var graph = EmbossGraphXy.Create(program, value, multi, context);
        if (graph is null)
            return null;
        graph.Dataname = firstName;
        graph.Desc = a.RawAssoc(item, "gdesc") ?? defaultDescription;
        string? title = a.RawAssoc(item, "gtitle");
        if (title is not null)
            graph.Title = title;
        graph.InitialTitle = graph.Title;
        string? subtitle = a.RawAssoc(item, "gsubtitle");
        if (subtitle is not null)
            graph.Subtitle = subtitle;
        graph.Xaxis = a.RawAssoc(item, "gxtitle") ?? firstName;
        string? ytitle = a.RawAssoc(item, "gytitle");
        if (ytitle is not null)
            graph.Yaxis = ytitle;
        string? outfile = a.RawAssoc(item, "goutfile");
        if (outfile is not null)
            graph.SetOutfile(outfile);
        string? dir = a.RawAssoc(item, "gdirectory");
        if (dir is not null)
            graph.SetOutdir(dir);
        graph.ExistingDirectories = EmbossReportArgs.ExistingDirectories(files);
        return graph;
    }
}
