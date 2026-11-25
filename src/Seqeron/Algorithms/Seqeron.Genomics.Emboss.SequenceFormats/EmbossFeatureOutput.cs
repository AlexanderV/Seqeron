// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The files of the EMBOSS 6.6.0 program <c>diffseq</c> (<c>diffseq.c</c>; B09 A5-8c3, F196): <see cref="EmbossPairwisePrograms.Diffseq(string, string, int, bool, CancellationToken)"/> for the differences, the seven-tag report
/// (<c>reportWriteDiffseq</c> and every other <c>-rformat</c>) written by the one <see cref="EmbossReportWriter"/>, and the GFF3 feature files of both sequences (<c>diffseq_WordMatchListConvDiffToFeat</c>).
/// </summary>
public static class EmbossFeatureOutput
{

    // acdSetFeatout + ajFeattabOutOpen + featoutUfoProcess: the format and the file name of one feature output from its value on the command line (the default is name.diffgff) and -offormat / -ofopenfile
    public static (int Format, string File) OpenFeatureOutput(string ufo, EmbossFeatureOutQualifiers? q, string defaultFormat = "")
    {
        string fmt = q?.Format ?? "";
        string formatstr, file;
        if (!string.IsNullOrEmpty(q?.OpenFile))
        {
            formatstr = fmt.Length > 0 ? fmt : defaultFormat;      // the ACD attribute offormat of the output (wordmatch: gff; diffseq has none)
            file = q.OpenFile;
        }
        else
        {
            var m = System.Text.RegularExpressions.Regex.Match(ufo, "^([A-Za-z0-9]*):+(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant);
            string rest = ufo;
            if (m.Success)
            {
                formatstr = m.Groups[1].Value.Length > 0 ? m.Groups[1].Value : "unknown";
                rest = m.Groups[2].Value;
            }
            else
                formatstr = fmt.Length > 0 ? fmt : "gff";
            file = System.Text.RegularExpressions.Regex.IsMatch(rest, "^([^:]+)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant) ? rest : "unknown." + formatstr;
        }

        int index = EmbossReportWriter.FeatureFormatIndex(formatstr);
        if (index < 0)
            index = formatstr.Length == 0 ? 0 : 1;               // featFindOutFormatC: a name that is not a format is written as gff3, an empty one stays "unknown"
        return (index, file);
    }
}
