// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854
using System.Text;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The feature output of an ACD item <c>featout</c> (<c>-outfeat</c> with <c>-offormat -ofopenfile -ofname -ofdirectory</c>; <c>acdSetFeatout</c>, <c>ajFeattabOutOpen</c>, <c>featoutUfoProcess</c>, <c>ajFeattableWrite</c>, <c>ajFeattabOutDel</c>: B09 WP63, F282):
/// the file name and format are found, a file that cannot be opened ends the program, a table in a format the sequence type does not have is an error of the library, a file that got no table is a warning.
/// </summary>
public sealed class EmbossFeatout
{
    private readonly EmbossFeatureReadContext _context;
    private readonly string _program;
    private readonly DateTime _date;
    private readonly EmbossFeatureWriteState _state = new();
    private readonly StringBuilder _text = new();
    private int _count;

    private EmbossFeatout(EmbossFeatureReadContext context, string program, int index, string formatStr, string file, DateTime date)
    {
        _context = context;
        _program = program;
        Index = index;
        FormatStr = formatStr;
        File = file;
        _date = date;
    }

    /// <summary>The path of the file the features are written to.</summary>
    public string File { get; }

    /// <summary>The index of the format in the table of the output formats.</summary>
    public int Index { get; }

    /// <summary>The name of the format as the program was given it (the messages print it).</summary>
    public string FormatStr { get; }

    /// <summary>The text written so far.</summary>
    public string Text => _text.ToString();

    /// <summary>The program stops: the messages are in <paramref name="context"/>, the death line is the caller's.</summary>
    /// <returns>The output; null when the program dies.</returns>
    public static EmbossFeatout? Open(EmbossAcdArgs a, EmbossFeatureReadContext context, string master, string program, string baseName, IReadOnlyDictionary<string, string> files, DateTime date)
    {
        string fmt = a.RawAssoc(master, "offormat") ?? "";
        string? openFile = a.RawAssoc(master, "ofopenfile");
        string name = a.RawAssoc(master, "ofname") is { } ofname ? EmbossOutPath.FmtWord(ofname) : baseName;
        string ext = fmt.Length > 0 ? fmt : "gff";                                        // ajFeatOutFormatDefault
        string ufo = a.Raw(master) ?? (name + "." + ext);
        string dir = EmbossOutPath.FixDirectory(a.RawAssoc(master, "ofdirectory"));
        if (ufo.Length == 0)
        {
            context.Error("Output UFO is required");
            return null;
        }

        // featoutUfoProcess
        string formatstr, file;
        if (!string.IsNullOrEmpty(openFile))
        {
            formatstr = fmt;
            file = openFile;
        }
        else
        {
            var m = Regex.Match(ufo, "^([A-Za-z0-9]*):+(.*)$", RegexOptions.CultureInvariant);
            string rest = ufo;
            if (m.Success)
            {
                formatstr = m.Groups[1].Value.Length > 0 ? m.Groups[1].Value : "unknown";
                rest = m.Groups[2].Value;
            }
            else
            {
                formatstr = fmt.Length > 0 ? fmt : "gff";
            }

            file = Regex.IsMatch(rest, "^([^:]+)$", RegexOptions.CultureInvariant) ? rest : "unknown." + formatstr;
        }

        int index = EmbossReportWriter.FeatureFormatIndex(formatstr);
        if (index < 0)
        {
            index = formatstr.Length == 0 ? 0 : 1;               // featFindOutFormatC: a name that is not a format is written as gff3, an empty one stays "unknown"
            context.Error($"unknown output feature format '{formatstr}' will write as {(index == 0 ? "unknown" : "gff3")} instead\n");
        }

        string path = EmbossOutPath.Combine(dir, file);
        if (!EmbossOutPath.DirectoryExists(EmbossReportArgs.ExistingDirectories(files), path))
        {
            context.Warn($"Failed to open output file '{file}' in directory '{dir}' error:2 'No such file or directory'");
            context.Error($"Unable to open features output '{dir}{ufo}'");
            return null;
        }

        return new EmbossFeatout(context, program, index, formatstr, path, date);
    }

    /// <summary><c>ajFeattableWrite</c> of a table.</summary>
    public void Write(EmbossReport table)
    {
        string? text = EmbossReportWriter.WriteFeatureTable(table, Index, _date, _program, _state, _count);
        if (text is null)
        {
            _context.Error($"Feature format '{FormatStr}' not supported for {(table.Protein ? "protein" : "nucleotide")} sequence features");
            return;
        }

        _text.Append(text);
        _count++;
    }

    /// <summary><c>ajFeattabOutDel</c>: a file that got no table is a warning.</summary>
    public void Close()
    {
        if (_count == 0)
            _context.Warn($"No features written to output file '{File}'");
    }
}
