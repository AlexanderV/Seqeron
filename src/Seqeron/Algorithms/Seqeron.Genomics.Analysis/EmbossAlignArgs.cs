using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The adapter between the ACD command line of the stream runner and the alignment writer (<see cref="EmbossAlignmentWriter"/>): <c>acdSetAlign</c> of an ACD item <c>align</c> (<c>-outfile</c>, its associated qualifiers
/// <c>-aformat -aextension -adirectory -aname -awidth -aaccshow -adesshow -ausashow -aglobal</c> and the default format of the program), the opening of the file and the writer (B09 WP63, F282).
/// </summary>
internal static class EmbossAlignArgs
{
    /// <summary>The alignment output of an open ACD item.</summary>
    internal sealed record Opened(EmbossAlignmentWriter Writer, string File, EmbossAlignFormat Format);

    /// <summary>
    /// <c>acdSetAlign</c>: the associated qualifiers are read (a bad value dies), the format must exist ("Unknown alignment format"), the file is opened (<c>-outfile</c> as given, else <c>BASE.EXT</c> of <c>-aname</c> / the name of the first sequence and
    /// <c>-aextension</c> / the program name; <c>-adirectory</c> in front of a name that has none).
    /// </summary>
    /// <returns>The output; null when the program dies (the messages are in <paramref name="context"/>).</returns>
    public static Opened? Open(EmbossAcdArgs a, EmbossFeatureReadContext context, string master, string program, string defaultFormat, string baseName, DateTime runDate, IReadOnlyDictionary<string, string> files,
        SubstitutionMatrix? matrix, double? gapOpen, double? gapExtend)
    {
        string formatText = a.RawAssoc(master, "aformat") ?? defaultFormat;
        int? width = a.IntAssoc(master, "awidth");
        bool? accshow = a.BoolAssoc(master, "aaccshow"), desshow = a.BoolAssoc(master, "adesshow"), usashow = a.BoolAssoc(master, "ausashow"), global = a.BoolAssoc(master, "aglobal");
        if (!EmbossAlignmentWriter.TryParseFormat(formatText, out var format))
        {
            context.Error($"Unknown alignment format '{formatText}'");
            context.Messages.Add($"Died: Alignment option -{master}: Validation failed");
            context.Fatal = true;
            return null;
        }

        string? given = a.Raw(master);
        string file;
        if (given is null)
        {
            string name = a.RawAssoc(master, "aname") is { Length: > 0 } n ? EmbossOutPath.FmtWord(n) : baseName;
            string ext = a.RawAssoc(master, "aextension") is { Length: > 0 } e ? e : program;
            file = name + "." + ext;
        }
        else
        {
            file = given;
        }

        if (file.Length == 0)
        {
            context.Error("Output file is required");
            return null;
        }

        string dir = EmbossOutPath.FixDirectory(a.RawAssoc(master, "adirectory"));
        if (dir.Length > 0 && !file.Contains('/', StringComparison.Ordinal))
            file = dir + file;
        if (!EmbossOutPath.DirectoryExists(EmbossReportArgs.ExistingDirectories(files), file))
        {
            context.Warn($"Failed to open output file '{file}' error:2 'No such file or directory'");
            context.Error($"Unable to open file '{file}' for output");
            return null;
        }

        // ajAlignValid: a width given below 10 is reset to 10 (no width: 50)
        int effectiveWidth = 0;
        if (width is { } w)
        {
            if (w < 10)
            {
                context.Warn($"Alignment width (-awidth={w.ToString(System.Globalization.CultureInfo.InvariantCulture)}) too narrow, reset to 10");
                w = 10;
            }

            effectiveWidth = w;
        }

        var options = new EmbossAlignOptions
        {
            Program = program,
            CommandLine = a.CommandLineEntries,
            RunDate = runDate,
            ReportFile = file,
            Width = effectiveWidth,
            Global = global ?? false,
            ShowAccession = accshow ?? false,
            ShowDescription = desshow ?? false,
            ShowUsa = usashow ?? false,
            Matrix = matrix,
            GapPenalty = gapOpen,
            ExtendPenalty = gapExtend,
        };
        return new Opened(new EmbossAlignmentWriter(format, options), file, format);
    }
}
