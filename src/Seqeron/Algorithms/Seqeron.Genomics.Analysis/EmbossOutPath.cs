using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The output file names of the sequence outputs of EMBOSS 6.6.0 (A5-8c3h, B09 WP59): <c>acdOutDirectory</c> / <c>ajDirnameFix</c> (the directory of <c>-osdirectory</c> gets a trailing slash) and
/// <c>ajFileNewOutNamePathS</c> (the directory is put in front of a file name that has no slash; <c>stdout</c> and <c>stderr</c> are never in a directory).
/// </summary>
internal static class EmbossOutPath
{
    /// <summary><c>acdOutDirectory</c>: empty stays empty, any other directory ends with a slash.</summary>
    internal static string FixDirectory(string? directory)
    {
        if (string.IsNullOrEmpty(directory))
            return "";
        return directory.EndsWith('/') ? directory : directory + "/";
    }

    /// <summary><c>ajFileNewOutNamePathS</c>: the path of the file that is opened.</summary>
    internal static string Combine(string directory, string name)
        => directory.Length == 0 || name is "stdout" or "stderr" || name.Contains('/', StringComparison.Ordinal) ? name : directory + name;

    /// <summary>
    /// Whether the directory of a file exists: <paramref name="existing"/> lists the directories (without trailing slash) of the place the program runs in; null means every directory exists.
    /// </summary>
    internal static bool DirectoryExists(IReadOnlyCollection<string>? existing, string path)
    {
        int slash = path.LastIndexOf('/');
        if (existing is null || slash < 0)
            return true;
        string dir = slash == 0 ? "/" : path[..slash];
        return dir is "." or "/" || existing.Contains(dir);
    }

    /// <summary>The messages the library prints when <c>ajSeqoutOpen</c> cannot open the output file of a non-single output (<c>acdSetSeqout</c> then dies: "Bad value for '-qualifier'").</summary>
    internal static void OpenFailed(EmbossFeatureReadContext context, string name, string directory, string usa)
    {
        context.Warn($"Failed to open output file '{name}' in directory '{directory}' error:2 'No such file or directory'");
        if (directory.Length > 0)
        {
            context.Error($"failed to open filename '{name}' in directory '{directory}'");
            context.Error($"Unable to write sequence to '{directory}{usa}'");
        }
        else
        {
            context.Error($"failed to open filename '{name}'");
            context.Error($"Unable to write sequence to '{usa}'");
        }
    }

    /// <summary>
    /// <c>seqFileReopen</c> of <c>-ossingle</c> in a directory that does not exist: the library warns "Failed to open output file" for every sequence and the writers that print through the file with
    /// <c>ajFmtPrintF</c> die by a segmentation fault, but the DAS and SAM writers write nothing and the program ends normally (measured on 6.6.0: every writer of the library).
    /// </summary>
    internal static bool SingleOpenSurvives(EmbossSeqOutFormat format) => format.Writer is "Dasdna" or "Dasseq" or "Sam";

    /// <summary>The warning of a file of <c>-ossingle</c> that cannot be opened.</summary>
    internal static void SingleOpenFailed(EmbossFeatureReadContext context, string file, string directory)
        => context.Warn($"Failed to open output file '{file}' in directory '{directory}' error:2 'No such file or directory'");

    /// <summary>
    /// <c>acdSetReport</c>: the report file of a program. Without a file name (<paramref name="given"/> null) the default is <c>BASE.EXT</c> (<c>acdOutFilename</c>: <c>-rname</c> with <c>ajStrFmtWord</c> applied, else the saved name of the first input - the name of its first
    /// sequence, lower case -, else <c>outfile</c>; <c>-rextension</c>, else the name of the program); <c>-rdirectory</c> is put in front of a name that has no directory (<c>ajFilenameReplacePathS</c>: unlike the sequence outputs also in front of <c>stdout</c>). A file that is
    /// in a directory that does not exist cannot be opened and the program dies.
    /// </summary>
    /// <returns>The path of the report file; null when the program dies (the messages are in <paramref name="context"/>).</returns>
    internal static string? OpenReport(EmbossFeatureReadContext context, string program, string? given, EmbossReportQualifiers? rq, string saved, IReadOnlyCollection<string>? existingDirectories)
    {
        string file;
        if (given is null)
        {
            string name = rq?.Name is { Length: > 0 } n ? FmtWord(n) : saved;
            if (name.Length == 0)
                name = "outfile";
            string ext = rq?.Extension is { Length: > 0 } e ? e : program;
            file = ext.Length > 0 ? name + "." + ext : name;
        }
        else
        {
            file = given;
        }

        if (file.Length == 0)
        {
            context.Error("Report file is required");
            return null;
        }

        string dir = FixDirectory(rq?.Directory);
        if (dir.Length > 0 && !file.Contains('/', StringComparison.Ordinal))
            file = dir + file;
        if (!DirectoryExists(existingDirectories, file))
        {
            context.Warn($"Failed to open output file '{file}' error:2 'No such file or directory'");
            context.Error($"Unable to open report file '{file}'");
            return null;
        }

        return file;
    }

    // ajStrFmtWord: a string with blanks loses the excess white space and its blanks become underscores; the part before the last colon is cut
    internal static string FmtWord(string name)
    {
        string t = name;
        if (t.Length == 0 || t.Any(char.IsWhiteSpace))
            t = Regex.Replace(t.Trim(), "\\s+", " ", RegexOptions.CultureInvariant).Replace(' ', '_');
        if (t.EndsWith(':'))
            t = t[..^1];
        int colon = t.LastIndexOf(':');
        return colon >= 0 ? t[(colon + 1)..] : t;
    }

    // acdInFileSave (the name of the first sequence, else the USA) + ajFilenameTrimAll + ajStrFmtLower, then ajStrFmtWord of acdOutFilename
    internal static string DefaultBaseName(string seqName, string usa)
    {
        string name = TrimAll(seqName.Length > 0 ? seqName : usa).ToLowerInvariant();
        int blank = name.IndexOfAny([' ', '\t', '\n', '\r']);
        return blank >= 0 ? name[..blank] : name;
    }

    // ajFilenameTrimAll
    private static string TrimAll(string name)
    {
        var m = Regex.Match(name, ":([A-Za-z0-9_-]+)$", RegexOptions.CultureInvariant);
        if (m.Success)
            return m.Groups[1].Value;
        m = Regex.Match(name, "([A-Za-z0-9_-]+)[.][A-Za-z0-9_-]+$", RegexOptions.CultureInvariant);
        if (m.Success)
            return m.Groups[1].Value;
        m = Regex.Match(name, "([A-Za-z0-9_-]+)[^A-Za-z0-9_-]*$", RegexOptions.CultureInvariant);
        return m.Success ? m.Groups[1].Value : "";
    }
}
