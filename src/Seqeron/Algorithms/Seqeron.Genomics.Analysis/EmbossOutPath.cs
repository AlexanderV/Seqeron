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
}
