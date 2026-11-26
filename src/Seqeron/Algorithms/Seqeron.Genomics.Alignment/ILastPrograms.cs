namespace Seqeron.Genomics.Alignment;

/// <summary>The outcome of a run of an external program on a working directory: standard output, standard error and the exit status.</summary>
/// <param name="StandardOutput">The text written to the standard output.</param>
/// <param name="StandardError">The text written to the standard error.</param>
/// <param name="ExitCode">The exit status of the program.</param>
public sealed record LastProcessResult(string StandardOutput, string StandardError, int ExitCode);

/// <summary>
/// The programs <c>lastdb</c> and <c>lastal</c> of LAST that <c>mafft --lastpair</c> and <c>--lastmultipair</c> start (<c>pairlocalalign -R</c> / <c>-r</c>). LAST is licensed under the GPL: it is
/// not part of this (MIT-licensed) assembly, which only knows this interface; an implementation is supplied through <see cref="MafftInstallation.Last"/> (the GPL-licensed port of LAST lives in a separate assembly).
/// Without it the script stops with the "Install LAST" message of the installation that has no LAST.
/// </summary>
public interface ILastPrograms
{
    /// <summary>Runs one program of LAST in a working directory.</summary>
    /// <param name="programPath">The path of the program as the caller starts it (for example <c>/opt/last/bin/lastal</c>): the program writes it in front of its error messages.</param>
    /// <param name="arguments">The arguments (without the program name).</param>
    /// <param name="workingDirectory">The files of the working directory (name -> text); the program reads and writes them.</param>
    /// <param name="standardInput">The text of the standard input.</param>
    LastProcessResult Run(string programPath, IReadOnlyList<string> arguments, IDictionary<string, string> workingDirectory, string standardInput = "");
}
