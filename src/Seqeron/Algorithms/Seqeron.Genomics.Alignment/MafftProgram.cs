#pragma warning disable IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The result of one run of the <c>mafft</c> script: the return code, the text on standard output (the alignment in the output format) and on standard error (the messages of the script and of the binaries it runs), and the files the run wrote (<c>--out</c>, <c>input.tree</c>, <c>input.hat2</c>, <c>input.map</c>, <c>input.notused</c>, <c>--progress</c> file).</summary>
public sealed record MafftProgramResult(int ExitCode, string Stdout, string Stderr, IReadOnlyDictionary<string, string> Files);

/// <summary>
/// MAFFT 7.505 as a command line: the shell script <c>mafft</c> (<c>core/mafft.tmpl</c>) ported option by option - the parsing of the arguments (every option, the checks of its argument, "Unknown option", the usage text), the derived variables
/// (<c>--auto</c>, <c>--similaritylevel</c>, <c>--allowshift</c>, the number of tree builds ...), the pipeline of binaries (<c>replaceu</c>, <c>disttbfast</c> / <c>tbfast</c> + <c>pairlocalalign</c> / <c>splittbfast</c>, <c>dndpre</c>, <c>dvtditr</c>,
/// <c>restoreu</c>, <c>f2cl</c>) on the ported programs, and the files the script copies next to the input. The arguments are those of the program, <b>the last argument is the input file</b>; every file the arguments name
/// (the input, <c>--add</c>, <c>--merge</c>, <c>--aamatrix</c>, ...) is given as text in <c>files</c> (the program reads the file system). The messages equal the installed Debian <c>mafft</c> 7.505-1 byte for byte
/// (the oracle generators <c>scripts/alignment/generate_mafft_command_oracle.py</c> ...; B09.md F315 ff.); the one text that depends on the machine - the number of physical cores written for <c>--thread -1</c> - is that of this machine.
/// </summary>
public static class MafftProgram
{
    /// <summary>The path the script writes in its messages for <c>$0</c> ("/usr/bin/mafft": the Debian installation).</summary>
    public const string ProgramPath = "/usr/bin/mafft";

    /// <summary>Runs <c>mafft</c> with the given arguments (options followed by the input file).</summary>
    /// <param name="arguments">The command-line arguments, one per element (the input file is the last one).</param>
    /// <param name="files">The text of the files the arguments name: file name -> content.</param>
    /// <param name="stdin">The text of standard input (only read when the input file is <c>-</c>).</param>
    /// <exception cref="NotSupportedException">The arguments use an option whose binary is not ported yet (the message names it).</exception>
    /// <param name="cancellationToken">Cancels the run at the next row / iteration of the dynamic programs.</param>
    public static MafftProgramResult Run(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string>? files = null, string? stdin = null, CancellationToken cancellationToken = default)
    {
        SequenceAligner.RequireProgramCall(arguments?.ToList(), files);
        using var cancelScope = SequenceAligner.BeginCancellationScope(cancellationToken);
        var run = new SequenceAligner.MafftScriptRun(arguments!, files ?? new Dictionary<string, string>(), stdin);
        return run.Execute();
    }

    /// <summary>The words of a command line as the shell splits them (blanks separate words; single and double quotes group; a backslash escapes the next character outside single quotes).</summary>
    public static IReadOnlyList<string> SplitCommandLine(string commandLine)
    {
        ArgumentNullException.ThrowIfNull(commandLine);
        var words = new List<string>();
        var cur = new StringBuilder();
        bool inWord = false;
        for (int i = 0; i < commandLine.Length; i++)
        {
            char c = commandLine[i];
            if (c == '\'')
            {
                inWord = true;
                int j = commandLine.IndexOf('\'', i + 1);
                if (j < 0)
                    throw new ArgumentException("Unterminated single quote in the command line.", nameof(commandLine));
                cur.Append(commandLine, i + 1, j - i - 1);
                i = j;
            }
            else if (c == '"')
            {
                inWord = true;
                int j = i + 1;
                while (j < commandLine.Length && commandLine[j] != '"')
                {
                    if (commandLine[j] == '\\' && j + 1 < commandLine.Length && commandLine[j + 1] is '"' or '\\' or '$' or '`')
                        j++;
                    cur.Append(commandLine[j]);
                    j++;
                }
                if (j >= commandLine.Length)
                    throw new ArgumentException("Unterminated double quote in the command line.", nameof(commandLine));
                i = j;
            }
            else if (c == '\\' && i + 1 < commandLine.Length)
            {
                inWord = true;
                cur.Append(commandLine[++i]);
            }
            else if (char.IsWhiteSpace(c))
            {
                if (inWord)
                {
                    words.Add(cur.ToString());
                    cur.Clear();
                    inWord = false;
                }
            }
            else
            {
                inWord = true;
                cur.Append(c);
            }
        }
        if (inWord)
            words.Add(cur.ToString());
        return words;
    }
}
