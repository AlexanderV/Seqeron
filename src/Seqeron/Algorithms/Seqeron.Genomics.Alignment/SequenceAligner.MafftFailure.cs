using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    /// <summary>
    /// A MAFFT program stops with <c>exit( 1 )</c> after it wrote its message to standard error (an unknown letter, a sequence that is too short for the 6-mer tables ...): <see cref="Text"/> is the whole standard error of the run
    /// up to and including that message. It is an <see cref="ArgumentException"/> for the callers of the library; the script layer (<see cref="MafftProgram"/>) turns it into the exit status of the program.
    /// </summary>
    internal sealed class MafftProgramFailureException : ArgumentException
    {
        internal MafftProgramFailureException(string text, int exitCode = 1, string? summary = null)
            : base(summary ?? Summary(text))
        {
            Text = text;
            ExitCode = exitCode;
        }

        /// <summary>The standard error of the run up to the end of the message of the failing program.</summary>
        internal string Text { get; }

        internal int ExitCode { get; }

        private static string Summary(string text)
        {
            var lines = text.Split('\n', StringSplitOptions.RemoveEmptyEntries);
            return lines.Length == 0 ? "MAFFT stopped." : "MAFFT stopped: " + lines[^1].Trim();
        }
    }

    internal sealed partial class MafftRun
    {
        /// <summary><c>seqcheck</c> of <c>mltaln9.c</c>: the box that the programs write for a letter outside the alphabet, <paramref name="site"/> and <paramref name="sequence"/> are 1-based.</summary>
        internal static string UnknownAlphabetBox(char c, int site, int sequence) =>
            "========================================================================= \n" +
            "========================================================================= \n" +
            "=== \n" +
            "=== Alphabet '" + c + "' is unknown.\n" +
            "=== Please check site " + site.ToString(CultureInfo.InvariantCulture) + " in sequence " + sequence.ToString(CultureInfo.InvariantCulture) + ".\n" +
            "=== \n" +
            "=== To make an alignment that has unusual characters (U, @, #, etc), try  \n" +
            "=== % mafft --anysymbol input > output\n" +
            "=== \n" +
            "========================================================================= \n" +
            "========================================================================= \n";

        /// <summary><c>seqcheck</c> + the message of the calling program: null when every letter is in the alphabet.</summary>
        internal MafftProgramFailureException? CheckAlphabet(IReadOnlyList<string> seq, string illegalWord = "Illegal")
        {
            for (int i = 0; i < seq.Count; i++)
            {
                string s = seq[i];
                for (int k = 0; k < s.Length; k++)
                {
                    char c = s[k];
                    if (c >= G.AminoN.Length || G.AminoN[c] == -1)
                        return new MafftProgramFailureException(UnknownAlphabetBox(c, k + 1, i + 1) + illegalWord + " character " + c + "\n");
                }
            }
            return null;
        }
    }

    /// <summary>The log of a run followed by the message of a failure, as one <see cref="MafftProgramFailureException"/>.</summary>
    internal static T GuardRun<T>(MafftRun run, StringBuilder log, Func<T> body)
    {
        try
        {
            return body();
        }
        catch (MafftProgramFailureException e)
        {
            log.Append(run.Log).Append(e.Text);
            throw new MafftProgramFailureException(log.ToString(), e.ExitCode);
        }
    }
}
