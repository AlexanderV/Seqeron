using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: a user-defined amino-acid matrix (--aamatrix; loadaamtx of io.c)

    internal sealed partial class MafftGlobals
    {
        [ThreadStatic]
        private static string? _ambientAaMatrix;

        /// <summary>The text of the file <c>_aamtx</c> of the run in progress on this thread (<c>mafft --aamatrix file</c>; the script copies the file there, the programs read it).</summary>
        internal string? AaMatrixText = _ambientAaMatrix;

        /// <summary>Makes <paramref name="text"/> the file <c>_aamtx</c> for every program that is started on this thread until the scope is disposed.</summary>
        internal static IDisposable UseAaMatrix(string? text)
        {
            var previous = _ambientAaMatrix;
            _ambientAaMatrix = text;
            return new Scope(() => _ambientAaMatrix = previous);
        }

        private sealed class Scope(Action restore) : IDisposable
        {
            public void Dispose() => restore();
        }

        private const string AaMatrixExample = "Format error in aa matrix\n# Example:\n# comment\n   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n"
            + "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0\nR -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3\n...\n"
            + "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4\nfrequency 0.07 0.05 0.04 0.05 0.02 .. \n# Example end\nOnly the lower half is loaded\nThe last line (frequency) is optional.\n";

        /// <summary>
        /// <c>loadaamtx</c> (<c>io.c</c>) with <c>BLOSUMmtx( -1 ... )</c>: the first line that does not start with '#' lists the 20 letters (any order, the first 20 letters of the line), then one row per letter with the lower half of the matrix,
        /// then optionally a line "frequency" with 20 numbers or "norescale". The lines are read with <c>fgets( 999 )</c> into buffers that keep their old contents, as the program does.
        /// </summary>
        internal (double[] LowerTriangle, double[] Frequencies, bool Rescale) LoadAaMatrix(string? text)
        {
            if (text == null)
                throw new MafftProgramFailureException("Cannot open the _aamtx file\n");
            const string aaorder = "ARNDCQEGHILKMFPSTWYV";
            int pos = 0;
            bool eof = false;
            void Fgets(char[] buf)
            {
                if (pos >= text.Length)
                {
                    eof = true;
                    return;
                }
                int n = 0;
                while (n < 998 && pos < text.Length)
                {
                    char c = text[pos++];
                    buf[n++] = c;
                    if (c == '\n')
                        break;
                }
                if (pos >= text.Length && (n == 0 || buf[n - 1] != '\n'))
                    eof = true;
                buf[n] = '\0';
            }
            var inorder = new char[1000];
            while (!eof)
            {
                Fgets(inorder);
                if (inorder[0] != '#')
                    break;
            }
            int p1 = 0;
            for (int p2 = 0; inorder[p2] != '\0'; p2++)
            {
                if (inorder[p2] is >= 'a' and <= 'z' or >= 'A' and <= 'Z')
                    inorder[p1++] = char.ToUpperInvariant(inorder[p2]);
            }
            inorder[20] = '\0';
            string order = Cstr(inorder, 0);
            var map = new int[20];
            for (int i = 0; i < 20; i++)
            {
                int k = order.IndexOf(aaorder[i]);
                if (k < 0)
                    throw new MafftProgramFailureException(aaorder[i] + ": not found in the first 20 letters.\n" + AaMatrixExample);
                map[i] = k;
            }
            var raw = new double[21][];
            for (int i = 0; i < 21; i++)
                raw[i] = new double[20];
            var line = new char[1000];
            int row = 0;
            while (!eof)
            {
                Fgets(line);
                if (line[0] == '#')
                    continue;
                int ptr = 0;
                for (int j = 0; j <= row; j++)
                {
                    while (!(line[ptr] is >= '0' and <= '9') && line[ptr] != '-' && line[ptr] != '.')
                    {
                        ptr++;
                        if (ptr >= line.Length)
                            throw new NotSupportedException("The matrix line has fewer numbers than the program expects: it reads behind its buffer (undefined behaviour).");
                    }
                    raw[row][j] = CAtof(Cstr(line, ptr));
                    ptr = IndexOf(line, ptr, ' ');
                    if (ptr < 0)
                    {
                        if (j < row)
                            throw new MafftProgramFailureException(AaMatrixExample);
                        break;
                    }
                }
                row++;
                if (row > 19)
                    break;
            }
            bool rescale = true;
            for (int i = 0; i < 20; i++)
                raw[20][i] = -1.0;
            string key = "";
            while (!eof)
            {
                Fgets(line);
                string[] tokens = Cstr(line, 0).Split(new[] { ' ', '\t', '\n', '\v', '\f', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                if (tokens.Length > 0)
                    key = tokens[0];
                if (key == "norescale")
                {
                    Report("no rescale\n");
                    rescale = false;
                    break;
                }
                if (key == "frequency")
                {
                    int ptr = 0;
                    for (int j = 0; j < 20; j++)
                    {
                        while (!(line[ptr] is >= '0' and <= '9') && line[ptr] != '-' && line[ptr] != '.')
                        {
                            ptr++;
                            if (ptr >= line.Length)
                                throw new NotSupportedException("The frequency line has fewer numbers than the program expects: it reads behind its buffer (undefined behaviour).");
                        }
                        raw[20][j] = CAtof(Cstr(line, ptr));
                        ptr = IndexOf(line, ptr, ' ');
                        if (ptr < 0)
                        {
                            if (j < 19)
                                throw new MafftProgramFailureException(AaMatrixExample);
                            break;
                        }
                    }
                    break;
                }
            }
            var tri = new double[210];
            int count = 0;
            for (int i = 0; i < 20; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    int ii, jj;
                    if (i != j)
                    {
                        ii = Math.Max(map[i], map[j]);
                        jj = Math.Min(map[i], map[j]);
                    }
                    else
                    {
                        ii = jj = map[i];
                    }
                    tri[count++] = raw[ii][jj];
                }
            }
            var freq = new double[20];
            double first = raw[20][map[0]];
            if (first != -1.0)
            {
                for (int i = 0; i < 20; i++)
                    freq[i] = raw[20][map[i]];
                double av = 0.0;
                for (int i = 0; i < 20; i++)
                    av += freq[i];
                for (int i = 0; i < 20; i++)
                    freq[i] /= av;
            }
            else
            {
                Array.Copy(MafftTables.BlosumFrequencies, freq, 20);
            }
            Report("inorder = " + order + "\n");
            return (tri, freq, rescale);
        }

        private static string Cstr(char[] buf, int start)
        {
            int e = start;
            while (e < buf.Length && buf[e] != '\0')
                e++;
            return new string(buf, start, e - start);
        }

        private static int IndexOf(char[] buf, int start, char c)
        {
            for (int i = start; i < buf.Length && buf[i] != '\0'; i++)
                if (buf[i] == c)
                    return i;
            return -1;
        }
    }

    #endregion
}
