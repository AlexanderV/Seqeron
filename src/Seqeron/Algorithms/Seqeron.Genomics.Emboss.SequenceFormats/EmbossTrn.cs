// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S3358 // the library functions keep their shape
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The translation tables of the EMBOSS 6.6.0 AJAX library (<c>ajtrn.c</c>, <c>AjPTrn</c>; B09 WP61, F280): <c>ajTrnNewI</c> reads the genetic code file <c>EGC.n</c> (<see cref="ReadFile"/> = <c>ajTrnReadFile</c>),
/// <c>trnComplete</c> fills the ambiguity codes of the table (the wobble of the third base and, for every amino acid, ONE "most ambiguous codon" built by OR-ing the bases of its codons: not the translation of the ambiguity code
/// by the genetic code), the translation of a string (<c>ajTrnSeqC</c>, <c>ajTrnSeqInitC</c>, <c>ajTrnSeqRevC</c>, <c>ajTrnSeqRevinitC</c>, <c>ajTrnSeqAltRevC</c>, <c>ajTrnSeqAltRevinitC</c>, <c>ajTrnSeqFrameC</c>,
/// <c>ajTrnSeqFrameInitC</c>), the dangling bases (<c>ajTrnSeqDangleC</c>) and the type of a codon (<c>ajTrnCodonstrTypeC</c>: start 1, stop -1). Read as ajtrn.c of EMBOSS-6.6.0 (raw.githubusercontent.com mirrors).
/// </summary>
public sealed partial class EmbossTrn
{
    // the conversion of a base character to the index of the translation arrays (ajtrn.c trnconv / trncomp: A C G T/U M R W S Y K V H D B N, anything else is 14)
    private static readonly int[] Conv = MakeConv(false);
    private static readonly int[] Comp = MakeConv(true);

    private static int[] MakeConv(bool complement)
    {
        var t = new int[128];
        Array.Fill(t, 14);
        // base, index of the base, index of its complement
        (char Base, int Fwd, int Rev)[] map =
        [
            ('A', 0, 3), ('C', 1, 2), ('G', 2, 1), ('T', 3, 0), ('U', 3, 0), ('M', 4, 9), ('R', 5, 8), ('W', 6, 6), ('S', 7, 7), ('Y', 8, 5), ('K', 9, 4), ('V', 10, 13), ('H', 11, 12), ('D', 12, 11), ('B', 13, 10),
        ];
        foreach (var (b, f, r) in map)
        {
            t[b] = complement ? r : f;
            t[char.ToLowerInvariant(b)] = complement ? r : f;
        }

        return t;
    }

    private static int C(char c) => c < 128 ? Conv[c] : 14;

    private static int K(char c) => c < 128 ? Comp[c] : 14;

    private readonly char[,,] _gc = new char[16, 16, 16];
    private readonly char[,,] _starts = new char[16, 16, 16];

    /// <summary>The title line of the genetic code file.</summary>
    public string Title { get; private set; } = "";

    /// <summary>The name of the genetic code file (<c>EGC.n</c>).</summary>
    public string FileName { get; private set; } = "";

    private static readonly Dictionary<int, EmbossTrn> Cache = [];

    /// <summary><c>ajTrnNewI</c>: the table of the file <c>EGC.n</c> (a table is read once and not changed afterwards).</summary>
    /// <exception cref="EmbossFatalException">There is no file <c>EGC.n</c> (<c>ajFatal</c>).</exception>
    public static EmbossTrn New(int table)
    {
        lock (Cache)
        {
            if (Cache.TryGetValue(table, out var t))
                return t;
            string? text = null;
            foreach (var (n, body) in DataFiles)
            {
                if (n == table)
                    text = body;
            }

            if (text is null)
                throw new EmbossFatalException($"Translation table file 'EGC.{table}' not found\n");
            var trn = new EmbossTrn { FileName = "EGC." + table.ToString(System.Globalization.CultureInfo.InvariantCulture) };
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 16; j++)
                    for (int k = 0; k < 16; k++)
                    {
                        trn._gc[i, j, k] = 'X';
                        trn._starts[i, j, k] = '-';
                    }

            trn.ReadFile(text);
            Cache[table] = trn;
            return trn;
        }
    }

    /// <summary><c>ajTrnName</c>: the name of the genetic code of the table <c>EGC.index</c> lists ("unknown" when it does not).</summary>
    public static string Name(int table)
    {
        foreach (string raw in IndexText.Split('\n'))
        {
            string line = raw.Trim();
            if (line.Length == 0 || line[0] == '#')
                continue;
            int sp = line.IndexOf(' ');
            string tok1 = sp < 0 ? line : line[..sp];
            string tok2 = sp < 0 ? "" : line[(sp + 1)..];
            if (tok1 == table.ToString(System.Globalization.CultureInfo.InvariantCulture))
                return tok2;
        }

        return "unknown";
    }

    // ajTrnReadFile (the text has no comment lines; trnNoComment also cuts a comment behind text)
    private void ReadFile(string text)
    {
        const string white = " \t\n\r!@#$%^&()_+=|\\~`{[}]:;\"'<,>.?/";
        var lines = new List<string>();
        foreach (string raw in text.Split('\n'))
        {
            string l = raw.Trim(' ', '\t', '\r', '\n', '\f', '\v');
            int hash = l.IndexOf('#');
            if (hash >= 0)
                l = l[..hash];
            lines.Add(l);
        }

        int pos = 0;
        while (pos < lines.Count)
        {
            string l = lines[pos++];
            if (l.Length > 0)
            {
                if (!l.Contains("Genetic Code", StringComparison.Ordinal))
                    throw new EmbossFatalException($"The file '{FileName}' is not a valid Genetic Code file.\nThe 'Genetic Code' line was not found.");
                break;
            }
        }

        while (pos < lines.Count)
        {
            string l = lines[pos++];
            if (l.Length > 0)
            {
                Title = l;
                break;
            }
        }

        var sb = new StringBuilder();
        while (pos < lines.Count)
        {
            string l = lines[pos++];
            if (l.Length > 0)
                sb.Append(l).Append(' ');
        }

        string[] tokens = sb.ToString().Split(white.ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
        // AAs <aa> Starts <starts> Base1 <b1> Base2 <b2> Base3 <b3>
        string aa = tokens[1], starts = tokens[3], base1 = tokens[5], base2 = tokens[7], base3 = tokens[9];
        for (int i = 0; i < starts.Length; i++)
            _starts[C(base1[i]), C(base2[i]), C(base3[i])] = starts[i];
        for (int i = 0; i < aa.Length; i++)
            _gc[C(base1[i]), C(base2[i]), C(base3[i])] = aa[i];
        Complete();
    }

    // trnComplete
    private void Complete()
    {
        const string bases = "ACGT";
        int[] code = [1, 2, 4, 8];
        int[] ambigcodes = [1, 2, 4, 8, 3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15];
        int[] trncodes = [14, 0, 1, 4, 2, 5, 7, 10, 3, 6, 8, 11, 9, 12, 13, 14];
        _ = bases;
        int[] codonval = new int[4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int kfirst = 0;
                while (kfirst < 4)
                {
                    codonval[2] = 0;
                    char aa = _gc[i, j, kfirst];
                    int newkfirst = 4;
                    for (int k = kfirst; k < 4; k++)
                    {
                        if (_gc[i, j, k] == aa)
                            codonval[2] |= code[k];
                        else if (newkfirst == 4)
                            newkfirst = k;
                    }

                    for (int k = 4; k < 15; k++)
                    {
                        int kk = ambigcodes[k];
                        if ((kk | codonval[2]) == codonval[2] && (kk & codonval[2]) != 0)
                            _gc[i, j, k] = aa;
                    }

                    kfirst = newkfirst;
                }
            }
        }

        char newaa = 'X';
        int ifirst = 0, jfirst = 0, kfirst2 = 0, jj = 0, kk2 = 0;
        var aalist = new StringBuilder();
        while (newaa != '\0')
        {
            newaa = '\0';
            char aa = _gc[ifirst, jfirst, kfirst2];
            codonval[0] = codonval[1] = codonval[2] = 0;
            for (int i = ifirst; i < 4; i++)
            {
                for (int j = jfirst; j < 4; j++)
                {
                    for (int k = kfirst2; k < 4; k++)
                    {
                        if (_gc[i, j, k] == aa)
                        {
                            codonval[0] |= code[i];
                            codonval[1] |= code[j];
                            codonval[2] |= code[k];
                        }
                        else if (newaa == '\0' && !aalist.ToString().Contains(_gc[i, j, k]))
                        {
                            newaa = _gc[i, j, k];
                            ifirst = i;
                            jj = j;
                            kk2 = k;
                            aalist.Append(newaa);
                        }
                    }

                    kfirst2 = 0;
                }

                jfirst = 0;
            }

            jfirst = jj;
            kfirst2 = kk2;
            _gc[trncodes[codonval[0]], trncodes[codonval[1]], trncodes[codonval[2]]] = aa;
        }
    }

    // ---- the translation of codons ------------------------------------------------------------------------------------------

    /// <summary><c>ajTrnCodonK</c> / <c>ajTrnCodonC</c>: the residue of three bases.</summary>
    public char Codon(char b1, char b2, char b3) => _gc[C(b1), C(b2), C(b3)];

    /// <summary><c>ajTrnCodonstrTypeC</c>: the residue of the codon at <paramref name="pos"/> and its type (1 = start, -1 = stop, 0).</summary>
    public int CodonType(string s, int pos, out char aa)
    {
        int a = C(s[pos]), b = C(s[pos + 1]), c = C(s[pos + 2]);
        aa = _gc[a, b, c];
        if (_starts[a, b, c] == 'M')
            return 1;
        return aa == '*' ? -1 : 0;
    }

    /// <summary><c>ajTrnSeqC</c>: the codons of the first <paramref name="len"/> bases of <paramref name="s"/> from <paramref name="offset"/>, appended.</summary>
    public void SeqC(string s, int offset, int len, StringBuilder pep)
    {
        int lenmod3 = len - len % 3;
        for (int i = 0; i < lenmod3; i += 3)
            pep.Append(_gc[C(s[offset + i]), C(s[offset + i + 1]), C(s[offset + i + 2])]);
    }

    /// <summary><c>ajTrnSeqInitC</c>: as <see cref="SeqC"/>, a start codon of the first position is an M.</summary>
    public void SeqInitC(string s, int offset, int len, StringBuilder pep)
    {
        int lenmod3 = len - len % 3;
        int istart = 0;
        if (lenmod3 > 2 && _starts[C(s[offset]), C(s[offset + 1]), C(s[offset + 2])] == 'M')
        {
            istart = 3;
            pep.Append('M');
        }

        for (int i = istart; i < lenmod3; i += 3)
            pep.Append(_gc[C(s[offset + i]), C(s[offset + i + 1]), C(s[offset + i + 2])]);
    }

    /// <summary><c>ajTrnSeqRevC</c>: the reverse complement of the codons from the end of the first <paramref name="len"/> bases, in the order of the codons read from the end.</summary>
    public void SeqRevC(string s, int offset, int len, StringBuilder pep)
    {
        int end = len / 3 * 3 - 1;
        for (int i = end; i > 1; i -= 3)
            pep.Append(_gc[K(s[offset + i]), K(s[offset + i - 1]), K(s[offset + i - 2])]);
    }

    /// <summary><c>ajTrnSeqRevinitC</c>.</summary>
    public void SeqRevinitC(string s, int offset, int len, StringBuilder pep)
    {
        int end = len / 3 * 3 - 1;
        int iend = end;
        int cp = offset + end;
        if (iend > 2 && _starts[K(s[cp]), K(s[cp - 1]), K(s[cp - 2])] == 'M')
        {
            cp -= 3;
            iend -= 3;
            pep.Append('M');
        }

        for (int i = iend; i > 1; i -= 3)
        {
            pep.Append(_gc[K(s[cp]), K(s[cp - 1]), K(s[cp - 2])]);
            cp -= 3;
        }
    }

    /// <summary><c>ajTrnSeqAltRevC</c>.</summary>
    public void SeqAltRevC(string s, int len, StringBuilder pep)
    {
        for (int i = len - 1; i > 1; i -= 3)
            pep.Append(_gc[K(s[i]), K(s[i - 1]), K(s[i - 2])]);
    }

    /// <summary><c>ajTrnSeqAltRevinitC</c>.</summary>
    public void SeqAltRevinitC(string s, int len, StringBuilder pep)
    {
        int iend = len - 1;
        if (iend > 1 && _starts[K(s[iend]), K(s[iend - 1]), K(s[iend - 2])] == 'M')
        {
            iend -= 3;
            pep.Append('M');
        }

        for (int i = iend; i > 1; i -= 3)
            pep.Append(_gc[K(s[i]), K(s[i - 1]), K(s[i - 2])]);
    }

    /// <summary><c>ajTrnSeqFrameC</c> (the frames 1..3, -1..-3 and the alternative reverse frames -4..-6; 4, 5, 6 are -1, -2, -3).</summary>
    /// <exception cref="EmbossFatalException">An invalid frame (<c>ajFatal</c>).</exception>
    public void SeqFrameC(string s, int len, int frame, StringBuilder pep)
    {
        if (frame > 3)
            frame = -frame + 3;
        if (frame is >= 1 and <= 3)
            SeqC(s, frame - 1, len - frame + 1, pep);
        else if (frame is >= -3 and <= -1)
            SeqRevC(s, -frame - 1, len + frame + 1, pep);
        else if (frame is >= -6 and <= -4)
            SeqAltRevC(s, len + frame + 4, pep);
        else
            throw new EmbossFatalException($"Invalid frame '{frame}' in ajTrnSeqFrameC()\n");
    }

    /// <summary><c>ajTrnSeqFrameInitC</c>.</summary>
    public void SeqFrameInitC(string s, int len, int frame, StringBuilder pep)
    {
        if (frame > 3)
            frame = -frame + 3;
        if (frame is >= 1 and <= 3)
            SeqInitC(s, frame - 1, len - frame + 1, pep);
        else if (frame is >= -3 and <= -1)
            SeqRevinitC(s, -frame - 1, len + frame + 1, pep);
        else if (frame is >= -6 and <= -4)
            SeqAltRevinitC(s, len + frame + 4, pep);
        else
            throw new EmbossFatalException($"Invalid frame '{frame}' in ajTrnSeqFrameC()\n");
    }

    // the character of a C string at an index (the terminating NUL past the end)
    private static char At(string s, int i) => i >= 0 && i < s.Length ? s[i] : '\0';

    /// <summary><c>ajTrnSeqDangleS</c>: translates the two bases that are left over at the end of the frame as the codon b1 b2 N (the reverse sense: the complements of the first two bases), a single base as X; returns the number of dangling bases.</summary>
    public int SeqDangle(string seq, int frame, StringBuilder pep)
    {
        int end = 0;
        int dangle;
        int len = seq.Length;
        if (frame > 3)
            frame = -frame + 3;
        if (frame > 0)
        {
            end = frame + (len - frame + 1) / 3 * 3 - 1;
            dangle = len - end;
        }
        else if (frame <= -4)
            dangle = (len + frame + 4) % 3;
        else
            dangle = -frame - 1;
        if (dangle == 2)
        {
            if (frame is >= 1 and <= 3)
                pep.Append(_gc[C(seq[end]), C(seq[end + 1]), C('\0')]);
            else
                pep.Append(_gc[K(At(seq, 1)), K(At(seq, 0)), K('\0')]);
        }
        else if (dangle == 1)
        {
            pep.Append('X');
        }

        return dangle;
    }

    /// <summary>The frame number that <c>ajTrnNewPep</c> appends to the name of the translation (<c>name_frame</c>): the frames 4..6 and -4..-6 are numbered 4, 5, 6; 0 appends nothing.</summary>
    public static int PepFrameNumber(int frame)
    {
        if (frame == 0)
            return 0;
        if (frame < -3)
            frame += 3;
        if (frame < 0)
            frame = -frame + 3;
        return frame;
    }

    /// <summary><c>ajTrnSeqOrig</c> (<paramref name="init"/> false) / <c>ajTrnSeqInit</c> (true): the translation of a frame with the dangling bases translated; the reverse frames 4, 5, 6 are reversed for display under the sequence.</summary>
    public string SeqOrig(string seq, int frame, bool init)
    {
        var trn = new StringBuilder();
        if (init)
            SeqFrameInitC(seq, seq.Length, frame, trn);
        else
            SeqFrameC(seq, seq.Length, frame, trn);
        SeqDangle(seq, frame, trn);
        if (frame > 3)
        {
            var chars = trn.ToString().ToCharArray();
            Array.Reverse(chars);
            return new string(chars);
        }

        return trn.ToString();
    }
}
