using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>Translate_DNA</c> of MUMmer 3.23 (<c>translate.cc</c>, tables of <c>translate.hh</c>): the translation of a DNA sequence in one of the six reading frames (1 - 3 forward, 4 - 6 reverse complement), with the 16 symbol IUPAC
/// code (an ambiguous codon becomes <c>X</c> or one of the ambiguity codes <c>B Z</c> ...) and the universal genetic code.
/// </summary>
public static class MummerTranslate
{
    private const int BadPepChar = -1, DnaXn = 14;

    /// <summary>Translates <paramref name="a"/> (1-based: <c>a[1 .. len]</c>, as <c>Read_String</c> stores it) in <paramref name="frame"/> (1 - 6).</summary>
    /// <param name="a">The DNA, 1-based; element 0 is padding.</param>
    /// <param name="len">The length of the DNA.</param>
    /// <param name="frame">The reading frame, 1 - 6.</param>
    /// <param name="warnings">Receives the warnings of the program (<c>WARNING: Forcing unrecognized DNA char to N</c>, once per such character of a forward frame).</param>
    /// <returns>The amino acids, 1-based (<c>result[1 .. n]</c>, a zero byte after the last), and their number <c>n</c>; n = -1 for a frame outside 1 - 6.</returns>
    /// <remarks>
    /// A reverse frame complements the symbol index with <c>compdna [transdna [c]]</c>; for a character that is not an IUPAC nucleotide (<c>transdna</c> = BAD_PEP_CHAR = -1) the program reads <c>compdna [-1]</c>, the
    /// word in front of the array, which is 0 (= <c>DNA_A</c>) in the installed binary: such a character acts as <c>a</c> in a reverse frame and raises no warning (measured: <c>prepro</c> translates
    /// <c>ccc e ccc</c> in frames 4 - 6 like <c>ccc t ccc</c>, in frames 1 - 3 it warns). The port reproduces that value.
    /// </remarks>
    public static (byte[] Amino, int Length) Translate(byte[] a, long len, int frame, StringBuilder warnings)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(warnings);
        var t = new byte[len / 3 + 4];
        int n = 0;
        if (frame is >= 1 and <= 3)
        {
            for (long p = frame; p <= len - 2; p += 3)
            {
                int index = (Forward(a[p], warnings) << 8) + (Forward(a[p + 1], warnings) << 4) + Forward(a[p + 2], warnings);
                t[++n] = (byte)MummerPromerTables.Universal[index];
            }

            return (t, n);
        }

        if (frame is >= 4 and <= 6)
        {
            long f = frame - 3;
            for (long p = len - f + 1; p >= 3; p -= 3)
            {
                int index = (Reverse(a[p]) << 8) + (Reverse(a[p - 1]) << 4) + Reverse(a[p - 2]);
                t[++n] = (byte)MummerPromerTables.Universal[index];
            }

            return (t, n);
        }

        return (t, -1);
    }

    private static int Forward(byte c, StringBuilder warnings)
    {
        int d = MummerPromerTables.TransDna[c];
        if (d == BadPepChar)
        {
            warnings.Append("WARNING: Forcing unrecognized DNA char to N\n");
            d = DnaXn;
        }

        return d;
    }

    private static int Reverse(byte c)
    {
        int d = MummerPromerTables.TransDna[c];
        return d < 0 ? 0 : MummerPromerTables.CompDna[d];   // compdna [-1] reads 0 in the installed binary
    }
}

/// <summary>
/// <c>prepro</c> of MUMmer 3.23 (<c>src/tigr/prepro.cc</c>): the sequences of a FASTA file translated in all six reading frames for the amino-acid matching step of <c>promer</c>. A stretch of amino acids of at most
/// <c>-m</c> letters between two stop codons (or a sequence end) is masked (<c>X</c> in the reference, <c>O</c> in the query), a stop codon is written <c>J</c>. The reference output has one header
/// (<c>&gt;allcontigs file</c>) and all frames of all sequences concatenated, every frame ended by the mask character; the query output keeps a header <c>&gt;id.frame</c> per frame.
/// </summary>
public static class PromerPrepro
{
    private const int CharsPerLine = 60;
    private const int DefaultMaskLength = 10;
    private const byte TranslateMask = (byte)'X', ReferenceMask = (byte)'X', QueryMask = (byte)'O', StopMask = (byte)'J', StopChar = (byte)'*';

    /// <summary><c>prepro [-m len] -r|-q file</c>: returns the text the program writes to standard output.</summary>
    /// <param name="fasta">The FASTA text.</param>
    /// <param name="isReference"><c>-r</c> (the reference file; false = <c>-q</c>).</param>
    /// <param name="maskLength">The maximum book-end masking length <c>-m</c> (default 10; a negative value is ignored with a warning).</param>
    /// <param name="fileName">The file name printed in the reference header.</param>
    /// <param name="messages">What the program writes to standard error.</param>
    public static string Run(string fasta, bool isReference, int maskLength, string fileName, out string messages)
    {
        ArgumentNullException.ThrowIfNull(fasta);
        ArgumentNullException.ThrowIfNull(fileName);
        var err = new StringBuilder();
        if (maskLength < 0)
        {
            err.Append("WARNING: Invalid maximum mask length ").Append(maskLength.ToString(CultureInfo.InvariantCulture)).Append(", ignored\n");
            maskLength = DefaultMaskLength;
        }

        byte maskChar = isReference ? ReferenceMask : QueryMask;
        var sb = new StringBuilder();
        if (isReference)
            sb.Append(">allcontigs ").Append(fileName).Append('\n');
        var reader = new MummerFastaReader(fasta);
        int ct = 0;
        while (reader.Read(out string id, out byte[] a, out long lenA, err))
        {
            for (int frame = 1; frame <= 6; frame++)
            {
                if (!isReference)
                    sb.Append('>').Append(id).Append('.').Append(frame.ToString(CultureInfo.InvariantCulture)).Append('\n');
                var (t, lent) = MummerTranslate.Translate(a, lenA, frame, err);
                t[++lent] = maskChar;

                long lastIndex = 0;
                long i;
                for (i = 1; i <= lent; i++)
                {
                    if (maskChar != TranslateMask && t[i] == TranslateMask)
                        t[i] = maskChar;
                    else if (t[i] == StopChar)
                    {
                        t[i] = StopMask;
                        if (i - lastIndex - 1 <= maskLength)
                            Mask(t, maskChar, lastIndex + 1, i - 1);
                        lastIndex = i;
                    }
                }

                if (lent - lastIndex - 1 <= maskLength)
                    Mask(t, maskChar, lastIndex + 1, i - 1);

                for (i = 1; i <= lent; i++)
                {
                    sb.Append((char)t[i]);
                    if (++ct == CharsPerLine)
                    {
                        ct = 0;
                        sb.Append('\n');
                    }
                }

                if (!isReference)
                {
                    if (ct != 0)
                        sb.Append('\n');
                    ct = 0;
                }
            }
        }

        if (ct != 0)
            sb.Append('\n');
        messages = err.ToString();
        return sb.ToString();
    }

    private static void Mask(byte[] a, byte mask, long x, long y)
    {
        for (; x <= y; x++)
            a[x] = mask;
    }
}
