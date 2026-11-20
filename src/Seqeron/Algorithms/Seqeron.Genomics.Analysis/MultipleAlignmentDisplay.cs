using Seqeron.Genomics.Alignment;
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>What <see cref="MultipleAlignmentDisplay.ShowAlign"/> shows of the sequences (EMBOSS <c>showalign -show</c>).</summary>
public enum ShowAlignContent
{
    /// <summary><c>A</c>: all residues; with similarcase the dissimilar ones are lower case.</summary>
    All = 'A',

    /// <summary><c>I</c>: identities to the reference; everything else is <c>.</c>.</summary>
    Identities = 'I',

    /// <summary><c>N</c> (default): non-identities; identities are <c>.</c>.</summary>
    NonIdentities = 'N',

    /// <summary><c>S</c>: similarities (positive matrix score); dissimilar residues are <c>.</c>.</summary>
    Similarities = 'S',

    /// <summary><c>D</c>: dissimilarities; similar residues are <c>.</c>.</summary>
    Dissimilarities = 'D',
}

/// <summary>Output order of the sequences (EMBOSS <c>showalign -order</c>).</summary>
public enum ShowAlignOrder
{
    /// <summary><c>I</c>: input order.</summary>
    Input = 'I',

    /// <summary><c>A</c>: alphabetical order of the names.</summary>
    Alphabetical = 'A',

    /// <summary><c>S</c>: decreasing similarity (matrix score) to the reference.</summary>
    Similarity = 'S',
}

/// <summary>Options of <see cref="MultipleAlignmentDisplay.ShowAlign"/>; the defaults are those of <c>showalign</c> 6.6.0.</summary>
public sealed record ShowAlignOptions
{
    /// <summary>Scoring matrix (<c>-matrix</c>); null = EDNAFULL for a nucleotide set, EBLOSUM62 for a protein set.</summary>
    public SubstitutionMatrix? Matrix { get; init; }

    /// <summary><c>-refseq</c>: a 1-based number or a sequence name; "0" (default) = the consensus is the reference.</summary>
    public string RefSeq { get; init; } = "0";

    /// <summary><c>-bottom</c> (default true): the consensus at the bottom instead of the top.</summary>
    public bool Bottom { get; init; } = true;

    /// <summary><c>-show</c>.</summary>
    public ShowAlignContent Show { get; init; } = ShowAlignContent.NonIdentities;

    /// <summary><c>-order</c>.</summary>
    public ShowAlignOrder Order { get; init; } = ShowAlignOrder.Input;

    /// <summary><c>-similarcase</c> (default true).</summary>
    public bool SimilarCase { get; init; } = true;

    /// <summary><c>-consensus</c> (default true): show the consensus line.</summary>
    public bool Consensus { get; init; } = true;

    /// <summary><c>-uppercase</c>: 1-based inclusive column ranges put in upper case.</summary>
    public IReadOnlyList<(int Start, int End)> Uppercase { get; init; } = [];

    /// <summary><c>-number</c> (default true): the position-number line.</summary>
    public bool Number { get; init; } = true;

    /// <summary><c>-ruler</c> (default true): the tick line.</summary>
    public bool Ruler { get; init; } = true;

    /// <summary><c>-width</c> (default 60, at least 1).</summary>
    public int Width { get; init; } = 60;

    /// <summary><c>-margin</c>: width of the name column; 0 = no names, -1 (default) = the longest name plus one.</summary>
    public int Margin { get; init; } = -1;

    /// <summary><c>-plurality</c> (default 50, in percent of the total weight) of the consensus.</summary>
    public float Plurality { get; init; } = 50f;

    /// <summary><c>-setcase</c>: weight at or below which the consensus is lower case; null = half the number of sequences.</summary>
    public float? SetCase { get; init; }

    /// <summary><c>-identity</c> (default 0, percent) of the consensus.</summary>
    public float Identity { get; init; }

    /// <summary><c>-gaps</c> (default true): a gap character can be the consensus.</summary>
    public bool Gaps { get; init; } = true;
}

/// <summary>
/// EMBOSS <c>showalign</c>: the "pretty" display of a multiple sequence alignment (EMBOSS 6.6.0 <c>emboss/showalign.c</c>, <c>nucleus/embcons.c</c> consensus, <c>embshow.c</c>
/// <c>embShowUpperRange</c>).
/// </summary>
public static class MultipleAlignmentDisplay
{
    /// <summary>
    /// <see cref="ShowAlign"/>: the largest display text it builds, in characters (blocks x lines per block x (margin + width + 1)): 2^27 characters = 256 MiB. A larger request - a huge margin, a width of one over a very long
    /// alignment - is refused with an <see cref="ArgumentException"/> before the text is built; the program streams its output and has no limit.
    /// </summary>
    public const long ShowAlignMaxOutputCharacters = 1L << 27;

    /// <summary>
    /// The text of <c>showalign</c> for an alignment (the HTML formatting and <c>-highlight</c> are not part of this port). The reference sequence (the consensus unless
    /// <see cref="ShowAlignOptions.RefSeq"/> names a row) is the one against which every other row is converted according to <see cref="ShowAlignOptions.Show"/>:
    /// identities / non-identities (case-insensitive equality, <c>.</c> marks the suppressed residues), similarities / dissimilarities (matrix score &gt; 0 / ≤ 0, an unknown or gap code scores 0), and
    /// with <see cref="ShowAlignOptions.SimilarCase"/> the case encodes similarity. The consensus is the EMBOSS <c>embConsCalc</c> consensus (the algorithm of <c>cons</c>, here with the gap option),
    /// of the alignment as read (the case of the rows is kept). The set is typed like the program does (nucleotide when the first row looks like one); rows shorter than the longest are padded with trailing gaps.
    /// </summary>
    /// <param name="names">Sequence names (first word is used; shown left-justified and truncated to the margin).</param>
    /// <param name="rows">The aligned sequences (letters, <c>-</c> <c>.</c> <c>~</c> gaps, <c>?</c> read as X).</param>
    /// <param name="options">The options, null = all defaults.</param>
    /// <returns>The display text, lines ending in <c>\n</c>, a blank line after each block.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">No row, different counts of names and rows, the reference is neither a number in range nor a name, or an invalid character.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Width &lt; 1, margin &lt; -1, a non-finite or out-of-range plurality / identity.</exception>
    public static string ShowAlign(IReadOnlyList<string> names, IReadOnlyList<string> rows, ShowAlignOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(names);
        ArgumentNullException.ThrowIfNull(rows);
        var o = options ?? new ShowAlignOptions();
        if (rows.Count == 0)
            throw new ArgumentException("At least one sequence is required.", nameof(rows));
        if (names.Count != rows.Count)
            throw new ArgumentException("One name per row is required.", nameof(names));
        ArgumentOutOfRangeException.ThrowIfLessThan(o.Width, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.Margin, -1);
        if (!(o.Plurality >= 0f && o.Plurality <= 100f))
            throw new ArgumentOutOfRangeException(nameof(options), "Plurality must be in [0, 100].");
        if (!(o.Identity >= 0f && o.Identity <= 100f))
            throw new ArgumentOutOfRangeException(nameof(options), "Identity must be in [0, 100].");
        foreach (var name in names)
            ArgumentNullException.ThrowIfNull(name);

        // the sequence set as the reader hands it over: padded (ajSeqsetFill), gaps read as '-', case kept
        int nseqs = rows.Count;
        if (rows.Any(r => r is null))
            throw new ArgumentException("Sequences cannot contain null elements.", nameof(rows));
        var read = rows.Select(EmbossSequenceReader.CleanGapped).ToList();          // the file reader keeps the letters and * . ~ ? # + - (A4-3)
        int mlen = read.Max(r => r.Length);
        string[] seqNames = names.Select(n => n.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries).FirstOrDefault() ?? string.Empty).ToArray();
        int nref = ReferenceIndex(o.RefSeq, seqNames);

        int margin = o.Margin;
        if (margin == -1)
        {
            int longest = o.Consensus || nref == -1 ? 9 : 0;
            foreach (var n in seqNames)
                longest = Math.Max(longest, n.Length);
            margin = longest + 1;
        }

        // A4-2: a width beyond the alignment shows one block (the same text), and the display text must stay below a documented size - a margin or a width that would make it larger is refused
        // before anything is built (the program streams its output and has no limit; the 32-bit position arithmetic of Output is done in 64 bit)
        int width = Math.Min(o.Width, Math.Max(1, mlen));
        long blocks = (mlen + (long)width - 1) / width;
        long linesPerBlock = nseqs + (o.Consensus ? 1 : 0) + (o.Number ? 1 : 0) + (o.Ruler ? 1 : 0) + 1;
        long estimate = blocks * linesPerBlock * ((long)margin + width + 1);
        if (estimate > ShowAlignMaxOutputCharacters)
            throw SequenceAligner.ResourceLimit("EMBOSS showalign", "builds the whole display text", $"{blocks} blocks x {linesPerBlock} lines x {(long)margin + width + 1} characters = {estimate} characters",
                ShowAlignMaxOutputCharacters, nameof(MultipleAlignmentDisplay) + "." + nameof(ShowAlignMaxOutputCharacters), "The program streams its output and has no such limit", "use a wider display (a larger width), a smaller margin or fewer rows", nameof(options));

        var padded = read.Select(r => r.PadRight(mlen, '-')).ToList();
        bool nucleotide = EmbossSequenceReader.IsNucleotide(padded[0]);
        var matrix = o.Matrix ?? (nucleotide ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62);
        var seqs = MotifFinder.NormalizedAlignmentChars(padded, ConsensusResidueType.Auto, nucleotide, nameof(rows), keepCase: true);
        var codes = matrix;

        // consensus (embConsCalc): plurality in weight units, identity as a count of sequences
        float total = nseqs;
        float fplural = total * o.Plurality / 100f;
        int ident = nseqs * (int)o.Identity / 100;
        char[] consensus = MotifFinder.EmbossConsensusCore(padded, ConsensusResidueType.Auto, fplural, ident, o.SetCase ?? total / 2f, null,
            o.Gaps, keepCase: true, writeAsNucleotide: false, matrix, nameof(rows)).ToCharArray();

        char[] reference = nref == -1 ? consensus : seqs[nref];
        var order = OrderRows(o.Order, seqNames, seqs, reference, nref, codes);
        foreach (int i in Enumerable.Range(0, nseqs).Where(i => i != nref))
            Convert(o.Show, reference, seqs[i], codes, o.SimilarCase);

        return Output(o, seqNames, seqs, nref, margin, consensus, order, mlen, width);
    }

    // showalign_Getrefseq: an exact name, else a unique name ignoring case, else a number 0..n (0 = the consensus; n > count is rejected)
    private static int ReferenceIndex(string refseq, string[] names)
    {
        ArgumentNullException.ThrowIfNull(refseq);
        int caseMatch = 0, found = 0;
        for (int i = 0; i < names.Length; i++)
        {
            if (string.Equals(names[i], refseq, StringComparison.Ordinal))
                return i;
            if (string.Equals(names[i], refseq, StringComparison.OrdinalIgnoreCase))
            {
                caseMatch = i;
                found++;
            }
        }

        if (found == 1)
            return caseMatch;
        if (int.TryParse(refseq, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int number))
        {
            if (number < 0 || number > names.Length)
                throw new ArgumentException($"Reference sequence number {number} is < 0 or > the number of sequences {names.Length}.", nameof(refseq));
            return number - 1;
        }

        throw new ArgumentException($"Reference sequence '{refseq}' is not a sequence name or a number.", nameof(refseq));
    }

    private static int[] OrderRows(ShowAlignOrder order, string[] names, char[][] seqs, char[] reference, int nref, SubstitutionMatrix codes)
    {
        var others = Enumerable.Range(0, seqs.Length).Where(i => i != nref).ToList();
        switch (order)
        {
            case ShowAlignOrder.Input:
                return others.ToArray();
            case ShowAlignOrder.Alphabetical:
                return others.OrderBy(i => names[i], StringComparer.Ordinal).ToArray();   // stable; qsort of the program is not guaranteed to be for equal names
            case ShowAlignOrder.Similarity:
                {
                    var score = others.ToDictionary(i => i, i =>
                    {
                        int sum = 0;
                        for (int k = 0; k < seqs[i].Length && k < reference.Length; k++)
                            sum += Score(codes, reference[k], seqs[i][k]);
                        return sum;
                    });
                    return others.OrderByDescending(i => score[i]).ToArray();
                }
            default:
                throw new ArgumentOutOfRangeException(nameof(order));
        }
    }

    // the EMBOSS programs read a pair that involves a character outside the matrix (a gap, an unknown letter) as score 0
    private static int Score(SubstitutionMatrix matrix, char a, char b) => matrix.TryGetScore(a, b, out int s) ? s : 0;

    private static char Up(char c) => c is >= 'a' and <= 'z' ? (char)(c - 32) : c;

    private static char Low(char c) => c is >= 'A' and <= 'Z' ? (char)(c + 32) : c;

    private static void Convert(ShowAlignContent show, char[] r, char[] s, SubstitutionMatrix codes, bool similarCase)
    {
        int lenref = r.Length, lenseq = s.Length;
        int i;
        switch (show)
        {
            case ShowAlignContent.All:
                if (!similarCase)
                    return;
                for (i = 0; i < lenref && i < lenseq; i++)
                    if (s[i] != '-')
                        s[i] = Score(codes, r[i], s[i]) <= 0 ? Low(s[i]) : Up(s[i]);
                for (; i < lenseq; i++)
                    if (s[i] != '-')
                        s[i] = Low(s[i]);
                break;
            case ShowAlignContent.Identities:
                for (i = 0; i < lenref && i < lenseq; i++)
                    if (Up(s[i]) != Up(r[i]) && s[i] != '-')
                        s[i] = '.';
                for (; i < lenseq; i++)
                    if (s[i] != '-')
                        s[i] = '.';
                break;
            case ShowAlignContent.NonIdentities:
                for (i = 0; i < lenref && i < lenseq; i++)
                {
                    if (s[i] == '-')
                        continue;
                    if (Up(s[i]) == Up(r[i]))
                        s[i] = '.';
                    else if (similarCase)
                        s[i] = Score(codes, r[i], s[i]) <= 0 ? Up(s[i]) : Low(s[i]);
                }

                if (lenseq > lenref && similarCase)
                    for (i = lenref; i < lenseq; i++)
                        if (s[i] != '-')
                            s[i] = Up(s[i]);
                break;
            case ShowAlignContent.Similarities:
                for (i = 0; i < lenref && i < lenseq; i++)
                {
                    if (s[i] == '-')
                        continue;
                    if (Score(codes, r[i], s[i]) <= 0)
                        s[i] = '.';
                    else if (similarCase)
                        s[i] = Up(s[i]) == Up(r[i]) ? Up(s[i]) : Low(s[i]);
                }

                for (i = lenref; i < lenseq; i++)
                    if (s[i] != '-')
                        s[i] = '.';
                break;
            case ShowAlignContent.Dissimilarities:
                for (i = 0; i < lenref && i < lenseq; i++)
                    if (s[i] != '-' && Score(codes, r[i], s[i]) > 0)
                        s[i] = '.';
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(show));
        }
    }

    private static string Output(ShowAlignOptions o, string[] names, char[][] seqs, int nref, int margin, char[] consensus, int[] order, int mlen, int width)
    {
        var sb = new StringBuilder();
        string refstr = new(nref == -1 ? consensus : SourceRef(seqs, nref));
        string tickline = Ticks(nref, refstr, out int firstpos, out int lastpos);
        string numline = Numbers(nref, refstr, firstpos, lastpos, margin);
        int begin = 0, end = mlen - 1;

        for (int pos = begin; pos <= end; pos += width)
        {
            if (o.Number)
            {
                int imargin = margin;
                int istart = pos + margin;
                while (istart > 0 && CharAt(numline, istart - 1) == ' ')
                {
                    istart--;
                    if (imargin > 0)
                        imargin--;
                }

                long iendWide = (long)pos + width + margin - 1;           // 64 bit: the sum of the three 32-bit quantities overflowed for a huge width or margin
                int iend;
                if (iendWide < (long)end + margin)
                {
                    iend = (int)iendWide;
                    while (iend < end && CharAt(numline, iend + 1) != ' ')
                        iend++;
                }
                else
                {
                    iend = end + margin;
                }

                sb.Append(' ', imargin).Append(Substring(numline, istart, iend)).Append('\n');
            }

            if (o.Ruler)
                sb.Append(' ', margin).Append(Substring(tickline, pos, pos + width - 1)).Append('\n');

            if (nref != -1)
                OutputSeq(sb, names[nref], seqs[nref], pos, end, width, margin, o.Uppercase);
            else if (o.Consensus && !o.Bottom)
                OutputSeq(sb, "Consensus", consensus, pos, end, width, margin, o.Uppercase);

            foreach (int i in order)
                OutputSeq(sb, names[i], seqs[i], pos, end, width, margin, o.Uppercase);

            if (o.Consensus && o.Bottom)
                OutputSeq(sb, "Consensus", consensus, pos, end, width, margin, o.Uppercase);
            sb.Append('\n');
        }

        return sb.ToString();
    }

    // the reference row is never converted, so seqs[nref] still holds the original letters
    private static char[] SourceRef(char[][] seqs, int nref) => seqs[nref];

    private static char CharAt(string s, int index) => index >= 0 && index < s.Length ? s[index] : '\0';

    private static string Substring(string s, int first, int last) => EmbossText.Substring(s, first, last);

    // one row of one block: the slice [pos, min(end, pos + width - 1)] of the row (EmbossText.Substring clipping), taken from the char[] directly - O(width) per block instead of a copy of the whole row
    private static void OutputSeq(StringBuilder sb, string name, char[] seq, int pos, int end, int width, int margin, IReadOnlyList<(int Start, int End)> uppercase)
    {
        long stopWide = Math.Min((long)end, (long)pos + width - 1);
        int stop = (int)stopWide;
        if (margin > 0)
            sb.Append(name.Length >= margin ? name[..margin] : name.PadRight(margin));
        int lineStart = sb.Length;
        if (seq.Length > 0)
        {
            int first = Math.Clamp(pos, 0, seq.Length - 1);
            int last = Math.Clamp(stop, first, seq.Length - 1);
            sb.Append(seq, first, last - first + 1);
        }

        int length = sb.Length - lineStart;
        if (uppercase.Count > 0)
        {
            foreach (var (start, lastUpper) in uppercase)
                for (int j = Math.Max(start - 1, pos); j <= lastUpper - 1 && j - pos < length; j++)
                    sb[lineStart + j - pos] = Up(sb[lineStart + j - pos]);
        }

        sb.Append('\n');
    }

    private static char TickChar(int i)
    {
        if (i % 10 == 0)
            return '|';
        return i % 5 == 0 ? ':' : '-';
    }

    // showalign_OutputTicks
    private static string Ticks(int nref, string refstr, out int firstpos, out int lastpos)
    {
        int iend = refstr.Length;
        var line = new StringBuilder();
        int ifirst = 0, ilast = 0;
        if (nref == -1)
        {
            for (int i = 1; i <= iend; i++)
                line.Append(TickChar(i));
            ifirst = 1;
            ilast = iend;
        }
        else
        {
            int i = 0, ipos = 0;
            foreach (char ch in refstr)
            {
                ipos++;
                if (ch == '-')
                {
                    line.Append('.');
                    continue;
                }

                i++;
                line.Append(TickChar(i));
                if (ifirst == 0)
                    ifirst = line.Length;
                ilast = ipos;
            }
        }

        int j = 0;
        if (ifirst > 1)
            for (int i = ifirst - 2; i > 0; i--)
            {
                j++;
                if (j % 10 == 0) line[i] = 'V';
                else if (j % 5 == 0) line[i] = 'v';
            }

        if (ilast < iend)
        {
            j = 0;
            for (int i = ilast; i < iend; i++)
            {
                j++;
                if (j % 10 == 0) line[i] = 'V';
                else if (j % 5 == 0) line[i] = 'v';
            }
        }

        firstpos = ifirst;
        lastpos = ilast;
        return line.ToString();
    }

    // showalign_OutputNums
    private static string Numbers(int nref, string refstr, int firstpos, int lastpos, int margin)
    {
        var line = new StringBuilder();
        int iend = refstr.Length;
        int iwidth;
        if (nref == -1)
        {
            iwidth = margin + 10;
            for (int i = 10; i < iend; i += 10)
            {
                line.Append(i.ToString(CultureInfo.InvariantCulture).PadLeft(iwidth));
                iwidth = 10;
            }

            iwidth = iend % 10;
        }
        else
        {
            iwidth = margin;
            int j = firstpos;
            for (int i = 1; i < firstpos; i++)
            {
                j--;
                iwidth++;
                if (j % 10 == 0)
                {
                    line.Append(("-" + j.ToString(CultureInfo.InvariantCulture)).PadLeft(iwidth));
                    iwidth = 0;
                }
            }

            j = 0;
            for (int i = firstpos - 1; i < lastpos; i++)
            {
                iwidth++;
                if (refstr[i] != '-')
                {
                    j++;
                    if (j % 10 == 0)
                    {
                        line.Append(j.ToString(CultureInfo.InvariantCulture).PadLeft(iwidth));
                        iwidth = 0;
                    }
                }
            }

            j = 0;
            for (int i = lastpos; i <= iend; i++)
            {
                iwidth++;
                j++;
                if (j % 10 == 0)
                {
                    line.Append(("+" + j.ToString(CultureInfo.InvariantCulture)).PadLeft(iwidth));
                    iwidth = 0;
                }
            }
        }

        if (iwidth > 0)
            line.Append(' ', iwidth);
        return line.ToString();
    }
}
