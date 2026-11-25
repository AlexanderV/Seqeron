// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

// B09 WP61 (F280): trimest, newseq, yank, seqcount, splitter, maskseq, maskambignuc, maskambigprot.
public static partial class EmbossSequencePrograms
{
    // ---------------------------------------------------------------------------------------------------------------- trimest

    /// <summary><c>trimest</c>: the poly-A tail (and the poly-T head: the sequence is then reverse-complemented) of ESTs removed, or only written in lower case (<c>-tolower</c>).</summary>
    internal static EmbossRunResult Trimest(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("trimest");
        var all = EmbossProgramBase.Open("trimest", a, files, ctx, "dna");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        int minlength, mismatches;
        bool reverse, fiveprime, tolower;
        EmbossOutseqOptions outOptions;
        try
        {
            minlength = a.IntLimited("minlength", 4, 1, int.MaxValue, ctx);
            mismatches = a.IntLimited("mismatches", 1, 0, int.MaxValue, ctx);
            reverse = a.Bool("reverse", true);
            tolower = a.Bool("tolower", false);
            fiveprime = a.Bool("fiveprime", true);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("trimest", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        int tail5 = 0;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            string desc = seq.Desc;
            string text = seq.SeqText;
            if (fiveprime)
                tail5 = EstTail(text, 5, minlength, mismatches);
            int tail3 = EstTail(text, 3, minlength, mismatches);
            string str = text;
            if (tail5 > tail3)
            {
                if (tolower)
                    str = LowerRange(str, 0, tail5 - 1);
                else
                    str = KeepRange(str, tail5, text.Length - 1);
                desc += " [poly-T tail removed]";
            }
            else if (tail3 > tail5)
            {
                if (tolower)
                    str = LowerRange(str, text.Length - tail3, text.Length);
                else
                    str = KeepRange(str, 0, text.Length - tail3 - 1);
                desc += " [poly-A tail removed]";
            }

            EmbossProgramBase.AssignSeq(seq, str);
            if (tail5 > tail3 && reverse)
            {
                seq.Rev = true;
                EmbossSeqRange.ReverseDo(seq);                                       // ajSeqReverseForce
                desc += " [reverse complement]";
            }

            seq.Desc = desc;
            stream.Write(seq);
        }

        return Finish(stream);
    }

    // ajStrKeepRange
    private static string KeepRange(string s, long pos1, long pos2)
    {
        int ibegin = EmbossInputSequence.PositionStart(s.Length, 0, pos1);
        int iend = EmbossInputSequence.PositionStart(s.Length, 0, pos2);
        if (iend == s.Length)
            iend--;
        return s.Length > 0 && iend >= ibegin ? s.Substring(ibegin, iend - ibegin + 1) : "";
    }

    // trimest_tolower / ajRangeStrToLower of one range [st, en] (0-based): the region lower case, cut and inserted again
    private static string LowerRange(string s, long st, long en)
    {
        string sub = EmbossProgramBase.SubS(s, st, en).ToLowerInvariant();
        s = CutRange(s, st, en);
        return Insert(s, st, sub);
    }

    // ajStrCutRange
    private static string CutRange(string s, long pos1, long pos2)
    {
        if (s.Length == 0)
            return s;
        int ibegin = EmbossInputSequence.PositionStart(s.Length, 0, pos1);
        int iend = EmbossInputSequence.PositionStart(s.Length, ibegin, pos2) + 1;
        return s.Remove(ibegin, iend - ibegin);
    }

    // ajStrInsertS: the position can be at the start or after the end
    private static string Insert(string s, long pos, string text)
    {
        int ibegin = EmbossInputSequence.PositionStart(s.Length + 1, 0, pos);
        return s.Insert(Math.Min(ibegin, s.Length), text);
    }

    // trimest_get_tail
    private static int EstTail(string s, int direction, int minlength, int mismatches)
    {
        char t;
        int inc, start, end;
        if (direction == 5)
        {
            t = 'T';
            inc = 1;
            start = 0;
            end = s.Length;
        }
        else
        {
            t = 'A';
            inc = -1;
            start = s.Length - 1;
            end = -1;
        }

        int mismatchcount = 0, polycount = 0, length = 0, result = 0;
        for (int i = start; i != end; i += inc, length++)
        {
            char c = char.ToUpperInvariant(s[i]);
            if (c == t)
            {
                polycount++;
                mismatchcount = 0;
            }
            else if (c != 'N')
            {
                polycount = 0;
                mismatchcount++;
            }

            if (polycount >= minlength)
                result = length + 1;
            if (mismatchcount > mismatches)
                break;
        }

        return result;
    }

    // ---------------------------------------------------------------------------------------------------------------- newseq

    /// <summary><c>newseq</c>: a sequence from a name, a description, a type and the residues (no input sequence).</summary>
    internal static EmbossRunResult Newseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("newseq");
        string name, desc, sequence, type;
        EmbossOutseqOptions outOptions;
        try
        {
            name = a.Raw("name") ?? throw a.Die("name");
            desc = a.Raw("description") ?? throw a.Die("description");
            type = a.Menu("type", "N", "N", "P");
            sequence = a.Raw("sequence") ?? throw a.Die("sequence");
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("newseq", ctx, outOptions, type == "P" ? "gapprotein" : "gapnucleotide");
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        var seq = EmbossProgramBase.NewSeq(name, desc, "", type);
        if (type == "N")
            seq.SetNuc();
        EmbossProgramBase.AssignSeq(seq, sequence);
        stream.Write(seq);
        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- yank, seqcount

    // an ACD outfile: the name "stdout" is the standard output (a file opened for append is a file of that name)
    private static EmbossRunResult OutfileResult(EmbossFeatureReadContext ctx, string outfile, string text, bool append = false)
    {
        if (outfile == "stdout" && !append)
            return new EmbossRunResult(text, EmbossRunResult.Join(ctx.Messages), 0, new Dictionary<string, string>());
        var filesOut = new Dictionary<string, string> { [outfile] = text };
        return new EmbossRunResult("", EmbossRunResult.Join(ctx.Messages), 0, filesOut);
    }

    /// <summary><c>yank</c>: the USA of the sequence written to a file.</summary>
    internal static EmbossRunResult Yank(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("yank");
        var all = EmbossProgramBase.Open("yank", a, files, ctx, "gapany", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Next()!;
        string outfile;
        bool newfile;
        try
        {
            newfile = a.Bool("newfile", false);
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "yank", "sequence", seq.Name);
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        // -newfile no: the file is opened for append (the ACD attribute append: @(!$(newfile)))
        string old = !newfile && files.TryGetValue(outfile, out string? existing) ? existing : "";
        return OutfileResult(ctx, outfile, old + UsaOf(seq) + "\n", !newfile);
    }

    // ajSeqGetUsaS / seqMakeUsa: format::file:name with the range [begin:end] or [end:begin:r] of the sequence
    public static string UsaOf(EmbossSeq t)
    {
        string usa;
        if (t.Db.Length > 0)
        {
            usa = $"{t.Db}-id:{t.Entryname}";
        }
        else if (t.Formatstr == "asis")
        {
            usa = "asis::" + (t.Reversed ? EmbossNuc.ReverseComplement(t.SeqText) : t.SeqText);
        }
        else
        {
            usa = t.Entryname.Length > 0 ? $"{t.Formatstr}::{t.Filename}:{t.Entryname}" : $"{t.Formatstr}::{t.Filename}";
        }

        var range = new StringBuilder("[");
        if (t.Reversed)
        {
            if (t.End != 0)
                range.Append((-t.End).ToString(CultureInfo.InvariantCulture));
            range.Append(':');
            if (t.Begin != 0)
                range.Append((-t.Begin).ToString(CultureInfo.InvariantCulture));
            range.Append(":r");
        }
        else
        {
            if (t.Begin != 0)
                range.Append(t.Begin.ToString(CultureInfo.InvariantCulture));
            range.Append(':');
            if (t.End != 0)
                range.Append(t.End.ToString(CultureInfo.InvariantCulture));
        }

        range.Append(']');
        return range.Length > 3 ? usa + range : usa;
    }

    /// <summary><c>seqcount</c>: the number of sequences.</summary>
    internal static EmbossRunResult Seqcount(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("seqcount");
        var all = EmbossProgramBase.Open("seqcount", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        string firstName = all.Peek().Name;
        uint count = 0;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out _))
            count++;
        string outfile;
        try
        {
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "seqcount", "sequence", firstName);
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        return OutfileResult(ctx, outfile, count.ToString(CultureInfo.InvariantCulture) + "\n");
    }

    // ---------------------------------------------------------------------------------------------------------------- splitter

    /// <summary><c>splitter</c>: sequences split into pieces of <c>-size</c> (with <c>-overlap</c>; <c>-addoverlap</c> adds the overlap to the pieces), named <c>name_begin-end</c>.</summary>
    internal static EmbossRunResult Splitter(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("splitter");
        if (a.Raw("feature") is not null && a.Bool("feature", false))
            throw new NotSupportedException("splitter -feature (the features of the pieces) is not ported.");
        var all = EmbossProgramBase.Open("splitter", a, files, ctx, "any");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        int size, overlap;
        bool addover;
        EmbossOutseqOptions outOptions;
        try
        {
            size = a.IntLimited("size", 10000, 1, int.MaxValue, ctx);
            overlap = a.IntLimited("overlap", 0, 0, size - 1, ctx);
            addover = a.Bool("addoverlap", false);
            a.Bool("feature", false);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("splitter", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            EmbossSeqRange.Trim(seq);
            int len = seq.Seq.Length, pos = 0;
            string text = seq.SeqText;
            if (!addover)
            {
                while (pos + size <= len - 1)
                {
                    Chunk(stream, seq, text, pos, pos + size - 1);
                    pos += size - overlap;
                }
            }
            else
            {
                while (pos + size + overlap < len - 1)
                {
                    Chunk(stream, seq, text, pos, pos + size + overlap - 1);
                    pos += size;
                }
            }

            Chunk(stream, seq, text, pos, len - 1);
        }

        return Finish(stream);
    }

    private static void Chunk(EmbossSeqoutStream stream, EmbossSeq seq, string text, int start, int end)
    {
        int begin = EmbossSeqRange.GetBegin(seq);
        string name = seq.Name + "_" + ((uint)(begin + start)).ToString(CultureInfo.InvariantCulture) + "-" + ((uint)(begin + end)).ToString(CultureInfo.InvariantCulture);
        // ajStrAssignSubC(&str, seq, start, end): negative positions count from the end
        long ibegin = start < 0 ? text.Length + start : start, iend = end < 0 ? text.Length + end : end;
        string str = iend >= ibegin && ibegin >= 0 && iend < text.Length ? text.Substring((int)ibegin, (int)(iend - ibegin + 1)) : "";
        var sub = EmbossProgramBase.NewSeq(name, seq.Desc, str, "");
        EmbossProgramBase.SeqType(sub);
        stream.Write(sub);
    }

    // ---------------------------------------------------------------------------------------------------------------- maskseq

    /// <summary><c>maskseq</c>: regions of a sequence replaced by a mask character (default N, X for protein) or written in lower case.</summary>
    internal static EmbossRunResult Maskseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("maskseq");
        var all = EmbossProgramBase.Open("maskseq", a, files, ctx, "any", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Next()!;
        IReadOnlyList<EmbossRangeElement> regions;
        bool tolower;
        string maskchar;
        EmbossOutseqOptions outOptions;
        try
        {
            regions = Ranges(a, "regions", ctx);
            tolower = a.Bool("tolower", false);
            maskchar = a.Raw("maskchar") ?? (EmbossSeqType.IsProt(seq.Type, seq.SeqText) ? "X" : "N");
            if (a.Raw("maskchar") is not null && maskchar.Length < 1)
                throw a.Die("maskchar", "Error: Too short (maskchar) - minimum length is 1 characters");
            if (maskchar.Length > 1)
                throw a.Die("maskchar", "Error: Too long (maskchar) - maximum length is 1 characters");
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("maskseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        int beg = EmbossSeqRange.GetBegin(seq) - 1, end = EmbossSeqRange.GetEnd(seq) - 1;
        // ajStrAssignSubS(&str, seq, beg, end)
        string text = seq.SeqText;
        int ib = EmbossInputSequence.PositionStart(text.Length, 0, beg);
        int ie = EmbossInputSequence.PositionStart(text.Length, ib, end);
        if (ie == text.Length)
            ie--;
        string str = text.Length > 0 && ie >= ib ? text.Substring(ib, ie - ib + 1) : "";
        // ajRangeSetOffset(regions, beg + 1): the positions as unsigned numbers
        var offsets = new List<(uint St, uint En)>();
        foreach (var r in regions)
            offsets.Add(((uint)r.Start - (uint)beg, (uint)r.End - (uint)beg));
        bool lower = tolower || maskchar.Length == 0 || maskchar == " ";
        foreach (var (st0, en0) in offsets)
        {
            if (st0 < 1 || en0 < 1 || st0 > en0)
                throw new InvalidOperationException("maskseq: a region outside the range of the sequence (the reference program then uses uninitialised values).");
            long st = st0 - 1, en = en0 - 1;
            if (lower)
            {
                str = LowerRange(str, st, en);
            }
            else
            {
                str = CutRange(str, st, en);
                for (long j = st; j <= en; ++j)
                    str = Insert(str, st, maskchar);
            }
        }

        EmbossProgramBase.AssignSeq(seq, str);
        stream.Write(seq);
        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- maskambignuc, maskambigprot

    /// <summary><c>maskambignuc</c> / <c>maskambigprot</c>: the ambiguity codes of the sequences replaced by N / X (done by the type of the sequence input).</summary>
    internal static EmbossRunResult Maskambig(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files, string program, string type)
    {
        var ctx = EmbossProgramBase.Context(program);
        if (a.Raw("feature") is not null && a.Bool("feature", false))
            throw new NotSupportedException(program + " -feature (the features of the sequences) is not ported.");
        var all = EmbossProgramBase.Open(program, a, files, ctx, type);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        EmbossOutseqOptions outOptions;
        try
        {
            a.Bool("feature", false);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream(program, ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
            stream.Write(seq);
        return Finish(stream);
    }

    internal static EmbossRunResult Maskambignuc(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files) => Maskambig(a, files, "maskambignuc", "gapnucleotidesimple");

    internal static EmbossRunResult Maskambigprot(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files) => Maskambig(a, files, "maskambigprot", "gapproteinsimple");

    // the default name of an output file (acdOutFilename of the saved name acdInFName): the name of the first sequence in lower case (acdInFileSave of the object name, trimmed by ajFilenameTrimAll), for an input that is no sequence
    // the base name of the first input file (up to the last dot); the extension is the program
    internal static string DefaultOutfile(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files, string program, string inputParameter, string? firstSequenceName = null, string? extension = null)
    {
        string input = a.Raw(inputParameter) ?? "";
        if (firstSequenceName is not null)
            return EmbossOutPath.DefaultBaseName(firstSequenceName, input) + "." + (extension ?? program);
        int colon = input.IndexOf("::", StringComparison.Ordinal);
        if (colon >= 0)
            input = input[(colon + 2)..];
        int bracket = input.IndexOf('[', StringComparison.Ordinal);
        if (bracket > 0)
            input = input[..bracket];
        string slashless = input[(input.LastIndexOf('/') + 1)..];
        int dot2 = slashless.LastIndexOf('.');
        string baseName = dot2 > 0 ? slashless[..dot2] : slashless;
        return baseName + "." + program;
    }

    // ---------------------------------------------------------------------------------------------------------------- text programs: trimspace, nospace, notab, noreturn

    // ajReadline: the lines of a text with the newline; a last line without one gets it
    public static List<string> ReadLines(string text)
    {
        var lines = new List<string>();
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            if (nl < 0)
            {
                lines.Add(text[pos..] + "\n");
                break;
            }

            lines.Add(text.Substring(pos, nl - pos + 1));
            pos = nl + 1;
        }

        return lines;
    }

    // acdSetInfile of the parameter infile (the file must exist) and acdSetOutfile (the default name <input base>.<program>)
    private static EmbossRunResult? TextOpen(EmbossFeatureReadContext ctx, EmbossAcdArgs a, IReadOnlyDictionary<string, string> files, string program, out string text, out string outfile)
    {
        text = "";
        outfile = "";
        string infile = a.Raw("infile") ?? "";
        if (!files.TryGetValue(infile, out string? t))
        {
            ctx.Error($"Unable to open file '{infile}' for input");
            return EmbossProgramBase.Died(ctx, new EmbossAcdDieException([$"Died: {program} terminated: Bad value for '-infile' with -auto defined"]));
        }

        text = t;
        return null;
    }

    private static EmbossRunResult TextProgram(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files, string program, Func<EmbossAcdArgs, EmbossFeatureReadContext, Func<string, string>?> options)
    {
        var ctx = EmbossProgramBase.Context(program);
        var died = TextOpen(ctx, a, files, program, out string text, out _);
        if (died is not null)
            return died;
        Func<string, string>? line;
        string outfile;
        try
        {
            line = options(a, ctx);
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, program, "infile");
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var sb = new StringBuilder();
        foreach (string l in ReadLines(text))
            sb.Append(line!(l));
        return OutfileResult(ctx, outfile, sb.ToString());
    }

    // ajStrRemoveWhiteExcess
    internal static string RemoveWhiteExcess(string s)
    {
        if (s.Length == 0)
            return s;
        var sb = new StringBuilder(s.Length);
        foreach (char c in s)
            sb.Append(c is '\t' or '\n' or '\r' ? ' ' : c);
        string t = sb.ToString().Trim(' ');
        var o = new StringBuilder(t.Length);
        char prev = '\0';
        foreach (char c in t)
        {
            if (c == ' ' && prev == ' ')
                continue;
            o.Append(c);
            prev = c;
        }

        return o.ToString();
    }

    /// <summary><c>trimspace</c>: leading and trailing white space removed and every run of blanks made one.</summary>
    internal static EmbossRunResult Trimspace(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
        => TextProgram(a, files, "trimspace", (_, _) => l => RemoveWhiteExcess(l) + "\n");

    /// <summary><c>nospace</c>: all white space (<c>-menu all</c>), the trailing white space (<c>end</c>) or the excess white space (<c>excess</c>) removed.</summary>
    internal static EmbossRunResult Nospace(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
        => TextProgram(a, files, "nospace", (args, _) =>
        {
            string menu = args.Menu("menu", "all", "all", "end", "excess");
            return menu switch
            {
                "all" => l => new string([.. l.Select(c => c is '\t' or '\n' or '\r' ? ' ' : c).Where(c => c != ' ')]) + "\n",
                "end" => l => l.TrimEnd('\t', ' ', '\n', '\r') + "\n",
                _ => l => RemoveWhiteExcess(l) + "\n",
            };
        });

    /// <summary><c>notab</c>: tabs replaced by blanks up to the next multiple of eight columns.</summary>
    internal static EmbossRunResult Notab(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
        => TextProgram(a, files, "notab", (_, _) => l =>
        {
            var outline = new StringBuilder();
            int ilen = 0;
            int tab;
            while ((tab = l.IndexOf('\t', StringComparison.Ordinal)) >= 0)
            {
                outline.Append(l, 0, tab);
                ilen += tab;
                int irest = 8 - (ilen - 8 * (ilen / 8));
                outline.Append(' ', irest);
                ilen += irest;
                l = l[(tab + 1)..];
            }

            return outline.Append(l).ToString();
        });

    /// <summary><c>noreturn</c>: the line ends of a file made those of <c>-system</c> unix, pc or mac (a CR inside a line is a line end).</summary>
    internal static EmbossRunResult Noreturn(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
        => TextProgram(a, files, "noreturn", (args, _) =>
        {
            string system = args.Menu("system", "unix", "unix", "pc", "mac");
            return l =>
            {
                // ajReadlineTrim: the newline and then a carriage return at the end are cut
                if (l.EndsWith('\n'))
                    l = l[..^1];
                if (l.EndsWith('\r'))
                    l = l[..^1];
                if (system != "mac")
                    l = l.Replace('\r', '\n');
                return system switch { "pc" => l + "\r\n", "mac" => l + "\r", _ => l + "\n" };
            };
        });

    // what a program wrote to its outputs, in the order it wrote it: (output number, text)
    public sealed class OutLog
    {
        public readonly List<(int Target, string Text)> Segments = [];

        public void Add(int target, string text) => Segments.Add((target, text));
    }

    public static void AssignUpper(EmbossSeq seq, string upper)
    {
        seq.Seq.Clear().Append(upper);
    }
}
