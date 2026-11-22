#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The sequence manipulation programs of EMBOSS 6.6.0 that read the sequences of a file and write sequences (B09 WP61, F280): <c>revseq</c>, <c>transeq</c>, <c>extractseq</c>, <c>cutseq</c>, <c>degapseq</c>,
/// <c>trimseq</c>, <c>skipseq</c>, <c>nthseq</c>, <c>notseq</c>, <c>pasteseq</c>. Every program is the statement-by-statement port of its <c>main</c> over the objects of this port (<see cref="EmbossSeqall"/>,
/// <see cref="EmbossSeq"/>, <see cref="EmbossSeqoutStream"/>, <see cref="EmbossTrn"/>): the arguments are resolved by <see cref="EmbossAcdArgs"/> and the output is what the program leaves.
/// </summary>
internal static partial class EmbossSequencePrograms
{
    private static EmbossRunResult Finish(EmbossSeqoutStream stream) => EmbossRunResult.From(stream.Close());

    // ---------------------------------------------------------------------------------------------------------------- revseq

    /// <summary><c>revseq</c>: the reverse complement of the sequences (<c>-noreverse</c> / <c>-nocomplement</c> for the parts; <c>-tag</c> puts "Reversed:" in front of the description when both are done).</summary>
    internal static EmbossRunResult Revseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("revseq");
        var all = EmbossProgramBase.Open("revseq", a, files, ctx, "gapdna");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        bool reverse, complement, tag;
        EmbossOutseqOptions outOptions;
        try
        {
            reverse = a.Bool("reverse", true);
            complement = a.Bool("complement", true);
            tag = a.Bool("tag", true);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("revseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            if (EmbossSeqType.IsProt(seq.Type, seq.SeqText))
                return EmbossProgramBase.Fatal(ctx, "Input sequence must be nucleic");
            if (reverse && complement)
            {
                seq.Rev = true;
                EmbossSeqRange.ReverseDo(seq);                       // ajSeqReverseForce: reverses and complements
                if (tag)
                    seq.Desc = seq.Desc.Length > 0 ? "Reversed: " + seq.Desc : "Reversed:";   // ajSeqTagRev
            }
            else if (reverse)
            {
                // ajSeqReverseOnly: the residues are reversed and the range is mirrored
                int ibegin = seq.Begin, iend = seq.End;
                if (ibegin != 0)
                    seq.End = -ibegin;
                if (iend != 0)
                    seq.Begin = -iend;
                var chars = seq.SeqText.ToCharArray();
                Array.Reverse(chars);
                seq.Seq.Clear().Append(chars);
            }
            else if (complement)
            {
                var sb = new StringBuilder(seq.Seq.Length);
                foreach (char c in seq.SeqText)
                    sb.Append(EmbossNuc.Comp(c));
                seq.Seq.Clear().Append(sb);                          // ajSeqComplement
            }

            stream.Write(seq);
        }

        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- transeq

    private static readonly string[] TableCodes = ["0", "1", "2", "3", "4", "5", "6", "9", "10", "11", "12", "13", "14", "15", "16", "21", "22", "23"];

    /// <summary>The ranges of an ACD <c>range</c> item (<c>ajRangeNewString</c>): the warning and the error of the library when the text is bad.</summary>
    internal static IReadOnlyList<EmbossRangeElement> Ranges(EmbossAcdArgs a, string name, EmbossFeatureReadContext ctx)
    {
        string? raw = a.Raw(name);
        if (raw is null)
            return [];
        try
        {
            return EmbossRangeSpec.Parse(raw);
        }
        catch (ArgumentException e)
        {
            throw a.Die(name, "Warning: " + System.Text.RegularExpressions.Regex.Replace(e.Message, @" \(Parameter '[^']*'\)$", ""), $"Error: Bad range specification '{raw}'");
        }
    }

    /// <summary><c>transeq</c>: the translation of the frames of the sequences (<c>-frame</c> 1 2 3 F -1 -2 -3 R 6, <c>-table</c>, <c>-regions</c>, <c>-trim</c>, <c>-clean</c>, <c>-alternative</c>, <c>-methionine</c>).</summary>
    internal static EmbossRunResult Transeq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("transeq");
        var all = EmbossProgramBase.Open("transeq", a, files, ctx, "nucleotide");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var frames = new bool[6];
        int table;
        IReadOnlyList<EmbossRangeElement> regions;
        bool trim, clean, alternative, methionine;
        EmbossOutseqOptions outOptions;
        try
        {
            foreach (string f in FrameList(a))
            {
                switch (f)
                {
                    case "1": frames[0] = true; break;
                    case "2": frames[1] = true; break;
                    case "3": frames[2] = true; break;
                    case "-1": frames[3] = true; break;
                    case "-2": frames[4] = true; break;
                    case "-3": frames[5] = true; break;
                    case "F" or "f": frames[0] = frames[1] = frames[2] = true; break;
                    case "R" or "r": frames[3] = frames[4] = frames[5] = true; break;
                    case "6": Array.Fill(frames, true); break;
                }
            }

            table = int.Parse(a.Menu("table", "0", TableCodes), CultureInfo.InvariantCulture);
            regions = Ranges(a, "regions", ctx);
            trim = a.Bool("trim", false);
            clean = a.Bool("clean", false);
            alternative = a.Bool("alternative", false);
            methionine = a.Bool("methionine", true);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var trn = EmbossTrn.New(table);
        var stream = new EmbossSeqoutStream("transeq", ctx, outOptions, "stopprotein");
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            EmbossSeqRange.Trim(seq);
            bool whole = regions.Count == 0;                                  // ajRangeIsWhole
            if (!whole)
                RangeExtract(seq, regions);                                   // ajRangeSeqExtract
            string residues = seq.SeqText;
            for (int i = 0; i < 6; i++)
            {
                if (!frames[i])
                    continue;
                int frame = i < 3 ? i + 1 : alternative ? -i - 1 : 2 - i;
                string pep = trn.SeqOrig(residues, frame, methionine && !whole);
                if (trim)
                {
                    int end = pep.Length;
                    while (end > 0 && pep[end - 1] is 'X' or '*')
                        end--;
                    pep = pep[..end];
                }

                if (clean)
                    pep = pep.Replace('*', 'X');
                int number = EmbossTrn.PepFrameNumber(frame);
                var o = EmbossProgramBase.NewSeq(number != 0 ? seq.Name + "_" + number.ToString(CultureInfo.InvariantCulture) : seq.Name, seq.Desc, pep, "P");
                stream.Write(o);
            }
        }

        return Finish(stream);
    }

    // the values of the multiple list -frame (acdSetListI: the delimiter is a comma, a value is one of the short codes)
    private static IEnumerable<string> FrameList(EmbossAcdArgs a)
    {
        string? raw = a.Raw("frame");
        if (raw is null)
            return ["1"];
        string[] codes = ["1", "2", "3", "F", "-1", "-2", "-3", "R", "6"];
        var result = new List<string>();
        foreach (string part in raw.Split(','))
        {
            string t = part.Trim();
            string? code = codes.FirstOrDefault(c => string.Equals(c, t, StringComparison.OrdinalIgnoreCase));
            if (code is null)
            {
                string first = raw.Contains(',', StringComparison.Ordinal) ? raw : t;
                throw a.Die("frame", $"Error: '{t}' is not a valid menu option", "Accepted short codes are: " + string.Join(",", codes), $"Error: Bad menu option '{first}'");
            }

            result.Add(code);
        }

        return result;
    }

    // ajRangeSeqExtract: the residues of the ranges in the order of the elements (ajStrAppendSubS clips), the sequence is assigned (range and offsets reset)
    internal static void RangeExtract(EmbossSeq seq, IReadOnlyList<EmbossRangeElement> regions)
    {
        var sb = new StringBuilder();
        string text = seq.SeqText;
        foreach (var r in regions)
        {
            if (r.Start < 1 || r.End < 1 || r.Start > r.End)
                continue;                                                      // ajRangeElementGetValues
            sb.Append(EmbossProgramBase.SubS(text, r.Start - 1, r.End - 1));
        }

        EmbossProgramBase.AssignSeq(seq, sb.ToString());
    }

    // ---------------------------------------------------------------------------------------------------------------- extractseq

    /// <summary><c>extractseq</c>: regions of a sequence, joined or one sequence per region (<c>-separate</c>, named <c>name_start_end</c>).</summary>
    internal static EmbossRunResult Extractseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("extractseq");
        var all = EmbossProgramBase.Open("extractseq", a, files, ctx, "any", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Next()!;
        IReadOnlyList<EmbossRangeElement> regions;
        bool separate;
        EmbossOutseqOptions outOptions;
        try
        {
            if (a.Raw("regions") is null)
                regions = [new EmbossRangeElement(EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), "")];     // default: $(sequence.begin)-$(sequence.end)
            else
                regions = Ranges(a, "regions", ctx);
            separate = a.Bool("separate", false);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("extractseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        string text = seq.SeqText;
        if (separate)
        {
            foreach (var r in regions)
            {
                if (r.Start < 1 || r.End < 1 || r.Start > r.End)
                    continue;
                string str = EmbossProgramBase.SubS(text, r.Start - 1, r.End - 1);
                bool nuc = EmbossSeqType.IsNuc(seq.Type, seq.SeqText);
                var o = EmbossProgramBase.NewSeq(seq.Name + "_" + r.Start.ToString(CultureInfo.InvariantCulture) + "_" + r.End.ToString(CultureInfo.InvariantCulture), seq.Desc, str, nuc ? "N" : "P");
                stream.Write(o);
            }
        }
        else
        {
            RangeExtract(seq, regions);
            stream.Write(seq);
        }

        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- cutseq

    /// <summary><c>cutseq</c>: the residues from <c>-from</c> to <c>-to</c> removed.</summary>
    internal static EmbossRunResult Cutseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("cutseq");
        var all = EmbossProgramBase.Open("cutseq", a, files, ctx, "gapany", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Next()!;
        int from, to;
        EmbossOutseqOptions outOptions;
        int b = EmbossSeqRange.GetBegin(seq), e = EmbossSeqRange.GetEnd(seq);
        try
        {
            from = a.IntLimited("from", b, b, e, ctx);
            to = a.IntLimited("to", e, from, e, ctx);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("cutseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        int beg = b - 1, end = e - 1;
        string str = EmbossProgramBase.SubS(seq.SeqText, beg, end);
        // ajStrCutRange(&str, from-beg, to-beg)
        {
            int ibegin = EmbossInputSequence.PositionStart(str.Length, 0, from - 1 - beg);
            int iend = EmbossInputSequence.PositionStart(str.Length, ibegin, to - 1 - beg) + 1;
            str = str.Length == 0 ? "" : str.Remove(ibegin, iend - ibegin);
        }

        EmbossProgramBase.AssignSeq(seq, str);
        stream.Write(seq);
        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- degapseq

    /// <summary><c>degapseq</c>: the gap characters removed (everything that is no letter and no star).</summary>
    internal static EmbossRunResult Degapseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("degapseq");
        var all = EmbossProgramBase.Open("degapseq", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        EmbossOutseqOptions outOptions;
        try
        {
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("degapseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            EmbossProgramBase.AssignSeq(seq, EmbossSeqType.RemoveGap(seq.SeqText));
            stream.Write(seq);
        }

        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- trimseq

    /// <summary><c>trimseq</c>: the ends of the sequences trimmed of ambiguity characters and gaps in a window.</summary>
    internal static EmbossRunResult Trimseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("trimseq");
        var all = EmbossProgramBase.Open("trimseq", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        int window;
        float percent;
        bool strict, star, left, right;
        EmbossOutseqOptions outOptions;
        try
        {
            window = a.IntLimited("window", 1, 1, all.Peek().Seq.Length, ctx);
            percent = a.Float("percent", 100.0f);
            strict = a.Bool("strict", false);
            star = a.Bool("star", false);
            left = a.Bool("left", true);
            right = a.Bool("right", true);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("trimseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            bool isnuc = EmbossSeqType.IsNuc(seq.Type, seq.SeqText);
            string text = seq.SeqText;
            int start = left ? TrimEnd(text, 1, isnuc, window, percent, strict, star) + 1 : 0;
            int end = right ? TrimEnd(text, 0, isnuc, window, percent, strict, star) - 1 : text.Length - 1;
            // ajStrKeepRange(&str, start, end)
            int ibegin = EmbossInputSequence.PositionStart(text.Length, 0, start);
            int iend = EmbossInputSequence.PositionStart(text.Length, 0, end);
            if (iend == text.Length)
                iend--;
            string kept = text.Length > 0 && iend >= ibegin ? text.Substring(ibegin, iend - ibegin + 1) : "";
            EmbossProgramBase.AssignSeq(seq, kept);
            stream.Write(seq);
        }

        return Finish(stream);
    }

    // trimseq_trim
    private static int TrimEnd(string s, int sense, bool isnuc, int window, float percent, bool strict, bool star)
    {
        var gang = new bool[256];
        Array.Fill(gang, true);
        void Parole(string good)
        {
            foreach (char c in good)
            {
                gang[char.ToLowerInvariant(c)] = false;
                gang[char.ToUpperInvariant(c)] = false;
            }
        }

        if (isnuc)
        {
            Parole("acgtu.-~ ");
            if (!strict)
                Parole("mrwsykvhdb");
        }
        else
        {
            Parole("arndcqeghilkmfpstwyv.-~ ");
            if (!strict)
                Parole("bz");
            if (star)
                Parole("*");
        }

        int len = s.Length, aa, z, inc, leroy, suspect;
        if (sense != 0 ? window > len : window >= len)
            throw new InvalidOperationException("trimseq: the window is not shorter than the sequence; the reference program then reads the memory around the sequence (the result is not defined).");
        if (sense != 0)
        {
            aa = 0;
            z = len - window;
            inc = 1;
            leroy = -1;
            suspect = -1;
        }
        else
        {
            aa = len - 1;
            z = window;
            inc = -1;
            leroy = len;
            suspect = len;
        }

        bool Bad(int i) => s[i] < 256 && gang[s[i]];
        bool Gap(int i) => s[i] is '.' or '-' or '~' or ' ';
        for (; aa != z; aa += inc)
        {
            if (Bad(aa) || Gap(aa))
                leroy = aa;
            else
                break;
        }

        for (; aa != z; aa += inc)
        {
            int count = 0;
            for (int look = 0; look < window && look > -window; look += inc)
            {
                if (Bad(aa + look))
                {
                    count++;
                    suspect = aa + look;
                }
            }

            float pc = 100.0f * count / window;
            if (pc < percent)
                break;
            if (sense != 0)
            {
                if (suspect > leroy)
                    leroy = suspect;
            }
            else if (suspect < leroy)
            {
                leroy = suspect;
            }
        }

        for (aa = leroy + inc; aa != z; aa += inc)
        {
            if (Gap(aa))
                leroy = aa;
            else
                break;
        }

        return leroy;
    }

    // ---------------------------------------------------------------------------------------------------------------- skipseq, nthseq

    /// <summary><c>skipseq</c>: the first <c>-skip</c> sequences are left out.</summary>
    internal static EmbossRunResult Skipseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("skipseq");
        var all = EmbossProgramBase.Open("skipseq", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        int skip;
        EmbossOutseqOptions outOptions;
        try
        {
            skip = a.Int("skip", 0);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("skipseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        int skipcnt = 0;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            skipcnt++;
            if (skipcnt > skip)
                stream.Write(seq);
        }

        return Finish(stream);
    }

    /// <summary><c>nthseq</c>: the <c>-number</c>th sequence.</summary>
    internal static EmbossRunResult Nthseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("nthseq");
        var all = EmbossProgramBase.Open("nthseq", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        int n;
        EmbossOutseqOptions outOptions;
        try
        {
            n = a.IntLimited("number", 1, 1, int.MaxValue, ctx);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("nthseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        int count = 0;
        EmbossSeq? seq = null;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out seq))
        {
            if (++count == n)
                break;
        }

        if (count != n)
            ctx.Error($"No such sequence - only {count} sequences were input.");
        else
            stream.Write(seq!);
        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- notseq

    /// <summary><c>notseq</c>: the sequences whose name or accession matches <c>-exclude</c> (wildcards; <c>@file</c> = a file of names) go to <c>-junkoutseq</c>, the others to <c>-outseq</c>.</summary>
    internal static EmbossRunResult Notseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("notseq");
        var all = EmbossProgramBase.Open("notseq", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        string exclude;
        EmbossOutseqOptions outOptions, junkOptions;
        try
        {
            exclude = a.Raw("exclude") ?? "";
            outOptions = a.Seqout("outseq");
            junkOptions = a.Seqout("junkoutseq") with { OutSeq = a.Raw("junkoutseq") ?? "/dev/null" };
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("notseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        var junk = new EmbossSeqoutStream("notseq", ctx, junkOptions, null);
        if (junk.OpenFailed)
            return EmbossRunResult.From(junk.DiedOutput("junkoutseq"));
        string pattern;
        if (!exclude.StartsWith('@'))
        {
            pattern = exclude;
        }
        else
        {
            string name = exclude.TrimStart('@');                                    // ajStrTrimC(&filename, "@")
            if (!files.TryGetValue(name, out string? text))
                return EmbossProgramBase.Fatal(ctx, $"Cannot open the file of sequence names: '{name}'");
            var sb = new StringBuilder();
            foreach (string raw in text.Split('\n'))
            {
                string line = raw.TrimEnd('\r');
                if (line.Length == 0 || line[0] is '#' or '!')
                    continue;
                sb.Append(line).Append(',');
            }

            pattern = sb.ToString();
        }

        string[] keys = pattern.Split([' ', '\t', '\n', '\r', ',', ';'], StringSplitOptions.RemoveEmptyEntries);
        bool Matches(string s) => keys.Any(k => EmbossWild.MatchCase(s, k));
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            if (Matches(seq.Name) || Matches(seq.Acc))
                junk.Write(seq);
            else
                stream.Write(seq);
        }

        var main = stream.Close();
        var junked = junk.Close();
        var result = EmbossRunResult.From(main with { Messages = junked.Messages });
        var filesOut = new Dictionary<string, string>(result.Files);
        if (junkOptions.OutSeq != "/dev/null" && junkOptions.OutSeq != "stdout")
            filesOut[junked.SequenceFileName] = junked.SequenceText;
        string stdout = result.Stdout + (junkOptions.OutSeq == "stdout" ? junked.SequenceText : "");
        return new EmbossRunResult(stdout, result.Stderr, result.ExitCode, filesOut);
    }

    // ---------------------------------------------------------------------------------------------------------------- pasteseq

    /// <summary><c>pasteseq</c>: <c>bsequence</c> inserted into <c>asequence</c> after <c>-pos</c>.</summary>
    internal static EmbossRunResult Pasteseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("pasteseq");
        var allA = EmbossProgramBase.Open("pasteseq", a, files, ctx, "gapany", "asequence", true);
        if (allA is null)
            return EmbossProgramBase.Died(ctx, "asequence");
        var seq = allA.Next()!;
        bool aProtein = EmbossSeqType.IsProt(seq.Type, seq.SeqText);
        var allB = EmbossProgramBase.Open("pasteseq", a, files, ctx, aProtein ? "gapstopprotein" : "gapnucleotide", "bsequence", true);
        if (allB is null)
            return EmbossProgramBase.Died(ctx, "bsequence");
        var insert = allB.Next()!;
        int pos;
        EmbossOutseqOptions outOptions;
        try
        {
            pos = a.IntLimited("pos", EmbossSeqRange.GetEnd(seq), 0, EmbossSeqRange.GetEnd(seq), ctx);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException ex)
        {
            return EmbossProgramBase.Died(ctx, ex);
        }

        var stream = new EmbossSeqoutStream("pasteseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        EmbossSeqRange.Trim(seq);
        EmbossSeqRange.Trim(insert);
        string str = seq.SeqText;
        int ibegin = EmbossInputSequence.PositionStart(str.Length + 1, 0, pos);              // ajStrInsertS: can be at the start or after the end
        string result = str.Insert(Math.Min(ibegin, str.Length), insert.SeqText);
        var o = EmbossProgramBase.NewSeq(seq.Name, seq.Desc, result, "");
        EmbossProgramBase.SeqType(o);
        stream.Write(o);
        return Finish(stream);
    }
}
