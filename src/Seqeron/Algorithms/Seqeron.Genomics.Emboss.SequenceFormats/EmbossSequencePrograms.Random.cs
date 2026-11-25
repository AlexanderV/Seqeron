// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The random sequence programs of EMBOSS 6.6.0 (B09 WP62, F281): <c>shuffleseq</c>, <c>msbar</c>, <c>makenucseq</c>, <c>makeprotseq</c>. Every program is the statement-by-statement port of its <c>main</c> over
/// <see cref="EmbossRandom"/> (<c>ajRandomSeed</c> / <c>ajRandomDouble</c> of <c>ajmath.c</c>, seeded from the clock of the run: <see cref="EmbossRunEnvironment"/>).
/// </summary>
public static partial class EmbossSequencePrograms
{
    private static EmbossRandom NewRandom() => new(EmbossRunEnvironment.Current.Seed());

    // (ajint) of a double as the compiled code does it: a value that does not fit is INT_MIN (cvttsd2si)
    private static int CInt(double d) => d is >= -2147483648.0 and < 2147483648.0 ? (int)d : int.MinValue;

    // ---------------------------------------------------------------------------------------------------------------- shuffleseq

    /// <summary><c>shuffleseq</c>: the residues of the sequences in random order (<c>-shuffle</c> times, every shuffle is a sequence of the output).</summary>
    internal static EmbossRunResult Shuffleseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("shuffleseq");
        var all = EmbossProgramBase.Open("shuffleseq", a, files, ctx, "any");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        int shuffles;
        EmbossOutseqOptions outOptions;
        try
        {
            shuffles = a.Int("shuffle", 1);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var random = NewRandom();
        random.Seed();
        var stream = new EmbossSeqoutStream("shuffleseq", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            for (int n = 0; n < shuffles; ++n)
            {
                string shuffled = random.Shuffle(seq.SeqText);
                EmbossProgramBase.AssignSeq(seq, shuffled);
                stream.Write(seq);
            }
        }

        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- msbar

    private const string MsbarValues = "0:None; 1:Any of the following; 2:Insertions; 3:Deletions; 4:Changes; 5:Duplications; 6:Moves";

    /// <summary><c>msbar</c>: point, codon and block mutations of the sequences (insertions, deletions, changes, duplications, moves), not equal to any of the sequences of <c>-othersequence</c> (10 attempts).</summary>
    internal static EmbossRunResult Msbar(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("msbar");
        var all = EmbossProgramBase.Open("msbar", a, files, ctx, "any");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        EmbossSeqall? other;
        try
        {
            other = EmbossSeqall.Open("msbar", a.Raw("othersequence") ?? "asis:N", files, a.Sequence("othersequence"), false, ctx, "any", false);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (other is null)
            return EmbossProgramBase.Died(ctx, "othersequence");
        int count, min, max;
        bool inframe;
        IReadOnlyList<string> pointlist, blocklist, codonlist;
        EmbossOutseqOptions outOptions;
        try
        {
            count = a.IntLimited("count", 1, 0, int.MaxValue, ctx);
            pointlist = a.List("point", "0", MsbarValues, 1, 4);
            blocklist = a.List("block", "0", MsbarValues, 1, 4);
            codonlist = a.List("codon", "0", MsbarValues, 1, 4);
            inframe = a.Bool("inframe", false);
            min = a.IntLimited("minimum", 1, 0, int.MaxValue, ctx);
            max = a.IntLimited("maximum", 10, min, int.MaxValue, ctx);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        const int attempts = 10;
        var random = NewRandom();
        var stream = new EmbossSeqoutStream("msbar", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        random.Seed();
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            EmbossSeqRange.Trim(seq);
            bool isnuc = EmbossSeqType.IsNuc(seq.Type, seq.SeqText);
            string str = "";
            int attempt;
            for (attempt = 0; attempt < attempts; attempt++)
            {
                str = seq.SeqText;
                for (int i = 0; i < count; i++)
                {
                    str = MsbarBlock(random, str, isnuc, blocklist, min, max, inframe);
                    str = MsbarCodon(random, str, isnuc, codonlist, inframe);
                    str = MsbarPoint(random, str, isnuc, pointlist);
                }

                if (MsbarUnlike(str, other, ctx))
                    break;
            }

            if (attempt >= attempts)
                ctx.Warn($"No unique mutation found after {attempts} attempts");
            EmbossProgramBase.AssignSeq(seq, str);
            stream.Write(seq);
        }

        return Finish(stream);
    }

    private static bool MsbarUnlike(string str, EmbossSeqall other, EmbossFeatureReadContext ctx)
    {
        while (EmbossSeqretPrograms.TryNext(other, ctx, out var next))
        {
            if (string.Equals(str, next.SeqText, StringComparison.OrdinalIgnoreCase))
                return false;
        }

        return true;
    }

    // the operation of a list entry: 1 = any of 2..6 (a random one)
    private static int MsbarOperation(EmbossRandom random, string code)
    {
        int opt = int.Parse(code, CultureInfo.InvariantCulture);
        if (opt == 1)
            opt = CInt(random.Double() * 5.0) + 2;
        return opt;
    }

    private static string MsbarBlock(EmbossRandom random, string str, bool isnuc, IReadOnlyList<string> blocklist, int min, int max, bool inframe)
    {
        int imax = max, imin = min;
        if (inframe && imin < 3)
            imin = 3;
        if (imax < imin)
            imax = imin;
        double xlen = unchecked((ulong)(str.Length - 1));                // ajStrGetLen(*str)-1 in size_t arithmetic
        foreach (string code in blocklist)
        {
            if (code == "0")
                return str;
            int opt = MsbarOperation(random, code);
            int rposstart, rposend;
            if (inframe)
            {
                double tda = random.Double() * xlen / 3.0;
                rposstart = CInt(tda);
                double xdiff = imax - imin;
                double tdb = random.Double() * xdiff / 3.0;
                rposend = rposstart + imin / 3 + CInt(tdb) - 1;
                rposstart *= 3;
                rposend *= 3;
                rposend--;
            }
            else
            {
                double tda = random.Double() * xlen;
                rposstart = CInt(tda);
                double xdiff = imax - imin;
                double tdb = random.Double() * xdiff;
                rposend = rposstart + imin + CInt(tdb) - 1;
            }

            switch (opt)
            {
                case 2:
                    str = MsbarInsert(random, str, isnuc, rposstart, rposend);
                    break;
                case 3:
                    str = MsbarCut(str, rposstart, rposend);
                    break;
                case 4:
                    str = MsbarCut(str, rposstart, rposend);
                    str = MsbarInsert(random, str, isnuc, rposstart, rposend);
                    break;
                case 5:
                    str = MsbarMove(str, rposstart, rposend, rposstart);
                    break;
                case 6:
                    {
                        int rpos2;
                        xlen = str.Length;
                        if (inframe)
                            rpos2 = CInt(random.Double() * xlen / 3.0) * 3;
                        else
                            rpos2 = CInt(random.Double() * xlen);

                        str = MsbarMove(str, rposstart, rposend, rpos2);
                        break;
                    }
            }
        }

        return str;
    }

    private static string MsbarCodon(EmbossRandom random, string str, bool isnuc, IReadOnlyList<string> codonlist, bool inframe)
    {
        foreach (string code in codonlist)
        {
            if (code == "0")
                return str;
            int opt = MsbarOperation(random, code);
            double xlen = str.Length;
            int rpos;
            if (inframe)
            {
                double td = random.Double() * xlen / 3.0;
                rpos = CInt(td) * 3;
            }
            else
            {
                double td = random.Double() * xlen;
                rpos = CInt(td);
            }

            switch (opt)
            {
                case 2:
                    str = MsbarInsert(random, str, isnuc, rpos, rpos + 2);
                    break;
                case 3:
                    str = MsbarCut(str, rpos, rpos + 2);
                    break;
                case 4:
                    str = MsbarCut(str, rpos, rpos + 2);
                    str = MsbarInsert(random, str, isnuc, rpos, rpos + 2);
                    break;
                case 5:
                    str = MsbarMove(str, rpos, rpos + 2, rpos);
                    break;
                case 6:
                    {
                        xlen = str.Length;
                        int rpos2;
                        if (inframe)
                            rpos2 = CInt(random.Double() * xlen / 3.0) * 3;
                        else
                            rpos2 = CInt(random.Double() * xlen);
                        str = MsbarMove(str, rpos, rpos + 2, rpos2);
                        break;
                    }
            }
        }

        return str;
    }

    private static string MsbarPoint(EmbossRandom random, string str, bool isnuc, IReadOnlyList<string> pointlist)
    {
        foreach (string code in pointlist)
        {
            if (code == "0")
                return str;
            int opt = MsbarOperation(random, code);
            double xlen = str.Length;
            int rpos = CInt(random.Double() * xlen);
            switch (opt)
            {
                case 2:
                    str = MsbarInsert(random, str, isnuc, rpos, rpos);
                    break;
                case 3:
                    str = MsbarCut(str, rpos, rpos);
                    break;
                case 4:
                    str = MsbarCut(str, rpos, rpos);
                    str = MsbarInsert(random, str, isnuc, rpos, rpos);
                    break;
                case 5:
                    str = MsbarMove(str, rpos, rpos, rpos);
                    break;
                case 6:
                    {
                        xlen = str.Length;
                        int rpos2 = CInt(random.Double() * xlen);
                        str = MsbarMove(str, rpos, rpos, rpos2);
                        break;
                    }
            }
        }

        return str;
    }

    // ajStrCutRange of an empty string: the length of the string becomes -1 (size_t) and the end character is written before the buffer (valgrind: invalid write)
    private static string MsbarCut(string str, long start, long end)
    {
        if (str.Length == 0)
            throw new InvalidOperationException("ajStrCutRange of an empty string writes before the buffer of the string: the reference program has no defined result.");
        return CutRange(str, start, end);
    }

    private static string MsbarInsert(EmbossRandom random, string str, bool isnuc, int start, int end)
    {
        const string nuc = "ACGT";
        const string prot = "ARNDCQEGHILKMFPSTWYV";
        int count = end - start + 1;
        if (count <= 0)
            throw new InvalidOperationException("msbar_Insert with no residues to insert calls ajStrInsertS on a NULL string (a count of 0: segmentation fault) or counts down through the whole range of an int (a negative count): the reference program dies or never ends.");
        var ins = new StringBuilder();
        while (count-- != 0)
        {
            string alphabet = isnuc ? nuc : prot;
            double xlen = alphabet.Length;
            int r = CInt(random.Double() * xlen);
            ins.Append(alphabet[r]);
            if (ins.Length > 100_000_000)
                throw new InvalidOperationException("The insertion is longer than the reference program can hold.");
        }

        return Insert(str, start, ins.ToString());
    }

    private static string MsbarMove(string str, long start, long end, long destination) => Insert(str, destination, KeepRange(str, start, end));

    // ---------------------------------------------------------------------------------------------------------------- makenucseq, makeprotseq

    // the message of ajFatal: a blank line, the program source and line, the text
    private static EmbossRunResult FatalAt(EmbossFeatureReadContext ctx, string source, int line, string text) =>
        EmbossProgramBase.Fatal(ctx, $"\n   EMBOSS An error in {source} at line {line}:\n{text}");

    // makenucseq_random_sequence / makeprotseq_random_sequence: the strings of the distribution picked by a random number
    private static string RandomSequence(EmbossRandom random, IReadOnlyList<string> seqchar, int scmax, int length)
    {
        var seq = new StringBuilder();
        int len = length;
        while (len-- > 0)
        {
            int idx = CInt(random.Double() * scmax);
            seq.Append(seqchar[idx]);
        }

        return seq.ToString();
    }

    /// <summary><c>makenucseq</c>: random nucleotide sequences (<c>-amount</c> of <c>-length</c>, an optional insert; with a codon usage file the sequence is made of codons with the frequencies of the file).</summary>
    internal static EmbossRunResult Makenucseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files) => MakeSeq(a, files, "makenucseq", true);

    /// <summary><c>makeprotseq</c>: random protein sequences (<c>-amount</c> of <c>-length</c>, an optional insert; with a pepstats file the composition of its first sequence).</summary>
    internal static EmbossRunResult Makeprotseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files) => MakeSeq(a, files, "makeprotseq", false);

    private static EmbossRunResult MakeSeq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files, string program, bool nuc)
    {
        var ctx = EmbossProgramBase.Context(program);
        string dataName = a.Raw(nuc ? "codonfile" : "pepstatsfile") ?? "";
        string? data = null;
        EmbossCodonUsage? codon = null;
        string insert;
        int start, length, amount;
        EmbossOutseqOptions outOptions;
        try
        {
            if (dataName.Length > 0)
            {
                if (nuc)
                {
                    codon = EmbossCodonUsage.Read(dataName, files);
                    if (codon is null)
                        throw a.Die("codonfile", $"Error: Unable to read codon usage '{dataName}'");
                }
                else
                {
                    if (!files.TryGetValue(dataName, out data))
                        throw a.Die("pepstatsfile", $"Error: Unable to open file '{dataName}' for input");
                }
            }

            amount = a.IntLimited("amount", 100, 1, int.MaxValue, ctx);
            length = a.IntLimited("length", 100, 1, int.MaxValue, ctx);
            a.Bool("useinsert", false);
            insert = a.String("insert", "");
            start = a.IntLimited("start", 1, 1, int.MaxValue, ctx);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        string source = program + ".c";
        var list = new List<string>();
        int extra = 0;
        var stream = new EmbossSeqoutStream(program, ctx, outOptions, nuc ? "nucleotide" : "protein");
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput());
        if (insert.Length > 0)
        {
            length -= insert.Length;
            start--;
            if (length <= 0)
                return FatalAt(ctx, source, nuc ? 77 : 80, "Sequence smaller than inserted part. Unable to create sequences.");
        }

        if (nuc)
        {
            if (codon is not null)
            {
                codon.AppendCodonList(list);
                extra = length % 3;
                length /= 3;
                if (extra != 0)
                {
                    length++;
                    extra = 0 - (3 - extra);
                }
            }
            else
            {
                foreach (char c in "ACGTacgt")
                    list.Add(c.ToString());
            }
        }
        else if (data is not null)
        {
            var lines = ReadLines(data).Select(l => l.TrimEnd('\n').TrimEnd('\r')).ToList();
            string first = lines.Count > 0 ? lines[0] : "";
            if (first.StartsWith("PEPSTATS", StringComparison.Ordinal))
                ParsePepstats(list, lines, ctx);
            else
            {
                ctx.Warn("Not pepstats file. Making completely random sequences.");
                foreach (char c in "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
                    list.Add(c.ToString());
            }
        }
        else
        {
            foreach (char c in "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
                list.Add(c.ToString());
        }

        int scmax = list.Count;
        if (scmax == 0)
            return FatalAt(ctx, source, nuc ? 114 : 123, "No strings in list. No characters to make the sequence.");
        var random = NewRandom();
        random.Seed();
        while (amount-- > 0)
        {
            string seqstr = RandomSequence(random, list, scmax, length);
            if (insert.Length > 0)
                seqstr = Insert(seqstr, start, insert);
            seqstr = seqstr.ToLowerInvariant();
            if (nuc)
            {
                seqstr = seqstr.Replace('u', 't');
                if (extra < 0)
                    seqstr = seqstr[..Math.Max(0, seqstr.Length + extra)];
            }

            var seq = EmbossProgramBase.NewSeq("", "", "", "");
            EmbossProgramBase.AssignSeq(seq, seqstr);
            if (nuc)
                seq.SetNuc();
            else
                seq.SetProt();
            stream.Write(seq);
        }

        return Finish(stream);
    }

    // makeprotseq_parse_pepstats: the lines after the line "Residue ..." until an empty line; the fifth column (Mole%) gives the number of copies of the residue (lower case) in the distribution (pushed at the head of the list)
    private static void ParsePepstats(List<string> list, List<string> lines, EmbossFeatureReadContext ctx)
    {
        int i = 1;
        while (i < lines.Count)
        {
            if (lines[i++].StartsWith("Residue", StringComparison.Ordinal))
                break;
        }

        var parse = new EmbossStrParse(ctx);
        for (; i < lines.Count; i++)
        {
            string line = lines[i];
            if (ParseCount(line, [' ', '\t']) == 0)
                break;
            string? tok = parse.Next(line);
            string ch = (tok ?? "").ToLowerInvariant();
            for (int count = 1; count < 5; count++)
                tok = parse.Next(null);
            double value = 0.0;
            if (tok is not null && EmbossStrtod.TryParse(tok, out double v))
                value = v;
            int copies = (int)(value * 100) + (value - (int)value >= 0.5 ? 1 : 0);
            for (int k = 0; k < copies; k++)
                list.Insert(0, ch);
        }
    }

    // ajStrParseCountC: the steps of the tokenizer (a run of non-delimiters, then the delimiters after it)
    private static int ParseCount(string s, char[] delimiters)
    {
        if (s.Length == 0)
            return 0;
        int pos = 0, count = 0;
        while (pos < s.Length)
        {
            count++;
            while (pos < s.Length && Array.IndexOf(delimiters, s[pos]) < 0)
                pos++;
            while (pos < s.Length && Array.IndexOf(delimiters, s[pos]) >= 0)
                pos++;
        }

        return count;
    }
}
