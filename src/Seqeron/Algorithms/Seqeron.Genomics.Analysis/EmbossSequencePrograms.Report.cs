#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The programs of EMBOSS 6.6.0 that write a report about the sequences they read: <c>infoseq</c>, <c>wordcount</c>, <c>seqxref</c> and the display programs (B09 WP62, F281).</summary>
internal static partial class EmbossSequencePrograms
{
    // ---------------------------------------------------------------------------------------------------------------- infoseq

    // the static state of infoseq_printheader / infoseq_print: the delimiter is left out of the first token of a line, a column that was full is followed by a blank
    private sealed class InfoseqLine
    {
        public bool HeaderNoDelim = true, HeaderColFull;
        public bool NoDelim = true, ColFull;
    }

    private static bool InfoseqPrintHeader(StringBuilder o, InfoseqLine st, bool html, bool instring, string text, int wid, bool columns, string delimiter)
    {
        if (!instring)
        {
            st.HeaderNoDelim = true;
            st.HeaderColFull = false;
        }

        if (html)
            o.Append("<th>").Append(text).Append("</th>");
        else if (columns)
        {
            if (st.HeaderColFull)
                o.Append(' ');
            o.Append(text.PadRight(wid));
            st.HeaderColFull = text.Length >= wid;
        }
        else
        {
            if (!st.HeaderNoDelim)
                o.Append(delimiter);
            o.Append(text);
        }

        st.HeaderNoDelim = false;
        return true;
    }

    private static bool InfoseqPrint(StringBuilder o, InfoseqLine st, bool html, bool instring, string str, bool usewid, int wid, bool columns, string delimiter)
    {
        if (!instring)
        {
            st.NoDelim = true;
            st.ColFull = false;
        }

        if (html)
            o.Append("<td>").Append(str).Append("</td>");
        else if (columns)
        {
            if (st.ColFull)
                o.Append(' ');
            if (usewid)
            {
                o.Append(str.PadRight(wid));
                st.ColFull = str.Length >= wid;
            }
            else
            {
                o.Append(str);
            }
        }
        else
        {
            if (!st.NoDelim)
                o.Append(delimiter);
            o.Append(str);
        }

        st.NoDelim = false;
        return true;
    }

    /// <summary><c>ajMeltGC</c>: the fraction of G, C and S (ambiguity codes count their share of G or C) of the first <paramref name="len"/> characters.</summary>
    internal static float MeltGC(string s, int len)
    {
        double count = 0.0;
        for (int i = 0; i < len; ++i)
        {
            char t = char.ToUpperInvariant(s[i]);
            if ("GCS".Contains(t, StringComparison.Ordinal))
                count++;
            else if ("RYMK".Contains(t, StringComparison.Ordinal) || "NX".Contains(t, StringComparison.Ordinal))
                count += 0.5;
            else if ("BV".Contains(t, StringComparison.Ordinal))
                count += 0.6666667;
            else if ("DH".Contains(t, StringComparison.Ordinal))
                count += 0.3333333;
        }

        return (float)(count / len);
    }

    /// <summary><c>infoseq</c>: the USA, database, name, accession, GI, version, type, length, GC content, organism and description of every sequence in columns (or delimited, or as an HTML table).</summary>
    internal static EmbossRunResult Infoseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("infoseq");
        var all = EmbossProgramBase.Open("infoseq", a, files, ctx, "gapany");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        bool columns, html, only, doheader, dousa, dodb, doname, doacc, dogi, dosv, dotype, dolength, dopgc, doorg, dodesc;
        string delimiter, outfile;
        try
        {
            columns = a.Bool("columns", true);
            delimiter = a.String("delimiter", "|");
            outfile = a.Raw("outfile") ?? "stdout";
            html = a.Bool("html", false);
            only = a.Bool("only", false);
            doheader = a.Bool("heading", true);
            dousa = a.Bool("usa", !only);
            dodb = a.Bool("database", !only);
            doname = a.Bool("name", !only);
            doacc = a.Bool("accession", !only);
            dogi = a.Bool("gi", false);
            dosv = a.Bool("seqversion", false);
            dotype = a.Bool("type", !only);
            dolength = a.Bool("length", !only);
            dopgc = a.Bool("pgc", !only);
            doorg = a.Bool("organism", !only);
            dodesc = a.Bool("description", !only);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (delimiter == "\\t")
            delimiter = "\t";
        var o = new StringBuilder();
        var st = new InfoseqLine();
        bool instring = false, firsttime = true;
        if (html)
            o.Append("<table border cellpadding=4 bgcolor=\"#FFFFF0\">\n");
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            EmbossSeqRange.Trim(seq);
            string text = seq.SeqText;
            bool type = EmbossSeqType.IsProt(seq.Type, text);
            if (firsttime)
            {
                if (doheader)
                {
                    if (html)
                        o.Append("<tr>");
                    if (dousa)
                        instring = InfoseqPrintHeader(o, st, html, instring, "USA", 25, columns, delimiter);
                    if (dodb)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Database", 10, columns, delimiter);
                    if (doname)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Name", 15, columns, delimiter);
                    if (doacc)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Accession", 15, columns, delimiter);
                    if (dogi)
                        instring = InfoseqPrintHeader(o, st, html, instring, "GI", 15, columns, delimiter);
                    if (dosv)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Version", 8, columns, delimiter);
                    if (dotype)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Type", 5, columns, delimiter);
                    if (dolength)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Length", 7, columns, delimiter);
                    if (!type && dopgc)
                        instring = InfoseqPrintHeader(o, st, html, instring, "%GC", 7, columns, delimiter);
                    if (doorg)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Organism", 20, columns, delimiter);
                    if (dodesc)
                        instring = InfoseqPrintHeader(o, st, html, instring, "Description", 12, columns, delimiter);
                    o.Append(html ? "</tr>\n" : "\n");
                    instring = false;
                }

                firsttime = false;
            }

            string usa = UsaOf(seq);
            if (usa.Length == 0)
                usa = "-";
            string db = seq.Db.Length == 0 ? "-" : seq.Db;
            string name = seq.Name.Length == 0 ? "-" : seq.Name;
            string acc = seq.Acc.Length == 0 ? "-" : seq.Acc;
            string gi = seq.Gi.Length == 0 ? "-" : seq.Gi;
            string sv = seq.Sv.Length == 0 ? "-" : seq.Sv;
            int length = text.Length;
            float pgc = 0.0f;
            if (dopgc && !type)
            {
                pgc = MeltGC(text, length);
                pgc *= 100;
            }

            if (html)
                o.Append("<tr>");
            if (dousa)
                instring = InfoseqPrint(o, st, html, instring, usa, true, 25, columns, delimiter);
            if (dodb)
                instring = InfoseqPrint(o, st, html, instring, db, true, 15, columns, delimiter);
            if (doname)
                instring = InfoseqPrint(o, st, html, instring, name, true, 15, columns, delimiter);
            if (doacc)
                instring = InfoseqPrint(o, st, html, instring, acc, true, 15, columns, delimiter);
            if (dogi)
                instring = InfoseqPrint(o, st, html, instring, gi, true, 15, columns, delimiter);
            if (dosv)
                instring = InfoseqPrint(o, st, html, instring, sv, true, 8, columns, delimiter);
            if (dotype)
                instring = InfoseqPrint(o, st, html, instring, type ? "P" : "N", true, 5, columns, delimiter);
            if (dolength)
                instring = InfoseqPrint(o, st, html, instring, length.ToString(CultureInfo.InvariantCulture), true, 7, columns, delimiter);
            if (!type && dopgc)
                instring = InfoseqPrint(o, st, html, instring, CFormat.Fixed(pgc, 2), true, 7, columns, delimiter);
            if (doorg)
                instring = InfoseqPrint(o, st, html, instring, seq.Tax, true, 20, columns, delimiter);
            if (dodesc)
                instring = InfoseqPrint(o, st, html, instring, seq.Desc, false, 0, columns, delimiter);
            o.Append(html ? "</tr>\n" : "\n");
            instring = false;
        }

        if (html)
            o.Append("</table>\n");
        return OutfileResult(ctx, outfile, o.ToString());
    }

    // ---------------------------------------------------------------------------------------------------------------- wordcount

    /// <summary><c>wordcount</c>: the words of <c>-wordsize</c> letters of all the sequences with their counts (at least <c>-mincount</c>), the most frequent first.</summary>
    internal static EmbossRunResult Wordcount(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("wordcount");
        var all = EmbossProgramBase.Open("wordcount", a, files, ctx, "any");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var first = all.Peek();
        int wordsize, mincount;
        string outfile;
        try
        {
            wordsize = a.IntLimited("wordsize", EmbossSeqType.IsProt(first.Type, first.SeqText) ? 2 : 4, 1, int.MaxValue, ctx);
            mincount = a.IntLimited("mincount", 1, 1, int.MaxValue, ctx);
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "wordcount", "sequence");
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        EmbossWordTable? table = null;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
            EmbossWordTable.GetTable(ref table, seq.SeqText, EmbossSeqType.IsNuc(seq.Type, seq.SeqText), wordsize);
        return OutfileResult(ctx, outfile, table is null || table.Length == 0 ? "" : table.Print(mincount));
    }

    // ---------------------------------------------------------------------------------------------------------------- seqxref

    // ajStrCmpCaseS: the comparison of the upper case texts (byte by byte, the shorter first)
    private static int CmpCase(string a, string b) => string.CompareOrdinal(a.ToUpperInvariant(), b.ToUpperInvariant()) switch { < 0 => -1, > 0 => 1, _ => 0 };

    /// <summary><c>seqxref</c>: the database cross-references of every sequence entry (<c>db:id</c>, sorted by database name, identifier ..., equal ones once), after <c>#USA: count</c>.</summary>
    internal static EmbossRunResult Seqxref(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("seqxref");
        var all = EmbossProgramBase.Open("seqxref", a, files, ctx, "any", features: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        string outfile;
        try
        {
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "seqxref", "sequence");
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var o = new StringBuilder();
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            var xrefs = seq.Xreflist.ToList();
            int nrefs = xrefs.Count;
            // ajSeqxreflistSort: the merge sort of qsort (stable) and the removal of an item equal to the one kept before it
            Comparison<EmbossSeqXref> cmp = (x, y) =>
            {
                int r = CmpCase(x.Db, y.Db);
                if (r != 0)
                    return r;
                r = CmpCase(x.Id, y.Id);
                if (r != 0)
                    return r;
                r = CmpCase(x.Secid, y.Secid);
                if (r != 0)
                    return r;
                r = CmpCase(x.Terid, y.Terid);
                return r != 0 ? r : CmpCase(x.Quatid, y.Quatid);
            };
            var sorted = xrefs.Select((x, i) => (x, i)).OrderBy(t => t, Comparer<(EmbossSeqXref x, int i)>.Create((l, r) => cmp(l.x, r.x))).Select(t => t.x).ToList();
            var unique = new List<EmbossSeqXref>();
            foreach (var x in sorted)
            {
                if (unique.Count > 0 && cmp(x, unique[^1]) == 0)
                    continue;
                unique.Add(x);
            }

            o.Append('#').Append(UsaOf(seq)).Append(": ").Append(nrefs.ToString(CultureInfo.InvariantCulture)).Append('\n');
            foreach (var x in unique)
                o.Append(x.Db).Append(':').Append(x.Id).Append('\n');
        }

        return OutfileResult(ctx, outfile, o.ToString());
    }


    // the write of an element of the array of protein positions: an index outside the array is a heap overflow (valgrind: invalid write), the result of the reference program is not defined
    private static void PutPos(int[] ppos, int pos, int len, int value)
    {
        if (pos < 0 || pos >= len)
            throw new InvalidOperationException("The program writes an element of its array of positions outside the array (valgrind: invalid write): the reference program has no defined result.");
        ppos[pos] = value;
    }

    // the read of an element of the array of protein positions in the order the C code makes it: an index outside the array (ajint ppos[len]) is a read of memory the program does not own (valgrind: invalid read, a value that is not defined)
    private static int PosAt(int[] ppos, int pos, int len)
    {
        if (pos < 0 || pos >= len)
            throw new InvalidOperationException("The program reads an element of its array of positions outside the array (valgrind: invalid read): the reference program has no defined result.");
        return ppos[pos];
    }

    // ---------------------------------------------------------------------------------------------------------------- showorf

    private const string TableValues = "0:Standard; 1:Standard (with alternative initiation codons); 2:Vertebrate Mitochondrial; 3:Yeast Mitochondrial; 4:Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma; 5:Invertebrate Mitochondrial; 6:Ciliate Macronuclear and Dasycladacean; 9:Echinoderm Mitochondrial; 10:Euplotid Nuclear; 11:Bacterial; 12:Alternative Yeast Nuclear; 13:Ascidian Mitochondrial; 14:Flatworm Mitochondrial; 15:Blepharisma Macronuclear; 16:Chlorophycean Mitochondrial; 21:Trematode Mitochondrial; 22:Scenedesmus obliquus; 23:Thraustochytrium Mitochondrial";

    // ajBaseAlphaToBin (Ebases.iub) and ajCodBase: the index of the first base a code stands for
    private static int CodBase(char c)
    {
        int v = char.ToUpperInvariant(c) switch
        {
            'A' => 1, 'B' => 14, 'C' => 2, 'D' => 13, 'G' => 4, 'H' => 11, 'K' => 12, 'M' => 3, 'N' => 15, 'R' => 5, 'S' => 6, 'T' => 8, 'U' => 8, 'V' => 7, 'W' => 9, 'X' => 15, 'Y' => 10, _ => 0,
        };
        return (v & 1) != 0 ? 0 : (v & 2) != 0 ? 1 : (v & 4) != 0 ? 2 : (v & 8) != 0 ? 3 : 0;
    }

    /// <summary><c>showorf</c>: a nucleotide sequence in blocks of <c>-width</c> with the translations of the chosen frames (the amino acid of the first base of every codon, a ruler and the numbers).</summary>
    internal static EmbossRunResult Showorf(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("showorf");
        var all = EmbossProgramBase.Open("showorf", a, files, ctx, "dna", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        IReadOnlyList<string> frames;
        string codestr, outfile;
        bool isrule, isp, isn;
        int width;
        try
        {
            frames = a.List("frames", "1,2,3,4,5,6", "0: None, 1: F1,2: F2,3: F3,4: R1,5: R2,6: R3", 0, 6, ",");
            codestr = a.List("table", "0", TableValues, 1, 1)[0];
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "showorf", "sequence");
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
            isrule = a.Bool("ruler", true);
            isp = a.Bool("plabel", true);
            isn = a.Bool("nlabel", true);
            width = a.IntLimited("width", 50, 10, int.MaxValue, ctx);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var seq = all.Peek();
        var mark = new bool[6];
        foreach (string f in frames)
        {
            int v = int.Parse(f.Trim(), CultureInfo.InvariantCulture);
            if (v == 0)
            {
                Array.Clear(mark);
                break;
            }

            mark[v - 1] = true;
        }

        var trn = EmbossTrn.New(int.Parse(codestr, CultureInfo.InvariantCulture));
        int beg = EmbossSeqRange.GetBegin(seq), end = EmbossSeqRange.GetEnd(seq);
        string substr = EmbossProgramBase.SubS(seq.SeqText, beg - 1, end - 1);
        int len = substr.Length;
        string revstr = EmbossNuc.ReverseComplement(substr);
        var pseqs = new char[6][];
        var ppos = new int[6][];
        for (int i = 0; i < 6; i++)
        {
            pseqs[i] = substr.ToCharArray();
            ppos[i] = new int[len + 4];
        }

        // showorf_SixTranslate / showorf_DoTrans
        for (int n = 0; n < 6; n++)
        {
            string p;
            int po;
            if (n < 3)
            {
                p = substr;
                po = n % 3;
            }
            else
            {
                p = revstr;
                po = len % 3;
                po -= n % 3;
                if (po < 0)
                    po += 3;
            }

            var q = pseqs[n];
            if (po > len)
                throw new InvalidOperationException("showorf_DoTrans pads the frame with po blanks beyond the end of a string shorter than po (it overwrites the end of the string): the reference program has no defined result.");
            for (int i = 0; i < po; ++i)
                q[i] = ' ';
            int c = 0;
            for (int i = po; i < len; ++i, ++c)
            {
                if (c % 3 != 0 || len - i < 3)
                {
                    q[i] = ' ';
                    continue;
                }

                string tri = p.Substring(i, 3).ToUpperInvariant();
                if (tri == "NNN")
                {
                    q[i] = 'X';
                    continue;
                }

                int b1 = CodBase(tri[0]), b2 = CodBase(tri[1]), b3 = CodBase(tri[2]);
                const string acgt = "ACGT";
                char aa = trn.Codon(acgt[b1], acgt[b2], acgt[b3]);
                q[i] = aa == '*' ? '*' : char.ToUpperInvariant(aa);
            }

            if (n > 2)
                Array.Reverse(q);
        }

        // showorf_MakeRuler
        var ruler = new char[len];
        var npos = new int[len];
        for (int i = 0; i < len; ++i)
        {
            npos[i] = i + beg;
            ruler[i] = (i + beg) % 10 == 0 ? '|' : '-';
        }

        const int POFF = 1000000;
        char P(int f, int pos) => pos >= 0 && pos < len ? pseqs[f][pos] : '\0';
        // showorf_CalcProteinPos
        for (int i = 0; i < 3; ++i)
        {
            int pos = 0, v = 1;
            while (P(i, pos) == ' ')
                PutPos(ppos[i], pos++, len, 0);
            while (pos < len)
            {
                if (P(i, pos) == '*')
                {
                    PutPos(ppos[i], pos, len, 0);
                    ++pos;
                    while (P(i, pos) == ' ')
                    {
                        PutPos(ppos[i], pos, len, 0);
                        ++pos;
                    }

                    v = 1;
                    continue;
                }

                if (P(i, pos) != ' ')
                {
                    PutPos(ppos[i], pos, len, v + POFF);
                    ++pos;
                    for (int j = 0; j < 2 && P(i, pos) != '\0'; ++j, ++pos)
                        PutPos(ppos[i], pos, len, v);
                    if (P(i, pos) == ' ')
                        v = 1;
                    else
                        ++v;
                }
                else
                {
                    PutPos(ppos[i], pos++, len, 0);
                }
            }
        }

        for (int i = 3; i < 6 && len > 0; ++i)
        {
            int pos = len - 1, v = 1;
            while (P(i, pos) == ' ')
                PutPos(ppos[i], pos--, len, 0);
            while (pos > -1)
            {
                if (P(i, pos) == '*')
                {
                    PutPos(ppos[i], pos, len, 0);
                    --pos;
                    while (P(i, pos) == ' ')
                    {
                        PutPos(ppos[i], pos, len, 0);
                        --pos;
                    }

                    v = 1;
                    continue;
                }

                if (P(i, pos) != ' ')
                {
                    PutPos(ppos[i], pos, len, v + POFF);
                    --pos;
                    for (int j = 0; j < 2 && P(i, pos) != '\0'; ++j, --pos)
                        PutPos(ppos[i], pos, len, v);
                    if (pos < 0)
                        continue;
                    if (P(i, pos) == ' ')
                        v = 1;
                    else
                        ++v;
                }
                else
                {
                    PutPos(ppos[i], pos--, len, 0);
                }
            }
        }

        var o = new StringBuilder();
        o.Append("SHOWORF of ").Append(seq.Name).Append(" from ").Append(beg.ToString(CultureInfo.InvariantCulture)).Append(" to ").Append((beg + len - 1).ToString(CultureInfo.InvariantCulture)).Append("\n\n");
        int at = 0;
        while (at < len)
        {
            if (at + width < len)
            {
                ShoworfBlock(o, ppos, len, npos, pseqs, substr, mark, ruler, isrule, isp, isn, at, at + width - 1);
                at += width;
                continue;
            }

            ShoworfBlock(o, ppos, len, npos, pseqs, substr, mark, ruler, isrule, isp, isn, at, len - 1);
            break;
        }

        return OutfileResult(ctx, outfile, o.ToString());
    }

    private static void ShoworfBlock(StringBuilder o, int[][] ppos, int len, int[] npos, char[][] pseqs, string substr, bool[] mark, char[] ruler, bool isrule, bool isp, bool isn, int start, int end)
    {
        string[] fr = ["F1", "F2", "F3", "R1", "R2", "R3"];
        const int POFF = 1000000;
        string Sub(char[] c) => new(c, start, end - start + 1);
        if (isrule)
            o.Append("           ").Append(new string(ruler, start, end - start + 1)).Append('\n');
        if (isn)
            o.Append(npos[start].ToString(CultureInfo.InvariantCulture).PadLeft(10)).Append(' ');
        else
            o.Append("           ");
        o.Append(substr.AsSpan(start, end - start + 1)).Append(' ');
        if (isn)
            o.Append(npos[end].ToString(CultureInfo.InvariantCulture));
        o.Append('\n');
        for (int i = 0; i < 3; ++i)
        {
            if (!mark[i])
                continue;
            o.Append(fr[i]).Append(' ');
            int b = 0, e = 0;
            if (isp)
            {
                int pos = start, v;
                while (pos <= end)
                {
                    if ((v = ppos[i][pos]) == 0)
                    {
                        ++pos;
                        continue;
                    }

                    if (v < POFF)
                    {
                        while (PosAt(ppos[i], pos, len) == v && pos <= end)
                            ++pos;
                        if (pos > end)
                            break;
                        continue;
                    }

                    b = v - POFF;
                    break;
                }

                pos = end;
                while (pos >= start)
                {
                    if ((v = ppos[i][pos]) == 0)
                    {
                        --pos;
                        continue;
                    }

                    if (v > POFF)
                    {
                        v -= POFF;
                    }
                    else
                    {
                        while (PosAt(ppos[i], pos, len) == v)
                        {
                            --pos;
                            if (pos < start)
                            {
                                v = 0;
                                break;
                            }
                        }
                    }

                    e = v;
                    break;
                }

                o.Append(b == 0 ? "        " : b.ToString(CultureInfo.InvariantCulture).PadLeft(7) + " ");
            }
            else
            {
                o.Append("        ");
            }

            o.Append(Sub(pseqs[i])).Append(' ');
            if (isp && e != 0)
                o.Append(e.ToString(CultureInfo.InvariantCulture));
            o.Append('\n');
        }

        ShoworfReverse(o, ppos, len, pseqs, mark, isp, start, end);
    }


    // the reverse frames of showorf_showTransb
    private static void ShoworfReverse(StringBuilder o, int[][] ppos, int len, char[][] pseqs, bool[] mark, bool isp, int start, int end)
    {
        string[] fr = ["F1", "F2", "F3", "R1", "R2", "R3"];
        const int POFF = 1000000;
        for (int i = 3; i < 6; ++i)
        {
            if (!mark[i])
                continue;
            o.Append(fr[i]).Append(' ');
            int b = 0, e = 0;
            if (isp)
            {
                int pos = start, v;
                while (pos <= end)
                {
                    if ((v = ppos[i][pos]) == 0)
                    {
                        ++pos;
                        continue;
                    }

                    if (v > POFF)
                    {
                        v -= POFF;
                    }
                    else
                    {
                        while (PosAt(ppos[i], pos, len) == v)
                        {
                            ++pos;
                            if (pos > end)
                            {
                                v = 0;
                                break;
                            }
                        }
                    }

                    b = v;
                    break;
                }

                pos = end;
                while (pos >= start)
                {
                    if ((v = ppos[i][pos]) == 0)
                    {
                        --pos;
                        continue;
                    }

                    if (v < POFF)
                    {
                        while (PosAt(ppos[i], pos, len) == v && pos >= start)
                            --pos;
                        if (pos < start)
                            break;
                        continue;
                    }

                    e = v - POFF;
                    break;
                }

                o.Append(b == 0 ? "        " : b.ToString(CultureInfo.InvariantCulture).PadLeft(7) + " ");
            }
            else
            {
                o.Append("        ");
            }

            o.Append(new string(pseqs[i], start, end - start + 1)).Append(' ');
            if (isp && e != 0)
                o.Append(e.ToString(CultureInfo.InvariantCulture));
            o.Append('\n');
        }

        o.Append('\n');
    }

    // ---------------------------------------------------------------------------------------------------------------- prettyseq

    /// <summary><c>prettyseq</c>: a nucleotide sequence in blocks of <c>-width</c> with the translation of the ranges (the amino acid under the first base of every codon of the ranges, written in lower case in the sequence).</summary>
    internal static EmbossRunResult Prettyseq(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("prettyseq");
        var all = EmbossProgramBase.Open("prettyseq", a, files, ctx, "dna", single: true);
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var seq = all.Peek();
        int beg = EmbossSeqRange.GetBegin(seq), end = EmbossSeqRange.GetEnd(seq);
        IReadOnlyList<EmbossRangeElement> range;
        string codestr, outfile;
        bool isrule, isp, isn;
        int width;
        try
        {
            if (a.Raw("range") is null)
                range = EmbossRangeSpec.Parse($"{beg}-{end}");
            else
                range = Ranges(a, "range", ctx);
            codestr = a.List("table", "0", TableValues, 1, 1)[0];
            isrule = a.Bool("ruler", true);
            isp = a.Bool("plabel", true);
            isn = a.Bool("nlabel", true);
            width = a.IntLimited("width", 60, 10, int.MaxValue, ctx);
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "prettyseq", "sequence");
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var trn = EmbossTrn.New(int.Parse(codestr, CultureInfo.InvariantCulture));
        string substr = EmbossProgramBase.SubS(seq.SeqText, beg - 1, end - 1).ToUpperInvariant();
        int len = substr.Length;
        // prettyseq_Translate: the ranges in sub-sequence positions (ajuint), the area of the ranges in lower case, the translation under the first base of every codon of lower case bases
        int nr = range.Count;
        var starts = new int[nr];
        var ends = new int[nr];
        for (int i = 0; i < nr; i++)
        {
            starts[i] = range[i].Start - beg;
            ends[i] = range[i].End - beg;
        }

        uint limit = (uint)len;
        uint st = 0, en = 0;
        for (int i = 0; i < nr; i++)
        {
            st = unchecked((uint)starts[i]);
            en = unchecked((uint)ends[i]);
            if (st >= limit || en >= limit)
            {
                var died = FatalAt(ctx, "prettyseq.c", 191, $"Range outside length of sequence [{unchecked((int)(st + (uint)beg)).ToString(CultureInfo.InvariantCulture)}-{(end + beg).ToString(CultureInfo.InvariantCulture)}]");
                // the output file was opened by the ACD before the program died
                return outfile == "stdout" ? died : died with { Files = new Dictionary<string, string> { [outfile] = "" } };
            }
        }

        var p = new char[len + 8];
        substr.CopyTo(0, p, 0, len);
        for (int i = 0; i < nr; i++)
        {
            // ajRangeElementGetValues refuses a start below 1 (a range that starts at the first base), an end below 1 and a start after the end: st and en stay what they were
            if (starts[i] >= 1 && ends[i] >= 1 && starts[i] <= ends[i])
            {
                st = (uint)starts[i];
                en = (uint)ends[i];
            }

            for (uint j = st; j <= en; ++j)
                p[j] = char.ToLowerInvariant(p[j]);
        }

        if (len == 1)
            throw new InvalidOperationException("prettyseq_Translate computes limit = 3 * ((len + 1) / 3) - 1 in unsigned arithmetic: a sequence of one residue gives 4294967295 and the loop runs over the whole memory (the reference program dies by a fault).");
        limit = (uint)(3 * ((len + 1) / 3) - 1);
        var q = new char[len + 8];
        substr.CopyTo(0, q, 0, len);
        uint iq;
        for (iq = 0; iq < limit; ++iq)
        {
            if (char.IsAsciiLetterUpper(p[iq]))
            {
                q[iq] = ' ';
                continue;
            }

            if (p[iq + 1] == '\0' || char.IsAsciiLetterUpper(p[iq + 1]) || p[iq + 2] == '\0' || char.IsAsciiLetterUpper(p[iq + 2]))
            {
                q[iq] = ' ';
                if (q[iq + 1] != '\0')
                    q[iq + 1] = ' ';
                if (q[iq + 2] != '\0')
                    q[iq + 2] = ' ';
                iq += 2;
                continue;
            }

            string tri = new([char.ToUpperInvariant(p[iq]), char.ToUpperInvariant(p[iq + 1]), char.ToUpperInvariant(p[iq + 2])]);
            char c;
            if (tri == "NNN")
            {
                c = 'X';
            }
            else
            {
                const string acgt = "ACGT";
                char aa = trn.Codon(acgt[CodBase(tri[0])], acgt[CodBase(tri[1])], acgt[CodBase(tri[2])]);
                c = aa == '*' ? '*' : char.ToUpperInvariant(aa);
            }

            q[iq] = c;
            q[iq + 1] = q[iq + 2] = ' ';
            iq += 2;
        }

        uint jj = (uint)len - limit;
        while (jj-- != 0)
            q[iq++] = ' ';
        string pro = new(q, 0, len);
        string text = new(p, 0, len);

        // prettyseq_makeRuler, prettyseq_calcProteinPos
        const int POFF = 1000000;
        var ruler = new char[len];
        var npos = new int[len];
        for (int i = 0; i < len; ++i)
        {
            npos[i] = i + beg;
            ruler[i] = (i + beg) % 10 == 0 ? '|' : '-';
        }

        var ppos = new int[len + 8];
        char Q(int pos) => pos >= 0 && pos < q.Length ? q[pos] : '\0';
        void Put(int pos, int value)
        {
            if (pos >= len)
                throw new InvalidOperationException("prettyseq_calcProteinPos writes past the end of the array of positions (heap overflow, valgrind: invalid write): the reference program has no defined result.");
            ppos[pos] = value;
        }

        int pp = 0, v0 = 1;
        while (Q(pp) == ' ')
            Put(pp++, 0);
        while (pp < len)
        {
            if (Q(pp) == '*')
            {
                Put(pp, 0);
                ++pp;
                while (pp < len && Q(pp) == ' ')
                {
                    Put(pp, 0);
                    ++pp;
                }

                v0 = 1;
                continue;
            }

            if (Q(pp) != ' ')
            {
                Put(pp, v0 + POFF);
                ++pp;
                for (int j = 0; j < 2 && Q(pp) != '\0'; ++j, ++pp)
                    Put(pp, v0);
                if (Q(pp) == ' ')
                    v0 = 1;
                else
                    ++v0;
            }
            else
            {
                Put(pp++, 0);
            }
        }

        var o = new StringBuilder();
        o.Append("PRETTYSEQ of ").Append(seq.Name).Append(" from ").Append(beg.ToString(CultureInfo.InvariantCulture)).Append(" to ").Append((beg + len - 1).ToString(CultureInfo.InvariantCulture)).Append("\n\n");
        int at = 0;
        while (at < len)
        {
            if (at + width < len)
            {
                PrettyBlock(o, ppos, len, npos, pro, text, ruler, isrule, isp, isn, at, at + width - 1);
                at += width;
                continue;
            }

            PrettyBlock(o, ppos, len, npos, pro, text, ruler, isrule, isp, isn, at, len - 1);
            break;
        }

        return OutfileResult(ctx, outfile, o.ToString());
    }

    private static void PrettyBlock(StringBuilder o, int[] ppos, int len, int[] npos, string pro, string substr, char[] ruler, bool isrule, bool isp, bool isn, int start, int end)
    {
        const int POFF = 1000000;
        if (isrule)
            o.Append("           ").Append(new string(ruler, start, end - start + 1)).Append('\n');
        if (isn)
            o.Append(npos[start].ToString(CultureInfo.InvariantCulture).PadLeft(10)).Append(' ');
        else
            o.Append("           ");
        o.Append(substr.AsSpan(start, end - start + 1)).Append(' ');
        if (isn)
            o.Append(npos[end].ToString(CultureInfo.InvariantCulture));
        o.Append('\n');
        int b = 0, e = 0;
        if (isp)
        {
            int pos = start, v;
            while (pos <= end)
            {
                if ((v = ppos[pos]) == 0)
                {
                    ++pos;
                    continue;
                }

                if (v < POFF)
                {
                    while (PosAt(ppos, pos, len) == v && pos <= end)
                        ++pos;
                    if (pos > end)
                        break;
                    continue;
                }

                b = v - POFF;
                break;
            }

            pos = end;
            while (pos >= start)
            {
                if ((v = ppos[pos]) == 0)
                {
                    --pos;
                    continue;
                }

                if (v > POFF)
                {
                    v -= POFF;
                }
                else
                {
                    while (PosAt(ppos, pos, len) == v)
                    {
                        --pos;
                        if (pos < start)
                        {
                            v = 0;
                            break;
                        }
                    }
                }

                e = v;
                break;
            }

            o.Append(b == 0 ? "           " : "   " + b.ToString(CultureInfo.InvariantCulture).PadLeft(7) + " ");
        }
        else
        {
            o.Append("           ");
        }

        o.Append(pro.AsSpan(start, end - start + 1)).Append(' ');
        if (isp && e != 0)
            o.Append(e.ToString(CultureInfo.InvariantCulture));
        o.Append("\n\n");
    }
}
