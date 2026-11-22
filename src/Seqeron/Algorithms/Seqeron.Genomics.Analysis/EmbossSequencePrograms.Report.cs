#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125
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
}
