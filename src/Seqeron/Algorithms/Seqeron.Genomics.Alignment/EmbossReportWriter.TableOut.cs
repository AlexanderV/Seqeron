#pragma warning disable S3358, S1121
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// The feature output formats of ajfeatwrite.c for a FULL feature table (EmbossSeqFeatureTable: sub-features, flags, remote ids, GFF3 tags): ajFeattableWrite and the writers of gff3, gff2, embl, genbank (ddbj), debug. A5-8c3e, B09 F241.
public static partial class EmbossReportWriter
{
    /// <summary>The output formats <see cref="WriteFeatureTable(EmbossSeqFeatureTable, string, DateTime, string, EmbossFeatureReadContext, EmbossFeatureWriteState)"/> writes.</summary>
    public static IReadOnlyList<string> FeatureTableFormats { get; } = ["unknown", "gff", "gff2", "gff3", "bed", "embl", "em", "genbank", "gb", "ddbj", "refseq", "refseqp", "pir", "nbrf", "swiss", "sw", "swissprot", "uniprot", "dasgff", "draw", "debug"];

    /// <summary>
    /// <c>ajFeattableWrite</c> of one table in an output feature format: a table that is to be reversed or trimmed is copied and trimmed first (<c>ajFeattableNewFtable</c> + <c>ajFeattableTrim</c>); a format the table type does not allow gives null with the library's
    /// error in the messages of <paramref name="context"/>.
    /// </summary>
    internal static string? WriteFeatureTable(EmbossSeqFeatureTable table, string format, DateTime date, string program, EmbossFeatureReadContext context, EmbossFeatureWriteState state, int tablesBefore = 0)
    {
        ArgumentNullException.ThrowIfNull(table);
        int index = FeatureFormatIndex(format);
        if (index < 0)
            throw new NotSupportedException($"The feature output format '{format}' is not an output format of the feature library.");
        var (_, writer, nucOk, protOk) = FeatureFormatTable[index];
        bool isNuc = table.IsNuc, isProt = table.IsProt;
        if (!((isNuc && nucOk) || (isProt && protOk)))
        {
            context.Error(isNuc ? $"Feature format '{format}' not supported for nucleotide sequence features" : isProt ? $"Feature format '{format}' not supported for protein sequence features" : $"Feature format '{format}' failed: unknown sequence feature type");
            return null;
        }

        var t = table;
        if (table.Rev || table.GetBegin() > 1 || table.GetEnd() < table.Length)
        {
            t = table.NewFtable();
            t.Trim();
        }

        var sb = new StringBuilder();
        switch (writer)
        {
            case "gff3":
                TableGff3(sb, t, date, context);
                break;
            case "gff2":
                TableGff2(sb, t, date, context, state);
                break;
            case "embl":
                TableEmbl(sb, t, true, context, state);
                break;
            case "genbank":
                TableEmbl(sb, t, false, context, state);
                break;
            case "bed":
                TableBed(sb, t);
                break;
            case "refseqp":
                TableRefseqp(sb, t, context, state);
                break;
            case "pir":
                TablePir(sb, t, context, state);
                break;
            case "swiss":
                TableSwiss(sb, t, context, state);
                break;
            case "dasgff":
                TableDasgff(sb, t, tablesBefore, context, state);
                break;
            case "draw":
                TableDraw(sb, t);
                break;
            default:
                sb.Append(t.Debug());
                break;
        }

        return sb.ToString();
    }

    private static char FrameChar(int frame) => frame is < 0 or > 3 ? '.' : ".012"[frame];

    private static char FrameCharNuc(int frame) => frame is <= 0 or > 3 ? '0' : ".012"[frame];

    private static string ScoreG4(float score) => score != 0f ? CFormat.General(score, 4) : ".";

    // ajStrFmtPercentEncodeC: every character of the set as %XX (upper-case hexadecimal)
    private static string PercentEncode(string s, string chars)
    {
        var sb = new StringBuilder();
        foreach (char c in s)
        {
            if (chars.Contains(c, StringComparison.Ordinal))
                sb.Append('%').Append(((int)c).ToString("X2", CultureInfo.InvariantCulture));
            else
                sb.Append(c);
        }

        return sb.ToString();
    }

    private static string RemoveWhiteAll(string s) => string.Concat(s.Where(c => !(c is ' ' or '\t' or '\n' or '\r' or '\f' or '\v')));

    // featTagAllLimit scans the list of allowed values for COMMAS, but the vocabulary stores them separated by semicolons: it never finds a value and never changes one
    private static void TagAllLimit(ref string value, string limits)
    {
        _ = value;
        _ = limits;
    }

    // ---------------------------------------------------------------------------------------------------------------- gff3

    private static void TableGff3(StringBuilder sb, EmbossSeqFeatureTable t, DateTime date, EmbossFeatureReadContext context)
    {
        sb.Append("##gff-version 3\n");
        sb.Append("##sequence-region ").Append(t.SeqId).Append(" 1 ").Append(((uint)t.Length).ToString(CultureInfo.InvariantCulture)).Append('\n');
        sb.Append("#!Date ").Append(EmbossDates.Gff(date)).Append('\n');
        sb.Append(t.Type == "P" ? "#!Type Protein\n" : "#!Type DNA\n");
        sb.Append("#!Source-version EMBOSS 6.6.0.0\n");
        foreach (var gf in t.Features)
        {
            if (t.Circular && gf.Start == 1 && gf.End == t.Length)
                gf.GffTags.Add(("Is_circular", "true"));
            DumpGff3(sb, gf, null, t, context);
        }
    }

    private static void Gff3Tag(StringBuilder sb, string tag, string value, ref bool firsttag, EmbossFeatureReadContext context)
    {
        if (tag.Length == 0)
        {
            context.Warn($"empty tag name for tag value '{value}'");
            tag = "Note";
        }

        sb.Append(firsttag ? "" : ";").Append(tag).Append('=').Append(value);
        firsttag = false;
    }

    // ajFeatGetGfftags: the GFF3 special tags of a feature, the values percent-encoded; the tags that can have several values are joined by commas
    private static Dictionary<string, string> Gfftags(EmbossSeqFeature f)
    {
        var ret = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var (tag, value) in f.GffTags)
        {
            bool multi = tag is "Alias" or "Parent" or "Note" or "Dbxref" or "Ontology_term";
            if (tag is not ("ID" or "Name" or "Alias" or "Parent" or "Target" or "Gap" or "Derives_from" or "Note" or "Dbxref" or "Ontology_term" or "Is_circular"))
                continue;
            string enc = PercentEncode(value, ";=&%,");
            if (multi && ret.TryGetValue(tag, out string? old) && old.Length > 0)
                ret[tag] = old + "," + enc;
            else
                ret[tag] = enc;
        }

        return ret;
    }

    private static void DumpGff3(StringBuilder sb, EmbossSeqFeature thys, EmbossSeqFeature? parent, EmbossSeqFeatureTable owner, EmbossFeatureReadContext context)
    {
        if (thys.Flags.HasFlag(EmbossLocationBits.Generated) && thys.Subfeatures is { Count: > 0 })
        {
            foreach (var sub in thys.Subfeatures)
                DumpGff3(sb, sub, thys, owner, context);
            return;
        }

        var vocab = EmbossFeatureVocab.ForFormat("gff3", thys.Protein)!;
        string outtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, thys.Type);
        char frame = thys.Type == "SO:0000316" ? FrameCharNuc(thys.Frame) : FrameChar(thys.Frame);
        var flagdata = new StringBuilder();
        void Flag(string text)
        {
            if (flagdata.Length > 0)
                flagdata.Append(',');
            flagdata.Append(text);
        }

        if (thys.Subfeatures is { Count: > 0 })
        {
            string subtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, thys.Subfeatures[0].Type);
            if (outtyp == subtyp)
            {
                Flag("type:" + outtyp);
                outtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, "biological_region");
            }
        }

        if (thys.Flags.HasFlag(EmbossLocationBits.StartBeforeSeq)) Flag("start_before:true");
        if (thys.Flags.HasFlag(EmbossLocationBits.EndAfterSeq)) Flag("end_after:true");
        if (thys.Flags.HasFlag(EmbossLocationBits.BetweenSeq)) Flag("between:true");
        if (thys.Start2 != 0) Flag("start2:" + thys.Start2.ToString(CultureInfo.InvariantCulture));
        if (thys.End2 != 0) Flag("end2:" + thys.End2.ToString(CultureInfo.InvariantCulture));
        if (!string.IsNullOrEmpty(thys.Label)) Flag("label:" + thys.Label);
        if (thys.Flags.HasFlag(EmbossLocationBits.StartUnsure)) Flag("start_unsure:true");
        if (thys.Flags.HasFlag(EmbossLocationBits.EndUnsure)) Flag("end_unsure:true");
        char strand = thys.Protein ? '.' : thys.Strand == '-' ? '-' : '+';
        uint tmpbeg = unchecked((uint)thys.Start), tmpend = unchecked((uint)thys.End);
        if (strand == '-' && tmpend < tmpbeg)
            (tmpbeg, tmpend) = (tmpend, tmpbeg);
        string tmpid = !string.IsNullOrEmpty(thys.Remote) ? thys.Remote : owner.SeqId;
        sb.Append(tmpid).Append('\t').Append(thys.Source).Append('\t').Append(outtyp).Append('\t').Append(unchecked((int)tmpbeg).ToString(CultureInfo.InvariantCulture)).Append('\t').Append(unchecked((int)tmpend).ToString(CultureInfo.InvariantCulture)).Append('\t');
        sb.Append(ScoreG4(thys.Score)).Append('\t').Append(strand).Append('\t').Append(frame).Append('\t');
        var gt = Gfftags(thys);
        string G(string key) => gt.TryGetValue(key, out string? v) ? v : "";
        bool firsttag = true;
        if (G("ID").Length > 0)
        {
            Gff3Tag(sb, "ID", G("ID"), ref firsttag, context);
        }
        else if (parent is null && G("Name").Length == 0)
        {
            sb.Append(firsttag ? "" : ";").Append("ID=").Append(owner.SeqId).Append('.').Append(thys.Group.ToString(CultureInfo.InvariantCulture));
            firsttag = false;
        }

        if (G("Parent").Length > 0)
        {
            Gff3Tag(sb, "Parent", G("Parent"), ref firsttag, context);
        }
        else if (parent is not null && !parent.Flags.HasFlag(EmbossLocationBits.Generated))
        {
            sb.Append(firsttag ? "" : ";").Append("Parent=").Append(owner.SeqId).Append('.').Append(parent.Group.ToString(CultureInfo.InvariantCulture));
            firsttag = false;
        }

        foreach (string key in new[] { "Name", "Alias", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular" })
        {
            if (G(key).Length > 0)
                Gff3Tag(sb, key, G(key), ref firsttag, context);
        }

        if (flagdata.Length > 0)
            Gff3Tag(sb, "featflags", flagdata.ToString(), ref firsttag, context);
        for (int i = 0; i < thys.Tags.Count; i++)
        {
            var (tag, value) = thys.Tags[i];
            string outtag = tag.ToLowerInvariant();                       // ajTagvalFmtLower: the tag of the feature itself is changed
            thys.Tags[i] = (outtag, value);
            string fmt = EmbossFeatureVocab.TagFormat(vocab.Tags, outtag);
            string val = value;
            string kind = fmt.Length >= 2 ? fmt[..2] : fmt;
            bool skip = false;
            switch (kind)
            {
                case "LI" or "QL":
                    TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(vocab.Tags, outtag));
                    break;
                case "TE" or "QW":
                    val = RemoveWhiteAll(val);
                    break;
                case "QT" or "":
                    break;
                case "QS" or "SP":
                    EmbossSeqFeatureTable.SpecialValue(outtag, ref val, context, "GFF3");     // ajFeattagSpecialGff3: the answer is not used, the value may be corrected
                    break;
                case "VO":
                    val = "true";
                    break;
                default:
                    context.Warn($"Unknown GFF3 feature tag type '{fmt}' for '{outtag}'");
                    skip = true;
                    break;
            }

            if (skip)
                continue;
            Gff3Tag(sb, outtag, PercentEncode(val, ";=%&,\t"), ref firsttag, context);
        }

        sb.Append('\n');
        if (thys.Subfeatures is { Count: > 0 })
        {
            foreach (var sub in thys.Subfeatures)
                DumpGff3(sb, sub, thys, owner, context);
        }
    }

    // ---------------------------------------------------------------------------------------------------------------- gff2

    private static void TableGff2(StringBuilder sb, EmbossSeqFeatureTable t, DateTime date, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        sb.Append("##gff-version 2.0\n");
        sb.Append("##date ").Append(EmbossDates.Gff(date)).Append('\n');
        if (t.Type == "N")
            sb.Append("##Type DNA ").Append(t.SeqId).Append('\n');
        else if (t.Type == "P")
            sb.Append("##Type Protein ").Append(t.SeqId).Append('\n');
        else
            sb.Append("##Type unknown <").Append(t.TypeNull ? "<null>" : t.Type).Append("> ").Append(t.SeqId).Append('\n');
        foreach (var gf in t.Features)
            DumpGff2(sb, gf, t, context, state);
    }

    private static void DumpGff2(StringBuilder sb, EmbossSeqFeature thys, EmbossSeqFeatureTable owner, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        var vocab = EmbossFeatureVocab.ForFormat("gff2", thys.Protein)!;
        string outtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, thys.Type);
        char strand = thys.Protein ? '.' : thys.Strand == '-' ? '-' : '+';
        uint a = unchecked((uint)thys.Start), b = unchecked((uint)thys.End);
        if (strand == '-' && b < a)
            (a, b) = (b, a);
        sb.Append(owner.SeqId).Append('\t').Append(thys.Source).Append('\t').Append(outtyp).Append('\t').Append(unchecked((int)a).ToString(CultureInfo.InvariantCulture)).Append('\t').Append(unchecked((int)b).ToString(CultureInfo.InvariantCulture)).Append('\t');
        sb.Append(ScoreG4(thys.Score)).Append('\t').Append(strand).Append('\t').Append(FrameChar(thys.Frame)).Append('\t');
        var flagdata = new StringBuilder();
        void Flag(string text)
        {
            if (flagdata.Length > 0)
                flagdata.Append(' ');
            flagdata.Append(text);
        }

        if (thys.Flags != 0) flagdata.Append("0x").Append(((uint)thys.Flags).ToString("x", CultureInfo.InvariantCulture));
        if (thys.Start2 != 0) Flag("start2:" + thys.Start2.ToString(CultureInfo.InvariantCulture));
        if (thys.End2 != 0) Flag("end2:" + thys.End2.ToString(CultureInfo.InvariantCulture));
        if (!string.IsNullOrEmpty(thys.Remote)) Flag("remoteid:" + thys.Remote);
        if (!string.IsNullOrEmpty(thys.Label)) Flag("label:" + thys.Label);
        sb.Append("Sequence \"").Append(owner.SeqId).Append('.').Append(thys.Group.ToString(CultureInfo.InvariantCulture)).Append('"');
        int i = 1;
        if (flagdata.Length > 0)
        {
            i++;
            sb.Append(" ; ").Append("FeatFlags \"").Append(flagdata).Append('"');
        }

        foreach (var (tag, value) in thys.Tags)
        {
            string outtag = EmbossFeatureVocab.TagNameOf(vocab.Tags, tag, out _);
            string fmt = EmbossFeatureVocab.TagFormat(vocab.Tags, outtag);
            if (i++ > 0)
                sb.Append(" ; ");
            string str = state.Leftover + outtag;                         // featoutStr is appended to: what an earlier writer left in it comes first
            string val = value;
            string kind = fmt.Length >= 2 ? fmt[..2] : fmt;
            switch (kind)
            {
                case "LI":
                    TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(vocab.Tags, outtag));
                    str += " " + val;
                    break;
                case "QL":
                    TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(vocab.Tags, outtag));
                    str += " " + QuoteGff2(val);
                    break;
                case "TE":
                    str += " " + RemoveWhiteAll(val);
                    break;
                case "QT":
                    str += " " + QuoteGff2(val);
                    break;
                case "QW":
                    str += " " + RemoveWhiteAll(QuoteGff2(val));
                    break;
                case "QS":
                    if (!EmbossSeqFeatureTable.SpecialValue(outtag, ref val, context, "GFF2"))
                        str = "note \"" + outtag + ": " + QuoteGff2(val) + "\"";
                    else
                        str += " " + QuoteGff2(val);
                    break;
                case "SP":
                    if (!EmbossSeqFeatureTable.SpecialValue(outtag, ref val, context, "GFF2"))
                        str = "note \"" + outtag + ": " + QuoteGff2(val) + "\"";
                    else
                        str += " " + val;
                    break;
                case "VO":
                    break;
                default:
                    context.Warn($"Unknown GFF 2.0 feature tag type '{fmt}' for '{outtag}'");
                    break;
            }

            sb.Append(str);
            state.Leftover = "";
        }

        sb.Append('\n');
    }

    // ---------------------------------------------------------------------------------------------------------------- embl / genbank

    private static void TableEmbl(StringBuilder sb, EmbossSeqFeatureTable t, bool isEmbl, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        sb.Append(isEmbl ? "FH   Key             Location/Qualifiers\nFH\n" : "FEATURES             Location/Qualifiers\n");
        var sorted = SortByGroup(t);
        foreach (var gf in sorted)
            WriteEmbl(sb, gf, t, isEmbl, context, state);
    }

    // ajListSort(thys->Features, &ajFeatCompByGroup): by group, then by exon (the signed difference of the unsigned numbers), the order of equal features kept
    private static List<EmbossSeqFeature> SortByGroup(EmbossSeqFeatureTable t)
    {
        var sorted = t.Features.OrderBy(f => f, Comparer<EmbossSeqFeature>.Create((a, b) =>
        {
            int val = unchecked((int)(a.Group - b.Group));
            return val != 0 ? val : unchecked(a.Exon - b.Exon);
        })).ToList();
        t.Features.Clear();
        t.Features.AddRange(sorted);
        return sorted;
    }

    private static string EmblPos(EmbossSeqFeature f)
    {
        var pos = new StringBuilder();
        if (f.Flags.HasFlag(EmbossLocationBits.RemoteId))
            pos.Append(f.Remote ?? "<null>").Append(':');
        if (f.Flags.HasFlag(EmbossLocationBits.Label))
            pos.Append(f.Label ?? "<null>");
        else if (f.Flags.HasFlag(EmbossLocationBits.StartBeforeSeq))
            pos.Append('<').Append(f.Start.ToString(CultureInfo.InvariantCulture));
        else if (f.Flags.HasFlag(EmbossLocationBits.StartTwo))
            pos.Append('(').Append(f.Start.ToString(CultureInfo.InvariantCulture)).Append('.').Append(f.Start2.ToString(CultureInfo.InvariantCulture)).Append(')');
        else
            pos.Append(f.Start.ToString(CultureInfo.InvariantCulture));
        if (!f.Flags.HasFlag(EmbossLocationBits.Point))
        {
            if (f.Flags.HasFlag(EmbossLocationBits.BetweenSeq))
                pos.Append('^').Append(f.End.ToString(CultureInfo.InvariantCulture));
            else if (f.Flags.HasFlag(EmbossLocationBits.EndAfterSeq))
                pos.Append("..>").Append(f.End.ToString(CultureInfo.InvariantCulture));
            else if (f.Flags.HasFlag(EmbossLocationBits.EndTwo))
                pos.Append("..(").Append(f.End2.ToString(CultureInfo.InvariantCulture)).Append('.').Append(f.End.ToString(CultureInfo.InvariantCulture)).Append(')');
            else
                pos.Append("..").Append(f.End.ToString(CultureInfo.InvariantCulture));
        }

        return pos.ToString();
    }

    // featWriteEmbl
    private static void WriteEmbl(StringBuilder sb, EmbossSeqFeature thys, EmbossSeqFeatureTable fttable, bool isEmbl, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        var location = new StringBuilder();
        if (thys.Flags.HasFlag(EmbossLocationBits.Multiple) && thys.Subfeatures is { Count: > 0 })
        {
            bool whole = thys.Flags.HasFlag(EmbossLocationBits.ComplementMain);
            if (whole)
                location.Append("complement(");
            location.Append(thys.Flags.HasFlag(EmbossLocationBits.Order) ? "order(" : "join(");
            int exon = 0;
            foreach (var gf in thys.Subfeatures)
            {
                if (exon++ > 0)
                    location.Append(',');
                string pos = EmblPos(gf);
                location.Append(gf.Strand == '-' && !whole ? "complement(" + pos + ")" : pos);
            }

            location.Append(')');
            if (whole)
                location.Append(')');
            DumpEmbl(sb, thys, location.ToString(), fttable.SeqId, isEmbl, context, state);
            return;
        }

        bool whole2 = thys.Flags.HasFlag(EmbossLocationBits.ComplementMain);
        if (whole2)
            location.Append("complement(");
        string p = EmblPos(thys);
        location.Append(thys.Strand == '-' && !whole2 ? "complement(" + p + ")" : p);
        if (whole2)
            location.Append(')');
        DumpEmbl(sb, thys, location.ToString(), fttable.SeqId, isEmbl, context, state);
        if (thys.Subfeatures is { Count: > 0 })
        {
            foreach (var sub in thys.Subfeatures)
                WriteEmbl(sb, sub, fttable, isEmbl, context, state);
        }
    }

    // featTagQuoteEmbl: the value in double quotes, an inner quote doubled
    private static string QuoteEmblValue(string v)
    {
        if (!v.Contains('"', StringComparison.Ordinal))
        {
            string q = v;
            if (q.Length == 0 || q[0] != '"')
                q = "\"" + q;
            if (q.Length == 1 || q[^1] != '"')
                q += "\"";
            return q;
        }

        return "\"" + v.Replace("\"", "\"\"", StringComparison.Ordinal) + "\"";
    }

    // featDumpEmbl
    private static void DumpEmbl(StringBuilder sb, EmbossSeqFeature feat, string location, string seqid, bool isEmbl, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        var vocab = EmbossFeatureVocab.ForFormat("embl", false)!;
        string tmptyp = EmbossFeatureVocab.TypeExternal(vocab.Types, feat.Type);
        string lead = isEmbl ? "FT" : "  ";
        int margin = isEmbl ? 80 : 79;
        string preftyploc = lead + "   " + (tmptyp.Length > 15 ? tmptyp[..15] : tmptyp).PadRight(15) + " ";
        string preftyptag = lead + "                   ";
        LocWrap(location, margin, preftyptag, preftyploc, sb);
        if (vocab.Types.FetchKey(feat.Type) is null)
        {
            // a type the EMBL vocabulary does not know: the term of the ontology is named in a note
            string? term = EmbossSoTerms.Name(feat.Type);
            if (!feat.Type.Contains(':', StringComparison.Ordinal))
                context.Error($"Failed to open filename '{feat.Type}'");           // ajOboinQryS takes a name that is no id for a file
            if (term is not null)
                TagWrapState("/note=\"*Type " + feat.Type + " " + term + "\"\n", margin, preftyptag, sb, state);
        }

        foreach (var (tag, value) in feat.GffTags)
        {
            string str;
            if (tag is "Dbxref" or "Ontology_term")
            {
                string tmptag = EmbossFeatureVocab.TagNameOf(vocab.Tags, "db_xref", out _);
                str = "/" + tmptag + "=" + QuoteEmblValue(value) + "\n";
            }
            else if (tag is "ID" or "Parent")
            {
                continue;
            }
            else
            {
                str = "/note=\"*" + tag + " " + value + "\"\n";
            }

            TagWrapState(str, margin, preftyptag, sb, state);
        }

        DumpEmblTags(sb, feat, vocab.Tags, margin, preftyptag, seqid, context, state);
    }

    // the tags of a feature in the EMBL layout (featDumpEmbl, featDumpRefseq, featDumpRefseqp share the loop; only the vocabulary differs)
    private static void DumpEmblTags(StringBuilder sb, EmbossSeqFeature feat, EmbossFeatureVocab.Table tags, int margin, string preftyptag, string seqid, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        foreach (var (tag, value) in feat.Tags)
        {
            string tmptag = EmbossFeatureVocab.TagNameOf(tags, tag, out _);
            string fmt = EmbossFeatureVocab.TagFormat(tags, tmptag);
            string str = "/" + tmptag;
            string val = value;
            string kind = fmt.Length >= 2 ? fmt[..2] : fmt;
            switch (kind)
            {
                case "LI":
                    TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(tags, tmptag));
                    str += "=" + val + "\n";
                    break;
                case "QL":
                    TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(tags, tmptag));
                    str += "=" + QuoteEmblValue(val) + "\n";
                    break;
                case "QS":
                    if (!EmbossSeqFeatureTable.SpecialValue(tmptag, ref val, context))
                    {
                        context.Warn($"{seqid}: Bad special tag value");
                        str = "/note=\"" + tmptag + ": " + QuoteEmblValue(val) + "\"";
                    }
                    else
                    {
                        str += "=" + QuoteEmblValue(val) + "\n";
                    }

                    break;
                case "SP":
                    if (!EmbossSeqFeatureTable.SpecialValue(tmptag, ref val, context))
                    {
                        context.Warn($"{seqid}: Bad special tag value");
                        str = "/note=\"" + tmptag + ": " + QuoteEmblValue(val) + "\"";
                    }
                    else
                    {
                        str += "=" + val + "\n";
                    }

                    break;
                case "TE":
                    str += "=" + RemoveWhiteAll(val) + "\n";
                    break;
                case "VO":
                    break;
                case "QT":
                    str += "=" + QuoteEmblValue(val) + "\n";
                    break;
                case "QW":
                    str += "=" + RemoveWhiteAll(QuoteEmblValue(val)) + "\n";
                    break;
                default:
                    context.Warn($"Unknown EMBL feature tag type '{fmt}' for '{tmptag}'");
                    break;
            }

            TagWrapState(str, margin, preftyptag, sb, state);
        }
    }

    // featTagEmblWrapC on the static featoutStr: the wrapped text is written and the collapsed string stays behind for the next writer
    private static void TagWrapState(string str, int margin, string prefix, StringBuilder sb, EmbossFeatureWriteState state)
    {
        TagWrap(str, margin, prefix, sb);
        state.Leftover = CollapseWhite(str);
    }
}
