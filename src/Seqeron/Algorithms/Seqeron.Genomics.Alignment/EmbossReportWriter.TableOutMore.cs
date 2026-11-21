#pragma warning disable S3358, S1121
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// The feature output formats of ajfeatwrite.c for a FULL feature table that are not in EmbossReportWriter.TableOut.cs: bed (ajFeattableWriteBed, featDumpBed), refseqp (feattableWriteRefseqp, featWriteRefseqp, featDumpRefseqp), refseq
// (the format table of the programs writes "refseq" with ajFeattableWriteGenbank; feattableWriteRefseq / featDumpRefseq are reached by no format and no program), pir (ajFeattableWritePir, featDumpPir, featTypePirOut), swiss (ajFeattableWriteSwiss,
// featDumpSwiss, featTagSwissWrapC), dasgff (ajFeattableWriteDasgff) and draw (ajFeattableWriteDraw). A5-8c3f, B09 F244.
public static partial class EmbossReportWriter
{
    private static string I(int v) => v.ToString(CultureInfo.InvariantCulture);

    private static string U(uint v) => v.ToString(CultureInfo.InvariantCulture);

    // ------------------------------------------------------------------------------------------------------------------ bed

    // ajFeattableWriteBed: ajFeattableGetScorerange never reports a range (the minimum and maximum it keeps are the 0.0 it starts with), so useScore is 0 and every score is 0
    private static void TableBed(StringBuilder sb, EmbossSeqFeatureTable t)
    {
        sb.Append("track name=").Append(t.SeqId).Append(" description=\"EMBOSS:6.6.0.0\" useScore=0\n");
        foreach (var thys in t.Features)
        {
            var vocab = EmbossFeatureVocab.ForFormat("gff2", thys.Protein)!;
            string outtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, thys.Type);
            char strand = thys.Protein ? '.' : thys.Strand == '-' ? '-' : '+';
            uint istart, iend;
            if (strand == '-' && unchecked((uint)thys.End) < unchecked((uint)thys.Start))
            {
                iend = unchecked((uint)thys.Start - 1);
                istart = unchecked((uint)thys.End - 1);
            }
            else
            {
                istart = unchecked((uint)thys.Start - 1);
                iend = unchecked((uint)thys.End - 1);
            }

            sb.Append(t.SeqId).Append('\t').Append(I(unchecked((int)istart))).Append('\t').Append(I(unchecked((int)iend))).Append('\t').Append(outtyp).Append("\t0\t").Append(strand).Append('\t')
              .Append(I(unchecked((int)istart))).Append('\t').Append(I(unchecked((int)iend))).Append("\t0");
            if (thys.Subfeatures is { Count: > 0 })
            {
                int nsub = 0;
                var lengths = new StringBuilder();
                var offsets = new StringBuilder();
                string punc = "";
                foreach (var gf in thys.Subfeatures)
                {
                    if (!string.IsNullOrEmpty(gf.Remote))
                        continue;
                    nsub++;
                    int ilen = unchecked((int)((uint)gf.End - (uint)gf.Start + 1));
                    int ioff = unchecked((int)((uint)gf.Start - istart - 1));
                    lengths.Append(punc).Append(I(ilen));
                    offsets.Append(punc).Append(I(ioff));
                    punc = ",";
                }

                sb.Append('\t').Append(U((uint)nsub)).Append('\t').Append(nsub > 0 ? lengths.ToString() : "<null>").Append('\t').Append(nsub > 0 ? offsets.ToString() : "<null>").Append('\n');
            }
            else
            {
                sb.Append("\t1\t0\t").Append(I(unchecked((int)((uint)thys.End - (uint)thys.Start + 1)))).Append('\n');
            }
        }
    }

    // ------------------------------------------------------------------------------------------------------------------ draw

    // ajFeattableWriteDraw
    private static void TableDraw(StringBuilder sb, EmbossSeqFeatureTable t)
    {
        sb.Append("Start ").Append(I(t.GetBegin())).Append('\n');
        sb.Append("End   ").Append(I(t.GetEnd())).Append('\n');
        sb.Append('\n').Append("group\n");
        foreach (var f in t.Features)
        {
            var vocab = f.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna;           // ajFeatGetType
            sb.Append("label\n");
            sb.Append("Range ").Append(I(f.Start)).Append(' ').Append(I(f.End)).Append(" | | 9 H\n");
            sb.Append(EmbossFeatureVocab.TypeExternal(vocab.Types, f.Type)).Append('\n');
            sb.Append("endlabel\n\n");
        }

        sb.Append("endgroup\n");
    }

    // ------------------------------------------------------------------------------------------------------------------ refseq / refseqp

    // the position text of a protein location (featWriteRefseqp): < before the start, ^ for between, ..> after the end
    private static string RefseqpPos(EmbossSeqFeature f)
    {
        var pos = new StringBuilder();
        pos.Append(f.Flags.HasFlag(EmbossLocationBits.StartBeforeSeq) ? "<" + I(f.Start) : I(f.Start));
        if (!f.Flags.HasFlag(EmbossLocationBits.Point))
        {
            if (f.Flags.HasFlag(EmbossLocationBits.BetweenSeq))
                pos.Append('^').Append(I(f.End));
            else if (f.Flags.HasFlag(EmbossLocationBits.EndAfterSeq))
                pos.Append("..>").Append(I(f.End));
            else
                pos.Append("..").Append(I(f.End));
        }

        return pos.ToString();
    }

    // feattableWriteRefseqp
    private static void TableRefseqp(StringBuilder sb, EmbossSeqFeatureTable t, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        if (t.IsNuc)
            return;                                  // feattableWriteRefseqp returns ajFalse for a table that may be nucleotide (the type unknown is both): nothing is written, not even the heading
        sb.Append("FEATURES             Location/Qualifiers\n");
        foreach (var gf in SortByGroup(t))
            WriteRefseqp(sb, gf, t.SeqId, context, state);
    }

    // featWriteRefseqp: the features of the sub-feature list are written after their parent (the ones of a joined location are inside it)
    private static void WriteRefseqp(StringBuilder sb, EmbossSeqFeature thys, string seqid, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
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
                location.Append(RefseqpPos(gf));
            }

            location.Append(')');
            if (whole)
                location.Append(')');
            DumpRefseqp(sb, thys, location.ToString(), seqid, context, state);
            return;
        }

        location.Append(RefseqpPos(thys));
        DumpRefseqp(sb, thys, location.ToString(), seqid, context, state);
        if (thys.Subfeatures is { Count: > 0 })
        {
            foreach (var sub in thys.Subfeatures)
                WriteRefseqp(sb, sub, seqid, context, state);
        }
    }

    // featDumpRefseqp: the GenBank layout with the vocabulary of the format (no ontology note, no GFF tags)
    private static void DumpRefseqp(StringBuilder sb, EmbossSeqFeature feat, string location, string seqid, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        var vocab = EmbossFeatureVocab.ForFormat("refseqp", true)!;
        string tmptyp = EmbossFeatureVocab.TypeExternal(vocab.Types, feat.Type);
        string preftyploc = "     " + (tmptyp.Length > 15 ? tmptyp[..15] : tmptyp).PadRight(15) + " ";           // "%s   %-15.15S " with "  "
        string preftyptag = "                     ";
        LocWrap(location, 79, preftyptag, preftyploc, sb);
        DumpEmblTags(sb, feat, vocab.Tags, 79, preftyptag, seqid, context, state);
    }

    // ------------------------------------------------------------------------------------------------------------------ pir

    // ajStrFmtTitle: lower case, then the first letter upper case
    private static string FmtTitle(string s)
    {
        char[] c = s.ToLowerInvariant().ToCharArray();
        for (int i = 0; i < c.Length; i++)
        {
            if (char.IsAsciiLetter(c[i]))
            {
                c[i] = char.ToUpperInvariant(c[i]);
                break;
            }
        }

        return new string(c);
    }

    // ajFeattableWritePir
    private static void TablePir(StringBuilder sb, EmbossSeqFeatureTable t, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        foreach (var gf in SortByGroup(t))
        {
            var location = new StringBuilder();
            if (gf.Subfeatures is { Count: > 0 })
            {
                foreach (var gfsub in gf.Subfeatures)
                {
                    if (location.Length > 0)
                        location.Append(',');
                    location.Append(I(gfsub.Start));
                    if (gfsub.End != gfsub.Start)
                        location.Append('-').Append(I(gfsub.End));
                }
            }
            else
            {
                location.Append(I(gf.Start));
                if (gf.End != gf.Start)
                    location.Append('-').Append(I(gf.End));
            }

            DumpPir(sb, t, gf, location.ToString(), state);
        }
    }

    // featDumpPir with featTypePirOut
    private static void DumpPir(StringBuilder sb, EmbossSeqFeatureTable owner, EmbossSeqFeature thys, string location, EmbossFeatureWriteState state)
    {
        var tags = EmbossFeatureVocab.ForFormat("pir", true)!;
        var copy = owner.NewFeat(thys);
        var outcomm = new StringBuilder();
        string pirtype = EmbossFeatureVocab.TypeExternal(tags.Types, copy.Type);
        string pirintype = EmbossFeatureVocab.TypeInternal(tags.Types, pirtype) ?? "";
        bool typmod = !string.Equals(copy.Type, pirintype, StringComparison.OrdinalIgnoreCase);
        if (typmod)
        {
            // ajFeatSetDescApp(copy, ajFeatTypeProt(copy->Type)): ", " and the protein name of the type are appended to the first note, or a note is made
            string desc = EmbossFeatureVocab.TypeExternal(EmbossFeatureVocab.Protein.Types, copy.Type);
            int note = copy.Tags.FindIndex(g => string.Equals(g.Name, "note", StringComparison.OrdinalIgnoreCase));
            if (note >= 0)
                copy.Tags[note] = (copy.Tags[note].Name, copy.Tags[note].Value + ", " + desc);
            else
                copy.Tags.Add((EmbossFeatureVocab.TagNameOf(copy.Protein ? EmbossFeatureVocab.Protein.Tags : EmbossFeatureVocab.Dna.Tags, "note", out _), desc));
        }

        string outtyp = FmtTitle(pirtype).Replace('_', ' ');
        sb.Append("F;").Append(location).Append('/').Append(outtyp).Append(':');
        foreach (var (name, value) in copy.Tags)
        {
            string outtag = EmbossFeatureVocab.TagNameOf(tags.Tags, name, out _);
            if (outtag.Equals("comment", StringComparison.OrdinalIgnoreCase))
            {
                outcomm.Append(" #").Append(value);
                continue;
            }

            state.Leftover += " " + value;               // every value type of the PIR vocabulary is written the same way, after the string the last writer left behind
            sb.Append(state.Leftover);
            state.Leftover = "";
        }

        sb.Append(outcomm).Append('\n');
    }

    // ------------------------------------------------------------------------------------------------------------------ swiss

    // featTagSwissWrapC: wraps at the last space that fits, breaks before a " /FTId=", and a word that fills the line makes the library loop for ever (ajStrFindlastC returns -1 for an unsigned position)
    private static void SwissWrap(string pval, int margin, string prefix, StringBuilder sb)
    {
        int left = prefix.Length;
        int width = margin - left;
        if (pval.Length <= left)
        {
            sb.Append(pval).Append('\n');
            return;
        }

        var ret = new StringBuilder(EmbossSeqFeatureTable.Sub(pval, 0, left - 1));
        string valstr = EmbossSeqFeatureTable.Sub(pval, left, -1);
        int len = valstr.Length;
        bool isftidstart = false;
        if (valstr.StartsWith("/FTId=", StringComparison.Ordinal))
        {
            string trimmed = ret.ToString().TrimEnd(' ');
            ret.Clear().Append(trimmed).Append('\n');
            isftidstart = true;
        }

        int k;
        for (int i = 0; i < len; i += k)
        {
            int last = i + width - 1;
            string tmp = EmbossSeqFeatureTable.Sub(valstr, i, len - 1);
            int iftid = tmp.IndexOf(" /FTId=", StringComparison.Ordinal);
            long j;
            if (iftid >= 0 && iftid + i <= last)
            {
                if (iftid > 0)
                {
                    tmp = EmbossSeqFeatureTable.Sub(valstr, i, i + iftid - 1);
                    j = iftid;
                }
                else
                {
                    tmp = EmbossSeqFeatureTable.Sub(valstr, i, i + iftid - 2);
                    j = 0;
                }
            }
            else if (last + 1 >= len)
            {
                tmp = EmbossSeqFeatureTable.Sub(valstr, i, len - 1);
                j = 0;
            }
            else if (valstr[last + 1] == ' ')
            {
                tmp = EmbossSeqFeatureTable.Sub(valstr, i, last);
                j = width;
            }
            else
            {
                tmp = EmbossSeqFeatureTable.Sub(valstr, i, last);
                j = tmp.LastIndexOf(' ');
                if (j < 0)
                    throw new InvalidOperationException("The reference program never stops on this input: featTagSwissWrapC finds no space in a line (ajStrFindlastC returns -1 into an unsigned position) and repeats the line for ever.");
            }

            if (j < 1)
            {
                j = tmp.Length;
                k = (int)j;
            }
            else
            {
                k = (int)j + 1;
            }

            if (i > 0 || isftidstart)
                ret.Append(prefix);
            ret.Append(tmp, 0, (int)Math.Min(j, tmp.Length)).Append('\n');
        }

        sb.Append(ret);
    }

    // ajFeattableWriteSwiss
    private static void TableSwiss(StringBuilder sb, EmbossSeqFeatureTable t, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        foreach (var gf in t.Features)
        {
            if (gf.Subfeatures is { Count: > 0 })
            {
                foreach (var gfsub in gf.Subfeatures)
                    DumpSwiss(sb, gfsub, gf, context, state);
            }
            else
            {
                DumpSwiss(sb, gf, gf, context, state);
            }
        }
    }

    // featDumpSwiss: the tags are those of the parent
    private static void DumpSwiss(StringBuilder sb, EmbossSeqFeature thys, EmbossSeqFeature gftop, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        var vocab = EmbossFeatureVocab.ForFormat("swiss", true)!;
        string outtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, thys.Type);
        string fromstr = thys.Flags.HasFlag(EmbossLocationBits.StartUnsure) ? (thys.Start != 0 ? "?" + I(thys.Start) : "?")
            : thys.Flags.HasFlag(EmbossLocationBits.StartBeforeSeq) ? "<" + I(thys.Start) : I(thys.Start);
        string tostr = thys.Flags.HasFlag(EmbossLocationBits.EndUnsure) ? (thys.End != 0 ? "?" + I(thys.End) : "?")
            : thys.Flags.HasFlag(EmbossLocationBits.EndAfterSeq) ? ">" + I(thys.End) : I(thys.End);
        var str = new StringBuilder("FT   ");
        str.Append((outtyp.Length > 8 ? outtyp[..8] : outtyp).PadRight(8)).Append(' ');
        str.Append((fromstr.Length > 6 ? fromstr[..6] : fromstr).PadLeft(6)).Append(' ').Append((tostr.Length > 6 ? tostr[..6] : tostr).PadLeft(6));
        int i = 0;
        bool wasnote = false;
        foreach (var (name, value) in gftop.Tags)
        {
            string outtag = EmbossFeatureVocab.TagNameOf(vocab.Tags, name, out _);
            string fmt = EmbossFeatureVocab.TagFormat(vocab.Tags, outtag);
            str.Append(i++ > 0 ? " " : "       ");
            string kind = fmt.Length >= 2 ? fmt[..2] : fmt;
            string val = value;
            switch (kind)
            {
                case "LI" or "QL":
                    TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(vocab.Tags, outtag));
                    str.Append(val).Append('.');
                    wasnote = false;
                    break;
                case "TA":
                    str.Append(outtag.Equals("ftid", StringComparison.OrdinalIgnoreCase) ? "/FTId=" + val + "." : "/" + outtag + "=" + val);
                    wasnote = false;
                    break;
                case "TE":
                    str.Append(val).Append('.');
                    wasnote = true;
                    break;
                case "BT":
                    if (wasnote)
                    {
                        str.Length = Math.Max(0, str.Length - 2);
                        str.Append(' ');
                    }

                    str.Append('(').Append(val).Append(").");
                    wasnote = false;
                    break;
                default:
                    context.Warn($"Unknown SWISS feature tag type '{fmt}' for '{outtag}'");
                    wasnote = false;
                    break;
            }
        }

        SwissWrap(str.ToString(), 75, "FT                                ", sb);
        state.Leftover = "";
    }

    // ------------------------------------------------------------------------------------------------------------------ dasgff

    // featCategory of ajfeat.c, as ajFeatTypeGetCategory reads it (case-insensitive types; one list of the promoter category lacks a comma and so names no type)
    private static readonly Lazy<Dictionary<string, string>> FeatCategoryTable = new(() =>
    {
        (string Name, string Types)[] cat =
        [
            ("other", "SO:2000061,SO:0000110,SO:0000723,SO:0000298,SO:0000724,SO:0000704,SO:0000410,SO:0000409,SO:0000297,SO:0000313,SO:0005836_rep_origin"),
            ("component", ""), ("supercomponent", ""),
            ("translation", "SO:0000316,SO:0000419,SO:0000418,SO:0000725,SO:0000204,SO:0000205,SO:0000139"),
            ("transcription", "SO:0000147,SO:0000188,SO:0000673,SO:0005836,SO:0000234,SO:0000551,SO:0000553,SO:0000185_precursor_RNA,SO:0000185,SO:0000252,SO:0000013,SO:0000275,SO:0000274,SO:0000253,SO:0000557,SO:0000555"),
            ("promoter", "SO:0000167,SO:0000172,SO:0000174,SO:0000178,SO:0000140,SO:0000173" + "SO:0000175,SO:0000176,SO:0000165,SO:0000141"),
            ("variation", ""), ("structural", "SO:0005850,SO:0000002,SO:0000331,SO:0000730"), ("similarity", ""),
            ("repeat", "SO:0000657,SO:0000726,SO:0000005,SO:0000286"), ("experimental", ""),
            ("immunoglobulin", "SO:0000001_C_region,SO:0000001_J_segment,SO:0000001_N_region,SO:0000001_S_region,SO:0000001_V_region,SO:0000001_V_segment,SO:0000458"),
        ];
        var table = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase) { [""] = "other" };
        foreach (var (name, types) in cat)
        {
            foreach (string token in types.Split(',', StringSplitOptions.RemoveEmptyEntries))
                table[token] = name;                 // ajTablePut replaces the value of a key that is there
        }

        return table;
    });

    // ajFeatGetTagC: the nth tag of a name; a note "*name: value" counts as a tag of that name
    private static string? FeatGetTag(EmbossSeqFeature f, string tname, int num)
    {
        bool isnote = tname == "note";
        int inum = 0;
        foreach (var (name, tagval) in f.Tags)
        {
            if (name.Equals(tname, StringComparison.OrdinalIgnoreCase))
            {
                inum++;
                if (num == inum)
                    return tagval;
            }
            else if (!isnote && name.Equals("note", StringComparison.OrdinalIgnoreCase) && tagval.Length > 0 && tagval[0] == '*' && tagval.Length > tname.Length + 1
                     && tagval.AsSpan(1).StartsWith(tname, StringComparison.OrdinalIgnoreCase) && tagval[tname.Length + 1] == ':')
            {
                inum++;
                if (num == inum)
                    return EmbossSeqFeatureTable.Sub(tagval, tname.Length + 3, -1);
            }
        }

        return null;
    }

    // ajFeattableWriteDasgff: the XML header is written for the first table of the output only
    private static void TableDasgff(StringBuilder sb, EmbossSeqFeatureTable t, int tablesBefore, EmbossFeatureReadContext context, EmbossFeatureWriteState state)
    {
        if (tablesBefore == 0)
            sb.Append("<?xml version=\"1.0\" standalone=\"no\"?>\n<DASGFF>\n  <GFF version=\"1.0\" href=\"url\">\n");
        sb.Append("    <SEGMENT id=\"").Append(t.SeqId).Append("\" start=\"").Append(I(unchecked((int)(1 + t.Offset + t.Start)))).Append("\" stop=\"").Append(I(unchecked((int)(1 + t.Offset + t.End)))).Append("\"\n");
        sb.Append("                version=\"0.0\">\n");
        foreach (var gf in t.Features)
        {
            var vocab = EmbossFeatureVocab.ForFormat("gff3", gf.Protein)!;
            string outtyp = EmbossFeatureVocab.TypeExternal(vocab.Types, gf.Type);
            sb.Append("      <FEATURE id=\"").Append(t.SeqId).Append('.').Append(gf.Source).Append('.').Append(U(gf.Group)).Append("\">\n");
            string? category = FeatGetTag(gf, "category", 1);
            if (string.IsNullOrEmpty(category))
                category = FeatCategoryTable.Value.TryGetValue(gf.Type, out string? c) ? c : FeatCategoryTable.Value[""];
            sb.Append("         <TYPE id=\"").Append(gf.Type).Append("\" category=\"").Append(category).Append("\" reference=\"no\" subparts=\"no\">").Append(outtyp).Append("</TYPE>\n");
            if (gf.Source.Length > 0)
                sb.Append("         <METHOD id=\"").Append(gf.Source).Append("\"></METHOD>\n");
            else
                sb.Append("         <METHOD></METHOD>\n");
            sb.Append("         <START>").Append(I(gf.Start)).Append("</START>\n");
            sb.Append("         <END>").Append(I(gf.End)).Append("</END>\n");
            if (!(Math.Abs((double)gf.Score) <= (double)1.192e-6F))
                sb.Append("         <SCORE>").Append(CFormat.Fixed(gf.Score, 6)).Append("</SCORE>\n");
            if (gf.Strand != '\0')
                sb.Append("         <ORIENTATION>").Append(gf.Strand).Append("</ORIENTATION>\n");
            else
                sb.Append("         <ORIENTATION>0</ORIENTATION>\n");
            char frame = gf.Type == "SO:0000316" ? FrameCharNuc(gf.Frame) : FrameChar(gf.Frame);
            if (gf.Frame == 0)
                sb.Append("         <PHASE>-</PHASE>\n");
            else
                sb.Append("         <PHASE>").Append(frame).Append("</PHASE>\n");
            foreach (var (name, value) in gf.Tags)
            {
                string outtag = EmbossFeatureVocab.TagNameOf(vocab.Tags, name, out _);
                string fmt = EmbossFeatureVocab.TagFormat(vocab.Tags, outtag);
                sb.Append("         <NOTE>");
                string str = outtag;
                string val = value;
                string kind = fmt.Length >= 2 ? fmt[..2] : fmt;
                switch (kind)
                {
                    case "LI" or "QL":
                        TagAllLimit(ref val, EmbossFeatureVocab.TagLimit(vocab.Tags, outtag));
                        str += ":" + val;
                        break;
                    case "TE" or "QW":
                        str += ":" + RemoveWhiteAll(val);
                        break;
                    case "QT":
                        if (outtag == "note" && val.Length > 0 && val[0] == '*')
                        {
                            val = val[1..];
                            int a = 0;
                            while (a < val.Length && EmbossSeqFeatureTable.IsCSpace(val[a]))
                                a++;
                            if (a >= val.Length)
                            {
                                str = "(null):(null)";             // ajStrExtractWord finds no word: the two strings it was to fill are NULL
                            }
                            else
                            {
                                int b = a;
                                while (b < val.Length && !EmbossSeqFeatureTable.IsCSpace(val[b]))
                                    b++;
                                int r = b;
                                while (r < val.Length && EmbossSeqFeatureTable.IsCSpace(val[r]))
                                    r++;
                                str = val[a..b] + ":" + (r < val.Length ? val[r..] : "");
                            }
                        }
                        else
                        {
                            str += ":" + val;
                        }

                        break;
                    case "QS" or "SP":
                        if (!EmbossSeqFeatureTable.SpecialValue(outtag, ref val, context, "GFF3"))
                            str = outtag + ":" + PercentEncode(val, ";=%&,\t");          // featTagDasgffDefault
                        else
                            str += ":" + val;
                        break;
                    case "VO":
                        break;
                    default:
                        context.Warn($"Unknown GFF3 feature tag type '{fmt}' for '{outtag}'");
                        break;
                }

                sb.Append(str).Append("</NOTE>\n");
                state.Leftover = str;
            }

            sb.Append("      </FEATURE>\n");
        }

        sb.Append("    </SEGMENT>\n");
    }
}
