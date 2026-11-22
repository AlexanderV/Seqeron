using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// The parts of ajfeat.c / ajfeatwrite.c that prepare a table for writing and print its internal state (ajFeattableNewFtable, ajFeattableTrim, ajFeattableReverse, ajFeattablePrint = the "debug" feature format). A5-8c3e, B09 F241.
internal sealed partial class EmbossSeqFeatureTable
{
    /// <summary>The table is to be reversed (<c>thys-&gt;Rev</c>: the <c>-freverse</c> of the input).</summary>
    public bool Rev;

    /// <summary>The table was trimmed (<c>thys-&gt;Trimmed</c>).</summary>
    public bool Trimmed;

    /// <summary><c>ajFeattableNewFtable</c>: a copy; the features are copied with <see cref="NewFeat"/> and added with <c>ajFeattableAdd</c>, which raises the length to the largest position.</summary>
    public EmbossSeqFeatureTable NewFtable()
    {
        var ret = new EmbossSeqFeatureTable(_context)
        {
            SeqId = SeqId,
            Type = Type,
            TypeNull = false,
            FormatName = FormatName ?? "",
            Format = Format,
            Groups = Groups,
            Start = Start,
            End = End,
            Length = Length,
            Offset = Offset,
            Rev = Rev,
            Reversed = Reversed,
            Trimmed = Trimmed,
            Circular = Circular,
        };
        foreach (var f in Features)
        {
            var copy = NewFeat(f);
            if (!copy.Flags.HasFlag(EmbossLocationBits.RemoteId) && !copy.Flags.HasFlag(EmbossLocationBits.Label))
            {
                ret.Length = Math.Max(ret.Length, copy.Start);
                ret.Length = Math.Max(ret.Length, copy.End);
            }

            ret.Features.Add(copy);
        }

        return ret;
    }

    /// <summary><c>ajFeattableAdd</c>: the feature is appended, the length raised to its largest position unless it has a remote id or a label.</summary>
    public void Add(EmbossSeqFeature feature)
    {
        if (!feature.Flags.HasFlag(EmbossLocationBits.RemoteId) && !feature.Flags.HasFlag(EmbossLocationBits.Label))
        {
            Length = Math.Max(Length, feature.Start);
            Length = Math.Max(Length, feature.End);
        }

        Features.Add(feature);
    }

    /// <summary><c>ajFeattableMerge</c>: a copy (<see cref="NewFeat"/>, every tag checked again; the group numbers are kept) of every feature of the other table is added.</summary>
    public int Merge(EmbossSeqFeatureTable source)
    {
        foreach (var f in source.Features)
            Add(NewFeat(f));
        return Features.Count;
    }

    /// <summary><c>ajFeattableIsNuc</c>: a table of unknown type is both.</summary>
    public bool IsNuc => Type != "P";

    /// <summary><c>ajFeattableIsProt</c>.</summary>
    public bool IsProt => Type != "N";

    /// <summary><c>ajFeattableSetRange</c>: the part of the sequence the features are written for.</summary>
    public void SetRange(int begin, int end)
    {
        Start = unchecked((uint)PositionII(Length, 1, begin));
        End = unchecked((uint)PositionII(Length, unchecked((int)Start), end));
    }

    /// <summary><c>ajFeattableSetReverse</c>: the table is to be reversed (the flag is toggled, the table counts as not reversed yet).</summary>
    public void SetReverse()
    {
        Rev = !Rev;
        Reversed = false;
    }

    /// <summary><c>ajFeattableGetBegin</c>.</summary>
    public int GetBegin() => Start == 0 ? 1 : PositionII(Length, 1, unchecked((int)Start));

    /// <summary><c>ajFeattableGetEnd</c>.</summary>
    public int GetEnd() => End == 0 ? Length : PositionII(Length, GetBegin(), unchecked((int)End));

    /// <summary><c>ajFeattableTrim</c>: false when the table was trimmed already.</summary>
    public bool Trim()
    {
        if (Trimmed)
            return false;
        int begin = PositionII(Length, 1, unchecked((int)Start));
        int end = End != 0 ? PositionII(Length, begin, unchecked((int)End)) : Length;
        bool doBegin = begin > 1, doEnd = end < Length;
        Features.RemoveAll(f => !TrimOffRange(f, 0, begin, end, doBegin, doEnd));
        Offset = 0;
        Length = unchecked((int)(1 + End - Start));
        if (Rev)
            ReverseTable();
        Trimmed = true;
        return true;
    }

    /// <summary><c>ajFeattableReverse</c>: only a nucleotide table that is flagged for reversal.</summary>
    public void ReverseTable()
    {
        if (IsProt || !Rev)
            return;
        foreach (var f in Features)
        {
            if (f.Flags.HasFlag(EmbossLocationBits.RemoteId) || f.Flags.HasFlag(EmbossLocationBits.Label))
                continue;
            ReverseFeature(f, Length);
        }

        Rev = false;
        Reversed = !Reversed;
    }

    /// <summary><c>ajFeattablePrint</c>: the debugging trace of the full internal content (the <c>debug</c> feature format).</summary>
    public string Debug()
    {
        var sb = new StringBuilder();
        sb.Append("  Format: ").Append(Format).Append(" (").Append(FormatName ?? "<null>").Append(')');
        sb.Append("  Start: ").Append(Start).Append('\n');
        sb.Append("  End: ").Append(End).Append('\n');
        sb.Append("  Len: ").Append(unchecked((uint)Length)).Append('\n');
        sb.Append("  Offset: ").Append(Offset).Append('\n');
        sb.Append("  Groups: ").Append(Groups).Append('\n');
        int i = 0;
        foreach (var f in Features)
        {
            i++;
            sb.Append("\n  Feature ").Append(i).Append('\n');
            sb.Append("    Source: '").Append(f.Source).Append("'\n");
            sb.Append("    Type: '").Append(f.Type).Append("'\n");
            sb.Append("    Score: ").Append(CFormat.Fixed(f.Score, 6)).Append('\n');
            sb.Append("    Protein: ").Append(f.Protein ? "Yes" : "No").Append('\n');
            sb.Append("    Strand: '").Append(f.Strand).Append("'\n");
            sb.Append("    Start: ").Append(f.Start).Append('\n');
            sb.Append("    End: ").Append(f.End).Append('\n');
            sb.Append("    Start2: ").Append(f.Start2).Append('\n');
            sb.Append("    End2: ").Append(f.End2).Append('\n');
            sb.Append("    Remote: '").Append(f.Remote ?? "<null>").Append("'\n");
            sb.Append("    Label: '").Append(f.Label ?? "<null>").Append("'\n");
            sb.Append("    Frame: ").Append(f.Frame).Append('\n');
            sb.Append("    Exon: ").Append(f.Exon).Append('\n');
            sb.Append("    Group: ").Append(f.Group).Append('\n');
            sb.Append("    Flags: ").Append(((uint)f.Flags).ToString("x", System.Globalization.CultureInfo.InvariantCulture)).Append('\n');
            FlagLines(sb, f.Flags, "      ", true);
            sb.Append("    Tags: ").Append(f.Tags.Count).Append(" tags\n");
            TagLines(sb, f);
            if (f.Subfeatures is { Count: > 0 })
            {
                sb.Append("    Subfeatures: ").Append(f.Subfeatures.Count).Append('\n');
                int k = 0;
                foreach (var sub in f.Subfeatures)
                {
                    k++;
                    sb.Append("\n    Subfeature ").Append(k).Append('\n');
                    sb.Append("      Type: '").Append(sub.Type).Append("'\n");
                    sb.Append("      Start: ").Append(sub.Start).Append('\n');
                    sb.Append("      End: ").Append(sub.End).Append('\n');
                    sb.Append("      Flags: ").Append(((uint)sub.Flags).ToString("x", System.Globalization.CultureInfo.InvariantCulture)).Append('\n');
                    FlagLines(sb, sub.Flags, "        ", false);
                    TagLines(sb, sub);
                }
            }
        }

        return sb.ToString();
    }

    private static void TagLines(StringBuilder sb, EmbossSeqFeature f)
    {
        int j = 0;
        foreach (var (name, value) in f.Tags)
            sb.Append("      Tag ").Append((++j).ToString(System.Globalization.CultureInfo.InvariantCulture).PadLeft(3)).Append(' ').Append(name).Append(" : '").Append(value).Append("'\n");
    }

    private static void FlagLines(StringBuilder sb, EmbossLocationBits flags, string indent, bool main)
    {
        void Line(EmbossLocationBits bit, string name)
        {
            if (flags.HasFlag(bit))
                sb.Append(indent).Append(name).Append('\n');
        }

        Line(EmbossLocationBits.StartBeforeSeq, "START_BEFORE_SEQ");
        Line(EmbossLocationBits.EndAfterSeq, "END_AFTER_SEQ");
        if (main)
        {
            Line(EmbossLocationBits.Generated, "GENERATED");
            Line(EmbossLocationBits.BetweenSeq, "BETWEEN_SEQ");
            Line(EmbossLocationBits.StartTwo, "START_TWO");
            Line(EmbossLocationBits.EndTwo, "END_TWO");
        }

        Line(EmbossLocationBits.Point, "POINT");
        Line(EmbossLocationBits.ComplementMain, "COMPLEMENT_MAIN");
        Line(EmbossLocationBits.Multiple, "MULTIPLE");
        Line(EmbossLocationBits.Order, "ORDER");
        if (main)
        {
            Line(EmbossLocationBits.RemoteId, "REMOTEID");
            Line(EmbossLocationBits.Label, "LABEL");
        }

        Line(EmbossLocationBits.StartUnsure, "START_UNSURE");
        Line(EmbossLocationBits.EndUnsure, "END_UNSURE");
    }
}
