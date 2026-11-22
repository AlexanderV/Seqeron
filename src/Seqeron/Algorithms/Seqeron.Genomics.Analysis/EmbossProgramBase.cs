using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>What the sequence programs of this port share (B09 WP61, F280): opening the sequence input and the dying of the ACD.</summary>
internal static class EmbossProgramBase
{
    /// <summary>The run state of a program.</summary>
    public static EmbossFeatureReadContext Context(string program) => new() { Program = program };

    /// <summary><c>acdSetSeqall</c> / <c>acdSetSeq</c> of the input parameter <paramref name="parameter"/> (the first sequence is read here, as the ACD does when it processes the item); null: the program dies.</summary>
    public static EmbossSeqall? Open(string program, EmbossAcdArgs a, IReadOnlyDictionary<string, string> files, EmbossFeatureReadContext context, string type, string parameter = "sequence", bool single = false, bool features = false)
    {
        context.InputType ??= type.ToLowerInvariant();                         // acdInTypeSeqSave: the first input wins
        return EmbossSeqall.Open(program, a.Raw(parameter) ?? "", files, a.Sequence(parameter), features, context, type, single);
    }

    /// <summary>The program dies on the input <paramref name="parameter"/>.</summary>
    public static EmbossRunResult Died(EmbossFeatureReadContext context, string parameter = "sequence") => EmbossRunResult.From(EmbossSeqretPrograms.Died(context, parameter));

    /// <summary>The program dies on a bad value of an item (the messages of the library, then the <c>Died:</c> line).</summary>
    public static EmbossRunResult Died(EmbossFeatureReadContext context, EmbossAcdDieException e)
    {
        context.Messages.AddRange(e.Messages);
        return EmbossRunResult.From(new EmbossSeqretOutput("", "", 0, [.. context.Messages]) { Died = true });
    }

    /// <summary>The program stops at once with a message of <c>ajFatal</c> (no <c>Died:</c> line).</summary>
    public static EmbossRunResult Fatal(EmbossFeatureReadContext context, string message)
    {
        context.Messages.Add(message);
        return EmbossRunResult.From(new EmbossSeqretOutput("", "", 0, [.. context.Messages]) { Died = true });
    }

    /// <summary><c>ajStrAppendSubS</c> of the range [<paramref name="first"/>, <paramref name="last"/>] (0-based, negative from the end, clipped as <c>ajCvtSposToPosStart</c>).</summary>
    public static string SubS(string s, long first, long last)
    {
        int ibegin = Alignment.EmbossInputSequence.PositionStart(s.Length, 0, first);
        int iend = Alignment.EmbossInputSequence.PositionStart(s.Length, ibegin, last);
        return s.Length == 0 ? "" : s.Substring(ibegin, iend - ibegin + 1);
    }

    /// <summary>The nucleic or protein sequence object a program makes (<c>ajSeqNew</c>, <c>ajSeqSetNuc</c> / <c>ajSeqSetProt</c>, name, description, residues).</summary>
    public static EmbossSeq NewSeq(string name, string desc, string residues, string type)
    {
        var s = new EmbossSeq { Name = name, Desc = desc, Type = type };
        s.Seq.Append(residues);
        return s;
    }

    /// <summary><c>ajSeqType</c>: the type of a sequence that has none, from its residues (nucleic first, else protein); <c>ajSeqSetNuc</c> / <c>ajSeqSetProt</c> change the residues (x becomes n).</summary>
    public static void SeqType(EmbossSeq seq)
    {
        if (seq.Type.Length > 0)
            return;
        if (EmbossSeqType.IsNuc("", seq.SeqText))
            seq.SetNuc();
        else if (EmbossSeqType.IsProt("", seq.SeqText))
            seq.SetProt();
    }

    /// <summary><c>ajSeqAssignSeqS</c>: new residues; the range and the offsets are reset.</summary>
    public static void AssignSeq(EmbossSeq seq, string residues)
    {
        seq.Seq.Clear().Append(residues);
        seq.Begin = seq.End = seq.Offset = seq.Offend = 0;
        seq.Rev = seq.Reversed = false;
    }
}
