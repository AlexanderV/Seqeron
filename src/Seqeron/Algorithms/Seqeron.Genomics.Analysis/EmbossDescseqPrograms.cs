using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>descseq</c> (A5-8c3g, B09 F250).</summary>
/// <param name="Sequence">The value of <c>-sequence</c>: a file name with an optional <c>format::</c> prefix, or <c>asis::SEQUENCE</c>; only the first sequence is read.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers (<c>-sformat1</c>, <c>-slower1</c> ...).</param>
/// <param name="Name">The new name (<c>-name</c>); empty = the name is kept.</param>
/// <param name="Description">The new description (<c>-description</c>); empty = the description is kept.</param>
/// <param name="Append"><c>-append</c>: the name and the description are appended to the existing ones.</param>
public sealed record EmbossDescseqRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, string Name = "", string Description = "", bool Append = false);

/// <summary><c>descseq</c> of EMBOSS 6.6.0 (<c>emboss/descseq.c</c>): the name and the description of a sequence replaced or appended to, written in FASTA.</summary>
public static class EmbossDescseqPrograms
{
    /// <summary>Runs <c>descseq</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossSeqretOutput Descseq(EmbossDescseqRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "descseq" };
        var seqall = EmbossSeqall.Open("descseq", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, false, context);
        if (seqall is null)
            return EmbossSeqretPrograms.Died(context);
        var seq = seqall.Next()!;
        if (r.Append)
        {
            if (r.Name.Length > 0)
                seq.Name += r.Name;
            if (r.Description.Length > 0)
                seq.Desc += r.Description;
        }
        else
        {
            if (r.Name.Length > 0)
                seq.Name = r.Name;
            if (r.Description.Length > 0)
                seq.Desc = r.Description;
        }

        var sb = new StringBuilder();
        int unnamed = 0;
        EmbossSeqretPrograms.WriteFasta(sb, seq, ref unnamed);
        return new EmbossSeqretOutput(sb.ToString(), "", 1, [.. context.Messages]);
    }
}
