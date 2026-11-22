using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>seqret</c> (A5-8c3g, B09 F248).</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA): a file name with an optional <c>format::</c> prefix.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers (<c>-sformat1</c>, <c>-sbegin1</c> ...).</param>
/// <param name="Feature"><c>-feature</c>: the features of the sequences are read and written.</param>
/// <param name="FirstOnly"><c>-firstonly</c>: read one sequence and stop.</param>
/// <param name="OutSeq"><c>-outseq</c>: the sequence output file name, optionally with a <c>format::</c> prefix.</param>
/// <param name="OutSeqFormat"><c>-osformat2</c>: the sequence output format (default: <c>fasta</c>).</param>
/// <param name="OutFeat"><c>-ofname2</c>: the feature output file name (default: the name of the sequence output file with the extension of the feature format).</param>
/// <param name="OutFeatFormat"><c>-offormat2</c>: the feature output format (default <c>gff3</c>).</param>
/// <param name="RunDate">The date of the run (GFF headers); null = now.</param>
public sealed record EmbossSeqretRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, bool Feature = false, bool FirstOnly = false,
    string OutSeq = "stdout", string? OutSeqFormat = null, string? OutFeat = null, string? OutFeatFormat = null, DateTime? RunDate = null);

/// <summary>The output of <c>seqret</c>: the sequence file, the feature file and the messages the library printed.</summary>
/// <param name="SequenceText">The text of the sequence output file.</param>
/// <param name="FeatureText">The text of the feature output file (empty when no feature was written).</param>
/// <param name="Sequences">The number of sequences written.</param>
/// <param name="Messages">What the library printed on its error stream, in order (the ones before the program died, too).</param>
public sealed record EmbossSeqretOutput(string SequenceText, string FeatureText, int Sequences, IReadOnlyList<string> Messages)
{
    /// <summary>The name of the feature output file.</summary>
    public string FeatureFileName { get; init; } = "";

    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}

/// <summary>
/// <c>seqret</c> of EMBOSS 6.6.0 (<c>emboss/seqret.c</c>): the sequences of a file read by the sequence readers of <c>ajseqread.c</c> (<see cref="EmbossSeqin"/>) and written by the sequence and feature writers.
/// </summary>
public static class EmbossSeqretPrograms
{
    /// <summary>Runs <c>seqret</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program does not stop or dies by a fault (the library dereferences a null pointer).</exception>
    public static EmbossSeqretOutput Seqret(EmbossSeqretRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "seqret" };
        var seqall = EmbossSeqall.Open("seqret", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, r.Feature, context);
        if (seqall is null)
            return Died(context);
        string featFormat = r.OutFeatFormat ?? "gff3";
        int featIndex = EmbossReportWriter.FeatureFormatIndex(featFormat);
        string featName = EmbossReportWriter.FeatureFormatName(featIndex < 0 ? 1 : featIndex);
        var seqOut = new StringBuilder();
        var featOut = new StringBuilder();
        var state = new EmbossFeatureWriteState();
        int tables = 0, unnamed = 0, written = 0;
        EmbossSeq? seq;
        while ((seq = seqall.Next()) is not null)
        {
            WriteFasta(seqOut, seq, ref unnamed);
            written++;
            if (r.Feature && seq.Fttable is { } table)
            {
                string? ftext = EmbossReportWriter.WriteFeatureTable(table, featName, r.RunDate ?? DateTime.Now, "seqret", context, state, tables);
                if (ftext is not null)
                {
                    featOut.Append(ftext);
                    tables++;
                }
            }

            if (r.FirstOnly)
                break;
        }

        if (r.Feature && tables == 0)
            context.Warn($"No features written to output file '{r.OutFeat ?? "o.gff"}'");
        return new EmbossSeqretOutput(seqOut.ToString(), featOut.ToString(), written, [.. context.Messages]) { FeatureFileName = r.OutFeat ?? "" };
    }

    internal static EmbossSeqretOutput Died(EmbossFeatureReadContext context, string parameter = "sequence")
    {
        context.Messages.Add($"Died: {context.Program} terminated: Bad value for '-{parameter}' with -auto defined");
        return new EmbossSeqretOutput("", "", 0, [.. context.Messages]) { Died = true };
    }

    // seqWriteFasta
    internal static void WriteFasta(StringBuilder sb, EmbossSeq s, ref int unnamed, string sequenceType = "gapany")
    {
        // ajSeqoutSetNameDefaultS: a sequence with no name is called EMBOSS_001, EMBOSS_002, ... in the order of the unnamed sequences written
        string name = s.Name.Length > 0 ? s.Name : "EMBOSS_" + (++unnamed).ToString("000", System.Globalization.CultureInfo.InvariantCulture);
        sb.Append('>').Append(name);
        if (s.Sv.Length > 0)
            sb.Append(' ').Append(s.Sv);
        else if (s.Acc.Length > 0)
            sb.Append(' ').Append(s.Acc);
        if (s.Desc.Length > 0)
            sb.Append(' ').Append(s.Desc);
        sb.Append('\n');
        string seq = s.SeqText;
        if (sequenceType == "any")
            seq = new string([.. seq.Where(c => char.IsAsciiLetter(c) || c == '*')]);          // seqClone: ajSeqTypeCheckS of the output type (the type of the input) removes the gaps and every other character
        for (int i = 0; i < seq.Length; i += 60)
            sb.Append(seq, i, Math.Min(60, seq.Length - i)).Append('\n');
    }
}
