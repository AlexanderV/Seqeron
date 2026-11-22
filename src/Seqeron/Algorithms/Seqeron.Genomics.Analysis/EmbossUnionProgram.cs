using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>union</c> (A5-8c3h, B09 F253).</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA).</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers (<c>-sformat1</c>, <c>-slower1</c> ...).</param>
/// <param name="FindOverlap"><c>-findoverlap</c>: the bases by which a sequence overlaps the end of the one before are joined once.</param>
/// <param name="OverlapFile"><c>-overlapfile</c> is given: the overlap of every sequence is written to it.</param>
/// <param name="OutSeq"><c>-outseq</c>: the output file name, optionally with a <c>format::</c> prefix.</param>
/// <param name="OutSeqFormat"><c>-osformat2</c>: the sequence output format (default <c>fasta</c>).</param>
/// <param name="RunDate">The date of the run; null = now.</param>
public sealed record EmbossUnionRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, bool FindOverlap = false, bool OverlapFile = false,
    string OutSeq = "stdout", string? OutSeqFormat = null, DateTime? RunDate = null);

/// <summary>The output of <c>union</c>: the sequence file, the overlap file (one number per sequence read) and the messages of the library.</summary>
/// <param name="SequenceText">The text of the sequence output file.</param>
/// <param name="OverlapText">The text of the overlap file (empty without <c>-overlapfile</c>).</param>
/// <param name="Messages">What the library printed on its error stream, in order.</param>
public sealed record EmbossUnionOutput(string SequenceText, string OverlapText, IReadOnlyList<string> Messages)
{
    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}

/// <summary><c>union</c> of EMBOSS 6.6.0 (<c>emboss/union.c</c>): the sequences of a file joined into one sequence that carries the name, description and type of the first.</summary>
public static class EmbossUnionProgram
{
    /// <summary>Runs <c>union</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program does not stop or dies by a fault.</exception>
    /// <exception cref="NotSupportedException">The output format is not written by this port.</exception>
    public static EmbossUnionOutput Union(EmbossUnionRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        var context = new EmbossFeatureReadContext { Program = "union" };
        var seqall = EmbossSeqall.Open("union", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, false, context);
        if (seqall is null)
        {
            context.Messages.Add("Died: union terminated: Bad value for '-sequence' with -auto defined");
            return new EmbossUnionOutput("", "", [.. context.Messages]) { Died = true };
        }

        var overlapOut = new StringBuilder();
        var unistr = new StringBuilder();
        EmbossSeq? uniseq = null;
        string prev = "";
        long overlap = 0;
        bool first = true;
        EmbossSeq? seq;
        while ((seq = seqall.Next()) is not null)
        {
            if (first)
                uniseq = seq;
            string residues = seq.SeqText;
            if (r.FindOverlap)
            {
                if (!first)
                    overlap = GetOverlap(prev, residues);
                if (r.OverlapFile)
                    overlapOut.Append(overlap.ToString(CultureInfo.InvariantCulture)).Append('\n');
            }

            // ajStrAppendSubS(&unistr, seq, overlap, len - 1): the positions are clipped into the string
            int len = residues.Length;
            int ibegin = Math.Min((int)overlap, len - 1);
            unistr.Append(residues, ibegin, len - ibegin);
            first = false;
            prev = residues;
        }

        var result = new StringBuilder();
        var outFormat = EmbossSeqWrite.Find(r.OutSeqFormat ?? "fasta");
        string outName = r.OutSeq;
        string formatName = r.OutSeqFormat ?? "fasta";
        var um = System.Text.RegularExpressions.Regex.Match(outName, "^([A-Za-z0-9-]*)::?(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant | System.Text.RegularExpressions.RegexOptions.Singleline);
        if (um.Success)
        {
            formatName = um.Groups[1].Value.Length == 0 ? "unknown" : um.Groups[1].Value;
            outName = um.Groups[2].Value;
            outFormat = EmbossSeqWrite.Find(formatName);
        }

        int written = 0;
        if (outFormat is null)
            context.Error($"unknown output sequence format '{formatName}'");
        else
        {
            var u = uniseq!;
            bool isNuc = u.Type == "N", isProt = u.Type == "P";
            if (!(isNuc && outFormat.Nucleotide) && !(isProt && outFormat.Protein))
            {
                string kind = isNuc ? "not supported for nucleotide sequences" : "not supported for protein sequences";
                if (!isNuc && !isProt)
                    kind = "failed: unknown sequence type";
                context.Error($"Sequence format '{formatName}' {kind}");
            }
            else
            {
                if (!EmbossSeqWrite.IsPorted(outFormat))
                    throw new NotSupportedException($"The sequence output format '{outFormat.Name}' is not written by this port.");
                string name = u.Name.Length > 0 ? u.Name : "EMBOSS_001";
                EmbossSeqWrite.Write(result, outFormat, name, u.Desc, u.Acc, u.Sv, u.Type, unistr.ToString(), 0, r.RunDate ?? DateTime.Now, u.FullDescMulti);
                written++;
            }
        }

        if (written == 0 && outFormat is not null)
            context.Warn($"No sequences written to output file '{outName}'");
        return new EmbossUnionOutput(result.ToString(), overlapOut.ToString(), [.. context.Messages]);
    }

    // union_GetOverlap: the longest end of the first sequence that is the start of the second
    private static long GetOverlap(string first, string second)
    {
        int i = Math.Min(first.Length, second.Length);
        for (; i >= 0; --i)
        {
            if (string.CompareOrdinal(first, first.Length - i, second, 0, i) == 0)
                return i;
        }

        return 0;
    }
}
