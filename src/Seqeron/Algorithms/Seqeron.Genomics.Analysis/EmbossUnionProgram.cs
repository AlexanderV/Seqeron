using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
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
    string OutSeq = "stdout", string? OutSeqFormat = null, DateTime? RunDate = null)
{
    /// <summary><c>-feature</c>: the features of the sequences are read and joined into the features of the output (written to the feature file).</summary>
    public bool Feature { get; init; }

    /// <summary><c>-source</c>: a feature of the type source (source union, tag origid = the name of the sequence) is made for every sequence joined.</summary>
    public bool Source { get; init; }

    /// <summary><c>-offormat2</c>: the format of the feature file (default gff3).</summary>
    public string? OutFeatFormat { get; init; }

    /// <summary><c>-ofname2</c>: the name of the feature file; null = the base name of <c>-outseq</c> (<c>ajFilenameTrimAll</c>) and the feature format as typed (<c>gff</c> by default).</summary>
    public string? OutFeatName { get; init; }
}

/// <summary>The output of <c>union</c>: the sequence file, the overlap file (one number per sequence read) and the messages of the library.</summary>
/// <param name="SequenceText">The text of the sequence output file.</param>
/// <param name="OverlapText">The text of the overlap file (empty without <c>-overlapfile</c>).</param>
/// <param name="Messages">What the library printed on its error stream, in order.</param>
public sealed record EmbossUnionOutput(string SequenceText, string OverlapText, IReadOnlyList<string> Messages)
{
    /// <summary>The text of the feature file (empty without <c>-feature</c>).</summary>
    public string FeatureText { get; init; } = "";

    /// <summary>The name of the feature file.</summary>
    public string FeatureFileName { get; init; } = "";

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
        var seqall = EmbossSeqall.Open("union", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, r.Feature, context);
        if (seqall is null)
        {
            if (!context.Fatal)
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
        EmbossSeqFeatureTable? newTable = null;
        long offset = 0;
        while ((seq = seqall.Next()) is not null)
        {
            if (first)
                _ = seq.Fttable?.NewFtable();                          // uniseq = ajSeqNewSeq(seq): a copy of the sequence with a copy of its feature table (the tags are checked again, which can warn again)
            if (first && r.Feature)
            {
                // ajFeattableNewSeq: a table of the type of the sequence, as long as the sequence is
                newTable = new EmbossSeqFeatureTable(context) { SeqId = seq.Name, Type = seq.Type == "P" ? "P" : "N", Length = seq.Seq.Length };
            }

            if (seq.Fttable is { } ranged && (seq.Begin != 0 || seq.End != 0))
                ranged.SetRange(seq.Begin, seq.End);                   // the table has the range of the sequence when ajSeqTrim cuts it (measured on union -feature -sbegin1: the features are shifted by twice the offset)
            EmbossSeqRange.Trim(seq);                                  // union.c: ajSeqTrim(seq)
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

            if (r.Feature)
            {
                var old = seq.Fttable?.NewFtable();
                EmbossSeqFeature? sourceFeature = null;
                if (r.Source)
                {
                    sourceFeature = newTable!.Type == "P"
                        ? newTable.NewProt("union", "source", unistr.Length - (int)overlap + 1, unistr.Length + residues.Length - (int)overlap, 0f, EmbossLocationBits.None, null)
                        : newTable.NewNuc("union", "source", unistr.Length - (int)overlap + 1, unistr.Length + residues.Length - (int)overlap, 0f, '+', 0, 0, 0, 0, null, null, EmbossLocationBits.None, null);
                    if (sourceFeature is not null)
                        newTable.TagAdd(sourceFeature, "origid", seq.Name);
                }

                if (old is not null)
                    CopyFeatures(old, newTable!, offset - overlap, sourceFeature);
                offset += residues.Length - overlap;
            }

            // ajStrAppendSubS(&unistr, seq, overlap, len - 1): the positions are clipped into the string
            int len = residues.Length;
            int ibegin = Math.Min((int)overlap, len - 1);
            unistr.Append(residues, ibegin, len - ibegin);
            first = false;
            prev = residues;
            _ = seq.Fttable?.NewFtable();                              // prev_seq = ajSeqNewSeq(seq)
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
                EmbossSeqWrite.Write(result, outFormat, name, u.Desc, u.Acc, u.Sv, u.Type, unistr.ToString(), 0, r.RunDate ?? DateTime.Now, u.FullDescMulti, u.Gi, u.Db, u.Setdb);
                written++;
                if (EmbossSeqWrite.Footer(outFormat) is { Length: > 0 } footer)
                    result.Append(footer);                              // seqCleanDasdna / seqCleanDasseq at the close of the file
            }
        }

        if (written == 0 && outFormat is not null)
            context.Warn($"No sequences written to output file '{outName}'");
        string featureText = "", featureFile = "";
        bool tableWritten = false;
        if (r.Feature && newTable is not null)
        {
            string featFormat = r.OutFeatFormat ?? "gff";
            int featIndex = EmbossReportWriter.FeatureFormatIndex(featFormat);
            featureFile = r.OutFeatName ?? TrimAll(outName) + "." + featFormat;
            if (written > 0)
            {
                if (outFormat!.Features)
                    throw new NotSupportedException($"The features of the output format '{outFormat.Name}' are written inside the sequence file, which this port does not do.");
                if (featIndex < 0)
                {
                    // ajFeattabOutOpen / featoutUfoProcess: an unknown feature format is written as gff3 (the message is printed twice, with the two names of the fallback)
                    context.Error($"unknown output feature format '{featFormat}' will write as gff3 instead\n");
                    context.Error($"unknown output feature format '{featFormat}' will write as gff instead\n");
                }

                string featName = EmbossReportWriter.FeatureFormatName(featIndex < 0 ? 1 : featIndex);
                string? text = EmbossReportWriter.WriteFeatureTable(newTable, featName, r.RunDate ?? DateTime.Now, "union", context, new EmbossFeatureWriteState(), 0);
                if (text is null)
                    context.Warn("ajSeqWriteSeq features output failed UFO: ''");
                else
                    featureText = text;
                tableWritten = text is not null;
            }

            if (!tableWritten)
                context.Warn($"No features written to output file '{featureFile}'");
        }

        return new EmbossUnionOutput(result.ToString(), overlapOut.ToString(), [.. context.Messages]) { FeatureText = featureText, FeatureFileName = featureFile };
    }

    // union_CopyFeatures: copies of the features, the positions of the feature (not of the parts of its location) moved by the offset
    private static void CopyFeatures(EmbossSeqFeatureTable oldTable, EmbossSeqFeatureTable newTable, long offset, EmbossSeqFeature? sourceFeature)
    {
        foreach (var gf in oldTable.Features)
        {
            string type = TypeOf(gf);
            var copy = newTable.NewFeat(gf);
            if (copy.Start != 0)
                copy.Start += (int)offset;
            if (copy.End != 0)
                copy.End += (int)offset;
            if (copy.Start2 != 0)
                copy.Start2 += (int)offset;
            if (copy.End2 != 0)
                copy.End2 += (int)offset;
            if (sourceFeature is not null && type == TypeOf(sourceFeature) && copy.GetStart() == sourceFeature.GetStart() && copy.GetEnd() == sourceFeature.GetEnd()
                && gf.Tags.Any(t => string.Equals(t.Name, "origid", StringComparison.OrdinalIgnoreCase)))
                continue;                                               // the source features of an earlier union are not duplicated
            // ajFeattableAddNew
            if (!copy.Flags.HasFlag(EmbossLocationBits.RemoteId) && !copy.Flags.HasFlag(EmbossLocationBits.Label))
            {
                // the positions and the length are unsigned in the library: a start before the sequence (-1) is the largest length
                newTable.Length = unchecked((int)Math.Max((uint)newTable.Length, (uint)copy.Start));
                newTable.Length = unchecked((int)Math.Max((uint)newTable.Length, (uint)copy.End));
            }

            newTable.Groups++;
            copy.Group = newTable.Groups;
            if (copy.Subfeatures is not null)
            {
                foreach (var sub in copy.Subfeatures)
                    sub.Group = newTable.Groups;
            }

            newTable.Features.Add(copy);
        }
    }

    private static string TypeOf(EmbossSeqFeature gf) => EmbossFeatureVocab.TypeExternal((gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types, gf.Type);

    // ajFilenameTrimAll
    private static string TrimAll(string name)
    {
        var m = System.Text.RegularExpressions.Regex.Match(name, ":([A-Za-z0-9_-]+)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant);
        if (m.Success)
            return m.Groups[1].Value;
        m = System.Text.RegularExpressions.Regex.Match(name, "([A-Za-z0-9_-]+)[.][A-Za-z0-9_-]+$", System.Text.RegularExpressions.RegexOptions.CultureInvariant);
        if (m.Success)
            return m.Groups[1].Value;
        m = System.Text.RegularExpressions.Regex.Match(name, "([A-Za-z0-9_-]+)[^A-Za-z0-9_-]*$", System.Text.RegularExpressions.RegexOptions.CultureInvariant);
        return m.Success ? m.Groups[1].Value : "";
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
