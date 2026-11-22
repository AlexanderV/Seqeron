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
    string OutSeq = "stdout", string? OutSeqFormat = null, string? OutFeat = null, string? OutFeatFormat = null, DateTime? RunDate = null)
{
    /// <summary><c>-ossingle</c>: every sequence is written to a file of its own (<c>name.extension</c>, lower case) instead of the file of <c>-outseq</c>. Null: the default of the program (no for seqret, yes for seqretsplit).</summary>
    public bool? OutSingle { get; init; }

    /// <summary><c>-osdirectory</c>: the directory of the files of <c>-ossingle</c>.</summary>
    public string? OutDirectory { get; init; }

    /// <summary><c>-osextension</c>: the extension of the files of <c>-ossingle</c> (default: the name of the output format as typed).</summary>
    public string? OutExtension { get; init; }

    /// <summary>The texts of the configuration files (<c>.embossrc</c>, <c>emboss.default</c>) that define the databases of the USAs <c>db:id</c> (<c>DBNAME name [ type: ... format: ... method: direct directory: ... filename: ... ]</c>, <c>SET name value</c>); the files of the databases are in <see cref="Files"/>.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }

    /// <summary><c>-osname</c>: accepted; the library uses the names of the sequences (the value has no effect on seqret, measured on 6.6.0).</summary>
    public string? OutName { get; init; }

    /// <summary><c>-osdbname</c> (<c>outseq-&gt;Setoutdb</c>): the database name that the writers put in front of the names (<c>db:name</c>: fasta, fastq, gde) or into the identifiers (<c>gnl|db|name</c>: ncbi, gifasta).</summary>
    public string? OutDbName { get; init; }

    /// <summary>The arguments of the command line of the program run, which the SAM writer prints on its <c>@PG</c> line (<c>ajUtilGetCmdline</c>); null: built from the request.</summary>
    public string? CommandLine { get; init; }
}

/// <summary>The output of <c>seqret</c>: the sequence file, the feature file and the messages the library printed.</summary>
/// <param name="SequenceText">The text of the sequence output file.</param>
/// <param name="FeatureText">The text of the feature output file (empty when no feature was written).</param>
/// <param name="Sequences">The number of sequences written.</param>
/// <param name="Messages">What the library printed on its error stream, in order (the ones before the program died, too).</param>
public sealed record EmbossSeqretOutput(string SequenceText, string FeatureText, int Sequences, IReadOnlyList<string> Messages)
{
    /// <summary>The files of <c>-ossingle</c> by path (<c>directory/name.extension</c>); empty when the output is one file.</summary>
    public IReadOnlyDictionary<string, string> SingleFiles { get; init; } = new Dictionary<string, string>();

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
    public static EmbossSeqretOutput Seqret(EmbossSeqretRequest r) => Run(r, "seqret", false);

    /// <summary>
    /// Runs <c>seqretsplit</c> (<c>emboss/seqretsplit.c</c>): <c>seqret</c> that writes every sequence to a file of its own by default (<c>-ossingle</c> is yes in its ACD) and gives a sequence whose name was used by an earlier one
    /// (any case) a new name <c>name.002</c>, <c>name.003</c> ... with the warning <c>Duplicate name 'x' changed to 'x.002'</c>.
    /// </summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program does not stop or dies by a fault.</exception>
    public static EmbossSeqretOutput Seqretsplit(EmbossSeqretRequest r) => Run(r, "seqretsplit", true);

    private static EmbossSeqretOutput Run(EmbossSeqretRequest r, string program, bool split)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        bool outSingle = r.OutSingle ?? split;
        var names = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        int nseqDup = 1;
        var context = new EmbossFeatureReadContext { Program = program };
        var seqall = EmbossSeqall.Open(program, r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, r.Feature, context, "gapany", false, r.Definitions);
        if (seqall is null)
            return Died(context);
        string featFormat = r.OutFeatFormat ?? "gff3";
        int featIndex = EmbossReportWriter.FeatureFormatIndex(featFormat);
        string featName = EmbossReportWriter.FeatureFormatName(featIndex < 0 ? 1 : featIndex);
        var seqOut = new StringBuilder();
        var featOut = new StringBuilder();
        var singles = new Dictionary<string, string>();
        var state = new EmbossFeatureWriteState();
        var runDate = r.RunDate ?? DateTime.Now;
        var writeRun = new EmbossSeqWriteRun { RunDate = runDate, Program = program, CommandLine = r.CommandLine ?? "", Setoutdb = r.OutDbName ?? "" };
        string? lastSingle = null;
        int tables = 0, unnamed = 0, written = 0, count = 0;
        // ajSeqoutUsaProcess: "format::file" names the output format (an empty format is "unknown")
        string outName = r.OutSeq;
        string formatName = r.OutSeqFormat ?? "fasta";
        var um = System.Text.RegularExpressions.Regex.Match(outName, "^([A-Za-z0-9-]*)::?(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant | System.Text.RegularExpressions.RegexOptions.Singleline);
        if (um.Success)
        {
            formatName = um.Groups[1].Value.Length == 0 ? "unknown" : um.Groups[1].Value;
            outName = um.Groups[2].Value;
        }

        writeRun.FileName = outName;

        var outFormat = EmbossSeqWrite.Find(formatName);
        EmbossSeq? seq;
        while (TryNext(seqall, context, out seq))
        {
            if (outFormat is null)
            {
                // ajSeqoutWriteSeq: the format is looked up when the first sequence is written (the message is printed once)
                if (count++ == 0)
                    context.Error($"unknown output sequence format '{formatName}'");
                if (r.FirstOnly)
                    break;
                continue;
            }

            bool isNuc = seq.Type == "N", isProt = seq.Type == "P";
            if (!(isNuc && outFormat.Nucleotide) && !(isProt && outFormat.Protein))
            {
                string kind = isNuc ? "not supported for nucleotide sequences" : "not supported for protein sequences";
                if (!isNuc && !isProt)
                    kind = "failed: unknown sequence type";
                context.Error($"Sequence format '{formatName}' {kind}");
                if (r.FirstOnly)
                    break;
                continue;
            }

            if (split)
                nseqDup = SplitName(names, seq, nseqDup, context);
            if (outFormat.Save)
            {
                // seqClone then seqWriteListAppend: the sequence (a copy cut to its range) is kept and the set is written when the output is closed
                // -ossingle: seqWriteListAppend writes the list at once, to a file that is not open (no file is made for a set format); the writers that print through the null file with a conversion that dereferences it
                // (the file name, the date, ...) kill the process by a segmentation fault, the others write nothing
                var copy = EmbossSeqRange.CopyForList(seq);
                EmbossSeqRange.Trim(copy);
                if (outSingle && EmbossSeqWrite.WritesToNullFileDies(outFormat, Math.Max(copy.Seq.Length, writeRun.Saved.Count == 0 ? 0 : writeRun.Saved.Max(x => x.Seq.Length))))
                    throw new InvalidOperationException("A set format with -ossingle writes the saved list to a file that is not open (the reference program dies by a fault).");
                if (seq.Type.Length > 0)
                    writeRun.OutType = seq.Type;
                EmbossSeqWrite.ListAppend(writeRun, new EmbossSeqSaved { Name = copy.Name, Desc = copy.Desc, Type = copy.Type, Seq = copy.SeqText, Weight = copy.Weight, Gene = seq.Genelist.Count > 0 ? seq.Genelist[0].Name : null }, seq.Entryname);
                count++;
                written++;
                if (r.FirstOnly)
                    break;
                continue;
            }

            if (!EmbossSeqWrite.IsPorted(outFormat))
                throw new NotSupportedException($"The sequence output format '{outFormat.Name}' is not written by this port.");
            if (r.Feature && outFormat.Features)
                throw new NotSupportedException($"The features of the output format '{outFormat.Name}' are written inside the sequence file, which this port does not do.");
            // ajSeqoutSetNameDefaultS: an unnamed sequence is called EMBOSS_001, EMBOSS_002 ... (EMBOSS with -ossingle: the file of every unnamed sequence has the same name)
            string name = seq.Name;
            if (name.Length == 0)
                name = outSingle ? "EMBOSS" : "EMBOSS_" + (++unnamed).ToString("000", System.Globalization.CultureInfo.InvariantCulture);
            string residues = EmbossSeqRange.Cut(seq, out int ibegin, out _, out _);                  // seqClone: the range of the sequence
            var target = outSingle ? new StringBuilder() : seqOut;
            var item = EmbossSeqWrite.Clone(seq, residues, ibegin, writeRun);
            item.Name = name;
            writeRun.Count = outSingle ? 0 : count;
            writeRun.HeaderCount = count;
            EmbossSeqWrite.Write(target, outFormat, item, writeRun);
            if (outSingle)
            {
                string file = (name + "." + (r.OutExtension ?? r.OutSeqFormat ?? "fasta")).ToLowerInvariant();
                string dir = r.OutDirectory ?? "";
                lastSingle = dir.Length > 0 ? Path.Combine(dir, file) : file;
                singles[lastSingle] = target.ToString();
            }

            count++;
            written++;
            if (r.Feature && seq.Fttable is { } table)
            {
                string? ftext = EmbossReportWriter.WriteFeatureTable(table, featName, runDate, "seqret", context, state, tables);
                if (ftext is not null)
                {
                    featOut.Append(ftext);
                    tables++;
                }
            }

            if (r.FirstOnly)
                break;
        }

        if (outFormat is { Save: true } && writeRun.Saved.Count > 0 && !outSingle)
        {
            writeRun.Count = count;                                    // ajSeqoutClose writes the saved list
            EmbossSeqWrite.WriteSaved(seqOut, outFormat, writeRun);
        }

        // seqCleanDasdna / seqCleanDasseq: the footer of a DAS file is written when the file that is open at the end of the run is closed (with -ossingle only the last file has it)
        if (written > 0 && EmbossSeqWrite.Footer(outFormat!) is { Length: > 0 } footer)
        {
            if (lastSingle is not null)
                singles[lastSingle] += footer;
            else
                seqOut.Append(footer);
        }

        if (context.Fatal)
            return new EmbossSeqretOutput(seqOut.ToString(), featOut.ToString(), written, [.. context.Messages]) { Died = true, FeatureFileName = r.OutFeat ?? "", SingleFiles = singles };
        if (count == 0 && !outSingle)                                  // with -ossingle there is no output file that stays empty: the files are opened for the sequences
            context.Warn($"No sequences written to output file '{outName}'");
        if (r.Feature && tables == 0)
            context.Warn($"No features written to output file '{r.OutFeat ?? "o.gff"}'");
        return new EmbossSeqretOutput(seqOut.ToString(), featOut.ToString(), written, [.. context.Messages]) { FeatureFileName = r.OutFeat ?? "", SingleFiles = singles };
    }

    // seqretsplit_Name: a name that is in the table (any case) gets name.002, name.003 ... (the counter of the run grows with every duplicate); the name of the sequence is then put into the table
    private static int SplitName(HashSet<string> table, EmbossSeq seq, int nseq, EmbossFeatureReadContext context)
    {
        if (table.Contains(seq.Name))
        {
            nseq++;
            bool renamed = false;
            for (int i = 2; i <= nseq; i++)
            {
                string newname = seq.Name + "." + i.ToString("000", System.Globalization.CultureInfo.InvariantCulture);
                if (!table.Contains(newname))
                {
                    context.Warn($"Duplicate name '{seq.Name}' changed to '{newname}'");
                    seq.Name = newname;
                    renamed = true;
                    break;
                }
            }

            if (!renamed)
                context.Warn($"Unable to set new name for duplicate sequence number {nseq} '{seq.Name}'");
        }

        table.Add(seq.Name);
        return nseq;
    }

    // ajSeqallNext; a fatal error of the library (ajFatal) ends the program with its message
    internal static bool TryNext(EmbossSeqall seqall, EmbossFeatureReadContext context, [System.Diagnostics.CodeAnalysis.NotNullWhen(true)] out EmbossSeq? seq)
    {
        try
        {
            seq = seqall.Next();
        }
        catch (EmbossFatalException e)
        {
            context.Messages.Add(e.Message);
            context.Fatal = true;
            seq = null;
        }

        return seq is not null;
    }

    internal static EmbossSeqretOutput Died(EmbossFeatureReadContext context, string parameter = "sequence")
    {
        if (!context.Fatal)
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
        string seq = EmbossSeqRange.Cut(s, out _, out _, out _);                              // seqClone: the range of the sequence
        if (sequenceType == "any")
            seq = new string([.. seq.Where(c => char.IsAsciiLetter(c) || c == '*')]);          // seqClone: ajSeqTypeCheckS of the output type (the type of the input) removes the gaps and every other character
        for (int i = 0; i < seq.Length; i += 60)
            sb.Append(seq, i, Math.Min(60, seq.Length - i)).Append('\n');
    }
}
