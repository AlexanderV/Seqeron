#pragma warning disable S107, S3776, S3358, S1121

using System.Text;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>coderet</c> (A5-8c3h, B09 F259): the parameters of its ACD file.</summary>
/// <param name="Sequence">The value of <c>-seqall</c> (the USA): nucleotide sequences with feature tables.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers.</param>
/// <param name="CdsOutSeq"><c>-cdsoutseq</c>: the file of the coding sequences ([format::]file); null = <c>NAME.cds</c> (NAME: the first sequence, lower case); empty = not wanted.</param>
/// <param name="MrnaOutSeq"><c>-mrnaoutseq</c>: null = <c>NAME.mrna</c>.</param>
/// <param name="TranslationOutSeq"><c>-translationoutseq</c>: null = <c>NAME.prot</c>.</param>
/// <param name="RestOutSeq"><c>-restoutseq</c>: null = <c>NAME.noncoding</c>.</param>
/// <param name="CdsFormat"><c>-osformat3</c> (default fasta).</param>
/// <param name="MrnaFormat"><c>-osformat4</c>.</param>
/// <param name="TranslationFormat"><c>-osformat5</c>.</param>
/// <param name="RestFormat"><c>-osformat6</c>.</param>
public sealed record EmbossCoderetRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, string? CdsOutSeq = null, string? MrnaOutSeq = null, string? TranslationOutSeq = null,
    string? RestOutSeq = null, string? CdsFormat = null, string? MrnaFormat = null, string? TranslationFormat = null, string? RestFormat = null)
{
    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }

    /// <summary>The date of the run (the header of some sequence formats); null = now.</summary>
    public DateTime? RunDate { get; init; }
}

/// <summary>The output of <c>coderet</c>.</summary>
/// <param name="Report">The text of the log file (<c>-outfile</c>): the counts of every sequence.</param>
/// <param name="Files">The sequence files by name.</param>
/// <param name="Messages">What the library printed on its error stream.</param>
public sealed record EmbossCoderetOutput(string Report, IReadOnlyDictionary<string, string> Files, IReadOnlyList<string> Messages)
{
    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}

/// <summary>
/// <c>coderet</c> of EMBOSS 6.6.0 (<c>emboss/coderet.c</c>): the coding sequences (<c>ajFeatGetSeqJoin</c> of the CDS features), the messenger RNAs (mRNA features), the translations (the <c>translation</c> tag of the CDS features, <c>ajFeatGetTranslation</c>) and the
/// non-coding regions (what <c>ajFeatLocMark</c> of the CDS, mRNA and exon features leaves in upper case) of the sequences with feature tables, named <c>ACCESSION_cds_N</c>, <c>_mrna_N</c>, <c>_pro_N</c> and <c>_noncoding_START</c>.
/// </summary>
public static class EmbossCoderetProgram
{
    private sealed class Out(string? value, string? format, string defaultName, string extension)
    {
        public bool Wanted = value is null || value.Length > 0;
        public string FileName = value is { Length: > 0 } ? value : defaultName + "." + extension;
        public string FormatName = format ?? "fasta";
        public readonly StringBuilder Text = new();
        public int Count;
        public EmbossSeqOutFormat? Format;
        public bool UnknownReported;
        public readonly EmbossSeqWriteRun Run = new() { Program = "coderet" };
    }

    /// <summary>Runs <c>coderet</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    /// <exception cref="NotSupportedException">An output format is not written by this port.</exception>
    public static EmbossCoderetOutput Coderet(EmbossCoderetRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "coderet" };
        var seqall = EmbossSeqall.Open("coderet", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, true, context, "dna", false, r.Definitions);
        if (seqall is null)
        {
            var died = EmbossSeqretPrograms.Died(context, "seqall");
            return new EmbossCoderetOutput("", new Dictionary<string, string>(), died.Messages) { Died = true };
        }

        // acdOutFilename: the output files without a name are called after the first sequence
        string baseName = DefaultBaseName(seqall.Peek().Name, r.Sequence);
        var cds = new Out(r.CdsOutSeq, r.CdsFormat, baseName, "cds");
        var mrna = new Out(r.MrnaOutSeq, r.MrnaFormat, baseName, "mrna");
        var prot = new Out(r.TranslationOutSeq, r.TranslationFormat, baseName, "prot");
        var rest = new Out(r.RestOutSeq, r.RestFormat, baseName, "noncoding");
        foreach (var o in new[] { cds, mrna, prot, rest })
        {
            if (!o.Wanted)
                continue;
            // ajSeqoutUsaProcess: "format::file" names the output format (an empty format is "unknown")
            var um = Regex.Match(o.FileName, "^([A-Za-z0-9-]*)::?(.*)$", RegexOptions.CultureInvariant | RegexOptions.Singleline);
            if (um.Success)
            {
                o.FormatName = um.Groups[1].Value.Length == 0 ? "unknown" : um.Groups[1].Value;
                o.FileName = um.Groups[2].Value;
            }

            o.Format = EmbossSeqWrite.Find(o.FormatName);
        }

        var date = r.RunDate ?? DateTime.Now;
        foreach (var o in new[] { cds, mrna, prot, rest })
        {
            o.Run.RunDate = date;
            o.Run.FileName = o.FileName;
        }

        var remote = new EmbossRemoteSeq(r.Files, r.Definitions, context);
        var log = new StringBuilder();
        if (cds.Wanted)
            log.Append("   CDS");
        if (mrna.Wanted)
            log.Append("  mRNA");
        if (rest.Wanted)
            log.Append(" non-c");
        if (prot.Wanted)
            log.Append(" Trans");
        log.Append(" Total Sequence\n");
        if (cds.Wanted)
            log.Append(" =====");
        if (mrna.Wanted)
            log.Append(" =====");
        if (rest.Wanted)
            log.Append(" =====");
        if (prot.Wanted)
            log.Append(" =====");
        log.Append(" ===== ========\n");

        while (EmbossSeqretPrograms.TryNext(seqall, context, out var seq))
        {
            int icds = 0, imrna = 0, irest = 0, itran = 0;
            if (seq.Fttable is { } feat)
            {
                StringBuilder? copyseq = null;
                if (rest.Wanted)
                    copyseq = new StringBuilder(seq.Seq.ToString().ToUpperInvariant());
                foreach (var gf in feat.Features)
                {
                    if (cds.Wanted && EmbossExtractfeatProgram.TypeMatchC(gf, "CDS"))
                    {
                        icds++;
                        PutSeq(seq, EmbossExtractfeatProgram.GetSeqJoin(gf, seq, remote), icds, "cds", false, cds, context, date);
                    }

                    if (mrna.Wanted && EmbossExtractfeatProgram.TypeMatchC(gf, "mRNA"))
                    {
                        imrna++;
                        PutSeq(seq, EmbossExtractfeatProgram.GetSeqJoin(gf, seq, remote), imrna, "mrna", false, mrna, context, date);
                    }

                    if (rest.Wanted)
                    {
                        if (EmbossExtractfeatProgram.TypeMatchC(gf, "CDS") || EmbossExtractfeatProgram.TypeMatchC(gf, "mRNA") || EmbossExtractfeatProgram.TypeMatchC(gf, "exon"))
                            LocMark(gf, feat, copyseq!);
                        else
                            irest++;
                    }

                    if (prot.Wanted && EmbossExtractfeatProgram.TypeMatchC(gf, "CDS") && GetTagTranslation(gf) is { } trn)
                    {
                        itran++;
                        PutSeq(seq, trn, itran, "pro", true, prot, context, date);
                    }
                }

                if (rest.Wanted)
                    PutRest(seq, copyseq!.ToString(), rest, context, date);
            }

            if (cds.Wanted)
                log.Append(icds.ToString().PadLeft(6));
            if (mrna.Wanted)
                log.Append(imrna.ToString().PadLeft(6));
            if (rest.Wanted)
                log.Append(irest.ToString().PadLeft(6));
            if (prot.Wanted)
                log.Append(itran.ToString().PadLeft(6));
            log.Append((icds + imrna + irest).ToString().PadLeft(6)).Append(' ').Append(seq.Name).Append('\n');
        }

        var files = new Dictionary<string, string>();
        if (!context.Fatal)
        {
            foreach (var o in new[] { cds, mrna, rest, prot })
            {
                if (!o.Wanted)
                    continue;
                EmbossSeqWrite.CloseSet(o.Text, o.Format, o.Run, o.Count);
                if (o.Count > 0 && o.Format is not null && EmbossSeqWrite.Footer(o.Format) is { Length: > 0 } footer)
                    o.Text.Append(footer);                              // seqCleanDasdna / seqCleanDasseq at the close of the file
                files[o.FileName] = o.Text.ToString();
                if (o.Count == 0)
                    context.Warn($"No sequences written to output file '{o.FileName}'");
            }
        }

        return new EmbossCoderetOutput(log.ToString(), files, [.. context.Messages]) { Died = context.Fatal };
    }

    // acdInFileSave (the name of the first sequence, else the USA) + ajFilenameTrimAll + ajStrFmtLower, then ajStrFmtWord of acdOutFilename
    private static string DefaultBaseName(string seqName, string usa)
    {
        string name = TrimAll(seqName.Length > 0 ? seqName : usa).ToLowerInvariant();
        int blank = name.IndexOfAny([' ', '\t', '\n', '\r']);
        return blank >= 0 ? name[..blank] : name;
    }

    // ajFilenameTrimAll
    private static string TrimAll(string name)
    {
        var m = Regex.Match(name, ":([A-Za-z0-9_-]+)$", RegexOptions.CultureInvariant);
        if (m.Success)
            return m.Groups[1].Value;
        m = Regex.Match(name, "([A-Za-z0-9_-]+)[.][A-Za-z0-9_-]+$", RegexOptions.CultureInvariant);
        if (m.Success)
            return m.Groups[1].Value;
        m = Regex.Match(name, "([A-Za-z0-9_-]+)[^A-Za-z0-9_-]*$", RegexOptions.CultureInvariant);
        return m.Success ? m.Groups[1].Value : "";
    }

    // ajFeatGetTagC(feature, "translation", 1): the first tag of that name (any case), or a note that starts with "*translation:"
    private static string? GetTagTranslation(EmbossSeqFeature gf)
    {
        const string tname = "translation";
        foreach (var (name, value) in gf.Tags)
        {
            if (string.Equals(name, tname, StringComparison.OrdinalIgnoreCase))
                return value;
            if (string.Equals(name, "note", StringComparison.OrdinalIgnoreCase) && value.Length > tname.Length + 3 && value[0] == '*' && value.AsSpan(1).StartsWith(tname, StringComparison.OrdinalIgnoreCase) && value[tname.Length + 1] == ':')
                return value[(tname.Length + 3)..];
        }

        return null;
    }

    // ajFeatLocMark: the sequence in the locations of the features of the group of the feature becomes lower case
    private static void LocMark(EmbossSeqFeature thys, EmbossSeqFeatureTable table, StringBuilder seq)
    {
        foreach (var gf in table.Features)
        {
            if (gf.Group != thys.Group)
                continue;
            if (gf.Flags.HasFlag(EmbossLocationBits.BetweenSeq) || gf.Flags.HasFlag(EmbossLocationBits.RemoteId))
                continue;
            // ajStrFmtLowerSub: the positions are clipped into the string
            int begin = EmbossInputSequence.PositionStart(seq.Length, 0, gf.GetStart() - 1);
            int end = EmbossInputSequence.PositionStart(seq.Length, begin, gf.GetEnd() - 1);
            for (int i = begin; i <= end && i < seq.Length; i++)
                seq[i] = char.ToLowerInvariant(seq[i]);
        }
    }

    // coderet_put_seq
    private static void PutSeq(EmbossSeq seq, string strseq, int n, string name, bool protein, Out o, EmbossFeatureReadContext context, DateTime date)
    {
        string fn = $"{seq.Acc}_{name}_{n}".ToLowerInvariant();
        Write(o, fn, protein ? "P" : "N", strseq, context, date);
    }

    // coderet_put_rest
    private static void PutRest(EmbossSeq seq, string copy, Out o, EmbossFeatureReadContext context, DateTime date)
    {
        string strseq = seq.Seq.ToString();
        int i = 0, j = 0;
        bool isup = false;
        foreach (char c in copy)
        {
            if (char.IsAsciiLetterLower(c))
            {
                if (isup)
                {
                    Write(o, $"{seq.Acc}_noncoding_{i + 1}".ToLowerInvariant(), "N", Sub(strseq, i, j - 1), context, date);
                    isup = false;
                }
            }
            else if (!isup)
            {
                i = j;
                isup = true;
            }

            j++;
        }

        if (isup)
            Write(o, $"{seq.Acc}_noncoding_{i + 1}".ToLowerInvariant(), "N", Sub(strseq, i, j), context, date);
    }

    // ajStrAssignSubS: the positions are clipped into the string
    private static string Sub(string s, long pos1, long pos2)
    {
        int ibegin = EmbossInputSequence.PositionStart(s.Length, 0, pos1);
        int iend = EmbossInputSequence.PositionStart(s.Length, ibegin, pos2);
        int ilen = iend - ibegin + 1;
        return ilen > 0 && ibegin < s.Length ? s.Substring(ibegin, Math.Min(ilen, s.Length - ibegin)) : "";
    }

    // ajSeqoutWriteSeq
    private static void Write(Out o, string name, string type, string residues, EmbossFeatureReadContext context, DateTime date)
    {
        if (o.Format is null)
        {
            if (!o.UnknownReported)
                context.Error($"unknown output sequence format '{o.FormatName}'");
            o.UnknownReported = true;
            return;
        }

        bool isNuc = type == "N";
        if (!(isNuc && o.Format.Nucleotide) && !(!isNuc && o.Format.Protein))
        {
            context.Error($"Sequence format '{o.FormatName}' {(isNuc ? "not supported for nucleotide sequences" : "not supported for protein sequences")}");
            return;
        }

        if (!EmbossSeqWrite.IsPorted(o.Format))
            throw new NotSupportedException($"The sequence output format '{o.Format.Name}' is not written by this port.");
        if (!EmbossSeqWrite.Save(o.Format, o.Run, name, "", type, residues))
            EmbossSeqWrite.Write(o.Text, o.Format, name, "", "", "", type, residues, o.Count, date, null, run: o.Run);
        o.Count++;
    }
}
