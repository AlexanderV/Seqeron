// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.FeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S107, S3776, S3358, S1121

using System.Text;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.FeatureTools;

using Seqeron.Genomics.Analysis;

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
    /// <summary><c>-osdbname3</c> .. <c>-osdbname6</c> (<c>outseq-&gt;Setoutdb</c>): the database name that the writers put in front of the names or into the identifiers, of the CDS, mRNA, translation and non-coding outputs.</summary>
    public string? CdsOutDbName { get; init; }

    /// <summary>See <see cref="CdsOutDbName"/>.</summary>
    public string? MrnaOutDbName { get; init; }

    /// <summary>See <see cref="CdsOutDbName"/>.</summary>
    public string? TranslationOutDbName { get; init; }

    /// <summary>See <see cref="CdsOutDbName"/>.</summary>
    public string? RestOutDbName { get; init; }

    /// <summary>The arguments of the command line of the run, which the SAM writer prints on its <c>@PG</c> line (<c>ajUtilGetCmdline</c>); null: empty.</summary>
    public string? CommandLine { get; init; }

    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }

    /// <summary>The date of the run (the header of some sequence formats); null = now.</summary>
    public DateTime? RunDate { get; init; }

    /// <summary>The associated qualifiers of the CDS output (<c>-osname3</c>, <c>-osextension3</c>, <c>-osdirectory3</c>, <c>-ossingle3</c>).</summary>
    public EmbossSeqoutOptions? CdsOptions { get; init; }

    /// <summary>The associated qualifiers of the mRNA output (<c>-osname4</c> ...).</summary>
    public EmbossSeqoutOptions? MrnaOptions { get; init; }

    /// <summary>The associated qualifiers of the translation output (<c>-osname5</c> ...).</summary>
    public EmbossSeqoutOptions? TranslationOptions { get; init; }

    /// <summary>The associated qualifiers of the non-coding output (<c>-osname6</c> ...).</summary>
    public EmbossSeqoutOptions? RestOptions { get; init; }

    /// <summary>The directories (without trailing slash) that exist where the program runs; a file whose directory is not among them cannot be opened. Null: every directory exists.</summary>
    public IReadOnlyCollection<string>? ExistingDirectories { get; init; }
}

/// <summary>The associated qualifiers of one sequence output of <c>coderet</c> (<c>-osname3</c> .. <c>-osname6</c>, <c>-osextension</c>, <c>-osdirectory</c>, <c>-ossingle</c>).</summary>
/// <param name="Name">The base of the default file name (<c>-osname</c>); an output given a file name does not use it.</param>
/// <param name="Extension">The extension of the default file name and of the files of <c>-ossingle</c> (<c>-osextension</c>; default <c>cds</c>, <c>mrna</c>, <c>prot</c>, <c>noncoding</c>).</param>
/// <param name="Directory">The directory of the file (<c>-osdirectory</c>).</param>
/// <param name="Single">One file per sequence (<c>-ossingle</c>), named <c>name.extension</c> in lower case.</param>
public sealed record EmbossSeqoutOptions(string? Name = null, string? Extension = null, string? Directory = null, bool Single = false);

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
    private sealed class Out(string? value, string? format, string defaultName, string extension, EmbossSeqoutOptions? options, string outputType, string parameter)
    {
        public bool Wanted = value is null || value.Length > 0;
        // acdOutFilename: the base name (-osname, else the name of the first sequence; its first word) and the extension (-osextension, else the one of the ACD file)
        public readonly string Extension = options?.Extension is { Length: > 0 } e ? e : extension;
        public string FileName = value is { Length: > 0 } ? value : EmbossOutPath.FmtWord(options?.Name is { Length: > 0 } n ? n : defaultName) + "." + (options?.Extension is { Length: > 0 } e2 ? e2 : extension);
        public readonly string Usa = value is { Length: > 0 } ? value : EmbossOutPath.FmtWord(options?.Name is { Length: > 0 } n2 ? n2 : defaultName) + "." + (options?.Extension is { Length: > 0 } e3 ? e3 : extension);
        public string FormatName = format ?? "fasta";
        public readonly string Directory = EmbossOutPath.FixDirectory(options?.Directory);
        public readonly bool Single = options?.Single ?? false;
        public readonly string OutputType = outputType;
        public readonly string Parameter = parameter;
        public readonly StringBuilder Text = new();
        public readonly Dictionary<string, string> Singles = [];
        public string? LastSingle;
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
        string baseName = EmbossOutPath.DefaultBaseName(seqall.Peek().Name, r.Sequence);
        var cds = new Out(r.CdsOutSeq, r.CdsFormat, baseName, "cds", r.CdsOptions, "nucleotide", "cdsoutseq");
        var mrna = new Out(r.MrnaOutSeq, r.MrnaFormat, baseName, "mrna", r.MrnaOptions, "nucleotide", "mrnaoutseq");
        var prot = new Out(r.TranslationOutSeq, r.TranslationFormat, baseName, "prot", r.TranslationOptions, "protein", "translationoutseq");
        var rest = new Out(r.RestOutSeq, r.RestFormat, baseName, "noncoding", r.RestOptions, "nucleotide", "restoutseq");
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

        // acdSetSeqoutall (ajSeqoutOpen) in the order of the ACD file: a file that cannot be opened ends the program
        foreach (var o in new[] { cds, mrna, prot, rest })
        {
            if (!o.Wanted || o.Single)
                continue;
            if (!EmbossOutPath.DirectoryExists(r.ExistingDirectories, EmbossOutPath.Combine(o.Directory, o.FileName)))
            {
                EmbossOutPath.OpenFailed(context, o.FileName, o.Directory, o.Usa);
                var died2 = EmbossSeqretPrograms.Died(context, o.Parameter);
                return new EmbossCoderetOutput("", new Dictionary<string, string>(), died2.Messages) { Died = true };
            }
        }

        var date = r.RunDate ?? DateTime.Now;
        foreach (var o in new[] { cds, mrna, prot, rest })
        {
            o.Run.RunDate = date;
            o.Run.FileName = o.FileName;
            o.Run.CommandLine = r.CommandLine ?? "";
            o.Run.Setoutdb = (o == cds ? r.CdsOutDbName : o == mrna ? r.MrnaOutDbName : o == prot ? r.TranslationOutDbName : r.RestOutDbName) ?? "";
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
                        PutSeq(seq, EmbossExtractfeatProgram.GetSeqJoin(gf, seq, remote), icds, "cds", false, cds, context, date, r);
                    }

                    if (mrna.Wanted && EmbossExtractfeatProgram.TypeMatchC(gf, "mRNA"))
                    {
                        imrna++;
                        PutSeq(seq, EmbossExtractfeatProgram.GetSeqJoin(gf, seq, remote), imrna, "mrna", false, mrna, context, date, r);
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
                        PutSeq(seq, trn, itran, "pro", true, prot, context, date, r);
                    }
                }

                if (rest.Wanted)
                    PutRest(seq, copyseq!.ToString(), rest, context, date, r);
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
                if (!o.Single)
                    EmbossSeqWrite.CloseSet(o.Text, o.Format, o.Run, o.Count);
                else if (o.Format is { Writer: "Bam" } && EmbossSeqWrite.FinishBam(o.Run) is { } bam)
                    o.Singles[o.LastSingle!] += bam;
                if (o.Count > 0 && o.Format is not null && EmbossSeqWrite.Footer(o.Format) is { Length: > 0 } footer)
                {
                    // seqCleanDasdna / seqCleanDasseq at the close of the file that is open at the end of the run (with -ossingle only the last file has it)
                    if (o.LastSingle is not null)
                        o.Singles[o.LastSingle] += footer;
                    else if (!o.Single)
                        o.Text.Append(footer);
                }

                if (!o.Single)
                    files[EmbossOutPath.Combine(o.Directory, o.FileName)] = o.Text.ToString();
                foreach (var kv in o.Singles)
                    files[kv.Key] = kv.Value;
                if (o.Count == 0 && !o.Single)
                    context.Warn($"No sequences written to output file '{o.FileName}'");
            }
        }

        return new EmbossCoderetOutput(log.ToString(), files, [.. context.Messages]) { Died = context.Fatal };
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
            // ajStrFmtLowerSub(Pseqstr, ajFeatGetStart(gf)-1, ajFeatGetEnd(gf)-1): ajFeatGetStart / ajFeatGetEnd return an unsigned int, so the positions are unsigned before they become the ajlong of the call:
            // the negative positions of the features of a reversed GFF3 entry (its table has no length: ajFeatReverse gives 1 + 0 - End) are huge numbers that are clipped to the last residue, not positions counted from the end
            int begin = EmbossInputSequence.PositionStart(seq.Length, 0, unchecked((uint)(gf.GetStart() - 1)));
            int end = EmbossInputSequence.PositionStart(seq.Length, begin, unchecked((uint)(gf.GetEnd() - 1)));
            for (int i = begin; i <= end && i < seq.Length; i++)
                seq[i] = char.ToLowerInvariant(seq[i]);
        }
    }

    // coderet_put_seq
    private static void PutSeq(EmbossSeq seq, string strseq, int n, string name, bool protein, Out o, EmbossFeatureReadContext context, DateTime date, EmbossCoderetRequest r)
    {
        string fn = $"{seq.Acc}_{name}_{n}".ToLowerInvariant();
        Write(o, fn, protein ? "P" : "N", strseq, context, date, r);
    }

    // coderet_put_rest
    private static void PutRest(EmbossSeq seq, string copy, Out o, EmbossFeatureReadContext context, DateTime date, EmbossCoderetRequest r)
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
                    Write(o, $"{seq.Acc}_noncoding_{i + 1}".ToLowerInvariant(), "N", Sub(strseq, i, j - 1), context, date, r);
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
            Write(o, $"{seq.Acc}_noncoding_{i + 1}".ToLowerInvariant(), "N", Sub(strseq, i, j), context, date, r);
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
    private static void Write(Out o, string name, string type, string residues, EmbossFeatureReadContext context, DateTime date, EmbossCoderetRequest r)
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
        if (o.Format.Save)
        {
            // -ossingle: seqWriteListAppend writes the list at once, to a file that is not open (see seqret)
            if (o.Single && EmbossSeqWrite.WritesToNullFileDies(o.Format, Math.Max(residues.Length, o.Run.Saved.Count == 0 ? 0 : o.Run.Saved.Max(x => x.Seq.Length))))
                throw new InvalidOperationException("A set format with -ossingle writes the saved list to a file that is not open (the reference program dies by a fault).");
            EmbossSeqWrite.Save(o.Format, o.Run, name, "", type, residues);
        }
        else
        {
            var target = o.Single ? new StringBuilder() : o.Text;
            if (o.Format.Writer == "Debug")
            {
                o.Run.Debug = new EmbossSeqDebugInfo
                {
                    Usa = o.Usa, InputFormat = "", OutputFormat = o.FormatName, Filename = o.FileName, Directory = o.Directory, Entryname = name, Single = o.Single, Features = false, OutputType = o.OutputType,
                    FileNameOpen = o.Single ? (name + "." + o.Extension).ToLowerInvariant() : o.FileName, Extension = o.Extension,
                };
            }

            EmbossSeqWrite.Write(target, o.Format, name, "", "", "", type, residues, o.Single ? 0 : o.Count, date, null, run: o.Run, headerCount: o.Count);
            if (o.Single)
            {
                string file = (name + "." + o.Extension).ToLowerInvariant();
                string path = EmbossOutPath.Combine(o.Directory, file);
                if (!EmbossOutPath.DirectoryExists(r.ExistingDirectories, path))
                {
                    if (!EmbossOutPath.SingleOpenSurvives(o.Format))
                        throw new InvalidOperationException("The directory of a file of -ossingle does not exist: the library warns and the reference program dies by a segmentation fault.");
                    EmbossOutPath.SingleOpenFailed(context, file, o.Directory);
                }
                else
                {
                    o.LastSingle = path;
                    o.Singles[o.LastSingle] = target.ToString();
                }
            }
        }

        o.Count++;
    }
}
