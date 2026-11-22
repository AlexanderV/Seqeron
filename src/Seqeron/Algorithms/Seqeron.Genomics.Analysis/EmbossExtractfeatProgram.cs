#pragma warning disable S107, S3776, S1121, S3358, S3398

using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>extractfeat</c> (A5-8c3h, B09 F258): the options of its ACD file.</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA).</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers.</param>
/// <param name="Before"><c>-before</c>: sequence before each feature (negative: the start is that far before the END of the feature).</param>
/// <param name="After"><c>-after</c>: sequence after each feature (negative: the end is that far after the START of the feature).</param>
/// <param name="Join"><c>-join</c>: the parts of a join are written as one sequence.</param>
/// <param name="FeatInName"><c>-featinname</c>: the type of the feature is added to the name.</param>
/// <param name="Describe"><c>-describe</c>: the tags (names with wildcards, separated by blanks, commas, semicolons or bars) added with their values to the description.</param>
/// <param name="Source"><c>-source</c>: sources of the features to extract (default <c>*</c>).</param>
/// <param name="Type"><c>-type</c>: types of the features to extract (default <c>*</c>).</param>
/// <param name="Sense"><c>-sense</c>: 0 any, 1 forward, -1 reverse.</param>
/// <param name="MinScore"><c>-minscore</c>.</param>
/// <param name="MaxScore"><c>-maxscore</c> (both 0: the score is not tested).</param>
/// <param name="Tag"><c>-tag</c>: tags a feature must have (default <c>*</c>).</param>
/// <param name="Value"><c>-value</c>: values the tag must have (default <c>*</c>).</param>
/// <param name="OutSeq"><c>-outseq</c>: the output file name with an optional <c>format::</c> prefix.</param>
/// <param name="OutSeqFormat"><c>-osformat</c>: the sequence output format (default fasta).</param>
public sealed record EmbossExtractfeatRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, int Before = 0, int After = 0, bool Join = false, bool FeatInName = false,
    string Describe = "", string Source = "*", string Type = "*", int Sense = 0, float MinScore = 0f, float MaxScore = 0f, string Tag = "*", string Value = "*", string OutSeq = "stdout", string? OutSeqFormat = null)
{
    /// <summary><c>-osdbname</c> (<c>outseq-&gt;Setoutdb</c>): the database name that the writers put in front of the names or into the identifiers.</summary>
    public string? OutDbName { get; init; }

    /// <summary>The arguments of the command line of the run, which the SAM writer prints on its <c>@PG</c> line (<c>ajUtilGetCmdline</c>); null: empty.</summary>
    public string? CommandLine { get; init; }

    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }

    /// <summary>The date of the run (the header of some sequence formats); null = now.</summary>
    public DateTime? RunDate { get; init; }

    /// <summary><c>-ossingle</c>: one file per sequence, named <c>name.extension</c> in lower case; null = false.</summary>
    public bool? OutSingle { get; init; }

    /// <summary><c>-osdirectory</c>: the directory of the files of <c>-ossingle</c>.</summary>
    public string? OutDirectory { get; init; }

    /// <summary><c>-osextension</c>: the extension of the files of <c>-ossingle</c> (default: the name of the output format).</summary>
    public string? OutExtension { get; init; }

    /// <summary>The directories (without trailing slash) that exist where the program runs; a file whose directory is not among them cannot be opened. Null: every directory exists.</summary>
    public IReadOnlyCollection<string>? ExistingDirectories { get; init; }
}

/// <summary>
/// <c>extractfeat</c> of EMBOSS 6.6.0 (<c>emboss/extractfeat.c</c>): the features of the sequences that pass the filters (source, type, sense, score, tag, value) are extracted as sequences (<c>ajFeatGetSeq</c> / <c>ajFeatGetSeqJoin</c>) with
/// the sequence before and after them, named <c>seq_start_end[_type]</c> and described by <c>[type] (tags) description</c>, and written in a sequence output format.
/// </summary>
public static class EmbossExtractfeatProgram
{
    /// <summary>Runs <c>extractfeat</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    /// <exception cref="NotSupportedException">The output format is not written by this port.</exception>
    public static EmbossSeqretOutput Extractfeat(EmbossExtractfeatRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "extractfeat" };
        var seqall = EmbossSeqall.Open("extractfeat", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, true, context, "any", false, r.Definitions);
        if (seqall is null)
            return EmbossSeqretPrograms.Died(context);
        float minscore = r.MinScore, maxscore = r.MaxScore;
        bool testscore = minscore != 0f || maxscore != 0f;
        if (minscore != 0f && maxscore == 0f && minscore > maxscore)
            maxscore = minscore;
        if (minscore == 0f && maxscore != 0f && minscore > maxscore)
            minscore = maxscore;

        string formatName = r.OutSeqFormat ?? "fasta";
        string outName = r.OutSeq;
        var um = System.Text.RegularExpressions.Regex.Match(outName, "^([A-Za-z0-9-]*)::?(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant | System.Text.RegularExpressions.RegexOptions.Singleline);
        if (um.Success)
        {
            formatName = um.Groups[1].Value.Length == 0 ? "unknown" : um.Groups[1].Value;
            outName = um.Groups[2].Value;
        }

        // ajSeqoutOpen (at the start of the run, by acdSetSeqout): -osdirectory is the directory of the output file, too
        string outDir = EmbossOutPath.FixDirectory(r.OutDirectory);
        bool outSingle = r.OutSingle ?? false;
        string mainPath = outSingle ? outName : EmbossOutPath.Combine(outDir, outName);
        if (!outSingle && !EmbossOutPath.DirectoryExists(r.ExistingDirectories, mainPath))
        {
            EmbossOutPath.OpenFailed(context, outName, outDir, r.OutSeq);
            return EmbossSeqretPrograms.Died(context, "outseq");
        }

        var outFormat = EmbossSeqWrite.Find(formatName);
        var output = new StringBuilder();
        var run = new Run(r, context, outFormat, formatName, output, r.RunDate ?? DateTime.Now, new EmbossRemoteSeq(r.Files, r.Definitions, context));
        run.Out.FileName = outName;
        run.Out.Program = "extractfeat";
        run.Out.CommandLine = r.CommandLine ?? "";
        run.Out.Setoutdb = r.OutDbName ?? "";
        while (EmbossSeqretPrograms.TryNext(seqall, context, out var seq))
        {
            // the features of the sequence (a copy), the ones that pass the filters
            var table = seq.Fttable?.NewFtable();
            if (table is not null)
            {
                bool tagsmatch = false;
                table.Features.RemoveAll(gf => !MatchFeature(gf, r.Source, r.Type, r.Sense, testscore, minscore, maxscore, r.Tag, r.Value, ref tagsmatch));
            }

            run.Extract(seq, table);
        }

        if (!run.Single)
            EmbossSeqWrite.CloseSet(output, outFormat, run.Out, run.Count);
        else if (outFormat is { Writer: "Bam" } && EmbossSeqWrite.FinishBam(run.Out) is { } bam)
            run.Singles[run.LastSingle!] += bam;
        if (run.Written > 0 && outFormat is not null && EmbossSeqWrite.Footer(outFormat) is { Length: > 0 } footer)
        {
            // seqCleanDasdna / seqCleanDasseq at the close of the file that is open at the end of the run (with -ossingle only the last file has it)
            if (run.LastSingle is not null)
                run.Singles[run.LastSingle] += footer;
            else if (!run.Single)
                output.Append(footer);
        }

        if (context.Fatal)
            return new EmbossSeqretOutput(output.ToString(), "", run.Written, [.. context.Messages]) { Died = true, SingleFiles = run.Singles, SequenceFileName = mainPath };
        if (run.Count == 0 && !run.Single)
            context.Warn($"No sequences written to output file '{outName}'");
        return new EmbossSeqretOutput(output.ToString(), "", run.Written, [.. context.Messages]) { SingleFiles = run.Singles, SequenceFileName = mainPath };
    }

    private sealed class Run(EmbossExtractfeatRequest r, EmbossFeatureReadContext context, EmbossSeqOutFormat? format, string formatName, StringBuilder output, DateTime date, EmbossRemoteSeq remote)
    {
        public int Written, Count;
        public readonly EmbossSeqWriteRun Out = new() { RunDate = date };
        public readonly bool Single = r.OutSingle ?? false;
        public readonly Dictionary<string, string> Singles = [];
        public string? LastSingle;
        private bool _unknownReported;

        // extractfeat_FeatSeqExtract
        public void Extract(EmbossSeq seq, EmbossSeqFeatureTable? table)
        {
            if (table is null || table.Features.Count == 0)
                return;
            foreach (var gftop in table.Features)
            {
                bool compall = false, sense = true;
                int firstpos, lastpos;
                string describeout = MatchPatternDescribe(gftop, r.Describe);
                if (gftop.Subfeatures is { Count: > 0 } subs && !r.Join)
                {
                    foreach (var gf in subs)
                    {
                        if (StrandOf(gf) == '-')
                            sense = false;
                        string type = TypeOf(gf);
                        string featseq = GetSeq(gf, seq, remote);
                        firstpos = gf.GetStart() - 1;
                        lastpos = gf.GetEnd() - 1;
                        WriteOut(seq, featseq, compall, sense, firstpos, lastpos, type, describeout);
                    }
                }
                else
                {
                    string featseq;
                    if (r.Join)
                    {
                        compall = gftop.Flags.HasFlag(EmbossLocationBits.ComplementMain);
                        if (!compall && StrandOf(gftop) == '-')
                            sense = false;
                        featseq = GetSeqJoin(gftop, seq, remote);
                    }
                    else
                    {
                        if (StrandOf(gftop) == '-')
                            sense = false;
                        featseq = GetSeq(gftop, seq, remote);
                    }

                    firstpos = gftop.GetStart() - 1;
                    lastpos = gftop.GetEnd() - 1;
                    WriteOut(seq, featseq, compall, sense, firstpos, lastpos, TypeOf(gftop), describeout);
                }
            }
        }

        // extractfeat_WriteOut
        private void WriteOut(EmbossSeq seq, string featstr, bool compall, bool sense, int firstpos, int lastpos, string type, string describestr)
        {
            bool forward = sense;
            if (compall)
                forward = false;
            if (featstr.Length == 0)
            {
                context.Warn($"feature {seq.Name}_{firstpos + 1}_{lastpos + 1} [{type}] not written out because it has zero length\n");
                return;
            }

            featstr = BeforeAfter(seq, featstr, firstpos, lastpos, r.Before, r.After, forward);
            string name = $"{seq.Name}_{firstpos + 1}_{lastpos + 1}";
            if (r.FeatInName)
                name += "_" + type;
            string desc = "[" + type + "] " + describestr + seq.Desc;
            // ajSeqoutWriteSeq
            if (format is null)
            {
                if (!_unknownReported)
                    context.Error($"unknown output sequence format '{formatName}'");
                _unknownReported = true;
                Count++;
                return;
            }

            bool isNuc = seq.Type == "N", isProt = seq.Type == "P";
            if (!(isNuc && format.Nucleotide) && !(isProt && format.Protein))
            {
                string kind = isNuc ? "not supported for nucleotide sequences" : "not supported for protein sequences";
                if (!isNuc && !isProt)
                    kind = "failed: unknown sequence type";
                context.Error($"Sequence format '{formatName}' {kind}");
                return;
            }

            if (!EmbossSeqWrite.IsPorted(format))
                throw new NotSupportedException($"The sequence output format '{format.Name}' is not written by this port.");
            if (format.Save)
            {
                // -ossingle: seqWriteListAppend writes the list at once, to a file that is not open (see seqret)
                if (Single && EmbossSeqWrite.WritesToNullFileDies(format, Math.Max(featstr.Length, Out.Saved.Count == 0 ? 0 : Out.Saved.Max(x => x.Seq.Length))))
                    throw new InvalidOperationException("A set format with -ossingle writes the saved list to a file that is not open (the reference program dies by a fault).");
                EmbossSeqWrite.Save(format, Out, name, desc, seq.Type, featstr);
            }
            else
            {
                var target = Single ? new StringBuilder() : output;
                if (format.Writer == "Debug")
                {
                    string ext = r.OutExtension ?? r.OutSeqFormat ?? "fasta";
                    string dirName = r.OutDirectory ?? "";
                    if (dirName.Length > 0 && !dirName.EndsWith('/'))
                        dirName += "/";
                    Out.Debug = new EmbossSeqDebugInfo
                    {
                        Usa = r.OutSeq, InputFormat = "", OutputFormat = formatName, Filename = Out.FileName, Directory = dirName, Entryname = "", Single = Single, Features = false, OutputType = "any",
                        FileNameOpen = Single ? (name + "." + ext).ToLowerInvariant() : Out.FileName, Extension = ext,
                    };
                }

                EmbossSeqWrite.Write(target, format, name, desc, "", "", seq.Type, featstr, Single ? 0 : Count, date, null, run: Out, headerCount: Count);
                if (Single)
                {
                    string file = (name + "." + (r.OutExtension ?? r.OutSeqFormat ?? "fasta")).ToLowerInvariant();
                    string dirFixed = EmbossOutPath.FixDirectory(r.OutDirectory);
                    string path = EmbossOutPath.Combine(dirFixed, file);
                    if (!EmbossOutPath.DirectoryExists(r.ExistingDirectories, path))
                    {
                        if (!EmbossOutPath.SingleOpenSurvives(format))
                            throw new InvalidOperationException("The directory of a file of -ossingle does not exist: the library warns and the reference program dies by a segmentation fault.");
                        EmbossOutPath.SingleOpenFailed(context, file, dirFixed);
                    }
                    else
                    {
                        LastSingle = path;
                        Singles[LastSingle] = target.ToString();
                    }
                }
            }

            Count++;
            Written++;
        }

        // extractfeat_BeforeAfter
        private static string BeforeAfter(EmbossSeq seq, string featstr0, int firstpos, int lastpos, int before, int after, bool sense)
        {
            var featstr = new StringBuilder(featstr0);
            int featlen = featstr.Length - 1;
            int start = 0, end = featlen;
            if (before < 0)
                start = end + before + 1;
            if (after < 0)
                end = start - after - 1;
            if (end < start)
            {
                // ajWarn("Extraction region end less than start for %S [%d-%d]")
                return featstr0;
            }

            if (before < 0 || after < 0)
            {
                // ajStrKeepRange(featstr, start < 0 ? 0 : start, end > featlen ? featlen : end)
                int b = EmbossInputSequence.PositionStart(featstr.Length, 0, start < 0 ? 0 : start);
                int e = EmbossInputSequence.PositionStart(featstr.Length, 0, end > featlen ? featlen : end);
                if (e == featstr.Length)
                    e--;
                string kept = e < b ? "" : featstr.ToString(b, e - b + 1);
                featstr.Clear().Append(kept);
            }

            if (start < 0)
            {
                int len = -start;
                if (sense)
                    GetRegionPad(seq, featstr, firstpos - len, firstpos - 1, sense, true);
                else
                    GetRegionPad(seq, featstr, lastpos + 1, lastpos + len, sense, true);
            }

            if (end > featlen)
            {
                int len = end - featlen;
                if (sense)
                    GetRegionPad(seq, featstr, lastpos + 1, lastpos + len, sense, false);
                else
                    GetRegionPad(seq, featstr, firstpos - len, firstpos - 1, sense, false);
            }

            if (before > 0)
            {
                if (sense)
                    GetRegionPad(seq, featstr, firstpos - before, firstpos - 1, sense, true);
                else
                    GetRegionPad(seq, featstr, lastpos + 1, lastpos + before, sense, true);
            }

            if (after > 0)
            {
                if (sense)
                    GetRegionPad(seq, featstr, lastpos + 1, lastpos + after, sense, false);
                else
                    GetRegionPad(seq, featstr, firstpos - after, firstpos - 1, sense, false);
            }

            return featstr.ToString();
        }

        // extractfeat_GetRegionPad
        private static void GetRegionPad(EmbossSeq seq, StringBuilder featstr, int start, int end, bool sense, bool beginning)
        {
            if (start > end)
                return;
            var result = new StringBuilder();
            char pad = seq.Type == "N" ? 'N' : 'X';
            int seqlen = seq.Seq.Length;
            if (start < 0)
            {
                result.Append(pad, -start);
                start = 0;
            }

            int tmp = end > seqlen - 1 ? seqlen - 1 : end;
            if (start <= seqlen && tmp >= 0)
                AppendSub(result, seq.Seq.ToString(), start, tmp);
            if (end > seqlen - 1)
                result.Append(pad, end - seqlen + 1);
            string text = result.ToString();
            if (!sense)
                text = ReverseComplement(text);
            if (beginning)
                featstr.Insert(0, text);
            else
                featstr.Append(text);
        }
    }

    // the strand of a feature as ajFeatGetStrand returns it: '.' for a protein, '-' or '+'
    private static char StrandOf(EmbossSeqFeature gf) => gf.Protein ? '.' : gf.Strand == '-' ? '-' : '+';

    // ajFeatGetType: the type name of the vocabulary of the sequence kind
    private static string TypeOf(EmbossSeqFeature gf) => EmbossFeatureVocab.TypeExternal((gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types, gf.Type);

    // ajStrAppendSubS: the positions are clipped into the string
    private static void AppendSub(StringBuilder sb, string s, long pos1, long pos2)
    {
        int ibegin = EmbossInputSequence.PositionStart(s.Length, 0, pos1);
        int iend = EmbossInputSequence.PositionStart(s.Length, ibegin, pos2);
        int ilen = iend - ibegin + 1;
        if (ilen > 0 && ibegin < s.Length)
            sb.Append(s, ibegin, Math.Min(ilen, s.Length - ibegin));
    }

    // ajSeqstrReverse: ajBaseAlphacharComp of every character, from the end
    internal static string ReverseComplement(string s)
    {
        const string fwd = "ACGTURYWSMKBDHVNXacgturywsmkbdhvnx", rev = "TGCAAYRWSKMVHDBNXtgcaayrwskmvhdbnx";
        return string.Create(s.Length, s, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
            {
                char c = src[src.Length - 1 - i];
                int k = fwd.IndexOf(c, StringComparison.Ordinal);
                dest[i] = k >= 0 ? rev[k] : c;
            }
        });
    }

    // ajFeatGetSeq: the residues of a feature; a feature between two residues has none; a remote part is read from another entry by the library
    private static string GetSeq(EmbossSeqFeature gf, EmbossSeq seq, EmbossRemoteSeq remote)
    {
        if (gf.Flags.HasFlag(EmbossLocationBits.BetweenSeq))
            return "";
        if (gf.Flags.HasFlag(EmbossLocationBits.RemoteId))
            return remote.Get(gf, seq);
        var tmp = new StringBuilder();
        AppendSub(tmp, seq.Seq.ToString(), gf.GetStart() - 1, gf.GetEnd() - 1);
        return gf.Strand == '-' ? ReverseComplement(tmp.ToString()) : tmp.ToString();
    }

    // ajFeatGetSeqJoin: the parts of a join one after the other; a join that is complemented as a whole is reverse-complemented at the end
    internal static string GetSeqJoin(EmbossSeqFeature gf, EmbossSeq seq, EmbossRemoteSeq remote)
    {
        if (gf.Subfeatures is not { Count: > 0 } subs)
            return GetSeq(gf, seq, remote);
        bool compjoin = gf.Flags.HasFlag(EmbossLocationBits.ComplementMain);
        var result = new StringBuilder();
        foreach (var part in subs)
        {
            if (part.Flags.HasFlag(EmbossLocationBits.BetweenSeq))
                continue;
            if (part.Flags.HasFlag(EmbossLocationBits.RemoteId))
            {
                result.Append(remote.Get(part, seq));
                continue;
            }

            var tmp = new StringBuilder();
            AppendSub(tmp, seq.Seq.ToString(), part.GetStart() - 1, part.GetEnd() - 1);
            result.Append(part.Strand == '-' && !compjoin ? ReverseComplement(tmp.ToString()) : tmp.ToString());
        }

        return compjoin ? ReverseComplement(result.ToString()) : result.ToString();
    }

    // embMiscMatchPatternDelimC: the text matches (any case, wildcards) one of the patterns of the list (separated by white space and the delimiters)
    private static bool MatchPatternDelim(string text, string pattern, string delims)
    {
        char[] sep = [.. " \t\n\r", .. delims];
        foreach (string key in pattern.Split(sep, StringSplitOptions.RemoveEmptyEntries))
        {
            if (EmbossWild.MatchCase(text, key))
                return true;
        }

        return false;
    }

    // extractfeat_MatchPatternTags
    private static bool MatchPatternTags(EmbossSeqFeature feat, string tpattern, string vpattern)
    {
        if (tpattern == "*" && vpattern == "*")
            return true;
        foreach (var (tagnam, tagval) in feat.Tags)
        {
            bool tval = MatchPatternDelim(tagnam, tpattern, ",;|");
            bool vval = tagval.Length == 0 ? vpattern == "*" : MatchPatternDelim(tagval, vpattern, ",;|");
            if (tval && vval)
                return true;
        }

        return false;
    }

    // extractfeat_MatchPatternDescribe: "(tag="value", tag2) "
    private static string MatchPatternDescribe(EmbossSeqFeature feat, string describe)
    {
        var sb = new StringBuilder();
        bool val = false;
        foreach (var (tagnam, tagval) in feat.Tags)
        {
            if (!MatchPatternDelim(tagnam, describe, ",;|"))
                continue;
            sb.Append(val ? ", " : "(");
            val = true;
            sb.Append(tagnam);
            if (tagval.Length > 0)
                sb.Append("=\"").Append(tagval).Append('"');
        }

        if (val)
            sb.Append(") ");
        return sb.ToString();
    }

    // extractfeat_MatchFeature (the flag of the tags of the last feature that is no join is kept for the parts of a join: tagsmatch)
    private static bool MatchFeature(EmbossSeqFeature gf, string source, string type, int sense, bool testscore, float minscore, float maxscore, string tag, string value, ref bool tagsmatch)
    {
        if (!gf.Flags.HasFlag(EmbossLocationBits.Multiple))
            tagsmatch = MatchPatternTags(gf, tag, value);
        if (!MatchPatternDelim(gf.Source, source, ",;|"))
            return false;
        if (type.Length > 0)
        {
            bool val = false;
            foreach (string key in type.Split([' ', '\t', '\n', '\r', ',', ';', '|'], StringSplitOptions.RemoveEmptyEntries))
            {
                if (TypeMatchWild(gf, key))
                {
                    val = true;
                    break;
                }
            }

            if (!val)
                return false;
        }

        char strand = StrandOf(gf);
        if (strand == '+' && sense == -1)
            return false;
        if (strand == '-' && sense == 1)
            return false;
        if (testscore && gf.Score < minscore)
            return false;
        if (testscore && gf.Score > maxscore)
            return false;
        return tagsmatch;
    }

    // ajFeatTypeMatchWildS
    internal static bool TypeMatchWild(EmbossSeqFeature gf, string str0)
    {
        string str = str0.ToLowerInvariant();
        var table = (gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types;
        bool ret = TableTypeTestWild(table, gf.Type, str);
        if (!ret)
        {
            if (str.Contains('*') || str.Contains('?'))
                str = str.TrimEnd('*');
            ret = TypeMatch(gf, table, str);
        }

        return ret;
    }

    /// <summary><c>ajFeatTypeMatchC</c>: the type of the feature is the name or one of the alternative names of the type.</summary>
    internal static bool TypeMatchC(EmbossSeqFeature gf, string txt)
    {
        var table = (gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types;
        return TypeMatch(gf, table, txt);
    }

    // ajFeatTypeMatchS
    private static bool TypeMatch(EmbossSeqFeature gf, EmbossFeatureVocab.Table table, string str)
    {
        if (gf.Type == str)
            return true;
        return InternalLimit(table, gf.Type) == InternalLimit(table, str) && InternalLimit(table, str) is not null;
    }

    // featTableTypeInternalLimit
    private static string? InternalLimit(EmbossFeatureVocab.Table table, string type)
    {
        string? retkey = table.FetchKey(type);
        string? savekey = null;
        for (int i = 0; retkey is not null && i < 5; i++)
        {
            savekey = retkey;
            string tmpstr = table.Fetch(retkey) ?? "";
            if (tmpstr.Length == 0 || tmpstr[^1] != ';')
            {
                retkey = table.FetchKey(tmpstr);
            }
            else
            {
                retkey = table.FetchKey(tmpstr);
                if (retkey == savekey)
                    retkey = null;
            }
        }

        return savekey;
    }

    // featTableTypeTestWild
    private static bool TableTypeTestWild(EmbossFeatureVocab.Table table, string type, string str)
    {
        if (EmbossWild.Match(type, str))
            return true;
        string? retkey = table.FetchKey(type);
        for (int i = 0; retkey is not null && i < 5; i++)
        {
            string tmpstr = table.Fetch(retkey) ?? "";
            if (tmpstr.Length == 0 || tmpstr[^1] != ';')
            {
                if (EmbossWild.Match(tmpstr, str))
                    return true;
                retkey = table.FetchKey(tmpstr);
            }
            else
            {
                int semi = tmpstr.IndexOf(';');
                string tmpkey = tmpstr[..Math.Max(semi, 0)];
                if (EmbossWild.Match(tmpkey, str))
                    return true;
                retkey = table.FetchKey(tmpstr);
            }
        }

        foreach (var (key, _) in table.Items)
        {
            if (EmbossWild.Match(key, str) && InternalLimit(table, key) == InternalLimit(table, type))
                return true;
        }

        return false;
    }
}
