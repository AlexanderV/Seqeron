// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.OpenReadingFrames, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.OpenReadingFrames;

using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Emboss.SequenceFormats.EmbossSequencePrograms;

/// <summary><c>checktrans</c> of EMBOSS 6.6.0 (<c>emboss/checktrans.c</c>): the STOP codons and the open reading frames of a protein translation (B09 WP63, F282).</summary>
public static partial class EmbossOrfPrograms
{
    /// <summary>
    /// <c>checktrans</c>: for every protein sequence the ORFs (the stretches between stops, at least <c>-orfml</c> long) in the report <c>-outfile</c>, as sequences <c>name_N</c> in <c>-outseq</c> and as the protein feature table of <c>-outfeat</c>;
    /// <c>-addlast</c> forces an end with a stop (the stop that was added is not counted: the flag of the program is never reset, so it holds for the sequences after the first one that got a stop).
    /// </summary>
    internal static EmbossRunResult Checktrans(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("checktrans");
        var all = EmbossProgramBase.Open("checktrans", a, files, ctx, "stopprotein");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        var first = all.Peek();
        int orfml;
        bool addlast;
        string outfile;
        EmbossOutseqOptions outOptions;
        string baseName = EmbossOutPath.DefaultBaseName(first.Name, a.Raw("sequence") ?? "");
        try
        {
            orfml = a.IntClamped("orfml", 100, 1, int.MaxValue, ctx);
            addlast = a.Bool("addlast", true);
            outfile = a.Raw("outfile") ?? baseName + ".checktrans";
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
            string odir = EmbossOutPath.FixDirectory(a.RawAssoc("outfile", "odirectory"));
            if (odir.Length > 0 && !outfile.Contains('/', StringComparison.Ordinal))
                outfile = odir + outfile;                                               // ajFilenameReplacePathS: the directory is put in front of a name that has none (also of stdout)
            if (!EmbossOutPath.DirectoryExists(EmbossReportArgs.ExistingDirectories(files), outfile))
                throw a.Die("outfile", $"Warning: Failed to open output file '{outfile}' error:2 'No such file or directory'", $"Error: Unable to open file '{outfile}' for output");
            outOptions = a.Seqout("outseq");
            if (a.Raw("outseq") is null)
                outOptions = outOptions with { OutSeq = (a.RawAssoc("outseq", "osname") is { Length: > 0 } on ? EmbossOutPath.FmtWord(on) : baseName) + "." + (a.RawAssoc("outseq", "osextension") is { Length: > 0 } oe ? oe : a.RawAssoc("outseq", "osformat") ?? "fasta") };
            outOptions = outOptions with { ExistingDirectories = EmbossReportArgs.ExistingDirectories(files) };
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("checktrans", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResults.From(stream.DiedOutput()) with { Files = outfile == "stdout" ? new Dictionary<string, string>() : new Dictionary<string, string> { [outfile] = "" } };
        var fo = EmbossFeatout.Open(a, ctx, "outfeat", "checktrans", baseName, files, EmbossRunEnvironment.Current.RunDate ?? DateTime.Now);
        if (fo is null)
        {
            // the files the items before it opened exist (empty)
            var died = EmbossProgramBase.Died(ctx, "outfeat");
            var opened = new Dictionary<string, string>();
            foreach (string f in new[] { outfile, outOptions.Single ? "" : EmbossOutPath.Combine(EmbossOutPath.FixDirectory(outOptions.Directory), ChecktransSeqName(outOptions.OutSeq)) })
                if (f.Length > 0 && f != "stdout")
                    opened[f] = "";
            return died with { Files = opened };
        }

        var log = new OutLog();
        bool addedAsterisk = false;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            int begin = EmbossSeqRange.GetBegin(seq), end = EmbossSeqRange.GetEnd(seq);
            string full = seq.SeqText;
            string sub = EmbossProgramBase.SubS(full, begin - 1, end - 1);
            // ajSeqGetSeqC(seq)[end-1]: the character at the end of the range
            if (end - 1 < 0 || end - 1 >= full.Length)
                throw new InvalidOperationException("checktrans reads the character at the end of an empty range: the reference program has no defined result.");
            if (addlast && full[end - 1] != '*')
            {
                sub += "*";
                addedAsterisk = true;
            }

            sub = sub.ToUpperInvariant();
            int len = sub.Length;
            log.Add(0, "\n\nCHECKTRANS of " + seq.Name + " from " + begin.ToString(CultureInfo.InvariantCulture) + " to " + (begin + len - 1).ToString(CultureInfo.InvariantCulture) + "\n\n");
            ChecktransFindOrfs(stream, log, sub, seq.Name, orfml, addedAsterisk);
            ChecktransFeatures(fo, log, sub, seq.Name, orfml);
        }

        var closed = stream.Close();
        log.Add(1, stream.Take());
        fo.Close();
        // the destinations: the report, the sequences and the features go to their files; what goes to the standard output keeps the order it was written in
        var outFiles = new Dictionary<string, string>();
        foreach (var (path, text) in closed.SingleFiles)
            outFiles[path] = text;                                                        // -ossingle: a file for every sequence
        var stdout = new StringBuilder();
        string[] dest = [outfile, outOptions.Single ? "" : EmbossOutPath.Combine(EmbossOutPath.FixDirectory(outOptions.Directory), ChecktransSeqName(outOptions.OutSeq)), fo.File];
        foreach (var (target, text) in log.Segments)
        {
            string d = dest[target];
            if (text.Length == 0)
                continue;
            if (d == "stdout")
                stdout.Append(text);
            else
                outFiles[d] = outFiles.GetValueOrDefault(d, "") + text;
        }

        // the files every output opened exist (empty) even when nothing was written to them
        foreach (string d in dest)
            if (d.Length > 0 && d != "stdout" && !outFiles.ContainsKey(d))
                outFiles[d] = "";
        return new EmbossRunResult(stdout.ToString(), EmbossRunResult.Join(ctx.Messages), 0, outFiles);
    }

    // the name of the file of a sequence output: the value without the format prefix
    private static string ChecktransSeqName(string outseq)
    {
        var m = System.Text.RegularExpressions.Regex.Match(outseq, "^([A-Za-z0-9-]*)::?(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant | System.Text.RegularExpressions.RegexOptions.Singleline);
        return m.Success ? m.Groups[2].Value : outseq;
    }

    private static void ChecktransFindOrfs(EmbossSeqoutStream stream, OutLog log, string p, string name, int minOrf, bool addedAsterisk)
    {
        int count = 1, lastStop = 0, len = p.Length;
        log.Add(0, "\tORF#\tPos\tLen\tORF Range\tSequence name\n\n");
        for (int i = 0; i < len; ++i)
        {
            if (p[i] == '*')
            {
                int orflength = i - lastStop;
                if (orflength >= minOrf)
                {
                    log.Add(0, string.Create(CultureInfo.InvariantCulture, $"\t{count}\t{i + 1}\t{orflength}\t{i - orflength + 1}-{i}\t{name}_{count}\n"));
                    var s = EmbossProgramBase.NewSeq(name + "_" + count.ToString(CultureInfo.InvariantCulture), "", p.Substring(i - orflength, orflength), "");
                    stream.Write(s);
                    log.Add(1, stream.Take());
                }

                lastStop = ++i;
                ++count;
                while (i < len && p[i] == '*')
                {
                    lastStop = ++i;
                    ++count;
                }
            }
        }

        if (addedAsterisk)
            --count;
        log.Add(0, "\n\tTotal STOPS: " + (count - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + "\n\n ");
    }

    private static void ChecktransFeatures(EmbossFeatout fo, OutLog log, string p, string name, int minOrf)
    {
        int lastStop = 0, len = p.Length;
        var features = new List<EmbossFeature>();
        for (int i = 0; i < len; ++i)
        {
            if (p[i] == '*')
            {
                int orflength = i - lastStop;
                if (orflength >= minOrf)
                    features.Add(new EmbossFeature("polypeptide_region", "checktrans", i - orflength + 1, i, 0.0, '\0', 0, []));
                lastStop = ++i;
                while (i < len && p[i] == '*')
                    lastStop = ++i;
            }
        }

        features.Sort((x, y) => x.Start.CompareTo(y.Start));            // ajFeatSortByStart
        var seq = new EmbossReportSequence(name, name, "", "", "", 1, 0, 0, false, false);
        int length = features.Count == 0 ? 0 : features.Max(f => Math.Max(f.Start, f.End));
        string before = fo.Text;
        fo.Write(new EmbossReport(seq, features, true, "P", TableLength: length, WriteWhenEmpty: true));
        log.Add(2, fo.Text[before.Length..]);
    }
}
