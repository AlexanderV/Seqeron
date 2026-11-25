// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.PairwiseAlignment;

using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Emboss.SequenceFormats.EmbossSequencePrograms;

/// <summary><c>wordfinder</c> of EMBOSS 6.6.0 (<c>emboss/wordfinder.c</c>): the local alignment of large sequences found by their word matches (B09 WP63, F282).</summary>
public static partial class EmbossPairwisePrograms
{
    // acdSetFloat with the attributes minimum / maximum: a value outside the limits is reset to the limit with the warning "floating point value out of range V less than (reset to) L"
    private static float FloatClamped(EmbossAcdArgs a, string name, float dflt, float min, float max, EmbossFeatureReadContext ctx)
    {
        float v = a.Float(name, dflt);
        if (v < min)
        {
            ctx.Warn($"floating point value out of range {CFormat.Fixed(v, 3)} less than (reset to) {CFormat.Fixed(min, 3)}");
            v = min;
        }

        if (v > max)
        {
            ctx.Warn($"floating point value out of range {CFormat.Fixed(v, 3)} more than (reset to) {CFormat.Fixed(max, 3)}");
            v = max;
        }

        return v;
    }

    // the matrix of -datafile: an EMBOSS data matrix (E...), or the text of a matrix file of the work directory
    private static SubstitutionMatrix? ReadMatrix(string name, IReadOnlyDictionary<string, string> files)
    {
        if (files.TryGetValue(name, out string? text))
            return ParseMatrixFile(name, text);
        return name.StartsWith('E') && SubstitutionMatrix.TryGet(name, out var m) ? m : null;
    }

    // ajMatrixfRead: '#' lines are comments, the first other line holds the symbols, every following line a symbol and its scores
    private static SubstitutionMatrix? ParseMatrixFile(string name, string text)
    {
        var lines = text.Split('\n').Select(l => l.TrimEnd('\r')).Where(l => l.Length > 0 && l[0] != '#').ToList();
        if (lines.Count < 2)
            return null;
        string[] header = lines[0].Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
        var scores = new Dictionary<(char, char), int>();
        foreach (string line in lines.Skip(1))
        {
            string[] t = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
            if (t.Length != header.Length + 1 || t[0].Length != 1)
                return null;
            for (int i = 0; i < header.Length; i++)
            {
                if (!int.TryParse(t[i + 1], NumberStyles.Integer, CultureInfo.InvariantCulture, out int v))
                    return null;
                scores[(char.ToUpperInvariant(t[0][0]), char.ToUpperInvariant(header[i][0]))] = v;
            }
        }

        string alphabet = string.Concat(header.Select(h => char.ToUpperInvariant(h[0])));
        var matrix = new int[alphabet.Length, alphabet.Length];
        for (int i = 0; i < alphabet.Length; i++)
            for (int j = 0; j < alphabet.Length; j++)
                matrix[i, j] = scores.GetValueOrDefault((alphabet[i], alphabet[j]));
        return new SubstitutionMatrix(name, alphabet, matrix);
    }

    private sealed class WordGroup
    {
        public int Offset, Count, Total;
    }

    /// <summary>
    /// <c>wordfinder</c>: every sequence of <c>-bsequence</c> is indexed by its words of <c>-wordlen</c> residues (<c>embWordGetTable</c>); every sequence of <c>-asequence</c> is looked up in it (<c>embWordBuildMatchTable</c>), the diagonal with the most word residues
    /// fixes the start points, and the banded Smith-Waterman of <c>-width</c> cells (<c>embAlignPathCalcSWFast</c>) along it gives the alignment, written in <c>-aformat</c> unless the limits exclude it (the lines of <c>-errfile</c>).
    /// </summary>
    internal static EmbossRunResult Wordfinder(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("wordfinder");
        ctx.InputType ??= "any";
        var qrySet = EmbossSeqset.Read("wordfinder", a.Raw("asequence") ?? "", files, a.Sequence("asequence"), ctx, "any");
        if (qrySet is null)
            return EmbossProgramBase.Died(ctx, "asequence");
        bool protein = !qrySet.IsNuc;
        var trgAll = EmbossProgramBase.Open("wordfinder", a, files, ctx, protein ? "stopprotein" : "nucleotide", "bsequence");
        if (trgAll is null)
            return EmbossProgramBase.Died(ctx, "bsequence");
        SubstitutionMatrix matrix;
        float gapopen, gapextend;
        int width, wordlen, limitmatch, limitalign, lowmatch, lowalign;
        EmbossAlignArgs.Opened? align;
        string errfile;
        string baseName = EmbossOutPath.DefaultBaseName(qrySet.Seqs.Count > 0 ? qrySet.Seqs[0].Name : "", a.Raw("asequence") ?? "");
        try
        {
            string matrixName = a.Raw("datafile") ?? (protein ? "EBLOSUM62" : "EDNAFULL");
            var read = ReadMatrix(matrixName, files);
            if (read is null)
            {
                ctx.Error($"Unable to read matrix '{matrixName}'");
                throw a.Die("datafile");
            }

            matrix = read;
            gapopen = FloatClamped(a, "gapopen", 30.0f, 0.0f, 1000.0f, ctx);
            gapextend = FloatClamped(a, "gapextend", 1.5f, 0.0f, 10.0f, ctx);
            width = a.IntClamped("width", 16, 1, int.MaxValue, ctx);
            wordlen = a.IntClamped("wordlen", 6, 3, int.MaxValue, ctx);
            limitmatch = a.IntClamped("limitmatch", 0, 0, int.MaxValue, ctx);
            limitalign = a.IntClamped("limitalign", 0, 0, int.MaxValue, ctx);
            lowmatch = a.IntClamped("lowmatch", 0, 0, int.MaxValue, ctx);
            lowalign = a.IntClamped("lowalign", 0, 0, int.MaxValue, ctx);
            string fmtName = a.RawAssoc("outfile", "aformat") ?? "simple";
            bool plain = EmbossAlignmentWriter.TryParseFormat(fmtName, out var parsed) && parsed is EmbossAlignFormat.Match or EmbossAlignFormat.Score;       // the match and score writers have no matrix and no gap penalties (ajAlignDefineCC)
            align = EmbossAlignArgs.Open(a, ctx, "outfile", "wordfinder", "simple", baseName, EmbossRunEnvironment.Current.RunDate ?? DateTime.Now, files, plain ? null : matrix, plain ? null : (double)gapopen, plain ? null : (double)gapextend);
            if (align is null)
                return EmbossProgramBase.Died(ctx, "outfile");
            errfile = a.Raw("errfile") ?? "wordfinder.error";
            string edir = EmbossOutPath.FixDirectory(a.RawAssoc("errfile", "odirectory"));
            if (edir.Length > 0 && !errfile.Contains('/', StringComparison.Ordinal))
                errfile = edir + errfile;
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        gapopen = EmbossPairwisePrograms.RoundFloat(gapopen, 8);
        gapextend = EmbossPairwisePrograms.RoundFloat(gapextend, 8);
        var table = new EmbossScoreTable(matrix);
        foreach (var q in qrySet.Seqs)
            EmbossSeqRange.Trim(q);
        var log = new OutLog();
        bool showsSequences = align.Format is not (EmbossAlignFormat.Match or EmbossAlignFormat.Score);
        int itrg = 0;
        while (EmbossSeqretPrograms.TryNext(trgAll, ctx, out var trg))
        {
            itrg++;
            int imatches = 0;
            EmbossSeqRange.Trim(trg);
            string ptrg = trg.SeqText;
            int trglen = ptrg.Length;
            // embWordGetTable: a sequence shorter than the word has no table
            if (trglen >= wordlen)
            {
                var matcher = new EmbossRepeatPrograms.PairWordMatcher(ptrg, wordlen, EmbossSeqType.IsNuc(trg.Type, ptrg));
                if (matcher.HasTable)
                {
                    foreach (var qry in qrySet.Seqs)
                    {
                        string qqry = qry.SeqText;
                        int qrylen = qqry.Length;
                        int matchscore = WordfinderStartPoints(matcher, qqry, wordlen, trglen, qrylen, out int qrystart, out int trgstart, out int qryend, out int trgend, ctx, qry);
                        if (matchscore == 0 || (limitmatch != 0 && matchscore > limitmatch) || matchscore < lowmatch)
                        {
                            if (matchscore != 0)
                                log.Add(1, string.Create(CultureInfo.InvariantCulture, $"Match limits ({lowmatch}..{limitmatch}) exclude '{qry.Name}' '{trg.Name}' {matchscore}\n"));
                            continue;
                        }

                        imatches++;
                        int lenq = qryend - qrystart + 1, lent = trgend - trgstart + 1;
                        var path = new float[(long)(lent + 1) * width];
                        var compass = new int[path.Length];
                        string subq = qqry.Substring(qrystart, lenq), subt = ptrg.Substring(trgstart, lent);
                        float score = EmbossPairwisePrograms.SwPathCalcFast(subq, subt, lenq, lent, 0, width, gapopen, gapextend, path, compass, table);
                        var (nqry, mtrg) = EmbossPairwisePrograms.SwWalkFast(path, compass, gapopen, gapextend, qqry, ptrg, lenq, lent, 0, width, ref qrystart, ref trgstart);
                        EmbossPairAlignment pair;
                        int alignLength;
                        if (!showsSequences)
                        {
                            pair = new EmbossPairAlignment(new EmbossAlignedSequence(trg.Name, mtrg.Replace('.', '-'), Usa: "<null>::<null>"), new EmbossAlignedSequence(qry.Name, nqry.Replace('.', '-'), Usa: "<null>::<null>"), score);
                            alignLength = Math.Max(mtrg.Length, nqry.Length);
                        }
                        else
                        {
                            pair = WordfinderReportLocal(qry, trg, nqry, mtrg, qrystart, trgstart, score);
                            alignLength = Math.Max(mtrg.Length, nqry.Length);
                        }

                        if ((limitalign != 0 && alignLength > limitalign) || alignLength < lowalign)
                        {
                            imatches--;
                            log.Add(1, string.Create(CultureInfo.InvariantCulture, $"Align limits ({lowalign}..{limitalign}) excludes '{qry.Name}' '{trg.Name}' {alignLength}\n"));
                        }
                        else
                        {
                            if (showsSequences)
                                pair = pair with { ExtraSubHeader = string.Create(CultureInfo.InvariantCulture, $"\nWordscore:{matchscore}\nAlignlength:{alignLength}") };
                            align.Writer.Write(pair);
                            log.Add(0, align.Writer.Take());
                        }
                    }
                }
            }

            if (imatches > 0)
                log.Add(1, string.Create(CultureInfo.InvariantCulture, $"Target {itrg} {trg.Name} matches {imatches}\n"));
        }

        align.Writer.Close();
        log.Add(0, align.Writer.Take());
        var outFiles = new Dictionary<string, string>();
        var stdout = new StringBuilder();
        string[] dest = [align.File, errfile];
        foreach (var (target, text) in log.Segments)
        {
            if (text.Length == 0)
                continue;
            if (dest[target] == "stdout")
                stdout.Append(text);
            else
                outFiles[dest[target]] = outFiles.GetValueOrDefault(dest[target], "") + text;
        }

        foreach (string d in dest)
            if (d != "stdout" && !outFiles.ContainsKey(d))
                outFiles[d] = "";                                                         // an output that was opened exists even when nothing was written to it
        return new EmbossRunResult(stdout.ToString(), EmbossRunResult.Join(ctx.Messages), 0, outFiles);
    }

    // embAlignReportLocal: the rows with '.' for gaps become '-', the residues outside the aligned part are the offsets (before) and offends (after), added to those of the range
    private static EmbossPairAlignment WordfinderReportLocal(EmbossSeq first, EmbossSeq second, string rowFirst, string rowSecond, int startFirst, int startSecond, float score)
    {
        static int Gaps(string row) => row.Count(c => c is '.' or '-' or ' ');
        int offend1 = first.Seq.Length - rowFirst.Length - startFirst + Gaps(rowFirst);
        int offend2 = second.Seq.Length - rowSecond.Length - startSecond + Gaps(rowSecond);
        // ajSeqNewRangeC(row, start + Offset, offend + Offend, reversed): a reversed sequence takes the two offsets the other way round
        var a = new EmbossAlignedSequence(first.Name, rowFirst.Replace('.', '-'), first.Reversed ? offend1 + first.Offend : startFirst + first.Offset, first.Reversed ? startFirst + first.Offset : offend1 + first.Offend, first.Reversed, first.Acc, first.Desc, UsaOf(first));
        var b = new EmbossAlignedSequence(second.Name, rowSecond.Replace('.', '-'), second.Reversed ? offend2 + second.Offend : startSecond + second.Offset, second.Reversed ? startSecond + second.Offset : offend2 + second.Offend, second.Reversed, second.Acc, second.Desc, UsaOf(second));
        return new EmbossPairAlignment(a, b, score);
    }

    // wordfinder_findstartpoints: the diagonal with the most word residues (the first one of equal totals in the order of the sorted match list)
    private static int WordfinderStartPoints(EmbossRepeatPrograms.PairWordMatcher matcher, string qqry, int wordlen, int trglen, int qrylen, out int qrystart, out int trgstart, out int qryend, out int trgend, EmbossFeatureReadContext ctx, EmbossSeq qry)
    {
        qrystart = trgstart = qryend = trgend = 0;
        if (qrylen < wordlen)
        {
            ctx.Warn($"ERROR: Sequence {UsaOf(qry)} length {qrylen.ToString(CultureInfo.InvariantCulture)} less than word length {wordlen.ToString(CultureInfo.InvariantCulture)}");
            return 0;
        }

        // wordOrderMatchTable: by length descending, then the start in the target and in the query
        var matches = matcher.Matches(qqry).OrderByDescending(m => m.Length).ThenBy(m => m.PositionInText).ThenBy(m => m.PositionInQuery).ToList();
        if (matches.Count == 0)
            return 0;
        var groups = new List<WordGroup>();
        foreach (var m in matches)
        {
            int offset = m.PositionInText - m.PositionInQuery;
            var g = groups.Find(x => x.Offset == offset);
            if (g is null)
                groups.Add(new WordGroup { Offset = offset, Count = 1, Total = m.Length });
            else
            {
                g.Total += m.Length;
                g.Count++;
            }
        }

        int max = -10;
        WordGroup? best = null;
        foreach (var g in groups)
        {
            if (g.Total > max)
            {
                max = g.Total;
                best = g;
            }
        }

        int off = best!.Offset;
        if (off > 0)
        {
            trgstart = off;
            qrystart = 0;
        }
        else
        {
            qrystart = 0 - off;
            trgstart = 0;
        }

        trgend = trgstart;
        qryend = qrystart;
        int trgmax = trglen - 1, qrymax = qrylen - 1;
        while (trgend < trgmax && qryend < qrymax)
        {
            trgend++;
            qryend++;
        }

        return max;
    }
}
