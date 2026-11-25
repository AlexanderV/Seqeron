// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.OpenReadingFrames, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.OpenReadingFrames;

using Seqeron.Genomics.Analysis;

/// <summary>A run of <c>getorf</c> (B09 WP61, F280).</summary>
/// <param name="Sequence">The value of <c>-sequence</c>: a file name with an optional <c>format::</c> prefix, <c>asis::SEQUENCE</c>, a list file or a query.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers (<c>-sbegin1</c> ...).</param>
/// <param name="Table"><c>-table</c>: the genetic code file (0..6, 9..16, 21..23).</param>
/// <param name="MinSize"><c>-minsize</c>: the minimum size of an ORF in nucleotides.</param>
/// <param name="MaxSize"><c>-maxsize</c>: the maximum size of an ORF in nucleotides.</param>
/// <param name="Find"><c>-find</c>: the type of output (0..6, the <see cref="OrfFindMode"/> values).</param>
/// <param name="Methionine"><c>-methionine</c>: a START codon at the beginning of a protein is an M.</param>
/// <param name="Circular"><c>-circular</c>: the sequence is circular.</param>
/// <param name="Reverse"><c>-reverse</c>: ORFs of the reverse complement are found as well.</param>
/// <param name="Flanking"><c>-flanking</c>: the number of nucleotides around a codon for the modes 4 .. 6.</param>
/// <param name="Output">The qualifiers of the sequence output <c>-outseq</c>.</param>
public sealed record EmbossGetorfRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, int Table = 0, int MinSize = 30, int MaxSize = 1000000, int Find = 0,
    bool Methionine = true, bool Circular = false, bool Reverse = true, int Flanking = 100, EmbossOutseqOptions? Output = null);

/// <summary>
/// <c>getorf</c> of EMBOSS 6.6.0 (<c>emboss/getorf.c</c>): the open reading frames of the sequences of a file. The scan is <see cref="Translator.FindOrfsGetorfProgram"/> (the model of <c>getorf_FindORFs</c> over the
/// EMBOSS genetic code files); this class is the program around it: the sequences with their range, the names <c>name_N</c>, the header <c>[start - end]</c> with the offset of the range, the
/// <c>(REVERSE SENSE)</c> and <c>(ORF crosses the breakpoint)</c> remarks, <c>Around codon at N.</c>, and the sequence output.
/// </summary>
public static class EmbossGetorfProgram
{
    private static readonly string[] TableCodes = ["0", "1", "2", "3", "4", "5", "6", "9", "10", "11", "12", "13", "14", "15", "16", "21", "22", "23"];
    private static readonly string[] FindCodes = ["0", "1", "2", "3", "4", "5", "6"];

    /// <summary>Runs <c>getorf</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossSeqretOutput Run(EmbossGetorfRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "getorf" };
        var seqall = EmbossSeqall.Open("getorf", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, false, context, "dna");
        if (seqall is null)
            return EmbossSeqretPrograms.Died(context);
        return Execute(context, seqall, r);
    }

    internal static EmbossRunResult RunArgs(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var context = new EmbossFeatureReadContext { Program = "getorf" };
        var seqall = EmbossSeqall.Open("getorf", a.Raw("sequence") ?? "", files, a.Sequence("sequence"), false, context, "dna");
        if (seqall is null)
            return EmbossRunResults.From(EmbossSeqretPrograms.Died(context));
        EmbossGetorfRequest r;
        try
        {
            r = new EmbossGetorfRequest("", files, null,
                int.Parse(a.Menu("table", "0", TableCodes), CultureInfo.InvariantCulture), a.Int("minsize", 30), a.Int("maxsize", 1000000), int.Parse(a.Menu("find", "0", FindCodes), CultureInfo.InvariantCulture),
                a.Bool("methionine", true), a.Bool("circular", false), a.Bool("reverse", true), a.Int("flanking", 100), a.Seqout("outseq"));
        }
        catch (EmbossAcdDieException e)
        {
            context.Messages.AddRange(e.Messages);
            return EmbossRunResults.From(new EmbossSeqretOutput("", "", 0, [.. context.Messages]) { Died = true });
        }

        return EmbossRunResults.From(Execute(context, seqall, r));
    }

    private static EmbossSeqretOutput Execute(EmbossFeatureReadContext context, EmbossSeqall seqall, EmbossGetorfRequest r)
    {
        var stream = new EmbossSeqoutStream("getorf", context, r.Output ?? new EmbossOutseqOptions(), "stopprotein");
        if (stream.OpenFailed)
            return stream.DiedOutput();
        var mode = (OrfFindMode)r.Find;
        bool nucleotideOutput = r.Find is 2 or 3 or 4 or 5 or 6;
        while (EmbossSeqretPrograms.TryNext(seqall, context, out var seq))
        {
            EmbossSeqRange.Trim(seq);                                  // ajSeqTrim
            int len = seq.Seq.Length;
            string residues = seq.SeqText;
            if (r.Circular)
                seq.Offset = seq.Offend = 0;                           // ajSeqAssignSeqS of the tripled sequence resets the offsets
            int orfNo = 1;
            foreach (var orf in EmbossGetorf.FindOrfsGetorfProgram(residues, mode, r.Table, r.MinSize, r.MaxSize, r.Flanking, r.Reverse, r.Methionine, r.Circular))
            {
                var (s, e) = orf.ToGetOrfRange(len);
                // ajSeqGetOffset: the reverse pass (ajSeqReverseForce) swapped the offsets of the two ends of the range
                int offset = orf.IsReverseStrand ? seq.Offend : seq.Offset;
                string header = "[" + (s + offset).ToString(CultureInfo.InvariantCulture) + " - " + (e + offset).ToString(CultureInfo.InvariantCulture) + "] ";
                if (orf.IsReverseStrand)
                    header += "(REVERSE SENSE) ";
                if (s > len || e > len)
                    header += "(ORF crosses the breakpoint) ";
                if (mode is OrfFindMode.FlankingStart or OrfFindMode.FlankingInitialStop or OrfFindMode.FlankingEndingStop)
                    header += "Around codon at " + (orf.CodonPosition + 1).ToString(CultureInfo.InvariantCulture) + ". ";
                var o = new EmbossSeq { Name = seq.Name + "_" + (orfNo++).ToString(CultureInfo.InvariantCulture), Desc = header + seq.Desc, Type = nucleotideOutput ? "N" : "P" };
                o.Seq.Append(orf.Sequence);
                stream.Write(o);
            }

            if (r.Reverse)
                (seq.Offset, seq.Offend) = (seq.Offend, seq.Offset);   // ajSeqReverseForce
        }

        context.Fatal |= false;
        return stream.Close();
    }
}
