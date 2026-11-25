// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.FeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.FeatureTools;

using Seqeron.Genomics.Analysis;

/// <summary>A run of <c>maskfeat</c> (A5-8c3g, B09 F249).</summary>
/// <param name="Sequence">The value of <c>-sequence</c>: a file name with an optional <c>format::</c> prefix, or <c>asis::SEQUENCE</c>.</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers (<c>-sformat1</c>, <c>-slower1</c> ...).</param>
/// <param name="Type">The feature types to mask (<c>-type</c>, default <c>repeat*</c>): names separated by blanks, commas, semicolons or bars, <c>*</c> and <c>?</c> wildcards.</param>
/// <param name="ToLower"><c>-tolower</c>: the masked regions are made lower case.</param>
/// <param name="MaskChar"><c>-maskchar</c>: the masking character (default X for protein, N for nucleic); no character or a blank means lower case.</param>
public sealed record EmbossMaskfeatRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, string Type = "repeat*", bool ToLower = false, string? MaskChar = null);

/// <summary>
/// <c>maskfeat</c> of EMBOSS 6.6.0 (<c>emboss/maskfeat.c</c>): the regions of the features of the wanted types are masked with a character or made lower case, and the sequences are written.
/// </summary>
public static class EmbossMaskfeatPrograms
{
    private const string WhiteSpace = " \t\n\r,;|";

    /// <summary>Runs <c>maskfeat</c>; the sequence output is FASTA.</summary>
    /// <exception cref="ArgumentException">The mask character is longer than one character ("Too long (maskchar)").</exception>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossSeqretOutput Maskfeat(EmbossMaskfeatRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "maskfeat" };
        var seqall = EmbossSeqall.Open("maskfeat", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, true, context, "any");
        if (seqall is null)
            return EmbossSeqretPrograms.Died(context, "sequence");
        var first = seqall.Peek();
        string defaultMask = first.Type == "P" ? "X" : "N";                   // acdprotein: the first sequence is not nucleotide
        string maskchar = r.MaskChar ?? defaultMask;
        if (maskchar.Length > 1)
        {
            context.Error("Too long (maskchar) - maximum length is 1 characters");
            return EmbossSeqretPrograms.Died(context, "maskchar");
        }

        bool lower = r.ToLower || maskchar.Length == 0 || maskchar == " ";
        var sb = new StringBuilder();
        int unnamed = 0, written = 0;
        EmbossSeq? seq;
        while ((seq = seqall.Next()) is not null)
        {
            Mask(seq, r.Type, maskchar, lower);
            EmbossSeqretPrograms.WriteFasta(sb, seq, ref unnamed, "any");
            written++;
        }

        return new EmbossSeqretOutput(sb.ToString(), "", written, [.. context.Messages]);
    }

    // maskfeat_FeatSeqMask
    private static void Mask(EmbossSeq seq, string type, string maskchar, bool lower)
    {
        var chars = new StringBuilder(seq.SeqText);
        if (seq.Fttable is { } table && table.Features.Count > 0)
        {
            foreach (var gf in table.Features)
            {
                var tokens = new StrTok(type, WhiteSpace);
                while (tokens.Next(out string key))
                {
                    string typeName = EmbossFeatureVocab.TypeExternal((gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types, gf.Type);
                    if (!MatchWild(typeName, key))
                        continue;
                    if (lower)
                        ToLower(chars, unchecked((int)(uint)(gf.GetStart() - 1)), unchecked((int)(uint)(gf.GetEnd() - 1)));
                    else
                        MaskRange(chars, (long)unchecked((uint)(gf.GetStart() - 1)), (long)unchecked((uint)(gf.GetEnd() - 1)), maskchar[0]);
                }
            }
        }

        seq.Seq.Clear().Append(chars);
    }

    // ajCvtSposToPosStart
    private static int PosStart(int length, int minimum, long position) => EmbossInputSequence.PositionStart(length, minimum, position);

    // ajStrMaskRange
    private static void MaskRange(StringBuilder s, long pos1, long pos2, char c)
    {
        int ibegin = PosStart(s.Length, 0, pos1);
        int iend = PosStart(s.Length, ibegin, pos2);
        if (iend < ibegin || s.Length == 0)
            return;
        for (int i = ibegin; i <= iend && i < s.Length; i++)
            s[i] = c;
    }

    // maskfeat_StrToLower: the range is cut out and its lower-case copy inserted at the begin position
    private static void ToLower(StringBuilder s, int begin, int end)
    {
        if (s.Length == 0)
            return;
        int ib = PosStart(s.Length, 0, begin);
        int ie = PosStart(s.Length, ib, end);
        string sub = s.ToString(ib, ie - ib + 1);
        string lowered = new(sub.Select(ch => char.IsAsciiLetterUpper(ch) ? (char)(ch + 32) : ch).ToArray());
        int cb = PosStart(s.Length, 0, begin);
        int ce = PosStart(s.Length, cb, end) + 1;
        s.Remove(cb, ce - cb);
        int pos = PosStart(s.Length + 1, 0, begin);
        s.Insert(pos, lowered);
    }

    // ajCharMatchWildC (ajCharCmpWild): '*' matches any run, '?' one character
    internal static bool MatchWild(string text, string pattern) => CmpWild(text, 0, pattern, 0) == 0;

    private static int CmpWild(string txt, int q, string pat, int p)
    {
        if (p >= pat.Length && q >= txt.Length)
            return 0;
        if (p >= pat.Length)
            return -1;
        while (p < pat.Length)
        {
            if (q >= txt.Length && pat[p] != '*')
                return 1;
            switch (pat[p])
            {
                case '?':
                    p++;
                    q++;
                    break;
                case '*':
                    p++;
                    if (p >= pat.Length)
                        return 0;
                    if (q >= txt.Length)
                        return CmpWild(txt, q, pat, p);
                    while (q < txt.Length)
                    {
                        if (CmpWild(txt, q, pat, p) == 0)
                            return 0;
                        q++;
                    }

                    return 1;
                default:
                    if (pat[p] != txt[q])
                        return pat[p] > txt[q] ? -1 : 1;
                    p++;
                    if (q < txt.Length)
                        q++;
                    break;
            }
        }

        if (p < pat.Length)
            return -1;
        if (q < txt.Length)
            return 1;
        return 0;
    }
}
