// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The <c>seqset</c> input of an EMBOSS program (<c>acdSetSeqset</c> / <c>ajSeqsetRead</c> of <c>ajseqread.c</c>, <c>ajSeqsetFromList</c>, <c>ajSeqsetReverse</c>; B09 WP62, F281): every sequence of the USA is read into a new object with the range of the set,
/// the set takes the type of its first sequence, and with <c>-sreverse</c> a nucleic set is reversed (every member, <c>ajSeqReverseDo</c>).
/// </summary>
public sealed class EmbossSeqset
{
    /// <summary>The members, in the order they were read.</summary>
    public List<EmbossSeq> Seqs { get; } = [];

    /// <summary>The range of the set (<c>-sbegin</c> / <c>-send</c>, or those of the USA).</summary>
    public int Begin, End;

    /// <summary><c>ajSeqsetIsNuc</c>: the set's type is N, or the first member holds nucleic residues.</summary>
    public bool IsNuc => Seqs.Count > 0 && EmbossSeqType.IsNucSet(Seqs[0].Type, Seqs[0].SeqText);

    /// <summary>
    /// <c>acdSetSeqset</c>: reads the set; null when it cannot be read (the program dies with "Unable to read sequence", the messages of the library are in <paramref name="context"/>).
    /// </summary>
    public static EmbossSeqset? Read(string program, string spec, IReadOnlyDictionary<string, string> files, EmbossSequenceQualifiers q, EmbossFeatureReadContext context, string type)
    {
        var all = EmbossSeqall.Open(program, spec, files, q with { Reverse = false }, false, context, type, false);
        if (all is null)
            return null;
        var set = new EmbossSeqset();
        EmbossSeq? seq;
        while ((seq = EmbossSeqretPrograms.TryNextMember(all, context)) is not null)
            set.Seqs.Add(seq);
        if (q.Reverse && set.Seqs.Count > 0 && set.IsNuc)
        {
            foreach (var m in set.Seqs)
            {
                m.Rev = true;
                EmbossSeqRange.ReverseDo(m);
            }
        }

        return set;
    }

    /// <summary>
    /// <c>acdSetSeqsetall</c> / <c>ajSeqsetallRead</c>: the sets of the input: a new set is started when the reader of a format that holds several sets finishes one (<c>seqin-&gt;Multidone</c>), the others read the whole input as one set.
    /// Null when nothing can be read.
    /// </summary>
    public static List<EmbossSeqset>? ReadAll(string program, string spec, IReadOnlyDictionary<string, string> files, EmbossSequenceQualifiers q, EmbossFeatureReadContext context, string type, bool features)
    {
        var all = EmbossSeqall.Open(program, spec, files, q with { Reverse = false }, features, context, type, false);
        if (all is null)
            return null;
        var sets = new List<EmbossSeqset>();
        var current = new EmbossSeqset();
        EmbossSeq? seq;
        while ((seq = EmbossSeqretPrograms.TryNextMember(all, context, false)) is not null)
        {
            current.Seqs.Add(seq);
            if (all.MultiDone)
            {
                sets.Add(current);
                current = new EmbossSeqset();
            }
        }

        if (current.Seqs.Count > 0)
            sets.Add(current);
        if (q.Reverse && sets.Count > 0 && sets[0].IsNuc)
        {
            foreach (var set in sets)
            {
                foreach (var m in set.Seqs)
                {
                    m.Rev = true;
                    EmbossSeqRange.ReverseDo(m);
                }
            }
        }

        return sets;
    }
}
