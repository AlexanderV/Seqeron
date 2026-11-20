using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>Result of <see cref="GenomicAnalyzer.MegaMerge"/>: the merged sequence, the program's report text and the minimal set of word matches the merge follows.</summary>
/// <param name="MergedSequence">The merged sequence (matched regions lower case, regions taken from one sequence over a mismatch or an unmatched end upper case; empty when nothing was merged).</param>
/// <param name="Report">The text of the <c>megamerger</c> report file (empty when nothing was merged).</param>
/// <param name="Matches">The minimal non-overlapping matches (<see cref="WordMatch.Sequence1"/> = 0, <see cref="WordMatch.Sequence2"/> = 1), ordered by position in the first sequence.</param>
public sealed record MegaMergeResult(string MergedSequence, string Report, IReadOnlyList<WordMatch> Matches);

public static partial class GenomicAnalyzer
{
    #region megamerger (merge of two large sequences by word-match anchors)

    internal sealed class MegaMatch(int start1, int start2, int length)
    {
        public int Start1 { get; set; } = start1;

        public int Start2 { get; set; } = start2;

        public int Length { get; set; } = length;
    }

    /// <summary>
    /// <see cref="MegaMerge(string, string, int, bool, string, string, CancellationToken)"/> for two sequences read with the qualifiers (A5-7; ACD type <c>DNA</c> for both): megamerger trims both sequences
    /// (<c>ajSeqTrim</c>: the range of <c>-sbegin1</c> ... of the possibly reverse-complemented sequence), merges the trimmed texts and names them as the report does (the sequence names; positions are those of the trimmed texts).
    /// </summary>
    /// <exception cref="ArgumentException">A sequence is empty after reading, or <c>-sprotein</c> is given.</exception>
    public static MegaMergeResult MegaMerge(EmbossSequence first, EmbossSequence second, int wordSize = 20, bool preferFirst = false, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);
        var v1 = EmbossInputSequence.Read(first.Residues, first.Qualifiers, EmbossSequenceType.Dna);
        var v2 = EmbossInputSequence.Read(second.Residues, second.Qualifiers, EmbossSequenceType.Dna);
        return MegaMerge(v1.Residues, v2.Residues, wordSize, preferFirst, first.Name, second.Name, cancellationToken);
    }

    /// <summary>
    /// EMBOSS <c>megamerger</c> (6.6.0; <c>emboss/megamerger.c</c>, <c>nucleus/embword.c</c> <c>embWordBuildMatchTable</c> / <c>embWordMatchMin</c> / <c>wordDeadZone</c>): merges two large overlapping DNA sequences.
    /// All exact word matches of <paramref name="wordSize"/> symbols are listed (the same maximal matches as <see cref="SeqMatchAll"/>), ordered longest first, and reduced to a
    /// minimal non-overlapping set: the longest match is kept, every match inside its dead zone is dropped and one that crosses it is truncated (and dropped when shorter than the word), until none is left;
    /// the kept matches are ordered by position. The merged sequence starts with the part of the sequence that precedes the first match (the first sequence when its overlap is at the start of the
    /// second, and vice versa; when neither overlap is at a start the longer prefix, written upper case), continues through the matches (taken from the first sequence, lower case)
    /// and, between two matches, takes the region from the sequence whose mismatch region is farther from its ends (<paramref name="preferFirst"/>: always the first sequence), upper case, and ends with the part after the last
    /// match (the longer remainder in upper case when neither overlap is at an end). The report names every overlap, mismatch region and choice exactly as the program writes it.
    /// </summary>
    /// <param name="sequence1">First sequence (letters, any case).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="wordSize">Word size (<c>-wordsize</c>, default 20, at least 2).</param>
    /// <param name="preferFirst">Use the first sequence over every mismatch (<c>-prefer</c>).</param>
    /// <param name="name1">Name of the first sequence in the report.</param>
    /// <param name="name2">Name of the second sequence in the report.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>The merge; with an empty <see cref="MegaMergeResult.MergedSequence"/> when the second sequence is shorter than the word or has no match.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="wordSize"/> &lt; 2.</exception>
    /// <exception cref="ArgumentException">The first sequence has no word table (shorter than the word, or every word starts in a run of N): the program stops with "No match found".</exception>
    /// <exception cref="InvalidOperationException">Two matches lie on one diagonal in a way the program's truncation does not handle (the program dies there).</exception>
    public static MegaMergeResult MegaMerge(
        string sequence1, string sequence2, int wordSize = 20, bool preferFirst = false, string name1 = "seq1", string name2 = "seq2", CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(name1);
        ArgumentNullException.ThrowIfNull(name2);
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 2);
        sequence1 = EmbossSequenceReader.Clean(sequence1);
        sequence2 = EmbossSequenceReader.Clean(sequence2);
        AlignmentCancellation.Poll(cancellationToken);

        var matcher = new PairWordMatcher(sequence1, wordSize);
        if (!matcher.HasTable)
            throw new ArgumentException($"The first sequence has no word of {wordSize} symbols (shorter than the word, or every word starts inside a run of the unknown symbol): megamerger stops with 'No match found'.", nameof(sequence1));
        var matches = sequence2.Length < wordSize
            ? []
            : matcher.Matches(sequence2).Select(m => new MegaMatch(m.PositionInText, m.PositionInQuery, m.Length)).ToList();
        AlignmentCancellation.Poll(cancellationToken);

        MinimalMatchSet(matches, wordSize, cancellationToken);
        if (matches.Count == 0)
            return new MegaMergeResult(string.Empty, string.Empty, []);

        var (merged, report) = MegaMergeText(matches, sequence1, sequence2, preferFirst, name1, name2);
        return new MegaMergeResult(merged, report, matches.Select(m => new WordMatch(0, 1, m.Start1, m.Start2, m.Length)).ToList());
    }

    // wordMatchCmp: longest first, then by start in the first sequence, then in the second
    internal static int MegaOrder(MegaMatch a, MegaMatch b)
    {
        if (a.Length != b.Length)
            return a.Length < b.Length ? 1 : -1;
        if (a.Start1 != b.Start1)
            return a.Start1 > b.Start1 ? 1 : -1;
        if (a.Start2 != b.Start2)
            return a.Start2 > b.Start2 ? 1 : -1;
        return 0;
    }

    // embWordMatchMin. The program keeps the not yet decided matches in a list sorted by wordMatchCmp, takes the first one into the minimal set, removes every match inside its dead zone, shortens the ones
    // that overlap it and sorts the list again when a match was shortened. Same result here (and the same exceptions in the same order: every remaining match is examined against every taken one),
    // without the O(k) shift of RemoveAt(0) / RemoveAt(i) and without sorting the whole list again (A4-6): the head is an index, one pass compacts the survivors in place, and the shortened matches
    // (whose key changed) are sorted on their own and merged back into the sorted survivors.
    internal static void MinimalMatchSet(List<MegaMatch> matchlist, int wordLength, CancellationToken cancellationToken = default)
    {
        var comparer = Comparer<MegaMatch>.Create(MegaOrder);
        matchlist.Sort(comparer);
        var minlist = new List<MegaMatch>();
        var shortened = new List<MegaMatch>();
        var merged = new List<MegaMatch>();
        int head = 0;
        while (head < matchlist.Count)
        {
            AlignmentCancellation.Poll(cancellationToken);
            var thismatch = matchlist[head++];
            minlist.Add(thismatch);
            int deadx1 = thismatch.Start1, deady1 = thismatch.Start2;
            int deadx2 = thismatch.Start1 + thismatch.Length - 1, deady2 = thismatch.Start2 + thismatch.Length - 1;
            int write = head;
            shortened.Clear();
            for (int read = head; read < matchlist.Count; read++)
            {
                var m = matchlist[read];
                int result = DeadZone(m, deadx1, deady1, deadx2, deady2, wordLength);
                if (result == 1)
                    continue;
                if (result == 2)
                    shortened.Add(m);
                else
                    matchlist[write++] = m;
            }

            matchlist.RemoveRange(write, matchlist.Count - write);
            if (shortened.Count == 0)
                continue;

            // the shortened matches have new keys: sort them and merge them with the survivors (both runs are sorted by wordMatchCmp)
            shortened.Sort(comparer);
            merged.Clear();
            int a = head, b = 0;
            while (a < matchlist.Count && b < shortened.Count)
                merged.Add(MegaOrder(shortened[b], matchlist[a]) < 0 ? shortened[b++] : matchlist[a++]);
            while (a < matchlist.Count)
                merged.Add(matchlist[a++]);
            while (b < shortened.Count)
                merged.Add(shortened[b++]);
            matchlist.RemoveRange(head, matchlist.Count - head);
            matchlist.AddRange(merged);
        }

        matchlist.Clear();
        minlist.Sort((a, b) => a.Start1 != b.Start1 ? a.Start1.CompareTo(b.Start1) : a.Start2.CompareTo(b.Start2));
        matchlist.AddRange(minlist);
    }

    // wordDeadZone: 0 live, 1 dead (remove), 2 truncated
    internal static int DeadZone(MegaMatch match, int deadx1, int deady1, int deadx2, int deady2, int minlength)
    {
        int startx = match.Start1, starty = match.Start2;
        int endx = match.Start1 + match.Length - 1, endy = match.Start2 + match.Length - 1;

        if (startx > deadx2 && starty > deady2)
            return 0;
        if (endx < deadx1 && endy < deady1)
            return 0;
        if (starty >= deady1 && endx <= deadx2)
            return 1;
        if (endy <= deady2 && startx >= deadx1)
            return 1;

        string where = $"startx={startx}, starty={starty}, endx={endx}, endy={endy}, deadx1={deadx1}, deady1={deady1}, deadx2={deadx2}, deady2={deady2}";
        if (endy < deady2)
        {
            if (startx - starty < deadx1 - deady1)
                match.Length = deady1 - starty;
            else if (startx - starty > deadx1 - deady1)
                match.Length = deadx1 - startx;
            else
                throw new InvalidOperationException("Found a match where match is on the same diagonal: " + where);
        }
        else if (starty > deady1)
        {
            if (startx - starty < deadx1 - deady1)
            {
                match.Length = endx - deadx2;
                match.Start1 = deadx2 + 1;
                match.Start2 += deadx2 - startx + 1;
            }
            else if (startx - starty > deadx1 - deady1)
            {
                match.Length = endy - deady2;
                match.Start1 += deady2 - starty + 1;
                match.Start2 = deady2 + 1;
            }
            else
            {
                throw new InvalidOperationException("Found a match where match is on the same diagonal: " + where);
            }
        }
        else
        {
            throw new InvalidOperationException("Found a match that was not caught by any case: " + where);
        }

        return match.Length < minlength ? 1 : 2;
    }

    // megamerger_Merge (positions are 0-based; unsigned arithmetic of the program where it matters)
    private static (string Merged, string Report) MegaMergeText(List<MegaMatch> matches, string sequence1, string sequence2, bool prefer, string name1, string name2)
    {
        var report = new StringBuilder();
        var seqstr = new StringBuilder();
        string s1 = sequence1.ToLowerInvariant(), s2 = sequence2.ToLowerInvariant();
        int len1 = sequence1.Length, len2 = sequence2.Length;
        uint prev1end = 0, prev2end = 0;
        string Sub(string s, long first, long last) => EmbossText.Substring(s, (int)first, (int)last);
        string Num(long v) => v.ToString(CultureInfo.InvariantCulture);

        report.Append($"# Report of megamerger of: {name1} and {name2}\n\n");
        int count = 0;
        MegaMatch last = matches[0];
        foreach (var p in matches)
        {
            last = p;
            if (count++ == 0)
            {
                report.Append($"{name1} overlap starts at {Num(p.Start1 + 1)}\n");
                report.Append($"{name2} overlap starts at {Num(p.Start2 + 1)}\n\n");
                if (p.Start1 == 0)
                {
                    if (p.Start2 > 0)
                    {
                        report.Append($"Using {name2} 1-{Num(p.Start2)} as the initial sequence\n\n");
                        seqstr.Clear().Append(Sub(s2, 0, p.Start2 - 1));
                    }
                }
                else if (p.Start2 == 0)
                {
                    report.Append($"Using {name1} 1-{Num(p.Start1)} as the initial sequence\n\n");
                    seqstr.Clear().Append(Sub(s1, 0, p.Start1 - 1));
                }
                else
                {
                    report.Append("WARNING!\n");
                    report.Append("Neither sequence's overlap is at the start of the sequence\n");
                    if (p.Start1 > p.Start2)
                    {
                        report.Append($"Using {name1} 1-{Num(p.Start1)} as the initial sequence\n\n");
                        seqstr.Clear().Append(Sub(s1, 0, p.Start1 - 1).ToUpperInvariant());
                    }
                    else
                    {
                        report.Append($"Using {name2} 1-{Num(p.Start2)} as the initial sequence\n\n");
                        seqstr.Clear().Append(Sub(s2, 0, p.Start2 - 1).ToUpperInvariant());
                    }
                }
            }
            else
            {
                report.Append("\nWARNING!\nMismatch region found:\n");
                report.Append(prev1end < p.Start1 ? $"Mismatch {name1} {Num(prev1end + 1)}-{Num(p.Start1)}\n" : $"Mismatch {name1} {Num(prev1end)}\n");
                report.Append(prev2end < p.Start2 ? $"Mismatch {name2} {Num(prev2end + 1)}-{Num(p.Start2)}\n" : $"Mismatch {name2} {Num(prev2end)}\n");
                if (prefer)
                {
                    seqstr.Append(Sub(s1, prev1end, p.Start1 - 1).ToUpperInvariant());
                }
                else
                {
                    uint mid1 = (uint)(prev1end + p.Start1 - 1) / 2, mid2 = (uint)(prev2end + p.Start2 - 1) / 2;
                    uint dist1 = Math.Min(mid1, unchecked((uint)len1 - mid1 - 1)), dist2 = Math.Min(mid2, unchecked((uint)len2 - mid2 - 1));
                    if (dist1 < dist2)
                    {
                        report.Append($"Mismatch is closer to the ends of {name1}, so use {name2} in the merged sequence\n\n");
                        if (prev2end < p.Start2)
                            seqstr.Append(Sub(s2, prev2end, p.Start2 - 1).ToUpperInvariant());
                    }
                    else
                    {
                        report.Append($"Mismatch is closer to the ends of {name2}, so use {name1} in the merged sequence\n\n");
                        if (prev1end < p.Start1)
                            seqstr.Append(Sub(s1, prev1end, p.Start1 - 1).ToUpperInvariant());
                    }
                }
            }

            report.Append($"Matching region {name1} {Num(p.Start1 + 1)}-{Num(p.Start1 + p.Length)} : {name2} {Num(p.Start2 + 1)}-{Num(p.Start2 + p.Length)}\n");
            report.Append($"Length of match: {Num(p.Length)}\n");
            seqstr.Append(Sub(s1, p.Start1, p.Start1 + p.Length - 1));
            prev1end = (uint)(p.Start1 + p.Length);
            prev2end = (uint)(p.Start2 + p.Length);
        }

        report.Append($"\n{name1} overlap ends at {Num(last.Start1 + last.Length)}\n");
        report.Append($"{name2} overlap ends at {Num(last.Start2 + last.Length)}\n\n");
        if (prev2end >= len2 && prev1end < len1)
        {
            report.Append($"Using {name1} {Num(prev1end + 1)}-{Num(len1)} as the final sequence\n\n");
            seqstr.Append(Sub(s1, prev1end, len1 - 1L));
        }
        else if (prev1end >= len1 && prev2end < len2)
        {
            report.Append($"Using {name2} {Num(prev2end + 1)}-{Num(len2)} as the final sequence\n\n");
            seqstr.Append(Sub(s2, prev2end, len2 - 1L));
        }
        else if (prev1end < len1 && prev2end < len2)
        {
            report.Append("WARNING!\n");
            report.Append("Neither sequence's overlap is at the end of the sequence\n");
            if (len1 - prev1end > len2 - prev2end)
            {
                report.Append($"Using {name1} {Num(prev1end + 1)}-{Num(len1)} as the final sequence\n\n");
                seqstr.Append(Sub(s1, prev1end, len1 - 1L).ToUpperInvariant());
            }
            else
            {
                report.Append($"Using {name2} {Num(prev2end + 1)}-{Num(len2)} as the final sequence\n\n");
                seqstr.Append(Sub(s2, prev2end, len2 - 1L).ToUpperInvariant());
            }
        }

        return (seqstr.ToString(), report.ToString());
    }

    #endregion
}
