// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.Motifs;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Emboss.Motifs;

using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.GenomicAnalyzer;

/// <summary>A hit of <see cref="EmbossMotifPrograms.FindProteinPatternHits(string, string, CancellationToken)"/>: the 0-based start and the length of the matching window.</summary>
/// <param name="Position">0-based start (EMBOSS <c>fuzzpro</c> prints Position + 1).</param>
/// <param name="Length">Window length (a range such as <c>x(2,4)</c> yields hits of several lengths).</param>
/// <param name="Mismatches">Substituted positions of a search with <c>-pmismatch</c> (fuzzpro's <c>*mismatch</c> note, 0 for an exact hit).</param>
public readonly record struct ProteinPatternHit(int Position, int Length, int Mismatches = 0)
{
    /// <summary>Exclusive end of the window (0-based).</summary>
    public int End => Position + Length;
}

/// <summary>A hit of <see cref="EmbossMotifPrograms.FindTranslatedPatternHits(string, string, string, int, CancellationToken)"/> (EMBOSS <c>fuzztran</c>).</summary>
/// <param name="Frame">Translation frame: 1, 2, 3 forward, -1, -2, -3 reverse (EMBOSS transeq convention).</param>
/// <param name="ProteinStart">1-based first residue in the translated frame (the program's <c>*start</c> note).</param>
/// <param name="ProteinEnd">1-based last residue (<c>*end</c>).</param>
/// <param name="Translation">The matched residues (<c>*translation</c>).</param>
/// <param name="Position">0-based start of the codons on the forward strand of the input.</param>
/// <param name="Length">Length in nucleotides (3 per residue).</param>
/// <param name="EmbossStart">The first coordinate the program prints (<c>fuzztran_SourceFeature</c>): forward frames <c>S + frame + 3·start − 4</c> (= <see cref="Position"/> + 1 for <c>S = 1</c>), reverse frames <c>S − |frame| − 3·start + 4</c>, where
/// <c>S</c> is the first position of the sequence used (1, or the <c>-sbegin</c> of the range). The reverse-frame numbers are what the program prints: negative for a whole sequence (it counts the protein position from the start of the sequence with
/// the wrong sign; EMBOSS 6.6.0 defect, reproduced on 3 000+ hits, A5-7 oracle), the gff report lists the smaller of the two first.</param>
/// <param name="EmbossEnd">The second coordinate the program prints: forward frames <c>S + frame + 3·end − 2</c> (= <see cref="Position"/> + <see cref="Length"/> for <c>S = 1</c>), reverse frames <c>S − |frame| − 3·end + 2</c>.</param>
/// <param name="Mismatches">Substituted residues of a search with <c>-pmismatch</c> (the <c>*mismatch</c> note, 0 for an exact hit).</param>
public readonly record struct TranslatedPatternHit(
    int Frame, int ProteinStart, int ProteinEnd, string Translation, int Position, int Length, int? EmbossStart, int? EmbossEnd, int Mismatches = 0)
{
    /// <summary>True for the reverse frames (-1..-3).</summary>
    public bool IsReverseStrand => Frame < 0;
}

public static partial class EmbossMotifPrograms
{
    #region fuzzpro / fuzztran (protein patterns)

    /// <summary>
    /// EMBOSS <c>fuzzpro</c> (6.6.0): every window of a protein sequence that matches a <see cref="ProteinPattern"/> (PROSITE-like: letters, <c>x</c>, <c>[ILV]</c>, <c>{P}</c>,
    /// <c>x(3)</c>, <c>x(2,4)</c>, <c>&lt;</c> / <c>&gt;</c>), exact (the binary has no <c>-pmismatch</c>).
    /// </summary>
    /// <remarks>
    /// Like the binary, every matching window length of a range is reported (overlapping hits included); a start-anchored (<c>&lt;</c>) range pattern reports only its longest match
    /// (the regular-expression search class of <c>embpat.c</c>). The sequence is case-insensitive; a single trailing stop <c>*</c> is not part of it (the binary reads it off). Hits are
    /// ordered by position, then length. Verified against fuzzpro 6.6.0 on 900 random syntax-tree patterns over sequences with X, B, Z, U, J, O, * and lower case
    /// (docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md).
    /// </remarks>
    /// <param name="protein">The protein sequence.</param>
    /// <param name="pattern">The pattern text (<see cref="ProteinPattern.Parse"/>).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="maxMismatches">Maximum substituted residues per hit (<c>-pmismatch</c>, A5-8c2): the search of <c>embPatBruteForce</c> / Baeza-Yates-Perleberg / Tarhio-Ukkonen-Bleasby; with a range only the mandatory part may mismatch.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The pattern is invalid.</exception>
    public static IReadOnlyList<ProteinPatternHit> FindProteinPatternHits(string protein, string pattern, CancellationToken cancellationToken = default, int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(protein);
        return FindProteinPatternHits(protein, ProteinPattern.Parse(pattern), cancellationToken, maxMismatches);
    }

    /// <summary>
    /// <see cref="FindProteinPatternHits(string, string, CancellationToken)"/> for a sequence read with the qualifiers <c>-sbegin -send -slower -supper -sprotein</c> (A5-7; ACD type <c>protein</c>): the search is on the range, and
    /// <see cref="ProteinPatternHit.Position"/> counts from the start of the original sequence (<c>-sbegin</c> offset added). A protein is never reverse-complemented; <c>-snucleotide</c> is refused.
    /// </summary>
    /// <exception cref="ArgumentException">The sequence is empty after reading, or <c>-snucleotide</c> is given.</exception>
    public static IReadOnlyList<ProteinPatternHit> FindProteinPatternHits(string protein, EmbossSequenceQualifiers? qualifiers, string pattern, CancellationToken cancellationToken = default, int maxMismatches = 0)
    {
        var view = EmbossInputSequence.Read(protein, qualifiers, EmbossSequenceType.Protein);
        var hits = FindProteinPatternHits(view.Residues, pattern, cancellationToken, maxMismatches);
        return view.Offset == 0 ? hits : [.. hits.Select(h => h with { Position = h.Position + view.Offset })];
    }

    /// <summary><see cref="FindProteinPatternHits(string, string, CancellationToken)"/> with an already parsed pattern.</summary>
    /// <param name="protein">The protein sequence.</param>
    /// <param name="pattern">The parsed pattern.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="maxMismatches">Maximum substituted residues per hit (<c>-pmismatch</c>, A5-8c2): the search of <c>embPatBruteForce</c> / Baeza-Yates-Perleberg / Tarhio-Ukkonen-Bleasby; with a range only the mandatory part may mismatch.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static IReadOnlyList<ProteinPatternHit> FindProteinPatternHits(string protein, ProteinPattern pattern, CancellationToken cancellationToken = default, int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(protein);
        ArgumentNullException.ThrowIfNull(pattern);
        protein = EmbossSequenceReader.Clean(protein);                // the fuzzpro reader drops the gap characters and every other non-letter but * (A4-3), then reads a trailing stop off
        return ProteinPatternHitsCore(protein.EndsWith('*') ? protein[..^1] : protein, pattern, cancellationToken, maxMismatches);
    }

    private static List<ProteinPatternHit> ProteinPatternHitsCore(string text, ProteinPattern pattern, CancellationToken cancellationToken, int maxMismatches = 0)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        var bits = new uint[text.Length];
        for (int i = 0; i < bits.Length; i++)
            bits[i] = ProteinPattern.ResidueBit(text[i]);

        if (maxMismatches > 0)
        {
            if (pattern.UsesTarhioUkkonenBleasby(out int realLength))
            {
                if (maxMismatches >= realLength)
                    throw new ArgumentException($"embPatCompileII: Mismatches ({maxMismatches}) must be less than the real pattern length ({realLength})", nameof(maxMismatches));
                // type 6: the skip-table search of embpat.c (it jumps over windows that hold a symbol that is no letter)
                string tubText = pattern.TubText(out int m);
                var upper = text.ToUpperInvariant();
                return [.. EmbossTubSearch.Search(tubText, upper, EmbossTubSearch.Init(tubText, m, maxMismatches), m, maxMismatches, pattern.AnchoredAtStart, pattern.AnchoredAtEnd)
                    .Select(h => new ProteinPatternHit(h.Start, h.Length, h.Mismatches))];
            }

            return ProteinPatternHitsApproximate(bits, pattern, maxMismatches, cancellationToken);
        }

        var hits = new HashSet<(int Position, int Length)>();
        foreach (var expansion in pattern.Expand())
        {
            AlignmentCancellation.Poll(cancellationToken);
            int length = expansion.Length;
            for (int pos = 0; pos + length <= bits.Length; pos++)
            {
                if ((pos & 0xFFFF) == 0)
                    AlignmentCancellation.Poll(cancellationToken);        // A4-1: a long sequence with one expansion has no other poll
                if (pattern.AnchoredAtStart && pos != 0)
                    break;
                if (pattern.AnchoredAtEnd && pos + length != bits.Length)
                    continue;
                bool ok = true;
                for (int j = 0; j < length && ok; j++)
                    ok = (bits[pos + j] & expansion[j]) != 0;
                if (ok)
                    hits.Add((pos, length));
            }
        }

        // an all-optional range pattern: the empty hit at every start of the sequence (the search loop stops at the last residue); anchored at the end ('>') the regular expression 'pattern$' matches the empty string AT the end of the
        // sequence, found only when no match began at the last residue (A5-8c3b, F217)
        if (pattern.HasRange && pattern.ExpandPaths(includeEmpty: true).Any(path => path.Sets.Length == 0))
        {
            if (pattern.AnchoredAtEnd)
            {
                if (!pattern.AnchoredAtStart && bits.Length > 0 && !hits.Contains((bits.Length - 1, 1)))
                    hits.Add((bits.Length, 0));
            }
            else
            {
                for (int pos = 0; pos < bits.Length && !(pattern.AnchoredAtStart && pos != 0); pos++)
                    hits.Add((pos, 0));
            }
        }

        // regex class (embPatMatchFindAll with nterm): a start-anchored range pattern returns only the longest match at the anchor
        if (pattern.AnchoredAtStart && pattern.HasRange && hits.Count > 1)
        {
            int longest = hits.Max(h => h.Length);
            hits.RemoveWhere(h => h.Length != longest);
        }

        return hits.OrderBy(h => h.Position).ThenBy(h => h.Length).Select(h => new ProteinPatternHit(h.Position, h.Length)).ToList();
    }

    // fuzzpro / fuzztran with -pmismatch n (A5-8c2, F195): embPatGetTypeII gives a protein pattern with mismatches the Baeza-Yates-Perleberg, Tarhio-Ukkonen-Bleasby or recursive brute-force search (embPatBruteForce -> patOUBrute); all
    // report a window with its Hamming distance <= n. With a range the recursive matcher lets a mismatch fall only on the mandatory part (the extra repeats must match exactly) and fails when the sequence ends while elements are left over.
    // The recursion pushes one hit per successful unrolling (embPatPushHit), so a window that several unrollings match is listed several times (each with its own distance); the list is in push order: start ascending, then the
    // unrollings in the order of the recursion (fewer repeats of the earlier element first), which is the order of ProteinPattern.ExpandPaths.
    private static List<ProteinPatternHit> ProteinPatternHitsApproximate(uint[] bits, ProteinPattern pattern, int maxMismatches, CancellationToken cancellationToken)
    {
        bool recursiveMatcher = pattern.HasRange;
        var paths = pattern.ExpandPaths(includeEmpty: true);          // the empty unrolling of an all-optional pattern is a hit at every start (A5-8c3b, F217)
        var hits = new List<ProteinPatternHit>();
        for (int pos = 0; pos < bits.Length; pos++)
        {
            if ((pos & 0xFFF) == 0)
                AlignmentCancellation.Poll(cancellationToken);
            if (pattern.AnchoredAtStart && pos != 0)
                break;
            foreach (var path in paths)
            {
                int length = path.Sets.Length;
                if (pos + length > bits.Length)
                    continue;
                if (pattern.AnchoredAtEnd && pos + length != bits.Length)
                    continue;
                int bad = 0;
                bool allowed = true;
                for (int j = 0; j < length && allowed; j++)
                {
                    if ((bits[pos + j] & path.Sets[j]) != 0)
                        continue;
                    allowed = ++bad <= maxMismatches && !(recursiveMatcher && Array.IndexOf(path.ExtraPositions, j) >= 0);
                }

                if (!allowed || (recursiveMatcher && path.EndsWithEmptyElement && pos + length == bits.Length))
                    continue;
                hits.Add(new ProteinPatternHit(pos, length, bad));
            }
        }

        // a pattern without a range is one unrolling (the Hamming engines list each window once)
        return hits;
    }

    /// <summary>
    /// EMBOSS <c>fuzztran</c> (6.6.0): translates a nucleotide sequence in the chosen frames with an EMBOSS genetic-code table and searches every translation with a
    /// <see cref="ProteinPattern"/> (<see cref="FindProteinPatternHits(string, ProteinPattern, CancellationToken)"/>, without reading a trailing stop off).
    /// </summary>
    /// <remarks>
    /// Frames <c>1</c>, <c>2</c>, <c>3</c> (forward), <c>-1</c>, <c>-2</c>, <c>-3</c> (reverse), <c>F</c> (forward three), <c>R</c> (reverse three) or <c>6</c> (all) — the program's
    /// <c>-frame</c> — are translated as <c>transeq</c> does (<see cref="EmbossTranslator.TranslateSixFrames(string, GeneticCode?, SixFrameConvention)"/> with
    /// <see cref="SixFrameConvention.EmbossTranseq"/>; only the complete codons are searched, the dangling-base residue of transeq is dropped) and the sequence is upper-cased first. <paramref name="table"/> is the program's <c>-table</c> menu (0 =
    /// standard; otherwise the NCBI table number). Hits are listed per frame in the order 1, 2, 3, -1, -2, -3 (those requested), by protein position then length. Each hit carries the
    /// forward-strand nucleotide location (<see cref="TranslatedPatternHit.Position"/>, <see cref="TranslatedPatternHit.Length"/>) and the numbers the binary prints
    /// (<see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/>: the program's 1-based coordinates, exact for the forward frames; for the reverse frames EMBOSS 6.6.0
    /// prints meaningless numbers from a defect of its <c>fuzztran_SourceFeature</c> and they are not reproduced). Verified against fuzztran 6.6.0 on 700 random cases (docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md) and, with every frame, on the A5-7 oracle.
    /// </remarks>
    /// <param name="nucleotides">The nucleotide sequence (IUPAC codes allowed).</param>
    /// <param name="pattern">The protein pattern text.</param>
    /// <param name="frames">The frame selection: 1, 2, 3, -1, -2, -3, F, R or 6 (case-insensitive).</param>
    /// <param name="table">EMBOSS genetic code (0 standard, otherwise an NCBI table number).</param>
    /// <param name="cancellationToken">Cancellation token, polled per frame and every 65 536 residues of a translation (A4-1).</param>
    /// <param name="maxMismatches">Maximum substituted residues per hit (<c>-pmismatch</c>, A5-8c2): the search of <c>embPatBruteForce</c> / Baeza-Yates-Perleberg / Tarhio-Ukkonen-Bleasby; with a range only the mandatory part may mismatch.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The pattern or the frame selection is invalid, or the table is unknown.</exception>
    public static IReadOnlyList<TranslatedPatternHit> FindTranslatedPatternHits(string nucleotides, string pattern, string frames = "1", int table = 0, CancellationToken cancellationToken = default, int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(nucleotides);
        string text = EmbossSequenceReader.Clean(nucleotides).ToUpperInvariant();      // the fuzztran reader drops the gap characters and every other non-letter but * (A4-3)
        return TranslatedPatternHitsCore(text, 1, text.Length, pattern, frames, table, cancellationToken, maxMismatches);
    }

    // fuzztran on the text it translates (the sequence as used, upper case): `begin` / `end` are the 1-based first and last position of that text in the sequence the program numbers (ajSeqallGetseqBegin / End)
    private static IReadOnlyList<TranslatedPatternHit> TranslatedPatternHitsCore(string text, int begin, int end, string pattern, string frames, int table, CancellationToken cancellationToken, int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(frames);
        var parsed = ProteinPattern.Parse(pattern);
        string selector = frames.Trim().ToUpperInvariant();
        int[] selected = selector switch
        {
            "F" => [1, 2, 3],
            "R" => [-1, -2, -3],
            "6" => [1, 2, 3, -1, -2, -3],
            "1" => [1],
            "2" => [2],
            "3" => [3],
            "-1" => [-1],
            "-2" => [-2],
            "-3" => [-3],
            _ => throw new ArgumentException("frames must be 1, 2, 3, -1, -2, -3, F, R or 6 (fuzztran -frame).", nameof(frames)),
        };
        var code = GeneticCode.GetByTableNumberOrThrow(table == 0 ? 1 : table, nameof(table));
        var translations = EmbossTranslator.TranslateSixFrames(text, code, SixFrameConvention.EmbossTranseq);
        int length = text.Length;
        var result = new List<TranslatedPatternHit>();
        foreach (int frame in selected)
        {
            AlignmentCancellation.Poll(cancellationToken);
            // fuzztran translates only complete codons (the dangling-base residue transeq appends is not searched) and searches the translation as it is (a stop '*' stays)
            string protein = translations[frame].Sequence;
            int complete = Math.Max(0, (length - (Math.Abs(frame) - 1)) / 3);
            if (protein.Length > complete)
                protein = protein[..complete];
            foreach (var hit in ProteinPatternHitsCore(protein, parsed, cancellationToken, maxMismatches))
            {
                int s = hit.Position + 1, e = hit.Position + hit.Length;
                int framenum = Math.Abs(frame);
                int position = frame > 0 ? framenum + 3 * s - 4 : length - (framenum - 1) - 3 * e;
                // fuzztran_SourceFeature(start = S): forward S + frame + 3 start - 4, reverse S - |frame| - 3 start + 4; the program passes the END of the sequence instead of its begin for frame -3 of -frame R / 6
                int anchor = frame == -3 && selector is "R" or "6" ? end : begin;
                int embossStart = frame > 0 ? anchor + framenum + 3 * s - 4 : anchor - framenum - 3 * s + 4;
                int embossEnd = frame > 0 ? anchor + framenum + 3 * e - 2 : anchor - framenum - 3 * e + 2;
                result.Add(new TranslatedPatternHit(frame, s, e, EmbossText.SubstringFromEnd(protein, hit.Position, hit.Position + hit.Length - 1), position, 3 * hit.Length, embossStart, embossEnd, hit.Mismatches));
            }
        }

        return result;
    }

    /// <summary>
    /// <see cref="FindTranslatedPatternHits(string, string, string, int, CancellationToken)"/> for a sequence read with the qualifiers <c>-sbegin -send -sreverse -slower -supper -snucleotide -sprotein</c> (A5-7; ACD type <c>dna</c>):
    /// the translation is of the sequence as the program cuts it (reverse-complemented for <c>-sreverse</c>, then the text from <c>ajSeqallGetseqBegin</c> to <c>ajSeqallGetseqEnd</c>: a begin or end given on the command line is used as given, not mirrored, so with
    /// <c>-sreverse</c> the range is applied to the reverse complement), and the hits are numbered from the first position of that text (<see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/>;
    /// <see cref="TranslatedPatternHit.Position"/> stays relative to the text translated).
    /// </summary>
    /// <param name="nucleotides">The nucleotide sequence.</param>
    /// <param name="qualifiers">The sequence qualifiers; null = none.</param>
    /// <param name="pattern">The protein pattern text.</param>
    /// <param name="frames">The frame selection: 1, 2, 3, -1, -2, -3, F, R or 6.</param>
    /// <param name="table">EMBOSS genetic code (0 standard).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="maxMismatches">Maximum substituted residues per hit (<c>-pmismatch</c>, A5-8c2): the search of <c>embPatBruteForce</c> / Baeza-Yates-Perleberg / Tarhio-Ukkonen-Bleasby; with a range only the mandatory part may mismatch.</param>
    /// <exception cref="ArgumentException">The sequence is empty after reading, or <c>-sprotein</c> is given (fuzztran needs a nucleotide sequence).</exception>
    public static IReadOnlyList<TranslatedPatternHit> FindTranslatedPatternHits(string nucleotides, EmbossSequenceQualifiers? qualifiers, string pattern, string frames = "1", int table = 0, CancellationToken cancellationToken = default, int maxMismatches = 0)
    {
        var view = EmbossInputSequence.Read(nucleotides, qualifiers, EmbossSequenceType.Dna);
        var q = qualifiers ?? EmbossSequenceQualifiers.None;
        // fuzztran does not trim the sequence: it cuts the reverse-complemented text at ajSeqallGetseqBegin / End, which take the begin / end AS GIVEN when they are set (not mirrored) and the sequence's own (mirrored) otherwise
        int n = view.Full.Length;
        int rawBegin = q.Begin != 0 ? q.Begin : view.Begin;
        int rawEnd = q.End != 0 ? q.End : view.End;
        int begin = 1, end = n;
        if (rawBegin != 0)
            begin = 1 + EmbossInputSequence.PositionStart(n, 0, rawBegin > 0 ? rawBegin - 1 : rawBegin);
        if (rawEnd != 0)
            end = 1 + EmbossInputSequence.PositionStart(n, begin - 1, rawEnd > 0 ? rawEnd - 1 : rawEnd);
        return TranslatedPatternHitsCore(SequenceCase.ToUpperAscii(view.Full[(begin - 1)..end]), begin, end, pattern, frames, table, cancellationToken, maxMismatches);
    }

    #endregion
}
