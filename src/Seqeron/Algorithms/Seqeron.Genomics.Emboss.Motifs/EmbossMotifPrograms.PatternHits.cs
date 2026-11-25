// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.Motifs;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Emboss.Motifs;

using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.GenomicAnalyzer;

// B09 A2-32 / A2-35: the fuzznuc pattern language (FuzznucPattern) on a sequence, and the string overloads of the motif finders.
// The scan itself is the canonical ApproximateMatcher (Hamming, IUPAC pattern codes) run on every fixed-length expansion of the pattern.
public static partial class EmbossMotifPrograms
{
    #region Pattern hits (fuzznuc -pattern language)

    /// <summary>
    /// Every hit of an EMBOSS <c>fuzznuc</c> pattern (<see cref="FuzznucPattern"/>: IUPAC letters, <c>[ACG]</c> classes, <c>{AG}</c> exclusions, <c>N(3)</c> repeats,
    /// <c>A(2,4)</c> length ranges, <c>&lt;</c> / <c>&gt;</c> anchors) in a DNA sequence, on the forward strand and optionally on the reverse-complement strand
    /// (<c>fuzznuc -complement</c>), with at most <paramref name="maxMismatches"/> substituted positions per hit (<c>-pmismatch</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// A hit is a window <c>S[p..p+L)</c>: some fixed-length expansion of the pattern (<see cref="FuzznucPattern.Expand"/>) of length L agrees with it at all but at most
    /// k positions (a base agrees when it is in the base set of the expansion's IUPAC letter; Hamming distance, substitutions only, no gaps). Like the binary, a range reports
    /// <b>every</b> window length that matches (<c>A(2,3)</c> on <c>AAAA</c> gives 2-mers at 0, 1, 2 and 3-mers at 0, 1), overlapping hits included; the same window reached through
    /// several expansions is one hit with the smallest mismatch count (the binary's <c>*mismatch</c> note). A reverse hit is a window that matches the
    /// <see cref="FuzznucPattern.ReverseComplement"/> of the pattern; its <see cref="PatternHit.Position"/> is the forward-strand start (fuzznuc / GFF convention) and the anchors swap
    /// (<c>&lt;ACG</c> on the reverse strand must end at the last base of the sequence). <c>&lt;</c> requires <c>Position = 0</c>, <c>&gt;</c> requires <c>Position + Length = n</c>.
    /// Hits are ordered by position, then length, forward before reverse.
    /// </para>
    /// <para>
    /// Two fuzznuc search classes are reproduced exactly (EMBOSS <c>embpat.c</c>, read and checked on the binary): with <paramref name="maxMismatches"/> &gt; 0 a pattern that has
    /// a range <c>X(n,m)</c> is searched by a recursive matcher (<c>patOUBrute</c>) that tolerates mismatches only in the mandatory part — the fixed elements and the first n
    /// repeats of each range; the extra repeats beyond n must match exactly (<see cref="FuzznucPattern.ExpandWithExtras"/>) — and that matcher fails when pattern elements are
    /// left over at the end of the scanned strand, so a hit whose last element is an empty optional element never ends on the last base; and an exact (k = 0) range pattern that
    /// starts with <c>&lt;</c> is searched by a regular expression that returns only the longest match at the anchor (on the reverse strand: only the longest
    /// window that ends at the last base). Fixed-length patterns, and every other combination, report all windows.
    /// </para>
    /// <para>
    /// The scan is the canonical <see cref="ApproximateMatcher"/> Hamming matcher over each expansion (O(n·m) per expansion, no index). Verified against EMBOSS fuzznuc 6.6.0 on
    /// 3 000+ random syntax-tree patterns (positions, ends, strands and mismatch notes), see docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md.
    /// </para>
    /// </remarks>
    /// <param name="sequence">Subject sequence.</param>
    /// <param name="pattern">The pattern text.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit (default 0 = exact).</param>
    /// <param name="bothStrands">Also search the reverse-complement strand.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="pattern"/> is null.</exception>
    /// <exception cref="ArgumentException">The pattern is invalid (see <see cref="FuzznucPattern.Parse"/>).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static IReadOnlyList<PatternHit> FindPatternHits(
        DnaSequence sequence, string pattern, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindPatternHits(sequence, FuzznucPattern.Parse(pattern), maxMismatches, bothStrands, cancellationToken);
    }

    /// <summary>
    /// <see cref="FindPatternHits(DnaSequence, string, int, bool, CancellationToken)"/> with an already parsed pattern.
    /// </summary>
    /// <param name="sequence">Subject sequence.</param>
    /// <param name="pattern">The parsed pattern.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit.</param>
    /// <param name="bothStrands">Also search the reverse-complement strand.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static IReadOnlyList<PatternHit> FindPatternHits(
        DnaSequence sequence, FuzznucPattern pattern, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return PatternHitsCore(sequence.Sequence, pattern, maxMismatches, bothStrands, cancellationToken);
    }

    /// <summary>
    /// <see cref="FindPatternHits(DnaSequence, string, int, bool, CancellationToken)"/> on a nucleotide <b>string</b> (N runs, IUPAC codes, lower case). By default the sequence side is
    /// literal as in the canonical motif finders — a symbol other than A/C/G/T matches no pattern symbol (not even <c>N</c>) and counts as a mismatch. With
    /// <paramref name="sequenceAmbiguity"/> = <c>true</c> the sequence is read exactly as the <c>fuzznuc</c> binary reads it (A3-12a, measured on 6.6.0, 1 350 oracle rows): a pattern IUPAC letter matches every
    /// sequence symbol whose base set lies within the letter's (pattern <c>R</c> matches A, G and R; sequence <c>R</c> matches pattern R, V, D and N but not A, M or Y), <c>N</c> matches every symbol (X included),
    /// a class <c>[ACG]</c> matches its plain bases only, an exclusion <c>{AG}</c> matches every symbol not listed (ambiguity codes and X included), U reads as T, and the reverse strand is the IUPAC reverse complement
    /// of the sequence searched with the same pattern (see <see cref="FuzznucPattern.Matches(int, char)"/>). Everything else (mismatches, ranges, anchors, search classes) is unchanged.
    /// </summary>
    /// <param name="sequence">Subject nucleotide string.</param>
    /// <param name="pattern">The pattern text.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit.</param>
    /// <param name="bothStrands">Also search the reverse-complement strand.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="sequenceAmbiguity">Read IUPAC ambiguity codes in the sequence as fuzznuc does (default: they match nothing).</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The pattern is invalid.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static IReadOnlyList<PatternHit> FindPatternHits(
        string sequence, string pattern, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default, bool sequenceAmbiguity = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindPatternHits(sequence, FuzznucPattern.Parse(pattern), maxMismatches, bothStrands, cancellationToken, sequenceAmbiguity);
    }

    /// <summary>
    /// <see cref="FindPatternHits(string, string, int, bool, CancellationToken, bool)"/> for a sequence read with the qualifiers <c>-sbegin -send -sreverse -slower -supper -snucleotide -sprotein</c> (A5-7; ACD type <c>dna</c>, IUPAC letters of the sequence
    /// read as fuzznuc reads them). The search is on the sequence as used (the range, reverse-complemented for <c>-sreverse</c>) and the hits are reported on the ORIGINAL strand, as the program's feature table is reversed with the sequence:
    /// <see cref="PatternHit.Position"/> is the 0-based forward-strand start in the original sequence and <see cref="PatternHit.IsReverseStrand"/> is flipped for a reversed sequence (a pattern found on the reverse complement of the original is a hit on its minus strand).
    /// </summary>
    /// <param name="sequence">The nucleotide sequence.</param>
    /// <param name="qualifiers">The sequence qualifiers; null = none.</param>
    /// <param name="pattern">The pattern text.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit (<c>-pmismatch</c>).</param>
    /// <param name="bothStrands">Also search the reverse-complement strand (<c>-complement</c>).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="programOrder">A range pattern with mismatches: return the hits in the order of the program's list, a window once per matching unrolling (F195), instead of one hit per window sorted by position (the report programs use it).</param>
    /// <exception cref="ArgumentException">The sequence is empty after reading, or <c>-sprotein</c> is given (fuzznuc needs a nucleotide sequence), or the pattern is invalid.</exception>
    public static IReadOnlyList<PatternHit> FindPatternHits(
        string sequence, EmbossSequenceQualifiers? qualifiers, string pattern, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default, bool programOrder = false)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        return FindPatternHits(sequence, qualifiers, FuzznucPattern.Parse(pattern), maxMismatches, bothStrands, cancellationToken, programOrder);
    }

    /// <summary><see cref="FindPatternHits(string, EmbossSequenceQualifiers?, string, int, bool, CancellationToken)"/> with an already parsed pattern.</summary>
    public static IReadOnlyList<PatternHit> FindPatternHits(
        string sequence, EmbossSequenceQualifiers? qualifiers, FuzznucPattern pattern, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default, bool programOrder = false)
    {
        var view = EmbossInputSequence.Read(sequence, qualifiers, EmbossSequenceType.Dna);
        var hits = FindPatternHits(view.Residues, pattern, maxMismatches, bothStrands, cancellationToken, sequenceAmbiguity: true, programOrder: programOrder);
        var mapped = new List<PatternHit>(hits.Count);
        foreach (var h in hits)
        {
            var (start, _) = view.ReportedRange(h.Position, h.Length);
            mapped.Add(h with { Position = start - 1, IsReverseStrand = view.Reversed ? !h.IsReverseStrand : h.IsReverseStrand });
        }

        return mapped;
    }

    /// <summary><see cref="FindPatternHits(string, string, int, bool, CancellationToken, bool)"/> with an already parsed pattern.</summary>
    /// <param name="sequence">Subject nucleotide string.</param>
    /// <param name="pattern">The parsed pattern.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit.</param>
    /// <param name="bothStrands">Also search the reverse-complement strand.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="sequenceAmbiguity">Read IUPAC ambiguity codes in the sequence as fuzznuc does (see the string overload; default: they match nothing).</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static IReadOnlyList<PatternHit> FindPatternHits(
        string sequence, FuzznucPattern pattern, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default, bool sequenceAmbiguity = false, bool programOrder = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return SequenceTextHits(sequence, pattern, maxMismatches, bothStrands, sequenceAmbiguity, cancellationToken, programOrder);
    }

    /// <summary>
    /// Several named patterns at once (<c>fuzznuc -pattern @file</c>; see <see cref="FuzznucPattern.ParsePatternFile"/>): the hits of each pattern, keyed by its name, in the order
    /// of <paramref name="patterns"/>; a pattern without a hit is present with an empty list. <paramref name="maxMismatches"/> applies to every pattern (<c>-pmismatch</c>).
    /// </summary>
    /// <param name="sequence">Subject sequence.</param>
    /// <param name="patterns">The named patterns; names must be distinct.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit.</param>
    /// <param name="bothStrands">Also search the reverse-complement strand.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Two patterns have the same name.</exception>
    public static IReadOnlyDictionary<string, IReadOnlyList<PatternHit>> FindPatternHits(
        DnaSequence sequence, IEnumerable<NamedPattern> patterns, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patterns);
        return NamedPatternHits(sequence.Sequence, patterns, maxMismatches, bothStrands, false, cancellationToken);
    }

    /// <summary><see cref="FindPatternHits(DnaSequence, IEnumerable{NamedPattern}, int, bool, CancellationToken)"/> on a nucleotide string (non-A/C/G/T symbols match nothing).</summary>
    /// <param name="sequence">Subject nucleotide string.</param>
    /// <param name="patterns">The named patterns.</param>
    /// <param name="maxMismatches">Maximum substitutions per hit.</param>
    /// <param name="bothStrands">Also search the reverse-complement strand.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="sequenceAmbiguity">Read IUPAC ambiguity codes in the sequence as fuzznuc does (see the string overload; default: they match nothing).</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Two patterns have the same name.</exception>
    public static IReadOnlyDictionary<string, IReadOnlyList<PatternHit>> FindPatternHits(
        string sequence, IEnumerable<NamedPattern> patterns, int maxMismatches = 0, bool bothStrands = false, CancellationToken cancellationToken = default, bool sequenceAmbiguity = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patterns);
        return NamedPatternHits(sequence, patterns, maxMismatches, bothStrands, sequenceAmbiguity, cancellationToken);
    }

    private static List<PatternHit> SequenceTextHits(string sequence, FuzznucPattern pattern, int maxMismatches, bool bothStrands, bool sequenceAmbiguity, CancellationToken cancellationToken, bool programOrder = false)
    {
        string upper = EmbossSequenceReader.Clean(sequence).ToUpperInvariant();      // the fuzznuc reader drops the gap characters and every other non-letter but * (A4-3)
        return sequenceAmbiguity
            ? PatternHitsLiteral(upper.Replace('U', 'T'), pattern, maxMismatches, bothStrands, cancellationToken, programOrder)
            : PatternHitsCore(QueryText(upper), pattern, maxMismatches, bothStrands, cancellationToken);
    }

    private static Dictionary<string, IReadOnlyList<PatternHit>> NamedPatternHits(
        string text, IEnumerable<NamedPattern> patterns, int maxMismatches, bool bothStrands, bool sequenceAmbiguity, CancellationToken cancellationToken)
    {
        var result = new Dictionary<string, IReadOnlyList<PatternHit>>(StringComparer.Ordinal);
        foreach (var named in patterns)
        {
            string name = named.Name;
            if (result.ContainsKey(name))
                throw new ArgumentException($"Two patterns are named '{name}'.", nameof(patterns));
            // a pattern read from a pattern file carries its own mismatches (<mismatch=n> of its header, else -pmismatch)
            result[name] = SequenceTextHits(text, named.Pattern, named.Mismatch ?? maxMismatches, bothStrands, sequenceAmbiguity, cancellationToken);
        }

        return result;
    }

    // `text` is upper case; every symbol other than A/C/G/T is already QueryBreak (matches nothing).
    private static List<PatternHit> PatternHitsCore(string text, FuzznucPattern pattern, int maxMismatches, bool bothStrands, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        AlignmentCancellation.Poll(cancellationToken);

        var best = new Dictionary<(int Position, int Length, bool Reverse), int>();
        var expansions = pattern.ExpandWithExtras();
        bool hasRange = pattern.Elements.Any(e => e.IsRange);
        // fuzznuc picks its search class from the pattern (embPatGetTypeII): a range pattern is searched by the recursive matcher
        // patOUBrute when mismatches are allowed, otherwise by a regular expression (every other pattern reports all windows).
        bool recursiveMatcher = hasRange && maxMismatches > 0;
        bool regexClass = hasRange && !recursiveMatcher;

        Scan(expansions, pattern.AnchoredAtStart, pattern.AnchoredAtEnd, reverse: false);
        if (bothStrands)
        {
            // fuzznuc reverse-complements the sequence and runs the same pattern on it: mirror each expansion (and its extra positions) instead.
            var mirrored = expansions.Select(e => new PatternExpansion(
                DnaSequence.GetReverseComplementString(e.Iupac), e.ExtraPositions.Select(x => e.Iupac.Length - 1 - x).ToArray(), e.EndsWithEmptyElement)).ToList();
            Scan(mirrored, pattern.AnchoredAtEnd, pattern.AnchoredAtStart, reverse: true);
        }

        return OrderedPatternHits(best);

        // anchoredStart / anchoredEnd are in forward coordinates of the scanned strand's expansion (swapped for the reverse strand)
        void Scan(IReadOnlyList<PatternExpansion> list, bool anchoredStart, bool anchoredEnd, bool reverse)
        {
            var strandHits = new Dictionary<(int Position, int Length), int>();
            foreach (var expansion in list)
            {
                int length = expansion.Iupac.Length;
                foreach (var match in ApproximateMatcher.FindWithMismatches(text, expansion.Iupac, maxMismatches, iupacPattern: true, cancellationToken))
                {
                    if (anchoredStart && match.Position != 0)
                        continue;
                    if (anchoredEnd && match.Position + length != text.Length)
                        continue;
                    if (recursiveMatcher)
                    {
                        // patOUBrute: the extra repeats of a range must match exactly ...
                        if (expansion.ExtraPositions.Count > 0 && match.MismatchPositions.Any(m => expansion.ExtraPositions.Contains(m)))
                            continue;
                        // ... and a pattern with elements left over when the scanned sequence has ended fails, so a path whose last element is empty
                        // cannot end on the last base of the scanned strand (forward: position n, reverse: forward position 0).
                        if (expansion.EndsWithEmptyElement && (reverse ? match.Position == 0 : match.Position + length == text.Length))
                            continue;
                    }

                    var key = (Position: match.Position, Length: length);
                    if (!strandHits.TryGetValue(key, out int known) || match.Distance < known)
                        strandHits[key] = match.Distance;
                }
            }

            // regex class (embPatMatchFindAll with nterm): a start-anchored (<) pattern returns only the longest match at the anchor
            if (regexClass && pattern.AnchoredAtStart && strandHits.Count > 1)
            {
                int longest = strandHits.Keys.Max(k => k.Length);
                foreach (var key in strandHits.Keys.Where(k => k.Length != longest).ToList())
                    strandHits.Remove(key);
            }

            foreach (var (key, mismatches) in strandHits)
                best[(key.Position, key.Length, reverse)] = mismatches;
        }
    }

    // A3-12a: the sequence read literally as fuzznuc reads it (FuzznucPattern.Matches). `text` is upper case, U already T. The reverse strand is the IUPAC reverse complement of the
    // sequence searched with the SAME pattern (fuzznuc reverse-complements the sequence), so the anchors and the two search classes apply to the scanned strand unchanged.
    private static List<PatternHit> PatternHitsLiteral(string text, FuzznucPattern pattern, int maxMismatches, bool bothStrands, CancellationToken cancellationToken, bool programOrder = false)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        AlignmentCancellation.Poll(cancellationToken);

        // an all-optional range pattern also has the empty hit (A5-8c3b, F217): the program tries every start of the strand, so the empty unrolling is a hit at each of them (never at the end of the strand)
        var paths = pattern.ExpandPaths(includeEmpty: true);
        bool hasRange = pattern.Elements.Any(e => e.IsRange);
        bool recursiveMatcher = hasRange && maxMismatches > 0;
        bool regexClass = hasRange && !recursiveMatcher;
        if (recursiveMatcher && programOrder)
            return RecursiveMatcherHits(text);
        var best = new Dictionary<(int Position, int Length, bool Reverse), int>();

        foreach (var (position, length, mismatches) in ScanStrand(text))
            best[(position, length, false)] = mismatches;
        if (bothStrands)
            foreach (var (position, length, mismatches) in ScanStrand(DnaSequence.GetReverseComplementString(text)))
                best[(text.Length - position - length, length, true)] = mismatches;
        return OrderedPatternHits(best);

        // patOUBrute (a range pattern with mismatches, A5-8c2b F195): embPatBruteForce tries every start in turn and the recursion pushes one hit per successful unrolling (embPatPushHit), so a window several unrollings match is listed
        // several times, each with its own distance, in the order start ascending then the unrollings in the order of the recursion (ExpandPaths). The program lists the first pass in that push order and the complement pass in the REVERSED
        // push order (the list is popped without being reversed), so the hits are returned in the program's order: the strand scanned first, then the other strand (positions mapped back), instead of sorted by position.
        List<PatternHit> RecursiveMatcherHits(string forwardText)
        {
            var result = new List<PatternHit>();
            foreach (var (position, length, mismatches) in PushOrder(forwardText))
                result.Add(new PatternHit(position, length, false, mismatches));
            if (bothStrands)
            {
                var reverse = PushOrder(DnaSequence.GetReverseComplementString(forwardText));
                for (int i = reverse.Count - 1; i >= 0; i--)
                    result.Add(new PatternHit(forwardText.Length - reverse[i].Position - reverse[i].Length, reverse[i].Length, true, reverse[i].Mismatches));
            }

            return result;
        }

        List<(int Position, int Length, int Mismatches)> PushOrder(string strand)
        {
            var pushed = new List<(int, int, int)>();
            for (int start = 0; start < strand.Length; start++)
            {
                if (pattern.AnchoredAtStart && start != 0)
                    break;
                AlignmentCancellation.Poll(cancellationToken);
                foreach (var path in paths)
                {
                    int length = path.Elements.Count;
                    if (start + length > strand.Length || (pattern.AnchoredAtEnd && start + length != strand.Length))
                        continue;
                    int bad = 0;
                    bool allowed = true;
                    for (int j = 0; j < length && allowed; j++)
                    {
                        if (pattern.Matches(path.Elements[j], strand[start + j]))
                            continue;
                        allowed = ++bad <= maxMismatches && !path.ExtraPositions.Contains(j);
                    }

                    if (allowed && !(path.EndsWithEmptyElement && start + length == strand.Length))
                        pushed.Add((start, length, bad));
                }
            }

            return pushed;
        }

        List<(int Position, int Length, int Mismatches)> ScanStrand(string strand)
        {
            var strandHits = new Dictionary<(int Position, int Length), int>();
            foreach (var path in paths)
            {
                AlignmentCancellation.Poll(cancellationToken);
                int length = path.Elements.Count;
                // the empty hit: at every start of the strand (the loops of the program stop at the last residue)
                for (int start = 0; start + length <= strand.Length && (length > 0 || start < strand.Length); start++)
                {
                    if ((pattern.AnchoredAtStart && start != 0) || (pattern.AnchoredAtEnd && start + length != strand.Length))
                        continue;
                    int bad = 0;
                    bool allowed = true;
                    for (int j = 0; j < length && allowed; j++)
                    {
                        if (pattern.Matches(path.Elements[j], strand[start + j]))
                            continue;
                        // patOUBrute: the extra repeats of a range must match exactly
                        allowed = ++bad <= maxMismatches && !(recursiveMatcher && path.ExtraPositions.Contains(j));
                    }

                    // patOUBrute fails when elements are left over at the end of the scanned strand
                    if (!allowed || (recursiveMatcher && path.EndsWithEmptyElement && start + length == strand.Length))
                        continue;
                    var key = (start, length);
                    if (!strandHits.TryGetValue(key, out int known) || bad < known)
                        strandHits[key] = bad;
                }
            }

            // regex class anchored at the end ('>'): embPatMatchFindAll looks for the matches of the regular expression 'pattern$' start by start while the search position is not the end of the string; the empty string matches AT the end,
            // and is found only when no match began at the last residue (the position after it is the end: the search stops there; A5-8c3b, F217)
            if (regexClass && pattern.AnchoredAtEnd && !pattern.AnchoredAtStart && strand.Length > 0 && paths.Any(p => p.Elements.Count == 0) && !strandHits.ContainsKey((strand.Length - 1, 1)))
                strandHits[(strand.Length, 0)] = 0;

            // regex class: a start-anchored (<) range pattern returns only the longest match at the anchor
            if (regexClass && pattern.AnchoredAtStart && strandHits.Count > 1)
            {
                int longest = strandHits.Keys.Max(k => k.Length);
                foreach (var key in strandHits.Keys.Where(k => k.Length != longest).ToList())
                    strandHits.Remove(key);
            }

            return strandHits.Select(kv => (kv.Key.Position, kv.Key.Length, kv.Value)).ToList();
        }
    }

    private static List<PatternHit> OrderedPatternHits(Dictionary<(int Position, int Length, bool Reverse), int> best)
        => best
            .OrderBy(kv => kv.Key.Position).ThenBy(kv => kv.Key.Length).ThenBy(kv => kv.Key.Reverse)
            .Select(kv => new PatternHit(kv.Key.Position, kv.Key.Length, kv.Key.Reverse, kv.Value))
            .ToList();

    #endregion

    #region String overloads of the motif finders (A2-35)

    /// <summary>
    /// <see cref="FindMotif(DnaSequence, string)"/> on a nucleotide <b>string</b> (N runs, IUPAC codes, lower case). Case-insensitive; a symbol other than A/C/G/T in the
    /// sequence matches nothing, so no window containing one is an occurrence (and a motif that contains one has none, as on a <see cref="DnaSequence"/>). For an A/C/G/T string
    /// the result equals the <see cref="DnaSequence"/> overload.
    /// </summary>
    /// <param name="sequence">Subject nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="motif">Exact motif; a null or empty motif yields an empty list.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IReadOnlyList<int> FindMotif(string sequence, string motif)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (string.IsNullOrEmpty(motif))
            return Array.Empty<int>();
        return KmerAnalyzer.FindKmerPositions(QueryText(sequence.ToUpperInvariant()), motif).ToList();
    }

    /// <summary>
    /// <see cref="FindKnownMotifs(DnaSequence, IEnumerable{string})"/> on a nucleotide <b>string</b>; the rule of <see cref="FindMotif(string, string)"/> (non-A/C/G/T never match).
    /// Many-motif sets use one linear scan per motif (the suffix-tree route needs a <see cref="DnaSequence"/>).
    /// </summary>
    /// <param name="sequence">Subject nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="motifs">Motifs to search.</param>
    /// <exception cref="ArgumentNullException">An argument is <c>null</c>.</exception>
    public static Dictionary<string, IReadOnlyList<int>> FindKnownMotifs(string sequence, IEnumerable<string> motifs)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(motifs);
        var found = SearchKnownMotifs(QueryText(sequence.ToUpperInvariant()), null, motifs, bothStrands: false, allowIupac: false);
        var result = new Dictionary<string, IReadOnlyList<int>>(found.Count);
        foreach (var (key, (forward, _)) in found)
            result[key] = forward;
        return result;
    }

    /// <summary>
    /// <see cref="FindKnownMotifHits(DnaSequence, IEnumerable{string}, bool, bool, int, CancellationToken)"/> on a nucleotide <b>string</b>: strands, IUPAC motifs and the mismatch budget
    /// work as there, but a sequence symbol other than A/C/G/T matches no motif symbol (not even <c>N</c> with <paramref name="allowIupac"/>) and costs a mismatch when
    /// <paramref name="maxMismatches"/> &gt; 0 — the sequence side is literal, as in the canonical motif finders.
    /// </summary>
    /// <param name="sequence">Subject nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="motifs">Motifs to search.</param>
    /// <param name="bothStrands">Also report occurrences on the reverse strand.</param>
    /// <param name="allowIupac">Interpret IUPAC ambiguity codes in the motifs.</param>
    /// <param name="maxMismatches">Maximum number of mismatching positions per hit.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">An argument is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static Dictionary<string, IReadOnlyList<KnownMotifHit>> FindKnownMotifHits(
        string sequence,
        IEnumerable<string> motifs,
        bool bothStrands = false,
        bool allowIupac = false,
        int maxMismatches = 0,
        CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(motifs);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        AlignmentCancellation.Poll(cancellationToken);

        string text = QueryText(sequence.ToUpperInvariant());
        if (maxMismatches > 0)
            return FindKnownMotifHitsApproximate(text, motifs, bothStrands, allowIupac, maxMismatches, cancellationToken);

        var found = SearchKnownMotifs(text, null, motifs, bothStrands, allowIupac, cancellationToken);
        var result = new Dictionary<string, IReadOnlyList<KnownMotifHit>>(found.Count);
        foreach (var (key, (forward, reverse)) in found)
            result[key] = MergeStrandHits(forward, reverse);
        return result;
    }

    #endregion
}

/// <summary>
/// A hit of <see cref="EmbossMotifPrograms.FindPatternHits(DnaSequence, string, int, bool, CancellationToken)"/>: the 0-based forward-strand start and the length of the matching window,
/// the strand it matched on and its Hamming distance to the pattern (fuzznuc's <c>*mismatch</c> note, 0 for an exact hit).
/// </summary>
/// <param name="Position">0-based start of the window on the forward strand.</param>
/// <param name="Length">Window length (a range such as <c>A(2,4)</c> yields hits of several lengths).</param>
/// <param name="IsReverseStrand">The window matches the reverse complement of the pattern.</param>
/// <param name="Mismatches">Substituted positions, ≤ the allowed maximum.</param>
public readonly record struct PatternHit(int Position, int Length, bool IsReverseStrand, int Mismatches = 0)
{
    /// <summary>Exclusive end of the window (0-based).</summary>
    public int End => Position + Length;
}
