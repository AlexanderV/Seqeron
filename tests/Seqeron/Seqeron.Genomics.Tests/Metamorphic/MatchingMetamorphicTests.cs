using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Matching area (exact pattern / motif matching).
///
/// Each test encodes a metamorphic relation (MR) — a property that relates the
/// outputs of several executions under input transformations, without needing a
/// hardcoded oracle. The relations are derived from the exact-matching
/// *definition*, not from the current implementation's output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PAT-EXACT-001 — exact pattern/motif matching (Matching).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 8.
/// Relations (row 8):
///   • SHIFT: prepend flank shifts positions by |flank|.
///   • COMP:  exact ⊆ hamming(maxDist=0)  (in fact set equality).
///   • INV:   duplicate (S+S) → doubled count, second block shifted by |S|.
///
/// Source (exact-matching definition):
///   An exact occurrence of pattern P in text T is every 0-based start position i
///   with T[i .. i+|P|-1] = P. ALL occurrences are reported, INCLUDING overlapping
///   ones (e.g. text "GATATATGCATATACTT", pattern "ATAT" → positions {1, 3, 9},
///   where 1 and 3 overlap). The window scan is exhaustive over i ∈ [0, |T|-|P|].
///   — docs/algorithms/Pattern_Matching/Exact_Pattern_Search.md §2.1–2.2, §7.1
///     [Gusfield 1997; Ukkonen 1995; Rosalind SUBS].
///
/// API surface under test:
///   • MotifFinder.FindExactMotif(DnaSequence, string) → IEnumerable&lt;int&gt;
///     (sorted 0-based start positions; uppercases the motif; empty for empty motif).
///     — src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs.
///   • ApproximateMatcher.FindWithMismatches(string seq, string pattern, int maxMismatches)
///     → IEnumerable&lt;ApproximateMatchResult&gt; (Hamming). With maxMismatches = 0 it is,
///     by the doc's own edge-case table, "equivalent to exact matching" — used for COMP.
///     — src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs.
///
/// Overlap policy (confirmed): BOTH surfaces report overlapping occurrences. The
/// exact matcher collects every suffix-tree leaf under the matched locus, and the
/// Hamming matcher scans every length-|P| window, so neither de-overlaps. This is
/// load-bearing for the SHIFT and INV relations below (they count overlapping hits).
///
/// Why these relations hold (from the definition above):
///   • SHIFT: occurrences depend only on the local substring T[i..i+|P|-1]. Prepending
///     a flank F that introduces no occurrence overlapping the F|S junction relabels
///     every interior position i → i+|F| and adds none, so the position SET maps
///     bijectively under +|F| and the count is preserved exactly.
///   • COMP: a Hamming-distance-0 window is, by definition, an exact match and vice
///     versa, so the two position sets are equal (not merely subset).
///   • INV: in S+S every occurrence of P in the first |S| characters recurs |S| later
///     in the second copy; choosing S so the junction spawns no extra occurrence pins
///     the count at exactly 2k with the second block shifted by |S|.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MatchingMetamorphicTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed so random inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260619);

    /// <summary>Generates a random DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>0-based exact-match start positions via the production matcher (sorted).</summary>
    private static List<int> ExactPositions(string text, string pattern) =>
        MotifFinder.FindExactMotif(new DnaSequence(text), pattern).OrderBy(p => p).ToList();

    /// <summary>0-based Hamming-distance-0 (≡ exact) start positions, sorted.</summary>
    private static List<int> HammingZeroPositions(string text, string pattern) =>
        ApproximateMatcher.FindWithMismatches(text, pattern, 0)
            .Select(r => r.Position).OrderBy(p => p).ToList();

    /// <summary>
    /// Naive reference scan of every length-|P| window — the textbook definition of an
    /// exact occurrence. Used as a definition-derived oracle, independent of either
    /// production matcher, so the relations are pinned to THEORY, not to code output.
    /// </summary>
    private static List<int> NaivePositions(string text, string pattern)
    {
        var result = new List<int>();
        if (pattern.Length == 0 || pattern.Length > text.Length) return result;
        for (int i = 0; i <= text.Length - pattern.Length; i++)
        {
            if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern))
                result.Add(i);
        }
        return result;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-EXACT-001 — exact pattern/motif matching
    // ═══════════════════════════════════════════════════════════════════

    #region MR1: SHIFT — prepending a flank shifts positions by |flank|

    /// <summary>
    /// MR1: for a flank F that creates no occurrence spanning the F|S junction,
    /// positions(F + S) == { p + |F| : p ∈ positions(S) }, and the count is preserved
    /// exactly. Occurrences depend only on the local window T[i..i+|P|-1]; a leading
    /// flank relabels every interior start i → i+|F| and (by construction) adds none.
    ///
    /// Junction safety: each flank is a homopolymer run of a base that does NOT begin
    /// the pattern, so no window straddling the boundary can equal the pattern. We also
    /// assert, via the naive reference scan, that no NEW occurrence appears at the seam.
    /// </summary>
    [Test]
    public void FindExactMotif_PrependNonDisruptingFlank_ShiftsPositionsByFlankLength()
    {
        // (sequence, pattern) cases — overlapping occurrences included deliberately.
        var cases = new[]
        {
            (seq: "ATGCATGCATGC", pat: "ATGC"),
            (seq: "GATATATGCATATACTT", pat: "ATAT"), // overlapping occurrences at 1,3,9
            (seq: "AAAAAA", pat: "AA"),               // dense overlapping occurrences
            (seq: RandomDna(120), pat: "GTAC"),
            (seq: RandomDna(200), pat: "CG"),
        };

        // Flanks: homopolymer runs whose base never starts the corresponding pattern,
        // so no junction-spanning window can equal the pattern.
        foreach (var (seq, pat) in cases)
        {
            char flankBase = "ACGT".First(b => b != pat[0]);
            foreach (int flankLen in new[] { 1, 3, 7, 25 })
            {
                string flank = new(flankBase, flankLen);
                string shifted = flank + seq;

                var basePositions = ExactPositions(seq, pat);
                var shiftedPositions = ExactPositions(shifted, pat);

                // No new occurrence introduced (definition-derived check on the seam too).
                NaivePositions(shifted, pat).Should().BeEquivalentTo(
                    basePositions.Select(p => p + flankLen),
                    because: $"flank '{flank}' starts with '{flankBase}' (≠ pattern start '{pat[0]}'), so the F|S junction creates no new occurrence of '{pat}'");

                shiftedPositions.Should().Equal(basePositions.Select(p => p + flankLen),
                    because: $"prepending a {flankLen}-base flank relabels every occurrence of '{pat}' in '{seq}' by +{flankLen}");

                shiftedPositions.Should().HaveCount(basePositions.Count,
                    because: $"a non-disrupting flank neither adds nor removes occurrences of '{pat}', so the count is preserved");
            }
        }
    }

    /// <summary>
    /// MR1-b: the SHIFT also holds end-to-end against the textbook naive scan — both the
    /// production matcher and the definition agree on the +|F| relabelling. This guards
    /// against the matcher silently de-overlapping (which would break the bijection).
    /// </summary>
    [Test]
    public void FindExactMotif_AgreesWithNaiveScan_AndPreservesOverlaps()
    {
        var cases = new[]
        {
            (seq: "GATATATGCATATACTT", pat: "ATAT"),
            (seq: "AAAAAAAA", pat: "AAA"),
            (seq: RandomDna(150), pat: "AT"),
            (seq: RandomDna(150), pat: "ACG"),
        };

        foreach (var (seq, pat) in cases)
        {
            ExactPositions(seq, pat).Should().Equal(NaivePositions(seq, pat),
                because: $"FindExactMotif reports every (overlapping) start window where '{pat}' occurs in '{seq}', matching the definition");
        }
    }

    #endregion

    #region MR2: COMP — exact ⊆ hamming(maxDist=0), and in fact set equality

    /// <summary>
    /// MR2: the exact-match position set equals the Hamming-distance-0 position set.
    /// A 0-mismatch window IS an exact occurrence and vice versa, so:
    ///   exact(S, P) ⊆ hamming0(S, P)  AND  hamming0(S, P) ⊆ exact(S, P).
    /// Both surfaces are asserted to coincide with the naive definition scan, pinning
    /// equality to theory rather than to either implementation.
    /// </summary>
    [Test]
    public void ExactMatch_EqualsHammingDistanceZero()
    {
        var cases = new[]
        {
            (seq: "ATGCATGCATGC", pat: "ATGC"),
            (seq: "GATATATGCATATACTT", pat: "ATAT"),
            (seq: "AAAAAA", pat: "AA"),
            (seq: "ACGTACGTACGTACGT", pat: "TTTTT"), // pattern absent — both sets must be empty
            (seq: RandomDna(180), pat: "GAT"),
            (seq: RandomDna(180), pat: "CGCG"),
        };

        foreach (var (seq, pat) in cases)
        {
            var exact = ExactPositions(seq, pat);
            var hamming0 = HammingZeroPositions(seq, pat);
            var reference = NaivePositions(seq, pat);

            exact.Should().BeSubsetOf(hamming0,
                because: $"every exact occurrence of '{pat}' is a window at Hamming distance 0, so it is a 0-mismatch match");
            hamming0.Should().BeSubsetOf(exact,
                because: $"a Hamming-distance-0 window of '{pat}' is by definition an exact occurrence");

            exact.Should().Equal(reference,
                because: $"the exact matcher must reproduce the definition's occurrence set for '{pat}' in '{seq}'");
            hamming0.Should().Equal(reference,
                because: $"maxMismatches = 0 is equivalent to exact matching, so hamming0 must equal the definition set for '{pat}'");
        }
    }

    /// <summary>
    /// MR2-b: every reported 0-mismatch result genuinely has Distance == 0 and its matched
    /// window equals the (uppercased) pattern — the Hamming surface does not smuggle in
    /// non-exact windows under the maxDist = 0 setting.
    /// </summary>
    [Test]
    public void HammingDistanceZero_ResultsAreGenuinelyExact()
    {
        foreach (var (seq, pat) in new[]
        {
            (seq: "ATGCATGCATGC", pat: "ATGC"),
            (seq: RandomDna(120), pat: "GTA"),
            (seq: RandomDna(120), pat: "TT"),
        })
        {
            foreach (var r in ApproximateMatcher.FindWithMismatches(seq, pat, 0))
            {
                r.Distance.Should().Be(0,
                    because: "maxMismatches = 0 admits only zero-mismatch windows");
                r.MatchedSequence.Should().Be(pat.ToUpperInvariant(),
                    because: $"a 0-mismatch window must equal the pattern '{pat}' exactly (case-folded)");
            }
        }
    }

    #endregion

    #region MR3: INV — duplicating the sequence doubles the occurrence count

    /// <summary>
    /// MR3: for a pattern occurring k times in S, the concatenation S+S contains exactly
    /// 2k occurrences when the S|S junction creates none, with the second block's positions
    /// shifted by |S|. Formally:
    ///   positions(S + S) == positions(S) ∪ { p + |S| : p ∈ positions(S) },  count = 2k.
    ///
    /// Junction safety: each S ends and begins with bases chosen so that no window crossing
    /// the |S| boundary equals the pattern. We additionally assert, via the naive reference
    /// scan, that the seam introduces no extra occurrence (so the count is EXACTLY 2k, not ≥).
    /// </summary>
    [Test]
    public void FindExactMotif_DuplicatedSequence_DoublesCount_WithSecondBlockShifted()
    {
        // Each S is chosen so the S|S junction cannot complete the pattern:
        // S ends with a base, S begins with a base, and the boundary substring ≠ pattern.
        var cases = new[]
        {
            (seq: "TATGCATGCAT", pat: "ATGC"),   // ends 'T', begins 'T'; no "ATGC" across seam
            (seq: "TACGTACGTA",  pat: "ACGT"),
            (seq: "CATATATAC",   pat: "ATAT"),    // internal overlapping occurrences
            (seq: "GATTACAG",    pat: "AT"),
        };

        foreach (var (seq, pat) in cases)
        {
            var single = ExactPositions(seq, pat);
            int k = single.Count;
            k.Should().BeGreaterThan(0,
                because: $"the test is only meaningful when '{pat}' occurs in '{seq}'");

            string doubled = seq + seq;
            var doubledPositions = ExactPositions(doubled, pat);

            var expected = single
                .Concat(single.Select(p => p + seq.Length))
                .OrderBy(p => p)
                .ToList();

            // Definition-derived guard: the seam introduces no extra occurrence.
            NaivePositions(doubled, pat).Should().Equal(expected,
                because: $"the S|S junction of '{seq}' cannot complete '{pat}', so S+S has exactly the two shifted blocks");

            doubledPositions.Should().Equal(expected,
                because: $"each of the {k} occurrences of '{pat}' recurs |S|={seq.Length} later in S+S");
            doubledPositions.Should().HaveCount(2 * k,
                because: $"with a non-disrupting junction, duplicating '{seq}' exactly doubles the occurrence count of '{pat}'");
        }
    }

    /// <summary>
    /// MR3-b: the doubling is a lower bound in GENERAL (≥ 2k) even when the junction may
    /// add seam-spanning occurrences. Concatenating S+S always contains the two original
    /// blocks of occurrences, so count(S+S) ≥ 2·count(S) for any S — a weaker but
    /// unconditionally-true companion relation. Verified on random sequences (fixed seed)
    /// where junction effects are not controlled.
    /// </summary>
    [Test]
    public void FindExactMotif_DuplicatedSequence_CountIsAtLeastDoubled_General()
    {
        for (int t = 0; t < 8; t++)
        {
            string seq = RandomDna(40 + Rng.Next(60));
            // A short pattern guarantees several occurrences in random DNA.
            string pat = RandomDna(2);

            int k = ExactPositions(seq, pat).Count;
            int doubled = ExactPositions(seq + seq, pat).Count;

            doubled.Should().BeGreaterThanOrEqualTo(2 * k,
                because: $"S+S contains both copies of every occurrence of '{pat}' in '{seq}', so the count is at least doubled (junction may add more)");
        }
    }

    #endregion
}
