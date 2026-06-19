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
/// Units: PAT-EXACT-001 — exact pattern/motif matching (Matching);
///        PAT-APPROX-001 — approximate (Hamming-distance) matching (Matching).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, rows 8 and 9.
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

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-APPROX-001 — approximate (Hamming-distance) matching
    // ═══════════════════════════════════════════════════════════════════
    //
    //  Relations (checklist row 9):
    //    • SYM:  hamming(a,b) = hamming(b,a)            (metric symmetry; corollary d(a,a)=0)
    //    • MON:  higher maxDist → ⊇ matches             (positions monotone in the threshold)
    //    • COMP: exact ⊆ approx(d=0)                    (in fact set equality)
    //
    //  Source (Hamming-distance definition):
    //    d_H(s,t) = |{ i : s[i] ≠ t[i], 0 ≤ i < n }| for equal-length s,t (case-insensitive).
    //    Approximate matching with at most k mismatches compares the pattern to every
    //    equal-length window of the text and reports windows with d_H ≤ k.
    //    — docs/algorithms/Pattern_Matching/Approximate_Matching_Hamming.md §2.2, §2.4
    //      (INV-01 non-negativity, INV-02 d=0 iff identical, INV-03 symmetry); Hamming 1950;
    //      Rosalind HAMM; Navarro 2001; Compeau & Pevzner ch.1 (ROSALIND BA1H/BA1I).
    //
    //  APIs under test:
    //    • ApproximateMatcher.HammingDistance(string, string)
    //        → int; throws ArgumentException on unequal lengths; case-insensitive.
    //    • ApproximateMatcher.FindWithMismatches(string seq, string pat, int maxMismatches)
    //        → IEnumerable<ApproximateMatchResult>; equal-length window scan, d_H ≤ k.
    //    • MotifFinder.FindExactMotif(DnaSequence, string)  (via ExactPositions, for COMP).
    //
    //  Differentiation from PAT-APPROX-002 (already done in MetamorphicTests.cs): that unit
    //  covers the EDIT-distance path — FindWithEdits/EditDistance maxEdits-monotonicity,
    //  EditDistance symmetry/non-negativity, exact ⊆ approximate under EDIT distance. THIS
    //  unit exercises the HAMMING/substitution path exclusively: HammingDistance(string,string)
    //  symmetry and the FindWithMismatches mismatch-threshold monotonicity — disjoint surfaces.

    #region MR4: SYM — Hamming distance is symmetric (and d(a,a)=0)

    /// <summary>
    /// MR4: <c>HammingDistance(a, b) == HammingDistance(b, a)</c> for every equal-length
    /// pair — positionwise mismatch counting is order-independent (doc INV-03). The metric
    /// corollary <c>d(a, a) == 0</c> (doc INV-02) is asserted alongside, and a hand-computed
    /// anchor (Rosalind HAMM sample, expected 7) pins the absolute value to theory.
    /// </summary>
    [Test]
    public void HammingDistance_SwapArguments_IsSymmetric()
    {
        // Fixed anchors with a hand-known distance, plus randomized equal-length pairs.
        var fixedPairs = new[]
        {
            (a: "ACGT", b: "ACGT", d: 0),
            (a: "ACGT", b: "TGCA", d: 4),
            (a: "AAAA", b: "AAAT", d: 1),
            // Rosalind HAMM sample — independently known distance = 7.
            (a: "GAGCCTACTAACGGGAT", b: "CATCGTAATGACGGCCT", d: 7),
        };

        foreach (var (a, b, d) in fixedPairs)
        {
            int ab = ApproximateMatcher.HammingDistance(a, b);
            int ba = ApproximateMatcher.HammingDistance(b, a);

            ab.Should().Be(d,
                because: $"d_H('{a}','{b}') counts the positions where the two strings differ, which is {d}");
            ab.Should().Be(ba,
                because: $"positionwise mismatch counting is order-independent, so d_H('{a}','{b}') = d_H('{b}','{a}') (doc INV-03)");
        }

        // Randomized equal-length pairs (fixed seed) — symmetry, non-negativity, and the
        // self-distance corollary must hold for every pair.
        for (int t = 0; t < 40; t++)
        {
            int len = 1 + Rng.Next(60);
            string a = RandomDna(len);
            string b = RandomDna(len);

            int ab = ApproximateMatcher.HammingDistance(a, b);
            int ba = ApproximateMatcher.HammingDistance(b, a);

            ab.Should().Be(ba,
                because: "Hamming distance is symmetric: swapping the operands cannot change the per-position mismatch count");
            ab.Should().BeGreaterThanOrEqualTo(0,
                because: "a Hamming distance is a count of mismatching positions, hence never negative (doc INV-01)");
            ApproximateMatcher.HammingDistance(a, a).Should().Be(0,
                because: "a string has zero mismatches against itself (doc INV-02, identity corollary)");
        }
    }

    /// <summary>
    /// MR4-b: symmetry is robust to the case-insensitive comparison — uppercasing one
    /// operand does not perturb the distance, and the swap remains symmetric. Guards
    /// against a case-folding asymmetry between the two argument positions.
    /// </summary>
    [Test]
    public void HammingDistance_CaseFolding_PreservesSymmetry()
    {
        foreach (var (a, b) in new[]
        {
            (a: "acgt", b: "ACGT"),
            (a: "AcGtAcGt", b: "acGTacGT"),
            (a: "GATTACA", b: "gattaca"),
        })
        {
            int ab = ApproximateMatcher.HammingDistance(a, b);
            int ba = ApproximateMatcher.HammingDistance(b, a);

            ab.Should().Be(ba,
                because: "case-insensitive comparison must treat both argument positions identically, so the distance stays symmetric");
            ApproximateMatcher.HammingDistance(a.ToUpperInvariant(), b.ToUpperInvariant())
                .Should().Be(ab,
                because: "uppercasing the inputs cannot change a case-insensitive Hamming distance");
        }
    }

    /// <summary>
    /// MR4-c: unequal-length operands are outside the Hamming metric's domain, so the
    /// distance is undefined and the implementation must reject them (doc §3.3). This
    /// fixes the equal-length precondition that the SYM relation relies on.
    /// </summary>
    [Test]
    public void HammingDistance_UnequalLengths_Throws()
    {
        Action act = () => ApproximateMatcher.HammingDistance("ACGT", "ACG");
        act.Should().Throw<ArgumentException>(
            because: "Hamming distance is defined only for equal-length strings, so a length mismatch is rejected");
    }

    #endregion

    #region MR5: MON — a higher mismatch threshold yields a superset of matches

    /// <summary>
    /// MR5: the approximate-occurrence position set is monotone non-decreasing in the
    /// mismatch threshold. A window with Hamming distance ≤ k also has distance ≤ k+1,
    /// so for every k:
    ///   positions(seq, pat, k) ⊆ positions(seq, pat, k+1),
    /// and |positions(·, k)| ≤ |positions(·, k+1)|. We build a chain k = 0,1,2,…,|pat|
    /// and assert subset + non-decreasing count along the whole chain. (HAMMING/mismatch
    /// matcher — distinct from the edit-distance MON of PAT-APPROX-002.)
    /// </summary>
    [Test]
    public void FindWithMismatches_HigherThreshold_YieldsSupersetOfPositions()
    {
        var cases = new[]
        {
            (seq: "ATGCATGCATGC", pat: "ATGC"),
            (seq: "GATATATGCATATACTT", pat: "ATAT"),
            (seq: "AAAAAAAA", pat: "AAA"),
            (seq: RandomDna(200), pat: "GTAC"),
            (seq: RandomDna(200), pat: "CGT"),
            (seq: RandomDna(150), pat: RandomDna(5)),
        };

        foreach (var (seq, pat) in cases)
        {
            // Chain of thresholds 0..|pat| (at |pat| every window qualifies).
            for (int k = 0; k < pat.Length; k++)
            {
                var lower = MismatchPositions(seq, pat, k);
                var higher = MismatchPositions(seq, pat, k + 1);

                lower.Should().BeSubsetOf(higher,
                    because: $"a window with ≤{k} mismatches of '{pat}' also has ≤{k + 1}, so raising maxMismatches never drops a match");
                higher.Count.Should().BeGreaterThanOrEqualTo(lower.Count,
                    because: $"the occurrence count of '{pat}' is non-decreasing in the mismatch threshold (k={k} → k={k + 1})");
            }
        }
    }

    /// <summary>
    /// MR5-b: at <c>maxMismatches = |pat|</c> every equal-length window qualifies, so the
    /// match count equals the number of windows, |seq| − |pat| + 1 — the saturation point
    /// of the monotone chain (a window can differ in at most |pat| positions). Anchored to
    /// an independent naive Hamming scan to pin the set, not just the count.
    /// </summary>
    [Test]
    public void FindWithMismatches_ThresholdAtPatternLength_MatchesEveryWindow()
    {
        foreach (var (seq, pat) in new[]
        {
            (seq: "ATGCATGC", pat: "ACGT"),
            (seq: RandomDna(80), pat: "GAT"),
            (seq: RandomDna(80), pat: "TTTT"),
        })
        {
            var saturated = MismatchPositions(seq, pat, pat.Length);
            var allWindows = Enumerable.Range(0, seq.Length - pat.Length + 1).ToList();

            saturated.Should().Equal(allWindows,
                because: $"with maxMismatches = |pat| = {pat.Length} every length-{pat.Length} window of '{seq}' qualifies, so all {allWindows.Count} start positions are reported");

            // Independent oracle: a naive Hamming scan at the same threshold agrees.
            saturated.Should().Equal(NaiveHammingPositions(seq, pat, pat.Length),
                because: "the production mismatch scan must reproduce the definition's occurrence set at the saturating threshold");
        }
    }

    /// <summary>
    /// MR5-c: at every threshold along the chain the production matcher agrees with an
    /// independent naive Hamming scan, AND each reported result's <c>Distance</c> is the
    /// true Hamming distance of its window and is ≤ the threshold — so the monotone chain
    /// is anchored to theory at each link, not merely internally consistent.
    /// </summary>
    [Test]
    public void FindWithMismatches_AtEachThreshold_AgreesWithNaiveHammingScan()
    {
        var cases = new[]
        {
            (seq: "GATATATGCATATACTT", pat: "ATAT"),
            (seq: RandomDna(120), pat: "ACG"),
            (seq: RandomDna(120), pat: RandomDna(4)),
        };

        foreach (var (seq, pat) in cases)
        {
            for (int k = 0; k <= pat.Length; k++)
            {
                MismatchPositions(seq, pat, k).Should().Equal(NaiveHammingPositions(seq, pat, k),
                    because: $"FindWithMismatches must report exactly the windows of '{pat}' in '{seq}' with Hamming distance ≤ {k}");

                foreach (var r in ApproximateMatcher.FindWithMismatches(seq, pat, k))
                {
                    r.Distance.Should().BeLessThanOrEqualTo(k,
                        because: $"a reported match must be within the requested threshold {k}");
                    r.Distance.Should().Be(
                        ApproximateMatcher.HammingDistance(pat, seq.Substring(r.Position, pat.Length)),
                        because: $"the reported Distance must equal the true Hamming distance of the window at position {r.Position}");
                }
            }
        }
    }

    #endregion

    #region MR6: COMP — exact ⊆ approx(d=0) (and in fact set equality)

    /// <summary>
    /// MR6: every exact-match position appears among the 0-mismatch approximate matches,
    /// and the two sets are equal — a 0-mismatch Hamming window IS an exact occurrence and
    /// vice versa. Both surfaces are anchored to the naive definition scan:
    ///   exact(S,P) ⊆ approx(S,P,0)  AND  approx(S,P,0) ⊆ exact(S,P),  i.e. set equality.
    /// </summary>
    [Test]
    public void ExactMatch_IsSubsetOf_ApproxDistanceZero_AndEqual()
    {
        var cases = new[]
        {
            (seq: "ATGCATGCATGC", pat: "ATGC"),
            (seq: "GATATATGCATATACTT", pat: "ATAT"),
            (seq: "AAAAAA", pat: "AA"),
            (seq: "ACGTACGTACGTACGT", pat: "TTTTT"), // absent — both sets empty
            (seq: RandomDna(180), pat: "GAT"),
            (seq: RandomDna(180), pat: RandomDna(4)),
        };

        foreach (var (seq, pat) in cases)
        {
            var exact = ExactPositions(seq, pat);
            var approx0 = MismatchPositions(seq, pat, 0);
            var reference = NaivePositions(seq, pat);

            exact.Should().BeSubsetOf(approx0,
                because: $"every exact occurrence of '{pat}' is a window at Hamming distance 0, hence a 0-mismatch match");
            approx0.Should().BeSubsetOf(exact,
                because: $"a 0-mismatch (Hamming-distance-0) window of '{pat}' is by definition an exact occurrence");
            approx0.Should().Equal(exact,
                because: $"exact matching and approx(d=0) describe the same occurrence set for '{pat}'");

            exact.Should().Equal(reference,
                because: $"the exact matcher must reproduce the definition's occurrence set for '{pat}' in '{seq}'");
            approx0.Should().Equal(reference,
                because: $"maxMismatches = 0 is equivalent to exact matching, so approx(d=0) must equal the definition set for '{pat}'");
        }
    }

    #endregion

    #region Helpers — Hamming (approximate) matching

    /// <summary>0-based start positions of approximate (Hamming) matches with ≤ k mismatches, sorted.</summary>
    private static List<int> MismatchPositions(string text, string pattern, int maxMismatches) =>
        ApproximateMatcher.FindWithMismatches(text, pattern, maxMismatches)
            .Select(r => r.Position).OrderBy(p => p).ToList();

    /// <summary>
    /// Naive reference scan of every equal-length window, reporting starts whose Hamming
    /// distance to the (uppercased) pattern is ≤ maxMismatches. The textbook definition of
    /// an approximate occurrence — independent of the production matcher, so the MON/COMP
    /// relations are pinned to THEORY rather than to code output.
    /// </summary>
    private static List<int> NaiveHammingPositions(string text, string pattern, int maxMismatches)
    {
        var result = new List<int>();
        if (pattern.Length == 0 || pattern.Length > text.Length) return result;
        string t = text.ToUpperInvariant();
        string p = pattern.ToUpperInvariant();
        for (int i = 0; i <= t.Length - p.Length; i++)
        {
            int d = 0;
            for (int j = 0; j < p.Length; j++)
                if (t[i + j] != p[j]) d++;
            if (d <= maxMismatches) result.Add(i);
        }
        return result;
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MOTIF-CONS-001 — IUPAC consensus from aligned rows (Matching).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 169.
    //
    // API under test (MotifFinder.GenerateConsensus):
    //   Per column, bases whose fraction exceeds the IUPAC inclusion threshold are merged into a
    //   degenerate IUPAC code.
    //
    // Relations (derived from the fraction-based column rule, NOT from output):
    //   • INV  (row order independent): each column tally is a multiset, so reordering the rows yields
    //          the same consensus.
    //   • INV  (duplicating the rows preserves consensus): replicating every row scales counts and the
    //          total equally, so every column fraction — and the consensus — is unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly string[] ConsensusRows =
    {
        "ACGTACGT",
        "ACGTACGT",
        "ACGAACGT",
        "TCGTACGA",
    };

    #region MOTIF-CONS-001 INV — consensus is independent of row order

    [Test]
    [Description("INV: each column's base tally is a multiset, so reversing (permuting) the row order yields the identical consensus.")]
    public void Consensus_RowOrder_Invariant()
    {
        string original = MotifFinder.GenerateConsensus(ConsensusRows);
        MotifFinder.GenerateConsensus(ConsensusRows.Reverse().ToList()).Should().Be(original,
            because: "the per-column IUPAC code depends only on the base multiset, not the row order");
    }

    #endregion

    #region MOTIF-CONS-001 INV — duplicating the rows preserves the consensus

    [Test]
    [Description("INV: replicating every row scales each column's counts and the total equally, so all fractions — and therefore the IUPAC consensus — are unchanged.")]
    public void Consensus_DuplicateRows_Unchanged()
    {
        string original = MotifFinder.GenerateConsensus(ConsensusRows);
        var duplicated = ConsensusRows.Concat(ConsensusRows).ToList();

        MotifFinder.GenerateConsensus(duplicated).Should().Be(original,
            because: "doubling every row leaves every per-column base fraction identical");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MOTIF-DISCOVER-001 — overrepresented k-mer discovery (Matching).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 170.
    //
    // API under test (MotifFinder.DiscoverMotifs):
    //   Reports k-mers whose occurrence count ≥ minCount, with their positions and O/E enrichment.
    //
    // Relations (derived from the count threshold, NOT from output):
    //   • MON  (lower support ⇒ superset): lowering minCount admits every previously-reported k-mer
    //          plus lower-count ones, so the discovered set grows.
    //   • SHIFT (prepend flank shifts positions): prepending a flank that adds no occurrence of a
    //          motif shifts that motif's positions by the flank length.
    // ───────────────────────────────────────────────────────────────────────────

    #region MOTIF-DISCOVER-001 MON — lower support yields a superset

    [Test]
    [Description("MON: lowering minCount keeps every higher-count k-mer and admits lower-count ones, so the discovered motif set at the lower support is a superset.")]
    public void DiscoverMotifs_LowerSupport_Superset()
    {
        var seq = new DnaSequence("ATGATGATG"); // ATG×3, TGA×2, GAT×2 at k=3
        var support3 = MotifFinder.DiscoverMotifs(seq, k: 3, minCount: 3).Select(m => m.Sequence).ToHashSet();
        var support2 = MotifFinder.DiscoverMotifs(seq, k: 3, minCount: 2).Select(m => m.Sequence).ToHashSet();

        support3.IsSubsetOf(support2).Should().BeTrue(because: "every k-mer with count ≥3 also has count ≥2");
        support2.Count.Should().BeGreaterThan(support3.Count, because: "TGA and GAT (count 2) are admitted only at the lower support");
    }

    #endregion

    #region MOTIF-DISCOVER-001 SHIFT — a prepended flank shifts motif positions

    [Test]
    [Description("SHIFT: prepending a flank that introduces no new ATG occurrence shifts every ATG position by the flank length.")]
    public void DiscoverMotifs_PrependFlank_ShiftsPositions()
    {
        const string seq = "ATGATGATG";
        var original = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: 3, minCount: 2)
            .First(m => m.Sequence == "ATG").Positions.ToList();

        foreach (var flank in new[] { "CC", "GCCG" }) // contain no ATG and form none at the junction
        {
            var shifted = MotifFinder.DiscoverMotifs(new DnaSequence(flank + seq), k: 3, minCount: 2)
                .First(m => m.Sequence == "ATG").Positions.ToList();
            shifted.Should().Equal(original.Select(p => p + flank.Length),
                because: $"the {flank.Length}-base flank relocates every ATG occurrence by {flank.Length}");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MOTIF-GENERATE-001 — position weight matrix generation (Matching).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 171.
    //
    // API under test (MotifFinder.CreatePwm):
    //   Builds a log-odds PWM from equal-length sequences using per-column base counts.
    //
    // Relations (derived from column counting, NOT from output):
    //   • INV  (row order independent): the PWM is a function of each column's base counts, so
    //          reordering the input sequences yields the identical matrix and consensus.
    // ───────────────────────────────────────────────────────────────────────────

    #region MOTIF-GENERATE-001 INV — the PWM is independent of sequence order

    [Test]
    [Description("INV: each PWM column is built from the column's base counts, so reversing (permuting) the sequence order yields the identical matrix and consensus.")]
    public void Pwm_RowOrder_Invariant()
    {
        var seqs = new[] { "ACGTA", "ACGAA", "TCGTA", "ACGTC" };
        var original = MotifFinder.CreatePwm(seqs);
        var reordered = MotifFinder.CreatePwm(seqs.Reverse().ToList());

        reordered.Length.Should().Be(original.Length);
        reordered.Consensus.Should().Be(original.Consensus, because: "the consensus is derived from order-independent column counts");
        for (int b = 0; b < 4; b++)
            for (int i = 0; i < original.Length; i++)
                reordered.Matrix[b, i].Should().BeApproximately(original.Matrix[b, i], 1e-12,
                    because: $"PWM cell [{b},{i}] depends only on the column's base counts");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MOTIF-REGULATORY-001 — regulatory-element scan (Matching).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 172.
    //
    // API under test (MotifFinder.FindRegulatoryElements):
    //   Scans for a fixed catalogue of known promoter/regulatory motifs and reports each occurrence.
    //
    // Relations (derived from positional motif matching, NOT from output):
    //   • SHIFT (prepend flank shifts positions): a 5' flank with no regulatory motif shifts every
    //          element's position by the flank length.
    //   • SUB  (broader sequence ⇒ ≥ matches): appending more sequence keeps every existing element
    //          (same name and position) and can only add new ones, so the match set grows.
    // ───────────────────────────────────────────────────────────────────────────

    private const string RegulatorySeq = "GGGGTATAAAGGGGAATAAAGGGG"; // TATA box + poly(A) signal

    private static HashSet<(string, int)> RegulatoryHits(string seq) =>
        MotifFinder.FindRegulatoryElements(new DnaSequence(seq)).Select(e => (e.Name, e.Position)).ToHashSet();

    #region MOTIF-REGULATORY-001 SHIFT — a prepended flank shifts the positions

    [Test]
    [Description("SHIFT: prepending a flank that contains no regulatory motif shifts every detected element's position by the flank length.")]
    public void Regulatory_PrependFlank_ShiftsPositions()
    {
        var original = MotifFinder.FindRegulatoryElements(new DnaSequence(RegulatorySeq))
            .Select(e => (e.Name, e.Position)).ToList();
        original.Should().NotBeEmpty();

        foreach (var flank in new[] { "ACAC", "CACACA" })
        {
            var shifted = RegulatoryHits(flank + RegulatorySeq);
            shifted.Should().BeEquivalentTo(original.Select(e => (e.Name, e.Position + flank.Length)),
                because: $"the {flank.Length}-base flank relocates every regulatory element by {flank.Length}");
        }
    }

    #endregion

    #region MOTIF-REGULATORY-001 SUB — appending sequence yields a superset of matches

    [Test]
    [Description("SUB: appending more sequence preserves every existing regulatory element (same name and position) and can only add new ones, so the match set is a superset.")]
    public void Regulatory_AppendSequence_Superset()
    {
        var baseHits = RegulatoryHits(RegulatorySeq);
        var extendedHits = RegulatoryHits(RegulatorySeq + "CCCTTGACACCC"); // adds a -35 box (TTGACA)

        baseHits.IsSubsetOf(extendedHits).Should().BeTrue(because: "appending does not move or remove existing 5' matches");
        extendedHits.Count.Should().BeGreaterThan(baseHits.Count, because: "the appended -35 box adds a new regulatory element");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MOTIF-SHARED-001 — words shared across sequences (Matching).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 173.
    //
    // API under test (MotifFinder.FindSharedMotifs):
    //   Reports k-mers occurring in at least minSequences of the input sequences.
    //
    // Relations (derived from the matching-sequence quorum, NOT from output):
    //   • INV  (input order independent): membership in the shared set depends on which sequences
    //          contain a word, not their order, so the shared word set is order-independent.
    //   • SUB  (fewer inputs ⇒ ⊇ shared set): the "shared by ALL inputs" set is antitone — dropping a
    //          sequence relaxes the all-agree requirement, so the shared set can only grow.
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly DnaSequence SharedS0 = new("ACGTCGTAA");
    private static readonly DnaSequence SharedS1 = new("TTCGTGGTT");
    private static readonly DnaSequence SharedS2 = new("AACGTACGA");

    private static HashSet<string> SharedWords(IEnumerable<DnaSequence> seqs, int k, int minSeq) =>
        MotifFinder.FindSharedMotifs(seqs, k, minSeq).Select(m => m.Sequence).ToHashSet();

    #region MOTIF-SHARED-001 INV — the shared set is independent of input order

    [Test]
    [Description("INV: a word is shared based on which sequences contain it, not their order, so reordering the inputs yields the same shared word set.")]
    public void SharedMotifs_InputOrder_Invariant()
    {
        var ordered = new[] { SharedS0, SharedS1, SharedS2 };
        var shuffled = new[] { SharedS2, SharedS0, SharedS1 };

        SharedWords(shuffled, 3, 2).Should().BeEquivalentTo(SharedWords(ordered, 3, 2),
            because: "the shared word set is a function of the sequence set, not its order");
    }

    #endregion

    #region MOTIF-SHARED-001 SUB — dropping an input grows the shared-by-all set

    [Test]
    [Description("SUB: the set of words shared by ALL inputs is antitone — removing a sequence relaxes the all-agree requirement, so the shared set of the smaller input collection is a superset.")]
    public void SharedMotifs_FewerInputs_Superset()
    {
        var all = SharedWords(new[] { SharedS0, SharedS1, SharedS2 }, 3, minSeq: 3); // shared by all three
        var fewer = SharedWords(new[] { SharedS0, SharedS2 }, 3, minSeq: 2);          // shared by these two

        all.IsSubsetOf(fewer).Should().BeTrue(because: "a word in all three sequences is in any two of them");
        fewer.Count.Should().BeGreaterThan(all.Count, because: "words shared only by S0 and S2 (e.g. ACG) appear once S1 is dropped");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PAT-APPROX-003 — best approximate match (Matching).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 174.
    //
    // API under test (ApproximateMatcher.FindBestMatch):
    //   Returns the equal-length window with the minimum Hamming distance to the pattern.
    //
    // Relations (derived from the minimum-distance definition, NOT from output):
    //   • INV  (exact match ⇒ 0): when the pattern occurs verbatim, the best distance is 0.
    //   • MON  (best ≤ any candidate): the returned distance is the minimum over all windows, so it
    //          is ≤ the Hamming distance of every equal-length window.
    // ───────────────────────────────────────────────────────────────────────────

    #region PAT-APPROX-003 INV — an exact occurrence gives distance 0

    [Test]
    [Description("INV: when the pattern occurs verbatim in the sequence, the best approximate match has Hamming distance 0.")]
    public void BestMatch_ExactOccurrence_DistanceZero()
    {
        var result = ApproximateMatcher.FindBestMatch("ACGTACGTACGT", "CGTA");
        result.Should().NotBeNull();
        result!.Value.Distance.Should().Be(0, because: "CGTA occurs exactly, so the minimum-distance window has distance 0");
    }

    #endregion

    #region PAT-APPROX-003 MON — the best distance is ≤ every window's distance

    [Test]
    [Description("MON: the best match is the minimum-distance window, so its distance is at most the Hamming distance of any equal-length window of the sequence.")]
    public void BestMatch_NoGreaterThanAnyWindow()
    {
        const string seq = "ACGTACGTACGT";
        const string pattern = "TTTT"; // no exact occurrence
        var best = ApproximateMatcher.FindBestMatch(seq, pattern);
        best.Should().NotBeNull();

        for (int i = 0; i <= seq.Length - pattern.Length; i++)
        {
            int windowDistance = ApproximateMatcher.HammingDistance(pattern, seq.Substring(i, pattern.Length));
            best!.Value.Distance.Should().BeLessThanOrEqualTo(windowDistance,
                because: $"the best distance is the minimum over windows, so it cannot exceed window {i}'s distance");
        }
    }

    #endregion
}
