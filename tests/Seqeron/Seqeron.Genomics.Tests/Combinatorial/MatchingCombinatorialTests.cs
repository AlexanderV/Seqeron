namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Matching area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing (every t-tuple of parameter values is exercised with a
/// real business assertion; small grids use the exhaustive <c>[Combinatorial]</c>
/// product, a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Matching")]
public class MatchingCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-EXACT-001 — Exact pattern matching (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 8.
    // Dimensions: patLen(3) × seqLen(3) × algorithm(2: SuffixTree/KMP).
    //             Full grid 3×3×2 = 18 cells.
    //
    // Model: exact string matching returns EVERY 0-based start position p with
    // text[p .. p+m) = pattern (m = pattern length). The result is an algorithm
    // INVARIANT: the suffix-tree search (production MotifFinder.FindExactMotif,
    // O(m+occ)) and the Knuth–Morris–Pratt automaton (O(n+m)) must return the
    // identical position multiset, and both must equal the brute-force ground
    // truth — Gusfield (1997) "Algorithms on Strings, Trees and Sequences".
    // Cross-ref docs/checklists/08_DIFFERENTIAL_TESTING.md row 8 ("Suffix tree /
    // KMP" vs "String.IndexOf loop" → "Same positions").
    //
    // The combinatorial point: correctness must hold for every (patLen, seqLen)
    // pairing under BOTH algorithms — short patterns with many overlapping hits,
    // long patterns near the text length, etc. The `algorithm` axis makes the
    // implementation-independence explicit at each grid cell.
    // ═══════════════════════════════════════════════════════════════════════

    public enum ExactAlgorithm { SuffixTree, Kmp }

    // A periodic base unit (≥ 12 nt so the longest tested pattern is a prefix of
    // it). Repeating it embeds each prefix pattern at many positions, exercising
    // multiplicity and overlapping matches.
    private const string BaseUnit = "ACGTTGCAACGGTACC";

    private static string BuildText(int length)
    {
        var sb = new System.Text.StringBuilder(length);
        while (sb.Length < length)
            sb.Append(BaseUnit);
        return sb.ToString(0, length);
    }

    private static List<int> RunAlgorithm(ExactAlgorithm algorithm, string text, string pattern) => algorithm switch
    {
        ExactAlgorithm.SuffixTree => MotifFinder.FindExactMotif(new DnaSequence(text), pattern).ToList(),
        ExactAlgorithm.Kmp => KmpSearch(text, pattern),
        _ => throw new ArgumentOutOfRangeException(nameof(algorithm)),
    };

    /// <summary>Independent KMP reference (Knuth–Morris–Pratt 1977): O(n+m), all start positions ascending.</summary>
    private static List<int> KmpSearch(string text, string pattern)
    {
        var result = new List<int>();
        int m = pattern.Length, n = text.Length;
        if (m == 0 || m > n) return result;

        var lps = new int[m];
        for (int i = 1, len = 0; i < m;)
        {
            if (pattern[i] == pattern[len]) lps[i++] = ++len;
            else if (len != 0) len = lps[len - 1];
            else lps[i++] = 0;
        }

        for (int i = 0, j = 0; i < n;)
        {
            if (text[i] == pattern[j])
            {
                i++; j++;
                if (j == m) { result.Add(i - j); j = lps[j - 1]; }
            }
            else if (j != 0) j = lps[j - 1];
            else i++;
        }
        return result;
    }

    /// <summary>Brute-force ground truth: scan every start offset.</summary>
    private static List<int> BruteForce(string text, string pattern)
    {
        var result = new List<int>();
        if (pattern.Length == 0 || pattern.Length > text.Length) return result;
        for (int i = 0; i + pattern.Length <= text.Length; i++)
            if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern))
                result.Add(i);
        return result;
    }

    /// <summary>
    /// Pairwise grid: every (patLen × seqLen × algorithm) cell. The selected
    /// algorithm must return exactly the brute-force occurrence set for a pattern
    /// known to recur, and exactly nothing for a pattern known to be absent.
    /// </summary>
    [Test, Combinatorial]
    public void PatExact_MatchesBruteForce_AcrossAlgorithms(
        [Values(3, 6, 12)] int patLen,
        [Values(30, 80, 240)] int seqLen,
        [Values(ExactAlgorithm.SuffixTree, ExactAlgorithm.Kmp)] ExactAlgorithm algorithm)
    {
        string text = BuildText(seqLen);
        string present = BaseUnit.Substring(0, patLen);   // a recurring prefix of the unit

        var expected = BruteForce(text, present);
        expected.Should().NotBeEmpty("the periodic text must contain the prefix pattern");

        var actual = RunAlgorithm(algorithm, text, present);
        actual.Should().Equal(expected,
            $"{algorithm} must return the brute-force positions of \"{present}\" (patLen={patLen}, seqLen={seqLen})");

        // Every reported position is a genuine match within bounds.
        foreach (int p in actual)
        {
            p.Should().BeInRange(0, seqLen - patLen);
            text.Substring(p, patLen).Should().Be(present);
        }

        // A pattern guaranteed absent (a homopolymer run longer than any in the
        // periodic unit) yields no matches under the same algorithm.
        string absent = AbsentPattern(patLen);
        var none = RunAlgorithm(algorithm, text, absent);
        none.Should().BeEmpty($"\"{absent}\" does not occur in the periodic text under {algorithm}");
    }

    /// <summary>A pattern guaranteed absent from the periodic BaseUnit text: a homopolymer run longer than any in the unit.</summary>
    private static string AbsentPattern(int patLen) => new string('A', patLen + 2); // unit has no run of ≥3 identical bases

    /// <summary>
    /// Interaction witness: the suffix-tree and KMP algorithms return the IDENTICAL
    /// position list for the same (text, pattern) — the implementation-independence
    /// that the `algorithm` axis asserts cell-by-cell, pinned here on overlapping
    /// matches where a naive scan and a failure-function automaton could diverge.
    /// </summary>
    [Test]
    public void PatExact_SuffixTreeAndKmp_AgreeOnOverlappingMatches()
    {
        string text = "AAAAAAAAAA";   // pattern "AAA" overlaps at positions 0..7
        string pattern = "AAA";

        var st = RunAlgorithm(ExactAlgorithm.SuffixTree, text, pattern);
        var kmp = RunAlgorithm(ExactAlgorithm.Kmp, text, pattern);

        st.Should().Equal(BruteForce(text, pattern));
        kmp.Should().Equal(BruteForce(text, pattern));
        st.Should().Equal(kmp);
        st.Should().Equal(new[] { 0, 1, 2, 3, 4, 5, 6, 7 });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-APPROX-001 — Approximate matching with mismatches (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 9.
    // Dimensions: patLen(3) × maxDist(3) × seqLen(3). Full grid 3×3×3 = 27 cells.
    //
    // Model (Compeau & Pevzner, Bioinformatics Algorithms ch.1, ROSALIND BA1H/BA1N):
    // FindWithMismatches returns every start position i whose equal-length window
    // text[i .. i+m) is within HAMMING distance maxDist of the pattern
    // (substitutions only — no indels). Each result carries that exact distance.
    //
    // The combinatorial point: `maxDist` interacts with `patLen` and `seqLen`.
    // The occurrence set must equal the brute-force Hamming oracle in every cell,
    // and must grow MONOTONICALLY in maxDist (a hit at threshold k is still a hit
    // at k+1) — an interaction across the maxDist axis.
    // ═══════════════════════════════════════════════════════════════════════

    private static int IndependentHamming(ReadOnlySpan<char> a, ReadOnlySpan<char> b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++)
            if (a[i] != b[i]) d++;
        return d;
    }

    private static List<int> BruteForceHamming(string text, string pattern, int maxDist)
    {
        var result = new List<int>();
        for (int i = 0; i + pattern.Length <= text.Length; i++)
            if (IndependentHamming(text.AsSpan(i, pattern.Length), pattern) <= maxDist)
                result.Add(i);
        return result;
    }

    /// <summary>
    /// Pairwise grid: every (patLen × maxDist × seqLen) cell. The Hamming
    /// approximate matcher reproduces the brute-force occurrence set exactly, and
    /// each reported result has a recomputed distance ≤ maxDist that equals the
    /// stored Distance, over an equal-length (substitution-only) window.
    /// </summary>
    [Test, Combinatorial]
    public void PatApprox_Mismatches_MatchBruteForce(
        [Values(4, 8, 15)] int patLen,
        [Values(0, 1, 3)] int maxDist,
        [Values(30, 80, 240)] int seqLen)
    {
        string text = BuildText(seqLen);
        string pattern = BaseUnit.Substring(0, patLen);

        var expected = BruteForceHamming(text, pattern, maxDist);
        var actual = ApproximateMatcher.FindWithMismatches(text, pattern, maxDist).ToList();

        actual.Select(r => r.Position).Should().Equal(expected,
            $"Hamming matcher must equal brute force (patLen={patLen}, maxDist={maxDist}, seqLen={seqLen})");

        foreach (var r in actual)
        {
            r.Position.Should().BeInRange(0, seqLen - patLen);
            r.MatchedSequence.Length.Should().Be(patLen, "Hamming matching uses equal-length windows only");
            r.Distance.Should().BeLessThanOrEqualTo(maxDist);
            r.Distance.Should().Be(IndependentHamming(r.MatchedSequence, pattern));
            r.MismatchType.Should().Be(MismatchType.Substitution);
        }
    }

    /// <summary>
    /// Interaction witness across the maxDist axis: the occurrence set is monotone
    /// non-decreasing in the mismatch threshold (positions found at k remain found
    /// at k+1), and an exact occurrence is present at every threshold ≥ 0.
    /// </summary>
    [Test, Combinatorial]
    public void PatApprox_Mismatches_MonotoneInThreshold(
        [Values(6, 12)] int patLen,
        [Values(60, 180)] int seqLen)
    {
        string text = BuildText(seqLen);
        string pattern = BaseUnit.Substring(0, patLen);

        var d0 = ApproximateMatcher.FindWithMismatches(text, pattern, 0).Select(r => r.Position).ToHashSet();
        var d1 = ApproximateMatcher.FindWithMismatches(text, pattern, 1).Select(r => r.Position).ToHashSet();
        var d2 = ApproximateMatcher.FindWithMismatches(text, pattern, 2).Select(r => r.Position).ToHashSet();

        d0.Should().NotBeEmpty("an exact occurrence of a prefix exists in the periodic text");
        d0.Should().BeSubsetOf(d1, "raising the threshold cannot drop a match");
        d1.Should().BeSubsetOf(d2);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-APPROX-002 — Approximate matching by edit distance (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 10.
    // Dimensions: patLen(3) × maxEdits(3) × seqLen(3) × substCost(2).
    //             Full grid 3×3×3×2 = 54 cells.
    //
    // Model (Levenshtein 1966; Sellers 1980 sequence-vs-text approximate match):
    // FindWithEdits scans every window text[i .. i+len) with len ∈ [m−k, m+k]
    // (k = maxEdits) and reports it when the edit distance to the pattern is ≤ k.
    // The shipped EditDistance is UNIT-cost: substitution = 1, indel = 1. The
    // substitution-cost knob is the classic weighted-Levenshtein generalisation
    // (Wagner–Fischer 1974): with substCost = 2 a substitution costs the same as a
    // delete+insert, so it is never strictly preferred.
    //
    // The combinatorial point: substCost INTERACTS with maxEdits. The reference
    // weighted distance is validated to coincide with the production matcher at
    // substCost = 1 (so it is not unchecked test code); raising substCost to 2 can
    // only SHRINK the match set at a fixed threshold (a substitution gets dearer).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Weighted Levenshtein: substitution cost = <paramref name="subCost"/>, indel cost = 1. At subCost=1 this is identical to the production EditDistance.</summary>
    private static int WeightedEdit(string a, string b, int subCost)
    {
        int m = a.Length, n = b.Length;
        var prev = new int[n + 1];
        var curr = new int[n + 1];
        for (int j = 0; j <= n; j++) prev[j] = j;
        for (int i = 1; i <= m; i++)
        {
            curr[0] = i;
            for (int j = 1; j <= n; j++)
            {
                int cost = a[i - 1] == b[j - 1] ? 0 : subCost;
                curr[j] = Math.Min(Math.Min(prev[j] + 1, curr[j - 1] + 1), prev[j - 1] + cost);
            }
            (prev, curr) = (curr, prev);
        }
        return prev[n];
    }

    /// <summary>Reference window scan replicating FindWithEdits' geometry (len ∈ [m−k, m+k]) under a given substitution cost.</summary>
    private static List<(int Pos, int Len, int Dist)> RefFindEdits(string text, string pattern, int maxEdits, int subCost)
    {
        var result = new List<(int, int, int)>();
        int minLen = Math.Max(1, pattern.Length - maxEdits);
        int maxLen = pattern.Length + maxEdits;
        for (int i = 0; i <= text.Length - minLen; i++)
            for (int len = minLen; len <= maxLen && i + len <= text.Length; len++)
            {
                int d = WeightedEdit(pattern, text.Substring(i, len), subCost);
                if (d <= maxEdits) result.Add((i, len, d));
            }
        return result;
    }

    /// <summary>
    /// Pairwise grid: every (patLen × maxEdits × seqLen × substCost) cell.
    ///   • substCost = 1 (the shipped unit cost): the production edit-distance
    ///     matcher reproduces the reference window/distance set EXACTLY.
    ///   • substCost = 2: a dearer substitution can only remove matches — the
    ///     weighted match windows are a subset of the unit-cost (production) ones,
    ///     and the weighted distance dominates the unit distance window-by-window.
    /// </summary>
    [Test, Combinatorial]
    public void PatApprox_Edits_UnitMatchesProduction_HigherSubstCostShrinks(
        [Values(4, 8, 14)] int patLen,
        [Values(0, 1, 2)] int maxEdits,
        [Values(40, 100, 240)] int seqLen,
        [Values(1, 2)] int substCost)
    {
        string text = BuildText(seqLen);
        string pattern = BaseUnit.Substring(0, patLen);

        var production = ApproximateMatcher.FindWithEdits(text, pattern, maxEdits)
            .Select(r => (r.Position, r.MatchedSequence.Length, r.Distance)).ToList();
        var unitRef = RefFindEdits(text, pattern, maxEdits, subCost: 1);

        // Production is the unit-cost matcher: it must equal the validated reference.
        production.Should().Equal(unitRef,
            $"unit-cost edit matcher must equal the reference (patLen={patLen}, maxEdits={maxEdits}, seqLen={seqLen})");

        foreach (var r in ApproximateMatcher.FindWithEdits(text, pattern, maxEdits))
        {
            r.Distance.Should().BeLessThanOrEqualTo(maxEdits);
            r.Distance.Should().Be(WeightedEdit(pattern, r.MatchedSequence, 1));
        }

        if (substCost == 2)
        {
            var weightedRef = RefFindEdits(text, pattern, maxEdits, subCost: 2);
            var unitWindows = unitRef.Select(t => (t.Pos, t.Len)).ToHashSet();

            weightedRef.Select(t => (t.Pos, t.Len)).Should()
                .BeSubsetOf(unitWindows, "a dearer substitution cannot create new matches");

            foreach (var (pos, len, dist2) in weightedRef)
            {
                int dist1 = WeightedEdit(pattern, text.Substring(pos, len), 1);
                dist2.Should().BeGreaterThanOrEqualTo(dist1, "weighted distance dominates unit distance");
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-IUPAC-001 — Degenerate (IUPAC) motif matching (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 11.
    // Dimensions: patLen(3) × ambiguityLevel(3) × seqLen(3). Full grid 3×3×3 = 27.
    //
    // Model (NC-IUB 1984): a degenerate motif matches window text[i .. i+m) iff at
    // every position j the IUPAC code motif[j] expands to a set CONTAINING the base
    // text[i+j]. FindDegenerateMotif reports every such i. Ambiguity relaxes
    // constraints monotonically: replacing a concrete base by a broader code (up to
    // 'N' = {A,C,G,T}) can only ADD matches.
    //
    // The combinatorial point: ambiguityLevel interacts with patLen and seqLen.
    // The match set must equal an independent IUPAC oracle (Core IupacHelper, a
    // different implementation than MotifFinder's inline switch) in every cell,
    // and grow as the motif becomes more degenerate; the all-'N' motif matches
    // every in-bounds position.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Independent IUPAC oracle via Core IupacHelper (distinct from MotifFinder's switch).</summary>
    private static List<int> BruteForceIupac(string text, string motif)
    {
        var result = new List<int>();
        for (int i = 0; i + motif.Length <= text.Length; i++)
        {
            bool all = true;
            for (int j = 0; j < motif.Length && all; j++)
                all = IupacHelper.MatchesIupac(text[i + j], motif[j]);
            if (all) result.Add(i);
        }
        return result;
    }

    /// <summary>Builds a motif from the concrete prefix with its first <paramref name="nAmbiguous"/> positions set to 'N' (nested relaxation).</summary>
    private static string DegenerateMotif(int patLen, int nAmbiguous)
    {
        char[] m = BaseUnit.Substring(0, patLen).ToCharArray();
        for (int j = 0; j < nAmbiguous; j++) m[j] = 'N';
        return new string(m);
    }

    /// <summary>
    /// Pairwise grid: every (patLen × ambiguityLevel × seqLen) cell. The degenerate
    /// matcher equals the independent IUPAC oracle; matches grow with ambiguity; a
    /// zero-ambiguity motif reduces to exact matching and an all-'N' motif matches
    /// every in-bounds offset.
    /// </summary>
    [Test, Combinatorial]
    public void PatIupac_MatchesOracle_MonotoneInAmbiguity(
        [Values(4, 8, 12)] int patLen,
        [Values(0, 1, 2)] int ambiguityLevel,   // 0 = none, 1 = half 'N', 2 = all 'N'
        [Values(30, 80, 240)] int seqLen)
    {
        string text = BuildText(seqLen);
        var dna = new DnaSequence(text);
        int nAmbiguous = ambiguityLevel switch { 0 => 0, 1 => patLen / 2, _ => patLen };
        string motif = DegenerateMotif(patLen, nAmbiguous);

        var expected = BruteForceIupac(text, motif);
        var actual = MotifFinder.FindDegenerateMotif(dna, motif).ToList();

        actual.Select(m => m.Position).Should().Equal(expected,
            $"degenerate matcher must equal the IUPAC oracle (patLen={patLen}, amb={ambiguityLevel}, seqLen={seqLen})");

        foreach (var m in actual)
        {
            m.Position.Should().BeInRange(0, seqLen - patLen);
            m.Score.Should().Be(1.0);
            for (int j = 0; j < patLen; j++)
                IupacHelper.MatchesIupac(m.MatchedSequence[j], motif[j]).Should().BeTrue();
        }

        if (ambiguityLevel == 0)
            actual.Select(m => m.Position).Should()
                .Equal(MotifFinder.FindExactMotif(dna, motif), "a non-degenerate motif is exact matching");
        if (ambiguityLevel == 2)
            actual.Should().HaveCount(seqLen - patLen + 1, "an all-'N' motif matches every in-bounds offset");
    }

    /// <summary>
    /// Interaction witness across the ambiguity axis: the match set is monotone
    /// non-decreasing as concrete positions are relaxed to 'N' (none ⊆ half ⊆ all).
    /// </summary>
    [Test, Combinatorial]
    public void PatIupac_MonotoneRelaxation(
        [Values(6, 10)] int patLen,
        [Values(60, 200)] int seqLen)
    {
        var dna = new DnaSequence(BuildText(seqLen));
        var none = MotifFinder.FindDegenerateMotif(dna, DegenerateMotif(patLen, 0)).Select(m => m.Position).ToHashSet();
        var half = MotifFinder.FindDegenerateMotif(dna, DegenerateMotif(patLen, patLen / 2)).Select(m => m.Position).ToHashSet();
        var all = MotifFinder.FindDegenerateMotif(dna, DegenerateMotif(patLen, patLen)).Select(m => m.Position).ToHashSet();

        none.Should().BeSubsetOf(half, "relaxing bases to 'N' cannot drop a match");
        half.Should().BeSubsetOf(all);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-PWM-001 — Position-weight-matrix scanning (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 12.
    // Dimensions: motifLen(3) × threshold(3) × seqLen(3) × pseudocount(2).
    //             Full grid 3×3×3×2 = 54 cells.
    //
    // Model (Stormo et al. 1982; Stormo 2000): a PWM stores per-position log-odds
    // log₂(f/0.25) with pseudocount smoothing. The window score is the SUM of the
    // per-position log-odds; ScanWithPwm reports windows whose score ≥ threshold.
    //
    // The combinatorial point: threshold and pseudocount interact with the motif.
    // The reported set must equal the brute-force score filter in every cell; it is
    // monotone non-increasing in `threshold`; and a smaller pseudocount sharpens
    // the matrix so the consensus window scores at least as high.
    // ═══════════════════════════════════════════════════════════════════════

    private static int BaseIndex(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => -1 };

    private static double PwmScore(PositionWeightMatrix pwm, string window)
    {
        double score = 0;
        for (int j = 0; j < pwm.Length; j++)
        {
            int b = BaseIndex(window[j]);
            if (b < 0) return double.NaN;
            score += pwm.Matrix[b, j];
        }
        return score;
    }

    /// <summary>Training set with a strong consensus = the BaseUnit prefix, plus minor variants giving the PWM structure.</summary>
    private static List<string> PwmTraining(int motifLen)
    {
        string consensus = BaseUnit.Substring(0, motifLen);
        var set = new List<string>();
        for (int n = 0; n < 8; n++) set.Add(consensus);          // dominant consensus
        char[] v = consensus.ToCharArray();
        v[0] = v[0] == 'A' ? 'C' : 'A';                          // one minor substitution at position 0
        set.Add(new string(v));
        return set;
    }

    /// <summary>
    /// Pairwise grid: every (motifLen × threshold × seqLen × pseudocount) cell. The
    /// PWM scanner equals the brute-force log-odds score filter; every reported
    /// score is ≥ threshold and equals the recomputed sum.
    /// </summary>
    [Test, Combinatorial]
    public void PatPwm_ScanEqualsScoreFilter(
        [Values(4, 6, 8)] int motifLen,
        [Values(-5.0, 0.0, 3.0)] double threshold,
        [Values(40, 100, 240)] int seqLen,
        [Values(0.25, 1.0)] double pseudocount)
    {
        var pwm = MotifFinder.CreatePwm(PwmTraining(motifLen), pseudocount);
        var dna = new DnaSequence(BuildText(seqLen));
        string text = dna.Sequence;

        var expected = new List<int>();
        for (int i = 0; i + motifLen <= text.Length; i++)
            if (PwmScore(pwm, text.Substring(i, motifLen)) >= threshold)
                expected.Add(i);

        var actual = MotifFinder.ScanWithPwm(dna, pwm, threshold).ToList();
        actual.Select(m => m.Position).Should().Equal(expected,
            $"PWM scan must equal the score filter (motifLen={motifLen}, thr={threshold}, seqLen={seqLen}, ps={pseudocount})");

        foreach (var m in actual)
        {
            m.Score.Should().BeGreaterThanOrEqualTo(threshold);
            m.Score.Should().BeApproximately(PwmScore(pwm, m.MatchedSequence), 1e-9);
        }
    }

    /// <summary>
    /// Interaction witnesses: raising the threshold yields a subset of matches, and
    /// a smaller pseudocount sharpens the matrix so the consensus window scores at
    /// least as high (the pseudocount × score interaction).
    /// </summary>
    [Test]
    public void PatPwm_ThresholdSubset_AndPseudocountSharpens()
    {
        var dna = new DnaSequence(BuildText(160));
        var pwm = MotifFinder.CreatePwm(PwmTraining(6), 0.25);

        var low = MotifFinder.ScanWithPwm(dna, pwm, 0.0).Select(m => m.Position).ToHashSet();
        var high = MotifFinder.ScanWithPwm(dna, pwm, 4.0).Select(m => m.Position).ToHashSet();
        high.Should().BeSubsetOf(low, "raising the threshold cannot add matches");

        string consensus = BaseUnit.Substring(0, 6);
        double sharp = PwmScore(MotifFinder.CreatePwm(PwmTraining(6), 0.25), consensus);
        double smooth = PwmScore(MotifFinder.CreatePwm(PwmTraining(6), 1.0), consensus);
        sharp.Should().BeGreaterThanOrEqualTo(smooth, "a smaller pseudocount sharpens consensus log-odds");
    }
}
