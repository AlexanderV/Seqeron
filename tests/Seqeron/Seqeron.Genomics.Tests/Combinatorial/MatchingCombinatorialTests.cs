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
        st.Should().Equal(0, 1, 2, 3, 4, 5, 6, 7);
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

    // ── Shared helpers for the motif/approximate units below ────────────────────────────────
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static Dictionary<string, List<int>> BruteKmerPositions(string s, int k)
    {
        var d = new Dictionary<string, List<int>>();
        for (int i = 0; i + k <= s.Length; i++)
        {
            string km = s.Substring(i, k);
            if (!d.TryGetValue(km, out var l))
            {
                l = new List<int>();
                d[km] = l;
            }
            l.Add(i);
        }
        return d;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MOTIF-DISCOVER-001 — Overrepresented-motif discovery (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 170.
    // Spec: tests/TestSpecs/MOTIF-DISCOVER-001.md (canonical MotifFinder.DiscoverMotifs). ADVANCED §10.
    // Dimensions: k(3) × support(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive ⊇ pairwise).
    //
    // Model (Compeau & Pevzner): DiscoverMotifs enumerates every length-k window, keeps k-mers whose
    // occurrence count ≥ minCount, and scores each by the O/E enrichment Count / E, where the i.i.d.
    // uniform background expects E = (N − k + 1) / 4ᵏ occurrences of any specific k-mer.
    //
    // The combinatorial point: k, the support floor (minCount) and the sequence length interact —
    // every returned motif's Count and Positions match an independent recount (INV-1), Count ≥ support
    // (INV-3), and Enrichment is exactly Count·4ᵏ/(N−k+1) (INV-2) at every cell.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void MotifDiscover_OverrepresentedKmers_AcrossKSupportLength(
        [Values(4, 6, 8)] int k,
        [Values(2, 3, 4)] int support,
        [Values(48, 96, 160)] int seqLen)
    {
        string text = BuildText(seqLen); // periodic ⇒ many recurrent k-mers
        var dna = new DnaSequence(text);

        var motifs = MotifFinder.DiscoverMotifs(dna, k, support).ToList();
        var brute = BruteKmerPositions(text, k);
        int n = text.Length;
        double expected = (n - k + 1.0) / Math.Pow(4, k);

        motifs.Select(m => m.Sequence).Should().BeEquivalentTo(
            brute.Where(kv => kv.Value.Count >= support).Select(kv => kv.Key),
            "exactly the k-mers occurring ≥ support times are reported (INV-3)");

        foreach (var m in motifs)
        {
            m.Sequence.Length.Should().Be(k);
            m.Count.Should().Be(brute[m.Sequence].Count, "Count is the occurrence count (INV-1)");
            m.Positions.Should().Equal(brute[m.Sequence], "Positions are the 0-based window starts (INV-1)");
            m.Count.Should().BeGreaterThanOrEqualTo(support, "INV-3");
            m.Enrichment.Should().BeApproximately(m.Count / expected, 1e-9, "Enrichment = Count/E (INV-2)");
            m.Enrichment.Should().BeGreaterThan(0, "INV-4");
        }
    }

    /// <summary>
    /// Interaction witness — the support floor is monotone (raising it can only drop motifs) and a
    /// planted tandem k-mer is enriched far above the chance expectation.
    /// </summary>
    [Test]
    public void MotifDiscover_SupportMonotone_AndEnrichmentReflectsOverrepresentation()
    {
        var dna = new DnaSequence(BuildText(160));
        var at2 = MotifFinder.DiscoverMotifs(dna, 6, 2).Select(m => m.Sequence).ToHashSet();
        var at4 = MotifFinder.DiscoverMotifs(dna, 6, 4).Select(m => m.Sequence).ToHashSet();
        at4.Should().BeSubsetOf(at2, "a higher support floor retains no more motifs");

        // A heavily repeated 6-mer is enriched ≫ 1 (observed ≫ expected under the uniform background).
        var planted = new DnaSequence(string.Concat(Enumerable.Repeat("GGGCCC", 20)));
        MotifFinder.DiscoverMotifs(planted, 6, 2).Should().Contain(m => m.Sequence == "GGGCCC" && m.Enrichment > 1.0);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MOTIF-SHARED-001 — Motifs shared across sequences (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 173.
    // Spec: tests/TestSpecs/MOTIF-SHARED-001.md (canonical MotifFinder.FindSharedMotifs). ADVANCED §10.
    // Dimensions: nSeqs(3) × k(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (RSAT oligo-analysis; Das & Dai): a length-k word is "shared" when it occurs (exactly, no
    // substitutions) in at least minSequences of the inputs. FindSharedMotifs counts each word once
    // per sequence and reports words meeting the quorum, with Prevalence = matching-seqs / total-seqs.
    //
    // Axis mapping (documented): the quorum minSequences is fixed at the default 2; the grid varies the
    // number of input sequences and k. The combinatorial point: every reported motif has length k
    // (INV-1) with distinct valid sequence indices (INV-2), the reported set is exactly the words
    // meeting the quorum (INV-3 vs an independent recount), and a word planted in EVERY sequence has
    // Prevalence 1.0 (INV-4).
    // ═══════════════════════════════════════════════════════════════════════

    private const string SharedCore = "TTGGCCAATTGGCCAA"; // 16 nt planted in every sequence

    [Test, Combinatorial]
    public void MotifShared_QuorumWords_AcrossSeqCountAndK(
        [Values(2, 3, 4)] int nSeqs,
        [Values(4, 6, 8)] int k)
    {
        // Each sequence = unique prefix + the shared core + unique suffix.
        var texts = Enumerable.Range(0, nSeqs)
            .Select(i => DiverseDna(12, (uint)(0x1000 + i)) + SharedCore + DiverseDna(12, (uint)(0x9000 + i)))
            .ToList();
        var seqs = texts.Select(t => new DnaSequence(t)).ToList();

        var shared = MotifFinder.FindSharedMotifs(seqs, k, 2).ToList();

        // Independent recount: for each word, the set of sequence indices containing it.
        var present = new Dictionary<string, HashSet<int>>();
        for (int s = 0; s < texts.Count; s++)
            foreach (var km in BruteKmerPositions(texts[s], k).Keys)
            {
                if (!present.TryGetValue(km, out var set))
                {
                    set = new HashSet<int>();
                    present[km] = set;
                }
                set.Add(s);
            }
        var quorumTruth = present.Where(kv => kv.Value.Count >= 2).Select(kv => kv.Key).ToHashSet();

        shared.Select(m => m.Sequence).Should().BeEquivalentTo(quorumTruth, "exactly the quorum words are reported (INV-3)");

        foreach (var m in shared)
        {
            m.Sequence.Length.Should().Be(k, "INV-1");
            m.SequenceIndices.Should().OnlyHaveUniqueItems("each sequence appears at most once (INV-2)");
            m.SequenceIndices.Should().OnlyContain(idx => idx >= 0 && idx < nSeqs, "valid indices (INV-2)");
            m.Prevalence.Should().BeApproximately((double)m.SequenceIndices.Count / nSeqs, 1e-12, "INV-4");
            m.Prevalence.Should().BeInRange(2.0 / nSeqs - 1e-12, 1.0 + 1e-12, "quorum ≥ 2 ⇒ prevalence in (0,1]");
        }

        // The planted core word appears in every sequence ⇒ prevalence 1.0.
        shared.Should().Contain(m => m.Sequence == SharedCore.Substring(0, k) && Math.Abs(m.Prevalence - 1.0) < 1e-12,
            "a word in all sequences has full prevalence");
    }

    /// <summary>
    /// Interaction witness — matching is exact (a single substitution makes a distinct word that is
    /// no longer shared) and the quorum gates reporting.
    /// </summary>
    [Test]
    public void MotifShared_ExactMatchingAndQuorum()
    {
        var seqs = new[] { "AAAGGGGCCCAAA", "TTTGGGGCCCTTT", "TTTGAGGCCCTTT" } // first two share GGGGCCC; third has GAGGCCC
            .Select(t => new DnaSequence(t)).ToList();

        var shared = MotifFinder.FindSharedMotifs(seqs, 7, 2).Select(m => m.Sequence).ToHashSet();
        shared.Should().Contain("GGGGCCC", "two sequences share the exact 7-mer");
        shared.Should().NotContain("GAGGCCC", "a 1-substitution variant occurs in only one sequence (INV-5)");

        // Raising the quorum to all 3 drops GGGGCCC (only 2 sequences carry it).
        MotifFinder.FindSharedMotifs(seqs, 7, 3).Select(m => m.Sequence)
            .Should().NotContain("GGGGCCC", "quorum 3 is not met by a 2-sequence word");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PAT-APPROX-003 — Frequent words & best match with mismatches (Matching)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 174.
    // Spec: tests/TestSpecs/PAT-APPROX-003.md (canonical FindFrequentKmersWithMismatches,
    //       CountApproximateOccurrences, FindBestMatch). ADVANCED §10.
    // Dimensions: patLen(3) × maxDist(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Compeau & Pevzner ch.1, ROSALIND BA1H/BA1I): Count_d(Text,Pattern) is the number of
    // windows within Hamming distance d (INV-2); Count_d ≥ Count_0 = exact count (INV-1).
    // FindFrequentKmersWithMismatches returns every k-mer maximising Count_d; FindBestMatch returns
    // the leftmost equal-length window of minimum Hamming distance (IsExact ⇔ distance 0).
    //
    // The combinatorial point: pattern length, mismatch budget and text length interact — Count_d
    // equals an independent Hamming recount and is monotone in d; the most-frequent-with-mismatches
    // tally cross-checks against Count_d; and the best match's distance is the global window minimum.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void PatApprox003_CountFrequentAndBestMatch_AcrossPatLenDistLength(
        [Values(4, 6, 8)] int patLen,
        [Values(0, 1, 2)] int maxDist,
        [Values(32, 48, 64)] int seqLen)
    {
        string text = BuildText(seqLen);
        string pattern = BaseUnit.Substring(0, patLen); // occurs exactly in the periodic text

        // Count_d == independent Hamming recount (INV-2), and ≥ exact count (INV-1).
        int countD = ApproximateMatcher.CountApproximateOccurrences(text, pattern, maxDist);
        countD.Should().Be(BruteForceHamming(text, pattern, maxDist).Count, "Count_d = #windows within d (INV-2)");
        int count0 = ApproximateMatcher.CountApproximateOccurrences(text, pattern, 0);
        countD.Should().BeGreaterThanOrEqualTo(count0, "Count_d ≥ Count_0 (INV-1)");

        // FindFrequentKmersWithMismatches: all returned k-mers tie at the max Count_d, and that count
        // equals an independent CountApproximateOccurrences for the same k-mer (cross-check, INV-3).
        var frequent = ApproximateMatcher.FindFrequentKmersWithMismatches(text, patLen, maxDist).ToList();
        frequent.Should().NotBeEmpty();
        int maxCount = frequent[0].Count;
        frequent.Should().OnlyContain(f => f.Count == maxCount, "all returned k-mers share the maximum Count_d (INV-3)");
        foreach (var (kmer, count) in frequent)
            ApproximateMatcher.CountApproximateOccurrences(text, kmer, maxDist).Should().Be(count,
                "the neighborhood tally equals Count_d for that k-mer");
        maxCount.Should().BeGreaterThanOrEqualTo(countD, "the maximum Count_d is ≥ any specific pattern's Count_d");

        // FindBestMatch: distance is the global minimum Hamming over equal-length windows (INV-4),
        // and the planted pattern is found exactly at its leftmost occurrence (INV-5).
        var best = ApproximateMatcher.FindBestMatch(text, pattern);
        best.Should().NotBeNull();
        int minDist = Enumerable.Range(0, text.Length - patLen + 1)
            .Min(i => IndependentHamming(pattern, text.AsSpan(i, patLen)));
        best!.Value.Distance.Should().Be(minDist, "best distance is the global window minimum (INV-4)");
        best.Value.IsExact.Should().Be(minDist == 0, "IsExact ⇔ distance 0 (INV-4)");
        best.Value.Position.Should().Be(text.IndexOf(pattern, StringComparison.Ordinal), "leftmost minimum (INV-5)");
    }

    /// <summary>
    /// Interaction witness — Count_d is monotone non-decreasing in the mismatch budget, and a pattern
    /// absent exactly can still have approximate occurrences once d is large enough.
    /// </summary>
    [Test]
    public void PatApprox003_CountMonotoneInDistance()
    {
        string text = BuildText(96);
        string pattern = BaseUnit.Substring(0, 6);

        int c0 = ApproximateMatcher.CountApproximateOccurrences(text, pattern, 0);
        int c1 = ApproximateMatcher.CountApproximateOccurrences(text, pattern, 1);
        int c2 = ApproximateMatcher.CountApproximateOccurrences(text, pattern, 2);
        c1.Should().BeGreaterThanOrEqualTo(c0);
        c2.Should().BeGreaterThanOrEqualTo(c1, "a larger mismatch budget admits ≥ as many windows");

        // A pattern one substitution off every occurrence: 0 exact, but > 0 at d ≥ 1.
        string near = string.Concat("T", BaseUnit.AsSpan(1, 5)); // differs from BaseUnit[0..6] only at position 0
        ApproximateMatcher.CountApproximateOccurrences(text, near, 0)
            .Should().BeLessThan(ApproximateMatcher.CountApproximateOccurrences(text, near, 1));
    }
}
