using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for exact/IUPAC/PWM pattern matching.
///
/// Test Units: PAT-EXACT-001, PAT-IUPAC-001, PAT-PWM-001 (Property Extensions), MOTIF-CONS-001, MOTIF-DISCOVER-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class PatternMatchingProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen = 8) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    #region PAT-EXACT-001: R: positions ∈ [0, len-patLen]; M: substring → count≥1; D: deterministic; P: total ≤ len-patLen+1

    /// <summary>
    /// INV-1: Every occurrence position is within valid bounds [0, seqLen - patLen].
    /// Evidence: A match at position p requires p + patLen ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_Positions_WithinBounds_Property()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return positions.All(p => p >= 0 && p + pattern.Length <= seq.Length)
                .Label($"All positions must be in [0, {seq.Length - pattern.Length}]");
        });
    }

    /// <summary>
    /// INV-2: Substring at each found position equals the pattern.
    /// Evidence: Exact match guarantees character-by-character equality.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_FoundSubstrings_EqualPattern_Property()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return positions.All(p => seq.Substring(p, pattern.Length) == pattern)
                .Label("Substring at each position must equal the pattern");
        });
    }

    /// <summary>
    /// INV-3: If the pattern is a known substring, at least one match exists.
    /// Evidence: FindExactMotif must find the pattern at its known position.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_KnownSubstring_AtLeastOneMatch()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            int start = Math.Min(2, seq.Length - 1);
            int len = Math.Min(3, seq.Length - start);
            string pattern = seq.Substring(start, len);
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return (positions.Count >= 1)
                .Label($"Pattern '{pattern}' is a substring at {start}, expected ≥1 match, got {positions.Count}");
        });
    }

    /// <summary>
    /// INV-4: Total matches ≤ seqLen - patLen + 1 (upper bound on overlapping matches).
    /// Evidence: Maximum matches occur when every sliding window matches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_TotalMatches_AtMost_LenMinusPatLenPlus1()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();
            int maxPossible = seq.Length - pattern.Length + 1;

            return (positions.Count <= maxPossible)
                .Label($"Matches={positions.Count} must be ≤ {maxPossible}");
        });
    }

    /// <summary>
    /// INV-5: Exact match is deterministic — same input always yields same positions.
    /// Evidence: FindExactMotif delegates to SuffixTree which is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions1 = MotifFinder.FindExactMotif(dna, pattern).ToList();
            var positions2 = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return positions1.SequenceEqual(positions2)
                .Label("FindExactMotif must be deterministic");
        });
    }

    /// <summary>
    /// A pattern not in the sequence yields no results.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ExactMatch_PatternNotPresent_Empty()
    {
        var dna = new DnaSequence("AAAAAAAAAA");
        var positions = MotifFinder.FindExactMotif(dna, "CCCC").ToList();
        Assert.That(positions, Is.Empty);
    }

    #endregion

    #region PAT-IUPAC-001: M: more degenerate code → ≥ matches; P: N matches all 4 bases; D: deterministic

    /// <summary>
    /// INV-1: All IUPAC match positions are within valid bounds [0, seqLen - motifLen].
    /// Evidence: A match at position p requires p + motifLen ≤ seqLen.
    /// Source: IUPAC-IUB nomenclature (Cornish-Bowden 1985).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IupacMatch_Positions_WithinBounds_Property()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var dna = new DnaSequence(seq);
            string motif = "RYN"; // purine, pyrimidine, any
            var matches = MotifFinder.FindDegenerateMotif(dna, motif).ToList();

            return matches.All(m => m.Position >= 0 && m.Position + motif.Length <= seq.Length)
                .Label($"All positions must be in [0, {seq.Length - motif.Length}]");
        });
    }

    /// <summary>
    /// INV-2: 'N' matches all 4 bases — a motif of all N's yields maximal matches.
    /// Evidence: IUPAC code N = {A, C, G, T}, so every valid window matches.
    /// Source: IUPAC-IUB nomenclature (Cornish-Bowden 1985).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IupacMatch_N_MatchesAllBases_MaximalMatches()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            int motifLen = Math.Min(3, seq.Length);
            string motif = new string('N', motifLen);
            var matches = MotifFinder.FindDegenerateMotif(dna, motif).ToList();
            int expectedCount = seq.Length - motifLen + 1;

            return (matches.Count == expectedCount)
                .Label($"N-motif of len {motifLen} should yield {expectedCount} matches, got {matches.Count}");
        });
    }

    /// <summary>
    /// INV-3: More degenerate motif → ≥ matches than less degenerate.
    /// Evidence: If code C₁ ⊆ C₂, then every C₁ match is also a C₂ match.
    /// Example: Specific base 'A' ⊆ purine 'R' ⊆ any 'N'.
    /// Source: Monotonicity of IUPAC pattern matching.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IupacMatch_MoreDegenerate_MoreOrEqualMatches()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            var dna = new DnaSequence(seq);
            // 'A' ⊂ 'R' (A or G) ⊂ 'N' (any)
            var matchesA = MotifFinder.FindDegenerateMotif(dna, "A").Count();
            var matchesR = MotifFinder.FindDegenerateMotif(dna, "R").Count();
            var matchesN = MotifFinder.FindDegenerateMotif(dna, "N").Count();

            return (matchesN >= matchesR && matchesR >= matchesA)
                .Label($"Expected N({matchesN}) ≥ R({matchesR}) ≥ A({matchesA})");
        });
    }

    /// <summary>
    /// INV-4: Matched sequence at each position conforms to the IUPAC pattern.
    /// Evidence: Each character in matched sequence must be in the IUPAC code's set.
    /// Source: IUPAC-IUB nomenclature (Cornish-Bowden 1985).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IupacMatch_MatchedSequence_ConformsToPattern()
    {
        var iupacSets = new Dictionary<char, HashSet<char>>
        {
            ['R'] = new() { 'A', 'G' },
            ['Y'] = new() { 'C', 'T' },
            ['N'] = new() { 'A', 'C', 'G', 'T' }
        };

        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var dna = new DnaSequence(seq);
            string motif = "RYN";
            var matches = MotifFinder.FindDegenerateMotif(dna, motif).ToList();

            return matches.All(m =>
                m.MatchedSequence.Length == motif.Length &&
                Enumerable.Range(0, motif.Length).All(j =>
                    iupacSets[motif[j]].Contains(m.MatchedSequence[j])))
                .Label("All matched sequences must conform to IUPAC pattern RYN");
        });
    }

    /// <summary>
    /// INV-5: IUPAC match is deterministic — same input always yields same result.
    /// Evidence: FindDegenerateMotif is a pure sliding-window algorithm.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IupacMatch_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var dna = new DnaSequence(seq);
            var matches1 = MotifFinder.FindDegenerateMotif(dna, "RYN").ToList();
            var matches2 = MotifFinder.FindDegenerateMotif(dna, "RYN").ToList();

            return (matches1.Count == matches2.Count &&
                    matches1.Zip(matches2).All(p => p.First.Position == p.Second.Position))
                .Label("FindDegenerateMotif must be deterministic");
        });
    }

    #endregion

    #region PAT-PWM-001: R: scores ∈ ℝ; M: lower threshold → ≥ matches; D: deterministic; P: consensus from PWM valid

    /// <summary>
    /// INV-1: All PWM scan scores are finite real numbers.
    /// Evidence: PWM scores are sums of log-odds ratios; with pseudocounts, no −∞ terms arise.
    /// Source: Stormo (2000) "DNA binding sites: representation and discovery".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PwmScan_Scores_AreFinite_Property()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            var training = new[] { seq[..6], seq[..6], seq[..6] };
            var pwm = MotifFinder.CreatePwm(training);
            var dna = new DnaSequence(seq);
            var matches = MotifFinder.ScanWithPwm(dna, pwm, threshold: -1000).ToList();

            return matches.All(m => double.IsFinite(m.Score))
                .Label("All PWM scores must be finite");
        });
    }

    /// <summary>
    /// INV-2: Lower threshold → ≥ matches (monotonicity).
    /// Evidence: ScanWithPwm returns all positions with score ≥ threshold;
    /// lowering threshold can only include more positions.
    /// Source: Threshold monotonicity of filter operations.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PwmScan_LowerThreshold_MoreOrEqualMatches()
    {
        var training = new[] { "ACGTAC", "ACGTGC", "ACGAAC" };
        var pwm = MotifFinder.CreatePwm(training);
        var dna = new DnaSequence("ACGTACGTGCACGAACACGTACTTTTGGGG");

        var matchesHigh = MotifFinder.ScanWithPwm(dna, pwm, threshold: 2.0).Count();
        var matchesMed = MotifFinder.ScanWithPwm(dna, pwm, threshold: 0.0).Count();
        var matchesLow = MotifFinder.ScanWithPwm(dna, pwm, threshold: -100.0).Count();

        Assert.That(matchesLow, Is.GreaterThanOrEqualTo(matchesMed),
            $"Low threshold matches ({matchesLow}) should be ≥ medium ({matchesMed})");
        Assert.That(matchesMed, Is.GreaterThanOrEqualTo(matchesHigh),
            $"Medium threshold matches ({matchesMed}) should be ≥ high ({matchesHigh})");
    }

    /// <summary>
    /// INV-3: PWM consensus from identical sequences scores highest at the embedded position.
    /// Evidence: When all training sequences are identical, the resulting PWM has maximum
    /// log-odds at positions matching that sequence.
    /// Source: Stormo (2000).
    /// </summary>
    [Test]
    [Category("Property")]
    public void PwmScan_ConsensusScoresHighest_AtKnownPosition()
    {
        string motif = "ACGTAC";
        var pwm = MotifFinder.CreatePwm(new[] { motif, motif, motif });
        string seq = "TTTTTT" + motif + "TTTTTT";
        var dna = new DnaSequence(seq);
        var matches = MotifFinder.ScanWithPwm(dna, pwm, threshold: -1000).ToList();

        Assert.That(matches, Is.Not.Empty, "PWM scan must find matches");
        var best = matches.OrderByDescending(m => m.Score).First();
        Assert.That(best.Position, Is.EqualTo(6),
            $"Consensus at position 6, but best match at {best.Position}");
    }

    /// <summary>
    /// INV-4: PWM scan match positions are within valid bounds.
    /// Evidence: Match at position p requires p + pwm.Length ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PwmScan_Positions_WithinBounds()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            string sub = seq[..Math.Min(6, seq.Length)];
            var pwm = MotifFinder.CreatePwm(new[] { sub, sub });
            var dna = new DnaSequence(seq);
            var matches = MotifFinder.ScanWithPwm(dna, pwm, threshold: -1000).ToList();

            return matches.All(m => m.Position >= 0 && m.Position + pwm.Length <= seq.Length)
                .Label($"All positions must be in [0, {seq.Length - pwm.Length}]");
        });
    }

    /// <summary>
    /// INV-5: PWM scan is deterministic.
    /// Evidence: ScanWithPwm is a pure sliding-window scoring algorithm.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PwmScan_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            string sub = seq[..Math.Min(6, seq.Length)];
            var pwm = MotifFinder.CreatePwm(new[] { sub, sub });
            var dna = new DnaSequence(seq);

            var m1 = MotifFinder.ScanWithPwm(dna, pwm, threshold: -100).ToList();
            var m2 = MotifFinder.ScanWithPwm(dna, pwm, threshold: -100).ToList();

            return (m1.Count == m2.Count &&
                    m1.Zip(m2).All(p => p.First.Position == p.Second.Position &&
                                        Math.Abs(p.First.Score - p.Second.Score) < 1e-10))
                .Label("ScanWithPwm must be deterministic");
        });
    }

    /// <summary>
    /// INV-6: PWM length matches input training sequence length.
    /// Evidence: CreatePwm builds a matrix with one column per position in the aligned sequences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pwm_LengthMatchesInput_Property()
    {
        var motifLenArb = Gen.Choose(4, 8).ToArbitrary();
        return Prop.ForAll(motifLenArb, motifLen =>
        {
            var training = Enumerable.Range(0, 3)
                .Select(_ => new string(Enumerable.Range(0, motifLen)
                    .Select(_ => "ACGT"[Random.Shared.Next(4)]).ToArray()))
                .ToArray();
            var pwm = MotifFinder.CreatePwm(training);
            return (pwm.Length == motifLen)
                .Label($"PWM length {pwm.Length} should equal training length {motifLen}");
        });
    }

    #endregion

    #region MOTIF-CONS-001: P: consensus length = alignment width; P: each column = majority residue; D: deterministic

    // CreateConsensusFromAlignment selects, per column, the most frequent base (ties broken
    // alphabetically A<C<G<T) — the classical per-position consensus (Rosalind CONS).

    /// <summary>Generates 2..5 aligned DNA sequences of equal length 4..10.</summary>
    private static Arbitrary<string[]> AlignedDnaArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int rows = 2 + rng.Next(4);
            int cols = 4 + rng.Next(7);
            const string bases = "ACGT";
            var rowsArr = new string[rows];
            for (int r = 0; r < rows; r++)
            {
                var c = new char[cols];
                for (int j = 0; j < cols; j++) c[j] = bases[rng.Next(4)];
                rowsArr[r] = new string(c);
            }
            return rowsArr;
        }).ToArbitrary();

    private static char ColumnMajority(string[] aln, int col)
    {
        var counts = new int[4]; // A,C,G,T
        foreach (var row in aln) counts["ACGT".IndexOf(row[col])]++;
        int best = 0;
        for (int b = 1; b < 4; b++) if (counts[b] > counts[best]) best = b;
        return "ACGT"[best];
    }

    /// <summary>
    /// INV-1 (P): The consensus length equals the alignment width.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignmentConsensus_Length_EqualsWidth()
    {
        return Prop.ForAll(AlignedDnaArbitrary(), aln =>
            (MotifFinder.CreateConsensusFromAlignment(aln).Length == aln[0].Length)
                .Label("consensus length ≠ alignment width"));
    }

    /// <summary>
    /// INV-2 (P): Each consensus position is the most frequent base of that column (ties broken
    /// alphabetically), verified against an independent per-column majority.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignmentConsensus_EachColumn_IsMajority()
    {
        return Prop.ForAll(AlignedDnaArbitrary(), aln =>
        {
            string consensus = MotifFinder.CreateConsensusFromAlignment(aln);
            bool ok = Enumerable.Range(0, aln[0].Length).All(c => consensus[c] == ColumnMajority(aln, c));
            return ok.Label("a consensus position was not the column majority");
        });
    }

    /// <summary>
    /// INV-3 (D + boundary): consensus is deterministic; empty input is empty; ragged or non-ACGT
    /// alignments are rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AlignmentConsensus_DeterministicAndBoundary()
    {
        var aln = new[] { "ACGT", "AGGT", "ACGA" };
        string first = MotifFinder.CreateConsensusFromAlignment(aln);
        string second = MotifFinder.CreateConsensusFromAlignment(aln);
        Assert.Multiple(() =>
        {
            Assert.That(second, Is.EqualTo(first), "deterministic");
            Assert.That(first, Is.EqualTo("ACGT"), "column majorities");
            Assert.That(MotifFinder.CreateConsensusFromAlignment(Array.Empty<string>()), Is.EqualTo(""));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateConsensusFromAlignment(new[] { "ACGT", "AC" }));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateConsensusFromAlignment(new[] { "ACGT", "ACGN" }));
        });
    }

    #endregion

    #region MOTIF-DISCOVER-001: R: motif length = k; M: lower support → ≥ motifs; D: deterministic

    // DiscoverMotifs reports k-mers occurring at least minCount times, with their positions and an
    // observed/expected enrichment (Compeau & Pevzner i.i.d. background).

    /// <summary>
    /// INV-1 (R): every discovered motif has length k, Count = #positions ≥ minCount, valid positions
    /// pointing at genuine occurrences, and positive enrichment.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DiscoverMotifs_AreWellFormed()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            const int k = 2;
            const int minCount = 2;
            var dna = new DnaSequence(seq);
            var motifs = MotifFinder.DiscoverMotifs(dna, k, minCount).ToList();
            bool ok = motifs.All(m =>
                m.Sequence.Length == k &&
                m.Count >= minCount && m.Count == m.Positions.Count &&
                m.Enrichment > 0 &&
                m.Positions.All(p => p >= 0 && p + k <= seq.Length && seq.Substring(p, k) == m.Sequence));
            return ok.Label("a discovered motif was malformed");
        });
    }

    /// <summary>
    /// INV-2 (M): a lower support threshold reports at least as many motifs — minCount 3 ⊆ minCount 2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DiscoverMotifs_LowerSupport_IsSuperset()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            var dna = new DnaSequence(seq);
            var loose = MotifFinder.DiscoverMotifs(dna, 2, 2).Select(m => m.Sequence).ToHashSet();
            var strict = MotifFinder.DiscoverMotifs(dna, 2, 3).Select(m => m.Sequence).ToHashSet();
            return strict.IsSubsetOf(loose).Label("minCount 3 result not ⊆ minCount 2 result");
        });
    }

    /// <summary>
    /// INV-3 (D): Motif discovery is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DiscoverMotifs_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(12), seq =>
        {
            var dna = new DnaSequence(seq);
            var a = MotifFinder.DiscoverMotifs(dna, 2, 2).Select(m => (m.Sequence, m.Count)).ToHashSet();
            var b = MotifFinder.DiscoverMotifs(dna, 2, 2).Select(m => (m.Sequence, m.Count)).ToHashSet();
            return a.SetEquals(b).Label("DiscoverMotifs must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (boundary): k &lt; 1 is rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DiscoverMotifs_InvalidK_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => MotifFinder.DiscoverMotifs(new DnaSequence("ACGTACGT"), 0).ToList());
    }

    #endregion
}
