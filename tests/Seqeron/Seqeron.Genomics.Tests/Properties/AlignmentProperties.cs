using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for global, local, semi-global, and multiple alignment.
/// Verifies structural invariants of alignment results.
///
/// Test Units: ALIGN-GLOBAL-001, ALIGN-LOCAL-001, ALIGN-SEMI-001, ALIGN-MULTI-001, ALIGN-STATS-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Alignment")]
public class AlignmentProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen = 8) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    #region ALIGN-GLOBAL-001: S: score(a,b)=score(b,a); R: aligned len ≥ max(len1,len2); P: identity → max score; P: aligned1.len = aligned2.len

    /// <summary>
    /// INV-1: Global alignment score is symmetric: score(A,B) == score(B,A).
    /// Evidence: Needleman-Wunsch recurrence is symmetric — S(a,b) = S(b,a)
    /// and gap penalties apply equally to both sequences.
    /// Source: Needleman–Wunsch algorithm (1970).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_ScoreSymmetry()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var ab = SequenceAligner.GlobalAlign(s1, s2);
            var ba = SequenceAligner.GlobalAlign(s2, s1);
            return (ab.Score == ba.Score)
                .Label($"score(A,B)={ab.Score} ≠ score(B,A)={ba.Score}");
        });
    }

    /// <summary>
    /// INV-2: Both aligned sequences have equal length.
    /// Evidence: Global alignment inserts gaps to make sequences equal length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_AlignedSequences_HaveEqualLength()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.GlobalAlign(s1, s2);
            return (result.AlignedSequence1.Length == result.AlignedSequence2.Length)
                .Label($"len1={result.AlignedSequence1.Length} ≠ len2={result.AlignedSequence2.Length}");
        });
    }

    /// <summary>
    /// INV-3: Aligned sequence length ≥ max(len1, len2).
    /// Evidence: A global alignment must cover both sequences entirely,
    /// so gaps can only increase the length beyond the longer input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_AlignedLength_AtLeastMaxInput()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.GlobalAlign(s1, s2);
            int maxLen = Math.Max(s1.Length, s2.Length);
            return (result.AlignedSequence1.Length >= maxLen)
                .Label($"Aligned len={result.AlignedSequence1.Length} < max(inputs)={maxLen}");
        });
    }

    /// <summary>
    /// INV-4: Identical sequences produce maximum score (identity → max score).
    /// Evidence: No gaps or mismatches needed for identical sequences.
    /// Score = matchScore × length with default SimpleDna matrix (1 per match).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_IdenticalSequences_MaxScore()
    {
        return Prop.ForAll(DnaArbitrary(6), seq =>
        {
            var result = SequenceAligner.GlobalAlign(seq, seq);
            int expectedScore = seq.Length; // SimpleDna: match=1
            return (result.Score == expectedScore)
                .Label($"Score={result.Score}, expected={expectedScore} for identical seqs of len {seq.Length}");
        });
    }

    /// <summary>
    /// INV-5: Global alignment is deterministic.
    /// Evidence: GlobalAlign is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GlobalAlign_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var r1 = SequenceAligner.GlobalAlign(s1, s2);
            var r2 = SequenceAligner.GlobalAlign(s1, s2);
            return (r1.Score == r2.Score &&
                    r1.AlignedSequence1 == r2.AlignedSequence1 &&
                    r1.AlignedSequence2 == r2.AlignedSequence2)
                .Label("GlobalAlign must be deterministic");
        });
    }

    /// <summary>
    /// CancellationToken overload produces the same result as the standard overload.
    /// Verifies the separate code path (non-pooled 2D array) is functionally equivalent.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GlobalAlign_CancellationOverload_SameResultAsStandard()
    {
        using var cts = new CancellationTokenSource();
        var standard = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT");
        var withToken = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT", null, cts.Token);

        Assert.Multiple(() =>
        {
            Assert.That(withToken.Score, Is.EqualTo(standard.Score),
                "CancellationToken overload must produce same score");
            Assert.That(withToken.AlignedSequence1, Is.EqualTo(standard.AlignedSequence1));
            Assert.That(withToken.AlignedSequence2, Is.EqualTo(standard.AlignedSequence2));
        });
    }

    #endregion

    #region ALIGN-LOCAL-001: R: score ≥ 0; P: aligned1.len = aligned2.len; M: identical substring → score ≥ matchScore×len

    /// <summary>
    /// INV-1: Local alignment score is non-negative.
    /// Evidence: Smith-Waterman matrix values are clamped to ≥ 0.
    /// Source: Smith &amp; Waterman (1981).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LocalAlign_Score_NonNegative()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.LocalAlign(s1, s2);
            return (result.Score >= 0)
                .Label($"Local alignment score={result.Score} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2: Both aligned sequences have equal length.
    /// Evidence: Local alignment pads with gaps to equalize aligned region length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LocalAlign_AlignedSequences_HaveEqualLength()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.LocalAlign(s1, s2);
            return (result.AlignedSequence1.Length == result.AlignedSequence2.Length)
                .Label($"len1={result.AlignedSequence1.Length} ≠ len2={result.AlignedSequence2.Length}");
        });
    }

    /// <summary>
    /// INV-3: When sequences share an identical substring, score ≥ matchScore × substring length.
    /// Evidence: Smith-Waterman must find at least the embedded common substring.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LocalAlign_IdenticalSubstring_ScoreAtLeastMatchTimesLength()
    {
        string common = "ACGTACGT";
        string s1 = "TTT" + common + "TTT";
        string s2 = "GGG" + common + "GGG";
        var result = SequenceAligner.LocalAlign(s1, s2);

        // SimpleDna: match = 1
        int minScore = common.Length;
        Assert.That(result.Score, Is.GreaterThanOrEqualTo(minScore),
            $"Score={result.Score} must be ≥ {minScore} for shared '{common}'");
    }

    /// <summary>
    /// INV-4: Local alignment is deterministic.
    /// Evidence: LocalAlign is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LocalAlign_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var r1 = SequenceAligner.LocalAlign(s1, s2);
            var r2 = SequenceAligner.LocalAlign(s1, s2);
            return (r1.Score == r2.Score)
                .Label("LocalAlign must be deterministic");
        });
    }

    #endregion

    #region ALIGN-SEMI-001: P: aligned1.len = aligned2.len; R: score ≥ 0; D: deterministic

    /// <summary>
    /// INV-1: Semi-global aligned sequences have equal length.
    /// Evidence: Semi-global alignment pads with gaps to equalize aligned region length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SemiGlobalAlign_AlignedSequences_HaveEqualLength()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(s1), new DnaSequence(s2));
            return (result.AlignedSequence1.Length == result.AlignedSequence2.Length)
                .Label($"len1={result.AlignedSequence1.Length} ≠ len2={result.AlignedSequence2.Length}");
        });
    }

    /// <summary>
    /// INV-2: Semi-global alignment score is finite.
    /// Evidence: SemiGlobalAlign produces a valid numeric score for any input pair.
    /// Note: score can be negative when mismatches dominate with default scoring.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SemiGlobalAlign_Score_IsFinite()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.SemiGlobalAlign(new DnaSequence(s1), new DnaSequence(s2));
            return double.IsFinite(result.Score)
                .Label($"Semi-global score={result.Score} must be finite");
        });
    }

    /// <summary>
    /// INV-3: Semi-global alignment is deterministic.
    /// Evidence: SemiGlobalAlign is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SemiGlobalAlign_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var d1 = new DnaSequence(s1);
            var d2 = new DnaSequence(s2);
            var r1 = SequenceAligner.SemiGlobalAlign(d1, d2);
            var r2 = SequenceAligner.SemiGlobalAlign(d1, d2);
            return (r1.Score == r2.Score &&
                    r1.AlignedSequence1 == r2.AlignedSequence1 &&
                    r1.AlignedSequence2 == r2.AlignedSequence2)
                .Label("SemiGlobalAlign must be deterministic");
        });
    }

    #endregion

    #region ALIGN-MULTI-001: P: all aligned sequences same length; R: score ≥ 0; D: deterministic

    /// <summary>An arbitrary set of 2–6 unrelated DNA sequences (≥ 6 bp each).</summary>
    private static Arbitrary<DnaSequence[]> DnaSetArbitrary(int minLen = 6) =>
        (from n in Gen.Choose(2, 6)
         from arr in (from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
                      select new DnaSequence(new string(chars))).ArrayOf(n)
         select arr).ToArbitrary();

    /// <summary>
    /// A set of homologous sequences: a common ancestor (≥ 40 bp) plus ≤ 3 substitutions each and no
    /// indels, so identity stays high and the sum-of-pairs score (SimpleDna: match +1, mismatch −1)
    /// remains non-negative. This is the biologically meaningful domain for the "score ≥ 0" claim —
    /// SP can legitimately go negative for unrelated sequences.
    /// </summary>
    private static Arbitrary<DnaSequence[]> HomologousDnaSetArbitrary() =>
        (from baseChars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 40)
         from k in Gen.Choose(2, 5)
         let baseStr = new string(baseChars)
         from edits in MutationPositionsGen(baseStr.Length).ArrayOf(k)
         select edits.Select(positions => new DnaSequence(ApplySubstitutions(baseStr, positions))).ToArray())
        .ToArbitrary();

    private static Gen<int[]> MutationPositionsGen(int len) =>
        from m in Gen.Choose(0, 3)
        from positions in Gen.Choose(0, len - 1).ArrayOf(m)
        select positions;

    private static string ApplySubstitutions(string s, int[] positions)
    {
        var arr = s.ToCharArray();
        foreach (int p in positions)
            arr[p] = arr[p] switch { 'A' => 'C', 'C' => 'G', 'G' => 'T', _ => 'A' };
        return new string(arr);
    }

    /// <summary>
    /// INV-1 (P): a multiple alignment is a rectangular block — every aligned row has the same length
    /// (gaps are inserted to equalize lengths). Holds for any input set.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MultipleAlign_AllSequences_HaveEqualLength()
    {
        return Prop.ForAll(DnaSetArbitrary(), seqs =>
        {
            var aligned = SequenceAligner.MultipleAlign(seqs).AlignedSequences;
            int len = aligned[0].Length;
            return aligned.All(s => s.Length == len)
                .Label($"aligned rows differ in length: [{string.Join(",", aligned.Select(s => s.Length))}]");
        });
    }

    /// <summary>
    /// INV-2 (R): the sum-of-pairs score is non-negative for homologous sequences (common ancestor
    /// with few substitutions). Identity-dominated columns keep SP ≥ 0 (SimpleDna match +1).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MultipleAlign_TotalScore_NonNegative_ForHomologs()
    {
        return Prop.ForAll(HomologousDnaSetArbitrary(), seqs =>
        {
            int score = SequenceAligner.MultipleAlign(seqs).TotalScore;
            return (score >= 0).Label($"SP score {score} must be ≥ 0 for homologous sequences");
        });
    }

    /// <summary>
    /// INV-3 (P): the consensus is derived column-by-column, so its length equals the aligned width.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MultipleAlign_ConsensusLength_EqualsAlignedLength()
    {
        return Prop.ForAll(DnaSetArbitrary(), seqs =>
        {
            var result = SequenceAligner.MultipleAlign(seqs);
            return (result.Consensus.Length == result.AlignedSequences[0].Length)
                .Label($"consensus {result.Consensus.Length} != aligned width {result.AlignedSequences[0].Length}");
        });
    }

    /// <summary>
    /// INV-4 (D): MultipleAlign is a pure function — identical inputs yield identical aligned rows,
    /// consensus and total score.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MultipleAlign_IsDeterministic()
    {
        return Prop.ForAll(DnaSetArbitrary(), seqs =>
        {
            var r1 = SequenceAligner.MultipleAlign(seqs);
            var r2 = SequenceAligner.MultipleAlign(seqs);
            bool same = r1.TotalScore == r2.TotalScore
                        && r1.Consensus == r2.Consensus
                        && r1.AlignedSequences.SequenceEqual(r2.AlignedSequences);
            return same.Label("MultipleAlign must be deterministic");
        });
    }

    #endregion

    #region Alignment Statistics

    /// <summary>
    /// Alignment statistics: identity + mismatch + gap percentages should account for full alignment.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Statistics_Matches_Plus_Mismatches_Plus_Gaps_EqualsLength()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGTT");
        var stats = SequenceAligner.CalculateStatistics(result);
        Assert.That(stats.Matches + stats.Mismatches + stats.Gaps, Is.EqualTo(stats.AlignmentLength));
    }

    /// <summary>
    /// Alignment statistics: percentage fields satisfy structural invariants.
    ///   - 0 ≤ Identity ≤ 100
    ///   - 0 ≤ Similarity ≤ 100
    ///   - 0 ≤ GapPercent ≤ 100
    ///   - Identity ≤ Similarity  (matches ≤ matches + mismatches)
    ///   - Similarity + GapPercent ≈ 100  (non-gap + gap = total)
    /// </summary>
    [Test]
    [Category("Property")]
    public void Statistics_PercentageFields_SatisfyInvariants()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGTT");
        var stats = SequenceAligner.CalculateStatistics(result);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Identity, Is.InRange(0.0, 100.0), "Identity in [0,100]");
            Assert.That(stats.Similarity, Is.InRange(0.0, 100.0), "Similarity in [0,100]");
            Assert.That(stats.GapPercent, Is.InRange(0.0, 100.0), "GapPercent in [0,100]");
            Assert.That(stats.Identity, Is.LessThanOrEqualTo(stats.Similarity),
                "Identity ≤ Similarity (matches ≤ matches + mismatches)");
            Assert.That(stats.Similarity + stats.GapPercent, Is.EqualTo(100.0).Within(0.001),
                "Similarity + GapPercent = 100 (all positions accounted for)");
        });
    }

    #endregion

    #region ALIGN-STATS-001: P: matches+mismatches+gaps = length; R: identity ∈ [0,100]; D: deterministic

    // CalculateStatistics tallies an alignment's columns. NOTE: Identity/Similarity/GapPercent are
    // EMBOSS-style PERCENTAGES in [0,100] (the checklist's [0,1] is the fraction = percentage/100).

    /// <summary>
    /// INV-1 (P + R): for any aligned pair, matches+mismatches+gaps equals the alignment length,
    /// identity equals matches/length·100 and lies in [0,100].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignmentStats_PartitionAndIdentity()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.GlobalAlign(s1, s2);
            var stats = SequenceAligner.CalculateStatistics(result);
            bool partition = stats.Matches + stats.Mismatches + stats.Gaps == stats.AlignmentLength;
            bool identityOk = stats.Identity is >= 0.0 and <= 100.0
                              && Math.Abs(stats.Identity - (double)stats.Matches / stats.AlignmentLength * 100) < 1e-6;
            return (partition && identityOk).Label($"partition={partition}, identity={stats.Identity}");
        });
    }

    /// <summary>
    /// INV-2 (positive control): identical sequences align to 100% identity with no gaps/mismatches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignmentStats_IdenticalSequences_FullIdentity()
    {
        return Prop.ForAll(DnaArbitrary(6), seq =>
        {
            var stats = SequenceAligner.CalculateStatistics(SequenceAligner.GlobalAlign(seq, seq));
            return (stats.Matches == stats.AlignmentLength && stats.Mismatches == 0 && stats.Gaps == 0
                    && Math.Abs(stats.Identity - 100.0) < 1e-9)
                .Label($"identical alignment not 100% identity: {stats.Identity}");
        });
    }

    /// <summary>
    /// INV-3 (D): Alignment statistics are deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlignmentStats_AreDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(6), DnaArbitrary(6), (s1, s2) =>
        {
            var result = SequenceAligner.GlobalAlign(s1, s2);
            var a = SequenceAligner.CalculateStatistics(result);
            var b = SequenceAligner.CalculateStatistics(result);
            return (a.Matches == b.Matches && a.Identity == b.Identity && a.AlignmentLength == b.AlignmentLength)
                .Label("CalculateStatistics must be deterministic");
        });
    }

    #endregion
}
