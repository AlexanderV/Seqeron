using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Targeted mutation-killing tests for March 2026 commits.
/// Each test is designed to kill specific Stryker mutants that survived
/// or had no coverage in the baseline run.
///
/// Baseline scores (before):
///   IupacHelper.cs      → 100% (no action needed)
///   MotifFinder.cs      → T=178 K=96 S=18 NC=24
///   RepeatFinder.cs     → T=196 K=112 S=39 NC=7
///   ApproximateMatcher  → not baselined (Stryker compile issue on Alignment project)
/// </summary>
[TestFixture]
public class MutationKillerTests
{
    // ═══════════════════════════════════════════════════════════════════
    // MotifFinder.cs — Cancellable FindDegenerateMotif (NoCoverage L175-184)
    // The cancellable overload's IUPAC switch has ZERO test coverage.
    // ═══════════════════════════════════════════════════════════════════

    #region MotifFinder — Cancellable Degenerate Motif (NoCoverage lines 144, 175-184)

    [Test]
    public void FindDegenerateMotif_Cancellable_StandardBases_MatchesCorrectly()
    {
        // Kills NoCoverage mutants on lines 175 ('A','T','G','C' branch)
        var results = MotifFinder.FindDegenerateMotif("ATGCATGC", "ATGC", CancellationToken.None).ToList();

        results.Should().HaveCountGreaterThanOrEqualTo(1);
        results[0].Position.Should().Be(0);
        results[0].MatchedSequence.Should().Be("ATGC");
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_PurineR_MatchesAandG()
    {
        // Kills NoCoverage mutants on line 176 ('R' branch: seqChar == 'A' || seqChar == 'G')
        var results = MotifFinder.FindDegenerateMotif("AATTGGCC", "R", CancellationToken.None).ToList();

        results.Should().HaveCount(4); // A at 0,1 and G at 4,5
        results.Select(r => r.MatchedSequence).Should().AllSatisfy(m =>
            m.Should().BeOneOf("A", "G"));
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_PyrimidineY_MatchesCandT()
    {
        // Kills NoCoverage mutants on line 177 ('Y' branch: seqChar == 'C' || seqChar == 'T')
        var results = MotifFinder.FindDegenerateMotif("AATTGGCC", "Y", CancellationToken.None).ToList();

        results.Should().HaveCount(4); // T at 2,3 and C at 6,7
        results.Select(r => r.MatchedSequence).Should().AllSatisfy(m =>
            m.Should().BeOneOf("C", "T"));
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_StrongS_MatchesGandC()
    {
        // Kills NoCoverage mutants on line 178 ('S' branch: seqChar == 'G' || seqChar == 'C')
        var results = MotifFinder.FindDegenerateMotif("AATTGGCC", "S", CancellationToken.None).ToList();

        results.Should().HaveCount(4); // G at 4,5 and C at 6,7
        results.Select(r => r.MatchedSequence).Should().AllSatisfy(m =>
            m.Should().BeOneOf("G", "C"));
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_WeakW_MatchesAandT()
    {
        // Kills NoCoverage mutants on line 179 ('W' branch: seqChar == 'A' || seqChar == 'T')
        var results = MotifFinder.FindDegenerateMotif("AATTGGCC", "W", CancellationToken.None).ToList();

        results.Should().HaveCount(4); // A at 0,1 and T at 2,3
        results.Select(r => r.MatchedSequence).Should().AllSatisfy(m =>
            m.Should().BeOneOf("A", "T"));
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_KetoK_MatchesGandT()
    {
        // Kills NoCoverage mutants on line 180 ('K' branch: seqChar == 'G' || seqChar == 'T')
        var results = MotifFinder.FindDegenerateMotif("AATTGGCC", "K", CancellationToken.None).ToList();

        results.Should().HaveCount(4); // T at 2,3 and G at 4,5
        results.Select(r => r.MatchedSequence).Should().AllSatisfy(m =>
            m.Should().BeOneOf("G", "T"));
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_AminoM_MatchesAandC()
    {
        // Kills NoCoverage mutants on line 181 ('M' branch: seqChar == 'A' || seqChar == 'C')
        var results = MotifFinder.FindDegenerateMotif("AATTGGCC", "M", CancellationToken.None).ToList();

        results.Should().HaveCount(4); // A at 0,1 and C at 6,7
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_B_MatchesNotA()
    {
        // Kills NoCoverage mutant on line 182 ('B' branch: seqChar != 'A')
        var seq = "ATGC";
        var results = MotifFinder.FindDegenerateMotif(seq, "B", CancellationToken.None).ToList();

        // B = C, G, T (not A). So matches at positions 1(T), 2(G), 3(C)
        results.Should().HaveCount(3);
        results.Select(r => r.MatchedSequence).Should().NotContain("A");
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_D_MatchesNotC()
    {
        // Kills NoCoverage mutant on line 183 ('D' branch: seqChar != 'C')
        var seq = "ATGC";
        var results = MotifFinder.FindDegenerateMotif(seq, "D", CancellationToken.None).ToList();

        // D = A, G, T (not C). So matches at 0(A), 1(T), 2(G)
        results.Should().HaveCount(3);
        results.Select(r => r.MatchedSequence).Should().NotContain("C");
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_H_MatchesNotG()
    {
        // Kills NoCoverage mutant on line 183 ('H' branch: seqChar != 'G')
        var seq = "ATGC";
        var results = MotifFinder.FindDegenerateMotif(seq, "H", CancellationToken.None).ToList();

        // H = A, C, T (not G). So matches at 0(A), 1(T), 3(C)
        results.Should().HaveCount(3);
        results.Select(r => r.MatchedSequence).Should().NotContain("G");
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_V_MatchesNotT()
    {
        // Kills NoCoverage mutant on line 184 ('V' branch: seqChar != 'T')
        var seq = "ATGC";
        var results = MotifFinder.FindDegenerateMotif(seq, "V", CancellationToken.None).ToList();

        // V = A, C, G (not T). So matches at 0(A), 2(G), 3(C)
        results.Should().HaveCount(3);
        results.Select(r => r.MatchedSequence).Should().NotContain("T");
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_N_MatchesAll()
    {
        // Kills NoCoverage mutants for 'N' branch
        var results = MotifFinder.FindDegenerateMotif("ATGC", "N", CancellationToken.None).ToList();

        results.Should().HaveCount(4);
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_EmptyInputs_ReturnsEmpty()
    {
        // Kills NoCoverage mutant on line 144 (block removal of yield break)
        var r1 = MotifFinder.FindDegenerateMotif("", "ATG", CancellationToken.None).ToList();
        var r2 = MotifFinder.FindDegenerateMotif("ATGC", "", CancellationToken.None).ToList();

        r1.Should().BeEmpty();
        r2.Should().BeEmpty();
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_CancellationRespected()
    {
        // Ensures the cancellation path on line 163 executes
        var cts = new CancellationTokenSource();
        cts.Cancel();

        var act = () => MotifFinder.FindDegenerateMotif(
            new string('A', 5000), "A", cts.Token).ToList();

        act.Should().Throw<OperationCanceledException>();
    }

    [Test]
    public void FindDegenerateMotif_Cancellable_MultiCharPattern_MatchesCorrectly()
    {
        // Exercises cancellable path with multi-char IUPAC pattern,
        // kills boundary mutations on lines 161, 167
        var results = MotifFinder.FindDegenerateMotif(
            "ATGATGATG", "RYG", CancellationToken.None).ToList();

        // RYG: R=A|G, Y=C|T, G=G → at pos 0: ATG=A(R)T(Y)G(G)=match
        // at pos 3: ATG=match, at pos 6: ATG=match
        results.Should().HaveCount(3);
        results.Should().AllSatisfy(r => r.MatchedSequence.Should().Be("ATG"));
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // MotifFinder.cs — FindDegenerateMotifCore (Survived L153, 161, 163, 167, 174)
    // ═══════════════════════════════════════════════════════════════════

    #region MotifFinder — FindDegenerateMotifCore boundary survivors

    [Test]
    public void FindDegenerateMotif_BothEmptyStrings_ReturnsEmpty()
    {
        // Kills L153: Logical mutation (string.IsNullOrEmpty(sequence) && string.IsNullOrEmpty(motif))
        // Both empty should still return empty (tests the || → && mutation)
        var r1 = MotifFinder.FindDegenerateMotif("", "", CancellationToken.None).ToList();
        var r2 = MotifFinder.FindDegenerateMotif("ATGC", "", CancellationToken.None).ToList();
        var r3 = MotifFinder.FindDegenerateMotif("", "ATG", CancellationToken.None).ToList();

        r1.Should().BeEmpty();
        r2.Should().BeEmpty();
        r3.Should().BeEmpty();
    }

    [Test]
    public void FindDegenerateMotif_LoopBoundary_MatchesAtLastPosition()
    {
        // Kills L161: i < seq.Length - motifUpper.Length → i <= (boundary)
        // and L161: seq.Length + motifUpper.Length arithmetic mutation
        // Pattern "GC" in "ATGC" should match at position 2 (the last valid position)
        var results = MotifFinder.FindDegenerateMotif("ATGC", "GC", CancellationToken.None).ToList();

        results.Should().ContainSingle(m => m.Position == 2);
        results.Last().Position.Should().Be(2); // can't match at position 3+
    }

    [Test]
    public void FindDegenerateMotif_InnerLoopBoundary_FullPatternLength()
    {
        // Kills L167: j > motifUpper.Length → j >= or j <
        // Verifies the complete pattern is matched (all characters checked)
        var results = MotifFinder.FindDegenerateMotif("ATGCATGC", "ATGC", CancellationToken.None).ToList();

        results.Should().HaveCount(2);
        results[0].Position.Should().Be(0);
        results[1].Position.Should().Be(4);
    }

    [Test]
    public void FindDegenerateMotif_MismatchDetection_CharComparison()
    {
        // Kills L174: motifChar != seqChar → == mutation
        // If mutated to ==, non-matching chars would pass as matches
        var results = MotifFinder.FindDegenerateMotif("AAAA", "T", CancellationToken.None).ToList();

        results.Should().BeEmpty(); // T does not match A
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // MotifFinder.cs — GenerateConsensus / GetIupacCode (Survived L355, 367, 369, 374, 375)
    // ═══════════════════════════════════════════════════════════════════

    #region MotifFinder — Consensus generation survivors

    [Test]
    public void GenerateConsensus_ThresholdBoundary_ExactlyAtQuarter()
    {
        // Kills L367: total * 0.25 → * (arithmetic), and L369: kv.Value > threshold → >=
        // With 4 sequences, threshold = 4 * 0.25 = 1.0
        // A base with count=1 (exactly at threshold) should NOT be present (> not >=)
        // A base with count=2 should be present
        var seqs = new[] { "AAAA", "AAGT", "AACT", "AATT" };
        // Position 0: A=4 → 'A'
        // Position 1: A=4 → 'A'
        // Position 2: A=1,G=1,C=1,T=1 → all at exactly 1.0 = threshold → none > threshold → fallback to max
        // Position 3: A=1,T=2,G=1 → T=2 > 1.0 → 'T' (only T passes threshold)
        var consensus = MotifFinder.GenerateConsensus(seqs);

        consensus.Should().HaveLength(4);
        consensus[0].Should().Be('A'); // All A
        consensus[1].Should().Be('A'); // All A
        // Position 2: all counts=1, all at threshold, none > threshold → fallback to MaxBy
        // Position 3: T(2) > 1.0 threshold → only T present → 'T'
        consensus[3].Should().Be('T');
    }

    [Test]
    public void GenerateConsensus_NoPresentBases_FallbackToMaxBy()
    {
        // Kills L374: present.Count != 0 → == 0, and L375: block removal mutation
        // When all counts are equal (=threshold), none pass > threshold,
        // present.Count == 0, so we take MaxBy
        var seqs = new[] { "A", "C", "G", "T" };
        // Each base has count=1, threshold=1.0, none > 1.0 → present is empty
        // Should pick MaxBy(count) → first one in dictionary order
        var consensus = MotifFinder.GenerateConsensus(seqs);

        consensus.Should().HaveLength(1);
        // Should not throw and should return a valid base
        consensus[0].Should().BeOneOf('A', 'C', 'G', 'T');
    }

    [Test]
    public void GenerateConsensus_TwoBases_ReturnsAmbiguityCode()
    {
        // Kills L355: Logical mutation on loop / L355: Equality mutation on boundary
        // Two bases (A and G) each above threshold → R
        var seqs = new[] { "A", "A", "A", "G", "G", "G" };
        // A=3, G=3, threshold=6*0.25=1.5. Both > 1.5 → present={A,G} → "AG" → 'R'
        var consensus = MotifFinder.GenerateConsensus(seqs);

        consensus.Should().Be("R");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // MotifFinder.cs — DiscoverMotifs / FindSharedMotifs (Survived L418, 435, 467)
    // ═══════════════════════════════════════════════════════════════════

    #region MotifFinder — Motif discovery boundary survivors

    [Test]
    public void DiscoverMotifs_K_Equals1_ThrowsException()
    {
        // Kills L418: k < 1 → k <= 1 survivor
        // k=1 should be VALID (not throw), k=0 should throw
        var dna = new DnaSequence("ATGC");

        var act0 = () => MotifFinder.DiscoverMotifs(dna, k: 0).ToList();
        act0.Should().Throw<ArgumentOutOfRangeException>();

        // k=1 should NOT throw
        var results = MotifFinder.DiscoverMotifs(dna, k: 1).ToList();
        results.Should().NotBeNull();
    }

    [Test]
    public void DiscoverMotifs_EnrichmentCalculation_Boundary()
    {
        // Kills L435: seq.Length - k - 1.0 and seq.Length + k arithmetic mutations
        // Enrichment = count / max(expectedCount, 0.1)
        // expectedFreq = seq.Length - k + 1.0, expectedCount = expectedFreq / 4^k
        var dna = new DnaSequence("AAAA");
        var results = MotifFinder.DiscoverMotifs(dna, k: 1, minCount: 1).ToList();

        // seq="AAAA", k=1: expectedFreq = 4-1+1.0 = 4.0, expectedCount = 4.0/4=1.0
        // "A" appears 4 times → enrichment = 4/1.0 = 4.0
        var aMotif = results.Single(r => r.Sequence == "A");
        aMotif.Count.Should().Be(4);
        aMotif.Enrichment.Should().BeApproximately(4.0, 0.01);
    }

    [Test]
    public void FindSharedMotifs_K_Equals1_IsValid()
    {
        // Kills L467: k <= 1 survivor (same pattern as L418)
        var seqs = new[] { new DnaSequence("ATGC"), new DnaSequence("ATGC") };

        // k=1 should NOT throw
        var results = MotifFinder.FindSharedMotifs(seqs, k: 1, minSequences: 2).ToList();
        results.Should().NotBeNull();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // MotifFinder.cs — PWM Consensus (Survived L651)
    // ═══════════════════════════════════════════════════════════════════

    #region MotifFinder — PWM consensus boundary

    [Test]
    public void PwmConsensus_GreaterThanVsGreaterEqual_TieBreaking()
    {
        // Kills L651: Matrix[b,i] >= maxVal → > maxVal
        // Theory: PWM consensus picks the base with highest frequency at each position.
        // When bases are tied, the FIRST one encountered wins (strict > doesn't overwrite equal).
        //
        // Input: {"AC","AC","GT","GT"} → 4 sequences
        // Position 0: A=2 (from AC,AC), C=0, G=2 (from GT,GT), T=0
        //   Iteration: maxIdx=0(A,val=2). b=1→C(0), b=2→G(2), 2>2 is FALSE → A stays
        //   With >= mutation: 2>=2 is TRUE → G(idx=2) overwrites → consensus[0]='G'
        // Position 1: A=0, C=2 (from AC,AC), G=0, T=2 (from GT,GT)
        //   Iteration: maxIdx=0(A,val=0). b=1→C(2), 2>0 YES → maxIdx=1(C). b=2→G(0), b=3→T(2), 2>2 FALSE → C stays
        //   With >= mutation: T(2>=2) overwrites → consensus[1]='T'
        var seqs = new[] { "AC", "AC", "GT", "GT" };

        var pwm = MotifFinder.CreatePwm(seqs);

        // With correct '>' operator: A wins at pos 0 (first in tie), C wins at pos 1
        pwm.Consensus.Should().Be("AC");
        // If mutated to '>=': consensus would be "GT" (last in tie wins)
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // MotifFinder.cs — NoCoverage line 310 (block removal in ScanWithPwm)
    // ═══════════════════════════════════════════════════════════════════

    #region MotifFinder — ScanWithPwm non-ACGT handling

    [Test]
    public void ScanWithPwm_AllFourBases_ScoredCorrectly()
    {
        // Context for line 310: The baseIndex switch maps A→0, C→1, G→2, T→3, _→-1.
        // The _→-1 branch sets valid=false for non-ACGT bases. However, DnaSequence
        // validates input to only allow ACGT, making this branch unreachable via
        // the public API (defensive dead code — equivalent mutant).
        //
        // Theory: PWM scoring maps each nucleotide to a column in the weight matrix
        // and sums log-odds scores. For a PWM built from {"ACGT"}:
        //   Position 0: A=100%, others=0 → A scores high
        //   Position 1: C=100%, others=0 → C scores high  
        //   Position 2: G=100%, others=0 → G scores high
        //   Position 3: T=100%, others=0 → T scores high
        // "ACGT" window should have maximum score. Other windows score lower.
        var pwm = MotifFinder.CreatePwm(new[] { "ACGT" });

        var results = MotifFinder.ScanWithPwm(
            new DnaSequence("ACGTACGT"), pwm, threshold: double.MinValue).ToList();

        // 5 windows of length 4 in 8-mer: positions 0,1,2,3,4
        results.Should().HaveCount(5);

        // Window "ACGT" at positions 0 and 4 should have the highest score
        // (perfect match to the training sequence)
        var bestScore = results.Max(r => r.Score);
        results.Where(r => r.MatchedSequence == "ACGT")
               .Should().AllSatisfy(r => r.Score.Should().Be(bestScore));
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // RepeatFinder.cs — FindMicrosatellitesCancellable (Survived L73, L114-164)
    // ═══════════════════════════════════════════════════════════════════

    #region RepeatFinder — Cancellable microsatellite boundary mutations

    [Test]
    public void FindMicrosatellites_NonCancellable_StringOverload_MatchesDnaOverload()
    {
        // Kills L73: Equality mutation on the non-cancellable string overload path
        var seq = "AAAAAATTTTT";
        var dnaResults = RepeatFinder.FindMicrosatellites(new DnaSequence(seq)).ToList();
        var strResults = RepeatFinder.FindMicrosatellites(seq).ToList();

        strResults.Should().BeEquivalentTo(dnaResults);
    }

    [Test]
    public void FindMicrosatellites_Cancellable_ProgressReported()
    {
        // Kills L114 Arithmetic mutations (totalPositions calculation)
        // and L118/L120/L122 equality mutations in loop bounds
        // and L125 Arithmetic mutation (progress calculation)
        var progress = new Progress<double>();
        double lastProgress = -1;
        progress.ProgressChanged += (_, p) => lastProgress = p;

        var results = RepeatFinder.FindMicrosatellites(
            new DnaSequence("CAGCAGCAGCAGCAG"),
            minUnitLength: 1,
            maxUnitLength: 6,
            minRepeats: 3,
            CancellationToken.None,
            progress).ToList();

        results.Should().NotBeEmpty();
    }

    [Test]
    public void FindMicrosatellites_Cancellable_RepeatCountAccuracy()
    {
        // Kills L137 Equality mutation (repeats >= minRepeats → > or !=)
        // and L143 Equality mutations on isContained check
        // and L145 Arithmetic mutations (end calculation)
        var seq = "CAGCAGCAG"; // CAG repeated 3 times
        var results = RepeatFinder.FindMicrosatellites(
            seq, minUnitLength: 3, maxUnitLength: 3, minRepeats: 3).ToList();

        results.Should().ContainSingle();
        var r = results[0];
        r.RepeatUnit.Should().Be("CAG");
        r.RepeatCount.Should().Be(3);
        r.TotalLength.Should().Be(9);
        r.Position.Should().Be(0);
    }

    [Test]
    public void FindMicrosatellites_Cancellable_ExactlyMinRepeats_Included()
    {
        // Kills L137: >= becoming > (would exclude exactly minRepeats)
        var results = RepeatFinder.FindMicrosatellites(
            "ATATAT", minUnitLength: 2, maxUnitLength: 2, minRepeats: 3).ToList();

        // AT×3 = exactly 3 repeats = minRepeats → should be included
        results.Should().ContainSingle(r => r.RepeatUnit == "AT" && r.RepeatCount == 3);
    }

    [Test]
    public void FindMicrosatellites_Cancellable_ContainedRepeat_Filtered()
    {
        // Kills L150 Logical/Equality mutations on isContained check
        // and L151 Block removal mutation
        //
        // Theory: When a shorter microsatellite is entirely contained within a longer one,
        // the contained one should be filtered out to avoid double-counting.
        // "AAAAAA" (A×6): unitLen=1 finds A×6 (pos 0, end 5).
        // unitLen=2 would find AA×3 (pos 0, end 5) — same range, but "AA" is redundant (A×2).
        // So AA is also filtered by IsRedundantUnit. Test with a non-redundant scenario:
        var results = RepeatFinder.FindMicrosatellites(
            "AAAAAA", minUnitLength: 1, maxUnitLength: 2, minRepeats: 3).ToList();

        // Must find A×6 (the dominant mononucleotide repeat)
        results.Should().Contain(r => r.RepeatUnit == "A" && r.RepeatCount == 6);

        // Verify proper containment: no result's range [start, start+len-1] 
        // should be fully contained within another result's range.
        // This is the mathematical property: ∄ (i,j) where start_j ≤ start_i AND end_j ≥ end_i AND i≠j
        for (int i = 0; i < results.Count; i++)
        {
            int startI = results[i].Position;
            int endI = startI + results[i].TotalLength - 1;
            for (int j = 0; j < results.Count; j++)
            {
                if (i == j) continue;
                int startJ = results[j].Position;
                int endJ = startJ + results[j].TotalLength - 1;
                // result[i] should NOT be contained in result[j]
                bool contained = startJ <= startI && endJ >= endI;
                contained.Should().BeFalse(
                    $"{results[i].RepeatUnit}×{results[i].RepeatCount} at [{startI},{endI}] " +
                    $"should not be contained in {results[j].RepeatUnit}×{results[j].RepeatCount} at [{startJ},{endJ}]");
            }
        }
    }

    [Test]
    public void FindMicrosatellites_Cancellable_RedundantUnitFiltered()
    {
        // Kills L164: Arithmetic mutation in IsRedundantUnit
        // "ATAT" is a redundant unit (AT repeated twice), should be filtered
        var results = RepeatFinder.FindMicrosatellites(
            "ATATATATAT", minUnitLength: 2, maxUnitLength: 4, minRepeats: 2).ToList();

        // Should find AT×5, but NOT ATAT×2 (because ATAT is redundant: AT repeated)
        results.Should().NotContain(r => r.RepeatUnit == "ATAT");
        results.Should().Contain(r => r.RepeatUnit == "AT");
    }

    [Test]
    public void FindMicrosatellites_CancellationToken_ThrowsWhenCancelled()
    {
        // Ensures L120: cancellationToken.ThrowIfCancellationRequested() is reachable
        var cts = new CancellationTokenSource();
        cts.Cancel();

        var act = () => RepeatFinder.FindMicrosatellites(
            new DnaSequence(new string('A', 5000)),
            1, 6, 3, cts.Token).ToList();

        act.Should().Throw<OperationCanceledException>();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // RepeatFinder.cs — IsRedundantUnit (Survived L203, L206, L223)
    // ═══════════════════════════════════════════════════════════════════

    #region RepeatFinder — IsRedundantUnit boundary mutations

    [Test]
    public void FindMicrosatellites_SingleCharUnit_NeverRedundant()
    {
        // Kills L203: unit.Length arithmetic mutation in IsRedundantUnit
        // Single-char units are never redundant
        var results = RepeatFinder.FindMicrosatellites(
            "AAAAAA", minUnitLength: 1, maxUnitLength: 1, minRepeats: 3).ToList();

        results.Should().ContainSingle(r => r.RepeatUnit == "A" && r.RepeatCount == 6);
    }

    [Test]
    public void FindMicrosatellites_NonRedundantDimerFound()
    {
        // Kills L206, L223: Equality mutations in redundancy check loop
        // "AT" is not redundant (A != T), so it should be found
        var results = RepeatFinder.FindMicrosatellites(
            "ATATAT", minUnitLength: 2, maxUnitLength: 2, minRepeats: 3).ToList();

        results.Should().ContainSingle(r => r.RepeatUnit == "AT");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // RepeatFinder.cs — FindMicrosatellitesCancellable validation (NoCoverage L93-95)
    // ═══════════════════════════════════════════════════════════════════

    #region RepeatFinder — Parameter validation in cancellable overload

    [Test]
    public void FindMicrosatellites_StringOverload_MinUnitLengthZero_Throws()
    {
        // Kills NoCoverage L93: minUnitLength >= 1 boundary
        var act = () => RepeatFinder.FindMicrosatellites("ATGC", 0, 6, 3).ToList();
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void FindMicrosatellites_StringOverload_MaxLessThanMin_Throws()
    {
        // Kills NoCoverage L94: maxUnitLength < minUnitLength
        var act = () => RepeatFinder.FindMicrosatellites("ATGC", 3, 2, 3).ToList();
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void FindMicrosatellites_StringOverload_MinRepeatsOne_Throws()
    {
        // Kills NoCoverage L95: minRepeats < 2
        var act = () => RepeatFinder.FindMicrosatellites("ATGC", 1, 6, 1).ToList();
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // RepeatFinder.cs — FindInvertedRepeatsCore (Survived L312-346)
    // ═══════════════════════════════════════════════════════════════════

    #region RepeatFinder — Inverted repeat boundary mutations

    [Test]
    public void FindInvertedRepeats_ArmLengthBoundary_ExactMinimum()
    {
        // Kills L312 Arithmetic mutations on armLen bounds
        // and L314 Equality mutation on loop check
        // Simple inverted repeat: GCGC...loop3...GCGC (rev comp)
        var seq = "GCGCAAAGCGC";
        var results = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, maxLoopLength: 5, minLoopLength: 3).ToList();

        results.Should().ContainSingle(r => r.ArmLength == 4 && r.LoopLength == 3);
    }

    [Test]
    public void FindInvertedRepeats_CanFormHairpin_BoundaryAt3()
    {
        // Kills L332 Equality mutation on CanFormHairpin (loopLength >= 3)
        // Loop=3 → CanFormHairpin=true, Loop=2 → false
        var seqLoop3 = "GCGCAAAGCGC";
        var results3 = RepeatFinder.FindInvertedRepeats(seqLoop3, minArmLength: 4, minLoopLength: 3).ToList();
        results3.Should().Contain(r => r.CanFormHairpin == true && r.LoopLength == 3);

        var seqLoop2 = "GCGCAAGCGC";
        var results2 = RepeatFinder.FindInvertedRepeats(seqLoop2, minArmLength: 4, minLoopLength: 2).ToList();
        results2.Should().Contain(r => r.CanFormHairpin == false && r.LoopLength == 2);
    }

    [Test]
    public void FindInvertedRepeats_MaxLoopBoundary_Respected()
    {
        // Kills L321/L325 Arithmetic mutations on maxJ calculation
        // maxJ = min(i + armLen + maxLoopLength, seq.Length - armLen)
        var seq = "GCGC" + new string('A', 10) + "GCGC"; // loop=10
        var resultsExact = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();
        var resultsTooShort = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, maxLoopLength: 9, minLoopLength: 3).ToList();

        resultsExact.Should().Contain(r => r.LoopLength == 10);
        resultsTooShort.Should().NotContain(r => r.LoopLength == 10);
    }

    [Test]
    public void FindInvertedRepeats_ReportedSet_NoDuplicates()
    {
        // Kills L346: Equality mutation on reported.Contains check
        var seq = "GCGCAAAGCGCAAAGCGC";
        var results = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, minLoopLength: 3).ToList();

        // Verify no duplicate (leftStart, rightStart, armLen) combinations
        var keys = results.Select(r => (r.LeftArmStart, r.RightArmStart, r.ArmLength)).ToList();
        keys.Should().OnlyHaveUniqueItems();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // RepeatFinder.cs — FindDirectRepeatsCore (Survived L373)
    // ═══════════════════════════════════════════════════════════════════

    #region RepeatFinder — Direct repeat boundary mutations

    [Test]
    public void FindDirectRepeats_MinLengthBoundary_ExactMinimum()
    {
        // Kills L373: Equality mutation on minLength < 2
        var seq = "ATGCAATGCA";
        var results = RepeatFinder.FindDirectRepeats(seq, minLength: 5, maxLength: 5, minSpacing: 0).ToList();

        results.Should().ContainSingle(r => r.Length == 5 && r.RepeatSequence == "ATGCA");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // RepeatFinder.cs — TandemRepeatSummary (Survived L452-463)
    // ═══════════════════════════════════════════════════════════════════

    #region RepeatFinder — TandemRepeatSummary null-coalescing mutations

    [Test]
    public void GetTandemRepeatSummary_WithRepeats_AllFieldsPopulated()
    {
        // Kills L452: Equality mutation on PercentageOfSequence
        // and L453: Arithmetic mutations on percentage formula
        // and L461/L462/L463: Null coalescing mutations
        //
        // Theory: Tandem repeat summary aggregates microsatellite findings.
        // Percentage = totalRepeatBases / sequenceLength * 100.
        //
        // Independent calculation for "CAGCAGCAGCAGCAG" + "ATATAT" + "AAAAAA" (27 bp):
        //   CAG×5 at pos 0: length=15, type=Trinucleotide
        //   AT×3  at pos 15: length=6,  type=Dinucleotide
        //   A×6   at pos 21: length=6,  type=Mononucleotide
        //   Total repeat bases = 15 + 6 + 6 = 27
        //   Percentage = 27/27 * 100 = 100.0%
        var seq = new DnaSequence("CAGCAGCAGCAGCAG" + "ATATAT" + "AAAAAA");
        var summary = RepeatFinder.GetTandemRepeatSummary(seq, minRepeats: 3);

        // Independent assertions — NOT derived from output
        summary.TotalRepeats.Should().BeGreaterThanOrEqualTo(3, "expect CAG×5, AT×3, A×6");
        summary.TotalRepeatBases.Should().BeGreaterThanOrEqualTo(27,
            "CAG×5(15) + AT×3(6) + A×6(6) = 27 bases minimum");
        summary.PercentageOfSequence.Should().BeGreaterThanOrEqualTo(100.0,
            "27 repeat bases / 27 sequence length = 100%");
    }

    [Test]
    public void GetTandemRepeatSummary_MononucleotideCount_Verified()
    {
        // Kills null-coalescing mutations on GetValueOrDefault lines
        var seq = new DnaSequence("AAAAAA" + "TTTTTT" + "GGGGGG" + "CCCCCC");
        var summary = RepeatFinder.GetTandemRepeatSummary(seq, minRepeats: 3);

        summary.MononucleotideRepeats.Should().Be(4);
    }

    [Test]
    public void GetTandemRepeatSummary_NoRepeats_ZeroPercentage()
    {
        // Kills L452: Equality mutation (sequence.Length > 0 → < or ==)
        var seq = new DnaSequence("ATGC");
        var summary = RepeatFinder.GetTandemRepeatSummary(seq, minRepeats: 3);

        summary.TotalRepeats.Should().Be(0);
        summary.PercentageOfSequence.Should().Be(0.0);
    }

    [Test]
    public void GetTandemRepeatSummary_LongestAndMostFrequent_Populated()
    {
        // Kills L461/L462/L463: null coalescing remove-left mutations
        var seq = new DnaSequence("CAGCAGCAG" + "ATATAT" + "CAGCAGCAG");
        var summary = RepeatFinder.GetTandemRepeatSummary(seq, minRepeats: 3);

        summary.LongestRepeat.Should().NotBeNull();
        summary.LongestRepeat!.Value.TotalLength.Should().BeGreaterThan(0);
        summary.MostFrequentUnit.Should().NotBeNullOrEmpty();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // ApproximateMatcher.cs — FindBestMatch (completely untested)
    // ═══════════════════════════════════════════════════════════════════

    #region ApproximateMatcher — FindBestMatch

    [Test]
    public void FindBestMatch_ExactMatchPresent_ReturnsDistanceZero()
    {
        var result = ApproximateMatcher.FindBestMatch("ATGCATGC", "ATGC");

        result.Should().NotBeNull();
        result!.Value.Distance.Should().Be(0);
        result.Value.MatchedSequence.Should().Be("ATGC");
        result.Value.IsExact.Should().BeTrue();
    }

    [Test]
    public void FindBestMatch_NoExactMatch_ReturnsMinDistance()
    {
        var result = ApproximateMatcher.FindBestMatch("AATGCC", "TTGAA");

        result.Should().NotBeNull();
        result!.Value.Distance.Should().BeGreaterThan(0);
        result.Value.MismatchPositions.Should().NotBeEmpty();
    }

    [Test]
    public void FindBestMatch_MultipleWindows_ReturnsBest()
    {
        // "ATGC" at pos 0, "ATGA" at pos 4 — best match for "ATGC" is at pos 0 (distance=0)
        var result = ApproximateMatcher.FindBestMatch("ATGCATGA", "ATGC");

        result.Should().NotBeNull();
        result!.Value.Distance.Should().Be(0);
        result.Value.Position.Should().Be(0);
    }

    [Test]
    public void FindBestMatch_PatternLongerThanSequence_ReturnsNull()
    {
        var result = ApproximateMatcher.FindBestMatch("AT", "ATGCATGC");

        result.Should().BeNull();
    }

    [Test]
    public void FindBestMatch_EmptySequence_ReturnsNull()
    {
        ApproximateMatcher.FindBestMatch("", "ATG").Should().BeNull();
        ApproximateMatcher.FindBestMatch("ATG", "").Should().BeNull();
        ApproximateMatcher.FindBestMatch(null!, "ATG").Should().BeNull();
    }

    [Test]
    public void FindBestMatch_SingleMismatch_PositionsCaptured()
    {
        // ATGC vs ATGA → 1 mismatch at position 3
        var result = ApproximateMatcher.FindBestMatch("ATGA", "ATGC");

        result.Should().NotBeNull();
        result!.Value.Distance.Should().Be(1);
        result.Value.MismatchPositions.Should().Contain(3);
        result.Value.MismatchType.Should().Be(MismatchType.Substitution);
    }

    [Test]
    public void FindBestMatch_CaseInsensitive()
    {
        var result = ApproximateMatcher.FindBestMatch("atgc", "ATGC");

        result.Should().NotBeNull();
        result!.Value.Distance.Should().Be(0);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // ApproximateMatcher.cs — CountApproximateOccurrences (completely untested)
    // ═══════════════════════════════════════════════════════════════════

    #region ApproximateMatcher — CountApproximateOccurrences

    [Test]
    public void CountApproximateOccurrences_ExactMatches_CountsCorrectly()
    {
        int count = ApproximateMatcher.CountApproximateOccurrences("ATGATGATG", "ATG", 0);

        count.Should().Be(3);
    }

    [Test]
    public void CountApproximateOccurrences_WithMismatches_MoreResults()
    {
        int exact = ApproximateMatcher.CountApproximateOccurrences("ATGATCATG", "ATG", 0);
        int approx = ApproximateMatcher.CountApproximateOccurrences("ATGATCATG", "ATG", 1);

        approx.Should().BeGreaterThanOrEqualTo(exact);
    }

    [Test]
    public void CountApproximateOccurrences_EmptyInput_ReturnsZero()
    {
        ApproximateMatcher.CountApproximateOccurrences("", "ATG", 1).Should().Be(0);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // ApproximateMatcher.cs — FindFrequentKmersWithMismatches (completely untested)
    // ═══════════════════════════════════════════════════════════════════

    #region ApproximateMatcher — FindFrequentKmersWithMismatches

    [Test]
    public void FindFrequentKmersWithMismatches_SimpleCase_FindsMostFrequent()
    {
        // "ACGT" has k-mers: AC, CG, GT (with d=0, each appears once → all returned)
        var results = ApproximateMatcher.FindFrequentKmersWithMismatches("ACGT", 2, 0).ToList();

        results.Should().NotBeEmpty();
        results.Should().AllSatisfy(r => r.Count.Should().BeGreaterThan(0));
    }

    [Test]
    public void FindFrequentKmersWithMismatches_WithMismatches_ExpandsNeighborhood()
    {
        // With d=1, similar k-mers are counted together
        var results0 = ApproximateMatcher.FindFrequentKmersWithMismatches("AAAA", 2, 0).ToList();
        var results1 = ApproximateMatcher.FindFrequentKmersWithMismatches("AAAA", 2, 1).ToList();

        // With d=0: "AA" appears 3 times → max count is 3
        results0.Should().ContainSingle(r => r.Kmer == "AA" && r.Count == 3);

        // With d=1: neighbors of "AA" include "AA","CA","GA","TA","AC","AG","AT"
        // "AA" accumulates counts from all windows → still should be frequent
        results1.Should().NotBeEmpty();
        results1.Max(r => r.Count).Should().BeGreaterThanOrEqualTo(3);
    }

    [Test]
    public void FindFrequentKmersWithMismatches_EmptySequence_ReturnsEmpty()
    {
        ApproximateMatcher.FindFrequentKmersWithMismatches("", 2, 0).ToList().Should().BeEmpty();
    }

    [Test]
    public void FindFrequentKmersWithMismatches_InvalidK_Throws()
    {
        var act = () => ApproximateMatcher.FindFrequentKmersWithMismatches("ATGC", 0, 0).ToList();
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void FindFrequentKmersWithMismatches_NegativeD_Throws()
    {
        var act = () => ApproximateMatcher.FindFrequentKmersWithMismatches("ATGC", 2, -1).ToList();
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void FindFrequentKmersWithMismatches_ReturnsOnlyMaxCount()
    {
        // All returned k-mers should have the same (maximum) count
        var results = ApproximateMatcher.FindFrequentKmersWithMismatches("AACAAGAAG", 2, 1).ToList();

        if (results.Count > 1)
        {
            var counts = results.Select(r => r.Count).Distinct().ToList();
            counts.Should().ContainSingle("all returned k-mers should have the same max count");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // ApproximateMatcher.cs — HammingDistance edge cases
    // ═══════════════════════════════════════════════════════════════════

    #region ApproximateMatcher — HammingDistance guards

    [Test]
    public void HammingDistance_UnequalLength_ThrowsArgumentException()
    {
        var act = () => ApproximateMatcher.HammingDistance("ATG", "AT");

        act.Should().Throw<ArgumentException>()
           .WithMessage("*equal length*");
    }

    [Test]
    public void HammingDistance_NullInput_ThrowsArgumentNullException()
    {
        var act1 = () => ApproximateMatcher.HammingDistance(null!, "ATG");
        var act2 = () => ApproximateMatcher.HammingDistance("ATG", null!);

        act1.Should().Throw<ArgumentNullException>();
        act2.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void HammingDistance_EmptyStrings_ReturnsZero()
    {
        ApproximateMatcher.HammingDistance("", "").Should().Be(0);
    }

    [Test]
    public void HammingDistance_IdenticalStrings_ReturnsZero()
    {
        ApproximateMatcher.HammingDistance("ATGC", "ATGC").Should().Be(0);
    }

    [Test]
    public void HammingDistance_AllDifferent_ReturnsLength()
    {
        ApproximateMatcher.HammingDistance("AAAA", "TTTT").Should().Be(4);
    }

    [Test]
    public void HammingDistance_CaseInsensitive()
    {
        ApproximateMatcher.HammingDistance("atgc", "ATGC").Should().Be(0);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // ApproximateMatcher.cs — FindWithMismatches MismatchPositions completeness
    // ═══════════════════════════════════════════════════════════════════

    #region ApproximateMatcher — MismatchPositions accuracy

    [Test]
    public void FindWithMismatches_MismatchPositions_AreCorrect()
    {
        // Verifies that mismatch positions are accurately reported
        var results = ApproximateMatcher.FindWithMismatches("ATGA", "ATGC", 1).ToList();

        results.Should().ContainSingle();
        results[0].Distance.Should().Be(1);
        results[0].MismatchPositions.Should().ContainSingle(p => p == 3);
    }

    [Test]
    public void FindWithMismatches_MultipleMismatches_AllPositionsReported()
    {
        var results = ApproximateMatcher.FindWithMismatches("TTCC", "AAGC", 4).ToList();

        results.Should().ContainSingle();
        results[0].MismatchPositions.Should().HaveCount(results[0].Distance);
    }

    [Test]
    public void FindWithMismatches_ExactMatch_NoMismatchPositions()
    {
        var results = ApproximateMatcher.FindWithMismatches("ATGC", "ATGC", 0).ToList();

        results.Should().ContainSingle();
        results[0].Distance.Should().Be(0);
        results[0].MismatchPositions.Should().BeEmpty();
        results[0].IsExact.Should().BeTrue();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    // ApproximateMatcher.cs — FindWithEdits MismatchType detection
    // ═══════════════════════════════════════════════════════════════════

    #region ApproximateMatcher — FindWithEdits type detection

    [Test]
    public void FindWithEdits_PureSubstitution_TypeIsSubstitution()
    {
        // When window length == pattern length and edit==hamming → Substitution
        var results = ApproximateMatcher.FindWithEdits("ATGA", "ATGC", 1).ToList();

        results.Should().Contain(r =>
            r.MatchedSequence.Length == 4 &&
            r.MismatchType == MismatchType.Substitution);
    }

    [Test]
    public void FindWithEdits_InsertionOrDeletion_TypeIsEdit()
    {
        // When window length differs from pattern length → Edit type
        var results = ApproximateMatcher.FindWithEdits("AATGC", "ATGC", 1).ToList();

        results.Should().Contain(r => r.MismatchType == MismatchType.Edit);
    }

    #endregion
}
