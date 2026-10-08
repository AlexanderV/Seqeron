namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Tests for RepeatFinder.FindDirectRepeats (Test Unit REP-DIRECT-001).
/// 
/// Direct repeats are identical sequences appearing multiple times in the same orientation.
/// Example: 5' TTACG------TTACG 3' where ------ is the spacing region.
/// 
/// Reporting convention: maximal repeated pairs (left- and right-maximal; Gusfield 1997 §7.12),
/// i.e. the forward-strand output of MUMmer repeat-match -f (Kurtz et al. 2004; mummer4 source
/// src/tigr/repeat-match.cc compiled and run), filtered by maxLength and Spacing ≥ minSpacing.
/// Only A/C/G/T match (MUMmer mummer -n). Expected values below were produced by repeat-match -f
/// and an independent brute-force reference of the definition (0 mismatches on 8000 random cases).
///
/// Sources:
/// - Gusfield (1997) Algorithms on Strings, Trees and Sequences §7.12 (maximal pairs)
/// - Kurtz et al. (2004) Genome Biol 5:R12 — MUMmer 3 repeat-match; Kurtz &amp; Schleiermacher (1999) REPuter
/// - Wikipedia: Direct repeat, Repeated sequence (DNA)
/// - Ussery et al. (2009): Computing for Comparative Microbial Genomics
/// - Richard (2021): PMC8145212 - Trinucleotide repeat expansions
/// </summary>
[TestFixture]
public class RepeatFinder_DirectRepeat_Tests
{
    #region MUST Tests - Core Algorithm

    /// <summary>
    /// M1: Core algorithm detects identical sequences at two positions.
    /// Evidence: Wikipedia - Direct repeat definition.
    /// </summary>
    [Test]
    public void FindDirectRepeats_SimpleRepeat_FindsMatchingPair()
    {
        // Arrange: "ACGTA" appears at position 0 and position 9 with spacing of 4
        var sequence = new DnaSequence("ACGTATTTTACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(results[0].RepeatSequence, Is.EqualTo("ACGTA"));
            Assert.That(results[0].FirstPosition, Is.EqualTo(0));
            Assert.That(results[0].SecondPosition, Is.EqualTo(9));
            Assert.That(results[0].Length, Is.EqualTo(5));
            Assert.That(results[0].Spacing, Is.EqualTo(4));
        });
    }

    /// <summary>
    /// M2: Adjacent repeats (minSpacing=0) should be detected.
    /// Evidence: Wikipedia - tandem direct repeats.
    /// </summary>
    [Test]
    public void FindDirectRepeats_AdjacentRepeats_WithZeroSpacing_Found()
    {
        // Arrange: "ACGTA" immediately followed by "ACGTA"
        var sequence = new DnaSequence("ACGTAACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(results[0].RepeatSequence, Is.EqualTo("ACGTA"));
            Assert.That(results[0].FirstPosition, Is.EqualTo(0));
            Assert.That(results[0].SecondPosition, Is.EqualTo(5));
            Assert.That(results[0].Spacing, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// M3: Sequence without repeated patterns returns empty.
    /// </summary>
    [Test]
    public void FindDirectRepeats_NoRepeats_ReturnsEmpty()
    {
        // Arrange: No 5+ bp pattern repeats
        var sequence = new DnaSequence("ACGTACGT");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M4: Empty input returns empty enumerable.
    /// </summary>
    [Test]
    public void FindDirectRepeats_EmptySequence_ReturnsEmpty()
    {
        // Act
        var results = RepeatFinder.FindDirectRepeats("", 5, 10, 1).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    #endregion

    #region MUST Tests - Parameter Validation

    /// <summary>
    /// M5: Null DnaSequence throws ArgumentNullException.
    /// </summary>
    [Test]
    public void FindDirectRepeats_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindDirectRepeats((DnaSequence)null!, 5, 10, 1).ToList());
    }

    /// <summary>
    /// M6: minLength less than 2 throws exception.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinLengthTooSmall_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindDirectRepeats(sequence, 1, 10, 1).ToList());
    }

    /// <summary>
    /// M7: maxLength less than minLength throws exception.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MaxLengthLessThanMinLength_ThrowsException()
    {
        var sequence = new DnaSequence("ACGTACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindDirectRepeats(sequence, 10, 5, 1).ToList());
    }

    #endregion

    #region MUST Tests - Invariants

    /// <summary>
    /// M8: Spacing = SecondPosition - FirstPosition - Length.
    /// </summary>
    [Test]
    public void FindDirectRepeats_SpacingCalculation_MatchesFormula()
    {
        // Arrange: Various repeats with known spacing
        var sequence = new DnaSequence("ACGTATTTTACGTA");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert: Verify invariant for all results
        foreach (var result in results)
        {
            int expectedSpacing = result.SecondPosition - result.FirstPosition - result.Length;
            Assert.That(result.Spacing, Is.EqualTo(expectedSpacing),
                $"Invariant violated for repeat at {result.FirstPosition}->{result.SecondPosition}");
        }
    }

    /// <summary>
    /// M9: FirstPosition is always less than SecondPosition.
    /// </summary>
    [Test]
    public void FindDirectRepeats_FirstPosition_AlwaysLessThanSecondPosition()
    {
        var sequence = new DnaSequence("ACGTAACGTATTTTACGTA");

        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        foreach (var result in results)
        {
            Assert.That(result.FirstPosition, Is.LessThan(result.SecondPosition),
                $"FirstPosition ({result.FirstPosition}) should be < SecondPosition ({result.SecondPosition})");
        }
    }

    /// <summary>
    /// M10: RepeatSequence equals actual substring at FirstPosition.
    /// </summary>
    [Test]
    public void FindDirectRepeats_RepeatSequence_MatchesSubstringAtFirstPosition()
    {
        var sequenceStr = "ACGTATTTTACGTA";
        var sequence = new DnaSequence(sequenceStr);

        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        foreach (var result in results)
        {
            string actual = sequenceStr.Substring(result.FirstPosition, result.Length);
            Assert.That(result.RepeatSequence, Is.EqualTo(actual),
                $"RepeatSequence should match substring at FirstPosition");
        }
    }

    #endregion

    #region MUST Tests - Filter Thresholds

    /// <summary>
    /// M11: Only repeats with Length &gt;= minLength are returned.
    /// "ACGTA" copies at 0 and 7 (maximal, length 5); "ACGT" at 14 pairs with 0 and 7 at length 4.
    /// repeat-match -f -n 4: (1,8,5), (1,15,4), (8,15,4) in 1-based coordinates.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinLength_RespectsThreshold()
    {
        const string seq = "ACGTATTACGTAGGACGTC";

        var min5 = RepeatFinder.FindDirectRepeats(seq, 5, 10, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
        var min4 = RepeatFinder.FindDirectRepeats(seq, 4, 10, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();

        Assert.That(min5, Is.EqualTo(new[] { (0, 7, 5) }));
        Assert.That(min4, Is.EqualTo(new[] { (0, 7, 5), (0, 14, 4), (7, 14, 4) }));
    }

    /// <summary>
    /// M12: Only repeats with Length &lt;= maxLength are returned. A maximal repeat longer than
    /// maxLength is NOT truncated into sub-windows — it is simply not reported.
    /// ACGTACGTAC+TTTT+ACGTACGTAC: maximal pairs (0,14,10), (0,18,6), (3,13,7) (brute force = repeat-match -f).
    /// </summary>
    [Test]
    public void FindDirectRepeats_MaxLength_RespectsThreshold()
    {
        var repeat = "ACGTACGTAC"; // 10bp
        var sequence = new DnaSequence(repeat + "TTTT" + repeat);

        var capped = RepeatFinder.FindDirectRepeats(sequence, 5, 8, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
        var full = RepeatFinder.FindDirectRepeats(sequence, 5, 50, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();

        Assert.That(capped, Is.EqualTo(new[] { (0, 18, 6), (3, 13, 7) }));
        Assert.That(full, Is.EqualTo(new[] { (0, 14, 10), (0, 18, 6), (3, 13, 7) }));
    }

    /// <summary>
    /// M13: Only repeats with Spacing >= minSpacing are returned.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinSpacing_RespectsThreshold()
    {
        // Arrange: Adjacent repeats with spacing=0
        var sequence = new DnaSequence("ACGTAACGTA");

        // Act: minSpacing=1 should exclude adjacent repeats
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 1).ToList();

        // Assert: Should be empty because only adjacent repeat exists
        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M14: Sequence shorter than 2×minLength returns empty.
    /// </summary>
    [Test]
    public void FindDirectRepeats_SequenceTooShort_ReturnsEmpty()
    {
        // Arrange: 8bp sequence, minLength=5 requires at least 10bp (5+5)
        var sequence = new DnaSequence("ACGTACGT");

        // Act
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 10, 0).ToList();

        // Assert
        Assert.That(results, Is.Empty);
    }

    #endregion

    #region SHOULD Tests

    /// <summary>
    /// S1: Three copies of a pattern with distinct flanks produce all three pairwise (maximal) pairs.
    /// Evidence: Wikipedia - "nucleotide sequences present in multiple copies"; Gusfield maximal pairs.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MultipleOccurrences_FindsAllPairs()
    {
        // "ACGTA" at 0, 7, 14; left flanks (start, T, G) and right flanks (T, G, C) all differ.
        var results = RepeatFinder.FindDirectRepeats(new DnaSequence("ACGTATTACGTAGGACGTACC"), 5, 5, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.RepeatSequence)).ToList();

        Assert.That(results, Is.EqualTo(new[] { (0, 7, "ACGTA"), (0, 14, "ACGTA"), (7, 14, "ACGTA") }));
    }

    /// <summary>
    /// S1b: Periodic copies ("ACGTATT" period 7) — copies 0/7 and 7/14 are sub-windows of ONE maximal
    /// overlapping repeat (0,7,12); only the (0,14,5) pair is left- and right-maximal.
    /// repeat-match -f -n 4 on ACGTATTACGTATTACGTA: (1,8,12), (1,15,5).
    /// </summary>
    [Test]
    public void FindDirectRepeats_PeriodicCopies_ReportsMaximalPairsOnly()
    {
        const string seq = "ACGTATTACGTATTACGTA";

        var spaced = RepeatFinder.FindDirectRepeats(seq, 4, 50, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
        var all = RepeatFinder.FindDirectRepeats(seq, 4, 50, int.MinValue)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Spacing)).ToList();

        Assert.That(spaced, Is.EqualTo(new[] { (0, 14, 5) }));
        Assert.That(all, Is.EqualTo(new[] { (0, 7, 12, -5), (0, 14, 5, 9) }));
    }

    /// <summary>
    /// S2: String overload produces consistent results with DnaSequence overload.
    /// </summary>
    [Test]
    public void FindDirectRepeats_StringOverload_MatchesDnaSequenceOverload()
    {
        const string seq = "ACGTAACGTA";
        var dnaSequence = new DnaSequence(seq);

        var stringResults = RepeatFinder.FindDirectRepeats(seq, 5, 10, 0).ToList();
        var dnaResults = RepeatFinder.FindDirectRepeats(dnaSequence, 5, 10, 0).ToList();

        Assert.That(stringResults, Has.Count.EqualTo(dnaResults.Count));
        for (int i = 0; i < stringResults.Count; i++)
        {
            Assert.That(stringResults[i], Is.EqualTo(dnaResults[i]));
        }
    }

    /// <summary>
    /// S3: Lowercase input is handled correctly (case-insensitive).
    /// </summary>
    [Test]
    public void FindDirectRepeats_LowercaseInput_HandledCorrectly()
    {
        // Arrange: lowercase sequence — same as "ACGTAACGTA" after normalization
        var results = RepeatFinder.FindDirectRepeats("acgtaacgta", 5, 10, 0).ToList();

        // Assert: Exact same result as uppercase input
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(results[0].RepeatSequence, Is.EqualTo("ACGTA"));
            Assert.That(results[0].FirstPosition, Is.EqualTo(0));
            Assert.That(results[0].SecondPosition, Is.EqualTo(5));
        });
    }

    /// <summary>
    /// S4: Repeats with large intervening region are detected.
    /// Evidence: Wikipedia - interspersed repeats.
    /// </summary>
    [Test]
    public void FindDirectRepeats_LongSpacing_Detected()
    {
        // Arrange: 20bp spacing between repeats, using alternating pattern
        var spacing = "ATGATGATGATGATGATGAT"; // 20bp alternating pattern
        var repeat = "CCCGGGCCC"; // 9bp unique repeat (not found in spacing)
        var sequence = new DnaSequence($"{repeat}{spacing}{repeat}");

        // Act: Use exact repeat length to avoid subpattern matches
        var results = RepeatFinder.FindDirectRepeats(sequence, 9, 9, 1).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].Spacing, Is.EqualTo(20));
        Assert.That(results[0].RepeatSequence, Is.EqualTo(repeat));
    }

    /// <summary>
    /// S5: CAG trinucleotide repeat (Huntington's disease related).
    /// Evidence: Richard (2021) - trinucleotide repeat disorders.
    /// </summary>
    [Test]
    public void FindDirectRepeats_BiologicalRepeat_TrinucleotideCAG()
    {
        // Arrange: CAG repeat with spacing (not tandem)
        var sequence = new DnaSequence("CAGCAGTTTTTTCAGCAG");

        // Act: Look for 6bp pattern (CAGCAG)
        var results = RepeatFinder.FindDirectRepeats(sequence, 6, 6, 1).ToList();

        // Assert
        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].RepeatSequence, Is.EqualTo("CAGCAG"));
    }

    #endregion

    #region COULD Tests

    /// <summary>
    /// C1: Homopolymer flanks — every left- and right-maximal pair is reported once at its full length;
    /// nested sub-windows (e.g. AAAA at 1..5 vs 11..15) are not. repeat-match -f -n 4 (1-based, spacing ≥ 1):
    /// (1,11,6), (1,12,5), (2,11,5), (1,13,4), (3,11,4). Previously 14 (i,j,len) windows were reported.
    /// </summary>
    [Test]
    public void FindDirectRepeats_OverlappingPatterns_OnlyMaximalPairs()
    {
        var sequence = new DnaSequence("AAAAAATTTTAAAAAA");

        var results = RepeatFinder.FindDirectRepeats(sequence, 4, 6, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();

        Assert.That(results, Is.EqualTo(new[] { (0, 10, 6), (0, 11, 5), (0, 12, 4), (1, 10, 5), (2, 10, 4) }));
    }

    /// <summary>
    /// C2: Performance baseline — O(n²) algorithm completes in bounded time.
    /// Evidence: DoD requirement for O(n²) complexity.
    /// </summary>
    [Test]
    public void FindDirectRepeats_LargeSequence_CompletesInReasonableTime()
    {
        // Arrange: 1000bp pseudo-random sequence (deterministic seed)
        var random = new Random(42);
        var bases = new[] { 'A', 'C', 'G', 'T' };
        var sb = new System.Text.StringBuilder(1000);
        for (int i = 0; i < 1000; i++)
            sb.Append(bases[random.Next(4)]);
        var sequence = new DnaSequence(sb.ToString());

        // Act: Should complete well within time bound for O(n²)
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var results = RepeatFinder.FindDirectRepeats(sequence, 5, 20, 1).ToList();
        sw.Stop();

        // Assert: valid output + anti-hang guard (generous bound catches O(2^n)/infinite-loop
        // blowups without flaking under parallel-suite CPU load; perf is measured by the [Explicit]
        // benchmark fixtures, not by a wall-clock budget in a correctness test).
        Assert.That(sw.ElapsedMilliseconds, Is.LessThan(30_000),
            $"O(n²) on 1000bp must not blow up, took {sw.ElapsedMilliseconds}ms");
        Assert.That(results, Is.Not.Null);

        // Verify invariants hold on all results
        foreach (var result in results)
        {
            Assert.That(result.Length, Is.GreaterThanOrEqualTo(5));
            Assert.That(result.Length, Is.LessThanOrEqualTo(20));
            Assert.That(result.FirstPosition, Is.LessThan(result.SecondPosition));
        }
    }

    #endregion

    #region Reference cross-check (MUMmer repeat-match -f) and edge cases

    /// <summary>
    /// Full maximal-pair lists (minSpacing = int.MinValue, i.e. overlaps admitted) equal the output of the
    /// compiled MUMmer 4 repeat-match -f -n L (converted to 0-based, sorted).
    /// </summary>
    [TestCase("ACGTACGTTTTTTTTTACGTACGT", 4,
        "0,4,4;0,16,8;0,20,4;3,15,5;7,8,8;7,9,7;7,10,6;7,11,5;7,12,4;15,19,5")]
    [TestCase("ACGTACGTACGT", 2, "0,4,8;0,8,4")]
    [TestCase("AAAAAATTTTAAAAAA", 5, "0,1,5;0,10,6;0,11,5;1,10,5;10,11,5")]
    [TestCase("ACGTATTTTACGTA", 5, "0,9,5")]
    public void FindDirectRepeats_AllMaximalPairs_MatchRepeatMatch(string seq, int minLength, string expected)
    {
        var actual = RepeatFinder.FindDirectRepeats(seq, minLength, int.MaxValue, int.MinValue)
            .Select(r => $"{r.FirstPosition},{r.SecondPosition},{r.Length}");
        Assert.That(string.Join(";", actual), Is.EqualTo(expected));
    }

    /// <summary>
    /// minSpacing filters maximal pairs by Spacing = j − i − L; negative values admit overlap but a
    /// position is never paired with itself (previously minSpacing = −L produced (i, i) self-pairs).
    /// </summary>
    [Test]
    public void FindDirectRepeats_NegativeMinSpacing_NoSelfPairs()
    {
        var results = RepeatFinder.FindDirectRepeats("ACGTACGTACGT", 4, 4, -4)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Spacing)).ToList();

        Assert.That(results, Is.EqualTo(new[] { (0, 8, 4, 4) }));
        Assert.That(RepeatFinder.FindDirectRepeats("ACGTACGTACGT", 2, 50, int.MaxValue), Is.Empty);
    }

    /// <summary>
    /// Non-ACGT symbols never match (MUMmer mummer -n "match only the characters a, c, g, or t").
    /// </summary>
    [Test]
    public void FindDirectRepeats_NonAcgt_NeverMatches()
    {
        Assert.That(RepeatFinder.FindDirectRepeats("NNNNNNNNNNNNNNNNNNNN", 5, 5, 1), Is.Empty);
        Assert.That(RepeatFinder.FindDirectRepeats("acgtannnnnacgta", 5, 50, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.RepeatSequence)),
            Is.EqualTo(new[] { (0, 10, "ACGTA") }));
        // N inside both copies of ACGNTACG splits them: N≠N, so (0,10) stops at 3 and (4,14) is left-maximal.
        Assert.That(RepeatFinder.FindDirectRepeats("ACGNTACGGGACGNTACG", 3, 50, 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)),
            Is.EqualTo(new[] { (0, 5, 3), (0, 10, 3), (0, 15, 3), (4, 14, 4), (5, 10, 3), (10, 15, 3) }));
    }

    /// <summary>
    /// Output is sorted by (FirstPosition, SecondPosition) and every position pair occurs once.
    /// </summary>
    [Test]
    public void FindDirectRepeats_Output_SortedAndUniquePerPositionPair()
    {
        var rng = new Random(7);
        var seq = new string(Enumerable.Range(0, 3000).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        var results = RepeatFinder.FindDirectRepeats(seq, 6, 50, 0).ToList();

        Assert.That(results.Select(r => (r.FirstPosition, r.SecondPosition)).Distinct().Count(), Is.EqualTo(results.Count));
        Assert.That(results, Is.Ordered.By(nameof(DirectRepeatResult.FirstPosition))
            .Then.By(nameof(DirectRepeatResult.SecondPosition)));
        foreach (var r in results)
        {
            Assert.That(seq.Substring(r.SecondPosition, r.Length), Is.EqualTo(r.RepeatSequence));
            Assert.That(r.FirstPosition == 0 || seq[r.FirstPosition - 1] != seq[r.SecondPosition - 1], "left-maximal");
            Assert.That(r.SecondPosition + r.Length == seq.Length || seq[r.FirstPosition + r.Length] != seq[r.SecondPosition + r.Length], "right-maximal");
        }
    }

    #endregion
}