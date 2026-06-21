namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Repeats area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
///
/// Construction convention: a known repeat structure is embedded inside
/// "A-free" neutral padding (a non-repetitive cycle of T/G/C). Because every
/// embedded repeat motif contains an 'A', the padding can never spuriously
/// reproduce it, so the detector's output is determined by the embedded signal
/// and the parameter cell.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Repeats")]
public class RepeatsCombinatorialTests
{
    /// <summary>Non-repetitive, A-free filler (period-4 cycle, so it forms no period-1..3 tandem repeat and no len-4..7 direct/palindrome motif containing 'A').</summary>
    private static string Pad(int n)
    {
        const string cycle = "TGTC";
        var sb = new System.Text.StringBuilder(n);
        while (sb.Length < n) sb.Append(cycle);
        return sb.ToString(0, Math.Max(0, n));
    }

    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-STR-001 — Microsatellite / short-tandem-repeat detection
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 13.
    // Dimensions: unitLen(3) × minReps(3) × seqLen(3). Full grid 3×3×3 = 27.
    //
    // Model (Benson 1999 TRF; STR definition): a microsatellite is a unit u of
    // length unitLen repeated consecutively ≥ minRepeats times. FindMicrosatellites
    // restricted to [unitLen, unitLen] reports the maximal run with its unit,
    // repeat count and total length.
    //
    // The combinatorial point: unitLen and minReps interact. A run of E copies is
    // reported iff E ≥ minReps, and the reported unit/length must be exact in every
    // (unitLen, seqLen) cell.
    // ═══════════════════════════════════════════════════════════════════════

    private const int EmbeddedReps = 6;

    private static (string Text, int Pos, string Unit) BuildStr(int unitLen, int seqLen, int reps)
    {
        string unit = unitLen switch { 1 => "A", 2 => "AC", _ => "ACG" };
        string run = string.Concat(Enumerable.Repeat(unit, reps));
        // 'T' spacers isolate the run: no embedded unit contains T, so the padding's
        // boundary base cannot rotate into the run and shift the detected unit.
        int padTotal = seqLen - run.Length - 2;
        int left = padTotal / 2;
        string text = Pad(left) + "T" + run + "T" + Pad(padTotal - left);
        return (text, left + 1, unit);
    }

    /// <summary>
    /// Pairwise grid: every (unitLen × minReps × seqLen) cell. A run of
    /// <see cref="EmbeddedReps"/> copies (≥ every tested minReps) is detected as
    /// exactly one microsatellite with the embedded unit, count and position.
    /// </summary>
    [Test, Combinatorial]
    public void RepStr_DetectsEmbeddedRun_ExactUnitAndCount(
        [Values(1, 2, 3)] int unitLen,
        [Values(3, 4, 5)] int minReps,
        [Values(30, 60, 120)] int seqLen)
    {
        var (text, pos, unit) = BuildStr(unitLen, seqLen, EmbeddedReps);

        var results = RepeatFinder.FindMicrosatellites(new DnaSequence(text), unitLen, unitLen, minReps).ToList();

        results.Should().ContainSingle("the padding is A-free so only the embedded run qualifies");
        var r = results[0];
        r.RepeatUnit.Should().Be(unit);
        r.RepeatCount.Should().Be(EmbeddedReps).And.BeGreaterThanOrEqualTo(minReps);
        r.Position.Should().Be(pos);
        r.TotalLength.Should().Be(EmbeddedReps * unitLen);
    }

    /// <summary>
    /// Interaction witness across the minReps axis: a run of E copies vanishes from
    /// the report once minReps exceeds E (the count threshold genuinely gates).
    /// </summary>
    [Test, Combinatorial]
    public void RepStr_RunBelowThreshold_NotReported([Values(1, 2, 3)] int unitLen)
    {
        var (text, _, _) = BuildStr(unitLen, 80, EmbeddedReps);
        RepeatFinder.FindMicrosatellites(new DnaSequence(text), unitLen, unitLen, EmbeddedReps + 1)
            .Should().BeEmpty($"a run of {EmbeddedReps} copies cannot satisfy minReps={EmbeddedReps + 1}");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-TANDEM-001 — Tandem repeats over a unit-length range
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 14.
    // Dimensions: minUnitLen(3) × maxUnitLen(3) × minReps(3) × seqLen(3).
    //             Full grid 3×3×3×3 = 81 cells.
    //
    // Model: FindMicrosatellites scans every unit length in [minUnitLen, maxUnitLen]
    // and reports runs of ≥ minReps. Contract: maxUnitLen ≥ minUnitLen (else
    // ArgumentOutOfRangeException). A run is reported iff its unit length lies in the
    // window AND its repeat count ≥ minReps.
    //
    // The combinatorial point: four interacting knobs. Invalid (maxUnitLen <
    // minUnitLen) cells must throw; valid cells must (a) only ever report units
    // within [minUnitLen, maxUnitLen] with count ≥ minReps and (b) surface the
    // embedded length-2 tandem exactly when 2 ∈ [minUnitLen, maxUnitLen] and minReps
    // ≤ embedded count.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void RepTandem_RangeAndThreshold(
        [Values(1, 2, 3)] int minUnitLen,
        [Values(2, 4, 6)] int maxUnitLen,
        [Values(3, 4, 5)] int minReps,
        [Values(40, 90, 160)] int seqLen)
    {
        // Embedded: the dinucleotide "AC" repeated 6 times in A-free padding.
        string run = string.Concat(Enumerable.Repeat("AC", EmbeddedReps));
        int padTotal = seqLen - run.Length;
        string text = Pad(padTotal / 2) + run + Pad(padTotal - padTotal / 2);
        var dna = new DnaSequence(text);

        if (maxUnitLen < minUnitLen)
        {
            Action act = () => RepeatFinder.FindMicrosatellites(dna, minUnitLen, maxUnitLen, minReps).ToList();
            act.Should().Throw<ArgumentOutOfRangeException>("maxUnitLen < minUnitLen violates the contract");
            return;
        }

        var results = RepeatFinder.FindMicrosatellites(dna, minUnitLen, maxUnitLen, minReps).ToList();

        foreach (var r in results)
        {
            r.RepeatUnit.Length.Should().BeInRange(minUnitLen, maxUnitLen);
            r.RepeatCount.Should().BeGreaterThanOrEqualTo(minReps);
            r.TotalLength.Should().Be(r.RepeatCount * r.RepeatUnit.Length);
        }

        // The per-result invariant above already guarantees no unit shorter than
        // minUnitLen leaks; here we pin that the embedded length-2 tandem is
        // surfaced exactly when the unit-length window admits it.
        bool embeddedVisible = minUnitLen <= 2 && maxUnitLen >= 2 && minReps <= EmbeddedReps;
        if (embeddedVisible)
            results.Should().Contain(r => r.RepeatUnit == "AC" && r.RepeatCount == EmbeddedReps,
                "the AC×6 tandem is in range and meets the count threshold");
        else
            results.Should().NotContain(r => r.RepeatUnit == "AC",
                "a length-2 unit must not be reported when minUnitLen excludes it");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-INV-001 — Inverted-repeat (hairpin) detection
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 15.
    // Dimensions: minArmLen(3) × maxGap(3) × seqLen(3). Full grid 3×3×3 = 27.
    //
    // Model: an inverted repeat is leftArm · loop · rightArm with
    // rightArm = reverseComplement(leftArm). FindInvertedRepeats reports such
    // structures with arm length ≥ minArmLen and loop length ≤ maxGap (the
    // maxLoopLength parameter).
    //
    // The combinatorial point: minArmLen and maxGap interact. EVERY result must
    // satisfy rightArm = revcomp(leftArm), armLen ≥ minArmLen and loop ≤ maxGap;
    // and the embedded arm-6/loop-4 hairpin appears exactly when maxGap ≥ 4 — the
    // loop ≤ maxGap invariant simultaneously proves its absence when maxGap < 4.
    // ═══════════════════════════════════════════════════════════════════════

    private const string IrArm = "ACGTGC";   // revcomp = GCACGT
    private const string IrLoop = "TTAA";    // loop length 4

    [Test, Combinatorial]
    public void RepInverted_ArmComplementAndLoopBound(
        [Values(3, 5, 6)] int minArmLen,
        [Values(2, 5, 10)] int maxGap,
        [Values(40, 80, 160)] int seqLen)
    {
        string core = IrArm + IrLoop + RevComp(IrArm);
        int padTotal = seqLen - core.Length;
        string text = Pad(padTotal / 2) + core + Pad(padTotal - padTotal / 2);
        var dna = new DnaSequence(text);

        var results = RepeatFinder.FindInvertedRepeats(dna, minArmLen, maxGap, minLoopLength: 0).ToList();

        foreach (var r in results)
        {
            r.RightArm.Should().Be(RevComp(r.LeftArm), "the right arm is the reverse complement of the left arm");
            r.ArmLength.Should().BeGreaterThanOrEqualTo(minArmLen);
            r.LoopLength.Should().BeLessThanOrEqualTo(maxGap);
        }

        if (maxGap >= IrLoop.Length)   // minArmLen ≤ 6 holds for every tested value
            results.Should().Contain(r => r.LeftArm == IrArm && r.RightArm == RevComp(IrArm) && r.LoopLength == IrLoop.Length,
                "the arm-6 / loop-4 hairpin fits within the gap bound");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-DIRECT-001 — Direct-repeat detection
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 16.
    // Dimensions: minLen(3) × maxGap(3) × seqLen(3). Full grid 3×3×3 = 27.
    //
    // Model: a direct repeat is the same substring occurring twice, R₁ at i and
    // R₂ at j (j > i + len), with spacing = j − i − len. FindDirectRepeats reports
    // copies of length ∈ [minLength, maxLength] separated by spacing ≥ minSpacing
    // (the gap knob).
    //
    // The combinatorial point: minLen and the spacing gate interact. Every result
    // is a genuine duplicate (R₁ = R₂) within bounds and spacing ≥ gate; the
    // embedded 7-mer duplicate with gap 5 appears iff minLen ≤ 7 and the gap gate
    // ≤ 5.
    // ═══════════════════════════════════════════════════════════════════════

    private const string DirectCopy = "ACGTGCA"; // length 7, contains 'A' → unique to embedding
    private const int DirectGap = 5;

    [Test, Combinatorial]
    public void RepDirect_DuplicateWithinBoundsAndGap(
        [Values(3, 5, 7)] int minLen,
        [Values(1, 3, 8)] int spacingGate,
        [Values(40, 90, 160)] int seqLen)
    {
        string core = DirectCopy + Pad(DirectGap) + DirectCopy;
        int padTotal = seqLen - core.Length;
        string text = Pad(padTotal / 2) + core + Pad(padTotal - padTotal / 2);
        var dna = new DnaSequence(text);

        var results = RepeatFinder.FindDirectRepeats(dna, minLen, maxLength: 20, minSpacing: spacingGate).ToList();

        foreach (var r in results)
        {
            text.Substring(r.FirstPosition, r.Length).Should().Be(r.RepeatSequence);
            text.Substring(r.SecondPosition, r.Length).Should().Be(r.RepeatSequence);
            r.Length.Should().BeGreaterThanOrEqualTo(minLen);
            r.SecondPosition.Should().BeGreaterThan(r.FirstPosition);
            r.Spacing.Should().Be(r.SecondPosition - r.FirstPosition - r.Length).And.BeGreaterThanOrEqualTo(spacingGate);
        }

        if (minLen <= DirectCopy.Length && spacingGate <= DirectGap)
            results.Should().Contain(r => r.RepeatSequence == DirectCopy && r.Spacing == DirectGap,
                "the embedded 7-mer duplicate with gap 5 satisfies this cell");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-PALIN-001 — DNA palindrome detection
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 17.
    // Dimensions: minLen(3) × maxLen(3) × seqLen(3). Full grid 3×3×3 = 27.
    //
    // Model: a DNA palindrome reads identically on both strands, i.e. equals its
    // own reverse complement (necessarily even length). FindPalindromes reports
    // even-length palindromes with length ∈ [minLength, maxLength]; minLength must
    // be even and ≥ 4 and maxLength ≥ minLength (else ArgumentOutOfRangeException).
    //
    // The combinatorial point: minLen and maxLen interact. Invalid (maxLen <
    // minLen) cells throw; valid cells report only self-reverse-complement even
    // windows within bounds, and the embedded EcoRI site GAATTC appears iff
    // minLen ≤ 6 ≤ maxLen.
    // ═══════════════════════════════════════════════════════════════════════

    private const string PalindromeSite = "GAATTC"; // EcoRI, self reverse-complement, length 6

    [Test, Combinatorial]
    public void RepPalindrome_SelfReverseComplementWithinBounds(
        [Values(4, 6, 8)] int minLen,
        [Values(6, 8, 12)] int maxLen,
        [Values(30, 60, 120)] int seqLen)
    {
        int padTotal = seqLen - PalindromeSite.Length;
        string text = Pad(padTotal / 2) + PalindromeSite + Pad(padTotal - padTotal / 2);
        var dna = new DnaSequence(text);

        if (maxLen < minLen)
        {
            Action act = () => RepeatFinder.FindPalindromes(dna, minLen, maxLen).ToList();
            act.Should().Throw<ArgumentOutOfRangeException>("maxLen < minLen violates the contract");
            return;
        }

        var results = RepeatFinder.FindPalindromes(dna, minLen, maxLen).ToList();

        foreach (var r in results)
        {
            (r.Length % 2).Should().Be(0, "DNA palindromes are even length");
            r.Length.Should().BeInRange(minLen, maxLen);
            r.Sequence.Should().Be(RevComp(r.Sequence), "a palindrome equals its own reverse complement");
            text.Substring(r.Position, r.Length).Should().Be(r.Sequence);
        }

        if (minLen <= PalindromeSite.Length && maxLen >= PalindromeSite.Length)
            results.Should().Contain(r => r.Sequence == PalindromeSite,
                "the EcoRI palindrome GAATTC is in range");
    }
}
