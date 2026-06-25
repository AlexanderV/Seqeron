using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology neoantigen candidate-peptide windowing unit — ONCO-NEO-001.
/// The unit under test is <see cref="OncologyAnalyzer.GenerateNeoantigenPeptides(string, char, int, int, int)"/>
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs, together with the
/// <see cref="OncologyAnalyzer.NeoantigenPeptide"/> record struct it emits.
///
/// This file is scoped STRICTLY to NEO-001 peptide-window generation. It does NOT touch MHC binding
/// classification (ONCO-MHC-001, row 110, <see cref="OncologyAnalyzer.ClassifyMhcBinding"/>) which is a
/// distinct learned-model unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts that the code NEVER fails in
/// an undisciplined way: no hang, no IndexOutOfRangeException, no negative-length Substring crash, and no
/// nonsense output. Every input must resolve to EITHER a well-defined, theory-correct set of windows OR a
/// *documented, intentional* outcome (an Argument*Exception for a malformed mutation spec).
/// For neoantigen windowing the headline hazards (all Malformed-Content edges of the windowing arithmetic)
/// are:
///   • a mutation AT A TERMINUS (position 1 or L): the spanning start range
///     s ∈ [max(0, p0−k+1), min(p0, L−k)] must clip to the protein bounds so NO window starts before 0 or
///     runs past L. The hazard is a negative start or a start+k > L feeding Substring → crash. The contract
///     is a TRUNCATED count of fully-valid windows (§2.2, §6.1).
///   • STOP-GAIN: a nonsense substitution to the stop sentinel '*'. The unit treats sequences as opaque
///     one-letter codes (no alphabet validation, §3.3) → '*' is a normal substituted residue; the windows
///     still tile correctly and carry '*' at the mutation offset. No crash, no special-casing collapse.
///   • NON-CODING / no-protein-consequence context: a mutation that yields NO full window — e.g. a protein
///     shorter than every requested length, or a single-residue range that fits no k-mer. The documented
///     result is an EMPTY list (§3.3 "if no length yields a window the result is empty"), never a throw,
///     never a partial/negative-length window.
///   • length &lt; 8 / k &gt; L: a requested length exceeding the protein length is SKIPPED (§3.3, §6.1) —
///     the negative Substring length (L−k &lt; 0) must never reach Substring. Shorter lengths still return.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-NEO-001 — Neoantigen candidate peptide window generation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 109.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content (невалідний контент).
///     Targets (checklist row 109): "mutation at protein terminus, stop-gain, non-coding, length&lt;8".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Target mapping: "mutation at protein terminus" = mutationPosition 1 or L (clip the start range);
/// "stop-gain" = substitution to '*' (premature stop, opaque-code handling, no crash on truncation);
/// "non-coding" = a mutation with no full-length peptide consequence ⇒ empty result; "length&lt;8" =
/// a peptide/protein context shorter than the requested length k ⇒ that length skipped (k &gt; L guard),
/// no negative-length Substring.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Neoantigen_Peptide_Generation.md (docs/algorithms/Oncology/Neoantigen_Peptide_Generation.md):
///   • Window spans the mutation iff start s ∈ [max(1, p−k+1), min(p, L−k+1)] (1-based)         (§2.2)
///   • INV-01: every peptide length ∈ [minLength, maxLength]                                     (§2.4)
///   • INV-02: every peptide spans the mutation: StartPosition + MutationOffset == p,
///       0 ≤ MutationOffset &lt; Length                                                          (§2.4)
///   • INV-03: mutant &amp; wild-type peptides equal length, differ at exactly the mutation offset (§2.4)
///   • INV-04: MutantPeptide[offset] == a (mutant residue); WildTypePeptide[offset] == P[p]      (§2.4)
///   • INV-05: interior mutation (≥ k−1 from both ends) ⇒ exactly k windows of length k          (§2.4)
///   • INV-06: output ordered by length ascending then start ascending                           (§2.4)
///   • §3.3 validation: null ⇒ ArgumentNullException; empty protein / non-substitution
///       (mutantResidue == WT) / minLength &lt; 1 / maxLength &lt; minLength ⇒ ArgumentException;
///       mutationPosition ∉ [1, L] ⇒ ArgumentOutOfRangeException; k &gt; L skipped; empty if none fit.
///   • §6.1: mutation at terminus ⇒ truncated count of fitting windows; requested length &gt; L skipped.
///   • §7.1 worked example: GenerateNeoantigenPeptides("MKTAYIAKQRSTVWLNDEFGH", 'C', 5) ⇒ 35 peptides
///       (5 per length 8..14); first 8-mer Mutant "MKTACIAK", WildType "MKTAYIAK", Start 1, Offset 4.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyNeoantigenFuzzTests
{
    private const string AminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    // ── Well-formed-peptide-set assertion helper ─────────────────────────────
    // Pins the full documented contract on EVERY returned window, for ANY input.
    // This is what stops a fuzz test from rubber-stamping a malformed window
    // (out-of-bounds substring, wrong offset, mismatched WT/mutant pair).
    private static void AssertWellFormedPeptides(
        IReadOnlyList<NeoantigenPeptide> peptides,
        string wildTypeProtein,
        char mutantResidue,
        int mutationPosition,
        int minLength,
        int maxLength)
    {
        peptides.Should().NotBeNull();
        int proteinLength = wildTypeProtein.Length;
        char wildTypeResidue = wildTypeProtein[mutationPosition - 1];

        foreach (var p in peptides)
        {
            // INV-01: length within the requested range.
            p.Length.Should().BeInRange(minLength, maxLength,
                "every peptide length ∈ [minLength, maxLength] (INV-01)");

            // Window must lie fully inside the protein (no clip leak / negative substring).
            p.StartPosition.Should().BeGreaterThanOrEqualTo(1, "1-based start, never < 1 (terminus clip)");
            (p.StartPosition + p.Length - 1).Should().BeLessThanOrEqualTo(proteinLength,
                "window must end inside the protein (no run past L)");

            // Both peptides are real substrings of the right protein at the right coordinates.
            p.MutantPeptide.Length.Should().Be(p.Length);
            p.WildTypePeptide.Length.Should().Be(p.Length);
            p.WildTypePeptide.Should().Be(wildTypeProtein.Substring(p.StartPosition - 1, p.Length),
                "WT peptide is the agretope at the same coordinates");

            // INV-02: spans the mutation.
            p.MutationOffset.Should().BeInRange(0, p.Length - 1, "0 ≤ offset < length (INV-02)");
            (p.StartPosition + p.MutationOffset).Should().Be(mutationPosition,
                "StartPosition + MutationOffset == mutationPosition (INV-02)");

            // INV-04: residues at the offset.
            p.MutantPeptide[p.MutationOffset].Should().Be(mutantResidue,
                "MutantPeptide[offset] == mutant residue (INV-04)");
            p.WildTypePeptide[p.MutationOffset].Should().Be(wildTypeResidue,
                "WildTypePeptide[offset] == wild-type residue (INV-04)");

            // INV-03: equal length, differ at EXACTLY the offset.
            int diffCount = 0;
            for (int i = 0; i < p.Length; i++)
            {
                if (p.MutantPeptide[i] != p.WildTypePeptide[i])
                {
                    diffCount++;
                    i.Should().Be(p.MutationOffset, "the only differing position is the mutation offset (INV-03)");
                }
            }
            diffCount.Should().Be(1, "mutant & WT differ at exactly one position (INV-03)");
        }

        // INV-06: ordered by length ascending then start ascending.
        for (int i = 1; i < peptides.Count; i++)
        {
            var prev = peptides[i - 1];
            var cur = peptides[i];
            bool ordered = prev.Length < cur.Length
                || (prev.Length == cur.Length && prev.StartPosition < cur.StartPosition);
            ordered.Should().BeTrue("output ordered by length asc then start asc (INV-06)");
        }
    }

    // Documented spanning-window count for length k around 1-based position p in a protein of length L:
    // |[max(1, p−k+1), min(p, L−k+1)]|, or 0 if k > L.
    private static int ExpectedWindowCount(int proteinLength, int mutationPosition, int k)
    {
        if (k > proteinLength)
        {
            return 0;
        }

        int firstStart = Math.Max(1, mutationPosition - k + 1);
        int lastStart = Math.Min(mutationPosition, proteinLength - k + 1);
        return Math.Max(0, lastStart - firstStart + 1);
    }

    private static char DifferentResidue(Random rng, char from)
    {
        char c;
        do
        {
            c = AminoAcids[rng.Next(AminoAcids.Length)];
        }
        while (c == from);
        return c;
    }

    #region ONCO-NEO-001 — Positive sanity (documented worked example & interior counts)

    [Test]
    public void GenerateNeoantigenPeptides_WorkedExample_MatchesDocumentedWindows()
    {
        // Docs §7.1: MKTAYIAKQRSTVWLNDEFGH (L=21), Y5C, default 8..14 ⇒ 35 peptides (5 per length).
        const string wt = "MKTAYIAKQRSTVWLNDEFGH";
        var peptides = GenerateNeoantigenPeptides(wt, 'C', 5);

        peptides.Should().HaveCount(35, "5 windows × 7 lengths (8..14) for an interior residue near the N-terminus");
        AssertWellFormedPeptides(peptides, wt, 'C', 5, 8, 14);

        // Exactly 5 per length (p=5 keeps the right clamp at 5 for every k≤14 since 21-k+1 ≥ 8 ≥ 5).
        for (int k = 8; k <= 14; k++)
        {
            peptides.Count(p => p.Length == k).Should().Be(5, $"length {k} yields 5 spanning windows for Y5C in L=21");
        }

        // First 8-mer is the documented one.
        var first8 = peptides.First(p => p.Length == 8);
        first8.StartPosition.Should().Be(1);
        first8.MutationOffset.Should().Be(4);
        first8.MutantPeptide.Should().Be("MKTACIAK");
        first8.WildTypePeptide.Should().Be("MKTAYIAK");
    }

    [Test]
    public void GenerateNeoantigenPeptides_InteriorMutation_YieldsExactlyKWindowsPerLength()
    {
        // INV-05: an interior mutation (≥ k−1 from both ends) yields exactly k windows of length k.
        var rng = new Random(109_001);
        for (int trial = 0; trial < 200; trial++)
        {
            int proteinLength = rng.Next(31, 60);
            var chars = new char[proteinLength];
            for (int i = 0; i < proteinLength; i++)
            {
                chars[i] = AminoAcids[rng.Next(AminoAcids.Length)];
            }
            string wt = new(chars);

            // Pick a deeply interior position (≥ 14 from both ends, so all k∈8..14 are interior).
            int pos = rng.Next(15, proteinLength - 13); // 1-based; comfortably interior.
            char mut = DifferentResidue(rng, wt[pos - 1]);

            var peptides = GenerateNeoantigenPeptides(wt, mut, pos);
            AssertWellFormedPeptides(peptides, wt, mut, pos, 8, 14);

            for (int k = 8; k <= 14; k++)
            {
                peptides.Count(p => p.Length == k).Should().Be(k,
                    $"interior mutation ⇒ exactly k={k} windows of length {k} (INV-05)");
            }
        }
    }

    [Test]
    public void GenerateNeoantigenPeptides_SingleLengthRange_ReturnsOnlyThatLength()
    {
        // §6.1: min==max ⇒ only that length.
        const string wt = "MKTAYIAKQRSTVWLNDEFGH";
        var peptides = GenerateNeoantigenPeptides(wt, 'C', 5, minLength: 9, maxLength: 9);

        peptides.Should().OnlyContain(p => p.Length == 9);
        peptides.Should().NotBeEmpty();
        AssertWellFormedPeptides(peptides, wt, 'C', 5, 9, 9);
    }

    #endregion

    #region ONCO-NEO-001 — MC: mutation at protein terminus (clip start range, no out-of-bounds)

    [Test]
    [CancelAfter(30_000)]
    public void GenerateNeoantigenPeptides_MutationAtEitherTerminus_ClipsWindows_NoCrash()
    {
        // The headline MC hazard: a mutation at position 1 or L. The start range must clip to the bounds,
        // producing a TRUNCATED count of fully-valid windows — never a negative start or start+k > L
        // feeding Substring (IndexOutOfRange / negative length).
        var rng = new Random(109_002);
        for (int trial = 0; trial < 300; trial++)
        {
            int proteinLength = rng.Next(8, 40);
            var chars = new char[proteinLength];
            for (int i = 0; i < proteinLength; i++)
            {
                chars[i] = AminoAcids[rng.Next(AminoAcids.Length)];
            }
            string wt = new(chars);

            // Alternate between the N-terminus (pos 1) and the C-terminus (pos L).
            int pos = (trial % 2 == 0) ? 1 : proteinLength;
            char mut = DifferentResidue(rng, wt[pos - 1]);

            var act = () => GenerateNeoantigenPeptides(wt, mut, pos);
            act.Should().NotThrow("terminus mutation must clip, never crash (§6.1)");
            var peptides = act();

            AssertWellFormedPeptides(peptides, wt, mut, pos, 8, 14);

            // The count matches the documented clipped range for each length.
            for (int k = 8; k <= 14; k++)
            {
                peptides.Count(p => p.Length == k).Should().Be(ExpectedWindowCount(proteinLength, pos, k),
                    $"clipped window count for k={k} at terminus pos={pos}, L={proteinLength}");
            }
        }
    }

    [Test]
    public void GenerateNeoantigenPeptides_NTerminusExactlyKLong_SingleWindowStartsAtOne()
    {
        // L == k and mutation at position 1: exactly one window (start 1), offset 0 — the tightest N-terminus
        // case. A negative firstStart (p0−k+1 = −7) must be clamped to 0.
        const string wt = "MKTAYIAK"; // L = 8.
        var peptides = GenerateNeoantigenPeptides(wt, 'W', 1, minLength: 8, maxLength: 8);

        peptides.Should().HaveCount(1);
        peptides[0].StartPosition.Should().Be(1);
        peptides[0].MutationOffset.Should().Be(0);
        peptides[0].MutantPeptide.Should().Be("WKTAYIAK");
        AssertWellFormedPeptides(peptides, wt, 'W', 1, 8, 8);
    }

    [Test]
    public void GenerateNeoantigenPeptides_CTerminusExactlyKLong_SingleWindowOffsetAtEnd()
    {
        // L == k and mutation at position L: one window, offset k−1 (last residue). lastStart = min(p0, L−k)
        // must clamp to 0, not run past the end.
        const string wt = "MKTAYIAK"; // L = 8.
        var peptides = GenerateNeoantigenPeptides(wt, 'W', 8, minLength: 8, maxLength: 8);

        peptides.Should().HaveCount(1);
        peptides[0].StartPosition.Should().Be(1);
        peptides[0].MutationOffset.Should().Be(7);
        peptides[0].MutantPeptide.Should().Be("MKTAYIAW");
        AssertWellFormedPeptides(peptides, wt, 'W', 8, 8, 8);
    }

    #endregion

    #region ONCO-NEO-001 — MC: stop-gain (substitution to '*', opaque one-letter code, no special collapse)

    [Test]
    [CancelAfter(30_000)]
    public void GenerateNeoantigenPeptides_StopGainSubstitution_TilesNormally_NoCrash()
    {
        // Stop-gain: the mutant residue is the stop sentinel '*'. §3.3: sequences are opaque one-letter codes
        // with NO alphabet validation; '*' is treated as any other substituted residue. The windows must
        // still tile correctly and carry '*' at the mutation offset — no crash on the "truncation", no
        // silent collapse to an empty set.
        var rng = new Random(109_003);
        for (int trial = 0; trial < 200; trial++)
        {
            int proteinLength = rng.Next(12, 50);
            var chars = new char[proteinLength];
            for (int i = 0; i < proteinLength; i++)
            {
                chars[i] = AminoAcids[rng.Next(AminoAcids.Length)];
            }
            string wt = new(chars);
            int pos = rng.Next(1, proteinLength + 1);
            // WT residue is from AminoAcids (never '*'), so '*' is always a genuine substitution.

            var act = () => GenerateNeoantigenPeptides(wt, '*', pos);
            act.Should().NotThrow("stop-gain '*' is an opaque substitution, no crash (§3.3)");
            var peptides = act();

            AssertWellFormedPeptides(peptides, wt, '*', pos, 8, 14);
            foreach (var p in peptides)
            {
                p.MutantPeptide[p.MutationOffset].Should().Be('*', "stop-gain residue carried at the offset");
            }
        }
    }

    #endregion

    #region ONCO-NEO-001 — MC: non-coding / no peptide consequence (empty result, never a throw)

    [Test]
    public void GenerateNeoantigenPeptides_ProteinShorterThanAllLengths_ReturnsEmpty_NoCrash()
    {
        // "Non-coding" / no peptide consequence: a protein shorter than every requested length. Each k > L is
        // skipped (§3.3) ⇒ documented EMPTY result, never a throw, never a negative-length Substring.
        var rng = new Random(109_004);
        for (int len = 1; len <= 7; len++)
        {
            var chars = new char[len];
            for (int i = 0; i < len; i++)
            {
                chars[i] = AminoAcids[rng.Next(AminoAcids.Length)];
            }
            string wt = new(chars);
            int pos = rng.Next(1, len + 1);
            char mut = DifferentResidue(rng, wt[pos - 1]);

            var act = () => GenerateNeoantigenPeptides(wt, mut, pos); // default 8..14, all > L.
            act.Should().NotThrow($"protein length {len} < min length 8 ⇒ empty, not a crash");
            act().Should().BeEmpty($"no 8..14-mer fits in a length-{len} protein (§3.3)");
        }
    }

    [Test]
    public void GenerateNeoantigenPeptides_SingleResidueProtein_DefaultLengths_Empty()
    {
        // The extreme non-coding edge: a 1-residue protein. The only valid position is 1; no k∈8..14 fits.
        var peptides = GenerateNeoantigenPeptides("M", 'K', 1);
        peptides.Should().BeEmpty();
    }

    #endregion

    #region ONCO-NEO-001 — MC: length<8 / k>L (skip the length, never a negative-length substring)

    [Test]
    public void GenerateNeoantigenPeptides_RangeStraddlingProteinLength_SkipsOversizedLengths()
    {
        // length<8 target: a protein of length 9 with the default 8..14 range. k=8,9 fit; k=10..14 (> L) must
        // be skipped (the L−k < 0 substring length must never reach Substring). Shorter lengths still return.
        const string wt = "MKTAYIAKQ"; // L = 9.
        var peptides = GenerateNeoantigenPeptides(wt, 'C', 5); // default 8..14.

        peptides.Should().OnlyContain(p => p.Length == 8 || p.Length == 9,
            "k=10..14 (> L=9) are skipped (§3.3, §6.1)");
        peptides.Should().Contain(p => p.Length == 8);
        peptides.Should().Contain(p => p.Length == 9);
        AssertWellFormedPeptides(peptides, wt, 'C', 5, 8, 14);
    }

    [Test]
    [CancelAfter(30_000)]
    public void GenerateNeoantigenPeptides_AllRequestedLengthsExceedProtein_Empty()
    {
        // Fuzz the k>L guard: random short proteins with a length range entirely above L ⇒ empty, no crash.
        var rng = new Random(109_005);
        for (int trial = 0; trial < 200; trial++)
        {
            int proteinLength = rng.Next(1, 10);
            var chars = new char[proteinLength];
            for (int i = 0; i < proteinLength; i++)
            {
                chars[i] = AminoAcids[rng.Next(AminoAcids.Length)];
            }
            string wt = new(chars);
            int pos = rng.Next(1, proteinLength + 1);
            char mut = DifferentResidue(rng, wt[pos - 1]);

            int min = proteinLength + 1 + rng.Next(0, 5);
            int max = min + rng.Next(0, 5);

            var act = () => GenerateNeoantigenPeptides(wt, mut, pos, min, max);
            act.Should().NotThrow("all k > L ⇒ each skipped, no negative-length Substring");
            act().Should().BeEmpty("no length in the range fits the protein (§3.3)");
        }
    }

    [Test]
    public void GenerateNeoantigenPeptides_ShortLengthBelowEight_StillTilesWhenFits()
    {
        // length<8 also means SUB-8 peptide lengths are LEGAL (minLength ≥ 1 per §3.1). A length-3 window on a
        // short protein must tile correctly — exercising the small-k path the default 8..11 never reaches.
        const string wt = "ACDEF"; // L = 5.
        var peptides = GenerateNeoantigenPeptides(wt, 'W', 3, minLength: 3, maxLength: 3);

        // pos 3, k 3 ⇒ starts s ∈ [max(1,1), min(3,3)] = [1,3] ⇒ 3 windows.
        peptides.Should().HaveCount(3);
        AssertWellFormedPeptides(peptides, wt, 'W', 3, 3, 3);
        peptides.Select(p => p.StartPosition).Should().Equal(1, 2, 3);
    }

    #endregion

    #region ONCO-NEO-001 — MC/BE: documented validation surface (malformed mutation spec ⇒ documented throw)

    [Test]
    public void GenerateNeoantigenPeptides_NullProtein_ThrowsArgumentNull()
    {
        var act = () => GenerateNeoantigenPeptides(null!, 'C', 1);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void GenerateNeoantigenPeptides_EmptyProtein_ThrowsArgument()
    {
        var act = () => GenerateNeoantigenPeptides(string.Empty, 'C', 1);
        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void GenerateNeoantigenPeptides_NonSubstitution_MutantEqualsWildType_ThrowsArgument()
    {
        // §3.3: mutantResidue == WT residue is not a missense ⇒ ArgumentException.
        const string wt = "MKTAYIAKQRST";
        var act = () => GenerateNeoantigenPeptides(wt, 'Y', 5); // position 5 is already 'Y'.
        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void GenerateNeoantigenPeptides_PositionOutOfRange_ThrowsArgumentOutOfRange()
    {
        const string wt = "MKTAYIAKQRST"; // L = 12.
        var rng = new Random(109_006);
        for (int trial = 0; trial < 50; trial++)
        {
            // Positions ≤ 0 or > L.
            int pos = rng.Next(0, 2) == 0 ? -rng.Next(0, 100) : 13 + rng.Next(0, 100);
            var act = () => GenerateNeoantigenPeptides(wt, 'C', pos);
            act.Should().Throw<ArgumentOutOfRangeException>($"position {pos} ∉ [1, 12]");
        }
    }

    [Test]
    public void GenerateNeoantigenPeptides_InvalidLengthRange_ThrowsArgument()
    {
        const string wt = "MKTAYIAKQRST";
        ((Action)(() => GenerateNeoantigenPeptides(wt, 'C', 5, minLength: 0, maxLength: 8)))
            .Should().Throw<ArgumentException>("minLength < 1 (§3.3)");
        ((Action)(() => GenerateNeoantigenPeptides(wt, 'C', 5, minLength: 9, maxLength: 8)))
            .Should().Throw<ArgumentException>("maxLength < minLength (§3.3)");
    }

    #endregion
}
