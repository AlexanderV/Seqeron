using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Codon area — codon optimization (CODON-OPT-001) and the
/// Codon Adaptation Index (CODON-CAI-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// KeyNotFoundException, NullReferenceException, OverflowException, …). Every
/// input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception (ArgumentException /
/// ArgumentNullException). A raw runtime exception or a hang on garbage input is
/// a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CODON-OPT-001 — codon optimization (Codon)
/// Checklist: docs/checklists/03_FUZZING.md, row 58.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — non-coding sequences (no recognizable codons,
///          input shorter than one codon), non-DNA / non-nucleotide characters
///          that produce un-translatable codons.
///   • BE = Boundary Exploitation — empty string, null, a single base (&lt; 3),
///          input length not a multiple of 3 (a trailing partial codon), and an
///          extremely long sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The codon-optimization contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Codon optimization replaces each codon with a synonymous codon preferred by
/// the target organism (to improve heterologous expression) WITHOUT changing the
/// encoded protein.[Sharp &amp; Li 1987; Plotkin &amp; Kudla 2011]
///   — docs/algorithms/Codon_Optimization/Sequence_Optimization.md §1, §2.2.
///
/// API entry: CodonOptimizer.OptimizeSequence(string codingSequence,
///   CodonUsageTable targetOrganism, OptimizationStrategy strategy = ...,
///   double gcTargetMin = 0.40, double gcTargetMax = 0.60,
///   double rareCodonThreshold = 0.15)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs).
///
/// The method is LENIENT by documented design — it does NOT validate the
/// alphabet; instead it normalizes, trims, and translates whatever it is given:
///   • null OR empty input → an all-empty OptimizationResult with zero-valued
///     metrics; NOT an exception (explicit early return, CodonOptimizer.cs
///     lines 247–250). Sequence_Optimization.md §3.3, §6.1.
///   • input is upper-cased and `T` is replaced by `U` (RNA notation), so DNA
///     and lowercase input round-trip identically (CodonOptimizer.cs line 252).
///   • input length NOT divisible by 3 → the trailing partial codon is TRIMMED
///     away before optimization (`rna = rna.Substring(0, (rna.Length/3)*3)`,
///     CodonOptimizer.cs lines 254–258). This is the KEY no-crash boundary: a
///     final partial codon must be ignored, NEVER an IndexOutOfRangeException.
///     Sequence_Optimization.md §6.1 ("Incomplete final codon → Trimmed away").
///   • a sequence shorter than one codon (length &lt; 3) → SplitIntoCodons yields
///     ZERO codons (the loop guard is `i + 2 < length`, CodonOptimizer.cs lines
///     687–695) → an empty optimized sequence and empty protein; no codon is ever
///     indexed out of range.
///   • a codon that is NOT in the standard genetic code (because it contains a
///     non-DNA character, or any symbol other than A/C/G/U) → TranslateCodon
///     returns the sentinel "X" (GetValueOrDefault default, CodonOptimizer.cs
///     line 699), NOT a KeyNotFoundException. SelectOptimalCodon then finds no
///     synonymous set for "X" and returns the codon UNCHANGED (lines 323–324).
///     So non-DNA input is carried through verbatim — never a crash, never a
///     KeyNotFound, and never a wrong-length result.
///
/// KEY INVARIANT (INV-01, Sequence_Optimization.md §2.4): the optimized sequence
/// encodes the SAME protein as the (normalized, trimmed) input — replacement
/// codons are drawn only from the synonymous set of the original amino acid and
/// stop codons are preserved. This fuzz suite pins INV-01 directly by translating
/// BOTH the input codons and the OptimizedSequence with the standard genetic code
/// and asserting equality, on every input that is exercised.
///
/// LENGTH INVARIANTS we also pin so a malformed input cannot silently corrupt the
/// output length:
///   • INV-02: `OptimizedSequence.Length % 3 == 0` (only complete codons survive).
///   • INV-03: `OriginalSequence.Length == OptimizedSequence.Length` (codon-for-
///     codon replacement preserves length; the stored OriginalSequence is the
///     trimmed, normalized RNA input).
/// — Sequence_Optimization.md §2.4 (INV-02, INV-03), §3.2.
///
/// Determinism note: every test uses a FIXED input (no shared static Rng) and
/// avoids the `HarmonizeExpression` strategy, which performs weighted RANDOM
/// codon selection (CodonOptimizer.cs lines 361–376, documented non-deterministic
/// in Sequence_Optimization.md §5.3). The deterministic strategies
/// (MaximizeCAI, BalancedOptimization, AvoidRareCodeons) are used throughout, so
/// every assertion is reproducible.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class CodonFuzzTests
{
    #region Helpers

    /// <summary>
    /// The standard genetic code (RNA codon → one-letter amino acid, '*' = stop),
    /// mirroring CodonOptimizer's internal table. Used to INDEPENDENTLY translate
    /// a sequence so the protein-preservation invariant (INV-01) can be checked
    /// against the optimizer's output without relying on the optimizer itself.
    /// Unknown / malformed codons map to the sentinel "X" (matching
    /// CodonOptimizer.TranslateCodon's GetValueOrDefault default).
    /// </summary>
    private static readonly Dictionary<string, string> StandardGeneticCode = new()
    {
        { "UUU", "F" }, { "UUC", "F" },
        { "UUA", "L" }, { "UUG", "L" }, { "CUU", "L" }, { "CUC", "L" }, { "CUA", "L" }, { "CUG", "L" },
        { "AUU", "I" }, { "AUC", "I" }, { "AUA", "I" },
        { "AUG", "M" },
        { "GUU", "V" }, { "GUC", "V" }, { "GUA", "V" }, { "GUG", "V" },
        { "UCU", "S" }, { "UCC", "S" }, { "UCA", "S" }, { "UCG", "S" }, { "AGU", "S" }, { "AGC", "S" },
        { "CCU", "P" }, { "CCC", "P" }, { "CCA", "P" }, { "CCG", "P" },
        { "ACU", "T" }, { "ACC", "T" }, { "ACA", "T" }, { "ACG", "T" },
        { "GCU", "A" }, { "GCC", "A" }, { "GCA", "A" }, { "GCG", "A" },
        { "UAU", "Y" }, { "UAC", "Y" },
        { "UAA", "*" }, { "UAG", "*" }, { "UGA", "*" },
        { "CAU", "H" }, { "CAC", "H" },
        { "CAA", "Q" }, { "CAG", "Q" },
        { "AAU", "N" }, { "AAC", "N" },
        { "AAA", "K" }, { "AAG", "K" },
        { "GAU", "D" }, { "GAC", "D" },
        { "GAA", "E" }, { "GAG", "E" },
        { "UGU", "C" }, { "UGC", "C" },
        { "UGG", "W" },
        { "CGU", "R" }, { "CGC", "R" }, { "CGA", "R" }, { "CGG", "R" }, { "AGA", "R" }, { "AGG", "R" },
        { "GGU", "G" }, { "GGC", "G" }, { "GGA", "G" }, { "GGG", "G" }
    };

    /// <summary>The deterministic target organism table used by every test.</summary>
    private static readonly CodonOptimizer.CodonUsageTable Target = CodonOptimizer.EColiK12;

    /// <summary>The deterministic strategies (HarmonizeExpression is RANDOM → excluded).</summary>
    private static readonly CodonOptimizer.OptimizationStrategy[] DeterministicStrategies =
    {
        CodonOptimizer.OptimizationStrategy.MaximizeCAI,
        CodonOptimizer.OptimizationStrategy.BalancedOptimization,
        CodonOptimizer.OptimizationStrategy.AvoidRareCodeons,
        CodonOptimizer.OptimizationStrategy.MinimizeSecondary,
    };

    /// <summary>
    /// Independently translates a normalized RNA sequence codon-by-codon, exactly
    /// as CodonOptimizer does (3-base step, partial trailing bases ignored, unknown
    /// codon → "X"). The result is the protein this sequence encodes.
    /// </summary>
    private static string Translate(string rna)
    {
        var protein = new System.Text.StringBuilder();
        for (int i = 0; i + 2 < rna.Length; i += 3)
            protein.Append(StandardGeneticCode.GetValueOrDefault(rna.Substring(i, 3), "X"));
        return protein.ToString();
    }

    /// <summary>Normalizes a raw coding string to the optimizer's RNA notation.</summary>
    private static string Normalize(string coding) =>
        coding.ToUpperInvariant().Replace('T', 'U');

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CODON-OPT-001 — codon optimization : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CODON-OPT-001 — codon optimization

    #region Positive sanity — a coding sequence optimizes, preserving the protein

    /// <summary>
    /// Positive sanity: a well-formed coding sequence (M-A-L-R-Stop) optimizes
    /// without error and PRESERVES the encoded protein across every deterministic
    /// strategy (INV-01), keeps the codon part a multiple of 3 (INV-02), and keeps
    /// OriginalSequence.Length == OptimizedSequence.Length (INV-03). This is the
    /// baseline that proves the fuzz targets below are measured against a working
    /// happy path, not a uniformly-broken method.
    /// </summary>
    [Test]
    public void OptimizeSequence_WellFormedCoding_PreservesProteinAndLength()
    {
        const string coding = "ATGGCTCTACGTTAA"; // M  A  L  R  *
        string expectedProtein = Translate(Normalize(coding));
        expectedProtein.Should().Be("MALR*", "the input encodes Met-Ala-Leu-Arg-Stop");

        foreach (var strategy in DeterministicStrategies)
        {
            var result = CodonOptimizer.OptimizeSequence(coding, Target, strategy);

            // INV-02 / INV-03 — the codon part stays a clean multiple of 3 and
            // length is preserved codon-for-codon.
            (result.OptimizedSequence.Length % 3).Should().Be(0,
                $"[{strategy}] only complete codons survive (INV-02)");
            result.OriginalSequence.Length.Should().Be(result.OptimizedSequence.Length,
                $"[{strategy}] codon-for-codon replacement preserves length (INV-03)");

            // INV-01 (KEY) — synonymous optimization must not change the protein.
            Translate(result.OptimizedSequence).Should().Be(expectedProtein,
                $"[{strategy}] codon optimization is synonymous: translate(optimized) == translate(input)");
            result.ProteinSequence.Should().Be(expectedProtein,
                $"[{strategy}] the reported protein matches the input translation");
        }
    }

    #endregion

    #region BE — Boundary: empty string / null (no codons at all)

    /// <summary>
    /// BE: empty string and null are the "no input" boundary. The documented
    /// contract returns an all-empty OptimizationResult with zero-valued metrics
    /// via an explicit early return (CodonOptimizer.cs lines 247–250) — NOT a
    /// NullReferenceException, NOT a crash. Pins that "no sequence" is a defined
    /// no-op, not an error.
    /// </summary>
    [TestCase(null, TestName = "OptimizeSequence_Null_IsEmptyResultNoThrow")]
    [TestCase("", TestName = "OptimizeSequence_Empty_IsEmptyResultNoThrow")]
    public void OptimizeSequence_NullOrEmpty_ReturnsEmptyResult(string? input)
    {
        CodonOptimizer.OptimizationResult result = default;
        var act = () => result = CodonOptimizer.OptimizeSequence(input!, Target);

        act.Should().NotThrow("null/empty is a defined no-op early return, not an error");
        result.OptimizedSequence.Should().BeEmpty("an empty input optimizes to nothing");
        result.OriginalSequence.Should().BeEmpty();
        result.ProteinSequence.Should().BeEmpty();
        result.ChangedCodons.Should().Be(0);
        result.OriginalCAI.Should().Be(0);
        result.OptimizedCAI.Should().Be(0);
    }

    #endregion

    #region BE/MC — Boundary: sequence shorter than one codon (non-coding, < 3)

    /// <summary>
    /// BE/MC: a sequence shorter than a single codon has NO complete codon. The
    /// optimizer trims to the largest multiple-of-3 prefix — which is empty — so
    /// SplitIntoCodons yields zero codons (loop guard `i + 2 < length`,
    /// CodonOptimizer.cs lines 687–695). The result is an empty optimized sequence
    /// and empty protein, with NO IndexOutOfRangeException from indexing a partial
    /// codon. Covers length 1 and length 2, valid bases and non-DNA, in both cases.
    /// </summary>
    [TestCase("A", TestName = "OptimizeSequence_SingleBase_NoCodon_IsEmpty")]
    [TestCase("AT", TestName = "OptimizeSequence_TwoBases_NoCodon_IsEmpty")]
    [TestCase("g", TestName = "OptimizeSequence_SingleLowercaseBase_NoCodon_IsEmpty")]
    [TestCase("Z", TestName = "OptimizeSequence_SingleNonDnaChar_NoCodon_IsEmpty")]
    [TestCase("ZZ", TestName = "OptimizeSequence_TwoNonDnaChars_NoCodon_IsEmpty")]
    public void OptimizeSequence_ShorterThanOneCodon_TrimsToEmpty(string input)
    {
        CodonOptimizer.OptimizationResult result = default;
        var act = () => result = CodonOptimizer.OptimizeSequence(input, Target);

        act.Should().NotThrow(
            "a sub-codon input must be trimmed to zero codons, never indexed out of range");
        result.OptimizedSequence.Should().BeEmpty(
            "no complete codon exists, so nothing is optimized");
        result.ProteinSequence.Should().BeEmpty("no codon translates to an amino acid");
        (result.OptimizedSequence.Length % 3).Should().Be(0, "INV-02 holds trivially at length 0");
    }

    #endregion

    #region BE — Boundary: length NOT divisible by 3 (trailing partial codon)

    /// <summary>
    /// BE (KEY no-crash boundary): when the input length is not a multiple of 3,
    /// the trailing partial codon is TRIMMED before optimization
    /// (CodonOptimizer.cs lines 254–258) — it must NEVER trigger an
    /// IndexOutOfRangeException on the final 1–2 leftover bases. The complete
    /// codons that remain are optimized and the protein over them is preserved
    /// (INV-01); the output length equals the trimmed input length and stays a
    /// multiple of 3 (INV-02, INV-03). Verified for both a +1 and a +2 remainder,
    /// across every deterministic strategy.
    /// </summary>
    [TestCase("ATGGCTA", TestName = "OptimizeSequence_LenMod3Is1_TrimsTrailingBase")]   // 7 = 2 codons + 1
    [TestCase("ATGGCTAT", TestName = "OptimizeSequence_LenMod3Is2_TrimsTrailingTwo")]   // 8 = 2 codons + 2
    public void OptimizeSequence_LengthNotDivisibleBy3_TrimsPartialCodonAndPreservesProtein(string input)
    {
        // The protein the optimizer should see = translation of the trimmed,
        // normalized prefix (the trailing partial codon contributes nothing).
        string expectedProtein = Translate(Normalize(input));

        foreach (var strategy in DeterministicStrategies)
        {
            CodonOptimizer.OptimizationResult result = default;
            var act = () => result = CodonOptimizer.OptimizeSequence(input, Target, strategy);

            act.Should().NotThrow(
                $"[{strategy}] a trailing partial codon must be trimmed, never cause IndexOutOfRange");

            (result.OptimizedSequence.Length % 3).Should().Be(0,
                $"[{strategy}] the partial codon is dropped, so the output is whole codons (INV-02)");
            result.OriginalSequence.Length.Should().Be(result.OptimizedSequence.Length,
                $"[{strategy}] trimmed input and optimized output have equal length (INV-03)");
            result.OriginalSequence.Length.Should().Be((Normalize(input).Length / 3) * 3,
                $"[{strategy}] OriginalSequence is the input trimmed to complete codons");

            Translate(result.OptimizedSequence).Should().Be(expectedProtein,
                $"[{strategy}] optimization over the complete codons preserves the protein (INV-01)");
        }
    }

    #endregion

    #region MC — Malformed Content: non-DNA chars in otherwise-codon-aligned input

    /// <summary>
    /// MC: non-DNA characters embedded in a length-multiple-of-3 sequence form
    /// codons that are NOT in the standard genetic code. The optimizer must map
    /// each such codon to the sentinel "X" (TranslateCodon's GetValueOrDefault
    /// default, CodonOptimizer.cs line 699) — NEVER a KeyNotFoundException — and,
    /// finding no synonymous set for "X", leave that codon UNCHANGED
    /// (SelectOptimalCodon lines 323–324). So the malformed codons pass through
    /// verbatim, the recognizable codons may be optimized, and the protein over
    /// the WHOLE sequence is still preserved (INV-01). Covers digits, gap, an
    /// embedded null byte, the ambiguity code N, and unicode (Greek, astral
    /// surrogate pair) — each as the middle codon of a 3-codon input.
    /// </summary>
    [TestCase("ATG123TAA", TestName = "OptimizeSequence_NonDna_Digits_NoCrashProteinPreserved")]
    [TestCase("ATG-GCTAA", TestName = "OptimizeSequence_NonDna_GapDash_NoCrashProteinPreserved")]
    [TestCase("ATG\0\0\0TAA", TestName = "OptimizeSequence_NonDna_NullBytes_NoCrashProteinPreserved")]
    [TestCase("ATGNNNTAA", TestName = "OptimizeSequence_NonDna_AmbiguityN_NoCrashProteinPreserved")]
    [TestCase("ATGαβγTAA", TestName = "OptimizeSequence_NonDna_GreekLetters_NoCrashProteinPreserved")]
    [TestCase("ATG😀xTAA", TestName = "OptimizeSequence_NonDna_AstralSurrogatePair_NoCrashProteinPreserved")]
    public void OptimizeSequence_NonDnaCharacters_DoNotCrashAndPreserveProtein(string input)
    {
        string expectedProtein = Translate(Normalize(input));

        foreach (var strategy in DeterministicStrategies)
        {
            CodonOptimizer.OptimizationResult result = default;
            var act = () => result = CodonOptimizer.OptimizeSequence(input, Target, strategy);

            act.Should().NotThrow(
                $"[{strategy}] a non-standard codon must map to the 'X' sentinel, never KeyNotFound/IndexOutOfRange");

            (result.OptimizedSequence.Length % 3).Should().Be(0,
                $"[{strategy}] the output is still whole codons (INV-02)");
            result.OriginalSequence.Length.Should().Be(result.OptimizedSequence.Length,
                $"[{strategy}] length is preserved codon-for-codon even with un-translatable codons (INV-03)");

            // INV-01 across the whole sequence: un-translatable codons pass through
            // unchanged (X), translatable ones stay synonymous → same protein.
            Translate(result.OptimizedSequence).Should().Be(expectedProtein,
                $"[{strategy}] translate(optimized) == translate(input); 'X' codons are carried through verbatim");
        }
    }

    /// <summary>
    /// MC: an entirely non-coding sequence — length is a multiple of 3 but no codon
    /// is a real codon (all 'X'). The optimizer must produce the input unchanged
    /// (every codon is left as-is), with a protein of all 'X', no crash, and the
    /// length invariants intact. Pins that "no real codons" is a no-op, not a
    /// KeyNotFound on a missing amino-acid entry.
    /// </summary>
    [Test]
    public void OptimizeSequence_FullyNonCoding_IsNoOpAllX()
    {
        const string input = "ZZZQQQJJJ"; // 3 codons, none in the genetic code

        var result = CodonOptimizer.OptimizeSequence(input, Target,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        result.ProteinSequence.Should().Be("XXX",
            "every non-standard codon translates to the 'X' sentinel");
        result.OptimizedSequence.Should().Be(Normalize(input),
            "no synonymous set exists for 'X', so every codon is carried through unchanged");
        result.ChangedCodons.Should().Be(0, "nothing is optimizable, so no codon changes");
    }

    #endregion

    #region BE/OVF — Boundary: extremely long sequence (no hang, invariants hold)

    /// <summary>
    /// BE/OVF: an extremely long valid coding sequence (300,000 codons) must
    /// optimize in linear time without hanging, without overflow, and with every
    /// invariant intact: protein preserved (INV-01), output a multiple of 3
    /// (INV-02), and length preserved (INV-03). A CancelAfter guards against a
    /// pathological hang. The input is the fixed repeat "AUG" (Met) so the
    /// expected protein is deterministic and known.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void OptimizeSequence_ExtremelyLong_StaysLinearAndPreservesProtein()
    {
        const int codonCount = 300_000;
        // Met (AUG) has no synonym, so this is a maximally cheap, fully-defined input.
        string coding = string.Concat(Enumerable.Repeat("ATG", codonCount));
        string expectedProtein = new string('M', codonCount);

        CodonOptimizer.OptimizationResult result = default;
        var act = () => result = CodonOptimizer.OptimizeSequence(
            coding, Target, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        act.Should().NotThrow("a long valid sequence must not overflow or hang");
        result.OptimizedSequence.Length.Should().Be(codonCount * 3,
            "every codon is preserved at length (INV-02/INV-03)");
        Translate(result.OptimizedSequence).Should().Be(expectedProtein,
            "even at scale the protein is preserved (INV-01)");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CODON-CAI-001 — codon adaptation index : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CODON-CAI-001 — codon adaptation index

    /// <summary>
    /// Fuzz tests for CODON-CAI-001 — the Codon Adaptation Index (CAI).
    /// Checklist: docs/checklists/03_FUZZING.md, row 59.
    /// Fuzz strategy exercised for THIS unit:
    ///   • BE = Boundary Exploitation — empty sequence, a sequence made up only of
    ///          stop codons (no evaluable codons), and an input whose length is not
    ///          a multiple of 3 (a trailing partial codon).
    /// — docs/checklists/03_FUZZING.md §Description (strategy codes).
    ///
    /// ───────────────────────────────────────────────────────────────────────────
    /// The CAI contract under test
    /// ───────────────────────────────────────────────────────────────────────────
    /// CAI measures how strongly a coding sequence favours codons preferred in a
    /// reference organism.[Sharp &amp; Li 1987] For each non-stop codon i encoding
    /// amino acid a, the relative adaptiveness is
    ///     w_i = f_i / max(f_j)   over the synonymous codons j of a
    /// and CAI is the geometric mean over the L evaluated codons, computed in the
    /// numerically-stable logarithmic form
    ///     CAI = (∏ w_i)^(1/L) = exp((1/L) · Σ ln w_i).
    /// — docs/algorithms/Codon_Optimization/CAI_Calculation.md §2.2.
    ///
    /// API entry: CodonOptimizer.CalculateCAI(string codingSequence,
    ///   CodonUsageTable table)
    ///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs lines
    ///   423–450).
    ///
    /// The method is LENIENT by documented design — it never throws on garbage; it
    /// normalises, trims, and scores whatever it is given:
    ///   • null OR empty input → returns 0 via an explicit early return
    ///     (CodonOptimizer.cs lines 425–426). CAI_Calculation.md §3.3, §6.1.
    ///   • DNA input is upper-cased and `T` is replaced by `U` before splitting
    ///     (line 428), so DNA / lowercase round-trip identically.
    ///   • input length NOT divisible by 3 → SplitIntoCodons drops the trailing
    ///     partial codon (loop guard `i + 2 &lt; length`, lines 687–695). This is the
    ///     KEY no-crash boundary: a final 1–2 leftover bases must be IGNORED, never
    ///     cause an IndexOutOfRangeException. CAI_Calculation.md §6.1.
    ///   • a sequence of ONLY stop codons → every codon translates to `*` and is
    ///     SKIPPED (line 440); the evaluated count L is 0, so the method returns 0
    ///     from `count > 0 ? … : 0` (line 449) — it must NOT evaluate exp(logSum/0)
    ///     (a 0/0 NaN) nor take ln(0). CAI_Calculation.md §2.4 INV-03, §6.1 (KEY).
    ///   • a codon that is not in the standard genetic code translates to the
    ///     sentinel `X`; CalculateRelativeAdaptiveness returns NaN for it
    ///     (AminoAcidToCodons lookup miss, lines 454–455) and CalculateCAI skips it
    ///     (line 443) — never a KeyNotFoundException, never a NaN result.
    ///
    /// KEY THEORY INVARIANTS this suite pins directly (CAI_Calculation.md §2.4):
    ///   • INV-01: 0 ≤ CAI ≤ 1 on every input that yields a score.
    ///   • A sequence built only from each amino acid's MOST-used synonym → every
    ///     w_i = 1 → CAI = 1 exactly (the worked optimal example, §7.1).
    ///   • INV-03: stop codons do not affect the result (excluded from L).
    /// One exact CAI is pinned against the documented formula (the §7.1 suboptimal
    /// worked example, ≈0.32) so the geometric-mean computation itself is verified,
    /// not just its range.
    ///
    /// Determinism note: CalculateCAI is a pure function of (sequence, table) with
    /// no randomness; every test uses a FIXED input and the deterministic
    /// EColiK12 reference table, so every assertion is reproducible.
    /// ───────────────────────────────────────────────────────────────────────────

    #region Positive sanity — optimal codons → CAI = 1, suboptimal → known value

    /// <summary>
    /// Positive sanity (KEY): a sequence built ONLY from each amino acid's most-used
    /// E. coli synonym makes every w_i = 1, so CAI = 1 exactly — the worked optimal
    /// example AUG·CUG·CCG·ACC from CAI_Calculation.md §7.1. This is the baseline
    /// proving the boundary targets below are measured against a working happy path,
    /// and it pins the upper bound of INV-01 (0 ≤ CAI ≤ 1) at its extreme. DNA
    /// notation (ATGCTGCCGACC) is used to also exercise the T→U normalisation.
    /// </summary>
    [Test]
    public void CalculateCAI_AllOptimalCodons_IsExactlyOne()
    {
        const string optimalDna = "ATGCTGCCGACC"; // M·L·P·T, each the top E. coli codon

        double cai = CodonOptimizer.CalculateCAI(optimalDna, Target);

        cai.Should().BeApproximately(1.0, 1e-12,
            "every codon is the most-used synonym, so w_i = 1 and the geometric mean is 1 (CAI_Calculation.md §7.1)");
        cai.Should().BeInRange(0.0, 1.0, "CAI is bounded by [0, 1] (INV-01)");
    }

    /// <summary>
    /// Positive sanity: the suboptimal worked example AUG·CUA·CCA·ACU from
    /// CAI_Calculation.md §7.1. Pins the geometric-mean computation against the
    /// documented formula exp((1/L)·Σ ln w_i) with the EColiK12 frequencies:
    ///   w = {1.00, 0.04/0.50, 0.19/0.53, 0.16/0.44} → CAI ≈ 0.3196 (doc rounds to
    /// ≈0.31). Confirms the value lies strictly inside (0, 1) — INV-01 at an
    /// interior point, not just the endpoints.
    /// </summary>
    [Test]
    public void CalculateCAI_SuboptimalCodons_MatchesDocumentedFormula()
    {
        const string suboptimalRna = "AUGCUACCAACU"; // M·L·P·T, weak E. coli codons

        // Reference value computed directly from w_i = f_i / max(f_j) and the
        // geometric mean exp((1/L)·Σ ln w_i) over the EColiK12 table.
        double[] w = { 1.00, 0.04 / 0.50, 0.19 / 0.53, 0.16 / 0.44 };
        double expected = Math.Exp(w.Select(x => Math.Log(x)).Sum() / w.Length);

        double cai = CodonOptimizer.CalculateCAI(suboptimalRna, Target);

        cai.Should().BeApproximately(expected, 1e-12,
            "CAI = exp((1/L)·Σ ln w_i) over the EColiK12 frequencies (CAI_Calculation.md §2.2, §7.1)");
        cai.Should().BeInRange(0.0, 1.0, "CAI is bounded by [0, 1] (INV-01)");
        cai.Should().BeApproximately(0.3196, 1e-3, "matches the §7.1 worked example (≈0.31)");
    }

    #endregion

    #region BE — Boundary: empty string / null (no codons at all)

    /// <summary>
    /// BE: empty string and null are the "no input" boundary. The documented
    /// contract returns 0 via an explicit early return (CodonOptimizer.cs lines
    /// 425–426) — NOT a NullReferenceException, and critically NOT exp(logSum/0)
    /// (a 0/0 NaN) from a geometric mean over zero codons. Pins that "no sequence"
    /// is a defined 0, never NaN and never a crash. CAI_Calculation.md §6.1.
    /// </summary>
    [TestCase(null, TestName = "CalculateCAI_Null_IsZeroNoThrow")]
    [TestCase("", TestName = "CalculateCAI_Empty_IsZeroNoThrow")]
    public void CalculateCAI_NullOrEmpty_ReturnsZero(string? input)
    {
        double cai = double.NaN;
        var act = () => cai = CodonOptimizer.CalculateCAI(input!, Target);

        act.Should().NotThrow("null/empty is a defined no-op early return, not an error");
        cai.Should().Be(0.0, "an empty sequence has no evaluable codons → CAI is 0");
        cai.Should().NotBe(double.NaN, "the empty boundary must never produce a 0/0 NaN");
    }

    #endregion

    #region BE — Boundary: sequence of ONLY stop codons (no evaluable codons)

    /// <summary>
    /// BE (KEY no-crash boundary): a sequence made up entirely of stop codons has
    /// NO evaluable codon — every codon translates to `*` and is skipped (line 440),
    /// so the evaluated count L is 0 and the method returns 0 from the
    /// `count > 0 ? … : 0` guard (line 449). It must NOT compute exp(logSum/0)
    /// (a 0/0 NaN), and must NOT take ln(0). This pins INV-03 (stop codons do not
    /// affect the result) at its extreme: a sequence of nothing but stops scores 0,
    /// not NaN, not a DivideByZero, not a crash. Covers each of the three stop
    /// codons alone and all three together, in both RNA and DNA notation.
    /// CAI_Calculation.md §2.4 INV-03, §6.1.
    /// </summary>
    [TestCase("UAAUAGUGA", TestName = "CalculateCAI_AllThreeStops_Rna_IsZero")]
    [TestCase("TAATAGTGA", TestName = "CalculateCAI_AllThreeStops_Dna_IsZero")]
    [TestCase("UAAUAAUAA", TestName = "CalculateCAI_OnlyUAA_IsZero")]
    [TestCase("UAGUAGUAG", TestName = "CalculateCAI_OnlyUAG_IsZero")]
    [TestCase("UGAUGAUGA", TestName = "CalculateCAI_OnlyUGA_IsZero")]
    public void CalculateCAI_OnlyStopCodons_ReturnsZeroNoNaN(string input)
    {
        double cai = double.NaN;
        var act = () => cai = CodonOptimizer.CalculateCAI(input, Target);

        act.Should().NotThrow(
            "a stops-only sequence leaves zero evaluated codons; the L=0 guard must fire, never ln(0)/div-by-zero");
        cai.Should().Be(0.0,
            "stop codons are excluded (INV-03); with no remaining codons CAI is the defined 0");
        double.IsNaN(cai).Should().BeFalse("the L=0 path must return 0, never exp(logSum/0) = NaN");
    }

    /// <summary>
    /// BE: a stop codon sandwiched between real codons must contribute NOTHING to
    /// the score (INV-03) — the CAI of M·*·L·* must equal the CAI of just M·L. Pins
    /// that the stop-skipping leaves a valid geometric mean over the real codons,
    /// not a NaN from a mis-counted L, and keeps the result in [0, 1].
    /// </summary>
    [Test]
    public void CalculateCAI_StopCodonsInterspersed_DoNotAffectScore()
    {
        const string withStops = "AUGUAACUGUGA"; // M · stop · L · stop
        const string withoutStops = "AUGCUG";      // M · L

        double caiWithStops = CodonOptimizer.CalculateCAI(withStops, Target);
        double caiWithoutStops = CodonOptimizer.CalculateCAI(withoutStops, Target);

        caiWithStops.Should().BeApproximately(caiWithoutStops, 1e-12,
            "stop codons are excluded from L, so they cannot change the score (INV-03)");
        caiWithStops.Should().BeInRange(0.0, 1.0, "CAI stays bounded by [0, 1] (INV-01)");
        double.IsNaN(caiWithStops).Should().BeFalse("a valid score must never be NaN");
    }

    #endregion

    #region BE — Boundary: length NOT divisible by 3 (trailing partial codon)

    /// <summary>
    /// BE (KEY no-crash boundary): when the input length is not a multiple of 3,
    /// SplitIntoCodons drops the trailing partial codon (loop guard `i + 2 &lt; length`,
    /// CodonOptimizer.cs lines 687–695) — the leftover 1–2 bases must NEVER trigger
    /// an IndexOutOfRangeException. The CAI scored over the complete codons that
    /// remain must EQUAL the CAI of the trimmed prefix alone, and stay in [0, 1].
    /// Here AUG·CUG (+ 1 or 2 trailing bases) must score identically to AUG·CUG,
    /// which is the all-optimal CAI = 1. Verified for both a +1 and a +2 remainder.
    /// CAI_Calculation.md §6.1, §2.3 ASM-02.
    /// </summary>
    [TestCase("AUGCUGA", TestName = "CalculateCAI_LenMod3Is1_TrimsTrailingBase")]   // 7 = 2 codons + 1
    [TestCase("AUGCUGAU", TestName = "CalculateCAI_LenMod3Is2_TrimsTrailingTwo")]   // 8 = 2 codons + 2
    public void CalculateCAI_LengthNotDivisibleBy3_TrimsPartialCodon(string input)
    {
        const string completePrefix = "AUGCUG"; // the 2 complete codons (M·L), both optimal
        double expected = CodonOptimizer.CalculateCAI(completePrefix, Target);

        double cai = double.NaN;
        var act = () => cai = CodonOptimizer.CalculateCAI(input, Target);

        act.Should().NotThrow(
            "a trailing partial codon must be trimmed, never cause IndexOutOfRange");
        cai.Should().BeApproximately(expected, 1e-12,
            "the partial codon is ignored, so CAI equals that of the complete-codon prefix (ASM-02)");
        cai.Should().BeInRange(0.0, 1.0, "CAI stays bounded by [0, 1] (INV-01)");
        cai.Should().BeApproximately(1.0, 1e-12,
            "AUG·CUG are the top E. coli codons, so the trimmed score is the optimal 1");
    }

    /// <summary>
    /// BE: the smallest non-trivial partial-codon inputs — lengths 1, 2, 4 and 5 —
    /// where after trimming there are either zero complete codons (lengths 1, 2) or
    /// one complete codon plus a partial (lengths 4, 5). None may crash; the
    /// zero-codon cases must return the defined 0 (not NaN), and the one-codon cases
    /// must return a score in [0, 1]. Pins the trim boundary right at the codon edge.
    /// </summary>
    [TestCase("A", 0.0, TestName = "CalculateCAI_Len1_NoCompleteCodon_IsZero")]
    [TestCase("AU", 0.0, TestName = "CalculateCAI_Len2_NoCompleteCodon_IsZero")]
    [TestCase("AUGA", 1.0, TestName = "CalculateCAI_Len4_OneCodonAUG_IsOne")]    // AUG (M, w=1) + 'A'
    [TestCase("AUGAU", 1.0, TestName = "CalculateCAI_Len5_OneCodonAUG_IsOne")]   // AUG (M, w=1) + 'AU'
    public void CalculateCAI_TinyPartialInputs_TrimAtCodonEdge(string input, double expected)
    {
        double cai = double.NaN;
        var act = () => cai = CodonOptimizer.CalculateCAI(input, Target);

        act.Should().NotThrow("sub-codon trailing bases are trimmed, never indexed out of range");
        cai.Should().BeApproximately(expected, 1e-12,
            "zero complete codons → defined 0; one optimal codon (AUG, Met) → w=1 → CAI 1");
        double.IsNaN(cai).Should().BeFalse("no boundary input may yield NaN");
    }

    #endregion

    #endregion
}
