using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Codon area — codon optimization (CODON-OPT-001), the
/// Codon Adaptation Index (CODON-CAI-001), rare-codon detection (CODON-RARE-001)
/// and the codon-usage table (CODON-USAGE-001).
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

    // ═══════════════════════════════════════════════════════════════════
    //  CODON-RARE-001 — rare codon detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CODON-RARE-001 — rare codon detection

    /// <summary>
    /// Fuzz tests for CODON-RARE-001 — rare codon detection.
    /// Checklist: docs/checklists/03_FUZZING.md, row 60.
    /// Fuzz strategy exercised for THIS unit:
    ///   • BE = Boundary Exploitation — empty sequence, an all-rare sequence, an
    ///          all-common sequence, threshold=0 and threshold=1 (the two extremes
    ///          of the [0,1] frequency cutoff).
    /// — docs/checklists/03_FUZZING.md §Description (strategy codes).
    ///
    /// ───────────────────────────────────────────────────────────────────────────
    /// The rare-codon-detection contract under test
    /// ───────────────────────────────────────────────────────────────────────────
    /// Rare codon detection flags codons whose reference frequency in a target
    /// organism is STRICTLY below a threshold (rare codons slow translation because
    /// they are decoded by low-abundance tRNAs).[Kane 1995; Shu et al. 2006]
    ///   — docs/algorithms/Codon_Optimization/Rare_Codon_Detection.md §1, §2.2.
    ///
    /// API entry: CodonOptimizer.FindRareCodons(string codingSequence,
    ///   CodonUsageTable table, double threshold = 0.15)
    ///   → IEnumerable&lt;(int Position, string Codon, string AminoAcid, double Frequency)&gt;
    ///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs lines
    ///   606–629).
    ///
    /// THRESHOLD SEMANTICS (verified): `threshold` is a per-amino-acid relative
    /// FREQUENCY in [0, 1] — the same value stored in CodonUsageTable.CodonFrequencies
    /// — NOT a percentile and NOT an RSCU. A codon is flagged iff
    ///     frequency(codon) &lt; threshold     (STRICT `<`, never `<=`)
    /// (CodonOptimizer.cs line 623; Rare_Codon_Detection.md §2.2, §2.4 INV-02). The
    /// frequency is read with GetValueOrDefault(codon, 0), so a codon absent from the
    /// table is treated as frequency 0 and is flagged whenever threshold &gt; 0
    /// (Rare_Codon_Detection.md §3.3, §6.1).
    ///
    /// The method is LENIENT by documented design — it never throws on garbage; it
    /// normalises, trims, and screens whatever it is given:
    ///   • null OR empty input → yields NO results via an explicit `yield break`
    ///     (CodonOptimizer.cs lines 614–615); NOT a NullReferenceException.
    ///     Rare_Codon_Detection.md §6.1 (Empty sequence → no results).
    ///   • DNA input is upper-cased and `T` is replaced by `U` before splitting
    ///     (line 617), so DNA / lowercase round-trip to RNA codons identically.
    ///   • input length NOT divisible by 3 → SplitIntoCodons drops the trailing
    ///     partial codon (loop guard `i + 2 &lt; length`, lines 687–695); a leftover
    ///     1–2 bases must be IGNORED, never cause an IndexOutOfRangeException.
    ///   • an unknown / non-standard codon → frequency defaults to 0 (flagged when
    ///     threshold &gt; 0) and translates to the sentinel `X`; never a
    ///     KeyNotFoundException. Rare_Codon_Detection.md §6.1.
    ///
    /// KEY THEORY INVARIANTS this suite pins directly (Rare_Codon_Detection.md §2.4):
    ///   • INV-01: every reported Position is a multiple of 3 (Position = index*3).
    ///   • INV-02: every reported Frequency is STRICTLY &lt; threshold.
    ///   • INV-03: every reported Codon has length 3.
    ///   • SOUNDNESS + COMPLETENESS: the flagged set is EXACTLY the codons whose
    ///     table frequency is &lt; threshold — no false positives, no false negatives.
    ///     The suite recomputes the expected flagged set independently from the
    ///     EColiK12 table and asserts set-equality, not just count.
    ///
    /// THRESHOLD EXTREMES (verified against the EColiK12 table):
    ///   • threshold = 0 → no frequency can be &lt; 0 (frequencies are ≥ 0, and even
    ///     unknown codons default to exactly 0, which is NOT &lt; 0) → NONE flagged.
    ///   • threshold = 1 → every codon with frequency &lt; 1 is flagged. In EColiK12
    ///     only AUG (Met) and UGG (Trp) have frequency EXACTLY 1.00, so by the strict
    ///     `<` they are NEVER flagged even at threshold 1; every other codon IS
    ///     flagged. This is the documented strict-comparison boundary (§2.4 INV-02,
    ///     §5.2): "a codon exactly at the threshold is not reported."
    ///
    /// Determinism note: FindRareCodons is a pure function of (sequence, table,
    /// threshold) with no randomness (Rare_Codon_Detection.md §2.4 INV-04); every
    /// test uses a FIXED input and the deterministic EColiK12 reference table, so
    /// every assertion is reproducible.
    /// ───────────────────────────────────────────────────────────────────────────

    #region Helpers (RARE)

    /// <summary>
    /// Independently computes the rare-codon contract from the EColiK12 table: for
    /// each complete codon of the NORMALIZED RNA sequence, flag it iff its table
    /// frequency (default 0 when absent) is STRICTLY less than the threshold. Returns
    /// the expected (position, codon, frequency) flagged set, so a test can assert
    /// FindRareCodons matches it exactly (soundness + completeness) without trusting
    /// the method under test.
    /// </summary>
    private static List<(int Position, string Codon, double Frequency)> ExpectedRare(
        string coding, CodonOptimizer.CodonUsageTable table, double threshold)
    {
        string rna = Normalize(coding);
        var expected = new List<(int, string, double)>();
        for (int i = 0; i + 2 < rna.Length; i += 3)
        {
            string codon = rna.Substring(i, 3);
            double freq = table.CodonFrequencies.GetValueOrDefault(codon, 0);
            if (freq < threshold)
                expected.Add((i, codon, freq));
        }
        return expected;
    }

    #endregion

    #region Positive sanity — a known rare codon is flagged, a common one is not

    /// <summary>
    /// Positive sanity (KEY baseline): the worked example AUG·AGA·AGG·CGA from
    /// Rare_Codon_Detection.md §7.1. At threshold 0.10 the three Arg codons AGA
    /// (0.04), AGG (0.02), CGA (0.06) are below threshold and AUG (Met, 1.00) is
    /// not, so EXACTLY the codons at positions 3, 6, 9 are flagged — and AUG at
    /// position 0 is NOT. This proves the boundary targets below are measured
    /// against a working happy path that both FLAGS a rare codon at a known
    /// position and SPARES a common one, with the documented tuple fields intact.
    /// </summary>
    [Test]
    public void FindRareCodons_KnownRareAndCommon_FlagsExactlyTheRare()
    {
        const string sequence = "AUGAGAAGGCGA"; // M · R(AGA) · R(AGG) · R(CGA)

        var rare = CodonOptimizer.FindRareCodons(sequence, Target, 0.10).ToList();

        // EXACTLY the three rare Arg codons, at the documented positions, with the
        // documented codon / amino-acid / frequency fields (§7.1).
        rare.Should().HaveCount(3, "AGA, AGG and CGA are below 0.10; AUG (1.00) is not");
        rare.Should().BeEquivalentTo(new[]
        {
            (Position: 3, Codon: "AGA", AminoAcid: "R", Frequency: 0.04),
            (Position: 6, Codon: "AGG", AminoAcid: "R", Frequency: 0.02),
            (Position: 9, Codon: "CGA", AminoAcid: "R", Frequency: 0.06),
        }, "the flagged set is exactly the codons with frequency < 0.10 (Rare_Codon_Detection.md §7.1)");

        // The common start codon AUG (freq 1.00) at position 0 is NOT flagged.
        rare.Should().NotContain(r => r.Codon == "AUG",
            "AUG has frequency 1.00, far above 0.10, so it is a common codon");

        // INV-01/02/03 on every reported item.
        rare.Should().OnlyContain(r => r.Position % 3 == 0, "INV-01: positions are multiples of 3");
        rare.Should().OnlyContain(r => r.Frequency < 0.10, "INV-02: every flagged frequency is strictly < threshold");
        rare.Should().OnlyContain(r => r.Codon.Length == 3, "INV-03: every flagged codon is a triplet");
    }

    #endregion

    #region BE — Boundary: empty / null sequence (no codons → no rare codons)

    /// <summary>
    /// BE: empty string and null are the "no input" boundary. The documented
    /// contract yields NO results via an explicit `yield break` (CodonOptimizer.cs
    /// lines 614–615) — NOT a NullReferenceException. Theory: no codons → no rare
    /// codons. Exercised across the default and both extreme thresholds so the empty
    /// boundary is a no-op regardless of cutoff. Rare_Codon_Detection.md §6.1.
    /// </summary>
    [TestCase(null, TestName = "FindRareCodons_Null_IsEmptyNoThrow")]
    [TestCase("", TestName = "FindRareCodons_Empty_IsEmptyNoThrow")]
    public void FindRareCodons_NullOrEmpty_YieldsNothing(string? input)
    {
        foreach (double threshold in new[] { 0.0, 0.15, 1.0 })
        {
            List<(int, string, string, double)> rare = null!;
            var act = () => rare = CodonOptimizer.FindRareCodons(input!, Target, threshold).ToList();

            act.Should().NotThrow(
                $"[threshold={threshold}] null/empty is a defined no-op (yield break), not an error");
            rare.Should().BeEmpty(
                $"[threshold={threshold}] no codons exist, so no codon can be flagged rare");
        }
    }

    #endregion

    #region BE — Boundary: all-rare sequence (every codon flagged)

    /// <summary>
    /// BE (all rare): a sequence built ENTIRELY from below-threshold codons must
    /// have EVERY codon flagged (soundness + completeness at the upper extreme of
    /// the rare set). Uses the three rarest E. coli Arg codons AGG (0.02), AGA
    /// (0.04), CGA (0.06), all &lt; 0.10. At threshold 0.10 all three of the three
    /// codons are reported, at the documented positions, and the flagged set equals
    /// the independently-recomputed expected set exactly. DNA notation also
    /// exercises the T→U normalisation.
    /// </summary>
    [Test]
    public void FindRareCodons_AllRare_FlagsEveryCodon()
    {
        const string allRareDna = "AGGAGACGA"; // R(0.02) · R(0.04) · R(0.06), all < 0.10
        const double threshold = 0.10;

        var rare = CodonOptimizer.FindRareCodons(allRareDna, Target, threshold).ToList();
        var expected = ExpectedRare(allRareDna, Target, threshold);

        rare.Should().HaveCount(3, "all three codons are below 0.10, so all are flagged");
        rare.Select(r => (r.Position, r.Codon, r.Frequency))
            .Should().BeEquivalentTo(expected,
                "the flagged set equals the independently-computed below-threshold set (completeness)");
        rare.Should().OnlyContain(r => r.Frequency < threshold, "INV-02: every flagged frequency is strictly < threshold");
        rare.Should().OnlyContain(r => r.Position % 3 == 0 && r.Codon.Length == 3, "INV-01/INV-03 hold for every item");
    }

    #endregion

    #region BE — Boundary: none-rare sequence (no codon flagged)

    /// <summary>
    /// BE (none rare): a sequence built ENTIRELY from above-threshold (common)
    /// codons must flag NOTHING (no false positives — soundness at the lower
    /// extreme). Uses each amino acid's most-used E. coli synonym AUG (1.00), CUG
    /// (0.50), CCG (0.53), ACC (0.44), all well above the default 0.15. The result
    /// must be empty, with no crash, at the default threshold.
    /// </summary>
    [Test]
    public void FindRareCodons_NoneRare_FlagsNothing()
    {
        const string allCommonDna = "ATGCTGCCGACC"; // M(1.00) · L(0.50) · P(0.53) · T(0.44)

        var rare = CodonOptimizer.FindRareCodons(allCommonDna, Target, 0.15).ToList();

        rare.Should().BeEmpty(
            "every codon's frequency exceeds the 0.15 threshold, so none is rare (no false positives)");
    }

    #endregion

    #region BE — Boundary: threshold = 0 (nothing is below zero → none flagged)

    /// <summary>
    /// BE (threshold = 0): the lower extreme of the [0,1] cutoff. No frequency can
    /// be strictly &lt; 0 — frequencies are non-negative, and even an unknown codon
    /// defaults to EXACTLY 0, which is not &lt; 0 — so NONE is flagged, regardless of
    /// how rare the input is. This pins the documented strict-comparison boundary
    /// (Rare_Codon_Detection.md §3.3, §6.1: unknown codons are flagged only when
    /// threshold &gt; 0). Exercised on the all-rare sequence AND a non-coding (unknown,
    /// freq-0) sequence to prove even freq-0 codons escape at threshold 0.
    /// </summary>
    [TestCase("AGGAGACGA", TestName = "FindRareCodons_Threshold0_AllRareSeq_FlagsNothing")]
    [TestCase("ZZZQQQJJJ", TestName = "FindRareCodons_Threshold0_UnknownFreq0Codons_FlagsNothing")]
    public void FindRareCodons_ThresholdZero_FlagsNothing(string input)
    {
        List<(int, string, string, double)> rare = null!;
        var act = () => rare = CodonOptimizer.FindRareCodons(input, Target, 0.0).ToList();

        act.Should().NotThrow("threshold 0 is a valid boundary, not an error");
        rare.Should().BeEmpty(
            "no frequency is strictly < 0, and freq-0 codons are NOT < 0, so nothing is flagged at threshold 0");
    }

    #endregion

    #region BE — Boundary: threshold = 1 (every sub-1 codon flagged; freq-1 spared)

    /// <summary>
    /// BE (threshold = 1): the upper extreme of the [0,1] cutoff. By the strict
    /// `<`, every codon with frequency &lt; 1 is flagged, while AUG (Met) and UGG (Trp)
    /// — the only EColiK12 codons with frequency EXACTLY 1.00 — are NOT flagged even
    /// here (Rare_Codon_Detection.md §2.4 INV-02, §5.2: "a codon exactly at the
    /// threshold is not reported"). Verified against an independently-recomputed
    /// expected set on a mixed sequence containing both a freq-1 codon (AUG) and
    /// several sub-1 codons, so completeness AND the strict-boundary exclusion are
    /// pinned together.
    /// </summary>
    [Test]
    public void FindRareCodons_ThresholdOne_FlagsAllExceptFrequencyOneCodons()
    {
        // M(AUG,1.00) · W(UGG,1.00) · L(CUG,0.50) · R(AGA,0.04) · A(GCC,0.27)
        const string mixed = "AUGUGGCUGAGAGCC";
        const double threshold = 1.0;

        var rare = CodonOptimizer.FindRareCodons(mixed, Target, threshold).ToList();
        var expected = ExpectedRare(mixed, Target, threshold);

        // Completeness + soundness: exactly the sub-1 codons, none of the freq-1 ones.
        rare.Select(r => (r.Position, r.Codon, r.Frequency))
            .Should().BeEquivalentTo(expected,
                "at threshold 1 the flagged set is exactly the codons with frequency < 1");
        rare.Should().NotContain(r => r.Codon == "AUG",
            "AUG has frequency EXACTLY 1.00, so the strict `<` spares it even at threshold 1 (INV-02)");
        rare.Should().NotContain(r => r.Codon == "UGG",
            "UGG has frequency EXACTLY 1.00, so it too is never flagged (strict `<`)");
        rare.Select(r => r.Codon).Should().Contain(new[] { "CUG", "AGA", "GCC" },
            "every sub-1 codon (CUG, AGA, GCC) is below threshold 1 and is flagged");
        rare.Should().OnlyContain(r => r.Frequency < threshold, "INV-02: every flagged frequency is strictly < 1");
    }

    /// <summary>
    /// BE (threshold = 1, every codon flagged): when NO codon has frequency 1, the
    /// threshold-1 boundary flags EVERY codon. CUG·CCG·ACC·GCC are all sub-1 common
    /// codons, so at threshold 1 all four are reported — the all-flagged extreme,
    /// confirming the strict `<` admits everything below 1. Asserted against the
    /// independently-recomputed expected set.
    /// </summary>
    [Test]
    public void FindRareCodons_ThresholdOne_NoFreqOneCodon_FlagsAll()
    {
        const string noFreqOne = "CUGCCGACCGCC"; // L(0.50) · P(0.53) · T(0.44) · A(0.27), all < 1
        const double threshold = 1.0;

        var rare = CodonOptimizer.FindRareCodons(noFreqOne, Target, threshold).ToList();
        var expected = ExpectedRare(noFreqOne, Target, threshold);

        rare.Should().HaveCount(4, "no codon is at frequency 1, so all four sub-1 codons are flagged");
        rare.Select(r => (r.Position, r.Codon, r.Frequency))
            .Should().BeEquivalentTo(expected,
                "every codon with frequency < 1 is flagged at threshold 1 (completeness)");
    }

    #endregion

    #region BE — Boundary: length not divisible by 3 (trailing partial codon)

    /// <summary>
    /// BE (no-crash boundary): a length not a multiple of 3 must trim the trailing
    /// partial codon (loop guard `i + 2 &lt; length`) — the leftover 1–2 bases must
    /// NEVER trigger an IndexOutOfRangeException. The flagged set over AUG·AGA (+ a
    /// trailing base or two) must equal the flagged set of just AUG·AGA: at
    /// threshold 0.10 only AGA (0.04) at position 3 is rare. Verified for both a +1
    /// and a +2 remainder. Rare_Codon_Detection.md §6.1.
    /// </summary>
    [TestCase("AUGAGAA", TestName = "FindRareCodons_LenMod3Is1_TrimsTrailingBase")]   // 7 = 2 codons + 1
    [TestCase("AUGAGAAU", TestName = "FindRareCodons_LenMod3Is2_TrimsTrailingTwo")]   // 8 = 2 codons + 2
    public void FindRareCodons_LengthNotDivisibleBy3_TrimsPartialCodon(string input)
    {
        const double threshold = 0.10;

        List<(int Position, string Codon, string AminoAcid, double Frequency)> rare = null!;
        var act = () => rare = CodonOptimizer.FindRareCodons(input, Target, threshold).ToList();

        act.Should().NotThrow("a trailing partial codon must be trimmed, never cause IndexOutOfRange");

        // Only the complete codons AUG·AGA are screened; AGA (0.04) is the lone rare one.
        rare.Should().ContainSingle("only AGA is below 0.10 among the complete codons AUG·AGA")
            .Which.Should().BeEquivalentTo(
                (Position: 3, Codon: "AGA", AminoAcid: "R", Frequency: 0.04),
                "the partial codon is ignored, so the flagged set is exactly that of the trimmed prefix");
        rare.Should().OnlyContain(r => r.Position % 3 == 0, "INV-01: positions stay multiples of 3");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CODON-USAGE-001 — codon usage table : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CODON-USAGE-001 — codon usage

    /// <summary>
    /// Fuzz tests for CODON-USAGE-001 — the codon-usage table (raw codon counts).
    /// Checklist: docs/checklists/03_FUZZING.md, row 61.
    /// Fuzz strategy exercised for THIS unit:
    ///   • BE = Boundary Exploitation — empty sequence, a sequence of length 1, a
    ///          sequence of length 2 (the two partial-codon boundaries where no
    ///          complete codon exists), and an extremely long sequence.
    /// — docs/checklists/03_FUZZING.md §Description (strategy codes).
    ///
    /// ───────────────────────────────────────────────────────────────────────────
    /// The codon-usage contract under test
    /// ───────────────────────────────────────────────────────────────────────────
    /// Codon usage tallies how often each codon appears in a coding sequence by
    /// splitting it into non-overlapping in-frame triplets and counting each codon:
    ///     Count(c) = |{ i : seq[3i : 3i+3] == c }|.[Plotkin &amp; Kudla 2011]
    /// — docs/algorithms/Codon_Optimization/Codon_Usage_Analysis.md §2.2.
    ///
    /// API entry: CodonOptimizer.CalculateCodonUsage(string codingSequence)
    ///   → Dictionary&lt;string,int&gt; of raw codon counts
    ///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs lines
    ///   634–652).
    ///
    /// The method is LENIENT by documented design — it never throws on garbage; it
    /// normalises, splits, and counts whatever it is given:
    ///   • null OR empty input → returns an EMPTY dictionary via an explicit early
    ///     return (CodonOptimizer.cs lines 638–639) — NOT a NullReferenceException,
    ///     and critically NOT a DivideByZero: CalculateCodonUsage only COUNTS, it
    ///     never normalises to frequencies, so there is no division at all on the
    ///     count path. Codon_Usage_Analysis.md §6.1 (Empty → {}).
    ///   • DNA input is upper-cased and `T` is replaced by `U` before splitting
    ///     (line 641), so DNA / lowercase round-trip to RNA codons identically.
    ///   • a sequence shorter than one codon (length 1 or 2) → SplitIntoCodons
    ///     yields ZERO codons (loop guard `i + 2 &lt; length`, lines 690–693) → an
    ///     EMPTY dictionary, with NO IndexOutOfRangeException from indexing a partial
    ///     codon. This is the KEY partial-codon boundary. Codon_Usage_Analysis.md
    ///     §6.1 (Incomplete trailing bases → Ignored), ASM-01.
    ///   • input length NOT a multiple of 3 → the trailing 1–2 bases are likewise
    ///     IGNORED; only complete triplets are counted.
    ///
    /// KEY THEORY INVARIANT this suite pins directly (Codon_Usage_Analysis.md §2.4
    /// INV-01; tests/TestSpecs/CODON-USAGE-001.md S3):
    ///   • INV-01: sum(counts.Values) == floor(normalizedLength / 3) — the counter
    ///     increments EXACTLY once per extracted complete codon, so the total tally
    ///     equals the number of complete codons. At the boundaries this means
    ///     empty / len-1 / len-2 all give a total of 0 (an empty table), and a long
    ///     sequence of K codons gives a total of exactly K. Every codon key is a
    ///     length-3 triplet and every count is ≥ 1 (absent codons are NOT stored as
    ///     zero entries — §5.2).
    ///
    /// Determinism note: CalculateCodonUsage is a pure function of the sequence with
    /// no randomness; every test uses a FIXED input, so every assertion is
    /// reproducible.
    /// ───────────────────────────────────────────────────────────────────────────

    #region Positive sanity — a known coding sequence yields the expected counts

    /// <summary>
    /// Positive sanity (KEY baseline): a well-formed coding sequence with known,
    /// repeated codons yields the exact expected per-codon counts, and the total
    /// equals floor(len/3) (INV-01). AUG·AUG·CUG·CUG·CUG·UAA = 6 codons → AUG×2,
    /// CUG×3, UAA×1. DNA notation (ATG…) also exercises the T→U normalisation, so a
    /// 'U'-keyed table is asserted from a 'T'-spelled input. This proves the boundary
    /// targets below are measured against a working happy path, not a uniformly
    /// broken counter.
    /// </summary>
    [Test]
    public void CalculateCodonUsage_KnownCoding_CountsExactlyAndTotalsFloorLenOver3()
    {
        const string codingDna = "ATGATGCTGCTGCTGTAA"; // 18 bases = 6 codons

        var usage = CodonOptimizer.CalculateCodonUsage(codingDna);

        usage.Should().BeEquivalentTo(new Dictionary<string, int>
        {
            ["AUG"] = 2, // ATG → AUG, twice
            ["CUG"] = 3, // CTG → CUG, three times
            ["UAA"] = 1, // TAA → UAA, once (the stop codon is still counted)
        }, "each complete triplet is counted once, with T→U normalisation");

        // INV-01: the tally totals the number of complete codons = floor(len/3).
        usage.Values.Sum().Should().Be(codingDna.Length / 3,
            "sum(counts) == floor(len/3): the counter increments exactly once per codon (INV-01)");
        usage.Values.Sum().Should().Be(6);
        usage.Keys.Should().OnlyContain(c => c.Length == 3, "every codon key is a triplet");
        usage.Values.Should().OnlyContain(v => v >= 1, "absent codons are not stored as zero entries");
    }

    #endregion

    #region BE — Boundary: empty / null sequence (no codons → empty table, no div-by-zero)

    /// <summary>
    /// BE: empty string and null are the "no input" boundary. The documented contract
    /// returns an EMPTY dictionary via an explicit early return (CodonOptimizer.cs
    /// lines 638–639) — NOT a NullReferenceException, and critically NOT a
    /// DivideByZeroException: the count path performs no normalisation, so there is no
    /// division over zero codons. Pins that "no sequence" is a defined empty table
    /// with total 0 (INV-01 at len 0), never a crash. Codon_Usage_Analysis.md §6.1.
    /// </summary>
    [TestCase(null, TestName = "CalculateCodonUsage_Null_IsEmptyNoThrow")]
    [TestCase("", TestName = "CalculateCodonUsage_Empty_IsEmptyNoThrow")]
    public void CalculateCodonUsage_NullOrEmpty_ReturnsEmptyTable(string? input)
    {
        Dictionary<string, int> usage = null!;
        var act = () => usage = CodonOptimizer.CalculateCodonUsage(input!);

        act.Should().NotThrow(
            "null/empty is a defined no-op early return — never a NullReference, never a divide-by-zero");
        usage.Should().BeEmpty("no codons exist, so the usage table is empty");
        usage.Values.Sum().Should().Be(0, "INV-01: floor(0/3) == 0 — an empty table totals 0");
    }

    #endregion

    #region BE — Boundary: length 1 and length 2 (partial codon, no complete codon)

    /// <summary>
    /// BE (KEY partial-codon boundary): a sequence of length 1 or length 2 has NO
    /// complete codon. SplitIntoCodons yields zero codons (loop guard
    /// `i + 2 &lt; length`, CodonOptimizer.cs lines 690–693), so the usage table is
    /// EMPTY — and the trailing 1–2 leftover bases must NEVER trigger an
    /// IndexOutOfRangeException from indexing a partial triplet. Covers length 1 and
    /// length 2, in valid DNA, lowercase, and non-DNA characters, so the trim holds
    /// regardless of alphabet. INV-01: floor(1/3) == floor(2/3) == 0.
    /// </summary>
    [TestCase("A", TestName = "CalculateCodonUsage_Len1_DnaBase_IsEmpty")]
    [TestCase("AT", TestName = "CalculateCodonUsage_Len2_DnaBases_IsEmpty")]
    [TestCase("g", TestName = "CalculateCodonUsage_Len1_Lowercase_IsEmpty")]
    [TestCase("at", TestName = "CalculateCodonUsage_Len2_Lowercase_IsEmpty")]
    [TestCase("Z", TestName = "CalculateCodonUsage_Len1_NonDna_IsEmpty")]
    [TestCase("ZZ", TestName = "CalculateCodonUsage_Len2_NonDna_IsEmpty")]
    public void CalculateCodonUsage_LengthOneOrTwo_NoCompleteCodon_IsEmpty(string input)
    {
        Dictionary<string, int> usage = null!;
        var act = () => usage = CodonOptimizer.CalculateCodonUsage(input);

        act.Should().NotThrow(
            "a sub-codon input must yield zero codons, never index a partial triplet out of range");
        usage.Should().BeEmpty("no complete codon exists at length 1 or 2, so nothing is counted");
        usage.Values.Sum().Should().Be(input.Length / 3,
            "INV-01: floor(len/3) == 0 for len 1 and 2 — the table totals 0");
    }

    /// <summary>
    /// BE: a sequence whose length is NOT a multiple of 3 must IGNORE the trailing
    /// partial codon (the +1/+2 leftover bases) — never index it out of range. The
    /// counts over the complete prefix must be exactly those of the trimmed prefix,
    /// and the total must equal floor(len/3) (INV-01). Covers a +1 and a +2 remainder
    /// just past the first complete codon (lengths 4 and 5 → exactly one codon).
    /// </summary>
    [TestCase("AUGA", TestName = "CalculateCodonUsage_Len4_OneCodonPlusOne_CountsOne")]   // AUG + 'A'
    [TestCase("AUGAU", TestName = "CalculateCodonUsage_Len5_OneCodonPlusTwo_CountsOne")]  // AUG + 'AU'
    public void CalculateCodonUsage_LengthNotMultipleOf3_TrimsTrailingPartialCodon(string input)
    {
        Dictionary<string, int> usage = null!;
        var act = () => usage = CodonOptimizer.CalculateCodonUsage(input);

        act.Should().NotThrow("a trailing partial codon must be ignored, never cause IndexOutOfRange");
        usage.Should().BeEquivalentTo(new Dictionary<string, int> { ["AUG"] = 1 },
            "only the single complete codon AUG is counted; the trailing 1–2 bases are dropped");
        usage.Values.Sum().Should().Be(input.Length / 3,
            "INV-01: sum(counts) == floor(len/3) — the partial codon contributes nothing");
    }

    #endregion

    #region BE/OVF — Boundary: extremely long sequence (no hang, INV-01 holds at scale)

    /// <summary>
    /// BE/OVF: an extremely long valid coding sequence (300,000 codons) must tally in
    /// linear time without hanging, without overflow, and with INV-01 intact at scale:
    /// the total count equals exactly the number of complete codons, floor(len/3).
    /// A CancelAfter guards against a pathological hang. The input is the fixed repeat
    /// "ATG" (→ AUG, Met) so the expected table is a single key with a known count,
    /// making the assertion deterministic and the total exact.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void CalculateCodonUsage_ExtremelyLong_StaysLinearAndTotalsFloorLenOver3()
    {
        const int codonCount = 300_000;
        string coding = string.Concat(Enumerable.Repeat("ATG", codonCount)); // all AUG (Met)

        Dictionary<string, int> usage = null!;
        var act = () => usage = CodonOptimizer.CalculateCodonUsage(coding);

        act.Should().NotThrow("a long valid sequence must not overflow or hang");
        usage.Should().ContainKey("AUG").WhoseValue.Should().Be(codonCount,
            "every one of the 300,000 ATG triplets normalises to AUG and is counted");
        usage.Should().HaveCount(1, "the sequence is a single repeated codon, so one key");
        usage.Values.Sum().Should().Be(coding.Length / 3,
            "INV-01 holds at scale: sum(counts) == floor(len/3) == 300,000");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CODON-ENC-001 — effective number of codons (Wright Nc) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CODON-ENC-001 — effective number of codons (ENC / Nc)

    /// <summary>
    /// Fuzz tests for CODON-ENC-001 — the Effective Number of Codons (ENC / Nc),
    /// Wright's 1990 estimator of synonymous codon-usage bias in a single gene.
    /// Checklist: docs/checklists/03_FUZZING.md, row 213.
    /// Fuzz strategy exercised for THIS unit:
    ///   • BE = Boundary Exploitation — the degenerate extremes of the [20, 61]
    ///          range: a SINGLE codon (no estimable degeneracy class), UNIFORM
    ///          synonymous usage (no bias → the upper re-adjustment to 61), an input
    ///          whose LENGTH is NOT a multiple of 3 (trailing partial codon), plus
    ///          empty / null and non-ACGT garbage.
    /// — docs/checklists/03_FUZZING.md §Description (strategy codes);
    ///   docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
    ///
    /// ───────────────────────────────────────────────────────────────────────────
    /// The ENC / Nc contract under test (Effective_Number_of_Codons.md)
    /// ───────────────────────────────────────────────────────────────────────────
    /// Nc answers "how many codons are effectively in use in this gene?" and ranges
    /// from 20 (extreme bias — exactly one codon per amino acid) to 61 (no bias —
    /// every synonymous codon used equally).[Wright 1990; Fuglsang 2004]
    ///   — Effective_Number_of_Codons.md §1, §2.4 INV-01.
    ///
    /// For an amino acid with k synonymous codons, total count n and frequencies
    /// p_i = n_i/n, the codon homozygosity is (Wright Eq. 1):
    ///     F̂ = ( n·Σ p_i² − 1 ) / ( n − 1 )
    /// and the gene-level value aggregates class averages F̂₂, F̂₃, F̂₄, F̂₆ over the
    /// standard-code degeneracy classes (Wright Eq. 3):
    ///     N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆
    /// where 2 is the Met+Trp single-codon contribution and 9/1/5/3 are the numbers
    /// of two-/three-/four-/six-fold degenerate amino acids in the standard genetic
    /// code. Stop codons are excluded.
    ///   — Effective_Number_of_Codons.md §2.2 (Eq. 1, Eq. 3), §4.2.
    ///
    /// API entry: CodonUsageAnalyzer.CalculateEnc(string)  → double in [20, 61]
    ///            CodonUsageAnalyzer.CalculateEnc(DnaSequence)
    ///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs,
    ///    CalculateEnc / CalculateEncCore).
    ///
    /// The string overload is LENIENT by documented design — it never throws on
    /// garbage; it upper-cases, reads codons in non-overlapping in-frame triplets
    /// from index 0, and scores whatever it is given:
    ///   • null OR empty string → returns 0 via an explicit early return
    ///     (CalculateEnc(string) `string.IsNullOrEmpty` guard); NOT a crash and NOT
    ///     a value in [20, 61]. Effective_Number_of_Codons.md §3.2, §3.3, §6.1.
    ///   • input length NOT divisible by 3 → CountCodonsCore reads only complete
    ///     triplets (loop guard `i + 3 &lt;= seq.Length`); the trailing 1–2 leftover
    ///     bases are IGNORED, never an IndexOutOfRangeException. §3.3, §6.1.
    ///   • a codon containing any non-ACGT character → IsValidCodon is false and the
    ///     codon is SKIPPED (consistent with CountCodons), never a KeyNotFound. §6.1.
    ///   • an amino acid with total count n ≤ 1 → skipped (F̂ undefined, denominator
    ///     n − 1); its degeneracy class falls back to the within-class average
    ///     (Eq. 4) or, if no class member is estimable, to the class's full codon
    ///     count. §2.3 ASM-02, §4.1.
    /// The DnaSequence overload throws ArgumentNullException for null (§3.3, §6.1).
    ///
    /// KEY THEORY INVARIANTS this suite pins directly (Effective_Number_of_Codons.md
    /// §2.4):
    ///   • INV-01: 20 ≤ Nc ≤ 61 on every scored input (a randomized boundary sweep
    ///     asserts this on hundreds of random sequences).
    ///   • INV-02: one codon per amino acid, each used ≥2× ⇒ Nc = 20 (maximum bias —
    ///     every F̂ = 1, sum = 9+1+5+3+2 = 20).
    ///   • INV-03: uniform synonymous usage (each synonym used ≥2×) ⇒ Nc re-adjusted
    ///     to exactly 61 (Wright's overshoot cap — the maximally-unbiased extreme).
    ///   • INV-04: deterministic — a pure function of the codon counts.
    /// Two exact, hand-checkable values are pinned against the documented Wright
    /// formula: the §7.1 worked example (Phe TTT×3, TTC×1 ⇒ Nc = 29.0) and a single
    /// 2-fold family at uniform usage (TTT×2, TTC×2 ⇒ Nc = 38.0), so the F̂ / Eq. 3
    /// computation itself is verified, not merely its range.
    ///
    /// SUBTLETY pinned here (verified independently from Wright Eq. 1, NOT echoed off
    /// the code): "uniform usage" maps to Nc = 61 ONLY when each synonym is used at
    /// least TWICE. When every codon appears EXACTLY ONCE, each 2-fold family has
    /// n = 2, p = (0.5, 0.5), Σp² = 0.5, so F̂ = (2·0.5 − 1)/(2 − 1) = 0; an F̂ of 0
    /// is unestimable (ClassContribution requires f &gt; 0) so EVERY class falls back to
    /// its full codon count and the aggregate collapses to the structural floor 20.
    /// Both the "×1 ⇒ 20" and "×2 ⇒ 61" cases are asserted so the boundary is pinned
    /// on the correct side of the n &gt; 1 requirement, not assumed.
    ///
    /// Determinism note (INV-04): CalculateEnc is a pure function of the sequence
    /// with no randomness. The randomized sweep uses a LOCALLY-seeded `new Random(seed)`
    /// (never a shared static Rng); each generated sequence is fully reproducible and
    /// the same input always yields the same Nc, which the sweep also asserts.
    /// ───────────────────────────────────────────────────────────────────────────

    #region Helpers (ENC)

    /// <summary>Lower / upper bounds of Wright's Nc (Effective_Number_of_Codons.md §2.4 INV-01).</summary>
    private const double EncMin = 20.0;
    private const double EncMax = 61.0;

    /// <summary>
    /// Builds a coding sequence in which EVERY non-stop sense codon of the standard
    /// genetic code appears exactly `copiesPerCodon` times — i.e. perfectly UNIFORM
    /// synonymous usage. With copiesPerCodon ≥ 2 every represented amino acid has
    /// n ≥ 2 per codon, F̂ is small and Nc overshoots → re-adjusted to 61 (INV-03).
    /// The codon order is fixed, so the result is deterministic.
    /// </summary>
    private static string UniformAllCodons(int copiesPerCodon)
    {
        var sense = StandardGeneticCode
            .Where(kv => kv.Value != "*")
            .Select(kv => kv.Key.Replace('U', 'T')) // DNA notation for CalculateEnc
            .OrderBy(c => c, StringComparer.Ordinal)
            .ToList();
        var sb = new System.Text.StringBuilder();
        foreach (var codon in sense)
            for (int i = 0; i < copiesPerCodon; i++)
                sb.Append(codon);
        return sb.ToString();
    }

    #endregion

    #region Positive sanity — exact hand-checkable Wright values from the doc

    /// <summary>
    /// Positive sanity (KEY): the §7.1 worked example. A gene with only Phe
    /// (TTT×3, TTC×1): n = 4, p = (0.75, 0.25), Σp² = 0.625, F̂ = (4·0.625 − 1)/3 = 0.5,
    /// N̂c(Phe) = 1/F̂ = 2. No other degeneracy class is estimable, so each contributes
    /// its full codon count (9, 1, 5, 3): Nc = 2 + 9/0.5 + 1 + 5 + 3 = 29.0. This is
    /// the documented numerical walk-through (Effective_Number_of_Codons.md §7.1),
    /// computed here independently from Wright Eq. 1/Eq. 3 — it both verifies the
    /// formula and proves the fuzz targets below are measured against a working happy
    /// path, not a uniformly-broken method.
    /// </summary>
    [Test]
    public void CalculateEnc_PheOnlyWorkedExample_Is29()
    {
        const string gene = "TTTTTTTTTTTC"; // TTT, TTT, TTT, TTC  (Phe×4, ratio 3:1)

        double nc = CodonUsageAnalyzer.CalculateEnc(gene);

        nc.Should().BeApproximately(29.0, 1e-9,
            "Wright Eq. 1/Eq. 3 on Phe (TTT×3, TTC×1): 2 + 9/0.5 + 1 + 5 + 3 = 29.0 (Effective_Number_of_Codons.md §7.1)");
        nc.Should().BeInRange(EncMin, EncMax, "Nc is bounded by [20, 61] (INV-01)");
    }

    /// <summary>
    /// Positive sanity (interior value): a single 2-fold family at UNIFORM usage with
    /// each synonym used twice — Phe TTT×2, TTC×2. n = 4, p = (0.5, 0.5), Σp² = 0.5,
    /// F̂ = (4·0.5 − 1)/3 = 1/3, so 9/F̂₂ = 27; no other class estimable (full counts
    /// 1, 5, 3): Nc = 2 + 27 + 1 + 5 + 3 = 38.0. An exact interior point of the
    /// [20, 61] range, hand-derived from Wright Eq. 1/Eq. 3 — it pins the F̂
    /// computation away from both clamps, so a wrong homozygosity formula could not
    /// pass by accidentally hitting a boundary.
    /// </summary>
    [Test]
    public void CalculateEnc_SinglePheFamilyUniform_Is38()
    {
        const string gene = "TTTTTCTTTTTC"; // TTT, TTC, TTT, TTC  (Phe×4, ratio 1:1)

        double nc = CodonUsageAnalyzer.CalculateEnc(gene);

        nc.Should().BeApproximately(38.0, 1e-9,
            "Wright Eq. 1/Eq. 3 on Phe (TTT×2, TTC×2): F̂₂ = 1/3 ⇒ 2 + 9/(1/3) + 1 + 5 + 3 = 38.0");
        nc.Should().BeInRange(EncMin, EncMax, "Nc is bounded by [20, 61] (INV-01)");
    }

    #endregion

    #region BE — Boundary: single codon (no estimable class → floor 20)

    /// <summary>
    /// BE (single codon, KEY): one codon is the minimal coding input. ATG (Met) is a
    /// single-codon amino acid (degeneracy 1) and is excluded from the F̂ loop; no
    /// degeneracy class has any estimable F̂, so every class contributes its FULL
    /// codon count and the aggregate is 2 + 9 + 1 + 5 + 3 = 20 — the maximum-bias
    /// floor (INV-01 lower bound, INV-02). A single 2-fold codon (e.g. one TTT) is
    /// also skipped because its amino acid has n = 1 ≤ 1 (F̂ undefined, ASM-02), so it
    /// too yields 20. Neither may crash, and neither may exceed [20, 61]. Covers Met,
    /// a 2-fold sense codon, a 6-fold sense codon and a 4-fold sense codon.
    /// — Effective_Number_of_Codons.md §2.3 ASM-02, §2.4 INV-01/INV-02, §4.1.
    /// </summary>
    [TestCase("ATG", TestName = "CalculateEnc_SingleCodon_Met_Is20")]
    [TestCase("TTT", TestName = "CalculateEnc_SingleCodon_Phe2Fold_Is20")]
    [TestCase("CTG", TestName = "CalculateEnc_SingleCodon_Leu6Fold_Is20")]
    [TestCase("GCC", TestName = "CalculateEnc_SingleCodon_Ala4Fold_Is20")]
    public void CalculateEnc_SingleCodon_IsFloor20(string singleCodon)
    {
        double nc = double.NaN;
        var act = () => nc = CodonUsageAnalyzer.CalculateEnc(singleCodon);

        act.Should().NotThrow("a single in-frame codon is a valid degenerate input, not an error");
        nc.Should().Be(20.0,
            "no degeneracy class is estimable from one codon, so every class contributes its full count: 2+9+1+5+3 = 20 (INV-02)");
        nc.Should().BeInRange(EncMin, EncMax, "Nc is bounded by [20, 61] (INV-01)");
        double.IsNaN(nc).Should().BeFalse("a single codon must never produce NaN");
    }

    /// <summary>
    /// BE (maximum bias, INV-02): one codon per amino acid, each used TWICE, makes
    /// F̂ = 1 for every represented family (Σp² = 1 ⇒ F̂ = (n·1 − 1)/(n − 1) = 1), so
    /// each N̂c(aa) = 1 and the sum is 9 + 1 + 5 + 3 + 2 = 20 — Wright's extreme-bias
    /// limit. Verifies the lower extreme of [20, 61] from a genuinely biased gene
    /// (not just the structural fallback). The sequence repeats the first synonym of
    /// each sense amino acid twice.
    /// </summary>
    [Test]
    public void CalculateEnc_OneCodonPerAminoAcid_IsExtremeBias20()
    {
        // First codon of each sense amino acid, each used twice → every F̂ = 1.
        var firstCodonPerAa = StandardGeneticCode
            .Where(kv => kv.Value != "*")
            .GroupBy(kv => kv.Value)
            .Select(g => g.OrderBy(kv => kv.Key, StringComparer.Ordinal).First().Key.Replace('U', 'T'))
            .OrderBy(c => c, StringComparer.Ordinal);
        string gene = string.Concat(firstCodonPerAa.Select(c => c + c)); // each twice

        double nc = CodonUsageAnalyzer.CalculateEnc(gene);

        nc.Should().BeApproximately(20.0, 1e-9,
            "one codon per amino acid (each ≥2×) ⇒ every F̂ = 1 ⇒ Nc = 9+1+5+3+2 = 20 (INV-02)");
    }

    #endregion

    #region BE — Boundary: uniform synonymous usage (no bias → re-adjusted to 61)

    /// <summary>
    /// BE (uniform usage, KEY — INV-03): perfectly even synonymous usage is the
    /// no-bias extreme. With EVERY sense codon used twice (n ≥ 4 per 2-fold family),
    /// each F̂ is small, the aggregate 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆ overshoots 61,
    /// and Wright's re-adjustment caps it at EXACTLY 61 — the maximally-unbiased
    /// value. Verified at 2, 3 and 5 copies per codon so the cap is robust, not a
    /// fluke of one count. — Effective_Number_of_Codons.md §2.4 INV-03, §6.1
    /// ("near-uniform short gene → re-adjusted to 61").
    /// </summary>
    [TestCase(2, TestName = "CalculateEnc_UniformUsage_TwoCopies_Is61")]
    [TestCase(3, TestName = "CalculateEnc_UniformUsage_ThreeCopies_Is61")]
    [TestCase(5, TestName = "CalculateEnc_UniformUsage_FiveCopies_Is61")]
    public void CalculateEnc_UniformSynonymousUsage_ReadjustedTo61(int copiesPerCodon)
    {
        string gene = UniformAllCodons(copiesPerCodon);

        double nc = CodonUsageAnalyzer.CalculateEnc(gene);

        nc.Should().BeApproximately(61.0, 1e-9,
            "uniform synonymous usage (each codon ≥2×) is maximally unbiased ⇒ Nc re-adjusted to 61 (INV-03)");
        nc.Should().BeInRange(EncMin, EncMax, "Nc never exceeds the [20, 61] range (INV-01)");
    }

    /// <summary>
    /// BE (uniform usage SUBTLETY, pinned independently from Wright Eq. 1): when every
    /// sense codon appears EXACTLY ONCE, each 2-fold family has n = 2, p = (0.5, 0.5),
    /// Σp² = 0.5, so F̂ = (2·0.5 − 1)/(2 − 1) = 0. An F̂ of 0 is unestimable
    /// (a class contribution requires F̂ &gt; 0), so EVERY degeneracy class falls back to
    /// its full codon count and the aggregate collapses to the structural floor 20 —
    /// NOT 61. This pins the boundary on the correct side of the n &gt; 1 / F̂ &gt; 0
    /// requirement (ASM-02): "uniform usage → 61" holds only when each synonym is used
    /// at least twice; with single copies the result is the floor, never a NaN or a
    /// divide-by-zero. The value 20 is hand-derived, not echoed off the code.
    /// </summary>
    [Test]
    public void CalculateEnc_AllCodonsExactlyOnce_CollapsesToFloor20NotNaN()
    {
        string gene = UniformAllCodons(1); // each sense codon once → every F̂ = 0

        double nc = double.NaN;
        var act = () => nc = CodonUsageAnalyzer.CalculateEnc(gene);

        act.Should().NotThrow("F̂ = 0 for every family must not divide by zero or throw");
        nc.Should().Be(20.0,
            "every F̂ = 0 is unestimable ⇒ all classes use their full counts ⇒ floor 20, not 61 (ASM-02)");
        double.IsNaN(nc).Should().BeFalse("an F̂ of 0 must never produce a 0-division NaN");
    }

    #endregion

    #region BE — Boundary: empty / null string (defined 0, never in [20,61])

    /// <summary>
    /// BE: empty string and null are the "no gene" boundary. The string overload
    /// returns 0 via an explicit `string.IsNullOrEmpty` early return — NOT a
    /// NullReferenceException and, critically, NOT a value in [20, 61] (an empty gene
    /// has no Wright rule). Pins the documented degenerate-input contract
    /// (Effective_Number_of_Codons.md §3.2, §3.3, §6.1).
    /// </summary>
    [TestCase(null, TestName = "CalculateEnc_NullString_IsZeroNoThrow")]
    [TestCase("", TestName = "CalculateEnc_EmptyString_IsZeroNoThrow")]
    public void CalculateEnc_NullOrEmptyString_ReturnsZero(string? input)
    {
        double nc = double.NaN;
        var act = () => nc = CodonUsageAnalyzer.CalculateEnc(input!);

        act.Should().NotThrow("null/empty string is a defined no-op early return, not an error");
        nc.Should().Be(0.0, "an empty gene has no codons ⇒ the documented degenerate value 0");
        double.IsNaN(nc).Should().BeFalse("the empty boundary must never produce NaN");
    }

    /// <summary>
    /// BE: the DnaSequence overload's documented contract is to THROW
    /// ArgumentNullException for a null reference (Effective_Number_of_Codons.md §3.3,
    /// §6.1) — a *documented, intentional* validation exception, the disciplined
    /// alternative to a raw NullReferenceException (ADVANCED_TESTING_CHECKLIST.md §8).
    /// </summary>
    [Test]
    public void CalculateEnc_NullDnaSequence_ThrowsArgumentNullException()
    {
        Seqeron.Genomics.Core.DnaSequence? nullSeq = null;
        var act = () => CodonUsageAnalyzer.CalculateEnc(nullSeq!);

        act.Should().Throw<ArgumentNullException>(
            "the DnaSequence overload validates null by contract (§3.3), not a raw NRE");
    }

    #endregion

    #region BE — Boundary: length NOT a multiple of 3 (trailing partial codon)

    /// <summary>
    /// BE (length not %3, KEY no-crash boundary): when the input length is not a
    /// multiple of 3, CountCodonsCore reads only complete in-frame triplets (loop
    /// guard `i + 3 &lt;= length`) — the trailing 1–2 leftover bases are IGNORED, never
    /// an IndexOutOfRangeException. The Nc over the complete codons that remain must
    /// EQUAL the Nc of the trimmed prefix alone. Here the prefix is Phe (TTT×3, TTC×1)
    /// whose documented Nc is 29.0 (§7.1); appending +1 or +2 trailing bases must not
    /// change it. Verified for both a +1 and a +2 remainder.
    /// — Effective_Number_of_Codons.md §3.3, §6.1.
    /// </summary>
    [TestCase("TTTTTTTTTTTCA", TestName = "CalculateEnc_LenMod3Is1_TrimsTrailingBase")]  // 13 = 4 codons + 1
    [TestCase("TTTTTTTTTTTCAT", TestName = "CalculateEnc_LenMod3Is2_TrimsTrailingTwo")] // 14 = 4 codons + 2
    public void CalculateEnc_LengthNotMultipleOf3_TrimsTrailingPartialCodon(string input)
    {
        const string completePrefix = "TTTTTTTTTTTC"; // Phe×4 (3:1) → Nc 29.0

        double nc = double.NaN;
        var act = () => nc = CodonUsageAnalyzer.CalculateEnc(input);

        act.Should().NotThrow(
            "a trailing partial codon must be ignored, never cause IndexOutOfRange");
        nc.Should().BeApproximately(CodonUsageAnalyzer.CalculateEnc(completePrefix), 1e-9,
            "the partial codon is dropped, so Nc equals that of the complete-codon prefix");
        nc.Should().BeApproximately(29.0, 1e-9, "the trimmed prefix is the §7.1 Phe gene ⇒ Nc 29.0");
        nc.Should().BeInRange(EncMin, EncMax, "Nc stays bounded by [20, 61] (INV-01)");
    }

    /// <summary>
    /// BE: a sub-codon input (length 1 or 2) has NO complete codon at all; the loop
    /// never executes, CountCodonsCore returns an empty table, no class is estimable,
    /// and the aggregate is the full-count floor 20 — no IndexOutOfRange on the 1–2
    /// leftover bases. Pins the trim boundary right at the codon edge.
    /// </summary>
    [TestCase("A", TestName = "CalculateEnc_LenOne_NoCompleteCodon_IsFloor20")]
    [TestCase("AT", TestName = "CalculateEnc_LenTwo_NoCompleteCodon_IsFloor20")]
    public void CalculateEnc_SubCodonLength_IsFloor20NoThrow(string input)
    {
        double nc = double.NaN;
        var act = () => nc = CodonUsageAnalyzer.CalculateEnc(input);

        act.Should().NotThrow("sub-codon trailing bases are ignored, never indexed out of range");
        nc.Should().Be(20.0, "no complete codon ⇒ no estimable class ⇒ full-count floor 20");
        nc.Should().BeInRange(EncMin, EncMax, "Nc stays bounded by [20, 61] (INV-01)");
    }

    #endregion

    #region MC/INJ — non-ACGT codons are skipped (never a crash)

    /// <summary>
    /// MC/INJ: codons containing any non-ACGT character are SKIPPED by IsValidCodon
    /// (consistent with CountCodons), never a KeyNotFoundException. Threading garbage
    /// codons (digits, gap, null byte, ambiguity N, unicode) between real Phe codons
    /// must leave the result EQUAL to the Phe-only gene's Nc (29.0) — the garbage
    /// codons contribute nothing — with no crash and the value in [20, 61].
    /// — Effective_Number_of_Codons.md §6.1 (non-ACGT codon → skipped).
    /// </summary>
    [TestCase("TTT123TTTNNNTTTGGGTTC", TestName = "CalculateEnc_NonAcgt_DigitsNAmbig_SkippedPheStays29")]
    [TestCase("TTT---TTT\0\0\0TTTTTC", TestName = "CalculateEnc_NonAcgt_GapNullByte_SkippedPheStays29")]
    [TestCase("TTTαβγTTTTTTαβγTTC", TestName = "CalculateEnc_NonAcgt_Unicode_SkippedPheStays29")]
    public void CalculateEnc_NonAcgtCodons_SkippedAndDoNotCrash(string input)
    {
        // GGG is a real codon (Gly), so the first case keeps one Gly — but Gly has
        // n = 1 there (single copy) so it too is skipped (ASM-02). The net evaluable
        // content is Phe TTT×3, TTC×1 in every case ⇒ Nc 29.0.
        double nc = double.NaN;
        var act = () => nc = CodonUsageAnalyzer.CalculateEnc(input);

        act.Should().NotThrow(
            "non-ACGT codons must be skipped (IsValidCodon false), never a KeyNotFound/IndexOutOfRange");
        nc.Should().BeApproximately(29.0, 1e-9,
            "garbage codons (and the lone n=1 Gly) contribute nothing; the evaluable Phe gene scores 29.0");
        nc.Should().BeInRange(EncMin, EncMax, "Nc stays bounded by [20, 61] (INV-01)");
    }

    #endregion

    #region BE — Randomized boundary sweep (INV-01 always holds; INV-04 determinism)

    /// <summary>
    /// BE (randomized sweep, KEY): hundreds of random sequences — random ACGT content,
    /// random lengths (including sub-codon, not-%3, and longer), plus a fraction
    /// salted with non-ACGT noise — must EVERY time yield a finite Nc that is either
    /// the degenerate 0 (empty input) OR a value in [20, 61] (INV-01), with NO crash,
    /// NO hang, NO NaN and NO Infinity. The Random is LOCALLY seeded (never a shared
    /// static Rng) so the whole sweep is reproducible; each sequence is scored twice
    /// to also pin determinism (INV-04). A CancelAfter guards against a pathological
    /// hang on any single input. This is the core fuzz bar for the unit.
    /// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"; Effective_Number_of_Codons.md §2.4.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void CalculateEnc_RandomSweep_AlwaysFiniteAndInRangeOrZero()
    {
        const int seed = 20213; // local, fixed → reproducible
        var rng = new Random(seed);
        const string acgt = "ACGT";
        // A small noise alphabet to occasionally produce invalid (skipped) codons.
        const string noisy = "ACGTNXacgt-0\0αβ ";

        for (int trial = 0; trial < 600; trial++)
        {
            int len = rng.Next(0, 240); // includes 0, 1, 2, not-%3 and longer
            bool addNoise = rng.Next(0, 4) == 0; // ~25% of inputs carry non-ACGT noise
            string alphabet = addNoise ? noisy : acgt;

            var sb = new System.Text.StringBuilder(len);
            for (int i = 0; i < len; i++)
                sb.Append(alphabet[rng.Next(alphabet.Length)]);
            string seq = sb.ToString();

            double nc = double.NaN;
            var act = () => nc = CodonUsageAnalyzer.CalculateEnc(seq);

            act.Should().NotThrow(
                $"trial {trial} (len {len}, noise {addNoise}) must not crash on random/garbage input");
            double.IsNaN(nc).Should().BeFalse($"trial {trial}: Nc must never be NaN");
            double.IsInfinity(nc).Should().BeFalse($"trial {trial}: Nc must never be Infinity");

            if (nc != 0.0)
                nc.Should().BeInRange(EncMin, EncMax,
                    $"trial {trial}: a scored Nc must lie in [20, 61] (INV-01)");

            // INV-04: deterministic — the same input always yields the same Nc.
            CodonUsageAnalyzer.CalculateEnc(seq).Should().Be(nc,
                $"trial {trial}: CalculateEnc is a pure function (INV-04)");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CODON-RSCU-001 — relative synonymous codon usage (RSCU) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CODON-RSCU-001 — relative synonymous codon usage (RSCU)

    /// <summary>
    /// Fuzz tests for CODON-RSCU-001 — Relative Synonymous Codon Usage (RSCU),
    /// Sharp, Tuohy &amp; Mosurski's 1986 within-family normalization of codon counts.
    /// Checklist: docs/checklists/03_FUZZING.md, row 214.
    /// Fuzz strategy exercised for THIS unit:
    ///   • BE = Boundary Exploitation — the checklist boundaries "single codon,
    ///          missing amino acid, empty": a SINGLE codon (the one synonym takes all
    ///          family usage ⇒ RSCU = n_i), an amino acid MISSING from the input
    ///          (family total 0 ⇒ no divide-by-zero, defined 0), and EMPTY / null
    ///          input; plus an input whose LENGTH is NOT a multiple of 3 (trailing
    ///          partial codon) and non-ACGT garbage.
    /// — docs/checklists/03_FUZZING.md §Description (strategy codes);
    ///   docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
    ///
    /// ───────────────────────────────────────────────────────────────────────────
    /// The RSCU contract under test (Relative_Synonymous_Codon_Usage.md)
    /// ───────────────────────────────────────────────────────────────────────────
    /// RSCU measures, per codon, how often it is used relative to the usage expected
    /// if all synonymous codons of the same amino acid were used equally.[Sharp,
    /// Tuohy &amp; Mosurski 1986] For amino acid i with degeneracy n_i (number of
    /// synonymous codons) and observed counts x_{i,j}, the value is (§2.2, Eq.):
    ///     RSCU_{i,j} = x_{i,j} / ((1/n_i)·Σ_k x_{i,k}) = (n_i · x_{i,j}) / Σ_k x_{i,k}
    /// 1 means no bias, &gt; 1 over-representation, &lt; 1 under-representation.
    ///   — Relative_Synonymous_Codon_Usage.md §1, §2.2.
    ///
    /// API entry: CodonUsageAnalyzer.CalculateRscu(string)  → Dictionary&lt;string,double&gt;
    ///            CodonUsageAnalyzer.CalculateRscu(DnaSequence)
    ///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs,
    ///    CalculateRscu / CalculateRscuCore).
    ///
    /// The string overload is LENIENT by documented design — it never throws on
    /// garbage; it upper-cases, reads codons in non-overlapping in-frame triplets
    /// from index 0, and scores whatever it is given:
    ///   • null OR empty string → returns an EMPTY dictionary via an explicit early
    ///     return (`string.IsNullOrEmpty` guard); NOT a crash, NOT a divide-by-zero.
    ///     Relative_Synonymous_Codon_Usage.md §3.3, §6.1.
    ///   • input length NOT divisible by 3 → CountCodonsCore reads only complete
    ///     triplets (loop guard `i + 3 &lt;= seq.Length`); the trailing 1–2 leftover
    ///     bases are IGNORED, never an IndexOutOfRangeException. §3.3, §6.1.
    ///   • a triplet containing any non-ACGT character → IsValidCodon is false and the
    ///     triplet is EXCLUDED from the counts, never a KeyNotFound. §3.3, §6.1.
    ///   • a family that NEVER occurs (total Σx = 0) → every codon of that family maps
    ///     to 0 (the 0/0 convention), NOT NaN and NOT a DivideByZeroException. §3.2,
    ///     §5.3, §6.1.
    /// The DnaSequence overload throws ArgumentNullException for null (§6.1).
    ///
    /// KEY THEORY INVARIANTS this suite pins directly (Relative_Synonymous_Codon_Usage.md
    /// §2.4), each derived INDEPENDENTLY from Sharp &amp; Li 1986 and the standard
    /// genetic-code degeneracy — never echoed off the implementation's own arrays:
    ///   • INV-01: RSCU = (n_i·x_{i,j})/Σx for a present family (the worked example
    ///     §7.1 is pinned to exact values).
    ///   • INV-02: equal usage within a family ⇒ every RSCU = 1 (no bias).
    ///   • INV-03: 0 ≤ RSCU ≤ n_i for every codon (a randomized sweep asserts this).
    ///   • INV-04: Σ_j RSCU_{i,j} = n_i for every PRESENT family (a randomized sweep
    ///     asserts the per-family sum to within float tolerance).
    ///   • INV-05: a single-codon family (Met ATG, Trp TGG) ⇒ RSCU = 1 when present.
    /// No value may ever be NaN or ±Infinity (asserted on every codon of every input).
    ///
    /// The degeneracy n_i and the synonymous families used to compute every expected
    /// value are derived here from the standard genetic code (the same RNA table
    /// <see cref="StandardGeneticCode"/> already used by the sibling units, translated
    /// to DNA), NOT from CodonUsageAnalyzer's internal table — so a test would fail if
    /// the implementation's family grouping disagreed with the genetic code.
    ///
    /// Determinism note: CalculateRscu is a pure function of the codon counts with no
    /// randomness. The randomized boundary sweep uses a LOCALLY-seeded `new Random(seed)`
    /// (never a shared static Rng); every generated sequence is fully reproducible.
    /// ───────────────────────────────────────────────────────────────────────────

    #region Helpers (RSCU)

    /// <summary>
    /// The standard genetic code in DNA notation (codon → one-letter amino acid,
    /// '*' = stop), derived INDEPENDENTLY of CodonUsageAnalyzer from the RNA-keyed
    /// <see cref="StandardGeneticCode"/> table by replacing U with T. Used to compute
    /// the synonymous families and degeneracies for the expected RSCU values so the
    /// tests never echo the implementation's own table.
    /// </summary>
    private static readonly Dictionary<string, string> DnaGeneticCode =
        StandardGeneticCode.ToDictionary(kv => kv.Key.Replace('U', 'T'), kv => kv.Value);

    /// <summary>Degeneracy n_i = number of synonymous codons sharing the codon's amino acid.</summary>
    private static int Degeneracy(string codon)
    {
        string aa = DnaGeneticCode[codon];
        return DnaGeneticCode.Count(kv => kv.Value == aa);
    }

    /// <summary>The synonymous family (all codons) of the amino acid encoded by <paramref name="codon"/>.</summary>
    private static List<string> Family(string codon)
    {
        string aa = DnaGeneticCode[codon];
        return DnaGeneticCode.Where(kv => kv.Value == aa).Select(kv => kv.Key).ToList();
    }

    #endregion

    #region Positive sanity — the exact worked example (§7.1)

    /// <summary>
    /// Positive sanity (KEY): the §7.1 worked example. Input CTGCTGCTGCTA → Leu
    /// (6-fold) with CTG×3, CTA×1, family total Σx = 4, n_i = 6. Computed
    /// independently from RSCU = (n_i·x)/Σx (Sharp &amp; Li 1986):
    ///   RSCU(CTG) = 6·3/4 = 4.5;  RSCU(CTA) = 6·1/4 = 1.5;
    ///   RSCU(TTA)=RSCU(TTG)=RSCU(CTT)=RSCU(CTC) = 0;  Σ over the family = 6.0 (= n_i).
    /// This both verifies the formula at exact, hand-checkable values and proves the
    /// boundary targets below are measured against a working happy path, not a
    /// uniformly-broken method. (Relative_Synonymous_Codon_Usage.md §7.1, INV-01/04.)
    /// </summary>
    [Test]
    public void CalculateRscu_LeuWorkedExample_MatchesDocumentedExactValues()
    {
        const string coding = "CTGCTGCTGCTA"; // Leu (L) 6-fold: CTG×3, CTA×1

        var rscu = CodonUsageAnalyzer.CalculateRscu(coding);

        rscu["CTG"].Should().BeApproximately(4.5, 1e-12, "RSCU(CTG) = 6·3/4 = 4.5 (§7.1)");
        rscu["CTA"].Should().BeApproximately(1.5, 1e-12, "RSCU(CTA) = 6·1/4 = 1.5 (§7.1)");
        rscu["TTA"].Should().BeApproximately(0.0, 1e-12, "unused Leu synonym ⇒ 0");
        rscu["TTG"].Should().BeApproximately(0.0, 1e-12, "unused Leu synonym ⇒ 0");
        rscu["CTT"].Should().BeApproximately(0.0, 1e-12, "unused Leu synonym ⇒ 0");
        rscu["CTC"].Should().BeApproximately(0.0, 1e-12, "unused Leu synonym ⇒ 0");

        // INV-04 — the present family sums to its degeneracy n_i = 6.
        Family("CTG").Sum(c => rscu[c]).Should().BeApproximately(6.0, 1e-12,
            "Σ_j RSCU over the present Leu family equals n_i = 6 (INV-04)");
    }

    /// <summary>
    /// Positive sanity (INV-02): perfectly EQUAL usage within a family ⇒ every RSCU = 1
    /// (the no-bias reference value). A 4-fold family (Ala) with each of GCT/GCC/GCA/GCG
    /// used once: x = Σx/n_i for all four, so RSCU = 1 for each. Computed independently
    /// from the definition. (Relative_Synonymous_Codon_Usage.md §2.4 INV-02.)
    /// </summary>
    [Test]
    public void CalculateRscu_UniformFamilyUsage_AllRscuEqualOne()
    {
        const string coding = "GCTGCCGCAGCG"; // Ala (A) 4-fold, each synonym once

        var rscu = CodonUsageAnalyzer.CalculateRscu(coding);

        foreach (var codon in Family("GCT"))
            rscu[codon].Should().BeApproximately(1.0, 1e-12,
                $"equal usage ⇒ RSCU({codon}) = 1, the no-bias value (INV-02)");
    }

    #endregion

    #region BE — Boundary: a SINGLE codon (one synonym takes all usage ⇒ RSCU = n_i)

    /// <summary>
    /// BE (checklist "single codon"): a sequence of ONE codon used once means that
    /// codon takes ALL of its family's usage, so RSCU = n_i for it and 0 for every
    /// other synonym — and Σ over the family = n_i (INV-04), RSCU ∈ [0, n_i] (INV-03).
    /// Verified across families of every degeneracy in the standard code: 6-fold (Leu
    /// CTG ⇒ 6), 4-fold (Ala GCT ⇒ 4), 3-fold (Ile ATT ⇒ 3), 2-fold (Phe TTT ⇒ 2).
    /// The expected n_i is derived independently from the genetic-code degeneracy, not
    /// the implementation. (Relative_Synonymous_Codon_Usage.md §2.4 INV-03, §6.1.)
    /// </summary>
    [TestCase("CTG", TestName = "CalculateRscu_SingleCodon_Leu6Fold_RscuEqualsDegeneracy")]
    [TestCase("GCT", TestName = "CalculateRscu_SingleCodon_Ala4Fold_RscuEqualsDegeneracy")]
    [TestCase("ATT", TestName = "CalculateRscu_SingleCodon_Ile3Fold_RscuEqualsDegeneracy")]
    [TestCase("TTT", TestName = "CalculateRscu_SingleCodon_Phe2Fold_RscuEqualsDegeneracy")]
    public void CalculateRscu_SingleCodon_RscuEqualsDegeneracy(string codon)
    {
        int n = Degeneracy(codon);
        n.Should().BeGreaterThan(1, "this case targets a multi-codon family");

        Dictionary<string, double> rscu = null!;
        var act = () => rscu = CodonUsageAnalyzer.CalculateRscu(codon);

        act.Should().NotThrow("a single in-frame codon is a valid minimal coding input");

        // The one present codon takes all family usage ⇒ RSCU = n_i (its upper bound).
        rscu[codon].Should().BeApproximately(n, 1e-12,
            $"a single use of {codon} takes all of its {n}-fold family usage ⇒ RSCU = n_i (INV-03 upper bound)");

        // Every other synonym is unused ⇒ 0.
        foreach (var syn in Family(codon).Where(c => c != codon))
            rscu[syn].Should().BeApproximately(0.0, 1e-12, $"{syn} is unused ⇒ RSCU 0");

        // INV-04 — the present family sums to n_i; INV-03 — every codon in [0, n_i].
        Family(codon).Sum(c => rscu[c]).Should().BeApproximately(n, 1e-12,
            "Σ_j RSCU over the present family equals n_i (INV-04)");
        foreach (var c in Family(codon))
            rscu[c].Should().BeInRange(0.0, n + 1e-12, $"RSCU({c}) ∈ [0, n_i] (INV-03)");
    }

    /// <summary>
    /// BE (INV-05): a single SINGLE-codon-family codon. Met (ATG) and Trp (TGG) have
    /// degeneracy 1, so when present their RSCU is x/x = 1 exactly — there is no
    /// synonym to be biased toward. Pins INV-05 at the n_i = 1 extreme.
    /// (Relative_Synonymous_Codon_Usage.md §2.4 INV-05, §6.1.)
    /// </summary>
    [TestCase("ATG", TestName = "CalculateRscu_SingleCodon_Met_RscuIsOne")] // n_i = 1
    [TestCase("TGG", TestName = "CalculateRscu_SingleCodon_Trp_RscuIsOne")] // n_i = 1
    public void CalculateRscu_SingleCodonFamilyPresent_RscuIsOne(string codon)
    {
        Degeneracy(codon).Should().Be(1, "Met / Trp are the single-codon amino acids");

        var rscu = CodonUsageAnalyzer.CalculateRscu(codon);

        rscu[codon].Should().BeApproximately(1.0, 1e-12,
            $"a single-codon family ⇒ RSCU({codon}) = x/x = 1 when present (INV-05)");
    }

    #endregion

    #region BE — Boundary: a MISSING amino acid (family total 0 ⇒ defined 0, no NaN)

    /// <summary>
    /// BE (checklist "missing amino acid"): an amino acid that NEVER occurs in the
    /// input has family total Σx = 0. The 0/0 case must resolve to a DEFINED 0 for
    /// every codon of that family — NOT NaN, NOT ±Infinity, NOT a DivideByZeroException
    /// (Relative_Synonymous_Codon_Usage.md §3.2, §5.3 "Intentionally simplified",
    /// §6.1). The input encodes ONLY Phe (TTT/TTC), so every OTHER amino acid is an
    /// absent family; we assert the whole Leu family (6-fold) and both single-codon
    /// families Met/Trp are an exact, finite 0, and that NO codon anywhere is NaN/Inf.
    /// </summary>
    [Test]
    public void CalculateRscu_MissingAminoAcid_AbsentFamilyIsDefinedZero()
    {
        const string onlyPhe = "TTTTTCTTTTTC"; // Phe (F) only: TTT×2, TTC×2

        Dictionary<string, double> rscu = null!;
        var act = () => rscu = CodonUsageAnalyzer.CalculateRscu(onlyPhe);

        act.Should().NotThrow("an absent family is the 0/0 case, resolved to 0 — never a div-by-zero");

        // Leu never occurs ⇒ every Leu codon is a finite, exact 0.
        foreach (var leu in Family("CTG"))
        {
            rscu[leu].Should().BeApproximately(0.0, 1e-12,
                $"{leu} (Leu) is absent ⇒ family total 0 ⇒ defined RSCU 0, not NaN");
            double.IsNaN(rscu[leu]).Should().BeFalse($"{leu}: absent family must never yield NaN");
            double.IsInfinity(rscu[leu]).Should().BeFalse($"{leu}: absent family must never yield ±Infinity");
        }

        // The single-codon families Met/Trp are absent here too ⇒ 0 (not the present-1).
        rscu["ATG"].Should().BeApproximately(0.0, 1e-12, "Met is absent ⇒ 0 (not the present-family 1)");
        rscu["TGG"].Should().BeApproximately(0.0, 1e-12, "Trp is absent ⇒ 0");

        // The PRESENT family (Phe) is well-defined and finite: uniform usage ⇒ 1 each.
        rscu["TTT"].Should().BeApproximately(1.0, 1e-12, "Phe used uniformly (×2 each) ⇒ RSCU 1 (INV-02)");
        rscu["TTC"].Should().BeApproximately(1.0, 1e-12, "Phe used uniformly (×2 each) ⇒ RSCU 1 (INV-02)");

        // Global no-NaN / no-Inf guarantee over EVERY codon in the result.
        rscu.Values.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v),
            "no RSCU value may ever be NaN or ±Infinity on any input");
    }

    #endregion

    #region BE — Boundary: empty / null input

    /// <summary>
    /// BE (checklist "empty"): the string overload returns an EMPTY dictionary for
    /// null OR empty input via the documented early return — NOT a crash, NOT a NaN,
    /// NOT a divide-by-zero over zero counts. (Relative_Synonymous_Codon_Usage.md
    /// §3.3, §6.1.)
    /// </summary>
    [TestCase(null, TestName = "CalculateRscu_Null_IsEmptyDictNoThrow")]
    [TestCase("", TestName = "CalculateRscu_Empty_IsEmptyDictNoThrow")]
    public void CalculateRscu_NullOrEmptyString_ReturnsEmptyDictionary(string? input)
    {
        Dictionary<string, double> rscu = null!;
        var act = () => rscu = CodonUsageAnalyzer.CalculateRscu(input!);

        act.Should().NotThrow("null/empty string is a defined no-op early return, not an error");
        rscu.Should().BeEmpty("no codons exist, so no RSCU value is produced");
    }

    /// <summary>
    /// BE: the DnaSequence overload guards null with ArgumentNullException — its
    /// documented, INTENTIONAL validation exception (a disciplined failure, the
    /// acceptable alternative to a defined value under the fuzz bar).
    /// (Relative_Synonymous_Codon_Usage.md §3.3, §6.1.)
    /// </summary>
    [Test]
    public void CalculateRscu_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => CodonUsageAnalyzer.CalculateRscu((DnaSequence)null!);

        act.Should().Throw<ArgumentNullException>(
            "a null DnaSequence is rejected by the documented input guard, not a NullReference");
    }

    #endregion

    #region BE/MC — Boundary: length not %3 and non-ACGT triplets (counts unaffected)

    /// <summary>
    /// BE: a trailing partial codon (length not a multiple of 3) is IGNORED — the
    /// complete in-frame triplets are read and the leftover 1–2 bases contribute
    /// nothing, never an IndexOutOfRangeException (CountCodonsCore loop guard
    /// `i + 3 &lt;= seq.Length`). CTGCTGCTGCTA + 1 or 2 trailing bases must yield the
    /// SAME RSCU as the worked example. (Relative_Synonymous_Codon_Usage.md §3.3, §6.1.)
    /// </summary>
    [TestCase("CTGCTGCTGCTAA", TestName = "CalculateRscu_LenMod3Is1_TrimsTrailingBase")]  // 13 = 4 codons + 1
    [TestCase("CTGCTGCTGCTAAC", TestName = "CalculateRscu_LenMod3Is2_TrimsTrailingTwo")] // 14 = 4 codons + 2
    public void CalculateRscu_LengthNotDivisibleBy3_TrimsPartialCodon(string input)
    {
        Dictionary<string, double> rscu = null!;
        var act = () => rscu = CodonUsageAnalyzer.CalculateRscu(input);

        act.Should().NotThrow("a trailing partial codon must be ignored, never indexed out of range");

        // Only CTG×3, CTA×1 are complete in-frame Leu codons ⇒ same exact values as §7.1.
        rscu["CTG"].Should().BeApproximately(4.5, 1e-12, "the partial codon contributes nothing");
        rscu["CTA"].Should().BeApproximately(1.5, 1e-12, "the partial codon contributes nothing");
        rscu.Values.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v),
            "no RSCU value may be NaN or ±Infinity");
    }

    /// <summary>
    /// MC: a triplet containing any non-ACGT character is EXCLUDED from the counts
    /// (IsValidCodon over {A,C,G,T}), never a KeyNotFoundException. The valid Leu
    /// codons around the garbage triplet still produce the §7.1 RSCU, and no value is
    /// NaN/Inf. Covers digits, the ambiguity code N, lowercase-invalid, unicode and an
    /// embedded null byte as the excluded triplet.
    /// (Relative_Synonymous_Codon_Usage.md §3.3, §6.1.)
    /// </summary>
    [TestCase("CTGCTG123CTGCTA", TestName = "CalculateRscu_NonAcgt_Digits_ExcludedNoCrash")]
    [TestCase("CTGCTGNNNCTGCTA", TestName = "CalculateRscu_NonAcgt_AmbiguityN_ExcludedNoCrash")]
    [TestCase("CTGCTG\0\0\0CTGCTA", TestName = "CalculateRscu_NonAcgt_NullBytes_ExcludedNoCrash")]
    [TestCase("CTGCTGαβγCTGCTA", TestName = "CalculateRscu_NonAcgt_GreekLetters_ExcludedNoCrash")]
    public void CalculateRscu_NonAcgtTriplet_ExcludedFromCounts(string input)
    {
        Dictionary<string, double> rscu = null!;
        var act = () => rscu = CodonUsageAnalyzer.CalculateRscu(input);

        act.Should().NotThrow("a non-ACGT triplet must be excluded from counts, never KeyNotFound/IndexOutOfRange");

        // The four valid Leu codons (CTG×3, CTA×1) are exactly the §7.1 example.
        rscu["CTG"].Should().BeApproximately(4.5, 1e-12, "valid Leu codons score as §7.1; the bad triplet is excluded");
        rscu["CTA"].Should().BeApproximately(1.5, 1e-12, "valid Leu codons score as §7.1; the bad triplet is excluded");
        rscu.Values.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v),
            "no RSCU value may be NaN or ±Infinity");
    }

    #endregion

    #region BE — Randomized boundary sweep (per-family sum = n_i, RSCU ∈ [0, n_i], finite)

    /// <summary>
    /// BE (randomized sweep, INV-03/INV-04): hundreds of random in-frame coding
    /// sequences must, on EVERY input, satisfy the core algebraic contract derived
    /// independently from Sharp &amp; Li 1986:
    ///   • every RSCU value is finite (no NaN, no ±Infinity);
    ///   • every RSCU value lies in [0, n_i] for its family (INV-03);
    ///   • for every PRESENT family (Σx &gt; 0) the per-family RSCU sum equals n_i to
    ///     within float tolerance (INV-04), and for every ABSENT family the sum is 0.
    /// The family / degeneracy used for the bounds and the per-family sum are computed
    /// from the standard genetic code, NOT from the implementation. A LOCALLY-seeded
    /// `new Random(seed)` keeps every trial reproducible; CancelAfter guards a hang.
    /// (Relative_Synonymous_Codon_Usage.md §2.4 INV-03/INV-04.)
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void CalculateRscu_RandomizedSweep_FamilySumIsDegeneracyAndBounded()
    {
        const string bases = "ACGT";
        var allCodons = DnaGeneticCode.Keys.ToList();
        // Group families once (by representative codon) for the per-family checks.
        var families = allCodons
            .GroupBy(c => DnaGeneticCode[c])
            .Select(g => (Codons: g.ToList(), N: g.Count()))
            .ToList();

        for (int trial = 0; trial < 400; trial++)
        {
            var rng = new Random(20140513 + trial); // locally seeded, reproducible
            int codonCount = rng.Next(1, 60);       // 1..59 in-frame codons (includes the single-codon edge)

            var sb = new System.Text.StringBuilder();
            for (int i = 0; i < codonCount; i++)
                for (int j = 0; j < 3; j++)
                    sb.Append(bases[rng.Next(bases.Length)]);
            string seq = sb.ToString();

            Dictionary<string, double> rscu = null!;
            var act = () => rscu = CodonUsageAnalyzer.CalculateRscu(seq);
            act.Should().NotThrow($"trial {trial}: random valid coding input must never crash");

            // Finiteness over every codon.
            rscu.Values.Should().OnlyContain(v => !double.IsNaN(v) && !double.IsInfinity(v),
                $"trial {trial}: no RSCU value may be NaN or ±Infinity");

            foreach (var (codons, n) in families)
            {
                double familySum = 0;
                foreach (var codon in codons)
                {
                    double value = rscu.GetValueOrDefault(codon, 0);
                    value.Should().BeInRange(0.0, n + 1e-9,
                        $"trial {trial}: RSCU({codon}) ∈ [0, n_i={n}] (INV-03)");
                    familySum += value;
                }

                // INV-04: present family ⇒ Σ = n_i; absent family ⇒ Σ = 0.
                if (familySum > 1e-9)
                    familySum.Should().BeApproximately(n, 1e-6,
                        $"trial {trial}: Σ_j RSCU over a present family equals n_i = {n} (INV-04)");
                else
                    familySum.Should().BeApproximately(0.0, 1e-9,
                        $"trial {trial}: an absent family sums to 0 (INV-04)");
            }
        }
    }

    #endregion

    #endregion
}
