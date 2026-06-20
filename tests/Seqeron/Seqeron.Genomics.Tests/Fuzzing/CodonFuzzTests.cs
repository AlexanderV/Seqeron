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
}
