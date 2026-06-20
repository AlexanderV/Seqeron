using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area NUCLEOTIDE SEQUENCE SUMMARY aggregator.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, DivideByZeroException, OverflowException, …). Every
/// input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception. A raw runtime exception or a
/// hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-SUMMARY-001 — Sequence Summary (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 128.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — empty, single base, very long, mixed case.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 128 targets.
///
/// ───────────────────────────────────────────────────────────────────────────
/// SCOPING vs row 121 (SEQ-COMPOSITION-001) and row 127 (SEQ-STATS-001)
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-SUMMARY-001 owns the NUCLEOTIDE aggregator
/// SequenceStatistics.SummarizeNucleotideSequence(string?) — a pure aggregation
/// that bundles length, GC content, Shannon entropy, linguistic complexity,
/// melting temperature, and an A/T/G/C/U/N composition map into one
/// SequenceSummary record (docs/algorithms/Statistics/Sequence_Summary.md §1).
/// It is DISTINCT from:
///   • row 121 SEQ-COMPOSITION-001 — CalculateNucleotideComposition scalar
///     counts/fractions/skews (fuzzed in SequenceCompositionStatFuzzTests.cs);
///   • row 127 SEQ-STATS-001 — the PROTEIN AminoAcidComposition aggregate record
///     (fuzzed in SequenceStatisticsAggregateFuzzTests.cs).
/// This file is scoped to the nucleotide SUMMARY aggregate and asserts EVERY
/// reported field of SequenceSummary.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The summary contract under test (Sequence_Summary.md §2.2, §3, §6; source)
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.SummarizeNucleotideSequence(string?)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 990–1020), returning a SequenceSummary record struct with fields
///    Length, GcContent, Entropy, Complexity, MeltingTemperature, Composition.
///
/// Field definitions (Sequence_Summary.md §2.2, §3.2; INV-01..INV-06):
///   • Length  = |S| raw character count = CalculateNucleotideComposition(S).Length
///               = S.Length (INCLUDES N and other chars, not just A/T/G/C/U).
///   • GcContent = GC FRACTION in [0,1] = (#G+#C)/(#A+#T+#G+#C+#U) over counted
///               bases, case-insensitive; 0 for empty input.
///   • Entropy = Shannon entropy H = −Σ p·log₂p over per-symbol frequencies (bits);
///               0 for empty input. Maximum is log₂k for k equiprobable symbols.
///   • Complexity = linguistic complexity (mean of vocabulary-usage ratios across
///               word sizes k=1..6), in [0,1]; 0 for empty input.
///   • MeltingTemperature = Wallace 2(A+T)+4(G+C) when |S| < 14, else GC/Marmur-Doty
///               64.9 + 41·(GC−16.4)/N; 0 for empty input
///               (ThermoConstants.WallaceMaxLength = 14, strict <).
///   • Composition = a 6-entry map {A,T,G,C,U,N} of the composition counts.
///
/// Empty/null handling (§3.3): null is normalised to empty (SequenceStatistics.cs
///   line 995), and an empty input yields a DEGENERATE summary — Length 0,
///   GcContent 0, Entropy 0, Complexity 0, MeltingTemperature 0, all six
///   composition counts 0. NO exception, NO DivideByZero.
///
/// Aggregation contract (§2.4 INV-01..INV-06): every summary field reproduces,
///   bit-for-bit, the value its canonical per-metric method returns on the same
///   input. The fuzz oracle below re-derives each field from those methods.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceSummaryFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The six composition symbols the summary map reports (§2.2).</summary>
    private static readonly char[] CompositionSymbols = { 'A', 'T', 'G', 'C', 'U', 'N' };

    /// <summary>Wallace-vs-GC branch boundary (ThermoConstants.WallaceMaxLength), strict &lt;.</summary>
    private const int WallaceMaxLength = 14;

    /// <summary>Generates a random string of arbitrary BMP code points (0x0000–0xFFFF),
    /// spanning control characters, the null byte, lone surrogate halves, unicode
    /// letters and digits — random-byte fuzz fodder for the aggregator.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>Generates a random nucleotide-ish string over a configurable alphabet.</summary>
    private static string RandomOver(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal well-formedness contract that must hold for ANY input:
    /// every reported field is finite and in its documented range, the composition
    /// map has exactly the six A/T/G/C/U/N keys with non-negative counts, and no
    /// field is NaN/Infinity (overflow / DivideByZero guard).
    /// — Sequence_Summary.md §2.4 INV-01..INV-07, §3.2.
    /// </summary>
    private static void AssertWellFormed(SequenceStatistics.SequenceSummary s)
    {
        s.Length.Should().BeGreaterThanOrEqualTo(0);

        // GC fraction is a probability in [0,1] (NOT a percentage) — INV-07.
        s.GcContent.Should().BeInRange(0.0, 1.0);
        double.IsFinite(s.GcContent).Should().BeTrue("GC content must be finite");

        // Shannon entropy is non-negative and finite (bits).
        s.Entropy.Should().BeGreaterThanOrEqualTo(0.0);
        double.IsFinite(s.Entropy).Should().BeTrue("entropy must be finite");

        // Linguistic complexity is a non-negative, finite vocabulary-usage ratio.
        // NOTE: the documented [0,1] bound (INV-07) holds for DNA/RNA fragments, where
        // the 4-symbol alphabet caps the observed/possible word ratio at 1. For arbitrary
        // unicode garbage the observed distinct-symbol count can exceed the 4^k cap, so the
        // ratio may exceed 1; the [0,1] bound is therefore asserted only on nucleotide
        // inputs (AssertDnaWellFormed), not in this universal helper.
        s.Complexity.Should().BeGreaterThanOrEqualTo(0.0);
        double.IsFinite(s.Complexity).Should().BeTrue("complexity must be finite");

        // Melting temperature is finite (can be negative for the GC formula on very
        // low-GC long sequences, so no lower bound is asserted here).
        double.IsFinite(s.MeltingTemperature).Should().BeTrue("Tm must be finite, not NaN/Inf");

        // Composition map: exactly the six documented symbols, all counts ≥ 0.
        s.Composition.Should().NotBeNull();
        s.Composition.Keys.Should().BeEquivalentTo(CompositionSymbols,
            "the summary reports exactly the A,T,G,C,U,N counts");
        foreach (var (_, count) in s.Composition)
            count.Should().BeGreaterThanOrEqualTo(0);

        // The five standard bases counted toward GC/AT can never exceed the raw length.
        int countedBases = s.Composition['A'] + s.Composition['T'] + s.Composition['G'] +
                           s.Composition['C'] + s.Composition['U'];
        (countedBases + s.Composition['N']).Should().BeLessThanOrEqualTo(s.Length,
            "A+T+G+C+U+N counts cannot exceed the raw character length");
    }

    /// <summary>
    /// Stronger well-formedness for DNA inputs: asserts the universal contract PLUS the
    /// documented Complexity ∈ [0,1] bound (INV-07). That bound is the vocabulary-usage
    /// ratio observed/min(4^k, n−k+1); it is ≤ 1 only when the alphabet is the 4-letter
    /// DNA alphabet {A,T,G,C}, because the 4^k cap in CalculateLinguisticComplexity
    /// assumes exactly four symbols. Inputs containing U or N (≥ 5 distinct symbols) can
    /// push the k=1 ratio above 1, so this helper is applied ONLY to {A,T,G,C} sequences;
    /// other nucleotide inputs use AssertWellFormed (which still pins finiteness and ≥ 0).
    /// </summary>
    private static void AssertDnaWellFormed(SequenceStatistics.SequenceSummary s)
    {
        AssertWellFormed(s);
        s.Complexity.Should().BeInRange(0.0, 1.0,
            "for the 4-letter DNA alphabet the vocabulary-usage ratio is bounded by 1 (INV-07)");
    }

    /// <summary>
    /// The aggregation oracle: re-derives every field directly from the canonical
    /// per-metric methods and asserts the summary reproduces each one bit-for-bit
    /// (INV-01..INV-06). This is the core of SEQ-SUMMARY-001 — the summary must be
    /// a faithful copy of its components, never a recomputation that could drift.
    /// </summary>
    private static void AssertMatchesComponentMetrics(string? input)
    {
        string seq = input ?? string.Empty;
        var comp = SequenceStatistics.CalculateNucleotideComposition(seq);
        double entropy = SequenceStatistics.CalculateShannonEntropy(seq);
        double complexity = SequenceStatistics.CalculateLinguisticComplexity(seq);
        double tm = SequenceStatistics.CalculateMeltingTemperature(
            seq, useWallaceRule: seq.Length < WallaceMaxLength);

        var s = SequenceStatistics.SummarizeNucleotideSequence(input);

        s.Length.Should().Be(comp.Length, "INV-01: Length copies composition Length");
        s.GcContent.Should().Be(comp.GcContent, "INV-02: GcContent copies composition GcContent");
        s.Entropy.Should().Be(entropy, "INV-03: Entropy = CalculateShannonEntropy");
        s.Complexity.Should().Be(complexity, "INV-04: Complexity = CalculateLinguisticComplexity");
        s.MeltingTemperature.Should().Be(tm, "INV-05: Tm = CalculateMeltingTemperature(len<14)");

        // INV-06: composition map equals the composition record's counts.
        s.Composition['A'].Should().Be(comp.CountA);
        s.Composition['T'].Should().Be(comp.CountT);
        s.Composition['G'].Should().Be(comp.CountG);
        s.Composition['C'].Should().Be(comp.CountC);
        s.Composition['U'].Should().Be(comp.CountU);
        s.Composition['N'].Should().Be(comp.CountN);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-SUMMARY-001 — nucleotide sequence summary : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact result

    /// <summary>
    /// Positive baseline (not a boundary): the worked example from the algorithm doc
    /// must reproduce EXACTLY. "ATGCATGC" → A=2,T=2,G=2,C=2:
    ///   Length 8; GcContent = 4/8 = 0.5; four equiprobable symbols → Entropy log₂4 = 2.0;
    ///   length 8 &lt; 14 → Wallace Tm = 2·(2+2) + 4·(2+2) = 8 + 16 = 24.0 °C;
    ///   Complexity = mean vocabulary-usage ratio = 0.8396825396825397 (hand-derived,
    ///   externally re-grounded — docs/Validation/FINDINGS_REGISTER.md A39).
    /// Confirms the suite asserts the BUSINESS contract, not just non-throwing.
    /// — Sequence_Summary.md §7.1; SequenceStatistics.cs lines 990–1020.
    /// </summary>
    [Test]
    public void Summary_DocWorkedExample_MatchesHandComputedFields()
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGC");

        s.Length.Should().Be(8);
        s.GcContent.Should().BeApproximately(0.5, Tolerance, "GC = (2+2)/8");
        s.Entropy.Should().BeApproximately(2.0, Tolerance, "four equiprobable symbols → log2 4");
        s.MeltingTemperature.Should().BeApproximately(24.0, Tolerance,
            "len 8 < 14 → Wallace 2*(A+T)+4*(G+C) = 2*4 + 4*4");
        s.Complexity.Should().BeApproximately(0.8396825396825397, 1e-10,
            "externally-derived vocabulary-usage-mean lock (FINDINGS_REGISTER A39)");

        s.Composition['A'].Should().Be(2);
        s.Composition['T'].Should().Be(2);
        s.Composition['G'].Should().Be(2);
        s.Composition['C'].Should().Be(2);
        s.Composition['U'].Should().Be(0);
        s.Composition['N'].Should().Be(0);

        AssertDnaWellFormed(s);
        AssertMatchesComponentMetrics("ATGCATGC");
    }

    /// <summary>
    /// Positive baseline for the GC/Marmur-Doty branch: a 16-mer (≥ 14) selects the
    /// GC formula. "ATGCATGCATGCATGC" → GC = 8/16 = 0.5; Marmur-Doty Tm =
    /// 64.9 + 41·(8 − 16.4)/16 = 43.375 (hand-derived, FINDINGS_REGISTER A39).
    /// Guards the branch-selection field of the aggregator at the long-input boundary.
    /// </summary>
    [Test]
    public void Summary_SixteenMer_UsesGcFormula_MatchesHandComputedTm()
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGCATGCATGC");

        s.Length.Should().Be(16);
        s.GcContent.Should().BeApproximately(0.5, Tolerance);
        // GC formula because length 16 ≥ 14: 64.9 + 41*(8 - 16.4)/16.
        s.MeltingTemperature.Should().BeApproximately(64.9 + 41.0 * (8 - 16.4) / 16, Tolerance,
            "len 16 ≥ 14 → Marmur-Doty GC formula, not Wallace");
        AssertDnaWellFormed(s);
        AssertMatchesComponentMetrics("ATGCATGCATGCATGC");
    }

    #endregion

    #region BE — Boundary: empty / null (degenerate summary, no DivideByZero)

    /// <summary>
    /// BE: the empty string is the lower size boundary. Documented as a DEGENERATE
    /// summary — Length 0, GcContent 0, Entropy 0, Complexity 0, Tm 0, all six
    /// composition counts 0; NO DivideByZero on any fraction denominator, no throw.
    /// — Sequence_Summary.md §3.3, §6.1; SequenceStatistics.cs lines 995, 1013–1019.
    /// </summary>
    [Test]
    public void Summary_EmptyString_IsDegenerateAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.SummarizeNucleotideSequence(string.Empty);

        act.Should().NotThrow("the empty string is a defined boundary, not an error");

        var s = act();
        s.Length.Should().Be(0);
        s.GcContent.Should().Be(0.0);
        s.Entropy.Should().Be(0.0);
        s.Complexity.Should().Be(0.0);
        s.MeltingTemperature.Should().Be(0.0);
        foreach (char sym in CompositionSymbols)
            s.Composition[sym].Should().Be(0, $"empty input ⇒ zero {sym} count");
        AssertWellFormed(s);
    }

    /// <summary>
    /// BE: null is normalised to empty (SequenceStatistics.cs line 995) — identical
    /// degenerate summary, no NullReferenceException. The method signature is
    /// explicitly nullable (string?), so null is a documented input.
    /// </summary>
    [Test]
    public void Summary_Null_IsDegenerateAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.SummarizeNucleotideSequence(null);

        act.Should().NotThrow("null is normalised to empty, not an error");

        var s = act();
        s.Length.Should().Be(0);
        s.GcContent.Should().Be(0.0);
        s.Entropy.Should().Be(0.0);
        s.Complexity.Should().Be(0.0);
        s.MeltingTemperature.Should().Be(0.0);
        s.Composition.Values.Should().OnlyContain(c => c == 0);
        AssertWellFormed(s);
    }

    #endregion

    #region BE — Boundary: single base (length 1)

    /// <summary>
    /// BE: a single G or C is the length-1 lower content boundary for GC = 1.0.
    /// Length 1; that base count 1, all others 0; GcContent = 1/1 = 1.0; a single
    /// symbol → Entropy −1·log₂1 = 0.0; Complexity for k=1 only = 1/1 = 1.0;
    /// length 1 &lt; 14 → Wallace Tm = 4·1 = 4.0. No crash, no DivideByZero.
    /// — Sequence_Summary.md §2.2; SequenceStatistics.cs (entropy/Wallace).
    /// </summary>
    [TestCase('G')]
    [TestCase('C')]
    public void Summary_SingleGcBase_GcIsOne_EntropyZero(char baseChar)
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence(baseChar.ToString());

        s.Length.Should().Be(1);
        s.GcContent.Should().BeApproximately(1.0, Tolerance, "the only base is G or C");
        s.Entropy.Should().BeApproximately(0.0, Tolerance, "one symbol carries no information");
        s.Complexity.Should().BeApproximately(1.0, Tolerance, "k=1 vocabulary fully used");
        s.MeltingTemperature.Should().BeApproximately(4.0, Tolerance, "Wallace 4*(G+C) on one GC base");
        s.Composition[baseChar].Should().Be(1);
        AssertDnaWellFormed(s);
        AssertMatchesComponentMetrics(baseChar.ToString());
    }

    /// <summary>
    /// BE: a single A or T is the length-1 boundary for GC = 0.0 with the base still
    /// counted: Length 1; GcContent 0.0; Entropy 0.0; Wallace Tm = 2·1 = 2.0
    /// (Wallace sums A+T only — SequenceStatistics.cs line 580).
    /// </summary>
    [TestCase('A')]
    [TestCase('T')]
    public void Summary_SingleAtBase_GcIsZero_WallaceTwo(char baseChar)
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence(baseChar.ToString());

        s.Length.Should().Be(1);
        s.GcContent.Should().BeApproximately(0.0, Tolerance, "A/T are not GC");
        s.Entropy.Should().BeApproximately(0.0, Tolerance);
        s.MeltingTemperature.Should().BeApproximately(2.0, Tolerance, "Wallace 2*(A+T) on one AT base");
        s.Composition[baseChar].Should().Be(1);
        AssertDnaWellFormed(s);
        AssertMatchesComponentMetrics(baseChar.ToString());
    }

    /// <summary>
    /// BE: a single 'U' (RNA) is counted (Length 1, CountU 1) with GcContent 0.0, but
    /// the Wallace Tm is 0.0 — Wallace sums only A+T and G+C (SequenceStatistics.cs
    /// line 580: `comp.CountA + comp.CountT`), so a lone U contributes nothing to Tm.
    /// This is the documented per-metric behaviour the summary faithfully copies, and a
    /// distinct boundary from A/T whose Wallace Tm is 2.0.
    /// </summary>
    [Test]
    public void Summary_SingleU_GcZero_WallaceZero()
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence("U");

        s.Length.Should().Be(1);
        s.Composition['U'].Should().Be(1);
        s.GcContent.Should().BeApproximately(0.0, Tolerance, "U is not GC");
        s.Entropy.Should().BeApproximately(0.0, Tolerance);
        s.MeltingTemperature.Should().Be(0.0, "Wallace counts only A/T (not U) → lone U gives Tm 0");
        AssertWellFormed(s);
        AssertMatchesComponentMetrics("U");
    }

    /// <summary>
    /// BE: a single 'N' is counted (Length 1, CountN 1) but is NOT a standard base, so
    /// the GC denominator (A+T+G+C+U) is 0 — the documented guard returns GcContent 0
    /// with NO DivideByZero. Entropy still sees one symbol (0.0). Tm: Wallace counts
    /// only A/T/G/C, so Tm = 0.0. Length (raw) includes the N though GC excludes it.
    /// — SequenceStatistics.cs lines 76 (GC guard), 82 (Length = raw).
    /// </summary>
    [Test]
    public void Summary_SingleN_CountedButNotGc_NoDivideByZero()
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence("N");

        s.Length.Should().Be(1, "raw length includes N");
        s.Composition['N'].Should().Be(1);
        s.GcContent.Should().Be(0.0, "no standard base ⇒ GC denominator 0, guarded to 0");
        s.MeltingTemperature.Should().Be(0.0, "N contributes no A/T/G/C to Wallace");
        AssertWellFormed(s);
        AssertMatchesComponentMetrics("N");
    }

    #endregion

    #region BE — Boundary: very long input (no overflow, terminates, finite)

    /// <summary>
    /// BE: a very long sequence must compute in a non-overflowing width, stay finite,
    /// and TERMINATE (the linguistic-complexity k-mer sets up to k=6 are the dominant
    /// cost). 1,000,000 bases of a repeating "ATGC" → A=T=G=C=250000; GcContent 0.5;
    /// Entropy 2.0; length ≥ 14 → GC/Marmur-Doty Tm finite. Length must be the exact
    /// raw count (no int overflow at 10⁶). [CancelAfter] guards against a hang.
    /// — Sequence_Summary.md §4.3 (O(n)); SequenceStatistics.cs.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Summary_VeryLongRepeat_FieldsFiniteAndExact_NoOverflow()
    {
        const int repeats = 250_000; // 1,000,000 bases total
        string seq = string.Concat(Enumerable.Repeat("ATGC", repeats));
        seq.Length.Should().Be(1_000_000);

        var s = SequenceStatistics.SummarizeNucleotideSequence(seq);

        s.Length.Should().Be(1_000_000, "raw length must not overflow at 10^6");
        s.Composition['A'].Should().Be(repeats);
        s.Composition['T'].Should().Be(repeats);
        s.Composition['G'].Should().Be(repeats);
        s.Composition['C'].Should().Be(repeats);
        s.GcContent.Should().BeApproximately(0.5, Tolerance);
        s.Entropy.Should().BeApproximately(2.0, Tolerance, "four equiprobable symbols");
        double.IsFinite(s.MeltingTemperature).Should().BeTrue("GC-formula Tm must be finite at 10^6");
        AssertDnaWellFormed(s);
    }

    /// <summary>
    /// BE: a very long single-symbol homopolymer is the worst case for GC=1, entropy=0
    /// and complexity boundaries simultaneously. 500,000 G's → Length 500000, all G,
    /// GcContent 1.0, Entropy 0.0 (one symbol), finite GC-formula Tm, no overflow/hang.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Summary_VeryLongHomopolymer_GcOne_EntropyZero_NoOverflow()
    {
        string seq = new string('G', 500_000);

        var s = SequenceStatistics.SummarizeNucleotideSequence(seq);

        s.Length.Should().Be(500_000);
        s.Composition['G'].Should().Be(500_000);
        s.GcContent.Should().BeApproximately(1.0, Tolerance);
        s.Entropy.Should().BeApproximately(0.0, Tolerance, "a single symbol has zero entropy");
        double.IsFinite(s.MeltingTemperature).Should().BeTrue();
        AssertDnaWellFormed(s);
    }

    #endregion

    #region BE — Boundary: mixed case (case-insensitive, identical to uppercase)

    /// <summary>
    /// BE: the summary is case-insensitive (every per-metric method uppercases
    /// internally). A mixed-case sequence MUST yield a summary IDENTICAL to its
    /// all-uppercase form — every numeric field and every composition count. Guards a
    /// case bug that would split 'a'/'A' or miss a lowercase base in GC/entropy/Tm.
    /// — Sequence_Summary.md §3.3, §6.1 (lowercase → identical to uppercase).
    /// </summary>
    [TestCase("atgcatgc")]
    [TestCase("AtGcAtGc")]
    [TestCase("gattaca")]
    [TestCase("GaTtAcAggTccaATGn")]
    [TestCase("auauauaucgcg")]
    public void Summary_MixedCase_EqualsUppercase(string seq)
    {
        var lower = SequenceStatistics.SummarizeNucleotideSequence(seq);
        var upper = SequenceStatistics.SummarizeNucleotideSequence(seq.ToUpperInvariant());

        lower.Length.Should().Be(upper.Length);
        lower.GcContent.Should().Be(upper.GcContent, "case must not affect GC content");
        lower.Entropy.Should().Be(upper.Entropy, "case must not affect entropy");
        lower.Complexity.Should().Be(upper.Complexity, "case must not affect complexity");
        lower.MeltingTemperature.Should().Be(upper.MeltingTemperature, "case must not affect Tm");
        lower.Composition.Should().BeEquivalentTo(upper.Composition, "case must not affect counts");
        AssertWellFormed(lower);
    }

    /// <summary>
    /// BE: a single lowercase 'g' must be recognised as G — GcContent 1.0, CountG 1 —
    /// proving lowercase is not silently dropped from the GC/Tm computation.
    /// </summary>
    [Test]
    public void Summary_SingleLowercaseG_RecognisedAsGc()
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence("g");

        s.Length.Should().Be(1);
        s.Composition['G'].Should().Be(1, "lowercase g is upper-cased before counting");
        s.GcContent.Should().BeApproximately(1.0, Tolerance);
        s.MeltingTemperature.Should().BeApproximately(4.0, Tolerance);
        AssertDnaWellFormed(s);
    }

    #endregion

    #region BE/RB — Random fuzz: never throws, always matches component metrics

    /// <summary>
    /// BE/RB: a large batch of arbitrary BMP strings (control chars, null byte, lone
    /// surrogate halves, unicode letters/digits, mixed lengths incl. 0 and 1) must
    /// NEVER throw and must ALWAYS produce a well-formed summary whose every field
    /// equals the corresponding per-metric method (the aggregation contract). This is
    /// the core fuzz guarantee: no DivideByZero, no IndexOutOfRange, no NullReference,
    /// no overflow/NaN on garbage. [CancelAfter] guards against a hang.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Summary_RandomGarbageStrings_NeverThrow_MatchComponentMetrics()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            SequenceStatistics.SequenceSummary s = default;
            var act = () => s = SequenceStatistics.SummarizeNucleotideSequence(input);

            act.Should().NotThrow($"garbage input (len {len}) must never crash the aggregator");
            AssertWellFormed(s);
            AssertMatchesComponentMetrics(input);
        }
    }

    /// <summary>
    /// BE: a large batch of random NUCLEOTIDE-alphabet strings (A/T/G/C/U/N, lower and
    /// upper case, random lengths down to 0/1) cross-checks the summary against an
    /// independent oracle: re-counted composition, GC fraction over standard bases, and
    /// the field-by-field component equality. Exercises many content shapes including
    /// the single-base, all-N and case boundaries simultaneously.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Summary_RandomNucleotideStrings_FieldsMatchOracle()
    {
        var rng = new Random(424242);
        const string alphabet = "ATGCUNatgcun";

        for (int iteration = 0; iteration < 1500; iteration++)
        {
            int len = rng.Next(0, 300);
            string seq = RandomOver(rng, alphabet, len);

            var s = SequenceStatistics.SummarizeNucleotideSequence(seq);

            // Raw length oracle.
            s.Length.Should().Be(len);

            // Composition oracle (case-folded counts).
            string up = seq.ToUpperInvariant();
            s.Composition['A'].Should().Be(up.Count(c => c == 'A'));
            s.Composition['T'].Should().Be(up.Count(c => c == 'T'));
            s.Composition['G'].Should().Be(up.Count(c => c == 'G'));
            s.Composition['C'].Should().Be(up.Count(c => c == 'C'));
            s.Composition['U'].Should().Be(up.Count(c => c == 'U'));
            s.Composition['N'].Should().Be(up.Count(c => c == 'N'));

            // GC-fraction oracle over standard bases (denominator excludes N).
            int gc = s.Composition['G'] + s.Composition['C'];
            int standard = gc + s.Composition['A'] + s.Composition['T'] + s.Composition['U'];
            double expectedGc = standard > 0 ? (double)gc / standard : 0.0;
            s.GcContent.Should().BeApproximately(expectedGc, Tolerance);

            AssertWellFormed(s);
            AssertMatchesComponentMetrics(seq);
        }
    }

    #endregion
}
