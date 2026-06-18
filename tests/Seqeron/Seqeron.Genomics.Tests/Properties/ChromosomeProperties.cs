using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for chromosome-level analysis.
///
/// Test Units: CHROM-TELO-001 (this file is shared; CHROM-CENT/KARYO/ANEU/SYNT
/// extend it later in sibling regions).
///
/// CHROM-TELO-001 covers <see cref="ChromosomeAnalyzer.AnalyzeTelomeres"/> and
/// <see cref="ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio"/>:
///   R  positions/lengths valid (INV-01/02, exact repeat-unit multiples);
///   P  detected telomere is a repeat tract (purity, motif containment);
///   M  more repeats → strictly longer measured region (length = 6·k);
///   D  identical inputs → identical result.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Chromosome")]
public class ChromosomeProperties
{
    #region CHROM-TELO-001

    // ---- Domain constants (transcribed independently from Telomere_Analysis.md §2.1) ----

    /// <summary>Canonical vertebrate 3'-end telomere repeat (Meyne 1989).</summary>
    private const string Vertebrate3Prime = "TTAGGG";

    /// <summary>Reverse complement of <see cref="Vertebrate3Prime"/>, expected at the 5' end.</summary>
    private const string Vertebrate5Prime = "CCCTAA";

    // ===================== Generators =====================

    /// <summary>
    /// Generates a 3'-telomeric sequence: a non-telomeric prefix of A's whose length is a
    /// multiple of 6 (keeps repeat-sized windows aligned and never extends the tract — the
    /// last 6-base window "AAAAAA" matches TTAGGG at only 1/6 &lt; 0.7), followed by exactly
    /// <c>k</c> pure TTAGGG units placed at the very end. Yields (sequence, k).
    /// </summary>
    private static Arbitrary<(string seq, int k)> Telomeric3PrimeArbitrary() =>
        (from prefixUnits in Gen.Choose(0, 30)
         from k in Gen.Choose(1, 400)
         let prefix = new string('A', 6 * prefixUnits)
         let tract = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k))
         select (prefix + tract, k)).ToArbitrary();

    /// <summary>
    /// Generates a 5'-telomeric sequence: exactly <c>k</c> pure CCCTAA units at the start,
    /// followed by a non-telomeric A-suffix. The first post-tract window "AAAAAA" matches
    /// CCCTAA at only 2/6 &lt; 0.7, so the tract is not extended. Yields (sequence, k).
    /// </summary>
    private static Arbitrary<(string seq, int k)> Telomeric5PrimeArbitrary() =>
        (from k in Gen.Choose(1, 400)
         from suffixLen in Gen.Choose(0, 120)
         let tract = string.Concat(Enumerable.Repeat(Vertebrate5Prime, k))
         select (tract + new string('A', suffixLen), k)).ToArbitrary();

    /// <summary>
    /// Generates non-telomeric DNA: only A and T. Neither motif can reach the 70% threshold —
    /// against TTAGGG the three G positions never match (max 3/6 = 0.5), and against CCCTAA the
    /// three C positions never match (max 3/6 = 0.5). Empty allowed.
    /// </summary>
    private static Arbitrary<string> NonTelomericArbitrary() =>
        Gen.Elements('A', 'T').ArrayOf().Select(a => new string(a)).ToArbitrary();

    /// <summary>Generates arbitrary DNA (A/C/G/T, possibly empty) for general invariants.</summary>
    private static Arbitrary<string> AnyDnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Select(a => new string(a)).ToArbitrary();

    /// <summary>Generates non-negative (tsRatio, refRatio &gt; 0, refLength) triples for the T/S oracle.</summary>
    private static Arbitrary<(double ts, double refRatio, double refLen)> TsRatioArbitrary() =>
        (from ts in Gen.Choose(0, 100_000).Select(i => i / 1000.0)
         from rr in Gen.Choose(1, 100_000).Select(i => i / 1000.0)
         from rl in Gen.Choose(0, 100_000).Select(i => (double)i)
         select (ts, rr, rl)).ToArbitrary();

    // ===================== Independent oracles =====================

    /// <summary>
    /// Independent T/S oracle, transcribed verbatim from Telomere_Analysis.md §2.2 /
    /// INV-04 — NOT routed through production.
    /// </summary>
    private static double ExpectedTsLength(double tsRatio, double referenceRatio, double referenceLength) =>
        referenceLength * tsRatio / referenceRatio;

    #endregion

    #region CHROM-TELO-001 — M: more repeats → longer region (central business invariant)

    /// <summary>
    /// INV-M (3' end): a tract of <c>k</c> pure TTAGGG units placed at the very end measures
    /// exactly <c>6·k</c> bases — the measured length is a linear function of repeat count,
    /// hence strictly increasing in k. The internal scan counts only complete 6-bp windows,
    /// and the A-prefix terminates the scan, so no over- or under-count occurs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_ThreePrimeLength_EqualsSixTimesRepeatCount()
    {
        return Prop.ForAll(Telomeric3PrimeArbitrary(), t =>
        {
            var (seq, k) = t;
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: 0);
            return (r.TelomereLength3Prime == 6 * k)
                .Label($"3' length {r.TelomereLength3Prime} != 6·{k} for tract of {k} units");
        });
    }

    /// <summary>
    /// INV-M (5' end): a tract of <c>k</c> pure CCCTAA units at the start measures exactly
    /// <c>6·k</c> bases (mirror of the 3' invariant against the reverse-complement motif).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_FivePrimeLength_EqualsSixTimesRepeatCount()
    {
        return Prop.ForAll(Telomeric5PrimeArbitrary(), t =>
        {
            var (seq, k) = t;
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: 0);
            return (r.TelomereLength5Prime == 6 * k)
                .Label($"5' length {r.TelomereLength5Prime} != 6·{k} for tract of {k} units");
        });
    }

    /// <summary>
    /// INV-M (strict monotonicity): adding one more repeat unit strictly lengthens the measured
    /// 3' tract. This is the literal "more repeats → longer region" claim, checked as a delta of
    /// exactly one repeat unit rather than relying on the 6·k identity alone.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_OneMoreRepeat_StrictlyLongerThreePrime()
    {
        return Prop.ForAll(Gen.Choose(1, 400).ToArbitrary(), k =>
        {
            string seqK = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
            string seqK1 = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k + 1));
            int lenK = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seqK, minTelomereLength: 0).TelomereLength3Prime;
            int lenK1 = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seqK1, minTelomereLength: 0).TelomereLength3Prime;
            return (lenK1 > lenK && lenK1 - lenK == 6)
                .Label($"k={k}: len {lenK} → {lenK1} (expected +6, strictly larger)");
        });
    }

    /// <summary>Exact anchors for the length=6·k law at fixed repeat counts.</summary>
    [TestCase(1, 6, TestName = "3' tract of 1 unit → length 6")]
    [TestCase(100, 600, TestName = "3' tract of 100 units → length 600")]
    [TestCase(200, 1200, TestName = "3' tract of 200 units → length 1200")]
    [Category("Property")]
    public void AnalyzeTelomeres_ThreePrimeLength_KnownRepeatCounts(int k, int expectedLength)
    {
        string seq = new string('A', 1000) + string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: 0);
        Assert.That(r.TelomereLength3Prime, Is.EqualTo(expectedLength),
            $"a 3' tract of {k} TTAGGG units must measure {expectedLength} bases");
    }

    #endregion

    #region CHROM-TELO-001 — P: detected telomere is a repeat tract

    /// <summary>
    /// INV-P (purity): a pure TTAGGG 3' tract yields perfect repeat purity (1.0) — every
    /// counted base matches the motif. Has3PrimeTelomere is true iff the measured length 6·k
    /// reaches the configured minTelomereLength.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_PureThreePrimeTract_PerfectPurityAndThresholdedFlag()
    {
        return Prop.ForAll(Telomeric3PrimeArbitrary(), t =>
        {
            var (seq, k) = t;
            const int minLen = 500;
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: minLen);
            bool purityPerfect = Math.Abs(r.RepeatPurity3Prime - 1.0) < 1e-12;
            bool flag = r.Has3PrimeTelomere == (6 * k >= minLen);
            return (purityPerfect && flag)
                .Label($"k={k}: purity={r.RepeatPurity3Prime}, has3'={r.Has3PrimeTelomere} (len {6 * k} vs {minLen})");
        });
    }

    /// <summary>
    /// INV-P (no false positives): non-telomeric DNA (only A/C, no window can reach 70%
    /// similarity to either motif) yields zero-length tracts and both presence flags false.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_NonTelomericDna_NoTelomereDetected()
    {
        return Prop.ForAll(NonTelomericArbitrary(), seq =>
        {
            if (seq.Length == 0) return true.ToProperty(); // empty is a separate edge case.
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            return (r.TelomereLength3Prime == 0 && r.TelomereLength5Prime == 0
                    && !r.Has3PrimeTelomere && !r.Has5PrimeTelomere)
                .Label($"non-telomeric '{seq}' yielded len5={r.TelomereLength5Prime}, len3={r.TelomereLength3Prime}");
        });
    }

    /// <summary>
    /// INV-P (custom motif, ASM-02): a caller-supplied repeat (Arabidopsis TTTAGGG, 7 bp) is
    /// detected at the 3' end with length = 7·k when passed as <c>telomereRepeat</c>, proving
    /// the motif is parameterized rather than hard-coded.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeTelomeres_CustomArabidopsisMotif_DetectedWithCorrectLength()
    {
        const string motif = "TTTAGGG"; // 7 bp
        const int k = 100;
        // Prefix length a multiple of 7 keeps windows aligned; "AAAAAAA" vs TTTAGGG = 1/7 < 0.7.
        string seq = new string('A', 7 * 10) + string.Concat(Enumerable.Repeat(motif, k));
        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, telomereRepeat: motif, minTelomereLength: 0);
        Assert.Multiple(() =>
        {
            Assert.That(r.TelomereLength3Prime, Is.EqualTo(7 * k), "7-bp motif tract must measure 7·k bases");
            Assert.That(r.RepeatPurity3Prime, Is.EqualTo(1.0).Within(1e-12), "pure custom tract → purity 1.0");
        });
    }

    #endregion

    #region CHROM-TELO-001 — R: positions/lengths valid (INV-01, INV-02)

    /// <summary>
    /// INV-01/INV-02 (R): for ANY DNA input, both measured lengths are ≥ 0, never exceed the
    /// scanned window <c>min(searchLength, |sequence|)</c>, are exact multiples of the repeat-unit
    /// length (only complete units counted), and both purities lie in [0, 1]. When a tract is
    /// non-empty its purity is ≥ 0.7 (every accepted window is ≥ 70% similar).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_LengthsAndPurities_AlwaysValid()
    {
        return Prop.ForAll(AnyDnaArbitrary(), seq =>
        {
            const int searchLength = 10000;
            const int unit = 6; // default TTAGGG
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, searchLength: searchLength);
            int cap = Math.Min(searchLength, seq.Length);

            bool nonNeg = r.TelomereLength5Prime >= 0 && r.TelomereLength3Prime >= 0;          // INV-01
            bool capped = r.TelomereLength5Prime <= cap && r.TelomereLength3Prime <= cap;
            bool multiples = r.TelomereLength5Prime % unit == 0 && r.TelomereLength3Prime % unit == 0;
            bool purityRange = r.RepeatPurity5Prime is >= 0.0 and <= 1.0                        // INV-02
                            && r.RepeatPurity3Prime is >= 0.0 and <= 1.0;
            bool accepted5 = r.TelomereLength5Prime == 0 || r.RepeatPurity5Prime >= 0.7 - 1e-12;
            bool accepted3 = r.TelomereLength3Prime == 0 || r.RepeatPurity3Prime >= 0.7 - 1e-12;

            return (nonNeg && capped && multiples && purityRange && accepted5 && accepted3)
                .Label($"invalid result for '{seq}': len5={r.TelomereLength5Prime}, len3={r.TelomereLength3Prime}, " +
                       $"pur5={r.RepeatPurity5Prime}, pur3={r.RepeatPurity3Prime}, cap={cap}");
        });
    }

    /// <summary>
    /// INV-03 (derived flags): Has5PrimeTelomere and Has3PrimeTelomere are exactly the
    /// comparison of the measured length to <c>minTelomereLength</c>, recomputed independently
    /// across random thresholds on a pure tract whose length is known.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_PresenceFlags_DerivedFromMinLength()
    {
        var arb = (from k in Gen.Choose(1, 400)
                   from minLen in Gen.Choose(0, 3000)
                   select (k, minLen)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            var (k, minLen) = t;
            string seq = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: minLen);
            bool expected3 = r.TelomereLength3Prime >= minLen;
            bool expected5 = r.TelomereLength5Prime >= minLen;
            return (r.Has3PrimeTelomere == expected3 && r.Has5PrimeTelomere == expected5)
                .Label($"k={k}, min={minLen}: has3'={r.Has3PrimeTelomere} (len {r.TelomereLength3Prime}), " +
                       $"has5'={r.Has5PrimeTelomere} (len {r.TelomereLength5Prime})");
        });
    }

    #endregion

    #region CHROM-TELO-001 — IsCriticallyShort rule & short-sequence edge cases

    /// <summary>
    /// IsCriticallyShort rule: independently recomputed as
    /// <c>(Has5 &amp;&amp; len5 &lt; criticalLength) || (Has3 &amp;&amp; len3 &lt; criticalLength)</c>
    /// over random thresholds on a 3' tract. A detected-but-short telomere is critical; a long or
    /// absent one is not.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_IsCriticallyShort_MatchesRule()
    {
        var arb = (from k in Gen.Choose(1, 600)
                   from minLen in Gen.Choose(0, 1500)
                   from critical in Gen.Choose(1, 4000)
                   select (k, minLen, critical)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            var (k, minLen, critical) = t;
            string seq = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: minLen, criticalLength: critical);
            bool expected = (r.Has5PrimeTelomere && r.TelomereLength5Prime < critical)
                         || (r.Has3PrimeTelomere && r.TelomereLength3Prime < critical);
            return (r.IsCriticallyShort == expected)
                .Label($"k={k}, min={minLen}, crit={critical}: critical={r.IsCriticallyShort}, expected={expected}");
        });
    }

    /// <summary>
    /// Edge case (M4 / §6.1): empty or null sequence returns an all-zero result with
    /// IsCriticallyShort = true (special case).
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeTelomeres_EmptyOrNull_AllZeroAndCriticallyShort()
    {
        Assert.Multiple(() =>
        {
            foreach (var seq in new[] { "", null })
            {
                var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq!);
                Assert.That(r.TelomereLength5Prime, Is.Zero);
                Assert.That(r.TelomereLength3Prime, Is.Zero);
                Assert.That(r.Has5PrimeTelomere, Is.False);
                Assert.That(r.Has3PrimeTelomere, Is.False);
                Assert.That(r.RepeatPurity5Prime, Is.Zero);
                Assert.That(r.RepeatPurity3Prime, Is.Zero);
                Assert.That(r.IsCriticallyShort, Is.True, "empty input is a special-case critically-short result");
            }
        });
    }

    /// <summary>
    /// Edge case (§6.1): a sequence shorter than the repeat unit cannot contain a complete
    /// window, so both measured lengths are 0 for any sub-6-base input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_ShorterThanRepeatUnit_ZeroLengths()
    {
        var shortGen = (from len in Gen.Choose(1, 5)
                        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                        select new string(chars)).ToArbitrary();
        return Prop.ForAll(shortGen, seq =>
        {
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            return (r.TelomereLength5Prime == 0 && r.TelomereLength3Prime == 0)
                .Label($"sub-unit '{seq}' gave len5={r.TelomereLength5Prime}, len3={r.TelomereLength3Prime}");
        });
    }

    #endregion

    #region CHROM-TELO-001 — INV-04: T/S ratio helper

    /// <summary>
    /// INV-04: <see cref="ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio"/> equals the
    /// independent proportional oracle <c>referenceLength · tsRatio / referenceRatio</c> within
    /// 1e-9, and is non-negative for non-negative inputs (Cawthon 2002).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateTelomereLengthFromTSRatio_MatchesProportionalOracle()
    {
        return Prop.ForAll(TsRatioArbitrary(), t =>
        {
            var (ts, rr, rl) = t;
            double actual = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(ts, rr, rl);
            double expected = ExpectedTsLength(ts, rr, rl);
            return (Math.Abs(actual - expected) < 1e-9 && actual >= 0.0)
                .Label($"T/S({ts},{rr},{rl}) = {actual}, expected {expected}");
        });
    }

    /// <summary>Documented T/S anchors (TestSpec table): exact base-pair estimates.</summary>
    [TestCase(1.0, 1.0, 7000, 7000, TestName = "T/S=1 → reference length 7000")]
    [TestCase(1.5, 1.0, 7000, 10500, TestName = "T/S=1.5 → 10500")]
    [TestCase(0.5, 1.0, 7000, 3500, TestName = "T/S=0.5 → 3500")]
    [TestCase(2.0, 1.0, 7000, 14000, TestName = "T/S=2.0 → 14000")]
    [TestCase(0.0, 1.0, 7000, 0, TestName = "T/S=0 → 0")]
    [Category("Property")]
    public void EstimateTelomereLengthFromTSRatio_CanonicalValues(
        double ts, double refRatio, double refLen, double expected)
    {
        double actual = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(ts, refRatio, refLen);
        Assert.That(actual, Is.EqualTo(expected).Within(1e-9));
    }

    #endregion

    #region CHROM-TELO-001 — D: determinism

    /// <summary>
    /// INV-D: identical inputs produce an identical <see cref="ChromosomeAnalyzer.TelomereResult"/>
    /// (all eight fields), checked via record-struct value equality across arbitrary DNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_IsDeterministic()
    {
        return Prop.ForAll(AnyDnaArbitrary(), seq =>
        {
            var a = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            var b = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            return a.Equals(b).Label($"non-deterministic TelomereResult for '{seq}'");
        });
    }

    #endregion
}
