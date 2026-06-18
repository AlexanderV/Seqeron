using System.Text;
using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for primer/probe design: melting temperature, hairpin, dimer, probe validation.
///
/// Test Units: PRIMER-TM-001, PRIMER-DESIGN-001, PRIMER-STRUCT-001, PROBE-DESIGN-001, PROBE-VALID-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class PrimerProbeProperties
{
    #region Generators &amp; Theory Oracle (PRIMER-TM-001)

    /// <summary>
    /// Generates random valid DNA primers (only A/C/G/T) within an explicit length window.
    /// </summary>
    private static Gen<string> ValidPrimerGen(int minLen, int maxLen) =>
        from len in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
        select new string(chars);

    /// <summary>Wraps <see cref="ValidPrimerGen"/> as an arbitrary.</summary>
    private static Arbitrary<string> ValidPrimerArbitrary(int minLen, int maxLen) =>
        ValidPrimerGen(minLen, maxLen).ToArbitrary();

    /// <summary>
    /// Generates arbitrary character soup: valid DNA bases mixed with case variants,
    /// IUPAC ambiguity codes, RNA uracil, and pure junk — including the empty string.
    /// Used to verify defined behaviors (Section 3 of the TestSpec).
    /// </summary>
    private static Arbitrary<string> MixedSequenceArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T', 'a', 'c', 'g', 't',
                     'N', 'n', 'U', 'u', 'R', 'Y', '-', ' ', 'x', '7')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Generates a length N and two G/C counts (lo ≤ hi) within [0, N], used to build
    /// equal-length sequences that differ only in GC composition.
    /// </summary>
    private static Arbitrary<(int n, int gcLo, int gcHi)> EqualLengthGcPairArbitrary() =>
        (from n in Gen.Choose(1, 60)
         from a in Gen.Choose(0, n)
         from b in Gen.Choose(0, n)
         select (n, Math.Min(a, b), Math.Max(a, b))).ToArbitrary();

    /// <summary>
    /// Generates a valid primer paired with two Na⁺ concentrations (mM) for salt-correction tests.
    /// </summary>
    private static Arbitrary<(string primer, int na1, int na2)> PrimerSaltPairArbitrary() =>
        (from p in ValidPrimerGen(8, 30)
         from na1 in Gen.Choose(1, 2000)
         from na2 in Gen.Choose(1, 2000)
         select (p, na1, na2)).ToArbitrary();

    /// <summary>
    /// Generates a clean ACGT primer together with a "dirty" copy that has non-ACGT
    /// characters sprinkled around every base. Tm must be identical (Section 3.1).
    /// </summary>
    private static Arbitrary<(string clean, string dirty)> CleanDirtyArbitrary() =>
        (from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 1)
         from junk in Gen.Elements('N', 'n', 'U', 'u', 'R', 'Y', '-', ' ', 'x', '9', '.')
                         .ArrayOf().Where(a => a.Length >= 1)
         let clean = new string(chars)
         select (clean, Sprinkle(clean, new string(junk)))).ToArbitrary();

    /// <summary>Builds a string with a junk character inserted before and after every clean base.</summary>
    private static string Sprinkle(string clean, string junk)
    {
        var sb = new StringBuilder();
        sb.Append(junk[0]);
        for (int i = 0; i < clean.Length; i++)
        {
            sb.Append(clean[i]);
            sb.Append(junk[i % junk.Length]);
        }
        return sb.ToString();
    }

    /// <summary>Counts valid DNA bases (case-insensitive), partitioned into A/T and G/C.</summary>
    private static (int at, int gc) CountAtGc(string s)
    {
        int at = 0, gc = 0;
        foreach (char ch in s.ToUpperInvariant())
        {
            if (ch is 'A' or 'T') at++;
            else if (ch is 'G' or 'C') gc++;
        }
        return (at, gc);
    }

    /// <summary>
    /// Independent theory oracle for Tm, derived directly from the published equations
    /// (Thein &amp; Wallace 1986; Marmur &amp; Doty 1962) with literal constants — deliberately
    /// NOT routed through <c>ThermoConstants</c>, so a wrong production constant is caught.
    /// </summary>
    private static double ExpectedTm(string s)
    {
        var (at, gc) = CountAtGc(s);
        int n = at + gc;
        if (n == 0) return 0;
        // Wallace rule for short oligos (&lt; 14 valid bases): Tm = 2(A+T) + 4(G+C).
        if (n < 14) return 2 * at + 4 * gc;
        // Marmur-Doty for longer primers: Tm = 64.9 + 41(GC − 16.4)/N, floored at 0.
        return Math.Max(0, 64.9 + 41.0 * (gc - 16.4) / n);
    }

    /// <summary>Upper bound on Tm for any input: the Marmur-Doty asymptote 64.9 + 41 = 105.9°C.</summary>
    private static double TmUpperBound =>
        ThermoConstants.MarmurDotyBase + ThermoConstants.MarmurDotyGcCoefficient;

    #endregion

    #region PRIMER-TM-001 — Formula Conformance (close to theory)

    /// <summary>
    /// INV-Wallace: For any primer with &lt; 14 valid bases, Tm equals the Wallace-rule value
    /// 2(A+T) + 4(G+C) computed independently from base counts.
    /// Evidence: Thein &amp; Wallace (1986). Validates regime selection, base counting, and
    /// the A/T = 2°C, G/C = 4°C contributions across hundreds of random short oligos.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_ShortPrimer_MatchesWallaceRule()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 13), primer =>
        {
            double actual = PrimerDesigner.CalculateMeltingTemperature(primer);
            return (Math.Abs(actual - ExpectedTm(primer)) < 1e-9)
                .Label($"Wallace Tm mismatch for '{primer}': {actual} vs {ExpectedTm(primer)}");
        });
    }

    /// <summary>
    /// INV-MarmurDoty: For any primer with ≥ 14 valid bases, Tm equals
    /// max(0, 64.9 + 41(GC − 16.4)/N) computed independently from base counts.
    /// Evidence: Marmur &amp; Doty (1962). Validates the long-primer regime and GC/length scaling.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_LongPrimer_MatchesMarmurDoty()
    {
        return Prop.ForAll(ValidPrimerArbitrary(14, 40), primer =>
        {
            double actual = PrimerDesigner.CalculateMeltingTemperature(primer);
            return (Math.Abs(actual - ExpectedTm(primer)) < 1e-9)
                .Label($"Marmur-Doty Tm mismatch for '{primer}': {actual} vs {ExpectedTm(primer)}");
        });
    }

    /// <summary>
    /// INV-Threshold: The Wallace ↔ Marmur-Doty switch happens exactly at 14 valid bases.
    /// A 13-mer uses Wallace (integer multiple of 2/4); a 14-mer uses Marmur-Doty.
    /// Evidence: TestSpec Section 2 (threshold &lt; 14). Pins the boundary that the
    /// random-length properties exercise only stochastically.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperature_RegimeThreshold_SwitchesAt14()
    {
        // 13 valid bases → Wallace: 7×A/T (2) + 6×G/C (4) = 14 + 24 = 38.
        double tm13 = PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACGTA");
        // 14 valid bases → Marmur-Doty: 64.9 + 41×(7−16.4)/14 ≈ 37.36 (non-integer ⇒ not Wallace).
        double tm14 = PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACGTAC");

        Assert.Multiple(() =>
        {
            Assert.That(tm13, Is.EqualTo(38.0).Within(1e-9), "13-mer must use the Wallace rule");
            Assert.That(tm14, Is.EqualTo(64.9 + 41.0 * (7 - 16.4) / 14).Within(1e-9),
                "14-mer must use the Marmur-Doty formula");
        });
    }

    #endregion

    #region PRIMER-TM-001 — Evidence-Based Canonical Values (literature anchors)

    /// <summary>
    /// Anchors Tm to absolute, literature-published values so that a self-consistent but
    /// wrong constant (shared between production and oracle) cannot pass unnoticed.
    /// Sources: Thein &amp; Wallace (1986) for short oligos; Marmur &amp; Doty (1962) for ≥14-mers.
    /// </summary>
    [TestCase("ATATATAT", 16.0, TestName = "Wallace: 8×A/T → 16°C")]
    [TestCase("GCGCGCGC", 32.0, TestName = "Wallace: 8×G/C → 32°C")]
    [TestCase("ACGTACGT", 24.0, TestName = "Wallace: 4 A/T + 4 G/C → 24°C")]
    [TestCase("ATATATATATATATATATAT", 31.28, TestName = "Marmur-Doty: 20bp 0% GC → 31.28°C")]
    [TestCase("ACGTACGTACGTACGTACGT", 51.78, TestName = "Marmur-Doty: 20bp 50% GC → 51.78°C")]
    [TestCase("GCGCGCGCGCGCGCGCGCGC", 72.28, TestName = "Marmur-Doty: 20bp 100% GC → 72.28°C")]
    [Category("Property")]
    public void MeltingTemperature_CanonicalLiteratureValues(string primer, double expected)
    {
        double actual = PrimerDesigner.CalculateMeltingTemperature(primer);
        Assert.That(actual, Is.EqualTo(expected).Within(0.01),
            $"Tm for '{primer}' must match the published value {expected}°C");
    }

    #endregion

    #region PRIMER-TM-001 — Range &amp; Validity Invariants

    /// <summary>
    /// INV-1 (R): Tm is always within the physically meaningful band [0, 105.9) for ANY input,
    /// including junk, mixed case, and the empty string. The upper bound is the Marmur-Doty
    /// asymptote (64.9 + 41), which an all-G/C primer approaches but never reaches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_AlwaysInPhysicalRange()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double tm = PrimerDesigner.CalculateMeltingTemperature(seq);
            return (tm >= 0.0 && tm < TmUpperBound)
                .Label($"Tm {tm} outside [0, {TmUpperBound}) for '{seq}'");
        });
    }

    /// <summary>
    /// INV-2 (P): Tm is strictly positive iff the sequence contains at least one valid DNA base,
    /// and exactly 0 otherwise. This separates a correct DNA-aware implementation from one that
    /// would emit a spurious temperature for sequences like "NNNN" or "" (Section 3.1/3.3).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_PositiveIffValidBasePresent()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var (at, gc) = CountAtGc(seq);
            bool hasValidBase = at + gc > 0;
            double tm = PrimerDesigner.CalculateMeltingTemperature(seq);
            return ((tm > 0.0) == hasValidBase)
                .Label($"Tm {tm} vs validBases={at + gc} for '{seq}'");
        });
    }

    #endregion

    #region PRIMER-TM-001 — Composition Monotonicity (the core thermodynamic claim)

    /// <summary>
    /// INV-3 (M): For two sequences of equal valid length, the one with more G/C bases has a
    /// strictly higher Tm; equal G/C content yields an equal Tm. This is the central biological
    /// invariant — G·C pairs (3 hydrogen bonds) stabilize the duplex more than A·T pairs (2 H-bonds),
    /// so higher GC ⇒ higher melting temperature. Equal length keeps both sequences in the same
    /// formula regime, making the comparison rigorous rather than confounded by the 14-base threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_HigherGc_StrictlyHigherTm_AtEqualLength()
    {
        return Prop.ForAll(EqualLengthGcPairArbitrary(), triple =>
        {
            var (n, gcLo, gcHi) = triple;
            string seqLo = new string('G', gcLo) + new string('A', n - gcLo);
            string seqHi = new string('G', gcHi) + new string('A', n - gcHi);
            double tmLo = PrimerDesigner.CalculateMeltingTemperature(seqLo);
            double tmHi = PrimerDesigner.CalculateMeltingTemperature(seqHi);

            bool ok = gcHi > gcLo ? tmHi > tmLo : Math.Abs(tmHi - tmLo) < 1e-9;
            return ok.Label($"n={n}, gc {gcLo}→{tmLo} vs {gcHi}→{tmHi}");
        });
    }

    #endregion

    #region PRIMER-TM-001 — Robustness Invariants (defined behaviors)

    /// <summary>
    /// INV-4 (P): Tm is case-insensitive — uppercase, lowercase, and mixed-case spellings of the
    /// same primer all yield the identical temperature (Section 3.2).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IsCaseInsensitive()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 40), primer =>
        {
            char[] m = primer.ToCharArray();
            for (int i = 0; i < m.Length; i += 2) m[i] = char.ToLowerInvariant(m[i]);
            string mixed = new string(m);

            double up = PrimerDesigner.CalculateMeltingTemperature(primer.ToUpperInvariant());
            double lo = PrimerDesigner.CalculateMeltingTemperature(primer.ToLowerInvariant());
            double mi = PrimerDesigner.CalculateMeltingTemperature(mixed);

            return (Math.Abs(up - lo) < 1e-9 && Math.Abs(up - mi) < 1e-9)
                .Label($"Case sensitivity: up={up}, lo={lo}, mixed={mi} for '{primer}'");
        });
    }

    /// <summary>
    /// INV-5 (P): Non-ACGT characters are ignored — sprinkling N's, U's, gaps, and other junk
    /// around a primer leaves its Tm unchanged, because only valid DNA bases are counted
    /// (Section 3.1). This protects against an implementation that lets ambiguity codes leak
    /// into the temperature.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IgnoresNonAcgtCharacters()
    {
        return Prop.ForAll(CleanDirtyArbitrary(), pair =>
        {
            var (clean, dirty) = pair;
            double tmClean = PrimerDesigner.CalculateMeltingTemperature(clean);
            double tmDirty = PrimerDesigner.CalculateMeltingTemperature(dirty);
            return (Math.Abs(tmClean - tmDirty) < 1e-9)
                .Label($"junk changed Tm: clean '{clean}'={tmClean}, dirty '{dirty}'={tmDirty}");
        });
    }

    /// <summary>
    /// INV-6 (D): Tm is deterministic — repeated evaluation of the same input is bit-identical.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double a = PrimerDesigner.CalculateMeltingTemperature(seq);
            double b = PrimerDesigner.CalculateMeltingTemperature(seq);
            return a.Equals(b).Label($"non-deterministic Tm for '{seq}': {a} vs {b}");
        });
    }

    #endregion

    #region PRIMER-TM-001 — Salt Correction Invariants

    /// <summary>
    /// INV-7 (P): Salt-corrected Tm equals the base Tm plus the published salt term
    /// 16.6 × log₁₀([Na⁺]/1000), rounded to one decimal. The correction is computed here
    /// directly from the Owczarzy (2004) formula — independent of production code.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_IsAdditiveCorrection()
    {
        return Prop.ForAll(PrimerSaltPairArbitrary(), t =>
        {
            var (primer, na, _) = t;
            double baseTm = PrimerDesigner.CalculateMeltingTemperature(primer);
            double expected = Math.Round(baseTm + 16.6 * Math.Log10(na / 1000.0), 1);
            double actual = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, na);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"salt Tm mismatch for '{primer}' @ {na}mM: {actual} vs {expected}");
        });
    }

    /// <summary>
    /// INV-8 (M): Higher salt stabilizes the duplex — increasing [Na⁺] never decreases the
    /// corrected Tm. Evidence: Owczarzy et al. (2004); higher cation concentration shields the
    /// phosphate backbone and raises melting temperature.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_IsMonotonicInSalt()
    {
        return Prop.ForAll(PrimerSaltPairArbitrary(), t =>
        {
            var (primer, na1, na2) = t;
            int lo = Math.Min(na1, na2), hi = Math.Max(na1, na2);
            double tmLo = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, lo);
            double tmHi = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, hi);
            return (tmHi >= tmLo)
                .Label($"salt monotonicity broken for '{primer}': {lo}mM→{tmLo}, {hi}mM→{tmHi}");
        });
    }

    /// <summary>
    /// INV-9 (P): At the 1 M (1000 mM) reference concentration the salt term is exactly 0
    /// (log₁₀(1) = 0), so the corrected Tm collapses to the base Tm (rounded). This pins the
    /// zero-crossing of the correction curve.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_ReferenceConcentration_EqualsBaseTm()
    {
        return Prop.ForAll(ValidPrimerArbitrary(8, 30), primer =>
        {
            double expected = Math.Round(PrimerDesigner.CalculateMeltingTemperature(primer), 1);
            double actual = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 1000);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"@1000mM '{primer}': {actual} vs base {expected}");
        });
    }

    /// <summary>
    /// INV-10 (M): Standard PCR salt (50 mM) is below the 1 M reference, so it destabilizes the
    /// duplex and yields a Tm strictly lower than the uncorrected base Tm for any real primer.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_StandardPcrSalt_LowersTm()
    {
        return Prop.ForAll(ValidPrimerArbitrary(8, 30), primer =>
        {
            double baseTm = PrimerDesigner.CalculateMeltingTemperature(primer);
            double salted = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 50);
            return (salted < baseTm)
                .Label($"50mM did not lower Tm for '{primer}': base={baseTm}, salted={salted}");
        });
    }

    /// <summary>
    /// INV-11 (Edge): Salt correction is undefined for an empty/null primer — both return 0,
    /// regardless of the requested concentration (Section 3.3).
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperatureWithSalt_EmptyOrNull_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureWithSalt("", 50), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureWithSalt(null!, 50), Is.EqualTo(0.0));
        });
    }

    #endregion

    // -- PRIMER-STRUCT-001 --

    /// <summary>
    /// Self-complementary primer has hairpin potential.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasHairpinPotential_SelfComplementary_ReturnsTrue()
    {
        // GCGC...GCGC is self-complementary
        string selfComp = "GCGCGCTTTTGCGCGC";
        bool result = PrimerDesigner.HasHairpinPotential(selfComp, minStemLength: 4, minLoopLength: 3);
        Assert.That(result, Is.True, "Self-complementary primer should have hairpin potential");
    }

    /// <summary>
    /// Primer dimer detection: identical forward/reverse with self-complementary 3' ends.
    /// Source: Wikipedia Primer-dimer — A₂₀ vs A₂₀, revcomp = T₂₀, A-T complementary.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasPrimerDimer_SelfComplementaryPair_ReturnsTrue()
    {
        // AAAAAAAA vs AAAAAAAA: revcomp(AAAAAAAA) = TTTTTTTT
        // 3' end (A) vs 5' revcomp (T) = A-T complementary
        bool result = PrimerDesigner.HasPrimerDimer("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", minComplementarity: 4);
        Assert.That(result, Is.True, "Poly-A self-dimer: A-T complementarity with revcomp");
    }

    /// <summary>
    /// 3' stability score is finite and negative for valid sequences.
    /// Source: SantaLucia (1998) — all NN ΔG values are negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Calculate3PrimeStability_ValidSequence_IsNegativeAndFinite()
    {
        double stability = PrimerDesigner.Calculate3PrimeStability("ACGTACGTACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(double.IsFinite(stability), Is.True);
            Assert.That(stability, Is.LessThan(0.0));
        });
    }

    #region Generators &amp; Theory Oracle (PRIMER-DESIGN-001)

    /// <summary>
    /// Generates a random ACGT template at least <paramref name="minLen"/> bases long.
    /// Only valid DNA bases are produced so <see cref="DnaSequence"/> construction never throws,
    /// and the sequence is realistic enough to yield candidates across the filter windows.
    /// </summary>
    private static Gen<string> TemplateGen(int minLen, int maxLen) =>
        from len in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
        select new string(chars);

    /// <summary>
    /// Generates a randomized parameter window (the design "config knobs") together with a
    /// template long enough to host both flanking search regions and a target between them.
    /// The length/GC/Tm windows are deliberately varied so the invariants are exercised across
    /// configurations rather than a single hard-coded setting (cf. ADVANCED_TESTING_CHECKLIST
    /// "multiple configuration knobs": length range + GC range + Tm range).
    /// </summary>
    private static Arbitrary<(string template, int targetStart, int targetEnd, PrimerParameters param)>
        DesignScenarioArbitrary() =>
        (from minLen in Gen.Choose(15, 20)
         from extra in Gen.Choose(0, 8)
         let maxLen = minLen + extra
         from minGc in Gen.Choose(20, 45)
         from gcSpan in Gen.Choose(15, 50)
         let maxGc = Math.Min(100, minGc + gcSpan)
         from minTm in Gen.Choose(40, 58)
         from tmSpan in Gen.Choose(8, 30)
         let maxTm = minTm + tmSpan
         // Template must hold: forward flank (>= maxLen) + target (>= 1) + reverse flank (>= maxLen).
         from fwdFlank in Gen.Choose(maxLen, maxLen + 40)
         from targetSpan in Gen.Choose(1, 30)
         from revFlank in Gen.Choose(maxLen, maxLen + 40)
         let total = fwdFlank + targetSpan + revFlank
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(total)
         let param = new PrimerParameters(
             MinLength: minLen, MaxLength: maxLen, OptimalLength: (minLen + maxLen) / 2,
             MinGcContent: minGc, MaxGcContent: maxGc,
             MinTm: minTm, MaxTm: maxTm, OptimalTm: (minTm + maxTm) / 2.0,
             MaxHomopolymer: 4, MaxDinucleotideRepeats: 4,
             Avoid3PrimeGC: false, Check3PrimeStability: true)
         select (new string(chars), fwdFlank, fwdFlank + targetSpan, param)).ToArbitrary();

    /// <summary>
    /// Generates a candidate-enumeration scenario: a template, a region inside it, an orientation,
    /// and a randomized parameter window. Used to test <see cref="PrimerDesigner.GeneratePrimerCandidates"/>.
    /// </summary>
    private static Arbitrary<(string template, int regionStart, int regionEnd, bool forward, PrimerParameters param)>
        CandidateScenarioArbitrary() =>
        (from minLen in Gen.Choose(15, 20)
         from extra in Gen.Choose(0, 8)
         let maxLen = minLen + extra
         from minGc in Gen.Choose(20, 45)
         from gcSpan in Gen.Choose(15, 50)
         let maxGc = Math.Min(100, minGc + gcSpan)
         from minTm in Gen.Choose(40, 58)
         from tmSpan in Gen.Choose(8, 30)
         let maxTm = minTm + tmSpan
         from regionLen in Gen.Choose(maxLen, maxLen + 60)
         from forward in Gen.Elements(true, false)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(regionLen)
         let param = new PrimerParameters(
             MinLength: minLen, MaxLength: maxLen, OptimalLength: (minLen + maxLen) / 2,
             MinGcContent: minGc, MaxGcContent: maxGc,
             MinTm: minTm, MaxTm: maxTm, OptimalTm: (minTm + maxTm) / 2.0,
             MaxHomopolymer: 4, MaxDinucleotideRepeats: 4,
             Avoid3PrimeGC: false, Check3PrimeStability: true)
         select (new string(chars), 0, regionLen, forward, param)).ToArbitrary();

    /// <summary>Independent GC% oracle: 100 × (#G + #C) / length, computed from raw base counts.</summary>
    private static double ExpectedGcPercent(string seq)
    {
        if (seq.Length == 0) return 0;
        var (_, gc) = CountAtGc(seq);
        return 100.0 * gc / seq.Length;
    }

    /// <summary>
    /// Independent re-implementation of the pair-validity predicate (INV-02), derived from the
    /// doc contract rather than reused production code: a pair is valid iff |Tm_f − Tm_r| ≤ 5
    /// AND the two primers do not form a primer-dimer.
    /// </summary>
    private static bool ExpectedPairValid(PrimerCandidate fwd, PrimerCandidate rev) =>
        Math.Abs(fwd.MeltingTemperature - rev.MeltingTemperature) <= 5.0
        && !PrimerDesigner.HasPrimerDimer(fwd.Sequence, rev.Sequence);

    #endregion

    #region PRIMER-DESIGN-001 — Candidate Filter Invariants (R: length, P: GC, R: Tm)

    /// <summary>
    /// INV-R-len: Every VALID candidate returned by <see cref="PrimerDesigner.GeneratePrimerCandidates"/>
    /// has Length within the configured [MinLength, MaxLength] window, and an out-of-window length is
    /// never accepted. Exercised across randomized templates AND randomized length windows so the
    /// bound is a true config-driven invariant, not a single-template coincidence.
    /// Source: Primer_Design.md §3 contract (length 18-25 default; window from parameters).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_ValidLength_InConfiguredWindow()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var candidates = PrimerDesigner
                .GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param)
                .ToList();

            bool allValidInWindow = candidates
                .Where(c => c.IsValid)
                .All(c => c.Length >= s.param.MinLength && c.Length <= s.param.MaxLength);

            return allValidInWindow.Label(
                $"valid candidate length outside [{s.param.MinLength},{s.param.MaxLength}] " +
                $"in region len {s.template.Length}");
        });
    }

    /// <summary>
    /// INV-P-GC: Every VALID candidate has GcContent (a percentage 0-100 in source) within the
    /// configured [MinGcContent, MaxGcContent] window; a candidate whose GC% falls outside the window
    /// is never marked valid. GcContent is recomputed independently from base counts (ExpectedGcPercent)
    /// to confirm the production value used for filtering is itself correct.
    /// Source: Primer_Design.md §3 / §4.2 (GC 40-60 default; window from parameters).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_ValidGcContent_InConfiguredWindow()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var candidates = PrimerDesigner
                .GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param)
                .ToList();

            foreach (var c in candidates)
            {
                // GcContent is rounded to 1 decimal in source; oracle within that tolerance.
                if (Math.Abs(c.GcContent - ExpectedGcPercent(c.Sequence)) > 0.06)
                    return false.Label($"GcContent {c.GcContent} != oracle {ExpectedGcPercent(c.Sequence)} for {c.Sequence}");

                if (c.IsValid &&
                    (c.GcContent < s.param.MinGcContent || c.GcContent > s.param.MaxGcContent))
                    return false.Label(
                        $"valid candidate GC {c.GcContent}% outside [{s.param.MinGcContent},{s.param.MaxGcContent}]");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-R-Tm: Every VALID candidate has MeltingTemperature within the configured [MinTm, MaxTm]
    /// window; a candidate whose Tm is outside the window is never marked valid.
    /// Source: Primer_Design.md §3 (Tm 57-63 default; window from parameters).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_ValidTm_InConfiguredWindow()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var candidates = PrimerDesigner
                .GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param)
                .ToList();

            bool allValidInWindow = candidates
                .Where(c => c.IsValid)
                .All(c => c.MeltingTemperature >= s.param.MinTm && c.MeltingTemperature <= s.param.MaxTm);

            return allValidInWindow.Label(
                $"valid candidate Tm outside [{s.param.MinTm},{s.param.MaxTm}]");
        });
    }

    /// <summary>
    /// INV-R-len (negation): an invalid-length primer is never marked valid by
    /// <see cref="PrimerDesigner.EvaluatePrimer"/>. Builds primers one base shorter than MinLength
    /// and one longer than MaxLength and asserts IsValid is false in both cases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EvaluatePrimer_LengthOutsideWindow_IsNeverValid()
    {
        var arb = (from minLen in Gen.Choose(16, 22)
                   from extra in Gen.Choose(2, 8)
                   let maxLen = minLen + extra
                   from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(maxLen + 6)
                   select (new string(chars), minLen, maxLen)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var (seq, minLen, maxLen) = t;
            var param = new PrimerParameters(minLen, maxLen, (minLen + maxLen) / 2,
                0, 100, 0, 200, 60, 99, 99, false, false);

            string tooShort = seq.Substring(0, minLen - 1);
            string tooLong = seq.Substring(0, maxLen + 1);

            bool shortInvalid = !PrimerDesigner.EvaluatePrimer(tooShort, 0, true, param).IsValid;
            bool longInvalid = !PrimerDesigner.EvaluatePrimer(tooLong, 0, true, param).IsValid;

            return (shortInvalid && longInvalid)
                .Label($"len-window violation accepted: short={shortInvalid}, long={longInvalid} for [{minLen},{maxLen}]");
        });
    }

    #endregion

    #region PRIMER-DESIGN-001 — Pair Invariants (INV-02 validity, INV-03 product size)

    /// <summary>
    /// INV-03 (product size): when <see cref="PrimerDesigner.DesignPrimers"/> returns both primers,
    /// ProductSize equals reverse.Position + reverse.Sequence.Length − forward.Position, recomputed
    /// independently from the returned candidates. Source: Primer_Design.md §2.4 INV-03.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_ProductSize_MatchesIndependentFormula()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var result = PrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            if (result.Forward is null || result.Reverse is null)
                return true.ToProperty(); // INV-03 only applies when both primers exist.

            int expected = result.Reverse.Position + result.Reverse.Sequence.Length - result.Forward.Position;
            return (result.ProductSize == expected)
                .Label($"ProductSize {result.ProductSize} != {expected}");
        });
    }

    /// <summary>
    /// INV-02 (pair validity): result.IsValid ⟺ (|Tm_f − Tm_r| ≤ 5 AND not HasPrimerDimer(fwd,rev)),
    /// with the predicate recomputed independently (ExpectedPairValid) from the returned candidates.
    /// Source: Primer_Design.md §2.4 INV-02.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_PairValidity_MatchesIndependentPredicate()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var result = PrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            if (result.Forward is null || result.Reverse is null)
                return result.IsValid == false ? true.ToProperty()  // INV-01: missing side ⇒ invalid.
                    : false.ToProperty().Label("missing primer but IsValid==true");

            bool expected = ExpectedPairValid(result.Forward, result.Reverse);
            return (result.IsValid == expected)
                .Label($"IsValid {result.IsValid} != predicate {expected} " +
                       $"(Tm {result.Forward.MeltingTemperature}/{result.Reverse.MeltingTemperature})");
        });
    }

    /// <summary>
    /// INV-R/P/R on the chosen pair: whenever DesignPrimers selects primers, EACH selected primer
    /// satisfies the configured length, GC, and Tm windows (selection only ever draws from valid
    /// candidates). Ties the three filter invariants to the public design entry point.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_SelectedPrimers_SatisfyAllWindows()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var result = PrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            foreach (var c in new[] { result.Forward, result.Reverse })
            {
                if (c is null) continue;
                bool ok = c.Length >= s.param.MinLength && c.Length <= s.param.MaxLength
                          && c.GcContent >= s.param.MinGcContent && c.GcContent <= s.param.MaxGcContent
                          && c.MeltingTemperature >= s.param.MinTm && c.MeltingTemperature <= s.param.MaxTm;
                if (!ok)
                    return false.Label(
                        $"selected primer len={c.Length} gc={c.GcContent} tm={c.MeltingTemperature} " +
                        $"violates windows len[{s.param.MinLength},{s.param.MaxLength}] " +
                        $"gc[{s.param.MinGcContent},{s.param.MaxGcContent}] tm[{s.param.MinTm},{s.param.MaxTm}]");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region PRIMER-DESIGN-001 — Determinism (D)

    /// <summary>
    /// INV-D (DesignPrimers): identical inputs yield identical pair results — same validity,
    /// product size, message, and selected primer sequences/positions.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_IsDeterministic()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var a = PrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);
            var b = PrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            bool same = a.IsValid == b.IsValid
                        && a.ProductSize == b.ProductSize
                        && a.Message == b.Message
                        && a.Forward?.Sequence == b.Forward?.Sequence
                        && a.Forward?.Position == b.Forward?.Position
                        && a.Reverse?.Sequence == b.Reverse?.Sequence
                        && a.Reverse?.Position == b.Reverse?.Position;

            return same.Label("DesignPrimers not deterministic for identical inputs");
        });
    }

    /// <summary>
    /// INV-D (GeneratePrimerCandidates): the candidate list is identical across repeated calls —
    /// same count, sequences, positions, validity flags, and scores.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_IsDeterministic()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var a = PrimerDesigner.GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param).ToList();
            var b = PrimerDesigner.GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param).ToList();

            if (a.Count != b.Count)
                return false.Label($"candidate count differs: {a.Count} vs {b.Count}");

            for (int i = 0; i < a.Count; i++)
            {
                if (a[i].Sequence != b[i].Sequence || a[i].Position != b[i].Position
                    || a[i].IsValid != b[i].IsValid || a[i].Score != b[i].Score)
                    return false.Label($"candidate {i} differs between runs");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region PRIMER-DESIGN-001 — Edge Cases &amp; Evidence-Based Anchors

    /// <summary>
    /// Edge (§6.1): an invalid target region throws ArgumentException — targetStart &lt; 0,
    /// targetEnd ≥ template.Length, or targetStart ≥ targetEnd. Source: explicit guard in DesignPrimers.
    /// </summary>
    [TestCase(-1, 50, TestName = "targetStart < 0 throws")]
    [TestCase(0, 500, TestName = "targetEnd >= length throws")]
    [TestCase(60, 60, TestName = "targetStart == targetEnd throws")]
    [TestCase(80, 40, TestName = "targetStart > targetEnd throws")]
    [Category("Property")]
    public void DesignPrimers_InvalidTargetRegion_ThrowsArgumentException(int start, int end)
    {
        var dna = new DnaSequence(string.Concat(Enumerable.Repeat("ACGT", 75))); // length 300
        Assert.That(() => PrimerDesigner.DesignPrimers(dna, start, end),
            NUnit.Framework.Throws.ArgumentException);
    }

    /// <summary>
    /// Evidence anchor (GcContent + Tm): a hand-constructed 22-mer with a known composition.
    /// "CGGTTCACTACGTCCGTTCTGG" has 13 G/C of 22 bases ⇒ GC% = 1300/22 = 59.09 → 59.1 (rounded),
    /// and (≥14 bases) Marmur-Doty Tm = 64.9 + 41×(13 − 16.4)/22 = 58.56 → 58.6 (rounded).
    /// Both fall inside the default windows (GC 40-60, Tm 57-63), pinning EvaluatePrimer's
    /// per-primer metrics to literally computed values independent of production constants.
    /// </summary>
    [Test]
    [Category("Property")]
    public void EvaluatePrimer_KnownMer_HasComputedGcAndTm()
    {
        var c = PrimerDesigner.EvaluatePrimer("CGGTTCACTACGTCCGTTCTGG", 0, true);
        Assert.Multiple(() =>
        {
            Assert.That(c.Length, Is.EqualTo(22));
            Assert.That(c.GcContent, Is.EqualTo(59.1).Within(0.05),
                "GC% = 100×13/22 = 59.09 → 59.1");
            Assert.That(c.MeltingTemperature, Is.EqualTo(58.6).Within(0.05),
                "Marmur-Doty Tm = 64.9 + 41×(13−16.4)/22 = 58.56 → 58.6");
        });
    }

    /// <summary>
    /// Evidence anchor (valid pair + product size): a hand-constructed 48 bp template engineered so
    /// the design produces a known-valid pair. The forward flank is the verified valid 22-mer
    /// "CGGTTCACTACGTCCGTTCTGG" (positions 0-21); a 4 bp target ("AAAA") occupies 22-25; the reverse
    /// flank (positions 26-47) is the reverse complement of the verified valid 22-mer
    /// "CAAGCCGGGGCTAATCCGTCAT", so reverse-complementing it back recovers that primer. Both Tm = 58.6
    /// (|ΔTm| = 0 ≤ 5) and they are not a dimer ⇒ a valid pair. ProductSize = reverse.Position(26) +
    /// reverse.Length(22) − forward.Position(0) = 48. This pins INV-03 and INV-02 to literal values.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignPrimers_KnownTemplate_ProductSizeEqualsSpan()
    {
        const string forwardPrimer = "CGGTTCACTACGTCCGTTCTGG";
        const string reversePrimer = "CAAGCCGGGGCTAATCCGTCAT";
        string reverseFlank = new DnaSequence(reversePrimer).ReverseComplement().Sequence;
        string template = forwardPrimer + "AAAA" + reverseFlank; // length 48

        var dna = new DnaSequence(template);
        var result = PrimerDesigner.DesignPrimers(dna, targetStart: 22, targetEnd: 26);

        Assert.That(result.Forward, Is.Not.Null, "expected a forward primer on this template");
        Assert.That(result.Reverse, Is.Not.Null, "expected a reverse primer on this template");

        int expected = result.Reverse!.Position + result.Reverse.Sequence.Length - result.Forward!.Position;
        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.True, "engineered pair must be compatible");
            Assert.That(result.Forward!.Sequence, Is.EqualTo(forwardPrimer));
            Assert.That(result.Reverse!.Sequence, Is.EqualTo(reversePrimer));
            Assert.That(result.ProductSize, Is.EqualTo(expected));
            Assert.That(result.ProductSize, Is.EqualTo(48),
                "ProductSize = 26 + 22 − 0 = 48");
        });
    }

    #endregion

    // -- PROBE-DESIGN-001 --

    /// <summary>
    /// Designed probes have GC content within valid range.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_GcContent_InRange()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var p in probes)
            Assert.That(p.GcContent, Is.InRange(0.0, 1.0));
    }

    /// <summary>
    /// Designed probes have finite melting temperature.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_Tm_IsFinite()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var p in probes)
            Assert.That(double.IsFinite(p.Tm), Is.True);
    }

    /// <summary>
    /// Probe score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_Score_InRange()
    {
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 10));
        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var p in probes)
            Assert.That(p.Score, Is.InRange(0.0, 1.0));
    }

    // -- PROBE-VALID-001 --

    /// <summary>
    /// ValidateProbe returns a valid ProbeValidation with specificity in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_Specificity_InRange()
    {
        string probe = "ACGTACGTACGTACGTACGTACGT";
        var refs = new[] { "ACGTACGTACGTACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTTTTTT" };
        var validation = ProbeDesigner.ValidateProbe(probe, refs);

        Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0));
    }

    /// <summary>
    /// ValidateProbe off-target hits is non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_OffTargetHits_NonNegative()
    {
        string probe = "ACGTACGTACGTACGTACGTACGT";
        var refs = new[] { "TTTTTTTTTTTTTTTTTTTTTTTTTTTT" };
        var validation = ProbeDesigner.ValidateProbe(probe, refs);

        Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0));
    }
}
