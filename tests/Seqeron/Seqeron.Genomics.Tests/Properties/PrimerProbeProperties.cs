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

    #region Generators &amp; Theory Oracle (PRIMER-STRUCT-001)

    /// <summary>
    /// Independent SantaLucia (1998) nearest-neighbor ΔG°37 oracle for the 3' terminal 5-mer.
    /// Constants are transcribed LITERALLY from Primer_Structure_Analysis.md §2.2/§4.2 and the
    /// SantaLucia (1998) unified table — deliberately NOT routed through production, so a wrong
    /// production constant is caught. ΔG = Σ(four NN dinucleotide steps) + two terminal initiation
    /// terms (+0.98 for a terminal G·C, +1.03 for a terminal A·T). Mirrors the doc anchors
    /// GCGCG = −6.86 (most stable) and TATAT = −0.86 (least stable).
    /// </summary>
    private static double ExpectedThreePrimeStability(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 5)
            return 0;

        // Unified NN ΔG°37 (kcal/mol), SantaLucia (1998) — written here independently.
        var nn = new Dictionary<string, double>
        {
            ["AA"] = -1.00, ["TT"] = -1.00,
            ["AT"] = -0.88,
            ["TA"] = -0.58,
            ["CA"] = -1.45, ["TG"] = -1.45,
            ["GT"] = -1.44, ["AC"] = -1.44,
            ["CT"] = -1.28, ["AG"] = -1.28,
            ["GA"] = -1.30, ["TC"] = -1.30,
            ["CG"] = -2.17,
            ["GC"] = -2.24,
            ["GG"] = -1.84, ["CC"] = -1.84,
        };

        string last5 = sequence.ToUpperInvariant()[^5..];
        double dg = 0;
        for (int i = 0; i < last5.Length - 1; i++)
            if (nn.TryGetValue(last5.Substring(i, 2), out double step))
                dg += step;

        // Terminal initiation: +0.98 kcal/mol per terminal G·C, +1.03 per terminal A·T.
        dg += last5[0] is 'G' or 'C' ? 0.98 : 1.03;
        dg += last5[^1] is 'G' or 'C' ? 0.98 : 1.03;
        return dg;
    }

    /// <summary>Generates valid ACGT sequences with length ≥ 5 (the 3'-stability operating regime).</summary>
    private static Arbitrary<string> StabilitySequenceArbitrary() =>
        ValidPrimerGen(5, 40).ToArbitrary();

    /// <summary>
    /// Generates a valid ACGT sequence together with hairpin stem/loop knobs. The sequence length
    /// is drawn small enough to frequently fall below the structural minimum 2·stem + loop, so
    /// INV-01 is exercised on both sides of the threshold.
    /// </summary>
    private static Arbitrary<(string seq, int stem, int loop)> HairpinKnobArbitrary() =>
        (from stem in Gen.Choose(2, 6)
         from loop in Gen.Choose(1, 6)
         from len in Gen.Choose(1, 40)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
         select (new string(chars), stem, loop)).ToArbitrary();

    /// <summary>
    /// Generates a pair of valid ACGT primers together with two complementarity thresholds, for
    /// HasPrimerDimer monotonicity. Primers are ≥ 8 bases so the 3' comparison window is full.
    /// </summary>
    private static Arbitrary<(string p1, string p2, int thrLo, int thrHi)> DimerThresholdArbitrary() =>
        (from p1 in ValidPrimerGen(8, 30)
         from p2 in ValidPrimerGen(8, 30)
         from a in Gen.Choose(1, 8)
         from b in Gen.Choose(1, 8)
         select (p1, p2, Math.Min(a, b), Math.Max(a, b))).ToArbitrary();

    /// <summary>
    /// Independent transcription of the production complement relation (A↔T, G↔C). Used to
    /// recompute the terminal complementary-pair count for HasPrimerDimer anchors.
    /// </summary>
    private static bool IsWatsonCrick(char a, char b) =>
        (a, b) is ('A', 'T') or ('T', 'A') or ('G', 'C') or ('C', 'G');

    /// <summary>Reverse complement built independently (no production helper) for dimer oracle work.</summary>
    private static string RevComp(string s)
    {
        var chars = s.ToUpperInvariant().ToCharArray();
        Array.Reverse(chars);
        for (int i = 0; i < chars.Length; i++)
            chars[i] = chars[i] switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', var c => c };
        return new string(chars);
    }

    #endregion

    #region PRIMER-STRUCT-001 — Calculate3PrimeStability (R: ΔG ≤ 0, NN oracle)

    /// <summary>
    /// INV-NN-oracle: for ANY ACGT sequence of length ≥ 5, production
    /// <see cref="PrimerDesigner.Calculate3PrimeStability"/> equals the independent SantaLucia
    /// nearest-neighbor oracle within 1e-9. This is the rigorous core: it pins every NN step and
    /// both initiation terms to literally transcribed published constants rather than the code's own.
    /// Source: SantaLucia (1998); Primer_Structure_Analysis.md §2.2/§4.2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_MatchesSantaLuciaOracle()
    {
        return Prop.ForAll(StabilitySequenceArbitrary(), seq =>
        {
            double actual = PrimerDesigner.Calculate3PrimeStability(seq);
            double expected = ExpectedThreePrimeStability(seq);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"ΔG mismatch for '{seq}': {actual} vs oracle {expected}");
        });
    }

    /// <summary>
    /// INV-01 (R): ΔG ≤ 0 for every sequence of length ≥ 5. Every NN step is negative and the
    /// most positive initiation pair (+1.03 twice) cannot outweigh the four steps, so the duplex
    /// stability is never positive (the anchor TATAT = −0.86 is the documented least-stable case).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_IsNonPositive_ForLengthAtLeast5()
    {
        return Prop.ForAll(StabilitySequenceArbitrary(), seq =>
        {
            double dg = PrimerDesigner.Calculate3PrimeStability(seq);
            return (dg <= 0.0 && double.IsFinite(dg))
                .Label($"ΔG {dg} not ≤ 0 (finite) for '{seq}'");
        });
    }

    /// <summary>
    /// INV-02 (Edge): ΔG = 0 for any sequence shorter than 5 nt (the method short-circuits).
    /// Source: Primer_Structure_Analysis.md §2.4 INV-02.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_IsZero_ForLengthBelow5()
    {
        var shortGen = (from len in Gen.Choose(0, 4)
                        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                        select new string(chars)).ToArbitrary();
        return Prop.ForAll(shortGen, seq =>
            (PrimerDesigner.Calculate3PrimeStability(seq) == 0.0)
                .Label($"ΔG not 0 for short '{seq}'"));
    }

    /// <summary>
    /// INV-M (terminal G·C ⇒ not less stable): replacing the terminal A/T of an A/T-flanked 5-mer
    /// with G/C makes the duplex no LESS stable (ΔG not greater). G·C contributes a stronger NN
    /// step and a lower initiation penalty (+0.98 vs +1.03), so more terminal G/C ⇒ ΔG decreases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_MoreTerminalGc_NotLessStable()
    {
        // Inner trinucleotide is fixed; only the two terminal bases vary A/T → G/C.
        var arb = (from inner in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(3)
                   select new string(inner)).ToArbitrary();
        return Prop.ForAll(arb, inner =>
        {
            double dgAt = PrimerDesigner.Calculate3PrimeStability("A" + inner + "T");
            double dgGc = PrimerDesigner.Calculate3PrimeStability("G" + inner + "C");
            return (dgGc <= dgAt + 1e-9)
                .Label($"GC-terminal not ≤ AT-terminal for inner '{inner}': {dgGc} vs {dgAt}");
        });
    }

    /// <summary>
    /// INV-D (Determinism): repeated evaluation of the same sequence is bit-identical.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_IsDeterministic()
    {
        return Prop.ForAll(StabilitySequenceArbitrary(), seq =>
        {
            double a = PrimerDesigner.Calculate3PrimeStability(seq);
            double b = PrimerDesigner.Calculate3PrimeStability(seq);
            return a.Equals(b).Label($"non-deterministic ΔG for '{seq}': {a} vs {b}");
        });
    }

    /// <summary>
    /// Evidence anchors: the two canonical 5-mers documented in §2.2 with literal published values —
    /// GCGCG = −6.86 kcal/mol (most stable) and TATAT = −0.86 (least stable). A shared-but-wrong
    /// constant in production+oracle could pass the oracle property; these absolute values cannot.
    /// </summary>
    [TestCase("GCGCG", -6.86, TestName = "SantaLucia: GCGCG → −6.86 (most stable)")]
    [TestCase("TATAT", -0.86, TestName = "SantaLucia: TATAT → −0.86 (least stable)")]
    [Category("Property")]
    public void Calculate3PrimeStability_CanonicalAnchors(string seq, double expected)
    {
        double actual = PrimerDesigner.Calculate3PrimeStability(seq);
        Assert.That(actual, Is.EqualTo(expected).Within(0.01),
            $"ΔG for '{seq}' must match the documented SantaLucia value {expected} kcal/mol");
    }

    #endregion

    #region PRIMER-STRUCT-001 — HasHairpinPotential (INV-01, monotonicity, case)

    /// <summary>
    /// INV-01 (R): HasHairpinPotential returns false whenever the sequence is shorter than the
    /// structural minimum 2·minStemLength + minLoopLength — a valid stem-loop cannot fit. Exercised
    /// over randomized sequences AND randomized stem/loop knobs. Source: §2.4 INV-01.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_BelowMinStructureSize_IsFalse()
    {
        return Prop.ForAll(HairpinKnobArbitrary(), t =>
        {
            var (seq, stem, loop) = t;
            if (seq.Length >= 2 * stem + loop)
                return true.ToProperty(); // property only constrains the sub-minimum regime.
            return (!PrimerDesigner.HasHairpinPotential(seq, stem, loop))
                .Label($"hairpin reported below min size: len {seq.Length} < 2·{stem}+{loop} for '{seq}'");
        });
    }

    /// <summary>
    /// INV-M (monotonicity in permissiveness): a shorter required stem is strictly more permissive,
    /// so if a hairpin is detected at stem = k it must still be detected at stem = k−1 (loop fixed).
    /// Detection sets are nested: lowering the stem threshold never turns true into false.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_LowerStem_NeverLosesDetection()
    {
        var arb = (from stem in Gen.Choose(3, 6)
                   from loop in Gen.Choose(1, 4)
                   from len in Gen.Choose(2 * stem + loop, 60)
                   from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                   select (new string(chars), stem, loop)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var (seq, stem, loop) = t;
            bool atK = PrimerDesigner.HasHairpinPotential(seq, stem, loop);
            bool atKMinus1 = PrimerDesigner.HasHairpinPotential(seq, stem - 1, loop);
            // detected at stricter stem ⇒ must be detected at more permissive (shorter) stem.
            return (!atK || atKMinus1)
                .Label($"monotonicity broken for '{seq}' loop={loop}: stem {stem}={atK}, stem {stem - 1}={atKMinus1}");
        });
    }

    /// <summary>
    /// INV-P (case-insensitivity): upper-, lower-, and mixed-case spellings of the same sequence
    /// give the identical hairpin verdict (the source upper-cases before comparison). Source: §3.3.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_IsCaseInsensitive()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 40), seq =>
        {
            char[] m = seq.ToCharArray();
            for (int i = 0; i < m.Length; i += 2) m[i] = char.ToLowerInvariant(m[i]);
            bool up = PrimerDesigner.HasHairpinPotential(seq.ToUpperInvariant());
            bool lo = PrimerDesigner.HasHairpinPotential(seq.ToLowerInvariant());
            bool mi = PrimerDesigner.HasHairpinPotential(new string(m));
            return (up == lo && up == mi)
                .Label($"case sensitivity for '{seq}': up={up}, lo={lo}, mixed={mi}");
        });
    }

    /// <summary>
    /// INV-D (Determinism): the hairpin verdict is identical across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_IsDeterministic()
    {
        return Prop.ForAll(HairpinKnobArbitrary(), t =>
        {
            var (seq, stem, loop) = t;
            bool a = PrimerDesigner.HasHairpinPotential(seq, stem, loop);
            bool b = PrimerDesigner.HasHairpinPotential(seq, stem, loop);
            return (a == b).Label($"non-deterministic hairpin for '{seq}'");
        });
    }

    /// <summary>
    /// Constructed positive anchor: a 5' stem, a loop, then the reverse complement of the stem as
    /// the 3' stem forms a genuine stem-loop. Stem "GGGG" (positions 0-3), loop "AAA" (3 nt ≥
    /// minLoop), 3' stem "CCCC" = revcomp("GGGG"). With minStem=4/minLoop=3 this must be detected.
    /// Source: §2.2 stem-loop model; §6.1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasHairpinPotential_EngineeredStemLoop_ReturnsTrue()
    {
        // 5' stem GGGG | loop AAA | 3' stem CCCC (= reverse complement of GGGG read 3'→5').
        const string hairpin = "GGGGAAACCCC";
        Assert.That(PrimerDesigner.HasHairpinPotential(hairpin, minStemLength: 4, minLoopLength: 3),
            Is.True, "engineered stem(GGGG)+loop(AAA)+revcomp-stem(CCCC) must form a hairpin");
    }

    /// <summary>
    /// Constructed negative anchor: a homopolymer has no self-complementary stem (A never pairs with
    /// A), so no hairpin is possible regardless of length. Source: §6.1 "No self-complementary regions".
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasHairpinPotential_Homopolymer_ReturnsFalse()
    {
        Assert.That(PrimerDesigner.HasHairpinPotential(new string('A', 40)),
            Is.False, "a poly-A run has no complementary stem and cannot form a hairpin");
    }

    #endregion

    #region PRIMER-STRUCT-001 — HasPrimerDimer (empty, monotonicity, constructed)

    /// <summary>
    /// INV-Edge: an empty or null primer on either side yields false — no terminal complementarity
    /// can be evaluated. Source: §3.3 / §6.1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasPrimerDimer_EmptyOrNull_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.HasPrimerDimer("", "ACGTACGT"), Is.False);
            Assert.That(PrimerDesigner.HasPrimerDimer("ACGTACGT", ""), Is.False);
            Assert.That(PrimerDesigner.HasPrimerDimer(null!, "ACGTACGT"), Is.False);
            Assert.That(PrimerDesigner.HasPrimerDimer("ACGTACGT", null!), Is.False);
        });
    }

    /// <summary>
    /// INV-M (monotonicity in threshold): a higher minComplementarity is stricter, so if a dimer is
    /// flagged at the high threshold it must also be flagged at the lower one (same primers). Raising
    /// the threshold can never INCREASE detections. Source: §2.2 (threshold on terminal window).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasPrimerDimer_HigherThreshold_NeverIncreasesDetection()
    {
        return Prop.ForAll(DimerThresholdArbitrary(), t =>
        {
            var (p1, p2, lo, hi) = t;
            bool atHi = PrimerDesigner.HasPrimerDimer(p1, p2, hi);
            bool atLo = PrimerDesigner.HasPrimerDimer(p1, p2, lo);
            // detected at the stricter (higher) threshold ⇒ must be detected at the looser (lower) one.
            return (!atHi || atLo)
                .Label($"threshold monotonicity broken for '{p1}'/'{p2}': thr {hi}={atHi}, thr {lo}={atLo}");
        });
    }

    /// <summary>
    /// INV-D (Determinism): the dimer verdict is identical across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasPrimerDimer_IsDeterministic()
    {
        return Prop.ForAll(DimerThresholdArbitrary(), t =>
        {
            var (p1, p2, _, thr) = t;
            bool a = PrimerDesigner.HasPrimerDimer(p1, p2, thr);
            bool b = PrimerDesigner.HasPrimerDimer(p1, p2, thr);
            return (a == b).Label($"non-deterministic dimer for '{p1}'/'{p2}'");
        });
    }

    /// <summary>
    /// Constructed positive anchor with independent recomputation: the production rule compares the
    /// last min(8,len) bases of primer1 against the FIRST min(8,len) bases of revcomp(primer2) and
    /// counts Watson-Crick pairs. We build a primer pair whose terminal windows are fully complementary
    /// and verify (a) production flags the dimer at minComplementarity=4, and (b) the independently
    /// recomputed complementary count equals 8 (the full window) — confirming the rule, not just the bit.
    /// Source: §2.2 / §5.2 (last min(8,len) bases vs start of revcomp of the other primer).
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasPrimerDimer_EngineeredComplementaryEnds_ReturnsTrue()
    {
        // primer2 fixed; choose primer1 so its 3' 8-mer == first 8 of revcomp(primer2) ⇒ all complementary.
        const string primer2 = "TTTTTTTTGGCCAATT";
        string rc2 = RevComp(primer2);                 // independent reverse complement
        string window2 = rc2.Substring(0, 8);          // first 8 bases of revcomp(primer2)
        // Make primer1 end in a stretch that is Watson-Crick complementary to window2 base-by-base.
        var tail = new char[8];
        for (int i = 0; i < 8; i++)
            tail[i] = window2[i] switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', _ => 'A' };
        string primer1 = "ACGTACGT" + new string(tail); // length 16, last 8 = engineered complement

        // Independent recomputation of the terminal complementary-pair count.
        string end1 = primer1.Substring(primer1.Length - 8);
        int comp = 0;
        for (int i = 0; i < 8; i++)
            if (IsWatsonCrick(end1[i], window2[i])) comp++;

        Assert.Multiple(() =>
        {
            Assert.That(comp, Is.EqualTo(8), "engineered terminal window must be fully complementary");
            Assert.That(PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity: 4),
                Is.True, "fully complementary 3' ends must be flagged as a primer-dimer");
        });
    }

    #endregion

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

    #region Generators &amp; Theory Oracle (PROBE-DESIGN-001)

    /// <summary>
    /// Independent GC fraction oracle: (#G + #C) / length on the uppercased sequence, computed
    /// from raw base counts — NOT routed through any production helper. Returns 0 for empty input.
    /// (The production code stores GcContent as a fraction in [0,1], cf. INV-02.)
    /// </summary>
    private static double ExpectedGcFraction(string seq)
    {
        if (seq.Length == 0) return 0;
        var (_, gc) = CountAtGc(seq);
        return (double)gc / seq.Length;
    }

    /// <summary>
    /// Generates a randomized probe-design scenario: a target long enough to host probe windows,
    /// plus a randomized parameter window (the design "config knobs" — length, GC and Tm bands).
    /// The windows are deliberately varied so the invariants are exercised across configurations
    /// rather than a single hard-coded preset (cf. ADVANCED_TESTING_CHECKLIST multi-knob probe
    /// mention: length range + GC range + Tm range). Lengths are kept short (probe windows
    /// ~12-26 bp) so hundreds of cases run quickly while still spanning the length filter.
    /// </summary>
    private static Arbitrary<(string target, ProbeDesigner.ProbeParameters param)>
        ProbeScenarioArbitrary() =>
        (from minLen in Gen.Choose(12, 18)
         from extra in Gen.Choose(0, 8)
         let maxLen = minLen + extra
         from minGcPct in Gen.Choose(25, 50)
         from gcSpanPct in Gen.Choose(10, 40)
         let maxGcPct = Math.Min(95, minGcPct + gcSpanPct)
         from minTm in Gen.Choose(30, 60)
         from tmSpan in Gen.Choose(5, 35)
         let maxTm = minTm + tmSpan
         // Target must be at least maxLen long so windows exist; add slack for multiple windows.
         from slack in Gen.Choose(0, 30)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(maxLen + slack)
         let param = new ProbeDesigner.ProbeParameters(
             MinLength: minLen, MaxLength: maxLen,
             MinTm: minTm, MaxTm: maxTm,
             MinGc: minGcPct / 100.0, MaxGc: maxGcPct / 100.0,
             MaxHomopolymer: 6,
             AvoidSecondaryStructure: false,
             MaxSelfComplementarity: 0.5)
         select (new string(chars), param)).ToArbitrary();

    /// <summary>
    /// Generates a target plus two NESTED Tm windows (B ⊆ A) over the same length/GC parameters,
    /// for the subset-monotonicity property. Both windows share identical length and GC bands so
    /// the only difference is the Tm acceptance band — window B is strictly no wider than window A.
    /// </summary>
    private static Arbitrary<(string target, ProbeDesigner.ProbeParameters wide, ProbeDesigner.ProbeParameters narrow)>
        NestedTmScenarioArbitrary() =>
        (from minLen in Gen.Choose(12, 18)
         from extra in Gen.Choose(0, 6)
         let maxLen = minLen + extra
         // Wide Tm window [aMin, aMax]; narrow window [bMin, bMax] nested inside it.
         from aMin in Gen.Choose(20, 50)
         from aSpan in Gen.Choose(30, 70)
         let aMax = aMin + aSpan
         from bLeftPad in Gen.Choose(0, aSpan / 2)
         from bRightPad in Gen.Choose(0, aSpan / 2)
         let bMin = aMin + bLeftPad
         let bMax = aMax - bRightPad
         from slack in Gen.Choose(10, 40)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(maxLen + slack)
         let wide = new ProbeDesigner.ProbeParameters(
             MinLength: minLen, MaxLength: maxLen,
             MinTm: aMin, MaxTm: aMax,
             MinGc: 0.30, MaxGc: 0.70,
             MaxHomopolymer: 6, AvoidSecondaryStructure: false, MaxSelfComplementarity: 0.5)
         let narrow = wide with { MinTm = bMin, MaxTm = bMax }
         select (new string(chars), wide, narrow)).ToArbitrary();

    #endregion

    #region PROBE-DESIGN-001 — Length Invariant (R: strict, by construction)

    /// <summary>
    /// INV-R-len (strict, by construction): EVERY returned <see cref="ProbeDesigner.Probe"/> has
    /// Sequence.Length within [MinLength, MaxLength], with coordinates consistent
    /// (End − Start + 1 == Sequence.Length, Start ≥ 0, End &lt; target.Length). This is a HARD
    /// guarantee: the candidate loop only emits windows whose length lies in the configured range
    /// (ProbeDesigner.cs §DesignProbesOptimized). Exercised across randomized targets AND randomized
    /// length windows so the bound is a true config-driven invariant. Source: §4.1 step 3; the
    /// length filter is structural, never a soft penalty.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_Length_StrictlyInConfiguredWindow()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 100000).ToList();

            foreach (var p in probes)
            {
                if (p.Sequence.Length < param.MinLength || p.Sequence.Length > param.MaxLength)
                    return false.Label(
                        $"length {p.Sequence.Length} outside [{param.MinLength},{param.MaxLength}]");
                if (p.End - p.Start + 1 != p.Sequence.Length)
                    return false.Label($"coords inconsistent: {p.Start}..{p.End} vs len {p.Sequence.Length}");
                if (p.Start < 0 || p.End >= target.Length)
                    return false.Label($"coords out of target bounds: {p.Start}..{p.End}, target len {target.Length}");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region PROBE-DESIGN-001 — GC Invariant (P: RELAXED band + warning equivalence)

    /// <summary>
    /// INV-P-GC (RELAXED band, NOT strict): GC is only a −0.3 soft penalty in
    /// <c>EvaluateProbeWithGc</c>; the sole HARD GC bound is the early-rejection filter
    /// <c>gc &lt; MinGc − 0.1 || gc &gt; MaxGc + 0.1</c>. Therefore the guaranteed invariant is
    /// GcContent ∈ [MinGc − 0.1, MaxGc + 0.1] (a fraction) — NOT the strict [MinGc, MaxGc].
    /// The production GcContent is also recomputed independently (ExpectedGcFraction) to confirm
    /// the stored fraction is itself correct (INV-02: gcCount/length). Source: ProbeDesigner.cs
    /// §DesignProbesOptimized early-rejection; §2.2 penalty table; doc §6 "outside range" is a warning.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_GcContent_WithinRelaxedBand_AndMatchesOracle()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 100000).ToList();

            foreach (var p in probes)
            {
                double oracleGc = ExpectedGcFraction(p.Sequence);
                if (Math.Abs(p.GcContent - oracleGc) > 1e-9)
                    return false.Label($"GcContent {p.GcContent} != oracle {oracleGc} for '{p.Sequence}'");

                // Relaxed hard band: [MinGc − 0.1, MaxGc + 0.1]. A small epsilon absorbs fp jitter.
                if (p.GcContent < param.MinGc - 0.1 - 1e-9 || p.GcContent > param.MaxGc + 0.1 + 1e-9)
                    return false.Label(
                        $"GcContent {p.GcContent} outside relaxed band " +
                        $"[{param.MinGc - 0.1},{param.MaxGc + 0.1}]");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-P-GC-warning (equivalence): a returned probe carries a GC warning IFF its GcContent is
    /// strictly outside the configured [MinGc, MaxGc] band. The GC penalty branch in
    /// <c>EvaluateProbeWithGc</c> adds exactly the "GC content … outside range" warning precisely
    /// when <c>gc &lt; MinGc || gc &gt; MaxGc</c>. We recompute the GC fraction independently and
    /// verify the warning⇔outside-strict-range biconditional. This is the honest replacement for a
    /// (WRONG) strict-band assertion. Source: ProbeDesigner.cs §EvaluateProbeWithGc.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_GcWarning_IffOutsideStrictBand()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 100000).ToList();

            foreach (var p in probes)
            {
                double gc = ExpectedGcFraction(p.Sequence);
                bool outsideStrict = gc < param.MinGc || gc > param.MaxGc;
                bool hasGcWarning = p.Warnings.Any(w => w.StartsWith("GC content"));

                if (outsideStrict != hasGcWarning)
                    return false.Label(
                        $"GC warning {hasGcWarning} != outsideStrict {outsideStrict} " +
                        $"(gc {gc}, band [{param.MinGc},{param.MaxGc}]) for '{p.Sequence}'");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region PROBE-DESIGN-001 — Tm Monotonicity (M: SUBSET, not naive count)

    /// <summary>
    /// INV-M-Tm (subset monotonicity): Tm is only a −0.3 soft penalty (never a hard filter);
    /// probes are dropped solely when total Score ≤ 0, then the top <c>maxProbes</c> by score are
    /// returned. So for the SAME target and parameters identical except the Tm band, if window
    /// B ⊆ window A (B narrower) then count(B) ≤ count(A): each probe's score under B is ≤ its
    /// score under A (the only difference is possibly an extra −0.3 Tm penalty), so the
    /// positive-score set under B is a subset of that under A. A LARGE maxProbes (100000) is used
    /// so neither result saturates at the cap — otherwise both saturate and the relation is vacuous.
    /// This is the HONEST monotonic claim, NOT the naive "narrower window ⇒ fewer windows scanned".
    /// Source: ProbeDesigner.cs §EvaluateProbeWithGc Tm penalty + §DesignProbesOptimized score gate.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_NarrowerTmWindow_NeverMoreProbes()
    {
        return Prop.ForAll(NestedTmScenarioArbitrary(), s =>
        {
            var (target, wide, narrow) = s;
            const int cap = 100000;
            int countWide = ProbeDesigner.DesignProbes(target, wide, maxProbes: cap).Count();
            int countNarrow = ProbeDesigner.DesignProbes(target, narrow, maxProbes: cap).Count();

            return (countNarrow <= countWide)
                .Label($"narrower Tm window produced MORE probes: " +
                       $"narrow [{narrow.MinTm},{narrow.MaxTm}]={countNarrow} > " +
                       $"wide [{wide.MinTm},{wide.MaxTm}]={countWide}");
        });
    }

    #endregion

    #region PROBE-DESIGN-001 — Determinism &amp; Edge Cases

    /// <summary>
    /// INV-D (Determinism): identical (target, params, maxProbes) ⇒ identical probe list —
    /// same sequences, starts, Tm, GcContent, Score, AND order. Two runs are compared element-wise.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_IsDeterministic()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var a = ProbeDesigner.DesignProbes(target, param, maxProbes: 50).ToList();
            var b = ProbeDesigner.DesignProbes(target, param, maxProbes: 50).ToList();

            if (a.Count != b.Count)
                return false.Label($"count differs: {a.Count} vs {b.Count}");

            for (int i = 0; i < a.Count; i++)
            {
                if (a[i].Sequence != b[i].Sequence || a[i].Start != b[i].Start ||
                    !a[i].Tm.Equals(b[i].Tm) || !a[i].GcContent.Equals(b[i].GcContent) ||
                    !a[i].Score.Equals(b[i].Score))
                    return false.Label($"probe #{i} differs between runs");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-Edge: null and empty targets yield no probes; a target shorter than MinLength yields
    /// no probes (the explicit early return in <see cref="ProbeDesigner.DesignProbes"/>).
    /// Source: ProbeDesigner.cs §DesignProbes guard; doc §6.1 edge cases.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_NullEmptyOrTooShort_ReturnsNoProbes()
    {
        var param = ProbeDesigner.Defaults.qPCR; // MinLength 20
        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.DesignProbes(null!, param).ToList(), Is.Empty, "null target");
            Assert.That(ProbeDesigner.DesignProbes("", param).ToList(), Is.Empty, "empty target");
            Assert.That(ProbeDesigner.DesignProbes(new string('A', 19), param).ToList(), Is.Empty,
                "target shorter than MinLength (20)");
        });
    }

    #endregion

    #region PROBE-DESIGN-001 — Evidence-Based Engineered Anchors

    /// <summary>
    /// Engineered length + GC anchor with hand-computed expectations. A 24-mer with exactly 50%
    /// GC, in a window forcing MinLength == MaxLength == 24 and a GC band that contains 0.50, must
    /// return a single probe whose Sequence is the whole target, whose length is 24, whose
    /// coordinates are 0..23, and whose GcContent equals the independently computed 12/24 = 0.5.
    /// Because GC is inside the strict band, NO GC warning is present.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_EngineeredFiftyPercentGc_FixedLengthWindow()
    {
        // 24 bases: 12 G/C (ACGT repeated 6× → 6 G + 6 C = 12) ⇒ GC fraction = 12/24 = 0.5.
        const string target = "ACGTACGTACGTACGTACGTACGT";
        var param = new ProbeDesigner.ProbeParameters(
            MinLength: 24, MaxLength: 24,
            MinTm: 0, MaxTm: 200,        // wide Tm so no Tm penalty drops the probe
            MinGc: 0.40, MaxGc: 0.60,    // 0.50 is strictly inside
            MaxHomopolymer: 6, AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 10).ToList();

        Assert.That(probes, Has.Count.EqualTo(1), "exactly one length-24 window exists");
        var p = probes[0];
        Assert.Multiple(() =>
        {
            Assert.That(p.Sequence, Is.EqualTo(target));
            Assert.That(p.Sequence.Length, Is.EqualTo(24));
            Assert.That(p.Start, Is.EqualTo(0));
            Assert.That(p.End, Is.EqualTo(23));
            Assert.That(p.GcContent, Is.EqualTo(0.5).Within(1e-9), "12 G/C of 24 bases = 0.5");
            Assert.That(p.Warnings.Any(w => w.StartsWith("GC content")), Is.False,
                "GC 0.5 is inside [0.40,0.60] ⇒ no GC warning");
        });
    }

    /// <summary>
    /// Engineered warning anchor: a low-GC target (20% GC) whose fraction (0.20) is still inside
    /// the RELAXED band [MinGc − 0.1, MaxGc + 0.1] = [0.20, 0.70] but strictly BELOW the configured
    /// MinGc = 0.30 ⇒ the returned probe survives the early-rejection filter yet MUST carry a GC
    /// warning. Demonstrates the relaxed-band + warning-equivalence contract on a concrete case.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_GcInsideRelaxedBandButBelowStrict_CarriesWarning()
    {
        // 20 bases, 4 of them G/C (G,C,G,C) ⇒ GC fraction = 4/20 = 0.20.
        const string target20 = "AAGCAAGCAAAAAAAAAAAA";
        var param = new ProbeDesigner.ProbeParameters(
            MinLength: 20, MaxLength: 20,
            MinTm: 0, MaxTm: 200,
            MinGc: 0.30, MaxGc: 0.60,   // strict band; relaxed band = [0.20, 0.70]
            MaxHomopolymer: 20, AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

        var probes = ProbeDesigner.DesignProbes(target20, param, maxProbes: 10).ToList();

        Assert.That(probes, Has.Count.EqualTo(1), "0.20 is the lower edge of the relaxed band ⇒ retained");
        var p = probes[0];
        Assert.Multiple(() =>
        {
            Assert.That(p.GcContent, Is.EqualTo(0.20).Within(1e-9), "4 G/C of 20 bases = 0.20");
            Assert.That(p.GcContent, Is.LessThan(param.MinGc), "0.20 < strict MinGc 0.30");
            Assert.That(p.Warnings.Any(w => w.StartsWith("GC content")), Is.True,
                "GC below strict MinGc ⇒ GC warning present");
        });
    }

    #endregion

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
