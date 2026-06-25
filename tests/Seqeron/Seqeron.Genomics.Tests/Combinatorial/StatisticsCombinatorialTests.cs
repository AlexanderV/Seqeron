namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Statistics area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Statistics")]
public class StatisticsCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-HYDRO-001 — Hydrophobicity (Kyte-Doolittle GRAVY + profile) (Statistics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 123.
    // Spec: tests/TestSpecs/SEQ-HYDRO-001.md (SequenceStatistics.CalculateHydrophobicity / CalculateHydrophobicityProfile).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Kyte & Doolittle (1982); Biopython gravy()/protein_scale (GRAVY = Σkd/N; profile = N−W+1
    // unweighted window means; window 9 surface / 19 transmembrane).
    //
    // Checklist axes scale(2) × windowSize(3) × seqLen(3). Only the Kyte-Doolittle scale is implemented;
    // per the campaign convention we map the axes onto the real knobs and document it:
    //   • scale     → the input letter case {Upper, Lower}: the scale lookup is case-insensitive, so the
    //     profile is INVARIANT to case (INV-4).
    //   • windowSize→ the sliding-window size ∈ {3, 9, 19} (surface vs transmembrane windows).
    //   • seqLen    → protein length ∈ {10, 20, 40} (W > N yields an empty profile).
    // Grid = 2 × 3 × 3 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the profile length is a JOINT function of window and sequence length
    // (N − W + 1, or empty when W > N), every window value is the unweighted Kyte-Doolittle mean, and the
    // whole profile is invariant to input case. Each cell is checked against the KD sliding-mean re-derived
    // from the sequence.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>The published Kyte-Doolittle hydropathy scale (the ground-truth theory).</summary>
    private static readonly IReadOnlyDictionary<char, double> KyteDoolittle = new Dictionary<char, double>
    {
        ['A'] = 1.8, ['R'] = -4.5, ['N'] = -3.5, ['D'] = -3.5, ['C'] = 2.5,
        ['E'] = -3.5, ['Q'] = -3.5, ['G'] = -0.4, ['H'] = -3.2, ['I'] = 4.5,
        ['L'] = 3.8, ['K'] = -3.9, ['M'] = 1.9, ['F'] = 2.8, ['P'] = -1.6,
        ['S'] = -0.8, ['T'] = -0.7, ['W'] = -0.9, ['Y'] = -1.3, ['V'] = 4.2,
    };

    private const string CanonicalResidues = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>
    /// For every (case, window size, sequence length) the hydropathy profile has length N − W + 1 (or is
    /// empty when W &gt; N), every window value equals the unweighted Kyte-Doolittle mean, and the profile
    /// is invariant to the input letter case.
    /// </summary>
    [Test, Combinatorial]
    public void HydrophobicityProfile_CaseWindowLengthGrid_MatchesKyteDoolittleMeans(
        [Values(true, false)] bool upperCase,
        [Values(3, 9, 19)] int windowSize,
        [Values(10, 20, 40)] int seqLen)
    {
        string sequence = BuildProtein(seqLen);
        string cased = upperCase ? sequence.ToUpperInvariant() : sequence.ToLowerInvariant();

        // Independent ground truth (Kyte-Doolittle sliding mean).
        var expected = new List<double>();
        for (int i = 0; i + windowSize <= seqLen; i++)
        {
            double sum = 0;
            for (int j = 0; j < windowSize; j++)
                sum += KyteDoolittle[sequence[i + j]];
            expected.Add(sum / windowSize);
        }

        var profile = SequenceStatistics.CalculateHydrophobicityProfile(cased, windowSize).ToList();

        profile.Should().HaveCount(Math.Max(0, seqLen - windowSize + 1), "[INV-3] profile length = N − W + 1 (0 when W > N)");
        for (int i = 0; i < expected.Count; i++)
            profile[i].Should().BeApproximately(expected[i], 1e-9, "each window is the unweighted KD mean");
    }

    /// <summary>
    /// Interaction witness (window × length, empty profile): a window wider than the sequence yields no
    /// profile values, while a window that fits yields N − W + 1. Source: Biopython range(N−W+1).
    /// </summary>
    [Test]
    public void HydrophobicityProfile_WindowExceedsLength_IsEmpty()
    {
        string sequence = BuildProtein(10);

        SequenceStatistics.CalculateHydrophobicityProfile(sequence, 19).Should().BeEmpty("W = 19 > N = 10");
        SequenceStatistics.CalculateHydrophobicityProfile(sequence, 9).Should().HaveCount(2, "N − W + 1 = 10 − 9 + 1");
    }

    /// <summary>
    /// Interaction witness (case invariance + GRAVY, INV-2/INV-4): GRAVY = Σkd/N is case-insensitive, and a
    /// single residue's GRAVY equals its Kyte-Doolittle value. Source: Expasy / Biopython gravy().
    /// </summary>
    [Test]
    public void Hydrophobicity_GravyCaseInsensitiveAndSingleResidue()
    {
        string sequence = BuildProtein(20);

        SequenceStatistics.CalculateHydrophobicity(sequence.ToLowerInvariant())
            .Should().BeApproximately(SequenceStatistics.CalculateHydrophobicity(sequence.ToUpperInvariant()), 1e-12,
                "[INV-4] GRAVY is case-insensitive");
        SequenceStatistics.CalculateHydrophobicity("I").Should().BeApproximately(4.5, 1e-12, "[INV-2] single residue → its kd value");
    }

    /// <summary>
    /// Witness (worked GRAVY): the GRAVY of the all-20-residue string is the mean of the Kyte-Doolittle
    /// values, re-derived independently. Source: Kyte & Doolittle (1982).
    /// </summary>
    [Test]
    public void Hydrophobicity_AllResidues_MatchesMeanOfScale()
    {
        double expected = CanonicalResidues.Sum(c => KyteDoolittle[c]) / CanonicalResidues.Length;

        SequenceStatistics.CalculateHydrophobicity(CanonicalResidues).Should().BeApproximately(expected, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-SECSTRUCT-001 — Chou-Fasman secondary-structure propensity profile (Statistics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 126.
    // Spec: tests/TestSpecs/SEQ-SECSTRUCT-001.md (SequenceStatistics.PredictSecondaryStructure).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Chou & Fasman (1978) conformational propensities (Pa helix, Pb sheet, Pt turn); per-window
    // mean over a sliding window, stepping one residue, N→C order; profile length N − W + 1.
    //
    // Checklist axes method(2) × windowSize(3) × seqLen(3). A single Chou-Fasman method is implemented
    // (returning all three propensities); per the campaign convention we map the axes onto the real knobs:
    //   • method     → input letter case {Upper, Lower}: the propensity lookup is case-insensitive, so the
    //     profile is INVARIANT to case.
    //   • windowSize → the sliding-window size ∈ {3, 7, 11}.
    //   • seqLen     → protein length ∈ {10, 20, 40} (W > N yields an empty profile).
    // Grid = 2 × 3 × 3 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the profile length is a JOINT function of window and sequence length, every
    // window emits the unweighted mean of all three Chou-Fasman propensities, and the profile is invariant
    // to case. Each cell is checked against the propensity means re-derived from the sequence.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>The published Chou-Fasman (Pa helix, Pb sheet, Pt turn) propensities (the ground-truth theory).</summary>
    private static readonly IReadOnlyDictionary<char, (double Helix, double Sheet, double Turn)> ChouFasman =
        new Dictionary<char, (double, double, double)>
        {
            ['A'] = (1.42, 0.83, 0.66), ['R'] = (0.98, 0.93, 0.95), ['N'] = (0.67, 0.89, 1.56),
            ['D'] = (1.01, 0.54, 1.46), ['C'] = (0.70, 1.19, 1.19), ['E'] = (1.51, 0.37, 0.74),
            ['Q'] = (1.11, 1.10, 0.98), ['G'] = (0.57, 0.75, 1.56), ['H'] = (1.00, 0.87, 0.95),
            ['I'] = (1.08, 1.60, 0.47), ['L'] = (1.21, 1.30, 0.59), ['K'] = (1.14, 0.74, 1.01),
            ['M'] = (1.45, 1.05, 0.60), ['F'] = (1.13, 1.38, 0.60), ['P'] = (0.57, 0.55, 1.52),
            ['S'] = (0.77, 0.75, 1.43), ['T'] = (0.83, 1.19, 0.96), ['W'] = (1.08, 1.37, 0.96),
            ['Y'] = (0.69, 1.47, 1.14), ['V'] = (1.06, 1.70, 0.50),
        };

    /// <summary>
    /// For every (case, window size, sequence length) the secondary-structure profile has length N − W + 1
    /// (empty when W &gt; N), every window emits the unweighted mean of the three Chou-Fasman propensities,
    /// and the profile is invariant to the input case.
    /// </summary>
    [Test, Combinatorial]
    public void PredictSecondaryStructure_CaseWindowLengthGrid_MatchesChouFasmanMeans(
        [Values(true, false)] bool upperCase,
        [Values(3, 7, 11)] int windowSize,
        [Values(10, 20, 40)] int seqLen)
    {
        string sequence = BuildProtein(seqLen);
        string cased = upperCase ? sequence.ToUpperInvariant() : sequence.ToLowerInvariant();

        var profile = SequenceStatistics.PredictSecondaryStructure(cased, windowSize).ToList();

        profile.Should().HaveCount(Math.Max(0, seqLen - windowSize + 1), "profile length = N − W + 1 (0 when W > N)");
        for (int i = 0; i < profile.Count; i++)
        {
            double h = 0, s = 0, t = 0;
            for (int j = 0; j < windowSize; j++)
            {
                var p = ChouFasman[sequence[i + j]];
                h += p.Helix; s += p.Sheet; t += p.Turn;
            }
            profile[i].Helix.Should().BeApproximately(h / windowSize, 1e-9, "window helix propensity mean");
            profile[i].Sheet.Should().BeApproximately(s / windowSize, 1e-9, "window sheet propensity mean");
            profile[i].Turn.Should().BeApproximately(t / windowSize, 1e-9, "window turn propensity mean");
        }
    }

    /// <summary>
    /// Interaction witness (window × length, empty profile): a window wider than the sequence yields no
    /// profile, while a fitting window yields N − W + 1 entries. Source: Chou & Fasman windowing.
    /// </summary>
    [Test]
    public void PredictSecondaryStructure_WindowExceedsLength_IsEmpty()
    {
        var sequence = BuildProtein(8);

        SequenceStatistics.PredictSecondaryStructure(sequence, 11).Should().BeEmpty("W = 11 > N = 8");
        SequenceStatistics.PredictSecondaryStructure(sequence, 7).Should().HaveCount(2, "N − W + 1 = 8 − 7 + 1");
    }

    /// <summary>
    /// Interaction witness (conformation discrimination): a helix-favouring stretch (poly-E, Pa = 1.51)
    /// gives a window helix propensity above 1 and well above its sheet propensity, while a sheet-favouring
    /// stretch (poly-V, Pb = 1.70) is the reverse. Source: Chou & Fasman (1978).
    /// </summary>
    [Test]
    public void PredictSecondaryStructure_FavouringStretches_DiscriminateConformation()
    {
        var helix = SequenceStatistics.PredictSecondaryStructure(new string('E', 7), 7).Single();
        helix.Helix.Should().BeApproximately(1.51, 1e-9).And.BeGreaterThan(helix.Sheet);

        var sheet = SequenceStatistics.PredictSecondaryStructure(new string('V', 7), 7).Single();
        sheet.Sheet.Should().BeApproximately(1.70, 1e-9).And.BeGreaterThan(sheet.Helix);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-THERMO-001 — DNA nearest-neighbor thermodynamics (Statistics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 129.
    // Spec: tests/TestSpecs/SEQ-THERMO-001.md (SequenceStatistics.CalculateThermodynamics).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Allawi & SantaLucia (1997) / SantaLucia (1998) unified NN model — ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000;
    // Tm = (1000·ΔH°)/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.987; method-5 Na⁺ salt correction.
    //
    // Checklist axes saltConc(3) × seqLen(3) × gcContent(3) map onto the real knobs: naConcentration ∈
    // {0.01, 0.05, 1.0} M, length ∈ {6, 12, 20}, GC fraction ∈ {0.0, 0.5, 1.0}. Grid = 3³ = 27 = the
    // checklist's "Full Combos" for this row.
    //
    // The combinatorial point: ΔH/ΔS/ΔG/Tm are a JOINT function of salt, length and GC, but the two
    // defining thermodynamic identities (the Gibbs relation, INV-02, and the Tm equation, INV-03) must hold
    // in EVERY cell. Each cell re-derives both relations from the reported ΔH/ΔS and checks them.
    // ═══════════════════════════════════════════════════════════════════════

    private const double GasConstant = 1.987;       // cal/(mol·K)
    private const double StrandConcentration = 2.5e-7; // C_T (default), F = 4

    /// <summary>
    /// For every (salt, length, GC) the reported thermodynamics satisfy the Gibbs relation
    /// ΔG = ΔH − 310.15·ΔS/1000 (INV-02) and the Tm equation Tm = 1000·ΔH/(ΔS + R·ln(C_T/4)) − 273.15
    /// (INV-03), within rounding tolerance.
    /// </summary>
    [Test, Combinatorial]
    public void CalculateThermodynamics_SaltLengthGcGrid_SatisfiesGibbsAndTmRelations(
        [Values(0.01, 0.05, 1.0)] double naConcentration,
        [Values(6, 12, 20)] int seqLen,
        [Values(0.0, 0.5, 1.0)] double gcFraction)
    {
        string dna = BuildDna(seqLen, gcFraction);

        var t = SequenceStatistics.CalculateThermodynamics(dna, naConcentration);

        double expectedG = t.DeltaH - 310.15 * t.DeltaS / 1000.0;
        t.DeltaG.Should().BeApproximately(expectedG, 0.05, "[INV-02] ΔG = ΔH − 310.15·ΔS/1000");

        double expectedTm = t.DeltaH * 1000.0 / (t.DeltaS + GasConstant * Math.Log(StrandConcentration / 4.0)) - 273.15;
        t.MeltingTemperature.Should().BeApproximately(expectedTm, 0.2, "[INV-03] Tm = 1000·ΔH/(ΔS + R·ln(C_T/4)) − 273.15");
    }

    /// <summary>
    /// Interaction witness (salt axis monotonicity): raising the Na⁺ concentration raises Tm (the method-5
    /// salt correction makes ΔS less negative). Source: SantaLucia (1998) salt correction.
    /// </summary>
    [Test]
    public void CalculateThermodynamics_SaltAxis_RaisesTm()
    {
        const string dna = "GATCGATCGATC";

        double low = SequenceStatistics.CalculateThermodynamics(dna, 0.01).MeltingTemperature;
        double mid = SequenceStatistics.CalculateThermodynamics(dna, 0.05).MeltingTemperature;
        double high = SequenceStatistics.CalculateThermodynamics(dna, 1.0).MeltingTemperature;

        mid.Should().BeGreaterThan(low);
        high.Should().BeGreaterThan(mid);
    }

    /// <summary>
    /// Interaction witness (GC axis monotonicity): a GC-rich duplex melts higher than an AT-rich one of the
    /// same length and salt (three hydrogen bonds vs two; stronger NN stacking). Source: Allawi & SantaLucia (1997).
    /// </summary>
    [Test]
    public void CalculateThermodynamics_GcAxis_RaisesTm()
    {
        double atTm = SequenceStatistics.CalculateThermodynamics(BuildDna(12, 0.0), 0.05).MeltingTemperature;
        double gcTm = SequenceStatistics.CalculateThermodynamics(BuildDna(12, 1.0), 0.05).MeltingTemperature;

        gcTm.Should().BeGreaterThan(atTm, "GC-rich duplexes are more stable");
    }

    /// <summary>
    /// Witness (INV-05/06): the NN model is case-insensitive and returns all-zero for inputs shorter than a
    /// dinucleotide. Source: NN model undefined for length &lt; 2.
    /// </summary>
    [Test]
    public void CalculateThermodynamics_CaseInsensitiveAndShortInput()
    {
        var upper = SequenceStatistics.CalculateThermodynamics("GATCGATC", 0.05);
        var lower = SequenceStatistics.CalculateThermodynamics("gatcgatc", 0.05);
        lower.Should().Be(upper, "[INV-05] case-insensitive");

        SequenceStatistics.CalculateThermodynamics("A", 0.05)
            .Should().Be(new SequenceStatistics.ThermodynamicProperties(0, 0, 0, 0), "[INV-06] length < 2 → zeros");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-TM-001 — Simple melting temperature (Wallace / Marmur-Doty) (Statistics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 130.
    // Spec: tests/TestSpecs/SEQ-TM-001.md (SequenceStatistics.CalculateMeltingTemperature).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Wallace rule Tm = 2·(A+T) + 4·(G+C) for short oligos (< 14 bp); Marmur-Doty
    // Tm = 64.9 + 41·(GC − 16.4)/N otherwise.
    //
    // Checklist axes method(2) × seqLen(3) × gcContent(3) map onto the real knobs: useWallaceRule ∈
    // {true, false}, length ∈ {8, 12, 20}, GC fraction ∈ {0.0, 0.5, 1.0}. Grid = 2 × 3 × 3 = 18.
    //
    // The combinatorial point: which formula is applied is a JOINT function of the method flag AND the
    // length — Wallace applies only when useWallaceRule is set AND length < 14, otherwise Marmur-Doty — and
    // the Tm value then depends on the GC content. Each cell re-derives the selected formula from the inputs.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// For every (method flag, length, GC) the melting temperature matches the formula the selection logic
    /// picks: Wallace (2·AT + 4·GC) when useWallaceRule AND length &lt; 14, else Marmur-Doty.
    /// </summary>
    [Test, Combinatorial]
    public void CalculateMeltingTemperature_MethodLengthGcGrid_MatchesSelectedFormula(
        [Values(true, false)] bool useWallaceRule,
        [Values(8, 12, 20)] int seqLen,
        [Values(0.0, 0.5, 1.0)] double gcFraction)
    {
        string dna = BuildDna(seqLen, gcFraction);
        int gc = dna.Count(c => c is 'G' or 'C');
        int at = dna.Length - gc;

        // Independent ground truth (selection + formula).
        bool wallace = useWallaceRule && seqLen < 14;
        double expected = wallace ? 2.0 * at + 4.0 * gc : 64.9 + 41.0 * (gc - 16.4) / seqLen;

        SequenceStatistics.CalculateMeltingTemperature(dna, useWallaceRule)
            .Should().BeApproximately(expected, 1e-9, "Wallace if requested AND short, else Marmur-Doty");
    }

    /// <summary>
    /// Interaction witness (method × length): for a short oligo the method flag selects Wallace vs
    /// Marmur-Doty (different Tm), but at ≥ 14 bp the flag is ignored and Marmur-Doty always applies — the
    /// formula choice flips on the method axis only for short oligos. Source: Wallace vs Marmur-Doty regimes.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_MethodAxis_OnlyMattersForShortOligos()
    {
        string shortOligo = BuildDna(8, 0.5);
        SequenceStatistics.CalculateMeltingTemperature(shortOligo, useWallaceRule: true)
            .Should().NotBe(SequenceStatistics.CalculateMeltingTemperature(shortOligo, useWallaceRule: false),
                "Wallace and Marmur-Doty differ for a short oligo");

        string longOligo = BuildDna(20, 0.5);
        SequenceStatistics.CalculateMeltingTemperature(longOligo, useWallaceRule: true)
            .Should().Be(SequenceStatistics.CalculateMeltingTemperature(longOligo, useWallaceRule: false),
                "at ≥ 14 bp the Wallace flag is ignored");
    }

    /// <summary>
    /// Interaction witness (GC axis monotonicity): a GC-rich oligo melts higher than an AT-rich one of the
    /// same length under both formulas. Source: Wallace / Marmur-Doty GC weighting.
    /// </summary>
    [Test]
    public void CalculateMeltingTemperature_GcAxis_RaisesTm()
    {
        SequenceStatistics.CalculateMeltingTemperature(BuildDna(12, 1.0))
            .Should().BeGreaterThan(SequenceStatistics.CalculateMeltingTemperature(BuildDna(12, 0.0)), "Wallace: 4·GC > 2·AT");
        SequenceStatistics.CalculateMeltingTemperature(BuildDna(20, 1.0))
            .Should().BeGreaterThan(SequenceStatistics.CalculateMeltingTemperature(BuildDna(20, 0.0)), "Marmur-Doty: +GC term");
    }

    /// <summary>Witness (INV-6): an empty sequence returns Tm 0.</summary>
    [Test]
    public void CalculateMeltingTemperature_EmptySequence_IsZero()
    {
        SequenceStatistics.CalculateMeltingTemperature("").Should().Be(0.0);
    }

    /// <summary>
    /// Builds a DNA sequence of <paramref name="length"/> bp at approximately <paramref name="gcFraction"/>
    /// GC content: AT-only (0.0, A/T alternating), balanced (0.5, GATC repeat), or GC-only (1.0, G/C alternating).
    /// </summary>
    private static string BuildDna(int length, double gcFraction)
    {
        string pattern = gcFraction switch
        {
            <= 0.0 => "AT",
            >= 1.0 => "GC",
            _ => "GATC",
        };
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = pattern[i % pattern.Length];
        return new string(chars);
    }

    /// <summary>Builds a protein of <paramref name="length"/> recognized residues by cycling the 20 canonical amino acids.</summary>
    private static string BuildProtein(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = CanonicalResidues[i % CanonicalResidues.Length];
        return new string(chars);
    }
}
