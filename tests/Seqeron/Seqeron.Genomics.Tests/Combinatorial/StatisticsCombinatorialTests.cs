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

    /// <summary>Builds a protein of <paramref name="length"/> recognized residues by cycling the 20 canonical amino acids.</summary>
    private static string BuildProtein(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = CanonicalResidues[i % CanonicalResidues.Length];
        return new string(chars);
    }
}
