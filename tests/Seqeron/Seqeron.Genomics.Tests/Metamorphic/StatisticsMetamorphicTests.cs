using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Statistics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMPOSITION-001 — nucleotide composition (Statistics).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 121.
///
/// API under test (SequenceStatistics.CalculateNucleotideComposition):
///   Counts A/T/G/C/U/N and derives GC%, AT%, GC-skew and AT-skew — all functions of the
///   per-base counts only.
///
/// Relations (derived from count-based composition, NOT from output):
///   • INV  (permutation invariant): composition depends only on the multiset of bases, so any
///          reordering of the sequence yields the same composition.
///   • P    (complement swaps counts): the DNA complement maps A↔T and C↔G, so the complemented
///          sequence's A/T and C/G counts are the originals swapped.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class StatisticsMetamorphicTests
{
    private static string Complement(string dna) =>
        new(dna.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    #region SEQ-COMPOSITION-001 INV — composition is permutation invariant

    [Test]
    [Description("INV: nucleotide composition depends only on the multiset of bases, so reordering the sequence (here reversing it) yields the same counts and derived statistics.")]
    public void NucleotideComposition_Permutation_Invariant()
    {
        const string seq = "AACGTACGTTGGCCAATAC";
        var original = SequenceStatistics.CalculateNucleotideComposition(seq);
        var reversed = SequenceStatistics.CalculateNucleotideComposition(new string(seq.Reverse().ToArray()));

        reversed.Should().Be(original,
            because: "counts, GC%, AT% and the skews are all functions of the per-base counts, which a permutation preserves");
    }

    #endregion

    #region SEQ-COMPOSITION-001 P — the complement swaps A↔T and C↔G counts

    [Test]
    [Description("P: the DNA complement maps A↔T and C↔G, so the complemented sequence's A/T and C/G counts are the originals swapped (and the length is unchanged).")]
    public void NucleotideComposition_Complement_SwapsCounts()
    {
        const string seq = "AACGTACGTTGGCCAATAC";
        var original = SequenceStatistics.CalculateNucleotideComposition(seq);
        var complemented = SequenceStatistics.CalculateNucleotideComposition(Complement(seq));

        complemented.CountA.Should().Be(original.CountT, because: "complement turns every T into A");
        complemented.CountT.Should().Be(original.CountA, because: "complement turns every A into T");
        complemented.CountC.Should().Be(original.CountG, because: "complement turns every G into C");
        complemented.CountG.Should().Be(original.CountC, because: "complement turns every C into G");
        complemented.Length.Should().Be(original.Length, because: "complement is length-preserving");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-DINUC-001 — dinucleotide frequencies (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 122.
    //
    // API under test (SequenceStatistics.CalculateDinucleotideFrequencies):
    //   f_XY = count(XY) / (N−1) over the sliding dinucleotide windows.
    //
    // Relations (derived from dinucleotide counting, NOT from output):
    //   • INV  (reverse-complement maps each dinucleotide to its revcomp): the count of XY in a
    //          sequence equals the count of revcomp(XY) in its reverse complement, so the
    //          frequency of XY equals the frequency of revcomp(XY) in the reverse complement.
    //   • SHIFT (prepend flank adds only boundary dinucleotides): the dinucleotide multiset of
    //          flank+seq is the multisets of flank and seq plus the single boundary dinucleotide.
    // ───────────────────────────────────────────────────────────────────────────

    private static string RevComp(string dna) =>
        new(dna.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static string RevCompDinuc(string xy) => RevComp(xy);

    // Integer dinucleotide counts reconstructed from frequencies (exact for clean A/C/G/T sequences).
    private static System.Collections.Generic.Dictionary<string, int> DinucCounts(string seq)
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        int positions = seq.Length - 1;
        return freq.ToDictionary(kv => kv.Key, kv => (int)System.Math.Round(kv.Value * positions));
    }

    #region SEQ-DINUC-001 INV — reverse complement maps each dinucleotide to its revcomp

    [Test]
    [Description("INV: the count of XY in a sequence equals the count of revcomp(XY) in its reverse complement, so f_seq[XY] = f_revcomp[revcomp(XY)].")]
    public void DinucleotideFrequencies_ReverseComplement_MapsToRevcomp()
    {
        const string seq = "ACGTACGTTGGCCAATAC";
        var freqSeq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        var freqRc = SequenceStatistics.CalculateDinucleotideFrequencies(RevComp(seq));

        foreach (var (dinuc, f) in freqSeq)
            freqRc.GetValueOrDefault(RevCompDinuc(dinuc)).Should().BeApproximately(f, 1e-12,
                because: $"dinucleotide {dinuc} maps to {RevCompDinuc(dinuc)} on the reverse-complement strand");
    }

    #endregion

    #region SEQ-DINUC-001 SHIFT — prepending a flank adds only boundary dinucleotides

    [Test]
    [Description("SHIFT: the dinucleotide multiset of flank+seq equals the multisets of the flank and the sequence plus the single boundary dinucleotide (flank-last, seq-first).")]
    public void DinucleotideFrequencies_PrependFlank_AddsOnlyBoundary()
    {
        const string seq = "ACGTACGT";

        foreach (var flank in new[] { "TT", "GCGC", "AATT" })
        {
            var combined = DinucCounts(flank + seq);

            // Expected multiset = flank dinucs ⊎ seq dinucs ⊎ {boundary}.
            var expected = new System.Collections.Generic.Dictionary<string, int>();
            void Add(string d, int n) => expected[d] = expected.GetValueOrDefault(d) + n;
            foreach (var (d, n) in DinucCounts(flank)) Add(d, n);
            foreach (var (d, n) in DinucCounts(seq)) Add(d, n);
            Add($"{flank[^1]}{seq[0]}", 1); // boundary dinucleotide

            combined.Should().BeEquivalentTo(expected,
                because: $"prepending '{flank}' adds the flank's dinucleotides and one boundary dinucleotide, preserving the sequence's own");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-HYDRO-001 — hydropathy (GRAVY) (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 123.
    //
    // API under test (SequenceStatistics.CalculateHydrophobicity / CalculateHydrophobicityProfile):
    //   GRAVY = mean Kyte-Doolittle hydropathy over residues; the profile is the sliding-window
    //   mean.
    //
    // Relations (derived from the mean over residues, NOT from output):
    //   • INV  (permutation changes the profile but not the mean): GRAVY is the residue-mean,
    //          invariant to order, while the windowed profile is order-dependent.
    //   • MON  (adding a hydrophobic residue ⇒ ≥ mean): appending the most-hydrophobic residue
    //          (Ile, 4.5) — above any sub-4.5 mean — raises the GRAVY.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-HYDRO-001 INV — permutation preserves the mean but changes the profile

    [Test]
    [Description("INV: GRAVY is the residue-mean (permutation invariant), while the sliding-window profile depends on residue order.")]
    public void Hydrophobicity_Permutation_PreservesMeanChangesProfile()
    {
        const string protein = "MKLVAGWTYSDE";
        string reversed = new string(protein.Reverse().ToArray());

        SequenceStatistics.CalculateHydrophobicity(reversed)
            .Should().BeApproximately(SequenceStatistics.CalculateHydrophobicity(protein), 1e-12,
                because: "GRAVY is the mean over residues, independent of their order");

        SequenceStatistics.CalculateHydrophobicityProfile(reversed, 9).ToList()
            .Should().NotEqual(SequenceStatistics.CalculateHydrophobicityProfile(protein, 9).ToList(),
                because: "the windowed profile depends on residue order, so reversing changes it");
    }

    #endregion

    #region SEQ-HYDRO-001 MON — adding a hydrophobic residue raises the mean

    [Test]
    [Description("MON: appending isoleucine (the most hydrophobic residue, 4.5) — above any sub-4.5 mean — monotonically raises the GRAVY.")]
    public void Hydrophobicity_AddHydrophobicResidue_HigherMean()
    {
        const string hydrophilicCore = "DEKR"; // charged residues → strongly negative GRAVY
        double previous = double.MinValue;
        foreach (int isoleucines in new[] { 0, 1, 3, 6, 12 })
        {
            double gravy = SequenceStatistics.CalculateHydrophobicity(hydrophilicCore + new string('I', isoleucines));
            gravy.Should().BeGreaterThan(previous, because: $"adding {isoleucines} isoleucines (4.5, above the mean) raises GRAVY");
            previous = gravy;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-MW-001 — protein molecular weight (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 124.
    //
    // API under test (SequenceStatistics.CalculateMolecularWeight):
    //   MW = Σ(residue masses) − (n − 1)·water (one water lost per peptide bond).
    //
    // Relations (derived from the peptide-bond model, NOT from output):
    //   • ADD  (MW(a+b) = MW(a) + MW(b) − water): concatenation forms one extra peptide bond,
    //          releasing exactly one water.
    //   • INV  (permutation invariant): MW is a function of the residue multiset only.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-MW-001 ADD — concatenation loses one water

    [Test]
    [Description("ADD: concatenating two peptides forms one extra peptide bond, so MW(a+b) = MW(a) + MW(b) − water.")]
    public void MolecularWeight_Concatenation_LosesOneWater()
    {
        // Derive the water mass from the API itself: MW(\"AA\") = 2·MW(\"A\") − water.
        double water = 2 * SequenceStatistics.CalculateMolecularWeight("A") - SequenceStatistics.CalculateMolecularWeight("AA");
        water.Should().BeGreaterThan(0, because: "a peptide bond releases a positive water mass");

        foreach (var (a, b) in new[] { ("MKL", "VAG"), ("ACDEF", "GHIKL"), ("W", "PQRST") })
        {
            double expected = SequenceStatistics.CalculateMolecularWeight(a) + SequenceStatistics.CalculateMolecularWeight(b) - water;
            SequenceStatistics.CalculateMolecularWeight(a + b)
                .Should().BeApproximately(expected, 1e-6,
                    because: $"joining '{a}' and '{b}' forms one peptide bond, removing one water");
        }
    }

    #endregion

    #region SEQ-MW-001 INV — permutation invariant

    [Test]
    [Description("INV: MW depends only on the residue multiset, so reordering the residues leaves it unchanged.")]
    public void MolecularWeight_Permutation_Invariant()
    {
        const string protein = "MKLVAGWTYSDE";
        SequenceStatistics.CalculateMolecularWeight(new string(protein.Reverse().ToArray()))
            .Should().BeApproximately(SequenceStatistics.CalculateMolecularWeight(protein), 1e-9,
                because: "MW sums residue masses and removes (n−1) waters — both invariant to order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-PI-001 — isoelectric point (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 125.
    //
    // API under test (SequenceStatistics.CalculateIsoelectricPoint):
    //   pI = the pH at which the net charge from the ionisable groups is zero — a function of
    //   the charged-residue counts.
    //
    // Relations (derived from the charge model, NOT from output):
    //   • INV  (permutation invariant): pI depends only on the residue counts, not their order.
    //   • MON  (more acidic residues ⇒ lower pI): adding acidic residues (Asp) shifts the
    //          net-charge curve down, lowering the pI.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-PI-001 INV — permutation invariant

    [Test]
    [Description("INV: pI is a function of the charged-residue counts, so reordering the residues leaves it unchanged.")]
    public void IsoelectricPoint_Permutation_Invariant()
    {
        const string protein = "MKLVAGWTYSDERH";
        SequenceStatistics.CalculateIsoelectricPoint(new string(protein.Reverse().ToArray()))
            .Should().BeApproximately(SequenceStatistics.CalculateIsoelectricPoint(protein), 1e-9,
                because: "pI depends only on the counts of ionisable residues, which a permutation preserves");
    }

    #endregion

    #region SEQ-PI-001 MON — more acidic residues lower the pI

    [Test]
    [Description("MON: adding acidic residues (Asp) shifts the net-charge curve down, monotonically lowering the isoelectric point.")]
    public void IsoelectricPoint_MoreAcidicResidues_LowerPi()
    {
        const string basicCore = "KRKR"; // basic residues → high pI
        double previous = double.MaxValue;
        foreach (int aspartates in new[] { 0, 1, 3, 6 })
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(basicCore + new string('D', aspartates));
            pi.Should().BeLessThan(previous, because: $"adding {aspartates} acidic Asp residues lowers the pI");
            previous = pi;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-SECSTRUCT-001 — secondary-structure propensity profile (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 126.
    //
    // API under test (SequenceStatistics.PredictSecondaryStructure):
    //   Sliding-window mean Chou-Fasman (helix, sheet, turn) propensities; one tuple per window.
    //
    // Relations (derived from the sliding window, NOT from output):
    //   • SHIFT (prepend flank shifts assignments): windows entirely within the original sequence
    //          appear as the suffix of the combined profile, shifted by the flank length.
    //   • INV  (deterministic): the profile is a pure function of the sequence.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-SECSTRUCT-001 SHIFT — a prepended flank shifts the window profile

    [Test]
    [Description("SHIFT: windows entirely within the original sequence appear unchanged as the suffix of the flank+seq profile, shifted by the flank length.")]
    public void SecondaryStructure_PrependFlank_ShiftsProfile()
    {
        const string seq = "MKLVAGWTYSDERHF";
        const int window = 7;
        var seqProfile = SequenceStatistics.PredictSecondaryStructure(seq, window).ToList();
        seqProfile.Should().NotBeEmpty();

        foreach (var flank in new[] { "AAAA", "WWWWWW" })
        {
            var combined = SequenceStatistics.PredictSecondaryStructure(flank + seq, window).ToList();
            combined.Skip(flank.Length).Take(seqProfile.Count)
                .Should().Equal(seqProfile,
                    because: $"windows fully inside the original sequence are unchanged, only shifted by the {flank.Length}-residue flank");
        }
    }

    #endregion

    #region SEQ-SECSTRUCT-001 INV — the profile is deterministic

    [Test]
    [Description("INV: PredictSecondaryStructure is a pure function — repeated calls give the identical profile.")]
    public void SecondaryStructure_SameSequence_SameProfile()
    {
        const string seq = "MKLVAGWTYSDERHF";
        SequenceStatistics.PredictSecondaryStructure(seq, 7).ToList()
            .Should().Equal(SequenceStatistics.PredictSecondaryStructure(seq, 7).ToList(),
                because: "the propensity profile has no hidden state");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-STATS-001 — sequence composition statistics (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 127.
    //
    // API under test (SequenceStatistics.CalculateNucleotideComposition):
    //   Per-base counts (A/T/G/C/U/N/Other) and the derived GC%, AT%, GC-skew and AT-skew.
    //   Every field is a function of the per-base counts alone.
    //
    // Relations (derived from count-based composition, NOT from output):
    //   • INV  (permutation invariant): the counts depend only on the multiset of bases, so any
    //          reordering — here a deterministic sort — leaves every count and derived statistic
    //          unchanged.
    //   • P    (concatenation sums counts): counting is additive over concatenation, so each
    //          per-base count and the length of a+b equal the sum of the parts' counts.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-STATS-001 INV — composition is permutation invariant

    [Test]
    [Description("INV: every composition statistic depends only on the multiset of bases, so a sort-permutation of a mixed DNA/RNA sequence (with U and N) leaves the full result unchanged.")]
    public void NucleotideStats_Permutation_Invariant()
    {
        const string seq = "AACGTUACGTTNGGCCAAUNAC";
        string sorted = new string(seq.OrderBy(c => c).ToArray());

        SequenceStatistics.CalculateNucleotideComposition(sorted)
            .Should().Be(SequenceStatistics.CalculateNucleotideComposition(seq),
                because: "counts, GC%/AT% and the skews are all functions of the per-base counts, which any permutation preserves");
    }

    #endregion

    #region SEQ-STATS-001 P — concatenation sums the per-base counts

    [Test]
    [Description("P: base counting is additive over concatenation, so CountA/T/G/C/U/N/Other and Length of a+b equal the sums of the parts'.")]
    public void NucleotideStats_Concatenation_SumsCounts()
    {
        foreach (var (a, b) in new[] { ("AACGT", "GGCCN"), ("AAUUGGCC", "TTTACGU"), ("N", "ACGTACGT") })
        {
            var ca = SequenceStatistics.CalculateNucleotideComposition(a);
            var cb = SequenceStatistics.CalculateNucleotideComposition(b);
            var cab = SequenceStatistics.CalculateNucleotideComposition(a + b);

            cab.Length.Should().Be(ca.Length + cb.Length, because: "length is additive over concatenation");
            cab.CountA.Should().Be(ca.CountA + cb.CountA, because: $"A count of '{a}'+'{b}' is the sum of the parts");
            cab.CountT.Should().Be(ca.CountT + cb.CountT, because: $"T count of '{a}'+'{b}' is the sum of the parts");
            cab.CountG.Should().Be(ca.CountG + cb.CountG, because: $"G count of '{a}'+'{b}' is the sum of the parts");
            cab.CountC.Should().Be(ca.CountC + cb.CountC, because: $"C count of '{a}'+'{b}' is the sum of the parts");
            cab.CountU.Should().Be(ca.CountU + cb.CountU, because: $"U count of '{a}'+'{b}' is the sum of the parts");
            cab.CountN.Should().Be(ca.CountN + cb.CountN, because: $"N count of '{a}'+'{b}' is the sum of the parts");
            cab.CountOther.Should().Be(ca.CountOther + cb.CountOther, because: $"Other count of '{a}'+'{b}' is the sum of the parts");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-SUMMARY-001 — aggregated sequence summary (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 128.
    //
    // API under test (SequenceStatistics.SummarizeNucleotideSequence):
    //   Aggregates Length, GcContent, Shannon Entropy, linguistic Complexity, melting
    //   temperature and the per-base Composition dictionary into one record.
    //
    // Relations (derived from the per-metric definitions, NOT from output):
    //   • INV  (permutation invariant for composition fields): Length, the Composition counts,
    //          GcContent, Entropy and Tm are all functions of the per-base counts only, so a
    //          permutation leaves them unchanged. (Complexity is order-dependent and is therefore
    //          NOT asserted invariant.)
    //   • SHIFT (length additive on concatenation): Length — and the per-base Composition counts —
    //          of a+b equal the sums of the parts'.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-SUMMARY-001 INV — composition fields are permutation invariant

    [Test]
    [Description("INV: the count-derived summary fields (Length, Composition, GcContent, Entropy, Tm) depend only on the base multiset, so a sort-permutation leaves them unchanged. Complexity is order-dependent and not asserted here.")]
    public void Summary_Permutation_PreservesCompositionFields()
    {
        const string seq = "AACGTACGTTGGCCAATACGT";
        string sorted = new string(seq.OrderBy(c => c).ToArray());

        var original = SequenceStatistics.SummarizeNucleotideSequence(seq);
        var permuted = SequenceStatistics.SummarizeNucleotideSequence(sorted);

        permuted.Length.Should().Be(original.Length, because: "length is unchanged by reordering");
        permuted.GcContent.Should().BeApproximately(original.GcContent, 1e-12, because: "GC% is a function of the G/C counts");
        permuted.Entropy.Should().BeApproximately(original.Entropy, 1e-12, because: "Shannon entropy is a function of the base frequencies");
        permuted.MeltingTemperature.Should().BeApproximately(original.MeltingTemperature, 1e-12,
            because: "Tm depends on the GC count and length, both permutation-invariant (length unchanged keeps the same formula branch)");
        permuted.Composition.Should().BeEquivalentTo(original.Composition, because: "per-base counts depend only on the multiset of bases");
    }

    #endregion

    #region SEQ-SUMMARY-001 SHIFT — length and composition are additive on concatenation

    [Test]
    [Description("SHIFT: Length and the per-base Composition counts are additive over concatenation, so the summary of a+b sums the parts' length and counts.")]
    public void Summary_Concatenation_AdditiveLengthAndCounts()
    {
        foreach (var (a, b) in new[] { ("AACGT", "GGCCN"), ("AAUUGGCC", "TTTACGU"), ("ACGTACGT", "N") })
        {
            var sa = SequenceStatistics.SummarizeNucleotideSequence(a);
            var sb = SequenceStatistics.SummarizeNucleotideSequence(b);
            var sab = SequenceStatistics.SummarizeNucleotideSequence(a + b);

            sab.Length.Should().Be(sa.Length + sb.Length, because: $"length of '{a}'+'{b}' is the sum of the parts");

            foreach (char baseChar in sab.Composition.Keys)
                sab.Composition[baseChar].Should().Be(
                    sa.Composition.GetValueOrDefault(baseChar) + sb.Composition.GetValueOrDefault(baseChar),
                    because: $"the {baseChar} count of '{a}'+'{b}' is the sum of the parts' counts");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-THERMO-001 — nearest-neighbour duplex thermodynamics (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 129.
    //
    // API under test (SequenceStatistics.CalculateThermodynamics):
    //   ΔH°, ΔS°, ΔG°₃₇ and Tm via the Allawi & SantaLucia (1997) / SantaLucia (1998) unified
    //   nearest-neighbour (NN) model: a sum of overlapping-dinucleotide stacking parameters plus
    //   helix-initiation terms at both termini.
    //
    // Relations (derived from the NN model, NOT from output):
    //   • MON  (more GC pairs ⇒ lower ΔG): GC stacks are more stabilising than AT stacks, so
    //          replacing AT steps with GC steps at fixed length makes ΔG°₃₇ progressively more
    //          negative (more stable).
    //   • INV  (reverse-complement symmetry): the NN parameters obey param(XY)=param(revcomp(XY))
    //          and initiation depends only on the GC/AT class of each terminus — both preserved by
    //          reverse complement — so a permutation to the complementary strand only relabels each
    //          NN context to its symmetric partner, leaving ΔH/ΔS/ΔG/Tm identical.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-THERMO-001 MON — more GC steps lower (stabilise) ΔG

    [Test]
    [Description("MON: replacing AT dinucleotide steps with GC steps at fixed length makes ΔG°₃₇ strictly more negative, since GC stacking is more stabilising than AT stacking.")]
    public void Thermodynamics_MoreGcSteps_LowerDeltaG()
    {
        // Fixed-length ladder: convert leading 'AT' units to 'GC' one at a time.
        double previous = double.MaxValue;
        for (int gcUnits = 0; gcUnits <= 6; gcUnits++)
        {
            string seq = string.Concat(Enumerable.Repeat("GC", gcUnits)) +
                         string.Concat(Enumerable.Repeat("AT", 6 - gcUnits));
            double dG = SequenceStatistics.CalculateThermodynamics(seq).DeltaG;
            dG.Should().BeLessThan(previous,
                because: $"replacing an AT step with a GC step (now {gcUnits} GC units) adds stacking stability, lowering ΔG");
            previous = dG;
        }
    }

    #endregion

    #region SEQ-THERMO-001 INV — reverse-complement leaves the thermodynamics unchanged

    [Test]
    [Description("INV: the NN parameters satisfy param(XY)=param(revcomp(XY)) and terminal initiation depends only on the GC/AT class of each end, so reverse-complementing the strand yields identical ΔH/ΔS/ΔG/Tm.")]
    public void Thermodynamics_ReverseComplement_Invariant()
    {
        foreach (var seq in new[] { "ATGCGATTACAGGCAT", "GCGCGCATAT", "ACGTACGTACGT", "TTTTAAAAGGGGCCCC" })
        {
            SequenceStatistics.CalculateThermodynamics(RevComp(seq))
                .Should().Be(SequenceStatistics.CalculateThermodynamics(seq),
                    because: $"reverse complement maps every NN step to its symmetric partner with identical parameters, so '{seq}' and its reverse complement are thermodynamically identical");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-TM-001 — simple melting temperature (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 130.
    //
    // API under test (SequenceStatistics.CalculateMeltingTemperature):
    //   Tm via the Wallace rule Tm = 2(A+T) + 4(G+C) for short oligos (length < 14) or the
    //   Marmur-Doty GC formula Tm = 81.5 + 0.41·%GC − 675/length otherwise.
    //
    // Relations (derived from both Tm formulas, NOT from output):
    //   • MON  (more GC ⇒ higher Tm): in both regimes GC pairs raise Tm — Wallace weights GC at
    //          4 vs AT at 2, and the Marmur-Doty term grows with the GC fraction — so replacing AT
    //          with GC at fixed length monotonically increases Tm.
    //   • INV  (case-insensitive): the implementation upper-cases its input, so lower-, upper- and
    //          mixed-case spellings of the same sequence give the same Tm.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-TM-001 MON — more GC raises Tm (both formula regimes)

    [Test]
    [Description("MON: at fixed length, replacing AT with GC raises Tm in both regimes — the Wallace rule (length < 14) and the Marmur-Doty GC formula (length ≥ 14).")]
    public void MeltingTemperature_MoreGc_HigherTm()
    {
        // (units, branch) pairs: 6 units → length 12 (Wallace), 10 units → length 20 (Marmur-Doty).
        foreach (int units in new[] { 6, 10 })
        {
            double previous = double.MinValue;
            for (int gcUnits = 0; gcUnits <= units; gcUnits++)
            {
                string seq = string.Concat(Enumerable.Repeat("GC", gcUnits)) +
                             string.Concat(Enumerable.Repeat("AT", units - gcUnits));
                double tm = SequenceStatistics.CalculateMeltingTemperature(seq);
                tm.Should().BeGreaterThan(previous,
                    because: $"a length-{seq.Length} oligo with {gcUnits} GC units melts higher than one with fewer");
                previous = tm;
            }
        }
    }

    #endregion

    #region SEQ-TM-001 INV — Tm is case-insensitive

    [Test]
    [Description("INV: the implementation upper-cases its input, so lower-, upper- and mixed-case spellings yield the same Tm in both length regimes.")]
    public void MeltingTemperature_CaseInsensitive_SameTm()
    {
        foreach (var seq in new[] { "acgtacgtag", "atgcgattacaggcatacgt" })
        {
            string upper = seq.ToUpperInvariant();
            string mixed = new string(seq.Select((c, i) => i % 2 == 0 ? char.ToUpperInvariant(c) : char.ToLowerInvariant(c)).ToArray());

            double expected = SequenceStatistics.CalculateMeltingTemperature(upper);
            SequenceStatistics.CalculateMeltingTemperature(seq).Should().Be(expected, because: "lower case is upper-cased before counting");
            SequenceStatistics.CalculateMeltingTemperature(mixed).Should().Be(expected, because: "mixed case is upper-cased before counting");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-CODON-FREQ-001 — codon usage frequencies (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 227.
    //
    // API under test (SequenceStatistics.CalculateCodonFrequencies):
    //   Reads consecutive non-overlapping triplets from the reading frame and reports
    //   frequency = count(codon) / total counted codons (Kazusa CUTG).
    //
    // Relations (derived from the count/total ratio, NOT from output):
    //   • INV   (codon-preserving shuffle): the table depends only on the codon multiset, so
    //           reordering whole codons leaves every frequency unchanged.
    //   • SCALE (triplicating the sequence): repeating a frame-aligned sequence k times multiplies
    //           every codon count and the total by k, leaving the frequencies unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region SEQ-CODON-FREQ-001 — Helpers

    // A frame-aligned coding sequence (Met-Lys-Leu-Gly-Phe-Glu-Arg-Pro), 24 nt = 8 codons.
    private const string CodonFreqSeq = "ATGAAACTAGGTTTTGAACGTCCC";

    private static System.Collections.Generic.List<string> SplitCodons(string seq)
    {
        var codons = new System.Collections.Generic.List<string>();
        for (int i = 0; i + 3 <= seq.Length; i += 3)
            codons.Add(seq.Substring(i, 3));
        return codons;
    }

    #endregion

    #region SEQ-CODON-FREQ-001 INV — reordering whole codons keeps the frequencies

    [Test]
    [Description("INV: the codon-frequency table depends only on the codon multiset, so permuting whole codons leaves every frequency unchanged.")]
    public void CodonFrequencies_CodonPreservingShuffle_Invariant()
    {
        var original = SequenceStatistics.CalculateCodonFrequencies(CodonFreqSeq);
        original.Values.Sum().Should().BeApproximately(1.0, 1e-12, because: "frequencies are count/total and partition the codons — a non-vacuous fixture");

        var codons = SplitCodons(CodonFreqSeq);
        var rng = new System.Random(20260620);
        for (int i = codons.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (codons[i], codons[j]) = (codons[j], codons[i]);
        }
        var shuffled = SequenceStatistics.CalculateCodonFrequencies(string.Concat(codons));

        shuffled.Should().BeEquivalentTo(original, because: "codon frequencies ignore the order of codons");
    }

    #endregion

    #region SEQ-CODON-FREQ-001 SCALE — repeating the sequence preserves the frequencies

    [Test]
    [Description("SCALE: repeating a frame-aligned sequence k times multiplies every codon count and the total by k, leaving the frequencies unchanged.")]
    public void CodonFrequencies_Triplication_PreservesFrequencies()
    {
        var original = SequenceStatistics.CalculateCodonFrequencies(CodonFreqSeq);

        foreach (int k in new[] { 2, 3, 5 })
        {
            string repeated = string.Concat(System.Linq.Enumerable.Repeat(CodonFreqSeq, k));
            SequenceStatistics.CalculateCodonFrequencies(repeated)
                .Should().BeEquivalentTo(original,
                    because: $"repeating the frame-aligned sequence {k}× scales every count and the total by {k}, a common factor that cancels");
        }
    }

    #endregion
}
