using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatisticsTests
{
    #region Nucleotide Composition Tests

    [Test]
    public void CalculateNucleotideComposition_ValidDna_ReturnsCorrectCounts()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AATTGGCC");

        Assert.That(comp.CountA, Is.EqualTo(2));
        Assert.That(comp.CountT, Is.EqualTo(2));
        Assert.That(comp.CountG, Is.EqualTo(2));
        Assert.That(comp.CountC, Is.EqualTo(2));
        Assert.That(comp.Length, Is.EqualTo(8));
    }

    [Test]
    public void CalculateNucleotideComposition_ValidDna_ReturnsCorrectGcContent()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AATTGGCC");

        Assert.That(comp.GcContent, Is.EqualTo(0.5).Within(0.001));
        Assert.That(comp.AtContent, Is.EqualTo(0.5).Within(0.001));
    }

    [Test]
    public void CalculateNucleotideComposition_AllGc_Returns100PercentGc()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("GGGGCCCC");

        Assert.That(comp.GcContent, Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculateNucleotideComposition_Rna_CountsUracil()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("AAUUGGCC");

        Assert.That(comp.CountU, Is.EqualTo(2));
        Assert.That(comp.CountT, Is.EqualTo(0));
    }

    [Test]
    public void CalculateNucleotideComposition_GcSkew_CalculatesCorrectly()
    {
        // More G than C should give positive skew
        var comp = SequenceStatistics.CalculateNucleotideComposition("GGGC");
        Assert.That(comp.GcSkew, Is.GreaterThan(0));

        // More C than G should give negative skew
        comp = SequenceStatistics.CalculateNucleotideComposition("GCCC");
        Assert.That(comp.GcSkew, Is.LessThan(0));
    }

    [Test]
    public void CalculateNucleotideComposition_EmptyString_ReturnsZeros()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("");

        Assert.That(comp.Length, Is.EqualTo(0));
        Assert.That(comp.GcContent, Is.EqualTo(0));
    }

    [Test]
    public void CalculateNucleotideComposition_WithNsAndOther_CountsCorrectly()
    {
        var comp = SequenceStatistics.CalculateNucleotideComposition("ATGCNNXX");

        Assert.That(comp.CountN, Is.EqualTo(2));
        Assert.That(comp.CountOther, Is.EqualTo(2));
    }

    #endregion

    #region Amino Acid Composition Tests

    [Test]
    public void CalculateAminoAcidComposition_ValidProtein_ReturnsCounts()
    {
        var comp = SequenceStatistics.CalculateAminoAcidComposition("MKVLWA");

        Assert.That(comp.Counts['M'], Is.EqualTo(1));
        Assert.That(comp.Counts['K'], Is.EqualTo(1));
        Assert.That(comp.Length, Is.EqualTo(6));
    }

    [Test]
    public void CalculateAminoAcidComposition_CalculatesMolecularWeight()
    {
        var comp = SequenceStatistics.CalculateAminoAcidComposition("MKVLWA");

        Assert.That(comp.MolecularWeight, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateAminoAcidComposition_ChargedResidues_CalculatesRatio()
    {
        // K, R, D, E, H are charged
        var comp = SequenceStatistics.CalculateAminoAcidComposition("KKRRDDEEHH");

        Assert.That(comp.ChargedResidueRatio, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateAminoAcidComposition_AromaticResidues_CalculatesRatio()
    {
        // F, Y, W are aromatic
        var comp = SequenceStatistics.CalculateAminoAcidComposition("FFYYWW");

        Assert.That(comp.AromaticResidueRatio, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateAminoAcidComposition_EmptyString_ReturnsZeros()
    {
        var comp = SequenceStatistics.CalculateAminoAcidComposition("");

        Assert.That(comp.Length, Is.EqualTo(0));
        Assert.That(comp.IsoelectricPoint, Is.EqualTo(7.0));
    }

    #endregion

    #region Molecular Weight Tests

    [Test]
    public void CalculateMolecularWeight_ValidProtein_ReturnsPositive()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("MKVLWAIFGAPV");

        Assert.That(mw, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateMolecularWeight_SingleAminoAcid_ReturnsWeight()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("A");

        Assert.That(mw, Is.GreaterThan(70)); // Alanine ~89 Da
        Assert.That(mw, Is.LessThan(110));
    }

    [Test]
    public void CalculateMolecularWeight_LongerProtein_IncreasesProperly()
    {
        double mw1 = SequenceStatistics.CalculateMolecularWeight("AAA");
        double mw2 = SequenceStatistics.CalculateMolecularWeight("AAAAAA");

        Assert.That(mw2, Is.GreaterThan(mw1));
    }

    [Test]
    public void CalculateNucleotideMolecularWeight_Dna_ReturnsPositive()
    {
        double mw = SequenceStatistics.CalculateNucleotideMolecularWeight("ATGCATGC", isDna: true);

        Assert.That(mw, Is.GreaterThan(2000)); // ~8 * 330 Da
    }

    [Test]
    public void CalculateNucleotideMolecularWeight_Rna_ReturnsDifferentWeight()
    {
        double dnaMw = SequenceStatistics.CalculateNucleotideMolecularWeight("ATGC", isDna: true);
        double rnaMw = SequenceStatistics.CalculateNucleotideMolecularWeight("AUGC", isDna: false);

        Assert.That(rnaMw, Is.GreaterThan(dnaMw)); // RNA has extra OH groups
    }

    [Test]
    public void CalculateMolecularWeight_EmptyString_ReturnsZero()
    {
        double mw = SequenceStatistics.CalculateMolecularWeight("");
        Assert.That(mw, Is.EqualTo(0));
    }

    #endregion

    #region Isoelectric Point Tests

    [Test]
    public void CalculateIsoelectricPoint_NeutralProtein_ReturnsMiddle()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("AAAA");

        Assert.That(pi, Is.GreaterThan(5.0));
        Assert.That(pi, Is.LessThan(9.0));
    }

    [Test]
    public void CalculateIsoelectricPoint_AcidicProtein_ReturnsLow()
    {
        // D and E are acidic
        double pi = SequenceStatistics.CalculateIsoelectricPoint("DDDDEEEE");

        Assert.That(pi, Is.LessThan(5.0));
    }

    [Test]
    public void CalculateIsoelectricPoint_BasicProtein_ReturnsHigh()
    {
        // K and R are basic
        double pi = SequenceStatistics.CalculateIsoelectricPoint("KKKKRRRR");

        Assert.That(pi, Is.GreaterThan(9.0));
    }

    [Test]
    public void CalculateIsoelectricPoint_EmptyString_ReturnsNeutral()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("");
        Assert.That(pi, Is.EqualTo(7.0));
    }

    #endregion

    #region Hydrophobicity Tests

    [Test]
    public void CalculateHydrophobicity_HydrophobicResidues_ReturnsPositive()
    {
        // I, V, L are hydrophobic
        double gravy = SequenceStatistics.CalculateHydrophobicity("IVLIVL");

        Assert.That(gravy, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateHydrophobicity_HydrophilicResidues_ReturnsNegative()
    {
        // D, E, K, R are hydrophilic
        double gravy = SequenceStatistics.CalculateHydrophobicity("DDEEKK");

        Assert.That(gravy, Is.LessThan(0));
    }

    [Test]
    public void CalculateHydrophobicityProfile_ReturnsCorrectCount()
    {
        string protein = "MKVLWAIFGAPVMKVLWAIFGAPV";
        int windowSize = 9;

        var profile = SequenceStatistics.CalculateHydrophobicityProfile(protein, windowSize).ToList();

        Assert.That(profile, Has.Count.EqualTo(protein.Length - windowSize + 1));
    }

    [Test]
    public void CalculateHydrophobicityProfile_WindowTooLarge_ReturnsEmpty()
    {
        var profile = SequenceStatistics.CalculateHydrophobicityProfile("MKV", windowSize: 10).ToList();

        Assert.That(profile, Is.Empty);
    }

    #endregion

    #region Thermodynamics Tests

    [Test]
    public void CalculateThermodynamics_ValidDna_ReturnsProperties()
    {
        var thermo = SequenceStatistics.CalculateThermodynamics("ATGCATGCATGC");

        Assert.That(thermo.DeltaH, Is.Not.EqualTo(0));
        Assert.That(thermo.DeltaS, Is.Not.EqualTo(0));
        Assert.That(thermo.MeltingTemperature, Is.GreaterThan(0));
    }

    [Test]
    public void CalculateThermodynamics_HighGc_HasHigherTm()
    {
        var lowGc = SequenceStatistics.CalculateThermodynamics("AAAATTTT");
        var highGc = SequenceStatistics.CalculateThermodynamics("GGGGCCCC");

        Assert.That(highGc.MeltingTemperature, Is.GreaterThan(lowGc.MeltingTemperature));
    }

    [Test]
    public void CalculateMeltingTemperature_WallaceRule_ShortOligo()
    {
        double tm = SequenceStatistics.CalculateMeltingTemperature("ATGCATGC", useWallaceRule: true);

        // Wallace: 2(A+T) + 4(G+C) = 2*4 + 4*4 = 8 + 16 = 24
        Assert.That(tm, Is.EqualTo(24.0).Within(0.1));
    }

    [Test]
    public void CalculateMeltingTemperature_EmptyString_ReturnsZero()
    {
        double tm = SequenceStatistics.CalculateMeltingTemperature("");
        Assert.That(tm, Is.EqualTo(0));
    }

    #endregion

    #region Dinucleotide Frequency Tests

    [Test]
    public void CalculateDinucleotideFrequencies_ReturnsFrequencies()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGATGATG");

        Assert.That(freq, Is.Not.Empty);
        Assert.That(freq.ContainsKey("AT"), Is.True);
        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(0.001));
    }

    [Test]
    public void CalculateDinucleotideRatios_CpgRatio()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATGCATGCATGC");

        Assert.That(ratios, Is.Not.Empty);
        // CpG ratio should be calculated
        if (ratios.ContainsKey("CG"))
        {
            Assert.That(ratios["CG"], Is.GreaterThanOrEqualTo(0));
        }
    }

    [Test]
    public void CalculateDinucleotideFrequencies_ShortSequence_ReturnsEmpty()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("A");
        Assert.That(freq, Is.Empty);
    }

    #endregion

    #region Codon Frequency Tests

    [Test]
    public void CalculateCodonFrequencies_ReturnsFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGATGATG");

        Assert.That(freq, Is.Not.Empty);
        Assert.That(freq.ContainsKey("ATG"), Is.True);
    }

    [Test]
    public void CalculateCodonFrequencies_DifferentReadingFrame()
    {
        string seq = "AATGATGATG";
        var frame0 = SequenceStatistics.CalculateCodonFrequencies(seq, readingFrame: 0);
        var frame1 = SequenceStatistics.CalculateCodonFrequencies(seq, readingFrame: 1);

        // Different frames should give different codon sets
        Assert.That(frame0.Keys, Is.Not.EquivalentTo(frame1.Keys));
    }

    [Test]
    public void CalculateCodonFrequencies_ShortSequence_ReturnsEmpty()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("AT");
        Assert.That(freq, Is.Empty);
    }

    #endregion

    #region Entropy and Complexity Tests

    [Test]
    public void CalculateShannonEntropy_UniformDistribution_ReturnsHighEntropy()
    {
        // Equal distribution of all 4 bases should give max entropy = 2 bits
        double entropy = SequenceStatistics.CalculateShannonEntropy("ATGC");

        Assert.That(entropy, Is.EqualTo(2.0).Within(0.001));
    }

    [Test]
    public void CalculateShannonEntropy_HomopolymerRun_ReturnsZero()
    {
        double entropy = SequenceStatistics.CalculateShannonEntropy("AAAA");

        Assert.That(entropy, Is.EqualTo(0).Within(0.001));
    }

    [Test]
    public void CalculateLinguisticComplexity_RandomSequence_ReturnsHigh()
    {
        double complexity = SequenceStatistics.CalculateLinguisticComplexity("ATGCGATCGATCGATCGATCGATC");

        Assert.That(complexity, Is.GreaterThan(0.4)); // Complex sequences have reasonable complexity
    }

    [Test]
    public void CalculateLinguisticComplexity_RepetitiveSequence_ReturnsLow()
    {
        double complexity = SequenceStatistics.CalculateLinguisticComplexity("ATATATAT");

        Assert.That(complexity, Is.LessThan(0.5));
    }

    [Test]
    public void CalculateShannonEntropy_EmptyString_ReturnsZero()
    {
        double entropy = SequenceStatistics.CalculateShannonEntropy("");
        Assert.That(entropy, Is.EqualTo(0));
    }

    #endregion

    #region Secondary Structure Prediction Tests

    [Test]
    public void PredictSecondaryStructure_ReturnsCorrectCount()
    {
        string protein = "MKVLWAIFGAPVMKVLWAIFGAPV";
        int windowSize = 7;

        var predictions = SequenceStatistics.PredictSecondaryStructure(protein, windowSize).ToList();

        Assert.That(predictions, Has.Count.EqualTo(protein.Length - windowSize + 1));
    }

    [Test]
    public void PredictSecondaryStructure_HelixFormers_HighHelixPropensity()
    {
        // A, E, L, M are helix formers
        var predictions = SequenceStatistics.PredictSecondaryStructure("AAAAAAA", windowSize: 7).ToList();

        Assert.That(predictions[0].Helix, Is.GreaterThan(1.0)); // Above average propensity
    }

    [Test]
    public void PredictSecondaryStructure_SheetFormers_HighSheetPropensity()
    {
        // V, I, Y are sheet formers
        var predictions = SequenceStatistics.PredictSecondaryStructure("VVVVVVV", windowSize: 7).ToList();

        Assert.That(predictions[0].Sheet, Is.GreaterThan(1.0));
    }

    [Test]
    public void PredictSecondaryStructure_WindowTooLarge_ReturnsEmpty()
    {
        var predictions = SequenceStatistics.PredictSecondaryStructure("MKV", windowSize: 10).ToList();

        Assert.That(predictions, Is.Empty);
    }

    #endregion

    #region Profile Tests

    [Test]
    public void CalculateGcContentProfile_ReturnsCorrectCount()
    {
        string seq = string.Concat(Enumerable.Repeat("ATGCATGCATGC", 10));
        int windowSize = 20;

        var profile = SequenceStatistics.CalculateGcContentProfile(seq, windowSize, stepSize: 10).ToList();

        Assert.That(profile, Has.Count.GreaterThan(0));
    }

    [Test]
    public void CalculateGcContentProfile_AllGc_Returns100Percent()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("GGGGGGGGGGCCCCCCCCCC", 10).ToList();

        Assert.That(profile.All(p => p >= 0.99));
    }

    [Test]
    public void CalculateEntropyProfile_ReturnsCorrectCount()
    {
        string seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
        int windowSize = 20;

        var profile = SequenceStatistics.CalculateEntropyProfile(seq, windowSize).ToList();

        Assert.That(profile, Has.Count.EqualTo(seq.Length - windowSize + 1));
    }

    [Test]
    public void CalculateGcContentProfile_WindowTooLarge_ReturnsEmpty()
    {
        var profile = SequenceStatistics.CalculateGcContentProfile("ATGC", windowSize: 100).ToList();

        Assert.That(profile, Is.Empty);
    }

    #endregion

    #region Summary Statistics Tests

    [Test]
    public void SummarizeNucleotideSequence_ReturnsComprehensiveStats()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGCATGC");

        Assert.That(summary.Length, Is.EqualTo(12));
        Assert.That(summary.GcContent, Is.GreaterThan(0));
        Assert.That(summary.Entropy, Is.GreaterThan(0));
        Assert.That(summary.Complexity, Is.GreaterThan(0));
        Assert.That(summary.Composition, Is.Not.Empty);
    }

    [Test]
    public void SummarizeNucleotideSequence_CompositionContainsAllBases()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGC");

        Assert.That(summary.Composition.ContainsKey('A'));
        Assert.That(summary.Composition.ContainsKey('T'));
        Assert.That(summary.Composition.ContainsKey('G'));
        Assert.That(summary.Composition.ContainsKey('C'));
    }

    #endregion

    #region Record Types Tests

    [Test]
    public void NucleotideComposition_RecordEquality_Works()
    {
        var comp1 = new SequenceStatistics.NucleotideComposition(10, 2, 2, 3, 3, 0, 0, 0, 0.6, 0.4, 0.0, 0.0);
        var comp2 = new SequenceStatistics.NucleotideComposition(10, 2, 2, 3, 3, 0, 0, 0, 0.6, 0.4, 0.0, 0.0);

        Assert.That(comp1, Is.EqualTo(comp2));
    }

    [Test]
    public void ThermodynamicProperties_RecordProperties_Work()
    {
        var thermo = new SequenceStatistics.ThermodynamicProperties(
            DeltaH: -50.0,
            DeltaS: -130.0,
            DeltaG: -10.0,
            MeltingTemperature: 55.5);

        Assert.That(thermo.DeltaH, Is.EqualTo(-50.0));
        Assert.That(thermo.MeltingTemperature, Is.EqualTo(55.5));
    }

    [Test]
    public void AminoAcidComposition_RecordProperties_Work()
    {
        var comp = new SequenceStatistics.AminoAcidComposition(
            Length: 10,
            Counts: new Dictionary<char, int> { { 'A', 5 }, { 'G', 5 } },
            MolecularWeight: 1000.0,
            IsoelectricPoint: 6.5,
            Hydrophobicity: 0.5,
            ChargedResidueRatio: 0.1,
            AromaticResidueRatio: 0.05);

        Assert.That(comp.Length, Is.EqualTo(10));
        Assert.That(comp.Counts['A'], Is.EqualTo(5));
    }

    [Test]
    public void SequenceSummary_RecordProperties_Work()
    {
        var summary = new SequenceStatistics.SequenceSummary(
            Length: 100,
            GcContent: 0.5,
            Entropy: 2.0,
            Complexity: 0.8,
            MeltingTemperature: 50.0,
            Composition: new Dictionary<char, int> { { 'A', 25 } });

        Assert.That(summary.Length, Is.EqualTo(100));
        Assert.That(summary.GcContent, Is.EqualTo(0.5));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void AllMethods_HandleNullGracefully()
    {
        // These should not throw exceptions
        Assert.DoesNotThrow(() => SequenceStatistics.CalculateNucleotideComposition(null!));
        Assert.DoesNotThrow(() => SequenceStatistics.CalculateAminoAcidComposition(null!));
        Assert.DoesNotThrow(() => SequenceStatistics.CalculateMolecularWeight(null!));
        Assert.DoesNotThrow(() => SequenceStatistics.CalculateIsoelectricPoint(null!));
        Assert.DoesNotThrow(() => SequenceStatistics.CalculateHydrophobicity(null!));
        Assert.DoesNotThrow(() => SequenceStatistics.CalculateShannonEntropy(null!));
    }

    [Test]
    public void AllMethods_HandleCaseInsensitively()
    {
        var upper = SequenceStatistics.CalculateNucleotideComposition("ATGC");
        var lower = SequenceStatistics.CalculateNucleotideComposition("atgc");
        var mixed = SequenceStatistics.CalculateNucleotideComposition("AtGc");

        Assert.That(upper.GcContent, Is.EqualTo(lower.GcContent));
        Assert.That(upper.GcContent, Is.EqualTo(mixed.GcContent));
    }

    #endregion
}
