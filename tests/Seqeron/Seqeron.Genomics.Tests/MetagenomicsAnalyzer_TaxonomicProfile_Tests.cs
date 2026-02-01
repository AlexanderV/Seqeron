using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for META-PROF-001: Taxonomic Profile Generation.
/// Canonical method: <see cref="MetagenomicsAnalyzer.GenerateTaxonomicProfile"/>
/// </summary>
/// <remarks>
/// Evidence sources:
/// - Wikipedia (Metagenomics): Taxonomic profiling, relative abundance normalization
/// - Shannon (1948): Shannon diversity index formula
/// - Simpson (1949): Simpson concentration index formula
/// - Segata et al. (2012): MetaPhlAn paper on taxonomic profiling
/// </remarks>
[TestFixture]
public class MetagenomicsAnalyzer_TaxonomicProfile_Tests
{
    #region Test Data Helpers

    private static MetagenomicsAnalyzer.TaxonomicClassification CreateClassification(
        string readId,
        string kingdom,
        string phylum = "",
        string genus = "",
        string species = "")
    {
        return new MetagenomicsAnalyzer.TaxonomicClassification(
            ReadId: readId,
            Kingdom: kingdom,
            Phylum: phylum,
            Class: "",
            Order: "",
            Family: "",
            Genus: genus,
            Species: species,
            Confidence: 0.9,
            MatchedKmers: 100,
            TotalKmers: 110);
    }

    private static MetagenomicsAnalyzer.TaxonomicClassification CreateUnclassified(string readId)
    {
        return new MetagenomicsAnalyzer.TaxonomicClassification(
            ReadId: readId,
            Kingdom: "Unclassified",
            Phylum: "",
            Class: "",
            Order: "",
            Family: "",
            Genus: "",
            Species: "",
            Confidence: 0,
            MatchedKmers: 0,
            TotalKmers: 110);
    }

    #endregion

    #region Abundance Tests

    /// <summary>
    /// M1: Empty input returns TotalReads=0, ClassifiedReads=0, empty abundances.
    /// Evidence: Standard robustness requirement.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_EmptyInput_ReturnsEmptyProfile()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>();

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.Multiple(() =>
        {
            Assert.That(profile.TotalReads, Is.EqualTo(0));
            Assert.That(profile.ClassifiedReads, Is.EqualTo(0));
            Assert.That(profile.KingdomAbundance, Is.Empty);
            Assert.That(profile.PhylumAbundance, Is.Empty);
            Assert.That(profile.GenusAbundance, Is.Empty);
            Assert.That(profile.SpeciesAbundance, Is.Empty);
            Assert.That(profile.ShannonDiversity, Is.EqualTo(0));
            Assert.That(profile.SimpsonDiversity, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// M2: Single classified read produces abundance=1.0 for that taxon.
    /// Evidence: Abundance formula (count / total = 1/1 = 1.0).
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_SingleClassification_AbundanceIsOne()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("read1", "Bacteria", "Proteobacteria", "Escherichia", "coli")
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.Multiple(() =>
        {
            Assert.That(profile.TotalReads, Is.EqualTo(1));
            Assert.That(profile.ClassifiedReads, Is.EqualTo(1));
            Assert.That(profile.KingdomAbundance["Bacteria"], Is.EqualTo(1.0));
            Assert.That(profile.PhylumAbundance["Proteobacteria"], Is.EqualTo(1.0));
            Assert.That(profile.GenusAbundance["Escherichia"], Is.EqualTo(1.0));
            Assert.That(profile.SpeciesAbundance["coli"], Is.EqualTo(1.0));
        });
    }

    /// <summary>
    /// M8: Sum of kingdom abundances ≈ 1.0.
    /// Evidence: Normalization invariant (abundances are fractions).
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_MultipleClassifications_AbundancesSumToOne()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Bacteria", "P1", "G1", "sp2"),
            CreateClassification("r3", "Archaea", "P2", "G2", "sp3"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        double kingdomSum = profile.KingdomAbundance.Values.Sum();
        double speciesSum = profile.SpeciesAbundance.Values.Sum();

        Assert.Multiple(() =>
        {
            Assert.That(kingdomSum, Is.EqualTo(1.0).Within(0.001));
            Assert.That(speciesSum, Is.EqualTo(1.0).Within(0.001));
        });
    }

    /// <summary>
    /// M14: Multiple taxa with equal counts produce equal abundances.
    /// Evidence: Fair distribution under uniformity.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_EqualCounts_ProduceEqualAbundances()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Archaea", "P2", "G2", "sp2"),
            CreateClassification("r3", "Eukarya", "P3", "G3", "sp3"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        double expectedAbundance = 1.0 / 3.0;

        Assert.Multiple(() =>
        {
            Assert.That(profile.KingdomAbundance["Bacteria"], Is.EqualTo(expectedAbundance).Within(0.001));
            Assert.That(profile.KingdomAbundance["Archaea"], Is.EqualTo(expectedAbundance).Within(0.001));
            Assert.That(profile.KingdomAbundance["Eukarya"], Is.EqualTo(expectedAbundance).Within(0.001));
        });
    }

    /// <summary>
    /// M13: Empty rank values are filtered from rank-specific abundances.
    /// Evidence: Implementation note — clean output requirement.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_EmptyRankValues_FilteredFromAbundance()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "", "", ""),  // Only kingdom populated
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.Multiple(() =>
        {
            Assert.That(profile.KingdomAbundance, Has.Count.EqualTo(1));
            Assert.That(profile.PhylumAbundance, Is.Empty);
            Assert.That(profile.GenusAbundance, Is.Empty);
            Assert.That(profile.SpeciesAbundance, Is.Empty);
        });
    }

    #endregion

    #region Unclassified Handling Tests

    /// <summary>
    /// M3, M6: Unclassified reads are excluded from ClassifiedReads count.
    /// Evidence: MetaPhlAn documentation — unclassified excluded from profiling.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_UnclassifiedReads_ExcludedFromClassifiedCount()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("read1", "Bacteria", "Proteobacteria", "Escherichia", "coli"),
            CreateUnclassified("read2"),
            CreateUnclassified("read3"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.Multiple(() =>
        {
            Assert.That(profile.TotalReads, Is.EqualTo(3));
            Assert.That(profile.ClassifiedReads, Is.EqualTo(1));
        });
    }

    /// <summary>
    /// M4: Unclassified reads are excluded from abundance denominator.
    /// Evidence: MetaPhlAn/Wikipedia — abundance = count / classified_total.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_UnclassifiedReads_ExcludedFromAbundanceDenominator()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateUnclassified("r2"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        // Single classified read should have abundance 1.0 (not 0.5)
        Assert.That(profile.KingdomAbundance["Bacteria"], Is.EqualTo(1.0));
    }

    /// <summary>
    /// Edge case: All reads unclassified produces empty abundances.
    /// Evidence: Logical extension of M3/M4.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_AllUnclassified_ReturnsEmptyAbundances()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateUnclassified("r1"),
            CreateUnclassified("r2"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.Multiple(() =>
        {
            Assert.That(profile.TotalReads, Is.EqualTo(2));
            Assert.That(profile.ClassifiedReads, Is.EqualTo(0));
            Assert.That(profile.KingdomAbundance, Is.Empty);
        });
    }

    #endregion

    #region Read Count Invariants

    /// <summary>
    /// M5: TotalReads equals input classification count.
    /// Evidence: Basic counting invariant.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_TotalReads_EqualsInputCount()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Archaea", "P2", "G2", "sp2"),
            CreateUnclassified("r3"),
            CreateClassification("r4", "Bacteria", "P1", "G1", "sp3"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.That(profile.TotalReads, Is.EqualTo(4));
    }

    /// <summary>
    /// M7: ClassifiedReads is less than or equal to TotalReads.
    /// Evidence: Logical bound invariant.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_ClassifiedReads_LessThanOrEqualToTotalReads()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateUnclassified("r2"),
            CreateClassification("r3", "Archaea", "P2", "G2", "sp2"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.That(profile.ClassifiedReads, Is.LessThanOrEqualTo(profile.TotalReads));
    }

    #endregion

    #region Diversity Metrics Tests

    /// <summary>
    /// M11: Single species produces Shannon = 0.
    /// Evidence: Shannon (1948) — entropy is zero when outcome is certain.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_SingleSpecies_ShannonIsZero()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Bacteria", "P1", "G1", "sp1"),  // Same species
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.That(profile.ShannonDiversity, Is.EqualTo(0).Within(0.001));
    }

    /// <summary>
    /// M12: Single species produces Simpson = 1.0.
    /// Evidence: Simpson (1949) — concentration is 1 when only one species.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_SingleSpecies_SimpsonIsOne()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Bacteria", "P1", "G1", "sp1"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.That(profile.SimpsonDiversity, Is.EqualTo(1.0).Within(0.001));
    }

    /// <summary>
    /// M9: Shannon diversity is non-negative.
    /// Evidence: Shannon formula — -Σp·ln(p) where ln(p) ≤ 0 for p ≤ 1.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_ShannonDiversity_IsNonNegative()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Archaea", "P2", "G2", "sp2"),
            CreateClassification("r3", "Eukarya", "P3", "G3", "sp3"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.That(profile.ShannonDiversity, Is.GreaterThanOrEqualTo(0));
    }

    /// <summary>
    /// M10: Simpson diversity is in range [0, 1].
    /// Evidence: Simpson formula — Σp² where each p ∈ [0,1].
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_SimpsonDiversity_InZeroOneRange()
    {
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "Bacteria", "P1", "G1", "sp1"),
            CreateClassification("r2", "Archaea", "P2", "G2", "sp2"),
            CreateClassification("r3", "Eukarya", "P3", "G3", "sp3"),
            CreateClassification("r4", "Bacteria", "P1", "G1", "sp4"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        Assert.That(profile.SimpsonDiversity, Is.InRange(0, 1));
    }

    /// <summary>
    /// S1: Uniform distribution produces expected Shannon value.
    /// Evidence: Shannon theory — H = ln(n) for uniform distribution of n species.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_UniformDistribution_ExpectedShannonValue()
    {
        // 4 species, each appearing once → uniform distribution
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "K1", "P1", "G1", "sp1"),
            CreateClassification("r2", "K2", "P2", "G2", "sp2"),
            CreateClassification("r3", "K3", "P3", "G3", "sp3"),
            CreateClassification("r4", "K4", "P4", "G4", "sp4"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        // Expected Shannon = ln(4) ≈ 1.386
        double expectedShannon = Math.Log(4);
        Assert.That(profile.ShannonDiversity, Is.EqualTo(expectedShannon).Within(0.01));
    }

    /// <summary>
    /// S2: Uniform distribution produces expected Simpson value.
    /// Evidence: Simpson theory — D = 1/n for uniform distribution of n species.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_UniformDistribution_ExpectedSimpsonValue()
    {
        // 4 species, each appearing once → uniform distribution
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>
        {
            CreateClassification("r1", "K1", "P1", "G1", "sp1"),
            CreateClassification("r2", "K2", "P2", "G2", "sp2"),
            CreateClassification("r3", "K3", "P3", "G3", "sp3"),
            CreateClassification("r4", "K4", "P4", "G4", "sp4"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        // Expected Simpson = 4 × (1/4)² = 4 × 1/16 = 1/4 = 0.25
        double expectedSimpson = 0.25;
        Assert.That(profile.SimpsonDiversity, Is.EqualTo(expectedSimpson).Within(0.01));
    }

    /// <summary>
    /// S4: High skew produces low Shannon, high Simpson.
    /// Evidence: Diversity theory — dominance reduces Shannon, increases Simpson.
    /// </summary>
    [Test]
    public void GenerateTaxonomicProfile_HighSkew_LowShannonHighSimpson()
    {
        // Dominant species (9 reads) vs rare species (1 read)
        var classifications = new List<MetagenomicsAnalyzer.TaxonomicClassification>();
        for (int i = 0; i < 9; i++)
        {
            classifications.Add(CreateClassification($"r{i}", "Bacteria", "P1", "G1", "dominant"));
        }
        classifications.Add(CreateClassification("r9", "Archaea", "P2", "G2", "rare"));

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        // Uniform 2-species would have Shannon = ln(2) ≈ 0.693
        // Skewed should be much lower
        Assert.Multiple(() =>
        {
            Assert.That(profile.ShannonDiversity, Is.LessThan(0.5));
            Assert.That(profile.SimpsonDiversity, Is.GreaterThan(0.8));
        });
    }

    #endregion
}
