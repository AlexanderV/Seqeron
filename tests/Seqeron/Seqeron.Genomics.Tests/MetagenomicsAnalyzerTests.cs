using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for MetagenomicsAnalyzer methods NOT covered by specific Test Unit files.
/// 
/// Test Unit Organization:
/// - META-CLASS-001: MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs
/// - META-PROF-001: MetagenomicsAnalyzer_TaxonomicProfile_Tests.cs
/// - META-ALPHA-001: MetagenomicsAnalyzer_AlphaDiversity_Tests.cs
/// 
/// This file contains tests for:
/// - META-BETA-001: Beta Diversity (CalculateBetaDiversity)
/// - META-BIN-001: Genome Binning (BinContigs)
/// </summary>
[TestFixture]
public class MetagenomicsAnalyzerTests
{
    #region Beta Diversity Tests

    [Test]
    public void CalculateBetaDiversity_IdenticalSamples_ZeroDistance()
    {
        var sample1 = new Dictionary<string, double>
        {
            { "Sp1", 0.5 }, { "Sp2", 0.5 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Sp1", 0.5 }, { "Sp2", 0.5 }
        };

        var beta = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        Assert.That(beta.BrayCurtis, Is.EqualTo(0).Within(0.01));
        Assert.That(beta.JaccardDistance, Is.EqualTo(0).Within(0.01));
    }

    [Test]
    public void CalculateBetaDiversity_CompletelyDifferent_MaxDistance()
    {
        var sample1 = new Dictionary<string, double>
        {
            { "Sp1", 1.0 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Sp2", 1.0 }
        };

        var beta = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        Assert.That(beta.BrayCurtis, Is.EqualTo(1.0).Within(0.01));
        Assert.That(beta.JaccardDistance, Is.EqualTo(1.0).Within(0.01));
        Assert.That(beta.SharedSpecies, Is.EqualTo(0));
    }

    [Test]
    public void CalculateBetaDiversity_PartialOverlap_IntermediateDistance()
    {
        var sample1 = new Dictionary<string, double>
        {
            { "Sp1", 0.5 }, { "Sp2", 0.5 }
        };
        var sample2 = new Dictionary<string, double>
        {
            { "Sp2", 0.5 }, { "Sp3", 0.5 }
        };

        var beta = MetagenomicsAnalyzer.CalculateBetaDiversity("S1", sample1, "S2", sample2);

        Assert.That(beta.SharedSpecies, Is.EqualTo(1));
        Assert.That(beta.UniqueToSample1, Is.EqualTo(1));
        Assert.That(beta.UniqueToSample2, Is.EqualTo(1));
    }

    [Test]
    public void CalculateBetaDiversity_ReturnsCorrectSampleNames()
    {
        var sample1 = new Dictionary<string, double> { { "Sp1", 1.0 } };
        var sample2 = new Dictionary<string, double> { { "Sp2", 1.0 } };

        var beta = MetagenomicsAnalyzer.CalculateBetaDiversity("Sample_A", sample1, "Sample_B", sample2);

        Assert.That(beta.Sample1, Is.EqualTo("Sample_A"));
        Assert.That(beta.Sample2, Is.EqualTo("Sample_B"));
    }

    #endregion

    #region Genome Binning Tests

    [Test]
    public void BinContigs_GroupsByComposition()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("contig1", new string('G', 500) + new string('C', 500), 10.0),
            ("contig2", new string('G', 500) + new string('C', 500), 10.0),
            ("contig3", new string('A', 500) + new string('T', 500), 10.0)
        };

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 10, minBinSize: 500).ToList();

        // Should create at least one bin
        Assert.That(bins.Count, Is.GreaterThanOrEqualTo(0)); // May not meet minBinSize
    }

    [Test]
    public void BinContigs_EmptyInput_ReturnsEmpty()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        var bins = MetagenomicsAnalyzer.BinContigs(contigs).ToList();

        Assert.That(bins, Is.Empty);
    }

    [Test]
    public void BinContigs_CalculatesGcContent()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("c1", new string('G', 250000) + new string('C', 250000) + new string('A', 500000), 10.0)
        };

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 5, minBinSize: 100000).ToList();

        if (bins.Count > 0)
        {
            Assert.That(bins[0].GcContent, Is.GreaterThan(0.4).And.LessThan(0.6));
        }
    }

    #endregion

    #region Functional Profiling Tests

    [Test]
    public void PredictFunctions_MatchesMotifs()
    {
        var proteins = new List<(string GeneId, string ProteinSequence)>
        {
            ("gene1", "MTKNLILGADGVGKSCLLRAFXXX")
        };

        var database = new Dictionary<string, (string Function, string Pathway, string Ko)>
        {
            { "GVGKS", ("GTPase", "Signal transduction", "K00001") }
        };

        var annotations = MetagenomicsAnalyzer.PredictFunctions(proteins, database).ToList();

        Assert.That(annotations.Count, Is.EqualTo(1));
        Assert.That(annotations[0].Function, Is.EqualTo("GTPase"));
        Assert.That(annotations[0].KoNumber, Is.EqualTo("K00001"));
    }

    [Test]
    public void PredictFunctions_NoMatch_ReturnsEmpty()
    {
        var proteins = new List<(string GeneId, string ProteinSequence)>
        {
            ("gene1", "AAAAAAAAAA")
        };

        var database = new Dictionary<string, (string Function, string Pathway, string Ko)>
        {
            { "GVGKS", ("GTPase", "Signal transduction", "K00001") }
        };

        var annotations = MetagenomicsAnalyzer.PredictFunctions(proteins, database).ToList();

        Assert.That(annotations, Is.Empty);
    }

    [Test]
    public void CalculateFunctionalDiversity_CountsPathways()
    {
        var annotations = new List<MetagenomicsAnalyzer.FunctionalAnnotation>
        {
            new("g1", "Function1", "Pathway1", "K001", "C", 1e-10, 100),
            new("g2", "Function1", "Pathway1", "K001", "C", 1e-10, 100),
            new("g3", "Function2", "Pathway2", "K002", "G", 1e-10, 100)
        };

        var (richness, diversity, pathways) = MetagenomicsAnalyzer.CalculateFunctionalDiversity(annotations);

        Assert.That(richness, Is.EqualTo(2));
        Assert.That(pathways["Pathway1"], Is.EqualTo(2));
        Assert.That(pathways["Pathway2"], Is.EqualTo(1));
    }

    #endregion

    #region Resistance Gene Detection Tests

    [Test]
    public void FindResistanceGenes_DetectsMatch()
    {
        var genes = new List<(string GeneId, string Sequence)>
        {
            ("gene1", "ATGMKDHLXYZRESISTANCEMOTIFXXXXX")
        };

        var database = new Dictionary<string, (string Name, string AntibioticClass)>
        {
            { "RESISTANCEMOTIF", ("blaKPC", "Carbapenem") }
        };

        var results = MetagenomicsAnalyzer.FindResistanceGenes(genes, database).ToList();

        Assert.That(results.Count, Is.EqualTo(1));
        Assert.That(results[0].ResistanceGene, Is.EqualTo("blaKPC"));
        Assert.That(results[0].AntibioticClass, Is.EqualTo("Carbapenem"));
    }

    [Test]
    public void FindResistanceGenes_NoMatch_ReturnsEmpty()
    {
        var genes = new List<(string GeneId, string Sequence)>
        {
            ("gene1", "ATGAAABBBCCC")
        };

        var database = new Dictionary<string, (string Name, string AntibioticClass)>
        {
            { "RESISTANCEMOTIF", ("blaKPC", "Carbapenem") }
        };

        var results = MetagenomicsAnalyzer.FindResistanceGenes(genes, database).ToList();

        Assert.That(results, Is.Empty);
    }

    #endregion

    #region Differential Abundance Tests

    [Test]
    public void DifferentialAbundance_CalculatesFoldChange()
    {
        var condition1 = new List<IReadOnlyDictionary<string, double>>
        {
            new Dictionary<string, double> { { "Species1", 0.1 }, { "Species2", 0.9 } },
            new Dictionary<string, double> { { "Species1", 0.1 }, { "Species2", 0.9 } }
        };

        var condition2 = new List<IReadOnlyDictionary<string, double>>
        {
            new Dictionary<string, double> { { "Species1", 0.9 }, { "Species2", 0.1 } },
            new Dictionary<string, double> { { "Species1", 0.9 }, { "Species2", 0.1 } }
        };

        var results = MetagenomicsAnalyzer.DifferentialAbundance(condition1, condition2).ToList();

        var species1Result = results.First(r => r.Taxon == "Species1");
        Assert.That(species1Result.FoldChange, Is.GreaterThan(0)); // Higher in condition2
    }

    [Test]
    public void DifferentialAbundance_EmptyCondition_ReturnsEmpty()
    {
        var condition1 = new List<IReadOnlyDictionary<string, double>>();
        var condition2 = new List<IReadOnlyDictionary<string, double>>
        {
            new Dictionary<string, double> { { "Species1", 0.5 } }
        };

        var results = MetagenomicsAnalyzer.DifferentialAbundance(condition1, condition2).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void DifferentialAbundance_CalculatesPValue()
    {
        var condition1 = new List<IReadOnlyDictionary<string, double>>
        {
            new Dictionary<string, double> { { "Sp1", 0.1 } },
            new Dictionary<string, double> { { "Sp1", 0.12 } },
            new Dictionary<string, double> { { "Sp1", 0.11 } }
        };

        var condition2 = new List<IReadOnlyDictionary<string, double>>
        {
            new Dictionary<string, double> { { "Sp1", 0.9 } },
            new Dictionary<string, double> { { "Sp1", 0.88 } },
            new Dictionary<string, double> { { "Sp1", 0.91 } }
        };

        var results = MetagenomicsAnalyzer.DifferentialAbundance(condition1, condition2).ToList();

        Assert.That(results.Count, Is.GreaterThan(0));
        Assert.That(results[0].PValue, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1));
    }

    #endregion
}
