using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ProteinMotifFinderTests
{
    // Signal Peptide tests moved to ProteinMotifFinder_DomainPrediction_Tests.cs (PROTMOTIF-DOMAIN-001)

    #region Transmembrane Prediction Tests

    [Test]
    public void PredictTransmembraneHelices_HydrophobicStretch_FindsHelix()
    {
        // 20+ hydrophobic amino acids
        string protein = "AAAA" + new string('L', 22) + "EEEE";
        var helices = PredictTransmembraneHelices(protein, windowSize: 19, threshold: 1.0).ToList();

        Assert.That(helices, Has.Count.EqualTo(1), "Single TM helix expected");
        Assert.That(helices[0].Start, Is.EqualTo(0));
        Assert.That(helices[0].End, Is.EqualTo(29));
        Assert.That(helices[0].Score, Is.EqualTo(3.8).Within(0.01));
    }

    [Test]
    public void PredictTransmembraneHelices_NoHydrophobic_ReturnsEmpty()
    {
        string protein = new string('E', 50);
        var helices = PredictTransmembraneHelices(protein).ToList();

        Assert.That(helices, Is.Empty);
    }

    [Test]
    public void PredictTransmembraneHelices_MultipleTM_FindsAll()
    {
        // Multi-pass membrane protein: 3 TM segments separated by loops
        string loop = "EEEEEEEEEEE";
        string tm = new string('L', 22);
        string protein = tm + loop + tm + loop + tm;

        var helices = PredictTransmembraneHelices(protein, threshold: 1.0).ToList();

        Assert.That(helices, Has.Count.EqualTo(3), "Three TM helices expected");
        Assert.Multiple(() =>
        {
            Assert.That(helices[0].Start, Is.EqualTo(0));
            Assert.That(helices[1].Start, Is.EqualTo(26));
            Assert.That(helices[2].Start, Is.EqualTo(59));
        });
    }

    [Test]
    public void PredictTransmembraneHelices_ShortSequence_ReturnsEmpty()
    {
        var helices = PredictTransmembraneHelices("LLLLL").ToList();
        Assert.That(helices, Is.Empty);
    }

    #endregion

    #region Disorder Prediction Tests

    [Test]
    public void PredictDisorderedRegions_DisorderProne_FindsRegions()
    {
        // Disorder-promoting: P, E, K, S flanked by ordered residues
        string protein = "LLLLVVVV" + "PEKSPEKSPPEKSPEKS" + "LLLLVVVV";
        var regions = PredictDisorderedRegions(protein, threshold: 0.4).ToList();

        Assert.That(regions, Has.Count.EqualTo(1), "One disordered region expected");
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(22));
        Assert.That(regions[0].Score, Is.EqualTo(0.6548).Within(0.01));
    }

    [Test]
    public void PredictDisorderedRegions_Ordered_ReturnsEmpty()
    {
        // Order-promoting: I, V, L, W, F
        string protein = new string('I', 30) + new string('V', 30);
        var regions = PredictDisorderedRegions(protein, threshold: 0.6).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void PredictDisorderedRegions_ShortSequence_ReturnsEmpty()
    {
        var regions = PredictDisorderedRegions("PPPP").ToList();
        Assert.That(regions, Is.Empty);
    }

    #endregion

    // Coiled-coil prediction tests moved to the canonical file
    // ProteinMotifFinder_PredictCoiledCoils_Tests.cs (PROTMOTIF-CC-001). The previous tests here
    // encoded the now-removed fabricated position-weight table and are superseded.

    // Low-complexity-region tests moved to the canonical file
    // ProteinMotifFinder_FindLowComplexityRegions_Tests.cs (PROTMOTIF-LC-001). The previous tests here
    // encoded an invented "dominant single-AA frequency" rule (now removed in favour of the SEG
    // Wootton & Federhen complexity measure) and are superseded.

    // Domain Finding tests moved to ProteinMotifFinder_DomainPrediction_Tests.cs (PROTMOTIF-DOMAIN-001)

    // Case sensitivity tests moved to ProteinMotifFinder_DomainPrediction_Tests.cs (PROTMOTIF-DOMAIN-001)

    #region Integration Tests

    [Test]
    public void FullWorkflow_AnalyzeProtein()
    {
        // A sample protein with multiple features
        string protein = "MKRLLLLLLLLLLLLLLLLLLASAG" + // Signal peptide
                        "NFTAAAA" +                      // N-glycosylation
                        "SARK" +                         // PKC site
                        "RGDAAA" +                       // Cell attachment
                        new string('L', 22) +            // TM helix
                        "EEEEE";

        var motifs = FindCommonMotifs(protein).ToList();
        // Signal-peptide prediction is covered by PROTMOTIF-SP-001
        // (ProteinMotifFinder_PredictSignalPeptide_Tests.cs).

        Assert.Multiple(() =>
        {
            Assert.That(motifs, Has.Count.EqualTo(64), "Total motifs including NES, SUMO, glycosylation, PKC, RGD, leucine zipper");
            Assert.That(motifs.Any(m => m.MotifName == "RGD"), Is.True, "RGD cell attachment motif");
            Assert.That(motifs.Any(m => m.MotifName == "ASN_GLYCOSYLATION"), Is.True, "N-glycosylation site");
        });
    }

    [Test]
    public void FullWorkflow_LargeProtein()
    {
        // Generate a deterministic 500-residue protein (seed=42)
        var random = new Random(42);
        var aas = "ACDEFGHIKLMNPQRSTVWY";
        var protein = new string(Enumerable.Range(0, 500)
            .Select(_ => aas[random.Next(aas.Length)]).ToArray());

        var motifs = FindCommonMotifs(protein).ToList();
        var tm = PredictTransmembraneHelices(protein).ToList();
        var disorder = PredictDisorderedRegions(protein).ToList();
        var domains = FindDomains(protein).ToList();
        // Signal-peptide prediction is covered by PROTMOTIF-SP-001
        // (ProteinMotifFinder_PredictSignalPeptide_Tests.cs).

        Assert.Multiple(() =>
        {
            Assert.That(motifs, Has.Count.EqualTo(23), "Deterministic motif count for seed=42");
            Assert.That(tm, Is.Empty, "No TM helices in random sequence");
            Assert.That(disorder, Has.Count.EqualTo(23), "Disordered regions in random sequence");
            Assert.That(domains, Is.Empty, "No domains in random sequence");
        });
    }

    #endregion
}

