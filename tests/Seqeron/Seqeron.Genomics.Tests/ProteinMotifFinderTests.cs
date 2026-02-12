using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ProteinMotifFinderTests
{
    #region Signal Peptide Tests

    [Test]
    public void PredictSignalPeptide_ClassicSignal_PredictsSite()
    {
        // Signal peptide: M + positive (RK) + hydrophobic + small residues + cleavage
        string protein = "MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF";
        var signal = PredictSignalPeptide(protein);

        Assert.That(signal, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(signal!.Value.CleavagePosition, Is.EqualTo(25), "Cleavage after LASAG region");
            Assert.That(signal.Value.Score, Is.EqualTo(0.96).Within(0.01), "Combined N/H/C region score");
            Assert.That(signal.Value.Probability, Is.EqualTo(1.0).Within(0.01));
        });
    }

    [Test]
    public void PredictSignalPeptide_NoSignal_ReturnsNull()
    {
        // All charged, no hydrophobic region
        string protein = "EEEEEEEEEEKKKKKKKKKKDDDDDRRRRR";
        var signal = PredictSignalPeptide(protein);

        Assert.That(signal, Is.Null);
    }

    [Test]
    public void PredictSignalPeptide_ShortSequence_ReturnsNull()
    {
        string protein = "MKKLLLL";
        var signal = PredictSignalPeptide(protein);

        Assert.That(signal, Is.Null);
    }

    [Test]
    public void PredictSignalPeptide_ReturnsRegions()
    {
        string protein = "MKRLLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF";
        var signal = PredictSignalPeptide(protein);

        Assert.That(signal, Is.Not.Null, "Signal peptide must be detected");
        Assert.Multiple(() =>
        {
            Assert.That(signal!.Value.CleavagePosition, Is.EqualTo(26));
            Assert.That(signal.Value.NRegion, Is.EqualTo("MKRLL"));
            Assert.That(signal.Value.HRegion, Has.Length.EqualTo(16));
            Assert.That(signal.Value.CRegion, Is.EqualTo("LASAG"));
            Assert.That(signal.Value.Score, Is.EqualTo(0.96).Within(0.01));
        });
    }

    #endregion

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

    #region Coiled-Coil Prediction Tests

    [Test]
    public void PredictCoiledCoils_HeptadPattern_FindsCoil()
    {
        // Ideal coiled-coil: L at a,d positions, E at e, K at g
        // abcdefg pattern repeated 6 times = 42 residues
        string coiledCoil = "";
        for (int i = 0; i < 6; i++)
        {
            coiledCoil += "LAEALEK"; // L-A-E-A-L-E-K pattern
        }

        var coils = PredictCoiledCoils(coiledCoil, threshold: 0.3).ToList();

        Assert.That(coils, Has.Count.EqualTo(5), "Five overlapping coiled-coil windows");
        Assert.That(coils[0].Start, Is.EqualTo(0));
        Assert.That(coils[0].Score, Is.EqualTo(0.525).Within(0.01));
    }

    [Test]
    public void PredictCoiledCoils_NoPattern_ReturnsEmpty()
    {
        string protein = new string('P', 50); // Proline breaks helices
        var coils = PredictCoiledCoils(protein, threshold: 0.5).ToList();

        Assert.That(coils, Is.Empty);
    }

    [Test]
    public void PredictCoiledCoils_ShortSequence_ReturnsEmpty()
    {
        var coils = PredictCoiledCoils("LAELAE").ToList();
        Assert.That(coils, Is.Empty);
    }

    #endregion

    #region Low Complexity Tests

    [Test]
    public void FindLowComplexityRegions_PolyAlanine_Finds()
    {
        string protein = "MKKK" + new string('A', 15) + "VVVV";
        var regions = FindLowComplexityRegions(protein, threshold: 0.5).ToList();

        Assert.That(regions, Has.Count.EqualTo(1), "Single low-complexity region");
        Assert.That(regions[0].DominantAa, Is.EqualTo('A'));
        Assert.That(regions[0].Frequency, Is.EqualTo(1.0).Within(0.01), "Pure poly-A frequency");
    }

    [Test]
    public void FindLowComplexityRegions_Diverse_ReturnsEmpty()
    {
        string protein = "ACDEFGHIKLMNPQRSTVWY"; // All different
        var regions = FindLowComplexityRegions(protein, threshold: 0.5).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void FindLowComplexityRegions_MultipleRegions_FindsAll()
    {
        string protein = new string('G', 12) + "MKLVFP" + new string('S', 12);
        var regions = FindLowComplexityRegions(protein, threshold: 0.5).ToList();

        Assert.That(regions, Has.Count.EqualTo(2), "Two low-complexity regions: poly-G and poly-S");
        Assert.That(regions[0].DominantAa, Is.EqualTo('G'));
        Assert.That(regions[1].DominantAa, Is.EqualTo('S'));
    }

    #endregion

    #region Domain Finding Tests

    [Test]
    public void FindDomains_ZincFinger_Finds()
    {
        // C-x(2,4)-C-x(3)-L-x(8)-H-x(3,5)-H
        string protein = "AAAACXXCXXXLXXXXXXXXHXXXHAAA";
        var domains = FindDomains(protein).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Zinc Finger C2H2 domain expected");
        Assert.That(domains[0].Name, Is.EqualTo("Zinc Finger C2H2"));
        Assert.That(domains[0].Start, Is.EqualTo(4));
        Assert.That(domains[0].End, Is.EqualTo(24));
    }

    [Test]
    public void FindDomains_PLloop_FindsKinase()
    {
        // [AG]-x(4)-G-K-[ST]
        string protein = "AAAAGXXXXGKSAAAA";
        var domains = FindDomains(protein).ToList();

        Assert.That(domains, Has.Count.EqualTo(1), "Kinase ATP-binding domain expected");
        Assert.That(domains[0].Name, Is.EqualTo("Protein Kinase ATP-binding"));
        Assert.That(domains[0].Start, Is.EqualTo(4));
        Assert.That(domains[0].End, Is.EqualTo(11));
    }

    [Test]
    public void FindDomains_EmptySequence_ReturnsEmpty()
    {
        var domains = FindDomains("").ToList();
        Assert.That(domains, Is.Empty);
    }

    #endregion

    #region Case Sensitivity Tests

    [Test]
    public void PredictSignalPeptide_HandlesLowercase()
    {
        // Evidence: Protein sequences should be case-insensitive
        string proteinLower = "mkrllllllllllllllllllasagdddeeefff";
        string proteinUpper = "MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF";

        var signalLower = PredictSignalPeptide(proteinLower);
        var signalUpper = PredictSignalPeptide(proteinUpper);

        // Both should produce consistent results (both null or both not null)
        bool lowerHasSignal = signalLower != null;
        bool upperHasSignal = signalUpper != null;
        Assert.That(lowerHasSignal, Is.EqualTo(upperHasSignal),
            "Lowercase and uppercase should produce same signal peptide detection result");
    }

    #endregion

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
        var signal = PredictSignalPeptide(protein);

        Assert.Multiple(() =>
        {
            Assert.That(motifs, Has.Count.EqualTo(63), "Total motifs including NES, SUMO, glycosylation, PKC, RGD, leucine zipper");
            Assert.That(signal, Is.Not.Null);
            Assert.That(signal!.Value.CleavagePosition, Is.EqualTo(25));
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
        var signal = PredictSignalPeptide(protein);
        var tm = PredictTransmembraneHelices(protein).ToList();
        var disorder = PredictDisorderedRegions(protein).ToList();
        var domains = FindDomains(protein).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(motifs, Has.Count.EqualTo(21), "Deterministic motif count for seed=42");
            Assert.That(signal, Is.Not.Null, "Signal peptide detected in random sequence");
            Assert.That(tm, Is.Empty, "No TM helices in random sequence");
            Assert.That(disorder, Has.Count.EqualTo(23), "Disordered regions in random sequence");
            Assert.That(domains, Is.Empty, "No domains in random sequence");
        });
    }

    #endregion
}

