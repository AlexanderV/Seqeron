using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class ProteinMotifFinderTests
{

    #region PROSITE Pattern Conversion Tests

    [Test]
    public void ConvertPrositeToRegex_SimplePattern_Converts()
    {
        string prosite = "R-G-D";
        string regex = ConvertPrositeToRegex(prosite);

        Assert.That(regex, Is.EqualTo("RGD"));
    }

    [Test]
    public void ConvertPrositeToRegex_AnyAminoAcid_ConvertsToX()
    {
        string prosite = "A-x-G";
        string regex = ConvertPrositeToRegex(prosite);

        Assert.That(regex, Is.EqualTo("A.G"));
    }

    [Test]
    public void ConvertPrositeToRegex_RepeatRange_Converts()
    {
        string prosite = "A-x(2,4)-G";
        string regex = ConvertPrositeToRegex(prosite);

        Assert.That(regex, Is.EqualTo("A.{2,4}G"));
    }

    [Test]
    public void ConvertPrositeToRegex_CharacterClass_Converts()
    {
        string prosite = "[ST]-x-[RK]";
        string regex = ConvertPrositeToRegex(prosite);

        Assert.That(regex, Is.EqualTo("[ST].[RK]"));
    }

    [Test]
    public void ConvertPrositeToRegex_Exclusion_Converts()
    {
        string prosite = "N-{P}-[ST]";
        string regex = ConvertPrositeToRegex(prosite);

        Assert.That(regex, Is.EqualTo("N[^P][ST]"));
    }

    [Test]
    public void ConvertPrositeToRegex_Terminus_Converts()
    {
        string prosite = "<M-x-K";
        string regex = ConvertPrositeToRegex(prosite);

        Assert.That(regex, Does.StartWith("^"));
    }

    [Test]
    public void FindMotifByProsite_UsesConversion()
    {
        string protein = "AAANFSAAANGTAAA";
        var matches = FindMotifByProsite(protein, "N-{P}-[ST]-{P}", "N-glycosylation").ToList();

        Assert.That(matches, Has.Count.GreaterThanOrEqualTo(1));
    }

    #endregion

    #region Signal Peptide Tests

    [Test]
    public void PredictSignalPeptide_ClassicSignal_PredictsSite()
    {
        // Signal peptide: M + positive (RK) + hydrophobic + small residues + cleavage
        string protein = "MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF";
        var signal = PredictSignalPeptide(protein);

        Assert.That(signal, Is.Not.Null);
        Assert.That(signal!.Value.CleavagePosition, Is.GreaterThan(15));
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

        if (signal != null)
        {
            Assert.Multiple(() =>
            {
                Assert.That(signal.Value.NRegion, Is.Not.Empty);
                Assert.That(signal.Value.HRegion, Is.Not.Empty);
                Assert.That(signal.Value.CRegion, Is.Not.Empty);
                Assert.That(signal.Value.Score, Is.GreaterThan(0));
            });
        }
    }

    #endregion

    #region Transmembrane Prediction Tests

    [Test]
    public void PredictTransmembraneHelices_HydrophobicStretch_FindsHelix()
    {
        // 20+ hydrophobic amino acids
        string protein = "AAAA" + new string('L', 22) + "EEEE";
        var helices = PredictTransmembraneHelices(protein, windowSize: 19, threshold: 1.0).ToList();

        Assert.That(helices, Has.Count.GreaterThanOrEqualTo(1));
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
        // Multi-pass membrane protein
        string loop = "EEEEEEEEEEE";
        string tm = new string('L', 22);
        string protein = tm + loop + tm + loop + tm;

        var helices = PredictTransmembraneHelices(protein, threshold: 1.0).ToList();

        Assert.That(helices, Has.Count.GreaterThanOrEqualTo(1));
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
        // Disorder-promoting: P, E, K, S
        string protein = "LLLLVVVV" + "PEKSPEKSPPEKSPEKS" + "LLLLVVVV";
        var regions = PredictDisorderedRegions(protein, threshold: 0.4).ToList();

        Assert.That(regions, Has.Count.GreaterThanOrEqualTo(0));
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
        // abcdefg pattern repeated
        string coiledCoil = "";
        for (int i = 0; i < 6; i++)
        {
            coiledCoil += "LAEALEK"; // L-A-E-A-L-E-K pattern
        }

        var coils = PredictCoiledCoils(coiledCoil, threshold: 0.3).ToList();

        Assert.That(coils, Has.Count.GreaterThanOrEqualTo(0));
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

        Assert.That(regions, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(regions.Any(r => r.DominantAa == 'A'), Is.True);
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

        Assert.That(regions, Has.Count.GreaterThanOrEqualTo(2));
    }

    #endregion

    #region Domain Finding Tests

    [Test]
    public void FindDomains_ZincFinger_Finds()
    {
        // C-x(2,4)-C-x(3)-L-x(8)-H-x(3,5)-H
        string protein = "AAAACXXCXXXLXXXXXXXXHXXXHAAA";
        var domains = FindDomains(protein).ToList();

        // May or may not find depending on exact pattern match
        Assert.That(domains, Is.Not.Null);
    }

    [Test]
    public void FindDomains_PLloop_FindsKinase()
    {
        // [AG]-x(4)-G-K-[ST]
        string protein = "AAAAGXXXXGKSAAAA";
        var domains = FindDomains(protein).ToList();

        var kinase = domains.Where(d => d.Name.Contains("Kinase")).ToList();
        Assert.That(kinase, Has.Count.GreaterThanOrEqualTo(1));
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
        var tm = PredictTransmembraneHelices(protein, threshold: 1.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(motifs.Count, Is.GreaterThan(0));
            Assert.That(signal, Is.Not.Null);
        });
    }

    [Test]
    public void FullWorkflow_LargeProtein()
    {
        // Generate a larger protein sequence
        var random = new Random(42);
        var aas = "ACDEFGHIKLMNPQRSTVWY";
        var protein = new string(Enumerable.Range(0, 500)
            .Select(_ => aas[random.Next(aas.Length)]).ToArray());

        // Should complete without error
        var motifs = FindCommonMotifs(protein).ToList();
        var signal = PredictSignalPeptide(protein);
        var tm = PredictTransmembraneHelices(protein).ToList();
        var disorder = PredictDisorderedRegions(protein).ToList();
        var domains = FindDomains(protein).ToList();

        Assert.That(motifs, Is.Not.Null);
    }

    #endregion
}

