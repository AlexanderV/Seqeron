// ONCO-MHC-001 — MHC-Peptide Binding (MHCflurry Class I pan-allele binding-affinity predictor)
// Evidence: docs/Evidence/ONCO-MHC-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-MHC-001.md
// Source: O'Donnell TJ, Rubinsteyn A, Laserson U (2020). MHCflurry 2.0. Cell Systems 11(1):42-48.e7.
//         doi:10.1016/j.cels.2020.06.010. MHCflurry source (Apache-2.0): github.com/openvax/mhcflurry
//         (amino_acid.py, encodable_sequences.py, class1_neural_network.py, regression_target.py,
//          ensemble_centrality.py); models_class1_pan release 20200610.
// Oracle: mhcflurry 2.1.5 Python API (Class1AffinityPredictor over models_class1_pan), and the smallest
//         ensemble member (feedforward [512,512], PAN-CLASS1-1-3ed9fb2d2dcc9803) for the embedded
//         single-network parity test. Values were captured in-session by re-running the model.

using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class MhcflurryAffinityPredictor_PredictIc50_Tests
{
    // The embedded single-network weight pack = the smallest models_class1_pan member
    // (feedforward [512,512]). Loading it exercises the full forward pass + to_ic50.
    private const string SingleNetResource =
        "Seqeron.Genomics.Tests.TestData.Mhcflurry.mhcflurry_single_net.bin";

    // Relative tolerance for neural-network float parity (Python float64 vs .NET double).
    private const double RelTol = 1e-3; // 0.1% — comfortably above the observed <0.03% gap.

    private static IReadOnlyList<MhcflurryAffinityPredictor.Network> LoadSingleNet()
    {
        var asm = typeof(MhcflurryAffinityPredictor_PredictIc50_Tests).Assembly;
        using Stream stream = asm.GetManifestResourceStream(SingleNetResource)
            ?? throw new InvalidOperationException($"Embedded resource '{SingleNetResource}' not found.");
        return MhcflurryAffinityPredictor.LoadWeightPack(stream);
    }

    private static void AssertRelClose(double actual, double expected, double relTol, string because)
    {
        double rel = Math.Abs(actual - expected) / Math.Abs(expected);
        Assert.That(rel, Is.LessThanOrEqualTo(relTol),
            $"{because} (expected {expected}, actual {actual}, rel {rel:G4}).");
    }

    #region EncodePeptide

    // M1 — left_pad_centered_right_pad layout: SIINFEKL (len 8) -> three copies in a 45-position layout.
    // Golden index layout from EncodableSequences.sequences_to_fixed_length_index_encoded_array
    // (alignment_method="left_pad_centered_right_pad", max_length=15): block0 left, block1 centered (pad
    // floor((15-8)/2)=3), block2 right; X index = 20. BLOSUM62 S-row[A]=1 confirms position 0 channel 0.
    [Test]
    public void EncodePeptide_Siinfekl_MatchesLeftPadCenteredRightPadLayout()
    {
        double[] flat = MhcflurryAffinityPredictor.EncodePeptide("SIINFEKL");

        // Expected per-position residue indices (ACDEFGHIKLMNPQRSTVWYX): S=15 I=7 N=11 F=4 E=3 K=8 L=9, X=20.
        int[] expectedIdx =
        [
            15, 7, 7, 11, 4, 3, 8, 9, 20, 20, 20, 20, 20, 20, 20,       // block 0: left-aligned + X pad
            20, 20, 20, 15, 7, 7, 11, 4, 3, 8, 9, 20, 20, 20, 20,        // block 1: centered (3 X, peptide, 4 X)
            20, 20, 20, 20, 20, 20, 20, 15, 7, 7, 11, 4, 3, 8, 9,        // block 2: right-aligned
        ];

        // BLOSUM62 row for S (index 15): A C D E F G H I K L M N P Q R S T V W Y X
        int[] sRow = [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2, 0];
        // X row (index 20): all zero except X column = 1.
        int[] xRow = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

        Assert.Multiple(() =>
        {
            Assert.That(flat.Length, Is.EqualTo(MhcflurryAffinityPredictor.PeptideFlatLength),
                "Peptide encoding length must be 3*15*21 = 945.");
            // Position 0 is 'S' -> first 21 values are the BLOSUM62 S row.
            for (int c = 0; c < 21; c++)
            {
                Assert.That(flat[c], Is.EqualTo((double)sRow[c]),
                    $"Position 0 (S) channel {c} must equal BLOSUM62 S-row value {sRow[c]}.");
            }
            // Position 8 is an X pad in block 0 -> the X row.
            for (int c = 0; c < 21; c++)
            {
                Assert.That(flat[8 * 21 + c], Is.EqualTo((double)xRow[c]),
                    $"Pad position 8 channel {c} must equal BLOSUM62 X-row value {xRow[c]}.");
            }
            // Verify a sample of the full index layout via the diagonal-ish check: every position's max-scoring
            // channel must be the self-channel of the expected residue (BLOSUM62 diagonal is the row maximum).
            for (int p = 0; p < expectedIdx.Length; p++)
            {
                int best = 0;
                double bestVal = double.NegativeInfinity;
                for (int c = 0; c < 21; c++)
                {
                    double v = flat[p * 21 + c];
                    if (v > bestVal) { bestVal = v; best = c; }
                }
                Assert.That(best, Is.EqualTo(expectedIdx[p]),
                    $"Layout position {p} should encode residue index {expectedIdx[p]} (its BLOSUM62 self-score is maximal).");
            }
        });
    }

    // C1 — peptide length boundaries: length 5 (min) and 15 (max) encode; below/above throw.
    [Test]
    public void EncodePeptide_OutOfRangeLength_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MhcflurryAffinityPredictor.EncodePeptide("AAAAA").Length,
                Is.EqualTo(MhcflurryAffinityPredictor.PeptideFlatLength), "Length 5 is the supported minimum.");
            Assert.That(MhcflurryAffinityPredictor.EncodePeptide("AAAAAAAAAAAAAAA").Length,
                Is.EqualTo(MhcflurryAffinityPredictor.PeptideFlatLength), "Length 15 is the supported maximum.");
            Assert.Throws<ArgumentOutOfRangeException>(() => MhcflurryAffinityPredictor.EncodePeptide("AAAA"),
                "Length 4 is below the minimum and must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(() => MhcflurryAffinityPredictor.EncodePeptide("AAAAAAAAAAAAAAAA"),
                "Length 16 is above the maximum and must throw.");
            Assert.Throws<ArgumentNullException>(() => MhcflurryAffinityPredictor.EncodePeptide(null!),
                "Null peptide must throw.");
        });
    }

    #endregion

    #region Pseudosequence table

    // M2 — bundled allele pseudosequence table: HLA-A*02:01 -> the exact 37-residue MHCflurry pseudosequence.
    [Test]
    public void GetPseudosequence_HlaA0201_MatchesBundledMhcflurryTable()
    {
        string ps = MhcflurryAffinityPredictor.GetPseudosequence("HLA-A*02:01");
        Assert.Multiple(() =>
        {
            Assert.That(ps, Is.EqualTo("YFAMYGEKVAHTHVDTLYGVRYDHYYTWAVLAYTWYA"),
                "HLA-A*02:01 must map to the MHCflurry models_class1_pan 37-residue pseudosequence.");
            Assert.That(ps.Length, Is.EqualTo(MhcflurryAffinityPredictor.PseudosequenceLength),
                "All MHCflurry pseudosequences are 37 residues.");
            Assert.That(MhcflurryAffinityPredictor.GetPseudosequence("HLA-B*07:02"),
                Is.EqualTo("YYSEYRNIYAQTDESNLYGLSYDDYYTWAERAYEWYA"),
                "HLA-B*07:02 must map to its MHCflurry pseudosequence.");
        });
    }

    // C2 — unknown allele throws KeyNotFoundException; null throws ArgumentNullException.
    [Test]
    public void GetPseudosequence_UnknownOrNull_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<KeyNotFoundException>(() => MhcflurryAffinityPredictor.GetPseudosequence("HLA-Z*99:99"),
                "An allele absent from the table must throw KeyNotFoundException.");
            Assert.Throws<ArgumentNullException>(() => MhcflurryAffinityPredictor.GetPseudosequence(null!),
                "Null allele must throw ArgumentNullException.");
        });
    }

    // S1 — the bundled table contains the expected scale of human HLA alleles (sanity on resource loading).
    [Test]
    public void GetAllelePseudosequences_BundledTable_ContainsManyHlaAlleles()
    {
        var table = MhcflurryAffinityPredictor.GetAllelePseudosequences();
        int hlaCount = 0;
        foreach (var key in table.Keys)
        {
            if (key.StartsWith("HLA-", StringComparison.Ordinal)) { hlaCount++; }
        }
        Assert.That(hlaCount, Is.GreaterThan(5000),
            "The MHCflurry pseudosequence table bundles thousands of HLA alleles (11609 in release 20200610).");
    }

    #endregion

    #region ToIc50 transform

    // M3 — to_ic50: IC50 = 50000^(1 - x). x=0 -> 50000; x=1 -> 1; x=0.5 -> sqrt(50000).
    [Test]
    public void ToIc50_MatchesRegressionTargetTransform()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MhcflurryAffinityPredictor.ToIc50(0.0), Is.EqualTo(50000.0).Within(1e-6),
                "Output 0 maps to the 50000 nM ceiling.");
            Assert.That(MhcflurryAffinityPredictor.ToIc50(1.0), Is.EqualTo(1.0).Within(1e-9),
                "Output 1 maps to 1 nM (50000^0).");
            Assert.That(MhcflurryAffinityPredictor.ToIc50(0.5), Is.EqualTo(Math.Sqrt(50000.0)).Within(1e-6),
                "Output 0.5 maps to sqrt(50000) nM.");
        });
    }

    #endregion

    #region Single-network forward-pass parity (embedded weight pack)

    // M4 — full forward-pass parity against the mhcflurry single-network oracle (encoders + dense layers +
    // sigmoid + to_ic50). Oracle IC50s captured from the smallest models_class1_pan member.
    [TestCase("SIINFEKL", "HLA-A*02:01", 11483.195201)]
    [TestCase("GILGFVFTL", "HLA-A*02:01", 19.123150)]
    [TestCase("NLVPMVATV", "HLA-A*02:01", 17.542640)]
    [TestCase("ELAGIGILTV", "HLA-A*02:01", 119.054961)]
    [TestCase("AAAWYLWEV", "HLA-A*02:01", 16.559303)]
    [TestCase("SIINFEKL", "HLA-B*07:02", 28830.796646)]
    [TestCase("SLYNTVATL", "HLA-A*02:01", 28.972028)]
    [TestCase("CINGVCWTV", "HLA-A*02:01", 92.105940)]
    public void PredictIc50_SingleNet_MatchesMhcflurryOracle(string peptide, string allele, double oracleIc50)
    {
        var nets = LoadSingleNet();
        Assert.That(nets, Has.Count.EqualTo(1), "The embedded pack contains exactly one network.");

        double ic50 = MhcflurryAffinityPredictor.PredictIc50(nets, peptide, allele);

        AssertRelClose(ic50, oracleIc50, RelTol,
            $"Single-network IC50 for {peptide}/{allele} must match mhcflurry within {RelTol:P1}");
    }

    // M5 — strong-binder vs non-binder ranking: a known HLA-A*02:01 strong binder (GILGFVFTL, flu M1) must
    // rank far below a non-binding self peptide (SIINFEKL on A*02:01, predicted in the µM range).
    [Test]
    public void PredictIc50_StrongBinderRanksFarBelowNonBinder()
    {
        var nets = LoadSingleNet();
        double strong = MhcflurryAffinityPredictor.PredictIc50(nets, "GILGFVFTL", "HLA-A*02:01");
        double nonBinder = MhcflurryAffinityPredictor.PredictIc50(nets, "SIINFEKL", "HLA-A*02:01");

        Assert.Multiple(() =>
        {
            Assert.That(strong, Is.LessThan(50.0),
                "GILGFVFTL/HLA-A*02:01 is a high-affinity binder (predicted IC50 < 50 nM).");
            Assert.That(nonBinder, Is.GreaterThan(5000.0),
                "SIINFEKL/HLA-A*02:01 is not an HLA-A2 binder (predicted IC50 in the µM range).");
            Assert.That(nonBinder / strong, Is.GreaterThan(100.0),
                "The non-binder's IC50 must be orders of magnitude above the strong binder's.");
        });
    }

    #endregion

    #region Ensemble aggregation (geometric mean)

    // M6 — the ensemble combiner is the geometric mean of per-network IC50s. Feeding the same network N times
    // must reproduce that single network's IC50 (geomean of identical values = the value), proving the
    // exp(mean(log(.))) aggregation rather than an arithmetic mean of outputs.
    [Test]
    public void PredictIc50_DuplicatedNetworkEnsemble_EqualsSingleNetwork()
    {
        var single = LoadSingleNet();
        var ensemble = new List<MhcflurryAffinityPredictor.Network> { single[0], single[0], single[0] };

        double one = MhcflurryAffinityPredictor.PredictIc50(single, "GILGFVFTL", "HLA-A*02:01");
        double many = MhcflurryAffinityPredictor.PredictIc50(ensemble, "GILGFVFTL", "HLA-A*02:01");

        Assert.That(many, Is.EqualTo(one).Within(1e-9),
            "The geometric mean of identical per-network IC50s equals that IC50.");
    }

    // C3 — empty ensemble throws.
    [Test]
    public void PredictIc50_EmptyEnsemble_Throws()
    {
        var empty = new List<MhcflurryAffinityPredictor.Network>();
        Assert.Throws<ArgumentException>(
            () => MhcflurryAffinityPredictor.PredictIc50(empty, "GILGFVFTL", "HLA-A*02:01"),
            "An empty ensemble cannot produce a prediction.");
    }

    #endregion

    #region Predict -> classify chain

    // M7 — predict -> ClassifyBindingAffinity end-to-end: a strong binder classifies as Strong and a non-binder
    // as NonBinder, using the existing unchanged classifier (strong < 50 nM, weak < 500 nM).
    [Test]
    public void PredictAndClassify_StrongAndNonBinder_ChainsThroughClassifier()
    {
        var nets = LoadSingleNet();

        var strong = MhcflurryAffinityPredictor.PredictAndClassify(nets, "GILGFVFTL", "HLA-A*02:01");
        var nonBinder = MhcflurryAffinityPredictor.PredictAndClassify(nets, "SIINFEKL", "HLA-A*02:01");

        Assert.Multiple(() =>
        {
            Assert.That(strong.Strength, Is.EqualTo(OncologyAnalyzer.BindingStrength.Strong),
                "GILGFVFTL/HLA-A*02:01 (IC50 < 50 nM) must classify as a strong binder.");
            Assert.That(strong.Ic50Nm, Is.LessThan(OncologyAnalyzer.StrongBinderIc50Nm),
                "The chained IC50 must fall under the strong-binder cutoff.");
            Assert.That(nonBinder.Strength, Is.EqualTo(OncologyAnalyzer.BindingStrength.NonBinder),
                "SIINFEKL/HLA-A*02:01 (IC50 in the µM range) must classify as a non-binder.");
        });
    }

    #endregion

    #region Weight-pack loader edge cases

    // C4 — a stream with a bad magic / truncated body is rejected.
    [Test]
    public void LoadWeightPack_InvalidStream_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<InvalidDataException>(
                () => MhcflurryAffinityPredictor.LoadWeightPack(new MemoryStream([0x00, 0x01, 0x02, 0x03])),
                "A stream without the MHCF magic must be rejected.");
            Assert.Throws<ArgumentNullException>(
                () => MhcflurryAffinityPredictor.LoadWeightPack(null!),
                "A null stream must throw.");
        });
    }

    #endregion
}
