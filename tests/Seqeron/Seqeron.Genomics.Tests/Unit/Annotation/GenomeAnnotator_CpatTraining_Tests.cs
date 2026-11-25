// ANNOT-CODING-001 - CPAT training: make_hexamer_tab, make_logitModel (R glm binomial), TG-ROC cutoff (A4)
// Evidence: docs/Evidence/ANNOT-CODING-001-Evidence.md (CPAT training oracle)
// TestSpec: tests/TestSpecs/ANNOT-CODING-001.md (8.3)
// Source: Wang L et al. (2013) NAR 41(6):e74; CPAT 3.0.5 make_hexamer_tab / make_logitModel (cpmodule/__init__.py,
//         FrameKmer.py, utils.make_logit = R glm(Label ~ mRNA + ORF + Fickett + Hexamer, binomial)); Greiner et al. 2000 (TG-ROC).
// The training sequences below are synthetic (seeded generator, no CPAT data). Expected values are the output of
// the real make_hexamer_tab / make_logitModel --stop TAG,TAA,TGA / cpat (3.0.5, R 4.x) run on the same files;
// R prints 17 significant digits (coefficients) / 15 (cpat probabilities).

namespace Seqeron.Genomics.Tests.Unit.Annotation;

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.GenomeAnnotator;

[TestFixture]
public class GenomeAnnotator_CpatTraining_Tests
{
    private static readonly string[] CodingTranscripts =
    {
        "CTGATCGGGTGAGGAGGGCCGATAGGCATTGTCCGAGTCTACGGTATATGATCCTCAAGACACATTAAGTTATACGCGGAACAAACGACAACTAATCCCC",
        "AGTATTGGATGCGTGCGTGGGTTGACTGCAACCGGAACGCGTGACTTTGAACGATTGATCGGGCTTTGCCCTTTTTGCTCGATCCTATTATGTG",
        "CGCGCAATTCAATTTCAACGTGGAAAATGTTTGTCACCGGTTTGGAGACGCTCCGGAGTGACAAATGGGCTGCGTTAGCTATTATGAGAAACGTTTCCCAGAGTGCGTCCTGTAGGTGATAACAGATTTCCGGCCCGTT",
        "CGATGGTAGGAGGTCACCATAGATGAGGACATGCCGAAATACGGCGAGCTCGACCTTCGGGACGCGAGCCGAGCCACCTATGTACCAGGAGAGTAGGAGGACACTGGTATGAT",
        "AGGATGGACTTGAGTCTCACCCCGGCCCGCAAGAAGGGGTTCCTAACTAGATGGACCTTTATAACGGTACCCCCGTCCTACGTCAAATAGGAGAGTATGCA",
        "TCTGACATGAACCGCCAGCACTATGAACTTTCGGTACCTCACAATTAGCTGTTCTGGCTCCGCAATTTTTGTGCAACGGAACTTTCTCCTCAACTCCACAGCGGCTGCCTCCTAAG",
        "ATCCTATGCGGTTGTTTAATGATGCAACACCTGCAAGAGGAATTCTTTAAGGAGGTGTCTGTTGCGGCGCCGCACTACTCTTAACGCAACCGGTCT",
        "GACACTACGTATTAGTTATGTGTGATCAGCCCCTATTGGGGGTATGGTCCCATAGGTGCGACGTGGCGATGAACAAGCGTTGCCACCTGACGATGAATACTAAGACATCAGAGGTAAAATAGCCATGTAGGCAA",
        "CGCCGATTATGGAATCCCCCCCGCCCGGCTCGAGCTTGGAGATGGCCTATCTCAACGCTGCTGACGGTTGATATA",
        "CGCTCCTAGAGCTAGCTTGACCCAATGATGGCGACTGACATAGTTATTGCCGTCTGGTCTTCGGTTGAACCCGGTTAGTCTTCTTTA",
        "GAAGACTGTTAATGAGATGCACCCACTCTGCAGCGGCGCACGTGATTGCATCTGCCCGGTTTCGAGGTCGTTTCAGAACACGTTGAAATGGTAAGGGATAA",
        "CTGCCGCATTTTGAGTGGTTGAATGTTGGAAAAAATCGACGCGGACTTTGAATCCTTCACCTTGTCGTTGATCACAGGGAACCAATTCATAAAGTACTTCCAACATAAGCGTCACTAGTTAGCGGCTAGTGGTAGGC",
        "GGTAATCGTTTATGCGAGAAGCCGGCCATGCGGACATCGCCCAATGGCTCGTACTTCTCATCAGATGAC",
        "AGAAGCCATAGTACTGCATGTCCTCAACTAGCATTACGCCCTACCTTCGTACCCGATTTGGGCGTCTCTACGTGAGCGGCCCGCTTATGCCGCCCCCTTAAG",
    };
    private static readonly string[] NoncodingTranscripts =
    {
        "CTATACATCGATGATTGTCGCTACCCTGGTAATCGGACACCCGTGGTTATTTCGGGAGTATGTGCGTTAGCTCAGAGATG",
        "GAGCATACATATGTGCAAGCTACTAGGTGTACCCTGTAGTAGATATACTAGTCTGTAGATATACAATC",
        "CTAGAGTGAAGTAGGTCTCTATCCAAGGACAAACTTTCCTTTGAGTTAGCATGGATTTGTTGCTTTCTGTATTCTAACACTGGC",
        "TGTTCAGTCGTTACAATCGCCACTGCACGCGTAGATACTCATGTCCTGTATAGAATCTCCAGTGCTTGAGCTCGGCGGTCCTATGGCTTGGACATTAGGACGGCTCGCTGTGCTGGCAGT",
        "CAGTAGACACCTTACTACCACCCTCATCTGTCATCGGTCCTGCGCTATGCTCAAGTGCTACGTGGTCTGCCTCATGGGTACACGTGGCGCTTTAGTGTTGTGTAGGGCAGTCAAAGTTCGATCAGACACGTATAAGCC",
        "AGTCTGGAATCCTCAACTCGTACGGCGGCTTGAGAATATATTGTATTCTAGGACCTGGGGGGCCAGGCAGTCTCTCACCTTTGCCCGTACTGTGATCAAGGATTATAGTTCGCTTAACTATTCGTGTTGGTAGTTGCTGACC",
        "TTAGGAGGATAACCTTTAGGCACCTTTAAGGCGTCTAAGAGGTGGGGCGGAAGTAGGGTGAGGAATTACCCCGCAAAGTACTTGATTAAAGGATTAGCCA",
        "TACAATGCGTTCTAGCGCTATACAAATCCAACAAGATTGTGGATTGCATGGATCCCTATGGGGAAACAGGAAGGGGTTCTCCGATTTCCACTAACCGTAATGCCACTA",
        "GGAAATATCGATCCGAATCCTAAACCGCTAGTAACACGAACAACAACCTGAAGGTTTATGAATAACTCGACCAGCATCAAGGACGGACGGATTGTCTAAGTA",
        "TCTTCCAAAACCACTGGATACCTCTTAACGGTATATGCGCCGCCTTACTACATTCGAAGTGACTTTATAAACCGGGCCGGCGGCGCCTTGCTCTGAG",
        "GTACGGGATTACTACTTAGCGTAAACAACCAGTGCGCACAGACTTTGGGTATCCTGTTGAAGGCAAAGAACGTAATTACGCCCCAAAAGGCTGCCAT",
        "TTTCCACACAATCCTGTCCATATCTGGCAAAGAGCAATTTAACCGGGTCGCTCCCAAGAACTGTGCCAGGCATGCTTCGGTGCATATATTGAGAGTGCGAGGGATGGGCTGCACTGCTCGATTGCAATTCATGTCTGTAGT",
        "TGCGTGGAAACTCAATACTTCGCCTAGTTGTTGATAACGCCCAGGATACTTTTGTTTTTCCGTGCAATTACGCATGCCCGTGCCTGCTG",
        "ACATAAGTGTGGCGAATACGGATGACTGAGACCTGTCTATGCTTAAGTAGAAGTTGGACAAACAAGGAAGATTTATGATAATTGACGATCGCATTTGCCAAGCGTATCCCGCAGACCTTGCCGCGGACCAAGCGTG",
    };
    private static readonly string[] CodingCds =
    {
        "ATGTTCCGAAAGTCAACCGATTGGATCACCCGCATGGACAATGAGTTCGCGGATAATGATAAGATCGAGAGACGCGACGAGGAGGGCACGGAGAGTACCGGTGAGGCCTAG",
        "ATGCGGCAGGATCGTCAGTACGGGATGGAGGAATCTGTACACTTGACGCGCGTGAGCCTCAAGGGCGTGGAGTCTCAGGGGGCAGATCGTTCCTTACTCGGAAATGGGGCTTTTTCACACCCTGGCAACTAG",
        "ATGAGGACGAATGTGACTATGCCTTCTACCCGTGACTGCTCATCGCACCGCCCCTGCTTCAGCAGGATGATGCTCGCACAAAAACGCCCGGAGTTCGGGCTTGCCGGCAATTGA",
        "ATGCGACGAACGGCGAATACAATCAGCCCCAAGAAAAGCACGCCGAACTCTGCCAGCATTCCCCCGACTTCGTTTCTACTCACCACGCCCCAATCCAGCAACTTCGAAAGCGACAACGAATTGTCGGGGTAG",
        "ATGCCGAGGGTTGAAGTTCCAATCTATTTGTCGTATTCAGTCATCCTTCCCAACAACCACATTGGGTTCCGTTTATTGTTCCCGATCATCATAATAGCGTGCTCCCTCGGATACACCGGATGGGTTTAA",
        "ATGTGCCGGGCGACCTTAGTACTAGTCACCAGCTTTAACGCACTAGGAGCCGCACTCTTAAGCGACAACTACGTAACGGTCAGATTGAGTTCTAACATCCTGGTGCATGACTCCAATGTCTGCTGA",
        "ATGTATCTTTCATCTAGAAGCCTTTCTGTTGCCGTTAGCAGAGCCACTAATGGCGAACTCGTTTGCAATAGTAAGACCCGTAACTGCAATCATTTGGCCATACACCTGTATCTGATGGAATCTTGTTGGTGCTAA",
        "ATGAGAACTCGGACAATGACAGGCTACTGCAGCGGTGGAAGCTGCTTGATCTTCCGACTCTTTTCCAGCCGGCTTCCCAGTCATGTAAACGACCGGAACCCAGCGCTTAATTTGTACTTGTAA",
        "ATGCAAAAGGTTCAAACTCGCTGGTCAATTAGCTGGGAGGGGACCAATGGGCGGACATTAGCGTACACCTCGACATCGCGTCCCTGA",
        "ATGTGTCAAACTGAAATGACGAGCAAGATACAGCAGGTGCACCTGAGGCCCAACATTAACGTCGACAGTCTCTGTACGTCACGGCTGTCCCGACTGGGGTTGCCCTCCGGTGCGCGCTAG",
        "ATGTCCAACGCGTCTCGGAATGGCCCTCTGTTGGCAGGGAAGAGTCAAAAACCGCCATTGCATAAGTGCATCAGGCCAACCTAG",
        "ATGGGAATGGTTACCTCAAGCGATCGACAGGCGACCTTCTCGAGCTGCTACGCAGTGGTACTGCACGCCAGCTCAGGTCTTCGGGGGCTGTGGGTGGAGTATAACGGGTACGGTTTCGCTACCAAGTGA",
        "ATGTTCCCGCGTAAGAACGAAATGCAACTAACATGTTGCACCGAAGTCCTCTATTCCATGGTGATGGCTCTATCCACATGGCACAAGGCGATAGCCCCAGCTGGGGCGAGATAG",
        "ATGCTCTTGACGATCCTGGGTCACCTAGTGTTTGACCAGACCCTCCACACACTGTTTTTTCCTTGGTGTACAACCGTCCAATAG",
    };
    private static readonly string[] NoncodingForTable =
    {
        "GGCACCTCATGCGCCACCCCTAAGCACTATTTCAGGATCCTTCATCTGAATCGATCTGACCCGGCACGAGAGTGCCTCGAAGCCTAATCGACTCGCGGTTATCGCCTT",
        "AGTCGACTTTTAAAGGTGGTCGATATGCACGGAAGGCAGACCTGTCAATGCACGCCGTAACCGATGGATCTAATCTCAACGTACATAAACTTGTTCATGAGCCCAACC",
        "ATACTTCACCTACTACTTAAATCGTCCTGCGCCATTGTAACTTCCCCTATGTTGCCCTAATAAAGCTTTCAACGCACGGGTGATCGCGAAAATGCCTGATGTACTTTACAATGTCTCCTTAGGACAACGGAACGACATGAACGGAAGCCCTCGCCT",
        "CAGAATCTCTCAAGGTAGCGACGTACTTAATTCATCTGTATACTTAACGCAAGTAGTATAGGGCCGCGCCAGCATTCACAACACGAGAGTGAGCTATTGTATCTTATTACTATCTAGTACGCAAATGTATAGGCTGGTATGCTCTGCTCAG",
        "CGCTTTTGGACGCCACGAAGCGGCCACTAATGCAAGCTTATTGGGGCAAAAAATAGGACGGGCGGAGTGCGCCTAGATCGATTTCAGGGTTCTTAGATTAATA",
        "CATGGGTCTGTCGGCAGATAAATCAAAGGGGAGCTTTATTTAAAATTGTGCAATACCACCGCCTGTGACCTCGGATATAAATTACCTACCTCTCAGCAGATGTGCCCTATAACTGCTACGTTGC",
        "GTCGCCGCACGGGCTTGGTGGCATCTTAGAAACCCATCTCATCATCACGAATGGGTGGTGGAAACCGTACGATCATGGTTCATTACAGTACTCACATCCGCTAGTGGGCCAATGGGACCT",
        "TTCATATAGGTAGTTCTTGTCTAGGTGCTTCGCCGATACCAGCTGCGGAACGGTACTACGGGTATGGCCAATGAACTCGCTGAAGGACGCCCCTGACGAAGTTACTAACAAGCT",
        "TACCGTCTCGCGTGTCCCGGTAGTGACAATGGTTGGGGGGAGCTAAATTGCTTGCCCCCGTAGAATGTCACGTTCCAGAAAAGATATAAGAGAATGTGAGGCTCTGAAGTACGATCACGCGATCGGCCTTTGCAT",
        "TGTGAATGTACCTACGCCATTTCATTATTGTCGGCCTGTCACGTGACGGCACGCTCGAAATCGCTTATCTCCACACGAATGGATCACGGATCTCAACGCTTTGTGACAATGATCAGATTAGATCGGTCACAAGGCGCCCCCCTTT",
        "GAATATTTACATTGCCTCGAAAGGATCCACTCACCTATATGACCCACTTTACTAGCTATCAATGACTAGAATGAAAGCCAGGTGAAGCGATCTTCATCTCATCACGAGCGAGTTCGCTATTCGGGATAA",
        "CGTATGCTGCACCAACCGTAATGAAGCCTAGGGGATACCTGTGTCCGCAGCTGAAACAGACCTCATTATGTGGTAAAAAAACAATTTAGTCGTTGCGGCCGACGCT",
        "GCCAAGCGACTAGCCTCTAATCACCGAATGACCGGTTGGGTCTGTATATGCTGCTTGCCCGCCGCCTTCCACAGGATCTTTACTAAGCTTGTGGGCCGCTGGACCGGTTCCGTCGTACT",
        "CAGCCTGTGGTATCAGAGGGGAAGGTCAAAATGTTGAAGAGGAATCCCTGGGAACAGTTAAGAACTAAAAAGAATCCCTTAGGGCCGATAACATCCCTCGTCTGCCCCCATAAAGTGCTACACGGCCATGTTGCTGCGATGATGTAGAAGGT",
    };
    private static readonly string[] HeldOut =
    {
        "TACGTATGATCCGCCTCCAGAATCCTGCTCTCGAAATAATTCTTACTCATTGCTGGTGCCTCTGTGAGACGCTGTTGCCACCCACCCCGACGTTGCTGTAGTCGCAAACGG",
        "TGCTGGTGCGCTATAAGAGATGGTTGGAACTCGACTGTAAATATCTGCTCGATTAGCACGTTTGGAGGCTCAAAGTAGGGCGAACATATAACCCTAGTCGTGAAACAGCCGTAACCTGCA",
        "AATGAGAGCAATGGTCGGTCTTAGGGGACGGGGAACCTGGAAAATCCGATAACGCCAGCACTAAAACGCGCATGCAGTAACGAATGTTATAA",
    };

    [Test]
    public void BuildHexamerTable_HandCheckable_CodingInFrameNoncodingAllFrames()
    {
        // coding ATGAAATAA: in-frame hexamers ATGAAA, AAATAA (1/2 each); noncoding AAAAAAA: AAAAAA at 0 and 1 (1.0).
        var (coding, noncoding) = CpatModel.BuildHexamerTable(new[] { "atgaaataa" }, new[] { "AAAAAAA" });
        Assert.That(coding, Has.Count.EqualTo(4096));
        Assert.That(noncoding, Has.Count.EqualTo(4096));
        Assert.That(coding["ATGAAA"], Is.EqualTo(0.5));
        Assert.That(coding["AAATAA"], Is.EqualTo(0.5));
        Assert.That(coding.Values.Sum(), Is.EqualTo(1.0));
        Assert.That(noncoding["AAAAAA"], Is.EqualTo(1.0));
        Assert.That(noncoding["ATGAAA"], Is.EqualTo(0.0));
    }

    [Test]
    public void BuildHexamerTable_HexamersWithN_AreNotCounted()
    {
        var (coding, _) = CpatModel.BuildHexamerTable(new[] { "ATGNAAGGGAAA" }, new[] { "ACGTAC" });
        // i=0 ATGNAA (N, skipped); i=3 NAAGGG (N, skipped); i=6 GGGAAA counted -> total 1.
        Assert.That(coding["GGGAAA"], Is.EqualTo(1.0));
        Assert.That(coding.Values.Sum(), Is.EqualTo(1.0));
    }

    [Test]
    public void BuildHexamerTable_MatchesMakeHexamerTabOracle()
    {
        var (coding, noncoding) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        Assert.That(coding["AGGACG"], Is.EqualTo(0.0019011406844106464));
        Assert.That(noncoding["AGGACG"], Is.EqualTo(0.001176470588235294));
        Assert.That(coding["GGGAAG"], Is.EqualTo(0.0019011406844106464));
        Assert.That(noncoding["GGGAAG"], Is.EqualTo(0.000588235294117647));
        Assert.That(coding["AAAAAA"], Is.EqualTo(0.0));
        Assert.That(noncoding["AAAAAA"], Is.EqualTo(0.0017647058823529412));
        Assert.That(coding.Values.Count(v => v > 0), Is.EqualTo(494));
        Assert.That(noncoding.Values.Count(v => v > 0), Is.EqualTo(1347));
    }

    [Test]
    public void ToHexamerTable_FormatAndRoundTrip()
    {
        var (coding, noncoding) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        var model = new CpatModel(coding, noncoding, 0, 0, 0, 0, 0, 0.5);
        var lines = model.ToHexamerTable();

        Assert.That(lines, Has.Count.EqualTo(4097));
        Assert.That(lines[0], Is.EqualTo("hexamer	coding	noncoding"));
        Assert.That(lines[1], Is.EqualTo("AAAAAA	0.0	0.0017647058823529412")); // Python str(float)
        Assert.That(lines[4096], Is.EqualTo("TTTTTT	" + PyFloat(coding["TTTTTT"]) + "	" + PyFloat(noncoding["TTTTTT"])));

        var (c2, n2) = CpatModel.ParseHexamerTable(lines);
        Assert.That(c2.Count, Is.EqualTo(4096));
        foreach (var kv in coding) Assert.That(c2[kv.Key], Is.EqualTo(kv.Value));
        foreach (var kv in noncoding) Assert.That(n2[kv.Key], Is.EqualTo(kv.Value));
    }

    [Test]
    public void ToHexamerTable_SmallValue_UsesPythonExponentForm()
    {
        var model = new CpatModel(new Dictionary<string, double> { ["AAAAAA"] = 1e-5, ["CCCCCC"] = 0.25 },
            new Dictionary<string, double> { ["AAAAAA"] = 1.0, ["CCCCCC"] = 3e-10, ["GGGGGG"] = 0.5 }, 0, 0, 0, 0, 0, 0.5);
        Assert.That(model.ToHexamerTable(), Is.EqualTo(new[]
        {
            "hexamer	coding	noncoding", "AAAAAA	1e-05	1.0", "CCCCCC	0.25	3e-10",
        }));
    }

    private static string PyFloat(double v)
    {
        string s = v.ToString("R", System.Globalization.CultureInfo.InvariantCulture);
        if (s.Contains('E')) return s.Replace('E', 'e');
        return s.Contains('.') ? s : s + ".0";
    }

    [Test]
    public void Train_MatchesRGlmAndCpatOracle()
    {
        var (coding, noncoding) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        var result = CpatModel.Train(CodingTranscripts, NoncodingTranscripts, coding, noncoding);
        var m = result.Model;

        // Rscript: coef(mylogit) of the make_logitModel fit, deviance, iter.
        Assert.That(m.Intercept, Is.EqualTo(-4.5133790317038329).Within(1e-9).Percent);
        Assert.That(m.MrnaLengthCoefficient, Is.EqualTo(-0.055765617893335499).Within(1e-9).Percent);
        Assert.That(m.OrfLengthCoefficient, Is.EqualTo(0.10892994611389771).Within(1e-9).Percent);
        Assert.That(m.FickettCoefficient, Is.EqualTo(5.0407879119821963).Within(1e-9).Percent);
        Assert.That(m.HexamerCoefficient, Is.EqualTo(3.8232944199889674).Within(1e-9).Percent);
        Assert.That(result.Deviance, Is.EqualTo(17.373176228178064).Within(1e-9).Percent);
        Assert.That(result.Iterations, Is.EqualTo(6));
        Assert.That(result.Converged, Is.True);
        Assert.That(result.CodingSamples, Is.EqualTo(13));      // 14 minus one without an ORF > 30 nt
        Assert.That(result.NoncodingSamples, Is.EqualTo(9));
        Assert.That(m.CodingHexamerFrequencies, Has.Count.EqualTo(4096));
    }

    [Test]
    public void Train_TgRocCutoff_MatchesIndependentOracle()
    {
        var (coding, noncoding) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        var result = CpatModel.Train(CodingTranscripts, NoncodingTranscripts, coding, noncoding);
        // R: argmin |Se-Sp| over the fitted probabilities (Se = P(p>=c|coding), Sp = P(p<c|noncoding)).
        Assert.That(result.Model.CodingProbabilityCutoff, Is.EqualTo(0.64722981501159438).Within(1e-12));
        Assert.That(result.Sensitivity, Is.EqualTo(10.0 / 13.0).Within(1e-15));
        Assert.That(result.Specificity, Is.EqualTo(7.0 / 9.0).Within(1e-15));
    }

    [Test]
    public void Train_HeldOutProbabilities_MatchCpat()
    {
        var (coding, noncoding) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        var model = CpatModel.Train(CodingTranscripts, NoncodingTranscripts, coding, noncoding).Model;
        var options = new CpatOptions { MinOrfLength = 30 };

        // cpat -g h.fa -d sm.logit.RData -x hex.tsv --min-orf 30 --top-orf 5 -> ORF_prob.best.tsv (15 significant digits).
        var s0 = PredictBestCodingOrf(HeldOut[0], model, options)!;
        Assert.That(s0.OrfLength, Is.EqualTo(96));
        Assert.That(s0.FickettScore, Is.EqualTo(1.1087));
        Assert.That(s0.HexamerScore, Is.EqualTo(-0.571428571428571).Within(1e-13));
        Assert.That(s0.CodingProbability, Is.EqualTo(0.959218646900687).Within(1e-9));

        var s2 = PredictBestCodingOrf(HeldOut[2], model, options)!;
        Assert.That(s2.OrfLength, Is.EqualTo(51));
        Assert.That(s2.CodingProbability, Is.EqualTo(0.192130655232657).Within(1e-9));

        // cpat lists S1 in no_ORF.txt.
        Assert.That(PredictBestCodingOrf(HeldOut[1], model, options), Is.Null);
    }

    [Test]
    public void Train_TrainedModelTableRoundTripsThroughParse()
    {
        var model = CpatModel.Train(CodingTranscripts, NoncodingTranscripts, CodingCds, null);
        var (c, n) = CpatModel.ParseHexamerTable(model.Model.ToHexamerTable());
        Assert.That(c, Is.EqualTo(model.Model.CodingHexamerFrequencies));
        Assert.That(n, Is.EqualTo(model.Model.NoncodingHexamerFrequencies));
    }

    [Test]
    public void Train_StopCodonOption_ChangesTheLongestOrf()
    {
        // The default stops are TAG, TAA, TGA (CPAT docs). A TAG-only stop set (what the CPAT CLI's
        // unstripped default "TAG, TAA, TGA" effectively yields) produces longer ORFs and a different fit.
        var (c, n) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        var all = CpatModel.Train(CodingTranscripts, NoncodingTranscripts, c, n);
        var tagOnly = CpatModel.Train(CodingTranscripts, NoncodingTranscripts, c, n,
            new CpatTrainingOptions { StopCodons = new[] { "TAG" } });
        Assert.That(tagOnly.Model.OrfLengthCoefficient, Is.Not.EqualTo(all.Model.OrfLengthCoefficient));
    }

    [Test]
    public void Train_ArgumentContract()
    {
        Assert.Throws<ArgumentNullException>(() => CpatModel.Train(null!, NoncodingTranscripts, CodingCds));
        Assert.Throws<ArgumentNullException>(() => CpatModel.BuildHexamerTable(null!, NoncodingForTable));
        Assert.Throws<ArgumentException>(() => CpatModel.BuildHexamerTable(new[] { "ATGXAA" }, new[] { "AAAAAA" }));
        Assert.Throws<ArgumentException>(() => CpatModel.BuildHexamerTable(new[] { "ATG" }, new[] { "AAAAAA" })); // no hexamer
        var (c, n) = CpatModel.BuildHexamerTable(CodingCds, NoncodingForTable);
        // no coding transcript has an ORF
        Assert.Throws<ArgumentException>(() => CpatModel.Train(new[] { "CCCCCCCCC" }, NoncodingTranscripts, c, n));
        Assert.Throws<ArgumentOutOfRangeException>(() => CpatModel.Train(CodingTranscripts, NoncodingTranscripts, c, n, new CpatTrainingOptions { MinOrfLength = -1 }));
    }
}
