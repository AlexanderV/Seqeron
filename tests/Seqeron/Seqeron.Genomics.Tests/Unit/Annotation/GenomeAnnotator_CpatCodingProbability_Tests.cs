// ANNOT-CODING-001 — CPAT coding probability (F13)
// Evidence: docs/Evidence/ANNOT-CODING-001-Evidence.md (§ CPAT end-to-end oracle)
// TestSpec: tests/TestSpecs/ANNOT-CODING-001.md
// Source: Wang L et al. (2013) CPAT, NAR 41(6):e74, https://doi.org/10.1093/nar/gkt006;
//         CPAT 3.0.5 (liguowang/cpat @ ea5aca3): cpmodule/__init__.py `cpat`, find_orfs.ORFFinder,
//         fickett.fickett_value, FrameKmer.kmer_ratio, utils.coding_prediction (R predict on the glm).
// Every expected row below is the ORF_prob.tsv / ORF_prob.best.tsv output of the real `cpat` command
// (CPAT 3.0.5 + R 4.x) run on the five transcripts with the synthetic hexamer table built by
// SyntheticTable() and CPAT's prebuilt Human_logitModel.RData (coefficients below = coef(mylogit)).
// R writes 15 significant digits, hence the 1e-13 relative tolerance.

namespace Seqeron.Genomics.Tests.Unit.Annotation;

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.GenomeAnnotator;

[TestFixture]
public class GenomeAnnotator_CpatCodingProbability_Tests
{
    // CPAT prebuilt Human model: coef(mylogit) of Label ~ mRNA + ORF + Fickett + Hexamer; cutoff 0.364
    // (Human_cutoff.txt, TG-ROC).
    private const double Intercept = -7.8284881840360567;
    private const double MrnaCoef = -5.4831277888113823e-05;
    private const double OrfCoef = 1.1187815405519514e-02;
    private const double FickettCoef = 3.4471695641643225;
    private const double HexamerCoef = 6.0682759241106110;
    private const double HumanCutoff = 0.364;

    private const string T1 = "TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTATGGTATGCCGTGAATTTTCTCATTTCCGCACATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTCGTATTAGAGGCGTTAGTGTGTAACAGGTAAGATACGCCAACGGAATTAACCTCTATGGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTACTGACGTCTGTGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGTAGAAGCGTCGAATCCTGATACA";
    private const string T2 = "CCATGTCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGGCGTGTATAAACAGGAAATAAACTGACGAATGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGC";
    private const string T3 = "TTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGG";
    private const string T4 = "GCGAAGGGTGCGTACGGTGTTTATATCATAGCTTGTAAGAGAGAACCACAAGCAGAAGAACGATGGTCGTTGACGATCTTTGAAGGATTGCATAAGGTCGGAACGATAGGTTTAAAAGGAATTGAATTGCCACCGTCCTGAGATTCTTCAGAATTACCCTCTTGGGGCCGTCACGCCTCATATACACCCCAATGGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGATGAGATCG";
    private const string T5 = "ATGCAGCTCGCGGAGCTAGCGACTGCCCTTGTGTCGGTGATCCCAAGGGTCCTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGGGTGAACTGAAAGGAGCATATTTGAGATCTTGAATGAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTAGGATGCAGAACTCCCGGGAACGCTTTTCAACTAGGAACCCCGTCATGGAAGTAGATCGCGTCGAAGGGAGTTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAG";

    private static readonly Dictionary<string, string> Transcripts = new()
    {
        ["t1"] = T1, ["t2"] = T2, ["t3"] = T3, ["t4"] = T4, ["t5"] = T5,
    };

    private sealed record Row(string Id, int Sn, char Strand, int Frame, int Start, int End, int Orf,
        double Fickett, double Hexamer, double Prob);

    // cpat -x synth_hex.tsv -d Human_logitModel.RData -g tx.fa  (defaults: sense only, top 5, min 75)
    private static readonly Row[] SenseRows =
    {
        new("t1", 1, '+', 2, 41, 226, 186, 1.0961, 0.0651496247708103, 0.168675901817624),
        new("t1", 2, '+', 3, 258, 368, 111, 0.8586, 0.00708755665998371, 0.0264635700760982),
        new("t2", 1, '+', 2, 107, 244, 138, 0.6273, 0.422637335920057, 0.172093961549197),
        new("t2", 2, '+', 3, 3, 98, 96, 0.8541, 0.331122475098671, 0.140123250478773),
        new("t4", 1, '+', 3, 192, 281, 90, 1.1267, 0.379887062364257, 0.34342209160935),
        new("t5", 1, '+', 1, 1, 126, 126, 0.7124, 0.453476340596222, 0.225827839785713),
        new("t5", 2, '+', 1, 127, 252, 126, 1.0753, 0.0977936238061278, 0.105323618942506),
        new("t5", 3, '+', 2, 254, 379, 126, 0.7782, 0.338781673973669, 0.15430906424874),
    };

    // cpat ... --antisense
    private static readonly Row[] AntisenseRows =
    {
        new("t1", 1, '+', 2, 41, 226, 186, 1.0961, 0.0651496247708103, 0.168675901817624),
        new("t1", 2, '+', 3, 258, 368, 111, 0.8586, 0.00708755665998371, 0.0264635700760982),
        new("t1", 3, '-', 2, 147, 61, 87, 0.9645, 0.288768966290966, 0.141936544139421),
        new("t2", 1, '+', 2, 107, 244, 138, 0.6273, 0.422637335920057, 0.172093961549197),
        new("t2", 2, '+', 3, 3, 98, 96, 0.8541, 0.331122475098671, 0.140123250478773),
        new("t2", 3, '-', 2, 165, 73, 93, 0.8701, 0.0609437912434101, 0.03130432592611),
        new("t4", 1, '-', 1, 268, 26, 243, 1.1479, 0.359893085478537, 0.734076720470705),
        new("t4", 2, '-', 2, 93, 1, 93, 0.8028, 0.258640487105822, 0.0782176700141873),
        new("t4", 3, '+', 3, 192, 281, 90, 1.1267, 0.379887062364257, 0.34342209160935),
        new("t5", 1, '-', 3, 245, 66, 180, 0.9403, 0.0557120404287372, 0.0948320657718226),
        new("t5", 2, '-', 2, 297, 166, 132, 0.7244, 0.36419154624484, 0.159047379295967),
        new("t5", 3, '+', 1, 1, 126, 126, 0.7124, 0.453476340596222, 0.225827839785713),
        new("t5", 4, '+', 1, 127, 252, 126, 1.0753, 0.0977936238061278, 0.105323618942506),
        new("t5", 5, '+', 2, 254, 379, 126, 0.7782, 0.338781673973669, 0.15430906424874),
    };

    /// <summary>
    /// Synthetic hexamer table (identical file fed to CPAT): hexamers in lexicographic A&lt;C&lt;G&lt;T order,
    /// index i → coding (i mod 7)/1000, noncoding (i mod 5)/1000 — exercises the ln, +1, −1 and
    /// both-zero branches of kmer_ratio.
    /// </summary>
    private static CpatModel SyntheticModel()
    {
        var lines = new List<string> { "hexamer\tcoding\tnoncoding" };
        const string bases = "ACGT";
        for (int i = 0; i < 4096; i++)
        {
            var chars = new char[6];
            for (int p = 0, v = i; p < 6; p++, v /= 4) chars[5 - p] = bases[v % 4];
            lines.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{new string(chars)}\t{(i % 7) / 1000.0:R}\t{(i % 5) / 1000.0:R}"));
        }

        var (coding, noncoding) = CpatModel.ParseHexamerTable(lines);
        return new CpatModel(coding, noncoding, Intercept, MrnaCoef, OrfCoef, FickettCoef, HexamerCoef, HumanCutoff);
    }

    private static void AssertRows(Row[] expected, bool antisense)
    {
        var model = SyntheticModel();
        var options = new CpatOptions { SearchAntisense = antisense };
        foreach (var (id, seq) in Transcripts)
        {
            var actual = PredictCodingProbability(seq, model, options);
            var rows = expected.Where(r => r.Id == id).ToArray();
            Assert.That(actual, Has.Count.EqualTo(rows.Length), $"{id}: candidate count");
            for (int k = 0; k < rows.Length; k++)
            {
                var (e, a) = (rows[k], actual[k]);
                Assert.Multiple(() =>
                {
                    Assert.That((a.Strand, a.Frame, a.Start, a.End, a.OrfLength, a.MrnaLength),
                        Is.EqualTo((e.Strand, e.Frame, e.Start, e.End, e.Orf, seq.Length)), $"{id}_ORF_{e.Sn}");
                    Assert.That(a.OrfSequence, Has.Length.EqualTo(e.Orf));
                    Assert.That(a.FickettScore, Is.EqualTo(e.Fickett).Within(1e-11).Percent, $"{id}_ORF_{e.Sn} Fickett");
                    Assert.That(a.HexamerScore, Is.EqualTo(e.Hexamer).Within(1e-11).Percent, $"{id}_ORF_{e.Sn} Hexamer");
                    Assert.That(a.CodingProbability, Is.EqualTo(e.Prob).Within(1e-11).Percent, $"{id}_ORF_{e.Sn} prob");
                    Assert.That(a.IsCoding, Is.EqualTo(e.Prob >= HumanCutoff));
                });
            }
        }
    }

    [Test]
    public void PredictCodingProbability_SenseStrand_MatchesCpat()
        => AssertRows(SenseRows, antisense: false);

    [Test]
    public void PredictCodingProbability_Antisense_MatchesCpat()
        => AssertRows(AntisenseRows, antisense: true);

    [Test]
    public void PredictCodingProbability_NoQualifyingOrf_ReturnsEmpty()
    {
        // CPAT lists t3 in no_ORF.txt.
        Assert.That(PredictCodingProbability(T3, SyntheticModel()), Is.Empty);
        Assert.That(PredictBestCodingOrf(T3, SyntheticModel()), Is.Null);
    }

    [Test]
    public void PredictCodingProbability_UnterminatedOrf_RunsToLastCompleteCodon()
    {
        // t2 ORF_1: ATG at 107 with no in-frame stop → ends at the transcript end (244), as CPAT.
        var orf = PredictCodingProbability(T2, SyntheticModel())[0];
        Assert.That((orf.Start, orf.End, orf.OrfSequence), Is.EqualTo((107, 244, T2.Substring(106))));
    }

    [Test]
    public void PredictBestCodingOrf_ByProbabilityAndByLength_MatchCpatBestFiles()
    {
        var model = SyntheticModel();
        var anti = new CpatOptions { SearchAntisense = true };
        Assert.Multiple(() =>
        {
            // ORF_prob.best.tsv, --best-orf p (default), --antisense
            Assert.That(Best(T1, model, anti, false), Is.EqualTo(('+', 41, 226)));
            Assert.That(Best(T2, model, anti, false), Is.EqualTo(('+', 107, 244)));
            Assert.That(Best(T4, model, anti, false), Is.EqualTo(('-', 268, 26)));
            Assert.That(Best(T5, model, anti, false), Is.EqualTo(('+', 1, 126)));
            // --best-orf l
            Assert.That(Best(T4, model, anti, true), Is.EqualTo(('-', 268, 26)));
            Assert.That(Best(T5, model, anti, true), Is.EqualTo(('-', 245, 66)));
        });
    }

    private static (char, int, int) Best(string seq, CpatModel model, CpatOptions options, bool byLength)
    {
        var b = PredictBestCodingOrf(seq, model, options, byLength)!;
        return (b.Strand, b.Start, b.End);
    }

    [Test]
    public void PredictCodingProbability_CustomCodonsAndTopCount_MatchCpat()
    {
        // cpat ... --start "ATG,GTG,TTG" --stop "TAA,TAG" --top-orf 2
        var options = new CpatOptions { StartCodons = new[] { "ATG", "GTG", "TTG" }, StopCodons = new[] { "TAA", "TAG" }, TopOrfCount = 2 };
        var model = SyntheticModel();
        var t2 = PredictCodingProbability(T2, model, options);
        var t5 = PredictCodingProbability(T5, model, options);
        Assert.Multiple(() =>
        {
            Assert.That(t2.Select(p => (p.Frame, p.Start, p.End, p.OrfLength)), Is.EqualTo(new[] { (2, 20, 244, 225), (3, 3, 98, 96) }));
            Assert.That(t2[0].CodingProbability, Is.EqualTo(0.505327372224478).Within(1e-11).Percent);
            Assert.That(t2[0].IsCoding, Is.True);
            Assert.That(t5.Select(p => (p.Frame, p.Start, p.End, p.OrfLength)), Is.EqualTo(new[] { (1, 1, 252, 252), (2, 188, 379, 192) }));
            Assert.That(t5[0].HexamerScore, Is.EqualTo(0.259886402174).Within(1e-10).Percent);
        });
    }

    [Test]
    public void PredictCodingProbability_EmptyStartList_AnyNonStopCodonOpens()
    {
        // ORFFinder(seq, 75).orf_candidates(start_coden=[], n_candidate=3) on t1 →
        // (+,2,16,226,210), (+,3,2,140,138), (+,3,233,368,135) (0-based starts).
        var options = new CpatOptions { StartCodons = Array.Empty<string>(), TopOrfCount = 3 };
        var r = PredictCodingProbability(T1, SyntheticModel(), options);
        Assert.That(r.Select(p => (p.Frame, p.Start, p.End, p.OrfLength)),
            Is.EqualTo(new[] { (2, 17, 226, 210), (3, 3, 140, 138), (3, 234, 368, 135) }));
    }

    [Test]
    public void PredictCodingProbability_Features_DelegateToCanonicalScores()
    {
        var model = SyntheticModel();
        foreach (var p in PredictCodingProbability(T5, model, new CpatOptions { SearchAntisense = true }))
        {
            Assert.That(p.FickettScore, Is.EqualTo(CalculateFickettScore(p.OrfSequence, FickettTableSet.Cpat)));
            Assert.That(p.HexamerScore, Is.EqualTo(CalculateCodingPotential(p.OrfSequence, model.CodingHexamerFrequencies, model.NoncodingHexamerFrequencies)));
            Assert.That(p.CodingProbability, Is.EqualTo(model.PredictProbability(p.MrnaLength, p.OrfLength, p.FickettScore, p.HexamerScore)));
        }
    }

    [Test]
    public void PredictCodingProbability_IsCaseInsensitive()
    {
        var model = SyntheticModel();
        var upper = PredictCodingProbability(T4, model, new CpatOptions { SearchAntisense = true });
        var lower = PredictCodingProbability(T4.ToLowerInvariant(), model, new CpatOptions { SearchAntisense = true });
        Assert.That(lower, Is.EqualTo(upper));
    }

    [Test]
    public void CpatModel_PredictProbability_IsLogisticResponse()
    {
        var model = SyntheticModel();
        double eta = Intercept + MrnaCoef * 1000 + OrfCoef * 300 + FickettCoef * 0.9 + HexamerCoef * 0.2;
        Assert.That(model.PredictProbability(1000, 300, 0.9, 0.2), Is.EqualTo(1 / (1 + Math.Exp(-eta))));
    }

    [Test]
    public void CpatModel_ParseHexamerTable_ReadsCpatFormat()
    {
        var (c, n) = CpatModel.ParseHexamerTable(new[] { "hexamer\tcoding\tnoncoding", "AAAAAA\t0.000647111\t0.002072893", "", "AAAAAC 1e-3   2E-3" });
        Assert.Multiple(() =>
        {
            Assert.That(c, Has.Count.EqualTo(2));
            Assert.That(c["AAAAAA"], Is.EqualTo(0.000647111));
            Assert.That(n["AAAAAA"], Is.EqualTo(0.002072893));
            Assert.That((c["AAAAAC"], n["AAAAAC"]), Is.EqualTo((0.001, 0.002)));
        });
        Assert.Throws<FormatException>(() => CpatModel.ParseHexamerTable(new[] { "AAAAAA\t0.1" }));
        Assert.Throws<FormatException>(() => CpatModel.ParseHexamerTable(new[] { "AAAAAA\tx\t0.1" }));
        Assert.Throws<ArgumentNullException>(() => CpatModel.ParseHexamerTable(null!));
    }

    [Test]
    public void InvalidArguments_Throw()
    {
        var model = SyntheticModel();
        var t = new Dictionary<string, double>();
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => PredictCodingProbability(null!, model));
            Assert.Throws<ArgumentNullException>(() => PredictCodingProbability(T1, null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => PredictCodingProbability(T1, model, new CpatOptions { MinOrfLength = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PredictCodingProbability(T1, model, new CpatOptions { TopOrfCount = 0 }));
            Assert.Throws<ArgumentNullException>(() => new CpatModel(null!, t, 0, 0, 0, 0, 0, 0.5));
            Assert.Throws<ArgumentException>(() => new CpatModel(t, t, double.NaN, 0, 0, 0, 0, 0.5));
            Assert.Throws<ArgumentException>(() => new CpatModel(t, t, 0, 0, 0, 0, 0, 1.5));
        });
    }
}
