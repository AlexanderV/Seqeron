// ONCO-ASCAT-001 — Battenberg determine_copynumber alternative solutions A–F, SDfrac and bootstrap CIs (B24 F41)
// Evidence: docs/Evidence/ONCO-ASCAT-001-Evidence.md (§ F41)
// TestSpec: tests/TestSpecs/ONCO-ASCAT-001.md (§15)
// Source: Wedge-lab/battenberg master (57a8f7e) R/fitcopynumber.R determine_copynumber (sub-clonal branch: all.edges,
//         sdl/sdtau, "Bootstrapping to obtain 95% confidence intervals": sample(BAFke, replace = T), sd(permFraction),
//         sort(permFraction)[25] / [975]); R/orderEdges.R; R 4.3.3 RNG.c (Mersenne-Twister, set.seed, R_unif_index).
// Expected rows: set.seed(seed); determine_copynumber(...) sourced verbatim in R 4.3.3, every column printed with %.17g.
// NA (R) = NaN here; a solution with NA copy numbers has null CNs.

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_BattenbergBootstrap_Tests
{
    private static readonly Dictionary<string, (double Rho, double Psit, double LogR, double BafSeg, double[] Snps, int Seed, int Perms)> Inputs = new()
    {
        ["s1"] = (0.8, 2.5, 0.222392, 0.672857, new[] { 0.6689, 0.6759, 0.6739, 0.6709, 0.6749 }, 7, 1000),
        ["s1n500"] = (0.8, 2.5, 0.222392, 0.672857, new[] { 0.6689, 0.6759, 0.6739, 0.6709, 0.6749 }, 7, 500),
        ["s1n20"] = (0.8, 2.5, 0.222392, 0.672857, new[] { 0.6689, 0.6759, 0.6739, 0.6709, 0.6749 }, 7, 20),
        ["s2"] = (0.8, 2.5, 0.222392, 0.672857, new[] { 0.6729 }, 3, 1000),
        ["s3"] = (0.55, 3.1, 0.01946, 0.796306, new[] { 0.7938, 0.7866, 0.8057 }, -12345, 1000),
        ["s4"] = (0.35, 1.9, 0.103462, 0.698502, new[] { 0.693, 0.7098, 0.7167, 0.7037, 0.6761 }, 2147483647, 1000),
        ["s5"] = (1.0, 2.0, 0.584963, 0.696667, new[] { 0.6927, 0.6997, 0.6977, 0.6947, 0.6987, 0.71, 0.69, 0.702, 0.688, 0.6955, 0.7011, 0.6932 }, 0, 1000),
    };

    // Per solution A–F: nMaj1, nMin1, frac1, nMaj2, nMin2, frac2, SDfrac, SDfrac_BS, frac1_0.025, frac1_0.975.
    private static readonly Dictionary<string, (double PVal, double Ntot, double[][] Solutions)> Expected = new()
    {
        ["s1"] = (PVal: 2.1020858614373501e-05, Ntot: 2.9999989778315501, Solutions: new[]
        {
            new double[] { 2.0, 0.0, 0.15605024544591201, 2.0, 1.0, 0.84394975455408805, 0.0064798217571942303, 0.0058290159705588597, 0.145295959445356, 0.16716042067841799 },
            new double[] { 2.0, 0.0, 0.57802512272295603, 2.0, 2.0, 0.42197487727704402, 0.00323991087859726, 0.0029698681831674002, 0.57214765100671106, 0.58358021033920904 },
            new double[] { 2.0, 1.0, 0.67904096985110596, 3.0, 1.0, 0.32095903014889399, 0.015228814911237801, 0.013369937771056301, 0.652662357648509, 0.70637329286798101 },
            new double[] { 1.0, 1.0, 0.33952048492555298, 3.0, 1.0, 0.66047951507444702, 0.00761440745561906, 0.0068019841476156797, 0.326331178824255, 0.35203461888855198 },
            new double[] { 2.0, 1.0, 0.83952048492555298, 4.0, 1.0, 0.16047951507444699, 0.00761440745561887, 0.0066334300041643503, 0.826331178824255, 0.85203461888855203 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
        ["s1n500"] = (PVal: 2.1020858614373501e-05, Ntot: 2.9999989778315501, Solutions: new[]
        {
            new double[] { 2.0, 0.0, 0.15605024544591201, 2.0, 1.0, 0.84394975455408805, 0.0064798217571942303, 0.0060085874798433301, 0.146296020271278, double.NaN },
            new double[] { 2.0, 0.0, 0.57802512272295603, 2.0, 2.0, 0.42197487727704402, 0.00323991087859726, 0.0028224199939346399, 0.573647742512293, double.NaN },
            new double[] { 2.0, 1.0, 0.67904096985110596, 3.0, 1.0, 0.32095903014889399, 0.015228814911237801, 0.0141561888468551, 0.65502922177791401, double.NaN },
            new double[] { 1.0, 1.0, 0.33952048492555298, 3.0, 1.0, 0.66047951507444702, 0.00761440745561906, 0.0068654071560284996, 0.327514610888957, double.NaN },
            new double[] { 2.0, 1.0, 0.83952048492555298, 4.0, 1.0, 0.16047951507444699, 0.00761440745561887, 0.0069727341025116896, 0.82751461088895695, double.NaN },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
        ["s1n20"] = (PVal: 2.1020858614373501e-05, Ntot: 2.9999989778315501, Solutions: new[]
        {
            new double[] { 2.0, 0.0, 0.15605024544591201, 2.0, 1.0, 0.84394975455408805, 0.0064798217571942303, 0.0084178392254478495, double.NaN, double.NaN },
            new double[] { 2.0, 0.0, 0.57802512272295603, 2.0, 2.0, 0.42197487727704402, 0.00323991087859726, 0.00323009097983953, double.NaN, double.NaN },
            new double[] { 2.0, 1.0, 0.67904096985110596, 3.0, 1.0, 0.32095903014889399, 0.015228814911237801, 0.0126758378895141, double.NaN, double.NaN },
            new double[] { 1.0, 1.0, 0.33952048492555298, 3.0, 1.0, 0.66047951507444702, 0.00761440745561906, 0.0069929824182679403, double.NaN, double.NaN },
            new double[] { 2.0, 1.0, 0.83952048492555298, 4.0, 1.0, 0.16047951507444699, 0.00761440745561887, 0.0075944341733923497, double.NaN, double.NaN },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
        ["s2"] = (PVal: 0, Ntot: 2.9999989778315501, Solutions: new[]
        {
            new double[] { 2.0, 0.0, 0.15605024544591201, 2.0, 1.0, 0.84394975455408805, double.NaN, 0.0, 0.15626393223361601, 0.15626393223361601 },
            new double[] { 2.0, 0.0, 0.57802512272295603, 2.0, 2.0, 0.42197487727704402, double.NaN, 0.0, 0.57813196611680795, 0.57813196611680795 },
            new double[] { 2.0, 1.0, 0.67904096985110596, 3.0, 1.0, 0.32095903014889399, double.NaN, 0.0, 0.67853867318862904, 0.67853867318862904 },
            new double[] { 1.0, 1.0, 0.33952048492555298, 3.0, 1.0, 0.66047951507444702, double.NaN, 0.0, 0.33926933659431502, 0.33926933659431502 },
            new double[] { 2.0, 1.0, 0.83952048492555298, 4.0, 1.0, 0.16047951507444699, double.NaN, 0.0, 0.83926933659431402, 0.83926933659431402 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
        ["s3"] = (PVal: 0.036949670267513897, Ntot: 3.16431994277464, Solutions: new[]
        {
            new double[] { 3.0, 0.0, 0.84149606261457699, 3.0, 1.0, 0.15850393738542301, 0.0335351207441906, 0.027099948479094901, 0.78233132237708902, 0.89740146907827201 },
            new double[] { 3.0, 0.0, 0.920748031307288, 3.0, 2.0, 0.0792519686927117, 0.0167675603720953, 0.013998323303221399, 0.89116566118854401, 0.94870073453913595 },
            new double[] { 4.0, 0.0, 0.58569741234696604, 4.0, 1.0, 0.41430258765303402, 0.042318128558144603, 0.035189054225982799, 0.51103714490442098, 0.65624471097972403 },
            new double[] { 4.0, 0.0, 0.79284870617348302, 4.0, 2.0, 0.20715129382651701, 0.021159064279072302, 0.0178105869527482, 0.75551857245220999, 0.82812235548986202 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
        ["s4"] = (PVal: 0.024885230920557801, Ntot: 2.4786006635438498, Solutions: new[]
        {
            new double[] { 2.0, 0.0, 0.554567052328231, 3.0, 0.0, 0.445432947671769, 0.14517451315980001, 0.12696742413780401, 0.28061684037572598, 0.75502869611456302 },
            new double[] { 1.0, 0.0, 0.277283526164115, 3.0, 0.0, 0.722716473835885, 0.072587256579900103, 0.064188520491898299, 0.14030842018786299, 0.38476312419974301 },
            new double[] { 2.0, 0.0, 0.777283526164116, 4.0, 0.0, 0.222716473835884, 0.072587256579899798, 0.066282075799958304, 0.63852438722455995, 0.89764697311867003 },
            new double[] { 3.0, 0.0, 0.76062937952524601, 3.0, 1.0, 0.23937062047475399, 0.070707155684898604, 0.061797589675289401, 0.65037812681791796, 0.886126385212617 },
            new double[] { 3.0, 0.0, 0.88031468976262295, 3.0, 2.0, 0.119685310237377, 0.035353577842449198, 0.032248357762904599, 0.81808624472784697, 0.943063192606308 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
        ["s5"] = (PVal: 2.1343911706822401e-09, Ntot: 3.0000010382213498, Solutions: new[]
        {
            new double[] { 2.0, 1.0, 0.70329308054184703, 3.0, 1.0, 0.29670691945815297, 0.0187275272879537, 0.018022906735212799, 0.66620363939436, 0.73567095563232898 },
            new double[] { 1.0, 1.0, 0.35164654027092301, 3.0, 1.0, 0.64835345972907699, 0.0093637636439768395, 0.0088868975252763202, 0.33203602802179499, 0.366858107188546 },
            new double[] { 2.0, 1.0, 0.85164654027092301, 4.0, 1.0, 0.14835345972907699, 0.0093637636439767892, 0.0090525980304564605, 0.83096887256947405, 0.86721365009388496 },
            new double[] { 2.0, 0.0, 0.12918797646508301, 2.0, 1.0, 0.87081202353491705, 0.0071004834457107898, 0.0068333264277050301, 0.117290252837667, 0.14438693557023 },
            new double[] { 2.0, 0.0, 0.56459398823254103, 2.0, 2.0, 0.43540601176745902, 0.0035502417228554001, 0.0034994784451735901, 0.558887461119985, 0.57232936361581399 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        }),
    };

    private static void AssertSolutions(IReadOnlyList<OncologyAnalyzer.BattenbergSubclonalSolution> actual, double[][] expected, string label)
    {
        Assert.That(actual, Has.Count.EqualTo(expected.Length), $"{label}: six solutions A–F.");
        for (int k = 0; k < expected.Length; k++)
        {
            var a = actual[k];
            var e = expected[k];
            string s = $"{label} solution {(char)('A' + k)}";
            int? Cn(double v) => double.IsNaN(v) ? null : (int)v;
            Assert.That((a.MajorCopyNumber1, a.MinorCopyNumber1, a.MajorCopyNumber2, a.MinorCopyNumber2),
                Is.EqualTo((Cn(e[0]), Cn(e[1]), Cn(e[3]), Cn(e[4]))), $"{s} nMaj1/nMin1/nMaj2/nMin2.");
            AssertNum(a.Fraction1, e[2], 1e-12, $"{s} frac1.");
            AssertNum(a.Fraction2, e[5], 1e-12, $"{s} frac2.");
            AssertNum(a.FractionSd, e[6], 1e-12, $"{s} SDfrac.");
            AssertNum(a.FractionBootstrapSd, e[7], 1e-12, $"{s} SDfrac_BS (R resamples reproduced).");
            AssertNum(a.Fraction1Lower, e[8], 1e-12, $"{s} frac1_0.025.");
            AssertNum(a.Fraction1Upper, e[9], 1e-12, $"{s} frac1_0.975.");
        }
    }

    private static void AssertNum(double actual, double expected, double tol, string message)
    {
        if (double.IsNaN(expected))
        {
            Assert.That(actual, Is.NaN, message + " (R NA)");
        }
        else
        {
            Assert.That(actual, Is.EqualTo(expected).Within(tol), message);
        }
    }

    private static OncologyAnalyzer.BattenbergSegmentCall Run(string name)
    {
        var (rho, psit, logR, bafSeg, snps, seed, perms) = Inputs[name];
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 1000L * snps.Length, logR, bafSeg, snps.Length);
        return OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(
            new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, snps) }, rho, psit, seed, permutations: perms)[0];
    }

    // s1: (2,0)+(2,1) with x = 0 ⇒ solution F needs nMin = −1 ⇒ NA, moved last. s1n500: 975th order statistic missing ⇒
    // frac1_0.975 NA; s1n20: both NA. s2: one SNP ⇒ SDfrac NA (sd of one value), bootstrap constant (SD_BS 0, bounds =
    // τ(SNP BAF) ≠ τ(l)). s3–s5: other genomes/edges; seeds −12345, 2³¹ − 1 and 0 exercise R's unsigned seeding.
    [TestCase("s1")]
    [TestCase("s1n500")]
    [TestCase("s1n20")]
    [TestCase("s2")]
    [TestCase("s3")]
    [TestCase("s4")]
    [TestCase("s5")]
    public void FitSubclonalCopyNumberWithBootstrap_MatchesBattenberg(string name)
    {
        var call = Run(name);
        var (pval, ntot, solutions) = Expected[name];
        Assert.Multiple(() =>
        {
            Assert.That(call.PValue, Is.EqualTo(pval).Within(1e-12 * Math.Max(pval, 1e-300)), "pval.");
            Assert.That(call.TotalCopyNumber, Is.EqualTo(ntot).Within(1e-12), "ntot.");
            Assert.That(call.Fit.IsSubclonal, Is.True);
            Assert.That(call.Baf, Is.EqualTo(Inputs[name].BafSeg).Within(1e-15), "BAF column = l.");
            AssertSolutions(call.Solutions, solutions, name);
            Assert.That((call.Fit.PrimaryState.MajorCopyNumber, call.Fit.PrimaryState.MinorCopyNumber, call.Fit.PrimaryState.CellFraction),
                Is.EqualTo((call.Solutions[0].MajorCopyNumber1!.Value, call.Solutions[0].MinorCopyNumber1!.Value, call.Solutions[0].Fraction1)),
                "Fit = solution A.");
        });
    }

    // Multi-segment end to end (F40 tracks): R set.seed once, then the RNG stream runs through the sub-clonal segments in
    // order (clonal segments draw nothing). e1: t1, ρ 0.7, ψ 2.6, seed 4711; e2: t2, ρ 0.85, ψ 3, seed 99. null = clonal.
    private static readonly double[][]?[] E1 =
    {
        null,
        new[]
        {
            new double[] { 2.0, 0.0, 0.47709449087559203, 2.0, 1.0, 0.52290550912440803, 0.020956394990033399, 0.021882553818262800, 0.42786173000326799, 0.51237317814035799 },
            new double[] { 2.0, 0.0, 0.73854724543779604, 2.0, 2.0, 0.26145275456220401, 0.010478197495016699, 0.010552093851411599, 0.71505178985032802, 0.75609582488656302 },
            new double[] { 3.0, 0.0, 0.085309869471425395, 3.0, 1.0, 0.91469013052857495, 0.029585498809459099, 0.029926746827075499, 0.016801363005953002, 0.13825288451664100 },
            new double[] { 3.0, 0.0, 0.54265493473571302, 3.0, 2.0, 0.45734506526428698, 0.014792749404729600, 0.015203296946707899, 0.50978959048207795, 0.56868633581708905 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        },
        new[]
        {
            new double[] { 1.0, 1.0, 0.37649750519839498, 2.0, 1.0, 0.62350249480160502, 0.084545436246698594, 0.068788195074560604, 0.56958083011405403, 0.83834293480680000 },
            new double[] { 0.0, 1.0, 0.18824875259919699, 2.0, 1.0, 0.81175124740080296, 0.042272718123349297, 0.033279011590574600, 0.28606236361796999, 0.41552414097986701 },
            new double[] { 1.0, 1.0, 0.688248752599197045, 3.0, 1.0, 0.31175124740080301, 0.042272718123349200, 0.034937582522478298, 0.781436876640085987, 0.91390969753780604 },
            new double[] { 1.0, 0.0, 0.43405738924477399, 1.0, 1.0, 0.56594261075522601, 0.040961526058712103, 0.047504763487642897, 0.14162418377875999, 0.32924918954892302 },
            new double[] { 1.0, 0.0, 0.71702869462238705, 1.0, 2.0, 0.282971305377613, 0.020480763029355999, 0.0232999982326237, 0.57692563543023201, 0.66584114988686005 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        },
    };

    private static readonly double[][]?[] E2 =
    {
        new[]
        {
            new double[] { 1.0, 1.0, 0.75731976285090297, 2.0, 1.0, 0.24268023714909701, 0.015061041969408301, 0.014678822762364801, 0.75770698279828796, 0.81580436169284898 },
            new double[] { 0.0, 1.0, 0.37865988142545098, 2.0, 1.0, 0.62134011857454896, 0.0075305209847041600, 0.0075680234658234199, 0.37830232085880700, 0.40793077466171801 },
            new double[] { 1.0, 1.0, 0.87865988142545104, 3.0, 1.0, 0.12134011857454900, 0.007530520984704140, 0.0072117411591522901, 0.87866760222262297, 0.90734436514755801 },
            new double[] { 2.0, 1.0, 0.37218482282030302, 2.0, 2.0, 0.62781517717969704, 0.0191481673448746989, 0.0196009688430185992, 0.29401265185485298, 0.36860086786568902 },
            new double[] { 2.0, 0.0, 0.18609241141015101, 2.0, 2.0, 0.81390758858984902, 0.0095740836724374796, 0.0099138104764448804, 0.14749891915848201, 0.18615274864469600 },
            new double[] { 2.0, 1.0, 0.68609241141015198, 2.0, 3.0, 0.31390758858984802, 0.0095740836724375108, 0.0093505892611737296, 0.64852167292974905, 0.68560634338449000 },
        },
        new[]
        {
            new double[] { 3.0, 0.0, 0.39163557478371203, 3.0, 1.0, 0.60836442521628797, 0.011898571826941300, 0.012333675825412600, 0.35898235699402198, 0.40742306607074202 },
            new double[] { 3.0, 0.0, 0.69581778739185596, 3.0, 2.0, 0.30418221260814399, 0.0059492859134706597, 0.0059142331131774896, 0.67995895525010197, 0.70263577463408700 },
            new double[] { 4.0, 0.0, 0.14455788536376901, 4.0, 1.0, 0.85544211463623099, 0.015644418513200399, 0.0154702257345826003, 0.10311127224337500, 0.16433940578928100 },
            new double[] { 4.0, 0.0, 0.57227894268188495, 4.0, 2.0, 0.42772105731811499, 0.0078222092566003192, 0.0078875839380796505, 0.55099979694716905, 0.58188842620113701 },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
            new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN },
        },
        null,
        null,
        new[]
        {
            new double[] { 2.0, 1.0, 0.79257756414276204, 3.0, 1.0, 0.20742243585723799, 0.024048213929423200, 0.023398839970892600, 0.73445587546712199, 0.82921530038811497 },
            new double[] { 1.0, 1.0, 0.39628878207138102, 3.0, 1.0, 0.60371121792861904, 0.0120241069647111006, 0.0118837785489453000, 0.36677802735613002, 0.41382941787706201 },
            new double[] { 2.0, 1.0, 0.89628878207138096, 4.0, 1.0, 0.10371121792861900, 0.012024106964710899, 0.0123875232515336005, 0.86647173790351695, 0.91645764257113305 },
            new double[] { 3.0, 1.0, 0.60885652850799299, 3.0, 2.0, 0.39114347149200701, 0.0158132168621975995, 0.0156492862692958999, 0.58572309274979795, 0.64682203102866498 },
            new double[] { 3.0, 0.0, 0.30442826425399699, 3.0, 2.0, 0.69557173574600295, 0.0079066084310987807, 0.0076008234034460999, 0.29234195869592100, 0.32246888112922101 },
            new double[] { 3.0, 1.0, 0.80442826425399605, 3.0, 3.0, 0.19557173574600401, 0.0079066084310988605, 0.0078415937866751200, 0.79284517555004597, 0.82326311601009605 },
        },
        null,
    };

    [TestCase("e1")]
    [TestCase("e2")]
    public void FitSubclonalCopyNumberWithBootstrap_MultiSegmentStream_MatchesBattenberg(string name)
    {
        bool first = name == "e1";
        var rows = OncologyAnalyzer.SegmentPhasedBaf(OncologyAnalyzer_BattenbergPhasedSegmentation_Tests.Track(first ? "t1" : "t2"));
        var logR = first
            ? OncologyAnalyzer_BattenbergPhasedSegmentation_Tests.LogRTrack(rows, new[] { 0.0, 0.15, -0.2 }, new[] { 100, 200, 300 })
            : OncologyAnalyzer_BattenbergPhasedSegmentation_Tests.LogRTrack(rows, new[] { -0.1, 0.3, 0.0, 0.05, 0.2 }, new[] { 500, 900, 1200, 1500, 2000 });
        if (first)
        {
            logR.Add(new OncologyAnalyzer.LogRProbe("1", 150500, double.PositiveInfinity));
        }

        var segments = OncologyAnalyzer.BuildBattenbergSegments(rows, logR);
        var calls = first
            ? OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(segments, 0.7, 2.6, seed: 4711)
            : OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(segments, 0.85, 3.0, seed: 99);
        var expected = first ? E1 : E2;
        Assert.That(calls, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                if (expected[k] is null)
                {
                    Assert.That(calls[k].Solutions, Is.Empty, $"segment {k + 1} clonal ⇒ NA solutions.");
                }
                else
                {
                    AssertSolutions(calls[k].Solutions, expected[k]!, $"{name} segment {k + 1}");
                }
            }
        });
    }

    // Statistical agreement (seed-independent): R 4.3.3, segment s4, 400 seeds (1..400) of set.seed + determine_copynumber:
    // mean (sd) SDfrac_A_BS 0.130353 (0.002746), frac1_A_0.025 0.280433 (0.007417), frac1_A_0.975 0.785209 (0.012439).
    // 100 C# seeds disjoint from R's must agree within 5 Monte-Carlo standard errors of the difference of means.
    [Test]
    public void FitSubclonalCopyNumberWithBootstrap_OtherSeeds_AgreeWithRWithinMonteCarloError()
    {
        var (rho, psit, logR, bafSeg, snps, _, _) = Inputs["s4"];
        var seg = new OncologyAnalyzer.SubclonalSegmentSnpBafs(
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, logR, bafSeg, snps.Length), snps);
        const int runs = 100;
        double sd = 0, lo = 0, hi = 0;
        for (int s = 0; s < runs; s++)
        {
            var a = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(new[] { seg }, rho, psit, seed: 100_001 + s)[0].Solutions[0];
            sd += a.FractionBootstrapSd / runs;
            lo += a.Fraction1Lower / runs;
            hi += a.Fraction1Upper / runs;
        }

        double Tol(double rSd) => 5 * rSd * Math.Sqrt((1.0 / runs) + (1.0 / 400));
        Assert.Multiple(() =>
        {
            Assert.That(sd, Is.EqualTo(0.130353).Within(Tol(0.002746)), "mean SDfrac_A_BS.");
            Assert.That(lo, Is.EqualTo(0.280433).Within(Tol(0.007417)), "mean frac1_A_0.025.");
            Assert.That(hi, Is.EqualTo(0.785209).Within(Tol(0.012439)), "mean frac1_A_0.975.");
            // The delta-method SDfrac (deterministic) is the same whatever the seed.
            Assert.That(Run("s4").Solutions[0].FractionSd, Is.EqualTo(0.14517451315980001).Within(1e-12));
        });
    }

    [Test]
    public void FitSubclonalCopyNumberWithBootstrap_SameSeedReproducible_ClonalHasNoSolutions_Guards()
    {
        var a = Run("s1");
        var b = Run("s1");
        var clonalSeg = new OncologyAnalyzer.SubclonalSegmentSnpBafs(
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1, 3, 0.222392, 0.647857, 3), new[] { 0.6379, 0.6579, 0.6484 });
        var clonal = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(new[] { clonalSeg }, 0.8, 2.5, seed: 1)[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.Solutions.Select(x => x.FractionBootstrapSd), Is.EqualTo(b.Solutions.Select(x => x.FractionBootstrapSd)), "deterministic per seed.");
            Assert.That(clonal.Fit.IsSubclonal, Is.False);
            Assert.That(clonal.PValue, Is.EqualTo(1.0), "within maxdist.");
            Assert.That(clonal.Solutions, Is.Empty, "clonal ⇒ Battenberg NA columns.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(new[] { clonalSeg }, 0.8, 2.5, 1, permutations: 0));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(null!, 0.8, 2.5, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(new[] { clonalSeg }, 0.0, 2.5, 1));
        });
    }
}
