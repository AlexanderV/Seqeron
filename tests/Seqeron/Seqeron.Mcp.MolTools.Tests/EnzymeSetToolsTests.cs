using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Review 2026-09 B08 D4 part 2 (WP12b, F61-F63): the optional <c>enzyme_set</c> argument
/// ("common" default | "rebase") of the enzyme-LISTING MCP tools, the default output locked
/// byte-identical to the pre-existing common-set path, and the JSON form of a two-cut enzyme's
/// second cut pair. Reference values are Biopython 1.88's (Bio.Restriction, REBASE emboss 404) over
/// the 754 enzymes with a defined cut (<c>fst5 is not None</c>): <c>RestrictionBatch.search</c>,
/// <c>compatible_end</c>, <c>is_blunt</c>, <c>.size</c>, <c>equischizomers</c>/<c>neoschizomers</c>,
/// BaeI <c>scd5 = 23</c>, <c>size + scd3 = 11 + 7 = 18</c>.
/// </summary>
[TestFixture]
public class EnzymeSetToolsTests
{
    private const string Probe = "GAATTCGGTCTCAAGCTTGGATCCGAGACCCTGCAGGCGGCCGCAAAACCGGGATATCAGATCT";

    // Biopython 1.88: sorted(str(k) for k, v in RestrictionBatch(defined).search(Seq(Probe), linear=True).items() if v)
    private static readonly string[] BiopythonCuttersLinear =
    {
        "AbaSI", "AciI", "AclWI", "AcoI", "AcsI", "AfiI", "AluBI", "AluI", "Alw26I", "AlwI", "AoxI", "ApoI",
        "Asi256I", "AspBHI", "AsuC2I", "BamHI", "BcnI", "BcoDI", "BfmI", "BglII", "BinI", "BisI", "BlsI", "Bme1390I",
        "BmeDI", "BmiI", "BmrFI", "BpuMI", "BsaI", "Bsc4I", "BscXI", "BseLI", "BseX3I", "Bsh1285I", "BshFI", "BsiEI",
        "BsiSI", "BsiYI", "BslI", "BsmAI", "BsnI", "Bso31I", "Bsp143I", "BspACI", "BspANI", "BspLI", "BspMAI",
        "BspPI", "BspTNI", "BssMI", "BstC8I", "BstKTI", "BstMAI", "BstMBI", "BstMCI", "BstSCI", "BstSFI", "BstX2I",
        "BstYI", "BstZI", "BsuRI", "BthCI", "Cac8I", "CauII", "CciNI", "CfrI", "ChaI", "CviJI", "CviKI_1", "CviRI",
        "DpnI", "DpnII", "EaeI", "EagI", "EclXI", "Eco31I", "Eco32I", "Eco52I", "EcoBLMcrX", "EcoHI", "EcoO157SI",
        "EcoRI", "EcoRV", "Fnu4HI", "Fsp4HI", "FspEI", "GdiII", "GluI", "HaeIII", "HapII", "Hin4I", "HindIII",
        "HpaII", "Hpy188I", "HpyCH4V", "Kzo9I", "LpnPI", "MalI", "MboI", "McrI", "MflI", "MluCI", "MspI", "MspJI",
        "MspR9I", "NciI", "NdeII", "NgoAVII", "NlaIV", "NotI", "Pfl8569I", "PkrI", "PspN4I", "PstI", "PsuI", "SatI",
        "Sau3AI", "SbfI", "ScrFI", "SdaI", "SetI", "SfcI", "SfeI", "SgeI", "SgrTI", "SimI", "SmlI", "SmoI",
        "Sse8387I", "Sse9I", "SsiI", "Sth132I", "Sth302II", "StyD4I", "TasI", "TauI", "TspEI", "XapI", "XhoII",
        "XmaIII", "YkrI",
    };

    // Same, linear=False.
    private static readonly string[] BiopythonCuttersCircular =
    {
        "AbaSI", "AciI", "AclWI", "AcoI", "AcsI", "AfiI", "AluBI", "AluI", "Alw26I", "AlwI", "AoxI", "ApoI",
        "Asi256I", "AspBHI", "AsuC2I", "BamHI", "Bce83I", "BcnI", "BcoDI", "BfmI", "BglII", "BinI", "BisI", "BlsI",
        "Bme1390I", "BmeDI", "BmiI", "BmrFI", "BpuEI", "BpuMI", "BsaI", "Bsc4I", "BscXI", "BseLI", "BseX3I",
        "Bsh1285I", "BshFI", "BsiEI", "BsiSI", "BsiYI", "BslI", "BsmAI", "BsnI", "Bso31I", "Bsp143I", "BspACI",
        "BspANI", "BspLI", "BspMAI", "BspPI", "BspTNI", "BssMI", "BstC8I", "BstKTI", "BstMAI", "BstMBI", "BstMCI",
        "BstSCI", "BstSFI", "BstX2I", "BstYI", "BstZI", "BsuRI", "BthCI", "Cac8I", "CauII", "CciNI", "CfrI", "ChaI",
        "CviJI", "CviKI_1", "CviRI", "DpnI", "DpnII", "EaeI", "EagI", "EclXI", "Eco31I", "Eco32I", "Eco52I",
        "EcoBLMcrX", "EcoHI", "EcoO157SI", "EcoRI", "EcoRV", "Fnu4HI", "Fsp4HI", "FspEI", "GdiII", "GluI", "HaeIII",
        "HapII", "Hin4I", "HindIII", "HpaII", "Hpy188I", "HpyCH4V", "Kzo9I", "LpnPI", "MalI", "MboI", "McrI", "MflI",
        "MluCI", "MspI", "MspJI", "MspR9I", "NciI", "NdeII", "NgoAVII", "NlaIV", "NotI", "Pfl8569I", "PkrI",
        "PspN4I", "PstI", "PsuI", "SatI", "Sau3AI", "SbfI", "ScrFI", "SdaI", "SetI", "SfcI", "SfeI", "SgeI", "SgrTI",
        "SimI", "SmlI", "SmoI", "Sse8387I", "Sse9I", "SsiI", "Sth132I", "Sth302II", "StyD4I", "TaqII", "TasI",
        "TauI", "TspEI", "XapI", "XhoII", "XmaIII", "YkrI",
    };

    // Biopython 1.88: sorted(str(e) for e in BamHI.compatible_end(RestrictionBatch(defined)))
    private static readonly string[] BiopythonBamHICompatible =
    {
        "AarI", "Acc36I", "AceIII", "Alw26I", "AspBHI", "BamHI", "BbsI", "BbvI", "BbvII", "BclI", "BcoDI", "BfuAI",
        "BglII", "BmsI", "BpiI", "BsaI", "BseXI", "BslFI", "BsmAI", "BsmBI", "BsmFI", "Bso31I", "Bsp143I", "BspMI",
        "BspTNI", "BssMI", "BstMAI", "BstMBI", "BstV1I", "BstV2I", "BstX2I", "BstYI", "BtgZI", "BveI", "DpnII",
        "Eco31I", "Esp3I", "FaqI", "FbaI", "FokI", "FspEI", "Ksp22I", "Kzo9I", "LpnPI", "Lsp1109I", "LweI", "MboI",
        "MflI", "MspJI", "NdeII", "PaqCI", "PsuI", "Sau3AI", "SfaNI", "SgeI", "SgrTI", "Sth132I", "StsI", "XhoII",
    };

    // Biopython 1.88: sorted(str(e) for e in EcoRI.compatible_end(RestrictionBatch(defined)))
    private static readonly string[] BiopythonEcoRICompatible =
    {
        "AarI", "Acc36I", "AceIII", "AcsI", "Alw26I", "ApoI", "AspBHI", "BbsI", "BbvI", "BbvII", "BcoDI", "BfuAI",
        "BmsI", "BpiI", "BsaI", "BseXI", "BslFI", "BsmAI", "BsmBI", "BsmFI", "Bso31I", "BspMI", "BspTNI", "BstMAI",
        "BstV1I", "BstV2I", "BtgZI", "BveI", "Eco31I", "EcoRI", "Esp3I", "FaqI", "FokI", "FspEI", "LpnPI",
        "Lsp1109I", "LweI", "MfeI", "MluCI", "MspJI", "MunI", "PaqCI", "SfaNI", "SgeI", "SgrTI", "Sse9I", "Sth132I",
        "StsI", "TasI", "TspEI", "XapI",
    };

    private static string Json<T>(T value) => JsonSerializer.Serialize(value, ModelContextProtocol.McpJsonUtilities.DefaultOptions);

    private static string[] Sorted(IEnumerable<string> names) => names.Distinct().Order(StringComparer.Ordinal).ToArray();

    #region Default ("common") output is byte-identical to the pre-existing path

    [Test]
    public void DefaultEnzymeSet_OutputIsByteIdenticalToTheCommonSetPath()
    {
        var dna = new DnaSequence(Probe);
        Assert.Multiple(() =>
        {
            Assert.That(Json(MolToolsTools.enzymes_by_cut_length(6)),
                Is.EqualTo(Json(new Models.EnzymeListResult(RestrictionAnalyzer.GetEnzymesByCutLength(6).ToList()))));
            Assert.That(Json(MolToolsTools.blunt_cutters()),
                Is.EqualTo(Json(new Models.EnzymeListResult(RestrictionAnalyzer.GetBluntCutters().ToList()))));
            Assert.That(Json(MolToolsTools.sticky_cutters()),
                Is.EqualTo(Json(new Models.EnzymeListResult(RestrictionAnalyzer.GetStickyCutters().ToList()))));
            Assert.That(Json(MolToolsTools.find_all_restriction_sites(Probe, MoleculeTopology.Circular)),
                Is.EqualTo(Json(new Models.RestrictionSiteListResult(RestrictionAnalyzer.FindAllSites(dna, MoleculeTopology.Circular).ToList()))));
            Assert.That(Json(MolToolsTools.restriction_map(Probe, System.Array.Empty<string>())),
                Is.EqualTo(Json(RestrictionAnalyzer.CreateMap(dna))));
            // Each pair also reports end2, the SECOND enzyme's own overhang window — equal to
            // compatibleEnd here, and different only for the 12 REBASE pairs that anneal through the
            // reverse complement (review 2026-09 B08 H30).
            Assert.That(Json(MolToolsTools.compatible_enzymes("BamHI")),
                Is.EqualTo(Json(new Models.CompatibleEnzymesResult(RestrictionAnalyzer.FindCompatibleEnzymes("BamHI")
                    .Select(t => new Models.EnzymeCompatibilityPair(
                        "BamHI", t.Enzyme, t.CompatibleEnd, OtherEnd(t.Enzyme, t.CompatibleEnd))).ToList()))));
            // The all-pairs listing now also reports totalPairs / truncated (F72); its pairs are unchanged.
            Assert.That(Json(MolToolsTools.compatible_enzymes()),
                Is.EqualTo(Json(new Models.CompatibleEnzymesResult(RestrictionAnalyzer.FindCompatibleEnzymes()
                    .Select(t => new Models.EnzymeCompatibilityPair(
                        t.Enzyme1, t.Enzyme2, t.CompatibleEnd, OtherEnd(t.Enzyme2, t.CompatibleEnd))).ToList(), 64, false))));
            Assert.That(MolToolsTools.get_enzyme("MspI").Equischizomers, Is.EqualTo(new[] { "HpaII" }));
            Assert.That(MolToolsTools.blunt_cutters().Enzymes, Has.Count.EqualTo(10));
            Assert.That(MolToolsTools.sticky_cutters().Enzymes, Has.Count.EqualTo(29));
            Assert.That(MolToolsTools.enzymes_by_cut_length(6).Enzymes, Has.Count.EqualTo(24));
        });
    }

    [Test]
    public void ExplicitCommon_AnyCase_EqualsTheDefault()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Json(MolToolsTools.blunt_cutters("COMMON")), Is.EqualTo(Json(MolToolsTools.blunt_cutters())));
            Assert.That(Json(MolToolsTools.enzymes_by_cut_length(4, "common")), Is.EqualTo(Json(MolToolsTools.enzymes_by_cut_length(4))));
            Assert.That(Json(MolToolsTools.get_enzyme("MspI", "Common")), Is.EqualTo(Json(MolToolsTools.get_enzyme("MspI"))));
            Assert.That(Json(MolToolsTools.find_all_restriction_sites(Probe, MoleculeTopology.Linear, null)),
                Is.EqualTo(Json(MolToolsTools.find_all_restriction_sites(Probe))));
        });
    }

    #endregion

    #region enzyme_set = "rebase" — Biopython-locked

    [Test]
    public void Rebase_ListingTools_MatchBiopythonCounts()
    {
        Assert.Multiple(() =>
        {
            // Biopython over the 754 defined-cut enzymes: is_blunt() 137, not blunt 617, .size == 6 → 385.
            Assert.That(MolToolsTools.blunt_cutters("rebase").Enzymes, Has.Count.EqualTo(137));
            Assert.That(MolToolsTools.sticky_cutters("rebase").Enzymes, Has.Count.EqualTo(617));
            Assert.That(MolToolsTools.enzymes_by_cut_length(6, "Rebase").Enzymes, Has.Count.EqualTo(385));
            Assert.That(MolToolsTools.blunt_cutters("rebase").Enzymes.All(e => e.IsBluntEnd), Is.True);
        });
    }

    [Test]
    public void Rebase_GetEnzyme_SchizomerListsMatchBiopython()
    {
        var bsaI = MolToolsTools.get_enzyme("BsaI", "rebase");
        var smaI = MolToolsTools.get_enzyme("SmaI", "rebase");

        Assert.Multiple(() =>
        {
            // BsaI.equischizomers() == [Bso31I, BspTNI, Eco31I], neoschizomers() == [];
            // SmaI.neoschizomers() == [Cfr9I, TspMI, XmaI], equischizomers() == [].
            Assert.That(Sorted(bsaI.Equischizomers!), Is.EqualTo(new[] { "Bso31I", "BspTNI", "Eco31I" }));
            Assert.That(bsaI.Neoschizomers, Is.Empty);
            Assert.That(Sorted(smaI.Neoschizomers!), Is.EqualTo(new[] { "Cfr9I", "TspMI", "XmaI" }));
            Assert.That(smaI.Equischizomers, Is.Empty);
            // Default (common) lists: no BsaI schizomer among the 39.
            Assert.That(MolToolsTools.get_enzyme("BsaI").Equischizomers, Is.Empty);
        });
    }

    [Test]
    public void Rebase_FindAllRestrictionSites_EnzymesWithASite_MatchBiopythonSearch()
    {
        var linear = MolToolsTools.find_all_restriction_sites(Probe, MoleculeTopology.Linear, "rebase");
        var circular = MolToolsTools.find_all_restriction_sites(Probe, MoleculeTopology.Circular, "rebase");

        Assert.Multiple(() =>
        {
            Assert.That(Sorted(linear.Sites.Select(s => s.Enzyme.Name)), Is.EqualTo(BiopythonCuttersLinear));
            Assert.That(Sorted(circular.Sites.Select(s => s.Enzyme.Name)), Is.EqualTo(BiopythonCuttersCircular));
        });
    }

    [Test]
    public void Rebase_RestrictionMap_NoNames_PartitionsTheFullSet()
    {
        var map = MolToolsTools.restriction_map(Probe, System.Array.Empty<string>(), MoleculeTopology.Linear, "rebase");

        Assert.Multiple(() =>
        {
            Assert.That(Sorted(map.SitesByEnzyme.Keys), Is.EqualTo(BiopythonCuttersLinear));
            Assert.That(map.SitesByEnzyme.Count + map.NonCutters.Count, Is.EqualTo(754));
            Assert.That(map.NonCutters, Does.Contain("BaeI"));
        });
    }

    [Test]
    public void Rebase_RestrictionMap_WithNames_IgnoresTheSet()
    {
        var names = new[] { "EcoRI", "BsaI" };
        Assert.That(Json(MolToolsTools.restriction_map(Probe, names, MoleculeTopology.Linear, "rebase")),
            Is.EqualTo(Json(MolToolsTools.restriction_map(Probe, names))));
    }

    [Test]
    public void Rebase_CompatibleEnzymes_SingleEnzyme_MatchesBiopythonCompatibleEnd()
    {
        var bamHI = MolToolsTools.compatible_enzymes("BamHI", "rebase");
        var ecoRI = MolToolsTools.compatible_enzymes("EcoRI", "rebase");

        Assert.Multiple(() =>
        {
            Assert.That(Sorted(bamHI.Pairs.Select(p => p.Enzyme2)), Is.EqualTo(BiopythonBamHICompatible));
            Assert.That(Sorted(ecoRI.Pairs.Select(p => p.Enzyme2)), Is.EqualTo(BiopythonEcoRICompatible));
            Assert.That(bamHI.Pairs.Single(p => p.Enzyme2 == "BclI").CompatibleEnd, Is.EqualTo("GATC"));
        });
    }

    [Test]
    public void Rebase_CompatibleEnzymes_AllPairsForm_IsBoundedAndPaged()
    {
        // B08 round 2 E10 (F72): 27 869 compatible pairs among the 754 REBASE enzymes (Biopython-locked on the
        // defined-overhang subset in RestrictionAnalyzer_RebaseSet_Tests); the response is a page, not the list.
        var first = MolToolsTools.compatible_enzymes(null, "rebase");
        var last = MolToolsTools.compatible_enzymes(null, "REBASE", max_pairs: 10000, offset: 20000);
        var tail = MolToolsTools.compatible_enzymes(null, "rebase", max_pairs: 10, offset: 27865);
        var small = MolToolsTools.compatible_enzymes(null, "rebase", max_pairs: 3);
        var all = RestrictionAnalyzer.FindCompatibleEnzymes(EnzymeSet.Rebase).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(first.Pairs, Has.Count.EqualTo(1000));
            Assert.That(first.TotalPairs, Is.EqualTo(27869));
            Assert.That(first.Truncated, Is.True);
            Assert.That(first.Pairs[0],
                Is.EqualTo(new Models.EnzymeCompatibilityPair("AanI", "Acc16I", "blunt", "blunt")),
                "both ends are reported; blunt is shared by definition (B08 H30)");
            Assert.That(small.Pairs, Has.Count.EqualTo(3));
            Assert.That(last.Pairs, Has.Count.EqualTo(7869));
            Assert.That(last.Truncated, Is.False);
            Assert.That(tail.Pairs, Has.Count.EqualTo(4), "page shorter than max_pairs at the end");
            Assert.That(tail.Pairs.Select(p => (p.Enzyme1, p.Enzyme2, p.CompatibleEnd)), Is.EqualTo(all.Skip(27865)));
            Assert.That(tail.Truncated, Is.False);
            Assert.That(MolToolsTools.compatible_enzymes(null, "rebase", offset: 30000).Pairs, Is.Empty);

            // F42: the EcoRI/MfeI 5' AATT pair is in the listing (paged by offset).
            var pages = Enumerable.Range(0, 28).SelectMany(i =>
                MolToolsTools.compatible_enzymes(null, "rebase", max_pairs: 1000, offset: i * 1000).Pairs);
            Assert.That(pages.Select(p => (p.Enzyme1, p.Enzyme2, p.CompatibleEnd)), Is.EqualTo(all));
            Assert.That(all, Does.Contain(("EcoRI", "MfeI", "AATT")));
        });
    }

    [Test]
    public void CompatibleEnzymes_AllPairs_Common_ReportsTotalsAndIsOtherwiseUnchanged()
    {
        var common = MolToolsTools.compatible_enzymes();

        Assert.Multiple(() =>
        {
            Assert.That(common.Pairs, Has.Count.EqualTo(64));
            Assert.That(common.TotalPairs, Is.EqualTo(64));
            Assert.That(common.Truncated, Is.False);
            Assert.That(MolToolsTools.compatible_enzymes("BamHI").TotalPairs, Is.Null, "per-enzyme form is not paged");
            Assert.That(MolToolsTools.compatible_enzymes("BamHI").Truncated, Is.Null);
        });
    }

    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(10001)]
    public void CompatibleEnzymes_MaxPairsOutOfRange_Throws(int maxPairs) =>
        Assert.Throws<System.ArgumentOutOfRangeException>(() => MolToolsTools.compatible_enzymes(null, "rebase", maxPairs));

    [Test]
    public void CompatibleEnzymes_NegativeOffset_Throws() =>
        Assert.Throws<System.ArgumentOutOfRangeException>(() => MolToolsTools.compatible_enzymes(null, "common", null, -1));

    [Test]
    public void CompatibleEnzymes_PagingParametersWithEnzymeName_AreRejected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<System.ArgumentException>(
                () => MolToolsTools.compatible_enzymes("EcoRI", "rebase", max_pairs: 5))!.ParamName, Is.EqualTo("max_pairs"));
            Assert.That(Assert.Throws<System.ArgumentException>(
                () => MolToolsTools.compatible_enzymes("EcoRI", "rebase", offset: 5))!.ParamName, Is.EqualTo("offset"));
        });
    }

    #endregion

    #region Invalid enzyme_set

    [TestCase("all")]
    [TestCase("")]
    [TestCase("rebase404")]
    public void UnknownEnzymeSet_IsAnArgumentError_OnEveryTool(string value)
    {
        Assert.Multiple(() =>
        {
            AssertSetError(() => MolToolsTools.get_enzyme("EcoRI", value));
            AssertSetError(() => MolToolsTools.enzymes_by_cut_length(6, value));
            AssertSetError(() => MolToolsTools.blunt_cutters(value));
            AssertSetError(() => MolToolsTools.sticky_cutters(value));
            AssertSetError(() => MolToolsTools.find_all_restriction_sites(Probe, MoleculeTopology.Linear, value));
            AssertSetError(() => MolToolsTools.restriction_map(Probe, System.Array.Empty<string>(), MoleculeTopology.Linear, value));
            AssertSetError(() => MolToolsTools.compatible_enzymes("EcoRI", value));
        });

        static void AssertSetError(TestDelegate call)
        {
            var ex = Assert.Throws<System.ArgumentException>(call);
            Assert.That(ex!.ParamName, Is.EqualTo("enzyme_set"));
        }
    }

    #endregion

    #region Name-taking tools over the full set; second cut in JSON

    [Test]
    public void EnzymesCompatible_PartnerOutsideTheCommonSet_ReturnsTheSharedEnd()
    {
        // BclI (TGATCA, 1/5) is a REBASE enzyme, not one of the 39. Biopython: BclI in BamHI.compatible_end(...).
        var result = MolToolsTools.enzymes_compatible("BamHI", "BclI");

        Assert.That((result.Compatible, result.CompatibleEnd), Is.EqualTo((true, "GATC")));
    }

    [Test]
    public void GetEnzyme_TwoCutEnzyme_SecondCutIsSerialised()
    {
        var baeI = Json(MolToolsTools.get_enzyme("BaeI"));
        var ecoRI = Json(MolToolsTools.get_enzyme("EcoRI"));

        Assert.Multiple(() =>
        {
            // Biopython BaeI: scd5 = 23, size + scd3 = 11 + 7 = 18.
            Assert.That(baeI, Does.Contain("\"secondCutForward\":23"));
            Assert.That(baeI, Does.Contain("\"secondCutReverse\":18"));
            Assert.That(baeI, Does.Not.Contain("\"secondCut\":"));
            // One-cut enzyme: the SDK omits null properties, so its JSON is unchanged (no secondCut* key).
            Assert.That(ecoRI, Does.Not.Contain("secondCut"));
        });
    }

    [Test]
    public void FindRestrictionSites_TwoCutEnzyme_SiteJsonCarriesNoEmptyTupleObject()
    {
        // Before WP12b every site record serialised RestrictionSite.CutOffsets (a value tuple) as
        // "cutOffsets":{}. The cut of each record is identified by isSecondCut + the enzyme's offsets.
        string s = new string('A', 20) + "ACAAAAGTACC" + new string('T', 20);
        var json = Json(MolToolsTools.find_restriction_sites(s, new[] { "BaeI" }));

        Assert.Multiple(() =>
        {
            Assert.That(json, Does.Not.Contain("{}"));
            Assert.That(json, Does.Not.Contain("cutOffsets"));
            Assert.That(json, Does.Contain("\"isSecondCut\":true"));
            Assert.That(json, Does.Contain("\"secondCutForward\":23"));
        });
    }

    #endregion

    /// <summary>
    /// B08 H30: <c>compatibleEnd</c> is the FIRST enzyme's own overhang window, not a string both
    /// enzymes leave, and the new <c>end2</c> makes the partner's window visible. 12 REBASE pairs
    /// anneal only through the reverse complement and all 12 are <c>GN</c> against <c>CN</c>; the
    /// tool description used to say the value "echoes the shared end".
    /// </summary>
    [Test]
    public void EnzymesCompatible_ReverseComplementOnlyPair_ReportsBothEnds()
    {
        var forward = MolToolsTools.enzymes_compatible("BsmI", "BsrI");
        var reversed = MolToolsTools.enzymes_compatible("BsrI", "BsmI");
        var blunt = MolToolsTools.enzymes_compatible("EcoRV", "SmaI");
        var incompatible = MolToolsTools.enzymes_compatible("EcoRI", "BamHI");

        Assert.Multiple(() =>
        {
            Assert.That(forward.Compatible, Is.True);
            Assert.That(forward.CompatibleEnd, Is.EqualTo("CN"), "BsmI's own end");
            Assert.That(forward.End2, Is.EqualTo("GN"), "BsrI's own end");
            Assert.That(reversed.CompatibleEnd, Is.EqualTo("GN"));
            Assert.That(reversed.End2, Is.EqualTo("CN"));

            Assert.That(blunt.CompatibleEnd, Is.EqualTo("blunt"));
            Assert.That(blunt.End2, Is.EqualTo("blunt"), "blunt is shared by definition");

            Assert.That(incompatible.Compatible, Is.False);
            Assert.That(incompatible.CompatibleEnd, Is.Null);
            Assert.That(incompatible.End2, Is.Null);
        });
    }

    /// <summary>
    /// B08 H32: <c>restriction_map</c>'s <c>totalSites</c>/<c>uniqueCutters</c> count cleavages, so a
    /// Type IIS enzyme is not a unique cutter merely because it has one forward-strand site — the
    /// tool description used to say "forward-strand site count", which reads wrong for every Golden
    /// Gate enzyme.
    /// </summary>
    [Test]
    public void RestrictionMap_TypeIISEnzyme_CountsBothStrandsCleavages()
    {
        const string molecule = "AAAAAAAAAAGGTCTCAAAAAAAAAAGAGACCAAAAAAAAAA";
        var map = MolToolsTools.restriction_map(molecule, ["BsaI"]);

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2));
            Assert.That(map.UniqueCutters, Is.Empty);
            Assert.That(map.SitesByEnzyme["BsaI"], Is.EqualTo(new[] { 10, 26 }));
        });
    }

    /// <summary>
    /// B08 H13: <c>get_enzyme</c> resolves the REBASE hyphen spelling of CviKI-1 (NEB R0710), which
    /// used to answer <c>enzyme: null</c>, and <c>find_restriction_sites</c> used to error on it.
    /// Reference: Bio.Restriction.CviKI_1 — site RGCY, blunt, and
    /// <c>search(Seq("AGCTAGCT"), linear=True) == [3, 7]</c> (1-based).
    /// </summary>
    [TestCase("CviKI-1")]
    [TestCase("cviki-1")]
    [TestCase("CviKI_1")]
    public void GetEnzyme_RebaseHyphenSpelling_Resolves(string spelling)
    {
        var looked = MolToolsTools.get_enzyme(spelling, "rebase");
        var digest = MolToolsTools.restriction_digest("AGCTAGCT", [spelling]);

        Assert.Multiple(() =>
        {
            Assert.That(looked.Enzyme, Is.Not.Null);
            Assert.That(looked.Enzyme!.Name, Is.EqualTo("CviKI_1"));
            Assert.That(looked.Enzyme.RecognitionSequence, Is.EqualTo("RGCY"));
            Assert.That(looked.OverhangSequence, Is.Empty, "Biopython CviKI_1.ovhgseq == ''");
            Assert.That(digest.Fragments.Select(f => f.Length), Is.EqualTo(new[] { 2, 4, 2 }),
                "cuts at the 0-based bonds 2 and 6 (Biopython [3, 7] 1-based)");
        });
    }

    /// <summary>
    /// Mirrors <c>MolToolsTools.OtherEnd</c>: the partner enzyme's own overhang window, which the
    /// wrapper reports as <c>end2</c> so a caller sees both ends of a compatible pair (B08 H30).
    /// </summary>
    private static string? OtherEnd(string otherEnzymeName, string? compatibleEnd)
        => compatibleEnd == "blunt" ? "blunt" : RestrictionAnalyzer.GetOverhangSequence(otherEnzymeName);
}
