// PROTMOTIF-COMMON-001 — Common Motif Finding
// Evidence: docs/Evidence/PROTMOTIF-COMMON-001-Evidence.md
// TestSpec: tests/TestSpecs/PROTMOTIF-COMMON-001.md
// Source: ExPASy PROSITE (PS00001/PS00005/PS00006/PS00016/PS00017) + ScanProsite syntax docs;
//         Sigrist CJA et al. (2013) Nucleic Acids Res 41(D1):D344-D347.
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// PROTMOTIF-COMMON-001: Common Motif Finding — canonical test file.
/// Validates <see cref="ProteinMotifFinder.FindCommonMotifs"/> aggregation semantics
/// (whole-dictionary scan, multi-pattern, multi-occurrence, identity propagation, invariants)
/// against official PROSITE entries (https://prosite.expasy.org/).
/// </summary>
[TestFixture]
public class ProteinMotifFinder_FindCommonMotifs_Tests
{
    #region FindCommonMotifs — single-pattern exact hits (MUST)

    // M1 — PS00001 N-{P}-[ST]-{P}: only window AAAA[NFTA]AAA satisfies all four elements.
    // https://prosite.expasy.org/PS00001
    [Test]
    public void FindCommonMotifs_NGlycosylationWindow_FindsExactSite()
    {
        const string protein = "AAAANFTAAAA";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "ASN_GLYCOSYLATION").ToList();

        Assert.That(hits, Has.Count.EqualTo(1),
            "PS00001 N-{P}-[ST]-{P} matches exactly the NFTA window once");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(4), "N of NFTA is at 0-based index 4");
            Assert.That(hits[0].End, Is.EqualTo(7), "4-element pattern spans indices 4..7");
            Assert.That(hits[0].Sequence, Is.EqualTo("NFTA"),
                "N(Asn), F(non-Pro), T(Ser/Thr), A(non-Pro) per PS00001");
        });
    }

    // M3 — PS00005 [ST]-x-[RK]: S(5)-A(6)-R(7). https://prosite.expasy.org/PS00005
    [Test]
    public void FindCommonMotifs_PkcWindow_FindsExactSite()
    {
        const string protein = "AAAAASARKAAA";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "PKC_PHOSPHO_SITE").ToList();

        Assert.That(hits, Has.Count.EqualTo(1), "PS00005 [ST]-x-[RK] matches SAR once");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(5), "S of SAR is at 0-based index 5");
            Assert.That(hits[0].End, Is.EqualTo(7), "3-element pattern spans 5..7");
            Assert.That(hits[0].Sequence, Is.EqualTo("SAR"), "[ST]=S, x=A, [RK]=R");
        });
    }

    // M4 — PS00006 [ST]-x(2)-[DE]: two windows SAAE(4..7) and SDED(9..12).
    // https://prosite.expasy.org/PS00006
    [Test]
    public void FindCommonMotifs_Ck2TwoWindows_FindsBothSites()
    {
        const string protein = "AAAASAAEASDEDAAA";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "CK2_PHOSPHO_SITE")
            .OrderBy(m => m.Start).ToList();

        Assert.That(hits, Has.Count.EqualTo(2),
            "PS00006 [ST]-x(2)-[DE] matches SAAE and SDED");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(4), "first S at index 4");
            Assert.That(hits[0].End, Is.EqualTo(7), "first window 4..7");
            Assert.That(hits[0].Sequence, Is.EqualTo("SAAE"), "[ST]=S, x(2)=AA, [DE]=E");
            Assert.That(hits[1].Start, Is.EqualTo(9), "second S at index 9");
            Assert.That(hits[1].End, Is.EqualTo(12), "second window 9..12");
            Assert.That(hits[1].Sequence, Is.EqualTo("SDED"), "[ST]=S, x(2)=DE, [DE]=D");
        });
    }

    // M5 — PS00017 [AG]-x(4)-G-K-[ST]: G(5)..S(12). https://prosite.expasy.org/PS00017
    [Test]
    public void FindCommonMotifs_PLoopWindow_FindsExactSite()
    {
        const string protein = "AAAAAGXXXXGKSAAAA";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "ATP_GTP_A").ToList();

        Assert.That(hits, Has.Count.EqualTo(1), "PS00017 P-loop matches once");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(5), "[AG]=G at index 5");
            Assert.That(hits[0].End, Is.EqualTo(12), "8-residue pattern spans 5..12");
            Assert.That(hits[0].Sequence, Is.EqualTo("GXXXXGKS"),
                "[AG]=G, x(4)=XXXX, G, K, [ST]=S");
        });
    }

    // M6 — PS00016 R-G-D: literal triplet at 2..4. https://prosite.expasy.org/PS00016
    [Test]
    public void FindCommonMotifs_RgdWindow_FindsExactSite()
    {
        const string protein = "AARGDKK";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "RGD").ToList();

        Assert.That(hits, Has.Count.EqualTo(1), "PS00016 R-G-D matches once");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(2), "R at index 2");
            Assert.That(hits[0].End, Is.EqualTo(4), "RGD spans 2..4");
            Assert.That(hits[0].Sequence, Is.EqualTo("RGD"), "R-G-D literal");
        });
    }

    #endregion

    #region FindCommonMotifs — exclusion / negation (MUST)

    // M2 — PS00001 {P}: Pro at the excluded second position rejects N-P-[ST] windows.
    // https://prosite.expasy.org/PS00001
    [Test]
    public void FindCommonMotifs_ProlineAtExcludedPosition_RejectsSite()
    {
        const string protein = "AAAANPSAAAAANPTAAA";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "ASN_GLYCOSYLATION").ToList();

        Assert.That(hits, Has.Count.EqualTo(0),
            "N-P-S and N-P-T are NOT N-glycosylation sites: {P} forbids Pro at position 2");
    }

    #endregion

    #region FindCommonMotifs — whole-dictionary aggregation (MUST)

    // M7 — whole-dictionary scan must surface two DIFFERENT pattern types from one sequence.
    // RGDNFTA contains exactly two PROSITE windows over the whole CommonMotifs library:
    // RGD (PS00016) at 0..2 and N-glycosylation (PS00001) at 3..6 — independently confirmed
    // by scanning every library pattern (no other entry matches this 7-mer).
    [Test]
    public void FindCommonMotifs_TwoDistinctPatterns_ReturnsBothTypes()
    {
        const string protein = "RGDNFTA";

        var hits = FindCommonMotifs(protein)
            .OrderBy(m => m.Start).ThenBy(m => m.MotifName).ToList();

        // Exact, not "contains": the entire library yields precisely these two hits.
        Assert.That(hits, Has.Count.EqualTo(2),
            "RGDNFTA matches exactly two library patterns across the whole CommonMotifs scan");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].MotifName, Is.EqualTo("RGD"), "first hit is RGD (PS00016)");
            Assert.That(hits[0].Pattern, Is.EqualTo("PS00016"));
            Assert.That(hits[0].Start, Is.EqualTo(0));
            Assert.That(hits[0].End, Is.EqualTo(2));
            Assert.That(hits[0].Sequence, Is.EqualTo("RGD"));

            Assert.That(hits[1].MotifName, Is.EqualTo("ASN_GLYCOSYLATION"),
                "second hit is N-glycosylation (PS00001) — proves the scan covers the whole library");
            Assert.That(hits[1].Pattern, Is.EqualTo("PS00001"));
            Assert.That(hits[1].Start, Is.EqualTo(3));
            Assert.That(hits[1].End, Is.EqualTo(6));
            Assert.That(hits[1].Sequence, Is.EqualTo("NFTA"));
        });
    }

    // M8 — two occurrences of one pattern must both be reported. RGDRGD -> (0..2),(3..5).
    [Test]
    public void FindCommonMotifs_TwoRgdOccurrences_ReturnsBoth()
    {
        const string protein = "RGDRGD";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "RGD")
            .OrderBy(m => m.Start).ToList();

        Assert.That(hits, Has.Count.EqualTo(2),
            "both RGD occurrences must be reported (PROSITE default reports all occurrences)");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(0), "first RGD at 0");
            Assert.That(hits[0].End, Is.EqualTo(2), "first RGD 0..2");
            Assert.That(hits[1].Start, Is.EqualTo(3), "second RGD at 3");
            Assert.That(hits[1].End, Is.EqualTo(5), "second RGD 3..5");
        });
    }

    // M8b — INV-03: genuinely OVERLAPPING occurrences of one pattern are both reported.
    // PS00005 [ST]-x-[RK] on "STRK": window S(0)-T(1)-R(2) = "STR" and T(1)-R(2)-K(3) = "TRK"
    // share residues 1..2 yet neither is contained in the other, so ScanProsite default
    // ("greedy, overlaps, no includes") reports BOTH. Independently hand-verified with the
    // lookahead translation (?=([ST].[RK])) over "STRK". Source: ScanProsite overlap doc
    // (https://prosite.expasy.org/scanprosite/scanprosite_doc.html).
    [Test]
    public void FindCommonMotifs_OverlappingPkcOccurrences_ReturnsBoth()
    {
        const string protein = "STRK";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "PKC_PHOSPHO_SITE")
            .OrderBy(m => m.Start).ToList();

        Assert.That(hits, Has.Count.EqualTo(2),
            "two overlapping [ST]-x-[RK] windows must both be reported (overlaps, no includes)");
        Assert.Multiple(() =>
        {
            Assert.That(hits[0].Start, Is.EqualTo(0), "first window S(0)..R(2)");
            Assert.That(hits[0].End, Is.EqualTo(2));
            Assert.That(hits[0].Sequence, Is.EqualTo("STR"));
            Assert.That(hits[1].Start, Is.EqualTo(1), "second window T(1)..K(3) overlaps the first");
            Assert.That(hits[1].End, Is.EqualTo(3));
            Assert.That(hits[1].Sequence, Is.EqualTo("TRK"));
        });
    }

    #endregion

    #region FindCommonMotifs — invariants (SHOULD)

    // S1 — INV-1: matched Sequence equals the substring at [Start..End].
    [Test]
    public void FindCommonMotifs_AllMatches_SatisfySubstringInvariant()
    {
        const string protein = "MNFTAKSARKGAAAAGKSRGDPPAY";

        var matches = FindCommonMotifs(protein).ToList();

        Assert.That(matches, Is.Not.Empty, "test sequence is designed to produce matches");
        Assert.Multiple(() =>
        {
            foreach (var m in matches)
            {
                Assert.That(m.Start, Is.InRange(0, protein.Length - 1),
                    $"{m.MotifName} Start within bounds");
                Assert.That(m.End, Is.InRange(m.Start, protein.Length - 1),
                    $"{m.MotifName} End >= Start and within bounds");
                Assert.That(m.Sequence, Is.EqualTo(protein.Substring(m.Start, m.End - m.Start + 1)),
                    $"{m.MotifName} Sequence must equal substring at its reported coordinates");
            }
        });
    }

    // S2 — INV-3: a match carries the exact Name and Accession of its CommonMotifs entry.
    [Test]
    public void FindCommonMotifs_RgdMatch_CarriesDictionaryIdentity()
    {
        const string protein = "AARGDKK";

        var hit = FindCommonMotifs(protein).Single(m => m.MotifName == "RGD");

        Assert.Multiple(() =>
        {
            Assert.That(hit.MotifName, Is.EqualTo("RGD"),
                "MotifName must be the PROSITE entry name from CommonMotifs[PS00016]");
            Assert.That(hit.Pattern, Is.EqualTo("PS00016"),
                "Pattern must be the PROSITE accession from CommonMotifs[PS00016]");
        });
    }

    // S3 — negative control: a sequence with no motif returns nothing.
    [Test]
    public void FindCommonMotifs_SequenceWithNoMotif_ReturnsEmpty()
    {
        const string protein = "AAAAAAAAAA";

        var matches = FindCommonMotifs(protein).ToList();

        Assert.That(matches, Is.Empty,
            "a poly-Ala sequence satisfies none of the library patterns");
    }

    #endregion

    #region FindCommonMotifs — edge cases (MUST) and extras (COULD)

    // M9 — null input.
    [Test]
    public void FindCommonMotifs_NullSequence_ReturnsEmpty()
    {
        var matches = FindCommonMotifs(null!).ToList();
        Assert.That(matches, Is.Empty, "null input yields no matches");
    }

    // M10 — empty input.
    [Test]
    public void FindCommonMotifs_EmptySequence_ReturnsEmpty()
    {
        var matches = FindCommonMotifs(string.Empty).ToList();
        Assert.That(matches, Is.Empty, "empty input yields no matches");
    }

    // C1 — INV-4: determinism.
    [Test]
    public void FindCommonMotifs_SameInput_IsDeterministic()
    {
        const string protein = "MNFTAKSARKGAAAAGKSRGD";

        var first = FindCommonMotifs(protein).ToList();
        var second = FindCommonMotifs(protein).ToList();

        Assert.That(
            second.Select(m => (m.MotifName, m.Start, m.End)),
            Is.EqualTo(first.Select(m => (m.MotifName, m.Start, m.End))),
            "two scans of the same sequence must produce identical ordered results");
    }

    // C2 — case-insensitivity: lowercase input is upper-cased before matching.
    [Test]
    public void FindCommonMotifs_LowercaseInput_FindsSite()
    {
        const string protein = "aargdkk";

        var hits = FindCommonMotifs(protein)
            .Where(m => m.MotifName == "RGD").ToList();

        Assert.That(hits, Has.Count.EqualTo(1),
            "lowercase rgd must still match RGD (input is upper-cased)");
        Assert.That(hits[0].Sequence, Is.EqualTo("RGD"),
            "matched sequence is reported upper-cased");
    }

    #endregion
}
