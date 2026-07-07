// ONCO-FUSION-002 — Known Fusion Database Lookup
// Evidence: docs/Evidence/ONCO-FUSION-002-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-FUSION-002.md
// Source: Bruford et al. (2021). HGNC recommendations for the designation of gene fusions.
//         Leukemia 35(11):3040-3043. https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
//         ("a double colon (::)"; "the 5′ partner gene should always be listed first"; "BCR::ABL1").

using Call = Seqeron.Genomics.Oncology.OncologyAnalyzer.FusionCall;
using Frame = Seqeron.Genomics.Oncology.OncologyAnalyzer.FusionReadingFrame;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_MatchKnownFusions_Tests
{
    private static Call MakeCall(string gene5p, string gene3p) =>
        new(gene5p, gene3p, JunctionReads: 5, DiscordantMates: 4, TotalSupport: 9, Frame.InFrame);

    #region GetFusionAnnotation — HGNC designation format

    // M1 — Bruford et al. (2021): 5' gene first, double-colon separator; BCR is 5', ABL1 is 3' → "BCR::ABL1".
    [Test]
    public void GetFusionAnnotation_BcrAbl1_FivePrimeFirstDoubleColon()
    {
        string designation = OncologyAnalyzer.GetFusionAnnotation("BCR", "ABL1");

        Assert.That(designation, Is.EqualTo("BCR::ABL1"),
            "HGNC: 5' partner first, joined by '::' — the worked example BCR::ABL1 (Bruford et al. 2021).");
    }

    // M2 — Direction matters: 5' partner is fixed by role, not sorted, so swapping yields a different designation.
    [Test]
    public void GetFusionAnnotation_Reciprocal_IsDirectional()
    {
        string forward = OncologyAnalyzer.GetFusionAnnotation("BCR", "ABL1");
        string reciprocal = OncologyAnalyzer.GetFusionAnnotation("ABL1", "BCR");

        Assert.Multiple(() =>
        {
            Assert.That(reciprocal, Is.EqualTo("ABL1::BCR"),
                "Reciprocal designation lists ABL1 as the 5' partner (Bruford et al. 2021).");
            Assert.That(reciprocal, Is.Not.EqualTo(forward),
                "A::B and B::A are different fusions because the 5' partner is always listed first.");
        });
    }

    // S2 — designation is a verbatim concatenation, so it preserves the input case of the symbols.
    [Test]
    public void GetFusionAnnotation_PreservesInputCase()
    {
        string designation = OncologyAnalyzer.GetFusionAnnotation("bcr", "abl1");

        Assert.That(designation, Is.EqualTo("bcr::abl1"),
            "The designation is verbatim gene5p + '::' + gene3p; case folding belongs to matching, not formatting.");
    }

    // M6 — null 5' partner is invalid input.
    [Test]
    public void GetFusionAnnotation_NullFivePrime_Throws()
    {
        Assert.That(() => OncologyAnalyzer.GetFusionAnnotation(null!, "ABL1"),
            NUnit.Framework.Throws.ArgumentException, "A null 5' partner symbol cannot form a designation.");
    }

    // M7 — empty 3' partner is invalid input.
    [Test]
    public void GetFusionAnnotation_EmptyThreePrime_Throws()
    {
        Assert.That(() => OncologyAnalyzer.GetFusionAnnotation("BCR", ""),
            NUnit.Framework.Throws.ArgumentException, "An empty 3' partner symbol cannot form a designation.");
    }

    // M7b — empty 5' partner is invalid input (the other half of the 5' validation branch).
    [Test]
    public void GetFusionAnnotation_EmptyFivePrime_Throws()
    {
        Assert.That(() => OncologyAnalyzer.GetFusionAnnotation("", "ABL1"),
            NUnit.Framework.Throws.ArgumentException, "An empty 5' partner symbol cannot form a designation.");
    }

    // M6b — null 3' partner is invalid input (the other half of the 3' validation branch).
    [Test]
    public void GetFusionAnnotation_NullThreePrime_Throws()
    {
        Assert.That(() => OncologyAnalyzer.GetFusionAnnotation("BCR", null!),
            NUnit.Framework.Throws.ArgumentException, "A null 3' partner symbol cannot form a designation.");
    }

    // M7c — whitespace-only symbols are rejected (contract: non-empty, non-whitespace).
    [Test]
    public void GetFusionAnnotation_WhitespaceSymbols_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.GetFusionAnnotation("   ", "ABL1"),
                NUnit.Framework.Throws.ArgumentException, "A whitespace-only 5' partner is not a valid symbol.");
            Assert.That(() => OncologyAnalyzer.GetFusionAnnotation("BCR", "\t"),
                NUnit.Framework.Throws.ArgumentException, "A whitespace-only 3' partner is not a valid symbol.");
        });
    }

    #endregion

    #region MatchKnownFusions — directional lookup against caller-supplied set

    // M3 — designation present in the caller's set → matched, annotation returned, keyed by 5'::3'.
    [Test]
    public void MatchKnownFusions_DesignationPresent_ReturnsAnnotation()
    {
        var known = new Dictionary<string, string>
        {
            ["EML4::ALK"] = "NSCLC driver",
        };

        OncologyAnalyzer.KnownFusionMatch match =
            OncologyAnalyzer.MatchKnownFusions(MakeCall("EML4", "ALK"), known);

        Assert.Multiple(() =>
        {
            Assert.That(match.IsKnown, Is.True, "EML4::ALK is present in the supplied set.");
            Assert.That(match.Designation, Is.EqualTo("EML4::ALK"), "Key built 5' first per HGNC.");
            Assert.That(match.Annotation, Is.EqualTo("NSCLC driver"), "The caller-supplied annotation is returned.");
        });
    }

    // M4 — only the reciprocal key is present; the lookup is directional → no match.
    [Test]
    public void MatchKnownFusions_OnlyReciprocalPresent_NoMatch()
    {
        var known = new Dictionary<string, string>
        {
            ["ALK::EML4"] = "reciprocal (different fusion)",
        };

        OncologyAnalyzer.KnownFusionMatch match =
            OncologyAnalyzer.MatchKnownFusions(MakeCall("EML4", "ALK"), known);

        Assert.Multiple(() =>
        {
            Assert.That(match.IsKnown, Is.False,
                "EML4::ALK ≠ ALK::EML4: the 5'-first rule makes the designation directional (Bruford et al. 2021).");
            Assert.That(match.Annotation, Is.Null, "An unmatched lookup carries no annotation.");
            Assert.That(match.Designation, Is.EqualTo("EML4::ALK"), "Designation reflects the query, not the set.");
        });
    }

    // M5 — query absent from the set → no match.
    [Test]
    public void MatchKnownFusions_DesignationAbsent_NoMatch()
    {
        var known = new Dictionary<string, string>
        {
            ["BCR::ABL1"] = "CML driver",
        };

        OncologyAnalyzer.KnownFusionMatch match =
            OncologyAnalyzer.MatchKnownFusions(MakeCall("EML4", "ALK"), known);

        Assert.That(match.IsKnown, Is.False, "EML4::ALK is not in a set containing only BCR::ABL1.");
    }

    // M8 — null known-fusion map is invalid input.
    [Test]
    public void MatchKnownFusions_NullKnownSet_Throws()
    {
        Assert.That(() => OncologyAnalyzer.MatchKnownFusions(MakeCall("EML4", "ALK"), null!),
            NUnit.Framework.Throws.ArgumentNullException, "The known-fusion set must be supplied.");
    }

    // M8b — a fusion with an invalid (empty) partner symbol propagates the designation-format contract.
    [Test]
    public void MatchKnownFusions_EmptyPartner_Throws()
    {
        var known = new Dictionary<string, string> { ["EML4::ALK"] = "NSCLC driver" };

        Assert.That(() => OncologyAnalyzer.MatchKnownFusions(MakeCall("", "ALK"), known),
            NUnit.Framework.Throws.ArgumentException,
            "An empty partner symbol cannot form the 5'::3' lookup key (delegates to GetFusionAnnotation).");
    }

    // S1 — case-insensitive symbol matching: lowercased query matches a stored uppercase designation.
    [Test]
    public void MatchKnownFusions_CaseInsensitiveSymbols_Match()
    {
        var known = new Dictionary<string, string>
        {
            ["EML4::ALK"] = "NSCLC driver",
        };

        OncologyAnalyzer.KnownFusionMatch match =
            OncologyAnalyzer.MatchKnownFusions(MakeCall("eml4", "alk"), known);

        Assert.Multiple(() =>
        {
            Assert.That(match.IsKnown, Is.True, "Symbol comparison is case-insensitive (INV-04).");
            Assert.That(match.Annotation, Is.EqualTo("NSCLC driver"), "The stored annotation is returned for the case-folded match.");
        });
    }

    // C1 — integration: a FusionCall produced by DetectFusions matches by its designation.
    [Test]
    public void MatchKnownFusions_FromDetectFusions_Match()
    {
        var candidates = new[] { new OncologyAnalyzer.FusionCandidate("EML4", "ALK", 3, 2, 4) };
        IReadOnlyList<Call> calls = OncologyAnalyzer.DetectFusions(candidates);
        var known = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase)
        {
            ["EML4::ALK"] = "NSCLC driver",
        };

        OncologyAnalyzer.KnownFusionMatch match = OncologyAnalyzer.MatchKnownFusions(calls.Single(), known);

        Assert.Multiple(() =>
        {
            Assert.That(match.Designation, Is.EqualTo("EML4::ALK"), "Designation derived from the detected fusion's partners.");
            Assert.That(match.IsKnown, Is.True, "The detected EML4::ALK fusion is present in the known set.");
        });
    }

    #endregion
}
