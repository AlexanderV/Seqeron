// DISORDER-PROPENSITY-001 — Disorder Propensity (TOP-IDP scale & Dunker classification)
// Evidence: docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md
// TestSpec: tests/TestSpecs/DISORDER-PROPENSITY-001.md
// Source: Campen A et al. (2008) Protein Pept Lett 15(9):956-963 (PMID 18991772, PMC2676888),
//         Table 2 — exact per-residue TOP-IDP values;
//         Dunker AK et al. (2001) J Mol Graph Model 19(1):26-59 (PMID 11381529),
//         order/disorder/ambiguous amino-acid classification.
//
// Expected values are copied from the published sources (TOP-IDP Table 2 and the Dunker
// classification sets), independently of the implementation's current output.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class DisorderPredictor_GetDisorderPropensity_Tests
{
    // TOP-IDP scale — Campen et al. (2008) Table 2 (PMC2676888), verbatim.
    private static readonly System.Collections.Generic.Dictionary<char, double> TopIdp = new()
    {
        ['A'] = 0.060, ['R'] = 0.180, ['N'] = 0.007, ['D'] = 0.192, ['C'] = 0.020,
        ['Q'] = 0.318, ['E'] = 0.736, ['G'] = 0.166, ['H'] = 0.303, ['I'] = -0.486,
        ['L'] = -0.326, ['K'] = 0.586, ['M'] = -0.397, ['F'] = -0.697, ['P'] = 0.987,
        ['S'] = 0.341, ['T'] = 0.059, ['W'] = -0.884, ['Y'] = -0.510, ['V'] = -0.121
    };

    // Dunker et al. (2001) classification sets, verbatim.
    private static readonly char[] DisorderPromoting = { 'A', 'E', 'G', 'K', 'P', 'Q', 'R', 'S' };
    private static readonly char[] OrderPromoting = { 'C', 'F', 'I', 'L', 'N', 'V', 'W', 'Y' };
    private static readonly char[] Ambiguous = { 'D', 'H', 'M', 'T' };

    #region GetDisorderPropensity — MUST

    // M1 — All 20 residues return their exact Campen (2008) Table 2 value.
    [Test]
    public void GetDisorderPropensity_AllTwentyAminoAcids_MatchScale()
    {
        Assert.Multiple(() =>
        {
            foreach (var (aa, expected) in TopIdp)
            {
                Assert.That(DisorderPredictor.GetDisorderPropensity(aa), Is.EqualTo(expected).Within(1e-10),
                    $"TOP-IDP({aa}) must equal {expected} — Campen et al. (2008) Table 2");
            }
        });
    }

    // M2 — Scale extrema: W is the global minimum, P the global maximum (Campen 2008 Table 2).
    [Test]
    public void GetDisorderPropensity_ExtremeAnchors_MinAtTrpMaxAtPro()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.GetDisorderPropensity('W'), Is.EqualTo(-0.884).Within(1e-10),
                "W is the most order-promoting residue (scale minimum) — Campen et al. (2008)");
            Assert.That(DisorderPredictor.GetDisorderPropensity('P'), Is.EqualTo(0.987).Within(1e-10),
                "P is the most disorder-promoting residue (scale maximum) — Campen et al. (2008)");

            double minValue = TopIdp.Values.Min();
            double maxValue = TopIdp.Values.Max();
            Assert.That(DisorderPredictor.GetDisorderPropensity('W'), Is.EqualTo(minValue).Within(1e-10),
                "W must be the minimum over all 20 residues");
            Assert.That(DisorderPredictor.GetDisorderPropensity('P'), Is.EqualTo(maxValue).Within(1e-10),
                "P must be the maximum over all 20 residues");
        });
    }

    #endregion

    #region IsDisorderPromoting — MUST

    // M3 — Disorder-promoting residues → true (Dunker 2001).
    [Test]
    [TestCase('A', TestName = "IsDisorderPromoting_Ala_True")]
    [TestCase('R', TestName = "IsDisorderPromoting_Arg_True")]
    [TestCase('G', TestName = "IsDisorderPromoting_Gly_True")]
    [TestCase('Q', TestName = "IsDisorderPromoting_Gln_True")]
    [TestCase('S', TestName = "IsDisorderPromoting_Ser_True")]
    [TestCase('P', TestName = "IsDisorderPromoting_Pro_True")]
    [TestCase('E', TestName = "IsDisorderPromoting_Glu_True")]
    [TestCase('K', TestName = "IsDisorderPromoting_Lys_True")]
    public void IsDisorderPromoting_DisorderPromotingResidues_ReturnsTrue(char aa)
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.True,
            $"'{aa}' must be disorder-promoting — Dunker et al. (2001)");
    }

    // M4 — Order-promoting residues → false (Dunker 2001).
    [Test]
    [TestCase('W', TestName = "IsDisorderPromoting_Trp_False")]
    [TestCase('C', TestName = "IsDisorderPromoting_Cys_False")]
    [TestCase('F', TestName = "IsDisorderPromoting_Phe_False")]
    [TestCase('I', TestName = "IsDisorderPromoting_Ile_False")]
    [TestCase('Y', TestName = "IsDisorderPromoting_Tyr_False")]
    [TestCase('V', TestName = "IsDisorderPromoting_Val_False")]
    [TestCase('L', TestName = "IsDisorderPromoting_Leu_False")]
    [TestCase('N', TestName = "IsDisorderPromoting_Asn_False")]
    public void IsDisorderPromoting_OrderPromotingResidues_ReturnsFalse(char aa)
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.False,
            $"'{aa}' is order-promoting and must not be classified disorder-promoting — Dunker et al. (2001)");
    }

    // M5 — Ambiguous residues → false (Dunker 2001: not in the disorder-promoting set).
    [Test]
    [TestCase('H', TestName = "IsDisorderPromoting_His_False")]
    [TestCase('M', TestName = "IsDisorderPromoting_Met_False")]
    [TestCase('T', TestName = "IsDisorderPromoting_Thr_False")]
    [TestCase('D', TestName = "IsDisorderPromoting_Asp_False")]
    public void IsDisorderPromoting_AmbiguousResidues_ReturnsFalse(char aa)
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.False,
            $"'{aa}' is ambiguous and must not be classified disorder-promoting — Dunker et al. (2001)");
    }

    // M10 — INV-3: IsDisorderPromoting(c) ⇔ c ∈ DisorderPromotingAminoAcids, for all 20 residues.
    [Test]
    public void IsDisorderPromoting_MatchesProperty_AllStandardResidues()
    {
        var set = DisorderPredictor.DisorderPromotingAminoAcids;
        Assert.Multiple(() =>
        {
            foreach (char aa in TopIdp.Keys)
            {
                bool expected = set.Contains(aa);
                Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.EqualTo(expected),
                    $"IsDisorderPromoting('{aa}') must equal membership in DisorderPromotingAminoAcids — INV-3");
            }
        });
    }

    #endregion

    #region Classification Properties — MUST

    // M6 — DisorderPromotingAminoAcids = {A,E,G,K,P,Q,R,S}, exactly 8 (Dunker 2001).
    [Test]
    public void DisorderPromotingAminoAcids_EqualsDunkerSet()
    {
        var actual = DisorderPredictor.DisorderPromotingAminoAcids;
        Assert.Multiple(() =>
        {
            Assert.That(actual.Count, Is.EqualTo(8), "Disorder-promoting set has 8 members — Dunker et al. (2001)");
            Assert.That(actual, Is.EquivalentTo(DisorderPromoting),
                "DisorderPromotingAminoAcids must equal {A,E,G,K,P,Q,R,S} — Dunker et al. (2001)");
        });
    }

    // M7 — OrderPromotingAminoAcids = {C,F,I,L,N,V,W,Y}, exactly 8 (Dunker 2001).
    [Test]
    public void OrderPromotingAminoAcids_EqualsDunkerSet()
    {
        var actual = DisorderPredictor.OrderPromotingAminoAcids;
        Assert.Multiple(() =>
        {
            Assert.That(actual.Count, Is.EqualTo(8), "Order-promoting set has 8 members — Dunker et al. (2001)");
            Assert.That(actual, Is.EquivalentTo(OrderPromoting),
                "OrderPromotingAminoAcids must equal {C,F,I,L,N,V,W,Y} — Dunker et al. (2001)");
        });
    }

    // M8 — AmbiguousAminoAcids = {D,H,M,T}, exactly 4 (Dunker 2001).
    [Test]
    public void AmbiguousAminoAcids_EqualsDunkerSet()
    {
        var actual = DisorderPredictor.AmbiguousAminoAcids;
        Assert.Multiple(() =>
        {
            Assert.That(actual.Count, Is.EqualTo(4), "Ambiguous set has 4 members — Dunker et al. (2001)");
            Assert.That(actual, Is.EquivalentTo(Ambiguous),
                "AmbiguousAminoAcids must equal {D,H,M,T} — Dunker et al. (2001)");
        });
    }

    // M9 — INV-4: the three sets are pairwise disjoint and cover all 20 standard residues (8+8+4).
    [Test]
    public void ClassificationSets_AreDisjointAndCoverAll20()
    {
        var disorder = DisorderPredictor.DisorderPromotingAminoAcids;
        var order = DisorderPredictor.OrderPromotingAminoAcids;
        var ambiguous = DisorderPredictor.AmbiguousAminoAcids;
        var union = disorder.Concat(order).Concat(ambiguous).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(union.Count, Is.EqualTo(20), "Three classes must total 20 residues — Dunker et al. (2001)");
            Assert.That(union.Distinct().Count(), Is.EqualTo(20), "Three classes must be pairwise disjoint");
            Assert.That(union, Is.EquivalentTo(TopIdp.Keys),
                "Union of the three classes must be exactly the 20 standard amino acids");
        });
    }

    // C1 — Property lists are returned in ascending (sorted) order.
    [Test]
    public void Properties_ReturnSortedAscending()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.DisorderPromotingAminoAcids, Is.Ordered.Ascending,
                "DisorderPromotingAminoAcids must be sorted ascending");
            Assert.That(DisorderPredictor.OrderPromotingAminoAcids, Is.Ordered.Ascending,
                "OrderPromotingAminoAcids must be sorted ascending");
            Assert.That(DisorderPredictor.AmbiguousAminoAcids, Is.Ordered.Ascending,
                "AmbiguousAminoAcids must be sorted ascending");
        });
    }

    #endregion

    #region Edge Cases — SHOULD

    // S1 — Unknown residue → 0.0 (implementation contract; not source-defined).
    [Test]
    public void GetDisorderPropensity_UnknownResidue_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.GetDisorderPropensity('X'), Is.EqualTo(0.0),
                "Unknown residue 'X' must return 0.0 (lookup default)");
            Assert.That(DisorderPredictor.GetDisorderPropensity('Z'), Is.EqualTo(0.0),
                "Unknown residue 'Z' must return 0.0 (lookup default)");
            Assert.That(DisorderPredictor.GetDisorderPropensity('B'), Is.EqualTo(0.0),
                "Unknown residue 'B' must return 0.0 (lookup default)");
            Assert.That(DisorderPredictor.GetDisorderPropensity('*'), Is.EqualTo(0.0),
                "Non-letter '*' must return 0.0 (lookup default)");
        });
    }

    // S2 — Lowercase input equals uppercase (INV-5 case-insensitivity).
    [Test]
    public void GetDisorderPropensity_LowercaseInput_SameAsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.GetDisorderPropensity('p'), Is.EqualTo(0.987).Within(1e-10),
                "Lowercase 'p' must return the same value as 'P'");
            Assert.That(DisorderPredictor.GetDisorderPropensity('w'), Is.EqualTo(-0.884).Within(1e-10),
                "Lowercase 'w' must return the same value as 'W'");
            Assert.That(DisorderPredictor.GetDisorderPropensity('e'), Is.EqualTo(0.736).Within(1e-10),
                "Lowercase 'e' must return the same value as 'E'");
        });
    }

    // S3 — IsDisorderPromoting is case-insensitive.
    [Test]
    public void IsDisorderPromoting_LowercaseInput_SameAsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.IsDisorderPromoting('p'), Is.True,
                "Lowercase 'p' must be disorder-promoting like 'P'");
            Assert.That(DisorderPredictor.IsDisorderPromoting('w'), Is.False,
                "Lowercase 'w' must be order-promoting like 'W'");
        });
    }

    #endregion
}
