using static Seqeron.Genomics.Annotation.VariantAnnotator;
using VariantType = Seqeron.Genomics.Annotation.VariantAnnotator.VariantType;
using Variant = Seqeron.Genomics.Annotation.VariantAnnotator.Variant;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// VARIANT-ANNOT-001 mutation killers: exact-value tests for the variant-classification, normalization,
/// consequence-ranking, conservation-scoring, regulatory-overlap and VCF parse/format helpers that the
/// canonical FunctionalImpact fixture left uncovered. Conforms to Ensembl VEP / SO consequence terms,
/// PhyloP/PhastCons/GERP conservation definitions and the VCFv4.x INFO format.
/// </summary>
[TestFixture]
public class VariantAnnotatorMutationTests
{
    private const double Tol = 1e-9;

    #region ClassifyVariant / NormalizeVariant

    [Test]
    public void ClassifyVariant_AllTypes()
    {
        Assert.That(ClassifyVariant("A", "G"), Is.EqualTo(VariantType.SNV));
        Assert.That(ClassifyVariant("AT", "GC"), Is.EqualTo(VariantType.MNV));
        Assert.That(ClassifyVariant("A", "ATG"), Is.EqualTo(VariantType.Insertion)); // ref is prefix of alt
        Assert.That(ClassifyVariant("ATG", "A"), Is.EqualTo(VariantType.Deletion));  // alt is prefix of ref
        Assert.That(ClassifyVariant("ATG", "GC"), Is.EqualTo(VariantType.Indel));
        Assert.That(ClassifyVariant("", "G"), Is.EqualTo(VariantType.Complex));
    }

    [Test]
    public void NormalizeVariant_TrimsCommonPrefixAndSuffix()
    {
        // GCAT/GCGT → trim suffix T, trim prefix GC (position +2) → A/G SNV at 102.
        var v = NormalizeVariant("chr1", 100, "GCAT", "GCGT");
        Assert.That(v.Position, Is.EqualTo(102));
        Assert.That(v.Reference, Is.EqualTo("A"));
        Assert.That(v.Alternate, Is.EqualTo("G"));
        Assert.That(v.Type, Is.EqualTo(VariantType.SNV));
    }

    #endregion

    #region Consequence rank / impact

    [Test]
    public void GetConsequenceRank_MoreSevereHasLowerRank()
    {
        Assert.That(GetConsequenceRank(ConsequenceType.StopGained),
            Is.LessThan(GetConsequenceRank(ConsequenceType.MissenseVariant)));
        Assert.That(GetConsequenceRank(ConsequenceType.MissenseVariant),
            Is.LessThan(GetConsequenceRank(ConsequenceType.SynonymousVariant)));
    }

    [Test]
    public void GetImpactLevel_SoTermToImpact()
    {
        Assert.That(GetImpactLevel(ConsequenceType.StopGained), Is.EqualTo(ImpactLevel.High));
        Assert.That(GetImpactLevel(ConsequenceType.MissenseVariant), Is.EqualTo(ImpactLevel.Moderate));
        Assert.That(GetImpactLevel(ConsequenceType.SynonymousVariant), Is.EqualTo(ImpactLevel.Low));
        Assert.That(GetImpactLevel(ConsequenceType.IntronVariant), Is.EqualTo(ImpactLevel.Modifier));
    }

    #endregion

    #region CalculateConservation / FindConservedElements

    [Test]
    public void CalculateConservation_ExactMetrics()
    {
        var rows = CalculateConservation(new[]
        {
            ("chr1", 10, (IReadOnlyList<char>)new[] { 'A', 'A', 'A', 'A' }),
            ("chr1", 11, (IReadOnlyList<char>)new[] { 'A', 'T', 'A', 'T' }),
            ("chr1", 12, (IReadOnlyList<char>)Array.Empty<char>()),
        }).ToList();

        // Fully conserved: fraction 1 ⇒ phyloP (1−0.5)·12 = 6, phastCons 1, GERP clamped to 6.
        Assert.That(rows[0].PhyloP, Is.EqualTo(6.0).Within(Tol));
        Assert.That(rows[0].PhastCons, Is.EqualTo(1.0).Within(Tol));
        Assert.That(rows[0].Gerp, Is.EqualTo(6.0).Within(Tol));
        Assert.That(rows[0].ConservedSpeciesCount, Is.EqualTo(4));
        // Half conserved: phyloP 0, phastCons 0.5, GERP (2−1)/(4−1)·6 = 2.
        Assert.That(rows[1].PhyloP, Is.EqualTo(0.0).Within(Tol));
        Assert.That(rows[1].PhastCons, Is.EqualTo(0.5).Within(Tol));
        Assert.That(rows[1].Gerp, Is.EqualTo(2.0).Within(Tol));
        // Empty ⇒ all zero.
        Assert.That(rows[2].PhastCons, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void FindConservedElements_ContiguousHighScoreRun()
    {
        var scores = Enumerable.Range(1, 25)
            .Select(p => new ConservationScore("chr1", p, 6.0, 1.0, 6.0, 4));
        var elements = FindConservedElements(scores, threshold: 0.8, minLength: 20).ToList();

        Assert.That(elements, Has.Count.EqualTo(1));
        Assert.That(elements[0].Start, Is.EqualTo(1));
        Assert.That(elements[0].End, Is.EqualTo(25));
        Assert.That(elements[0].Score, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion

    #region Regulatory overlap / TF binding

    [Test]
    public void AnnotateRegulatoryElements_ReportsOverlappingRegionOnly()
    {
        var v = new Variant("chr1", 100, "A", "G", VariantType.SNV);
        var regions = new (string, int, int, string, string?, double?, IReadOnlyList<string>)[]
        {
            ("chr1", 90, 110, "Promoter", null, null, new[] { "SP1" }),   // overlaps
            ("chr1", 200, 300, "Enhancer", null, null, new[] { "GATA1" }), // does not overlap
            ("chr2", 95, 105, "Promoter", null, null, new[] { "SP1" }),    // wrong chromosome
        };
        var anns = AnnotateRegulatoryElements(v, regions).ToList();
        Assert.That(anns, Has.Count.EqualTo(1));
        Assert.That(anns[0].FeatureType, Is.EqualTo("Promoter"));
        Assert.That(anns[0].Start, Is.EqualTo(90));
    }

    [Test]
    public void PredictTfBindingChange_NonSnvYieldsNothing()
    {
        var indel = new Variant("chr1", 100, "AT", "A", VariantType.Deletion);
        var result = PredictTfBindingChange(
            indel,
            new[] { ("CTCF", "CCGCGNGGNGGCAG", 5.0) },
            "ACGTACGTACGTACGTACGTACGTACGT");
        Assert.That(result.ToList(), Is.Empty);
    }

    #endregion

    #region VCF parse / format

    [Test]
    public void ParseVcfVariant_ClassifiesAndCarriesIdQuality()
    {
        var v = ParseVcfVariant("chr1", 100, "rs123", "A", "G", quality: 30.0);
        Assert.That(v.Type, Is.EqualTo(VariantType.SNV));
        Assert.That(v.Id, Is.EqualTo("rs123"));
        Assert.That(v.Quality, Is.EqualTo(30.0).Within(Tol));
    }

    [Test]
    public void FormatAsVcfInfo_ExactInfoString()
    {
        var ann = new VariantAnnotation(
            new Variant("chr1", 100, "A", "G", VariantType.SNV),
            TranscriptId: "ENST1", GeneId: "ENSG1", GeneName: "BRCA1",
            Consequence: ConsequenceType.MissenseVariant, Impact: ImpactLevel.Moderate,
            CodonChange: "c.1A>G", AminoAcidChange: "p.A1V",
            ProteinPosition: 1, CdsPosition: 1,
            SiftScore: 0.05, PolyphenScore: 0.9, CaddScore: null,
            ExistingVariation: null, PopulationFrequencies: null);

        Assert.That(FormatAsVcfInfo(ann), Is.EqualTo(
            "GENE=BRCA1;TRANSCRIPT=ENST1;CONSEQUENCE=MissenseVariant;IMPACT=Moderate;" +
            "HGVSP=p.A1V;HGVSC=c.1A>G;SIFT=0.050;POLYPHEN=0.900"));
    }

    #endregion
}
