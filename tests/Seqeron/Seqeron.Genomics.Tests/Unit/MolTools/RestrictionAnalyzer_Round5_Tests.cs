namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 5 on the restriction slice (B08 A5R-2, A5R-3, A5R-5 .. A5R-8 and the surviving mutants A5R-T1 .. T8).
/// Each test failed on the code the auditor read (commit 644056cd) or, for the mutation rows, on the mutant the auditor built;
/// expected values come from the documented rules, Biopython 1.88 and EMBOSS semantics, not from the member under test.
/// TestSpec: TestSpecs/RESTR-FIND-001.md, RESTR-DIGEST-001.md, RESTR-FILTER-001.md ("Round 5").
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Round5_Tests
{
    private static readonly DnaSequence EcoRiLocus = new("AAAGAATTCAAA");

    // ------------------------------------------------------------------------------------------------
    // A5R-2: the scan start near int.MaxValue
    // ------------------------------------------------------------------------------------------------

    [TestCase(int.MaxValue)]
    [TestCase(int.MaxValue - 1)]
    [TestCase(int.MaxValue - 2)]
    public void FindPatternOccurrences_HugeStartIndex_IsAnEmptyResult_NotAnOverflow(int startIndex)
    {
        // Documented: occurrences starting before startIndex are skipped. i + pattern.Length overflowed int and the loop read out of range.
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGT", "AC", startIndex), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGT", "ACGT", startIndex), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGT", "A", startIndex), Is.Empty);
        });
    }

    [Test]
    public void FindPatternOccurrences_StartIndexAtTheLastWindow_StillFindsIt()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGTAC", "AC", 4), Is.EqualTo(new[] { 4 }));
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGTAC", "AC", 5), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGTAC", "ACGTAC", 0), Is.EqualTo(new[] { 0 }));
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ACGTAC", "ACGTACG", 0), Is.Empty, "pattern longer than the sequence");
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-3: one meaning for the enzyme-name array, eager, in every form
    // ------------------------------------------------------------------------------------------------

    private static readonly IReadOnlyList<Methylase> Dam = MethylaseSets.Dam;

    private static IEnumerable<TestCaseData> NameArrayForms()
    {
        // (name, call taking the names array)
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.FindSites(EcoRiLocus, n))).SetName("FindSites");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.FindSites(EcoRiLocus, MoleculeTopology.Circular, n))).SetName("FindSites_topology");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.FindSites(EcoRiLocus, MoleculeTopology.Linear, Dam, n))).SetName("FindSites_methylation");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.FindSites("AAAGAATTCAAA", MoleculeTopology.Linear, Dam, n))).SetName("FindSites_raw_methylation");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.Digest(EcoRiLocus, n))).SetName("Digest");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.Digest(EcoRiLocus, MoleculeTopology.Circular, n))).SetName("Digest_topology");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.Digest(EcoRiLocus, MoleculeTopology.Linear, Dam, n))).SetName("Digest_methylation");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.Digest("AAAGAATTCAAA", MoleculeTopology.Linear, Dam, n))).SetName("Digest_raw_methylation");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.PartialDigest(EcoRiLocus, n))).SetName("PartialDigest");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.PartialDigest(EcoRiLocus, MoleculeTopology.Circular, (int?)null, n))).SetName("PartialDigest_topology");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.PartialDigest(EcoRiLocus, MoleculeTopology.Linear, null, Dam, n))).SetName("PartialDigest_methylation");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.GetDigestSummary(EcoRiLocus, n))).SetName("GetDigestSummary");
        yield return new TestCaseData(new Func<string[], object>(n => RestrictionAnalyzer.GetDigestSummary(EcoRiLocus, MoleculeTopology.Circular, Dam, n))).SetName("GetDigestSummary_methylation");
    }

    [TestCaseSource(nameof(NameArrayForms))]
    public void EnzymeNameArray_NullIsAMissingArgument_EmptyIsAWrongValue_UnknownAndNullNamesAreRefusedAtTheCall(Func<string[], object> call)
    {
        // The BCL split, applied to the array and to its elements; no enumeration (ToList) anywhere: all of it is eager.
        Assert.Throws<ArgumentNullException>(() => call(null!), "null array");
        Assert.Throws<ArgumentException>(() => call(Array.Empty<string>()), "empty array");
        Assert.Throws<ArgumentException>(() => call(new[] { "NoSuchEnzyme1", "EcoRI" }), "unknown name");
        Assert.Throws<ArgumentNullException>(() => call(new string[] { null!, "EcoRI" }), "null element");
        Assert.Throws<ArgumentException>(() => call(new[] { "", "EcoRI" }), "empty element");
    }

    [Test]
    public void EnzymeNameArray_ValidNamesStillWork_InTheFormsThatWereRewritten()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindSites(EcoRiLocus, MoleculeTopology.Linear, Dam, "EcoRI", "BamHI").Select(s => s.Enzyme.Name).Distinct(),
                Is.EqualTo(new[] { "EcoRI" }), "per-name order, EcoRI only cuts here");
            Assert.That(RestrictionAnalyzer.Digest(EcoRiLocus, MoleculeTopology.Linear, Dam, "EcoRI").Select(f => f.Length), Is.EqualTo(new[] { 4, 8 }));
            Assert.That(RestrictionAnalyzer.PartialDigest(EcoRiLocus, MoleculeTopology.Linear, null, Dam, "EcoRI"), Is.Not.Empty);
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-5: an undefined topology is refused by every entry point (it meant linear in FindSites / CreateMap and
    // circular in Digest / GetDigestSummary / PartialDigest, and the summary echoed it back)
    // ------------------------------------------------------------------------------------------------

    private const MoleculeTopology Undefined = (MoleculeTopology)7;

    private static IEnumerable<TestCaseData> TopologyForms()
    {
        var ecori = new[] { RestrictionAnalyzer.GetEnzyme("EcoRI")! };
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindSites(EcoRiLocus, "EcoRI", Undefined))).SetName("FindSites_name");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindSites("AAAGAATTCAAA", "EcoRI", Undefined))).SetName("FindSites_raw_name");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindSites(EcoRiLocus, ecori[0], Undefined))).SetName("FindSites_enzyme");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindSites(EcoRiLocus, Undefined, "EcoRI"))).SetName("FindSites_names");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindSites(EcoRiLocus, Undefined, ecori))).SetName("FindSites_enzymes");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindSites("AAAGAATTCAAA", Undefined, ecori, null))).SetName("FindSites_raw_enzymes");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindAllSites(EcoRiLocus, Undefined))).SetName("FindAllSites");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.FindAllSites("AAAGAATTCAAA", Undefined, EnzymeSet.Common, null))).SetName("FindAllSites_raw");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.Digest(EcoRiLocus, Undefined, "EcoRI"))).SetName("Digest_names");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.Digest(EcoRiLocus, Undefined, ecori))).SetName("Digest_enzymes");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.Digest("AAAGAATTCAAA", Undefined, null, "EcoRI"))).SetName("Digest_raw_names");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.GetDigestSummary(EcoRiLocus, Undefined, "EcoRI"))).SetName("GetDigestSummary_names");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.GetDigestSummary(EcoRiLocus, Undefined, ecori))).SetName("GetDigestSummary_enzymes");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.PartialDigest(EcoRiLocus, Undefined, (int?)null, "EcoRI"))).SetName("PartialDigest_names");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.PartialDigest(EcoRiLocus, Undefined, (int?)null, ecori))).SetName("PartialDigest_enzymes");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.CreateMap(EcoRiLocus, Undefined, "EcoRI"))).SetName("CreateMap_names");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.CreateMap(EcoRiLocus, Undefined, ecori))).SetName("CreateMap_enzymes");
        yield return new TestCaseData(new Func<object>(() => RestrictionAnalyzer.CreateMap("AAAGAATTCAAA", Undefined, null, "EcoRI"))).SetName("CreateMap_raw_names");
    }

    [TestCaseSource(nameof(TopologyForms))]
    public void UndefinedTopology_IsRefusedAtTheCall(Func<object> call)
    {
        var e = Assert.Throws<ArgumentOutOfRangeException>(() => call());
        Assert.That(e!.ParamName, Is.EqualTo("topology"));
    }

    [Test]
    public void DefinedTopologies_AreUnaffected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.Digest(EcoRiLocus, MoleculeTopology.Linear, "EcoRI").Count(), Is.EqualTo(2));
            Assert.That(RestrictionAnalyzer.Digest(EcoRiLocus, MoleculeTopology.Circular, "EcoRI").Count(), Is.EqualTo(1));
            Assert.That(RestrictionAnalyzer.GetDigestSummary(EcoRiLocus, MoleculeTopology.Circular, "EcoRI").Topology, Is.EqualTo(MoleculeTopology.Circular));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-6: a malformed methylase is rejected at the call whatever the sequence
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void RawStringOverloads_RejectAMalformedMethylase_ForNullAndEmptySequencesToo()
    {
        var bad = new[] { new Methylase("bad", "GATC") }; // no methylated-base marker
        var ecori = new[] { RestrictionAnalyzer.GetEnzyme("EcoRI")! };

        foreach (string? sequence in new string?[] { null, "", "AAAGAATTCAAA" })
        {
            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindAllSites(sequence, MoleculeTopology.Linear, EnzymeSet.Common, bad), $"FindAllSites '{sequence}'");
                Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(sequence, MoleculeTopology.Linear, ecori, bad), $"FindSites(enzymes) '{sequence}'");
                Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(sequence, MoleculeTopology.Linear, bad, "EcoRI"), $"FindSites(names) '{sequence}'");
                Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(sequence ?? "", "EcoRI", MoleculeTopology.Linear, bad), $"FindSites(name) '{sequence}' (the sibling that always validated)");
            });
        }
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-7: cut offsets near int.MaxValue on a circular molecule
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void CustomEnzyme_HugeCutOffsets_OnACircle_ReduceInLongArithmetic()
    {
        // GAATTC at 3 of a 12-mer. Forward cut: (3 + 2147483647) mod 12 = 10; the palindromic site's reverse record:
        // (3 + 6 - 2147483647) mod 12 = 2 (exact arithmetic, checked outside the library).
        var enzyme = new RestrictionEnzyme("Far", "GAATTC", int.MaxValue, int.MaxValue, Organism: "");

        var sites = RestrictionAnalyzer.FindSites("AAAGAATTCAAA", MoleculeTopology.Circular, new[] { enzyme }, null).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Single(s => s.IsForwardStrand).CutPosition, Is.EqualTo(10));
            Assert.That(sites.Single(s => !s.IsForwardStrand).CutPosition, Is.EqualTo(2));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-8: the raw-string entry points upper-case ASCII; output sequences are upper case
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void RawStringOutputs_AreUpperCase_ForASoftMaskedMolecule()
    {
        const string softMasked = "nnnGAATTCaaa";

        var digest = RestrictionAnalyzer.Digest(softMasked, MoleculeTopology.Linear, null, "EcoRI").ToList();
        var sites = RestrictionAnalyzer.FindSites("nnngaattcaaa", "EcoRI", MoleculeTopology.Linear).ToList();
        var partial = RestrictionAnalyzer.PartialDigest(softMasked, MoleculeTopology.Linear, null, null, "EcoRI").ToList();
        var map = RestrictionAnalyzer.CreateMap(softMasked, MoleculeTopology.Linear, null, "EcoRI");

        Assert.Multiple(() =>
        {
            Assert.That(digest.Select(f => f.Sequence), Is.EqualTo(new[] { "NNNG", "AATTCAAA" }), "the soft-mask is not carried over");
            Assert.That(digest.Select(f => f.StartPosition), Is.EqualTo(new[] { 0, 4 }), "coordinates and lengths are unaffected");
            Assert.That(sites.Select(s => s.RecognizedSequence).Distinct(), Is.EqualTo(new[] { "GAATTC" }));
            Assert.That(partial.Select(f => f.Fragment.Sequence), Has.All.Matches<string>(x => x == x.ToUpperInvariant()));
            Assert.That(map.Sites.Select(s => s.RecognizedSequence).Distinct(), Is.EqualTo(new[] { "GAATTC" }));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T1: a Type IIB site counts once (AlfI = GCA N6 TGC cuts 10/12 before and 24/22 after the site)
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void CreateMap_OneTypeIibSite_IsOneCleavage_AndAUniqueCutter()
    {
        string seq = new string('A', 40) + "GCAAAAAAATGC" + new string('A', 40);

        var map = RestrictionAnalyzer.CreateMap(new DnaSequence(seq), "AlfI");

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(1), "the two cut records of one two-sided site count once (F108)");
            Assert.That(map.UniqueCutters, Is.EqualTo(new[] { "AlfI" }));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T2 / T3: a methylated base inside a window that crosses the origin of a circle blocks the site
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void Methylation_BlocksAWindowThatCrossesTheCircularOrigin()
    {
        // GATC occupies 11, 0, 1, 2; Dam methylates the A of GATC (position 0), which lies inside MboI's window.
        const string circle = "ATCAAAAAAAAG";

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(circle), "MboI", MoleculeTopology.Circular, MethylaseSets.Dam), Is.Empty,
                "blocked by Dam");
            Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence(circle), "MboI", MoleculeTopology.Circular, null).Count(), Is.GreaterThan(0),
                "without methylation the wrapped GATC is a site (control)");
            Assert.That(RestrictionAnalyzer.FindAllSites((string?)circle, MoleculeTopology.Circular, EnzymeSet.Common, MethylaseSets.Dam)
                .Where(s => s.Enzyme.Name == "MboI"), Is.Empty, "the raw-string FindAllSites passes the circular flag to the mask");
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T4: names as supplied (documented on DigestSummary.EnzymesUsed)
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void RawGetDigestSummary_EnzymesUsed_AreTheNamesAsSupplied()
    {
        var summary = RestrictionAnalyzer.GetDigestSummary("AAAGAATTCAAA", MoleculeTopology.Linear, null, "ecori");

        Assert.That(summary.EnzymesUsed, Is.EqualTo(new[] { "ecori" }));
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T5: DigestFragment value equality compares the enzyme lists element-wise
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void DigestFragment_Equality_ComparesTheEnzymeListsByContent()
    {
        DigestFragment Make(params string[] left) =>
            new("ACGT", 0, 4, "A", null, 1) { LeftEnzymes = left };

        var first = Make("A", "B");
        var sameContent = Make("A", "B");
        Assert.Multiple(() =>
        {
            Assert.That(first, Is.EqualTo(sameContent), "equal content, different list instances");
            Assert.That(first.GetHashCode(), Is.EqualTo(sameContent.GetHashCode()));
            Assert.That(first, Is.Not.EqualTo(Make("A", "C")), "same count, different second enzyme");
            Assert.That(first, Is.Not.EqualTo(Make("A")));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T8: the fragment-count limit is inclusive
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void PartialDigest_CountExactlyAtTheLimit_IsAccepted_OneStepOverIsRefused()
    {
        // 1 000 EcoRI cuts on a circle: k * (m + 1) fragments for m < k - 1 -> m = 99 gives exactly 100 000, m = 100 gives 101 000.
        var circle = new DnaSequence(string.Concat(Enumerable.Repeat("GAATTCAAAA", 1000)));
        Assert.That(RestrictionAnalyzer.MaxPartialDigestFragments, Is.EqualTo(100_000));

        // The limit is checked when enumeration starts; one element is enough and nothing 100 000 long is materialised.
        Assert.DoesNotThrow(() => RestrictionAnalyzer.PartialDigest(circle, MoleculeTopology.Circular, 99, "EcoRI").Take(1).ToList());
        var e = Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.PartialDigest(circle, MoleculeTopology.Circular, 100, "EcoRI").Take(1).ToList());
        Assert.That(e!.ParamName, Is.EqualTo("maxMissedCuts"));
    }
}
