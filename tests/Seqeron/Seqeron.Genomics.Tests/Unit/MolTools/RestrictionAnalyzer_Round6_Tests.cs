using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 6 (the confirmation audit of round 5) on the restriction slice: B08 A6-2 .. A6-5.
/// Each test failed on the code the auditor read (bec761ac) or, for the mutation rows, on the mutant the auditor built;
/// expected values come from the documented rules and exact arithmetic, not from the member under test.
/// TestSpec: TestSpecs/RESTR-FIND-001.md, RESTR-DIGEST-001.md ("Round 6").
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Round6_Tests
{
    private const string EcoRiLocus = "AAAGAATTCAAA";

    // ------------------------------------------------------------------------------------------------
    // A6-2: U is not one of the 15 IUPAC nucleotide codes, so it matches nothing - on either strand. The reverse strand was
    // read through the RNA-aware complement (U -> A), so GAAUUC was seen there as GAATTC.
    // ------------------------------------------------------------------------------------------------

    [TestCase("AAAGAAUUCAAA")]
    [TestCase("AAAGAATUCAAA")]
    [TestCase("AAAGAAUTCAAA")]
    [TestCase("aaagaauucaaa")]
    [TestCase("UAAAGAATTCAAA", Description = "a control: the U is outside the site, so the site is found")]
    public void RawSequence_AUracilInsideASite_IsNotASite_OnEitherStrand(string sequence)
    {
        bool uInsideSite = sequence.ToUpperInvariant().Contains("GAAU") || sequence.ToUpperInvariant().Contains("GAATU") || sequence.ToUpperInvariant().Contains("GAAUT");
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI", topology).ToList();
            var digest = RestrictionAnalyzer.Digest(sequence, topology, null, "EcoRI").ToList();

            if (uInsideSite)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(sites, Is.Empty, $"{topology}: no EcoRI record");
                    Assert.That(digest, Has.Count.EqualTo(1), $"{topology}: nothing is cut");
                });
            }
            else
            {
                Assert.That(sites, Is.Not.Empty, $"{topology}: the control site is found");
            }
        }
    }

    [Test]
    public void RawSequence_AUracil_IsAlsoNotReadAcrossTheCircularOrigin()
    {
        // GAATTC split across the origin with a U in the wrapped part: GAA | UUC.
        Assert.That(RestrictionAnalyzer.FindSites("UUCAAAAAAAAAGAA", "EcoRI", MoleculeTopology.Circular), Is.Empty);
        // The control with T: found (forward and reverse record of the palindromic site).
        Assert.That(RestrictionAnalyzer.FindSites("TTCAAAAAAAAAGAA", "EcoRI", MoleculeTopology.Circular), Is.Not.Empty);
    }

    [Test]
    public void Methylation_AUracil_DoesNotCreateAMethylaseContextOnTheReverseStrand()
    {
        // Dcm = CmCWGG / CmCTGG. 'CCUGG' is not a Dcm site (U matches nothing), but its reverse complement with U read as A,
        // CCAGG, is: it marked the two G's at forward 3 and 4 and blocked the HaeIII site GGCC at 3..6.
        var sites = RestrictionAnalyzer.FindSites("CCUGGCCAAAA", MoleculeTopology.Linear, MethylaseSets.Dcm, "HaeIII").ToList();

        Assert.That(sites.Select(s => s.Position), Does.Contain(3), "HaeIII (GGCC at 3..6) is not blocked by a context that is not there");

        // The control: a genuine Dcm context (CCTGG at 2..6 marks its two C's, forward 2 and 3) does block a site that covers them -
        // HaeIII GGCC at 0..3 - so the mask is live and the assertion above is not vacuous.
        Assert.That(RestrictionAnalyzer.FindSites("GGCCTGGAAAA", MoleculeTopology.Linear, MethylaseSets.Dcm, "HaeIII"), Is.Empty);
        Assert.That(RestrictionAnalyzer.FindSites("GGCCTGGAAAA", MoleculeTopology.Linear, null, "HaeIII"), Is.Not.Empty);
    }

    [Test]
    public void Overhang_OfASiteWhoseWindowHoldsAUracil_KeepsTheUracil()
    {
        // BsaI (GGTCTC 1/5) cuts at 9 (top) and 13 (bottom), so its 4-nt 5' overhang is the molecule's bases 9..12 = UCGA. The right
        // fragment carries it as it reads; the left fragment carries it on the bottom strand, i.e. reverse-complemented. A base that is
        // not one of the 15 codes has no complement: the U stays a U (it was reported as A).
        var digest = RestrictionAnalyzer.Digest("AAGGTCTCAUCGAAAAAAA", MoleculeTopology.Linear, null, "BsaI").ToList();

        Assert.That(digest, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(digest[1].LeftOverhang, Is.EqualTo("UCGA"));
            Assert.That(digest[0].RightOverhang, Is.EqualTo("TCGU"));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A6-3: the stagger of a caller-supplied enzyme is bounded and computed in long
    // ------------------------------------------------------------------------------------------------

    public static IEnumerable<TestCaseData> UnboundedStaggers()
    {
        yield return new TestCaseData(0, int.MinValue).SetName("0_MinValue");
        yield return new TestCaseData(int.MinValue, 0).SetName("MinValue_0");
        yield return new TestCaseData(int.MaxValue, int.MinValue).SetName("MaxValue_MinValue");
        yield return new TestCaseData(int.MaxValue, 5).SetName("MaxValue_5");
        yield return new TestCaseData(5, int.MinValue).SetName("5_MinValue");
        yield return new TestCaseData(50_000_000, 5).SetName("50M_5");
        yield return new TestCaseData(0, RestrictionAnalyzer.MaxCutStagger + 1).SetName("limit_plus_one_forward_lower");
        yield return new TestCaseData(RestrictionAnalyzer.MaxCutStagger + 1, 0).SetName("limit_plus_one_forward_higher");
    }

    [TestCaseSource(nameof(UnboundedStaggers))]
    public void CustomEnzyme_StaggerBeyondTheLimit_IsRefusedAtTheCall_ByEveryMember(int forward, int reverse)
    {
        var enzyme = new RestrictionEnzyme("X", "GAATTC", forward, reverse, Organism: "");
        var dna = new DnaSequence("AAAGAATTCAAAGAATTCAAAGGGGGGGGGGGG");
        var ecoRi = RestrictionAnalyzer.GetEnzyme("EcoRI")!;

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Circular), "FindSites circular");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Linear), "FindSites linear");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, new[] { enzyme }), "Digest circular");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, new[] { enzyme }), "Digest linear");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.Digest(dna.Sequence, MoleculeTopology.Circular, new[] { enzyme }, null), "Digest raw");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetDigestSummary(dna, MoleculeTopology.Circular, new[] { enzyme }), "GetDigestSummary");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.CreateMap(dna, MoleculeTopology.Circular, new[] { enzyme }), "CreateMap");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Circular, (int?)null, new[] { enzyme }), "PartialDigest");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.GetOverhangSequence(enzyme), "GetOverhangSequence");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(enzyme, ecoRi), "AreCompatible first");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.AreCompatible(ecoRi, enzyme), "AreCompatible second");
        });
    }

    [Test]
    public void CustomEnzyme_StaggerAtTheLimit_IsAccepted_EvenAtTheEdgesOfTheIntRange()
    {
        int limit = RestrictionAnalyzer.MaxCutStagger;

        // Exact stagger = limit: a window of `limit` unknown bases, however far from the site the cuts lie.
        foreach (var (forward, reverse) in new[] { (0, limit), (limit, 0), (int.MaxValue - limit, int.MaxValue), (int.MinValue, int.MinValue + limit) })
        {
            var enzyme = new RestrictionEnzyme("Edge", "GAATTC", forward, reverse, Organism: "");

            string? overhang = RestrictionAnalyzer.GetOverhangSequence(enzyme);

            Assert.That(overhang, Has.Length.EqualTo(limit), $"({forward}, {reverse})");
        }
    }

    [Test]
    public void EveryRebaseEnzyme_IsWithinTheLimit_AndTheLargestStaggerIsNine()
    {
        var staggers = RestrictionAnalyzer.RebaseEnzymes.Values
            .Select(e => Math.Abs((long)e.CutPositionForward - e.CutPositionReverse))
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(staggers.Max(), Is.EqualTo(9), "TspRI / TscAI");
            Assert.That(RestrictionAnalyzer.MaxCutStagger, Is.GreaterThan(9));
            foreach (var enzyme in RestrictionAnalyzer.RebaseEnzymes.Values)
                Assert.DoesNotThrow(() => RestrictionAnalyzer.GetOverhangSequence(enzyme), enzyme.Name);
        });
    }

    [Test]
    public void SecondCutStagger_IsComparedInLong_NotWrappedInInt()
    {
        // first stagger 1; second (int.MaxValue, int.MinValue) has exact stagger 4294967295 whose int subtraction wraps to -1 ...
        var first = new RestrictionEnzyme("TwoSided", "GAATTC", 0, 1, Organism: "")
        {
            SecondCut = (int.MaxValue, int.MinValue)
        };

        // ... and so used to equal the first pair's stagger of -1 and pass the "same end" check.
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(new DnaSequence("AAAGAATTCAAA"), first, MoleculeTopology.Circular));
    }

    // ------------------------------------------------------------------------------------------------
    // A6-4: every public member that takes a topology refuses an undefined one; every params string[] member shares one contract
    // (found by reflection, so a new overload cannot be added without the guard)
    // ------------------------------------------------------------------------------------------------

    private static readonly RestrictionEnzyme[] EcoRiList = [RestrictionAnalyzer.GetEnzyme("EcoRI")!];

    private static object?[] ArgumentsFor(MethodInfo method, MoleculeTopology topology, string[]? names)
    {
        return method.GetParameters().Select<ParameterInfo, object?>(p =>
        {
            var t = p.ParameterType;
            if (t == typeof(MoleculeTopology)) return topology;
            if (t == typeof(DnaSequence)) return new DnaSequence(EcoRiLocus);
            if (t == typeof(string[])) return names;
            if (t == typeof(string)) return p.Name == "sequence" ? EcoRiLocus : "EcoRI";
            if (t == typeof(RestrictionEnzyme)) return EcoRiList[0];
            if (t == typeof(IEnumerable<RestrictionEnzyme>)) return EcoRiList;
            if (t == typeof(IReadOnlyList<Methylase>)) return null;
            if (t == typeof(int?)) return null;
            if (t == typeof(EnzymeSet)) return EnzymeSet.Common;
            throw new NotSupportedException($"{method.Name}: teach ArgumentsFor about parameter '{p.Name}' of type {t.Name}");
        }).ToArray();
    }

    private static Exception? CallAndCatch(MethodInfo method, object?[] args)
    {
        try
        {
            method.Invoke(null, args);
            return null;
        }
        catch (TargetInvocationException e)
        {
            return e.InnerException;
        }
    }

    private static IEnumerable<MethodInfo> PublicMembers() =>
        typeof(RestrictionAnalyzer).GetMethods(BindingFlags.Public | BindingFlags.Static).Where(m => !m.IsGenericMethod);

    [Test]
    public void EveryPublicMemberThatTakesATopology_RefusesAnUndefinedOne_AtTheCall()
    {
        var members = PublicMembers().Where(m => m.GetParameters().Any(p => p.ParameterType == typeof(MoleculeTopology))).ToList();
        Assert.That(members, Has.Count.GreaterThanOrEqualTo(41), "the public surface that takes a topology (41 at round 6)");

        var failures = new List<string>();
        foreach (var m in members)
        {
            var e = CallAndCatch(m, ArgumentsFor(m, (MoleculeTopology)7, ["EcoRI"]));
            if (e is not ArgumentOutOfRangeException { ParamName: "topology" })
                failures.Add($"{m} -> {(e is null ? "accepted" : e.GetType().Name)}");
        }

        Assert.That(failures, Is.Empty, string.Join(Environment.NewLine, failures));
    }

    [Test]
    public void EveryPublicMemberThatTakesATopology_StillWorks_ForTheDefinedOnes()
    {
        var members = PublicMembers().Where(m => m.GetParameters().Any(p => p.ParameterType == typeof(MoleculeTopology))).ToList();

        var failures = new List<string>();
        foreach (var m in members)
            foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
            {
                var e = CallAndCatch(m, ArgumentsFor(m, topology, ["EcoRI"]));
                if (e is not null)
                    failures.Add($"{m} {topology} -> {e.GetType().Name}: {e.Message}");
            }

        Assert.That(failures, Is.Empty, string.Join(Environment.NewLine, failures));
    }

    [Test]
    public void EveryParamsEnzymeNameMember_SharesOneContract_AtTheCall()
    {
        var members = PublicMembers().Where(m => m.GetParameters().Any(p => p.ParameterType == typeof(string[]) && p.Name == "enzymeNames")).ToList();
        Assert.That(members, Has.Count.GreaterThanOrEqualTo(21), "the public surface that takes enzyme names (21 at round 6)");

        var failures = new List<string>();
        foreach (var m in members)
        {
            bool emptyMeansEveryCommonEnzyme = m.Name == "CreateMap";
            var cases = new (string Label, string[]? Names, Type? Expected)[]
            {
                ("null array", null, typeof(ArgumentNullException)),
                ("empty array", [], emptyMeansEveryCommonEnzyme ? null : typeof(ArgumentException)),
                ("{null}", [null!], typeof(ArgumentNullException)),
                ("{\"\"}", [""], typeof(ArgumentException)),
                ("{\"Bogus\"}", ["Bogus"], typeof(ArgumentException)),
                ("{\"EcoRI\",\"Bogus\"}", ["EcoRI", "Bogus"], typeof(ArgumentException)),
            };

            foreach (var (label, names, expected) in cases)
            {
                var e = CallAndCatch(m, ArgumentsFor(m, MoleculeTopology.Linear, names));
                bool ok = expected is null ? e is null : e is not null && e.GetType() == expected;
                if (!ok)
                    failures.Add($"{m} with {label} -> {(e is null ? "accepted" : e.GetType().Name)}, expected {expected?.Name ?? "no exception"}");
            }
        }

        Assert.That(failures, Is.Empty, string.Join(Environment.NewLine, failures));
    }

    // ------------------------------------------------------------------------------------------------
    // A6-5: the circular reverse-record cut is reduced in long arithmetic for offsets of either sign
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void CustomEnzyme_MinValueOffsets_OnACircle_ReduceInLongArithmetic()
    {
        // GAATTC at 3 of a 12-mer. Forward cut: (3 + int.MinValue) mod 12 = 7; the palindromic site's reverse record:
        // (3 + 6 - int.MinValue) mod 12 = 5 (2^31 = 8 mod 12; checked outside the library). Same enzyme in int arithmetic wraps to 1.
        var enzyme = new RestrictionEnzyme("Far", "GAATTC", int.MinValue, int.MinValue, Organism: "");

        var sites = RestrictionAnalyzer.FindSites(EcoRiLocus, MoleculeTopology.Circular, new[] { enzyme }, null).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Single(s => s.IsForwardStrand).CutPosition, Is.EqualTo(7));
            Assert.That(sites.Single(s => !s.IsForwardStrand).CutPosition, Is.EqualTo(5));
        });
    }
}
