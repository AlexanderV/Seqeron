using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier (fuzz) coverage for the restriction-analysis argument contracts that B08 audit rounds 4-5 unified
/// and the older fuzz suites do not reach:
///   - an undefined <see cref="MoleculeTopology"/> value is refused at the call by every member that takes one
///     (a cast <c>(MoleculeTopology)2</c> used to fall through to the circular branch of a two-way ternary);
///   - the <c>params string[]</c> name arrays share one contract (null array / empty array / null name / empty name /
///     unknown name), resolved eagerly - no member defers it to first enumeration.
///
/// The exception split is the BCL's <c>ArgumentException.ThrowIfNullOrEmpty</c>: null is ArgumentNullException, empty or
/// unknown is ArgumentException. Test Units: RESTR-FIND-001, RESTR-DIGEST-001.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("MolTools")]
[Category("RESTR-FIND-001")]
[Category("RESTR-DIGEST-001")]
public class B08RestrictionAuditFuzzTests
{
    private const string Seq = "AAAGAATTCAAAGGATCCAAAGAATTCAAAGGATCCAAA";

    private static readonly DnaSequence Dna = new(Seq);
    private static readonly RestrictionEnzyme[] EcoRI = [RestrictionAnalyzer.GetEnzyme("EcoRI")!];
    private static readonly IReadOnlyList<Methylase> Dam = MethylaseSets.Dam;

    private static Arbitrary<int> UndefinedTopologies() =>
        Gen.OneOf(Gen.Choose(int.MinValue, -1), Gen.Choose(2, int.MaxValue)).ToArbitrary();

    #region Undefined topology

    /// <summary>
    /// Every public member that takes a topology throws ArgumentOutOfRangeException (parameter <c>topology</c>) for an
    /// integer outside the enum, at the call - no enumeration - for any such value.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EveryTopologyTakingMember_RejectsAnUndefinedTopology_AtTheCall()
    {
        return Prop.ForAll(UndefinedTopologies(), raw =>
        {
            var t = (MoleculeTopology)raw;
            var calls = new (string Name, Action Call)[]
            {
                ("FindSites(Dna, name, t)", () => RestrictionAnalyzer.FindSites(Dna, "EcoRI", t)),
                ("FindSites(string, name, t)", () => RestrictionAnalyzer.FindSites(Seq, "EcoRI", t)),
                ("FindSites(Dna, enzyme, t)", () => RestrictionAnalyzer.FindSites(Dna, EcoRI[0], t)),
                ("FindSites(Dna, t, names)", () => RestrictionAnalyzer.FindSites(Dna, t, "EcoRI")),
                ("FindSites(Dna, t, enzymes)", () => RestrictionAnalyzer.FindSites(Dna, t, (IEnumerable<RestrictionEnzyme>)EcoRI)),
                ("FindAllSites(Dna, t)", () => RestrictionAnalyzer.FindAllSites(Dna, t)),
                ("FindAllSites(Dna, t, set)", () => RestrictionAnalyzer.FindAllSites(Dna, t, EnzymeSet.Rebase)),
                ("FindAllSites(string, t, set, m)", () => RestrictionAnalyzer.FindAllSites(Seq, t, EnzymeSet.Common, null)),
                ("Digest(Dna, t, names)", () => RestrictionAnalyzer.Digest(Dna, t, "EcoRI")),
                ("Digest(Dna, t, enzymes)", () => RestrictionAnalyzer.Digest(Dna, t, (IEnumerable<RestrictionEnzyme>)EcoRI)),
                ("GetDigestSummary(Dna, t, names)", () => RestrictionAnalyzer.GetDigestSummary(Dna, t, "EcoRI")),
                ("GetDigestSummary(Dna, t, enzymes)", () => RestrictionAnalyzer.GetDigestSummary(Dna, t, (IEnumerable<RestrictionEnzyme>)EcoRI)),
                ("CreateMap(Dna, t, names)", () => RestrictionAnalyzer.CreateMap(Dna, t, "EcoRI")),
                ("CreateMap(Dna, t, enzymes)", () => RestrictionAnalyzer.CreateMap(Dna, t, (IEnumerable<RestrictionEnzyme>)EcoRI)),
                ("PartialDigest(Dna, t, names)", () => RestrictionAnalyzer.PartialDigest(Dna, t, "EcoRI")),
                ("PartialDigest(Dna, t, max, names)", () => RestrictionAnalyzer.PartialDigest(Dna, t, 1, "EcoRI")),
                ("PartialDigest(Dna, t, max, enzymes)", () => RestrictionAnalyzer.PartialDigest(Dna, t, 1, (IEnumerable<RestrictionEnzyme>)EcoRI)),
                ("FindSites(Dna, name, t, m)", () => RestrictionAnalyzer.FindSites(Dna, "EcoRI", t, Dam)),
                ("FindSites(Dna, enzyme, t, m)", () => RestrictionAnalyzer.FindSites(Dna, EcoRI[0], t, Dam)),
                ("FindSites(Dna, t, m, names)", () => RestrictionAnalyzer.FindSites(Dna, t, Dam, "EcoRI")),
                ("FindSites(Dna, t, enzymes, m)", () => RestrictionAnalyzer.FindSites(Dna, t, (IEnumerable<RestrictionEnzyme>)EcoRI, Dam)),
                ("FindSites(string, name, t, m)", () => RestrictionAnalyzer.FindSites(Seq, "EcoRI", t, Dam)),
                ("Digest(Dna, t, m, names)", () => RestrictionAnalyzer.Digest(Dna, t, Dam, "EcoRI")),
                ("Digest(Dna, t, enzymes, m)", () => RestrictionAnalyzer.Digest(Dna, t, (IEnumerable<RestrictionEnzyme>)EcoRI, Dam)),
                ("GetDigestSummary(Dna, t, m, names)", () => RestrictionAnalyzer.GetDigestSummary(Dna, t, Dam, "EcoRI")),
                ("CreateMap(Dna, t, m, names)", () => RestrictionAnalyzer.CreateMap(Dna, t, Dam, "EcoRI")),
                ("PartialDigest(Dna, t, max, m, names)", () => RestrictionAnalyzer.PartialDigest(Dna, t, 1, Dam, "EcoRI")),
                ("FindSites(raw, t, enzymes, m)", () => RestrictionAnalyzer.FindSites(Seq, t, (IEnumerable<RestrictionEnzyme>)EcoRI, Dam)),
                ("Digest(raw, t, m, names)", () => RestrictionAnalyzer.Digest(Seq, t, Dam, "EcoRI")),
                ("Digest(raw, t, enzymes, m)", () => RestrictionAnalyzer.Digest(Seq, t, (IEnumerable<RestrictionEnzyme>)EcoRI, Dam)),
                ("GetDigestSummary(raw, t, m, names)", () => RestrictionAnalyzer.GetDigestSummary(Seq, t, Dam, "EcoRI")),
                ("CreateMap(raw, t, m, names)", () => RestrictionAnalyzer.CreateMap(Seq, t, Dam, "EcoRI")),
                ("PartialDigest(raw, t, max, m, names)", () => RestrictionAnalyzer.PartialDigest(Seq, t, 1, Dam, "EcoRI")),
            };

            foreach (var (name, call) in calls)
            {
                try
                {
                    call();
                    return false.Label($"{name} accepted topology {raw}");
                }
                catch (ArgumentOutOfRangeException ex) when (ex.ParamName == "topology")
                {
                }
                catch (Exception ex)
                {
                    return false.Label($"{name} threw {ex.GetType().Name} ({ex.Message}) for topology {raw}");
                }
            }

            return true.ToProperty();
        });
    }

    /// <summary>The two defined values still work through the same members (the guard is not a blanket refusal).</summary>
    [TestCase(MoleculeTopology.Linear)]
    [TestCase(MoleculeTopology.Circular)]
    public void DefinedTopologies_AreAccepted_ByTheSameMembers(MoleculeTopology t)
    {
        RestrictionAnalyzer.FindSites(Dna, "EcoRI", t).Should().NotBeEmpty();
        RestrictionAnalyzer.Digest(Dna, t, "EcoRI").Should().NotBeEmpty();
        RestrictionAnalyzer.GetDigestSummary(Dna, t, "EcoRI").Topology.Should().Be(t);
        RestrictionAnalyzer.CreateMap(Dna, t, "EcoRI").TotalSites.Should().BeGreaterThan(0);
        RestrictionAnalyzer.PartialDigest(Dna, t, "EcoRI").Should().NotBeEmpty();
    }

    #endregion

    #region One enzyme-name contract for every params string[] entry point

    private static IEnumerable<(string Name, Func<string[], object> Call)> NameArrayMembers()
    {
        yield return ("FindSites(Dna, names)", names => RestrictionAnalyzer.FindSites(Dna, names));
        yield return ("FindSites(Dna, t, names)", names => RestrictionAnalyzer.FindSites(Dna, MoleculeTopology.Circular, names));
        yield return ("Digest(Dna, names)", names => RestrictionAnalyzer.Digest(Dna, names));
        yield return ("Digest(Dna, t, names)", names => RestrictionAnalyzer.Digest(Dna, MoleculeTopology.Circular, names));
        yield return ("GetDigestSummary(Dna, names)", names => RestrictionAnalyzer.GetDigestSummary(Dna, names));
        yield return ("GetDigestSummary(Dna, t, names)", names => RestrictionAnalyzer.GetDigestSummary(Dna, MoleculeTopology.Circular, names));
        yield return ("PartialDigest(Dna, names)", names => RestrictionAnalyzer.PartialDigest(Dna, names));
        yield return ("PartialDigest(Dna, t, names)", names => RestrictionAnalyzer.PartialDigest(Dna, MoleculeTopology.Circular, names));
        yield return ("PartialDigest(Dna, t, max, names)", names => RestrictionAnalyzer.PartialDigest(Dna, MoleculeTopology.Linear, 1, names));
        yield return ("FindSites(Dna, t, m, names)", names => RestrictionAnalyzer.FindSites(Dna, MoleculeTopology.Linear, Dam, names));
        yield return ("Digest(Dna, t, m, names)", names => RestrictionAnalyzer.Digest(Dna, MoleculeTopology.Linear, Dam, names));
        yield return ("GetDigestSummary(Dna, t, m, names)", names => RestrictionAnalyzer.GetDigestSummary(Dna, MoleculeTopology.Linear, Dam, names));
        yield return ("PartialDigest(Dna, t, max, m, names)", names => RestrictionAnalyzer.PartialDigest(Dna, MoleculeTopology.Linear, 1, Dam, names));
        yield return ("FindSites(raw, t, m, names)", names => RestrictionAnalyzer.FindSites(Seq, MoleculeTopology.Linear, Dam, names));
        yield return ("Digest(raw, t, m, names)", names => RestrictionAnalyzer.Digest(Seq, MoleculeTopology.Linear, Dam, names));
        yield return ("GetDigestSummary(raw, t, m, names)", names => RestrictionAnalyzer.GetDigestSummary(Seq, MoleculeTopology.Linear, Dam, names));
        yield return ("PartialDigest(raw, t, max, m, names)", names => RestrictionAnalyzer.PartialDigest(Seq, MoleculeTopology.Linear, 1, Dam, names));
    }

    public static IEnumerable<string> NameArrayMemberNames => NameArrayMembers().Select(m => m.Name);

    private static Func<string[], object> Member(string name) => NameArrayMembers().First(m => m.Name == name).Call;

    /// <summary>A null array is ArgumentNullException and an empty array is ArgumentException, at the call.</summary>
    [TestCaseSource(nameof(NameArrayMemberNames))]
    public void NullArray_IsArgumentNull_EmptyArray_IsArgument_AtTheCall(string member)
    {
        var call = Member(member);

        var nullArray = () => call(null!);
        var emptyArray = () => call([]);

        nullArray.Should().Throw<ArgumentNullException>();
        emptyArray.Should().Throw<ArgumentException>().Which.Should().NotBeOfType<ArgumentNullException>();
    }

    /// <summary>
    /// A null element is ArgumentNullException; an empty, whitespace-only, unknown or REBASE-but-undefined name is
    /// ArgumentException; a name hiding behind a valid one is still found (every element is resolved), at the call.
    /// </summary>
    [TestCaseSource(nameof(NameArrayMemberNames))]
    public void BadElement_IsRefused_AtTheCall_WhereverItSits(string member)
    {
        var call = Member(member);

        foreach (var bad in new[] { null, "", " ", "NoSuchEnzyme", "EcoRI2", "EcoRI\0" })
            foreach (var arrangement in new Func<string?, string?[]>[] { b => [b], b => ["EcoRI", b], b => [b, "EcoRI"] })
            {
                var names = arrangement(bad)!;
                var act = () => call(names!);

                var thrown = act.Should().Throw<ArgumentException>($"names = [{string.Join(", ", names.Select(n => n is null ? "null" : $"'{n}'"))}]").Which;
                if (bad is null)
                    thrown.Should().BeOfType<ArgumentNullException>();
                else
                    thrown.Should().NotBeOfType<ArgumentNullException>();
            }
    }

    /// <summary>
    /// Surrounding whitespace is ignored, as in <c>MethylaseSets.FromName</c>, and REBASE's hyphen resolves to the table's
    /// underscore (<c>CviKI-1</c> / <c>CviKI_1</c>), in any casing - for every name-array member, not just GetEnzyme.
    /// </summary>
    [TestCaseSource(nameof(NameArrayMemberNames))]
    public void WhitespaceAndHyphenSpellings_Resolve_ThroughEveryMember(string member)
    {
        var call = Member(member);

        foreach (var spelling in new[] { " EcoRI", "ecori\t", "\n EcOrI  ", "CviKI-1", "cviki_1", " CVIKI-1 " })
        {
            var act = () => call([spelling]);

            act.Should().NotThrow($"'{spelling}' names a REBASE enzyme");
        }
    }

    /// <summary>Names are case-insensitive: any casing of a known enzyme resolves, and the result equals the canonical spelling's.</summary>
    [FsCheck.NUnit.Property]
    public Property EnzymeNameLookup_IsCaseInsensitive_ForEveryCasing()
    {
        var casings = Gen.Elements("EcoRI", "ecori", "ECORI", "eCoRi", "EcOri", "ecoRI").ToArbitrary();
        return Prop.ForAll(casings, name =>
        {
            var expected = RestrictionAnalyzer.FindSites(Dna, "EcoRI").Select(s => (s.Position, s.IsForwardStrand, s.CutPosition)).ToList();
            var actual = RestrictionAnalyzer.FindSites(Dna, name).Select(s => (s.Position, s.IsForwardStrand, s.CutPosition)).ToList();

            return expected.SequenceEqual(actual).Label(name);
        });
    }

    #endregion
}
