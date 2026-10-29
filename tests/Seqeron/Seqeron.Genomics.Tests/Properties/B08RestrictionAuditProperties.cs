using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based (heavy-tier) coverage for the restriction-analysis behaviours that B08 audit rounds 4-5 fixed or
/// added: circular coordinates, the partial-digest combinatorics, IUPAC matching of the recognition pattern against an
/// ambiguous sequence, the methylation mask, and raw-string case handling.
///
/// Each relation is a definition, not an observed output:
///   - a circular molecule has no origin, so rotating its sequence rotates every coordinate and changes nothing else;
///   - a linear molecule with k cuts has (k+1)(k+2)/2 partial-digest fragments and a circular one k^2 + 1
///     (every contiguous run of segments between cuts, plus the uncut circle);
///   - methylation can only block a site, never create one;
///   - recognition-pattern N matches any sequence symbol, every other pattern code matches only a concrete base
///     (EMBOSS restrict / Biopython FormattedSeq semantics), so an ambiguity code in the sequence blocks a concrete
///     site code and is accepted inside an N run.
///
/// Test Units: RESTR-FIND-001, RESTR-DIGEST-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
[Category("RESTR-FIND-001")]
[Category("RESTR-DIGEST-001")]
public class B08RestrictionAuditProperties
{
    private static readonly char[] Bases = ['A', 'C', 'G', 'T'];

    // Concrete sites of the enzymes below plus the two interrupted ones (BglI GCCNNNNNGGC, SfiI GGCCNNNNNGGCC) and the
    // Dam / Dcm contexts (GATC, CCAGG, CCTGG).
    private static readonly string[] PlantedSites =
        ["GAATTC", "GGATCC", "AAGCTT", "CCCGGG", "GGTCTC", "GCCACGTGGGC", "GGCCACGTGGGCC", "GATC", "CCAGG", "CCTGG", "ATCGAT"];

    private static readonly string[] Enzymes = ["EcoRI", "BamHI", "HindIII", "SmaI", "BsaI", "Sau3AI", "BglI", "SfiI", "ClaI"];

    private static Gen<string> DnaWithSites(int min, int max) =>
        from len in Gen.Choose(min, max)
        from bases in Gen.Elements(Bases).ArrayOf(len)
        from count in Gen.Choose(0, 4)
        from sites in Gen.Elements(PlantedSites).ArrayOf(count)
        from offsets in Gen.Choose(0, len - 1).ArrayOf(count)
        select Plant(bases, sites, offsets);

    private static string Plant(char[] bases, string[] sites, int[] offsets)
    {
        for (int i = 0; i < sites.Length; i++)
            if (offsets[i] + sites[i].Length <= bases.Length)
                sites[i].CopyTo(0, bases, offsets[i], sites[i].Length);
        return new string(bases);
    }

    private static int Mod(int value, int n) => ((value % n) + n) % n;

    private static string Rotate(string s, int r) => s[r..] + s[..r];

    #region Circular molecules have no origin

    private static string SiteSignature(RestrictionSite s, int shift, int n) =>
        $"{s.Enzyme.Name}|{s.IsForwardStrand}|{s.IsSecondCut}|{Mod(s.Position - shift, n)}|{Mod(s.CutPosition - shift, n)}";

    /// <summary>
    /// Rotating a circular sequence by r moves every site, on both strands and across the origin, by exactly -r
    /// (mod n): the multiset of (enzyme, strand, position, cut) is the original's shifted. Covers the origin-spanning
    /// sites, the Type IIS cut that wraps (BsaI) and the interrupted recognition sites (BglI, SfiI).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CircularSites_AreEquivariant_UnderRotation()
    {
        var arb = (from s in DnaWithSites(50, 160) from r in Gen.Choose(0, s.Length - 1) select (s, r)).ToArbitrary();
        return Prop.ForAll(arb, input =>
        {
            var (seq, r) = input;
            int n = seq.Length;
            var original = new DnaSequence(seq);
            var rotated = new DnaSequence(Rotate(seq, r));

            foreach (string enzyme in Enzymes)
            {
                var before = RestrictionAnalyzer.FindSites(original, enzyme, MoleculeTopology.Circular)
                    .Select(s => SiteSignature(s, 0, n)).OrderBy(x => x, StringComparer.Ordinal).ToList();
                // Original coordinate p appears at p - r in the rotated sequence, i.e. rotated = original shifted by r.
                var after = RestrictionAnalyzer.FindSites(rotated, enzyme, MoleculeTopology.Circular)
                    .Select(s => SiteSignature(s, -r, n)).OrderBy(x => x, StringComparer.Ordinal).ToList();

                if (!before.SequenceEqual(after))
                    return false.Label($"{enzyme}, n={n}, r={r}: [{string.Join("; ", before)}] vs [{string.Join("; ", after)}]");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// A circular digest is rotation-invariant: the multiset of fragment lengths is the same and sums to n, with exactly
    /// as many fragments as distinct cuts (one uncut circle when there are none).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CircularDigest_FragmentLengths_AreRotationInvariant_AndSumToN()
    {
        var arb = (from s in DnaWithSites(50, 160) from r in Gen.Choose(0, s.Length - 1) select (s, r)).ToArbitrary();
        return Prop.ForAll(arb, input =>
        {
            var (seq, r) = input;
            var original = new DnaSequence(seq);
            var rotated = new DnaSequence(Rotate(seq, r));

            foreach (string[] enzymes in new[] { new[] { "EcoRI" }, new[] { "Sau3AI" }, new[] { "BsaI", "HindIII" }, new[] { "BglI", "SmaI", "ClaI" } })
            {
                var a = RestrictionAnalyzer.Digest(original, MoleculeTopology.Circular, enzymes).Select(f => f.Length).OrderBy(x => x).ToList();
                var b = RestrictionAnalyzer.Digest(rotated, MoleculeTopology.Circular, enzymes).Select(f => f.Length).OrderBy(x => x).ToList();

                if (!a.SequenceEqual(b))
                    return false.Label($"{string.Join("+", enzymes)} n={seq.Length} r={r}: [{string.Join(",", a)}] vs [{string.Join(",", b)}]");
                if (a.Sum() != seq.Length)
                    return false.Label($"{string.Join("+", enzymes)}: fragments sum {a.Sum()} != {seq.Length}");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// The circular restriction map is rotation-invariant too: total cleavages, the unique cutters and the non-cutters
    /// (the "binds once but cuts once" rule reads the reverse-strand reading of a site that wraps past the origin).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CircularMap_Totals_AndUniqueCutters_AreRotationInvariant()
    {
        var arb = (from s in DnaWithSites(50, 160) from r in Gen.Choose(0, s.Length - 1) select (s, r)).ToArbitrary();
        return Prop.ForAll(arb, input =>
        {
            var (seq, r) = input;
            var a = RestrictionAnalyzer.CreateMap(new DnaSequence(seq), MoleculeTopology.Circular, Enzymes);
            var b = RestrictionAnalyzer.CreateMap(new DnaSequence(Rotate(seq, r)), MoleculeTopology.Circular, Enzymes);

            return (a.TotalSites == b.TotalSites).Label($"TotalSites {a.TotalSites} vs {b.TotalSites} (n={seq.Length}, r={r})")
                .And(a.UniqueCutters.OrderBy(x => x, StringComparer.Ordinal).SequenceEqual(b.UniqueCutters.OrderBy(x => x, StringComparer.Ordinal))
                    .Label($"UniqueCutters [{string.Join(",", a.UniqueCutters)}] vs [{string.Join(",", b.UniqueCutters)}] (n={seq.Length}, r={r})"))
                .And(a.NonCutters.OrderBy(x => x, StringComparer.Ordinal).SequenceEqual(b.NonCutters.OrderBy(x => x, StringComparer.Ordinal)).Label("NonCutters"));
        });
    }

    /// <summary>
    /// One EcoRI site that straddles the origin of a circular molecule is a unique cutter: the molecule is cut once, so the
    /// reverse-strand reading of the same wrapped window must not be counted as a second cleavage.
    /// </summary>
    [TestCase(0)]
    [TestCase(1)]
    [TestCase(3)]
    [TestCase(5)]
    public void SingleWrappedEcoRISite_IsAUniqueCutter(int wrappedBases)
    {
        // GAATTC placed so that its first wrappedBases bases sit at the end of the sequence and the rest at the start (0 = no wrap, the control).
        const string body = "ACCATAGCACATCTACCTGACCTTGCATCAGTACATCTGACAG";
        string site = "GAATTC";
        string seq = site[wrappedBases..] + body + site[..wrappedBases];

        var map = RestrictionAnalyzer.CreateMap(new DnaSequence(seq), MoleculeTopology.Circular, "EcoRI");

        map.UniqueCutters.Should().Contain("EcoRI");
        RestrictionAnalyzer.Digest(new DnaSequence(seq), MoleculeTopology.Circular, "EcoRI").Should().ContainSingle();
    }

    /// <summary>
    /// A site found in the linear molecule is also found in the circular one (same enzyme, strand and position):
    /// closing the molecule can add sites across the origin but never removes one.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CircularSites_AreASupersetOf_LinearSites()
    {
        return Prop.ForAll(DnaWithSites(50, 160).ToArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            foreach (string enzyme in Enzymes)
            {
                var circular = RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Circular)
                    .Select(s => (s.IsForwardStrand, s.IsSecondCut, s.Position)).ToHashSet();
                var missing = RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Linear)
                    .Select(s => (s.IsForwardStrand, s.IsSecondCut, s.Position)).Where(x => !circular.Contains(x)).ToList();

                if (missing.Count > 0)
                    return false.Label($"{enzyme}: linear site(s) {string.Join(", ", missing)} missing from the circular result");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region Partial digest combinatorics

    /// <summary>
    /// With k cuts a linear molecule has (k+1)(k+2)/2 partial-digest fragments - k+1-m of them with m missed cuts, each
    /// the exact substring it claims to be - and a circular one k^2 + 1. The m = 0 fragments are the complete digest.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PartialDigest_Counts_FollowTheCombinatorics()
    {
        return Prop.ForAll(DnaWithSites(50, 160).ToArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            foreach (string enzyme in new[] { "EcoRI", "BamHI", "HindIII", "Sau3AI" })
            {
                var complete = RestrictionAnalyzer.Digest(dna, MoleculeTopology.Linear, enzyme).ToList();
                int k = complete.Count - 1;
                var partial = RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Linear, enzyme).ToList();

                if (partial.Count != (k + 1) * (k + 2) / 2)
                    return false.Label($"{enzyme} linear k={k}: {partial.Count} fragments, expected {(k + 1) * (k + 2) / 2}");
                for (int m = 0; m <= k; m++)
                    if (partial.Count(f => f.MissedCuts == m) != k + 1 - m)
                        return false.Label($"{enzyme} linear k={k}: {partial.Count(f => f.MissedCuts == m)} fragments with {m} missed cuts, expected {k + 1 - m}");
                if (partial.Any(f => seq.Substring(f.Fragment.StartPosition, f.Fragment.Length) != f.Fragment.Sequence))
                    return false.Label($"{enzyme} linear: a fragment is not the substring it reports");
                var cutFree = partial.Where(f => f.MissedCuts == 0).Select(f => f.Fragment.Sequence).OrderBy(x => x, StringComparer.Ordinal);
                if (!cutFree.SequenceEqual(complete.Select(f => f.Sequence).OrderBy(x => x, StringComparer.Ordinal)))
                    return false.Label($"{enzyme} linear: the zero-missed-cut fragments are not the complete digest");

                // A palindromic site is reported once per strand, so the number of cut events is read off the complete circular
                // digest: k >= 1 cuts give k fragments, none gives the single uncut circle.
                bool anySite = RestrictionAnalyzer.FindSites(dna, enzyme, MoleculeTopology.Circular).Any();
                int cuts = anySite ? RestrictionAnalyzer.Digest(dna, MoleculeTopology.Circular, enzyme).Count() : 0;
                var circular = RestrictionAnalyzer.PartialDigest(dna, MoleculeTopology.Circular, enzyme).ToList();
                if (circular.Count != cuts * cuts + 1)
                    return false.Label($"{enzyme} circular k={cuts}: {circular.Count} fragments, expected {cuts * cuts + 1}");
            }

            return true.ToProperty();
        });
    }

    /// <summary>The missed-cut limit keeps exactly the fragments with at most that many missed cuts.</summary>
    [FsCheck.NUnit.Property]
    public Property PartialDigest_MissedCutLimit_IsAFilter()
    {
        var arb = (from s in DnaWithSites(50, 160) from m in Gen.Choose(0, 5) select (s, m)).ToArbitrary();
        return Prop.ForAll(arb, input =>
        {
            var (seq, limit) = input;
            var dna = new DnaSequence(seq);
            foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
            {
                var all = RestrictionAnalyzer.PartialDigest(dna, topology, "Sau3AI").ToList();
                var limited = RestrictionAnalyzer.PartialDigest(dna, topology, limit, "Sau3AI").ToList();
                var expected = all.Where(f => f.MissedCuts <= limit).Select(f => (f.Fragment.StartPosition, f.Fragment.Length, f.MissedCuts)).OrderBy(x => x).ToList();

                if (!limited.Select(f => (f.Fragment.StartPosition, f.Fragment.Length, f.MissedCuts)).OrderBy(x => x).SequenceEqual(expected))
                    return false.Label($"{topology} limit {limit}: {limited.Count} fragments vs {expected.Count} expected");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region IUPAC symbols in the sequence

    private const string AllIupac = "ACGTRYSWKMBDHVN";

    /// <summary>
    /// Recognition-pattern N matches every one of the 15 IUPAC symbols in the sequence (and no other character): BglI
    /// (GCCNNNNNGGC) is found with any symbol at any of the five N positions, and lost when that position holds a
    /// non-IUPAC character.
    /// </summary>
    [Test]
    public void NRegion_AcceptsEveryIupacSymbol_AndOnlyThose()
    {
        const string flank = "ATATATAT";
        for (int i = 0; i < 5; i++)
        {
            foreach (char symbol in AllIupac)
            {
                string site = "GCC" + new string('A', i) + symbol + new string('A', 4 - i) + "GGC";
                // The interrupted site is non-palindromic in its N run, so it is reported once per strand at one position.
                RestrictionAnalyzer.FindSites(flank + site + flank, "BglI").Select(s => s.Position).Distinct().Should().ContainSingle(
                    $"BglI with '{symbol}' at N position {i}");
            }

            foreach (char junk in "X-*.?1")
            {
                string site = "GCC" + new string('A', i) + junk + new string('A', 4 - i) + "GGC";
                RestrictionAnalyzer.FindSites(flank + site + flank, "BglI").Should().BeEmpty(
                    $"'{junk}' is not an IUPAC symbol");
            }
        }
    }

    /// <summary>
    /// A concrete recognition code matches only the concrete base: an ambiguity symbol in the sequence at any of EcoRI's six
    /// positions removes the site (the symbol is not narrowed to a base the pattern would accept).
    /// </summary>
    [Test]
    public void ConcreteSiteCode_IsBlockedByEveryAmbiguitySymbol_AtEveryPosition()
    {
        const string flank = "ATATATAT";
        const string site = "GAATTC";
        for (int i = 0; i < site.Length; i++)
            foreach (char symbol in AllIupac.Where(c => c is not ('A' or 'C' or 'G' or 'T')))
            {
                string mutated = flank + site[..i] + symbol + site[(i + 1)..] + flank;

                RestrictionAnalyzer.FindSites(mutated, "EcoRI").Should().BeEmpty($"'{symbol}' at position {i} of GAATTC");
            }
    }

    /// <summary>
    /// Raw-string overloads accept either case: lower-casing the input (IUPAC symbols included) changes neither the
    /// sites nor the reported recognised sequence, which is upper case.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RawFindSites_IsCaseInsensitive_AndReportsUpperCase()
    {
        var iupac = Gen.Elements(AllIupac.ToCharArray()).ArrayOf().Where(a => a.Length is >= 20 and <= 120).Select(a => new string(a));
        var arb = Gen.OneOf(iupac, DnaWithSites(30, 120)).ToArbitrary();
        return Prop.ForAll(arb, seq =>
        {
            foreach (string enzyme in Enzymes)
            {
                var upper = RestrictionAnalyzer.FindSites(seq, enzyme).Select(s => (s.Position, s.IsForwardStrand, s.CutPosition, s.RecognizedSequence)).ToList();
                var lower = RestrictionAnalyzer.FindSites(seq.ToLowerInvariant(), enzyme).Select(s => (s.Position, s.IsForwardStrand, s.CutPosition, s.RecognizedSequence)).ToList();

                if (!upper.SequenceEqual(lower))
                    return false.Label($"{enzyme}: upper {upper.Count} sites vs lower {lower.Count}");
                if (lower.Any(s => s.RecognizedSequence.Any(char.IsLower)))
                    return false.Label($"{enzyme}: a recognised sequence is not upper case");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// The raw-string overloads that take a methylase set or a topology share the case rule: lower-casing the input changes
    /// neither the sites nor the digest fragments, and the fragment sequences come back upper-case.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RawMethylationAndTopologyOverloads_AreCaseInsensitive_AndReturnUpperCase()
    {
        return Prop.ForAll(DnaWithSites(40, 140).ToArbitrary(), seq =>
        {
            string lower = seq.ToLowerInvariant();
            foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
            {
                var sitesUpper = RestrictionAnalyzer.FindSites(seq, topology, MethylaseSets.Dam, "Sau3AI", "EcoRI")
                    .Select(s => (s.Enzyme.Name, s.Position, s.IsForwardStrand, s.CutPosition, s.RecognizedSequence)).ToList();
                var sitesLower = RestrictionAnalyzer.FindSites(lower, topology, MethylaseSets.Dam, "Sau3AI", "EcoRI")
                    .Select(s => (s.Enzyme.Name, s.Position, s.IsForwardStrand, s.CutPosition, s.RecognizedSequence)).ToList();
                if (!sitesUpper.SequenceEqual(sitesLower))
                    return false.Label($"{topology}: FindSites differs between cases");

                var fragsUpper = RestrictionAnalyzer.Digest(seq, topology, MethylaseSets.Dam, "Sau3AI", "EcoRI").Select(f => (f.StartPosition, f.Sequence)).ToList();
                var fragsLower = RestrictionAnalyzer.Digest(lower, topology, MethylaseSets.Dam, "Sau3AI", "EcoRI").Select(f => (f.StartPosition, f.Sequence)).ToList();
                if (!fragsUpper.SequenceEqual(fragsLower))
                    return false.Label($"{topology}: Digest differs between cases");
                if (fragsLower.Any(f => f.Sequence.Any(char.IsLower)))
                    return false.Label($"{topology}: a digest fragment is not upper case");

                var map = RestrictionAnalyzer.CreateMap(lower, topology, MethylaseSets.Dam, "Sau3AI", "EcoRI");
                var mapUpper = RestrictionAnalyzer.CreateMap(seq, topology, MethylaseSets.Dam, "Sau3AI", "EcoRI");
                if (map.TotalSites != mapUpper.TotalSites || !map.UniqueCutters.SequenceEqual(mapUpper.UniqueCutters))
                    return false.Label($"{topology}: CreateMap differs between cases");
            }

            return true.ToProperty();
        });
    }

    #endregion

    #region Methylation only ever blocks

    private static IEnumerable<(string Name, IReadOnlyList<Methylase> Set)> MethylaseCases() =>
    [
        ("Dam", MethylaseSets.Dam),
        ("Dcm", MethylaseSets.Dcm),
        ("DamDcm", MethylaseSets.DamDcm),
        ("Emboss", MethylaseSets.Emboss),
    ];

    /// <summary>
    /// For every methylase set, topology and enzyme, the sites found under methylation are a subset of those found without
    /// it (same strand, position and cut), and an empty set blocks nothing. Methylation cannot create a site.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Methylation_OnlyRemovesSites()
    {
        return Prop.ForAll(DnaWithSites(40, 140).ToArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
                foreach (string enzyme in Enzymes)
                {
                    var plain = RestrictionAnalyzer.FindSites(dna, enzyme, topology)
                        .Select(s => (s.IsForwardStrand, s.IsSecondCut, s.Position, s.CutPosition)).ToHashSet();
                    var none = RestrictionAnalyzer.FindSites(dna, enzyme, topology, Array.Empty<Methylase>())
                        .Select(s => (s.IsForwardStrand, s.IsSecondCut, s.Position, s.CutPosition)).ToHashSet();
                    if (!plain.SetEquals(none))
                        return false.Label($"{enzyme} {topology}: an empty methylase set changed the sites");

                    foreach (var (name, set) in MethylaseCases())
                    {
                        var blocked = RestrictionAnalyzer.FindSites(dna, enzyme, topology, set)
                            .Select(s => (s.IsForwardStrand, s.IsSecondCut, s.Position, s.CutPosition)).ToList();
                        if (blocked.Any(s => !plain.Contains(s)))
                            return false.Label($"{enzyme} {topology} {name}: methylation created a site");
                    }
                }

            return true.ToProperty();
        });
    }

    /// <summary>The methylated digest still tiles the molecule: fragment lengths sum to n for every set and topology.</summary>
    [FsCheck.NUnit.Property]
    public Property MethylatedDigest_StillTilesTheMolecule()
    {
        return Prop.ForAll(DnaWithSites(40, 140).ToArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
                foreach (var (name, set) in MethylaseCases())
                {
                    int total = RestrictionAnalyzer.Digest(dna, topology, set, "Sau3AI", "BamHI", "EcoRI").Sum(f => f.Length);
                    if (total != seq.Length)
                        return false.Label($"{topology} {name}: fragments sum {total} != {seq.Length}");
                }

            return true.ToProperty();
        });
    }

    #endregion
}
