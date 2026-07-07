namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Variants area — VARIANT-SNP-001 (SNP Detection).
/// The unit under test is the single-nucleotide substitution detector:
/// <see cref="VariantCaller.FindSnpsDirect"/> (positional Hamming-mismatch
/// enumeration over the common prefix — the canonical entry point) and
/// <see cref="VariantCaller.FindSnps"/> (align-then-filter-to-SNP delegate over
/// <c>CallVariants</c>), together with the SNP classifiers
/// <see cref="VariantCaller.ClassifyMutation"/> (transition / transversion) and
/// <see cref="VariantCaller.CalculateTiTvRatio"/>; implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference / ArgumentOutOfRange / DivideByZero). Every
/// input must resolve to EITHER a well-defined, theory-correct value OR a
/// *documented, intentional* outcome (here, an ArgumentNullException for null
/// FindSnps inputs; an empty result for empty/equal FindSnpsDirect inputs). —
/// docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: VARIANT-SNP-001 — SNP Detection (Variants)
/// Checklist: docs/checklists/03_FUZZING.md, row 189.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// MAPPING of the generic checklist BE targets ("ref=alt, multi-allelic, zero
/// depth") onto THIS unit's documented contract
/// (docs/algorithms/Variants/SNP_Detection.md). This unit is an exact positional
/// / alignment substitution scan, not a pileup caller, so the pileup-flavoured
/// targets are translated to their substitution-scan equivalents:
///   • "ref=alt" → a position where reference[i] == query[i]. A SNP is by
///       definition a substitution (REF≠ALT); an equal column is a MATCH, never a
///       variant ⇒ NO SNP must be emitted there — no false positive (INV-03,
///       §2.2, §6.1).
///   • "multi-allelic" → a single locus carrying MORE than one alternate base.
///       This pairwise detector represents one query allele per position; the
///       documented multi-allelic representation is therefore "one Variant per
///       differing position", and across several queries each alt is reported
///       independently with no dropped / merged / corrupted allele. Each emitted
///       record is still a clean single-base REF→single-base ALT (INV-02..INV-04).
///   • "zero depth" → a locus with NO read support. In the substitution-scan model
///       this is "nothing to compare": empty / null inputs and the all-match
///       (Hamming-distance-0) case. Result is the documented EMPTY SNP set with no
///       crash, and the derived Ti/Tv ratio must NOT divide by zero — it returns 0
///       when there are no transversions (CalculateTiTvRatio, §Mutation Classification).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (docs/algorithms/Variants/SNP_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • A SNP is a single-base substitution: REF≠ALT, both length 1, at one
///       position (§2.1, §2.2).                                          (INV-03)
///   • Positional model: for i ∈ [0, min(|r|,|q|)), r[i]≠q[i] ⇒ a SNP with
///       Position==i, REF==r[i], ALT==q[i], QueryPosition==i.            (INV-04)
///   • Every emitted variant has Type==VariantType.SNP.                  (INV-02)
///   • Identical equal-length inputs ⇒ zero SNPs (Hamming distance 0).   (INV-01)
///   • SNP count over equal-length inputs == Hamming distance.           (INV-05)
///   • FindSnpsDirect compares only the common prefix min(|r|,|q|).      (INV-06)
///   • FindSnps: null reference/query ⇒ ArgumentNullException; empty ⇒ empty.
///       FindSnpsDirect: null/empty input ⇒ empty.                       (§3.3, §6.1)
///   • Ti/Tv classification: purine↔purine or pyrimidine↔pyrimidine ⇒ Transition;
///       purine↔pyrimidine ⇒ Transversion; case-insensitive (§Mutation
///       Classification, refs [2][4]).
///   • Ti/Tv ratio over a SNP set: transitions / transversions, returning 0 when
///       there are no transversions (no DivideByZero).
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class VariantSnpFuzzTests
{
    // Purine / pyrimidine sets for the independent (oracle) Ti/Tv computation.
    private static bool IsPurine(char b) => char.ToUpperInvariant(b) is 'A' or 'G';
    private static bool IsPyrimidine(char b) => char.ToUpperInvariant(b) is 'C' or 'T';

    // ── Well-formed-SNP assertion helper ─────────────────────────────────────
    // Pins the documented per-SNP contract on EVERY emitted record, no matter how
    // degenerate the input. This is what stops a fuzz test from rubber-stamping
    // nonsense output green — a phantom ref==alt "SNP", an out-of-range position, a
    // multi-base allele, or a mis-classified transition/transversion would all
    // fail here:
    //   • Type is exactly SNP (never an indel) (INV-02).
    //   • REF and ALT are each a single, non-gap base and REF≠ALT — a real
    //     substitution, never a length-0 / ref==alt record (INV-03).
    //   • Position and QueryPosition are finite, non-negative coordinates in range
    //     (INV-04, INV-06).
    //   • When REF and ALT are both canonical A/C/G/T, the reported
    //     ClassifyMutation matches the purine/pyrimidine oracle.
    private static void AssertWellFormedSnp(Variant v, int refLen, int queryLen)
    {
        v.Type.Should().Be(VariantType.SNP, "this unit reports only substitutions (INV-02)");

        v.ReferenceAllele.Should().HaveLength(1, "a SNP REF is a single base (INV-03)");
        v.AlternateAllele.Should().HaveLength(1, "a SNP ALT is a single base (INV-03)");
        v.ReferenceAllele.Should().NotBe("-", "a SNP is a substitution, not a gap (INV-03)");
        v.AlternateAllele.Should().NotBe("-", "a SNP is a substitution, not a gap (INV-03)");
        v.ReferenceAllele.Should().NotBe(v.AlternateAllele,
            "a SNP requires REF≠ALT — an equal column is a match, never a SNP (INV-03)");

        v.Position.Should().BeGreaterThanOrEqualTo(0, "Position is a 0-based reference coordinate (INV-04)");
        v.Position.Should().BeLessThan(Math.Max(refLen, 1), "Position ∈ [0, reference.Length) for a substituted base (INV-04/06)");
        v.QueryPosition.Should().BeGreaterThanOrEqualTo(0, "QueryPosition is a 0-based query coordinate");
        v.QueryPosition.Should().BeLessThan(Math.Max(queryLen, 1), "QueryPosition addresses a real query base");

        // Ti/Tv oracle for canonical bases.
        char r = v.ReferenceAllele[0];
        char a = v.AlternateAllele[0];
        bool canonical = (IsPurine(r) || IsPyrimidine(r)) && (IsPurine(a) || IsPyrimidine(a));
        if (canonical)
        {
            var expected = IsPurine(r) == IsPurine(a) ? MutationType.Transition : MutationType.Transversion;
            VariantCaller.ClassifyMutation(v).Should().Be(expected,
                "Ti/Tv must follow purine/pyrimidine class (refs [2][4])");
        }
    }

    private static void AssertAllWellFormedSnps(IReadOnlyList<Variant> snps, int refLen, int queryLen)
    {
        foreach (var v in snps)
            AssertWellFormedSnp(v, refLen, queryLen);
    }

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-SNP-001 — SNP Detection (positive sanity)
    // ═════════════════════════════════════════════════════════════════════════

    // ── POSITIVE sanity: a ref/query pair with a hand-computed SNP set.
    //    ref ATGCATGC vs query ATTCATTC ⇒ two substitutions:
    //      index 2: G→T  (purine→pyrimidine ⇒ Transversion)
    //      index 6: G→T  (purine→pyrimidine ⇒ Transversion)
    //    Every other column is a match. (docs §7.1, §2.2) ──
    [Test]
    public void FindSnpsDirect_KnownSubstitutions_HandComputedCalls()
    {
        var snps = VariantCaller.FindSnpsDirect("ATGCATGC", "ATTCATTC").ToList();

        snps.Should().HaveCount(2, "exactly two columns differ ⇒ Hamming distance 2 (INV-05)");

        snps[0].Position.Should().Be(2);
        snps[0].ReferenceAllele.Should().Be("G");
        snps[0].AlternateAllele.Should().Be("T");
        snps[0].QueryPosition.Should().Be(2, "FindSnpsDirect indexes query == reference (INV-04)");
        VariantCaller.ClassifyMutation(snps[0]).Should().Be(MutationType.Transversion, "G→T is purine→pyrimidine");

        snps[1].Position.Should().Be(6);
        snps[1].ReferenceAllele.Should().Be("G");
        snps[1].AlternateAllele.Should().Be("T");

        AssertAllWellFormedSnps(snps, refLen: 8, queryLen: 8);
    }

    // A documented worked example from §7.1 verbatim: single G→T at index 2.
    [Test]
    public void FindSnpsDirect_DocWorkedExample_SingleSnp()
    {
        var snps = VariantCaller.FindSnpsDirect("ATGC", "ATTC").ToList();

        snps.Should().ContainSingle();
        var v = snps[0];
        v.Position.Should().Be(2);
        v.ReferenceAllele.Should().Be("G");
        v.AlternateAllele.Should().Be("T");
        v.Type.Should().Be(VariantType.SNP);
    }

    // A transition (A→G, purine→purine) is classified Transition; case-insensitive.
    [Test]
    public void ClassifyMutation_TransitionAndTransversion_HandComputed()
    {
        // A→G transition.
        var ti = VariantCaller.FindSnpsDirect("AAAA", "AGAA").Single();
        ti.ReferenceAllele.Should().Be("A");
        ti.AlternateAllele.Should().Be("G");
        VariantCaller.ClassifyMutation(ti).Should().Be(MutationType.Transition, "A↔G is purine↔purine");

        // C→A transversion, lowercase inputs ⇒ classification must be case-insensitive.
        var tv = VariantCaller.FindSnpsDirect("cccc", "caca").ToList();
        tv.Should().HaveCount(2);
        VariantCaller.ClassifyMutation(tv[0]).Should().Be(MutationType.Transversion, "c↔a is pyrimidine↔purine (case-insensitive)");
    }

    // FindSnps (the align-then-filter delegate) reports the substitution and drops
    // indels. Repeat-free inputs ⇒ unique alignment ⇒ deterministic SNP position.
    [Test]
    public void FindSnps_AlignmentBased_DropsIndelsKeepsSnp()
    {
        // docs §7.1: ref ATGCATGC vs query ATGAATGC ⇒ one SNP at pos 3 (C→A).
        var snps = VariantCaller.FindSnps(new DnaSequence("ATGCATGC"), new DnaSequence("ATGAATGC")).ToList();

        snps.Should().ContainSingle("one substitution, no indels");
        snps.Should().OnlyContain(v => v.Type == VariantType.SNP, "FindSnps filters to SNP only (INV-02)");
        snps[0].ReferenceAllele.Should().Be("C");
        snps[0].AlternateAllele.Should().Be("A");
        snps[0].Position.Should().Be(3);
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-SNP-001 — BE: "ref=alt" (equal column ⇒ NOT a SNP)
    // ═════════════════════════════════════════════════════════════════════════

    // Identical equal-length inputs: every column is ref==alt ⇒ zero SNPs. A false
    // SNP on an equal column would be a substitution with REF==ALT (INV-03 break).
    [Test]
    public void FindSnpsDirect_IdenticalSequences_NoSnp()
    {
        VariantCaller.FindSnpsDirect("GATTACA", "GATTACA").Should().BeEmpty(
            "identical sequences differ nowhere ⇒ Hamming distance 0 (INV-01)");
        VariantCaller.FindSnps(new DnaSequence("GATTACA"), new DnaSequence("GATTACA")).Should().BeEmpty(
            "no substitution between identical sequences (INV-01)");
    }

    // A single all-same base run: ref==alt at every position, no SNP, no false call.
    [Test]
    public void FindSnpsDirect_AllSameBase_NoSnp()
    {
        VariantCaller.FindSnpsDirect("AAAAAAAA", "AAAAAAAA").Should().BeEmpty("ref==alt everywhere ⇒ no SNP (INV-03)");
    }

    // Mixed: some columns equal, some different. The equal columns must contribute
    // NO record; only the differing columns are SNPs, and never with REF==ALT.
    [Test]
    public void FindSnpsDirect_MixedEqualAndDiffering_OnlyDifferingReported()
    {
        // index 0,2,4 equal; index 1,3 differ.
        var snps = VariantCaller.FindSnpsDirect("ACGTA", "ATGGA").ToList();

        snps.Should().HaveCount(2, "only the two differing columns are SNPs (equal columns are matches)");
        snps.Select(s => s.Position).Should().Equal(new[] { 1, 3 }, "the equal columns 0,2,4 contribute no record (INV-03)");
        snps.Should().OnlyContain(v => v.ReferenceAllele != v.AlternateAllele,
            "no emitted SNP has REF==ALT — that would be a phantom variant (INV-03)");
        AssertAllWellFormedSnps(snps, refLen: 5, queryLen: 5);
    }

    // Fuzz: query equals reference for random sequences ⇒ ALWAYS zero SNPs, never a
    // false ref==alt positive, never a crash.
    [Test]
    [CancelAfter(30_000)]
    public void FindSnpsDirect_QueryEqualsReference_NeverFalseSnp([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";

        for (int t = 0; t < 40; t++)
        {
            string s = RandomDna(rng, bases, 0, 30);

            List<Variant> direct = null!;
            var act = () => direct = VariantCaller.FindSnpsDirect(s, s).ToList();
            act.Should().NotThrow("identical inputs are valid — must never crash the scan");
            direct.Should().BeEmpty("ref==alt at every position ⇒ no SNP (INV-01/03)");

            // And via the alignment-based delegate (DnaSequence cannot be empty).
            if (s.Length > 0)
            {
                var ds = new DnaSequence(s);
                VariantCaller.FindSnps(ds, ds).Should().BeEmpty("identical sequences ⇒ no SNP (INV-01)");
            }
        }
    }

    // Reinforce the no-false-positive guarantee generally: for any random pair,
    // EVERY emitted SNP has REF≠ALT (a real substitution), never an equal column.
    [Test]
    [CancelAfter(30_000)]
    public void FindSnpsDirect_RandomPairs_NoSnpHasRefEqualsAlt([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGTN";

        for (int t = 0; t < 40; t++)
        {
            string r = RandomDna(rng, bases, 0, 28);
            string q = RandomDna(rng, bases, 0, 28);

            List<Variant> snps = null!;
            var act = () => snps = VariantCaller.FindSnpsDirect(r, q).ToList();
            act.Should().NotThrow("any A/C/G/T/N pair is valid input — must never crash");

            int n = Math.Min(r.Length, q.Length);
            // Independent oracle: the SNP set is exactly the mismatch set over the prefix.
            var expected = Enumerable.Range(0, n).Where(i => r[i] != q[i]).ToList();
            snps.Select(s => s.Position).Should().Equal(expected, "SNP set == Hamming mismatch set over the common prefix (INV-04/05/06)");
            snps.Should().OnlyContain(v => v.ReferenceAllele != v.AlternateAllele,
                "no SNP has REF==ALT (INV-03)");
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-SNP-001 — BE: "multi-allelic" (one Variant per differing position)
    // ═════════════════════════════════════════════════════════════════════════

    // A pairwise scan represents one query allele per locus. The documented
    // multi-allelic representation is "one Variant per differing position": several
    // distinct substitutions in one pair must each be reported independently with
    // its own REF/ALT — no dropped, merged, or corrupted allele.
    [Test]
    public void FindSnpsDirect_MultipleDistinctSubstitutions_EachReportedSeparately()
    {
        // ref ACGT vs query GTAC ⇒ all four columns differ with four distinct ALTs.
        var snps = VariantCaller.FindSnpsDirect("ACGT", "GTAC").ToList();

        snps.Should().HaveCount(4, "four differing columns ⇒ four independent SNPs (no merge/drop)");
        var expectedCalls = new[]
        {
            (Position: 0, Ref: "A", Alt: "G"),
            (Position: 1, Ref: "C", Alt: "T"),
            (Position: 2, Ref: "G", Alt: "A"),
            (Position: 3, Ref: "T", Alt: "C"),
        };
        snps.Select(s => (Position: s.Position, Ref: s.ReferenceAllele, Alt: s.AlternateAllele))
            .Should().Equal(expectedCalls, "each locus keeps its own REF→ALT pair (no allele corruption)");
        AssertAllWellFormedSnps(snps, refLen: 4, queryLen: 4);
    }

    // Multi-allelic ACROSS queries: the SAME reference locus carries DIFFERENT alt
    // bases in different queries. Each query yields the correct single alt at that
    // locus — the per-locus alt is never confused between samples.
    [Test]
    public void FindSnpsDirect_SameLocusDifferentAlts_EachAltDistinct()
    {
        const string reference = "AAAA";

        // Three "samples" each substituting a different base at locus 1.
        var alts = new[] { ("ACAA", 'C'), ("AGAA", 'G'), ("ATAA", 'T') };
        foreach (var (query, expectedAlt) in alts)
        {
            var snps = VariantCaller.FindSnpsDirect(reference, query).ToList();
            snps.Should().ContainSingle("one differing locus per sample");
            snps[0].Position.Should().Be(1);
            snps[0].ReferenceAllele.Should().Be("A");
            snps[0].AlternateAllele.Should().Be(expectedAlt.ToString(),
                "the alt at the shared locus is this sample's base, not another sample's (no merge)");
        }
    }

    // A locus carrying every non-reference alt, gathered from multiple queries, must
    // produce the full alt set with no dropped allele — the multi-allelic union.
    [Test]
    public void FindSnpsDirect_MultiAllelicUnionAcrossQueries_AllAltsPresent()
    {
        const string reference = "C";
        var observedAlts = new[] { "A", "G", "T" }
            .Select(b => VariantCaller.FindSnpsDirect(reference, b).Single().AlternateAllele)
            .ToHashSet();

        observedAlts.Should().BeEquivalentTo(new[] { "A", "G", "T" },
            "all three non-reference alleles at the locus are reported, none dropped or merged");
    }

    // Fuzz: many random pairs over the same reference. Every emitted SNP is a clean
    // single-base REF→single-base ALT (INV-03), and the count equals the number of
    // differing prefix columns — no allele is dropped or merged under fuzzing.
    [Test]
    [CancelAfter(30_000)]
    public void FindSnpsDirect_RandomMultiAllelic_NoAlleleDroppedOrMerged([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";
        string reference = RandomDna(rng, bases, 4, 20);

        for (int t = 0; t < 30; t++)
        {
            // Equal-length query so every column is a substitution candidate.
            var qChars = reference.ToCharArray();
            for (int i = 0; i < qChars.Length; i++)
                if (rng.Next(2) == 0) qChars[i] = bases[rng.Next(bases.Length)];
            string query = new string(qChars);

            List<Variant> snps = null!;
            var act = () => snps = VariantCaller.FindSnpsDirect(reference, query).ToList();
            act.Should().NotThrow();

            int expected = Enumerable.Range(0, reference.Length).Count(i => reference[i] != query[i]);
            snps.Should().HaveCount(expected, "every differing locus is an independent SNP — none dropped/merged");
            AssertAllWellFormedSnps(snps, refLen: reference.Length, queryLen: query.Length);
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-SNP-001 — BE: "zero depth" (no support ⇒ no call, no DivideByZero)
    // ═════════════════════════════════════════════════════════════════════════

    // Empty / null inputs are the substitution-scan analogue of "zero depth": no
    // bases to compare ⇒ documented empty SNP set (FindSnpsDirect) or
    // ArgumentNullException (FindSnps), never a crash.
    [Test]
    public void FindSnpsDirect_EmptyOrNullInput_EmptyNoCrash()
    {
        VariantCaller.FindSnpsDirect("", "").Should().BeEmpty("nothing to compare ⇒ no SNP (§6.1)");
        VariantCaller.FindSnpsDirect("", "ACGT").Should().BeEmpty("empty reference ⇒ no SNP (§6.1)");
        VariantCaller.FindSnpsDirect("ACGT", "").Should().BeEmpty("empty query ⇒ no SNP (§6.1)");
        VariantCaller.FindSnpsDirect(null!, "ACGT").Should().BeEmpty("null reference ⇒ guarded empty (§3.3)");
        VariantCaller.FindSnpsDirect("ACGT", null!).Should().BeEmpty("null query ⇒ guarded empty (§3.3)");
    }

    [Test]
    public void FindSnps_NullReferenceOrQuery_ThrowsArgumentNullException()
    {
        var nullRef = () => VariantCaller.FindSnps(null!, new DnaSequence("ACGT")).ToList();
        var nullQry = () => VariantCaller.FindSnps(new DnaSequence("ACGT"), null!).ToList();
        nullRef.Should().Throw<ArgumentNullException>("null reference is invalid input (§3.3)");
        nullQry.Should().Throw<ArgumentNullException>("null query is invalid input (§3.3)");
    }

    // The "zero depth" DivideByZero risk lives in CalculateTiTvRatio: ratio =
    // transitions / transversions. With zero transversions (or an empty SNP set) it
    // must NOT throw — it returns 0 by the documented contract.
    [Test]
    public void CalculateTiTvRatio_EmptyOrNoTransversions_ReturnsZeroNoDivideByZero()
    {
        // Empty SNP set ("zero depth").
        var emptyAct = () => VariantCaller.CalculateTiTvRatio(Enumerable.Empty<Variant>());
        emptyAct.Should().NotThrow<DivideByZeroException>();
        emptyAct().Should().Be(0, "no SNPs ⇒ ratio 0, never DivideByZero");

        // Only transitions, zero transversions ⇒ still guarded.
        var transitionsOnly = VariantCaller.FindSnpsDirect("AAAA", "GGGG").ToList(); // all A→G transitions
        transitionsOnly.Should().OnlyContain(v => VariantCaller.ClassifyMutation(v) == MutationType.Transition);
        var ratioAct = () => VariantCaller.CalculateTiTvRatio(transitionsOnly);
        ratioAct.Should().NotThrow<DivideByZeroException>("zero transversions must not divide by zero");
        ratioAct().Should().Be(0, "guarded: transversions==0 ⇒ ratio 0 (§Mutation Classification)");
    }

    // CalculateStatistics over identical / empty-equivalent inputs ("zero depth"):
    // zero SNPs, zero Ti/Tv, and the density computation must not divide by a zero
    // reference length.
    [Test]
    public void CalculateStatistics_IdenticalSequences_ZeroSnpsNoDivideByZero()
    {
        var seq = new DnaSequence("ACGTACGT");
        VariantStatistics stats = default;
        var act = () => stats = VariantCaller.CalculateStatistics(seq, seq);
        act.Should().NotThrow<DivideByZeroException>("identical inputs must not divide by zero");
        stats.Snps.Should().Be(0, "identical sequences carry no SNP (INV-01)");
        stats.TiTvRatio.Should().Be(0, "no SNPs ⇒ Ti/Tv 0, no DivideByZero");
    }

    // Fuzz the zero-depth / Ti/Tv path: random SNP sets (possibly empty, possibly
    // all-transition) must never throw DivideByZero and must always yield a finite,
    // non-negative ratio matching an independent oracle.
    [Test]
    [CancelAfter(30_000)]
    public void CalculateTiTvRatio_RandomSnpSets_FiniteNoDivideByZero([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";

        for (int t = 0; t < 40; t++)
        {
            string r = RandomDna(rng, bases, 0, 24);
            string q = RandomDna(rng, bases, 0, 24);
            var snps = VariantCaller.FindSnpsDirect(r, q).ToList();

            double ratio = 0;
            var act = () => ratio = VariantCaller.CalculateTiTvRatio(snps);
            act.Should().NotThrow<DivideByZeroException>("Ti/Tv must never divide by zero (zero-depth guard)");

            ratio.Should().BeGreaterThanOrEqualTo(0, "the ratio is non-negative");
            double.IsNaN(ratio).Should().BeFalse("the ratio is finite");
            double.IsInfinity(ratio).Should().BeFalse("the ratio is finite");

            // Independent oracle.
            int ti = snps.Count(v => IsPurine(v.ReferenceAllele[0]) == IsPurine(v.AlternateAllele[0]));
            int tv = snps.Count - ti;
            double expected = tv > 0 ? (double)ti / tv : 0;
            ratio.Should().BeApproximately(expected, 1e-9, "ratio matches transitions/transversions oracle");
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-SNP-001 — BE: unequal-length / common-prefix boundary
    // ═════════════════════════════════════════════════════════════════════════

    // FindSnpsDirect compares only the common prefix min(|r|,|q|); the trailing
    // region is NOT examined (it is indel territory, VARIANT-INDEL-001). A length
    // difference alone produces no SNP, and no IndexOutOfRange on the shorter side.
    [Test]
    public void FindSnpsDirect_UnequalLengths_ComparesCommonPrefixOnly()
    {
        // Common prefix "ATG" matches; query is shorter ⇒ no SNP, no crash.
        VariantCaller.FindSnpsDirect("ATGCATGC", "ATG").Should().BeEmpty(
            "matching common prefix ⇒ no SNP; trailing bases not scanned (INV-06)");

        // One mismatch inside the common prefix, with a length difference after it.
        var snps = VariantCaller.FindSnpsDirect("ATTCATGC", "ATG").ToList();
        snps.Should().ContainSingle("only the in-prefix mismatch at index 2 is a SNP (INV-06)");
        snps[0].Position.Should().Be(2);
    }

    [Test]
    [CancelAfter(30_000)]
    public void FindSnpsDirect_RandomUnequalLengths_NeverIndexError([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGTN";

        for (int t = 0; t < 40; t++)
        {
            string r = RandomDna(rng, bases, 0, 30);
            string q = RandomDna(rng, bases, 0, 30);

            List<Variant> snps = null!;
            var act = () => snps = VariantCaller.FindSnpsDirect(r, q).ToList();
            act.Should().NotThrow("unequal-length inputs must never throw IndexOutOfRange (INV-06)");

            int n = Math.Min(r.Length, q.Length);
            snps.Should().OnlyContain(v => v.Position < Math.Max(n, 1),
                "no SNP is reported past the common prefix (INV-06)");
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region Helpers
    // ═════════════════════════════════════════════════════════════════════════

    private static string RandomDna(Random rng, string bases, int minLen, int maxLen)
    {
        int len = rng.Next(minLen, maxLen + 1);
        var chars = new char[len];
        for (int i = 0; i < len; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion
}
