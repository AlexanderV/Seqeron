using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology driver-mutation-detection area — ONCO-DRIVER-001.
/// The unit under test is the deterministic Vogelstein 20/20-rule driver detector
/// <see cref="OncologyAnalyzer.IdentifyDriverMutations"/> (and its building blocks
/// <see cref="OncologyAnalyzer.ClassifyGene"/>, <see cref="OncologyAnalyzer.ScoreDriverPotential"/>
/// and the caller-supplied-hotspot membership test
/// <see cref="OncologyAnalyzer.MatchCancerHotspots"/>), implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / KeyNotFound / NullReference / Overflow). Every input must
/// resolve to EITHER a well-defined, theory-correct value OR a *documented,
/// intentional* outcome (here, an <see cref="ArgumentNullException"/> for a null
/// mutation enumerable or null hotspot set). For 20/20-rule driver detection the
/// headline hazards are:
///   • a DivideByZeroException computing a fraction count/N when N == 0 (the
///     documented contract is an empty list ⇒ Ambiguous with both fractions 0,
///     INV-? / §6.1 — NOT a throw, NOT a NaN);
///   • a false positive: an all-passenger spectrum (no gene passes the 20%
///     criterion, no hotspot) must yield ZERO drivers, never invent one
///     (INV-01: drivers ⊆ input);
///   • a fraction escaping [0, 1] — every fraction is a count over N so it must
///     stay in [0, 1], and ScoreDriverPotential = max(f_TSG, f_OG) ∈ [0, 1]
///     (INV-04);
///   • double-count corruption on duplicate hotspots / duplicate recurrent
///     missense — repeats must be counted by the documented rule deterministically
///     and the driver subset must remain a faithful (order-preserving) subset of
///     the input, NOT a corrupted multiset;
///   • a KeyNotFoundException / NullReferenceException on an unknown gene (a gene
///     absent from any driver criterion and absent from the hotspot set must be
///     silently skipped, never crash).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-DRIVER-001 — Driver mutation detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 89.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///   • MC = Malformed Content — невалідний контент (duplicate hotspots, genes
///     that meet no criterion, genes/positions absent from the hotspot catalog).
///     Targets (checklist row 89): "empty mutation list, all-passenger,
///     duplicate hotspots, unknown gene".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Driver_Mutation_Detection.md (docs/algorithms/Oncology/Driver_Mutation_Detection.md):
///   • f_OG = (# missense at recurrent positions) / N; a position is recurrent
///       when it carries ≥ 2 missense (RecurrentPositionMinCount = 2)       (§2.2, INV-05)
///   • f_TSG = (# truncating: Nonsense | Frameshift | SpliceSite) / N       (§2.2, §4.2)
///   • Oncogene  ⟺ f_OG > 0.20 (strict);  TumorSuppressor ⟺ f_TSG > 0.20    (§2.2, INV-02/03)
///       — DriverGeneFractionThreshold = 0.20, comparison is strict '>'.
///   • Dual pass (both > 0.20): the dominant fraction decides; an exact tie
///       ⇒ Ambiguous                                                        (§4.1, §5.4)
///   • Empty mutation list ⇒ Ambiguous, both fractions 0                    (§3.3, §6.1)
///   • ScoreDriverPotential = max(f_TSG, f_OG) ∈ [0, 1]                     (§3.2, INV-04)
///   • IdentifyDriverMutations returns the INPUT SUBSET, in INPUT ORDER, of
///       mutations whose gene is Oncogene/TumorSuppressor OR whose
///       (gene, position) is a known hotspot                                (§3.2, §4.1, INV-01)
///   • Hotspot catalogs are CALLER-SUPPLIED — there is NO hardcoded gene set;
///       null knownHotspots is treated as empty                            (§5.2, §3.1)
///   • Gene symbols matched with ORDINAL (case-sensitive) comparison        (§3.3)
///   • Null mutations / null knownHotspots ⇒ ArgumentNullException          (§3.3, §6.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyDriverFuzzTests
{
    private const double Threshold = DriverGeneFractionThreshold; // 0.20

    private static GeneMutation Mis(string gene, int pos) =>
        new(gene, pos, MutationConsequence.Missense);

    private static GeneMutation Trunc(string gene, int pos) =>
        new(gene, pos, MutationConsequence.Nonsense);

    private static GeneMutation Other(string gene, int pos) =>
        new(gene, pos, MutationConsequence.Other);

    private static IReadOnlySet<(string, int)> Hotspots(params (string, int)[] hs) =>
        new HashSet<(string, int)>(hs);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY classification: both
    // criterion fractions must be FINITE and inside [0, 1] (no DivideByZero NaN
    // from N == 0, no overshoot), MutationCount must be non-negative, and the
    // assigned Role must be exactly the deterministic 20/20-rule decision implied
    // by those two fractions. This is what stops a fuzz test from rubber-stamping
    // a NaN / out-of-range / mislabelled result green.
    private static void AssertWellFormedClassification(DriverGeneClassification c)
    {
        double.IsNaN(c.TruncatingFraction).Should().BeFalse("f_TSG must never be NaN");
        double.IsNaN(c.RecurrentMissenseFraction).Should().BeFalse("f_OG must never be NaN");
        double.IsInfinity(c.TruncatingFraction).Should().BeFalse("f_TSG must be finite");
        double.IsInfinity(c.RecurrentMissenseFraction).Should().BeFalse("f_OG must be finite");

        c.TruncatingFraction.Should().BeInRange(0.0, 1.0, "f_TSG = count/N ∈ [0, 1] (§2.2)");
        c.RecurrentMissenseFraction.Should().BeInRange(0.0, 1.0, "f_OG = count/N ∈ [0, 1] (§2.2)");
        c.MutationCount.Should().BeGreaterThanOrEqualTo(0, "the denominator N is a count");

        // Role must be the exact deterministic 20/20-rule decision (§2.2, §4.1, §5.4).
        bool isTsg = c.TruncatingFraction > Threshold;
        bool isOg = c.RecurrentMissenseFraction > Threshold;
        DriverGeneRole expected;
        if (isTsg && isOg)
        {
            expected = c.TruncatingFraction > c.RecurrentMissenseFraction ? DriverGeneRole.TumorSuppressor
                     : c.RecurrentMissenseFraction > c.TruncatingFraction ? DriverGeneRole.Oncogene
                     : DriverGeneRole.Ambiguous;
        }
        else if (isTsg)
        {
            expected = DriverGeneRole.TumorSuppressor;
        }
        else if (isOg)
        {
            expected = DriverGeneRole.Oncogene;
        }
        else
        {
            expected = DriverGeneRole.Ambiguous;
        }

        c.Role.Should().Be(expected, "Role must follow the strict-'>' 20/20 rule from the fractions");
    }

    // Confirms drivers ⊆ input as a faithful, order-preserving subsequence (INV-01).
    private static void AssertOrderedSubsetOfInput(
        IReadOnlyList<GeneMutation> drivers, IReadOnlyList<GeneMutation> input)
    {
        drivers.Count.Should().BeLessThanOrEqualTo(input.Count, "drivers ⊆ input (INV-01)");

        int j = 0;
        foreach (var d in drivers)
        {
            while (j < input.Count && !input[j].Equals(d))
            {
                j++;
            }

            j.Should().BeLessThan(input.Count, "each driver must appear in input order (INV-01)");
            j++;
        }
    }

    #region ONCO-DRIVER-001 — positive sanity (a real driver IS a driver; a passenger is NOT)

    // KRAS codon 12 is the canonical activating hotspot; ten missense all at G12 →
    // recurrent-missense fraction 1.0 > 0.20 → Oncogene → every mutation is a driver.
    [Test]
    public void ClassifyGene_RecurrentHotspotOncogene_IsOncogene_HandComputed()
    {
        var kras = Enumerable.Range(0, 10).Select(_ => Mis("KRAS", 12)).ToList();

        var c = ClassifyGene(kras);

        c.Role.Should().Be(DriverGeneRole.Oncogene);
        c.RecurrentMissenseFraction.Should().Be(1.0, "all 10 missense sit at the recurrent codon 12");
        c.TruncatingFraction.Should().Be(0.0);
        c.MutationCount.Should().Be(10);
        ScoreDriverPotential(kras).Should().Be(1.0);
    }

    // TP53 inactivated by dispersed truncating mutations → f_TSG = 1.0 > 0.20 → TumorSuppressor.
    [Test]
    public void ClassifyGene_DispersedTruncating_IsTumorSuppressor_HandComputed()
    {
        var tp53 = new[]
        {
            Trunc("TP53", 213),
            new GeneMutation("TP53", 245, MutationConsequence.Frameshift),
            new GeneMutation("TP53", 196, MutationConsequence.SpliceSite),
            Trunc("TP53", 306),
        };

        var c = ClassifyGene(tp53);

        c.Role.Should().Be(DriverGeneRole.TumorSuppressor);
        c.TruncatingFraction.Should().Be(1.0, "all 4 mutations are truncating");
        c.RecurrentMissenseFraction.Should().Be(0.0);
        ScoreDriverPotential(tp53).Should().Be(1.0);
    }

    // A clearly passenger gene: scattered single missense (no recurrence) + an Other →
    // f_OG = 0, f_TSG = 0 → Ambiguous → NOT a driver gene.
    [Test]
    public void ClassifyGene_ScatteredPassengers_IsAmbiguous_NotDriver()
    {
        var passenger = new[] { Mis("OR4F5", 10), Mis("OR4F5", 99), Other("OR4F5", 250) };

        var c = ClassifyGene(passenger);

        c.Role.Should().Be(DriverGeneRole.Ambiguous, "no recurrent missense, no truncating");
        c.RecurrentMissenseFraction.Should().Be(0.0, "two missense at DIFFERENT positions are not recurrent");
        c.TruncatingFraction.Should().Be(0.0);
    }

    // Mixed cohort: a true oncogene (KRAS G12 ×3) + scattered passengers (PASS) →
    // only the KRAS mutations come back, in input order.
    [Test]
    public void IdentifyDriverMutations_DriverPlusPassengers_ReturnsOnlyDrivers_InOrder()
    {
        var input = new[]
        {
            Mis("PASS", 5),     // passenger gene, single position
            Mis("KRAS", 12),    // driver gene
            Other("PASS", 80),  // passenger gene
            Mis("KRAS", 12),    // driver gene (recurrent)
            Mis("KRAS", 12),    // driver gene (recurrent)
        };

        var drivers = IdentifyDriverMutations(input);

        drivers.Should().HaveCount(3);
        drivers.Should().OnlyContain(m => m.Gene == "KRAS");
        AssertOrderedSubsetOfInput(drivers, input);
    }

    // Recurrence is hand-checked: 2 missense at codon 600 (recurrent) + 1 elsewhere →
    // recurrent-missense count = 2, N = 3, f_OG = 2/3 > 0.20 → Oncogene.
    [Test]
    public void ClassifyGene_RecurrenceCount_MatchesHandComputedFraction()
    {
        var braf = new[] { Mis("BRAF", 600), Mis("BRAF", 600), Mis("BRAF", 469) };

        var c = ClassifyGene(braf);

        c.RecurrentMissenseFraction.Should().BeApproximately(2.0 / 3.0, 1e-12,
            "the 2 recurrent-position missense / 3 total; the single 469 is not recurrent");
        c.Role.Should().Be(DriverGeneRole.Oncogene);
    }

    #endregion

    #region ONCO-DRIVER-001 — BE: empty mutation list (Ambiguous, fractions 0, no DivideByZero)

    [Test]
    public void ClassifyGene_EmptyList_Ambiguous_ZeroFractions_NoDivideByZero()
    {
        var c = ClassifyGene(Array.Empty<GeneMutation>());

        c.Role.Should().Be(DriverGeneRole.Ambiguous, "no evidence to classify (§6.1)");
        c.TruncatingFraction.Should().Be(0.0, "no DivideByZero / NaN at N = 0");
        c.RecurrentMissenseFraction.Should().Be(0.0);
        c.MutationCount.Should().Be(0);
        AssertWellFormedClassification(c);
    }

    [Test]
    public void ScoreDriverPotential_EmptyList_Zero_NoNaN()
    {
        double score = ScoreDriverPotential(Array.Empty<GeneMutation>());

        double.IsNaN(score).Should().BeFalse("max(0, 0) must not be NaN at N = 0");
        score.Should().Be(0.0);
    }

    [Test]
    public void IdentifyDriverMutations_EmptyList_ReturnsEmpty_NoThrow()
    {
        var drivers = IdentifyDriverMutations(Array.Empty<GeneMutation>());

        drivers.Should().BeEmpty();
    }

    [Test]
    public void IdentifyDriverMutations_EmptyList_NonEmptyHotspots_StillEmpty()
    {
        // A populated catalog cannot conjure drivers from no input (INV-01).
        var drivers = IdentifyDriverMutations(
            Array.Empty<GeneMutation>(), Hotspots(("KRAS", 12), ("BRAF", 600)));

        drivers.Should().BeEmpty();
    }

    [Test]
    public void NullMutations_ThrowArgumentNull()
    {
        FluentActions.Invoking(() => ClassifyGene(null!))
            .Should().Throw<ArgumentNullException>();
        FluentActions.Invoking(() => ScoreDriverPotential(null!))
            .Should().Throw<ArgumentNullException>();
        FluentActions.Invoking(() => IdentifyDriverMutations(null!))
            .Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void MatchCancerHotspots_NullHotspots_ThrowsArgumentNull()
    {
        FluentActions.Invoking(() => MatchCancerHotspots(Mis("KRAS", 12), null!))
            .Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region ONCO-DRIVER-001 — MC: all-passenger (no false positives, drivers ⊆ input)

    // No gene meets either criterion and no hotspots are supplied → ZERO drivers.
    [Test]
    public void IdentifyDriverMutations_AllPassenger_NoDrivers_NoFalsePositive()
    {
        var input = new[]
        {
            Mis("GENE_A", 10),   // single position, not recurrent
            Mis("GENE_A", 200),  // different position
            Other("GENE_B", 33),
            Mis("GENE_C", 7),
        };

        var drivers = IdentifyDriverMutations(input);

        drivers.Should().BeEmpty("no gene passes the 20% rule and no hotspots supplied");
    }

    // Exactly-at-threshold must NOT classify (strict '>' per "more than 20%", §5.4 #1).
    // 1 truncating out of 5 = 0.20 exactly → NOT a TumorSuppressor → no drivers.
    [Test]
    public void IdentifyDriverMutations_FractionExactlyTwentyPercent_NotDriver()
    {
        var input = new[]
        {
            Trunc("BOUND", 1),  // 1 truncating ...
            Other("BOUND", 2),
            Other("BOUND", 3),
            Other("BOUND", 4),
            Other("BOUND", 5),  // ... out of 5 → f_TSG = 0.20 exactly
        };

        ClassifyGene(input).TruncatingFraction.Should().BeApproximately(0.20, 1e-12);
        IdentifyDriverMutations(input).Should().BeEmpty("0.20 is NOT > 0.20 (strict)");
    }

    // Random all-passenger fuzz: every mutation is a single, distinct-position missense
    // in a gene of its own, so no recurrence anywhere → never a driver, ever.
    [Test]
    [CancelAfter(20_000)]
    public void IdentifyDriverMutations_RandomNonRecurrentSingletons_NeverDriver()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 40);
            var input = new List<GeneMutation>(n);
            for (int i = 0; i < n; i++)
            {
                // Unique gene per mutation ⇒ N = 1 per gene ⇒ no recurrence, f_TSG = 0.
                input.Add(Mis($"G{seed}_{i}", rng.Next(1, 1000)));
            }

            var drivers = IdentifyDriverMutations(input);

            drivers.Should().BeEmpty($"seed {seed}: singleton missense across distinct genes are passengers");
            AssertOrderedSubsetOfInput(drivers, input);
        }
    }

    #endregion

    #region ONCO-DRIVER-001 — MC: duplicate hotspots (deterministic dedup/recurrence, no corruption)

    // A duplicated known hotspot in the input must surface EVERY occurrence as a driver
    // (drivers preserve the input multiset/order), not be collapsed or double-counted.
    [Test]
    public void IdentifyDriverMutations_DuplicateHotspotMutations_AllOccurrencesReturned()
    {
        var hotspots = Hotspots(("IDH1", 132));
        var input = new[]
        {
            Mis("IDH1", 132),  // hotspot
            Mis("IDH1", 132),  // duplicate hotspot
            Mis("IDH1", 132),  // duplicate hotspot
        };

        // Note: IDH1 here is ALSO an oncogene (recurrent missense), so the test is
        // robust whether the driver call comes from the gene rule or the hotspot rule.
        var drivers = IdentifyDriverMutations(input, hotspots);

        drivers.Should().HaveCount(3, "each occurrence is a distinct observed mutation (no dedup of input)");
        AssertOrderedSubsetOfInput(drivers, input);
    }

    // Duplicate hotspot at a NON-driver gene: the hotspot rule alone must catch every
    // occurrence even though the gene fails the 20/20 rule (single Other mutations).
    [Test]
    public void IdentifyDriverMutations_DuplicateHotspot_NonDriverGene_HotspotRuleCatchesAll()
    {
        var hotspots = Hotspots(("WEIRD", 50));
        var input = new[]
        {
            Other("WEIRD", 50),  // hotspot position, non-truncating, gene is Ambiguous
            Other("WEIRD", 50),  // duplicate hotspot
            Other("WEIRD", 99),  // same gene, NOT a hotspot position → not a driver
        };

        ClassifyGene(input).Role.Should().Be(DriverGeneRole.Ambiguous, "no 20/20 signal");

        var drivers = IdentifyDriverMutations(input, hotspots);

        drivers.Should().HaveCount(2, "both position-50 occurrences match the hotspot; position 99 does not");
        drivers.Should().OnlyContain(m => m.ProteinPosition == 50);
    }

    // Duplicate missense at the SAME position is the recurrence signal: counting must be
    // deterministic and the recurrent-missense count must equal the hand-computed value.
    [Test]
    public void ClassifyGene_DuplicateRecurrentMissense_DeterministicHandComputedCount()
    {
        // 3 missense at codon 12, 1 missense at codon 13 (singleton), 1 Other.
        var input = new[]
        {
            Mis("KRAS", 12), Mis("KRAS", 12), Mis("KRAS", 12),
            Mis("KRAS", 13),     // singleton missense → NOT recurrent
            Other("KRAS", 61),
        };

        var c = ClassifyGene(input);

        // Recurrent missense = 3 (codon 12 only); N = 5 → f_OG = 0.6.
        c.RecurrentMissenseFraction.Should().BeApproximately(0.6, 1e-12);
        c.Role.Should().Be(DriverGeneRole.Oncogene);
    }

    // Determinism: classification and identification are invariant to call repetition
    // AND to input ordering of duplicates (the rule is a set/count operation, §4.1).
    [Test]
    public void IdentifyDriverMutations_DuplicateHotspots_Deterministic_OrderInsensitiveCount()
    {
        var hotspots = Hotspots(("EGFR", 858));
        var a = new[] { Mis("EGFR", 858), Mis("EGFR", 858), Mis("EGFR", 790) };
        var b = new[] { Mis("EGFR", 790), Mis("EGFR", 858), Mis("EGFR", 858) };

        var da1 = IdentifyDriverMutations(a, hotspots);
        var da2 = IdentifyDriverMutations(a, hotspots);

        da1.Should().Equal(da2, "repeated calls are deterministic");

        // Same multiset, different order: same DRIVER COUNT and same classification fractions.
        var db = IdentifyDriverMutations(b, hotspots);
        db.Count.Should().Be(da1.Count, "reordering duplicates must not change the driver count");

        var ca = ClassifyGene(a);
        var cb = ClassifyGene(b);
        cb.RecurrentMissenseFraction.Should().BeApproximately(ca.RecurrentMissenseFraction, 1e-12);
        cb.TruncatingFraction.Should().BeApproximately(ca.TruncatingFraction, 1e-12);
    }

    [Test]
    [CancelAfter(20_000)]
    public void IdentifyDriverMutations_RandomDuplicatedHotspots_NoCorruption()
    {
        var hotspots = Hotspots(("KRAS", 12), ("BRAF", 600), ("IDH1", 132));

        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 30);
            var pool = new (string, int)[] { ("KRAS", 12), ("BRAF", 600), ("PASS", 7) };
            var input = new List<GeneMutation>(n);
            for (int i = 0; i < n; i++)
            {
                var (g, p) = pool[rng.Next(pool.Length)];
                input.Add(Mis(g, p)); // heavy duplication of hotspot positions
            }

            var drivers = IdentifyDriverMutations(input, hotspots);

            AssertOrderedSubsetOfInput(drivers, input);

            // Every KRAS/BRAF mutation IS a driver (hotspot rule) — count them by hand.
            int expectedHotspotHits = input.Count(m =>
                hotspots.Contains((m.Gene, m.ProteinPosition)));
            drivers.Count.Should().BeGreaterThanOrEqualTo(expectedHotspotHits,
                $"seed {seed}: every hotspot occurrence must be a driver, none dropped");

            // PASS appears only as singletons-per-position across the pool? No: it may
            // recur. Verify PASS drivers (if any) are exactly those whose gene was
            // classified as a driver — never a spurious extra.
            foreach (var d in drivers)
            {
                bool byHotspot = hotspots.Contains((d.Gene, d.ProteinPosition));
                bool byGene = ClassifyGene(input.Where(m => m.Gene == d.Gene)).Role
                    != DriverGeneRole.Ambiguous;
                (byHotspot || byGene).Should().BeTrue(
                    $"seed {seed}: a driver must be justified by gene rule or hotspot rule");
            }
        }
    }

    #endregion

    #region ONCO-DRIVER-001 — MC: unknown gene (no KeyNotFound/NullReference; skipped/unclassified)

    // A gene that is in NO hotspot catalog and meets NO 20/20 criterion is "unknown" to
    // the detector: it must be silently skipped, never throw KeyNotFound/NullReference.
    [Test]
    public void IdentifyDriverMutations_UnknownGene_NotInHotspotSet_SkippedNoThrow()
    {
        var hotspots = Hotspots(("KRAS", 12)); // catalog does NOT mention the unknown gene
        var input = new[]
        {
            Mis("ZZZ_UNKNOWN", 999),  // unknown gene, unknown position
            Other("ZZZ_UNKNOWN", 1),
        };

        IReadOnlyList<GeneMutation> drivers = null!;
        FluentActions.Invoking(() => drivers = IdentifyDriverMutations(input, hotspots))
            .Should().NotThrow("unknown genes are skipped, not looked up unsafely");

        drivers.Should().BeEmpty("an unknown gene with no driver signal yields no drivers");
    }

    // Unknown gene mixed with a real driver: the driver still surfaces, the unknown is
    // ignored, and nothing throws on the unknown's absence from the catalog.
    [Test]
    public void IdentifyDriverMutations_UnknownGeneAlongsideDriver_OnlyDriverReturned()
    {
        var hotspots = Hotspots(("KRAS", 12));
        var input = new[]
        {
            Mis("UNK1", 5),
            Mis("KRAS", 12),
            Mis("UNK2", 77),
        };

        var drivers = IdentifyDriverMutations(input, hotspots);

        drivers.Should().ContainSingle().Which.Gene.Should().Be("KRAS");
    }

    // MatchCancerHotspots on an unknown (gene, position) must simply return false.
    [Test]
    public void MatchCancerHotspots_UnknownGene_ReturnsFalse_NoThrow()
    {
        var hotspots = Hotspots(("KRAS", 12), ("BRAF", 600));

        MatchCancerHotspots(Mis("KRAS", 12), hotspots).Should().BeTrue();
        MatchCancerHotspots(Mis("KRAS", 13), hotspots).Should().BeFalse("right gene, wrong position");
        MatchCancerHotspots(Mis("ZZZ", 12), hotspots).Should().BeFalse("unknown gene");
    }

    // Ordinal (case-sensitive) gene matching: "kras" is a DIFFERENT, unknown gene from
    // "KRAS" (§3.3) — it must not match the catalog and must not classify off "KRAS".
    [Test]
    public void IdentifyDriverMutations_CaseMismatchGene_TreatedAsUnknown()
    {
        var hotspots = Hotspots(("KRAS", 12));
        var input = new[] { Mis("kras", 12), Mis("kras", 12) }; // lowercase ⇒ different gene

        var drivers = IdentifyDriverMutations(input, hotspots);

        // The hotspot key ("KRAS",12) does NOT match ("kras",12) under Ordinal.
        drivers.Should().HaveCount(2,
            "lowercase 'kras' is its own gene; its 2 recurrent missense make IT an oncogene");
        drivers.Should().OnlyContain(m => m.Gene == "kras");

        // But it must NOT be matched via the uppercase hotspot catalog.
        MatchCancerHotspots(Mis("kras", 12), hotspots).Should().BeFalse("ordinal case-sensitive (§3.3)");
    }

    [Test]
    [CancelAfter(20_000)]
    public void IdentifyDriverMutations_RandomUnknownGenes_NeverThrow_WellFormed()
    {
        var hotspots = Hotspots(("KRAS", 12), ("TP53", 273));

        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 50);
            var consequences = Enum.GetValues<MutationConsequence>();
            var input = new List<GeneMutation>(n);
            for (int i = 0; i < n; i++)
            {
                // Random, mostly-unknown gene symbols + extreme/boundary positions.
                string gene = $"UNK_{rng.Next(0, 100)}";
                int pos = rng.Next(0, 5) switch
                {
                    0 => 0,
                    1 => -1,
                    2 => int.MaxValue,
                    3 => int.MinValue,
                    _ => rng.Next(1, 10_000),
                };
                input.Add(new GeneMutation(gene, pos, consequences[rng.Next(consequences.Length)]));
            }

            IReadOnlyList<GeneMutation> drivers = null!;
            FluentActions.Invoking(() => drivers = IdentifyDriverMutations(input, hotspots))
                .Should().NotThrow($"seed {seed}: unknown genes/extreme positions must not crash");

            AssertOrderedSubsetOfInput(drivers, input);

            // Per-gene classifications must all be well-formed (no NaN/out-of-range).
            foreach (var byGene in input.GroupBy(m => m.Gene))
            {
                AssertWellFormedClassification(ClassifyGene(byGene));
            }
        }
    }

    #endregion
}
