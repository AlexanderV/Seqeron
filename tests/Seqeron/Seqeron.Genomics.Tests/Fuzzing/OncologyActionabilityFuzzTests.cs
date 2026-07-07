using Level = Seqeron.Genomics.Oncology.OncologyAnalyzer.OncoKbLevel;
using Assoc = Seqeron.Genomics.Oncology.OncologyAnalyzer.TherapyAssociation;
using Input = Seqeron.Genomics.Oncology.OncologyAnalyzer.VariantActionabilityInput;
using Assessment = Seqeron.Genomics.Oncology.OncologyAnalyzer.ActionabilityAssessment;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology clinical-actionability area — ONCO-ACTION-001.
/// The unit under test is the OncoKB therapeutic-levels-of-evidence tiering engine
/// <see cref="OncologyAnalyzer.ClassifyActionabilityLevel"/> (single-variant highest
/// combined level), its batch wrapper <see cref="OncologyAnalyzer.AssessActionability"/>
/// (per-variant highest sensitive / resistance / combined level), the presentation
/// ordering <see cref="OncologyAnalyzer.GetTherapyRecommendations"/>, the comparator
/// <see cref="OncologyAnalyzer.CompareLevels"/>, and the grouping predicate
/// <see cref="OncologyAnalyzer.IsStandardCare"/>, implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs. This is the
/// surface that reduces a somatic variant's caller-supplied biomarker–drug associations
/// to their highest OncoKB therapeutic level of evidence — the actionability tiering an
/// oncology pipeline drives to decide whether a variant is therapeutically actionable.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary inputs to a unit and asserts that the code
/// NEVER fails in an undisciplined way: no hang, no nonsense output, and no
/// *unhandled* runtime fault — specifically NO NullReferenceException on an empty
/// association list, NO non-deterministic / order-dependent conflict resolution
/// when two associations assign DIFFERENT levels to one variant, and NO
/// KeyNotFoundException / NullReferenceException on an unknown drug name (a therapy
/// not in the caller's knowledgebase). Every input must resolve to EITHER a
/// well-defined, theory-correct level OR a *documented, intentional* validation
/// outcome (an <see cref="ArgumentNullException"/> for a null variants enumerable or
/// a null associations list — §3.3, §6.1). A silently-wrong tier on degenerate
/// input is just as much a bug as a crash. — docs/ADVANCED_TESTING_CHECKLIST.md §8
/// "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-ACTION-001 — Clinical actionability assessment (OncoKB levels)
/// Checklist: docs/checklists/03_FUZZING.md, row 118.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty).
///     Targets (checklist row 118): "no evidence, conflicting tiers, unknown drug".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// docs/algorithms/Oncology/Clinical_Actionability_Assessment.md (ONCO-ACTION-001):
///   • OncoKB defines seven therapeutic levels (§2.2): sensitivity 1, 2, 3A, 3B, 4
///     and resistance R1, R2 (plus None = no leveled association).
///   • Combined order (INV-01, §2.2): R1 > 1 > 2 > 3A > 3B > 4 > R2. The highest
///     combined level is the maximum association under this order; conflicting
///     levels resolve DETERMINISTICALLY to the maximum (highest-tier-wins), NOT to
///     an order-dependent value. {1, R1} -> R1 (§6.1). This precedence is encoded
///     in the integer order of the OncoKbLevel enum (None lowest … R1 highest), so
///     CompareLevels is an integer comparison (§4.2).
///   • Highest sensitive level (INV-02): 1 > 2 > 3A > 3B > 4, ignoring R1/R2.
///   • Highest resistance level (INV-03): R1 > R2, ignoring sensitivity levels.
///   • No-evidence (INV-04, §5.4, §6.1): a variant with ZERO leveled associations
///     yields None on every axis; IsActionable = false. An empty association list is
///     VALID (never null), no NullReference.
///   • Unknown drug (§3.3, §6.2): the classifier is "case- and content-agnostic
///     about drug names; it ranks only the supplied OncoKbLevel values." A drug not
///     in any knowledgebase is just an arbitrary string on a TherapyAssociation — it
///     is ranked by its level, never looked up, so NO KeyNotFound / NullReference.
///   • Standard-care grouping (§4.2): {1, 2, R1}; investigational/hypothetical
///     {3A, 3B, 4, R2}.
///   • Validation (§3.3, §6.1): AssessActionability(null) and
///     ClassifyActionabilityLevel / GetTherapyRecommendations on a variant whose
///     Associations is null throw ArgumentNullException; VariantActionabilityInput
///     rejects a null associations list at construction.
///   • Determinism / order (INV-05): one assessment per input variant, input order
///     preserved; the result depends only on the SET of levels, not their order.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class OncologyActionabilityFuzzTests
{
    // All seven leveled values plus None, in DESCENDING combined order, copied
    // verbatim from the §2.2 / INV-01 precedence — NOT read back from code output.
    private static readonly Level[] CombinedDescending =
    {
        Level.R1, Level.Level1, Level.Level2, Level.Level3A, Level.Level3B, Level.Level4, Level.R2
    };

    private static readonly Level[] AllLevelsIncludingNone =
    {
        Level.None, Level.Level1, Level.Level2, Level.Level3A, Level.Level3B,
        Level.Level4, Level.R1, Level.R2
    };

    private static readonly HashSet<Level> SensitivitySet = new()
    {
        Level.Level1, Level.Level2, Level.Level3A, Level.Level3B, Level.Level4
    };

    private static readonly HashSet<Level> ResistanceSet = new() { Level.R1, Level.R2 };

    private static readonly HashSet<Level> StandardCareSet = new() { Level.Level1, Level.Level2, Level.R1 };

    // Build a variant carrying exactly the given leveled drug associations (one
    // synthetic drug name per association). Drug names are arbitrary — the contract
    // ranks the LEVEL, never the drug.
    private static Input Variant(params Level[] levels) =>
        new("GENE", "p.X1Y", levels.Select((lvl, i) => new Assoc($"Drug{i}", lvl)).ToArray());

    private static Input VariantWithDrugs(params (string drug, Level level)[] assocs) =>
        new("GENE", "p.X1Y", assocs.Select(a => new Assoc(a.drug, a.level)).ToArray());

    // The hand-computed maximum of a level multiset under the combined order, used to
    // independently verify highest-tier-wins resolution. Mirrors the §2.2 ordering;
    // does NOT call the SUT.
    private static Level ExpectedCombined(IEnumerable<Level> levels)
    {
        foreach (var lvl in CombinedDescending)
        {
            if (levels.Contains(lvl))
            {
                return lvl;
            }
        }

        return Level.None;
    }

    private static Level ExpectedRestricted(IEnumerable<Level> levels, HashSet<Level> allowed)
        => ExpectedCombined(levels.Where(allowed.Contains));

    // A well-formed assessment: every reported level is a DEFINED enum value, the
    // sensitive axis is sensitivity-or-None, the resistance axis is
    // resistance-or-None, the combined axis is the documented max over both axes, and
    // IsActionable is exactly (combined != None) (INV-01..INV-04).
    private static void AssertWellFormed(Assessment a)
    {
        Enum.IsDefined(a.HighestSensitiveLevel).Should().BeTrue();
        Enum.IsDefined(a.HighestResistanceLevel).Should().BeTrue();
        Enum.IsDefined(a.HighestCombinedLevel).Should().BeTrue();

        (a.HighestSensitiveLevel == Level.None || SensitivitySet.Contains(a.HighestSensitiveLevel))
            .Should().BeTrue("the sensitive axis only reports sensitivity levels or None (INV-02)");
        (a.HighestResistanceLevel == Level.None || ResistanceSet.Contains(a.HighestResistanceLevel))
            .Should().BeTrue("the resistance axis only reports resistance levels or None (INV-03)");

        // Combined is the max of the two single-axis maxima under the combined order.
        var expectedCombined = OncologyAnalyzer.CompareLevels(a.HighestSensitiveLevel, a.HighestResistanceLevel) >= 0
            ? a.HighestSensitiveLevel
            : a.HighestResistanceLevel;
        a.HighestCombinedLevel.Should().Be(expectedCombined,
            "combined = max over both axes under R1 > 1 > 2 > 3A > 3B > 4 > R2 (INV-01)");

        a.IsActionable.Should().Be(a.HighestCombinedLevel != Level.None,
            "IsActionable iff there is a leveled association (INV-04)");
    }

    // ════════════════════════════════════════════════════════════════════════
    #region ONCO-ACTION-001 — clinical actionability tiering (OncoKB levels)
    // ════════════════════════════════════════════════════════════════════════

    // ───────────────────────────── POSITIVE SANITY ──────────────────────────
    // The business anchor: on good input the tiering is CORRECT. A clear Tier-I
    // (Level 1) match is assigned the top sensitivity tier and surfaced as the
    // top recommendation; conflicting entries resolve to the documented winner.

    // BRAF p.V600E worked example (§7.1): {Dabrafenib=Level1, OtherDrug=Level3A}
    // -> combined Level1 (1 > 3A), and Dabrafenib is the top recommendation.
    [Test]
    public void Positive_ClearLevel1Match_AssignsTopTierAndTopRecommendation()
    {
        var variant = VariantWithDrugs(("Dabrafenib", Level.Level1), ("OtherDrug", Level.Level3A));

        OncologyAnalyzer.ClassifyActionabilityLevel(variant).Should().Be(Level.Level1,
            "Level 1 outranks Level 3A under the combined order (§7.1 worked example)");

        var assessment = OncologyAnalyzer.AssessActionability(new[] { variant }).Single();
        AssertWellFormed(assessment);
        assessment.HighestCombinedLevel.Should().Be(Level.Level1);
        assessment.HighestSensitiveLevel.Should().Be(Level.Level1);
        assessment.HighestResistanceLevel.Should().Be(Level.None, "no resistance association is present");
        assessment.IsActionable.Should().BeTrue();
        OncologyAnalyzer.IsStandardCare(Level.Level1).Should().BeTrue("Level 1 is standard care (§4.2)");

        var recs = OncologyAnalyzer.GetTherapyRecommendations(variant);
        recs.Should().HaveCount(2);
        recs[0].Drug.Should().Be("Dabrafenib", "most-actionable association is surfaced first");
        recs[0].Level.Should().Be(Level.Level1);
    }

    // ──────────────────────── BE: no evidence (empty list) ───────────────────
    // A variant with ZERO leveled associations: documented None on every axis,
    // NotActionable, NO NullReference on the empty list (INV-04, §5.4, §6.1).

    [Test]
    public void Be_NoEvidence_EmptyAssociations_YieldsNoneOnEveryAxis()
    {
        var variant = Variant(); // empty association list — valid, never null

        OncologyAnalyzer.ClassifyActionabilityLevel(variant).Should().Be(Level.None,
            "no leveled association ⇒ None (INV-04)");

        var assessment = OncologyAnalyzer.AssessActionability(new[] { variant }).Single();
        AssertWellFormed(assessment);
        assessment.HighestSensitiveLevel.Should().Be(Level.None);
        assessment.HighestResistanceLevel.Should().Be(Level.None);
        assessment.HighestCombinedLevel.Should().Be(Level.None);
        assessment.IsActionable.Should().BeFalse("a variant with no evidence is not actionable");

        OncologyAnalyzer.GetTherapyRecommendations(variant).Should().BeEmpty(
            "no associations ⇒ empty recommendation list, never null");
    }

    // Batch of MANY empty-evidence variants: each maps to None, input order
    // preserved, one assessment per input (INV-05). No crash on a degenerate batch.
    [Test]
    [CancelAfter(30000)]
    public void Be_NoEvidence_BatchOfEmptyVariants_AllNoneInputOrder()
    {
        var variants = Enumerable.Range(0, 200)
            .Select(i => new Input($"GENE{i}", $"p.X{i}Y", Array.Empty<Assoc>()))
            .ToArray();

        var results = OncologyAnalyzer.AssessActionability(variants);

        results.Should().HaveCount(variants.Length, "one assessment per input variant (INV-05)");
        for (int i = 0; i < results.Count; i++)
        {
            results[i].Variant.Gene.Should().Be($"GENE{i}", "input order is preserved (INV-05)");
            AssertWellFormed(results[i]);
            results[i].IsActionable.Should().BeFalse();
            results[i].HighestCombinedLevel.Should().Be(Level.None);
        }
    }

    // ──────────────── BE: conflicting tiers (different levels) ───────────────
    // Two (or more) associations assign DIFFERENT levels to one variant. The
    // documented resolution is highest-tier-wins under the combined order, and it
    // must be DETERMINISTIC — independent of the order the associations appear in.

    // Canonical conflict from §6.1: {Level1, R1} -> combined R1 (R1 > 1), while the
    // single-axis maxima split cleanly: sensitive Level1, resistance R1.
    [Test]
    public void Be_ConflictingTiers_Level1AndR1_ResolvesToR1Combined()
    {
        var variant = Variant(Level.Level1, Level.R1);

        var assessment = OncologyAnalyzer.AssessActionability(new[] { variant }).Single();
        AssertWellFormed(assessment);
        assessment.HighestCombinedLevel.Should().Be(Level.R1, "{1, R1} ⇒ R1 (R1 > 1) — §6.1");
        assessment.HighestSensitiveLevel.Should().Be(Level.Level1, "sensitive axis ignores R1 (INV-02)");
        assessment.HighestResistanceLevel.Should().Be(Level.R1, "resistance axis ignores Level 1 (INV-03)");
    }

    // Conflict resolution is ORDER-INDEPENDENT: every permutation of a fixed level
    // multiset yields the SAME combined/sensitive/resistance maxima. This is the
    // core BE guard against "order-dependent garbage" tie resolution.
    [Test]
    [CancelAfter(60000)]
    public void Be_ConflictingTiers_ResolutionIsOrderIndependent()
    {
        var rng = new Random(118_2024);

        for (int trial = 0; trial < 400; trial++)
        {
            // A random multiset of 2..6 conflicting leveled associations.
            int k = 2 + rng.Next(5);
            var levels = Enumerable.Range(0, k)
                .Select(_ => CombinedDescending[rng.Next(CombinedDescending.Length)])
                .ToArray();

            var expectedCombined = ExpectedCombined(levels);
            var expectedSensitive = ExpectedRestricted(levels, SensitivitySet);
            var expectedResistance = ExpectedRestricted(levels, ResistanceSet);

            // Shuffle the SAME multiset; the assessment must be identical.
            var shuffled = levels.OrderBy(_ => rng.Next()).ToArray();

            var a = OncologyAnalyzer.AssessActionability(new[] { Variant(levels) }).Single();
            var b = OncologyAnalyzer.AssessActionability(new[] { Variant(shuffled) }).Single();

            AssertWellFormed(a);
            AssertWellFormed(b);

            a.HighestCombinedLevel.Should().Be(expectedCombined,
                "highest-tier-wins is the documented resolution (INV-01)");
            a.HighestSensitiveLevel.Should().Be(expectedSensitive);
            a.HighestResistanceLevel.Should().Be(expectedResistance);

            b.HighestCombinedLevel.Should().Be(a.HighestCombinedLevel,
                "conflict resolution must not depend on association order");
            b.HighestSensitiveLevel.Should().Be(a.HighestSensitiveLevel);
            b.HighestResistanceLevel.Should().Be(a.HighestResistanceLevel);
        }
    }

    // Adding a STRICTLY-LOWER conflicting association never changes the combined
    // winner (monotonicity of the max); adding a strictly-higher one always wins.
    [Test]
    public void Be_ConflictingTiers_LowerEntryNeverWins_HigherAlwaysWins()
    {
        // Pairs (higher, lower) adjacent in the combined order.
        for (int i = 0; i + 1 < CombinedDescending.Length; i++)
        {
            var higher = CombinedDescending[i];
            var lower = CombinedDescending[i + 1];

            OncologyAnalyzer.ClassifyActionabilityLevel(Variant(higher, lower)).Should().Be(higher,
                $"{higher} outranks {lower}; adding the lower tier cannot win");
            OncologyAnalyzer.ClassifyActionabilityLevel(Variant(lower, higher)).Should().Be(higher,
                "resolution is order-independent");
        }
    }

    // ──────────────────────── BE: unknown drug ───────────────────────────────
    // A drug name not in any knowledgebase is just an arbitrary string. The
    // classifier is content-agnostic about drug names (§3.3): it ranks the LEVEL,
    // never looks up the drug, so NO KeyNotFound / NullReference regardless of the
    // string — including empty, whitespace, null, unicode, and very long names.

    [Test]
    public void Be_UnknownDrug_DegenerateNames_RankedByLevelNoLookupFault()
    {
        var rng = new Random(118_5599);
        var weirdNames = new[]
        {
            "",                       // empty drug name
            "   ",                    // whitespace only
            "Юнікод-Препарат",        // non-ASCII unicode
            "drug\twith\nctl",        // control chars
            new string('X', 4096),    // very long
            "NotInAnyKnowledgebase",  // plausible-but-unknown therapy
            "💊",                     // emoji
        };

        foreach (var name in weirdNames)
        {
            var level = CombinedDescending[rng.Next(CombinedDescending.Length)];
            var variant = VariantWithDrugs((name, level));

            // No KeyNotFound / NullReference; the level alone drives the result.
            OncologyAnalyzer.ClassifyActionabilityLevel(variant).Should().Be(level,
                "drug name is irrelevant to ranking; only the level matters (§3.3)");

            var assessment = OncologyAnalyzer.AssessActionability(new[] { variant }).Single();
            AssertWellFormed(assessment);
            assessment.HighestCombinedLevel.Should().Be(level);

            var recs = OncologyAnalyzer.GetTherapyRecommendations(variant);
            recs.Should().ContainSingle();
            recs[0].Drug.Should().Be(name, "the unknown drug name is preserved verbatim, not dropped");
            recs[0].Level.Should().Be(level);
        }
    }

    // A null drug name is a valid TherapyAssociation field (record struct, no
    // validation); it must still rank by level with no NullReference.
    [Test]
    public void Be_UnknownDrug_NullName_StillRanksByLevel()
    {
        var variant = VariantWithDrugs((null!, Level.Level2), ("Real", Level.Level3A));

        OncologyAnalyzer.ClassifyActionabilityLevel(variant).Should().Be(Level.Level2,
            "a null drug name does not crash ranking; Level 2 > Level 3A");

        var recs = OncologyAnalyzer.GetTherapyRecommendations(variant);
        recs[0].Level.Should().Be(Level.Level2, "highest-leveled association is surfaced first");
    }

    // Mixed known + many unknown drugs: the highest level still wins regardless of
    // how many unknown-drug associations surround it.
    [Test]
    [CancelAfter(30000)]
    public void Be_UnknownDrug_ManyUnknownsDoNotMaskTheTopTier()
    {
        var rng = new Random(118_7777);

        for (int trial = 0; trial < 300; trial++)
        {
            // One guaranteed top-tier (R1) known drug among many unknown low-tier ones.
            var assocs = new List<(string, Level)> { ("KnownTopTier", Level.R1) };
            int noise = rng.Next(1, 12);
            for (int i = 0; i < noise; i++)
            {
                assocs.Add(($"Unknown_{Guid.NewGuid()}", Level.Level4));
            }

            var shuffled = assocs.OrderBy(_ => rng.Next()).ToArray();
            var variant = VariantWithDrugs(shuffled);

            OncologyAnalyzer.ClassifyActionabilityLevel(variant).Should().Be(Level.R1,
                "R1 is the combined maximum; unknown-drug noise cannot mask it");

            var recs = OncologyAnalyzer.GetTherapyRecommendations(variant);
            recs[0].Level.Should().Be(Level.R1, "most-actionable association ordered first");
            recs.Should().BeInDescendingOrder(r => (int)r.Level,
                "recommendations are ordered by descending combined level");
        }
    }

    // ──────────────────── BE: validation boundary (null) ─────────────────────
    // The documented intentional faults (§3.3, §6.1): these are NOT bugs — they are
    // the contract. A null variants enumerable / null associations throws
    // ArgumentNullException, distinct from the silent crashes BE hunts for.

    [Test]
    public void Be_NullVariantsBatch_ThrowsArgumentNullException()
    {
        Action act = () => OncologyAnalyzer.AssessActionability(null!);
        act.Should().Throw<ArgumentNullException>("a null batch is rejected (§3.3)");
    }

    [Test]
    public void Be_NullAssociations_RejectedAtConstruction()
    {
        Action act = () => new Input("GENE", "p.X1Y", null!);
        act.Should().Throw<ArgumentNullException>(
            "VariantActionabilityInput rejects a null associations list (§3.3)");
    }

    // ──────────────── BE: comparator / grouping totality ─────────────────────
    // CompareLevels must be a total order consistent with the documented combined
    // precedence; None ranks below every leveled value; IsStandardCare matches the
    // documented {1, 2, R1} set exactly. Exhaustive over all enum values.

    [Test]
    public void Be_CompareLevels_TotalOrderConsistentWithCombinedPrecedence()
    {
        // Strictly descending under the documented order.
        for (int i = 0; i + 1 < CombinedDescending.Length; i++)
        {
            OncologyAnalyzer.CompareLevels(CombinedDescending[i], CombinedDescending[i + 1])
                .Should().BePositive($"{CombinedDescending[i]} > {CombinedDescending[i + 1]}");
        }

        foreach (var lvl in CombinedDescending)
        {
            OncologyAnalyzer.CompareLevels(lvl, Level.None).Should().BePositive(
                "None ranks below every leveled value");
            OncologyAnalyzer.CompareLevels(lvl, lvl).Should().Be(0, "reflexive equality");
        }

        // Antisymmetry across all pairs.
        foreach (var x in AllLevelsIncludingNone)
        {
            foreach (var y in AllLevelsIncludingNone)
            {
                Math.Sign(OncologyAnalyzer.CompareLevels(x, y))
                    .Should().Be(-Math.Sign(OncologyAnalyzer.CompareLevels(y, x)),
                        $"compare({x},{y}) must be the negation of compare({y},{x})");
            }
        }
    }

    [Test]
    public void Be_IsStandardCare_ExactlyDocumentedSet()
    {
        foreach (var lvl in AllLevelsIncludingNone)
        {
            OncologyAnalyzer.IsStandardCare(lvl).Should().Be(StandardCareSet.Contains(lvl),
                $"standard care is exactly {{1, 2, R1}}; {lvl} membership must match (§4.2)");
        }
    }

    #endregion
}
