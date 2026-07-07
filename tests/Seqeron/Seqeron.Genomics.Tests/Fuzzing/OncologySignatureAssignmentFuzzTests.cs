using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for Oncology mutational-PROCESS classification — ONCO-SIG-004.
/// The unit under test is <see cref="OncologyAnalyzer.ClassifyMutationalProcess"/>: it takes per-signature
/// exposures (COSMIC SBS label + non-negative activity), normalises them to relative contributions
/// wᵢ = eᵢ / Σe, drops sub-cutoff signatures (deconstructSigs τ = 0.06, strict &lt;), maps the survivors to
/// their COSMIC mutational process via <see cref="OncologyAnalyzer.GetMutationalProcess"/>, sums the
/// surviving contributions per process, ranks the active processes by DESCENDING contribution (ties broken
/// by the process enum, ascending), and reports the largest-contribution process as the DominantProcess.
/// Implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This is SIG-004 (the exposure→process ASSIGNMENT / RANKING step), distinct from its three siblings in
/// the mutational-signature family:
///   • SIG-001 — the SBS-96 context CATALOGUE builder (OncologySignatureContextFuzzTests),
///   • SIG-002 — the NNLS exposure FITTING (OncologySignatureFittingFuzzTests),
///   • SIG-003 — the exposure BOOTSTRAP confidence intervals (OncologySignatureBootstrapFuzzTests).
/// SIG-004 consumes the exposures the others produce and assigns/ranks them into processes; it does no
/// fitting, no bootstrap, and no catalogue construction.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed exposure vectors and asserts the code NEVER fails in an
/// undisciplined way: no DivideByZero / NaN when normalising an ALL-ZERO exposure vector by a zero total
/// (the headline BE hazard here), no non-deterministic / order-dependent tie-break when two processes have
/// EXACTLY equal contribution, no proportion outside [0,1], no active set whose order disagrees with the
/// documented descending-contribution-then-enum rule, and the documented Argument* exceptions for invalid
/// input. Every input resolves to EITHER a well-defined, theory-correct classification OR a documented,
/// intentional outcome (empty active set + Unknown dominant for a zero total; ArgumentNullException /
/// ArgumentException / ArgumentOutOfRangeException for null / negative / NaN / out-of-range input).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-SIG-004 — mutational-process classification (exposure ranking/assignment), Oncology
/// Checklist: docs/checklists/03_FUZZING.md, row 99.
/// Fuzz strategy exercised: BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty).
///   Targets (checklist row 99): "tied exposures, all-zero exposures".
///     • TIED EXPOSURES — two or more processes with EXACTLY equal aggregated contribution. The documented
///       tie-break is DETERMINISTIC: descending contribution, then by process enum (ascending) — so the
///       lowest-valued enum among the tied processes is the dominant one (INV-04, §4.1 step 5, §3.2). The
///       result must be byte-for-byte identical regardless of the INPUT ORDER of the exposure list (a stable,
///       order-independent ranking), never "whatever the dictionary iteration happened to yield".
///     • ALL-ZERO EXPOSURES — every exposure is 0 (or the list is empty), so the total T = Σe = 0 and the
///       normalisation wᵢ = eᵢ / T is undefined. The documented contract: empty active set, dominant =
///       Unknown, with NO DivideByZero and NO NaN leaking out of the 0/0 normalisation (INV-05, §6.1).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Mutational_Process_Classification.md (docs/algorithms/Oncology/Mutational_Process_Classification.md):
///   • wᵢ = eᵢ / Σe; present iff wᵢ ≥ τ (τ = 0.06, strict &lt; excludes; 0.06 retained).            (§2.2, INV-02)
///   • C(P) = Σ surviving member contributions (additive weights).                                 (§2.2, INV-03)
///   • Active set = {P : C(P) > 0}; dominant = argmax_P C(P).                                       (§2.2, INV-04)
///   • Active processes ordered by DESCENDING contribution, then by process enum (ties).           (§3.2, §4.1)
///   • T = 0 (all-zero / empty) ⇒ empty active set, Unknown dominant, no ÷0.                        (§6.1, INV-05)
///   • Each surviving wᵢ ∈ [0,1]; Σ surviving ≤ 1 (sub-cutoff mass dropped).                        (INV-01)
///   • Signature-label matching is case-insensitive; unmapped labels ⇒ Unknown, no process.        (§3.3)
///   • null exposures / null label ⇒ ArgumentNullException; negative / NaN exposure ⇒ ArgumentException;
///     cutoff NaN or outside [0,1) ⇒ ArgumentOutOfRangeException.                                   (§3.3)
///   • Worked example (§7.1): SBS2 50, SBS13 30, SBS1 15, SBS4 5 ⇒ APOBEC 0.80, Aging 0.15
///     (SBS4 w=0.05 &lt; 0.06 dropped); DominantProcess == Apobec.
///   • COSMIC map: Aging{SBS1,SBS5}, APOBEC{SBS2,SBS13}, Tobacco{SBS4}, UV{SBS7a-d},
///     MMRd{SBS6,SBS15,SBS20,SBS26}; enum order Unknown,Aging,Apobec,Tobacco,UV,MMRd.              (§4.2)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologySignatureAssignmentFuzzTests
{
    private const double Eps = 1e-9;

    // The five mapped processes (all non-Unknown enum members) and one representative SBS label each, so a
    // fuzzer can build exposure vectors targeting specific processes with known contributions.
    private static readonly (MutationalProcess Process, string Label)[] MappedSignatures =
    {
        (MutationalProcess.Aging, "SBS1"),
        (MutationalProcess.Apobec, "SBS2"),
        (MutationalProcess.TobaccoSmoking, "SBS4"),
        (MutationalProcess.UltravioletLight, "SBS7a"),
        (MutationalProcess.MismatchRepairDeficiency, "SBS6"),
    };

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural contract on EVERY classification: every active contribution is finite
    // (no NaN/Inf from a 0/0 normalisation) and in (0,1]; the active set is sorted by the documented
    // descending-contribution-then-enum rule; the dominant process is exactly the first active process (or
    // Unknown when the set is empty); a non-empty active set never contains Unknown. This is what stops a
    // fuzz test from rubber-stamping a NaN contribution or a mis-ordered ranking green.
    private static void AssertWellFormedClassification(MutationalProcessClassification c)
    {
        c.ActiveProcesses.Should().NotBeNull();

        foreach (ProcessActivity a in c.ActiveProcesses)
        {
            double.IsFinite(a.Contribution).Should().BeTrue(
                "an aggregated contribution must never be NaN/Infinity (no ÷0 on a zero total, INV-05)");
            a.Contribution.Should().BeGreaterThan(0.0, "the active set holds only processes with C(P) > 0 (INV-04)");
            a.Contribution.Should().BeLessThanOrEqualTo(1.0 + Eps, "each contribution is a proportion ≤ 1 (INV-01)");
            a.Process.Should().NotBe(MutationalProcess.Unknown,
                "an Unknown-aetiology signature contributes to no recognised process (§3.3)");
        }

        // Σ surviving contributions ≤ 1 (sub-cutoff mass is dropped) — INV-01.
        c.ActiveProcesses.Sum(a => a.Contribution).Should().BeLessThanOrEqualTo(1.0 + Eps,
            "surviving contributions sum to ≤ 1 (sub-cutoff mass dropped, INV-01)");

        // Documented ordering: descending contribution, then ascending process enum on ties (§3.2, §4.1).
        // The source's tie-break (OrderByDescending(Contribution).ThenBy(Process)) keys on EXACT double
        // equality, so the helper mirrors that exactly: a strictly-greater predecessor is fine; only on an
        // exact-equal contribution must the enum order be ascending.
        for (int i = 1; i < c.ActiveProcesses.Count; i++)
        {
            ProcessActivity prev = c.ActiveProcesses[i - 1];
            ProcessActivity cur = c.ActiveProcesses[i];
            if (prev.Contribution == cur.Contribution)
            {
                ((int)prev.Process).Should().BeLessThan((int)cur.Process,
                    "exactly-tied contributions are broken by ascending process enum (deterministic, §4.1)");
            }
            else
            {
                prev.Contribution.Should().BeGreaterThan(cur.Contribution,
                    "active processes are ordered by descending contribution (§3.2)");
            }
        }

        // Dominant = first active process (largest contribution), or Unknown when none active (INV-04).
        if (c.ActiveProcesses.Count == 0)
        {
            c.DominantProcess.Should().Be(MutationalProcess.Unknown, "no active process ⇒ Unknown dominant (INV-05)");
        }
        else
        {
            c.DominantProcess.Should().Be(c.ActiveProcesses[0].Process,
                "the dominant process is the largest-contribution (first) active process (INV-04)");
        }
    }

    // Randomly shuffles an exposure list (Fisher-Yates) — used to prove the classification is INPUT-ORDER
    // INDEPENDENT (deterministic ranking, not dictionary-iteration garbage).
    private static List<(string, double)> Shuffle(Random rng, IReadOnlyList<(string, double)> items)
    {
        var copy = items.ToList();
        for (int i = copy.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (copy[i], copy[j]) = (copy[j], copy[i]);
        }

        return copy;
    }

    #region ONCO-SIG-004 — positive sanity (documented ranking / assignment / proportions)

    // The documented worked example (§7.1): a clearly dominant APOBEC signal ranks first with the correct
    // aggregated proportion, Aging second, and the sub-cutoff SBS4 (0.05 < 0.06) is dropped.
    [Test]
    public void ClassifyMutationalProcess_WorkedExample_RanksApobecDominant_WithCorrectProportions()
    {
        var exposures = new (string, double)[] { ("SBS2", 50), ("SBS13", 30), ("SBS1", 15), ("SBS4", 5) };

        var result = ClassifyMutationalProcess(exposures);

        AssertWellFormedClassification(result);
        result.DominantProcess.Should().Be(MutationalProcess.Apobec, "APOBEC has the largest contribution (§7.1)");
        result.ActiveProcesses.Should().HaveCount(2, "APOBEC and Aging survive; sub-cutoff SBS4 is dropped (§7.1)");
        result.ActiveProcesses[0].Process.Should().Be(MutationalProcess.Apobec);
        result.ActiveProcesses[0].Contribution.Should().BeApproximately(0.80, 1e-9, "SBS2+SBS13 = (50+30)/100 (§7.1)");
        result.ActiveProcesses[1].Process.Should().Be(MutationalProcess.Aging);
        result.ActiveProcesses[1].Contribution.Should().BeApproximately(0.15, 1e-9, "SBS1 = 15/100 (§7.1)");
    }

    // A single dominant signature ranks/assigns to its sole process with proportion 1.0.
    [Test]
    public void ClassifyMutationalProcess_SingleSignature_AssignsFullProportionToItsProcess()
    {
        var exposures = new (string, double)[] { ("SBS4", 42.0) };

        var result = ClassifyMutationalProcess(exposures);

        AssertWellFormedClassification(result);
        result.DominantProcess.Should().Be(MutationalProcess.TobaccoSmoking);
        result.ActiveProcesses.Should().ContainSingle();
        result.ActiveProcesses[0].Contribution.Should().BeApproximately(1.0, 1e-12,
            "a lone signature carries 100% of the normalised contribution (INV-01)");
    }

    #endregion

    #region ONCO-SIG-004 — BE: tied exposures (deterministic enum tie-break, order-independent)

    // Two processes with EXACTLY equal aggregated contribution: APOBEC (SBS2) vs Aging (SBS1), each 0.5.
    // The documented tie-break is ascending process enum ⇒ Aging (enum 1) beats APOBEC (enum 2), so Aging
    // is the deterministic dominant — NOT whichever the dictionary happened to enumerate first.
    [Test]
    public void ClassifyMutationalProcess_TiedExposures_TieBreaksToLowestEnum_Deterministic()
    {
        var exposures = new (string, double)[] { ("SBS2", 10.0), ("SBS1", 10.0) };

        var result = ClassifyMutationalProcess(exposures);

        AssertWellFormedClassification(result);
        result.ActiveProcesses.Should().HaveCount(2);
        result.ActiveProcesses[0].Contribution.Should().BeApproximately(0.5, 1e-12);
        result.ActiveProcesses[1].Contribution.Should().BeApproximately(0.5, 1e-12);
        result.DominantProcess.Should().Be(MutationalProcess.Aging,
            "on an exact tie the lowest process enum (Aging=1) wins, not APOBEC=2 (deterministic §4.1)");
        ((int)result.ActiveProcesses[0].Process).Should().BeLessThan((int)result.ActiveProcesses[1].Process,
            "tied active processes are emitted in ascending enum order");
    }

    // The tie-break must be INPUT-ORDER INDEPENDENT. The same tied exposure set, shuffled many ways, must
    // ALWAYS yield byte-for-byte the same active-set order and the same dominant — never order-dependent
    // garbage from dictionary iteration.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyMutationalProcess_TiedExposures_OrderIndependent_DeterministicAcrossShuffles()
    {
        // Three-way exact tie across Aging, APOBEC and Tobacco — each exactly 1/3.
        var baseExposures = new (string, double)[] { ("SBS1", 1.0), ("SBS2", 1.0), ("SBS4", 1.0) };
        var reference = ClassifyMutationalProcess(baseExposures);

        AssertWellFormedClassification(reference);
        reference.DominantProcess.Should().Be(MutationalProcess.Aging,
            "the three-way tie resolves to the lowest enum, Aging (§4.1)");

        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            List<(string, double)> permuted = Shuffle(rng, baseExposures);

            var result = ClassifyMutationalProcess(permuted);

            AssertWellFormedClassification(result);
            result.DominantProcess.Should().Be(reference.DominantProcess,
                $"seed {seed}: tie-break is input-order independent (deterministic ranking)");
            result.ActiveProcesses.Select(a => a.Process)
                .Should().Equal(reference.ActiveProcesses.Select(a => a.Process),
                    $"seed {seed}: the full active-set order is reproduced regardless of input order");
        }
    }

    // Generalised fuzz: build a vector of equal exposures across a RANDOM SUBSET of the mapped processes
    // (an N-way exact tie), shuffle the input, and assert (a) every contribution is exactly 1/N, (b) the
    // active set is in ascending enum order, and (c) the dominant is the lowest enum in the subset — for
    // any subset, any input order.
    [Test]
    [CancelAfter(30_000)]
    public void ClassifyMutationalProcess_NWayTie_Fuzz_AlwaysLowestEnumDominant_AscendingOrder()
    {
        for (int seed = 0; seed < 600; seed++)
        {
            var rng = new Random(seed);

            // Pick a random non-empty subset of the five mapped processes.
            var subset = MappedSignatures.Where(_ => rng.Next(2) == 0).ToList();
            if (subset.Count == 0)
            {
                subset.Add(MappedSignatures[rng.Next(MappedSignatures.Length)]);
            }

            // All equal exposures (a common random magnitude) ⇒ an exact N-way tie at 1/N each.
            double magnitude = rng.NextDouble() * 100.0 + 1.0;
            var exposures = subset.Select(s => (s.Label, magnitude)).ToList();
            List<(string, double)> permuted = Shuffle(rng, exposures);

            var result = ClassifyMutationalProcess(permuted);

            AssertWellFormedClassification(result);
            int n = subset.Count;
            result.ActiveProcesses.Should().HaveCount(n, $"seed {seed}: each distinct process is active");
            foreach (ProcessActivity a in result.ActiveProcesses)
            {
                a.Contribution.Should().BeApproximately(1.0 / n, 1e-12,
                    $"seed {seed}: an N-way exact tie gives each process exactly 1/N");
            }

            MutationalProcess expectedDominant = subset.Select(s => s.Process).Min();
            result.DominantProcess.Should().Be(expectedDominant,
                $"seed {seed}: the lowest enum in the tied subset is the deterministic dominant");
            result.ActiveProcesses.Select(a => (int)a.Process)
                .Should().BeInAscendingOrder($"seed {seed}: an all-tied active set is in ascending enum order");
        }
    }

    #endregion

    #region ONCO-SIG-004 — BE: all-zero / empty exposures (no ÷0, no NaN, Unknown dominant)

    // The all-zero exposure vector over MAPPED labels: total T = 0 ⇒ the wᵢ = eᵢ/T normalisation is 0/0.
    // The documented result is an EMPTY active set and an Unknown dominant — never a DivideByZero or a NaN
    // contribution (INV-05, §6.1).
    [Test]
    public void ClassifyMutationalProcess_AllZeroExposures_EmptyActiveSet_UnknownDominant_NoNaN()
    {
        var exposures = MappedSignatures.Select(s => (s.Label, 0.0)).ToArray();

        MutationalProcessClassification result = default;
        FluentActions.Invoking(() => result = ClassifyMutationalProcess(exposures))
            .Should().NotThrow("an all-zero exposure vector is a valid degenerate input — no DivideByZero (§6.1)");

        AssertWellFormedClassification(result);
        result.ActiveProcesses.Should().BeEmpty("T = 0 ⇒ no active processes; normalisation is undefined (INV-05)");
        result.DominantProcess.Should().Be(MutationalProcess.Unknown,
            "a zero total yields the Unknown dominant, never a NaN-contribution process (§6.1)");
    }

    // The empty exposure list (degenerate boundary): T = 0 likewise ⇒ empty active set, Unknown dominant.
    [Test]
    public void ClassifyMutationalProcess_EmptyExposureList_EmptyActiveSet_UnknownDominant()
    {
        var result = ClassifyMutationalProcess(Array.Empty<(string, double)>());

        AssertWellFormedClassification(result);
        result.ActiveProcesses.Should().BeEmpty("an empty list has total 0 (INV-05)");
        result.DominantProcess.Should().Be(MutationalProcess.Unknown);
    }

    // Fuzz: all-zero exposure vectors of RANDOM length over random (mapped, unmapped, mixed-case) labels are
    // ALWAYS the empty/Unknown result with no ÷0 and no NaN, whatever the labels or the vector length.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyMutationalProcess_AllZeroExposures_RandomLabels_Fuzz_AlwaysEmptyNoNaN()
    {
        string[] labelPool = { "SBS1", "SBS2", "SBS4", "sbs7a", "SBS6", "SBS13", "SBS99", "GARBAGE", "" };

        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 12);
            var exposures = new (string, double)[n];
            for (int i = 0; i < n; i++)
            {
                exposures[i] = (labelPool[rng.Next(labelPool.Length)], 0.0); // every exposure exactly zero
            }

            MutationalProcessClassification result = default;
            FluentActions.Invoking(() => result = ClassifyMutationalProcess(exposures))
                .Should().NotThrow($"seed {seed}: an all-zero vector is a valid degenerate input (no ÷0)");

            AssertWellFormedClassification(result);
            result.ActiveProcesses.Should().BeEmpty($"seed {seed}: T = 0 ⇒ empty active set (INV-05)");
            result.DominantProcess.Should().Be(MutationalProcess.Unknown,
                $"seed {seed}: a zero total ⇒ Unknown dominant, never a NaN (§6.1)");
        }
    }

    // BE corner: a single tiny-but-positive exposure surrounded by zeros. The total is positive (so
    // normalisation is DEFINED), the lone signature normalises to exactly 1.0 (≥ cutoff), and it is the sole
    // active/dominant process — proving the "zero-total ⇒ empty" branch is gated on T == 0, not on "looks
    // small". No NaN even for a denormal-scale exposure.
    [Test]
    public void ClassifyMutationalProcess_SingleTinyPositiveAmongZeros_NormalisesToOne_NoNaN()
    {
        var exposures = new (string, double)[] { ("SBS1", 0.0), ("SBS2", 1e-300), ("SBS4", 0.0) };

        var result = ClassifyMutationalProcess(exposures);

        AssertWellFormedClassification(result);
        result.ActiveProcesses.Should().ContainSingle("only the positive exposure contributes; T > 0");
        result.DominantProcess.Should().Be(MutationalProcess.Apobec);
        result.ActiveProcesses[0].Contribution.Should().BeApproximately(1.0, 1e-12,
            "the lone positive exposure normalises to 1.0 (e/e), never NaN");
    }

    #endregion

    #region ONCO-SIG-004 — BE: combined / general robustness sweep

    // Broad fuzz over arbitrary non-negative exposure vectors (mapped + unmapped + mixed-case labels, random
    // magnitudes including many exact zeros and ties), each fed in a RANDOM input order. Every classification
    // must be well-formed (finite contributions, INV-01 proportions, documented ordering, dominant = first /
    // Unknown) AND identical to the same input reshuffled — the core BE determinism+no-NaN robustness sweep.
    [Test]
    [CancelAfter(30_000)]
    public void ClassifyMutationalProcess_ArbitraryExposures_Fuzz_WellFormed_AndOrderIndependent()
    {
        string[] labelPool = { "SBS1", "SBS5", "SBS2", "SBS13", "SBS4", "sbs7a", "SBS7b", "SBS6", "SBS20", "SBS999", "junk" };

        for (int seed = 0; seed < 800; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 14);
            var exposures = new (string, double)[n];
            for (int i = 0; i < n; i++)
            {
                // Mostly small integer magnitudes (so exact ties and exact zeros occur often), with a 1-in-4
                // chance of an exact zero — the BE regime.
                double mag = rng.Next(4) == 0 ? 0.0 : rng.Next(0, 6);
                exposures[i] = (labelPool[rng.Next(labelPool.Length)], mag);
            }

            MutationalProcessClassification result = default;
            FluentActions.Invoking(() => result = ClassifyMutationalProcess(exposures))
                .Should().NotThrow($"seed {seed}: any non-negative exposure vector is a valid input");

            AssertWellFormedClassification(result);

            // Re-run on a shuffled copy: the documented ranking is input-order independent.
            List<(string, double)> permuted = Shuffle(rng, exposures);
            var permResult = ClassifyMutationalProcess(permuted);
            AssertWellFormedClassification(permResult);
            permResult.DominantProcess.Should().Be(result.DominantProcess,
                $"seed {seed}: the dominant process is input-order independent");
            permResult.ActiveProcesses.Select(a => a.Process)
                .Should().Equal(result.ActiveProcesses.Select(a => a.Process),
                    $"seed {seed}: the active-set order is reproduced regardless of input order");
            for (int i = 0; i < result.ActiveProcesses.Count; i++)
            {
                permResult.ActiveProcesses[i].Contribution.Should().BeApproximately(
                    result.ActiveProcesses[i].Contribution, 1e-12,
                    $"seed {seed}: per-process contributions are input-order independent");
            }
        }
    }

    #endregion

    #region ONCO-SIG-004 — BE: documented validation guards (null / negative / NaN / out-of-range)

    // null exposures ⇒ ArgumentNullException (§3.3).
    [Test]
    public void ClassifyMutationalProcess_NullExposures_ThrowsArgumentNull()
    {
        FluentActions.Invoking(() => ClassifyMutationalProcess(null!))
            .Should().Throw<ArgumentNullException>("a null exposure list is a documented guard (§3.3)");
    }

    // A negative or NaN exposure ⇒ ArgumentException; an out-of-range / NaN cutoff ⇒ ArgumentOutOfRangeException.
    // Fuzz: random invalid inputs always raise the documented Argument* family, never an undisciplined crash.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyMutationalProcess_RandomInvalidInput_AlwaysDocumentedArgumentException()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);
            int defect = rng.Next(4);
            Action act;
            switch (defect)
            {
                case 0: // negative exposure
                    act = () => ClassifyMutationalProcess(new (string, double)[]
                    {
                        ("SBS1", rng.Next(1, 5)), ("SBS2", -(rng.NextDouble() * 10.0 + 0.001)),
                    });
                    break;

                case 1: // NaN exposure
                    act = () => ClassifyMutationalProcess(new (string, double)[]
                    {
                        ("SBS1", rng.Next(1, 5)), ("SBS2", double.NaN),
                    });
                    break;

                case 2: // cutoff out of [0,1)
                    double badCutoff = rng.Next(2) == 0 ? 1.0 + rng.NextDouble() : -(rng.NextDouble() + 0.001);
                    act = () => ClassifyMutationalProcess(
                        new (string, double)[] { ("SBS1", 10.0) }, badCutoff);
                    break;

                default: // NaN cutoff
                    act = () => ClassifyMutationalProcess(
                        new (string, double)[] { ("SBS1", 10.0) }, double.NaN);
                    break;
            }

            if (defect <= 1)
            {
                FluentActions.Invoking(act).Should().Throw<ArgumentException>(
                    $"seed {seed}: a negative/NaN exposure is a documented ArgumentException (§3.3)");
            }
            else
            {
                FluentActions.Invoking(act).Should().Throw<ArgumentOutOfRangeException>(
                    $"seed {seed}: a cutoff outside [0,1) is a documented ArgumentOutOfRangeException (§3.3)");
            }
        }
    }

    // The cutoff boundary under fuzz: a contribution EXACTLY at the cutoff is retained (strict <); just below
    // is dropped. Build a two-signature vector whose minor signature sits exactly at a random cutoff τ and
    // confirm it survives; nudged below τ it disappears — proving the documented strict-< presence rule.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyMutationalProcess_ContributionAtCutoff_Retained_JustBelow_Dropped_Fuzz()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            double tau = rng.NextDouble() * 0.4 + 0.05; // cutoff in [0.05, 0.45)

            // Minor SBS2 at contribution exactly tau, major SBS1 at (1 - tau): exposures tau and (1-tau).
            var atCutoff = new (string, double)[] { ("SBS2", tau), ("SBS1", 1.0 - tau) };
            var atResult = ClassifyMutationalProcess(atCutoff, tau);
            AssertWellFormedClassification(atResult);
            atResult.ActiveProcesses.Select(a => a.Process).Should().Contain(MutationalProcess.Apobec,
                $"seed {seed}: a contribution exactly at the cutoff is retained (strict <, INV-02)");

            // Nudge the minor signature just below the cutoff ⇒ it is dropped (only Aging survives).
            double below = tau - 1e-6;
            var belowCutoff = new (string, double)[] { ("SBS2", below), ("SBS1", 1.0 - below) };
            var belowResult = ClassifyMutationalProcess(belowCutoff, tau);
            AssertWellFormedClassification(belowResult);
            belowResult.ActiveProcesses.Select(a => a.Process).Should().NotContain(MutationalProcess.Apobec,
                $"seed {seed}: a contribution just below the cutoff is dropped (INV-02)");
        }
    }

    #endregion
}
