using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology clonal-hematopoiesis (CHIP) filtering area — ONCO-CHIP-001.
/// The unit under test is the rule-based CHIP candidate flag + matched-WBC subtraction filter for plasma
/// cfDNA liquid biopsy. Implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs as the public entry points of the
/// <c>Clonal Hematopoiesis (CHIP) Filtering (ONCO-CHIP-001)</c> region:
///   • <see cref="OncologyAnalyzer.IsCanonicalChipGene(string, IReadOnlyCollection{string})"/> —
///       case-insensitive panel membership; null/empty gene ⇒ false;
///   • <see cref="OncologyAnalyzer.IdentifyCHIPVariants(IEnumerable{OncologyAnalyzer.ChipVariant}, IReadOnlyCollection{string}, double)"/> —
///       the gene + VAF candidate-CHIP flag: keep v ⟺ gene(v) ∈ G ∧ VAF(v) ≥ τ (default τ = 0.02), in input order;
///   • <see cref="OncologyAnalyzer.FilterCHIP(IEnumerable{OncologyAnalyzer.ChipVariant}, IEnumerable{OncologyAnalyzer.ChipVariant}, IReadOnlyCollection{string}, double, int)"/> —
///       drop v ⟺ (locus(v) ∈ matched-WBC loci) OR (gene+VAF heuristic), retain the complement, in input order.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts that the code NEVER fails
/// in an undisciplined way: no hang, no NullReferenceException on empty/degenerate collections, no nonsense
/// output. Every input must resolve to EITHER a well-defined, theory-correct CHIP result OR a *documented,
/// intentional* outcome (an <see cref="ArgumentNullException"/> for a null variant/WBC sequence, an
/// <see cref="ArgumentOutOfRangeException"/> for minVaf ∉ (0, 1] or minWbcAltReads &lt; 1). Because the CHIP
/// flag is a pure set/threshold rule the headline hazards are:
///   • a NullReferenceException / crash on an EMPTY gene panel or EMPTY variant list — empty is a DOCUMENTED
///     valid input (§3.3 "Empty inputs are valid (empty result …)"), so an empty panel ⇒ nothing flags and
///     an empty variant list ⇒ empty result, never a throw;
///   • an off-by-one at the VAF band edge — the threshold is inclusive at τ (≥ 0.02), so VAF EXACTLY τ MUST
///     flag (INV-02, §6.1 "VAF exactly 0.02 ⇒ flagged CHIP") and the next-representable double BELOW τ MUST
///     NOT (no false call just under the edge);
///   • a false call OUTSIDE the rule — a non-panel gene at any VAF, or a panel gene below τ, must NOT flag
///     (driver-gene requirement [1]); the result is always a SUBSET of the input in INPUT ORDER (INV-03);
///   • a wrong all-CHIP fold — when every variant is an in-band driver, IdentifyCHIPVariants returns ALL of
///     them and FilterCHIP retains NONE (heuristic removes them), deterministically, with no overflow.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CHIP-001 — clonal-hematopoiesis (CHIP) filtering (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 113.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 113): "empty gene list, VAF at band edges, all-CHIP".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Mapping of the BE targets onto the documented contract:
///   • "empty gene list" ⇒ a panel with ZERO driver genes (caller passes an empty IReadOnlyCollection). The
///       panel is gene-agnostic and empty inputs are valid (§3.3): NO gene can be in an empty panel, so
///       IsCanonicalChipGene ⇒ false for every gene, IdentifyCHIPVariants ⇒ empty, and FilterCHIP keeps
///       every cfDNA variant absent from matched WBC. No NullReferenceException on the empty collection.
///       (Also covers the dual empty: empty VARIANT list ⇒ empty result, no crash.)
///   • "VAF at band edges" ⇒ VAF EXACTLY at the lower band cutoff τ = 0.02 ⇒ flagged CHIP (inclusive ≥,
///       INV-02, §6.1); the next-representable double just BELOW τ ⇒ NOT flagged (no off-by-one). The CHIP
///       rule has a SINGLE documented edge — a one-sided inclusive lower threshold (§2.2 "VAF(v) ≥ τ"); there
///       is NO documented upper band edge, so a panel gene at VAF = 1.0 still flags (a high-VAF CH clone), and
///       this file asserts that one-sidedness explicitly rather than inventing an upper cutoff.
///   • "all-CHIP" ⇒ every variant is a driver gene at an in-band VAF (≥ τ). IdentifyCHIPVariants returns ALL
///       n of them in input order (INV-01); FilterCHIP with empty WBC retains NONE (rule (b) removes them);
///       deterministic across repeated calls; no overflow on a large all-CHIP batch.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Clonal_Hematopoiesis_Filtering.md (docs/algorithms/Oncology/Clonal_Hematopoiesis_Filtering.md):
///   • CHIP candidate flag: isCHIP(v) ⟺ gene(v) ∈ G ∧ VAF(v) ≥ τ, τ = 0.02, threshold INCLUSIVE (§2.2, INV-01, INV-02).
///   • Canonical default panel G = {DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D} (§2.2, §4.2); caller may override.
///   • Gene comparison case-insensitive; null/empty gene is NOT a CHIP gene (§3.3, §6.1).
///   • FilterCHIP: remove v ⟺ locus(v) ∈ matched-WBC loci (≥ minWbcAltReads alt reads, default 1) OR isCHIP(v);
///       retain the complement; output ⊆ input, input order preserved (§2.2, INV-03, INV-04, INV-05).
///   • null variants / null WBC ⇒ ArgumentNullException; minVaf ∉ (0, 1] ⇒ ArgumentOutOfRangeException;
///       minWbcAltReads &lt; 1 ⇒ ArgumentOutOfRangeException (§3.3). Empty inputs are valid (§3.3, §6.1).
///   • Edge cases (§6.1): VAF exactly 0.02 ⇒ flagged; non-CHIP gene high VAF ⇒ not CHIP; variant in matched
///       WBC (any gene) ⇒ removed; empty WBC ⇒ heuristic-only; null/empty gene ⇒ not CHIP.
///   • Worked example (§7.1): cfDNA {DNMT3A@0.06 (driver), EGFR@0.30 (tumour)}, empty WBC ⇒ FilterCHIP ⇒ {EGFR}.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyChipFuzzTests
{
    // ── Well-formed-CHIP-result assertion helper ─────────────────────────────
    // Pins the documented structural contract on EVERY accepted CHIP result, whether from
    // IdentifyCHIPVariants (candidate flag) or FilterCHIP (subtraction):
    //   • the result is a SUBSET of the input (no fabricated variants) (INV-03);
    //   • input ORDER is preserved among the kept variants (INV-03);
    //   • no NULL elements / no count overflow (count is in [0, n]).
    // This is what stops a fuzz test from rubber-stamping a result that re-orders, duplicates, or invents
    // variants the rule never saw.
    private static void AssertSubsetInInputOrder(
        IReadOnlyList<ChipVariant> input, IReadOnlyList<ChipVariant> result)
    {
        result.Should().NotBeNull();
        result.Count.Should().BeInRange(0, input.Count, "result ⊆ input ⇒ |result| ∈ [0, n] (INV-03)");

        // Every kept variant came from the input multiset.
        foreach (ChipVariant v in result)
        {
            input.Should().Contain(v, "the result must be a SUBSET of the input — no fabricated variants (INV-03)");
        }

        // Input order is preserved: the result is a subsequence of the input.
        int i = 0;
        foreach (ChipVariant v in input)
        {
            if (i < result.Count && result[i].Equals(v))
            {
                i++;
            }
        }
        i.Should().Be(result.Count,
            "the result must appear in INPUT ORDER (a subsequence of the input) (INV-03)");
    }

    // The documented canonical default panel (Steensma 2015 / Genovese 2014). Read from the source so the
    // tests use REAL gene symbols, not invented ones.
    private static readonly string[] CanonicalGenes =
        DefaultChipGenes.ToArray();

    // A non-CHIP gene (a classic solid-tumour driver, deliberately NOT in the canonical panel).
    private const string NonChipGene = "EGFR";

    // Build a cfDNA variant. Locus is varied per index so loci stay distinct unless we want a match.
    private static ChipVariant Variant(string gene, double vaf, int index = 0, int altReads = 0)
        => new(
            Chromosome: "chr" + ((index % 22) + 1),
            Position: 1_000 + index,
            ReferenceAllele: "C",
            AlternateAllele: "T",
            Gene: gene,
            Vaf: vaf,
            AltReads: altReads);

    #region ONCO-CHIP-001 — Clonal hematopoiesis (CHIP) filtering

    // ── BE: empty gene list ──────────────────────────────────────────────────
    // An empty driver panel is a DOCUMENTED valid input (the algorithm is gene-panel-agnostic, §3.3).
    // No gene can be a member of the empty set, so nothing flags — and the empty collection must NOT
    // trigger a NullReferenceException / crash on iteration.

    [Test]
    public void IsCanonicalChipGene_EmptyPanel_NeverMatches()
    {
        var empty = Array.Empty<string>();
        foreach (string gene in CanonicalGenes.Append(NonChipGene).Append(""))
        {
            IsCanonicalChipGene(gene, empty).Should().BeFalse(
                "no gene is a member of an EMPTY panel — empty list must not crash or match (§3.3)");
        }
    }

    [Test]
    public void IdentifyCHIPVariants_EmptyPanel_FlagsNothing_NoCrash()
    {
        var rng = new Random(11301);
        var variants = new List<ChipVariant>();
        for (int i = 0; i < 200; i++)
        {
            // Every variant is a canonical driver at a comfortably in-band VAF; only the empty panel
            // should suppress them.
            string gene = CanonicalGenes[rng.Next(CanonicalGenes.Length)];
            variants.Add(Variant(gene, 0.02 + rng.NextDouble() * 0.5, i));
        }

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(variants, Array.Empty<string>());

        flagged.Should().BeEmpty(
            "an EMPTY panel ⇒ no gene is a driver ⇒ NOTHING flags CHIP, no NullReference on the empty list (§3.3)");
        AssertSubsetInInputOrder(variants, flagged);
    }

    [Test]
    public void FilterCHIP_EmptyPanel_NoWbc_RetainsEveryVariant()
    {
        var rng = new Random(11302);
        var variants = new List<ChipVariant>();
        for (int i = 0; i < 150; i++)
        {
            string gene = CanonicalGenes[rng.Next(CanonicalGenes.Length)];
            variants.Add(Variant(gene, 0.05 + rng.NextDouble() * 0.4, i));
        }

        // Empty panel ⇒ heuristic (b) can never fire; empty WBC ⇒ subtraction (a) can never fire.
        IReadOnlyList<ChipVariant> retained = FilterCHIP(
            variants, Array.Empty<ChipVariant>(), chipGenes: Array.Empty<string>());

        retained.Should().BeEquivalentTo(variants, o => o.WithStrictOrdering(),
            "with an empty panel AND empty WBC, NOTHING is removed — every cfDNA variant is retained (§2.2, INV-05)");
        AssertSubsetInInputOrder(variants, retained);
    }

    [Test]
    public void EmptyVariantList_IsValid_EmptyResult_NoCrash()
    {
        var empty = Array.Empty<ChipVariant>();

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(empty);
        IReadOnlyList<ChipVariant> retained = FilterCHIP(empty, empty);

        flagged.Should().BeEmpty("empty variant input ⇒ empty result, never a crash (§3.3)");
        retained.Should().BeEmpty("empty variant input ⇒ empty result, never a crash (§3.3)");
    }

    // ── BE: VAF at band edges ────────────────────────────────────────────────
    // The CHIP rule has ONE documented edge: an inclusive lower threshold τ = 0.02. VAF EXACTLY τ flags;
    // the next-representable double just below τ does not. There is NO documented upper edge (a high-VAF
    // CH clone still flags), which this section asserts explicitly.

    [Test]
    public void IdentifyCHIPVariants_VafExactlyAtThreshold_IsFlagged_Inclusive()
    {
        // DNMT3A is the most prevalent canonical CHIP driver (§2.2).
        var atEdge = Variant("DNMT3A", OncologyAnalyzer.ChipVafThreshold, index: 0);

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(new[] { atEdge });

        flagged.Should().ContainSingle().Which.Should().Be(atEdge,
            "VAF EXACTLY 0.02 is inclusive ⇒ flagged CHIP (INV-02, §6.1 'VAF exactly 0.02 ⇒ flagged')");
    }

    [Test]
    public void IdentifyCHIPVariants_VafJustBelowThreshold_IsNotFlagged_NoOffByOne()
    {
        double justBelow = BitDecrement(OncologyAnalyzer.ChipVafThreshold);
        justBelow.Should().BeLessThan(OncologyAnalyzer.ChipVafThreshold);

        var belowEdge = Variant("TET2", justBelow, index: 0);

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(new[] { belowEdge });

        flagged.Should().BeEmpty(
            "the next double BELOW 0.02 is under the inclusive cutoff ⇒ NOT flagged (no off-by-one at the lower edge)");
    }

    [Test]
    public void IdentifyCHIPVariants_PanelGeneAtVafOne_StillFlags_NoUpperBand()
    {
        // The documented rule is one-sided (≥ τ only); a driver at VAF = 1.0 must still flag.
        var top = Variant("ASXL1", 1.0, index: 0);

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(new[] { top });

        flagged.Should().ContainSingle().Which.Should().Be(top,
            "the CHIP VAF rule is a ONE-SIDED inclusive lower threshold — there is no upper band edge (§2.2)");
    }

    [Test]
    public void IdentifyCHIPVariants_CustomThreshold_BandEdgesRespected()
    {
        var rng = new Random(11303);
        // Sweep custom thresholds across the valid (0, 1] range; for each, the variant exactly at τ flags
        // and the variant just below τ does not.
        for (int t = 0; t < 60; t++)
        {
            double tau = 0.0009765625 + rng.NextDouble() * (1.0 - 0.0009765625);
            string gene = CanonicalGenes[rng.Next(CanonicalGenes.Length)];

            var atEdge = Variant(gene, tau, index: 1);
            var belowEdge = Variant(gene, BitDecrement(tau), index: 2);

            IdentifyCHIPVariants(new[] { atEdge }, minVaf: tau)
                .Should().ContainSingle("VAF = τ is inclusive at any custom τ (INV-02)");
            IdentifyCHIPVariants(new[] { belowEdge }, minVaf: tau)
                .Should().BeEmpty("VAF just below τ never flags at any custom τ (no off-by-one)");
        }
    }

    [Test]
    public void FilterCHIP_PanelGeneAtThreshold_IsRemoved_BelowThreshold_IsRetained()
    {
        var inBand = Variant("JAK2", OncologyAnalyzer.ChipVafThreshold, index: 0);
        var outBand = Variant("JAK2", BitDecrement(OncologyAnalyzer.ChipVafThreshold), index: 1);

        IReadOnlyList<ChipVariant> retained = FilterCHIP(
            new[] { inBand, outBand }, Array.Empty<ChipVariant>());

        retained.Should().ContainSingle().Which.Should().Be(outBand,
            "driver at VAF = τ meets the heuristic ⇒ removed; just below τ ⇒ retained as candidate tumour (§2.2)");
    }

    // ── BE: all-CHIP ─────────────────────────────────────────────────────────
    // Every variant is a canonical driver at an in-band VAF. IdentifyCHIPVariants returns ALL of them in
    // input order; FilterCHIP (empty WBC) retains NONE; both are deterministic and overflow-free at scale.

    [Test]
    [CancelAfter(30_000)]
    public void IdentifyCHIPVariants_AllChip_FlagsEverything_InOrder()
    {
        var rng = new Random(11304);
        var variants = new List<ChipVariant>();
        for (int i = 0; i < 5_000; i++)
        {
            string gene = CanonicalGenes[rng.Next(CanonicalGenes.Length)];
            // In-band: VAF in [τ, 1].
            double vaf = OncologyAnalyzer.ChipVafThreshold
                + rng.NextDouble() * (1.0 - OncologyAnalyzer.ChipVafThreshold);
            variants.Add(Variant(gene, vaf, i));
        }

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(variants);

        flagged.Should().HaveCount(variants.Count,
            "every variant is an in-band canonical driver ⇒ ALL flag CHIP (INV-01)");
        flagged.Should().BeEquivalentTo(variants, o => o.WithStrictOrdering(),
            "all-CHIP ⇒ the output is the input verbatim, in input order (INV-01, INV-03)");
        AssertSubsetInInputOrder(variants, flagged);
    }

    [Test]
    [CancelAfter(30_000)]
    public void FilterCHIP_AllChip_NoWbc_RemovesEverything()
    {
        var rng = new Random(11305);
        var variants = new List<ChipVariant>();
        for (int i = 0; i < 5_000; i++)
        {
            string gene = CanonicalGenes[rng.Next(CanonicalGenes.Length)];
            double vaf = OncologyAnalyzer.ChipVafThreshold
                + rng.NextDouble() * (1.0 - OncologyAnalyzer.ChipVafThreshold);
            variants.Add(Variant(gene, vaf, i));
        }

        IReadOnlyList<ChipVariant> retained = FilterCHIP(variants, Array.Empty<ChipVariant>());

        retained.Should().BeEmpty(
            "all-CHIP with empty WBC ⇒ rule (b) removes every variant ⇒ no candidate tumour remains (§2.2)");
    }

    [Test]
    public void AllChip_Deterministic_AcrossRepeatedCalls()
    {
        var rng = new Random(11306);
        var variants = new List<ChipVariant>();
        for (int i = 0; i < 300; i++)
        {
            string gene = CanonicalGenes[rng.Next(CanonicalGenes.Length)];
            variants.Add(Variant(gene, 0.02 + rng.NextDouble() * 0.9, i));
        }

        IReadOnlyList<ChipVariant> a = IdentifyCHIPVariants(variants);
        IReadOnlyList<ChipVariant> b = IdentifyCHIPVariants(variants);

        b.Should().BeEquivalentTo(a, o => o.WithStrictOrdering(),
            "the CHIP flag is a deterministic pure rule — identical input ⇒ identical output");
    }

    // ── Positive sanity / business contract (not rubber-stamp green) ──────────
    // A driver gene in-band flags; the same gene out-of-band does not; a non-CHIP gene in-band does not.

    [Test]
    public void Sanity_DriverInBand_Flags_OutOfBand_And_NonChip_DoNot()
    {
        var driverInBand = Variant("DNMT3A", 0.06, index: 0);                 // CH driver, in band ⇒ CHIP
        var driverOutOfBand = Variant("DNMT3A", 0.005, index: 1);             // CH driver, below τ ⇒ not CHIP
        var nonChipInBand = Variant(NonChipGene, 0.30, index: 2);            // solid-tumour driver, high VAF ⇒ not CHIP

        var batch = new[] { driverInBand, driverOutOfBand, nonChipInBand };

        IReadOnlyList<ChipVariant> flagged = IdentifyCHIPVariants(batch);

        flagged.Should().ContainSingle().Which.Should().Be(driverInBand,
            "only a driver gene at VAF ≥ τ is CHIP; below-band driver and high-VAF non-driver are not (INV-01)");
    }

    [Test]
    public void Sanity_FilterCHIP_WorkedExample_KeepsTumourCandidate()
    {
        // §7.1 worked example: cfDNA {DNMT3A@0.06 (CH driver), EGFR@0.30 (tumour)}, empty WBC ⇒ {EGFR}.
        var cfDna = new[]
        {
            new ChipVariant("2", 25234374, "C", "T", "DNMT3A", 0.06),
            new ChipVariant("7", 55259515, "T", "G", NonChipGene, 0.30),
        };

        IReadOnlyList<ChipVariant> tumourCandidates =
            FilterCHIP(cfDna, Array.Empty<ChipVariant>());

        tumourCandidates.Should().ContainSingle().Which.Gene.Should().Be(NonChipGene,
            "the DNMT3A CH driver is filtered out; only the EGFR tumour candidate is retained (§7.1)");
    }

    [Test]
    public void Sanity_MatchedWbc_RemovesVariant_RegardlessOfGene()
    {
        // §6.1: a non-CHIP gene present in matched WBC is removed (matched-WBC origin test, rule (a)).
        var cf = new ChipVariant("7", 55259515, "T", "G", NonChipGene, 0.30, AltReads: 0);
        var wbc = new ChipVariant("7", 55259515, "T", "G", NonChipGene, 0.40, AltReads: 5);

        IReadOnlyList<ChipVariant> retained = FilterCHIP(new[] { cf }, new[] { wbc });

        retained.Should().BeEmpty(
            "a cfDNA locus present in matched WBC (≥ minWbcAltReads) is removed regardless of gene (INV-04, rule (a))");
    }

    [Test]
    public void IsCanonicalChipGene_NullOrEmptyGene_IsNeverChip()
    {
        IsCanonicalChipGene(null).Should().BeFalse("null gene is not a CHIP gene (§3.3)");
        IsCanonicalChipGene("").Should().BeFalse("empty gene is not a CHIP gene (§3.3)");
    }

    [Test]
    public void IsCanonicalChipGene_CaseInsensitive_OnCanonicalGenes()
    {
        foreach (string gene in CanonicalGenes)
        {
            IsCanonicalChipGene(gene.ToLowerInvariant()).Should().BeTrue(
                "panel membership is case-insensitive (HGNC symbols) (§3.3)");
        }
    }

    // ── Documented guard rails (degenerate parameter inputs must throw, not crash) ─

    [Test]
    public void NullVariants_Throws_ArgumentNullException()
    {
        Action identify = () => IdentifyCHIPVariants(null!);
        Action filterA = () => FilterCHIP(null!, Array.Empty<ChipVariant>());
        Action filterB = () => FilterCHIP(Array.Empty<ChipVariant>(), null!);

        identify.Should().Throw<ArgumentNullException>("null variants ⇒ ArgumentNullException (§3.3)");
        filterA.Should().Throw<ArgumentNullException>("null variants ⇒ ArgumentNullException (§3.3)");
        filterB.Should().Throw<ArgumentNullException>("null WBC ⇒ ArgumentNullException (§3.3)");
    }

    [Test]
    public void MinVafOutOfRange_Throws_ArgumentOutOfRangeException()
    {
        var rng = new Random(11307);
        foreach (double bad in new[] { 0.0, -0.0001, -1.0, 1.0000001, 2.0, double.NaN })
        {
            Action identify = () => IdentifyCHIPVariants(Array.Empty<ChipVariant>(), minVaf: bad);
            Action filter = () => FilterCHIP(
                Array.Empty<ChipVariant>(), Array.Empty<ChipVariant>(), minVaf: bad);
            identify.Should().Throw<ArgumentOutOfRangeException>(
                "minVaf ∉ (0, 1] ⇒ ArgumentOutOfRangeException (§3.3), bad = {0}", bad);
            filter.Should().Throw<ArgumentOutOfRangeException>(
                "minVaf ∉ (0, 1] ⇒ ArgumentOutOfRangeException (§3.3), bad = {0}", bad);
        }

        // The boundary value 1.0 is INSIDE (0, 1] and must NOT throw.
        Action atUpper = () => IdentifyCHIPVariants(Array.Empty<ChipVariant>(), minVaf: 1.0);
        atUpper.Should().NotThrow("minVaf = 1.0 is the inclusive upper bound of (0, 1] (§3.3)");

        // A random in-range minVaf never throws.
        double inRange = BitDecrement(1.0) * rng.NextDouble() + double.Epsilon;
        Action ok = () => IdentifyCHIPVariants(Array.Empty<ChipVariant>(), minVaf: Math.Clamp(inRange, double.Epsilon, 1.0));
        ok.Should().NotThrow();
    }

    [Test]
    public void MinWbcAltReadsBelowOne_Throws_ArgumentOutOfRangeException()
    {
        foreach (int bad in new[] { 0, -1, int.MinValue })
        {
            Action filter = () => FilterCHIP(
                Array.Empty<ChipVariant>(), Array.Empty<ChipVariant>(), minWbcAltReads: bad);
            filter.Should().Throw<ArgumentOutOfRangeException>(
                "minWbcAltReads < 1 ⇒ ArgumentOutOfRangeException (§3.3), bad = {0}", bad);
        }
    }

    #endregion

    // Next-representable double below x (no upper edge needed here; we only step below thresholds).
    private static double BitDecrement(double x) => Math.BitDecrement(x);
}
