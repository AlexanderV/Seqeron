using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the MiRNA area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-SEED-001 — miRNA seed region / target seed matching (MiRNA).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 74.
///
/// API under test (MiRnaAnalyzer.GetSeedSequence / CreateMiRna / FindTargetSites):
///   The seed is miRNA nucleotides 2–8 (0-based substring [1,8)). Target prediction scans an
///   mRNA for the reverse complement of the seed (the 6-mer core, optionally extended to
///   7mer/8mer); a seed match is therefore necessary for a predicted site.
///
/// Relations (derived from the seed definition, NOT from output):
///   • INV  (3' end changes don't affect seed): the seed is fixed to positions 2–8, so mutating
///          any nucleotide from position 9 onward leaves the seed unchanged.
///   • INV  (seed extraction deterministic): extraction is a pure, case-normalising function,
///          and CreateMiRna additionally maps T→U, so DNA and RNA spellings give the same seed.
///   • COMP (seed match ⊂ target prediction): a planted seed reverse-complement yields a target
///          site whose span contains the 6-mer seed-core match, and an mRNA with no seed match
///          yields no sites — target prediction is gated on the seed.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MiRnaMetamorphicTests
{
    // let-7a-5p (22 nt). Seed (positions 2–8) = "GAGGUAG".
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    #region MIRNA-SEED-001 INV — 3' end changes leave the seed unchanged

    [Test]
    [Description("INV: the seed is positions 2–8, so mutating any nucleotide from position 9 onward (the 3' end) does not change the extracted seed.")]
    public void Seed_3PrimeEndChanges_DoNotAffectSeed()
    {
        string baseSeed = MiRnaAnalyzer.GetSeedSequence(Let7a);
        baseSeed.Should().Be("GAGGUAG", because: "let-7a nucleotides 2–8 spell GAGGUAG");

        // Keep the first 8 nucleotides (which contain the seed) and rewrite the 3' tail.
        foreach (var tail in new[] { "CCCCCCCCCCCCCC", "AAAAAAAAAAAAAA", "GUGUGUGUGUGUGU" })
        {
            string mutated = Let7a.Substring(0, 8) + tail;
            MiRnaAnalyzer.GetSeedSequence(mutated).Should().Be(baseSeed,
                because: "the seed window (positions 2–8) is untouched by 3'-end edits");
            MiRnaAnalyzer.CreateMiRna("m", mutated).SeedSequence.Should().Be(baseSeed,
                because: "CreateMiRna derives the seed from the same fixed window");
        }
    }

    #endregion

    #region MIRNA-SEED-001 INV — seed extraction is deterministic and normalising

    [Test]
    [Description("INV: seed extraction is a pure, case-normalising function; CreateMiRna also maps T→U so DNA and RNA spellings of the same miRNA yield the same seed.")]
    public void Seed_Extraction_IsDeterministicAndNormalising()
    {
        MiRnaAnalyzer.GetSeedSequence(Let7a)
            .Should().Be(MiRnaAnalyzer.GetSeedSequence(Let7a), because: "extraction is deterministic");

        MiRnaAnalyzer.GetSeedSequence(Let7a.ToLowerInvariant())
            .Should().Be(MiRnaAnalyzer.GetSeedSequence(Let7a), because: "extraction upper-cases its input");

        // DNA spelling (T for U) of let-7a must yield the same RNA seed via CreateMiRna.
        string dna = Let7a.Replace('U', 'T');
        MiRnaAnalyzer.CreateMiRna("dna", dna).SeedSequence
            .Should().Be(MiRnaAnalyzer.CreateMiRna("rna", Let7a).SeedSequence,
                because: "CreateMiRna transcribes T→U before extracting the seed");
    }

    #endregion

    #region MIRNA-SEED-001 COMP — a seed match is necessary and contained in a target site

    [Test]
    [Description("COMP: a planted seed reverse-complement produces a target site whose span contains the 6-mer seed core, and an mRNA with no seed match produces no sites.")]
    public void Seed_TargetPrediction_RequiresAndContainsSeedMatch()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence); // RC of positions 2–8 (7 nt)
        string sixmerCore = seedRc.Substring(1, 6);                              // RC of positions 2–7

        // Plant the full seed-RC with an upstream position-8 match and a downstream 'A'
        // (opposite miRNA position 1) so an 8mer site forms and clears the default score gate.
        string mrna = "GGGGG" + seedRc + "A" + "GGGGG";

        var sites = MiRnaAnalyzer.FindTargetSites(mrna, mirna).ToList();

        sites.Should().NotBeEmpty(because: "a seed reverse-complement match yields a predicted target site");
        sites.Should().OnlyContain(s => mrna.Substring(s.Start, s.End - s.Start + 1).Contains(sixmerCore),
            because: "every predicted site's span contains the 6-mer seed-core match it was built around");

        // No seed match anywhere ⇒ no prediction: target prediction is gated on the seed.
        MiRnaAnalyzer.FindTargetSites(new string('A', 40), mirna, minScore: 0.0)
            .Should().BeEmpty(because: "without a seed reverse-complement match there is nothing to predict");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MIRNA-TARGET-001 — miRNA target-site scoring (MiRNA).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 75.
    //
    // API under test (MiRnaAnalyzer.FindTargetSites):
    //   Classifies each seed match by completeness — 8mer (seed + pos-8 + A1) > 7mer-m8 >
    //   7mer-A1 > 6mer — and assigns a base efficacy score (1.0/0.52/0.32/0.15) adjusted by
    //   the duplex; only sites scoring ≥ minScore are returned. The duplex (and hence score)
    //   is computed from the site and the following miRNA-length window only.
    //
    // Relations (derived from the site-type efficacy ladder, NOT from output):
    //   • MON  (more seed complementarity ⇒ higher score): completing a seed match from 6mer →
    //          7mer-m8 → 8mer (at the same locus) raises the site score.
    //   • SUB  (stringent ⇒ subset): the returned set is filtered by score ≥ minScore, so a
    //          higher threshold yields a subset of the sites found at a lower threshold.
    //   • INV  (distant 3'UTR change ⇒ same score): a site's score depends only on its own
    //          miRNA-length window, so editing the mRNA beyond that window leaves it unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region MIRNA-TARGET-001 MON — more seed complementarity raises the score

    [Test]
    [Description("MON: completing the seed match at one locus from 6mer → 7mer-m8 → 8mer raises the predicted site score.")]
    public void TargetScore_MoreSeedComplementarity_HigherScore()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence);
        string sixmerCore = seedRc.Substring(1, 6);

        // Same seed core, increasing completeness: 6mer (core only), 7mer-m8 (+pos8), 8mer (+A1).
        string mrna6 = "GG" + sixmerCore + "GG";
        string mrna7 = "GGGG" + seedRc + "GGGG";
        string mrna8 = "GGGG" + seedRc + "A" + "GGGG";

        double score6 = MiRnaAnalyzer.FindTargetSites(mrna6, mirna, minScore: 0.0).Max(s => s.Score);
        double score7 = MiRnaAnalyzer.FindTargetSites(mrna7, mirna, minScore: 0.0).Max(s => s.Score);
        double score8 = MiRnaAnalyzer.FindTargetSites(mrna8, mirna, minScore: 0.0).Max(s => s.Score);

        score8.Should().BeGreaterThan(score7, because: "an 8mer (perfect seed + A1) is a stronger site than a 7mer-m8");
        score7.Should().BeGreaterThan(score6, because: "a 7mer-m8 (seed + position 8) is a stronger site than a bare 6mer");
    }

    #endregion

    #region MIRNA-TARGET-001 SUB — a higher score threshold yields a subset of sites

    [Test]
    [Description("SUB: sites are filtered by score ≥ minScore, so raising the threshold returns a subset — here the weak 6mer drops while the strong 8mer survives.")]
    public void TargetSites_StricterThreshold_YieldsSubset()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence);
        string sixmerCore = seedRc.Substring(1, 6);

        // A strong 8mer near the start and a weak isolated 6mer downstream.
        string mrna = "GGGG" + seedRc + "A" + "GGGGGGGGGG" + sixmerCore + "GG";

        var lenient = SiteKeys(MiRnaAnalyzer.FindTargetSites(mrna, mirna, minScore: 0.0));
        var stringent = SiteKeys(MiRnaAnalyzer.FindTargetSites(mrna, mirna, minScore: 0.3));

        stringent.IsSubsetOf(lenient).Should().BeTrue(
            because: "raising minScore can only remove sites, never add them");
        lenient.Count.Should().BeGreaterThan(stringent.Count,
            because: "the 6mer (base efficacy 0.15 < 0.3) is dropped at the stricter threshold while the 8mer survives");
    }

    #endregion

    #region MIRNA-TARGET-001 INV — distant 3'UTR edits leave a site's score unchanged

    [Test]
    [Description("INV: a site's score is computed from its own miRNA-length window, so editing the mRNA beyond that window leaves the site's score (and type/span) unchanged.")]
    public void TargetScore_DistantEdit_DoesNotChangeSite()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence);

        // 8mer site at Start=4; its duplex window spans positions 4..25 (miRNA length 22).
        string baseMrna = "GGGG" + seedRc + "A" + new string('G', 20); // length 32
        var baseSite = MiRnaAnalyzer.FindTargetSites(baseMrna, mirna, minScore: 0.0).Single(s => s.Start == 4);

        // Edit only positions 28..31 — well beyond the site's window.
        char[] chars = baseMrna.ToCharArray();
        for (int i = 28; i < 32; i++) chars[i] = 'A';
        string editedMrna = new string(chars);

        var editedSite = MiRnaAnalyzer.FindTargetSites(editedMrna, mirna, minScore: 0.0).Single(s => s.Start == 4);

        editedSite.Score.Should().Be(baseSite.Score,
            because: "an edit outside the site's miRNA-length window cannot change its score");
        editedSite.Type.Should().Be(baseSite.Type, because: "the site classification is unchanged");
        editedSite.End.Should().Be(baseSite.End, because: "the site span is unchanged");
    }

    #endregion

    #region MIRNA-TARGET-001 — Helpers

    private static System.Collections.Generic.HashSet<(int Start, int End, MiRnaAnalyzer.TargetSiteType Type)> SiteKeys(
        System.Collections.Generic.IEnumerable<MiRnaAnalyzer.TargetSite> sites) =>
        sites.Select(s => (s.Start, s.End, s.Type)).ToHashSet();

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MIRNA-PRECURSOR-001 — pre-miRNA hairpin prediction (MiRNA).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 76.
    //
    // API under test (MiRnaAnalyzer.FindPreMiRnaHairpins):
    //   Treats a window as a pre-miRNA hairpin when its two ends pair into a stem of ≥18 bp
    //   closing a 3–25 nt loop. The dot-bracket structure marks the stem as '('…')' and the
    //   loop as '.', and the folding free energy sums Turner stem-stacking terms (negative)
    //   plus loop/terminal corrections. The stem is detected from the arms (ends) only.
    //
    // Relations (derived from the stem/loop model, NOT from output):
    //   • MON  (extend stem ⇒ more stable): lengthening the complementary stem adds stacking
    //          terms, lowering (more negative) the precursor free energy.
    //   • INV  (loop sequence ⇒ same classification): with the arms and loop length fixed,
    //          changing the loop's interior bases leaves the dot-bracket structure (and the
    //          loop-boundary-only energy) unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    // A pure G:C arm so every added stem step is a strongly stabilising stack.
    private const string GcArm = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"; // 30 nt

    /// <summary>A hairpin = (L-base G:C arm) + loop + reverse-complement of the arm.</summary>
    private static string Hairpin(int armLength, string loop) =>
        GcArm.Substring(0, armLength) + loop + MiRnaAnalyzer.GetReverseComplement(GcArm.Substring(0, armLength));

    /// <summary>The full-length precursor (the designed hairpin spanning the whole sequence).</summary>
    private static MiRnaAnalyzer.PreMiRna FullHairpin(string sequence) =>
        MiRnaAnalyzer.FindPreMiRnaHairpins(sequence).Single(p => p.Start == 0 && p.End == sequence.Length - 1);

    #region MIRNA-PRECURSOR-001 MON — a longer stem yields a more stable precursor

    [Test]
    [Description("MON: lengthening the complementary stem adds stacking terms, so the precursor free energy strictly decreases (more stable).")]
    public void Precursor_LongerStem_MoreStable()
    {
        // Fixed 19-nt loop with non-pairing 'A' boundaries; only the stem length changes.
        string loop = "A" + new string('C', 17) + "A";

        double previousEnergy = double.MaxValue;
        foreach (int armLength in new[] { 18, 20, 22, 24 })
        {
            double energy = FullHairpin(Hairpin(armLength, loop)).FreeEnergy;

            energy.Should().BeLessThan(previousEnergy,
                because: $"a {armLength}-bp stem adds stacking pairs beyond the shorter stem, lowering ΔG");
            previousEnergy = energy;
        }
    }

    #endregion

    #region MIRNA-PRECURSOR-001 INV — loop content does not change the structure classification

    [Test]
    [Description("INV: with the arms and loop length fixed, changing the loop's interior bases leaves the dot-bracket structure (and the loop-boundary-only energy) unchanged.")]
    public void Precursor_LoopContent_DoesNotChangeClassification()
    {
        const int armLength = 20;

        // Same 19-nt loop length and 'A' boundaries; only the interior 17 bases differ.
        var reference = FullHairpin(Hairpin(armLength, "A" + new string('C', 17) + "A"));

        foreach (var interior in new[]
                 {
                     new string('G', 17),
                     new string('U', 17),
                     "ACGUACGUACGUACGUA",
                 })
        {
            var variant = FullHairpin(Hairpin(armLength, "A" + interior + "A"));

            variant.Structure.Should().Be(reference.Structure,
                because: "the stem is detected from the arms only, so loop interior content does not change the dot-bracket");
            variant.FreeEnergy.Should().Be(reference.FreeEnergy,
                because: "stem stacking, loop initiation and the fixed loop-boundary terminal terms are all independent of the loop interior");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: MIRNA-PAIR-001 — miRNA-to-target seed pairing/alignment (MiRNA).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 225.
    //
    // API under test (MiRnaAnalyzer.FindTargetSites):
    //   Locates seed reverse-complement matches along the mRNA and reports the pairing as target
    //   sites (Start/End, type, alignment, score). The site span and score depend only on the local
    //   window at the match.
    //
    // Relations (derived from the local seed-match definition, NOT from output):
    //   • INV  (deterministic): a pure function — repeated calls give identical sites.
    //   • SHIFT (prepend flank shifts alignment): prepending a seed-free 5' flank to the mRNA shifts
    //          every site's Start/End by the flank length and leaves its type/score/alignment intact.
    //          A poly-A flank carries no seed match (target prediction is seed-gated), so it adds no
    //          sites — the relation is an exact translation.
    // ───────────────────────────────────────────────────────────────────────────

    #region MIRNA-PAIR-001 — Helpers

    private static (int Start, int End, MiRnaAnalyzer.TargetSiteType Type, double Score, int SeedMatchLength, string Target, string Alignment)
        PairKey(MiRnaAnalyzer.TargetSite s) => (s.Start, s.End, s.Type, s.Score, s.SeedMatchLength, s.TargetSequence, s.Alignment);

    #endregion

    #region MIRNA-PAIR-001 INV — pairing is deterministic

    [Test]
    [Description("INV: FindTargetSites is a pure function — repeated calls on the same inputs give identical sites.")]
    public void Pairing_SameInput_SameSites()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence);
        string mrna = "GGGG" + seedRc + "GGGGGGGG";

        var first = MiRnaAnalyzer.FindTargetSites(mrna, mirna, minScore: 0.0).Select(PairKey).ToList();
        var again = MiRnaAnalyzer.FindTargetSites(mrna, mirna, minScore: 0.0).Select(PairKey).ToList();

        first.Should().NotBeEmpty(because: "the embedded seed reverse-complement yields at least one pairing — a non-vacuous fixture");
        again.Should().Equal(first, because: "target pairing has no hidden state");
    }

    #endregion

    #region MIRNA-PAIR-001 SHIFT — prepending a seed-free flank translates the pairing

    [Test]
    [Description("SHIFT: prepending a poly-A (seed-free) 5' flank shifts every site's coordinates by the flank length while preserving its type, score and alignment.")]
    public void Pairing_PrependFlank_ShiftsSites()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string seedRc = MiRnaAnalyzer.GetReverseComplement(mirna.SeedSequence);
        string mrna = "GGGG" + seedRc + "GGGGGGGG";

        var baseline = MiRnaAnalyzer.FindTargetSites(mrna, mirna, minScore: 0.0).ToList();
        baseline.Should().NotBeEmpty();
        baseline.Should().OnlyContain(s => s.Start >= 1, because: "no site sits at index 0, so a 5' flank cannot change any site's seed classification");

        foreach (int offset in new[] { 4, 12 })
        {
            string flank = new string('A', offset); // poly-A carries no seed match
            var shifted = MiRnaAnalyzer.FindTargetSites(flank + mrna, mirna, minScore: 0.0).Select(PairKey).ToList();

            var expected = baseline.Select(s =>
                (s.Start + offset, s.End + offset, s.Type, s.Score, s.SeedMatchLength, s.TargetSequence, s.Alignment));

            shifted.Should().BeEquivalentTo(expected,
                because: $"a seed-free {offset}-nt 5' flank adds no sites and translates every pairing by {offset}");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-CONTEXT-001 — TargetScan context++ local-AU scoring (MiRNA)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (Agarwal et al. 2015 eLife 4:e05005; targetscan_70_context_scores.pl;
    //   tests/TestSpecs/MIRNA-CONTEXT-001.md):
    //   ScoreTargetSiteContextPlusPlus sums per-feature contributions of the TargetScan context++
    //   model; one feature is the LOCAL AU content of the 30-nt flanks around the seed match. Its
    //   coefficient is negative, so AU-rich context makes the contribution (and the overall
    //   score) MORE NEGATIVE — a stronger predicted repression. Two metamorphic relations
    //   (checklist row 252):
    //
    //   • MON (stronger local AU context → more negative score): raising the AU content of the
    //     flanks around a fixed seed match makes the local-AU contribution (and the partial context
    //     score) monotonically more negative.
    //   • INV (same site → same score): scoring is deterministic — the same (mRNA, miRNA, site)
    //     yields the identical contribution breakdown every time.
    //
    // API under test: MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus (ContextPlusPlusScore).

    #region MIRNA-CONTEXT-001 — context++ local-AU scoring

    // The canonical 8mer match for let-7a: the reverse complement of miRNA positions 1–8.
    private static string Let7a8merMatch => MiRnaAnalyzer.GetReverseComplement(Let7a.Substring(0, 8));

    private static MiRnaAnalyzer.TargetSite Find8merAt(string mrna, MiRnaAnalyzer.MiRna mirna, int start) =>
        MiRnaAnalyzer.FindTargetSites(mrna, mirna, minScore: 0.0)
            .First(s => s.Start == start && s.Type == MiRnaAnalyzer.TargetSiteType.Seed8mer);

    [Test]
    [Description("MON: the local-AU coefficient is negative, so raising the AU content of the flanks around a fixed seed match makes the local-AU contribution (and the partial context score) monotonically more negative.")]
    public void Context_StrongerLocalAu_GivesMoreNegativeScore()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string site8 = Let7a8merMatch;

        // Flanks of increasing AU content: 0% (GC) → 50% → 100% (AU). 30 nt each side.
        string[] flanks =
        {
            string.Concat(Enumerable.Repeat("GC", 15)), // 0% AU
            string.Concat(Enumerable.Repeat("GCAU", 8))[..30], // ~50% AU
            string.Concat(Enumerable.Repeat("AU", 15)), // 100% AU
        };

        double previousAu = double.PositiveInfinity;
        double previousTotal = double.PositiveInfinity;
        foreach (string flank in flanks)
        {
            string mrna = flank + site8 + flank;
            var site = Find8merAt(mrna, mirna, flank.Length);
            var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, mirna, site);

            ctx.LocalAuContribution.Should().BeLessThan(previousAu + 1e-9,
                because: "raising the flank AU content lowers (makes more negative) the local-AU contribution");
            ctx.ContextScorePartial.Should().BeLessThan(previousTotal + 1e-9,
                because: "a more AU-rich (more accessible) context cannot raise the repression score");
            previousAu = ctx.LocalAuContribution;
            previousTotal = ctx.ContextScorePartial;
        }

        // Non-vacuity: the AU-rich endpoint is strictly more negative than the GC endpoint.
        string gcMrna = flanks[0] + site8 + flanks[0];
        string auMrna = flanks[2] + site8 + flanks[2];
        var gc = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(gcMrna, mirna, Find8merAt(gcMrna, mirna, flanks[0].Length));
        var au = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(auMrna, mirna, Find8merAt(auMrna, mirna, flanks[2].Length));
        au.LocalAuContribution.Should().BeLessThan(gc.LocalAuContribution - 1e-6,
            because: "an AU-rich context yields a strictly more negative local-AU contribution than a GC-rich one");
    }

    [Test]
    [Description("INV: context++ scoring is deterministic — the same (mRNA, miRNA, site) yields the identical contribution breakdown on every call.")]
    public void Context_SameSite_GivesSameScore()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string flank = string.Concat(Enumerable.Repeat("GCAU", 8))[..30];
        string mrna = flank + Let7a8merMatch + flank;
        var site = Find8merAt(mrna, mirna, flank.Length);

        var a = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, mirna, site);
        var b = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, mirna, site);

        a.ContextScorePartial.Should().Be(b.ContextScorePartial, because: "context++ scoring is deterministic");
        a.LocalAuContribution.Should().Be(b.LocalAuContribution);
        a.SaContribution.Should().Be(b.SaContribution);
        a.SiteType.Should().Be(b.SiteType);
        a.OmittedFeatures.Should().Equal(b.OmittedFeatures,
            because: "the entire contribution breakdown is reproduced exactly for the same site");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-PCT-001 — probability of conserved targeting (PCT) (MiRNA)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (Friedman et al. 2009 Genome Res 19:92; targetscan_70_BL_PCT.pl;
    //   tests/TestSpecs/MIRNA-PCT-001.md):
    //   PCT maps a Friedman branch-length score (Bls — total branch length of the subtree of
    //   species in which the site is conserved) to a probability via a logistic curve
    //   PCT(Bls) = B0 + B1/(1 + e^(−B2·Bls + B3)), truncated at 0 (B2 > 0 ⇒ increasing in Bls). When
    //   no conservation evidence is supplied, the context++ scorer reports PCT = 0 (an omitted
    //   feature). Two metamorphic relations (checklist row 253):
    //
    //   • MON (deeper conservation → higher PCT): with B2 > 0 the sigmoid is increasing, so a longer
    //     branch length yields a higher (never lower) PCT.
    //   • INV (no conservation → PCT 0): with no conservation input, the branch-length score and PCT
    //     are 0 and PCT is reported as an omitted feature.
    //
    // API under test: MiRnaAnalyzer.PctFromBranchLength / .ScoreTargetSiteContextPlusPlus.

    #region MIRNA-PCT-001 — probability of conserved targeting

    [Test]
    [Description("MON: with B2 > 0 the PCT logistic is increasing in the branch-length score, so deeper conservation (a longer Bls) yields a higher PCT.")]
    public void Pct_DeeperConservation_RaisesPct()
    {
        // Published-regime parameters (B0 + B1 ≤ 1, B2 > 0): an increasing sigmoid in [0, ~0.85].
        var p = new MiRnaAnalyzer.PctSigmoidParameters(B0: 0.05, B1: 0.80, B2: 1.50, B3: 2.00);

        double previous = double.NegativeInfinity;
        bool sawStrictIncrease = false;
        foreach (double bls in new[] { 0.0, 0.5, 1.0, 2.0, 4.0, 8.0 })
        {
            double pct = MiRnaAnalyzer.PctFromBranchLength(bls, p);
            pct.Should().BeInRange(-1e-9, 1.0 + 1e-9, because: "PCT is a probability in [0,1]");
            pct.Should().BeGreaterThanOrEqualTo(previous - 1e-12,
                because: $"a longer branch length ({bls}) cannot lower the PCT (the sigmoid is increasing)");
            if (pct > previous + 1e-6) sawStrictIncrease = true;
            previous = pct;
        }

        sawStrictIncrease.Should().BeTrue(because: "deeper conservation genuinely raises PCT — the relation is non-vacuous");
    }

    [Test]
    [Description("INV: with no conservation evidence supplied, the context++ scorer reports a zero branch-length score and zero PCT, and lists PCT as an omitted feature.")]
    public void Pct_NoConservation_IsZero()
    {
        var mirna = MiRnaAnalyzer.CreateMiRna("let-7a", Let7a);
        string flank = string.Concat(Enumerable.Repeat("GCAU", 8))[..30];
        string mrna = flank + Let7a8merMatch + flank;
        var site = Find8merAt(mrna, mirna, flank.Length);

        // No ContextPlusPlusInputs.Conservation supplied (default inputs).
        var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, mirna, site);

        ctx.BranchLengthScore.Should().Be(0.0, because: "no conserved species ⇒ branch-length score 0");
        ctx.Pct.Should().Be(0.0, because: "no conservation evidence ⇒ PCT 0");
        ctx.PctContribution.Should().Be(0.0, because: "an omitted PCT feature contributes 0");
        ctx.OmittedFeatures.Any(f => f.StartsWith("PCT")).Should().BeTrue(
            because: "PCT is reported as an omitted feature when no conservation input is supplied");
    }

    #endregion
}
