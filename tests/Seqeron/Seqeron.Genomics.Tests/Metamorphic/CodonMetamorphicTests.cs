using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Codon area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CODON-OPT-001 — codon optimization (Codon).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 58.
///
/// API under test (CodonOptimizer.OptimizeSequence):
///   Replaces each codon with a synonymous codon chosen per the target organism's usage table
///   and strategy (MaximizeCAI = the argmax-frequency synonymous codon; AvoidRareCodeons =
///   replace only codons below the rare threshold). Stop codons are preserved.
///
/// Relations (derived from synonymous substitution, NOT from output):
///   • INV (optimized ⇒ same protein): every substitution is synonymous, so translating the
///          optimized sequence reproduces the original protein.
///   • INV (already optimal ⇒ no change): MaximizeCAI is idempotent — re-optimizing an
///          already-optimized sequence changes nothing, since every codon is already the
///          argmax for its amino acid.
///   • MON (more biased table ⇒ more codon changes): under AvoidRareCodeons, a codon is
///          replaced only when its target-table frequency is below the rare threshold; a table
///          that concentrates mass away from the input's codons (pushing them under the
///          threshold) replaces more of them, so the number of changed codons is non-decreasing
///          in that bias.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class CodonMetamorphicTests
{
    #region Helpers

    // A varied coding sequence (Met-Lys-Leu-Gly-Phe-Glu-Arg-Pro), no internal stop codon.
    private const string CodingSequence = "ATGAAACTAGGTTTTGAACGTCCC";

    /// <summary>Independent translation via the standard genetic code (RNA codons), used as the protein oracle.</summary>
    private static string TranslateRna(string rna)
    {
        var code = CodonOptimizer.EColiK12.CodonToAminoAcid;
        var sb = new StringBuilder();
        for (int i = 0; i + 3 <= rna.Length; i += 3)
            sb.Append(code.GetValueOrDefault(rna.Substring(i, 3), "X"));
        return sb.ToString();
    }

    #endregion

    #region INV — optimization preserves the protein

    [Test]
    [Description("INV: every codon substitution is synonymous, so the optimized sequence translates to the same protein as the original.")]
    public void OptimizeSequence_PreservesProtein()
    {
        var result = CodonOptimizer.OptimizeSequence(
            CodingSequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        string originalProtein = TranslateRna(result.OriginalSequence);
        string optimizedProtein = TranslateRna(result.OptimizedSequence);

        optimizedProtein.Should().Be(originalProtein,
            because: "synonymous codon substitution cannot change the encoded amino acids");
        result.ProteinSequence.Should().Be(originalProtein,
            because: "the reported protein is the translation of the (unchanged-in-meaning) coding sequence");
    }

    #endregion

    #region INV — re-optimizing an already-optimal sequence changes nothing

    [Test]
    [Description("INV: MaximizeCAI is idempotent — optimizing an already-optimized sequence yields zero further changes and an identical sequence.")]
    public void OptimizeSequence_AlreadyOptimal_NoChange()
    {
        var first = CodonOptimizer.OptimizeSequence(
            CodingSequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        var second = CodonOptimizer.OptimizeSequence(
            first.OptimizedSequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        second.ChangedCodons.Should().Be(0,
            because: "after MaximizeCAI every codon is already the argmax for its amino acid, so a second pass changes nothing");
        second.OptimizedSequence.Should().Be(first.OptimizedSequence,
            because: "an already-optimal sequence is a fixed point of the optimizer");
    }

    #endregion

    #region MON — a table biased away from the input's codons changes more codons

    [Test]
    [Description("MON: under AvoidRareCodeons, as the target table concentrates mass away from the input's codon (pushing it under the rare threshold), the number of changed codons is non-decreasing.")]
    public void OptimizeSequence_MoreBiasedTable_ChangesMoreCodons()
    {
        // Input: six Leucine codons, all "CTA" (RNA CUA).
        string input = string.Concat(Enumerable.Repeat("CTA", 6));

        // A family of tables increasingly biased AWAY from CUA (its frequency falls past the
        // 0.15 rare threshold), always offering CUG as a common replacement.
        int previousChanges = -1;
        foreach (double cuaFreq in new[] { 0.30, 0.20, 0.10, 0.02 })
        {
            var table = new CodonOptimizer.CodonUsageTable(
                OrganismName: "test",
                CodonFrequencies: new Dictionary<string, double> { ["CUA"] = cuaFreq, ["CUG"] = 0.80 },
                CodonToAminoAcid: CodonOptimizer.EColiK12.CodonToAminoAcid);

            int changes = CodonOptimizer.OptimizeSequence(
                input, table, CodonOptimizer.OptimizationStrategy.AvoidRareCodeons).ChangedCodons;

            changes.Should().BeGreaterThanOrEqualTo(previousChanges,
                because: "as the table grows more biased against the input codon, at least as many codons fall under the rare threshold and get replaced");
            previousChanges = changes;
        }

        // The most-biased table (CUA rare) replaces every codon; the least-biased (CUA common) replaces none.
        var common = new CodonOptimizer.CodonUsageTable("common",
            new Dictionary<string, double> { ["CUA"] = 0.30, ["CUG"] = 0.80 }, CodonOptimizer.EColiK12.CodonToAminoAcid);
        var biased = new CodonOptimizer.CodonUsageTable("biased",
            new Dictionary<string, double> { ["CUA"] = 0.02, ["CUG"] = 0.80 }, CodonOptimizer.EColiK12.CodonToAminoAcid);

        CodonOptimizer.OptimizeSequence(input, common, CodonOptimizer.OptimizationStrategy.AvoidRareCodeons)
            .ChangedCodons.Should().Be(0, because: "CUA is above the rare threshold, so it is not replaced");
        CodonOptimizer.OptimizeSequence(input, biased, CodonOptimizer.OptimizationStrategy.AvoidRareCodeons)
            .ChangedCodons.Should().Be(6, because: "CUA is rare under the biased table, so all six Leucine codons are replaced by CUG");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CODON-CAI-001 — Codon Adaptation Index (Codon).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 59.
    //
    // API under test (CodonOptimizer.CalculateCAI):
    //   CAI (Sharp & Li 1987) is the geometric mean of the per-codon relative
    //   adaptiveness over the L synonymous (non-stop) codons of a CDS:
    //       CAI = (∏ w_i)^(1/L) = exp( (1/L) Σ ln w_i ),   w_i = f(codon_i) / max_{syn} f
    //   where f is the target organism's codon-usage fraction within the amino-acid
    //   family and max_{syn} f is the fraction of the most frequent synonym. Hence
    //   0 < w_i ≤ 1, every w_i for an organism-optimal codon (and every single-codon
    //   amino acid) equals 1, and stop codons are excluded.
    //
    // Relations (derived from the geometric-mean definition, NOT from output):
    //   • INV  (same seq ⇒ same CAI): CAI is a deterministic pure function of the
    //          (sequence, table) pair — repeated evaluation is identical. Stronger
    //          form: the geometric mean depends only on the MULTISET of per-codon w_i,
    //          so any codon-order permutation (e.g. reversing the codon order) leaves
    //          CAI unchanged even though the sequence — and the encoded protein —
    //          changes.
    //   • COMP (all-optimal codons ⇒ CAI = 1): if every codon is the argmax synonym
    //          for its amino acid then every w_i = 1 and CAI = exp(0) = 1. MaximizeCAI
    //          optimization produces exactly such a sequence, so the optimizer composed
    //          with CalculateCAI must yield 1 for any input.
    //   • MON  (replace rare with optimal ⇒ CAI increases): replacing a sub-optimal
    //          codon by its optimal synonym raises that factor from w<1 to w=1 while
    //          preserving the amino acid; the geometric mean is strictly increasing in
    //          each factor, so the CAI of the (protein-identical) sequence strictly rises.
    // ───────────────────────────────────────────────────────────────────────────

    #region CODON-CAI-001 — Helpers

    /// <summary>Splits an RNA/DNA coding string into its complete 3-letter codons.</summary>
    private static List<string> Codons(string seq)
    {
        var codons = new List<string>();
        string rna = seq.ToUpperInvariant().Replace('T', 'U');
        for (int i = 0; i + 3 <= rna.Length; i += 3)
            codons.Add(rna.Substring(i, 3));
        return codons;
    }

    #endregion

    #region CODON-CAI-001 INV — CAI is deterministic and codon-order invariant

    [Test]
    [Description("INV: CalculateCAI is a pure function — re-evaluating the same sequence against the same table yields the identical CAI.")]
    public void CalculateCAI_SameSequence_SameValue()
    {
        double first = CodonOptimizer.CalculateCAI(CodingSequence, CodonOptimizer.EColiK12);
        double second = CodonOptimizer.CalculateCAI(CodingSequence, CodonOptimizer.EColiK12);

        second.Should().Be(first,
            because: "CAI is a deterministic function of (sequence, codon-usage table) with no hidden state");
    }

    [Test]
    [Description("INV: CAI is the geometric mean over the codon multiset, so reordering the codons (here, reversing them) leaves CAI unchanged though the sequence and its protein differ.")]
    public void CalculateCAI_CodonOrderPermutation_PreservesCAI()
    {
        // A varied coding sequence spanning several amino-acid families and adaptiveness levels.
        const string seq = "AUGCUACGUCCCGAAAGCACCGGG";

        var reversedCodons = Codons(seq);
        reversedCodons.Reverse();
        string permuted = string.Concat(reversedCodons);

        // Sanity: the transformation genuinely changed the sequence (so the test is not vacuous).
        permuted.Should().NotBe(seq.ToUpperInvariant().Replace('T', 'U'),
            because: "reversing the codon order must actually rearrange the sequence for the invariance to be meaningful");

        double original = CodonOptimizer.CalculateCAI(seq, CodonOptimizer.EColiK12);
        double reordered = CodonOptimizer.CalculateCAI(permuted, CodonOptimizer.EColiK12);

        reordered.Should().BeApproximately(original, 1e-9,
            because: "the geometric mean ∏ w_i depends only on the multiset of per-codon adaptiveness values, not on their order");
    }

    #endregion

    #region CODON-CAI-001 COMP — an all-optimal sequence has CAI = 1

    [Test]
    [Description("COMP: MaximizeCAI selects the argmax synonym for every amino acid (w_i = 1 each), so the optimized sequence has CAI = 1 — verified both via the reported OptimizedCAI and an independent CalculateCAI call.")]
    public void CalculateCAI_AllOptimalCodons_EqualsOne()
    {
        // Each input below mixes optimal, sub-optimal and rare codons; after MaximizeCAI
        // every synonymous family collapses to its most-frequent codon, forcing CAI = 1.
        foreach (string input in new[]
                 {
                     CodingSequence,                       // varied 8-codon CDS
                     "CUACUACUACUA",                       // all-rare Leucine
                     "AGAAGGCGACGGGGAGGG",                 // rare Arg + Gly mix
                     "AUGAAACUAGGUUUUGAACGUCCCUAA",        // CDS with a trailing stop codon
                 })
        {
            var result = CodonOptimizer.OptimizeSequence(
                input, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

            result.OptimizedCAI.Should().BeApproximately(1.0, 1e-9,
                because: "after MaximizeCAI every codon is the argmax synonym, so every w_i = 1 and CAI = exp(0) = 1");

            CodonOptimizer.CalculateCAI(result.OptimizedSequence, CodonOptimizer.EColiK12)
                .Should().BeApproximately(1.0, 1e-9,
                    because: "an independent CAI recomputation of the optimized sequence must also be 1");
        }
    }

    #endregion

    #region CODON-CAI-001 MON — replacing rare codons with optimal synonyms raises CAI

    [Test]
    [Description("MON: replacing a rare codon by its optimal synonym (a synonymous, protein-preserving substitution) strictly increases CAI; doing so repeatedly is strictly monotone up to CAI = 1.")]
    public void CalculateCAI_ReplaceRareWithOptimal_IncreasesCAI()
    {
        // Six Leucine codons. CUA is the rarest Leu codon in E. coli (w = 0.04/0.50 = 0.08);
        // CUG is the optimal Leu codon (w = 0.50/0.50 = 1). Each CUA→CUG swap keeps the
        // protein (poly-Leucine) identical but raises one adaptiveness factor from 0.08 to 1.
        const int n = 6;
        double previous = -1.0;

        for (int optimal = 0; optimal <= n; optimal++)
        {
            // `optimal` codons are the best synonym (CUG); the rest remain the rare one (CUA).
            string seq = string.Concat(Enumerable.Repeat("CUG", optimal))
                       + string.Concat(Enumerable.Repeat("CUA", n - optimal));

            // Protein invariance: every codon still encodes Leucine.
            TranslateRna(seq).Should().Be(new string('L', n),
                because: "CUA→CUG is a synonymous substitution, so the encoded protein never changes");

            double cai = CodonOptimizer.CalculateCAI(seq, CodonOptimizer.EColiK12);

            cai.Should().BeGreaterThan(previous,
                because: "each rare→optimal substitution raises one factor (0.08 → 1), strictly increasing the geometric mean");
            previous = cai;
        }

        // Endpoints follow from the definition: all-rare = the shared w = 0.08; all-optimal = 1.
        CodonOptimizer.CalculateCAI(string.Concat(Enumerable.Repeat("CUA", n)), CodonOptimizer.EColiK12)
            .Should().BeApproximately(0.08, 1e-6, because: "every codon shares w = 0.04/0.50, so the geometric mean is 0.08");
        CodonOptimizer.CalculateCAI(string.Concat(Enumerable.Repeat("CUG", n)), CodonOptimizer.EColiK12)
            .Should().BeApproximately(1.0, 1e-9, because: "every codon is optimal (w = 1), so CAI = 1");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CODON-RARE-001 — rare-codon detection (Codon).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 60.
    //
    // API under test (CodonOptimizer.FindRareCodons):
    //   Scans a CDS codon-by-codon and reports every codon whose target-organism usage
    //   fraction is strictly below a threshold τ:  flagged(τ) = { i : f(codon_i) < τ }.
    //   Each hit carries (Position = 3·i, Codon, AminoAcid, Frequency).
    //
    // Relations (derived from the f < τ predicate, NOT from output):
    //   • MON (raising τ ⇒ superset): the predicate f < τ is monotone in τ, so
    //          τ₁ ≤ τ₂ ⇒ flagged(τ₁) ⊆ flagged(τ₂); the flagged set never shrinks as
    //          the threshold rises, and grows exactly by the codons whose frequency lies
    //          in [τ₁, τ₂). (NOTE: this is the opposite direction to the loosely-worded
    //          checklist hint — the implementation flags codons BELOW the threshold, so
    //          a HIGHER threshold flags more.)
    //   • INV (no rare codons ⇒ empty): if no codon satisfies f < τ — e.g. τ = 0 (no
    //          frequency is negative) or a sequence built only from codons whose usage is
    //          ≥ τ — the result is empty.
    // ───────────────────────────────────────────────────────────────────────────

    #region CODON-RARE-001 — Helpers

    /// <summary>The set of codon indices (Position / 3) flagged as rare at a given threshold.</summary>
    private static HashSet<int> RarePositions(string seq, double threshold) =>
        CodonOptimizer.FindRareCodons(seq, CodonOptimizer.EColiK12, threshold)
            .Select(r => r.Position / 3)
            .ToHashSet();

    #endregion

    #region CODON-RARE-001 MON — raising the threshold yields a superset of flagged codons

    [Test]
    [Description("MON: f < τ is monotone in τ, so a higher threshold flags a superset of the codons flagged by a lower threshold, growing exactly by codons whose frequency falls in the opened band.")]
    public void FindRareCodons_HigherThreshold_YieldsSuperset()
    {
        // Codons chosen to span the E. coli usage range so each threshold step opens a new band:
        //   CUA Leu 0.04 | ACU Thr 0.16 | AUC Ile 0.42 | CUG Leu 0.50
        const string seq = "CUAACUAUCCUG";

        double[] thresholds = { 0.0, 0.05, 0.20, 0.45, 0.51, 1.0 };

        HashSet<int>? previous = null;
        foreach (double tau in thresholds)
        {
            var flagged = RarePositions(seq, tau);

            if (previous is not null)
                flagged.IsSupersetOf(previous).Should().BeTrue(
                    because: $"raising the threshold to {tau} can only add codons (f < τ is monotone in τ), never remove them");

            previous = flagged;
        }

        // Endpoints and the intermediate bands follow directly from the f < τ predicate.
        RarePositions(seq, 0.0).Should().BeEmpty(because: "no usage fraction is < 0");
        RarePositions(seq, 0.05).Should().BeEquivalentTo(new[] { 0 }, because: "only CUA (0.04) is below 0.05");
        RarePositions(seq, 0.20).Should().BeEquivalentTo(new[] { 0, 1 }, because: "CUA (0.04) and ACU (0.16) are below 0.20");
        RarePositions(seq, 0.45).Should().BeEquivalentTo(new[] { 0, 1, 2 }, because: "CUA, ACU and AUC (0.42) are below 0.45");
        RarePositions(seq, 0.51).Should().BeEquivalentTo(new[] { 0, 1, 2, 3 }, because: "all four codons (max 0.50) are below 0.51");
    }

    #endregion

    #region CODON-RARE-001 INV — no codon below threshold ⇒ empty result

    [Test]
    [Description("INV: when no codon satisfies f < τ the result is empty — both for τ = 0 (no frequency is negative) and for a sequence built only from codons whose usage is at or above τ.")]
    public void FindRareCodons_NoCodonBelowThreshold_ReturnsEmpty()
    {
        // (a) τ = 0: nothing can be strictly below zero, regardless of how rare the codons are.
        CodonOptimizer.FindRareCodons("CUACUACUA", CodonOptimizer.EColiK12, threshold: 0.0)
            .Should().BeEmpty(because: "no usage fraction is < 0, so even an all-rare sequence flags nothing at τ = 0");

        // (b) All-common sequence under the default threshold: CUG (0.50), CCG (0.53), ACC (0.44)
        //     are each well above 0.15, so none is rare.
        CodonOptimizer.FindRareCodons("CUGCCGACC", CodonOptimizer.EColiK12, threshold: 0.15)
            .Should().BeEmpty(because: "every codon's usage exceeds the threshold, so no codon is flagged as rare");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CODON-USAGE-001 — per-amino-acid codon usage table (Codon).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 61.
    //
    // API under test (CodonOptimizer.CreateCodonTableFromSequence):
    //   Builds an organism-style usage table from a reference CDS: it counts each codon,
    //   groups the observed codons by amino acid, and stores the WITHIN-FAMILY relative
    //   frequency  f(codon) = count(codon) / Σ_{syn present} count  for every codon that
    //   occurs. Thus the per-amino-acid frequencies form a probability distribution.
    //
    // Relations (derived from the count→fraction normalisation, NOT from output):
    //   • INV (duplicate ⇒ same ratios): every codon count scales by the same factor when
    //          the sequence is concatenated with itself (or its codons reordered), and the
    //          normalisation count/Σcount cancels that factor — so the usage ratios are
    //          invariant to sequence duplication and to codon-order permutation.
    //   • COMP (sum per AA = 1): within each amino-acid family the stored fractions are
    //          count/(family total), so they sum to exactly 1 over the codons present.
    // ───────────────────────────────────────────────────────────────────────────

    #region CODON-USAGE-001 — Helpers

    // A CDS where two amino-acid families each appear with a 2:1 codon split, so the
    // usage ratios are non-trivial: Leu CTG×2/CTA×1, Ala GCC×2/GCA×1, plus Lys AAA×1.
    private const string UsageReference = "CTGCTGCTAGCCGCCGCAAAA";

    /// <summary>Sums the table's stored frequencies within each amino-acid family that has any codon present.</summary>
    private static Dictionary<string, double> FamilySums(CodonOptimizer.CodonUsageTable table)
    {
        var byAa = new Dictionary<string, double>();
        foreach (var (codon, freq) in table.CodonFrequencies)
        {
            string aa = table.CodonToAminoAcid.GetValueOrDefault(codon, "X");
            byAa[aa] = byAa.GetValueOrDefault(aa, 0) + freq;
        }
        return byAa;
    }

    #endregion

    #region CODON-USAGE-001 INV — usage ratios are invariant to duplication and codon order

    [Test]
    [Description("INV: duplicating the reference sequence scales every codon count by 2; the count/Σcount normalisation cancels the factor, so the per-codon usage ratios are identical.")]
    public void CreateCodonTable_DuplicatedSequence_PreservesRatios()
    {
        var single = CodonOptimizer.CreateCodonTableFromSequence(UsageReference, "single");
        var doubled = CodonOptimizer.CreateCodonTableFromSequence(UsageReference + UsageReference, "doubled");

        doubled.CodonFrequencies.Keys.Should().BeEquivalentTo(single.CodonFrequencies.Keys,
            because: "duplication adds no new codons, only doubles existing counts");

        foreach (var (codon, freq) in single.CodonFrequencies)
            doubled.CodonFrequencies[codon].Should().BeApproximately(freq, 1e-12,
                because: $"the within-family fraction of {codon} is unchanged when all counts double");
    }

    [Test]
    [Description("INV: usage ratios depend only on the codon multiset, so reordering the codons (here reversing them) yields an identical table.")]
    public void CreateCodonTable_CodonOrderPermutation_PreservesRatios()
    {
        var reversedCodons = Codons(UsageReference);
        reversedCodons.Reverse();
        string permuted = string.Concat(reversedCodons);

        permuted.Should().NotBe(UsageReference.ToUpperInvariant().Replace('T', 'U'),
            because: "reversing codon order must actually rearrange the sequence for the invariance to be meaningful");

        var original = CodonOptimizer.CreateCodonTableFromSequence(UsageReference, "orig");
        var reordered = CodonOptimizer.CreateCodonTableFromSequence(permuted, "perm");

        reordered.CodonFrequencies.Keys.Should().BeEquivalentTo(original.CodonFrequencies.Keys);
        foreach (var (codon, freq) in original.CodonFrequencies)
            reordered.CodonFrequencies[codon].Should().BeApproximately(freq, 1e-12,
                because: "the codon multiset — and hence every within-family fraction — is unchanged by reordering");
    }

    #endregion

    #region CODON-USAGE-001 COMP — within-family frequencies sum to 1

    [Test]
    [Description("COMP: the stored frequencies are count/(family total), so within each amino-acid family that has any codon present they sum to exactly 1.")]
    public void CreateCodonTable_PerAminoAcidFrequencies_SumToOne()
    {
        var table = CodonOptimizer.CreateCodonTableFromSequence(UsageReference, "ref");

        var sums = FamilySums(table);
        sums.Should().NotBeEmpty(because: "the reference sequence contains several amino-acid families");

        foreach (var (aa, sum) in sums)
            sum.Should().BeApproximately(1.0, 1e-12,
                because: $"the codon fractions for amino acid '{aa}' partition that family's counts and must sum to 1");

        // Spot-check the engineered 2:1 splits to prove the fractions are the real ratios,
        // not an accidental 1.0 from single-codon families.
        table.CodonFrequencies["CUG"].Should().BeApproximately(2.0 / 3.0, 1e-12, because: "Leucine is CTG×2 vs CTA×1");
        table.CodonFrequencies["CUA"].Should().BeApproximately(1.0 / 3.0, 1e-12, because: "Leucine is CTG×2 vs CTA×1");
        table.CodonFrequencies["GCC"].Should().BeApproximately(2.0 / 3.0, 1e-12, because: "Alanine is GCC×2 vs GCA×1");
        table.CodonFrequencies["GCA"].Should().BeApproximately(1.0 / 3.0, 1e-12, because: "Alanine is GCC×2 vs GCA×1");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CODON-ENC-001 — Effective Number of Codons (Wright's Nc; Codon).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 213.
    //
    // API under test (CodonUsageAnalyzer.CalculateEnc):
    //   Wright (1990) Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆, where each F̂ is the average per-amino-acid
    //   codon homozygosity F̂ = (n·Σpᵢ² − 1)/(n − 1) for its degeneracy class. Nc ∈ [20, 61]:
    //   20 = extreme bias (one codon per amino acid), 61 = no bias (all synonyms equal).
    //
    // Relations (derived from the homozygosity formula, NOT from output):
    //   • MON  (more biased usage ⇒ lower ENC): concentrating usage onto fewer synonymous codons
    //          raises Σpᵢ², hence F̂, hence lowers each class contribution aaCount/F̂ — so Nc falls.
    //   • INV  (codon order independent): Nc is a function of codon COUNTS only, so permuting the
    //          codons of a sequence leaves it unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region CODON-ENC-001 — Helpers

    // One amino-acid family per degeneracy class (2/3/4/6) so all of F̂₂, F̂₃, F̂₄, F̂₆ are estimable
    // and Wright's fixed numerators (9, 1, 5, 3) all apply.
    private static readonly (char Aa, string[] Codons)[] EncFamilies =
    {
        ('F', new[] { "TTT", "TTC" }),                              // 2-fold
        ('I', new[] { "ATT", "ATC", "ATA" }),                       // 3-fold
        ('V', new[] { "GTT", "GTC", "GTA", "GTG" }),                // 4-fold
        ('L', new[] { "CTT", "CTC", "CTA", "CTG", "TTA", "TTG" }),  // 6-fold
    };

    private const int EncUnit = 12; // base per-codon multiplicity (keeps every family count an integer)

    // Builds a coding sequence; countFor(degeneracy, codonIndex) gives the copies of each codon.
    private static string BuildCoding(Func<int, int, int> countFor)
    {
        var sb = new StringBuilder();
        foreach (var (_, codons) in EncFamilies)
            for (int i = 0; i < codons.Length; i++)
                for (int c = 0; c < countFor(codons.Length, i); c++)
                    sb.Append(codons[i]);
        return sb.ToString();
    }

    // Uniform usage: every synonymous codon equally frequent (least bias).
    private static string UniformCoding() => BuildCoding((_, _) => EncUnit);
    // Moderate 2:1 preference for the first codon (intermediate bias).
    private static string ModerateCoding() => BuildCoding((_, i) => i == 0 ? 2 * EncUnit : EncUnit);
    // Single-codon usage: all mass on the first codon (extreme bias) — same n per family as uniform.
    private static string SingleCodonCoding() => BuildCoding((deg, i) => i == 0 ? deg * EncUnit : 0);

    #endregion

    #region CODON-ENC-001 MON — more biased usage lowers ENC

    [Test]
    [Description("MON: concentrating codon usage raises homozygosity and lowers Wright's Nc — uniform > moderate-bias > single-codon, bounded by the [20, 61] limits.")]
    public void Enc_MoreBiasedUsage_LowerEnc()
    {
        double uniform = CodonUsageAnalyzer.CalculateEnc(UniformCoding());
        double moderate = CodonUsageAnalyzer.CalculateEnc(ModerateCoding());
        double single = CodonUsageAnalyzer.CalculateEnc(SingleCodonCoding());

        uniform.Should().BeGreaterThan(moderate, because: "moving from equal usage to a 2:1 preference increases bias, lowering Nc");
        moderate.Should().BeGreaterThan(single, because: "collapsing each amino acid onto one codon is maximally biased, lowering Nc further");

        uniform.Should().BeApproximately(61.0, 1e-9, because: "unbiased usage reaches the upper limit of 61 sense codons");
        single.Should().BeApproximately(20.0, 1e-9, because: "one codon per amino acid is the extreme-bias limit of 20");
    }

    #endregion

    #region CODON-ENC-001 INV — ENC is independent of codon order

    [Test]
    [Description("INV: Nc is a function of codon counts only, so permuting the codons of a sequence leaves it unchanged.")]
    public void Enc_CodonOrder_Invariant()
    {
        string coding = ModerateCoding();
        double original = CodonUsageAnalyzer.CalculateEnc(coding);

        var codons = Codons(coding); // reuse the CODON-CAI-001 splitter (returns RNA codons)
        var rng = new Random(20260620);
        for (int i = codons.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (codons[i], codons[j]) = (codons[j], codons[i]);
        }
        // Codons() returns RNA (U); map back to DNA (T) for CalculateEnc, which counts DNA codons.
        string shuffled = string.Concat(codons).Replace('U', 'T');

        CodonUsageAnalyzer.CalculateEnc(shuffled).Should().BeApproximately(original, 1e-9,
            because: "Nc depends only on the multiset of codons, not their order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CODON-RSCU-001 — Relative Synonymous Codon Usage (Codon).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 214.
    //
    // API under test (CodonUsageAnalyzer.CalculateRscu):
    //   For codon j of an amino acid with k synonymous codons and counts x,
    //   RSCU_j = x_j / ((1/k)·Σx) = k·x_j/Σx (Sharp, Tuohy & Mosurski 1986). RSCU = 1 means no bias.
    //
    // Relations (derived from the normalisation, NOT from output):
    //   • P    (per-AA RSCU mean = 1): by construction ΣRSCU over a present amino acid's k synonymous
    //          codons equals k, so the per-amino-acid mean is exactly 1 — regardless of how biased the
    //          individual codons are.
    //   • INV  (codon order independent): RSCU is a function of codon COUNTS only, so permuting the
    //          codons of a sequence leaves the whole RSCU table unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region CODON-RSCU-001 — Helpers

    // Amino acid for a DNA codon, via the standard-code table (keyed by RNA codons).
    private static char AminoAcidOf(string dnaCodon) =>
        CodonOptimizer.EColiK12.CodonToAminoAcid.GetValueOrDefault(dnaCodon.Replace('T', 'U'), "X")[0];

    #endregion

    #region CODON-RSCU-001 P — every present amino acid has mean RSCU = 1

    [Test]
    [Description("P: ΣRSCU over a present amino acid's synonymous codons equals its degeneracy, so the per-amino-acid mean is exactly 1, however biased the individual codons are.")]
    public void Rscu_PerAminoAcidMean_IsOne()
    {
        // A 2:1-biased coding sequence over the F/I/V/L families (degeneracy 2/3/4/6).
        var rscu = CodonUsageAnalyzer.CalculateRscu(ModerateCoding());

        var present = new HashSet<char>();
        foreach (var group in rscu.GroupBy(kv => AminoAcidOf(kv.Key)))
        {
            var values = group.Select(kv => kv.Value).ToList();
            if (values.Any(v => v > 0)) // amino acid is present in the sequence
            {
                present.Add(group.Key);
                values.Average().Should().BeApproximately(1.0, 1e-9,
                    because: $"the RSCU values of amino acid '{group.Key}' average to 1 by the RSCU normalisation");
            }
        }

        present.Should().Contain(new[] { 'F', 'I', 'V', 'L' }, because: "the fixture encodes those four amino acids");

        // Non-vacuous: the usage is biased, so individual RSCU values genuinely depart from 1.
        rscu.Values.Should().Contain(v => v > 1.0 + 1e-9, because: "the over-used preferred codons have RSCU > 1");
        rscu.Values.Should().Contain(v => v > 0 && v < 1.0 - 1e-9, because: "the under-used synonyms have RSCU < 1");
    }

    #endregion

    #region CODON-RSCU-001 INV — RSCU is independent of codon order

    [Test]
    [Description("INV: RSCU is a function of codon counts only, so permuting the codons of a sequence leaves the whole RSCU table unchanged.")]
    public void Rscu_CodonOrder_Invariant()
    {
        string coding = ModerateCoding();
        var original = CodonUsageAnalyzer.CalculateRscu(coding);

        var codons = Codons(coding);
        var rng = new Random(20260620);
        for (int i = codons.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (codons[i], codons[j]) = (codons[j], codons[i]);
        }
        string shuffled = string.Concat(codons).Replace('U', 'T');

        var permuted = CodonUsageAnalyzer.CalculateRscu(shuffled);
        permuted.Should().BeEquivalentTo(original,
            because: "the RSCU table depends only on the multiset of codons, not their order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CODON-STATS-001 — codon-count statistics (Codon).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 215.
    //
    // API under test (CodonUsageAnalyzer.CountCodons / GetStatistics.TotalCodons):
    //   Splits a frame-aligned coding sequence into consecutive 3-mers and tallies each codon.
    //
    // Relations (derived from the in-frame tally, NOT from output):
    //   • INV  (order independent): the codon-count multiset depends only on which whole codons are
    //          present, so permuting the codons of a sequence leaves the counts unchanged.
    //   • ADD  (counts additive on concatenation): concatenating two frame-aligned sequences keeps
    //          both reading frames intact, so CountCodons(a+b) = CountCodons(a) + CountCodons(b) and
    //          TotalCodons is additive.
    // ───────────────────────────────────────────────────────────────────────────

    #region CODON-STATS-001 INV — codon counts are independent of codon order

    [Test]
    [Description("INV: the codon-count multiset depends only on which whole codons are present, so permuting the codons leaves the counts unchanged.")]
    public void CodonCounts_CodonOrder_Invariant()
    {
        string coding = ModerateCoding();
        var original = CodonUsageAnalyzer.CountCodons(coding);

        var codons = Codons(coding);
        var rng = new Random(20260620);
        for (int i = codons.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (codons[i], codons[j]) = (codons[j], codons[i]);
        }
        string shuffled = string.Concat(codons).Replace('U', 'T');

        CodonUsageAnalyzer.CountCodons(shuffled).Should().BeEquivalentTo(original,
            because: "codon counts depend only on the multiset of codons, not their order");
    }

    #endregion

    #region CODON-STATS-001 ADD — codon counts are additive over frame-aligned concatenation

    [Test]
    [Description("ADD: concatenating two frame-aligned sequences keeps both reading frames intact, so CountCodons(a+b) = CountCodons(a) + CountCodons(b) and TotalCodons is additive.")]
    public void CodonCounts_Concatenation_AreAdditive()
    {
        string a = CodingSequence;         // 24 nt = 8 codons (multiple of 3)
        string b = ModerateCoding();       // built from whole codons (multiple of 3)
        (a.Length % 3).Should().Be(0);
        (b.Length % 3).Should().Be(0);

        var countA = CodonUsageAnalyzer.CountCodons(a);
        var countB = CodonUsageAnalyzer.CountCodons(b);
        var countAB = CodonUsageAnalyzer.CountCodons(a + b);

        foreach (var codon in countA.Keys.Union(countB.Keys))
            countAB.GetValueOrDefault(codon, 0)
                .Should().Be(countA.GetValueOrDefault(codon, 0) + countB.GetValueOrDefault(codon, 0),
                    because: $"codon {codon}'s count in the concatenation is the sum of its counts in the parts");

        countAB.Keys.Should().BeEquivalentTo(countA.Keys.Union(countB.Keys),
            because: "a frame-aligned concatenation introduces no codons absent from both parts");

        CodonUsageAnalyzer.GetStatistics(a + b).TotalCodons
            .Should().Be(CodonUsageAnalyzer.GetStatistics(a).TotalCodons + CodonUsageAnalyzer.GetStatistics(b).TotalCodons,
                because: "the total codon count is additive over a frame-aligned concatenation");
    }

    #endregion
}
