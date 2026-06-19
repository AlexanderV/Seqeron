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
}
