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
}
