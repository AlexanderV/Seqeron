using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Translation area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: TRANS-CODON-001 — genetic-code codon table (Translation).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 62.
///
/// API under test (GeneticCode.Standard / GeneticCode.Translate):
///   A genetic code is a TOTAL function C : Σ³ → (20 amino acids ∪ {stop}) over the RNA
///   alphabet Σ = {A,C,G,U}. The standard code (NCBI table 1) assigns each of the 4³ = 64
///   codons a residue, of which exactly three are stop codons. Translate normalises its
///   input (upper-case, T→U) before the lookup, so DNA and RNA spellings are equivalent.
///
/// Relations (derived from the code being a fixed total map on 64 codons, NOT from output):
///   • COMP (covers all 64 codons): every codon in Σ³ is assigned a residue — the table has
///          exactly 64 entries, Translate succeeds for all 64, exactly 3 are stop codons and
///          61 are sense codons, and the 20 canonical amino acids are all represented.
///   • INV  (same code table always): the code is a stable singleton — repeated access to
///          GeneticCode.Standard is the same object with identical content, and Translate is
///          a pure, order-independent function (per-codon residue does not depend on call order).
///   • INV  (translation idempotent): Translate's normalisation is idempotent, so every
///          spelling of a codon that normalises to the same RNA triplet (mixed case, T-for-U)
///          yields the identical residue — re-feeding an already-normalised codon changes nothing.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class TranslationMetamorphicTests
{
    #region Helpers

    private static readonly char[] RnaBases = { 'U', 'C', 'A', 'G' };

    /// <summary>Enumerates all 4³ = 64 RNA codons over {U,C,A,G}.</summary>
    private static IEnumerable<string> AllRnaCodons()
    {
        foreach (char a in RnaBases)
            foreach (char b in RnaBases)
                foreach (char c in RnaBases)
                    yield return $"{a}{b}{c}";
    }

    // The 20 standard proteinogenic amino acids (single-letter codes).
    private static readonly HashSet<char> StandardAminoAcids = new("ACDEFGHIKLMNPQRSTVWY");

    #endregion

    #region COMP — the standard code assigns a residue to all 64 codons

    [Test]
    [Description("COMP: the standard genetic code is a total map on the 64 RNA codons — Translate succeeds for every codon in {A,C,G,U}³ and the table has exactly 64 entries.")]
    public void StandardCode_CoversAll64Codons()
    {
        var codons = AllRnaCodons().ToList();
        codons.Should().HaveCount(64, because: "there are 4³ = 64 codons over the four RNA bases");
        codons.Distinct().Should().HaveCount(64, because: "the enumeration must not repeat a codon");

        foreach (string codon in codons)
        {
            char residue = GeneticCode.Standard.Translate(codon);
            (StandardAminoAcids.Contains(residue) || residue == '*').Should().BeTrue(
                because: $"codon {codon} must map to a standard amino acid or a stop, never an undefined symbol");
        }

        GeneticCode.Standard.CodonTable.Should().HaveCount(64,
            because: "a complete genetic code assigns exactly one residue to each of the 64 codons");
        GeneticCode.Standard.CodonTable.Keys.Should().BeEquivalentTo(codons,
            because: "the table's domain is exactly the set of all 64 RNA codons");
    }

    [Test]
    [Description("COMP: of the 64 codons exactly 3 are stop codons and 61 are sense codons, and all 20 canonical amino acids are encoded (the standard code's redundancy structure).")]
    public void StandardCode_HasThreeStopsAnd61SenseCodonsCoveringAll20AminoAcids()
    {
        var residues = AllRnaCodons().Select(c => GeneticCode.Standard.Translate(c)).ToList();

        residues.Count(r => r == '*').Should().Be(3, because: "the standard code has three stop codons (UAA, UAG, UGA)");
        residues.Count(r => r != '*').Should().Be(61, because: "the remaining 61 codons are sense codons");

        residues.Where(r => r != '*').Distinct().Should().BeEquivalentTo(StandardAminoAcids,
            because: "the 61 sense codons encode all 20 proteinogenic amino acids (a degenerate, surjective code)");
    }

    #endregion

    #region INV — the code is a fixed singleton and translation is order-independent

    [Test]
    [Description("INV: GeneticCode.Standard is a stable singleton — repeated access returns the same object with identical content.")]
    public void StandardCode_RepeatedAccess_IsSameStableTable()
    {
        var first = GeneticCode.Standard;
        var second = GeneticCode.Standard;

        ReferenceEquals(first, second).Should().BeTrue(because: "Standard is a cached singleton, not rebuilt per access");
        first.CodonTable.Should().BeEquivalentTo(second.CodonTable, because: "the codon→residue mapping never changes");
        first.StopCodons.Should().BeEquivalentTo(new[] { "UAA", "UAG", "UGA" },
            because: "the standard stop set is fixed regardless of when it is read");
    }

    [Test]
    [Description("INV: Translate is a pure function — translating a set of codons gives the same per-codon residue irrespective of the order in which they are translated.")]
    public void Translate_OrderIndependent_PerCodonResidueStable()
    {
        var codons = AllRnaCodons().ToList();

        var forward = codons.ToDictionary(c => c, c => GeneticCode.Standard.Translate(c));

        var shuffled = codons.AsEnumerable().Reverse().ToList();
        var backward = shuffled.ToDictionary(c => c, c => GeneticCode.Standard.Translate(c));

        backward.Should().BeEquivalentTo(forward,
            because: "translation has no hidden state — each codon's residue is independent of evaluation order");
    }

    #endregion

    #region INV — normalisation is idempotent: equivalent spellings translate identically

    [Test]
    [Description("INV: Translate normalises (upper-case, T→U) before lookup, and that normalisation is idempotent, so DNA spelling, lower-case and the canonical RNA triplet all yield the same residue.")]
    public void Translate_EquivalentSpellings_YieldIdenticalResidue()
    {
        foreach (string rna in AllRnaCodons())
        {
            char canonical = GeneticCode.Standard.Translate(rna);

            string dna = rna.Replace('U', 'T');                 // DNA spelling
            string lower = rna.ToLowerInvariant();              // lower-case RNA
            string lowerDna = dna.ToLowerInvariant();           // lower-case DNA
            string mixed = rna[0] + dna.Substring(1).ToLowerInvariant(); // mixed case + T

            GeneticCode.Standard.Translate(dna).Should().Be(canonical, because: $"{dna} transcribes to {rna}");
            GeneticCode.Standard.Translate(lower).Should().Be(canonical, because: "case is normalised before lookup");
            GeneticCode.Standard.Translate(lowerDna).Should().Be(canonical, because: "case and T→U are both normalised");
            GeneticCode.Standard.Translate(mixed).Should().Be(canonical, because: "mixed spellings normalise to the same triplet");

            // Idempotence: re-translating the already-canonical RNA form is unchanged.
            GeneticCode.Standard.Translate(rna).Should().Be(canonical, because: "Translate is deterministic on its normalised domain");
        }
    }

    #endregion
}
