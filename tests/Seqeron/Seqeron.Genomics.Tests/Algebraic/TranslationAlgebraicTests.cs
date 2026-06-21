using System.Collections.Generic;
using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Translation area (the genetic code, DNA→protein).
///
/// Algebraic testing pins the completeness/surjectivity of the genetic-code
/// mapping (all 64 codons defined, every amino acid reachable) and the
/// length-conservation and fixed-point laws of translation.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 62, 63.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Translation")]
public class TranslationAlgebraicTests
{
    private const string Bases = "ACGU";
    private const string AminoAcids = "ACDEFGHIKLMNPQRSTVWY"; // the 20 standard amino acids

    private static IEnumerable<string> All64Codons()
    {
        foreach (char a in Bases)
            foreach (char b in Bases)
                foreach (char c in Bases)
                    yield return $"{a}{b}{c}";
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-CODON-001 — Genetic code (Translation)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 62.
    //
    // Model: the standard genetic code is a total function from the 64 triplets
    //        over {A,C,G,U} onto the 20 amino acids plus the stop signal; it is
    //        complete (every codon defined) and surjective onto the amino acids
    //        (every amino acid encoded by ≥ 1 codon — the degeneracy of the code).
    //   — docs/algorithms/Translation; GeneticCode.Standard.
    //
    // Laws under test (checklist row 62):
    //   • ID   — the code is complete: all 64 codons map to a residue or stop.
    //   • DIST — each amino acid is mapped from ≥ 1 codon (surjective/degenerate).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>ID: every one of the 64 codons is defined by the standard code.</summary>
    [Test]
    public void GeneticCode_Identity_AllSixtyFourCodonsDefined()
    {
        GeneticCode.Standard.CodonTable.Should().HaveCount(64);
        foreach (var codon in All64Codons())
        {
            GeneticCode.Standard.CodonTable.Should().ContainKey(codon);
            // Translate agrees with the table and never throws.
            GeneticCode.Standard.Translate(codon).Should().Be(GeneticCode.Standard.CodonTable[codon]);
        }
    }

    /// <summary>DIST: each of the 20 amino acids is encoded by at least one codon.</summary>
    [Test]
    public void GeneticCode_Distributive_EveryAminoAcidEncoded()
    {
        foreach (char aa in AminoAcids)
            GeneticCode.Standard.GetCodonsForAminoAcid(aa).Should().NotBeEmpty($"amino acid '{aa}'");
    }

    /// <summary>DIST: the codons partition into the amino acids and stop they encode —
    /// exactly the 64 codons are covered with no codon unmapped.</summary>
    [Test]
    public void GeneticCode_Distributive_CodonsPartitionByProduct()
    {
        var covered = AminoAcids.SelectMany(aa => GeneticCode.Standard.GetCodonsForAminoAcid(aa))
            .Concat(GeneticCode.Standard.StopCodons)
            .Distinct()
            .ToHashSet();
        covered.Should().BeEquivalentTo(All64Codons());
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-PROT-001 — DNA→protein translation (Translation)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 63.
    //
    // Model: translation reads non-overlapping triplets, so a length-n DNA yields
    //        at most ⌊n/3⌋ residues; ATG codes methionine; a stop codon terminates
    //        translation when reading to the first stop.
    //   — docs/algorithms/Translation; Translator.Translate.
    //
    // Laws under test (checklist row 63):
    //   • DIST — |protein| ≤ |dna|/3 (triplet length conservation).
    //   • ID   — ATG → M.
    //   • ID   — a stop codon terminates translation (toFirstStop).
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<string> DnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length > 0).Select(a => new string(a)).ToArbitrary();

    /// <summary>DIST: the protein length never exceeds ⌊|dna|/3⌋, with or without
    /// stop-codon truncation.</summary>
    [FsCheck.NUnit.Property]
    public Property Translate_Distributive_ProteinLengthBounded()
    {
        return Prop.ForAll(DnaArbitrary(), dna =>
        {
            int full = Translator.Translate(dna, GeneticCode.Standard, 0, false).Sequence.Length;
            int toStop = Translator.Translate(dna, GeneticCode.Standard, 0, true).Sequence.Length;
            int bound = dna.Length / 3;
            return (full <= bound && toStop <= full)
                .Label($"|prot|={full}, toStop={toStop}, bound={bound} for len={dna.Length}");
        });
    }

    /// <summary>ID: the start codon ATG translates to methionine (M).</summary>
    [Test]
    public void Translate_Identity_AtgIsMethionine()
    {
        Translator.Translate("ATG", GeneticCode.Standard, 0, false).Sequence.Should().Be("M");
        GeneticCode.Standard.Translate("AUG").Should().Be('M');
    }

    /// <summary>ID: each stop codon is '*', and reading to the first stop terminates
    /// translation right after the preceding residue.</summary>
    [Test]
    public void Translate_Identity_StopCodonTerminates()
    {
        foreach (var stop in new[] { "TAA", "TAG", "TGA" })
            GeneticCode.Standard.Translate(stop).Should().Be('*');

        // ATG (M) AAA (K) TAA (stop) GGG (G) — toFirstStop halts at TAA, giving "MK".
        Translator.Translate("ATGAAATAAGGG", GeneticCode.Standard, 0, true).Sequence.Should().Be("MK");
    }
}
