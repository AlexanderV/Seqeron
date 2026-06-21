using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Annotation area (ORF detection, GFF round-trip).
///
/// Algebraic testing pins the structural equations every annotation must obey:
/// the identity behaviour on inputs with no signal, the reading-frame
/// conservation law (lengths are multiples of the codon size), and the
/// parse∘serialize round-trip isomorphism of the GFF representation.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 28, 31.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Annotation")]
public class AnnotationAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-ORF-001 — Open reading frame detection (Annotation)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 28.
    //
    // Model: an ORF runs in-frame from a start codon to the next in-frame stop
    //        codon; its nucleotide span is therefore always a whole number of
    //        codons. This implementation accepts the bacterial start set
    //        {ATG, GTG, TTG} (GenomeAnnotator.StartCodons), so the precise
    //        identity law is "no START CODON present → no ORFs" — strictly
    //        stronger than the checklist's "no ATG" phrasing, which we honour by
    //        choosing a sequence containing NONE of the three start codons.
    //   — docs/algorithms/Annotation/ORF_Detection.md; GenomeAnnotator.FindOrfs.
    //
    // Laws under test (checklist row 28):
    //   • ID   — no start codon → no ORFs (with requireStartCodon = true).
    //   • DIST — reading-frame conservation: every ORF's nucleotide length
    //            (End − Start and Sequence.Length) is divisible by 3, and the
    //            translated protein has exactly length/3 residues.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: a poly-C sequence contains no start codon (CCC ∉ {ATG,GTG,TTG}) and no
    /// stop codon, so ORF detection returns the empty set.
    /// </summary>
    [Test]
    public void Orf_Identity_NoStartCodonYieldsNoOrfs()
    {
        string polyC = new string('C', 60);
        GenomeAnnotator.FindOrfs(polyC, minLength: 1, searchBothStrands: true, requireStartCodon: true)
            .Should().BeEmpty();
    }

    /// <summary>
    /// DIST: every detected ORF spans a whole number of codons — (End − Start) and
    /// Sequence.Length are multiples of 3, and ProteinSequence.Length = Sequence.Length / 3.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orf_Distributive_LengthDivisibleByThree()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 9), seq =>
        {
            var orfs = GenomeAnnotator
                .FindOrfs(seq, minLength: 1, searchBothStrands: true, requireStartCodon: true)
                .ToList();
            bool ok = orfs.All(o =>
                o.Sequence.Length % 3 == 0
                && (o.End - o.Start) % 3 == 0
                && o.Sequence.Length == o.End - o.Start
                && o.ProteinSequence.Length == o.Sequence.Length / 3);
            return ok.Label($"an ORF length was not a codon multiple in \"{seq}\"");
        });
    }

    /// <summary>
    /// DIST witness: ATG-AAA-TAA is a single 9-nt ORF (start→stop), length divisible
    /// by 3 and translating to a 3-residue product (M, K, stop).
    /// </summary>
    [Test]
    public void Orf_Distributive_WorkedExample()
    {
        var orfs = GenomeAnnotator
            .FindOrfs("ATGAAATAA", minLength: 1, searchBothStrands: false, requireStartCodon: true)
            .ToList();
        orfs.Should().NotBeEmpty();
        orfs.Should().OnlyContain(o => o.Sequence.Length % 3 == 0);
        orfs.Should().Contain(o => o.Sequence == "ATGAAATAA" && o.Sequence.Length == 9);
    }
}
