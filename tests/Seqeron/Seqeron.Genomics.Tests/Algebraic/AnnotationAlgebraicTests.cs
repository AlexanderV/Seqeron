using FsCheck;
using FsCheck.Fluent;

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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-GFF-001 — GFF3 annotation I/O (Annotation)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 31.
    //
    // Model: the lightweight GFF3 helper exports GeneAnnotation records to GFF3
    //        lines (ToGff3) and parses GFF3 lines back to GenomicFeature records
    //        (ParseGff3). The round-trip preserves the feature identity, type,
    //        strand, coordinates and attributes, with the documented coordinate
    //        convention: ToGff3 writes column-4 start as Start+1 (0-based →
    //        1-based GFF3 space), which ParseGff3 reads back unchanged. The score
    //        column is exported as '.', so score is intentionally NOT round-tripped.
    //   — docs/algorithms/Annotation/GFF3_IO.md; GenomeAnnotator.ToGff3 / ParseGff3.
    //
    // Laws under test (checklist row 31):
    //   • RT — ParseGff3(ToGff3(annotations)) reproduces each annotation's id,
    //          type, strand, end, the 1-based start (Start+1) and its attributes.
    //   • ID — empty annotation set → empty GFF (only the ##gff-version directive,
    //          which parses to zero features).
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<string> AlnumArbitrary() =>
        Gen.Elements("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789".ToCharArray())
            .ArrayOf()
            .Where(a => a.Length is > 0 and <= 8)
            .Select(a => new string(a))
            .ToArbitrary();

    private static Arbitrary<GenomeAnnotator.GeneAnnotation> GeneAnnotationArbitrary()
    {
        var alnum = AlnumArbitrary().Generator;
        return (from id in alnum
                from product in alnum
                from type in Gen.Elements("gene", "mRNA", "exon", "CDS", "tRNA")
                from start in Gen.Choose(0, 100_000)
                from len in Gen.Choose(1, 5_000)
                from strand in Gen.Elements('+', '-', '.')
                from akey in alnum
                from aval in alnum
                from hasAttr in Gen.Elements(true, false)
                let attrs = hasAttr && akey is not ("ID" or "product" or "translation")
                    ? (IReadOnlyDictionary<string, string>)new Dictionary<string, string> { [akey] = aval }
                    : new Dictionary<string, string>()
                select new GenomeAnnotator.GeneAnnotation(
                    GeneId: id,
                    Start: start,
                    End: start + len,
                    Strand: strand,
                    Type: type,
                    Product: product,
                    Attributes: attrs))
            .ToArbitrary();
    }

    /// <summary>
    /// RT: ParseGff3(ToGff3(annotations)) reproduces every annotation. Identity,
    /// type, strand and end survive verbatim; the start survives in the documented
    /// 1-based convention (parsed Start = original Start + 1); ID/product/extra
    /// attributes survive in column 9.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Gff_RoundTrip_ParseOfSerializeIsIdentity()
    {
        return Prop.ForAll(GeneAnnotationArbitrary().Generator.NonEmptyListOf().ToArbitrary(), annotations =>
        {
            var lines = GenomeAnnotator.ToGff3(annotations);
            var parsed = GenomeAnnotator.ParseGff3(lines).ToList();

            if (parsed.Count != annotations.Count)
                return false.Label($"count {parsed.Count} != {annotations.Count}");

            for (int i = 0; i < annotations.Count; i++)
            {
                var a = annotations[i];
                var f = parsed[i];
                bool ok = f.FeatureId == a.GeneId
                    && f.Type == a.Type
                    && f.Start == a.Start + 1
                    && f.End == a.End
                    && f.Strand == a.Strand
                    && f.Attributes.GetValueOrDefault("ID") == a.GeneId
                    && f.Attributes.GetValueOrDefault("product") == a.Product
                    && a.Attributes.All(kv => f.Attributes.GetValueOrDefault(kv.Key) == kv.Value);
                if (!ok)
                    return false.Label($"record {i} did not round-trip: a={a}, f={f}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// ID: serializing an empty annotation set yields only the ##gff-version
    /// directive, which parses back to zero features.
    /// </summary>
    [Test]
    public void Gff_Identity_EmptyAnnotationsYieldEmptyGff()
    {
        var lines = GenomeAnnotator.ToGff3(System.Array.Empty<GenomeAnnotator.GeneAnnotation>()).ToList();
        GenomeAnnotator.ParseGff3(lines).Should().BeEmpty();
    }
}
