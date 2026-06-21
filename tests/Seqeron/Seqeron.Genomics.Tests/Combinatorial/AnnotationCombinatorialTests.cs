namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Annotation area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Annotation")]
public class AnnotationCombinatorialTests
{
    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ANNOT-ORF-001 — Open-reading-frame detection (Annotation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 28.
    // Spec: tests/TestSpecs/ANNOT-ORF-001.md (canonical GenomeAnnotator.FindOrfs).
    // Dimensions: minLen(3) × frame(3: fwd/rev/both) × startCodon(2: ATG/alt). Grid 3×3×2 = 18.
    //
    // Model (Wikipedia ORF; Rosalind): an ORF is a maximal run start→stop in one of the
    // six frames; FindOrfs reports its amino-acid length as the codons from the start up
    // to (excluding) the stop, filters by minLength, recognises the canonical ATG plus the
    // alternative bacterial starts GTG/TTG, and — with requireStartCodon — emits only
    // runs that begin at a start codon AND terminate at a stop. searchBothStrands governs
    // whether the reverse complement is scanned; reverse hits carry IsReverseComplement.
    //
    // The combinatorial point: minLength, strand and start codon interact. A length-9 ORF
    // is detected for minLength ≤ 9 and filtered above it; it is found on whichever strand
    // it sits (with the correct IsReverseComplement flag); and both ATG- and GTG-started
    // ORFs are recognised. The grid embeds one clean ORF (no internal start/stop) in
    // start/stop-free poly-C filler so detection is exact and unambiguous.
    // ═══════════════════════════════════════════════════════════════════════

    public enum Orientation { Forward, Reverse, Both }

    private const int OrfAaLength = 9;                       // ATG/GTG + 8×AAA, stop excluded
    private static string MakeOrf(string start) => start + string.Concat(Enumerable.Repeat("AAA", 8)) + "TAA";
    private static readonly string Filler = new('C', 12);   // no start/stop codon on either strand

    [Test, Combinatorial]
    public void AnnotOrf_Detection_AcrossMinLenStrandAndStartCodon(
        [Values(5, 9, 15)] int minLength,
        [Values(Orientation.Forward, Orientation.Reverse, Orientation.Both)] Orientation orientation,
        [Values("ATG", "GTG")] string startCodon)
    {
        string orf = MakeOrf(startCodon);
        string template = orientation switch
        {
            Orientation.Forward => Filler + orf + Filler,
            Orientation.Reverse => Filler + RevComp(orf) + Filler,
            _ => Filler + orf + Filler + RevComp(orf) + Filler,
        };

        var orfs = GenomeAnnotator
            .FindOrfs(template, minLength, searchBothStrands: true, requireStartCodon: true)
            .ToList();

        bool detected = minLength <= OrfAaLength;
        var forward = orfs.Where(o => o.Sequence == orf && !o.IsReverseComplement).ToList();
        var reverse = orfs.Where(o => o.Sequence == orf && o.IsReverseComplement).ToList();

        int expFwd = orientation != Orientation.Reverse && detected ? 1 : 0;
        int expRev = orientation != Orientation.Forward && detected ? 1 : 0;
        forward.Should().HaveCount(expFwd, $"forward ORF presence under {orientation}, minLength {minLength}");
        reverse.Should().HaveCount(expRev, $"reverse ORF presence under {orientation}, minLength {minLength}");

        // Every reported match obeys the ORF invariants (start/stop/in-frame).
        foreach (var o in forward.Concat(reverse))
        {
            o.Sequence[..3].Should().BeOneOf("ATG", "GTG", "TTG");
            o.Sequence[^3..].Should().BeOneOf("TAA", "TAG", "TGA");
            (o.Sequence.Length % 3).Should().Be(0);
        }
    }

    /// <summary>
    /// Interaction witness: searchBothStrands is the strand gate. A reverse-strand ORF is
    /// invisible to a forward-only search but found once both strands are scanned.
    /// </summary>
    [Test]
    public void AnnotOrf_ForwardOnlySearch_IgnoresReverseStrandOrf()
    {
        string orf = MakeOrf("ATG");
        string template = Filler + RevComp(orf) + Filler; // ORF only on the reverse strand

        GenomeAnnotator.FindOrfs(template, minLength: 5, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == orf).Should().BeFalse("forward-only search cannot see a reverse ORF");
        GenomeAnnotator.FindOrfs(template, minLength: 5, searchBothStrands: true, requireStartCodon: true)
            .Any(o => o.Sequence == orf && o.IsReverseComplement).Should().BeTrue("both-strand search finds it");
    }

    /// <summary>
    /// Interaction witness: minLength is an inclusive amino-acid threshold — an ORF exactly
    /// at the threshold is kept, one amino acid longer than the threshold-minus-one is the
    /// boundary. — spec M06/M06b.
    /// </summary>
    [Test]
    public void AnnotOrf_MinLength_IsInclusiveThreshold()
    {
        string orf = MakeOrf("ATG"); // 9 aa
        string template = Filler + orf + Filler;

        GenomeAnnotator.FindOrfs(template, minLength: OrfAaLength, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == orf).Should().BeTrue("length == minLength is included");
        GenomeAnnotator.FindOrfs(template, minLength: OrfAaLength + 1, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == orf).Should().BeFalse("length < minLength is excluded");
    }

    /// <summary>
    /// Interaction witness: the alternative start TTG is also recognised, and requiring a
    /// start codon suppresses a start-less sense run that ends in a stop.
    /// </summary>
    [Test]
    public void AnnotOrf_StartCodonPolicy_AltStartRecognised_RequireStartFilters()
    {
        string ttgOrf = MakeOrf("TTG");
        string template = Filler + ttgOrf + Filler;
        GenomeAnnotator.FindOrfs(template, minLength: 5, searchBothStrands: false, requireStartCodon: true)
            .Any(o => o.Sequence == ttgOrf).Should().BeTrue("TTG is an alternative start codon");

        // A frame with no start codon: AAA…AAA TAA. requireStartCodon must yield nothing.
        string startless = string.Concat(Enumerable.Repeat("AAA", 9)) + "TAA";
        GenomeAnnotator.FindOrfs(startless, minLength: 1, searchBothStrands: false, requireStartCodon: true)
            .Should().BeEmpty("with requireStartCodon a run lacking a start codon is not an ORF");
    }

    /// <summary>
    /// Worked example: the textbook minimal ORF ATG-AAA-TAA is one frame-1 forward ORF
    /// spanning [0,9) translating to M-K (stop trimmed from the protein).
    /// </summary>
    [Test]
    public void AnnotOrf_MinimalOrf_WorkedExample()
    {
        var orfs = GenomeAnnotator.FindOrfs("ATGAAATAA", minLength: 1, searchBothStrands: false, requireStartCodon: true).ToList();

        orfs.Should().ContainSingle();
        var o = orfs[0];
        o.Start.Should().Be(0);
        o.End.Should().Be(9);
        o.Frame.Should().Be(1);
        o.IsReverseComplement.Should().BeFalse();
        o.Sequence.Should().Be("ATGAAATAA");
        o.ProteinSequence.Should().StartWith("MK");
    }
}
