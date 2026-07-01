// GENOMIC-ORF-001 — Open Reading Frame (ORF) Detection
// Evidence: docs/Evidence/GENOMIC-ORF-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-ORF-001.md
// Source: Rosalind "Open Reading Frames" (https://rosalind.info/problems/orf/, 2026-06-14);
//         NCBI Genetic Codes transl_table=1 (Standard).

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_FindOpenReadingFrames_Tests
{
    private static readonly string[] StopCodons = { "TAA", "TAG", "TGA" };

    // Translate an ORF DNA span to its protein candidate, excluding the terminating stop codon.
    private static string ProteinOf(OrfInfo orf) =>
        Translator.Translate(orf.Sequence, GeneticCode.Standard, frame: 0, toFirstStop: true).Sequence;

    #region FindOpenReadingFrames

    // M1 — Single forward ORF: ATG + AAA + AAA + TAA. Derived from the Standard genetic code:
    // span includes the stop (12 nt), protein candidate = MKK.
    [Test]
    public void FindOpenReadingFrames_SingleForwardOrf_ReturnsExactOrf()
    {
        var dna = new DnaSequence("ATGAAAAAATAA");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "Exactly one ATG→stop ORF exists in the only frame containing ATG.");
        Assert.Multiple(() =>
        {
            Assert.That(orfs[0].Sequence, Is.EqualTo("ATGAAAAAATAA"), "ORF span runs start codon through stop codon inclusive.");
            Assert.That(orfs[0].Position, Is.EqualTo(0), "ATG begins at offset 0.");
            Assert.That(orfs[0].Frame, Is.EqualTo(1), "Frame 1 corresponds to offset 0.");
            Assert.That(orfs[0].IsReverseComplement, Is.False, "ORF is on the forward strand.");
            Assert.That(ProteinOf(orfs[0]), Is.EqualTo("MKK"), "AAA AAA translate to KK after the start M; stop excluded.");
        });
    }

    // M2 — Rosalind sample dataset: six-frame search returns exactly four distinct proteins.
    // Source: Rosalind ORF problem, verbatim expected output.
    [Test]
    public void FindOpenReadingFrames_RosalindSampleDataset_ReturnsFourDistinctProteins()
    {
        var dna = new DnaSequence(
            "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG");
        var expected = new HashSet<string>
        {
            "MLLGSFRLIPKETLIQVAGSSPCNLS",
            "M",
            "MGMTPRLGLESLLE",
            "MTPRLGLESLLE",
        };

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();
        var proteins = orfs.Select(ProteinOf).ToHashSet();

        Assert.That(proteins, Is.EquivalentTo(expected),
            "Six-frame ORF detection must reproduce the exact set of protein candidates from the Rosalind worked example.");
    }

    // M3 — Nested ORFs sharing a stop: ATG GGG ATG CCC TAA. Both ATGs (offset 0 and 6) reach the
    // same TAA and must both be reported (Rosalind canonical semantics).
    [Test]
    public void FindOpenReadingFrames_NestedOrfsSharingStop_BothReported()
    {
        var dna = new DnaSequence("ATGGGGATGCCCTAA");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1)
            .Where(o => !o.IsReverseComplement && o.Frame == 1)
            .OrderBy(o => o.Position)
            .ToList();

        Assert.That(orfs, Has.Count.EqualTo(2), "Both the outer and the inner ATG terminate at the shared stop.");
        Assert.Multiple(() =>
        {
            Assert.That(orfs[0].Position, Is.EqualTo(0), "Outer ORF starts at the first ATG.");
            Assert.That(orfs[0].Sequence, Is.EqualTo("ATGGGGATGCCCTAA"), "Outer ORF spans both ATGs to the stop.");
            Assert.That(orfs[1].Position, Is.EqualTo(6), "Inner ORF starts at the second ATG.");
            Assert.That(orfs[1].Sequence, Is.EqualTo("ATGCCCTAA"), "Inner ORF spans the second ATG to the shared stop.");
        });
    }

    // M4 — ATG with no downstream in-frame stop is not a complete ORF (Rosalind: translate until a stop).
    [Test]
    public void FindOpenReadingFrames_AtgWithoutInFrameStop_ReturnsEmpty()
    {
        // ATG followed by lysine codons and no stop in this frame (or its RC).
        var dna = new DnaSequence("ATGAAAAAAAAA");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1)
            .Where(o => !o.IsReverseComplement)
            .ToList();

        Assert.That(orfs, Is.Empty, "An ATG with no in-frame stop is an incomplete ORF and must not be reported.");
    }

    // M5 — Reverse-complement-only ORF. Forward strand has no ATG; the reverse complement of
    // "TTAGGGGGGCAT" is "ATGCCCCCCTAA", a complete ORF.
    [Test]
    public void FindOpenReadingFrames_ReverseStrandOrf_DetectedOnReverseComplement()
    {
        var dna = new DnaSequence("TTAGGGGGGCAT");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "The single complete ORF lies on the reverse complement.");
        Assert.Multiple(() =>
        {
            Assert.That(orfs[0].IsReverseComplement, Is.True, "ORF is found on the reverse complement strand.");
            Assert.That(orfs[0].Sequence, Is.EqualTo("ATGCCCCCCTAA"), "Reverse-complement ORF span start→stop inclusive.");
            Assert.That(ProteinOf(orfs[0]), Is.EqualTo("MPP"), "CCC CCC translate to PP after the start M.");
        });
    }

    // M6 — minLength excludes a too-short ORF: ATGAAATAA is 9 nt; minLength 12 excludes it.
    [Test]
    public void FindOpenReadingFrames_BelowMinLength_Excluded()
    {
        var dna = new DnaSequence("ATGAAATAA");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 12).ToList();

        Assert.That(orfs, Is.Empty, "A 9 nt ORF is below the 12 nt threshold and must be excluded.");
    }

    // M6b — minLength inclusive: the same 9 nt ORF is kept at minLength 9.
    [Test]
    public void FindOpenReadingFrames_ExactlyMinLength_Included()
    {
        var dna = new DnaSequence("ATGAAATAA");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 9)
            .Where(o => !o.IsReverseComplement && o.Position == 0)
            .ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "minLength is an inclusive lower bound; a length-9 ORF passes minLength 9.");
        Assert.That(orfs[0].Sequence, Is.EqualTo("ATGAAATAA"), "Exact 9 nt ORF span.");
    }

    // M7/M8/M9 — Invariants over the Rosalind dataset: every ORF starts ATG, ends with a stop,
    // and has length divisible by 3.
    [Test]
    public void FindOpenReadingFrames_AllOrfs_SatisfyStructuralInvariants()
    {
        var dna = new DnaSequence(
            "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();

        Assert.That(orfs, Is.Not.Empty, "Sanity: the dataset contains ORFs.");
        Assert.Multiple(() =>
        {
            Assert.That(orfs.All(o => o.Sequence.StartsWith("ATG", StringComparison.Ordinal)), Is.True,
                "INV-01: every ORF begins with the ATG start codon.");
            Assert.That(orfs.All(o => StopCodons.Contains(o.Sequence[^3..])), Is.True,
                "INV-02: every ORF ends with a stop codon (TAA/TAG/TGA).");
            Assert.That(orfs.All(o => o.Length % 3 == 0), Is.True,
                "INV-03: every ORF length is divisible by 3 (whole codons start→stop).");
        });
    }

    // M10 — All three stop codons are recognized.
    [TestCase("ATGTAA")]
    [TestCase("ATGTAG")]
    [TestCase("ATGTGA")]
    public void FindOpenReadingFrames_EachStopCodon_TerminatesOrf(string sequence)
    {
        var dna = new DnaSequence(sequence);

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1)
            .Where(o => !o.IsReverseComplement && o.Position == 0 && o.Frame == 1)
            .ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), $"ATG immediately followed by {sequence[^3..]} is a complete 6 nt ORF.");
        Assert.That(orfs[0].Sequence, Is.EqualTo(sequence), "ORF span is the start codon plus the stop codon.");
    }

    // S1 — Lowercase input is handled (DnaSequence normalizes case): same ORF as M1.
    [Test]
    public void FindOpenReadingFrames_LowercaseInput_SameAsUppercase()
    {
        var lower = new DnaSequence("atgaaaaaataa");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(lower, minLength: 1).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "Lowercase input must be detected identically to uppercase.");
        Assert.That(orfs[0].Sequence, Is.EqualTo("ATGAAAAAATAA"), "Normalized ORF span matches the uppercase case.");
    }

    // S2 — ORFs in different frames are each reported with the correct frame number.
    // Frame 1: ATGAAATAA at offset 0. Frame 2: ATG at offset 1 (G ATGAAATGA...) reaching a stop.
    [Test]
    public void FindOpenReadingFrames_DifferentFrames_ReportedWithCorrectFrame()
    {
        // Offset 0: ATG AAA TGA  -> frame 1 ORF (6+3). Offset 1: TGA... no.
        // Use two clearly separate frames: frame1 ORF then a frame2 ORF.
        // "ATGAAATGA" = frame1: ATG AAA TGA -> ORF ATGAAATGA (protein MK).
        // Prefix one base to shift: "C" + "ATGCCCTAA" places an ORF at offset 1 (frame 2).
        var frame1 = new DnaSequence("ATGAAATGA");
        var frame2 = new DnaSequence("CATGCCCTAA");

        var f1 = GenomicAnalyzer.FindOpenReadingFrames(frame1, minLength: 1)
            .Where(o => !o.IsReverseComplement && o.Position == 0).ToList();
        var f2 = GenomicAnalyzer.FindOpenReadingFrames(frame2, minLength: 1)
            .Where(o => !o.IsReverseComplement && o.Position == 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(f1, Has.Count.EqualTo(1), "Frame-1 ORF at offset 0 is detected.");
            Assert.That(f1[0].Frame, Is.EqualTo(1), "Offset 0 maps to frame 1.");
            Assert.That(f2, Has.Count.EqualTo(1), "Frame-2 ORF at offset 1 is detected.");
            Assert.That(f2[0].Frame, Is.EqualTo(2), "Offset 1 maps to frame 2.");
        });
    }

    // C1 — Sequence too short to contain a stop codon yields no ORF.
    [Test]
    public void FindOpenReadingFrames_TooShortForStop_ReturnsEmpty()
    {
        var dna = new DnaSequence("ATG");

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength: 1).ToList();

        Assert.That(orfs, Is.Empty, "An ATG with no room for an in-frame stop is not a complete ORF.");
    }

    // Edge — null sequence throws.
    [Test]
    public void FindOpenReadingFrames_NullSequence_Throws()
    {
        Assert.That(() => GenomicAnalyzer.FindOpenReadingFrames(null!, minLength: 1).ToList(),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null input is a precondition violation.");
    }

    #endregion
}
