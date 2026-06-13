// TRANS-SIXFRAME-001 — Six-Frame Translation and ORF finding
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-SIXFRAME-001.md
// Source: Cock PJA et al. (2009) Biopython, Bioinformatics 25(11):1422-1423 (Bio/SeqUtils six_frame_translations);
//         Rice P et al. (2000) EMBOSS transeq/getorf; NCBI The Genetic Codes (table 1).

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class Translator_SixFrames_Tests
{
    // 39-nt evidence dataset (Evidence §Test Datasets). Expected proteins computed
    // by the Biopython six-frame algorithm under NCBI standard table 1.
    private const string Dna39 = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";

    #region TranslateSixFrames

    // M1 — INV-1: exactly six frames keyed +1,+2,+3,-1,-2,-3 (EMBOSS transeq -frame 6).
    [Test]
    public void TranslateSixFrames_Returns_SixFramesKeyedPlusMinus()
    {
        var dna = new DnaSequence(Dna39);

        var frames = Translator.TranslateSixFrames(dna);

        Assert.Multiple(() =>
        {
            Assert.That(frames, Has.Count.EqualTo(6),
                "A double-stranded sequence has exactly six reading frames (3 forward + 3 reverse).");
            Assert.That(frames.Keys.OrderBy(k => k),
                Is.EqualTo(new[] { -3, -2, -1, 1, 2, 3 }),
                "Frames are keyed +1,+2,+3 (forward) and -1,-2,-3 (reverse complement).");
        });
    }

    // M2 — Forward frame proteins of the 39-nt dataset (NCBI table 1).
    [Test]
    public void TranslateSixFrames_ForwardFrames_MatchEvidenceProteins()
    {
        var dna = new DnaSequence(Dna39);

        var frames = Translator.TranslateSixFrames(dna);

        Assert.Multiple(() =>
        {
            Assert.That(frames[1].Sequence, Is.EqualTo("MAIVMGR*KGAR*"),
                "Frame +1 = translation of the input at offset 0 (Biopython forward loop).");
            Assert.That(frames[2].Sequence, Is.EqualTo("WPL*WAAERVPD"),
                "Frame +2 = translation of the input at offset 1.");
            Assert.That(frames[3].Sequence, Is.EqualTo("GHCNGPLKGCPI"),
                "Frame +3 = translation of the input at offset 2.");
        });
    }

    // M3 — Reverse frame proteins (Biopython convention: revcomp at offset 0/1/2).
    [Test]
    public void TranslateSixFrames_ReverseFrames_MatchEvidenceProteins()
    {
        var dna = new DnaSequence(Dna39);

        var frames = Translator.TranslateSixFrames(dna);

        Assert.Multiple(() =>
        {
            Assert.That(frames[-1].Sequence, Is.EqualTo("LSGTLSAAHYNGH"),
                "Frame -1 = translation of the reverse complement at offset 0 (Biopython reverse loop).");
            Assert.That(frames[-2].Sequence, Is.EqualTo("YRAPFQRPITMA"),
                "Frame -2 = translation of the reverse complement at offset 1.");
            Assert.That(frames[-3].Sequence, Is.EqualTo("IGHPFSGPLQWP"),
                "Frame -3 = translation of the reverse complement at offset 2.");
        });
    }

    // M4 — INV-2: forward frames equal Translate at offsets 0/1/2.
    [Test]
    public void TranslateSixFrames_ForwardFrames_EqualTranslateAtOffsets()
    {
        var dna = new DnaSequence(Dna39);

        var frames = Translator.TranslateSixFrames(dna);

        Assert.Multiple(() =>
        {
            Assert.That(frames[1].Sequence, Is.EqualTo(Translator.Translate(dna, frame: 0).Sequence),
                "Frame +1 must equal direct translation at offset 0.");
            Assert.That(frames[2].Sequence, Is.EqualTo(Translator.Translate(dna, frame: 1).Sequence),
                "Frame +2 must equal direct translation at offset 1.");
            Assert.That(frames[3].Sequence, Is.EqualTo(Translator.Translate(dna, frame: 2).Sequence),
                "Frame +3 must equal direct translation at offset 2.");
        });
    }

    // M5 — INV-3: reverse frames equal translation of reverse complement at offsets 0/1/2.
    [Test]
    public void TranslateSixFrames_ReverseFrames_EqualReverseComplementOffsets()
    {
        var dna = new DnaSequence(Dna39);
        var revComp = dna.ReverseComplement();

        var frames = Translator.TranslateSixFrames(dna);

        Assert.Multiple(() =>
        {
            Assert.That(frames[-1].Sequence, Is.EqualTo(Translator.Translate(revComp, frame: 0).Sequence),
                "Frame -1 must equal translation of the reverse complement at offset 0.");
            Assert.That(frames[-2].Sequence, Is.EqualTo(Translator.Translate(revComp, frame: 1).Sequence),
                "Frame -2 must equal translation of the reverse complement at offset 1.");
            Assert.That(frames[-3].Sequence, Is.EqualTo(Translator.Translate(revComp, frame: 2).Sequence),
                "Frame -3 must equal translation of the reverse complement at offset 2.");
        });
    }

    // M6 — INV-4: trailing partial codon ignored (Biopython fragment_length truncation).
    [Test]
    public void TranslateSixFrames_PartialTrailingCodon_IsIgnored()
    {
        // ATG AAA TAG + trailing "GC" (11 nt). Frame +1 reads 3 full codons; "GC" dropped.
        var dna = new DnaSequence("ATGAAATAGGC");

        var frames = Translator.TranslateSixFrames(dna);

        Assert.That(frames[1].Sequence, Is.EqualTo("MK*"),
            "Frame +1 consumes only complete codons; the trailing 2 nt are ignored.");
    }

    // M7 — null input throws.
    [Test]
    public void TranslateSixFrames_NullInput_ThrowsArgumentNullException()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => Translator.TranslateSixFrames(null!),
            "Null sequence is invalid input.");
    }

    // M8 — empty sequence yields six empty frames.
    [Test]
    public void TranslateSixFrames_EmptySequence_ReturnsSixEmptyFrames()
    {
        var dna = new DnaSequence("");

        var frames = Translator.TranslateSixFrames(dna);

        Assert.Multiple(() =>
        {
            Assert.That(frames, Has.Count.EqualTo(6),
                "An empty sequence still produces all six (empty) frames.");
            Assert.That(frames.Values.All(p => p.Sequence.Length == 0), Is.True,
                "No complete codon exists, so every frame is the empty protein.");
        });
    }

    // C2 — TranslateSixFrames renders internal stop codons as '*' (no early termination).
    [Test]
    public void TranslateSixFrames_InternalStop_IsRenderedNotTerminated()
    {
        // ATG TAA GCT = M * A : the residue after the stop must still appear.
        var dna = new DnaSequence("ATGTAAGCT");

        var frames = Translator.TranslateSixFrames(dna);

        Assert.That(frames[1].Sequence, Is.EqualTo("M*A"),
            "Six-frame translation does not stop at an internal stop codon; it renders '*'.");
    }

    #endregion

    #region FindOrfs

    // M9 — INV-5: forward START->STOP ORF with exact positions and protein (EMBOSS getorf -find 1).
    [Test]
    public void FindOrfs_ForwardStartToStop_ReturnsExactPositionsAndProtein()
    {
        // GGG ATG AAA CCC TAA GGG : ATG at index 3, TAA at indices 12-14.
        var dna = new DnaSequence("GGGATGAAACCCTAAGGG");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "Exactly one START->STOP ORF exists on the forward strand.");
        var orf = orfs[0];
        Assert.Multiple(() =>
        {
            Assert.That(orf.StartPosition, Is.EqualTo(3), "Start = first base of the ATG start codon.");
            Assert.That(orf.EndPosition, Is.EqualTo(14), "End = last base of the TAA stop codon (inclusive).");
            Assert.That(orf.Frame, Is.EqualTo(1), "ORF is in forward frame +1.");
            Assert.That(orf.Protein.Sequence, Is.EqualTo("MKP"),
                "Protein includes the start residue (M) and excludes the stop codon.");
        });
    }

    // M10 — no START codon => no ORF (START->STOP model).
    [Test]
    public void FindOrfs_NoStartCodon_ReturnsEmpty()
    {
        // No ATG/TTG/CTG anywhere.
        var dna = new DnaSequence("AAACCCGGGAAACCCGGG");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Is.Empty, "With no START codon, the START->STOP model emits no ORF.");
    }

    // M11 — ORF shorter than minLength is filtered (getorf -minsize).
    [Test]
    public void FindOrfs_OrfBelowMinLength_IsFiltered()
    {
        // ATG AAA TAA : protein "MK" (2 aa). minLength 3 filters it out.
        var dna = new DnaSequence("ATGAAATAA");

        var orfs = Translator.FindOrfs(dna, minLength: 3, searchBothStrands: false).ToList();

        Assert.That(orfs, Is.Empty, "An ORF whose protein is shorter than minLength is discarded.");
    }

    // M12 — null input throws.
    [Test]
    public void FindOrfs_NullInput_ThrowsArgumentNullException()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => Translator.FindOrfs(null!).ToList(),
            "Null sequence is invalid input.");
    }

    // M13 — INV-6: derived length properties of the M9 ORF.
    [Test]
    public void FindOrfs_OrfResult_LengthDerivations_AreCorrect()
    {
        var dna = new DnaSequence("GGGATGAAACCCTAAGGG");

        var orf = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).Single();

        Assert.Multiple(() =>
        {
            Assert.That(orf.NucleotideLength, Is.EqualTo(12),
                "NucleotideLength = EndPosition - StartPosition + 1 = 14 - 3 + 1.");
            Assert.That(orf.AminoAcidLength, Is.EqualTo(3),
                "AminoAcidLength = protein length = len(\"MKP\").");
        });
    }

    // S1 — both strands: ORF present only on the reverse strand is found with negative frame.
    [Test]
    public void FindOrfs_BothStrands_FindsReverseStrandOrf()
    {
        // Forward strand has no START codon; its reverse complement is GGGATGAAACCCTAAGGG
        // which contains the START->STOP ORF (frame -1).
        var dna = new DnaSequence("CCCTTAGGGTTTCATCCC");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: true).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "Exactly one ORF, located on the reverse strand.");
        var orf = orfs[0];
        Assert.Multiple(() =>
        {
            Assert.That(orf.Frame, Is.EqualTo(-1), "Reverse-strand ORF carries a negative frame label.");
            Assert.That(orf.Protein.Sequence, Is.EqualTo("MKP"),
                "Reverse-strand ORF protein matches the START->STOP region of the reverse complement.");
            Assert.That(orf.StartPosition, Is.EqualTo(3), "Position is in the reverse-complement coordinate frame.");
            Assert.That(orf.EndPosition, Is.EqualTo(14), "Inclusive end of the stop codon in the reverse complement.");
        });
    }

    // S2 — forward-only search must not return the reverse-strand ORF.
    [Test]
    public void FindOrfs_ForwardOnly_DoesNotReturnReverseStrandOrf()
    {
        var dna = new DnaSequence("CCCTTAGGGTTTCATCCC");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Is.Empty,
            "With searchBothStrands=false, the reverse-strand ORF is not searched.");
    }

    // C1 — alternative start codon TTG initiates an ORF and is translated as its residue (L).
    [Test]
    public void FindOrfs_AlternativeStartCodonTtg_InitiatesOrf()
    {
        // GG TTG AAA GGG TAA CC : TTG start at index 2 (frame 3), TAA stop at indices 11-13.
        var dna = new DnaSequence("GGTTGAAAGGGTAACC");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "TTG is a START codon in NCBI standard table 1.");
        var orf = orfs[0];
        Assert.Multiple(() =>
        {
            Assert.That(orf.StartPosition, Is.EqualTo(2), "ORF starts at the TTG start codon.");
            Assert.That(orf.Frame, Is.EqualTo(3), "TTG at index 2 lies in forward frame +3.");
            Assert.That(orf.Protein.Sequence, Is.EqualTo("LKG"),
                "TTG is translated by its actual residue (Leu) at the initiator position in this implementation.");
        });
    }

    #endregion
}
