// TRANS-SIXFRAME-001 — Six-Frame Translation and ORF finding
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-SIXFRAME-001.md
// Source: Cock PJA et al. (2009) Biopython, Bioinformatics 25(11):1422-1423 (Bio/SeqUtils six_frame_translations);
//         Rice P et al. (2000) EMBOSS transeq/getorf; NCBI The Genetic Codes (table 1).

namespace Seqeron.Genomics.Tests.Unit.Core;

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
    // getorf prints [4 - 12] (1-based, STOP excluded); this API includes the STOP in the end
    // coordinate (INSDC feature table: CDS "location includes stop codon"), hence End = 11 + 3 = 14.
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

    // M11b — INV-5 / doc §6.1: an ORF that begins at a START codon but reaches the end of the
    // sequence with no in-frame STOP is emitted as an open (incomplete) ORF ending at the last base.
    // EMBOSS getorf -find 1 incomplete-ORF handling (WriteORF(start, pos+2): getorf reports [1 - 12]).
    // ATG AAA CCC GGG -> "MKPG", Start=0, End=11 (last base of the last complete codon, inclusive).
    [Test]
    public void FindOrfs_OrfRunsToSequenceEndWithoutStop_EmitsOpenOrf()
    {
        var dna = new DnaSequence("ATGAAACCCGGG");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "An ORF with no downstream STOP runs to the sequence end.");
        var orf = orfs[0];
        Assert.Multiple(() =>
        {
            Assert.That(orf.StartPosition, Is.EqualTo(0), "Start = first base of the ATG start codon.");
            Assert.That(orf.EndPosition, Is.EqualTo(11), "End = last base of the last complete codon (open ORF, no STOP).");
            Assert.That(orf.Frame, Is.EqualTo(1), "ORF is in forward frame +1.");
            Assert.That(orf.Protein.Sequence, Is.EqualTo("MKPG"),
                "Open ORF protein covers START to end with no terminating STOP residue.");
        });
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

    // C1 — alternative start codon TTG initiates an ORF; the initiator is reported as Met.
    // Corrected 2026-09 (B02 review): this test previously asserted "LKG" (TTG read as Leu).
    // EMBOSS getorf -methionine (default Y, getorf.acd: "Change initial START codons to
    // Methionine"; getorf.c appends 'M' for the START), Biopython translate(cds=True) and NCBI
    // The Genetic Codes ("the initiator codon ... is by default translated as methionine") all
    // give "MKG". Biopython: Seq("TTGAAAGGGTAA").translate(cds=True) == "MKG".
    [Test]
    public void FindOrfs_AlternativeStartCodonTtg_InitiatesOrfWithMethionine()
    {
        // GG TTG AAA GGG TAA CC : TTG start at index 2 (frame 3), TAA stop at indices 11-13.
        var dna = new DnaSequence("GGTTGAAAGGGTAACC");

        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "TTG is a START codon in NCBI standard table 1.");
        var orf = orfs[0];
        Assert.Multiple(() =>
        {
            Assert.That(orf.StartPosition, Is.EqualTo(2), "ORF starts at the TTG start codon.");
            Assert.That(orf.EndPosition, Is.EqualTo(13), "End = last base of the TAA stop codon.");
            Assert.That(orf.Frame, Is.EqualTo(3), "TTG at index 2 lies in forward frame +3.");
            Assert.That(orf.Protein.Sequence, Is.EqualTo("MKG"),
                "The initiator is read as Met whatever the start codon (getorf -methionine; Biopython cds=True).");
        });
    }

    // C1b — CTG (standard table 1 start) initiator also reported as Met; internal CTG stays Leu.
    // Biopython: Seq("CTGCTGTAA").translate(cds=True) == "ML".
    [Test]
    public void FindOrfs_AlternativeStartCodonCtg_InitiatorIsMethionine_InternalCtgIsLeucine()
    {
        var orf = Translator.FindOrfs(new DnaSequence("CTGCTGTAA"), minLength: 1, searchBothStrands: false).Single();

        Assert.Multiple(() =>
        {
            Assert.That(orf.Protein.Sequence, Is.EqualTo("ML"));
            Assert.That(orf.EndPosition, Is.EqualTo(8));
        });
    }

    // M11c — open ORF with a trailing partial codon ends at the last base of the last COMPLETE
    // codon (EMBOSS getorf WriteORF(start, pos+2): getorf reports [1 - 9] for this input).
    // Corrected 2026-09 (B02 review): the implementation previously reported End = Length-1 = 10
    // (NucleotideLength 11, not a whole number of codons).
    [Test]
    public void FindOrfs_OpenOrfWithTrailingPartialCodon_EndsAtLastCompleteCodon()
    {
        var dna = new DnaSequence("ATGAAACCCGG");

        var orf = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).Single();

        Assert.Multiple(() =>
        {
            Assert.That(orf.StartPosition, Is.EqualTo(0));
            Assert.That(orf.EndPosition, Is.EqualTo(8), "Last base of the last complete codon (CCC).");
            Assert.That(orf.NucleotideLength, Is.EqualTo(9), "Open ORF = 3 × protein length.");
            Assert.That(orf.Protein.Sequence, Is.EqualTo("MKP"));
        });
    }

    // M11d — same rule on the reverse strand: revcomp(CCGGGTTTCAT) = ATGAAACCCGG.
    [Test]
    public void FindOrfs_ReverseStrandOpenOrf_EndsAtLastCompleteCodon()
    {
        var orfs = Translator.FindOrfs(new DnaSequence("CCGGGTTTCAT"), minLength: 1, searchBothStrands: true).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(orfs[0].Frame, Is.EqualTo(-1));
            Assert.That(orfs[0].StartPosition, Is.EqualTo(0));
            Assert.That(orfs[0].EndPosition, Is.EqualTo(8));
            Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MKP"));
        });
    }

    // D1 — differential case locked to a Python port of EMBOSS getorf -find 1 (getorf.c
    // getorf_FindORFs, -methionine Y) with Biopython NCBI table 12 (CTG = Ser, but a START),
    // mapped to this API's convention (reverse strand in reverse-complement coordinates,
    // End includes the STOP when terminated). 35 nt: both ORFs run off the end (trailing
    // partial codons dropped) and the reverse-frame ORF starts at CTG -> 'M', not 'S'.
    [Test]
    public void FindOrfs_Table12_MatchesGetorfReference()
    {
        var dna = new DnaSequence("AGATTTTCATATTATGCAGAAAATCTACTTCGCCT");

        var orfs = Translator.FindOrfs(dna, GeneticCode.GetByTableNumber(12), minLength: 1)
            .OrderBy(o => o.Frame).Select(o => (o.Frame, o.StartPosition, o.EndPosition, o.Protein.Sequence))
            .ToList();

        Assert.That(orfs, Is.EqualTo(new[]
        {
            (-2, 16, 33, "MHNMKI"),
            (2, 13, 33, "MQKIYFA"),
        }));
    }

    // D2 — NCBI tables 27, 28, 31: every stop codon also codes for an amino acid (Biopython
    // CodonTable: stop_codons ⊂ forward_table), so no codon unambiguously ends an ORF. Rejected,
    // as Biopython rejects translate(to_stop=True) for these tables.
    [TestCase(27)]
    [TestCase(28)]
    [TestCase(31)]
    public void FindOrfs_DualCodingStopTable_Throws(int table)
    {
        var dna = new DnaSequence("ATGTGATAAGGG");

        Assert.Throws<System.ArgumentException>(
            () => Translator.FindOrfs(dna, GeneticCode.GetByTableNumber(table), minLength: 1));
    }

    #endregion
}
