namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Translation area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Translation")]
public class TranslationCombinatorialTests
{
    // Contains internal stops (TGA, TAA) and varied codons so frame/table/stop axes all bite.
    private const string Seq = "ATGAAATGAGGGTGATTTAAACCC"; // length 24

    private static string IndepTranslate(string seq, GeneticCode gc, int frame, bool toFirstStop)
    {
        string rna = seq.Replace('T', 'U');
        var sb = new System.Text.StringBuilder();
        for (int i = frame; i + 3 <= rna.Length; i += 3)
        {
            char aa = gc.Translate(rna.Substring(i, 3));
            if (toFirstStop && aa == '*') break;
            sb.Append(aa);
        }
        return sb.ToString();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-PROT-001 — Protein translation (Translation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 63.
    // Spec: tests/TestSpecs/TRANS-PROT-001.md (canonical Translator.Translate).
    // Dimensions: frame(3) × tableId(4) × stopHandling(2: stop/readthrough). Grid 3×4×2 = 24.
    //
    // Model (central dogma): translation reads codons from a chosen reading frame using a chosen
    // NCBI genetic-code table; "readthrough" emits a '*' for each stop and continues, while
    // "to first stop" truncates the protein before the first stop. Tables differ (e.g. UGA is a
    // stop in table 1 but Tryptophan in table 2). The checklist's table 5 is not implemented;
    // the fourth available table is 11 (bacterial/plastid).
    //
    // The combinatorial point: frame offset, code table and stop handling compose — the output
    // equals an independent codon-walk for every (frame, table, stop) cell; readthrough yields
    // one residue per complete codon, while to-first-stop never contains a '*'.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void TransProt_ComposesFrameTableAndStopHandling(
        [Values(0, 1, 2)] int frame,
        [Values(1, 2, 3, 11)] int tableId,
        [Values(true, false)] bool toFirstStop)
    {
        var gc = GeneticCode.GetByTableNumber(tableId);

        string actual = Translator.Translate(Seq, gc, frame, toFirstStop).Sequence;

        actual.Should().Be(IndepTranslate(Seq, gc, frame, toFirstStop),
            "translation equals an independent codon walk for this frame/table/stop");

        if (!toFirstStop)
            actual.Length.Should().Be((Seq.Length - frame) / 3, "readthrough emits one residue per complete codon");
        else
            actual.Should().NotContain("*", "to-first-stop truncates before any stop");
    }

    /// <summary>
    /// Interaction witness: the code table changes the amino acid for a codon — UGA is a stop in
    /// the standard table but Tryptophan in the vertebrate-mitochondrial table.
    /// </summary>
    [Test]
    public void TransProt_TableChangesCodonMeaning()
    {
        const string s = "ATGTGAAAA"; // M, TGA, K  (frame 0)
        Translator.Translate(s, GeneticCode.GetByTableNumber(1), 0, toFirstStop: false).Sequence
            .Should().Be("M*K", "UGA is a stop in the standard code");
        Translator.Translate(s, GeneticCode.GetByTableNumber(2), 0, toFirstStop: false).Sequence
            .Should().Be("MWK", "UGA is Tryptophan in the vertebrate mitochondrial code");
    }

    /// <summary>
    /// Interaction witness: the reading frame selects different codons, generally yielding a
    /// different protein.
    /// </summary>
    [Test]
    public void TransProt_FrameShiftsTheReading()
    {
        var gc = GeneticCode.Standard;
        string f0 = Translator.Translate(Seq, gc, 0, toFirstStop: false).Sequence;
        string f1 = Translator.Translate(Seq, gc, 1, toFirstStop: false).Sequence;
        f0.Should().NotBe(f1, "a frame shift re-reads the codons");
    }

    /// <summary>
    /// Interaction witness: with an internal stop, to-first-stop is a strict prefix of the
    /// readthrough protein (truncated at the first '*').
    /// </summary>
    [Test]
    public void TransProt_StopHandling_TruncatesAtFirstStop()
    {
        var gc = GeneticCode.Standard;
        string readthrough = Translator.Translate(Seq, gc, 0, toFirstStop: false).Sequence;
        string truncated = Translator.Translate(Seq, gc, 0, toFirstStop: true).Sequence;

        readthrough.Should().Contain("*", "the sequence has an in-frame stop");
        truncated.Should().Be(readthrough[..readthrough.IndexOf('*')], "to-first-stop is the pre-stop prefix");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-SIXFRAME-001 — Six-frame translation (Translation)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 223.
    // Spec: tests/TestSpecs/TRANS-SIXFRAME-001.md (canonical Translator.TranslateSixFrames). ADVANCED §10.
    // Dimensions: tableId(4) × seqLen(3). Grid 4×3 = 12 (full, exhaustive ⊇ pairwise).
    //
    // Model (Biopython six_frame_translations): the three forward frames (offsets 0/1/2) and the three
    // reverse-complement frames are each translated codon-by-codon with the chosen genetic code; trailing
    // partial codons are dropped.
    //
    // The combinatorial point: across genetic-code table and length, all six frames {±1,±2,±3} are
    // present and each frame's protein equals an independent codon-by-codon translation (forward frames
    // on the sequence, reverse frames on its reverse complement).
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void TransSixFrame_AllFramesMatchIndependentTranslation_AcrossTableAndLength(
        [Values(1, 2, 3, 11)] int tableId,
        [Values(12, 24, 36)] int seqLen)
    {
        var gc = GeneticCode.GetByTableNumber(tableId);
        var dna = new DnaSequence(Seq2(seqLen));
        string revComp = dna.ReverseComplement().Sequence;

        var frames = Translator.TranslateSixFrames(dna, gc);

        frames.Keys.Should().BeEquivalentTo(new[] { 1, 2, 3, -1, -2, -3 }, "all six reading frames are returned");
        for (int offset = 0; offset < 3; offset++)
        {
            frames[offset + 1].Sequence.Should().Be(IndepTranslate(dna.Sequence, gc, offset, false),
                "forward frame +{0} matches an independent translation", offset + 1);
            frames[-(offset + 1)].Sequence.Should().Be(IndepTranslate(revComp, gc, offset, false),
                "reverse frame −{0} translates the reverse complement", offset + 1);
        }
    }

    /// <summary>A deterministic ACGT sequence of the requested length (≥ 3).</summary>
    private static string Seq2(int n) => string.Concat(Enumerable.Repeat("ATGAAACCCGGGTTT", n / 15 + 1))[..n];

    /// <summary>
    /// Interaction witness — the genetic-code table changes a frame's translation: UGA is a stop under
    /// the standard code (table 1) but tryptophan (W) under the vertebrate-mitochondrial code (table 2).
    /// </summary>
    [Test]
    public void TransSixFrame_TableChangesTranslation()
    {
        var dna = new DnaSequence("TGATGATGA"); // frame +1 codons: TGA TGA TGA
        string standard = Translator.TranslateSixFrames(dna, GeneticCode.GetByTableNumber(1))[1].Sequence;
        string mito = Translator.TranslateSixFrames(dna, GeneticCode.GetByTableNumber(2))[1].Sequence;

        standard.Should().Be("***", "UGA is a stop under the standard code");
        mito.Should().Be("WWW", "UGA codes tryptophan under the vertebrate-mitochondrial code");
    }
}
