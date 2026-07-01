using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: TRANS-PROT-001 — sequence-to-protein translation (Translation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 63.
    //
    // API under test (Translator.Translate):
    //   Reads a DNA/RNA sequence codon-by-codon from a given reading frame and maps each
    //   codon through the genetic code. With toFirstStop = true it terminates at the first
    //   stop codon; otherwise stop codons are rendered as '*' and translation continues.
    //   Trailing 1–2 nucleotides that cannot form a codon are dropped.
    //
    // Relations (derived from the codon→residue mapping over a frame, NOT from output):
    //   • INV  (synonymous codon swap ⇒ same protein): the code is degenerate, so replacing
    //          any codon by a synonym of the SAME amino acid leaves the translated protein
    //          byte-for-byte identical.
    //   • COMP (stop codon ⇒ truncation): with toFirstStop the protein is exactly the full
    //          translation cut at its first '*', and is therefore invariant to whatever is
    //          appended after that stop codon.
    //   • SHIFT/INV (frame shift): reading at frame f re-partitions the codons, so
    //          Translate(seq, frame f) ≡ Translate(seq with its first f bases dropped, frame 0);
    //          because the partition changes, the three frames generally give DIFFERENT proteins.
    // ───────────────────────────────────────────────────────────────────────────

    #region TRANS-PROT-001 INV — synonymous codon substitution preserves the protein

    [Test]
    [Description("INV: the genetic code is degenerate, so swapping each codon for a synonym of the same amino acid yields an identical protein.")]
    public void Translate_SynonymousCodonSwap_PreservesProtein()
    {
        // Both spellings encode Met-Leu-Lys-Arg-Gly-Ser using DIFFERENT synonymous codons:
        //   variantA: AUG CUG AAA CGU GGU AGC
        //   variantB: AUG UUA AAG CGC GGA UCU
        const string variantA = "ATGCTGAAACGTGGTAGC";
        const string variantB = "ATGTTAAAGCGCGGATCT";

        variantA.Should().NotBe(variantB, because: "the two coding sequences must differ at the nucleotide level for the test to be meaningful");

        string proteinA = Translator.Translate(variantA).Sequence;
        string proteinB = Translator.Translate(variantB).Sequence;

        proteinB.Should().Be(proteinA,
            because: "synonymous codons encode the same amino acid, so codon swaps cannot change the protein");
        proteinA.Should().Be("MLKRGS", because: "the chosen codons spell Met-Leu-Lys-Arg-Gly-Ser");
    }

    #endregion

    #region TRANS-PROT-001 COMP — toFirstStop truncates at the first stop codon

    [Test]
    [Description("COMP: with toFirstStop the protein equals the full '*'-rendered translation cut at its first stop, and is invariant to any codons appended after that stop.")]
    public void Translate_ToFirstStop_TruncatesAtFirstStopAndIgnoresDownstream()
    {
        // AUG AAA UAA CCC GGG — a stop (UAA) sits after Met-Lys, followed by more codons.
        const string seq = "ATGAAATAACCCGGG";

        string full = Translator.Translate(seq, toFirstStop: false).Sequence;
        string truncated = Translator.Translate(seq, toFirstStop: true).Sequence;

        full.Should().Be("MK*PG", because: "without early termination every codon is rendered, the stop as '*'");

        int firstStop = full.IndexOf('*');
        truncated.Should().Be(full.Substring(0, firstStop),
            because: "toFirstStop terminates translation at the first stop codon");

        // Appending arbitrary codons after the stop cannot change the truncated protein.
        foreach (string tail in new[] { "", "TTTGGGAAA", "ATGATGATG", "TAGTAA" })
            Translator.Translate(seq + tail, toFirstStop: true).Sequence.Should().Be(truncated,
                because: "nothing downstream of the first stop is translated when toFirstStop is set");
    }

    #endregion

    #region TRANS-PROT-001 SHIFT/INV — frame shift re-partitions codons

    [Test]
    [Description("SHIFT: translating at frame f equals translating the f-base-dropped sequence at frame 0 — the exact re-partition identity that defines a reading frame.")]
    public void Translate_FrameShift_EqualsOffsetSequenceAtFrameZero()
    {
        const string seq = "ATGCATGCATGCATGCATGC";

        for (int frame = 0; frame <= 2; frame++)
        {
            string framed = Translator.Translate(seq, frame: frame).Sequence;
            string offset = Translator.Translate(seq.Substring(frame), frame: 0).Sequence;

            framed.Should().Be(offset,
                because: $"frame {frame} starts the codon partition at base {frame}, identical to dropping the first {frame} bases and reading at frame 0");
        }
    }

    [Test]
    [Description("SHIFT: because each frame imposes a different codon partition, the three reading frames of a designed sequence give pairwise-distinct proteins.")]
    public void Translate_DifferentFrames_GiveDifferentProteins()
    {
        // A sequence engineered so the three frames read disjoint codon sets.
        const string seq = "ATGCATGCATGCATGCATGC";

        string f0 = Translator.Translate(seq, frame: 0).Sequence;
        string f1 = Translator.Translate(seq, frame: 1).Sequence;
        string f2 = Translator.Translate(seq, frame: 2).Sequence;

        f0.Should().NotBe(f1, because: "shifting the frame by one base changes every codon, so the protein differs");
        f0.Should().NotBe(f2, because: "shifting the frame by two bases changes every codon, so the protein differs");
        f1.Should().NotBe(f2, because: "frames 1 and 2 impose different codon partitions");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: TRANS-SIXFRAME-001 — six-frame translation (Translation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 223.
    //
    // API under test (Translator.TranslateSixFrames):
    //   Translates a DNA sequence in all six reading frames: forward frames +1/+2/+3 at offsets
    //   0/1/2, and reverse frames −1/−2/−3 = forward frames of the reverse complement.
    //
    // Relations (derived from the frame definition, NOT from output):
    //   • P   (exactly 6 frames): the result has exactly the six keys {+1,+2,+3,−1,−2,−3}; frame 0
    //         does not exist.
    //   • INV (frames 4–6 = translation of revcomp): each reverse frame −k of a sequence equals the
    //         forward frame +k of its reverse complement (and symmetrically for the forward frames),
    //         since reverse-complementing twice is the identity.
    // ───────────────────────────────────────────────────────────────────────────

    #region TRANS-SIXFRAME-001 — Helpers

    // A sequence whose length is not a multiple of 3, to also exercise the trailing-codon handling.
    private const string SixFrameSeq = "ATGGCCTAAGGGTTTACGTACCGA"; // 24 nt; revcomp is non-palindromic

    private static string Frame(IReadOnlyDictionary<int, ProteinSequence> frames, int k) => frames[k].Sequence;

    #endregion

    #region TRANS-SIXFRAME-001 P — there are exactly six frames

    [Test]
    [Description("P: TranslateSixFrames returns exactly the six reading frames {+1,+2,+3,−1,−2,−3} and no frame 0.")]
    public void SixFrames_HasExactlySixFramesNoZero()
    {
        var frames = Translator.TranslateSixFrames(new DnaSequence(SixFrameSeq));

        frames.Keys.Should().BeEquivalentTo(new[] { 1, 2, 3, -1, -2, -3 },
            because: "three forward and three reverse reading frames span every codon partition of both strands");
        frames.Keys.Should().NotContain(0, because: "frame 0 is undefined — frames are 1-based per strand");
    }

    #endregion

    #region TRANS-SIXFRAME-001 INV — reverse frames are the forward frames of the reverse complement

    [Test]
    [Description("INV: each reverse frame −k equals the forward frame +k of the reverse complement (and vice versa), because revcomp∘revcomp is the identity.")]
    public void SixFrames_ReverseFrames_EqualForwardFramesOfRevComp()
    {
        var dna = new DnaSequence(SixFrameSeq);
        var revComp = dna.ReverseComplement();

        var frames = Translator.TranslateSixFrames(dna);
        var revFrames = Translator.TranslateSixFrames(revComp);

        for (int k = 1; k <= 3; k++)
        {
            Frame(frames, -k).Should().Be(Frame(revFrames, k),
                because: $"reverse frame −{k} of the sequence is forward frame +{k} of its reverse complement");
            Frame(frames, k).Should().Be(Frame(revFrames, -k),
                because: $"forward frame +{k} of the sequence is reverse frame −{k} of its reverse complement (revcomp is an involution)");
        }

        // Non-vacuous: forward and reverse translations genuinely differ for this sequence.
        Frame(frames, 1).Should().NotBe(Frame(frames, -1),
            because: "the sequence is not a reverse-complement palindrome, so its strands translate differently");
    }

    #endregion
}
