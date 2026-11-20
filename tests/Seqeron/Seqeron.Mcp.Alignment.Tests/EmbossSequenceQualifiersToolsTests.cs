using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>B09 A5-7: the <c>sequenceQualifiers</c> parameter of the EMBOSS pairwise tools. Values are what the installed EMBOSS 6.6.0 prints (also locked by the 180-row oracles of the library).</summary>
[TestFixture]
public class EmbossSequenceQualifiersToolsTests
{
    private const string A = "ACGTACGTTAGCCGATTAGGCATCGATCGGATCG";
    private const string B = "CGGATCGATCAGGACTATCGGCTAACGTACGTTAGC";

    [Test]
    public void NeedleAlign_RangeAndReverse_ReportsTheOriginalNumbering()
    {
        // needle a.fa b.fa -sbegin1 3 -send1 20 -sreverse2 -slower1 -aformat3 pair
        var r = AlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "-sbegin1 3 -send1 20 -sreverse2 -slower1");
        Assert.That(r.Formatted, Does.Contain("a                  3 ------gtacgttagccgattagg-------------     20"));
        Assert.That(r.Formatted, Does.Contain("b                 36 GCTAACGTACGTTAGCCGA-TAGTCCTGATCGATCCG      1"));
        Assert.That(r.Score, Is.EqualTo(66));
    }

    [Test]
    public void NeedleAlign_WithoutQualifiers_IsUnchanged()
    {
        var plain = AlignmentTools.NeedleAlign(A, B);
        var empty = AlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "");
        Assert.That(empty.Score, Is.EqualTo(plain.Score));
        Assert.That(empty.AlignedSequence1, Is.EqualTo(plain.AlignedSequence1));
    }

    [Test]
    public void ForcedTypeContradictingTheInput_IsRefused()
    {
        // needle: the second sequence of a nucleotide first sequence must be nucleotide ("Sequence is not nucleic")
        Assert.Throws<ArgumentException>(() => AlignmentTools.NeedleAlign("ACGTACGT", "ACDEFGAIK", sequenceQualifiers: "-sbegin1 1"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WaterAlign(A, B, sequenceQualifiers: "-gapopen 3"));
    }

    [Test]
    public void Dotmatcher_Qualifiers_ScansTheWholeLengthOverTheCutText()
    {
        // dotmatcher da.fa db.fa -windowsize 5 -threshold 15 -matrixfile EDNAFULL -sbegin1 10: the first segment of the plain run (32 32 40 40) stays: -sbegin shifts the range to the front
        var r = AlignmentTools.Dotmatcher("ACGTTGCATGCATTAGCGATCGATCGGATTACAGGCATCG", "TTGCATGCATTAGCGATCGATAAGGATTACAGGGGCATCG", 5, 15, "EDNAFULL", sequenceQualifiers: "-sbegin1 10");
        Assert.That(r.Segments, Is.Not.Empty);
        Assert.That(r.Segments.Any(s => s.Start1 == 32 && s.Start2 == 32), Is.True);
    }

    [Test]
    public void Merger_RangesAreIgnored_AsTheProgramDoes()
    {
        // merger a.fa b.fa -send1 -39 -send2 34 ... : merger does not trim, the whole sequences are merged (installed 6.6.0, oracle row of emboss_seqopts_oracle.tsv)
        const string a = "TCCGCCAGGTACACCCCTTGTGGGAAGAGTGCCTAGACC", b = "TTGTGGGAAGAGTGCCTAGACCGCCCATCTTGCG";
        var plain = AlignmentTools.Merger(a, b);
        var ranged = AlignmentTools.Merger(a, b, sequenceQualifiers: "-send1 -39 -send2 34");
        Assert.That(ranged.MergedSequence, Is.EqualTo("tccgccaggtacaccccttgtgggaagagtgcctagaccgcccatcttgcg"));
        Assert.That(ranged.MergedSequence, Is.EqualTo(plain.MergedSequence));
        Assert.That(ranged.Conflicts, Is.Empty);
    }

    [Test]
    public void Merger_ReverseAndCase_AreApplied()
    {
        // merger -sreverse1 -slower1 -sbegin2 52 -send2 24 (installed 6.6.0): the reverse complement of a is merged; the first conflict is a 24 'G' / b 2 'c' using 'G'
        var r = AlignmentTools.Merger("TCCTTTCCAAACCCCCGCGCTCTTGTTAAGGTCAAATG", "GCGCTCTTGTTAAGGTCAAATGCTACTGATTGACAG", sequenceQualifiers: "-sreverse1 -slower1 -sbegin2 52 -send2 24");
        Assert.That(r.MergedSequence, Is.EqualTo("catttgaccttaacaagagcgcgGgGGTttgTTaaggTcaaatgctactgattgacag"));
        Assert.That(r.Conflicts[0], Is.EqualTo(new MergerConflictDto(24, "G", 2, "c", "G")));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Merger("ACGT", "ACGT", sequenceQualifiers: "-sbegin1"));
    }
}
