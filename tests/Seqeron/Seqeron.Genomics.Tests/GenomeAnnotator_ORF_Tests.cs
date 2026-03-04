using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Tests.Helpers;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for ORF Detection (ANNOT-ORF-001).
/// Tests GenomeAnnotator.FindOrfs and FindLongestOrfsPerFrame.
///
/// Evidence: Wikipedia (Open reading frame), Rosalind ORF, NCBI ORF Finder,
/// Deonier et al. (2005), Claverie (1997).
/// </summary>
[TestFixture]
public class GenomeAnnotator_ORF_Tests
{
    #region Test Helpers

    /// <summary>
    /// Creates a valid ORF DNA sequence: startCodon + coding + stopCodon.
    /// </summary>
    private static string CreateOrf(string startCodon, int codingCodons, string stopCodon)
    {
        // Each coding codon is 3 nucleotides; use AAA (Lysine) for simplicity
        var sb = new StringBuilder(startCodon);
        for (int i = 0; i < codingCodons; i++)
        {
            sb.Append("AAA");
        }
        sb.Append(stopCodon);
        return sb.ToString();
    }

    private static string GetReverseComplement(string sequence)
        => TestSequenceHelpers.GetReverseComplement(sequence);

    private static readonly HashSet<string> ValidStartCodons = new(StringComparer.OrdinalIgnoreCase)
    {
        "ATG", "GTG", "TTG"
    };

    private static readonly HashSet<string> ValidStopCodons = new(StringComparer.OrdinalIgnoreCase)
    {
        "TAA", "TAG", "TGA"
    };

    #endregion

    #region M01-M05: Basic ORF Detection

    /// <summary>
    /// M01: Simple ORF with ATG start and TAA stop is detected.
    /// Evidence: Wikipedia, Rosalind ORF definition.
    /// Implementation note: Protein includes translated stop codon (*).
    /// </summary>
    [Test]
    public void FindOrfs_SimpleAtgTaaOrf_DetectsOrf()
    {
        // ATG + 50 codons + TAA = 52 characters in protein (M + 50 aa + *)
        string orf = CreateOrf("ATG", 50, "TAA");

        var orfs = GenomeAnnotator.FindOrfs(orf, minLength: 10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1), "Should find exactly one ORF");
            Assert.That(orfs[0].Start, Is.EqualTo(0), "ORF should start at position 0");
            Assert.That(orfs[0].End, Is.EqualTo(orf.Length), "ORF should end at sequence length");
            Assert.That(orfs[0].ProteinSequence, Does.StartWith("M"), "Protein should start with M");
            Assert.That(orfs[0].ProteinSequence, Does.EndWith("*"), "Protein should end with stop codon");
            Assert.That(orfs[0].ProteinSequence.TrimEnd('*').Length, Is.EqualTo(51), "Protein length (excluding stop) should be 51 aa");
        });
    }

    /// <summary>
    /// M02: Empty sequence returns empty collection.
    /// Evidence: Standard edge case handling.
    /// </summary>
    [Test]
    public void FindOrfs_EmptySequence_ReturnsEmpty()
    {
        var orfs = GenomeAnnotator.FindOrfs("", minLength: 1).ToList();

        Assert.That(orfs, Is.Empty);
    }

    /// <summary>
    /// M03: Sequence without start codon returns empty when requireStartCodon=true.
    /// Evidence: Algorithm definition.
    /// </summary>
    [Test]
    public void FindOrfs_NoStartCodon_RequireStart_ReturnsEmpty()
    {
        string noStart = "GGGGGGGGGGTAAGGGGGGGGG";

        var orfs = GenomeAnnotator.FindOrfs(noStart, minLength: 1, requireStartCodon: true).ToList();

        Assert.That(orfs, Is.Empty);
    }

    /// <summary>
    /// M04: Sequence without stop codon returns empty when requireStartCodon=true.
    /// Evidence: Algorithm definition - ORF requires stop codon.
    /// </summary>
    [Test]
    public void FindOrfs_NoStopCodon_RequireStart_ReturnsEmpty()
    {
        string noStop = "ATGAAAAAAAAAAAAAAAAAAA"; // ATG but no stop

        var orfs = GenomeAnnotator.FindOrfs(noStop, minLength: 1, requireStartCodon: true).ToList();

        Assert.That(orfs, Is.Empty);
    }

    /// <summary>
    /// M05: Alternative start codons GTG and TTG are detected.
    /// Evidence: Wikipedia, NCBI ORF Finder supports alternative initiation codons.
    /// </summary>
    [TestCase("GTG", Description = "GTG alternative start")]
    [TestCase("TTG", Description = "TTG alternative start")]
    public void FindOrfs_AlternativeStartCodons_Detected(string startCodon)
    {
        string orf = CreateOrf(startCodon, 50, "TAA");

        var orfs = GenomeAnnotator.FindOrfs(orf, minLength: 10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.GreaterThan(0), $"Should detect ORF with {startCodon} start");
            Assert.That(orfs[0].Sequence.StartsWith(startCodon, StringComparison.OrdinalIgnoreCase),
                $"ORF should start with {startCodon}");
        });
    }

    #endregion

    #region M06: Minimum Length Filtering

    /// <summary>
    /// M06: ORFs below minLength are excluded.
    /// Evidence: NCBI ORF Finder, Claverie (1997).
    /// </summary>
    [Test]
    public void FindOrfs_BelowMinLength_Excluded()
    {
        // ATG + 5 codons + TAA = 6 amino acids
        string shortOrf = CreateOrf("ATG", 5, "TAA");

        var orfs = GenomeAnnotator.FindOrfs(shortOrf, minLength: 10).ToList();

        Assert.That(orfs, Is.Empty, "ORF with 6 aa should be excluded when minLength=10");
    }

    /// <summary>
    /// M06b: ORFs at exactly minLength are included.
    /// </summary>
    [Test]
    public void FindOrfs_ExactlyMinLength_Included()
    {
        // ATG + 9 codons + TAA = 10 amino acids
        string orf = CreateOrf("ATG", 9, "TAA");

        var orfs = GenomeAnnotator.FindOrfs(orf, minLength: 10).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1), "ORF with exactly 10 aa should be included");
    }

    #endregion

    #region M07-M09: Six-Frame Search

    /// <summary>
    /// M07: ORFs found across multiple reading frames from a single sequence.
    /// Evidence: Wikipedia six-frame translation, Rosalind.
    /// </summary>
    [Test]
    public void FindOrfs_SixFrameSearch_FindsOrfsInMultipleFrames()
    {
        // Single sequence with ORFs in frame 1 and frame 2:
        // Frame 1: ATG at position 0 → (0 % 3) + 1 = frame 1, 51 nt ORF
        // Frame 2: ATG at position 52 → (52 % 3) + 1 = frame 2, 51 nt ORF
        string sequence = CreateOrf("ATG", 15, "TAA") + "G" + CreateOrf("ATG", 15, "TAA");

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 5, searchBothStrands: false).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(2), "Should find exactly 2 ORFs");
            Assert.That(orfs.Select(o => o.Frame).Distinct().Count(), Is.EqualTo(2),
                "ORFs should span 2 different frames");
            Assert.That(orfs.Any(o => o.Frame == 1), Is.True, "Should find ORF in frame 1");
            Assert.That(orfs.Any(o => o.Frame == 2), Is.True, "Should find ORF in frame 2");
        });
    }

    /// <summary>
    /// M08: ORFs on reverse complement strand are detected.
    /// Evidence: Rosalind ORF problem.
    /// </summary>
    [Test]
    public void FindOrfs_ReverseStrand_FindsOrfs()
    {
        // Create ORF on reverse strand: reverse complement of ATG...TAA
        string forwardOrf = CreateOrf("ATG", 30, "TAA");
        string reverseOrfSequence = GetReverseComplement(forwardOrf);

        var orfs = GenomeAnnotator.FindOrfs(reverseOrfSequence, minLength: 10, searchBothStrands: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1), "Should find exactly 1 ORF on reverse strand");
            Assert.That(orfs[0].IsReverseComplement, Is.True,
                "ORF should be from reverse complement strand");
            Assert.That(orfs[0].ProteinSequence, Does.StartWith("M"), "Protein should start with M");
            Assert.That(orfs[0].ProteinSequence, Does.EndWith("*"), "Protein should end with stop");
        });
    }

    /// <summary>
    /// M08b: When searchBothStrands=false, reverse strand is not searched.
    /// </summary>
    [Test]
    public void FindOrfs_ForwardOnly_DoesNotSearchReverse()
    {
        string forwardOrf = CreateOrf("ATG", 30, "TAA");
        string reverseOrfSequence = GetReverseComplement(forwardOrf);

        var orfs = GenomeAnnotator.FindOrfs(reverseOrfSequence, minLength: 10, searchBothStrands: false).ToList();

        Assert.That(orfs.Any(o => o.IsReverseComplement), Is.False);
    }

    /// <summary>
    /// M09: Each frame is reported separately with correct frame number.
    /// Evidence: Wikipedia - three frames per strand.
    /// </summary>
    [TestCase(0, 1, Description = "ATG at position 0 → frame 1")]
    [TestCase(1, 2, Description = "ATG at position 1 → frame 2")]
    [TestCase(2, 3, Description = "ATG at position 2 → frame 3")]
    public void FindOrfs_FrameNumber_CorrectlyAssigned(int offset, int expectedFrame)
    {
        string padding = new string('G', offset);
        string orf = CreateOrf("ATG", 30, "TAA");
        string sequence = padding + orf;

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 10, searchBothStrands: false).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1),
                $"Should find exactly 1 ORF with offset {offset}");
            Assert.That(orfs[0].Frame, Is.EqualTo(expectedFrame),
                $"ORF at offset {offset} should be in frame {expectedFrame}");
            Assert.That(orfs[0].Start, Is.EqualTo(offset),
                $"ORF should start at position {offset}");
        });
    }

    #endregion

    #region M10-M11: Multiple and Overlapping ORFs

    /// <summary>
    /// M10: Multiple non-overlapping ORFs in same sequence are all returned.
    /// </summary>
    [Test]
    public void FindOrfs_MultipleOrfs_AllReturned()
    {
        // Two separate ORFs in frame 1 (GGG spacer preserves frame alignment)
        string orf1 = CreateOrf("ATG", 20, "TAA"); // 66 nt, positions 0-65
        string orf2 = CreateOrf("ATG", 20, "TAG"); // 66 nt, positions 69-134
        string sequence = orf1 + "GGG" + orf2;

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 10, searchBothStrands: false).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(2), "Should find exactly 2 ORFs");
            Assert.That(orfs[0].Start, Is.EqualTo(0), "First ORF starts at position 0");
            Assert.That(orfs[1].Start, Is.EqualTo(69), "Second ORF starts at position 69");
        });
    }

    /// <summary>
    /// M11: Overlapping ORFs (nested start codons) are both reported if qualifying.
    /// </summary>
    [Test]
    public void FindOrfs_NestedOrfs_BothReportedIfQualifying()
    {
        // Outer: ATG at 0 + 30 A's + ATG at 33 + 30 A's + TAA = 69 nt
        // Inner ORF shares the same stop codon
        string sequence = "ATG" + "AAA".PadRight(30, 'A').Substring(0, 30) +
                          "ATG" + "AAA".PadRight(30, 'A').Substring(0, 30) + "TAA";

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 5, searchBothStrands: false).ToList();
        var frame1Orfs = orfs.Where(o => o.Frame == 1 && !o.IsReverseComplement)
                             .OrderBy(o => o.Start).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(frame1Orfs, Has.Count.EqualTo(2), "Should find both outer and nested ORF");
            Assert.That(frame1Orfs[0].Start, Is.EqualTo(0), "Outer ORF starts at position 0");
            Assert.That(frame1Orfs[1].Start, Is.EqualTo(33), "Inner ORF starts at position 33");
            Assert.That(frame1Orfs[0].End, Is.EqualTo(frame1Orfs[1].End),
                "Both ORFs share the same stop codon");
        });
    }

    #endregion

    #region M12: Rosalind Reference Dataset

    /// <summary>
    /// M12: Rosalind sample dataset produces expected proteins.
    /// Evidence: Rosalind ORF problem sample dataset.
    /// </summary>
    [Test]
    public void FindOrfs_RosalindDataset_FindsExpectedProteins()
    {
        // Rosalind sample dataset
        string rosalindSequence = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG";

        // Expected distinct proteins from Rosalind
        var expectedProteins = new HashSet<string>
        {
            "MLLGSFRLIPKETLIQVAGSSPCNLS",
            "M",
            "MGMTPRLGLESLLE",
            "MTPRLGLESLLE"
        };

        var orfs = GenomeAnnotator.FindOrfs(rosalindSequence, minLength: 1, searchBothStrands: true).ToList();
        var foundProteins = orfs.Select(o => o.ProteinSequence.TrimEnd('*')).Distinct().ToHashSet();

        // Check that all expected proteins are found
        foreach (var expected in expectedProteins)
        {
            Assert.That(foundProteins.Contains(expected), Is.True,
                $"Expected protein '{expected}' not found");
        }
    }

    #endregion

    #region M13-M15: ORF Invariants

    /// <summary>
    /// M13: All returned ORFs start with a valid start codon.
    /// Invariant: ORF must begin with ATG, GTG, or TTG.
    /// </summary>
    [Test]
    public void FindOrfs_Invariant_StartsWithStartCodon()
    {
        string sequence = CreateOrf("ATG", 50, "TAA") + "GG" + CreateOrf("GTG", 30, "TAG");

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 10, requireStartCodon: true).ToList();

        Assert.Multiple(() =>
        {
            foreach (var orf in orfs)
            {
                string firstCodon = orf.Sequence.Substring(0, 3).ToUpperInvariant();
                Assert.That(ValidStartCodons.Contains(firstCodon), Is.True,
                    $"ORF at {orf.Start} should start with valid start codon, got '{firstCodon}'");
            }
        });
    }

    /// <summary>
    /// M14: All returned ORFs end with a valid stop codon.
    /// Invariant: ORF must end with TAA, TAG, or TGA.
    /// </summary>
    [Test]
    public void FindOrfs_Invariant_EndsWithStopCodon()
    {
        string sequence = CreateOrf("ATG", 50, "TAA") + "GG" + CreateOrf("ATG", 30, "TGA");

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 10, requireStartCodon: true).ToList();

        Assert.Multiple(() =>
        {
            foreach (var orf in orfs)
            {
                string lastCodon = orf.Sequence.Substring(orf.Sequence.Length - 3, 3).ToUpperInvariant();
                Assert.That(ValidStopCodons.Contains(lastCodon), Is.True,
                    $"ORF at {orf.Start} should end with valid stop codon, got '{lastCodon}'");
            }
        });
    }

    /// <summary>
    /// M15: All ORF nucleotide lengths are divisible by 3.
    /// Invariant: Frame integrity.
    /// </summary>
    [Test]
    public void FindOrfs_Invariant_LengthDivisibleBy3()
    {
        string sequence = CreateOrf("ATG", 50, "TAA") + "G" + CreateOrf("ATG", 30, "TAG");

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 10).ToList();

        Assert.Multiple(() =>
        {
            foreach (var orf in orfs)
            {
                Assert.That(orf.Sequence.Length % 3, Is.EqualTo(0),
                    $"ORF nucleotide length {orf.Sequence.Length} should be divisible by 3");
            }
        });
    }

    #endregion

    #region M16-M17: FindLongestOrfsPerFrame

    /// <summary>
    /// M16: FindLongestOrfsPerFrame returns all 6 frame keys.
    /// </summary>
    [Test]
    public void FindLongestOrfsPerFrame_BothStrands_Returns6Keys()
    {
        string sequence = CreateOrf("ATG", 50, "TAA");

        var result = GenomeAnnotator.FindLongestOrfsPerFrame(sequence, searchBothStrands: true);

        Assert.Multiple(() =>
        {
            Assert.That(result.ContainsKey(1), Is.True, "Should contain frame +1");
            Assert.That(result.ContainsKey(2), Is.True, "Should contain frame +2");
            Assert.That(result.ContainsKey(3), Is.True, "Should contain frame +3");
            Assert.That(result.ContainsKey(-1), Is.True, "Should contain frame -1");
            Assert.That(result.ContainsKey(-2), Is.True, "Should contain frame -2");
            Assert.That(result.ContainsKey(-3), Is.True, "Should contain frame -3");
        });
    }

    /// <summary>
    /// M16b: Forward-only search returns only 3 keys.
    /// </summary>
    [Test]
    public void FindLongestOrfsPerFrame_ForwardOnly_Returns3Keys()
    {
        string sequence = CreateOrf("ATG", 50, "TAA");

        var result = GenomeAnnotator.FindLongestOrfsPerFrame(sequence, searchBothStrands: false);

        Assert.Multiple(() =>
        {
            Assert.That(result.ContainsKey(1), Is.True);
            Assert.That(result.ContainsKey(2), Is.True);
            Assert.That(result.ContainsKey(3), Is.True);
            Assert.That(result.ContainsKey(-1), Is.False);
            Assert.That(result.ContainsKey(-2), Is.False);
            Assert.That(result.ContainsKey(-3), Is.False);
        });
    }

    /// <summary>
    /// M17: Returns the longest ORF per frame when multiple exist.
    /// Implementation note: Protein includes translated stop codon (*).
    /// </summary>
    [Test]
    public void FindLongestOrfsPerFrame_ReturnsLongestPerFrame()
    {
        // Create sequence with short ORF then long ORF in same frame
        string shortOrf = CreateOrf("ATG", 10, "TAA");
        string longOrf = CreateOrf("ATG", 50, "TAA");
        string sequence = shortOrf + "GGG" + longOrf; // Both in frame 1

        var result = GenomeAnnotator.FindLongestOrfsPerFrame(sequence, searchBothStrands: false);

        // Protein length is 52: M + 50 aa + * (stop)
        Assert.That(result[1]?.ProteinSequence.TrimEnd('*').Length, Is.EqualTo(51),
            "Should return the longer ORF (51 aa excluding stop) for frame 1");
    }

    #endregion

    #region S01-S02: Case Handling

    /// <summary>
    /// S01: Lowercase input is handled correctly.
    /// </summary>
    [Test]
    public void FindOrfs_LowercaseInput_HandledCorrectly()
    {
        string uppercaseOrf = CreateOrf("ATG", 30, "TAA");
        string lowercaseOrf = uppercaseOrf.ToLowerInvariant();

        var uppercaseOrfs = GenomeAnnotator.FindOrfs(uppercaseOrf, minLength: 10).ToList();
        var lowercaseOrfs = GenomeAnnotator.FindOrfs(lowercaseOrf, minLength: 10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(lowercaseOrfs, Has.Count.EqualTo(1), "Should detect exactly 1 ORF in lowercase");
            Assert.That(lowercaseOrfs.Count, Is.EqualTo(uppercaseOrfs.Count),
                "Same ORF count as uppercase");
            Assert.That(lowercaseOrfs[0].ProteinSequence, Is.EqualTo(uppercaseOrfs[0].ProteinSequence),
                "Same protein as uppercase");
        });
    }

    /// <summary>
    /// S02: Mixed case input is handled correctly.
    /// </summary>
    [Test]
    public void FindOrfs_MixedCaseInput_HandledCorrectly()
    {
        // ATG + AAA + AAA + TAA = 4 codons → protein MKK*
        string mixedCase = "AtGaAaAaAtAa";

        var orfs = GenomeAnnotator.FindOrfs(mixedCase, minLength: 1, searchBothStrands: false).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1), "Should find exactly 1 ORF in mixed case");
            Assert.That(orfs[0].Frame, Is.EqualTo(1), "ORF should be in frame 1");
            Assert.That(orfs[0].ProteinSequence.TrimEnd('*'), Is.EqualTo("MKK"),
                "Protein should be MKK (Met + 2×Lys)");
        });
    }

    #endregion

    #region Stop Codon Variants (supplements M14)

    /// <summary>
    /// All three stop codons (TAA, TAG, TGA) are recognized.
    /// </summary>
    [TestCase("TAA", Description = "TAA (ochre) stop codon")]
    [TestCase("TAG", Description = "TAG (amber) stop codon")]
    [TestCase("TGA", Description = "TGA (opal) stop codon")]
    public void FindOrfs_AllStopCodons_Recognized(string stopCodon)
    {
        string orf = CreateOrf("ATG", 30, stopCodon);

        var orfs = GenomeAnnotator.FindOrfs(orf, minLength: 10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1), $"Should detect exactly 1 ORF with {stopCodon} stop");
            Assert.That(orfs[0].Sequence.EndsWith(stopCodon, StringComparison.OrdinalIgnoreCase),
                $"ORF should end with {stopCodon}");
        });
    }

    #endregion

    #region S03: Very Long Sequence (10kb+)

    /// <summary>
    /// S03: Very long sequence (10kb+) is handled correctly.
    /// Evidence: NCBI ORF Finder processes large sequences; performance baseline.
    /// </summary>
    [Test]
    public void FindOrfs_VeryLongSequence_CompletesCorrectly()
    {
        // 10,200 bp ORF: ATG + 3,398 coding codons (AAA) + TAA
        // Plus 2,000 nt non-ORF filler (GCGC repeats, no start codons)
        var sb = new StringBuilder();
        sb.Append(CreateOrf("ATG", 3398, "TAA")); // 3 + 3398*3 + 3 = 10,200 nt
        for (int i = 0; i < 500; i++)
            sb.Append("GCGC");

        string sequence = sb.ToString(); // 12,200 bp total

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 100, searchBothStrands: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1), "Should find exactly 1 ORF in 10kb+ sequence");
            Assert.That(orfs[0].Sequence.Length, Is.EqualTo(10200),
                "ORF should be 10,200 nt (ATG + 3398 codons + TAA)");
            Assert.That(orfs[0].ProteinSequence, Does.StartWith("M"),
                "Protein should start with M");
        });
    }

    #endregion

    #region S04: Sequences with N Characters

    /// <summary>
    /// S04a: Codon containing N does not match start codon.
    /// Evidence: NCBI C++ Toolkit orf.cpp — N is treated as unknown;
    /// codons containing N do not match ATG/GTG/TTG.
    /// </summary>
    [Test]
    public void FindOrfs_NInStartCodon_NotRecognizedAsStart()
    {
        // NTG should not match ATG (N is unknown nucleotide)
        string sequence = "NTG" + new string('A', 30 * 3) + "TAA";

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 1, searchBothStrands: false, requireStartCodon: true).ToList();

        Assert.That(orfs, Is.Empty, "Codon with N should not match start codon");
    }

    /// <summary>
    /// S04b: Codon containing N does not match stop codon.
    /// Evidence: NCBI C++ Toolkit orf.cpp — N prevents codon matching.
    /// An N-containing codon (e.g., NAA) is neither start nor stop.
    /// </summary>
    [Test]
    public void FindOrfs_NInStopCodon_NotRecognizedAsStop()
    {
        // ATG + coding + NAA (not a stop) + more coding + TAA (real stop)
        string sequence = "ATG" + new string('A', 10 * 3) + "NAA" + new string('A', 10 * 3) + "TAA";

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 1, searchBothStrands: false).ToList();
        var frame1Orfs = orfs.Where(o => o.Frame == 1 && !o.IsReverseComplement).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(frame1Orfs, Has.Count.EqualTo(1),
                "Should find exactly 1 ORF continuing through N-containing codon");
            Assert.That(frame1Orfs[0].Sequence, Does.EndWith("TAA"),
                "ORF should end at real stop codon TAA, not NAA");
        });
    }

    /// <summary>
    /// S04c: N in middle of ORF coding region does not break the ORF.
    /// Evidence: NCBI C++ Toolkit — individual N chars are not treated as stops.
    /// </summary>
    [Test]
    public void FindOrfs_NInCodingRegion_OrfContinues()
    {
        // ATG + AAA(x5) + NCC(not a stop) + AAA(x5) + TAA
        string sequence = "ATG" + "AAAAAAAAAAAAAAA" + "NCC" + "AAAAAAAAAAAAAAA" + "TAA";

        var orfs = GenomeAnnotator.FindOrfs(sequence, minLength: 1, searchBothStrands: false).ToList();
        var frame1Orfs = orfs.Where(o => o.Frame == 1 && !o.IsReverseComplement).ToList();

        Assert.That(frame1Orfs, Has.Count.EqualTo(1),
            "Should find exactly 1 ORF — N in coding region should not break it");
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Sequence too short to contain valid ORF returns empty.
    /// </summary>
    [Test]
    public void FindOrfs_VeryShortSequence_ReturnsEmpty()
    {
        var orfs = GenomeAnnotator.FindOrfs("AT", minLength: 1).ToList();

        Assert.That(orfs, Is.Empty, "Sequence shorter than 3bp cannot contain ORF");
    }

    /// <summary>
    /// Sequence with only start codon (no room for stop) returns empty.
    /// </summary>
    [Test]
    public void FindOrfs_OnlyStartCodon_ReturnsEmpty()
    {
        var orfs = GenomeAnnotator.FindOrfs("ATG", minLength: 1, requireStartCodon: true).ToList();

        Assert.That(orfs, Is.Empty, "ATG alone cannot form complete ORF");
    }

    /// <summary>
    /// Null sequence parameter is handled.
    /// </summary>
    [Test]
    public void FindOrfs_NullSequence_ReturnsEmpty()
    {
        // Implementation should handle null gracefully
        var orfs = GenomeAnnotator.FindOrfs(null!, minLength: 1).ToList();

        Assert.That(orfs, Is.Empty);
    }

    #endregion
}
