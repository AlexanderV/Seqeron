using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Annotation area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION (and its documentation), not from the
/// implementation's observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-ORF-001 — open reading frame detection (Annotation).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 28.
///
/// API under test (GenomeAnnotator.FindOrfs):
///   FindOrfs(dna, minLength = 100, searchBothStrands = true, requireStartCodon = true)
///     • Forward strand: in each of the 3 frames, scan codons; on a start codon
///       (ATG/GTG/TTG) push a pending start; on a stop codon (TAA/TAG/TGA) emit every
///       pending start s as the ORF [s, t+3) with aaLength = (t − s)/3, kept iff
///       aaLength ≥ minLength, then clear the pending starts.
///     • Reverse strand (searchBothStrands): the SAME scan runs on the reverse complement,
///       then each ORF is remapped to forward coordinates via
///       Start = |dna| − End_rev,  End = |dna| − Start_rev.
///   An ORF's identity is (Start, End, Frame, IsReverseComplement, Sequence, ProteinSequence);
///   Start/End are original-sequence coordinates (0-based, End exclusive).
///
/// Relations (derived from the definition / ORF_Detection.md §2.2–§2.4, NOT from output):
///   • MON   — minLength is ONLY a filter (aaLength ≥ minLength). Lowering it keeps every
///             ORF that already passed and may admit more, so the ORF set is a SUPERSET and
///             the count is non-decreasing along a decreasing-minLength chain.
///   • SHIFT — Start/End are forward-strand coordinates for both strands. Prepending an
///             IN-FRAME, non-coding flank F with |F| ≡ 0 (mod 3) preserves every reading
///             frame and creates no start/stop, so every ORF's Start and End advance by
///             exactly |F| with Frame, strand, Sequence and ProteinSequence preserved.
///             (Poly-C is start/stop-free; its reverse complement poly-G is too, so neither
///             strand gains an ORF; |F| ≡ 0 mod 3 keeps the frame labels.)
///   • INV   — a forward-strand ORF lying ENTIRELY upstream of an insertion point depends
///             only on codons before that point, so inserting ANY region at/after its End
///             leaves it byte-for-byte unchanged (same coordinates, sequence, protein).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class AnnotationMetamorphicTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed so random inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260619);

    /// <summary>Generates a random DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>Run-order-independent identity of an ORF.</summary>
    private static (int Start, int End, int Frame, bool Rc, string Seq, string Prot) OrfId(GenomeAnnotator.OpenReadingFrame o)
        => (o.Start, o.End, o.Frame, o.IsReverseComplement, o.Sequence, o.ProteinSequence);

    private static HashSet<(int, int, int, bool, string, string)> OrfSet(
        string dna, int minLength, bool bothStrands = true)
        => GenomeAnnotator.FindOrfs(dna, minLength, bothStrands).Select(OrfId).ToHashSet();

    /// <summary>
    /// Bodies carrying explicit forward-frame-0 ORFs (start codon … in-frame codons … stop),
    /// plus a fixed-seed random body. Codons avoid in-frame stops so the ORFs are unambiguous.
    /// </summary>
    private static IEnumerable<string> OrfBodies()
    {
        // ATG AAA CCC AAA TAA (aaLen 4) · spacer · ATG CCC AAA CCC TGA (aaLen 4)
        yield return "ATGAAACCCAAATAA" + "GGG" + "ATGCCCAAACCCTGA";
        // A single longer ORF: ATG + 8×AAA + TAA (aaLen 9)
        yield return "ATG" + string.Concat(Enumerable.Repeat("AAA", 8)) + "TAA";
        // Two ORFs sharing structure with extra GTG/TTG starts (alternate start codons)
        yield return "GTGAAATTTAAACCCTAA" + "TTT" + "TTGCCCAAATAA";
        // Fixed-seed random body (relations must hold for arbitrary input too)
        yield return RandomDna(120);
    }

    /// <summary>A small minimum length (in aa) so the short test ORFs qualify.</summary>
    private const int MinAa = 3;

    #endregion

    #region MON — lowering minLength yields a superset of ORFs (count non-decreasing)

    [Test]
    [Description("MON: along a DECREASING minLength chain the ORF set grows monotonically — each higher-threshold set is a subset of every lower-threshold one, count non-decreasing.")]
    public void FindOrfs_LoweringMinLength_YieldsSuperset_CountNonDecreasing()
    {
        int[] decreasingMinLength = { 10, 6, 4, 3, 2, 1 };

        foreach (var body in OrfBodies())
        {
            HashSet<(int, int, int, bool, string, string)>? previous = null;
            int previousCount = -1;

            foreach (int minLen in decreasingMinLength)
            {
                var orfs = OrfSet(body, minLen);

                if (previous is not null)
                {
                    orfs.IsSupersetOf(previous).Should().BeTrue(
                        because: $"minLength is only the filter aaLength ≥ minLength, so lowering it to {minLen} keeps every ORF that already passed and may add more");
                    orfs.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                        because: $"a lower minLength ({minLen}) admits-or-keeps each ORF, never removes one — the count is non-decreasing");
                }

                previous = orfs;
                previousCount = orfs.Count;
            }
        }
    }

    [Test]
    [Description("MON: every ORF found under a stricter minLength is also found under a looser one (membership preserved when the threshold drops).")]
    public void FindOrfs_StrictThreshold_SubsetOfLooseThreshold()
    {
        foreach (var body in OrfBodies())
        {
            var strict = OrfSet(body, minLength: 8);
            var loose = OrfSet(body, minLength: 2);

            loose.IsSupersetOf(strict).Should().BeTrue(
                because: "an ORF whose aaLength ≥ 8 also satisfies aaLength ≥ 2, and no other rule changed, so the strict ORF set is a subset of the loose one");
        }
    }

    #endregion

    #region SHIFT — prepending an in-frame non-coding flank shifts every ORF by |F|

    [Test]
    [Description("SHIFT: prepending a poly-C flank of length ≡ 0 (mod 3) advances every ORF's Start and End by exactly |F|, preserving Frame, strand, Sequence and ProteinSequence.")]
    public void FindOrfs_PrependInFrameNonCodingFlank_ShiftsAllOrfsByFlankLength()
    {
        foreach (var body in OrfBodies())
        {
            var baseOrfs = GenomeAnnotator.FindOrfs(body, MinAa, searchBothStrands: true).ToList();
            baseOrfs.Should().NotBeEmpty(because: "each body embeds at least one ORF at the small minimum length");

            foreach (int flankLen in new[] { 3, 6, 12 })
            {
                string flank = new string('C', flankLen); // poly-C: no start/stop; revComp poly-G: no start/stop
                var shifted = OrfSet(flank + body, MinAa);

                var expected = baseOrfs
                    .Select(o => (o.Start + flankLen, o.End + flankLen, o.Frame, o.IsReverseComplement, o.Sequence, o.ProteinSequence))
                    .ToHashSet();

                shifted.SetEquals(expected).Should().BeTrue(
                    because: $"an in-frame (|F| = {flankLen} ≡ 0 mod 3) non-coding prefix advances every ORF by exactly {flankLen} on both strands " +
                             "while preserving frame, strand, sequence and protein — and creates no new ORF (poly-C / poly-G carry no start or stop codon)");
            }
        }
    }

    [Test]
    [Description("SHIFT anchor: a known forward ORF starting at index 0 moves to exactly |F| when an in-frame poly-C flank is prepended.")]
    public void FindOrfs_PrependFlank_KnownForwardOrf_ShiftsExactly()
    {
        const string body = "ATG" + "AAACCCAAACCC" + "TAA"; // ATG at 0; stop TAA at 15; End 18
        var baseOrf = GenomeAnnotator.FindOrfs(body, MinAa, searchBothStrands: false)
            .Single(o => o.Frame == 1 && o.Start == 0);
        baseOrf.End.Should().Be(18, because: "the ORF spans ATG…TAA over 18 nucleotides");

        const int flankLen = 9;
        var shiftedOrf = GenomeAnnotator.FindOrfs(new string('C', flankLen) + body, MinAa, searchBothStrands: false)
            .Single(o => o.Frame == 1 && o.IsReverseComplement == false && o.Sequence == baseOrf.Sequence);

        shiftedOrf.Start.Should().Be(0 + flankLen, because: "the in-frame poly-C prefix shifts the ORF start by exactly the flank length");
        shiftedOrf.End.Should().Be(18 + flankLen, because: "the ORF end shifts by the same flank length");
    }

    #endregion

    #region INV — inserting a region downstream does not change ORFs entirely upstream of it

    [Test]
    [Description("INV: inserting ANY region at a codon boundary leaves every forward-strand ORF that lies entirely upstream of the insertion point byte-for-byte unchanged.")]
    public void FindOrfs_InsertDownstream_UpstreamForwardOrfsUnchanged()
    {
        // prefix ends right after a complete forward ORF (ATG…TAA); insertion goes at |prefix|.
        const string prefix = "ATGAAACCCAAATAA";           // forward ORF [0,15)
        const string suffix = "GGGATGCCCAAACCCTGA";          // a downstream ORF (will move)
        int insertPos = prefix.Length;

        // Forward-only search isolates the "upstream" direction on the coding strand.
        var baseUpstream = GenomeAnnotator.FindOrfs(prefix + suffix, MinAa, searchBothStrands: false)
            .Where(o => o.End <= insertPos).Select(OrfId).ToHashSet();
        baseUpstream.Should().NotBeEmpty(because: "the prefix contains a complete forward ORF upstream of the insertion point");

        foreach (var insert in new[] { "CCC", "CCCCCCCCC", RandomDna(9), RandomDna(12), "ATGTTTGGGTAA" })
        {
            var changedUpstream = GenomeAnnotator.FindOrfs(prefix + insert + suffix, MinAa, searchBothStrands: false)
                .Where(o => o.End <= insertPos).Select(OrfId).ToHashSet();

            changedUpstream.SetEquals(baseUpstream).Should().BeTrue(
                because: "an ORF ending at or before the insertion point reads only codons before it, which the downstream insertion never touches — " +
                         "so its coordinates, sequence and protein are preserved exactly regardless of what is inserted");
        }
    }

    #endregion
}
