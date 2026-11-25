using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// B09 A2-30 .. A2-36: the MCP parameters / tools added with the both-strand repeats, canonical tandem rotation, fuzznuc pattern language, common-region modes and strands,
/// k-mismatch common region and IUPAC input. Expected values are those locked by the reference-checked library tests (MUMmer 3.23 repeat-match / mummer, EMBOSS fuzznuc 6.6.0)
/// or derived by hand from the definitions.
/// </summary>
[TestFixture]
public class StrandPatternAndRepeatToolsTests
{
    [Test]
    public void FindRepeats_BothStrands_ReturnsClassesWithReversePositions()
    {
        // ACGGTTC ... GAACCGT: the second copy is the reverse complement (library test ReverseComplementCopy_IsFoundWithBothPositionLists)
        var r = AnalysisTools.FindRepeats("ACGGTTCAAAAGAACCGT", 7, maximalOnly: true, bothStrands: true).Items.Single();
        Assert.That((r.Sequence, r.Positions, r.ReversePositions, r.Count), Is.EqualTo(("ACGGTTC", new[] { 0 }, new[] { 11 }, 2)));
        // a palindromic site is a repeat on its own (repeat-match prints the pair i = k)
        var g = AnalysisTools.FindRepeats("ACGAATTCCA", 6, maximalOnly: true, bothStrands: true).Items.Single(i => i.Sequence == "GAATTC");
        Assert.That((g.Positions, g.ReversePositions), Is.EqualTo((new[] { 2 }, new[] { 2 })));
        // forward-only items carry no reversePositions
        Assert.That(AnalysisTools.FindRepeats("AAAA", 2).Items.All(i => i.ReversePositions is null), Is.True);
    }

    [Test]
    public void FindRepeats_IupacInput_NeverRepeatsAcrossN()
    {
        Assert.That(AnalysisTools.FindRepeats("NNNNNNNN", 1).Items, Is.Empty);
        Assert.That(AnalysisTools.FindRepeats("ACGTNACGTN", 1, maximalOnly: true).Items.Select(i => i.Sequence), Is.EqualTo(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRepeats("ACGU", 2));
    }

    [Test]
    public void FindAllLongestRepeats_ListsEveryTie()
    {
        // ATGATG: the longest repeat ATG at 0 and 3
        var one = AnalysisTools.FindAllLongestRepeats("ATGATG");
        Assert.That((one.Length, one.Items.Single().Sequence, one.Items.Single().Positions), Is.EqualTo((3, "ATG", new[] { 0, 3 })));
        // ACGTTGCA ... both ACG and TTG...: two ties of length 3: ACGxxxACG + GGGxxxGGG -> every longest, by first occurrence
        var ties = AnalysisTools.FindAllLongestRepeats("ACGCCCTGATGGGACGTGATGT");
        Assert.That(ties.Items.Select(i => i.Sequence), Is.EqualTo(ties.Items.Select(i => i.Sequence).Distinct()));
        Assert.That(ties.Items.All(i => i.Length == ties.Length), Is.True);
        Assert.That(ties.Items.Select(i => i.Positions[0]), Is.Ordered);
        // nothing repeats
        Assert.That(AnalysisTools.FindAllLongestRepeats("ACGT").Items, Is.Empty);
        Assert.That(AnalysisTools.FindAllLongestRepeats("ACGT").Length, Is.EqualTo(0));
        // both strands: the longest class is the reverse-complement pair
        var both = AnalysisTools.FindAllLongestRepeats("ACGGTTCAAAAGAACCGT", bothStrands: true);
        Assert.That((both.Length, both.Items.Single().Sequence, both.Items.Single().ReversePositions), Is.EqualTo((7, "ACGGTTC", new[] { 11 })));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindAllLongestRepeats(""));
    }

    [Test]
    public void FindTandemRepeats_CanonicalRotation_ReportsTheLeastRotation()
    {
        // TGATGATGAT: TGA x3 at 0; the least rotation of TGA is ATG, the original spelling is kept in originalUnit
        var plain = AnalysisTools.FindTandemRepeats("TGATGATGAT", 3, 2).Items.Single();
        Assert.That((plain.Unit, plain.OriginalUnit), Is.EqualTo(("TGA", (string?)null)));
        var canonical = AnalysisTools.FindTandemRepeats("TGATGATGAT", 3, 2, canonicalRotation: true).Items.Single();
        Assert.That((canonical.Unit, canonical.OriginalUnit, canonical.Position, canonical.Repetitions, canonical.TotalLength), Is.EqualTo(("ATG", "TGA", 0, 3, 9)));
    }

    [Test]
    public void FindCommonRegions_ModesAndStrands()
    {
        // the library test (MUMmer mummer -maxmatch / -mum -b): MEMs of the reference ACGGATCCGTTTTACGATCCGTAAGT against GGATCCGAAAACGGATCCGTTTT
        const string reference = "ACGGATCCGTTTTACGATCCGTAAGT", query = "GGATCCGAAAACGGATCCGTTTT";
        var mems = AnalysisTools.FindCommonRegions(reference, query, 4, "MaximalExactMatches", "Both").Items;
        // mummer -maxmatch -b -l 4 (1-based pos1 pos2 len; reverse positions in the reverse-complemented query: 1 4 13 / 16 7 7 / 2 17 7 / 16 19 5; with -c 20 / 17 / 7 / 5 = forward position of the last base)
        Assert.That(mems.Where(m => !m.IsReverseComplement).Select(m => (m.PositionInFirst, m.PositionInSecond, m.Length)).OrderBy(m => m),
            Is.EqualTo(new[] { (2, 0, 7), (15, 1, 6), (0, 10, 13), (15, 13, 7) }.OrderBy(m => m)));
        Assert.That(mems.Where(m => m.IsReverseComplement).Select(m => (m.PositionInFirst, m.PositionInSecond, m.Length)).OrderBy(m => m),
            Is.EqualTo(new[] { (0, 7, 13), (15, 10, 7), (1, 0, 7), (15, 0, 5) }.OrderBy(m => m)));
        Assert.That(AnalysisTools.FindCommonRegions(reference, query, 4).Items.All(m => !m.IsReverseComplement), Is.True);
        Assert.That(AnalysisTools.FindCommonRegions("ACGTACGTNNNNGGATCC", "TTACGTACGTAACCGGATCCTT", 5, "MaximalUnique").Items.Select(m => m.Sequence),
            Is.EquivalentTo(new[] { "GGATCC", "ACGTACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindCommonRegions(reference, query, 4, "nope"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindCommonRegions(reference, query, 4, "MaximalUnique", "sideways"));
    }

    [Test]
    public void FindLongestCommonRegionMismatches_KMismatchLcs()
    {
        // brute-force oracle row (common_region_kmismatch_oracle.tsv): CACAAAAC / AACCCCAACAGCAC, k = 1 -> 5:0:3:1 (length 5, starts 0 and 3, 1 mismatch)
        var r = AnalysisTools.FindLongestCommonRegionMismatches("CACAAAAC", "AACCCCAACAGCAC", 1);
        Assert.That((r.Length, r.PositionInFirst, r.PositionInSecond, r.Mismatches), Is.EqualTo((5, 0, 3, 1)));
        Assert.That(r.Sequence1, Is.EqualTo("CACAA").And.Length.EqualTo(r.Sequence2.Length));
        Assert.That(AnalysisTools.FindLongestCommonRegionMismatches("CACAAAAC", "AACCCCAACAGCAC", 0).Length, Is.EqualTo(3));
        Assert.That(AnalysisTools.FindLongestCommonRegionMismatches("TTTACGGTTCATTT", "GGGTGAACCGTAGG", 0, "ReverseComplement").IsReverseComplement, Is.True);
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLongestCommonRegionMismatches("ACGT", "ACGT", -1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongestCommonRegionMismatches("", "ACGT", 0));
    }
}
