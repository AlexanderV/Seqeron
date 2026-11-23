using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the audit of work package 1b.9c showed in <see cref="DifferenceCover"/> (the mutants D14 and D18 of <c>analysis/mutation/1b7d_dcs_differencecover_listed.json</c>): the pruned covers are the covers of the prototype
/// (E2), so that the rank files are, and the number of rounds of the pruning and the choice of the best round decide which cover that is (a seventh round for the periods above 1,024 finds a smaller cover for 1,169, 1,282 and 1,288). The figures were computed by a second implementation of the construction (the same
/// SplitMix64 stream and the same greedy removal, written again outside the production code), not read off the code under test.
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SaCoverGap1b9cTests
{
    [TestCase(10, 4, 21)]
    [TestCase(50, 9, 147)]
    [TestCase(100, 14, 225)]
    [TestCase(200, 20, 1404)]
    [TestCase(256, 23, 808)]
    [TestCase(257, 23, 1887)]
    [TestCase(300, 25, 2517)]
    [TestCase(323, 26, 2819)]
    [TestCase(341, 27, 1729)]
    [TestCase(380, 29, 3526)]
    [TestCase(443, 31, 4737)]
    [TestCase(672, 39, 8605)]
    [TestCase(754, 41, 6704)]
    [TestCase(1024, 47, 8176)]
    [TestCase(1025, 47, 16683)]
    [TestCase(1169, 52, 10605)]
    [TestCase(1282, 54, 23030)]
    [TestCase(1288, 54, 13410)]
    [TestCase(1500, 57, 29207)]
    [TestCase(2000, 66, 45020)]
    [TestCase(4096, 95, 69664)]
    public void ThePrunedCover_IsThePrototypesCover_OfTheSameSizeAndTheSameResidues(int period, int residues, long sum)
    {
        // D14 (the large periods get more rounds: a 13th round finds a smaller cover for 323, 341, 380, 443, 672 and 754) and D18 (the best round is the largest: the unpruned start is kept); the fresh mutant F54 (a seventh round above 1,024)
        DifferenceCover cover = DifferenceCover.Create(period);
        Assert.That(cover.Kind, Is.EqualTo(CoverKind.Pruned));
        Assert.That(cover.Verify(), Is.True);
        Assert.That(cover.Residues, Has.Length.EqualTo(residues), $"|D| of X = {period}");
        Assert.That(cover.Residues.Sum(r => (long)r), Is.EqualTo(sum), $"the residues of X = {period}");
    }

    [Test]
    public void ThePrunedCoverOfAHundred_IsTheRunOfTenAndFourMultiples()
    {
        Assert.That(DifferenceCover.Create(100).Residues, Is.EqualTo(new[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 60, 70 }));
        Assert.That(DifferenceCover.Create(10).Residues, Is.EqualTo(new[] { 3, 4, 6, 8 }));
    }
}
