using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// Every member of <see cref="ISuffixTree"/> on every string of the pre-registered small sizes (sigma 2 up to length 10, sigma 3 up to 7, sigma 4 up to 5) through the facade, against the
/// definitions and the Ukkonen tree; the other two child-table layouts on smaller sets. Heavy; explicit, category <c>Exhaustive</c>.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.6)")]
public class FacadeExhaustiveTests
{
    private static void All(int sigma, int maxLength, SuffixTreeBackend backend)
    {
        ulong seed = 1;
        foreach (string text in ExhaustiveStrings.All(sigma, maxLength))
        {
            ISuffixTree tree = backend.FromString(text);
            try
            {
                FacadeChecks.CheckEverything(tree, UkkonenTree.Build(text), null, text, seed++, $"[{backend.Name}] sigma {sigma}", maxPatternLength: 5);
            }
            finally
            {
                (tree as IDisposable)?.Dispose();
            }
        }
    }

    [Test]
    public void AllBinaryStrings_UpTo10() => All(2, 10, SuffixTreeBackends.Esa);

    [Test]
    public void AllTernaryStrings_UpTo7() => All(3, 7, SuffixTreeBackends.Esa);

    [Test]
    public void AllFourLetterStrings_UpTo5() => All(4, 5, SuffixTreeBackends.Esa);

    [Test]
    public void AllBinaryStrings_UpTo8_PackedLayout() => All(2, 8, SuffixTreeBackends.EsaPacked);

    [Test]
    public void AllBinaryStrings_UpTo8_ReferenceLayout() => All(2, 8, SuffixTreeBackends.EsaReference3);

    [Test]
    public void AllTernaryStrings_UpTo5_PackedAndReferenceLayouts()
    {
        All(3, 5, SuffixTreeBackends.EsaPacked);
        All(3, 5, SuffixTreeBackends.EsaReference3);
    }

    [Test]
    public void AllBinaryStrings_UpTo9_LogicalHashAndPrintTreeAgainstUkkonen()
    {
        foreach (string text in ExhaustiveStrings.All(2, 9))
        {
            ISuffixTree tree = SuffixTreeBackends.Esa.FromString(text);
            UkkonenTree reference = UkkonenTree.Build(text);
            Assert.That(SuffixTreeSerializer.CalculateLogicalHash(tree), Is.EqualTo(SuffixTreeSerializer.CalculateLogicalHash(reference)), text);
            FacadeChecks.CheckPrintTree(tree, reference, text, "exhaustive");
        }
    }
}
