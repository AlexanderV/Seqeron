using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// Every operation of <see cref="EsaSuffixIndex"/> on every locus of every string of the pre-registered small sizes (sigma 2 up to length 10, sigma 3 up to 7, sigma 4 up to 5), against the definitional
/// oracle (<see cref="EsaContractChecker"/>), with the default child layout; the two other layouts on binary strings up to 8. Heavy; explicit, category <c>Exhaustive</c>.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.5)")]
public class EsaSuffixIndexExhaustiveTests
{
    private static void All(int sigma, int maxLength, ChildLayoutKind layout)
    {
        foreach (string text in ExhaustiveStrings.All(sigma, maxLength))
        {
            var (index, _) = EsaSuffixIndexTests.Make(text, layout);
            EsaContractChecker.CheckAll(index, new IndexContractOracle(text), $"[{layout}] sigma {sigma} '{text}'", CheckOptions.Thorough);
        }
    }

    [Test]
    public void AllBinaryStrings_UpTo10() => All(2, 10, ChildLayoutKind.Tagged);

    [Test]
    public void AllTernaryStrings_UpTo7() => All(3, 7, ChildLayoutKind.Tagged);

    [Test]
    public void AllFourLetterStrings_UpTo5() => All(4, 5, ChildLayoutKind.Tagged);

    [Test]
    public void AllBinaryStrings_UpTo8_PackedLayout() => All(2, 8, ChildLayoutKind.Packed);

    [Test]
    public void AllBinaryStrings_UpTo8_ReferenceLayout() => All(2, 8, ChildLayoutKind.Reference3);

    [Test]
    public void AllTernaryStrings_UpTo6_PackedAndReferenceLayouts()
    {
        All(3, 6, ChildLayoutKind.Packed);
        All(3, 6, ChildLayoutKind.Reference3);
    }
}
