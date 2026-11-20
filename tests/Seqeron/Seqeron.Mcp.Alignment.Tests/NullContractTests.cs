using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// Audit A5-3 (B09 F169): the MCP tools over the secondary entry points refuse a null element / a nameless sequence / a matrix beyond the char-indexed table with a typed
/// <see cref="ArgumentException"/> (the shared message style of the library), not a <see cref="NullReferenceException"/> / <see cref="IndexOutOfRangeException"/>.
/// </summary>
[TestFixture]
public class NullContractTests
{
    [Test]
    public void CommandTools_NullArgumentElement_IsAnArgumentNullException()
    {
        string[] args = ["-i", null!];
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => AlignmentTools.KalignCommand(args))!.ParamName, Is.EqualTo("arguments"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => AlignmentTools.PoaCommand(args))!.ParamName, Is.EqualTo("arguments"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => AlignmentTools.ProbConsCommand(args))!.ParamName, Is.EqualTo("arguments"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.KalignCommand(null!));
            Assert.Throws<ArgumentException>(() => AlignmentTools.PoaCommand(null!));
            Assert.Throws<ArgumentException>(() => AlignmentTools.ProbConsCommand(null!));
        });
    }

    [Test]
    public void MultipleAlignPoa_NullName_AndWideMatrix_AreArgumentExceptions()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignPoa(["ACGT", "ACGA"], names: ["a", null!]))!.Message, Does.Contain("Sequence 1 has no name"));
            Assert.That(Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignPoa(["ACGT", "ACGA"], matrix: new string('A', 129) + "\n"))!.Message, Does.Contain("at most 128"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignPoa(["ACGT", "ACGA"], matrix: "Aé\n"));
        });
    }

    [Test]
    public void TCoffeeStageTools_NullElementAndForeignCharacters_AreArgumentExceptions()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => AlignmentTools.TCoffeeKtupSimilarity(["ACGTACGT", null!]))!.ParamName, Is.EqualTo("sequences"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeeLalignTopAlignments("ACGT*ACGT", "ACGTACGT"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeeLalignTopAlignments("ACGTACGT", "ACGéACGT"));
        });
    }

    [Test]
    public void MultipleAlignTools_NullSequenceList_AndNullStructureOutput_AreArgumentExceptions()
    {
        // found by the reflection fuzzer of audit A5-3: ParseDna(null) and a null element of structureAlignerOutputs were NullReferenceExceptions
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignClustalW(null!))!.ParamName, Is.EqualTo("sequences"));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignTCoffee(["ACGTACGT", "ACGTACGA"], structureAlignerOutputs: [null!]));
        });
    }
}
