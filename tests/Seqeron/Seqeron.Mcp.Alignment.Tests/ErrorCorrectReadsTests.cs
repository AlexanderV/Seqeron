using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class ErrorCorrectReadsTests
{
    [Test]
    public void ErrorCorrectReads_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.ErrorCorrectReads(new[] { "ACGT" }, 4, 2));
        // k-mer size below 1 is rejected by the underlying spectrum method.
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.ErrorCorrectReads(new[] { "ACGT" }, 0, 2));
        // Empty input returns empty output.
        Assert.That(AlignmentTools.ErrorCorrectReads(Array.Empty<string>(), 4, 2).Corrected, Is.Empty);
    }

    [Test]
    public void ErrorCorrectReads_Binding_InvokesSuccessfully()
    {
        // Three copies of ACGTACGT make its 4-mers trusted (freq >= 2); the fourth read
        // "ACGTTCGT" has a single error at the middle base that the k-mer spectrum method
        // corrects back to "ACGTACGT" (Musket/Quake two-sided substitution).
        var r = AlignmentTools.ErrorCorrectReads(
            new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT" }, 4, 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Corrected, Has.Length.EqualTo(4));
            Assert.That(r.Corrected[3], Is.EqualTo("ACGTACGT"));
            Assert.That(r.Corrected[0], Is.EqualTo("ACGTACGT"));
        });
    }

    [Test]
    public void ErrorCorrectReads_MusketOptions_DelegateToCore()
    {
        // Three reverse-strand copies of GATTACAGGCTTC + one forward read with C->G at index 6 (k=3,
        // cut-off 2): forward-only counting miscorrects index 0; Musket's canonical counting fixes index 6.
        var reads = new[] { "GAAGCCTGTAATC", "GAAGCCTGTAATC", "GAAGCCTGTAATC", "GATTACGGGCTTC" };
        var plain = AlignmentTools.ErrorCorrectReads(reads, 3, 2);
        var musket = AlignmentTools.ErrorCorrectReads(reads, 3, 2, musket: true);
        var canonical = AlignmentTools.ErrorCorrectReads(reads, 3, 2, countReverseComplement: true);
        // Two adjacent errors need Musket's aggressive stage with -maxerr >= 2.
        var twoErrors = new[] { "AAAAAAAA", "AAAAAAAA", "AAAAAAAA", "AACCAAAA" };
        Assert.Multiple(() =>
        {
            Assert.That(plain.Corrected[3], Is.EqualTo("AATTACGGGCTTC"));
            Assert.That(musket.Corrected[3], Is.EqualTo("GATTACAGGCTTC"));
            Assert.That(canonical.Corrected[3], Is.EqualTo("GATTACAGGCTTC"));
            Assert.That(AlignmentTools.ErrorCorrectReads(twoErrors, 4, 2, musket: true).Corrected[3], Is.EqualTo("AAAAAAAA"));
            Assert.That(AlignmentTools.ErrorCorrectReads(twoErrors, 4, 2, musket: true, maxErrorsPerKmer: 1).Corrected[3], Is.EqualTo("AACCAAAA"));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.ErrorCorrectReads(twoErrors, 4, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.ErrorCorrectReads(twoErrors, 4, 2, maxIterations: 0));
        });
    }
}
