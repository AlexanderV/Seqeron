using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindDiscordantPairsTests
{
    // Defaults: expected 400, sd 50, cutoffSd 3.0 -> concordant span window [250, 550].
    private static List<ReadPairInputDto> Pairs() =>
        new()
        {
            // Concordant: same chr, FR (+/-), span 400 in window -> NOT discordant.
            new ReadPairInputDto("concordant", "chr1", 100, '+', "chr1", 500, '-', 400),
            // Interchromosomal -> discordant (translocation signature).
            new ReadPairInputDto("interchrom", "chr1", 100, '+', "chr2", 500, '-', 400),
            // FR but span 800 > 550 -> discordant (deletion signature).
            new ReadPairInputDto("largespan",  "chr1", 100, '+', "chr1", 900, '-', 800),
            // Same-strand FF -> discordant (inversion signature).
            new ReadPairInputDto("inversion",   "chr1", 100, '+', "chr1", 500, '+', 400),
        };

    [Test]
    public void FindDiscordantPairs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindDiscordantPairs(Pairs()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindDiscordantPairs(null!));
    }

    [Test]
    public void FindDiscordantPairs_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.FindDiscordantPairs flags interchromosomal, out-of-span, and
        // abnormal-orientation pairs. The concordant FR pair with an in-window span is excluded.
        var result = AnnotationTools.FindDiscordantPairs(Pairs());

        Assert.That(result.Pairs.Select(p => p.ReadId),
            Is.EquivalentTo(new[] { "interchrom", "largespan", "inversion" }));
        Assert.That(result.Pairs.All(p => p.IsDiscordant), Is.True);
    }

    [Test]
    public void FindDiscordantPairs_AllConcordant_ReturnsEmpty()
    {
        var pairs = new List<ReadPairInputDto>
        {
            new ReadPairInputDto("c1", "chr1", 100, '+', "chr1", 500, '-', 400),
            new ReadPairInputDto("c2", "chr2", 100, '+', "chr2", 480, '-', 380),
        };
        var result = AnnotationTools.FindDiscordantPairs(pairs);
        Assert.That(result.Pairs, Is.Empty);
    }
}
