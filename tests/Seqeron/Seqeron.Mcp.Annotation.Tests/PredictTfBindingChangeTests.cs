using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictTfBindingChangeTests
{
    // Reference context, length 41, with the AP-1 motif TGACTCA placed at index 20..26.
    // The variant sits at contextOffset 20 (default), i.e. the motif's first base (T),
    // and mutates it to A, breaking one of seven motif positions.
    private const string Motif = "TGACTCA";
    private static readonly string Context =
        new string('N', 20) + Motif + new string('N', 41 - 20 - Motif.Length);

    private static AnnotatorVariantDto SnvAt20() =>
        // T -> A at the motif start; type must be SNV for the algorithm to score.
        new("chr1", 100, "T", "A", "SNV", null, null);

    private static List<TfMotifInputDto> Ap1() =>
        new() { new TfMotifInputDto("AP-1", Motif, 0.9) };

    [Test]
    public void PredictTfBindingChange_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictTfBindingChange(SnvAt20(), Ap1(), Context));

        Assert.Throws<ArgumentNullException>(
            () => AnnotationTools.PredictTfBindingChange(null!, Ap1(), Context));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.PredictTfBindingChange(SnvAt20(), Ap1(), ""));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.PredictTfBindingChange(SnvAt20(), new List<TfMotifInputDto>(), Context));

        // Non-SNV variant is rejected by the wrapper guard.
        var deletion = new AnnotatorVariantDto("chr1", 100, "TA", "T", "Deletion", null, null);
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.PredictTfBindingChange(deletion, Ap1(), Context));
    }

    [Test]
    public void PredictTfBindingChange_Binding_InvokesSuccessfully()
    {
        // VariantAnnotator.ScoreMotif is the fraction of matching motif positions (best window).
        // Ref: perfect match at index 20 => 7/7 = 1.0.
        // Alt: T->A at motif position 0 breaks 1 of 7 => best 6/7 ≈ 0.857142857.
        // |diff| = 1/7 ≈ 0.142857 > 0.1, so the change is reported.
        var result = AnnotationTools.PredictTfBindingChange(SnvAt20(), Ap1(), Context);

        Assert.That(result.Changes, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Changes[0].TfName, Is.EqualTo("AP-1"));
            Assert.That(result.Changes[0].RefScore, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.Changes[0].AltScore, Is.EqualTo(6.0 / 7.0).Within(1e-9));
            Assert.That(result.Changes[0].ScoreDifference, Is.EqualTo(1.0 / 7.0).Within(1e-9));
        });
    }

    [Test]
    public void PredictTfBindingChange_SmallChange_IsFilteredOut()
    {
        // A silent (identical) alternate base leaves refScore == altScore, |diff| = 0 <= 0.1,
        // so the algorithm yields nothing.
        var silent = new AnnotatorVariantDto("chr1", 100, "T", "T", "SNV", null, null);
        var result = AnnotationTools.PredictTfBindingChange(silent, Ap1(), Context);
        Assert.That(result.Changes, Is.Empty);
    }
}
