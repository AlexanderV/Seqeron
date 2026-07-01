using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class EpigeneticAgeTests
{
    private static Dictionary<string, double> Methylation() =>
        new() { ["cg1"] = 0.5, ["cg2"] = 0.25 };

    private static Dictionary<string, double> Coefficients() =>
        new() { ["cg1"] = 1.0, ["cg2"] = 2.0 };

    [Test]
    public void EpigeneticAge_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.EpigeneticAge(Methylation(), Coefficients(), 0.0));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.EpigeneticAge(new Dictionary<string, double>(), Coefficients()));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.EpigeneticAge(Methylation(), new Dictionary<string, double>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.EpigeneticAge(null!, Coefficients()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.EpigeneticAge(Methylation(), null!));
    }

    [Test]
    public void EpigeneticAge_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.CalculateEpigeneticAge:
        //   linearPredictor = intercept + Σ coef_i·β_i = 0.5 + 1.0*0.5 + 2.0*0.25 = 1.5.
        //   x >= 0 -> antiTransform = (1+20)*x + 20 = 21*1.5 + 20 = 51.5.
        var age = AnnotationTools.EpigeneticAge(Methylation(), Coefficients(), intercept: 0.5);
        Assert.That(age.Age, Is.EqualTo(51.5).Within(1e-9));
    }

    [Test]
    public void EpigeneticAge_NegativePredictor_UsesExponentialBranch()
    {
        // linearPredictor = -1.0 + 1.0*0.0 = -1.0 (< 0).
        //   antiTransform = (1+20)*exp(-1) - 1 = 6.72546826...
        var meth = new Dictionary<string, double> { ["cg1"] = 0.0 };
        var coef = new Dictionary<string, double> { ["cg1"] = 1.0 };
        var age = AnnotationTools.EpigeneticAge(meth, coef, intercept: -1.0);
        Assert.That(age.Age, Is.EqualTo(6.7254682646002895).Within(1e-9));
    }

    [Test]
    public void EpigeneticAge_UnmatchedCpG_IsIgnored()
    {
        // A CpG present in methylation but absent from the coefficient table contributes nothing.
        var meth = new Dictionary<string, double> { ["cg1"] = 0.5, ["unknown"] = 1.0 };
        var coef = new Dictionary<string, double> { ["cg1"] = 1.0 };
        // linearPredictor = 0.0 + 1.0*0.5 = 0.5 -> 21*0.5 + 20 = 30.5.
        var age = AnnotationTools.EpigeneticAge(meth, coef);
        Assert.That(age.Age, Is.EqualTo(30.5).Within(1e-9));
    }
}
