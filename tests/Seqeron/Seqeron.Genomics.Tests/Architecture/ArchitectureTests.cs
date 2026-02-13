using ArchUnitNET.Domain;
using ArchUnitNET.Loader;
using ArchUnitNET.NUnit;
using static ArchUnitNET.Fluent.ArchRuleDefinition;

namespace Seqeron.Genomics.Tests.Architecture;

/// <summary>
/// Architecture guard tests using ArchUnitNET.
/// These tests enforce structural rules across the codebase and prevent
/// architectural drift (e.g., Core depending on Analysis, circular dependencies).
///
/// Rules are checked at the IL level — they catch violations even across project boundaries.
/// </summary>
[TestFixture]
[Category("Architecture")]
public class ArchitectureTests
{
    private static readonly ArchUnitNET.Domain.Architecture Arch = new ArchLoader()
        .LoadAssemblies(
            typeof(DnaSequence).Assembly,             // Seqeron.Genomics.Core
            typeof(DisorderPredictor).Assembly,       // Seqeron.Genomics.Analysis
            typeof(SequenceAligner).Assembly,         // Seqeron.Genomics.Alignment
            typeof(FastaParser).Assembly)             // Seqeron.Genomics.IO
        .Build();

    /// <summary>
    /// Core must not depend on Analysis — Core is the innermost layer.
    /// </summary>
    [Test]
    public void Core_ShouldNotDependOn_Analysis()
    {
        Types().That().ResideInNamespace("Seqeron.Genomics.Core")
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Analysis"))
            .Check(Arch);
    }

    /// <summary>
    /// Core must not depend on IO — Core should be pure domain logic.
    /// </summary>
    [Test]
    public void Core_ShouldNotDependOn_IO()
    {
        Types().That().ResideInNamespace("Seqeron.Genomics.Core")
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.IO"))
            .Check(Arch);
    }

    /// <summary>
    /// Core must not depend on Alignment — Alignment is a higher-level module.
    /// </summary>
    [Test]
    public void Core_ShouldNotDependOn_Alignment()
    {
        Types().That().ResideInNamespace("Seqeron.Genomics.Core")
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Alignment"))
            .Check(Arch);
    }

    /// <summary>
    /// IO must not depend on Analysis — IO handles parsing, not algorithm logic.
    /// </summary>
    [Test]
    public void IO_ShouldNotDependOn_Analysis()
    {
        Types().That().ResideInNamespace("Seqeron.Genomics.IO")
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Analysis"))
            .Check(Arch);
    }

    /// <summary>
    /// Static analyzer/predictor/finder classes should remain static (abstract+sealed in IL).
    /// This prevents accidental introduction of instance state.
    /// </summary>
    [Test]
    public void AnalyzerClasses_ShouldBeStatic()
    {
        Classes().That().HaveNameEndingWith("Predictor")
            .Or().HaveNameEndingWith("Finder").And().ResideInNamespace("Seqeron.Genomics.Analysis")
            .Should().BeAbstract()
            .Check(Arch);
    }
}
