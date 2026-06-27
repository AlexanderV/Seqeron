using System.Reflection;
using ArchUnitNET.Domain;
using ArchUnitNET.Fluent.Slices;
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
///
/// Coverage maps to docs/checklists/07_ARCHITECTURE_TESTING.md (22 module-dependency rules).
/// Rules apply to modules (projects), not to individual algorithms; the algorithm-to-module
/// mapping in the checklist documents which rule guards each of the 258 algorithms.
/// </summary>
[TestFixture]
[Category("Architecture")]
public class ArchitectureTests
{
    // All Seqeron.Genomics module assemblies — loaded by ArchUnitNET for dependency rules
    // and cycle detection (full graph, not a subset), and reused by the reflection-based
    // immutability rule.
    private static readonly System.Reflection.Assembly[] SeqeronAssemblies =
    [
        typeof(DnaSequence).Assembly,                  // Seqeron.Genomics.Core
        typeof(ScoringMatrix).Assembly,                // Seqeron.Genomics.Infrastructure
        typeof(DisorderPredictor).Assembly,            // Seqeron.Genomics.Analysis
        typeof(SequenceAligner).Assembly,              // Seqeron.Genomics.Alignment
        typeof(FastaParser).Assembly,                  // Seqeron.Genomics.IO
        typeof(VariantCaller).Assembly,                // Seqeron.Genomics.Annotation
        typeof(CrisprDesigner).Assembly,               // Seqeron.Genomics.MolTools
        typeof(PhylogeneticAnalyzer).Assembly,         // Seqeron.Genomics.Phylogenetics
        typeof(PopulationGeneticsAnalyzer).Assembly,   // Seqeron.Genomics.Population
        typeof(GenomeAssemblyAnalyzer).Assembly,       // Seqeron.Genomics.Chromosome
        typeof(PanGenomeAnalyzer).Assembly,            // Seqeron.Genomics.Metagenomics
        typeof(ImmuneAnalyzer).Assembly,               // Seqeron.Genomics.Oncology
        typeof(ReportGenerator).Assembly,              // Seqeron.Genomics.Reports
    ];

    private static readonly ArchUnitNET.Domain.Architecture Arch =
        new ArchLoader().LoadAssemblies(SeqeronAssemblies).Build();

    private const string Core = "Seqeron.Genomics.Core";
    private const string Io = "Seqeron.Genomics.IO";

    // ---------------------------------------------------------------------
    // Rules 1-5 — existing layer boundaries (Core is the innermost layer)
    // ---------------------------------------------------------------------

    /// <summary>Rule 1: Core must not depend on Analysis — Core is the innermost layer.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Analysis()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Analysis"))
            .Check(Arch);
    }

    /// <summary>Rule 2: Core must not depend on IO — Core should be pure domain logic.</summary>
    [Test]
    public void Core_ShouldNotDependOn_IO()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace(Io))
            .Check(Arch);
    }

    /// <summary>Rule 3: Core must not depend on Alignment — Alignment is a higher-level module.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Alignment()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Alignment"))
            .Check(Arch);
    }

    /// <summary>Rule 4: IO must not depend on Analysis — IO handles parsing, not algorithm logic.</summary>
    [Test]
    public void IO_ShouldNotDependOn_Analysis()
    {
        Types().That().ResideInNamespace(Io)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Analysis"))
            .Check(Arch);
    }

    /// <summary>
    /// Rule 5: Static analyzer/predictor/finder classes should remain static
    /// (abstract+sealed in IL). This prevents accidental introduction of instance state.
    /// </summary>
    [Test]
    public void AnalyzerClasses_ShouldBeStatic()
    {
        Classes().That().HaveNameEndingWith("Predictor")
            .Or().HaveNameEndingWith("Finder").And().ResideInNamespace("Seqeron.Genomics.Analysis")
            .Should().BeAbstract()
            .Check(Arch);
    }

    // ---------------------------------------------------------------------
    // Rules 6-13 — Core must not depend on any higher-level module
    // ---------------------------------------------------------------------

    /// <summary>Rule 6: Core must not depend on Annotation.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Annotation()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Annotation"))
            .Check(Arch);
    }

    /// <summary>Rule 7: Core must not depend on MolTools.</summary>
    [Test]
    public void Core_ShouldNotDependOn_MolTools()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.MolTools"))
            .Check(Arch);
    }

    /// <summary>Rule 8: Core must not depend on Phylogenetics.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Phylogenetics()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Phylogenetics"))
            .Check(Arch);
    }

    /// <summary>Rule 9: Core must not depend on Population.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Population()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Population"))
            .Check(Arch);
    }

    /// <summary>Rule 10: Core must not depend on Chromosome.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Chromosome()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Chromosome"))
            .Check(Arch);
    }

    /// <summary>Rule 11: Core must not depend on Metagenomics.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Metagenomics()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Metagenomics"))
            .Check(Arch);
    }

    /// <summary>Rule 12: Core must not depend on Oncology.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Oncology()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Oncology"))
            .Check(Arch);
    }

    /// <summary>Rule 13: Core must not depend on Reports.</summary>
    [Test]
    public void Core_ShouldNotDependOn_Reports()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Reports"))
            .Check(Arch);
    }

    // ---------------------------------------------------------------------
    // Rules 14-16 — IO must not depend on higher-level analysis modules
    // ---------------------------------------------------------------------

    /// <summary>Rule 14: IO must not depend on Alignment.</summary>
    [Test]
    public void IO_ShouldNotDependOn_Alignment()
    {
        Types().That().ResideInNamespace(Io)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Alignment"))
            .Check(Arch);
    }

    /// <summary>Rule 15: IO must not depend on Phylogenetics.</summary>
    [Test]
    public void IO_ShouldNotDependOn_Phylogenetics()
    {
        Types().That().ResideInNamespace(Io)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Phylogenetics"))
            .Check(Arch);
    }

    /// <summary>Rule 16: IO must not depend on Oncology.</summary>
    [Test]
    public void IO_ShouldNotDependOn_Oncology()
    {
        Types().That().ResideInNamespace(Io)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespace("Seqeron.Genomics.Oncology"))
            .Check(Arch);
    }

    // ---------------------------------------------------------------------
    // Rules 17-19 — structural invariants
    // ---------------------------------------------------------------------

    /// <summary>
    /// Rule 17: No circular dependencies between modules. Each top-level Seqeron.Genomics.*
    /// module is treated as a slice; the dependency graph between slices must be acyclic.
    /// </summary>
    [Test]
    public void Modules_ShouldBeFreeOfCycles()
    {
        SliceRuleDefinition.Slices().Matching("Seqeron.Genomics.(*)")
            .Should().BeFreeOfCycles()
            .Check(Arch);
    }

    /// <summary>
    /// Rule 18: Core must not use System.IO directly — file/stream access belongs in the IO
    /// module. Core stays a pure, side-effect-free domain layer.
    /// </summary>
    [Test]
    public void Core_ShouldNotDependOn_SystemIO()
    {
        Types().That().ResideInNamespace(Core)
            .Should().NotDependOnAny(
                Types().That().ResideInNamespaceMatching(@"^System\.IO"))
            .Check(Arch);
    }

    /// <summary>
    /// Rule 19: Result/DTO types are immutable — they expose data via getters or init-only
    /// setters (records), never settable public setters. This keeps results safe to share and
    /// prevents callers from mutating analysis output.
    /// </summary>
    /// <remarks>
    /// Implemented with reflection rather than ArchUnitNET: in v0.13.2 the
    /// <c>NotHavePublicSetter</c> predicate flags record init-only setters as public setters,
    /// which would wrongly reject the immutable record DTOs this rule is meant to permit.
    /// Reflection distinguishes <c>init</c> (allowed) from a writable <c>set</c> (forbidden)
    /// via the <c>IsExternalInit</c> required modifier the compiler emits on init accessors.
    /// </remarks>
    [Test]
    public void ResultAndDtoTypes_ShouldBeImmutable()
    {
        const string isExternalInit = "System.Runtime.CompilerServices.IsExternalInit";

        var offenders = new List<string>();

        foreach (var type in SeqeronAssemblies
                     .SelectMany(a => a.GetExportedTypes())
                     .Where(t => (t.Name.EndsWith("Result") || t.Name.EndsWith("Dto"))
                                 && (t.IsClass || (t.IsValueType && !t.IsEnum))))
        {
            foreach (var property in type.GetProperties(BindingFlags.Public | BindingFlags.Instance))
            {
                var setter = property.SetMethod;
                if (setter is null || !setter.IsPublic)
                {
                    continue; // get-only or non-public setter — immutable from the outside
                }

                var isInitOnly = setter.ReturnParameter
                    .GetRequiredCustomModifiers()
                    .Any(m => m.FullName == isExternalInit);

                if (!isInitOnly)
                {
                    offenders.Add($"{type.FullName}.{property.Name}");
                }
            }
        }

        Assert.That(offenders, Is.Empty,
            "Result/DTO types must be immutable (records or get-only). Mutable public setters found on: "
            + string.Join(", ", offenders));
    }

    // ---------------------------------------------------------------------
    // Rules 20-22 — placement & naming-convention guards (anti-drift)
    //
    // These close the gap exposed when the algorithm roster grew to 258 and
    // several algorithm classes were found to have moved modules (e.g. codon
    // tools → MolTools, miRNA/splicing/epigenetics → Annotation) without any
    // test catching it. They enforce *where* algorithm types may live so a
    // future relocation into the wrong module fails the build.
    // ---------------------------------------------------------------------

    /// <summary>
    /// Algorithm-class name suffixes. A class whose name ends with one of these is an
    /// algorithm entry point (parser, predictor, analyzer, …), never a Core primitive.
    /// </summary>
    private static readonly string[] AlgorithmSuffixes =
    [
        "Parser", "Predictor", "Finder", "Designer", "Analyzer",
        "Optimizer", "Caller", "Aligner", "Assembler",
    ];

    /// <summary>
    /// Rule 20: File-format parsers live only in the IO module. Parsing is an IO concern;
    /// a <c>*Parser</c> appearing in any other module signals that file-format handling has
    /// leaked into algorithm or domain code.
    /// </summary>
    [Test]
    public void Parsers_ShouldResideIn_IO()
    {
        Classes().That().HaveNameEndingWith("Parser")
            .Should().ResideInNamespace(Io)
            .Check(Arch);
    }

    /// <summary>
    /// Rule 21: Core is the innermost domain layer and holds primitives only (sequences,
    /// genetic code, extensions) — no algorithm entry points. No Core type may carry an
    /// algorithm-class suffix; such a class belongs in an algorithm module instead.
    /// </summary>
    [Test]
    public void Core_ShouldContainNoAlgorithmClasses()
    {
        var coreAssembly = typeof(DnaSequence).Assembly;

        var offenders = coreAssembly.GetExportedTypes()
            .Where(t => t.IsClass)
            .Where(t => AlgorithmSuffixes.Any(s => t.Name.EndsWith(s, StringComparison.Ordinal)))
            .Select(t => t.FullName)
            .ToList();

        Assert.That(offenders, Is.Empty,
            "Core must contain primitives only; algorithm classes belong in an algorithm module. "
            + "Misplaced in Core: " + string.Join(", ", offenders));
    }

    /// <summary>
    /// Rule 22: A type's namespace must match the module (assembly) that physically contains
    /// it — project boundary equals namespace boundary. This catches a class file being moved
    /// into the wrong project (or its namespace edited) independently of the dependency rules.
    /// </summary>
    [Test]
    public void Types_ShouldResideInNamespaceMatchingTheirAssembly()
    {
        var offenders = new List<string>();

        foreach (var assembly in SeqeronAssemblies)
        {
            var root = assembly.GetName().Name!; // e.g. "Seqeron.Genomics.Analysis"

            foreach (var type in assembly.GetExportedTypes())
            {
                if (type.IsNested)
                {
                    continue; // nested types inherit their declaring type's namespace
                }

                var ns = type.Namespace;
                if (ns is null || (ns != root && !ns.StartsWith(root + ".", StringComparison.Ordinal)))
                {
                    offenders.Add($"{type.FullName} (in {root})");
                }
            }
        }

        Assert.That(offenders, Is.Empty,
            "Every type must reside in a namespace matching its assembly (project = namespace). "
            + "Mismatches: " + string.Join(", ", offenders));
    }
}
