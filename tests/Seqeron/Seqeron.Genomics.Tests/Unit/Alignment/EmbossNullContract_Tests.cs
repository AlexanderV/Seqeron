using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / GENOMIC-TANDEM-001 (audit A6-4, B09 F179): the EMBOSS input sets / sequences and the report layer refuse a null element, a null <c>Name</c> / <c>Residues</c>, a hit without the program's coordinates and an undefined format with a typed
/// <see cref="ArgumentException"/>-family exception naming the parameter. Before: <see cref="NullReferenceException"/> (<c>EmbossSetInput.Any</c>, <c>EmbossReportWriter.Write</c>), <see cref="InvalidOperationException"/> ("Nullable object must have a value",
/// <c>EmbossReports.Fuzztran</c>) and a silent GFF for <c>(EmbossReportFormat)9999</c>.
/// </summary>
[TestFixture]
public class EmbossNullContract_Tests
{
    private static readonly EmbossSequence Good = new("g", "ACGTACGTAGCTAGCTAGGATCGATCGATCGAT");
    private static readonly EmbossSequence NoResidues = new("n", null!);
    private static readonly EmbossSequence NoName = new(null!, "ACGTACGTAGCT");

    private static void ParamOf<T>(TestDelegate call, string paramName, string messagePart) where T : ArgumentException
    {
        var e = Assert.Throws<T>(call)!;
        Assert.That(e.ParamName, Is.EqualTo(paramName));
        Assert.That(e.Message, Does.Contain(messagePart));
    }

    [Test]
    public void NeedleAll_And_Supermatcher_RefuseANullMemberOfASet_WithItsIndex()
    {
        ParamOf<ArgumentNullException>(() => SequenceAligner.NeedleAll([Good, null!], [Good]), "setA", "Sequence 1 is null.");
        ParamOf<ArgumentNullException>(() => SequenceAligner.NeedleAll([Good], [null!]), "setB", "Sequence 0 is null.");
        ParamOf<ArgumentNullException>(() => SequenceAligner.Supermatcher([null!], [Good]), "queries", "Sequence 0 is null.");
        ParamOf<ArgumentNullException>(() => SequenceAligner.Supermatcher([Good], [Good, Good, null!]), "targets", "Sequence 2 is null.");
        Assert.That(Assert.Throws<ArgumentNullException>(() => SequenceAligner.NeedleAll(null!, [Good]))!.ParamName, Is.EqualTo("setA"));
        Assert.That(Assert.Throws<ArgumentNullException>(() => SequenceAligner.Supermatcher([Good], null!))!.ParamName, Is.EqualTo("targets"));
        // the qualifier-aware reading takes the same road
        var withQualifiers = Good with { Qualifiers = new EmbossSequenceQualifiers { Begin = 2 } };
        ParamOf<ArgumentNullException>(() => SequenceAligner.NeedleAll([withQualifiers, null!], [Good]), "setA", "Sequence 1 is null.");
    }

    [Test]
    public void SetMembersWithoutResiduesOrName_AreRefused()
    {
        ParamOf<ArgumentNullException>(() => SequenceAligner.NeedleAll([Good, NoResidues], [Good]), "setA", "Sequence 1 has no Residues.");
        ParamOf<ArgumentNullException>(() => SequenceAligner.Supermatcher([Good], [NoName]), "targets", "Sequence 0 has no Name.");
        ParamOf<ArgumentNullException>(() => SequenceAligner.NeedleProgram(NoResidues, Good), "first", "no Residues");
        ParamOf<ArgumentNullException>(() => SequenceAligner.WaterProgram(Good, NoName), "second", "no Name");
        ParamOf<ArgumentNullException>(() => SequenceAligner.StretcherProgram(NoResidues, Good), "first", "no Residues");
        ParamOf<ArgumentNullException>(() => SequenceAligner.MatcherProgram(Good, NoResidues), "second", "no Residues");
        ParamOf<ArgumentNullException>(() => SequenceAligner.MergerProgram(NoName, Good), "first", "no Name");
        ParamOf<ArgumentNullException>(() => GenomicAnalyzer.MegaMerge(Good, NoResidues), "second", "no Residues");
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.NeedleProgram(null!, Good));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.MegaMerge(Good, null!));
    }

    [Test]
    public void SeqMatchAll_NullMember_NamesTheParameterAndTheIndex()
    {
        ParamOf<ArgumentNullException>(() => GenomicAnalyzer.SeqMatchAll(["ACGTACGTACGT", null!]), "sequences", "Sequence 1 is null.");
        ParamOf<ArgumentNullException>(() => GenomicAnalyzer.SeqMatchAll([null!, "ACGTACGTACGT"], EmbossSequenceQualifiers.None), "sequences", "Sequence 0 is null.");
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.SeqMatchAll(null!));
    }

    [Test]
    public void Fuzztran_HitWithoutEmbossCoordinates_IsAnArgumentException_NamingTheHit()
    {
        var real = GenomicAnalyzer.FindTranslatedPatternHits("ATGGCCTTAGGATCCATGAAATAG", "M", "F")[0];
        Assert.That(real.EmbossStart, Is.Not.Null);
        ParamOf<ArgumentException>(() => EmbossReports.Fuzztran("s", "M", [real, default]), "hits", "Hit 1 has no EmbossStart");
        ParamOf<ArgumentException>(() => EmbossReports.Fuzztran("s", "M", [real with { EmbossEnd = null }]), "hits", "Hit 0 has no EmbossStart / EmbossEnd");
        ParamOf<ArgumentException>(() => EmbossReports.Fuzztran("s", "M", [real with { Translation = null! }]), "hits", "Hit 0 has no Translation");
        Assert.That(EmbossReports.Fuzztran("s", "M", [real]).Features, Has.Count.EqualTo(1));
        Assert.Throws<ArgumentNullException>(() => EmbossReports.Fuzztran((string)null!, "M", [real]));
        Assert.Throws<ArgumentNullException>(() => EmbossReports.Etandem(null!, "ACGT", []));
    }

    [Test]
    public void ReportWriter_RefusesAnUndefinedFormat_AndNullFieldsOfTheReport()
    {
        var report = EmbossReports.Fuzznuc("s", "ACGTACGTACGT", "ACG", [new PatternHit(0, 3, false, 0)]);
        Assert.That(EmbossReportWriter.Write(report, EmbossReportFormat.Gff, new DateTime(2026, 1, 1, 0, 0, 0, DateTimeKind.Utc)), Does.StartWith("##gff-version 3"));
        Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => EmbossReportWriter.Write(report, (EmbossReportFormat)9999, DateTime.Now))!.ParamName, Is.EqualTo("format"));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossReportWriter.Write(report, (EmbossReportFormat)(-1), DateTime.Now));
        Assert.Throws<ArgumentNullException>(() => EmbossReportWriter.Write(null!, EmbossReportFormat.Gff, DateTime.Now));

        var feature = report.Features[0];
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Features = null! }, EmbossReportFormat.Gff, DateTime.Now), "report", "Features list is null");
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Features = [feature, null!] }, EmbossReportFormat.Excel, DateTime.Now), "report", "Feature 1 is null");
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Sequence = null! }, EmbossReportFormat.Excel, DateTime.Now), "report", "sequence");
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Features = [feature with { Tags = null! }] }, EmbossReportFormat.Gff, DateTime.Now), "report", "Feature 0");
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Features = [feature with { Tags = [null!] }] }, EmbossReportFormat.Excel, DateTime.Now), "report", "Feature 0");
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Features = [feature with { Type = null! }] }, EmbossReportFormat.Gff, DateTime.Now), "report", "Feature 0");
        ParamOf<ArgumentException>(() => EmbossReportWriter.Write(report with { Features = [feature with { Source = null! }] }, EmbossReportFormat.Gff, DateTime.Now), "report", "Feature 0");
    }
}
