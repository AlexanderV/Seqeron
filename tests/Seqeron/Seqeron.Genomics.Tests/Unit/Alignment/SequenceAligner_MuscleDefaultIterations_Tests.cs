using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A3-7, B09 F134): the program default of MUSCLE's <c>-maxiters</c> is 8 (muscle3 3.8.1551 <c>params.cpp</c> line 109, <c>unsigned g_uMaxIters = 8;</c>); the program's <c>-help</c> says
/// "default 16" (stale) and <see cref="MuscleOptions.MaxIterations"/> is 2 for the callers of <c>MultipleAlignMuscle</c>. <see cref="MuscleOptions.ProgramDefaults"/> is the option set of a plain
/// <c>muscle -in x.fa</c>. The oracle (scripts/alignment/generate_muscle3_default_iters_oracle.py, installed muscle3 3.8.1551) runs two 40-sequence protein families with no <c>-maxiters</c>, 8, 16 and 2.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MuscleDefaultIterations_Tests
{
    private sealed record Row(string Id, string Variant, string[] Sequences, string[] Rows);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.muscle3_default_iters_oracle.tsv";
        using var stream = typeof(SequenceAligner_MuscleDefaultIterations_Tests).Assembly.GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        reader.ReadLine();
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], f[2].Split(';'), f[3].Split(';')));
        }

        return rows;
    });

    private static string[] Run(Row r, MuscleOptions o) =>
        SequenceAligner.MultipleAlignMuscle(r.Sequences, o with { SequenceType = MuscleSequenceType.Protein }).AlignedSequences.ToArray();

    private static Row Get(string id, string variant) => Oracle.Value.First(r => r.Id == id && r.Variant == variant);

    [Test]
    public void Oracle_HasTwoFamiliesAndFourVariants()
    {
        Assert.That(Oracle.Value, Has.Count.EqualTo(8));
        Assert.That(Oracle.Value.All(r => r.Sequences.Length == 40 && r.Rows.Length == 40), Is.True);
    }

    [TestCase("D09")]
    [TestCase("D27")]
    public void ProgramDefaults_EqualAPlainMuscleRun(string id)
        => Assert.That(Run(Get(id, "default"), MuscleOptions.ProgramDefaults), Is.EqualTo(Get(id, "default").Rows), "muscle3 -in x.fa");

    [TestCase("D09", "8")]
    [TestCase("D09", "16")]
    [TestCase("D09", "2")]
    [TestCase("D27", "8")]
    [TestCase("D27", "16")]
    [TestCase("D27", "2")]
    public void MaxIterations_EqualsTheProgramsMaxiters(string id, string variant)
        => Assert.That(Run(Get(id, variant), new MuscleOptions { MaxIterations = int.Parse(variant) }), Is.EqualTo(Get(id, variant).Rows), $"muscle3 -maxiters {variant}");

    [Test]
    public void TheProgramDefaultIs8_NotTheHelpTextsSixteen_AndNotTheOptionsDefaultOfTwo()
    {
        var d09Default = Get("D09", "default").Rows;
        Assert.Multiple(() =>
        {
            Assert.That(d09Default, Is.EqualTo(Get("D09", "8").Rows), "no -maxiters = -maxiters 8");
            Assert.That(d09Default, Is.Not.EqualTo(Get("D09", "16").Rows), "the -help text ('default 16') is stale: the rows differ from -maxiters 16");
            Assert.That(d09Default, Is.Not.EqualTo(Get("D09", "2").Rows), "and from the draft + tree refinement of -maxiters 2");
            Assert.That(MuscleOptions.ProgramDefaults.MaxIterations, Is.EqualTo(8));
            Assert.That(new MuscleOptions().MaxIterations, Is.EqualTo(2), "the options default stays 2 for the existing callers; ProgramDefaults is the named profile");
        });
    }
}
