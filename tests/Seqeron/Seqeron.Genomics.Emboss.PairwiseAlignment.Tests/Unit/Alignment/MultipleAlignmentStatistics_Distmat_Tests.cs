// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.PairwiseAlignment;
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-STATS-001 (audit A2-15, B09 F90): EMBOSS 6.6.0 <c>distmat</c> corrections and options of <see cref="EmbossAlignmentStatistics.CalculateDistanceMatrix(IReadOnlyList{string}, MsaDistanceOptions)"/>:
/// Kimura (2-parameter and protein), Tamura, Tajima-Nei, Jin-Nei Gamma (<c>-calculatea</c>, <c>-parametera</c>), <c>-ambiguous</c>, <c>-position</c> 12 / 1 / 2 / 3, <c>-gapweight</c> -
/// 1 100 random alignments run through the installed <c>distmat</c> (<c>scripts/alignment/generate_distmat_oracle.py</c>; source: EMBOSS <c>emboss/distmat.c</c> and <c>ajax/core/ajbase.c</c>).
/// The program prints 100 D as a C <c>float</c> with <c>%.2f</c>; the port reproduces the single-precision arithmetic, so the printed strings must be EQUAL (nan / inf as printed, either sign of a nan).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-STATS-001")]
public class MultipleAlignmentStatistics_Distmat_Tests
{
    private sealed record Row(int Index, string[] Rows, string Type, int Method, double GapWeight, int Position, bool Ambiguous, bool CalcA, double ParamA, string Expected);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.distmat_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(list.Count, f[0].Split(','), f[1], int.Parse(f[2]), double.Parse(f[3], CultureInfo.InvariantCulture), int.Parse(f[4]), f[5] == "1", f[6] == "1",
                double.Parse(f[7], CultureInfo.InvariantCulture), f[8]));
        }

        return list;
    });

    private static MsaDistanceOptions Options(Row r) => new()
    {
        Method = (MsaDistanceMethod)r.Method,
        GapWeight = r.GapWeight,
        Position = r.Position,
        Ambiguous = r.Ambiguous,
        CalculateA = r.CalcA,
        ParameterA = r.ParamA,
    };

    /// <summary>distmat's print: D = D * (float)100, then <c>%.2f</c> (glibc: nan / -nan / inf / -inf).</summary>
    private static string Printed(double d)
    {
        float v = (float)d * 100f;
        if (float.IsNaN(v))
            return "nan";
        if (float.IsInfinity(v))
            return v > 0 ? "inf" : "-inf";
        return v.ToString("F2", CultureInfo.InvariantCulture);
    }

    private static string Normalised(string token) => token == "-nan" ? "nan" : token;

    public static IEnumerable<TestCaseData> Cases() => Oracle.Value.Select(r => new TestCaseData(r.Index).SetArgDisplayNames($"D{r.Index}_{r.Type}_m{r.Method}_p{r.Position}{(r.Ambiguous ? "_amb" : "")}"));

    [Test]
    public void Oracle_Has1100Sets_AllMethodsAndOptionsExercised()
    {
        var o = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(o, Has.Count.EqualTo(1100));
            for (int m = 0; m <= 5; m++)
                Assert.That(o.Count(r => r.Type == "dna" && r.Method == m), Is.GreaterThan(30), $"nucmethod {m}");
            for (int m = 0; m <= 2; m++)
                Assert.That(o.Count(r => r.Type == "protein" && r.Method == m), Is.GreaterThan(30), $"protmethod {m}");
            foreach (int pos in new[] { 123, 12, 1, 2, 3 })
                Assert.That(o.Count(r => r.Type == "dna" && r.Position == pos), Is.GreaterThan(20), $"position {pos}");
            Assert.That(o.Count(r => r.Ambiguous && r.Type == "dna"), Is.GreaterThan(20));
            Assert.That(o.Count(r => r.Ambiguous && r.Type == "protein"), Is.GreaterThan(20));
            Assert.That(o.Count(r => r.CalcA), Is.GreaterThan(10));
            Assert.That(o.Count(r => r.Expected == "ERROR"), Is.GreaterThan(10), "positions the program rejects");
            Assert.That(o.Count(r => r.Expected == "TIMEOUT"), Is.GreaterThan(3), "-ambiguous on nucleotides with a gap does not terminate");
            Assert.That(o.Sum(r => r.Expected.Split(',').Count(t => t.Length > 0 && char.IsDigit(t[^1]) )), Is.GreaterThan(3000), "finite printed values");
        });
    }

    [TestCaseSource(nameof(Cases))]
    public void CalculateDistanceMatrix_PrintedValuesEqualDistmat(int index)
    {
        var r = Oracle.Value[index];
        string why = $"D{index} {r.Type} method {r.Method} position {r.Position} gapweight {r.GapWeight} ambiguous {r.Ambiguous} calculatea {r.CalcA} a {r.ParamA}: {string.Join(',', r.Rows)}";
        var o = Options(r);
        if (r.Expected == "ERROR")
        {
            Assert.Throws<ArgumentException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(r.Rows, o), why);
            return;
        }

        if (r.Expected == "TIMEOUT")
        {
            Assert.Throws<NotSupportedException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(r.Rows, o), why);
            return;
        }

        // distmat passes mlen (not the selected length) to Tajima-Nei and to the protein Kimura method: with a position other than 123 it reads heap memory beyond its buffers. The program ran and
        // printed values (recorded), but they depend on the allocator; the port refuses instead of guessing.
        // (EMBOSS decides the residue type from the first row: a protein set whose first row only has nucleotide codes runs the nucleotide Kimura)
        bool nucleotide = r.Rows[0].All(c => "ACGTUBDHKMNRSVWXY?.~-".Contains(char.ToUpperInvariant(c)));
        bool readsBeyond = o.Position != 123 && (r.Method == 4 || (r.Method == 2 && !nucleotide));
        if (readsBeyond)
        {
            Assert.That(r.Expected, Does.Not.StartWith("ERROR").And.Not.StartWith("TIMEOUT"), why);
            Assert.Throws<NotSupportedException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(r.Rows, o), why);
            return;
        }

        var d = EmbossAlignmentStatistics.CalculateDistanceMatrix(r.Rows, o);
        var expected = r.Expected.Split(',');
        int n = r.Rows.Length, idx = 0;
        var got = new List<string>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++, idx++)
            {
                got.Add(Printed(d[i, j]));
                Assert.That(d[j, i], Is.EqualTo(d[i, j]).Or.NaN, why);
            }

        Assert.That(got, Is.EqualTo(expected.Select(Normalised)), why);
    }

    [Test]
    public void DefaultOverload_IsTheUncorrectedJukesCantorOfBefore()
    {
        string[] rows = ["ACGTACGTAC", "ACGTTCGTAC", "ACGAACGTCC"];
        var a = EmbossAlignmentStatistics.CalculateDistanceMatrix(rows, MsaDistanceMethod.JukesCantor);
        var b = EmbossAlignmentStatistics.CalculateDistanceMatrix(rows, new MsaDistanceOptions { Method = MsaDistanceMethod.JukesCantor });
        Assert.That(a, Is.EqualTo(b));
        Assert.That(EmbossAlignmentStatistics.CalculateDistanceMatrix(rows)[0, 1], Is.EqualTo(0.1).Within(1e-6));
    }

    [Test]
    public void NucleotideOnlyMethods_RejectProteins_AndBadPositionsAreRejected()
    {
        string[] prot = ["MKVLAAGIV", "MKVLTAGIV", "MRVLAAGIL"];
        foreach (var m in new[] { MsaDistanceMethod.Tamura, MsaDistanceMethod.TajimaNei, MsaDistanceMethod.JinNei })
            Assert.Throws<ArgumentException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(prot, new MsaDistanceOptions { Method = m }));
        Assert.Throws<ArgumentException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(["ACGT", "ACGA"], new MsaDistanceOptions { Position = 23 }));
        Assert.Throws<ArgumentException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(prot, new MsaDistanceOptions { Position = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossAlignmentStatistics.CalculateDistanceMatrix(prot, new MsaDistanceOptions { Method = (MsaDistanceMethod)9 }));
    }
}
