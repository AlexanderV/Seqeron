// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Repeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text.Json;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / GENOMIC-TANDEM-001 (audit A4-4, B09 F156): the option boundaries of the minimap2 chaining, mgaps clustering and etandem ports against the installed programs
/// (<c>option_boundary_oracle.json</c>, scripts/alignment/generate_option_boundary_oracle.py): every option value that the program accepts is accepted and gives the program's result (negative or zero bandwidth, gap limit,
/// counts, scores, separations; NaN / infinite / huge scales and fractions through the C float-to-integer conversion), the values it rejects or dies on are refused with a clear exception.
/// </summary>
[TestFixture]
public class OptionBoundaries_Oracle_Tests
{

    [Test]
    public void Etandem_ThresholdIntMinValue_DiesLikeTheProgram_WhenARepeatQualifies_AndBelowIsExact()
    {
        int crash = 0, ok = 0, emptyAtMin = 0;
        foreach (var row in Oracle.Value.RootElement.GetProperty("etandem").EnumerateArray())
        {
            string sequence = row.GetProperty("sequence").GetString()!;
            int min = row.GetProperty("minRepeat").GetInt32(), max = row.GetProperty("maxRepeat").GetInt32(), threshold = row.GetProperty("threshold").GetInt32();
            if (row.GetProperty("rc").GetInt32() != 0)
            {
                // etandem 6.6.0: "Floating point exception" (SIGFPE, status 136): the report root (best score threshold - 1 = int.MaxValue) is the top entry and its repeat size is 0
                var ex = Assert.Throws<InvalidOperationException>(() => EmbossRepeatPrograms.FindTandemRepeatsEtandem(sequence, min, max, threshold), $"{sequence} threshold {threshold}")!;
                Assert.That(ex.Message, Does.Contain("-2147483648").And.Contain("SIGFPE"));
                crash++;
                continue;
            }

            var hits = EmbossRepeatPrograms.FindTandemRepeatsEtandem(sequence, min, max, threshold)
                .Select(h => string.Join(",", h.Score, h.Start, h.End, h.RepeatSize, h.RepeatCount, CFormat.Fixed(h.Identity, 1), h.Consensus));
            Assert.That(string.Join("|", hits), Is.EqualTo(row.GetProperty("rows").GetString()), $"{sequence} threshold {threshold}");
            ok++;
            if (threshold == int.MinValue && !hits.Any()) emptyAtMin++;
        }

        Assert.That(crash, Is.EqualTo(5));
        Assert.That(ok, Is.EqualTo(7));
        Assert.That(emptyAtMin, Is.EqualTo(2), "a sequence without a qualifying repeat prints nothing: no crash");
    }

    private static readonly Lazy<JsonDocument> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.option_boundary_oracle.json";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        return JsonDocument.Parse(stream);
    });
}
