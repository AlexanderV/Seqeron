// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14, B09 F92): <see cref="ClustalWAligner.GuessClustalWSequenceType"/> = ClustalW 2.1's residue-type guess when <c>-TYPE</c> is omitted
/// (<c>FileReader::readFile</c>: the flag is set from the SECOND sequence only; <c>Sequence::checkDNAFlag</c>: DNA when >= 85 % of its non-gap characters are A C G T U N).
/// Oracle: the unit (<c>bp</c> / <c>aa</c>) of the "Sequence n" lines printed by the installed <c>clustalw</c> 2.1 on 500 random sets
/// (<c>scripts/alignment/generate_clustalw_type_oracle.py</c>, <c>clustalw_type_oracle.tsv</c>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWTypeGuess_Tests
{
    private static List<(string[] Sequences, ClustalWSequenceType Expected)> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_type_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<(string[], ClustalWSequenceType)>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            rows.Add((f[0].Split(';'), f[1] == "dna" ? ClustalWSequenceType.Dna : ClustalWSequenceType.Protein));
        }

        return rows;
    }

    [Test]
    public void Guess_EqualsTheTypeTheInstalledClustalWPrints_On500Sets()
    {
        var rows = ReadOracle();
        var failures = rows.Where(r => ClustalWAligner.GuessClustalWSequenceType(r.Sequences) != r.Expected)
            .Select(r => $"{string.Join(';', r.Sequences)} expected {r.Expected}").ToList();
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(480));
            Assert.That(rows.Count(r => r.Expected == ClustalWSequenceType.Dna), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Expected == ClustalWSequenceType.Protein), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Sequences.Length == 1), Is.GreaterThan(20), "a single sequence is protein (the flag is only set at index 1)");
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        });
    }

    [Test]
    public void Guess_IsDecidedBySecondSequenceOnly_AndSingleSequenceIsProtein()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ClustalWAligner.GuessClustalWSequenceType(["MKVLAAGIV", "ACGTACGTAC", "MKVLAAGIV"]), Is.EqualTo(ClustalWSequenceType.Dna), "ClustalW printed 'bp' for all three (checked)");
            Assert.That(ClustalWAligner.GuessClustalWSequenceType(["ACGTACGTAC", "ACGTRYGTAC"]), Is.EqualTo(ClustalWSequenceType.Protein), "8 of 10 = 80 % < 85 %");
            Assert.That(ClustalWAligner.GuessClustalWSequenceType(["ACGTACGTAC"]), Is.EqualTo(ClustalWSequenceType.Protein));
            Assert.That(ClustalWAligner.GuessClustalWSequenceType(["ACGTACGTAC", "acgtacgtac"]), Is.EqualTo(ClustalWSequenceType.Dna));
        });
    }
}
