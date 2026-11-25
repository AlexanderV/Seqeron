// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.SequenceFormats;
using System.IO.Compression;
using NUnit.Framework;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3e (B09 F241): the type and tag tables of the ten feature vocabularies of EMBOSS 6.6.0 as the ported loader (<c>featVocabReadTags</c> / <c>featVocabReadTypes</c>, the aliasing of the <c>name = id</c> lines included) builds
/// them from the data files equal the tables of the installed library (TestData/Alignment/emboss_feature_vocab_oracle.tsv.gz, scripts/alignment/emboss_feature_vocab_dump.c: <c>ajFeatVocabGetTypes</c> / <c>ajFeatVocabGetTags</c>, 2 406 entries).
/// </summary>
[TestFixture]
public class EmbossFeatureVocab_Oracle_Tests
{
    [Test]
    public void VocabularyTables_EqualTheLibrary_AllEntries()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_feature_vocab_oracle.tsv.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var expected = new HashSet<string>(StringComparer.Ordinal);
        var vocabularies = new HashSet<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            expected.Add(line);
            vocabularies.Add(line.Split('\t')[0]);
        }

        var actual = new HashSet<string>(StringComparer.Ordinal);
        foreach (string v in vocabularies)
        {
            var vocab = EmbossFeatureVocab.Get(v);
            foreach (var (k, val) in vocab.Types.Items)
                actual.Add($"{v}\ttypes\t{k}\t{val}");
            foreach (var (k, val) in vocab.Tags.Items)
                actual.Add($"{v}\ttags\t{k}\t{val}");
        }

        var missing = expected.Except(actual).OrderBy(x => x, StringComparer.Ordinal).Take(8).ToList();
        var extra = actual.Except(expected).OrderBy(x => x, StringComparer.Ordinal).Take(8).ToList();
        Assert.That(vocabularies, Has.Count.EqualTo(10));
        Assert.That(missing, Is.Empty, "in the library, not in the port:\n" + string.Join("\n", missing));
        Assert.That(extra, Is.Empty, "in the port, not in the library:\n" + string.Join("\n", extra));
        Assert.That(actual, Has.Count.EqualTo(expected.Count));
    }
}
