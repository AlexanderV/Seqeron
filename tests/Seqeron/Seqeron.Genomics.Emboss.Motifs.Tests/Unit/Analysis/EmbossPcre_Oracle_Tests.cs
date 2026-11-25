// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.Motifs;
using System.IO.Compression;
using System.Reflection;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B09 WP67 (F306): the regular expression engine of the library of EMBOSS 6.6.0 (<see cref="EmbossPcre"/>: <c>ajRegComp</c> / <c>ajRegExecC</c> = <c>pcre_compile</c> / <c>pcre_exec</c> of the bundled PCRE 7, no options) equals the library on the rows of
/// TestData/Alignment/emboss_pcre_oracle.tsv.gz (scripts/alignment/generate_emboss_pcre_oracle.py: expressions that the PROSITE conversion can write and expressions of the PCRE syntax in general, each searched in several subjects):
/// the expressions the library refuses are refused with its message and position, the others find the same first match (offset, length).
/// </summary>
[TestFixture]
public class EmbossPcre_Oracle_Tests
{
    private static string Unescape(string s)
    {
        var sb = new StringBuilder();
        int i = 0;
        while (i < s.Length)
        {
            char c = s[i++];
            if (c != '\\')
            {
                sb.Append(c);
                continue;
            }

            c = s[i++];
            switch (c)
            {
                case '\\': sb.Append('\\'); break;
                case 't': sb.Append('\t'); break;
                case 'n': sb.Append('\n'); break;
                case 'r': sb.Append('\r'); break;
                default:
                    sb.Append((char)Convert.ToInt32(s.Substring(i, 2), 16));
                    i += 2;
                    break;
            }
        }

        return sb.ToString();
    }

    [Test]
    public void Engine_EqualsTheLibrary_AllRows()
    {
        const string resource = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_pcre_oracle.tsv.gz";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, matches = 0, errors = 0, unsupported = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var f = line.Split('\t');
            string pattern = Unescape(f[0]), subject = Unescape(f[1]), expected = f[2];
            rows++;
            string got;
            try
            {
                var rx = EmbossPcre.Compile(pattern);
                var m = rx.Exec(subject);
                got = m is null ? "-" : $"M {m.Value.Offset} {m.Value.Length}";
            }
            catch (EmbossPcreException e)
            {
                got = $"E {e.Position} {e.Message}";
            }
            catch (NotSupportedException)
            {
                unsupported++;
                continue;
            }

            if (got != expected)
                failures.Add($"/{f[0]}/ on '{f[1]}': {got} expected {expected}");
            else if (expected[0] == 'M')
                matches++;
            else if (expected[0] == 'E')
                errors++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(25)));
        Assert.That(matches, Is.GreaterThan(200));
        Assert.That(errors, Is.GreaterThan(200));
        Assert.That(unsupported, Is.LessThan(rows / 10), "the constructs that are not ported are rare in the rows");
    }
}
