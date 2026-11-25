// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using static Seqeron.Genomics.MolTools.CodonOptimizer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier fuzz tests for the entry points added by the B02 finisher
/// (docs/Validation/review-2026-09/B02.md, F31–F61, WP24/WP25). Random UTF-16 garbage (incl. U+017F 'ſ',
/// U+0131 'ı', U+212A Kelvin, surrogate halves, '\0'), random IUPAC/X nucleotide text and random option
/// values are fed to every new surface; each call must return a value satisfying its documented contract or
/// throw the documented exception (ArgumentException family) — never hang and never leak another runtime
/// exception (IndexOutOfRange, KeyNotFound, NullReference, InvalidOperation, …).
///
/// Test Units: TRANS-PROT-001, TRANS-CODON-001, TRANS-SIXFRAME-001, CODON-ENC-001, CODON-RSCU-001,
/// CODON-STATS-001, CODON-OPT-001, CODON-RARE-001, B02-SWEEP.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Translation")]
public class B02FinisherFuzzTests
{

    /// <summary>F44–F46, F53, F58: getorf modes / circular / transeq options on IUPAC + garbage input and random options.</summary>
    [Test]
    public void SixFrameAndGetOrf_RandomInputsAndOptions_DocumentedOutcomesOnly()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = rng.Next(4) == 0 ? Garbage(rng, rng.Next(0, 30)) : Nucleotides(rng, rng.Next(0, 60));
            var code = RandomCode(rng);
            var mode = (OrfFindMode)rng.Next(-1, 8);
            int minSize = rng.Next(-2, 40), maxSize = rng.Next(-2, 200), flanking = rng.Next(-2, 15);
            bool circular = rng.Next(2) == 0, both = rng.Next(2) == 0, met = rng.Next(2) == 0;
            Documented(() =>
            {
                var records = EmbossGetorf.FindOrfs(s, mode, code, minSize, maxSize, flanking, both, met, circular).ToList();
                int len = s.Length;
                foreach (var r in records)
                {
                    Assert.That(r.StartPosition, Is.InRange(0, Math.Max(0, len - 1)), $"getorf {trial}: start");
                    Assert.That(r.EndPosition, Is.GreaterThanOrEqualTo(r.StartPosition - 1), $"getorf {trial}: end");
                    if (!circular) Assert.That(r.EndPosition, Is.LessThan(len), $"getorf {trial}: linear end");
                    _ = r.ToGetOrfRange(len);
                }
            }, $"getorf {trial} mode {(int)mode} table {code.TableNumber}");
            Documented(() =>
            {
                foreach (var o in Translator.FindOrfs(s, code, rng.Next(-1, 10), both))
                    _ = o.ToGetOrfRange(s.Length);
            }, $"legacy orfs {trial}");

            var convention = (SixFrameConvention)rng.Next(-1, 4);
            Documented(() =>
            {
                var frames = EmbossTranslator.TranslateSixFrames(s, code, convention);
                Assert.That(frames.Keys, Is.EquivalentTo(new[] { 1, 2, 3, -1, -2, -3 }));
            }, $"six-frame {trial}");

            var regions = rng.Next(3) == 0
                ? null
                : Enumerable.Range(0, rng.Next(0, 3)).Select(_ => (rng.Next(-1, 40), rng.Next(-1, 80))).ToList();
            var options = new TranseqOptions { Trim = rng.Next(2) == 0, Clean = rng.Next(2) == 0, Methionine = rng.Next(2) == 0, Regions = regions };
            Documented(() =>
            {
                var frames = EmbossTranslator.TranslateSixFrames(s, code, convention, options);
                foreach (var p in frames.Values)
                {
                    if (options.Clean) Assert.That(p.Sequence, Does.Not.Contain("*"), $"transeq {trial}: clean");
                    if (options.Trim) Assert.That(p.Sequence.TrimEnd('X', '*'), Is.EqualTo(p.Sequence), $"transeq {trial}: trim");
                }
            }, $"transeq {trial}");
        }
    }

    private static GeneticCode RandomCode(Random rng) =>
        GeneticCode.GetByTableNumber(GeneticCode.SupportedTableNumbers[rng.Next(GeneticCode.SupportedTableNumbers.Count)]);

    // Runs the call; only the ArgumentException family is a documented failure.
    private static void Documented(Action call, string context)
    {
        try { call(); }
        catch (ArgumentException) { }
        catch (Exception e)
        {
            Assert.Fail($"{context}: undocumented {e.GetType().Name}: {e.Message}");
        }
    }

    private static string Nucleotides(Random rng, int length, string alphabet = "ACGTUacgtuRYSWKMBDHVNXx") =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private const int Seed = 20261014;

    private static string Garbage(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
        {
            chars[i] = rng.Next(6) switch
            {
                0 => (char)rng.Next(0, 0x10000),
                1 => (char)rng.Next(0, 128),
                2 => "ſKı\u0000￿𐀀"[rng.Next(7)],
                _ => "ACGTUNacgtunRYSWKMBDHVXx-.*"[rng.Next(27)],
            };
        }
        return new string(chars);
    }
}
