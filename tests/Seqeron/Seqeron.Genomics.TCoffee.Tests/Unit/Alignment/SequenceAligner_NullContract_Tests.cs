// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / ALIGN-STATS-001 (audit A5-3, B09 F169): the SECONDARY public entry points of the alignment classes refuse a null list / null element / out-of-range value with a typed
/// <see cref="ArgumentException"/> family exception naming the parameter (the library's contract; the main entry points were fixed in F81 / F82). Before: <see cref="NullReferenceException"/> /
/// <see cref="IndexOutOfRangeException"/> / <see cref="KeyNotFoundException"/> found by a reflection fuzzer over every public static method.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_NullContract_Tests
{

    private static void ThrowsNull(TestDelegate call, string paramName) =>
        Assert.That(Assert.Throws<ArgumentNullException>(call)!.ParamName, Is.EqualTo(paramName));

    // ------------------------------------------------------------------ (2) the program runners: a null element of the argument list

    [Test]
    public void KalignPoaProbConsRunners_NullArgumentOrNullFileText_IsAnArgumentNullException()
    {
        string?[] args = ["-i", null];
        var files = new Dictionary<string, string> { ["in.fa"] = null! };
        Assert.Multiple(() =>
        {
            ThrowsNull(() => ProbConsProgram.Parse(args!), "arguments");
            ThrowsNull(() => ProbConsProgram.Run(args!), "arguments");
            ThrowsNull(() => ProbConsProgram.Run(["in.fa"], files), "files");
            // the programs' own argument handling is untouched: an empty list is the usage screen, an unknown option is the program's error text
        });
    }

    // ------------------------------------------------------------------ (4) T-Coffee k-tuple similarity and the primary library: a null list / element

    [Test]
    public void TCoffeeKtupSimilarityAndPrimaryLibrary_NullListOrElement_IsAnArgumentNullException()
    {
        Assert.Multiple(() =>
        {
            ThrowsNull(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeKtupSimilarity(null!), "sequences");
            ThrowsNull(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeKtupSimilarity(["ACGTACGT", null!]), "sequences");
            ThrowsNull(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries(null!, [TCoffeePairMethod.SlowPair]), "sequences");
            ThrowsNull(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries(["ACGTACGT", null!], [TCoffeePairMethod.SlowPair]), "sequences");
            ThrowsNull(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries(["ACGTACGT", "ACGTACGA"], null!), "methods");
            Assert.That(Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries(["ACGTACGT", "ACGTACGA"], [TCoffeePairMethod.SlowPair]), Is.Not.Empty);
        });
    }

    // ------------------------------------------------------------------ (5) the T-Coffee pair stages: a character that is not a letter

    [TestCase("*")]
    [TestCase(" ")]
    [TestCase("-")]
    [TestCase("é")]
    [TestCase("中")]
    public void TCoffeePairStages_ACharacterThatIsNotALetter_IsAnArgumentException(string odd)
    {
        // measured: t_coffee 13.41 drops '*', '-' and blanks when it reads a sequence (rows of the cleaned sequence) and dies with "ERROR: COREDUMP" for U+00E9
        string bad = "ACGTAC" + odd + "GTTACG", good = "ACGTACGTTACG";
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeSlowPairAlign(bad, good))!.ParamName, Is.EqualTo("sequenceA"));
            Assert.That(Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeBiphasicPairAlign(good, bad))!.ParamName, Is.EqualTo("sequenceB"));
            Assert.That(Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeLalignTopAlignments(bad, good))!.ParamName, Is.EqualTo("sequenceA"));
            Assert.That(Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeSlowPairAlign(bad, good))!.Message, Does.Contain("index 6").And.Contain($"U+{(int)odd[0]:X4}"));
        });
        // the whole-program entry point refuses the same input (its contract: strings of letters)
        Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram([bad, good], TCoffeeProgramOptions.FromCommandLine("-method slow_pair")));
    }}
