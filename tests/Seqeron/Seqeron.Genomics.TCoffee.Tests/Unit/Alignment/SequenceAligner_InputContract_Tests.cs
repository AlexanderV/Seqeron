// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-38, A2-39, A2-41, A2-42, A2-43; B09 F81, F82, F83, F84, F85): the input contract of the multiple
/// aligners and the overload-resolution contract of the scoring argument.
///
/// A2-38: null collection / null element = <see cref="ArgumentNullException"/> on EVERY overload (the ClustalW and MUSCLE DNA overloads used to throw
///   <see cref="NullReferenceException"/>); malformed text of a user tree / a T-Coffee method parameter = <see cref="ArgumentException"/>, never <see cref="FormatException"/>.
/// A2-39: an empty member is refused by the three program ports with ONE message (the programs refuse it), the Seqeron aligners return an all-gap row.
/// A2-41: the T-Coffee proba pair-HMM and the MUSCLE profile DP refuse too large inputs with an <see cref="ArgumentException"/> (no int overflow / OOM).
/// A2-42: an exception inside the parallel pair stage reaches the caller with its original type for every k, never wrapped in <see cref="AggregateException"/>.
/// A2-43: <c>null</c> for the scoring argument compiles (no CS0121) and selects the default scoring.
/// </summary>
[TestFixture]
public class SequenceAligner_InputContract_Tests
{
    private static readonly ProteinSequence[] WithNullProtein = { new("MKVLA"), null! };
    // ------------------------------------------------------------------ A2-38 null elements

    private static readonly string[] WithNullString = { "ACGTACGT", null! };
    private static readonly DnaSequence[] WithNullDna = { new("ACGTACGT"), null! };

    [Test]
    public void NullElement_IsArgumentNullException_OnEveryClustalWOverload()
    {
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(WithNullString));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(WithNullProtein));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(WithNullString, new ClustalWProteinOptions()));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW((IEnumerable<string>)null!));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW((IEnumerable<DnaSequence>)null!));
    }

    [Test]
    public void NullElement_IsArgumentNullException_OnTCoffeeProgram_AndTheSeqeronAligners()
    {
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram(WithNullString));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignExact(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignAnchored(WithNullDna));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlign(WithNullProtein, new MatrixScoring(SubstitutionMatrix.Blosum62, -2, -1)));
    }

    [Test]
    public void NullElement_Message_NamesTheIndex()
    {
        var ex = Assert.Throws<ArgumentNullException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(WithNullString));
        Assert.That(ex!.Message, Does.Contain("Sequence 1 is null"));
        Assert.That(ex.ParamName, Is.EqualTo("sequences"));
    }

    [Test]
    public void MalformedTCoffeeMethodAndWeightValues_AreArgumentException_NotFormatException()
    {
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-method slow_pair@EP@GOP@abc"));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-method slow_pair@EP@GEP@1.5"));
        Assert.Throws<ArgumentException>(() => TCoffeeProgramOptions.FromCommandLine("-seq_weight @1,x,3"));
    }

    [Test]
    public void ClustalWUserMatrixParse_StaysFormatException_ByDocumentedDesign()
    {
        Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("not a matrix", "text"));
    }

    // ------------------------------------------------------------------ A2-39 empty members

    [Test]
    public void EmptyMember_ThreeProgramPorts_ShareOneMessage()
    {
        var withEmpty = new[] { "ACGTACGT", "", "ACGTACCT" };
        var c = Assert.Throws<ArgumentException>(() => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(withEmpty));
        var m = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMuscle(withEmpty));
        var t = Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram(withEmpty));
        Assert.That(c!.Message, Does.StartWith("ClustalW cannot align an empty sequence (sequence 1 is empty)"));
        Assert.That(m!.Message, Does.StartWith("MUSCLE cannot align an empty sequence (sequence 1 is empty)"));
        Assert.That(t!.Message, Does.StartWith("t_coffee cannot align an empty sequence (sequence 1 is empty)"));
        foreach (var ex in new[] { c, m, t })
            Assert.That(ex.Message, Does.Contain("MultipleAlign, MultipleAlignProgressive").And.Contain("return an all-gap row"));
    }

    [Test]
    public void EmptyMember_SeqeronAligners_ReturnAnAllGapRow()
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence(""), new DnaSequence("ACGTACCT") };
        foreach (var r in new[]
        {
            SequenceAligner.MultipleAlign(seqs), Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs),
            Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(seqs), Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(seqs),
        })
        {
            Assert.That(r.AlignedSequences[1], Is.EqualTo(new string('-', r.AlignedSequences[0].Length)));
        }
    }

    // ------------------------------------------------------------------ A2-41 size guards

    [Test]
    public void TCoffeeProba_RefusesAnOversizedPair_BeforeAllocating()
    {
        // 9 000 x 9 000 protein: (9 001)^2 x 5 = 405 M floats > 2^28; nothing is allocated (the guard runs first)
        string big = new('A', 9000);
        var ex = Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.TcProbaPairEntries(big, big, TcSeqKind.Protein));
        Assert.That(ex!.Message, Does.Contain("proba_pair").And.Contain("slow_pair"));
        Assert.That(Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeProbaMaxMatrixFloats, Is.EqualTo(1L << 28));
    }

    // ------------------------------------------------------------------ A2-43 null scoring overloads

    [Test]
    public void NullScoring_Compiles_AndSelectsTheDefault_OnEveryPairwiseFamily()
    {
        const string a = "ACGTACGTTAGC", b = "ACGTCGTTAGC";
        // each line below was CS0121 (ambiguous between ScoringMatrix? and MatrixScoring) before OverloadResolutionPriority
        Assert.That(SequenceAligner.GlobalAlign(a, b, null).Score, Is.EqualTo(SequenceAligner.GlobalAlign(a, b).Score));
        Assert.That(SequenceAligner.GlobalAlignAffine(a, b, null).Score, Is.EqualTo(SequenceAligner.GlobalAlignAffine(a, b).Score));
        Assert.That(SequenceAligner.LocalAlign(a, b, null).Score, Is.EqualTo(SequenceAligner.LocalAlign(a, b).Score));
        Assert.That(SequenceAligner.LocalAlignAffine(a, b, null).Score, Is.EqualTo(SequenceAligner.LocalAlignAffine(a, b).Score));
        Assert.That(SequenceAligner.SemiGlobalAlign(a, b, null).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b).Score));
        Assert.That(SequenceAligner.SemiGlobalAlignAffine(a, b, null).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b).Score));
        Assert.That(SequenceAligner.SemiGlobalAlign(a, b, null, EndGaps.All).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, SequenceAligner.SimpleDna, EndGaps.All).Score));
        Assert.That(SequenceAligner.LocalAlignTopN(a, b, null, 3, 1), Is.EqualTo(SequenceAligner.LocalAlignTopN(a, b, SequenceAligner.SimpleDna, 3, 1)).And.Count.EqualTo(3));
        Assert.That(SequenceAligner.LocalAlignAffineTopN(a, b, null, 3, 1), Is.EqualTo(SequenceAligner.LocalAlignAffineTopN(a, b, SequenceAligner.SimpleDna, 3, 1)).And.Count.EqualTo(3));
        Assert.That(SequenceAligner.LocalAlignCoOptimal(a, b, null, 5), Is.EqualTo(SequenceAligner.LocalAlignCoOptimal(a, b, SequenceAligner.SimpleDna, 5)));
        Assert.That(SequenceAligner.GlobalAlignCoOptimal(a, b, null, 5), Is.EqualTo(SequenceAligner.GlobalAlignCoOptimal(a, b, SequenceAligner.SimpleDna, 5)));
        Assert.That(SequenceAligner.GlobalAlignLinearSpace(a, b, null).Score, Is.EqualTo(SequenceAligner.GlobalAlign(a, b).Score));
        Assert.That(SequenceAligner.LocalAlignTopNLinearSpace(a, b, null, 3, 1), Is.EqualTo(SequenceAligner.LocalAlignTopN(a, b, SequenceAligner.SimpleDna, 3, 1)));
        var dna = new[] { new DnaSequence(a), new DnaSequence(b) };
        Assert.That(SequenceAligner.MultipleAlign(dna, null, null).AlignedSequences, Has.Length.EqualTo(2));
        Assert.That(Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(dna, null, null).AlignedSequences, Has.Length.EqualTo(2));
        Assert.That(Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(new[] { a, b }, null, null).AlignedSequences, Has.Length.EqualTo(2));
        Assert.That(SequenceAligner.MultipleAlignExact(dna, null, null).AlignedSequences, Has.Length.EqualTo(2));
    }

    // ------------------------------------------------------------------ A3-4 empty LIST
    // Installed programs on an empty FASTA file: clustalw "No sequences in file. No alignment!", muscle3 "No sequences in input file", t_coffee exit 1, kalign exit 1,
    // poa exit 1, probcons exit 1 "No sequences read."; mafft exit 0 with an empty output. The ports follow: one message for the six refusing programs, MAFFT returns an empty result.
    [Test]
    public void EmptyList_SixProgramPorts_ShareOneMessage_MafftReturnsEmpty()
    {
        var empty = Array.Empty<string>();
        var errors = new (string Program, Action Call)[]
        {
            ("ClustalW", () => Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(empty)),
            ("MUSCLE", () => SequenceAligner.MultipleAlignMuscle(empty)),
            ("t_coffee", () => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram(empty)),
            ("ProbCons", () => Seqeron.Genomics.TCoffee.ProbConsAligner.MultipleAlignProbCons(empty)),
        };
        Assert.Multiple(() =>
        {
            foreach (var (program, call) in errors)
            {
                var ex = Assert.Throws<ArgumentException>(() => call(), program);
                Assert.That(ex!.Message, Does.StartWith(program + " cannot align an empty list of sequences: the program refuses an empty input ("), program);
                Assert.That(ex.Message, Does.Contain("MultipleAlignMafft and the Seqeron aligners").And.Contain("return an empty result."), program);
            }

            Assert.That(SequenceAligner.MultipleAlignMafft(empty).AlignedSequences, Is.Empty, "mafft: exit 0, empty output");
        });
    }}
