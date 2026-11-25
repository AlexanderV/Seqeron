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
/// ALIGN-MULTI-001 (audit A2-14b / c / d, B09 F104-F106): the ClustalW 2.1 command line (<see cref="ClustalWCommandLine"/>, <c>CommandLineParser::checkParam / parseParams / setOptionalParam</c>) and the program as a function
/// (<see cref="ClustalWProgram"/>): option spelling, enumerated values matched like <c>findMatch</c>, ranges of the numeric options, the else-if chain of the verbs, files in and out. The oracle tests (output / structure / commands)
/// drive their runs through this parser, so every flag they use is verified against the program as well.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class ClustalWCommandLine_Tests
{
    [Test]
    public void Parse_AlignmentParameters_FollowTheProgramsRules()
    {
        var c = ClustalWCommandLine.Parse(["x.fa", "-type=dna", "-GAPOPEN=12.5", "/GapExt=3", "-MAXDIV=45", "-NOWEIGHTS", "-negative", "-quicktree", "-ktuple=3", "-pairgap=9", "-WINDOW=7", "-topdiags=6",
            "-pwgapopen=11", "-PWGAPEXT=0.5", "-transweight=0.25", "-clustering=upgma", "-iteration=tree", "-numiter=2", "-nopgap", "-nohgap", "-novgap", "-gapdist=6", "-hgapresidues=GPS", "-matrix=blosum", "-pwmatrix=pam"]);
        Assert.Multiple(() =>
        {
            Assert.That(c.Verb, Is.EqualTo(ClustalWVerb.Align));
            Assert.That(c.InFile, Is.EqualTo("x.fa"));
            Assert.That(c.Type, Is.EqualTo(ClustalWSequenceType.Dna));
            Assert.That(c.Dna.GapOpen, Is.EqualTo(12.5));
            Assert.That(c.Dna.GapExtend, Is.EqualTo(3.0));
            Assert.That(c.Dna.DivergenceCutoff, Is.EqualTo(45));
            Assert.That(c.Dna.NoWeights && c.Dna.NegativeMatrix && c.Dna.PairwiseMode == ClustalWPairwiseMode.Quick, Is.True);
            Assert.That(c.Dna.Quick!.KTuple, Is.EqualTo(3));
            Assert.That(c.Dna.Quick.PairGap, Is.EqualTo(9));
            Assert.That(c.Dna.Quick.Window, Is.EqualTo(7));
            Assert.That(c.Dna.Quick.TopDiagonals, Is.EqualTo(6));
            Assert.That(c.Dna.PairwiseGapOpen, Is.EqualTo(11.0));
            Assert.That(c.Dna.TransitionWeight, Is.EqualTo(0.25f));
            Assert.That(c.Dna.Clustering, Is.EqualTo(ClustalWClustering.Upgma));
            Assert.That(c.Dna.Iteration, Is.EqualTo(ClustalWIteration.Tree));
            Assert.That(c.Dna.NumIterations, Is.EqualTo(2));
            Assert.That(c.Protein.NoPreferredGapPenalties && c.Protein.NoHydrophilicGapPenalties && c.Protein.VariableRegionGapPenalties, Is.True);
            Assert.That(c.Protein.GapSeparationDistance, Is.EqualTo(6));
            Assert.That(c.Protein.HydrophilicResidues, Is.EqualTo("GPS"));
            Assert.That(c.Protein.Matrix, Is.EqualTo(ClustalWProfileMatrix.Blosum));
            Assert.That(c.Protein.Pairwise.Matrix, Is.EqualTo(ClustalWPairwiseMatrix.Pam));
            Assert.That(c.Protein.Pairwise.Quick.KTuple, Is.EqualTo(1), "a k-tuple of 3 is above the protein limit of 2: the program ignores it for protein");
        });
    }

    [Test]
    public void Parse_EnumeratedValues_AreMatchedByTheirFirstDistinguishingLetters()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ClustalWCommandLine.FindMatch("g", ["gcg", "gde", "pir"]), Is.EqualTo(-1), "'g' alone is ambiguous");
            Assert.That(ClustalWCommandLine.FindMatch("gc", ["gcg", "gde", "pir"]), Is.EqualTo(0));
            Assert.That(ClustalWCommandLine.FindMatch("p", ["gcg", "gde", "pir"]), Is.EqualTo(2));
            Assert.That(ClustalWCommandLine.FindMatch("x", ["gcg", "gde", "pir"]), Is.EqualTo(-1));
            var c = ClustalWCommandLine.Parse(["-INFILE=x", "-OUTPUT=phy", "-OUTORDER=in", "-CASE=up", "-SEQNOS=on", "-SEQNO_RANGE=on", "-RANGE=3:9", "-SECSTROUT=ma"]);
            Assert.That(c.Run.Output.Format, Is.EqualTo(ClustalWOutputFormat.Phylip));
            Assert.That(c.Run.OutputOrder, Is.EqualTo(ClustalWOutputOrder.Input));
            Assert.That(c.Run.Output.Lowercase, Is.False);
            Assert.That(c.Run.Output.SequenceNumbers && c.Run.Output.SequenceNumberRange, Is.True);
            Assert.That(c.Run.Output.Range, Is.EqualTo((3, 9)));
            Assert.That(c.Run.Output.Structure, Is.EqualTo(ClustalWStructureOutput.Mask));
        });
    }

    [Test]
    public void Parse_StructureParameters_AreRangeChecked()
    {
        var c = ClustalWCommandLine.Parse(["-INFILE=x", "-HELIXGAP=7", "-STRANDGAP=10", "-LOOPGAP=0", "-TERMINALGAP=5", "-HELIXENDIN=1", "-HELIXENDOUT=4", "-STRANDENDIN=2", "-STRANDENDOUT=3", "-NOSECSTR2"]);
        var s = c.Run.Structure;
        Assert.Multiple(() =>
        {
            Assert.That(s.HelixPenalty, Is.EqualTo(7));
            Assert.That(s.StrandPenalty, Is.EqualTo(4), "10 is outside 1..9: ignored");
            Assert.That(s.LoopPenalty, Is.EqualTo(1), "0 is outside 1..9: ignored");
            Assert.That(s.TerminalPenalty, Is.EqualTo(5));
            Assert.That(s.HelixEndInside, Is.EqualTo(1));
            Assert.That(s.HelixEndOutside, Is.EqualTo(0), "4 is outside 0..3: ignored");
            Assert.That(s.StrandEndInside, Is.EqualTo(2));
            Assert.That(s.StrandEndOutside, Is.EqualTo(3));
            Assert.That(s.UseStructure1 && !s.UseStructure2, Is.True);
        });
    }

    [Test]
    public void Parse_Verbs_FollowTheElseIfChainOfTheProgram()
    {
        var files = new Dictionary<string, string> { ["t.dnd"] = "(a:0.1,b:0.1);" };
        Assert.Multiple(() =>
        {
            Assert.That(ClustalWCommandLine.Parse(["x.fa"]).Verb, Is.EqualTo(ClustalWVerb.Align), "nothing asked: align");
            Assert.That(ClustalWCommandLine.Parse(["x.fa", "-NEWTREE=t.dnd"]).Verb, Is.EqualTo(ClustalWVerb.GuideTree));
            Assert.That(ClustalWCommandLine.Parse(["x.fa", "-ALIGN", "-NEWTREE=t.dnd"]).Verb, Is.EqualTo(ClustalWVerb.Align), "-ALIGN wins over -NEWTREE");
            var use = ClustalWCommandLine.Parse(["x.fa", "-USETREE=t.dnd"], files);
            Assert.That(use.Verb, Is.EqualTo(ClustalWVerb.AlignUsingTree));
            Assert.That(use.Run.UseTree, Is.EqualTo("(a:0.1,b:0.1);"));
            Assert.That(ClustalWCommandLine.Parse(["-PROFILE1=a", "-PROFILE2=b", "-PROFILE"]).Verb, Is.EqualTo(ClustalWVerb.Profile));
            Assert.That(ClustalWCommandLine.Parse(["-PROFILE1=a", "-PROFILE2=b", "-SEQUENCES"]).Verb, Is.EqualTo(ClustalWVerb.Sequences));
            Assert.That(ClustalWCommandLine.Parse(["x.fa", "-TREE"]).Verb, Is.EqualTo(ClustalWVerb.Tree));
            Assert.That(ClustalWCommandLine.Parse(["x.fa", "-BOOTSTRAP=100"]).Verb, Is.EqualTo(ClustalWVerb.Bootstrap));
            Assert.That(ClustalWCommandLine.Parse(["x.fa", "-CONVERT"]).Verb, Is.EqualTo(ClustalWVerb.Convert));
            var tw = ClustalWCommandLine.Parse(["-PROFILE1=a", "-PROFILE2=b", "-PROFILE", "-NEWTREE1=x"]);
            Assert.That(tw.Run.GuideTreeOnly, Is.True, "-NEWTREE1 stops after the profile trees");
        });
    }

    [Test]
    public void Parse_Errors_AreTheProgramsRefusals()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["x.fa", "-NOSUCHOPTION"]));
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["x.fa", "-GAPOPEN"]), "a value is required");
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["x.fa", "-GAPOPEN=abc"]));
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["x.fa", "-NUMITER=0"]));
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["x.fa", "-OUTPUT=zzz"]));
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["-ALIGN"]), "no input file");
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["x.fa", "-USETREE=missing.dnd"]), "the file's text is not given");
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["-PROFILE2=b"]), "profile 2 without profile 1");
            Assert.Throws<ArgumentException>(() => ClustalWCommandLine.Parse(["-PROFILE1=a", "-PROFILE2=b", "-SEQUENCES", "-NEWTREE1=x"]), "-NEWTREE1 is invalid with -SEQUENCES");
        });
    }

    [Test]
    public void Run_AlignWithGcgOutput_WritesTheMsfAndTheGuideTreeFiles()
    {
        // clustalw -INFILE=ex.fa -TYPE=DNA -ALIGN -OUTPUT=GCG -OUTORDER=INPUT writes ex.msf and ex.dnd
        const string fasta = ">human\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA\n>chimp\nACGTACGATTGACCAGTAGCCATCGATGGATTACA\n>mouse\nACGTAGGATTGACCGTAGCCATCGATGGTTACAC\n>rat\nACGTAGGATTGACCGTAGCCATCGAAGGTTACAC\n";
        var r = ClustalWProgram.Run(["-INFILE=ex.fa", "-TYPE=DNA", "-ALIGN", "-OUTPUT=GCG", "-OUTORDER=INPUT"], new Dictionary<string, string> { ["ex.fa"] = fasta });
        Assert.Multiple(() =>
        {
            Assert.That(r.Files.Keys, Is.EquivalentTo(new[] { "ex.msf", "ex.dnd" }));
            Assert.That(r.Files["ex.dnd"], Is.EqualTo("(\n(\nhuman:0.05798,\nchimp:-0.00084)\n:0.05966,\nmouse:0.00000,\nrat:0.02941);\n"));
            Assert.That(r.Files["ex.msf"], Does.StartWith("PileUp\n\n\n\n   MSF:   36  Type: N    Check:  6005   .. \n"));
        });
    }

    [Test]
    public void Run_ProfileCommand_WritesBothProfileTrees_AndUnsupportedVerbsAreRefused()
    {
        const string p1 = ">a1\nACGTACGTTAGCTAGCTAGG-ATCGATCGAT\n>a2\nACGTACGATAGCTAGCTAGGGATCGTTCGAT\n>a3\nACGAACGTTAGCTTGCTAGGGATCGATCAAT\n";
        const string p2 = ">b1\nTCGTACGTTAGCTAGCAAGGATCGATCGAT\n>b2\nACGTACGTTAGGTAGCTAGGTTCGATCGAA\n>b3\nACGTACGTTAGGTAGCTAGGTTCGATGGAA\n";
        var r = ClustalWProgram.Run(["-PROFILE1=p1.fa", "-PROFILE2=p2.fa", "-PROFILE", "-TYPE=DNA", "-OUTPUT=FASTA", "-OUTFILE=o.fa"], new Dictionary<string, string> { ["p1.fa"] = p1, ["p2.fa"] = p2 });
        Assert.Multiple(() =>
        {
            // clustalw writes p1.dnd and p2.dnd next to the profile files (an explicit -NEWTREE1 / -NEWTREE2 name is ignored by the program)
            Assert.That(r.Files["p1.dnd"], Is.EqualTo("(\na1:0.00000,\na2:0.07000,\na3:0.10000);\n"));
            Assert.That(r.Files["p2.dnd"], Is.EqualTo("(\nb1:0.00000,\nb2:0.07000,\nb3:0.10000);\n"));
            Assert.That(r.Files["o.fa"], Does.StartWith(">a1\nACGTACGTTAGCTAGCTAGG-ATCGATCGAT\n"));
            Assert.Throws<NotSupportedException>(() => ClustalWProgram.Run(["x.fa", "-TREE"], new Dictionary<string, string> { ["x.fa"] = ">a\nAC\n>b\nAC\n" }));
            Assert.Throws<ArgumentException>(() => ClustalWProgram.Run(["x.fa"], new Dictionary<string, string>()));
            Assert.Throws<NotSupportedException>(() => ClustalWProgram.Run(["x.pir"], new Dictionary<string, string> { ["x.pir"] = ">P1;a\n\nAC*\n" }));
        });
    }
}
