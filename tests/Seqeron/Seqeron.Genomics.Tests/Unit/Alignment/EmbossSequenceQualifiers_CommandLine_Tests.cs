using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A6-2 (B09 F177): the forms of the sequence-qualifier text that the EMBOSS 6.6.0 programs accept (values measured with the installed programs: the "Commandline" echo of the program header and <c>EMBOSS_ACDLOG</c>;
/// the whole space is locked by the 2 080 + 640 row oracles of <see cref="EmbossAcdCommandLine_Oracle_Tests"/> and <see cref="EmbossSequenceQualifiersText_Oracle_Tests"/>).
/// </summary>
[TestFixture]
public class EmbossSequenceQualifiers_CommandLine_Tests
{
    private static (int Begin, int End, bool Reverse) B(EmbossSequenceQualifiers q) => (q.Begin, q.End, q.Reverse);

    [Test]
    public void NumberedForm_IsAccepted_ByAOneSequenceProgram()
    {
        // fuzznuc a.fa -pattern ACG -sbegin1 5: hits start at 5, exactly as -sbegin 5 (the numbered qualifier of the only sequence)
        Assert.That(EmbossSequenceQualifiers.ParseInputs("fuzznuc", "-sbegin1 5")[0].Begin, Is.EqualTo(5));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("fuzznuc", "-sbegin 5")[0].Begin, Is.EqualTo(5));
        Assert.That(EmbossSequenceQualifiers.Parse("-sbegin1 5").Begin, Is.EqualTo(5));
        var hits = GenomicAnalyzer.FindPatternHits("ACGTACGTACGT", EmbossSequenceQualifiers.ParseInputs("fuzznuc", "-sbegin1 5")[0], "ACG", 0, false);
        Assert.That(hits.Select(h => h.Position + 1), Is.EqualTo(new[] { 5, 9 }));
    }

    [Test]
    public void EqualsForm_And_NegatedForm_And_Values()
    {
        var q = EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1=5 -send2=-3 -nosreverse1 -sreverse2 N -slower2=yes -supper1 1 -SCIRCULAR1")[0];
        Assert.That((q.Begin, q.Reverse, q.Upper, q.Circular), Is.EqualTo((5, false, true, true)));
        var q2 = EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1=5 -send2=-3 -nosreverse1 -sreverse2 N -slower2=yes")[1];
        Assert.That((q2.End, q2.Reverse, q2.Lower), Is.EqualTo((-3, false, true)));
        // '--' / '/' prefixes, any case, 'name=value' without a dash
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "--SBEGIN1 4 /send1 9 send2=7").Select(x => (x.Begin, x.End)), Is.EqualTo(new[] { (4, 9), (0, 7) }));
    }

    [Test]
    public void UnnumberedQualifier_AfterTheFiles_BelongsToTheSecondSequence()
    {
        // needle a.fa b.fa -sbegin 3 -send 6 (Commandline echo: -sbegin2 3 -send2 6)
        var q = EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin 3 -send 6");
        Assert.That(q[0].IsDefault, Is.True);
        Assert.That((q[1].Begin, q[1].End), Is.EqualTo((3, 6)));
        foreach (string program in new[] { "water", "fuzznuc", "dotmatcher", "merger", "megamerger" })
        {
            var r = EmbossSequenceQualifiers.ParseInputs(program, "-sbegin 3");
            Assert.That(r[^1].Begin, Is.EqualTo(3), program);
            Assert.That(r.Count == 1 || r[0].IsDefault, Is.True, program);
        }
    }

    [Test]
    public void NumberedQualifierOfAnotherInput_EndsTheTie_AndLaterUnnumberedOnesApplyToAll()
    {
        // needle a.fa b.fa -sbegin1 5 -sbegin 7: the echo is "-sbegin1 5 -sbegin 7" and both sequences begin at 7 (the second qualifier is found without a master and applies to every input)
        var q = EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1 5 -sbegin 7");
        Assert.That((q[0].Begin, q[1].Begin), Is.EqualTo((7, 7)));
        q = EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin 3 -sbegin1 5 -send 9");
        Assert.That(((q[0].Begin, q[0].End), (q[1].Begin, q[1].End)), Is.EqualTo(((5, 9), (3, 9))));
    }

    [Test]
    public void Abbreviations_AreUnique_AmongAllQualifiersOfTheProgram()
    {
        // needle a.fa b.fa -sbeg1 3 = -sbegin1 3; -sre2 = -sreverse2; -sbe 3 (no number) = -sbegin on every input
        Assert.That(B(EmbossSequenceQualifiers.ParseInputs("needle", "-sbeg1 3")[0]), Is.EqualTo((3, 0, false)));
        var q = EmbossSequenceQualifiers.ParseInputs("needle", "-sre2");
        Assert.That((q[0].Reverse, q[1].Reverse), Is.EqualTo((false, true)));
        q = EmbossSequenceQualifiers.ParseInputs("needle", "-sbe 3");
        Assert.That((q[0].Begin, q[1].Begin), Is.EqualTo((3, 3)));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-nosre1")[0].Reverse, Is.False);

        // -s is ambiguous everywhere; -se is send in needle but ambiguous with -sequence in a program whose sequence parameter is named 'sequence'
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-s 3"))!.Message, Does.Contain("ambiguous qualifier").And.Contain("Unknown qualifier -s"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("fuzznuc", "-se 5"))!.Message, Does.Contain("ambiguous qualifier 'se' (send,sequence)"));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-se 5").Select(x => x.End), Is.EqualTo(new[] { 5, 5 }));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("fuzznuc", "-sen 5")[0].End, Is.EqualTo(5));
    }

    [Test]
    public void Refusals_LikeTheProgram()
    {
        // "Value required for '-sbegin' before ..." / unknown / "'no' prefix used with value" / "Invalid integer value" / "Invalid Y/N value"
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin"))!.Message, Does.StartWith("Value required for '-sbegin'"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin -send 5"))!.Message, Does.Contain("before '-send'"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sxyz 5"))!.Message, Does.Contain("Unknown qualifier -sxyz"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-nosreverse1=Y"))!.Message, Does.Contain("'no' prefix used with value"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-nosbegin1"))!.Message, Does.Contain("'no' prefix invalid"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1 abc"))!.Message, Does.Contain("Invalid integer value 'abc'"));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sreverse1=maybe"))!.Message, Does.Contain("Invalid Y/N value 'maybe'"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin3 5"));              // there is no third sequence
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("fuzznuc", "-sbegin2 5"));            // nor a second one
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "stray"));                   // a word is a parameter of the program
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-gapopen 5"));              // a qualifier of the program, not of a sequence
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-st"));                     // unique prefix of -stdout
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("nosuchprogram", "-sbegin 3"));
    }

    [Test]
    public void Values_FollowStrtol_And_AjStrToBool()
    {
        // strtol base 10: sign, digits, saturated to int, a 64-bit overflow is refused ('+3', '007', '010' = 10, 2147483648 = int.MaxValue measured as the last residue on needle 6.6.0)
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1 +3")[0].Begin, Is.EqualTo(3));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1 010")[0].Begin, Is.EqualTo(10));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1 2147483648")[0].Begin, Is.EqualTo(int.MaxValue));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin1 -9223372036854775808")[0].Begin, Is.EqualTo(int.MinValue));
        foreach (string bad in new[] { "9223372036854775808", "1e3", "0x10", "1.5", "5-", "--5", "" })
            Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", $"-sbegin1={bad}"), bad);
        // ajStrToBool: y yes t true 1 / n no f false / a number (true unless zero), nothing else; a separate value token only when it is y yes n no t true f false 0 1
        foreach (string yes in new[] { "y", "YES", "t", "True", "1", "2", "0.5", "0x1" })
            Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", $"-slower1={yes}")[0].Lower, Is.True, yes);
        foreach (string no in new[] { "n", "No", "f", "FALSE", "0", "00", "+0", "-0.0", "+nan" })
            Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", $"-slower1={no}")[0].Lower, Is.False, no);
        foreach (string bad in new[] { "ye", "tru", "fals", "yes1", "maybe", "1e999", "1e2", "inf", "nan" })
            Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", $"-slower1={bad}"), bad);
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sreverse1 N -slower1 F")[0].IsDefault, Is.True);
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-sreverse1 tru"));           // 'tru' is a stray word
    }

    [Test]
    public void IgnoredQualifiers_AreAccepted_AndNonZeroOffsetIsRefused()
    {
        var q = EmbossSequenceQualifiers.ParseInputs("needle", "-sformat1 fasta -sask1 -squick1 -sid1 x -sdbname2 y -iquery1 z -ufo1 u -fformat1 f -fopenfile2 g -sbegin1 3");
        Assert.That(B(q[0]), Is.EqualTo((3, 0, false)));
        Assert.That(q[1].IsDefault, Is.True);
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-ioffset1 0")[0].IsDefault, Is.True);
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-ioffset1 3"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceQualifiers.ParseInputs("needle", "-ioffset1 x"));
    }

    [Test]
    public void ExplicitMasterForm_And_NoQualifiers()
    {
        // -sbegin_asequence 5 (qualifier_master) is accepted by the program and ties the qualifier to that parameter
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin_asequence 5").Select(x => x.Begin), Is.EqualTo(new[] { 5, 0 }));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "-sbegin_bsequence 5").Select(x => x.Begin), Is.EqualTo(new[] { 0, 5 }));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", null).All(x => x.IsDefault), Is.True);
        Assert.That(EmbossSequenceQualifiers.ParseInputs("needle", "  ").Count, Is.EqualTo(2));
        Assert.That(EmbossSequenceQualifiers.ParseInputs("seqmatchall", "-sbegin 3").Count, Is.EqualTo(1));
    }
}
