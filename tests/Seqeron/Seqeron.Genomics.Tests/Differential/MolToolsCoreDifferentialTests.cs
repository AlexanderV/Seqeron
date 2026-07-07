// 08_DIFFERENTIAL_TESTING rows 18, 21, 26, 27 (MolTools core ops). Production methods vs INDEPENDENT
// oracles: a forward/reverse PAM scan + spec-reconstruction property, the closed-form Wallace/Marmur-Doty
// Tm formulas, an IndexOf restriction-site scan, and a manual cut-and-split digest. Reverse complement in
// the oracles is from a literal table, not DnaSequence.GetReverseComplementString.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class MolToolsCoreDifferentialTests
{
    private const double Tol = 1e-9;

    private static readonly Dictionary<char, char> Comp = new() { ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G' };
    private static string RevComp(string s)
    {
        var a = s.Select(c => Comp[c]).ToArray();
        Array.Reverse(a);
        return new string(a);
    }

    // ---- Row 18: CRISPR-PAM-001 — FindPamSites vs independent NGG scan + spec reconstruction ----

    [Test]
    [Category("CRISPR-PAM-001")]
    public void FindPamSites_SpCas9_MatchesIndependentScanAndSpec()
    {
        const string seq = "ACCGGTACCGGTAAGGCCTTGGCCAACCGGTTACGTACGTAGGCCTTAAGGTTCCGGAACC";
        int n = seq.Length;
        var sites = CrisprDesigner.FindPamSites(new DnaSequence(seq), CrisprSystemType.SpCas9).ToList();

        // Independent forward/reverse counts: NGG with a full 20 nt guide 5' of the PAM.
        int FwdCount(string s)
        {
            int c = 0;
            for (int i = 0; i + 3 <= s.Length; i++)
                if (s[i + 1] == 'G' && s[i + 2] == 'G' && i - 20 >= 0) c++;
            return c;
        }
        Assert.That(sites.Count(p => p.IsForwardStrand), Is.EqualTo(FwdCount(seq)), "forward count");
        Assert.That(sites.Count(p => !p.IsForwardStrand), Is.EqualTo(FwdCount(RevComp(seq))), "reverse count");

        // Each site is reconstructible from the spec, independently of the search internals.
        foreach (var p in sites)
        {
            Assert.That(p.PamSequence, Is.EqualTo(seq.Substring(p.Position, 3)), "PAM on forward coords");
            if (p.IsForwardStrand)
            {
                Assert.That(p.PamSequence[1] == 'G' && p.PamSequence[2] == 'G', Is.True, "NGG");
                Assert.That(p.TargetStart, Is.EqualTo(p.Position - 20));
                Assert.That(p.TargetSequence, Is.EqualTo(seq.Substring(p.Position - 20, 20)));
            }
            else
            {
                // Reverse PAM reads "CCN" on the forward strand; target is the revcomp of the 20 nt 3' of it.
                Assert.That(p.PamSequence.StartsWith("CC", StringComparison.Ordinal), Is.True, "revcomp NGG = CCN");
                Assert.That(p.TargetSequence, Is.EqualTo(RevComp(seq.Substring(p.Position + 3, 20))));
            }
        }
    }

    // ---- Row 21: PRIMER-TM-001 — CalculateMeltingTemperature vs closed-form Wallace / Marmur-Doty ----

    private static double TmOracle(string primer)
    {
        var s = primer.ToUpperInvariant();
        int at = s.Count(c => c == 'A' || c == 'T');
        int gc = s.Count(c => c == 'G' || c == 'C');
        int vl = at + gc;
        if (vl == 0) return 0;
        if (vl < 14) return 2.0 * at + 4.0 * gc;                  // Wallace rule
        return Math.Max(0, 64.9 + 41.0 * (gc - 16.4) / vl);       // Marmur-Doty
    }

    [Test]
    [Category("PRIMER-TM-001")]
    [TestCase("ACGT")]
    [TestCase("GCGCGCGC")]
    [TestCase("acgtACGT")]                       // mixed case
    [TestCase("ACGTNNN")]                        // N excluded from both counts
    [TestCase("ATGCATGCATGCATGCATGC")]           // 20-mer -> Marmur-Doty branch
    [TestCase("GGGGCCCCGGGGCCCCGGGG")]           // GC-rich long
    [TestCase("ATATATATATATATATATAT")]           // AT-rich long
    [TestCase("")]                               // empty -> 0
    public void MeltingTemperature_MatchesClosedFormOracle(string primer)
    {
        Assert.That(PrimerDesigner.CalculateMeltingTemperature(primer), Is.EqualTo(TmOracle(primer)).Within(Tol));
    }

    // ---- Row 26: RESTR-FIND-001 — FindSites vs IndexOf scan + cut arithmetic from enzyme metadata ----

    [Test]
    [Category("RESTR-FIND-001")]
    public void FindSites_EcoRI_MatchesIndexOfScan()
    {
        const string seq = "AAGAATTCGGGAATTCAA"; // two EcoRI (GAATTC) sites at 2 and 10
        int n = seq.Length;
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI")!;
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(seq), "EcoRI").ToList();

        // Independent forward occurrences via IndexOf.
        var fwd = new List<int>();
        for (int p = seq.IndexOf("GAATTC", StringComparison.Ordinal); p >= 0; p = seq.IndexOf("GAATTC", p + 1, StringComparison.Ordinal))
            fwd.Add(p);

        var fwdSites = sites.Where(s => s.IsForwardStrand).Select(s => (s.Position, s.CutPosition)).ToList();
        Assert.That(fwdSites, Is.EqualTo(fwd.Select(p => (p, p + enzyme.CutPositionForward)).ToList()), "forward sites + cuts");

        // EcoRI is palindromic: a reverse-strand site sits at every forward position, cut at pos+CutReverse.
        var revSites = sites.Where(s => !s.IsForwardStrand).Select(s => (s.Position, s.CutPosition)).OrderBy(x => x.Position).ToList();
        Assert.That(revSites, Is.EqualTo(fwd.Select(p => (p, p + enzyme.CutPositionReverse)).ToList()), "reverse sites + cuts");
    }

    // ---- Row 27: RESTR-DIGEST-001 — Digest vs manual cut-collection + split ----

    [Test]
    [Category("RESTR-DIGEST-001")]
    public void Digest_EcoRI_MatchesManualSplit()
    {
        const string seq = "AAGAATTCGGGAATTCAA";
        int n = seq.Length;
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI")!;
        var fragments = RestrictionAnalyzer.Digest(new DnaSequence(seq), "EcoRI")
            .Select(f => (f.Sequence, f.StartPosition, f.Length)).ToList();

        // Independent: forward-strand cut positions, then split [0, cuts..., n].
        var cuts = new SortedSet<int>();
        for (int p = seq.IndexOf("GAATTC", StringComparison.Ordinal); p >= 0; p = seq.IndexOf("GAATTC", p + 1, StringComparison.Ordinal))
            cuts.Add(p + enzyme.CutPositionForward);

        var boundaries = new List<int> { 0 };
        boundaries.AddRange(cuts);
        boundaries.Add(n);

        var expected = new List<(string, int, int)>();
        for (int i = 0; i < boundaries.Count - 1; i++)
        {
            int start = boundaries[i], len = boundaries[i + 1] - start;
            if (len > 0) expected.Add((seq.Substring(start, len), start, len));
        }

        Assert.That(fragments, Is.EqualTo(expected));
    }
}
