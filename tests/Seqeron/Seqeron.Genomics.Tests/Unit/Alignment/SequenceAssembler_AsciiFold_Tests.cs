// ASSEMBLY-DBG-001 / ASSEMBLY-MERGE-001 / ASSEMBLY-SCAFFOLD-001 / ASSEMBLY-CONSENSUS-001 — B10 R1-WP1 (A15-7):
// the private char.ToUpperInvariant folds of SequenceAssembler were replaced by the canonical ASCII-only
// Infrastructure/SequenceCase.ToUpperAscii (consistent with B01 F21 / B02: U+017F 'ſ' is not 'S', EMBOSS and
// Biopython compare residue characters exactly and treat a non-ASCII character as an unknown residue), and the
// per-position Substring k-mer loops of BuildDeBruijnGraph / SpellStrandPairsOnce were replaced by the canonical
// KmerEnumerator + span dictionary lookups (allocation-free, result unchanged).
// Expected values are derived by hand from the cited rules (see each test).

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_AsciiFold_Tests
{
    // Biopython dumb_consensus: a residue is emitted only when it alone holds the maximum count and
    // count / numAtoms >= threshold. Case folding is ASCII-only: 'S', 'ſ' (U+017F), 's' -> S, ſ, S, so
    // S has 2/3 < 0.7 -> 'N'; with threshold 0.5: 2/3 >= 0.5 -> 'S'. (The old char.ToUpperInvariant fold
    // counted 3/3 for 'S' and returned "S" for both thresholds.)
    [Test]
    public void ComputeConsensus_LongS_IsNotCountedAsAsciiS()
    {
        var column = new[] { "S", "ſ", "s" };

        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.ComputeConsensus(column, 0.7, 'N'), Is.EqualTo("N"));
            Assert.That(SequenceAssembler.ComputeConsensus(column, 0.5, 'N'), Is.EqualTo("S"));
            Assert.That(SequenceAssembler.ComputeConsensus(column, 0.7, 'N', requireMultiple: false, ignoreCase: false),
                Is.EqualTo("N"), "case-sensitive: S x1, ſ x1, s x1 is a three-way tie");
        });
    }

    // EMBOSS merger_bestquality over "Aſ" / "ST", overlap 1: p = 'ſ', q = 'S' differ (ASCII fold), so the quality of the
    // 5-window around the overlap decides: a = (ſ=0) + (before: A=2) = 2 -> 2/11; b = (S=1: IUPAC strong/weak) + (T=2)
    // = 3 -> 3/11; b > a -> the second contig's base wins: "A" + "S" + "T" = "AST". (Old fold: ſ == S -> keeps 'ſ' -> "AſT".)
    [Test]
    public void MergeContigs_ResolveConflicts_LongSDiffersFromS_BestQualityDecides()
    {
        Assert.That(SequenceAssembler.MergeContigs("Aſ", "ST", 1, resolveOverlapConflicts: true), Is.EqualTo("AST"));
    }

    // ABySS createConsensus agreement is "equal ignoring ASCII case, or either base is N": 'S' vs 'ſ' disagree,
    // so the -1 gap (overlap 1) is not collapsed and the AGP unknown-size gap of 100 N is written instead.
    [Test]
    public void Scaffold_NegativeGap_LongSOverlapDisagrees_WritesUnknownGap()
    {
        var scaffold = SequenceAssembler.Scaffold(new[] { "AS", "ſT" }, new[] { (0, 1, -1) }).Single();

        Assert.That(scaffold, Is.EqualTo("AS" + new string('N', 100) + "ſT"));
    }

    // Control: the same overlap with an ASCII lower-case 's' agrees (case-insensitive) and is collapsed.
    [Test]
    public void Scaffold_NegativeGap_LowerCaseOverlapAgrees_Collapsed()
    {
        var scaffold = SequenceAssembler.Scaffold(new[] { "AS", "sT" }, new[] { (0, 1, -1) }).Single();

        Assert.That(scaffold, Is.EqualTo("AST"));
    }

    // Langmead DBG: a k-mer = edge prefix -> suffix over (k-1)-mer nodes, multi-edges kept; non-ASCII letters are
    // ordinary characters. Reads AİC and AİG (k = 3) -> node "Aİ" has edges to "İC" and "İG", in read order.
    [Test]
    public void BuildDeBruijnGraph_NonAsciiLetter_KeptVerbatim_SharedPrefixNodeAccumulates()
    {
        var graph = SequenceAssembler.BuildDeBruijnGraph(new[] { "AİC", "AİG" }, 3);

        Assert.Multiple(() =>
        {
            Assert.That(graph.Keys, Is.EqualTo(new[] { "Aİ" }));
            Assert.That(graph["Aİ"], Is.EqualTo(new[] { "İC", "İG" }));
        });
    }

    // Velvet read rule (velvetifySequence + writeNucleotide): a non-ACGT symbol, including the non-ASCII İ, is read as A before
    // the k-mers (and their reverse complements) are taken. Read ACİGT (k = 3) is the read ACAGT: k-mers ACA, CAG, AGT, and
    // their reverse complements TGT, CTG, ACT; the double-stranded graph has the single unitig ACAGT (its reverse complement
    // ACTGT is the same pair, reported once as the ordinally smaller string); the read has k-mers in it.
    [Test]
    public void AssembleDeBruijn_BothStrands_NonAsciiLetter_IsReadAsA_StrandPairEmittedOnce()
    {
        var result = SequenceAssembler.AssembleDeBruijn(
            new[] { "ACİGT" }, new SequenceAssembler.AssemblyParameters(KmerSize: 3, MinContigLength: 1), includeReverseComplement: true, clipTips: false);

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs, Is.EqualTo(new[] { "ACAGT" }));
            Assert.That(result.AssembledReads, Is.EqualTo(1));
        });
    }
}
