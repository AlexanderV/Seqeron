// PROBE-DESIGN-001 — Primer3 hybridization-probe picker, thermodynamic probe screens, probe Tm,
// molecular weight and nearest-neighbour extinction coefficient.
// Reference: primer3-py 2.3.1 design_primers(PRIMER_TASK = pick_hyb_probe_only), calc_tm;
//            Biopython 1.88 Bio.SeqUtils.molecular_weight; Cantor, Warshaw & Shapiro (1970) ε260 table.
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class ProbeDesigner_Primer3Probe_Tests
{
    // random.seed(1); 120 nt over ACGT (see B07.md F17).
    private const string Template =
        "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACC";

    [Test]
    public void ExtinctionCoefficientNearestNeighbor_MatchesCantorTable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("ACGT"), Is.EqualTo(40300.0));
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("A"), Is.EqualTo(15400.0));
            // AA AA: 2·27400 − 15400 = 39400
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("aaa"), Is.EqualTo(39400.0));
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("ACGU", isDna: false), Is.EqualTo(41300.0));
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("ACGU"), Is.NaN);
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor(""), Is.EqualTo(0.0));
        });
    }

    // Witness probes (length and G+C inside the preset window) whose Tm on the preset's own scale lies inside its Tm
    // window; expected Tm = primer3-py 2.3.1 calc_tm(seq, mv, dv 0, dntp 0, dna, max_nn_length).
    private static IEnumerable<TestCaseData> PresetWitnesses()
    {
        yield return new TestCaseData("Microarray", "AGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAAT", 86.25468810488348);
        yield return new TestCaseData("qPCR", "AGGAGCTTCACATCTGGCGCCGTGTGCCT", 68.58490905493761);
        yield return new TestCaseData("FISH", string.Concat(Enumerable.Repeat("ACGT", 50)), 77.4029020719779);
        yield return new TestCaseData("NorthernBlot", string.Concat(Enumerable.Repeat("ACGT", 25)), 74.4029020719779);
        yield return new TestCaseData("SouthernBlot", string.Concat(Enumerable.Repeat("AACGT", 40)), 73.30290207197791);
    }
}
