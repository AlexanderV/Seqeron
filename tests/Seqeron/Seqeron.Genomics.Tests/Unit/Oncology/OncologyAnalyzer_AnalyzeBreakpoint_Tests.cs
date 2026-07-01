// ONCO-FUSION-003 — Fusion Breakpoint Analysis
// Evidence: docs/Evidence/ONCO-FUSION-003-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-FUSION-003.md
// Source: Uhrig et al. (2021). Genome Research 31(3):448-460 (Arriba output spec:
//         https://github.com/suhrig/arriba/wiki/05-Output-files — reading_frame, site).
//         Murphy & Elemento (2016). AGFusion (model.py: cds_5prime + cds_3prime; translate;
//         protein[0:find("*")]; out-of-frame trims to whole codons).
//         Badger & Olsen (1999). Mol Biol Evol 16(4):512-524 (reading frame in triplets).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Site = Seqeron.Genomics.Oncology.OncologyAnalyzer.BreakpointSite;
using Frame = Seqeron.Genomics.Oncology.OncologyAnalyzer.BreakpointFrameStatus;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_AnalyzeBreakpoint_Tests
{
    private static OncologyAnalyzer.FusionBreakpoint Bp(
        Site site5, Site site3, int fivePrimeCodingBases, int threePrimeStartPhase) =>
        new("EML4", "ALK", site5, site3, fivePrimeCodingBases, threePrimeStartPhase);

    #region AnalyzeBreakpoint — reading-frame consequence

    // M1 — (9 - 0) % 3 == 0, both CDS → in-frame (AGFusion frame rule; triplet reading).
    [Test]
    public void AnalyzeBreakpoint_InFramePhase0_InFrame()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 9, 0));

        Assert.Multiple(() =>
        {
            Assert.That(result.FrameStatus, Is.EqualTo(Frame.InFrame),
                "(9 - 0) mod 3 == 0: the 3' partner stays in phase, so the junction is in-frame.");
            Assert.That(result.BreakpointInCoding, Is.True, "Both breakpoints are CDS → coding-to-coding junction.");
        });
    }

    // M2 — (10 - 1) % 3 == 0 → in-frame.
    [Test]
    public void AnalyzeBreakpoint_InFramePhase1_InFrame()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 10, 1));

        Assert.That(result.FrameStatus, Is.EqualTo(Frame.InFrame),
            "(10 - 1) mod 3 == 0 keeps the 3' partner in phase (AGFusion frame rule).");
    }

    // M3 — (11 - 2) % 3 == 0 → in-frame.
    [Test]
    public void AnalyzeBreakpoint_InFramePhase2_InFrame()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 11, 2));

        Assert.That(result.FrameStatus, Is.EqualTo(Frame.InFrame),
            "(11 - 2) mod 3 == 0 keeps the 3' partner in phase (AGFusion frame rule).");
    }

    // M4 — (10 - 0) % 3 == 1 != 0 → out-of-frame.
    [Test]
    public void AnalyzeBreakpoint_OutOfFrame_OutOfFrame()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 10, 0));

        Assert.That(result.FrameStatus, Is.EqualTo(Frame.OutOfFrame),
            "(10 - 0) mod 3 == 1 != 0: the reading frame is shifted, so the junction is out-of-frame.");
    }

    // M5 — (9 - 1) % 3 == 2 != 0 → out-of-frame.
    [Test]
    public void AnalyzeBreakpoint_OutOfFramePhase1_OutOfFrame()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 9, 1));

        Assert.That(result.FrameStatus, Is.EqualTo(Frame.OutOfFrame),
            "(9 - 1) mod 3 == 2 != 0: out-of-frame (AGFusion frame rule).");
    }

    // M6 — 3' breakpoint in 5'UTR → no coding-to-coding junction → frame NotPredicted (Arriba '.').
    [Test]
    public void AnalyzeBreakpoint_ThreePrimeUtr_NotPredicted()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.FivePrimeUtr, 9, 0));

        Assert.Multiple(() =>
        {
            Assert.That(result.FrameStatus, Is.EqualTo(Frame.NotPredicted),
                "A non-CDS 3' breakpoint cannot join two frames → Arriba reading_frame = '.'.");
            Assert.That(result.BreakpointInCoding, Is.False, "Only one breakpoint is CDS.");
        });
    }

    // M7 — 5' breakpoint in 3'UTR → frame NotPredicted (Arriba '.').
    [Test]
    public void AnalyzeBreakpoint_FivePrimeUtr_NotPredicted()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.ThreePrimeUtr, Site.Cds, 9, 0));

        Assert.That(result.FrameStatus, Is.EqualTo(Frame.NotPredicted),
            "A non-CDS 5' breakpoint cannot join two frames → Arriba reading_frame = '.'.");
    }

    // C1 — the analysis carries the partner symbols through unchanged.
    [Test]
    public void AnalyzeBreakpoint_PreservesPartners()
    {
        var result = OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 9, 0));

        Assert.Multiple(() =>
        {
            Assert.That(result.Gene5Prime, Is.EqualTo("EML4"), "5' partner is carried through.");
            Assert.That(result.Gene3Prime, Is.EqualTo("ALK"), "3' partner is carried through.");
            Assert.That(result.Site5Prime, Is.EqualTo(Site.Cds), "5' site is carried through.");
            Assert.That(result.Site3Prime, Is.EqualTo(Site.Cds), "3' site is carried through.");
        });
    }

    // S2 — phase outside {0,1,2} for a CDS-CDS junction is invalid (delegates to IsInFrame validation).
    [Test]
    public void AnalyzeBreakpoint_InvalidPhase_Throws()
    {
        Assert.That(() => OncologyAnalyzer.AnalyzeBreakpoint(Bp(Site.Cds, Site.Cds, 9, 3)),
            NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(),
            "A coding-start phase must be 0, 1, or 2 (a codon offset).");
    }

    #endregion

    #region PredictFusionProtein — chimeric CDS, translation, first-stop truncation

    // M8 — chimeric ATGAAA|GATGGT → ATG AAA GAT GGT → MKDG; in-frame, no stop.
    [Test]
    public void PredictFusionProtein_InFrameNoStop_TranslatesFullPeptide()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));

        Assert.Multiple(() =>
        {
            Assert.That(p.Peptide, Is.EqualTo("MKDG"),
                "ATG=M, AAA=K, GAT=D, GGT=G (NCBI Table 1); AGFusion translate of the chimeric CDS.");
            Assert.That(p.Effect, Is.EqualTo(Frame.InFrame), "6 - 0 mod 3 == 0 → in-frame.");
            Assert.That(p.HasPrematureStop, Is.False, "No stop codon before the end of the ORF.");
        });
    }

    // M9 — chimeric ATGAAA|GATTAAGGT → ATG AAA GAT TAA(stop) → MKD; premature stop.
    [Test]
    public void PredictFusionProtein_PrematureStop_TruncatesAtStop()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATTAAGGT"));

        Assert.Multiple(() =>
        {
            Assert.That(p.Peptide, Is.EqualTo("MKD"),
                "TAA is a stop codon; AGFusion truncates at the first stop: protein[0:find('*')].");
            Assert.That(p.HasPrematureStop, Is.True, "A stop codon was reached before the ORF end.");
        });
    }

    // S1 — the premature-stop case also flips HasPrematureStop, the StopCodon signal of Arriba.
    [Test]
    public void PredictFusionProtein_PrematureStop_FlagsStopCodon()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var noStop = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));
        var withStop = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATTAAGGT"));

        Assert.Multiple(() =>
        {
            Assert.That(noStop.HasPrematureStop, Is.False, "No stop → flag false (Arriba: not stop-codon).");
            Assert.That(withStop.HasPrematureStop, Is.True, "Stop reached → flag true (Arriba reading_frame stop-codon).");
        });
    }

    // M10 — 5' prefix off the codon boundary, 3' starts at native phase 0: ATGA(4 ≡ phase 1) + AAGGT(phase 0).
    // The 3' partner is read shifted relative to its native frame, so this is OUT-OF-FRAME under the Arriba
    // reading_frame model ("whether the 3' gene is fused in-frame or out-of-frame"): (4 - 0) mod 3 == 1 ≠ 0.
    // (AGFusion would label the *contiguous ORF* "in-frame (with mutation)" because len(chimeric) is a multiple
    // of 3, but the repo models Arriba's two-way 3'-gene-frame call, not AGFusion's three-way ORF-continuity call.)
    [Test]
    public void PredictFusionProtein_MidCodonJunctionPhaseMismatch_OutOfFrame()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 4, threePrimeStartPhase: 0);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGA", "AAGGT"));

        Assert.Multiple(() =>
        {
            Assert.That(p.ChimericCds, Is.EqualTo("ATGAAAGGT"),
                "5' prefix [0:4]='ATGA' joined to 3' suffix [0:]='AAGGT' (AGFusion concat).");
            Assert.That(p.Effect, Is.EqualTo(Frame.OutOfFrame),
                "(4 - 0) mod 3 == 1 ≠ 0: the 3' partner is read frameshifted → out-of-frame (Arriba reading_frame).");
            Assert.That(p.Peptide, Is.EqualTo("MKG"),
                "ATG=M, AAA=K, GGT=G: the 9-base chimeric CDS still translates cleanly even though the 3' gene is frameshifted.");
        });
    }

    // M10b — genuinely in-frame mid-codon junction: 5' contributes 4 bases (phase 1) and the 3' suffix begins at
    // its native phase 1, so the frames are compatible: (4 - 1) mod 3 == 0 → in-frame. 3' CDS 'TAAGGT' sliced at
    // offset 1 gives the suffix 'AAGGT'; chimeric ATGA + AAGGT = ATGAAAGGT → ATG AAA GGT → MKG. (AGFusion frame
    // rule: len(cds5)=4 (1/3 frac) and len(cds3)=5 (2/3 frac) complement → in-frame.)
    [Test]
    public void PredictFusionProtein_MidCodonJunctionPhaseMatch_InFrame()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 4, threePrimeStartPhase: 1);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGA", "TAAGGT"));

        Assert.Multiple(() =>
        {
            Assert.That(p.ChimericCds, Is.EqualTo("ATGAAAGGT"),
                "5' prefix [0:4]='ATGA' joined to 3' suffix [1:]='AAGGT' (AGFusion concat).");
            Assert.That(p.Effect, Is.EqualTo(Frame.InFrame),
                "(4 - 1) mod 3 == 0: the 3' partner is read in its native frame → in-frame.");
            Assert.That(p.Peptide, Is.EqualTo("MKG"),
                "ATG=M, AAA=K, GGT=G: the junction completes a codon across the breakpoint.");
        });
    }

    // M11 — out-of-frame: ATGAA(5) + GATGGT(6) → ATGAAGATGGT(11), trim to 9 → ATG AAG ATG → MKM.
    [Test]
    public void PredictFusionProtein_OutOfFrame_TrimsAndShifts()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 5, threePrimeStartPhase: 0);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAA", "GATGGT"));

        Assert.Multiple(() =>
        {
            Assert.That(p.Effect, Is.EqualTo(Frame.OutOfFrame), "5 - 0 mod 3 == 2 != 0 → out-of-frame.");
            Assert.That(p.ChimericCds, Is.EqualTo("ATGAAGATGGT"), "5' prefix 'ATGAA' ++ 3' suffix 'GATGGT'.");
            Assert.That(p.Peptide, Is.EqualTo("MKM"),
                "AGFusion trims to whole codons (9): ATG=M, AAG=K, ATG=M; trailing 'GT' is dropped.");
        });
    }

    // M12 — chimeric CDS composition is the 5' prefix concatenated with the 3' suffix.
    [Test]
    public void PredictFusionProtein_ChimericCds_ConcatenatesPrefixSuffix()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATGGT"));

        Assert.That(p.ChimericCds, Is.EqualTo("ATGAAAGATGGT"),
            "ChimericCds = transcript1[0:6] + transcript2[0:] = 'ATGAAA' + 'GATGGT' (AGFusion).");
    }

    // S3 — empty 3' suffix (junction3 == 3' CDS length): peptide is the translation of the 5' prefix only.
    [Test]
    public void PredictFusionProtein_EmptyThreePrimeSuffix_PeptideFromFivePrime()
    {
        // junction3 = 0 here but the 3' CDS is empty, so the suffix is empty.
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 6, threePrimeStartPhase: 0);

        var p = OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", ""));

        Assert.Multiple(() =>
        {
            Assert.That(p.ChimericCds, Is.EqualTo("ATGAAA"), "Empty 3' suffix → chimeric CDS is the 5' prefix.");
            Assert.That(p.Peptide, Is.EqualTo("MK"), "ATG=M, AAA=K from the 5' prefix alone.");
        });
    }

    // Validation — null 3' CDS is invalid input.
    [Test]
    public void PredictFusionProtein_NullCds_Throws()
    {
        var bp = Bp(Site.Cds, Site.Cds, 6, 0);

        Assert.That(() => OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", null!)),
            NUnit.Framework.Throws.ArgumentNullException, "A CDS sequence must be supplied.");
    }

    // Validation — a 5' prefix length beyond the 5' CDS is out of range.
    [Test]
    public void PredictFusionProtein_OffsetOutOfRange_Throws()
    {
        var bp = Bp(Site.Cds, Site.Cds, fivePrimeCodingBases: 99, threePrimeStartPhase: 0);

        Assert.That(() => OncologyAnalyzer.PredictFusionProtein(bp, ("ATGAAA", "GATGGT")),
            NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(),
            "The 5' prefix length cannot exceed the 5' CDS length.");
    }

    #endregion
}
