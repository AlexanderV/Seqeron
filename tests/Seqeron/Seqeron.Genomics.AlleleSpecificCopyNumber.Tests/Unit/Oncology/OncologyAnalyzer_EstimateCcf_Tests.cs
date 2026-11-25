// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
// ONCO-CCF-001 — Cancer Cell Fraction Estimation and CCF Clustering
// Evidence: docs/Evidence/ONCO-CCF-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CCF-001.md
// Source: McGranahan N et al. (2016). Science 351(6280):1463-1469. https://doi.org/10.1126/science.aaf1490
//         Tarabichi M et al. (2021). Nat. Methods 18:144-155 (Box 1). PMC7867630
//         Zheng J et al. (2022). Bioinformatics 38(15):3677-3683. https://doi.org/10.1093/bioinformatics/btac367
//         Lloyd SP (1982). IEEE Trans. Inf. Theory 28(2):129-137. https://doi.org/10.1109/TIT.1982.1056489
//
// Expected CCF values are computed independently from CCF = VAF·(rho·N_T + 2(1-rho))/(rho·m) — NOT
// copied from the implementation. Cluster centroids/assignments are derived from Lloyd's k-means on
// the sorted values; the clonal cluster is the highest centroid (Tarabichi 2021).

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_EstimateCcf_Tests
{


    // B24 review 2026-09 (ONCO-PURITY-001 dedup): EstimateCcf routes through the canonical CNAqc
    // purity/copy-number correction AdjustVAFForPurity (n_mut = VAF·(ρ·N_T + 2(1−ρ))/ρ, McGranahan 2016),
    // so RawCcf = n_mut / m exactly. Reference (Python, McGranahan 2016 formula):
    // VAF 0.3, ρ 0.7, N_T 3, m 2 ⇒ n_mut = 0.3·2.7/0.7 = 1.157142857…, CCF = 0.578571428…
    [Test]
    public void EstimateCcf_RawCcf_EqualsCanonicalPurityCorrectionOverMultiplicity()
    {
        OncologyAnalyzer.CcfEstimate estimate = OncologyAnalyzer.EstimateCcf(0.3, 0.7, 3, 2);

        Assert.Multiple(() =>
        {
            Assert.That(estimate.RawCcf, Is.EqualTo(OncologyAnalyzer.AdjustVAFForPurity(0.3, 0.7, 3) / 2));
            Assert.That(estimate.RawCcf, Is.EqualTo(0.3 * 2.7 / 0.7 / 2.0).Within(1e-12));
            Assert.That(AlleleSpecificCopyNumberAnalyzer.DeriveMultiplicity(0.3, 0.7, 3, 2), Is.EqualTo(1),
                "n_mut = 1.157… rounds to multiplicity 1 (McGranahan 2016)");
        });
    }
}
