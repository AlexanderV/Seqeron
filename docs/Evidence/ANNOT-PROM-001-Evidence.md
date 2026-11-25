# Evidence Artifact: ANNOT-PROM-001

**Test Unit ID:** ANNOT-PROM-001
**Algorithm:** Bacterial σ70 promoter boxes (`GenomeAnnotator.FindPromoterMotifs` → `MotifFinder.FindSigma70Promoters`)
**Date Collected:** 2026-11-25 (review-2026-09 B11, F7)

---

## Online Sources (what was actually opened)

### Harley & Reynolds (1987) — Nucleic Acids Res 15(5):2343-2361
**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC340638 (WebFetch; abstract only — the body is page images, the PDF link returned a download stub)
**Authority rank:** 1.
1. 263 E. coli promoters with known start points; "all bases in the -35 (TTGACA) and -10 (TATAAT) hexamers were highly conserved".
2. "92% of promoters had inter-region spacing of 17 +/- 1 bp"; consensus for spacers 16, 17, 18 did not differ.
3. The 15–21 bp spacer range is not stated in the abstract; it is the default of the canonical `FindSigma70Promoters`
   (B05 F30, MOTIF-REGULATORY-001 Evidence) and is exposed as a parameter here.

### Hawley & McClure (1983) — Nucleic Acids Res 11(8):2237-2255
**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC325881 (WebFetch; abstract) — 112-promoter consensus compilation, H&R's initial reference list.

### Wikipedia "Promoter (genetics)" (WebFetch, 2026-11-25)
−10 consensus TATAAT, −35 consensus TTGACA, "optimal spacing … is 17 bp". The per-position percentage table
(−35 T69 T79 G61 A56 C54 A54; −10 T77 A76 T60 A61 A56 T82) carries **no citation** on the page; WebSearch did not
locate a primary source for it. The previous implementation's score (sum of matched-position percentages / total)
and its prefix/suffix variant library (TTGAC, TGACA, TTGA, TATAA, ATAAT, TATA) are not described in any source.

---

## Stage A verdict — ❌ FAIL (definition)

Consensus hexamers and 17-bp optimum correct. Defects: (1) a σ70 promoter is a −35/−10 pair with a 15–21 bp
spacer (H&R), but the method reported unpaired boxes; (2) the variant library (exact 4/5-nt prefixes/suffixes, with
GACA/TAAT missing) and (3) the score formula are unsourced; mismatch-tolerant boxes (H&R: natural promoters rarely
match the consensus) were not found.

## Stage B — implementation vs the validated definition

**Fix (F7):** delegate to `MotifFinder.FindSigma70Promoters` (canonical, B05); report each box of every pair once;
score = 1 − d/6 (consensus identity, from the Hamming distance the canonical method reports). The previous inline
`Substring` scan (DUP_MAP §18 `GenomeAnnotator.cs:798`) is removed.

**Oracle:** independent Python brute force of the definition (no Seqeron code) — for each i with
Hamming(S[i,i+6), TTGACA) ≤ mm35 and each spacer 15…21 with Hamming(S[j,j+6), TATAAT) ≤ mm10 add both boxes with
score (6 − d)/6; sort by (position, −35 first).

```python
C35, C10 = "TTGACA", "TATAAT"
def ham(a,b): return sum(x!=y for x,y in zip(a,b))
def boxes(seq, mm35=2, mm10=2, smin=15, smax=21):
    seq=seq.upper(); n=len(seq); out=set()
    for i in range(n-5):
        d35=ham(seq[i:i+6],C35)
        if d35>mm35: continue
        for s in range(smin,smax+1):
            j=i+6+s
            if j+6>n: break
            d10=ham(seq[j:j+6],C10)
            if d10>mm10: continue
            out.add((i,0,"-35 box",seq[i:i+6],(6-d35)/6)); out.add((j,1,"-10 box",seq[j:j+6],(6-d10)/6))
    return [(p,t,h,sc) for p,o,t,h,sc in sorted(out)]
```

**Cross-check (C# console harness → Annotation project vs oracle JSON):** 400 random sequences (seed 20261125,
length 0–120, uniform and AT-biased, 1 in 5 with a planted consensus pair at spacer 14–22, 1 in 7 lowercase;
limits (2,2), (1,1), (0,0), (3,2)) → **1,131 boxes, 0 mismatches**, bit-identical scores; the single-argument
overload equals the (2,2,15,21) overload on every (2,2) case.

**Locked oracle values:**
| Input | Result |
|-------|--------|
| lacUV5 control region, (1,0) | (68, −35, TTTACA, 5/6), (92, −10, TATAAT, 1) — published −35/−10, spacer 18 |
| lacUV5, defaults | (54, TAGGCA, 4/6), (68, TTTACA, 5/6), (77, TATGCT, 4/6), (92, TATAAT, 1) |
| `CCGATATTTTAAAAGTTTTTTTAAAGACTAACTATAAATTTGACATTATG` | (7 TTTAAA 4/6), (8 TTAAAA 4/6), (19 TTTAAA 4/6), (32 TATAAA 5/6), (34 TAAATT 4/6), (43 CATTAT 4/6) |
| `TTGACA C15 TATAAT TATAAT`, defaults | (0 TTGACA 1), (21 TATAAT 1), (24 AATTAT 4/6), (27 TATAAT 1) |
| `GG TGGCCA C17 TATAAT GG`, (3,2) | (1 GTGGCC 3/6), (2 TGGCCA 4/6), (25 TATAAT 1) |
| Spacer 14 / 22 | empty |
