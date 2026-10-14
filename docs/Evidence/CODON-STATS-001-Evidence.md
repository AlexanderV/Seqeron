# Evidence Artifact: CODON-STATS-001

**Test Unit ID:** CODON-STATS-001
**Algorithm:** Codon Usage Statistics (GetStatistics: codon counts, RSCU, ENC, GC1/GC2/GC3, GC3s, total codons; CalculateCai; E. coli / human reference tables)
**Date Collected:** 2026-06-13

---

## Online Sources

### Sharp & Li (1987) — Codon Adaptation Index (primary paper)

**URL:** https://doi.org/10.1093/nar/15.3.1281 (PubMed: https://pubmed.ncbi.nlm.nih.gov/3547335/)
**Accessed:** 2026-06-13 (retrieved via WebSearch "Sharp Li 1987 codon adaptation index CAI formula"; the full-text PDF at academic.oup.com/nar/article-pdf/15/3/1281 returns a tokened redirect and could not be opened directly, so the w-value table was taken from the Biopython reproduction below and the formula from Wikipedia/seqinr/CodonW, all retrieved this session)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **CAI definition:** CAI is the geometric mean of the relative-adaptiveness (w) values of the codons of a gene. Range 0–1 (1 = uses only the most-frequent synonym of every amino acid).
2. **Relative adaptiveness:** `w_i = f_i / f_max`, the ratio of the usage of codon i to the usage of the most-abundant synonymous codon for the same amino acid.
3. **Reference set:** w values for E. coli were derived from 27 very highly expressed genes (ribosomal proteins, outer-membrane proteins, elongation factors); reported in Table 1 ("Values of RSCU and w for codons in very highly expressed genes from E.coli and yeast").

### Wikipedia — Codon Adaptation Index (formula transcription, cites Sharp & Li 1987)

**URL:** https://en.wikipedia.org/wiki/Codon_Adaptation_Index
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 4 (Wikipedia citing the primary Sharp & Li 1987)

**Key Extracted Points:**

1. **Weight formula (verbatim):** `w_i = f_i / max(f_j)` where `max(f_j)` is "the frequency of the most frequent synonymous codon for that amino acid."
2. **CAI formula (verbatim):** `CAI = (∏ w_i)^(1/L)`, equivalently `CAI = exp[(1/L) Σ ln(w_i)]`, with L the number of codons.

### seqinr `cai` reference documentation (R reference implementation)

**URL:** https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation, seqinr R package)

**Key Extracted Points:**

1. **w:** "The relative adaptiveness (w) of each codon is the ratio of the usage of each codon, to that of the most abundant codon for the same amino acid."
2. **CAI:** "The CAI index is defined as the geometric mean of these relative adaptiveness values"; computed via natural-log summation.
3. **Excluded codons:** "Non-synonymous codons and termination codons are excluded from the calculation (genetic code dependent)" — i.e. single-codon amino acids Met (ATG) and Trp (TGG) plus the stop codons.

### CodonW codon-usage indices documentation (Peden) — CAI and GC3s

**URL:** https://codonw.sourceforge.net/Indices.html
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation, CodonW)

**Key Extracted Points:**

1. **CAI exclusions (verbatim):** "Non-synonymous codons and termination codons (dependent on genetic code) are excluded."

### Peden (1999) — *Analysis of Codon Usage* thesis (definitive CodonW reference)

**URL:** https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf
**Accessed:** 2026-06-13 (WebFetch returned the binary PDF; extracted with `pdftotext`. Definition located at §1.8.2.1.3, line 2063, and the user-manual option 6, line 3257)
**Authority rank:** 1–3 (M.Sc./Ph.D. thesis defining the CodonW indices)

**Key Extracted Points:**

1. **GC3s definition (verbatim, §1.8.2.1.3):** "The index GC3s, is the frequency of G or C nucleotides present at the third position of synonymous codons (i.e. excluding Met, Trp and termination codons)."
2. **GC3s computation (verbatim, manual option 6):** "This option calculates the fraction of codons synonymous at the third codon position, and having guanine or cytosine at that third position."
3. **ENc range / clamp (verbatim):** "20 (when effectively only a single codon is used for each amino-acid) and 61 (when codons are used randomly). If ENc is greater than 61 ... it is corrected to 61."

### EMBOSS `cusp` documentation — GC position statistics

**URL:** https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation, EMBOSS)

**Key Extracted Points:**

1. **GC fields (verbatim labels):** "#Coding GC" (overall), "#1st letter GC" (codon position 1), "#2nd letter GC" (position 2), "#3rd letter GC" (position 3).
2. **Fraction field:** "the proportion of usage of the codon among its redundant set" (used to confirm RSCU-style per-amino-acid normalization).
3. **Source code (review 2026-09, EMBOSS 6.6 `ajcod.c`, opened locally):** `ajCodSetTripletsS` skips triplets that are not A/C/G/T(U) without shifting the frame and decrements the codon count; `ajCodWrite` sums G/C per position over all 64 codons (`AJCODSTART` = 64, i.e. **termination codons included**) and divides by `CodonCount`. "Coding GC" = mean of the three positions.

### CodonW 1.4.4 source + binary (review 2026-09)

**Opened:** `codon_us.c` (`gc_out`, `how_synon`, `base_sil_us_out`) and `README_indices.txt` of CodonW 1.4.4, compiled locally; Peden 1999 thesis PDF text-extracted (§1.8.2.1.3, verbatim quote above).

1. `how_synon`: a codon's degeneracy is the number of codons encoding the same amino acid **in the selected genetic code** (`-code`); `gc_out` skips stop codons (`ca == 11`) and codons with degeneracy 1 when accumulating GC3s/L_sym. So GC3s is genetic-code dependent (NCBI 2: Met = {ATA, ATG} and Trp = {TGA, TGG} are synonymous; AGA/AGG are stops).
2. CodonW's `GC` and GC1/GC2/GC3 (`-base`) exclude stop codons (denominator L_aa); EMBOSS cusp includes them. Seqeron's Gc1/Gc2/Gc3/OverallGc follow cusp by default (identical to CodonW for stop-free genes); `GetStatistics(…, GeneticCode, CodonGcConvention.CodonW)` (review 2026-10-14, F37) follows CodonW. `gc_out` source (verbatim): `if (pcu->ca[id] == 11) continue; /* skip if a stop codon */` before accumulating `base_1/2/3`, `base_tot` and `totalaa`; GC = `(base_tot[2]+base_tot[4])/(totalaa*3)`; GC1 = `(base_1[2]+base_1[4])/totalaa`; GC3s = `(bases[2]+bases[4])/tot_s` over codons with `ds > 1`; no output when `tot_s == 0 || totalaa == 0` ("appear to be too short"). Stops are those of `-code` (`pcu->ca`).
3. Cross-check (review 2026-09): 708 genes (random DNA/RNA/lower-case/N/R/Y, codon-sampled, few-codon) × 8 CodonW codes (NCBI 1,2,3,4,5,6,9,10) = 5664 comparisons — a Python port of `gc_out` equals CodonW's GC3s/GC/L_sym/L_aa/GC1-3 (3 dp) on every gene with output, and C# `GetStatistics(…, GeneticCode)` equals the port (GC3s) and the cusp formula (TotalCodons, GC1-3, OverallGc) to 1e-9 — 0 mismatches. Biopython 1.88 `GC123` equals (OverallGc, Gc1, Gc2, Gc3) on all 450 clean ACGT whole-codon genes.
4. CodonW-convention cross-check (review 2026-10-14, F37): CodonW 1.4.4 rebuilt from Ubuntu `codonw_1.4.4.orig.tar.gz`, run as `codonw genes.fa -base -enc -gc -gc3s -L_aa -code c -silent -nomenu -nowarn -machine` (c = 0..7) on 304 genes (seed 20261014; terminal/internal stops, N/R/Y triplets, lower-case, RNA, trailing bases). Every gene × code with CodonW output: C# `CodonW` convention GC (OverallGc), GC1, GC2, GC3 and GC3s equal CodonW at 3 dp — 0 mismatches (2220 of the 2432 cases differ from the cusp default because of stops). Locked genes: `ATGAAACTGTTCGCCAAATAA` (GC 0.389, GC1/2/3 0.333/0.167/0.667, GC3s 0.600), `ATGGCTTAAGCCNNNGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA` (code 0: GC 0.600, GC1/2/3 0.600/0.500/0.700, GC3s 0.625; code 1: 0.576, 0.545/0.545/0.636, 0.636), random r53 (code 0: GC 0.521, GC1/2/3 0.562/0.521/0.479; code 1: 0.529, 0.587/0.500/0.500).
5. EMBOSS 6.6.0.0 `cusp` binary re-verification (2026-10-14): `ATGGCTTAAGCCGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA` → Coding GC 47.62 %, 1st 42.86 %, 2nd 42.86 %, 3rd 57.14 % (14 codons, stops included) = Seqeron default. Note: the binary counts an ambiguous triplet as a codon (`NNN` listed as `AAA`, so `ATGGCAGCCTAAtggNNNgcg` gives 1st-letter GC 42.86 % over 7 codons); this `ajCodBase` artefact contradicts the documented `ajCodSetTripletsS` contract and is not reproduced (ambiguous triplets skipped, CODON-USAGE-001 F9).

### Biopython `SharpEcoliIndex` (reference implementation of Sharp & Li 1987 w values)

**URL:** https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py
**Accessed:** 2026-06-13 (fetched via `curl` from the `biopython-179` tag; the file was removed from `master`)
**Authority rank:** 3 (reference implementation, Biopython v1.79)

**Key Extracted Points:**

1. **Provenance (file header, verbatim):** "Codon adaption indxes, including Sharp and Li (1987) E. coli index ... from Sharp & Li, Nucleic Acids Res. 1987."
2. **Sample w values:** `GCT=1, GCC=0.122, GCA=0.586, GCG=0.424` (Ala); `CTG=1, CTC=0.037, CTA=0.007, CTT=0.042, TTA=0.02, TTG=0.02` (Leu); `CGT=1, CGC=0.356, CGA=0.004, CGG=0.004, AGA=0.004, AGG=0.002` (Arg); `TTT=0.296, TTC=1` (Phe); `ATG=1`, `TGG=1`. (Full 64-codon table transcribed into `EColiOptimalCodons`; stop codons added as 0.0.)

### Kazusa Codon Usage Database — *Homo sapiens* [gbpri]

**URL:** https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 5 (curated database; underlying method Nakamura et al. 2000)

**Key Extracted Points:**

1. **Dataset scale:** 93,487 coding sequences, 40,662,582 codons.
2. **Sample per-thousand frequencies:** CTG 39.6, GTG 28.1, GCC 27.7, AAG 31.9, CAG 34.2, GAG 39.6, ATG 22.0, TGG 13.2. (Full 64-codon per-thousand table used to derive `HumanOptimalCodons` RSCU; see Dataset below.)

### CodonW GC3s + ENc plot reference (Wright 1990 relationship, viral codon study)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7596632/
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 4 (peer-reviewed paper applying CodonW; confirms exclusion set)

**Key Extracted Points:**

1. **59-codon synonymous set (verbatim):** RSCU/GC3s are computed over "the 59 synonymous codons (excluding the codons of AUG, UGG and the three stop codons)" — confirms the Met/Trp/stop exclusion in Peden's GC3s definition.

---

### CodonW 1.4.4 Fop and CBI — source + binary (review 2026-10, audit A3-13)

**Source:** Ubuntu `codonw_1.4.4.orig.tar.gz` (http://archive.ubuntu.com/ubuntu/pool/universe/c/codonw/), compiled locally; `codon_us.c` `fop_out`, `cbi_out`, `how_synon`, `how_synon_aa`, `initilize_point`; `codonW.h` `FOP_STRUCT fop[]`; `commline.c`.
**Accessed:** 2026-10-14

**Key extracted points:**

1. `has_opt_info`: for codons x = 1..64, `if (pcu->ca[x] == 11 || *(ds + x) == 1) continue;` (stop or single-codon amino acid), `if (fop_cod[x] == 3) has_opt_info[ca[x]]++` — only amino acids with an optimal codon are scored; for the modified Fop a rare codon (value 1) also marks its amino acid.
2. Fop: over codons of scored amino acids `case 3: opt += n; case 2: std += n; case 1: nonopt += n;` → `ffop = opt / (opt + nonopt + std)`; modified (`factor_in_rare`, asked interactively only when not `-silent`): `(opt - nonopt) / (opt + nonopt + std)` — no clamp in the code although the header comment says negative values are adjusted to zero; total 0 → `ffop = 0.0`. Printed `%5.3f`.
3. CBI: `case 3: opt += n; tot_cod += n; exp_cod += naa[ca[x]] / da[ca[x]]; case 2: case 1: tot_cod += n;` → `fcbi = (opt - exp_cod) / (tot_cod - exp_cod)`, `if (tot_cod - exp_cod)` else 0.0 (rare codons count as common). `da` = codons per amino acid of the selected code (`how_synon_aa`).
4. `initilize_point`: `pfop = &fop[fop_species]; pcbi = &fop[fop_species];` — the same `-f_type` set feeds Fop and CBI. User sets: `-fop_file` / `-cbi_file`, 64 digits 1/2/3 in CodonW codon order (UUU, UCU, UAU, UGU, UUC, …).
5. Built-in sets (`fop[]`, `-f_type 0..7`): E. coli (Ikemura 1985, default), B. subtilis (Sharp et al. 1990; the only set with rare codons), D. discoideum (Sharp & Devine 1989), A. nidulans (Lloyd & Sharp 1991), S. cerevisiae (Sharp & Cowe 1991), D. melanogaster (Shields et al. 1988), C. elegans (Stenico, Lloyd & Sharp 1994), N. crassa (Lloyd & Sharp 1993); 24/19/15/19/22/20/20/22 optimal codons.
6. CodonW `-code 0..7` = NCBI tables 1, 2, 3, 4, 5, 6, 10, 9 (codonW.h `cu[]`).

**Binary lock (TestData/Codon/codonw_fop_cbi.tsv):** 70 genes (20 uniform random, 10 E. coli-optimal-biased, 8 non-optimal-biased, 6 GC-rich, 6 AT-rich, 4 with all 64 codons incl. in-frame stops, 16 edge cases: Met/Trp only, stops only, single codon, N/R triplets, lower case, RNA, trailing bases) × (8 built-in sets + 3 user `-fop_file`/`-cbi_file` sets, one with rare codons and one listing only stop/single-codon and alt-code codons) × `-code 0..7`, plus modified Fop (stdin answer `y`) for B. subtilis and user1: 7 280 Fop and 6 160 CBI printed values. C# `TryCalculateFop` / `TryCalculateCbi`: all within ±0.0005 of the printed value; every not-calculable result corresponds to a printed 0.000; 190 modified-Fop values are negative.

<details><summary>Generator (python3 -I gen.py SCRATCH codonw OUT.tsv)</summary>

```python
# WP23 generator: genes + CodonW 1.4.4 -fop/-cbi runs -> TSV fixture
import random, subprocess, os, sys, itertools
SCR=sys.argv[1]; CW=sys.argv[2]; OUT=sys.argv[3]
os.makedirs(SCR,exist_ok=True)
rnd=random.Random(20261014)
B="TCAG"
codons=[a+b+c for a in B for b in B for c in B]
ecoli_opt=set("AAA AAC ACC ACT AGC ATC CAC CAG CCG CGC CGT CTG GAA GAC GCG GCT GGC GGT GTT TAC TCC TCT TGC TTC".split())
genes=[]
def rgene(n, weights=None, stop=True):
    cs=[c for c in codons if c not in ("TAA","TAG","TGA")]
    w=[weights(c) for c in cs] if weights else None
    s="ATG"+"".join(rnd.choices(cs,weights=w,k=n-1))
    return s+("TAA" if stop else "")
for i in range(20): genes.append((f"unif{i}", rgene(rnd.randint(20,400))))
for i in range(10): genes.append((f"ecopt{i}", rgene(rnd.randint(30,400), lambda c: 6 if c in ecoli_opt else 1)))
for i in range(8): genes.append((f"econon{i}", rgene(rnd.randint(30,300), lambda c: 1 if c in ecoli_opt else 5)))
for i in range(6): genes.append((f"gcrich{i}", rgene(rnd.randint(40,300), lambda c: 1+4*sum(x in "GC" for x in c))))
for i in range(6): genes.append((f"atrich{i}", rgene(rnd.randint(40,300), lambda c: 1+4*sum(x in "AT" for x in c))))
# all 64 codons incl. internal stops (alt codes read TGA/AGA etc.)
for i in range(4): genes.append((f"all64_{i}", "".join(rnd.choices(codons,k=rnd.randint(60,250)))))
# edge cases
genes += [("short1","ATGAAATAA"),("mettrp","ATGTGGATGTGGTAA"),("stoponly","TAATAGTGA"),
  ("one_opt","CTG"),("one_nonopt","CTA"),("phe_opt_only","TTCTTCTTC"),("lys_mixed","AAAAAGAAG"),
  ("amb","ATGNNNCTGCTRAAACTGNAGTTCTAA"),("lower","atgctgctaaaagggggcggttaa"),("rna","AUGCUGCUAAAAGGGGGCGGUUAA"),
  ("partial","ATGCTGAAAGG"),("atarich","ATAATAATGTGATGAAGAAGAAGATAATG"),("ctn","CTTCTCCTACTGACTACC"),
  ("leu6","CTTCTCCTACTGTTATTG"),("ile3","ATTATCATAATAATT"),("gly4","GGTGGCGGAGGGGGG")]
with open(f"{SCR}/genes.fa","w") as f:
    for n,s in genes: f.write(f">{n}\n{s}\n")
NAMES=["EscherichiaColi","BacillusSubtilis","DictyosteliumDiscoideum","AspergillusNidulans","SaccharomycesCerevisiae","DrosophilaMelanogaster","CaenorhabditisElegans","NeurosporaCrassa"]
NCBI=[1,2,3,4,5,6,10,9]
def idx(c):  # CodonW codon order (1-based)
    return B.index(c[0])*16+B.index(c[2])*4+B.index(c[1])+1
def run(args, tag, stdin=None):
    out=f"{SCR}/{tag}.out"; blk=f"{SCR}/{tag}.blk"
    for p in (out,blk):
        if os.path.exists(p): os.remove(p)
    subprocess.run([CW,f"{SCR}/genes.fa",out,blk]+args,input=stdin,capture_output=True,text=True,check=False,cwd=SCR)
    lines=open(out).read().split("\n")
    hdr=lines[0].split(); res={}
    for l in lines[1:]:
        if not l.strip(): continue
        p=l.split(); res[p[0]]=dict(zip(hdr[1:],p[1:]))
    return res
rows=[]
for ft in range(8):
    for cc in range(8):
        r=run(["-fop","-cbi","-f_type",str(ft),"-code",str(cc),"-nomenu","-silent","-nowarn","-machine"],f"r{ft}_{cc}")
        for n,_ in genes: rows.append((NAMES[ft],NCBI[cc],n,r[n]["Fop"],r[n]["CBI"],"0"))
# modified Fop (interactive 'y' answer): B. subtilis, all codes
for cc in range(8):
    r=run(["-fop","-f_type","1","-code",str(cc),"-nomenu","-nowarn","-machine"],f"m1_{cc}",stdin="y\n"*20)
    for n,_ in genes: rows.append((NAMES[1],NCBI[cc],n,r[n]["Fop"],"NA","1"))
# user-supplied sets via -fop_file/-cbi_file
users={"user1":(sorted("TTT CTA ATA GTA TCG CCC ACA GCA TAT CAT CAA AAT AAG GAT GAG TGT AGA GGG TGA ATG TGG".split()),sorted("TTC CTG GCC GGC CGC AGG AAA".split())),
       "user2":(["CTG","GGC"],[]),"user3":(["TGA","AGA","ATA","ATG","TGG"],[])}
U=[]
for un,(opt,rare) in users.items():
    v=[2]*65
    for c in opt: v[idx(c)]=3
    for c in rare: v[idx(c)]=1
    fn=f"{SCR}/{un}.fop"; open(fn,"w").write("\n".join(str(v[x]) for x in range(1,65))+"\n")
    U.append((un,opt,rare))
    for cc in range(8):
        r=run(["-fop","-cbi","-fop_file",fn,"-cbi_file",fn,"-code",str(cc),"-nomenu","-silent","-nowarn","-machine"],f"{un}_{cc}")
        for n,_ in genes: rows.append((un,NCBI[cc],n,r[n]["Fop"],r[n]["CBI"],"0"))
        if rare:
            r=run(["-fop","-fop_file",fn,"-code",str(cc),"-nomenu","-nowarn","-machine"],f"{un}m_{cc}",stdin="y\n"*20)
            for n,_ in genes: rows.append((un,NCBI[cc],n,r[n]["Fop"],"NA","1"))
with open(OUT,"w") as f:
    f.write("# CODON-STATS-001 (FIN-B02 WP23, audit A3-13): CodonW 1.4.4 (Ubuntu codonw_1.4.4.orig.tar.gz, compiled locally) Fop/CBI printed values (%5.3f).\n")
    f.write("# `codonw genes.fa out blk -fop -cbi -f_type F -code C -nomenu -silent -nowarn -machine`; user sets via -fop_file/-cbi_file (64 values 1/2/3 in CodonW codon order);\n")
    f.write("# modified Fop (opt-rare)/total = same run without -silent answering 'y' to CodonW's prompt (MOD=1, CBI=NA). CodonW -code 0..7 = NCBI 1,2,3,4,5,6,10,9.\n")
    f.write("# G<TAB>name<TAB>gene | U<TAB>set<TAB>optimal codons<TAB>rare codons | R<TAB>set<TAB>NCBI table<TAB>gene<TAB>Fop<TAB>CBI<TAB>MOD. Generator: docs/Evidence/CODON-STATS-001-Evidence.md.\n")
    for n,s in genes: f.write(f"G\t{n}\t{s}\n")
    for un,opt,rare in U: f.write(f"U\t{un}\t{','.join(opt)}\t{','.join(rare)}\n")
    for r in rows: f.write("R\t"+"\t".join(map(str,r))+"\n")
print(len(genes),"genes",len(rows),"rows")
```

</details>

## Documented Corner Cases and Failure Modes

### From seqinr / CodonW

1. **Single-codon amino acids and stop codons:** excluded from CAI (and from GC3s). A sequence containing only Met/Trp/stop codons has no scorable codons → CAI is undefined (this implementation returns 0).
2. **Zero-frequency codons:** seqinr/EMBOSS substitute a small value (0.01, Bulmer 1988) to avoid `ln(0)`. *(Resolved 2026-09, CODON-CAI-001 F13: the implementation now does the same — w < 0.0001 → 0.01, CodonW `cai_out` — instead of skipping such codons.)*

### From Peden thesis

1. **GC3s denominator empty:** a sequence with no synonymous codons (only Met/Trp/stop) has an empty GC3s denominator → GC3s reported as 0.

---

## Test Datasets

### Dataset: E. coli relative adaptiveness (Sharp & Li 1987 via Biopython SharpEcoliIndex)

**Source:** Sharp & Li (1987) NAR 15:1281–1295, Table 1; reproduced in Biopython v1.79 `SharpEcoliIndex`.

| Codon (AA) | w |
|-----------|----|
| CTG (Leu) | 1.000 |
| CTC (Leu) | 0.037 |
| CTA (Leu) | 0.007 |
| GCT (Ala) | 1.000 |
| GCC (Ala) | 0.122 |
| GCA (Ala) | 0.586 |
| CGT (Arg) | 1.000 |
| AGG (Arg) | 0.002 |
| TTC (Phe) | 1.000 |
| TTT (Phe) | 0.296 |

### Dataset: Human RSCU (derived from Kazusa H. sapiens per-thousand frequencies)

**Source:** Kazusa Codon Usage Database, Homo sapiens [gbpri] (accessed 2026-06-13). RSCU_j = n·x_j / Σ_k x_k.

| Codon (AA) | per-1000 | RSCU |
|-----------|----------|------|
| CTG (Leu, n=6) | 39.6 | 2.3713 |
| GCC (Ala, n=4) | 27.7 | 1.5988 |
| ATG (Met, n=1) | 22.0 | 1.0000 |
| TGG (Trp, n=1) | 13.2 | 1.0000 |

### Dataset: Worked CAI / GC examples (derived from the formulas above)

| Input (DNA) | Quantity | Expected | Derivation |
|-------------|----------|----------|------------|
| `CTGATCGTTGCTCGTAAA` | CAI vs E. coli w | 1.0 | all six codons have w=1 → geo mean 1 |
| `GCTGCC` | CAI vs E. coli w | 0.34928498393146 | √(1×0.122) (Ala GCT w=1, GCC w=0.122) |
| `CTAATAGTC` | CAI vs E. coli w | 0.01114947479545 | ∛(0.007×0.003×0.066) |
| `ATGTGGTAA` | CAI vs E. coli w | 0 | Met+Trp+stop only → no scorable codon |
| `GCCGCA` | GC3s | 50.0 | 2 Ala codons; 3rd = C(GC), A(not) → 1/2 |
| `ATGGCA` | GC3s | 0.0 | Met excluded; Ala GCA 3rd = A → 0/1 |
| `ATGGCA` | GC3 (all) | 50.0 | positions 3 = G(ATG), A → 1/2 (shows GC3 ≠ GC3s) |
| `CTGGTTAAA` | GC1 / GC2 / GC3 | 66.666… / 0.0 / 33.333… | per-position counts over 3 codons |

---

## Assumptions

1. **ASSUMPTION: GC3s reported as a percentage.** CodonW reports GC3s as a fraction in [0,1]; this implementation reports it as a percentage (×100) for consistency with the existing GC1/GC2/GC3 fields (which follow EMBOSS cusp percentage style). Non-correctness-affecting unit/labeling choice; documented in the algorithm doc. The synonymous-codon *subset* used in the numerator/denominator is exactly per Peden.
2. **~~ASSUMPTION~~ (resolved 2026-09, CODON-CAI-001 F13 — now floored to 0.01 as CodonW): zero-w codons were skipped rather than floored to 0.01.** Sharp & Li / Bulmer floor missing codons to 0.01; this implementation skips codons whose relative adaptiveness is 0. For the supplied reference tables no synonymous codon has w=0, so CAI on real CDS is unaffected; only a gene using a codon entirely absent from the reference differs. Documented as a deviation.

---

## Recommendations for Test Coverage

1. **MUST Test:** CAI of an all-optimal E. coli sequence equals 1.0 — Evidence: Sharp & Li 1987 (w=1 for the preferred synonym).
2. **MUST Test:** CAI equals the geometric mean of w for a small hand-derived sequence (`GCTGCC` → √0.122) — Evidence: Sharp & Li 1987 formula.
3. **MUST Test:** CAI excludes Met, Trp and stop codons (`ATGTGGTAA` → 0) — Evidence: seqinr / CodonW exclusion rule.
4. **MUST Test:** GC3s counts only synonymous third positions, excluding Met/Trp/stop (`ATGGCA` GC3s=0 while GC3=50) — Evidence: Peden §1.8.2.1.3.
5. **MUST Test:** GC1/GC2/GC3 per-position percentages (`CTGGTTAAA`) — Evidence: EMBOSS cusp position GC.
6. **MUST Test:** `EColiOptimalCodons` reproduces Sharp & Li 1987 w values (CTG=1, GCC=0.122, …) — Evidence: Biopython SharpEcoliIndex.
7. **MUST Test:** `HumanOptimalCodons` reproduces Kazusa-derived RSCU (CTG≈2.3713, GCC≈1.5988, ATG=1) — Evidence: Kazusa H. sapiens.
8. **SHOULD Test:** null `DnaSequence` / null reference throw `ArgumentNullException`; empty string returns zeroed statistics / CAI 0 — Rationale: documented input contract.
9. **COULD Test:** `OverallGc` equals (GC1+GC2+GC3)/3 — Rationale: derived convenience property.

---

## References

1. Sharp PM, Li W-H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15(3):1281–1295. https://doi.org/10.1093/nar/15.3.1281
2. Peden JF. (1999). Analysis of Codon Usage (Ph.D. thesis, University of Nottingham; CodonW reference). https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf
3. CodonW — Codon usage indices. https://codonw.sourceforge.net/Indices.html
4. Charif D, Lobry JR. seqinr `cai` function documentation. https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html
5. EMBOSS `cusp` application documentation. https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp
6. Biopython v1.79, Bio.SeqUtils.CodonUsageIndices (`SharpEcoliIndex`). https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py
7. Nakamura Y, Gojobori T, Ikemura T. (2000). Codon usage tabulated from international DNA sequence databases. Nucleic Acids Research 28(1):292. Kazusa Codon Usage Database, Homo sapiens [gbpri]. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
8. Sharp PM, Tuohy TMF, Mosurski KR. (1986). Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Research 14(13):5125–5143. https://doi.org/10.1093/nar/14.13.5125
9. Wright F. (1990). The 'effective number of codons' used in a gene. Gene 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9
10. Ikemura T. (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. J. Mol. Biol. 151(3):389–409. https://doi.org/10.1016/0022-2836(81)90003-6
11. Bennetzen JL, Hall BD. (1982). Codon selection in yeast. J. Biol. Chem. 257(6):3026–3031. https://doi.org/10.1016/S0021-9258(19)81068-2

---

## Change History

- **2026-06-13**: Initial documentation (CODON-STATS-001).
- **2026-09-28**: Review 2026-09 (CODON-STATS-001 F18): CodonW/EMBOSS source opened; GC3s is genetic-code dependent (`GetStatistics(…, GeneticCode)`); cusp stop-inclusion recorded; CodonW/Biopython cross-check added.
- **2026-10-14**: Review 2026-09 (B02 F36/F37, audit A1-11/A1-14): `CodonGcConvention.CodonW` option (stops excluded, CodonW `gc_out`) verified on the CodonW 1.4.4 binary; cusp default re-verified on the EMBOSS 6.6 binary; `IsEncCalculable`.
- **2026-10-14**: Review 2026-10 (B02 FIN WP23, audit A3-13): CodonW 1.4.4 Fop and CBI (`fop_out`/`cbi_out`) implemented and locked on the binary (7 280 Fop + 6 160 CBI values, 8 built-in + 3 user sets × codes 0–7, incl. modified Fop).
