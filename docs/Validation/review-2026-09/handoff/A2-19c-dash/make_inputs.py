#!/usr/bin/env python3
"""Deterministic inputs of the hand-off A2-19c-dash: well-known protein sequences that have homologs in the PDB (so that the DASH server finds templates), the sets of `sets.tsv` and the command lines.
Written to ./inputs (sets.tsv: id TAB options TAB fasta-file; fasta files).  Usage: python3 make_inputs.py [outdir]"""
import os, random, sys

PROT = {
    'ubiquitin': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
    'hen_lysozyme': 'KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL',
    'whale_myoglobin': 'MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG',
    'horse_cytc': 'GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE',
    'bpti': 'RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCGGA',
    'gfp': 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK',
    'calmodulin': 'ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK',
    't4_lysozyme': 'MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAAINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL',
    'ecoli_trx': 'SDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA',
    'barnase': 'AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR',
    'rnase_a': 'KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV',
    'hba': 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR',
    'hbb': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
    'snase': 'ATSTKKLHKEPATLIKAIDGDTVKLMYKGQPMTFRLLLVDTPETKHPKKGVEKYGPEASAFTKKMVENAKKIEVEFDKGQRTDKYGRGLAYIYADGKMVNEALVRQGLAKVAYVYKPNNTHEQHLRKSEAQAKKEKLNIWSEDNADSGQ',
}
AA = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(rng, s, rate, ind=0.02):
    out = []
    for c in s:
        x = rng.random()
        if x < ind / 2:
            continue
        if x < ind:
            out.append(c)
            out.append(rng.choice(AA))
            continue
        out.append(rng.choice(AA) if rng.random() < rate else c)
    return ''.join(out)


def fasta(items):
    return ''.join('>%s\n%s\n' % (n, s) for n, s in items)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else 'inputs'
    os.makedirs(out, exist_ok=True)
    rng = random.Random(7450515)
    names = sorted(PROT)
    sets = []
    for k, n in enumerate(names):
        sets.append(('single-%s' % n, [(n, PROT[n])]))
    for k in range(10):
        a = rng.sample(names, rng.randint(2, 5))
        sets.append(('mix-%02d' % k, [(n, PROT[n]) for n in a]))
    for k in range(8):
        n = rng.choice(names)
        fam = [(n + '_wt', PROT[n])] + [('%s_v%d' % (n, i), mutate(rng, PROT[n], rng.choice([0.05, 0.15, 0.3]))) for i in range(rng.randint(2, 5))]
        sets.append(('family-%02d' % k, fam))
    sets.append(('globin-pair', [('hba', PROT['hba']), ('hbb', PROT['hbb']), ('mb', PROT['whale_myoglobin'])]))
    sets.append(('lysozymes', [('hen', PROT['hen_lysozyme']), ('t4', PROT['t4_lysozyme']), ('t4_mut', mutate(rng, PROT['t4_lysozyme'], 0.1))]))
    sets.append(('nucleotide', [('dna1', 'ACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATGCATGCATGCATCGATCGTAGCTAGCTAGC'), ('dna2', 'ACGTACGTACGAAGCTAGCTAGCTAGCATCGATCGTTCGATGCATGCATGCATCGATCGTAGCTAGCTAGC'), ('dna3', 'ACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATGCATGCATGCATCGATCGTAGCTAGCTAGA')]))
    sets.append(('unknown-letters', [('u1', PROT['ubiquitin'][:30] + 'XXUU' + PROT['ubiquitin'][30:]), ('u2', PROT['ubiquitin'])]))
    sets.append(('no-homolog', [('r1', ''.join(rng.choice(AA) for _ in range(120))), ('r2', ''.join(rng.choice(AA) for _ in range(90))), ('r3', ''.join(rng.choice(AA) for _ in range(150)))]))
    rows = []
    for sid, items in sets:
        fn = sid + '.fa'
        open(os.path.join(out, fn), 'w').write(fasta(items))
        rows.append((sid, fn))
    # the command lines (options are applied to every set of the group)
    opts = [
        ('dash', '--dash'), ('dash-quiet', '--dash --quiet'), ('dash-retree1', '--dash --retree 1'), ('dash-localpair', '--dash --localpair'), ('dash-globalpair', '--dash --globalpair'),
        ('dash-maxiterate0', '--dash --maxiterate 0'), ('dash-maxiterate10', '--dash --maxiterate 10'), ('dash-anysymbol', '--dash --anysymbol'), ('dash-originalseqonly', '--dash --originalseqonly'),
        ('dash-excludehomologs', '--dash --excludehomologs'), ('dash-reorder-clustal', '--dash --reorder --clustalout'), ('dash-treeout', '--dash --treeout --distout'), ('dash-op-ep', '--dash --op 2 --ep 0.2'),
        ('dash-server', '--dash --dashserver @SERVER@'), ('dash-server-noslash', '--dash --dashserver @SERVERNOSLASH@'), ('dash-merge', '--dash --merge @MERGE@'), ('dash-amino', '--dash --amino'), ('dash-nuc', '--dash --nuc'),
    ]
    with open(os.path.join(out, 'sets.tsv'), 'w') as f:
        for sid, fn in rows:
            f.write('%s\t%s\n' % (sid, fn))
    with open(os.path.join(out, 'options.tsv'), 'w') as f:
        for oid, o in opts:
            f.write('%s\t%s\n' % (oid, o))
    print(len(rows), 'sets,', len(opts), 'option lines')


if __name__ == '__main__':
    main()
