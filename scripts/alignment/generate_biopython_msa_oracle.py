#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-27b, B09 F126) oracle: Biopython 1.88 Bio.Align writers and readers of the multiple-alignment text formats fasta, a2m, clustal, phylip, nexus, stockholm.

Seqeron API: AlignmentFormats.Write(alignments, format, options) / AlignmentFormats.Parse(text, format) (Seqeron.Genomics.Alignment).
Reference (run here): Biopython 1.88 `Bio.Align.write(alignments, handle, fmt)` (bytes of the file) and `Bio.Align.parse(handle, fmt)` (the alignments read back: ids, descriptions, printed rows, annotations).
Output TSV (tests/.../TestData/Alignment/biopython_msa_oracle.tsv):  id  format  options  alignments_json  file  parsed_json
  alignments_json = [{records:[{id, desc, row, ann:{}, letter:{}, dbx:[]}], ann:{}, col:{}}]; options = JSON (interleave); file / parsed_json = escaped (\\n \\t \\\\).
  parsed_json = what Bio.Align.parse reads from `file`: the same structure (rows with all-gap columns removed; annotations as read).
Usage: python3 scripts/alignment/generate_biopython_msa_oracle.py OUTFILE [casesPerFormat] [seed]
"""
import io, json, random, sys
import Bio
from Bio import Align
from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

FORMATS = ['fasta', 'a2m', 'clustal', 'phylip', 'nexus', 'stockholm']


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


IDS = ['seq1', 'seq2', 'sp|P12345|ABC_HUMAN', 'gi|12|ref|NM_1.2|', 'name with space', 'x:y;z', "q'uote", 'A', 'Homo_sapiens_chr1_region_long_identifier', 'a(b)[c],d', 'tab_1', '12345', 'ID.copy', 'seq1']
DESCS = ['', '', 'alpha beta', 'hypothetical protein [Homo sapiens]', 'a description with several words', 'x']


def rand_ids(rnd, n, stockholm, friendly=False):
    if friendly:
        return ['s%d' % (i + 1) for i in range(n)] if rnd.random() < 0.5 else ['Seq_%d_%s' % (i, rnd.choice('xyz')) for i in range(n)]
    ids = []
    for i in range(n):
        s = rnd.choice(IDS) if rnd.random() < 0.7 else 'r%d' % i
        if stockholm: s = s.replace(' ', '_') or 'r%d' % i
        ids.append(s)
    if rnd.random() < 0.05 and not stockholm: ids[0] = ''
    return ids


def make_rows(rnd, n, length, alpha, gap_rate):
    rows = []
    for _ in range(n):
        row = []
        i = 0
        while i < length:
            if rnd.random() < gap_rate:
                k = rnd.randint(1, 6); row.extend('-' * k); i += k
            else:
                row.append(rnd.choice(alpha)); i += 1
        rows.append(''.join(row[:length]))
    return rows


def make_alignment(rnd, fmt, friendly=False):
    n = rnd.choice([1, 2, 2, 3, 3, 4, 5, 8])
    length = rnd.choice([1, 4, 11, 25, 49, 50, 51, 70, 71, 100, 140]) if not (fmt == 'nexus' and rnd.random() < 0.15) else rnd.choice([1001, 1100, 1500])
    mol = rnd.choice(['DNA', 'RNA', 'protein', 'DNA'])
    alpha = {'DNA': 'ACGTACGTN', 'RNA': 'ACGUACGUN', 'protein': 'ACDEFGHIKLMNPQRSTVWYX'}[mol]
    rows = make_rows(rnd, n, length, alpha, rnd.choice([0.0, 0.05, 0.15, 0.3]))
    # the printed-alignment parser drops columns that are gaps in every row: go through Bio like the readers do
    seqs, coords = Alignment.parse_printed_alignment([r.encode() for r in rows])
    ids = rand_ids(rnd, n, fmt == 'stockholm', friendly)
    recs = []
    for i in range(n):
        ann = {'molecule_type': mol} if fmt == 'nexus' else {}
        rec = SeqRecord(Seq(seqs[i]), id=ids[i], description=rnd.choice(DESCS) if fmt != 'phylip' else '', annotations=ann)
        recs.append(rec)
    al = Alignment(recs, coords)
    shape = al.shape
    cols = shape[1]
    if fmt == 'stockholm':
        # the Stockholm writer needs ids without a duplicate-free requirement; keep as is
        ann = {}
        if rnd.random() < 0.7: ann['identifier'] = rnd.choice(['PF00001', 'RF00010', 'my_family'])
        if rnd.random() < 0.5: ann['accession'] = rnd.choice(['PF00001.12', 'RF00010.3'])
        if rnd.random() < 0.5: ann['definition'] = rnd.choice(['7 transmembrane receptor', 'Some family', 'Long definition text ' * 3])
        if rnd.random() < 0.4: ann['author'] = rnd.choice([['Sonnhammer ELL'], ['Griffiths-Jones SR', 'Bateman A'], ['A', 'B', 'C']])
        if rnd.random() < 0.3: ann['type'] = rnd.choice(['Family', 'Domain', 'Motif'])
        if rnd.random() < 0.3: ann['trusted cutoff'] = '25.00 25.00'
        if rnd.random() < 0.3: ann['source of seed'] = 'Prosite'
        if rnd.random() < 0.3: ann['comment'] = rnd.choice(['A short comment.', 'A much longer comment text that has to be wrapped over several lines by the writer because it is longer than seventy-nine characters in total, truly.', 'Hyphen-ated words and ---- dashes, well-known two-line wrapped text for the fill of Python textwrap: pre-defined, out-of-line.'])
        if rnd.random() < 0.2: ann['database references'] = [{'reference': 'PDB; 1abc A; 1-100;', 'comment': 'x'}, {'reference': 'SO; 0000000; family'}]
        if rnd.random() < 0.2: ann['nested domains'] = [{'accession': 'PF00002', 'location': 'seq1/1-50'}]
        if rnd.random() < 0.15: ann['references'] = [{'number': 1, 'medline': '12345678', 'title': 'A title of a paper', 'author': 'Doe J, Roe R', 'location': 'J Mol Biol 1:2-3.', 'comment': 'Some reference comment'}]
        if rnd.random() < 0.15: ann['wikipedia'] = ['Some_article', 'Other']
        al.annotations = ann
        col = {}
        if rnd.random() < 0.5: col['consensus sequence'] = ''.join(rnd.choice('ACGT.') for _ in range(cols))
        if rnd.random() < 0.3: col['reference coordinate annotation'] = ''.join(rnd.choice('x.') for _ in range(cols))
        if rnd.random() < 0.2: col['consensus secondary structure'] = ''.join(rnd.choice('HEC.') for _ in range(cols))
        al.column_annotations = col
        for i, rec in enumerate(al.sequences):
            if rnd.random() < 0.4: rec.annotations['organism'] = rnd.choice(['Homo sapiens', 'Mus musculus'])
            if rnd.random() < 0.3: rec.annotations['accession'] = 'P%05d' % rnd.randint(0, 99999)
            if rnd.random() < 0.3: rec.dbxrefs = ['PDB; 1abc A;', 'SP; Q9 x;'][:rnd.randint(1, 2)]
            if rnd.random() < 0.25:
                n_res = len(rec.seq)
                rec.letter_annotations['secondary structure'] = ''.join(rnd.choice('HEC') for _ in range(n_res))
    elif fmt == 'a2m':
        al.column_annotations = {'state': ''.join(rnd.choice('DDDI') for _ in range(cols))}
    elif fmt == 'clustal':
        if rnd.random() < 0.6:
            al.column_annotations = {'clustal_consensus': ''.join(rnd.choice('*:. ') for _ in range(cols))}
    return al


def to_json(al, include_ann=True):
    recs = []
    for r in al.sequences:
        recs.append({'id': r.id, 'desc': r.description, 'row': None, 'ann': dict(r.annotations), 'letter': dict(r.letter_annotations), 'dbx': list(r.dbxrefs)})
    for i, row in enumerate(al):
        recs[i]['row'] = row
    out = {'records': recs}
    out['ann'] = dict(getattr(al, 'annotations', {}) or {})
    out['col'] = dict(getattr(al, 'column_annotations', {}) or {})
    return out


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    for fmt in FORMATS:
        made = 0
        while made < per:
            count = 1 if fmt in ('a2m', 'nexus') else rnd.choice([1, 1, 1, 2, 3])
            friendly = rnd.random() < 0.75
            als = [make_alignment(rnd, fmt, friendly) for _ in range(count)]
            options = {}
            if fmt == 'nexus': options['interleave'] = rnd.choice([None, None, True, False])
            if fmt == 'clustal' and count > 1 and False: pass
            buf = io.StringIO()
            try:
                if fmt == 'nexus':
                    Align.write(als, buf, fmt, interleave=options['interleave'])
                else:
                    Align.write(als, buf, fmt)
            except Exception as e:
                continue
            text = buf.getvalue()
            try:
                parsed = [to_json(a) for a in Align.parse(io.StringIO(text), fmt)]
            except Exception as e:
                parsed = None
            # Biopython's printed-alignment parser reads uninitialised memory when the rows of one file have different lengths (several alignments of different width in one fasta / clustal file): no expectation then
            if parsed is not None and fmt in ('fasta', 'a2m', 'clustal') and len({len(a_) for al_ in als for a_ in al_}) > 1:
                parsed = None
            ok = parsed is not None
            rows.append('\t'.join([str(len(rows)), fmt, json.dumps(options), json.dumps([to_json(a) for a in als]), esc(text), esc(json.dumps(parsed)) if ok else 'null']))
            made += 1
    with open(outfile, 'w') as f:
        f.write('id\tformat\toptions\talignments\tfile\tparsed\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 120, int(sys.argv[3]) if len(sys.argv) > 3 else 20261126)
