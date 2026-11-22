#!/usr/bin/env python3
"""A5-8c3h (B09 F261, WP58) oracle: EMBOSS 6.6.0 `union -feature` (emboss/union.c: union_CopyFeatures, the source features, the feature file)
    union -auto -sequence in.dat [-sformat1 F] -feature [-source] [-findoverlap] [-overlapfile ov.txt] -outseq o.fa [-osformat2 F] [-offormat2 FF] [-ofname2 NAME]
Inputs: EMBL entries made here (fragments of one sequence that overlap by random lengths, with CDS / exon / mRNA / gene features with joins, complements, orders; sometimes a feature that is itself a source feature with an origid tag), the entries
with feature tables of the reader oracle (embl, genbank, gff3, gff2, swiss, nbrf, refseqp), damaged files. JSONL.GZ rows: id, fmt, how, text, args, epoch, rc, files (every file the program wrote), err, hang, memerr
Usage: python3 scripts/alignment/generate_union_feature_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0, valgrind, gcc)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402
import generate_coderet_oracle as C  # noqa: E402

FMTS = ['own', 'own', 'own', 'own', 'own', 'embl', 'genbank', 'gff3', 'gff2', 'swiss', 'nbrf', 'refseqp', 'own']
FFORMATS = ['gff3', 'gff3', 'gff', 'gff2', 'embl', 'genbank', 'swiss', 'pir', 'tagseq', 'debug']
OSFORMATS = ['fasta', 'fasta', 'fasta', 'fasta', 'gcg', 'raw', 'text', 'ig', 'codata', 'fitch', 'mase', 'strider', 'gde', 'acedb', 'unknown', 'embl', 'em', 'emblnew', 'emold', 'emblold', 'genbank', 'gb', 'ddbj', 'refseq', 'genpept', 'refseqp', 'swiss', 'sw', 'swissprot', 'uniprot', 'swissnew', 'swold', 'fastq', 'fastq-sanger', 'fastq-illumina', 'fastq-solexa', 'sam', 'gff', 'gff2', 'gff3', 'scf', 'msf', 'clustal', 'aln', 'selex', 'phylip', 'phylipnon', 'phylip3', 'hennig86', 'mega', 'meganon', 'nexus', 'paup', 'nexusnon', 'paupnon', 'jackknifer', 'jackknifernon', 'treecon']


def fragments(rnd):
    n = rnd.randint(2, 4)
    total = ''.join(rnd.choice('acgt') for _ in range(rnd.randint(60, 200)))
    cuts = sorted(rnd.sample(range(10, len(total) - 10), n - 1)) if len(total) > 25 else []
    parts, start = [], 0
    for c in cuts + [len(total)]:
        ov = rnd.choice([0, 0, rnd.randint(1, 8)]) if parts else 0
        parts.append(total[max(0, start - ov):c]); start = c
    out = ''
    for i, seq in enumerate(parts):
        name = rnd.choice(['Y%d' % (i + 1), 'frag%d' % (i + 1), 'AB%06d' % rnd.randint(1, 999999)])
        entry = C.embl_entry(rnd, name)
        # replace the sequence by the fragment: keep the header and features, rewrite SQ
        head = entry[:entry.index('SQ   ')]
        n_ = len(seq)
        head = re.sub(r'(\d+) BP\.', '%d BP.' % n_, head, count=1)
        sq = ['SQ   Sequence %d BP; %d A; %d C; %d G; %d T; 0 other;' % (n_, seq.count('a'), seq.count('c'), seq.count('g'), seq.count('t'))]
        for j in range(0, n_, 60):
            chunk = seq[j:j + 60]
            sq.append('     %-66s%9d' % (' '.join(chunk[k:k + 10] for k in range(0, len(chunk), 10)), min(n_, j + 60)))
        # keep only the features that fit in the fragment
        lines = head.split('\n')
        keep, skip = [], False
        for ln in lines:
            if ln.startswith('FT   ') and not ln.startswith('FT                   '):
                nums = [int(x) for x in re.findall(r'\d+', ln[21:])]
                skip = bool(nums) and max(nums) > n_
            if ln.startswith('FT') and skip: continue
            keep.append(ln)
        out += '\n'.join(keep) + '\n' + '\n'.join(sq) + '\n//\n'
    return out


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            fmt = FMTS[len(rows) % len(FMTS)]
            if fmt == 'own':
                text, how = fragments(rnd), 'fragments'
                fmt = 'embl'
            else:
                text, how = S.make_input(rnd, gen, fmt)
            if rnd.random() < 0.06: text, how = S.damaged(rnd, text), how + '+damaged'
            args = ['-sformat1', fmt] if rnd.random() < 0.4 else []
            args.append('-feature')
            if rnd.random() < 0.6: args.append('-source')
            if rnd.random() < 0.5: args.append('-findoverlap')
            if rnd.random() < 0.3: args += ['-overlapfile', 'ov.txt']
            outseq = rnd.choice(['o.fa', 'o.fa', 'q.seq', 'plain', 'a.b.c', 'fasta::z.w'])
            args += ['-outseq', outseq]
            if rnd.random() < 0.3: args += ['-osformat2', rnd.choice(OSFORMATS)]
            if rnd.random() < 0.55: args += ['-offormat2', rnd.choice(FFORMATS)]
            if rnd.random() < 0.3: args += ['-ofname2', 'f.out']
            if rnd.random() < 0.04: args += ['-sbegin1', str(rnd.randint(1, 10))]
            if rnd.random() < 0.04: args += ['-send1', str(rnd.randint(20, 60))]
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            cmd = ['-auto', '-sequence', 'in.dat'] + args
            hang = False
            try:
                import resource
                def limit():
                    resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
                p = subprocess.run(['union'] + cmd, cwd=work, capture_output=True, timeout=5, preexec_fn=limit, env=env)
                rc, err = p.returncode, p.stderr.decode('latin-1')
                if rc == -25: hang = True
            except subprocess.TimeoutExpired:
                rc, err, hang = -9, '', True
                subprocess.run(['pkill', '-9', 'union'])
            files = {}
            if not hang:
                for nme in sorted(os.listdir(work)):
                    if nme == 'in.dat': continue
                    files[nme] = open(os.path.join(work, nme), encoding='latin-1', newline='').read()
            memerr = False
            if not hang and shutil.which('valgrind') and (rc != 0 or len(rows) % 4 == 0):
                for nme in os.listdir(work):
                    if nme != 'in.dat': os.remove(os.path.join(work, nme))
                try:
                    vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'union'] + cmd, cwd=work, capture_output=True, timeout=120)
                    memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                except subprocess.TimeoutExpired:
                    memerr = True
            if err.startswith('Concatenate multiple sequences into a single sequence\n'): err = err[len('Concatenate multiple sequences into a single sequence\n'):]
            rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'text': text, 'args': args, 'epoch': epoch, 'rc': rc, 'files': files, 'err': err, 'hang': hang, 'memerr': memerr})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'died' if r['rc'] else 'files' if r['files'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20265807)
