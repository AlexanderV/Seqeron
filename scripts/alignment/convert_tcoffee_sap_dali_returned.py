#!/usr/bin/env python3
"""B09 WP68 (F311): converts the returned runs of the Docker package docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr (runs 1, 3, 4, 5, 6: real SAP 1.1.3, DaliLite v5 on the PATH, DaliLite v3.3, t_coffee 13.41.0.28bdc39 with stderr recorded)
into the test data tests/.../TestData/Alignment/tcoffee_sap_dali_runs_oracle.tsv.gz. One row per (run, section, set):
  run  section(s1 = sap_real, s2 = -mode 3dcoffee with the real SAP, s3 = dali_pair with DaliLite v5 on the PATH, s4 = dali_pair with DaliLite v3.3)  idx  flags  seqs  rows  extra
  extra (compact JSON) = {"pdb": [...], "sap": {<pair file>: text}  (s1 only: the raw SAP output), "err": [normalised stderr lines of t_coffee that name the stop / the kill / the warning], "tm": {pair: TM-align output} (s3 / s4 sets with TMalign_pair),
                          "calls": [{"args", "rc", "out", "err"}]  (s4 only: every call of DaliLite)}
Usage: python3 -I convert_tcoffee_sap_dali_returned.py OUT.tsv.gz RUN1_DIR RUN3_DIR RUN4_DIR RUN5_DIR RUN6_DIR     (each = the unzipped archive directory that holds out/)"""
import csv, gzip, json, re, sys

csv.field_size_limit(10 ** 9)
KEEP = re.compile(r'ERROR|FATAL|COM:|TERMINATION STATUS|Killed|You have not provided|Could Not Be Found|SAP failed|WARNING', re.I)


def load(path):
    with open(path, newline='') as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        return [dict(zip(h, x)) for x in r]


def normalise(line):
    line = re.sub(r'/var/tmp/tco/\S+?/', '<tmp>/', line)
    line = re.sub(r'/tmp/\S+', '<path>', line)
    line = re.sub(r'<tmp>/\d+', '<tmp>/N', line)
    return re.sub(r'\b\d{2,}\b', 'N', line).strip()


def err_lines(text):
    out = []
    for l in text.split('\n'):
        if KEEP.search(l) and not l.startswith('pid') and 'MESSAGES' not in l:
            n = normalise(l)
            if n and n not in out:
                out.append(n)
    return out


def calls(extra):
    res = []
    for fn, t in sorted(extra.get('dali', {}).items()):
        head, _, rest = t.partition('--- STDOUT\n')
        out, _, err = rest.partition('--- STDERR\n')
        args = re.search(r'^ARGS: (.*)$', head, re.M).group(1)
        rc = int(re.search(r'^RC: (-?\d+)', head, re.M).group(1))
        res.append({'args': args, 'rc': rc, 'out': re.sub(r'\d{5,}', 'N', out), 'err': re.sub(r'\d{5,}', 'N', err)})
    return res


def main(outfile, *dirs):
    runs = ['1', '3', '4', '5', '6']
    lines = []
    for run, d in zip(runs, dirs):
        for section in ('s1', 's2', 's3', 's4'):
            try:
                rows = load(f'{d}/out/{section}/tcoffee_external_oracle.tsv')
            except FileNotFoundError:
                continue
            for i, r in enumerate(rows):
                e = json.loads(r['extra'])
                x = {'pdb': e['pdb']}
                if section == 's1':
                    x['sap'] = e.get('sap', {})
                x['err'] = err_lines(e.get('stderr', ''))
                if section in ('s3', 's4') and 'TMalign_pair' in r['flags']:
                    x['tm'] = e['out'].get('TMalign_pair', {})        # what the generator's own TM-align wrote for each pair (the input of the TMalign_pair part of the set)
                if section == 's4':
                    x['calls'] = calls(e)
                lines.append('\t'.join([run, section, str(i), r['flags'], r['seqs'], r['rows'], json.dumps(x, separators=(',', ':'))]))
    prov = ['# Seqeron review B09 WP68: runs 1, 3, 4, 5, 6 of the Docker package A2-12b (Ubuntu 24.04 linux/amd64 on the user\'s machine; t_coffee 13.41.0.28bdc39, SAP 1.1.3 sha256 d77350dd..., DaliLite v5 / v3.3) with the same sets (seeds of generate_tcoffee_external_oracle.py).',
            '# columns: run section idx flags seqs rows extra   (see scripts/alignment/convert_tcoffee_sap_dali_returned.py)']
    data = ('\n'.join(prov) + '\n' + '\n'.join(lines) + '\n').encode()
    with open(outfile, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', compresslevel=9, mtime=0, filename='') as f:
        f.write(data)
    print(len(lines), 'rows', len(data), 'bytes')


if __name__ == '__main__':
    main(*sys.argv[1:])
