#!/usr/bin/env python3
"""B09 WP68 (F311): converts the returned runs of the Docker package docs/Validation/review-2026-09/handoff/A2-12b-dali-sap-stderr (runs 1, 3, 4, 5, 6: real SAP 1.1.3, DaliLite v5 on the PATH, DaliLite v3.3, t_coffee 13.41.0.28bdc39 with stderr recorded)
into the test data tests/.../TestData/Alignment/tcoffee_sap_dali_runs_oracle.tsv.gz. One row per (run, section, set):
  run  section(s1 = sap_real, s2 = -mode 3dcoffee with the real SAP, s3 = dali_pair with DaliLite v5 on the PATH, s4 = dali_pair with DaliLite v3.3)  idx  flags  seqs  rows  extra
  extra (compact JSON) = {"pdb": [...], "sap": {<pair file>: text}  (s1 only: the raw SAP output), "err": [normalised stderr lines of t_coffee that name the stop / the kill / the warning], "tm": {pair: TM-align output} (s3 / s4 sets with TMalign_pair),
                          "calls": [{"args", "rc", "out", "err"}]  (s4 only: every call of DaliLite)}
B09 WP71 (F322): run 10 (A2-12b-s5) = section s4f: the real t_coffee 13.41 with the real DaliLite 3.3 on the FULL-BACKBONE structures (60 sets, 272 logged `DaliLite -pairwise pdb1.pdb pdb2.pdb` calls, all return code 0). One row per set, extra =
  {"pdb": [...], "calls": [{"args", "rc", "cwd" (the temporary directory of the call without its random number), "out" (the standard output of DaliLite = what the glue's `>tmp1` holds), "err", "saved" (the files DaliLite removed and the directory listing / remaining files, as logged by the wrapper)}]}
  with the calls in the order of their start (the pair jobs of t_coffee: the usable pairs x < y in lexicographic order, one worker). The standard error of t_coffee was not recorded in this run (TC_STDERR unset).
B09 WP76 (F351): run 11 (A2-12b-s5-stderr) = section s4g: the 15 sets that stopped in run 10, re-run with TC_STDERR=1. One row per stopped set; idx = the index of the same set (same sequences and flags) in section s4f of run 10; extra = {"pdb": [...], "stop": the block t_coffee wrote at the stop
  (from the line of asterisks with ERROR to the closing line, the temporary paths replaced by <p>), "err": the normalised key lines}. The 15 blocks are identical.
Usage: python3 -I convert_tcoffee_sap_dali_returned.py OUT.tsv.gz RUN1_DIR RUN3_DIR RUN4_DIR RUN5_DIR RUN6_DIR [RUN10_DIR [RUN11_DIR]]    (each = the unzipped archive directory that holds out/)"""
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


def calls_full(extra):
    """Run 10: every call of DaliLite with its standard output, standard error and the .saved log (nothing normalised: the byte-for-byte text of the program)."""
    names = [n for n in extra.get('dali', {}) if not n.endswith('.saved')]
    names.sort(key=lambda n: (int(n.split('_')[2].split('.')[0]), int(n.split('_')[1])))      # the nanosecond clock of the wrapper: the order the pair jobs ran
    res = []
    for fn in names:
        t = extra['dali'][fn]
        head, _, rest = t.partition('--- STDOUT\n')
        out, _, err = rest.partition('--- STDERR\n')
        args = re.search(r'^ARGS: (.*)$', head, re.M).group(1)
        rc = int(re.search(r'^RC: (-?\d+)', head, re.M).group(1))
        cwd = re.sub(r'\d+', 'N', re.search(r'^CWD: (.*)$', head, re.M).group(1))
        res.append({'args': args, 'rc': rc, 'cwd': cwd, 'out': out, 'err': err, 'saved': extra['dali'][fn + '.saved']})
    return res


def main(outfile, *dirs):
    runs = ['1', '3', '4', '5', '6', '10', '11']
    lines = []
    run10_rows = load(f'{dirs[5]}/out/s4/tcoffee_external_oracle.tsv') if len(dirs) > 6 else []
    for run, d in zip(runs, dirs):
        for section in (('s4f',) if run == '10' else ('s4g',) if run == '11' else ('s1', 's2', 's3', 's4')):
            try:
                rows = load(f'{d}/out/{"s4" if section in ("s4f", "s4g") else section}/tcoffee_external_oracle.tsv')
            except FileNotFoundError:
                continue
            for i, r in enumerate(rows):
                e = json.loads(r['extra'])
                x = {'pdb': e['pdb']}
                if section == 's4g':
                    twins = [k for k, q in enumerate(run10_rows) if q['seqs'] == r['seqs'] and q['flags'] == r['flags'] and q['rows'].startswith('FAIL:')]
                    assert len(twins) == 1 and r['rows'] == run10_rows[twins[0]]['rows'], 'the stopped set of run 11 must be one stopped set of run 10'
                    text = e['stderr']
                    start = text.index('******************ERROR')
                    x['stop'] = re.sub(r'/tmp/\S+', '<p>', text[start:])
                    x['err'] = err_lines(text)
                    lines.append('\t'.join([run, section, str(twins[0]), r['flags'], r['seqs'], r['rows'], json.dumps(x, separators=(',', ':'))]))
                    continue
                if section == 's1':
                    x['sap'] = e.get('sap', {})
                x['err'] = err_lines(e.get('stderr', ''))
                if section == 's4f':
                    del x['err']                                          # TC_STDERR was not set in run 10
                    x['calls'] = calls_full(e)
                if section in ('s3', 's4', 's4f') and 'TMalign_pair' in r['flags']:
                    x['tm'] = e['out'].get('TMalign_pair', {})        # what the generator's own TM-align wrote for each pair (the input of the TMalign_pair part of the set)
                if section == 's4':
                    x['calls'] = calls(e)
                lines.append('\t'.join([run, section, str(i), r['flags'], r['seqs'], r['rows'], json.dumps(x, separators=(',', ':'))]))
    prov = ['# Seqeron review B09 WP68 / WP71: runs 1, 3, 4, 5, 6, 10 (section s4f) and 11 (section s4g) of the Docker package A2-12b (Ubuntu 24.04 linux/amd64 on the user\'s machine; t_coffee 13.41.0.28bdc39, SAP 1.1.3 sha256 d77350dd..., DaliLite v5 / v3.3) with the same sets (seeds of generate_tcoffee_external_oracle.py).',
            '# columns: run section idx flags seqs rows extra   (see scripts/alignment/convert_tcoffee_sap_dali_returned.py)']
    data = ('\n'.join(prov) + '\n' + '\n'.join(lines) + '\n').encode()
    with open(outfile, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', compresslevel=9, mtime=0, filename='') as f:
        f.write(data)
    print(len(lines), 'rows', len(data), 'bytes')


if __name__ == '__main__':
    main(*sys.argv[1:])
