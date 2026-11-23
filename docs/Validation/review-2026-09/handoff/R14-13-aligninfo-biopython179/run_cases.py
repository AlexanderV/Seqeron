#!/usr/bin/env python3
"""Runs every call of cases.jsonl on Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79 (the installed one; the Dockerfile pins it) and writes results.jsonl: one line per case, {"id", "results": [{"ok": value} | {"error": [type, message]}]}.
Values: strings; replacement_dictionary as a list of [a, b, number]; pos_specific_score_matrix as {"axis": [...], "rows": [[letter, [[residue, count]...]]...], "str": str(pssm)}; information_content as a number with the ic_vector.
Floats are written with repr (shortest round trip). Usage: python3 -I run_cases.py cases.jsonl results.jsonl"""
import io, json, sys, warnings

warnings.simplefilter('ignore')
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo, print_info_content
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def build(rec):
    r = SeqRecord(Seq(rec['seq']), id=rec['id'])
    if 'weight' in rec:
        r.annotations['weight'] = rec['weight']
    return r


def pssm_json(p):
    return {'axis': [p.get_residue(i) for i in range(len(p.pssm))], 'rows': [[letter, [[k, v] for k, v in d.items()]] for letter, d in p.pssm], 'str': str(p)}


def call(summary, alignment, c):
    m, kw = c['method'], dict(c['kwargs'])
    if m == 'dumb_consensus':
        return str(summary.dumb_consensus(**kw))
    if m == 'gap_consensus':
        return str(summary.gap_consensus(**kw))
    if m == 'replacement_dictionary':
        d = summary.replacement_dictionary(**kw)
        return [[a, b, v] for (a, b), v in d.items()]
    if m == 'pos_specific_score_matrix':
        if kw.get('axis_seq') == 'dumb':
            kw['axis_seq'] = summary.dumb_consensus()
        elif 'axis_seq' in kw:
            kw['axis_seq'] = Seq(kw['axis_seq'])
        return pssm_json(summary.pos_specific_score_matrix(**kw))
    if m == 'information_content':
        v = summary.information_content(**kw)
        return {'value': v, 'ic_vector': list(summary.ic_vector)}
    if m == 'ic_vector_after':
        summary.information_content(**kw)
        return list(summary.ic_vector)
    if m == 'print_info_content':
        s = SummaryInfo(alignment)
        if kw:
            s.information_content(**kw)
        buf = io.StringIO()
        print_info_content(s, buf)
        return buf.getvalue()
    if m == 'get_column':
        return str(summary.get_column(**kw))
    raise ValueError(m)


def main():
    cases, out = sys.argv[1], sys.argv[2]
    with open(cases) as fi, open(out, 'w') as fo:
        for line in fi:
            case = json.loads(line)
            alignment = MultipleSeqAlignment([build(r) for r in case['records']])
            results = []
            for c in case['calls']:
                summary = SummaryInfo(alignment)           # a fresh summary per call: ic_vector starts empty
                try:
                    results.append({'ok': call(summary, alignment, c)})
                except Exception as e:                      # the type and the message are part of the behaviour
                    results.append({'error': [type(e).__name__, str(e)]})
            fo.write(json.dumps({'id': case['id'], 'results': results}) + '\n')
    print('done')


if __name__ == '__main__':
    main()
