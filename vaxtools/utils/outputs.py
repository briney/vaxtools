#!/usr/bin/env python
# filename: outputs.py


#
# Copyright (c) 2015 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import os

import numpy as np
import pandas as pd

from abutils.utils.pipeline import make_dir


def schief_csv_output(pairs, output_file, sep=',', legacy_abstar=True):
    make_dir(os.path.dirname(output_file))
    header = _get_schief_output_header(sep)
    output = [header, ]
    for p in sorted(pairs, key=lambda x: _get_name(x)):
        name = _get_name(p)
        line = [name, ]
        line += _get_pair_metadata(p)
        line += _schief_output_line(p.heavy, legacy_abstar)
        line += _schief_output_line(p.light, legacy_abstar)
        output.append(sep.join([str(l) for l in line]))
    open(output_file, 'w').write('\n'.join(output))


def _get_name(p):
    name = ''
    if p.heavy is not None:
        if 'seq_id' in p.heavy:
            name = p.heavy['seq_id']
    elif p.light is not None:
        if 'seq_id' in p.light:
            name = p.light['seq_id']
    return name if name != '' else p.name


def _get_pair_metadata(p):
    if p.heavy is not None:
        seq = p.heavy
    else:
        seq = p.light
    experiment = seq.get('experiment', '')
    group = seq.get('group', '')
    subject = seq.get('subject', '')
    timepoint = seq.get('timepoint', '')
    # experiment = seq['experiment'] if 'experiment' in seq else ''
    # group = seq['group'] if 'group' in seq else ''
    # subject = seq['subject'] if 'subject' in seq else ''
    # timepoint = seq['timepoint'] if 'timepoint' in seq else ''
    return [experiment, group, subject, timepoint]


def _get_fr_identity(seq, res='nt'):
    len_field = 'region_len_nt' if res == 'nt' else 'region_len_aa'
    mut_field = 'region_muts_nt' if res == 'nt' else 'region_muts_aa'
    regions = ['fr1', 'fr2', 'fr3']
    length = sum([seq[len_field][region] for region in regions])
    muts = sum([len(seq[mut_field][region]['muts']) for region in regions])
    return 100. * muts / length


def _get_schief_output_header(sep):
    fields = ['Sequence ID', 'Experiment', 'Group', 'Subject', 'Timepoint', 'VH gene', 'DH gene', 'JH gene', 'CDR3 length',
              'Junction AA', 'Junction NT seq', '% VH mutation (NT)', '% FR mutation (NT)',
              '% VH mutation (AA)', '% FR mutation (AA)', 'VH insertions', 'VH deletions',
              'VDJ AA seq', 'VDJ NT seq', 'Insertion count', 'Insertion lengths',
              'Deletion count', 'Deletion lengths', 'VDJ cysteine count', 'CDR3 cysteine count',
              'VL gene', 'JL gene', 'CDR3 length', 'Junction AA', 'Junction NT seq',
              '% VL mutation (NT)', '% FR mutation (NT)', '% VL mutation (AA)', '% FR mutation (AA)',
              'VL insertions', 'VL deletions', 'VDJ AA seq', 'VDJ NT seq',
              'Insertion count', 'Insertion lengths', 'Deletion count', 'Deletion lengths',
              'VDJ cysteine count', 'CDR3 cysteine count']
    return sep.join(fields)


def _schief_output_line(seq, legacy):
    if seq is None:
        return [''] * 20
    line = []
    # line.append(seq['experiment'] if 'experiment' in seq else '')
    # line.append(seq['group'] if 'group' in seq else '')
    # line.append(seq['subject'] if 'subject' in seq else '')
    # line.append(seq['timepoint'] if 'timepoint' in seq else '')
    line.append(seq['v_gene']['gene'])
    if seq['chain'] == 'heavy':
        line.append(seq['d_gene']['gene'] if 'd_gene' in seq else '')
    line.append(seq['j_gene']['gene'])
    line.append(seq['cdr3_len'])
    line.append(seq['junc_aa'])
    line.append(seq['junc_nt'])
    line.append(100. - seq['nt_identity']['v'])
    line.append(_get_fr_identity(seq, res='nt'))
    line.append(100. - seq['aa_identity']['v'])
    line.append(_get_fr_identity(seq, res='aa'))
    line.append('yes' if 'v_ins' in seq else '')
    line.append('yes' if 'v_del' in seq else '')
    line.append(seq['vdj_aa'])
    line.append(seq['vdj_nt'])
    if 'v_ins' in seq:
        len_field = 'len' if legacy else 'length'
        line.append(len(seq['v_ins']))
        line.append('[' + ' '.join([str(i[len_field]) for i in seq['v_ins']]) + ']')
    else:
        line.append('0')
        line.append('')
    if 'v_del' in seq:
        len_field = 'len' if legacy else 'length'
        line.append(len(seq['v_del']))
        line.append('[' + ' '.join([str(i[len_field]) for i in seq['v_del']]) + ']')
    else:
        line.append('0')
        line.append('')
    line.append(seq['vdj_aa'].upper().count('C'))
    line.append(seq['cdr3_aa'].upper().count('C'))
    return line
