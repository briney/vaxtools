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


import numpy as np
import pandas as pd

from vaxtools.utils.vrc01 import vrc01_class_mutation_count


def schief_csv_output(pairs, output_file, sep=','):
    header = _get_schief_output_header(sep)
    output = [header, ]
    for p in pairs:
        line = [p.name, ]
        line += _schief_output_line(p.heavy)
        line += _schief_output_line(p.light)
        output.append(sep.join([str(l) for l in line]))
    open(output_file, 'w').write('\n'.join(output))


def _get_fr_identity(seq, res='nt'):
    len_field = 'region_len_nt' if res == 'nt' else 'region_len_aa'
    mut_field = 'region_muts_nt' if res == 'nt' else 'region_muts_aa'
    regions = ['fr1', 'fr2', 'fr3']
    length = sum([seq[len_field][region] for region in regions])
    muts = sum([len(seq[mut_field][region]['muts']) for region in regions])
    return 100. * muts / length


def _get_schief_output_header(sep):
    fields = ['Sequence ID', 'VH gene', 'DH gene', 'JH gene', 'CDR3 length',
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


def _schief_output_line(seq):
    if seq is None:
        return [''] * 20
    line = []
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
    line.append('')
    line.append('')
    line.append(seq['vdj_aa'])
    line.append(seq['vdj_nt'])
    if 'v_ins' in seq:
        line.append(len(seq['v_ins']))
        line.append('[' + ' '.join([str(i['len']) for i in seq['v_ins']]) + ']')
    else:
        line.append('0')
        line.append('')
    if 'v_del' in seq:
        line.append(len(seq['v_del']))
        line.append('[' + ' '.join([str(i['len']) for i in seq['v_del']]) + ']')
    else:
        line.append('0')
        line.append('')
    line.append(seq['vdj_aa'].upper().count('C'))
    line.append(seq['cdr3_aa'].upper().count('C'))
    return line


def vrc01_summary_output(pairs, output_dir):
	all_stats = {}
	subjects = list(set([p.subject for p in pairs]))
	for subject in subjects:
		stats = {}
		all_pairs = [p for p in pairs if p.subject == subject]
		just_pairs = [p for p in all_pairs if p.is_pair]
		all_heavys = [p for p in all_pairs if p.heavy is not None]
		all_lights = [p for p in all_pairs if p.light is not None]


		# 1. the fraction of Abs in each animal that are (a) VH1-2 (b) 5-aa LCDR3 (c) both VH1-2 and 5-aa LCDR3 (d) not a or b (e) not c
		vh12_seqs = [p for p in all_heavys if p.heavy['v_gene']['gene'] == 'IGHV1-2']
		vh12_pairs = [p for p in just_pairs if p.heavy['v_gene']['gene'] == 'IGHV1-2']
		lcdr3_seqs = [p for p in all_lights if p.light['cdr3_len'] == 5]
		lcdr3_pairs = [p for p in just_pairs if p.light['cdr3_len'] == 5]
		vrc01like_pairs = [p for p in just_pairs if p.light['cdr3_len'] == 5 and p.heavy['v_gene']['gene'] == 'IGHV1-2']
		nonvrc01like_pairs = [p for p in just_pairs if not all([p.light['cdr3_len'] == 5, p.heavy['v_gene']['gene'] == 'IGHV1-2'])]
		stats['Fraction of VH1-2 heavy chains'] = 1. * len(vh12_seqs) / len(all_heavys)
		stats['Fraction of VH1-2 pairs'] = 1. * len(vh12_pairs) / len(just_pairs)
		stats['Fraction of 5AA LCDR3 light chains'] = 1. * len(lcdr3_seqs) / len(all_lights)
		stats['Fraction of 5AA LCDR3 pairs'] = 1. * len(lcdr3_pairs) / len(just_pairs)
		stats['Fraction of VRC01-like pairs'] = 1. * len(vrc01like_pairs) / len(just_pairs)
		stats['Fraction of non-VRC01-like pairs'] = 1. * len(nonvrc01like_pairs) / len(just_pairs)


		# 2. the %aa mutation (mean and sdev) in V genes (separately for H and L) for Abs in the a-e cases above.
		stats['AA mutation in VH1-2 heavy chains (mean)'] = np.mean([100. - p.heavy['aa_identity']['v'] for p in vh12_seqs])
		stats['AA mutation in VH1-2 heavy chains (stdev)'] = np.std([100. - p.heavy['aa_identity']['v'] for p in vh12_seqs])
		stats['AA mutation in 5AA LCDR3 light chains (mean)'] = np.mean([100. - p.light['aa_identity']['v'] for p in lcdr3_seqs])
		stats['AA mutation in 5AA LCDr3 light chains (stdev)'] = np.std([100. - p.light['aa_identity']['v'] for p in lcdr3_seqs])


		# 3. the number of VRC01-class mutations vs total number of mutations in H chain for cases a and c on a per animal basis.
		# Let's please agree on which VRC01-class Abs we are using to define this.
		# Please use the same list as in the Jardine, Ota, Sok Science 2015 paper: 12a12, 3BNC60, PGV04, PGV20, VRC-CH31, and VRC01.
		vh12_shared_mut_counts, vh12_total_mut_counts = vrc01_class_mutation_count([p.heavy for p in vh12_seqs])
		vrc01like_shared_mut_counts, vrc01like_total_mut_counts = vrc01_class_mutation_count([p.heavy for p in vrc01like_pairs])


		# 4. the number of VRC01-class mutations in L chain for case b. **make this lower priority because it is difficult**
		# This is tricky since different V genes are used, but we should try.
		# This will be very important in boosting experiments that include the 276 glycan (which we have already done and Devin has already sorted) so let's work out the method.

		# TODO


		# 5. the light chain V genes used in cases a-e, along with the LCDR1 length for those V genes


	# 6. the LCDR3 sequences for cases a-e

	# 7. the frequency distribution of insertion lengths in Abs in a-e

	# 8. ditto for deletion lengths

	# 9. the fraction of Abs in a-e that have a deletion in LCDR1, and the sizes of those deletions if any.

	# 10. The standard Abstar output that he usually gives me for all H-L pairs if available or just H or L individual when pairs not yet available.

	# 11. the % of aa mutations that are due to single vs double nt mutations. (for cases a-e)

	# 12. the % of silent vs non-silent mutations in cases a-e.



