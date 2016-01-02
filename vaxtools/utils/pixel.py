#!/usr/bin/env python
# filename: pixel.py


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


import math
import os
import subprocess as sp

from Bio import AlignIO
from Bio.Align import AlignInfo

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from abtools.utils.alignment import mafft

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')


def make_pixel(seqs, ffile, temp_dir='/tmp', consentroid=None):
	nuc_vals = {'-': 0,
				'A': 1,
				'C': 2,
				'G': 3,
				'T': 4}
	cmap = ListedColormap(['w', '#E12427', '#3B7FB6', '#63BE7A', '#E1E383'])
	need_black = False
	data = []
	aligned_seqs = _pixel_msa(consentroid, seqs, temp_dir)
	for seq in aligned_seqs:
		seq_data = []
		s_id = seq[0]
		sequence = seq[1]
		if s_id == 'consentroid':
			sequence = 'XXXXXX--' + sequence
		else:
			sequence = '--------' + sequence
		for res in sequence:
			seq_data.append(nuc_vals.get(res.upper(), 5))
		if 5 in seq_data:
			need_black = True
		data.append(seq_data)
	if need_black:
		cmap = ListedColormap(['w', '#E12427', '#3B7FB6', '#63BE7A', '#E1E383', 'k'])
	mag = (magnitude(len(data[0])) + magnitude(len(data))) / 2
	x_dim = len(data[0]) / 10**mag
	y_dim = len(data) / 10**mag
	plt.figure(figsize=(x_dim, y_dim), dpi=100)
	plt.imshow(data, cmap=cmap, interpolation='nearest')
	plt.axis('off')
	plt.savefig(os.path.expanduser(ffile), bbox_inches='tight', dpi=400)
	plt.close()


def _pixel_msa(consentroid, seqs, temp_dir):
	fasta = ''
	if consentroid is not None:
		seqs = [('consentroid', consentroid)] + seqs
	aln = mafft(seqs)
	aln_seqs = []
	for record in aln:
		aln_seqs.append((record.id, str(record.seq)))
	return aln_seqs


def magnitude(x):
	return int(math.log10(x))
