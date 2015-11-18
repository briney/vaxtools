#!/usr/bin/env python
# filename: scpcr_demultiplexing.py


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


from __future__ import print_function

import os
from StringIO import StringIO
import tempfile

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

from vaxtools.utils.sequence import Sequence


def mafft(sequences=None, alignment_file=None, fasta=None, fmt='fasta', as_file=False):
	'''
	Performs multiple sequence alignment with MAFFT

	MAFFT must be installed for this to work

	Input: sequences to be aligned, or a FASTA file of sequences to be aligned
	Returns: A biopython AlignIO object, or path to the alignment file (if ::as_file::)

	::sequences:: can be one of four different things:
		1) a FASTA-formatted string of sequences
		2) a list of biopythong SeqRecord objects
		3) a list of VaxTools Sequence objects
		4) a list of lists/tuples, of the format (seq_id, sequence)

	::fasta:: can be used to provide a FASTA-formatted file of sequences
	instead of providing ::sequences::

	If a name for the alignment file is not provided (via ::alignment_file::),
	a NamedTemporaryFile will be used

	Options for alignment output format (::fmt::) are "fasta" and "clustal".
	'''
	if sequences:
		fasta_string = _get_fasta_string(sequences)
		fasta_file = tempfile.NamedTemporaryFile()
		fasta_file.write(fasta_string)
		ffile = fasta_file.name
	elif fasta:
		ffile = fasta
	if alignment_file is None:
		alignment_file = tempfile.NamedTemporaryFile().name
	aln_format = ''
	if fmt = 'clustal':
		aln_format = '--clustalout '
	mafft_cline = 'mafft --thread -1 {}{} > {}'.format(aln_format, ffile, alignment_file)
	mafft = sp.Popen(str(mafft_cline),
					 stdout=sp.PIPE,
					 stderr=sp.PIPE,
					 universal_newlines=True,
					 shell=True)
	stdout, stderr = mafft.communicate()
	os.unlink(ffile)
	if as_file:
		return alignment_file
	aln = AlignIO.read(open(alignment_file), fmt)
	os.unlink(alignment_file)
	return aln


def muscle(sequences=None, alignment_file=None, fasta=None, fmt='fasta', as_file=False):
	'''
	Performs multiple sequence alignment with MUSCLE

	MUSCLE must be installed for this to work

	Input: sequences to be aligned, or a FASTA file of sequences to be aligned
	Returns: A biopython AlignIO object, or path to the alignment file (if ::as_file::)

	::sequences:: can be one of four different things:
		1) a FASTA-formatted string of sequences
		2) a list of biopythong SeqRecord objects
		3) a list of VaxTools Sequence objects
		4) a list of lists/tuples, of the format (seq_id, sequence)

	::fasta:: can be used to provide a FASTA-formatted file of sequences
	instead of providing ::sequences::

	If a name for the alignment file is not provided (via ::alignment_file::),
	a NamedTemporaryFile will be used

	Options for alignment output format (::fmt::) are "fasta" and "clustal".
	'''
	if sequences:
		fasta_string = _get_fasta_string(sequences)
	elif fasta:
		fasta_string = open(fasta, 'r').read()
	aln_format = ''
	if fmt == 'clustal':
		aln_format = ' -clwstrict'
	muscle_cline = 'muscle {}'.format(aln_format)
	muscle = sp.Popen(str(muscle_cline),
					  stdin=sp.PIPE,
					  stdout=sp.PIPE,
					  stderr=sp.PIPE,
					  universal_newlines=True,
					  shell=True)
	alignment = muscle.communicate(input=fasta_string)[0]
	aln = AlignIO.read(StringIO(alignment), 'clustal')
	if as_file:
		if not alignment_file:
			alignment_file = tempfile.NamedTemporaryFile().name
		AlignIO.write(aln, alignment_file, fmt)
		return alignment_file
	return aln


def _get_fasta_string(sequences):
	if type(sequences) == str:
		return sequences
	elif type(sequences[0]) == SeqRecord:
		return '\n'.join(['>{}\n{}'.format(seq.id, str(seq.seq).upper()) for seq in sequences])
	elif type(sequences[0]) == Sequence:
		return '\n'.join(['>{}\n{}'.format(seq.id, seq.seq) for seq in sequences])
	elif type(sequences[0]) in [list, tuple]:
		return '\n'.join(['>{}\n{}'.format(seq[0], seq[1]) for seq in sequences])
