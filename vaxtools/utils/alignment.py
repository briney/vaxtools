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
import subprocess as sp
import tempfile

from skbio.alignment import StripedSmithWaterman

import nwalign as nw

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

from vaxtools.utils.sequence import Sequence



# -------------------------------------
#
#     MULTIPLE SEQUENCE ALIGNMENT
#
# -------------------------------------



def mafft(sequences=None, alignment_file=None, fasta=None, fmt='fasta', threads=-1, as_file=False):
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

	::threads:: is the number of threads that MAFFT should use.
	Default is -1, which uses all available cores.

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
	if fmt == 'clustal':
		aln_format = '--clustalout '
	mafft_cline = 'mafft --thread {} {}{} > {}'.format(threads, aln_format, ffile, alignment_file)
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


def consensus(aln, name=None, threshold=0.51, ambiguous='N'):
	summary_align = AlignInfo.SummaryInfo(aln)
	consensus = summary_align.gap_consensus(threshold=threshold, ambiguous=ambiguous)
	if name is None:
		name = uuid.uuid4()
	consensus_string = str(consensus).replace('-', '')
	return (name, consensus_string.upper())


def _get_fasta_string(sequences):
	if type(sequences) == str:
		return sequences
	elif type(sequences[0]) == SeqRecord:
		return '\n'.join(['>{}\n{}'.format(seq.id, str(seq.seq).upper()) for seq in sequences])
	# elif type(sequences[0]) == Sequence:
	# 	return '\n'.join(['>{}\n{}'.format(seq.id, seq.seq) for seq in sequences])
	elif type(sequences[0]) in [list, tuple]:
		return '\n'.join(['>{}\n{}'.format(seq[0], seq[1]) for seq in sequences])



# ----------------------------
#
#     PAIRWISE ALIGNMENT
#
# ----------------------------



def local_alignment(query, target=None, targets=None, match=3, mismatch=-2, matrix=None,
	gap_open_penalty=5, gap_extend_penalty=2, aa=False):
	'''
	Wrapper for SSWAlignment, which performs fast Striped Smith-Waterman local pairwise alignment

	Input: query and target sequences
	Returns: a single SSWAlignment object or a list of multiple SSWAlignment objects

	Sequences can be one of four things:
		1) a nucleotide or amino acid sequence, as a string
		2) a Biopython SeqRecord object
		3) a VaxTools Sequence object
		4) an iterable of the format (seq_id, sequence)

	::query:: a single sequence

	::target:: can be one of two things:
		1) a single sequence, as a string
		2) an iterable containing one or more sequences as strings

	default scoring parameters:
		match = 3
		mismatch = -2
		gap_open = 5
		gap_extend = 2

	For protein sequences, set ::aa:: to True and optionally provide a scoring matrix.
	::matrix:: can be one of two things:
		1) the name of a built-in matrix (current options are 'blosum62' and 'pam250')
		2) a 2D dict containing match scores for each residue pair (either aa or nt)
	'''
	if aa and not matrix:
		print('\nERROR: You must supply a scoring matrix for amino acid alignments\n')
		sys.exit(1)
	if not target and not targets:
		print('\nERROR: You must supply a target sequence (or sequences).\n')
		sys.exit(1)
	if target:
		targets = [target, ]
	alignments = []
	for t in targets:
		alignment = SSWAlignment(query=query,
								 target=t,
								 match=match,
								 mismatch=mismatch,
								 matrix=matrix,
								 gap_open=gap_open_penalty,
								 gap_extend=gap_extend_penalty,
								 aa=aa)
		alignments.append(alignment)
	if len(alignments) == 1:
		return alignments[0]
	return alignments


def global_alignment(query, target=None, targets=None, match=3, mismatch=-2, gap_open=-5, gap_extend=-2,
		score_match=None, score_mismatch=None, score_gap_open=None,
		score_gap_extend=None, matrix=None, aa=False):
	'''

	'''
	import nwalign as nw
	if not target and not targets:
		print('\nERROR: You must supply a target sequence (or sequences).\n')
		sys.exit(1)
	if target:
		targets = [target, ]
	alignments = []
	for t in targets:
		alignment = NWAlignment(query=query,
								target=t,
								match=match,
								mismatch=mismatch,
								gap_open=gap_open,
								gap_extend=gap_extend,
								score_match=score_match,
								score_mismatch=score_mismatch,
								score_gap_open=score_gap_open,
								score_gap_extend=score_gap_extend,
								matrix=matrix,
								aa=aa)
		alignments.append(alignment)
	if len(alignments) == 1:
		return alignments[0]
	return alignments


class BaseAlignment(object):
	"""docstring for BaseAlignment"""
	def __init__(self, query, target, matrix,
		match, mismatch, gap_open, gap_extend, aa):
		super(BaseAlignment, self).__init__()
		self.query = self._process_sequence(query)
		self.target = self._process_sequence(target)
		self.raw_query = query
		self.raw_target = target
		self._matrix = matrix
		self._match = int(match)
		self._mismatch = int(mismatch)
		self._gap_open = int(gap_open)
		self._gap_extend = int(gap_extend)
		self._aa = bool(aa)

	def __repr__(self):
		if len(self.aligned_query) > 20:
			qstring = '{}...{}'.format(self.aligned_query[:10], self.aligned_query[-10:])
			mstring = '{}...{}'.format(self.alignment_midline[:10], self.alignment_midline[-10:])
			tstring = '{}...{}'.format(self.aligned_target[:10], self.aligned_target[-10:])
		else:
			qstring = self.aligned_query
			mstring = self.alignment_midline
			tstring = self.aligned_target
		return_string = '\n\n'
		return_string += 'Pairwise Alignment\n'
		return_string += '------------------\n\n'
		return_string += 'query:  {}\n'.format(qstring)
		return_string += '        {}\n'.format(mstring)
		return_string += 'target: {}\n\n'.format(tstring)
		return_string += 'score: {}\n'.format(str(self.score))
		return_string += 'type: {}\n'.format(self.alignment_type)
		return_string += 'length: {}'.format(str(len(self.aligned_query)))
		print(return_string)
		return ''


	def __str__(self):
		return_string = ''
		return_string += '{}\n'.format(self.aligned_query)
		return_string += '{}\n'.format(self.alignment_midline)
		return_string += '{}\n'.format(self.aligned_target)
		return return_string

	def __len__(self):
		return len(self.aligned_query)

	def __eq__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score == other.score

	def __lt__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score < other.score

	def __le__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score <= other.score

	def __gt__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score > other.score

	def __ge__(self, other):
		if not hasattr(other, 'score'):
			if type(other) in [int, float]:
				return self.score == other
			return False
		return self.score >= other.score

	@staticmethod
	def _process_sequence(sequence):
		if type(sequence) == Sequence:
			return sequence
		return Sequence(sequence)

	def _alignment_midline(self):
		midline = ''
		for q, t in zip(self.aligned_query, self.aligned_target):
			if q == t:
				midline += '|'
			else:
				midline += ' '
		return midline


class SSWAlignment(BaseAlignment):
	"""
	Stucture for performing and analyzing a Smith-Waterman alignment
	using the skbio StripedSmithWaterman method.

	Inputs are fairly straight-forward, with a few notable things:
		1) Sequences (::query:: and ::target::) can be one of several things:
			-- just the raw sequence, as a string
			-- an iterable, of the format (sequence ID, sequence)
			-- a Biopython SeqRecord object
			-- a VaxTools Sequence object
		2) ::match:: should be a POSITIVE integer, ::mismatch:: should be a NEGATIVE integer
		3) both gap penalties (gap_open and gap_extend) should be POSITIVE integers
		4) more complex scoring matrices can be specified by name with ::matrix::
			built-in matrices are 'blosum62' and 'pam250'
		5) if you'd like to align using one set of scoring parameters and score using
			a different set, provide all of the 'score_*' parameters

	Exposed properties and convenience methods are the same as NWAlignment objects, so local
	and global alignments can be handled the same way. In fact, since comparisons are made
	based on score, local and global alignments can be directly compared with constructions like
	local_aln == global_aln and local_aln > global_aln.
	"""
	def __init__(self, query, target, match=3, mismatch=-2, matrix=None,
		gap_open=5, gap_extend=2, aa=False):
		super(SSWAlignment, self).__init__(query, target, matrix,
			match, mismatch, gap_open, gap_extend, aa)

		self.alignment_type = 'local'
		self._alignment = self._align()
		self.aligned_query = self._alignment.aligned_query_sequence
		self.aligned_target = self._alignment.aligned_target_sequence
		self.alignment_midline = self._alignment_midline()
		self.score = self._alignment.optimal_alignment_score
		self.query_begin = self._alignment.query_begin
		self.query_end = self._alignment.query_end
		self.target_begin = self._alignment.target_begin
		self.target_end = self._alignment.target_end_optimal

	def _align(self):
		aligner = StripedSmithWaterman(self.query.sequence,
									   match_score=self._match,
									   mismatch_score=self._mismatch,
									   gap_open_penalty=self._gap_open,
									   gap_extend_penalty=self._gap_extend,
									   substitution_matrix=self._matrix,
									   protein=self._aa)
		return aligner(self.target.sequence)


class NWAlignment(BaseAlignment):
	"""
	Stucture for performing and analyzing a Needleman-Wunch alignment
	using the nwalign package.

	Inputs are fairly straight-forward, with a few notable things:
		1) Sequences (::query:: and ::target::) can be one of several things:
			-- just the raw sequence, as a string
			-- an iterable, of the format (sequence ID, sequence)
			-- a Biopython SeqRecord object
			-- a VaxTools Sequence object
		2) ::match:: (and ::score_match::) should be POSITIVE integers. ::mismatch:: (and
			::score_mismatch::) should be NEGATIVE integers or 0.
		2) Gap penalties (gap_open, gap_extend) should be NEGATIVE integers or 0.
		3) more complex scoring matrices can be specified by name with ::matrix::
			built-in matrices are 'blosum62' and 'pam250'
		4) if you'd like to align using one set of scoring parameters and score using
			a different set, provide all of the 'score_*' parameters

	Exposed properties and convenience methods are the same as SWAlignment objects, so local
	and global alignments can be handled the same way. In fact, since comparisons are made
	based on score, local and global alignments can be directly compared with constructions like
	local_aln == global_aln and local_aln > global_aln.
	"""

	def __init__(self, query, target, match=3, mismatch=-2,
		gap_open=-5, gap_extend=-2,
		score_match=None, score_mismatch=None,
		score_gap_open=None, score_gap_extend=None,
		matrix=None, aa=False):
		super(NWAlignment, self).__init__(query, target, matrix,
			match, mismatch, gap_open, gap_extend, aa)
		self.alignment_type = 'global'
		self._score_match = int(score_match) if score_match else None
		self._score_mismatch = int(score_mismatch) if score_mismatch else None
		self._score_gap_open = int(score_gap_open) if score_gap_open else None
		self._score_gap_extend = int(score_gap_extend) if score_gap_extend else None
		self._alignment = self._align()
		self.aligned_query = self._alignment[0]
		self.aligned_target = self._alignment[1]
		self.alignment_midline = self._alignment_midline()
		self.score = self._score_alignment()

	def _align(self):
		matrix = self._matrix
		if matrix is None:
			matrix = self._build_matrix(match=self._match,
										mismatch=self._mismatch)
		elif matrix in ['blosum62', ]:
			matrix_dict = _get_builtin_matrix_dict(matrix)
			matrix = self._build_matrix(matrix=matrix_dict)
		elif type(matrix) == dict:
			matrix = self._build_matrix(matrix=matrix)
		aln = nw.global_align(self.query.sequence,
							  self.target.sequence,
							  gap_open=self._gap_open,
							  gap_extend=self._gap_extend,
							  matrix=matrix)
		os.unlink(matrix)
		return aln

	def _score_alignment(self):
		matrix = self._matrix
		if all([self._score_match, self._score_mismatch]):
			matrix = self._build_matrix(match=self._score_match,
										mismatch=self._score_mismatch)
		elif matrix is None:
			matrix = self._build_matrix(match=self._match,
										mismatch=self._mismatch)
		gap_open = self._score_gap_open if self._score_gap_open else self._gap_open
		gap_extend = self._score_gap_extend if self._score_gap_extend else self._gap_extend
		aln = nw.score_alignment(self.aligned_query,
								self.aligned_target,
								gap_open=gap_open,
								gap_extend=gap_extend,
								matrix=matrix)
		os.unlink(matrix)
		return aln

	def _build_matrix(self, match=None, mismatch=None, matrix=None):
		if matrix:
			return self._build_matrix_from_dict(matrix)
		return self._build_matrix_from_params(match, mismatch)

	@staticmethod
	def _build_matrix_from_dict(matrix):
		matrix_file = tempfile.NamedTemporaryFile(delete=False)
		residues = sorted(matrix.keys())
		header = '   ' + '  '.join(residues)
		matlist = [header, ]
		for r1 in residues:
			resline = [r1, ]
			for r2 in residues:
				s = str(matrix[r1][r2])
				sstring = ' {}'.format(s) if len(s) == 1 else s
			matlist.append(' '.join(resline))
		matrix_file.write('\n'.join(matlist))
		return matrix_file.name

	@staticmethod
	def _build_matrix_from_params(match, mismatch):
		mstring = ' {}'.format(match) if len(str(match)) == 1 else str(match)
		mmstring = ' {}'.format(mismatch) if len(str(mismatch)) == 1 else str(mismatch)
		matrix_file = tempfile.NamedTemporaryFile(delete=False)
		residues = ['A', 'C', 'D', 'E', 'F',
					'G', 'H', 'I', 'K', 'L',
					'M', 'N', 'P', 'Q', 'R',
					'S', 'T', 'V', 'W', 'Y', '*']
		header = '   ' + '  '.join(residues)
		matlist = [header, ]
		for r1 in residues:
			resline = [r1, ]
			for r2 in residues:
				resline.append(mstring if r1 == r2 else mmstring)
			matlist.append(' '.join(resline))
		matrix_file.write('\n'.join(matlist))
		return matrix_file.name


def _get_builtin_matrix_dict(matrix_name):
	matrices = {'blosum62': {{}},

				'pam250': {{}}
				}
