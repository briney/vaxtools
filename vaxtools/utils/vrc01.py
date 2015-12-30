#!/usr/bin/env python
# filename: vrc01.py


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


from abtools.utils.sequence import Sequence
from abtools.utils.alignment import muscle


def vrc01_class_mutation_count(seqs, chain='heavy', aa=True):
	seqs = [Sequence(s, aa=aa) for s in seqs]
	shared = []
	total = []
	# get VRC01-class sequences
	vrc01_seqs = _get_vrc01_class_sequences(chain=chain)
	vrc01_names = [s.id for s in vrc01_seqs]
	# get glVRC01 sequence
	glvrc01 = _get_vrc01_germline_sequence()
	glvrc01_name = glvrc01.id
	vrc01_seqs.append(glvrc01)
	# identify VRC01-class mutations
	for s in seqs:
		alignment_seqs = [s] + vrc01_seqs
		aln = muscle(alignment_seqs)
		aln_seq = [seq for seq in aln if seq.id == s.id][0]
		aln_gl = [seq for seq in aln if seq.id == glvrc01_name][0]
		aln_vrc01s = [seq for seq in aln if seq.id in vrc01_names]
		total.append(sum([s != g for s, g in zip(str(aln_seq.seq), str(aln_gl.seq))]))
		all_shared = {}
		for vrc01 in aln_vrc01s:
			shared = []
			for q, g, v in zip(str(aln_seq.seq), str(aln_gl.seq), str(vrc01.seq)):
				if q == v and q != g:
					shared.append(True)
				else:
					shared.append(False)
			all_shared[vrc01.id] = shared
		any_shared = 0
		for pos in zip(*all_shared.values()):
			if any(pos):
				any_shared += 1
		shared_list.append(any_shared)
	return shared, total


def _get_vrc01_germline_sequence():
	gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS')
	return Sequence(gl_vrc01)


def _get_vrc01_class_sequences(chain='heavy'):
	heavy = [('minVRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS'),
			 ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCARQKFYTGGQGWYFDLWGRGTLIVVSS'),
			 ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCARAQKRGRSEWAYAHWGQGTPVVVSS'),
			 ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCARQRSDFWDFDVWGSGTQVTVSS'),
			 ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCARDGSGDDTSWHLDPWGQGTLVIVSA'),
			 ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCARRMRSQDREWDFQHWGQGTRIIVSS')]
	light = []
	seqs = heavy if chain == 'heavy' else light
	return [Sequence(s) for s in seqs]
