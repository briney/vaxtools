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


from collections import Counter
import os

import numpy as np
import pandas as pd

from abtools.utils.sequence import Sequence
from abtools.utils.alignment import muscle
from abtools.utils.pipeline import make_dir

from vaxtools.utils.mutations import single_nucleotide_aa_mutations, multiple_nucleotide_aa_mutations, silent_mutations, nonsilent_mutations
from vaxtools.utils.outputs import schief_csv_output


def vrc01_summary_output(pairs, output_dir):

	import warnings
	warnings.filterwarnings("ignore")

	make_dir(output_dir)

	# 1. the fraction of Abs in each animal that are (a) VH1-2 (b) 5-aa LCDR3 (c) both VH1-2 and 5-aa LCDR3 (d) not a or b (e) not c

	# 2. the %aa mutation (mean and sdev) in V genes (separately for H and L) for Abs in the a-e cases above.

	# 3. the number of VRC01-class mutations vs total number of mutations in H chain for cases a and c on a per animal basis.
	# Let's please agree on which VRC01-class Abs we are using to define this.
	# Please use the same list as in the Jardine, Ota, Sok Science 2015 paper: 12a12, 3BNC60, PGV04, PGV20, VRC-CH31, and VRC01.

	# 4. the number of VRC01-class mutations in L chain for case b. **make this lower priority because it is difficult**
	# This is tricky since different V genes are used, but we should try.
	# This will be very important in boosting experiments that include the 276 glycan (which we have already done and Devin has already sorted) so let's work out the method.

	# 9. the fraction of Abs in a-e that have a deletion in LCDR1, and the sizes of those deletions if any.

	# 11. the % of aa mutations that are due to single vs double nt mutations. (for cases a-e)

	# 12. the % of silent vs non-silent mutations in cases a-e.

	vrc01_summary_output_part1(pairs, output_dir)


	# 5. the light chain V genes used in cases a-e, along with the LCDR1 length for those V genes

	vrc01_summary_output_part2(pairs, output_dir)


	# 6. the LCDR3 sequences for cases a-e

	vrc01_summary_output_part3(pairs, output_dir)


	# 7. the frequency distribution of insertion lengths in Abs in a-e

	# 8. ditto for deletion lengths

	vrc01_summary_output_part4(pairs, output_dir)


	# 10. The standard Abstar output that he usually gives me for all H-L pairs if available or just H or L individual when pairs not yet available.

	vrc01_summary_output_part5(pairs, output_dir)



def vrc01_summary_output_part1(pairs, output_dir):
	all_stats = {}
	samples = list(set([p.sample for p in pairs]))
	for sample in samples:
		stats = {}
		all_pairs = [p for p in pairs if p.sample == sample]
		just_pairs = [p for p in all_pairs if p.is_pair]
		all_heavys = [p for p in all_pairs if p.heavy is not None]
		all_lights = [p for p in all_pairs if p.light is not None]

		# 1. the fraction of Abs in each animal that are (a) VH1-2 (b) 5-aa LCDR3 (c) both VH1-2 and 5-aa LCDR3 (d) not a or b (e) not c
		vh12_seqs = [p for p in all_heavys if p.heavy['v_gene']['gene'] == 'IGHV1-2']  # (a)
		vh12_pairs = [p for p in just_pairs if p.heavy['v_gene']['gene'] == 'IGHV1-2']  # (a)
		lcdr3_seqs = [p for p in all_lights if p.light['cdr3_len'] == 5]  # (b)
		lcdr3_pairs = [p for p in just_pairs if p.light['cdr3_len'] == 5]  # (b)
		vrc01like_pairs = [p for p in just_pairs if p.light['cdr3_len'] == 5 and p.heavy['v_gene']['gene'] == 'IGHV1-2']  # (c)
		strict_nonvrc01like_pairs = [p for p in just_pairs if all([p.light['cdr3_len'] != 5, p.heavy['v_gene']['gene'] != 'IGHV1-2'])]  # (d)
		nonvrc01like_pairs = [p for p in just_pairs if not all([p.light['cdr3_len'] == 5, p.heavy['v_gene']['gene'] == 'IGHV1-2'])]  # (e)
		stats['Fraction of VH1-2 heavy chains'] = 1. * len(vh12_seqs) / len(all_heavys) if len(all_heavys) > 0 else 0.0
		stats['Fraction of VH1-2 pairs'] = 1. * len(vh12_pairs) / len(just_pairs) if len(just_pairs) > 0 else 0.0
		stats['Fraction of 5AA LCDR3 light chains'] = 1. * len(lcdr3_seqs) / len(all_lights) if len(all_lights) > 0 else 0.0
		stats['Fraction of 5AA LCDR3 pairs'] = 1. * len(lcdr3_pairs) / len(just_pairs) if len(just_pairs) > 0 else 0.0
		stats['Fraction of VRC01-like pairs'] = 1. * len(vrc01like_pairs) / len(just_pairs) if len(just_pairs) > 0 else 0.0
		stats['Fraction of strict non-VRC01-like pairs'] = 1. * len(strict_nonvrc01like_pairs) / len(just_pairs) if len(just_pairs) > 0 else 0.0
		stats['Fraction of non-VRC01-like pairs'] = 1. * len(nonvrc01like_pairs) / len(just_pairs) if len(just_pairs) > 0 else 0.0

		# 2. the %aa mutation (mean and sdev) in V genes (separately for H and L) for Abs in the a-e cases above.
		stats['% AA mutation in VH1-2 heavy chains (mean)'] = np.mean([100. - p.heavy['aa_identity']['v'] for p in vh12_seqs])
		stats['% AA mutation in VH1-2 heavy chains (stdev)'] = np.std([100. - p.heavy['aa_identity']['v'] for p in vh12_seqs])
		stats['% AA mutation in 5AA LCDR3 light chains (mean)'] = np.mean([100. - p.light['aa_identity']['v'] for p in lcdr3_seqs])
		stats['% AA mutation in 5AA LCDr3 light chains (stdev)'] = np.std([100. - p.light['aa_identity']['v'] for p in lcdr3_seqs])
		stats['% AA mutation in VRC01-like paired heavy chains (mean)'] = np.mean([100. - p.heavy['aa_identity']['v'] for p in vrc01like_pairs])
		stats['% AA mutation in VRC01-like paired heavy chains (stdev)'] = np.std([100. - p.heavy['aa_identity']['v'] for p in vrc01like_pairs])
		stats['% AA mutation in VRC01-like paired light chains (mean)'] = np.mean([100. - p.light['aa_identity']['v'] for p in vrc01like_pairs])
		stats['% AA mutation in VRC01-like paired light chains (stdev)'] = np.std([100. - p.light['aa_identity']['v'] for p in vrc01like_pairs])
		stats['% AA mutation in strict nonVRC01-like paired heavy chains (mean)'] = np.mean([100. - p.heavy['aa_identity']['v'] for p in strict_nonvrc01like_pairs])
		stats['% AA mutation in strict nonVRC01-like paired heavy chains (stdev)'] = np.std([100. - p.heavy['aa_identity']['v'] for p in strict_nonvrc01like_pairs])
		stats['% AA mutation in strict nonVRC01-like paired light chains (mean)'] = np.mean([100. - p.light['aa_identity']['v'] for p in strict_nonvrc01like_pairs])
		stats['% AA mutation in strict nonVRC01-like paired light chains (stdev)'] = np.std([100. - p.light['aa_identity']['v'] for p in strict_nonvrc01like_pairs])
		stats['% AA mutation in nonVRC01-like paired heavy chains (mean)'] = np.mean([100. - p.heavy['aa_identity']['v'] for p in nonvrc01like_pairs])
		stats['% AA mutation in nonVRC01-like paired heavy chains (stdev)'] = np.std([100. - p.heavy['aa_identity']['v'] for p in nonvrc01like_pairs])
		stats['% AA mutation in nonVRC01-like paired light chains (mean)'] = np.mean([100. - p.light['aa_identity']['v'] for p in nonvrc01like_pairs])
		stats['% AA mutation in nonVRC01-like paired light chains (stdev)'] = np.std([100. - p.light['aa_identity']['v'] for p in nonvrc01like_pairs])

		# 3. the number of VRC01-class mutations vs total number of mutations in H chain for cases a and c on a per animal basis.
		vh12_shared_mut_counts, vh12_total_mut_counts = vrc01_class_mutation_count([p.heavy for p in vh12_seqs])
		vrc01like_shared_mut_counts, vrc01like_total_mut_counts = vrc01_class_mutation_count([p.heavy for p in vrc01like_pairs])
		stats['VRC01-class mutations in VH1-2 heavy chains (mean)'] = np.mean(vh12_shared_mut_counts)
		stats['Total mutations in VH1-2 heavy chains (mean)'] = np.mean(vh12_total_mut_counts)
		stats['VRC01-class mutations in VRC01-like paired heavy chains (mean)'] = np.mean(vrc01like_shared_mut_counts)
		stats['Total mutations in VRC01-like paired heavy chains (mean)'] = np.mean(vrc01like_total_mut_counts)

		# 4. the number of VRC01-class mutations in L chain for case b. **make this lower priority because it is difficult**
		# TODO

		# 9. the fraction of Abs in a-e that have a deletion in LCDR1, and the sizes of those deletions if any.
		vh12_lcdr1_flanks_lpairs = [p.light['fr1_nt'][-9:] + p.light['cdr1_nt'] + p.light['fr2_nt'][:9] for p in vh12_pairs]
		vh12_lcdr1_del_fraction_lpairs = 1. * len([p for p in vh12_lcdr1_flanks_lpairs if '-' in p]) / len(vh12_lcdr1_flanks_lpairs) if len(vh12_lcdr1_flanks_lpairs) > 0 else 0.0
		vh12_lcdr1_del_sizes_lpairs = [p.count('-') for p in vh12_lcdr1_flanks_lpairs if '-' in p]
		lcdr3_lcdr1_flanks_lseqs = [p.light['fr1_nt'][-9:] + p.light['cdr1_nt'] + p.light['fr2_nt'][:9] for p in lcdr3_seqs]
		lcdr3_lcdr1_del_fraction_lseqs = 1. * len([p for p in lcdr3_lcdr1_flanks_lseqs if '-' in p]) / len(lcdr3_lcdr1_flanks_lseqs) if len(lcdr3_lcdr1_flanks_lseqs) > 0 else 0.0
		lcdr3_lcdr1_del_sizes_lseqs = [p.count('-') for p in lcdr3_lcdr1_flanks_lseqs if '-' in p]
		lcdr3_lcdr1_flanks_lpairs = [p.light['fr1_nt'][-9:] + p.light['cdr1_nt'] + p.light['fr2_nt'][:9] for p in lcdr3_pairs]
		lcdr3_lcdr1_del_fraction_lpairs = 1. * len([p for p in lcdr3_lcdr1_flanks_lpairs if '-' in p]) / len(lcdr3_lcdr1_flanks_lpairs) if len(lcdr3_lcdr1_flanks_lpairs) > 0 else 0.0
		lcdr3_lcdr1_del_sizes_lpairs = [p.count('-') for p in lcdr3_lcdr1_flanks_lpairs if '-' in p]
		vrc01like_lcdr1_flanks_lpairs = [p.light['fr1_nt'][-9:] + p.light['cdr1_nt'] + p.light['fr2_nt'][:9] for p in vrc01like_pairs]
		vrc01like_lcdr1_del_fraction_lpairs = 1. * len([p for p in vrc01like_lcdr1_flanks_lpairs if '-' in p]) / len(vrc01like_lcdr1_flanks_lpairs) if len(vrc01like_lcdr1_flanks_lpairs) > 0 else 0.0
		vrc01like_lcdr1_del_sizes_lpairs = [p.count('-') for p in vrc01like_lcdr1_flanks_lpairs if '-' in p]
		nonvrc01like_lcdr1_flanks_lpairs = [p.light['fr1_nt'][-9:] + p.light['cdr1_nt'] + p.light['fr2_nt'][:9] for p in nonvrc01like_pairs]
		nonvrc01like_lcdr1_del_fraction_lpairs = 1. * len([p for p in nonvrc01like_lcdr1_flanks_lpairs if '-' in p]) / len(nonvrc01like_lcdr1_flanks_lpairs) if len(nonvrc01like_lcdr1_flanks_lpairs) > 0 else 0.0
		nonvrc01like_lcdr1_del_sizes_lpairs = [p.count('-') for p in nonvrc01like_lcdr1_flanks_lpairs if '-' in p]
		strict_nonvrc01like_lcdr1_flanks_lpairs = [p.light['fr1_nt'][-9:] + p.light['cdr1_nt'] + p.light['fr2_nt'][:9] for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_lcdr1_del_fraction_lpairs = 1. * len([p for p in strict_nonvrc01like_lcdr1_flanks_lpairs if '-' in p]) / len(strict_nonvrc01like_lcdr1_flanks_lpairs) if len(strict_nonvrc01like_lcdr1_flanks_lpairs) > 0 else 0.0
		strict_nonvrc01like_lcdr1_del_sizes_lpairs = [p.count('-') for p in strict_nonvrc01like_lcdr1_flanks_lpairs if '-' in p]
		stats['Fraction of LCDR1 deletions in VH1-2 paired light chains'] = vh12_lcdr1_del_fraction_lpairs
		stats['LCDR1 deletion sizes in VH1-2 paired light chains'] = ' '.join([str(d) for d in vh12_lcdr1_del_sizes_lpairs])
		stats['Fraction of LCDR1 deletions in 5AA LCDR3 light chains'] = lcdr3_lcdr1_del_fraction_lseqs
		stats['LCDR1 deletion sizes in 5AA LCDR3 light chains'] = ' '.join([str(d) for d in lcdr3_lcdr1_del_sizes_lseqs])
		stats['Fraction of LCDR1 deletions in 5AA LCDR3 paired light chains'] = lcdr3_lcdr1_del_fraction_lpairs
		stats['LCDR1 deletion sizes in 5AA LCDR3 paired light chains'] = ' '.join([str(d) for d in lcdr3_lcdr1_del_sizes_lpairs])
		stats['Fraction of LCDR1 deletions in VRC01-like paired light chains'] = vrc01like_lcdr1_del_fraction_lpairs
		stats['LCDR1 deletion sizes in VRC01-like paired light chains'] = ' '.join([str(d) for d in vrc01like_lcdr1_del_sizes_lpairs])
		stats['Fraction of LCDR1 deletions in nonVRC01-like paired light chains'] = nonvrc01like_lcdr1_del_fraction_lpairs
		stats['LCDR1 deletion sizes in nonVRC01-like paired light chains'] = ' '.join([str(d) for d in nonvrc01like_lcdr1_del_sizes_lpairs])
		stats['Fraction of LCDR1 deletions in strict nonVRC01-like paired light chains'] = strict_nonvrc01like_lcdr1_del_fraction_lpairs
		stats['LCDR1 deletion sizes in strict nonVRC01-like paired light chains'] = ' '.join([str(d) for d in strict_nonvrc01like_lcdr1_del_sizes_lpairs])

		# 11. the % of aa mutations that are due to single vs double nt mutations. (for cases a-e)
		vh12_single_nt_muts_hseqs = [single_nucleotide_aa_mutations(p.heavy) for p in vh12_seqs]
		vh12_multiple_nt_muts_hseqs = [multiple_nucleotide_aa_mutations(p.heavy) for p in vh12_seqs]
		vh12_total_nt_muts_hseqs = [s + m for s, m in zip(vh12_single_nt_muts_hseqs, vh12_multiple_nt_muts_hseqs)]
		vh12_single_nt_muts_hpairs = [single_nucleotide_aa_mutations(p.heavy) for p in vh12_pairs]
		vh12_multiple_nt_muts_hpairs = [multiple_nucleotide_aa_mutations(p.heavy) for p in vh12_pairs]
		vh12_total_nt_muts_hpairs = [s + m for s, m in zip(vh12_single_nt_muts_hpairs, vh12_multiple_nt_muts_hpairs)]
		vh12_single_nt_muts_lpairs = [single_nucleotide_aa_mutations(p.light) for p in vh12_pairs]
		vh12_multiple_nt_muts_lpairs = [multiple_nucleotide_aa_mutations(p.light) for p in vh12_pairs]
		vh12_total_nt_muts_lpairs = [s + m for s, m in zip(vh12_single_nt_muts_lpairs, vh12_multiple_nt_muts_lpairs)]
		stats['Fraction of AA mutations from single NT changes in VH1-2 heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_single_nt_muts_hseqs, vh12_total_nt_muts_hseqs)])
		stats['Fraction of AA mutations from multiple NT changes in VH1-2 heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_multiple_nt_muts_hseqs, vh12_total_nt_muts_hseqs)])
		stats['Fraction of AA mutations from single NT changes in VH1-2 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_single_nt_muts_hpairs, vh12_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from multiple NT changes in VH1-2 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_multiple_nt_muts_hpairs, vh12_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from single NT changes in VH1-2 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_single_nt_muts_lpairs, vh12_total_nt_muts_lpairs)])
		stats['Fraction of AA mutations from multiple NT changes in VH1-2 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_multiple_nt_muts_lpairs, vh12_total_nt_muts_lpairs)])

		lcdr3_single_nt_muts_lseqs = [single_nucleotide_aa_mutations(p.light) for p in lcdr3_seqs]
		lcdr3_multiple_nt_muts_lseqs = [multiple_nucleotide_aa_mutations(p.light) for p in lcdr3_seqs]
		lcdr3_total_nt_muts_lseqs = [s + m for s, m in zip(lcdr3_single_nt_muts_lseqs, lcdr3_multiple_nt_muts_lseqs)]
		lcdr3_single_nt_muts_hpairs = [single_nucleotide_aa_mutations(p.heavy) for p in lcdr3_pairs]
		lcdr3_multiple_nt_muts_hpairs = [multiple_nucleotide_aa_mutations(p.heavy) for p in lcdr3_pairs]
		lcdr3_total_nt_muts_hpairs = [s + m for s, m in zip(lcdr3_single_nt_muts_hpairs, lcdr3_multiple_nt_muts_hpairs)]
		lcdr3_single_nt_muts_lpairs = [single_nucleotide_aa_mutations(p.light) for p in lcdr3_pairs]
		lcdr3_multiple_nt_muts_lpairs = [multiple_nucleotide_aa_mutations(p.light) for p in lcdr3_pairs]
		lcdr3_total_nt_muts_lpairs = [s + m for s, m in zip(lcdr3_single_nt_muts_lpairs, lcdr3_multiple_nt_muts_lpairs)]
		stats['Fraction of AA mutations from single NT changes in 5AA LCDR3 light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_single_nt_muts_lseqs, lcdr3_total_nt_muts_lseqs)])
		stats['Fraction of AA mutations from multiple NT changes in 5AA LCDR3 light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_multiple_nt_muts_lseqs, lcdr3_total_nt_muts_lseqs)])
		stats['Fraction of AA mutations from single NT changes in 5AA LCDR3 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_single_nt_muts_lpairs, lcdr3_total_nt_muts_lpairs)])
		stats['Fraction of AA mutations from multiple NT changes in 5AA LCDR3 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_multiple_nt_muts_lpairs, lcdr3_total_nt_muts_lpairs)])
		stats['Fraction of AA mutations from single NT changes in 5AA LCDR3 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_single_nt_muts_hpairs, lcdr3_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from multiple NT changes in 5AA LCDR3 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_multiple_nt_muts_hpairs, lcdr3_total_nt_muts_hpairs)])

		vrc01like_single_nt_muts_hpairs = [single_nucleotide_aa_mutations(p.heavy) for p in vrc01like_pairs]
		vrc01like_multiple_nt_muts_hpairs = [multiple_nucleotide_aa_mutations(p.heavy) for p in vrc01like_pairs]
		vrc01like_total_nt_muts_hpairs = [s + m for s, m in zip(vrc01like_single_nt_muts_hpairs, vrc01like_multiple_nt_muts_hpairs)]
		vrc01like_single_nt_muts_lpairs = [single_nucleotide_aa_mutations(p.light) for p in vrc01like_pairs]
		vrc01like_multiple_nt_muts_lpairs = [multiple_nucleotide_aa_mutations(p.light) for p in vrc01like_pairs]
		vrc01like_total_nt_muts_lpairs = [s + m for s, m in zip(vrc01like_single_nt_muts_lpairs, vrc01like_multiple_nt_muts_lpairs)]
		stats['Fraction of AA mutations from single NT changes in VRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_single_nt_muts_hpairs, vrc01like_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from multiple NT changes in VRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_multiple_nt_muts_hpairs, vrc01like_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from single NT changes in VRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_single_nt_muts_lpairs, vrc01like_total_nt_muts_lpairs)])
		stats['Fraction of AA mutations from multiple NT changes in VRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_multiple_nt_muts_lpairs, vrc01like_total_nt_muts_lpairs)])

		nonvrc01like_single_nt_muts_hpairs = [single_nucleotide_aa_mutations(p.heavy) for p in nonvrc01like_pairs]
		nonvrc01like_multiple_nt_muts_hpairs = [multiple_nucleotide_aa_mutations(p.heavy) for p in nonvrc01like_pairs]
		nonvrc01like_total_nt_muts_hpairs = [s + m for s, m in zip(nonvrc01like_single_nt_muts_hpairs, nonvrc01like_multiple_nt_muts_hpairs)]
		nonvrc01like_single_nt_muts_lpairs = [single_nucleotide_aa_mutations(p.light) for p in nonvrc01like_pairs]
		nonvrc01like_multiple_nt_muts_lpairs = [multiple_nucleotide_aa_mutations(p.light) for p in nonvrc01like_pairs]
		nonvrc01like_total_nt_muts_lpairs = [s + m for s, m in zip(nonvrc01like_single_nt_muts_lpairs, nonvrc01like_multiple_nt_muts_lpairs)]
		stats['Fraction of AA mutations from single NT changes in nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_single_nt_muts_hpairs, nonvrc01like_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from multiple NT changes in nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_multiple_nt_muts_hpairs, nonvrc01like_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from single NT changes in nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_single_nt_muts_lpairs, nonvrc01like_total_nt_muts_lpairs)])
		stats['Fraction of AA mutations from multiple NT changes in nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_multiple_nt_muts_lpairs, nonvrc01like_total_nt_muts_lpairs)])

		strict_nonvrc01like_single_nt_muts_hpairs = [single_nucleotide_aa_mutations(p.heavy) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_multiple_nt_muts_hpairs = [multiple_nucleotide_aa_mutations(p.heavy) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_total_nt_muts_hpairs = [s + m for s, m in zip(strict_nonvrc01like_single_nt_muts_hpairs, strict_nonvrc01like_multiple_nt_muts_hpairs)]
		strict_nonvrc01like_single_nt_muts_lpairs = [single_nucleotide_aa_mutations(p.light) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_multiple_nt_muts_lpairs = [multiple_nucleotide_aa_mutations(p.light) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_total_nt_muts_lpairs = [s + m for s, m in zip(strict_nonvrc01like_single_nt_muts_lpairs, strict_nonvrc01like_multiple_nt_muts_lpairs)]
		stats['Fraction of AA mutations from single NT changes in strict nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_single_nt_muts_hpairs, strict_nonvrc01like_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from multiple NT changes in strict nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_multiple_nt_muts_hpairs, strict_nonvrc01like_total_nt_muts_hpairs)])
		stats['Fraction of AA mutations from single NT changes in strict nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_single_nt_muts_lpairs, strict_nonvrc01like_total_nt_muts_lpairs)])
		stats['Fraction of AA mutations from multiple NT changes in strict nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_multiple_nt_muts_lpairs, strict_nonvrc01like_total_nt_muts_lpairs)])

		# 12. the % of silent vs non-silent mutations in cases a-e.
		vh12_silent_muts_hseqs = [silent_mutations(p.heavy) for p in vh12_seqs]
		vh12_nonsilent_muts_hseqs = [nonsilent_mutations(p.heavy) for p in vh12_seqs]
		vh12_total_muts_hseqs = [s + n for s, n in zip(vh12_silent_muts_hseqs, vh12_nonsilent_muts_hseqs)]
		vh12_silent_muts_hpairs = [silent_mutations(p.heavy) for p in vh12_pairs]
		vh12_nonsilent_muts_hpairs = [nonsilent_mutations(p.heavy) for p in vh12_pairs]
		vh12_total_muts_hpairs = [s + n for s, n in zip(vh12_silent_muts_hpairs, vh12_nonsilent_muts_hpairs)]
		vh12_silent_muts_lpairs = [silent_mutations(p.light) for p in vh12_pairs]
		vh12_nonsilent_muts_lpairs = [nonsilent_mutations(p.light) for p in vh12_pairs]
		vh12_total_muts_lpairs = [s + n for s, n in zip(vh12_silent_muts_lpairs, vh12_nonsilent_muts_lpairs)]
		stats['Fraction of silent NT mutations in VH1-2 heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_silent_muts_hseqs, vh12_total_muts_hseqs)])
		stats['Fraction of non-silent NT mutations in VH1-2 heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_nonsilent_muts_hseqs, vh12_total_muts_hseqs)])
		stats['Fraction of silent NT mutations in VH1-2 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_silent_muts_hpairs, vh12_total_muts_hpairs)])
		stats['Fraction of non-silent NT mutations in VH1-2 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_nonsilent_muts_hpairs, vh12_total_muts_hpairs)])
		stats['Fraction of silent NT mutations in VH1-2 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_silent_muts_lpairs, vh12_total_muts_lpairs)])
		stats['Fraction of non-silent NT mutations in VH1-2 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vh12_nonsilent_muts_lpairs, vh12_total_muts_lpairs)])

		lcdr3_silent_muts_lseqs = [silent_mutations(p.light) for p in lcdr3_seqs]
		lcdr3_nonsilent_muts_lseqs = [nonsilent_mutations(p.light) for p in lcdr3_seqs]
		lcdr3_total_muts_lseqs = [s + n for s, n in zip(lcdr3_silent_muts_lseqs, lcdr3_nonsilent_muts_lseqs)]
		lcdr3_silent_muts_hpairs = [silent_mutations(p.heavy) for p in lcdr3_pairs]
		lcdr3_nonsilent_muts_hpairs = [nonsilent_mutations(p.heavy) for p in lcdr3_pairs]
		lcdr3_total_muts_hpairs = [s + n for s, n in zip(lcdr3_silent_muts_hpairs, lcdr3_nonsilent_muts_hpairs)]
		lcdr3_silent_muts_lpairs = [silent_mutations(p.light) for p in lcdr3_pairs]
		lcdr3_nonsilent_muts_lpairs = [nonsilent_mutations(p.light) for p in lcdr3_pairs]
		lcdr3_total_muts_lpairs = [s + n for s, n in zip(lcdr3_silent_muts_lpairs, lcdr3_nonsilent_muts_lpairs)]
		stats['Fraction of silent NT mutations in 5AA LCDR3 light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_silent_muts_lseqs, lcdr3_total_muts_lseqs)])
		stats['Fraction of non-silent NT mutations in 5AA LCDR3 light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_nonsilent_muts_lseqs, lcdr3_total_muts_lseqs)])
		stats['Fraction of silent NT mutations in 5AA LCDR3 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_silent_muts_hpairs, lcdr3_total_muts_hpairs)])
		stats['Fraction of non-silent NT mutations in 5AA LCDR3 paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_nonsilent_muts_hpairs, lcdr3_total_muts_hpairs)])
		stats['Fraction of silent NT mutations in 5AA LCDR3 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_silent_muts_lpairs, lcdr3_total_muts_lpairs)])
		stats['Fraction of non-silent NT mutations in 5AA LCDR3 paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(lcdr3_nonsilent_muts_lpairs, lcdr3_total_muts_lpairs)])

		vrc01like_silent_muts_hpairs = [silent_mutations(p.heavy) for p in vrc01like_pairs]
		vrc01like_nonsilent_muts_hpairs = [nonsilent_mutations(p.heavy) for p in vrc01like_pairs]
		vrc01like_total_muts_hpairs = [s + n for s, n in zip(vrc01like_silent_muts_hpairs, vrc01like_nonsilent_muts_hpairs)]
		vrc01like_silent_muts_lpairs = [silent_mutations(p.light) for p in vrc01like_pairs]
		vrc01like_nonsilent_muts_lpairs = [nonsilent_mutations(p.light) for p in vrc01like_pairs]
		vrc01like_total_muts_lpairs = [s + n for s, n in zip(vrc01like_silent_muts_lpairs, vrc01like_nonsilent_muts_lpairs)]
		stats['Fraction of silent NT mutations in VRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_silent_muts_hpairs, vrc01like_total_muts_hpairs)])
		stats['Fraction of non-silent NT mutations in VRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_nonsilent_muts_hpairs, vrc01like_total_muts_hpairs)])
		stats['Fraction of silent NT mutations in VRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_silent_muts_lpairs, vrc01like_total_muts_lpairs)])
		stats['Fraction of non-silent NT mutations in VRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(vrc01like_nonsilent_muts_lpairs, vrc01like_total_muts_lpairs)])

		nonvrc01like_silent_muts_hpairs = [silent_mutations(p.heavy) for p in nonvrc01like_pairs]
		nonvrc01like_nonsilent_muts_hpairs = [nonsilent_mutations(p.heavy) for p in nonvrc01like_pairs]
		nonvrc01like_total_muts_hpairs = [s + n for s, n in zip(nonvrc01like_silent_muts_hpairs, nonvrc01like_nonsilent_muts_hpairs)]
		nonvrc01like_silent_muts_lpairs = [silent_mutations(p.light) for p in nonvrc01like_pairs]
		nonvrc01like_nonsilent_muts_lpairs = [nonsilent_mutations(p.light) for p in nonvrc01like_pairs]
		nonvrc01like_total_muts_lpairs = [s + n for s, n in zip(nonvrc01like_silent_muts_lpairs, nonvrc01like_nonsilent_muts_lpairs)]
		stats['Fraction of silent NT mutations in nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_silent_muts_hpairs, nonvrc01like_total_muts_hpairs)])
		stats['Fraction of non-silent NT mutations in nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_nonsilent_muts_hpairs, nonvrc01like_total_muts_hpairs)])
		stats['Fraction of silent NT mutations in nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_silent_muts_lpairs, nonvrc01like_total_muts_lpairs)])
		stats['Fraction of non-silent NT mutations in nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(nonvrc01like_nonsilent_muts_lpairs, nonvrc01like_total_muts_lpairs)])

		strict_nonvrc01like_silent_muts_hpairs = [silent_mutations(p.heavy) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_nonsilent_muts_hpairs = [nonsilent_mutations(p.heavy) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_total_muts_hpairs = [s + n for s, n in zip(strict_nonvrc01like_silent_muts_hpairs, strict_nonvrc01like_nonsilent_muts_hpairs)]
		strict_nonvrc01like_silent_muts_lpairs = [silent_mutations(p.light) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_nonsilent_muts_lpairs = [nonsilent_mutations(p.light) for p in strict_nonvrc01like_pairs]
		strict_nonvrc01like_total_muts_lpairs = [s + n for s, n in zip(strict_nonvrc01like_silent_muts_lpairs, strict_nonvrc01like_nonsilent_muts_lpairs)]
		stats['Fraction of silent NT mutations in strict nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_silent_muts_hpairs, strict_nonvrc01like_total_muts_hpairs)])
		stats['Fraction of non-silent NT mutations in strict nonVRC01-like paired heavy chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_nonsilent_muts_hpairs, strict_nonvrc01like_total_muts_hpairs)])
		stats['Fraction of silent NT mutations in strict nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_silent_muts_lpairs, strict_nonvrc01like_total_muts_lpairs)])
		stats['Fraction of non-silent NT mutations in strict nonVRC01-like paired light chains'] = np.mean([1. * s / t if t > 0 else 0.0 for s, t in zip(strict_nonvrc01like_nonsilent_muts_lpairs, strict_nonvrc01like_total_muts_lpairs)])

		group = all_pairs[0].group
		all_stats['{}_{}'.format(group, sample)] = stats
	df = pd.DataFrame(all_stats).fillna(0)
	stats_csv = df.T.to_csv()
	open(os.path.join(output_dir, 'summary_output_part1.csv'), 'w').write(stats_csv)


def vrc01_summary_output_part2(pairs, output_dir):

	# 5. the light chain V genes used in cases a-e, along with the LCDR1 length for those V genes
	vl_gene_frequency_dir = os.path.join(output_dir, 'VL_gene_frequencies')
	make_dir(vl_gene_frequency_dir)
	samples = list(set([p.sample for p in pairs]))

	just_pairs = [p for p in pairs if p.is_pair]
	all_lights = [p for p in pairs if p.light is not None]

	vh12_lpairs = [p.light for p in just_pairs if p.heavy['v_gene']['gene'] == 'IGHV1-2']
	lcdr3_lseqs = [p.light for p in all_lights if p.light['cdr3_len'] == 5]
	lcdr3_lpairs = [p.light for p in just_pairs if p.light['cdr3_len'] == 5]
	vrc01like_lpairs = [p.light for p in just_pairs if p.light['cdr3_len'] == 5 and p.heavy['v_gene']['gene'] == 'IGHV1-2']
	strict_nonvrc01like_lpairs = [p.light for p in just_pairs if all([p.light['cdr3_len'] != 5, p.heavy['v_gene']['gene'] != 'IGHV1-2'])]
	nonvrc01like_lpairs = [p.light for p in just_pairs if not all([p.light['cdr3_len'] == 5, p.heavy['v_gene']['gene'] == 'IGHV1-2'])]

	sequences = {'VH1-2 paired light chain VL gene distributions': vh12_lpairs,
				 '5AA LCDR3 light chain VL gene distributions': lcdr3_lseqs,
				 '5AA LCDR3 paired light chain VL gene distributions': lcdr3_lpairs,
				 'VRC01-like paired light chain VL gene distributions': vrc01like_lpairs,
				 'strict nonVRC01-like paired light chain VL gene distributions': strict_nonvrc01like_lpairs,
				 'nonVRC01-like paired light chain VL gene distributions': nonvrc01like_lpairs}

	for sname in sequences.keys():
		data = {}
		seqs = sequences[sname]
		vl_lengths = {}
		for s in seqs:
			if s['v_gene']['gene'] not in vl_lengths:
				vl_lengths[s['v_gene']['gene']] = len(s['cdr1_germ_aa'].replace('-', ''))
		for sample in samples:
			sample_seqs = [s for s in seqs if s['sample'] == sample]
			if not sample_seqs:
				continue
			group = sample_seqs[0]['group']
			vl_genes = [s['v_gene']['gene'] for s in sample_seqs if s['chain'] in ['kappa', 'lambda']]
			vl_counts = Counter(vl_genes)
			data['{}_{}'.format(group, sample)] = vl_counts

		df = pd.DataFrame(data)
		df = df / df.sum()
		df = df.fillna(0)

		lengths = pd.Series([vl_lengths[v] for v in df.index], index=df.index)
		df['LCDR1 length'] = lengths

		outfile = os.path.join(vl_gene_frequency_dir, sname.replace(' ', '_') + '.csv')
		open(outfile, 'w').write(df.to_csv(sep=','))



def vrc01_summary_output_part3(pairs, output_dir):

	# 6. the LCDR3 sequences for cases a-e
	lcdr3_dir = os.path.join(output_dir, 'LCDR3_sequences')
	make_dir(lcdr3_dir)
	groups = sorted(list(set([p.group for p in pairs])))
	for group in groups:
		stats = {}
		all_pairs = sorted([p for p in pairs if p.group == group], key=lambda x: x.sample)
		just_pairs = [p for p in all_pairs if p.is_pair]
		all_heavys = [p for p in all_pairs if p.heavy is not None]
		all_lights = [p for p in all_pairs if p.light is not None]

		vh12_seqs = [p for p in all_heavys if p.heavy['v_gene']['gene'] == 'IGHV1-2']  # (a)
		vh12_pairs = [p for p in just_pairs if p.heavy['v_gene']['gene'] == 'IGHV1-2']  # (a)
		lcdr3_seqs = [p for p in all_lights if p.light['cdr3_len'] == 5]  # (b)
		lcdr3_pairs = [p for p in just_pairs if p.light['cdr3_len'] == 5]  # (b)
		vrc01like_pairs = [p for p in just_pairs if p.light['cdr3_len'] == 5 and p.heavy['v_gene']['gene'] == 'IGHV1-2']  # (c)
		strict_nonvrc01like_pairs = [p for p in just_pairs if all([p.light['cdr3_len'] != 5, p.heavy['v_gene']['gene'] != 'IGHV1-2'])]  # (d)
		nonvrc01like_pairs = [p for p in just_pairs if not all([p.light['cdr3_len'] == 5, p.heavy['v_gene']['gene'] == 'IGHV1-2'])]  # (e)

		vh12_pair_cdr3_csvs = ['{},{},{}'.format(p.sample, p.name, p.light['cdr3_aa']) for p in vh12_pairs]
		lcdr3_seq_cdr3_csvs = ['{},{},{}'.format(p.sample, p.name, p.light['cdr3_aa']) for p in lcdr3_seqs]
		lcdr3_pair_cdr3_csvs = ['{},{},{}'.format(p.sample, p.name, p.light['cdr3_aa']) for p in lcdr3_pairs]
		vrc01like_pair_cdr3_csvs = ['{},{},{}'.format(p.sample, p.name, p.light['cdr3_aa']) for p in vrc01like_pairs]
		strict_nonvrc01like_pair_cdr3_csvs = ['{},{},{}'.format(p.sample, p.name, p.light['cdr3_aa']) for p in strict_nonvrc01like_pairs]
		nonvrc01like_pair_cdr3_csvs = ['{},{},{}'.format(p.sample, p.name, p.light['cdr3_aa']) for p in nonvrc01like_pairs]

		open(os.path.join(lcdr3_dir, '{}_VH1-2_pair_LCDR3_sequences.fasta'.format(group)), 'w').write('\n'.join(vh12_pair_cdr3_csvs))
		open(os.path.join(lcdr3_dir, '{}_5AA_LCDR3_seq_LCDR3_sequences.fasta'.format(group)), 'w').write('\n'.join(lcdr3_seq_cdr3_csvs))
		open(os.path.join(lcdr3_dir, '{}_5AA_LCDR3_pair_LCDR3_sequences.fasta'.format(group)), 'w').write('\n'.join(lcdr3_pair_cdr3_csvs))
		open(os.path.join(lcdr3_dir, '{}_VRC01-like_pair_LCDR3_sequences.fasta'.format(group)), 'w').write('\n'.join(vrc01like_pair_cdr3_csvs))
		open(os.path.join(lcdr3_dir, '{}_nonVRC01-like_pair_LCDR3_sequences.fasta'.format(group)), 'w').write('\n'.join(nonvrc01like_pair_cdr3_csvs))
		open(os.path.join(lcdr3_dir, '{}_strict_nonVRC01-like_pair_LCDR3_sequences.fasta'.format(group)), 'w').write('\n'.join(strict_nonvrc01like_pair_cdr3_csvs))



def vrc01_summary_output_part4(pairs, output_dir):

	# 7. the frequency distribution of insertion lengths in Abs in a-e
	# 8. ditto for deletion lengths

	indel_dir = os.path.join(output_dir, 'indels')
	make_dir(indel_dir)
	samples = list(set([p.sample for p in pairs]))

	just_pairs = [p for p in pairs if p.is_pair]
	all_heavys = [p for p in pairs if p.heavy is not None]
	all_lights = [p for p in pairs if p.light is not None]

	vh12_hseqs = [p.heavy for p in all_heavys if p.heavy['v_gene']['gene'] == 'IGHV1-2']
	vh12_hpairs = [p.heavy for p in just_pairs if p.heavy['v_gene']['gene'] == 'IGHV1-2']
	vh12_lpairs = [p.light for p in just_pairs if p.heavy['v_gene']['gene'] == 'IGHV1-2']
	lcdr3_lseqs = [p.light for p in all_lights if p.light['cdr3_len'] == 5]
	lcdr3_hpairs = [p.heavy for p in just_pairs if p.heavy['cdr3_len'] == 5]
	lcdr3_lpairs = [p.light for p in just_pairs if p.light['cdr3_len'] == 5]
	vrc01like_hpairs = [p.heavy for p in just_pairs if p.light['cdr3_len'] == 5 and p.heavy['v_gene']['gene'] == 'IGHV1-2']
	vrc01like_lpairs = [p.light for p in just_pairs if p.light['cdr3_len'] == 5 and p.heavy['v_gene']['gene'] == 'IGHV1-2']
	strict_nonvrc01like_hpairs = [p.heavy for p in just_pairs if all([p.light['cdr3_len'] != 5, p.heavy['v_gene']['gene'] != 'IGHV1-2'])]
	strict_nonvrc01like_lpairs = [p.light for p in just_pairs if all([p.light['cdr3_len'] != 5, p.heavy['v_gene']['gene'] != 'IGHV1-2'])]
	nonvrc01like_hpairs = [p.heavy for p in just_pairs if not all([p.light['cdr3_len'] == 5, p.heavy['v_gene']['gene'] == 'IGHV1-2'])]
	nonvrc01like_lpairs = [p.light for p in just_pairs if not all([p.light['cdr3_len'] == 5, p.heavy['v_gene']['gene'] == 'IGHV1-2'])]

	sequences = {'VH1-2 heavy chain indel distributions': vh12_hseqs,
				 'VH1-2 paired heavy chain indel distributions': vh12_hpairs,
				 'VH1-2 paired light chain indel distributions': vh12_lpairs,
				 '5AA LCDR3 light chain indel distributions': lcdr3_lseqs,
				 '5AA LCDR3 paired heavy chain indel distributions': lcdr3_hpairs,
				 '5AA LCDR3 paired light chain indel distributions': lcdr3_lpairs,
				 'VRC01-like paired heavy chain indel distributions': vrc01like_hpairs,
				 'VRC01-like paired light chain indel distributions': vrc01like_lpairs,
				 'strict nonVRC01-like paired heavy chain indel distributions': strict_nonvrc01like_hpairs,
				 'strict nonVRC01-like paired light chain indel distributions': strict_nonvrc01like_lpairs,
				 'nonVRC01-like paired heavy chain indel distributions': nonvrc01like_hpairs,
				 'nonVRC01-like paired light chain indel distributions': nonvrc01like_lpairs}

	for sname in sequences.keys():
		data = {}
		seqs = sequences[sname]
		for sample in samples:
			sample_seqs = [s for s in seqs if s['sample'] == sample]
			if not sample_seqs:
				continue
			group = sample_seqs[0]['group']
			ins_counts = np.bincount([indel['len'] for s in sample_seqs if 'v_ins' in s for indel in s['v_ins']])
			ins_lengths = range(len(ins_counts))
			ins_dist = {l: c for l, c in zip(ins_lengths, ins_counts) if c > 0}
			del_counts = np.bincount([indel['len'] for s in sample_seqs if 'v_del' in s for indel in s['v_del']])
			del_lengths = range(len(del_counts))
			del_dist = {l: c for l, c in zip(del_lengths, del_counts) if c > 0}

			data['{}_{}_insertions'.format(group, sample)] = ins_dist
			data['{}_{}_deletions'.format(group, sample)] = del_dist

		df = pd.DataFrame(data)
		df = df / df.sum()
		df = df.fillna(0)
		outfile = os.path.join(indel_dir, sname.replace(' ', '_') + '.csv')
		open(outfile, 'w').write(df.to_csv(sep=','))



def vrc01_summary_output_part5(pairs, output_dir):

	# 10. The standard Abstar output that he usually gives me for all H-L pairs if available or just H or L individual when pairs not yet available.

	output_file = os.path.join(output_dir, 'abstar_output.csv')
	schief_csv_output(pairs, output_file)



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
			_shared = []
			for q, g, v in zip(str(aln_seq.seq), str(aln_gl.seq), str(vrc01.seq)):
				if q == v and q != g:
					_shared.append(True)
				else:
					_shared.append(False)
			all_shared[vrc01.id] = _shared
		any_shared = 0
		for pos in zip(*all_shared.values()):
			if any(pos):
				any_shared += 1
		shared.append(any_shared)
	return shared, total


def _get_vrc01_germline_sequence():
	gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS')
	return Sequence(gl_vrc01)


def _get_vrc01_class_sequences(chain='heavy'):
	heavy = [('VRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS'),
			 ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCARQKFYTGGQGWYFDLWGRGTLIVVSS'),
			 ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCARAQKRGRSEWAYAHWGQGTPVVVSS'),
			 ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCARQRSDFWDFDVWGSGTQVTVSS'),
			 ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCARDGSGDDTSWHLDPWGQGTLVIVSA'),
			 ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCARRMRSQDREWDFQHWGQGTRIIVSS')]
	light = []
	seqs = heavy if chain == 'heavy' else light
	return [Sequence(s) for s in seqs]


def _get_minvrc01_sequence():
	minvrc01 = ('minVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGCTLNWVRQAPGQGLEWMGWIKPRFGAVNYARKFQGRVTMTRDVYSDTAYMELSRLRSDDTAVYYCARGKNCDYNWDFQHWGQGTLVTVSS')
	return Sequence(minvrc01)
