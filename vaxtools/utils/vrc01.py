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


from __future__ import print_function

from collections import Counter, OrderedDict
from itertools import izip_longest
import os
import random

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg', warn=False)
import seaborn as sns
import matplotlib.pyplot as plt

from abtools.sequence import Sequence
from abtools.alignment import muscle
from abtools.pipeline import make_dir

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
        stats = OrderedDict()
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
        stats['Subject'] = all_pairs[0].subject if all_pairs[0].subject is not None else ''
        stats['Experiment'] = all_pairs[0].experiment if all_pairs[0].experiment is not None else ''
        stats['Group'] = all_pairs[0].group if all_pairs[0].group is not None else ''
        stats['Timepoint'] = all_pairs[0].timepoint if all_pairs[0].timepoint is not None else ''
        stats['Number of VH1-2 heavy chains'] = len(vh12_seqs)
        stats['Number of VH1-2 pairs'] = len(vh12_pairs)
        stats['Number of 5AA LCDR3 light chains'] = len(lcdr3_seqs)
        stats['Number of 5AA LCDR3 pairs'] = len(lcdr3_pairs)
        stats['Number of VRC01-like pairs'] = len(vrc01like_pairs)
        stats['Number of strict non-VRC01-like pairs'] = len(strict_nonvrc01like_pairs)
        stats['Number of non-VRC01-like pairs'] = len(nonvrc01like_pairs)
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
        stats['Fraction of LCDR1 deletions in 5AA LCDR3 light chains'] = lcdr3_lcdr1_del_fraction_lseqs
        stats['Fraction of LCDR1 deletions in 5AA LCDR3 paired light chains'] = lcdr3_lcdr1_del_fraction_lpairs
        stats['Fraction of LCDR1 deletions in VRC01-like paired light chains'] = vrc01like_lcdr1_del_fraction_lpairs
        stats['Fraction of LCDR1 deletions in nonVRC01-like paired light chains'] = nonvrc01like_lcdr1_del_fraction_lpairs
        stats['Fraction of LCDR1 deletions in strict nonVRC01-like paired light chains'] = strict_nonvrc01like_lcdr1_del_fraction_lpairs
        stats['LCDR1 deletion sizes in VH1-2 paired light chains'] = ' '.join([str(d) for d in vh12_lcdr1_del_sizes_lpairs])
        stats['LCDR1 deletion sizes in 5AA LCDR3 light chains'] = ' '.join([str(d) for d in lcdr3_lcdr1_del_sizes_lseqs])
        stats['LCDR1 deletion sizes in 5AA LCDR3 paired light chains'] = ' '.join([str(d) for d in lcdr3_lcdr1_del_sizes_lpairs])
        stats['LCDR1 deletion sizes in VRC01-like paired light chains'] = ' '.join([str(d) for d in vrc01like_lcdr1_del_sizes_lpairs])
        stats['LCDR1 deletion sizes in nonVRC01-like paired light chains'] = ' '.join([str(d) for d in nonvrc01like_lcdr1_del_sizes_lpairs])
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
    df = pd.DataFrame(all_stats).fillna(0).T
    summary_dfs = []
    summary_cols = [c for c in df.columns.values.tolist() if 'LCDR1 deletion sizes' not in c]
    summary_cols = [c for c in summary_cols if c not in ['Experiment', 'Group', 'Timepoint']]
    summary_dfs.append(pd.DataFrame({col: '' for col in summary_cols}, index=[""]))
    summary_dfs.append(pd.DataFrame({col: df[col].min() if col != 'Subject' else '' for col in summary_cols}, index=["min"]))
    summary_dfs.append(pd.DataFrame({col: df[col].max() if col != 'Subject' else '' for col in summary_cols}, index=["max"]))
    summary_dfs.append(pd.DataFrame({col: df[col].mean() if col != 'Subject' else '' for col in summary_cols}, index=["mean"]))
    summary_dfs.append(pd.DataFrame({col: df[col].std() if col != 'Subject' else '' for col in summary_cols}, index=["std"]))
    summary_dfs.append(pd.DataFrame({col: df[col].sem() if col != 'Subject' else '' for col in summary_cols}, index=["sem"]))
    for _df in summary_dfs:
        df = df.append(_df)
    df.loc['min', 'Subject'] = 'min'
    df.loc['max', 'Subject'] = 'max'
    df.loc['mean', 'Subject'] = 'mean'
    df.loc['std', 'Subject'] = 'std'
    df.loc['sem', 'Subject'] = 'sem'
    df = df[stats.keys()]
    stats_csv = df.to_csv(sep=',', index=False)
    # stats_csv = df.to_csv(sep=',', index=True)
    open(os.path.join(output_dir, 'summary_output.csv'), 'w').write(stats_csv)


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
        names = []
        for s in seqs:
            if s['v_gene']['gene'] not in vl_lengths:
                vl_lengths[s['v_gene']['gene']] = len(s['cdr1_germ_aa'].replace('-', ''))
        for sample in samples:
            sample_seqs = [s for s in seqs if _get_sample(s) == sample]
            if not sample_seqs:
                continue
            # if 'group' in sample_seqs[0]:
            #     group = sample_seqs[0]['group']
            # else:
            #     group = 'None'
            vl_genes = [s['v_gene']['gene'] for s in sample_seqs if s['chain'] in ['kappa', 'lambda']]
            vl_counts = Counter(vl_genes)
            norm_vl_counts = {k: 1. * v / sum(vl_counts.values()) for k, v in vl_counts.items()}
            norm_vl_counts['total sequences'] = sum(vl_counts.values())
            # name = '{}_{}'.format(group, sample)
            name = sample
            names.append(name)
            data[name] = norm_vl_counts
        df = pd.DataFrame(data)
        df = df.fillna(0)
        lengths = pd.Series([vl_lengths[v] if v != 'total sequences' else '' for v in df.index], index=df.index)
        df['LCDR1 length'] = lengths
        cols = ['LCDR1 length'] + sorted(names)
        df = df[cols]
        summary_cols = sorted(names)
        df['min'] = df[summary_cols].min(axis=1)
        df['max'] = df[summary_cols].max(axis=1)
        df['mean'] = df[summary_cols].mean(axis=1)
        df['std'] = df[summary_cols].std(axis=1)
        outfile = os.path.join(vl_gene_frequency_dir, sname.replace(' ', '_') + '.csv')
        open(outfile, 'w').write(df.to_csv(sep=','))


def _get_sample(s):
    slist = []
    if 'experiment' in s.keys():
        slist.append(s['experiment'])
    if 'group' in s.keys():
        slist.append(s['group'])
    if 'subject' in s.keys():
        slist.append(s['subject'])
    if 'timepoint' in s.keys():
        slist.append(s['timepoint'])
    return '|'.join(slist)


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
        open(os.path.join(lcdr3_dir, '{}_VH1-2_pair_LCDR3_sequences.csv'.format(group)), 'w').write('\n'.join(vh12_pair_cdr3_csvs))
        open(os.path.join(lcdr3_dir, '{}_5AA_LCDR3_seq_LCDR3_sequences.csv'.format(group)), 'w').write('\n'.join(lcdr3_seq_cdr3_csvs))
        open(os.path.join(lcdr3_dir, '{}_5AA_LCDR3_pair_LCDR3_sequences.csv'.format(group)), 'w').write('\n'.join(lcdr3_pair_cdr3_csvs))
        open(os.path.join(lcdr3_dir, '{}_VRC01-like_pair_LCDR3_sequences.csv'.format(group)), 'w').write('\n'.join(vrc01like_pair_cdr3_csvs))
        open(os.path.join(lcdr3_dir, '{}_nonVRC01-like_pair_LCDR3_sequences.csv'.format(group)), 'w').write('\n'.join(nonvrc01like_pair_cdr3_csvs))
        open(os.path.join(lcdr3_dir, '{}_strict_nonVRC01-like_pair_LCDR3_sequences.csv'.format(group)), 'w').write('\n'.join(strict_nonvrc01like_pair_cdr3_csvs))


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
    lcdr3_hpairs = [p.heavy for p in just_pairs if p.light['cdr3_len'] == 5]
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
            sample_seqs = [s for s in seqs if _get_sample(s) == sample]
            # sample_seqs = [s for s in seqs if s['sample'] == sample]
            ins_seqs = [s for s in sample_seqs if 'v_ins' in s]
            del_seqs = [s for s in sample_seqs if 'v_del' in s]
            if not sample_seqs:
                continue
            # try:
            #     group = sample_seqs[0]['group']
            # except KeyError:
            #     group = 'None'
            ins_dist = Counter([indel['len'] for s in ins_seqs for indel in s['v_ins']])
            del_dist = Counter([indel['len'] for s in del_seqs for indel in s['v_del']])
            norm_ins_dist = {k: 1. * v / sum(ins_dist.values()) for k, v in ins_dist.items()}
            norm_del_dist = {k: 1. * v / sum(del_dist.values()) for k, v in del_dist.items()}
            norm_ins_dist['total sequences'] = len(sample_seqs)
            norm_ins_dist['indel sequences'] = len(ins_seqs)
            norm_del_dist['total sequences'] = len(sample_seqs)
            norm_del_dist['indel sequences'] = len(del_seqs)
            # data['{}_{}_insertions'.format(group, sample)] = norm_ins_dist
            # data['{}_{}_deletions'.format(group, sample)] = norm_del_dist
            data['{}_insertions'.format(sample)] = norm_ins_dist
            data['{}_deletions'.format(sample)] = norm_del_dist
        df = pd.DataFrame(data)
        df = df.fillna(0)
        outfile = os.path.join(indel_dir, sname.replace(' ', '_') + '.csv')
        open(outfile, 'w').write(df.to_csv(sep=','))


def vrc01_summary_output_part5(pairs, output_dir):
    # 10. The standard Abstar output that he usually gives me for all H-L pairs if available or just H or L individual when pairs not yet available.
    output_file = os.path.join(output_dir, 'abstar_output.csv')
    schief_csv_output(pairs, output_file)


def vrc01_class_mutation_count(seqs, vrc01_class=True, minvrc01=True, min12a21=True,
                               vgene_only=True, chain='heavy', print_alignments=False):
    '''
    seqs should be an iterable of anything that abtools.utils.sequence.Sequence can handle
    '''
    input_seqs = [Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]
    vrc01_seqs = []
    shared = []
    total = []
    # get VRC01-class sequences
    if vrc01_class:
        vrc01_seqs += get_vrc01_class_sequences(chain=chain, vgene_only=vgene_only)
    if minvrc01:
        vrc01_seqs.append(get_minvrc01_sequence(vgene_only=vgene_only))
    if min12a21:
        vrc01_seqs.append(get_min12a21_sequence(vgene_only=vgene_only))
    vrc01_names = [s.id for s in vrc01_seqs]
    # get glVRC01 sequence
    glvrc01 = get_vrc01_germline_sequence(vgene_only=vgene_only)
    glvrc01_name = glvrc01.id
    # identify VRC01-class mutations
    for s in input_seqs:
        alignment_seqs = [s] + vrc01_seqs + [glvrc01]
        aln = muscle(alignment_seqs)
        if print_alignments:
        	print(aln)
        aln_seq = [seq for seq in aln if seq.id == s.id][0]
        aln_gl = [seq for seq in aln if seq.id == glvrc01_name][0]
        aln_vrc01s = [seq for seq in aln if seq.id in vrc01_names]
        if vgene_only:
            total.append(sum([_s != g for _s, g in zip(str(aln_seq.seq), str(aln_gl.seq)) if g != '-']))
        else:
            total.append(sum([s != g for s, g in zip(str(aln_seq.seq), str(aln_gl.seq))]))
        all_shared = {}
        for vrc01 in aln_vrc01s:
            _shared = []
            for q, g, v in zip(str(aln_seq.seq), str(aln_gl.seq), str(vrc01.seq)):
                if vgene_only and g == '-' and v == '-':
                    _shared.append(False)
                elif q == v and q != g:
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


def vrc01_class_mutation_positions(seqs, vrc01_class=True, minvrc01=True, min12a21=True,
                                   vgene_only=False, chain='heavy', aa=True, drop_gaps=True):
    data = []
    hiv_seqs = []
    input_seqs = [Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]
    input_names = [s.id for s in input_seqs]
    if vrc01_class:
        hiv_seqs += get_vrc01_class_sequences(vgene_only=vgene_only)
    if minvrc01:
        hiv_seqs.append(get_minvrc01_sequence(vgene_only=vgene_only))
    if min12a21:
        hiv_seqs.append(get_min12a21_sequence(vgene_only=vgene_only))
    all_hiv_names = [s.id for s in hiv_seqs]
    seqs_for_alignment = input_seqs + hiv_seqs
    seqs_for_alignment.append(get_vrc01_germline_sequence(vgene_only=vgene_only))
    aln = muscle(seqs_for_alignment)
    aln_seqs = [seq for seq in aln if seq.id in input_names]
    aln_gl = [seq for seq in aln if seq.id == 'glVRC01'][0]
    aln_mins = [seq for seq in aln if seq.id in ['minVRC01', 'min12A21']]
    aln_hiv = [seq for seq in aln if seq.id in vrc01_class_names]
    for seq in aln_seqs:
        seq_data = []
        for i, (s, g) in enumerate(zip(str(seq.seq), str(aln_gl.seq))):
            if g == '-' and s == '-':
                continue
            min_residues = [seq[i] for seq in aln_min]
            vrc01_residues = [seq[i] for seq in aln_hiv]
            if s == '-':
                seq_data.append(0)
            elif s == g:
                seq_data.append(0)
            elif s != g and s in min_residues:
                seq_data.append(2)
            elif s != g and s in vrc01_muts:
                seq_data.append(3)
            elif s != g and s not in vrc01_muts:
                seq_data.append(1)
            else:
                seq_data.append(0)
        data.append(seq_data)
    return np.array(data)


def pixel_plot(data, cmap, figfile=None, pad=2, labelsize=14, maxy_denom=30, maxx_denom=10):
    '''
    ::data:: is the output from vrc01_class_mutation_positions()
    '''
    max_y, max_x = data.shape
    f, ax = plt.subplots(figsize=(max_x / float(maxx_denom), max_y / float(maxy_denom)))
    plt.pcolor(data, cmap=cmap)
    ax.set_xlim([0, max_x])
    ax.set_ylim([0, max_y])
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_color('k')
    ax.xaxis.set_minor_locator(minorLocator)
    ticks = [26, 34, 50, 58, 99, 114]
    minor_ticks = [13, 30, 42, 54, 78.5, 106.5, 118]
    minor_labels = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    ax.set_xticks(ticks)
    ax.xaxis.set_tick_params(which='major', direction='out',
                             length=6, color='k', width=1, top='off', labelsize=0)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(minor_labels, minor=True, y=-0.02)
    ax.xaxis.set_tick_params(which='minor', direction='out',
                             length=12, color='white', width=1, top='off',
                             labelsize=labelsize, pad=pad)
    ticks = [26, 34, 50, 58, 99, 114]
    plt.xticks(ticks, [' '] * len(ticks))
    ax.xaxis.set_tick_params(which='major', direction='out',
                             length=6, color='k', width=1.5, top='off')
    ax.tick_params(axis='y', labelsize=0)
    plt.tight_layout()
    if figfile is None:
        plt.show()
    else:
        plt.savefig(figfile)


def get_vrc01_germline_sequence(vgene_only=True):
    if vgene_only:
        gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR')
    else:
        gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS')
    return Sequence(gl_vrc01)


def get_vrc01_class_sequences(chain='heavy', vgene_only=True):
    if vgene_only:
        heavy = [('VRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTR'),
                 ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCAR'),
                 ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCAR'),
                 ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCAR'),
                 ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCAR'),
                 ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCAR')]
        light = []
    else:
        heavy = [('VRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS'),
                 ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCARQKFYTGGQGWYFDLWGRGTLIVVSS'),
                 ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCARAQKRGRSEWAYAHWGQGTPVVVSS'),
                 ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCARQRSDFWDFDVWGSGTQVTVSS'),
                 ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCARDGSGDDTSWHLDPWGQGTLVIVSA'),
                 ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCARRMRSQDREWDFQHWGQGTRIIVSS')]
        light = []
    seqs = heavy if chain == 'heavy' else light
    return [Sequence(s) for s in seqs]


def get_minvrc01_sequence(vgene_only=True):
    if vgene_only:
        minvrc01 = ('minVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGCTLNWVRQAPGQGLEWMGWIKPRFGAVNYARKFQGRVTMTRDVYSDTAYMELSRLRSDDTAVYYCAR')
    else:
        minvrc01 = ('minVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGCTLNWVRQAPGQGLEWMGWIKPRFGAVNYARKFQGRVTMTRDVYSDTAYMELSRLRSDDTAVYYCARGKNCDYNWDFQHWGQGTLVTVSS')
    return Sequence(minvrc01)


def get_min12a21_sequence(vgene_only=True):
    if vgene_only:
        min12a21 = ('min12A21', 'QVQLVQSGAEVKKPGASVRVSCKASGYTFTNYILHWWRQAPGQGLEWMGWIKPVFGAVNYARQFQGRVTMTRDIYREIAYMELSRLRSDDTAVYYCAR')
    else:
        min12a21 = ('min12A21', 'QVQLVQSGAEVKKPGASVRVSCKASGYTFTNYILHWWRQAPGQGLEWMGWIKPVFGAVNYARQFQGRVTMTRDIYREIAYMELSRLRSDDTAVYYCARDESGDDLKWHLHPWGQGTQVIVSP')
    return Sequence(min12a21)


def shared_mutation_kde_plot(xs, ys, cmaps, figfile=None, figsize=(8, 8), n_levels=10, alpha=0.6,
                             y_label='VRC01-class mutations', x_label='Total amino acid mutations', labelsize=14,
                             scatter=True, x_jitter=0.15, y_jitter=0.15, scatter_color='grey', scatter_alpha=0.4):
    sns.set_style('whitegrid')
    f, ax = plt.subplots(figsize=figsize)

    for x, y, cmap in izip_longest(xs, ys, cmaps):
        x = np.array(x)
        y = np.array(y)
        ax = sns.kdeplot(x, y, cmap=cmap, shade=True,
                         shade_lowest=False, alpha=alpha, n_levels=n_levels)
        ax = sns.kdeplot(x, y, cmap=cmap, shade=False,
                         shade_lowest=False, alpha=alpha, n_levels=n_levels)
        if scatter:
            random.seed(1234)
            plt.scatter([_x + random.uniform(-x_jitter, x_jitter) for _x in x],
                        [_y + random.uniform(-y_jitter, y_jitter) for _y in y],
                        alpha=scatter_alpha, color=scatter_color)
    axes = plt.gca()
    axes.set_xlim([0, max(x) + 1])
    axes.set_ylim([0, max(y) + 1])
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    plt.ylabel(y_label, size=labelsize)
    plt.xlabel(x_label, size=labelsize)
    if figfile is None:
        plt.show()
    else:
        plt.savefig(figfile)


def shared_mutation_2dhist_plot(x, y, cmap, figfile=None, figsize=None, figsize_denom=4, n_levels=10, alpha=1.0, x_lim=None, y_lim=None,
                                y_label='VRC01-class mutations', x_label='Total amino acid mutations', labelsize=14,
                                tick_labelsize=12, pad=-4):
    # adjust the inputs to make them iterable
    if type(x[0]) in [int, float, np.int64, np.float64]:
        x = [x]
        y = [y]
        cmap = [cmap]
    sns.set_style('whitegrid')
    f, ax = plt.subplots(figsize=figsize)
    # make the plots
    for _x, _y, _cmap in zip(x, y, cmap):
        bin_x = max(_x) + 2 if x_lim is None else x_lim[1] + 2
        bin_y = max(_y) + 2 if y_lim is None else y_lim[1] + 2
        bins = [[i - 0.5 for i in range(bin_y)], [i - 0.5 for i in range(bin_x)]]
        data, x_edges, y_edges = np.histogram2d(_y, _x, bins=bins)
        data = data[::-1]
        if figsize is None:
            figsize = (float(bin_x) / figsize_denom, float(bin_y) / figsize_denom)
        mask = np.array([[val == 0 for val in subl] for subl in data])
        ax = sns.heatmap(data, cmap=_cmap, square=True, cbar=False, mask=mask,
                         linewidths=1, linecolor='w', alpha=alpha)
    if x_lim is None:
        x_lim = [0, len(data[0])]
    else:
        x_lim = [x_lim[0], x_lim[1]]
    if y_lim is None:
        y_lim = [0, len(data)]
    else:
        y_lim = [y_lim[0], y_lim[1]]
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    # format ticks and spines
    x_ticks = [t for t in range(x_lim[0], x_lim[1] + 1)]
    ax.set_xticks([t + 0.5 for t in x_ticks])
    x_ticklabels = [str(l) if l % 2 == 0 else '' for l in x_ticks]
    ax.set_xticklabels(x_ticklabels, rotation=0)
    y_ticks = [t for t in range(y_lim[0], y_lim[1] + 1)]
    ax.set_yticks([t + 0.5 for t in y_ticks])
    y_ticklabels = [str(l) if l % 2 == 0 else '' for l in y_ticks]
    ax.set_yticklabels(y_ticklabels, rotation=0)
    for position in ['right', 'left', 'top', 'bottom']:
        ax.spines[position].set_visible(True)
        ax.spines[position].set_color('k')
    # ax.spines['right'].set_visible(True)
    # ax.spines['left'].set_visible(True)
    # ax.spines['top'].set_visible(True)
    # ax.spines['bottom'].set_visible(True)
    # ax.spines['right'].set_color('k')
    # ax.spines['left'].set_color('k')
    # ax.spines['top'].set_color('k')
    # ax.spines['bottom'].set_color('k')
    ax.tick_params(axis='x', labelsize=tick_labelsize)
    ax.tick_params(axis='y', labelsize=tick_labelsize)
    plt.ylabel(y_label, size=labelsize)
    plt.xlabel(x_label, size=labelsize)
    ax.xaxis.set_tick_params(which='major', direction='out', length=3,
                             color='k', width=1, top='off', pad=pad)
    ax.yaxis.set_tick_params(which='major', direction='out', length=3,
                             color='k', width=1, top='off', pad=pad)

    plt.tight_layout()
    if figfile is None:
        plt.show()
    else:
        plt.savefig(figfile)
