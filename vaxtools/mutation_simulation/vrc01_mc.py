#!/usr/bin/env python
# filename: vrc01_mc.py


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
import sys
from datetime import datetime

import numpy as np
from scipy import stats

from abtools.alignment import global_alignment
from abtools.utils.progbar import progress_bar

from vaxtools.mutation_simulation.probability import get_mutations
from vaxtools.utils.vrc01 import get_vrc01_class_sequences, get_vrc01_germline_sequence


def simulate(probabilities, n_mutations=1, n_sequences=50000):
    '''
    Simulates mutation of antibody sequences, given a list of mutations
    and their probabilities.

    Inputs
    ------
    probabilities: a dictionary containing mutations as keys and probabilities
        as values.
    n_mutations: Number of mutations in each sequence. If n_mutations is
        greater than 1, mutations are selected from the pool of supplied
        mutations without replacement. Default is 1.
    n_sequences: Number of mutated sequences to generate. Default is 50,000.

    Returns
    -------
    Mean number of VRC01-like mutations (float) and 95% confidence interval
    (tuple of floats).
    '''
    vrc01_muts = get_vrc01_class_mutations()
    sim_muts = []
    muts, probs = list(zip(*list(probabilities.items())))
    start = datetime.now()
    for i in range(n_sequences):
        m = np.random.choice(muts, size=n_mutations, replace=False, p=probs)
        sim_muts.append(m)
        if (i + 1) % 100 == 0:
            progress_bar(i + 1, n_sequences, start)
    mut_counts = [len([x for x in sublist if x in vrc01_muts]) for sublist in sim_muts]
    # calculate mean and 95% confidence interval
    n, min_max, mean, var, skew, kurt = stats.describe(mut_counts)
    std = np.sqrt(var)
    R = stats.norm.interval(0.95, loc=mean, scale=std / np.sqrt(len(mut_counts)))
    return mean, R


def get_vrc01_class_mutations():
    vrc01_class = [s.sequence for s in get_vrc01_class_sequences()]
    glvrc01 = get_vrc01_germline_sequence().sequence
    return list(set(get_mutations(vrc01_class, glvrc01)))



# list intersection = http://stackoverflow.com/questions/642763/python-intersection-of-two-lists
