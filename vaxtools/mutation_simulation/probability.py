#!/usr/bin/env python
# filename: probability.py


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

from collections import Counter
from multiprocessing import Pool
import os
import sys

from abtools.alignment import global_alignment, local_alignment
from abtools.jobs import monitor_mp_jobs


def mutation_probabilities(collections, standard, outfile=None, query=None,
        seq_field='vdj_aa', chunksize=100):
    '''
    Calculates the probability of each mutation from the standard sequence,
    given one or more collections.

    Inputs
    ------
    collections: One or more pymongo collection objects, as an iterable.
    standard: The target amino acid sequence, as a string.
    query: Query parameters, as a dict. Will be passed directly to
        collection.find()
    seq_field: The MongoDB field to be used for comparison to the target.
        Default is 'vdj_aa'.
    chunksize: Number of sequences to be submitted for each alignment job.
        Default is 100.

    Returns
    -------
    A dictionary of normalized mutation probabilities, of the format:
        {'12A': 0.01, '15F': 0.12, ...}
    Mutation names are a concatenation of the mutation position (1-based
    indexing of the standard sequence) and the mutated residue.
    '''
    async_results = []
    p = Pool()
    for collection in collections:
        print('\n' + collection)
        print('-' * len(collection))
        print('querying for sequences...')
        sequences = get_sequences(collection, query=query, seq_field=seq_field)
        print('performing alignments:')
        for chunk in chunker(sequences, chunksize):
            async_results.append(p.apply_async(get_mutations, args=(chunk, standard)))
        monitor_mp_jobs(async_results)
    mutations = []
    for ar in async_results:
        mutations.extend(ar.get())
    print('\ncalculating mutation probabilities...')
    mcounts = Counter(mutations)
    total = sum(mcounts.values())
    norm_counts = {k: float(v) / total for k, v in mcounts.items()}
    prob_string = '\n'.join(['{}\t{}'.format(k, v) for k, v in list(norm_counts.items())])
    if outfile is not None:
        open(outfile, 'w').write(prob_string)
    return norm_counts


def chunker(l, n):
    'Generator that produces n-length chunks from iterable l.'
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_sequences(db, collection, query=None, seq_field='vdj_aa'):
    c = db[collection]
    if query is None:
        query = {}
    proj = {'_id': 0, seq_field: 1}
    results = c.find(query, proj)
    [r[seq_field] for r in results]


def get_mutations(seqs, standard):
    mutations = []
    for seq in seqs:
        aln = global_alignment(seq, target=standard,
            matrix='blosum62', gap_open=-15, gap_extend=-1)
        mutations.extend(parse_mutations(aln))
    return mutations


def parse_mutations(aln):
    muts = []
    tpos = 0
    for q, t in zip(aln.aligned_query, aln.aligned_target):
        if t == '-':
            continue
        tpos += 1
        if any([q == '-', q == t]):
            continue
        mut = '{}{}'.format(tpos, q)
        muts.append(mut)
    return muts
