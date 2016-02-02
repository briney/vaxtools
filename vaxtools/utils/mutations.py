#!/usr/bin/env python
# filename: mutations.py


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


from vaxtools.utils.codon_tables import translate


def single_nucleotide_aa_mutations(seq, v=True, j=True):
    '''
    Analyses codons to determine the number of AA mutations that requre
    single nucleotide changes.

    By default, both V and J genes are analysed. Either can be disabled
    by setting ::v:: or ::j:: to False.
    '''
    single = 0
    if v:
        vquery = seq['codons']['v']
        vgerm = seq['codons']['v_germ']
        for q, g in zip(vquery, vgerm):
            diffs = _codon_nt_changes(q, g)
            if diffs == 1:
                single += 1
    if j:
        jquery = seq['codons']['j']
        jgerm = seq['codons']['j_germ']
        for q, g in zip(jquery, jgerm):
            diffs = _codon_nt_changes(q, g)
            if diffs == 1:
                single += 1
    return single


def multiple_nucleotide_aa_mutations(seq, v=True, j=True):
    '''
    Analyses codons to determine the number of AA mutations that requre
    multiple nucleotide changes.

    By default, both V and J genes are analysed. Either can be disabled
    by setting ::v:: or ::j:: to False.
    '''
    multiple = 0
    if v:
        vquery = seq['codons']['v']
        vgerm = seq['codons']['v_germ']
        for q, g in zip(vquery, vgerm):
            diffs = _codon_nt_changes(q, g)
            if diffs > 1:
                multiple += 1
    if j:
        jquery = seq['codons']['j']
        jgerm = seq['codons']['j_germ']
        for q, g in zip(jquery, jgerm):
            diffs = _codon_nt_changes(q, g)
            if diffs > 1:
                multiple += 1
    return multiple


def _codon_nt_changes(q, g):
    if any(['-' in q, '-' in g]):
        return 0
    qaa = translate(q)
    gaa = translate(g)
    if q == g:
        return 0
    diffs = 0
    for _q, _g in zip(q, g):
        if _q != _g:
            diffs += 1
    return diffs



def silent_mutations(seq, v=True, j=True):
    silent = 0
    if v:
        vquery = seq['codons']['v']
        vgerm = seq['codons']['v_germ']
        for q, g in zip(vquery, vgerm):
            silent += _is_silent_mutation(q, g)
    if j:
        jquery = seq['codons']['j']
        jgerm = seq['codons']['j_germ']
        for q, g in zip(jquery, jgerm):
            silent += _is_silent_mutation(q, g)
    return silent


def nonsilent_mutations(seq, v=True, j=True):
    nonsilent = 0
    if v:
        vquery = seq['codons']['v']
        vgerm = seq['codons']['v_germ']
        for q, g in zip(vquery, vgerm):
            nonsilent += _is_nonsilent_mutation(q, g)
    if j:
        jquery = seq['codons']['j']
        jgerm = seq['codons']['j_germ']
        for q, g in zip(jquery, jgerm):
            nonsilent += _is_nonsilent_mutation(q, g)
    return nonsilent


def _is_silent_mutation(q, g):
    if q == g:
        return 0
    if any(['-' in q, '-' in g]):
        return 0
    qaa = translate(q)
    gaa = translate(g)
    if qaa != gaa:
        return 0
    if q != g:
        return 1


def _is_nonsilent_mutation(q, g):
    if q == g:
        return 0
    if any(['-' in q, '-' in g]):
        return 0
    qaa = translate(q)
    gaa = translate(g)
    if qaa == gaa:
        return 0
    return 1
