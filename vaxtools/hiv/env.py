#!/usr/bin/env python
# filename: env.py


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
import multiprocessing as mp
import os

from Bio import SeqIO

from abtools.alignment import global_alignment
from abtools.sequence import Sequence



CLADES = ['A', 'AC', 'AD', 'AE', 'AG', 'B', 'BC', 'BF', 'C', 'D', 'E', 'F', 'G', 'OTHER']
HXB2 = ('HXB2', 'MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLNCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYTLTSCNTSVISQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFMDNAKTIIVQLNTSVEINCTRPSNNTIKRIRIQRGPGRAFVTMGKIGDMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGKGNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHMTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFDITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL')


class Env(object):
    """docstring for Env"""
    def __init__(self, sequence, reference=HXB2):
        super(Env, self).__init__()
        self._sequence = Sequence(sequence)
        self._reference = Sequence(reference)
        self._aligned = None
        self._aligned_reference = None
        self._id = None
        self._reference_name = None
        self._clade = None


    def __getitem__(self, key):
        if isinstance(key, int):
            pos = self._get_aligned_position(key)
            return self.aligned[pos]
        elif isinstance(key, slice):
            start = self._get_aligned_position(key.start)
            stop = self._get_aligned_position(key.stop)
            return self.aligned[start:stop:key.step]


    @property
    def sequence(self):
        return self._sequence


    @property
    def id(self):
        if self._id is None:
            self._id = self.sequence.id
        return self._id


    @property
    def reference(self):
        return self._reference


    @property
    def reference_name(self):
        if self._reference_name is None:
            self._reference_name = self.reference.id
        return self._reference_name


    @property
    def aligned(self):
        if self._aligned is None:
            self._align_to_reference()
        return self._aligned


    @property
    def aligned_reference(self):
        if self._aligned_reference is None:
            self._align_to_reference()
        return self._aligned_reference


    @property
    def clade(self):
        if self._clade is None:
            raw_clade = ''.join(i for i in self.id.split('.')[0] if i.isalpha())
            self._clade = raw_clade if raw_clade in CLADES else 'OTHER'
        return self._clade


    def glycan(self, position):
        end = position + 3
        codon = self[position:end]
        if '-' in codon:
            codon = codon.replace('-', '')
            while len(codon) < 3:
                end += 1
                codon = self[position:end].replace('-', '')
        try:
            if all([codon[0].upper() == 'N', codon[1].upper() != 'P', codon[2].upper() in ['S', 'T']]):
                return True
            return False
        except IndexError:
            return False


    def region(positions):
        pass


    def all_glycans(self):
        glycans = []
        for i in range(1, len(self.reference.sequence) - 3):
            if self.glycan(i):
                glycans.append(i)
        return glycans


    def _get_aligned_position(self, pos):
        if pos is None:
            return pos
        raw = 0
        ref = 0
        riter = iter(self.aligned_reference)
        while ref < pos:
            try:
                p = riter.next()
                raw += 1
                if p != '-':
                    ref += 1
            except StopIteration:
                break
        return raw - 1


    def _align_to_reference(self):
        aln = global_alignment(self.sequence, self.reference, aa=True)
        self._aligned = aln.aligned_query
        self._aligned_reference = aln.aligned_target


def get_clades():
    return CLADES


def get_env_sequences(clade=None):
    if clade is None:
        clade = CLADES
    elif type(clade) in [str, unicode]:
        clade = [clade, ]

    envs = []
    seq_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'env_sequences')

    for c in clade:
        seq_file = os.path.join(seq_dir, '{}.fasta'.format(c.upper()))
        envs.extend([Sequence([s.id, str(s.seq).upper()]) for s in SeqIO.parse(open(seq_file, 'r'), 'fasta')])
    return envs


def get_all_hxb2_glycan_positions():
    env = Env(HXB2)
    return env.all_glycans()


def get_all_glycan_positions(envs=None, clade=None, minimum_frequency=0.01):
    if envs is None:
        envs = [Env(seq) for seq in get_env_sequences(clade)]
    glycans = []
    for e in envs:
        glycans.extend(_find_all_glycans(e))
    glycount = Counter(glycans)
    glycan_positions = [g for g, c in glycount.items() if float(c) / len(envs) >= minimum_frequency]
    return sorted(glycan_positions)


def force_reference_alignment(envs):
    for e in envs:
        e._align_to_reference()


def parse_envs_into_clades(env_file, output_directory, clades=CLADES):
    '''
    Parses Env sequences from the Los Alamos database into separate files for
    each clade.

    Inputs:
      env_file - the combined FASTA file of Env sequences, downloaded from LANL
      output_directory - directory into which the parsed files should be saved
      clades - a list of clades (strings). By default, CLADES will be used. All
               sequences that don't match a clade will be designated 'OTHER'.
    '''
    envs_by_clade = {c: [] for c in clades}
    envs_by_clade['OTHER'] = []
    ehandle = open(env_file, 'r')
    raw_clade_names = []
    for s in SeqIO.parse(ehandle, 'fasta'):
        raw_clade = ''.join(i for i in s.id.split('.')[0] if i.isalpha())
        clade = raw_clade if raw_clade in clades else 'OTHER'
        envs_by_clade[clade].append(s)
    ehandle.close()
    for c in envs_by_clade.keys():
        ofile = os.path.join(output_directory, '{}.fasta'.format(c))
        fastas = ['>{}\n{}'.format(s.id, str(s.seq)) for s in envs_by_clade[c]]
        open(ofile, 'w').write('\n'.join(fastas))


def _find_all_glycans(e):
    glycans = []
    ref = 0
    for i, (a, r) in enumerate(zip(e.aligned, e.aligned_reference)):
        if r != '-':
            ref += 1
        if a == 'N':
            try:
                codon = a
                adder = 1
                while len(codon) < 3:
                    res = e.aligned[i + adder]
                    if res != '-':
                        codon += res
                    adder += 1
                if all([codon[1] != 'P', codon[2] in ['S', 'T']]):
                    glycans.append(ref)
            except IndexError:
                continue
    return glycans


def _force_alignment(e):
    e._align_to_reference()
