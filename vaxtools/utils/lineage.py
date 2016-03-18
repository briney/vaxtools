#!/usr/bin/env python
# filename: lineage.py

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

import ete2

from abstar import run as run_abstar

from abtools.alignment import muscle
from abtools.phylogeny.utils import tree
from abtools.sequence import Sequence
from abtools.utils.properties import lazy_property

from vaxtools.utils.pair import Pair


class Lineage(object):
    '''
    Methods for manipulating an antibody lineage.


    INPUTS
    ------
    pairs: a list of one or more vaxtools.utils.pair.Pair objects


    PROPERTIES
    ----------
    name: the Clonify ID of the lineage (if Clonify was used to assign lineages).
        if ['clonify']['id'] does not exist in any of the heavy chains, name will
        be None.

    just_pairs: a list containing all lineage Pair objects that have both heavy
        and light chains.

    heavies: a list of all lineage Pair objects with a heavy chain (whether or not
        they also have a light chain)

    lights: a list of all lineage Pair objects with a light chain (whether or not
        they also have a heavy chain)

    uca: returns a Pair objecct of the unmutated common ancestor for the lineage.
        The uca is computed by taking the germline V(D)J regions, plus the N-addition
        region(s) from the least mutated heavy or light chain. If the lineage contains
        both heavy and light chains, the returned Pair will have ucas for both chains.

    subject: if any of the Pair objects contains a 'subject' property, this will be returned.
        If all 'subject' properties are the same, the return value will be a string. If there are
        multiple different 'subject' properties, the return value will be a list. If no Pairs
        have a 'subject' property, None will be returned.

    experiment: if any of the Pair objects contains a 'experiment' property, this will be returned.
        If all 'experiment' properties are the same, the return value will be a string. If there are
        multiple different 'experiment' properties, the return value will be a list. If no Pairs
        have a 'experiment' property, None will be returned.

    group: if any of the Pair objects contains a 'group' property, this will be returned.
        If all 'group' properties are the same, the return value will be a string. If there are
        multiple different 'group' properties, the return value will be a list. If no Pairs
        have a 'group' property, None will be returned.

    timepoint: if any of the Pair objects contains a 'timepoint' property, this will be returned.
        If all 'timepoint' properties are the same, the return value will be a string. If there are
        multiple different 'timepoint' properties, the return value will be a list. If no Pairs
        have a 'timepoint' property, None will be returned.


    METHODS
    -------

    Lineage.phylogeny()

        inputs
        ------
        project_dir: directory for alignment, tree and figure files (required)

        aln_file: if an alignment file has already been computed, passing the file path
            will force phylogeny() to use this alignment instead of computing a new one.

        tree_file: if a tree file has already been computed, passing the file path
            will force phylogeny() to use this tree instead of computing a new one.

        aa: if True, build an alignment based on amino acid sequences. Default is False.

        root: provide a sequence (string) to root the tree. Default is to use the UCA.

        colors: a dict that maps sequences to colors (as hex values). If used without
            orders or order_function, keys should be Pair names and values should be hex
            strings for all Pairs to be colored. Any Pair not present in the dict will
            be colored black. If used with orders or order_function, the keys should be
            orders and the values should be hex strings. As before, any order values
            not present in the dict will be colored black.

    '''
    def __init__(self, pairs, name=None):
        super(Lineage, self).__init__()
        self.pairs = pairs
        self.rmp_threshold = 0.9
        self.rmp_alt_allowed_mismatches = 1


    def __contains__(self, item):
        if item in self.pair_dict.keys():
            return True
        return False

    def __getitem__(self, key):
        return self.pair_dict.get(key, None)

    def __setitem__(self, key, val):
        self.pair_dict[key] = val


    @lazy_property
    def pair_dict(self):
        return {p.name: p for p in self.pairs}

    @lazy_property
    def name(self):
        clonify_ids = [p.heavy['clonify']['id'] for p in self.heavies if 'clonify' in p.heavy]
        if len(clonify_ids) > 0:
            return clonify_ids[0]
        return None

    @lazy_property
    def just_pairs(self):
        return [p for p in self.pairs if p.is_pair]

    @lazy_property
    def heavies(self):
        return [p for p in self.pairs if p.heavy is not None]

    @lazy_property
    def lights(self):
        return [p for p in self.pairs if p.light is not None]

    @lazy_property
    def uca(self):
        return self._calculate_uca()

    @lazy_property
    def rmp(self):
        return self._calculate_recalled_memory_precursor()

    @lazy_property
    def subject(self):
        return self._get_metadata('subject')

    @lazy_property
    def group(self):
        return self._get_metadata('group')

    @lazy_property
    def experiment(self):
        return self._get_metadata('experiment')

    @lazy_property
    def timepoints(self):
        if self.heavies:
            _timepoints = [p.heavy._mongo.get('timepoint', None) for p in self.heavies]
            return list(set([g for g in _timepoints if g is not None]))
        if self.lights:
            _timepoints = [p.light._mongo.get('timepoint', None) for p in self.lights]
            return list(set([g for g in _timepoints if g is not None]))
        return []


    def dot_alignment(self, seq_field='vdj_nt', name_field='seq_id', chain='heavy'):
        if chain == 'heavy':
            sequences = [p.heavy for p in self.heavies if seq_field in p.heavy]
            sequences.append(self.uca.heavy)
            seqs = [(s[name_field], s[seq_field]) for s in sequences]
        else:
            sequences = [p.light for p in self.lights if seq_field in p.light]
            sequences.append(self.uca.light)
            seqs = [(s[name_field], s[seq_field]) for s in sequences]
        aln = muscle(seqs)
        g_aln = [a for a in aln if a.id == 'UCA'][0]
        dots = [('UCA', str(g_aln.seq)), ]
        for seq in [a for a in aln if a.id != 'UCA']:
            s_aln = ''
            for g, q in zip(str(g_aln.seq), str(seq.seq)):
                if g == q == '-':
                    s_aln += '-'
                elif g == q:
                    s_aln += '.'
                else:
                    s_aln += q
            dots.append((seq.id, s_aln))
        name_len = max([len(d[0]) for d in dots]) + 2
        dot_aln = []
        for d in dots:
            spaces = name_len - len(d[0])
            dot_aln.append(d[0] + ' ' * spaces + d[1])
        return '\n'.join(dot_aln)


    def phylogeny(self, project_dir, aln_file=None, tree_file=None, aa=False,
            root=None, colors=None, color_function=None, orders=None, order_function=None,
            chain='heavy', filter_function=None, just_pairs=False,
            scale=None, branch_vert_margin=None, fontsize=12, show_names=True,
            mirror=False, min_order_fraction=0.1):
        seq_field = 'vdj_aa' if aa else 'vdj_nt'
        project_dir = os.path.abspath(project_dir)
        orientation = 1 if mirror else 0
        # get sequences for phylogeny
        if chain == 'heavy':
            seq_pool = self.just_pairs if just_pairs else self.heavies
            seqs = [p.heavy for p in seq_pool]
            if root is None:
                root = self.uca.heavy[seq_field]
            seqs += [Sequence(root, id='root')]
        else:
            seq_pool = self.just_pairs if just_pairs else self.lights
            seqs = [p.light for p in seq_pool]
            if root is None:
                root = self.uca.light[seq_field]
            seqs += [Sequence(root, id='root')]
        # filter sequences
        if filter_function is not None:
            seqs = [s for s in seqs if filter_function(s)]
        # setup colors
        if colors is None:
            if color_function is not None:
                colors = {seq['seq_id']: color_function(seq) for seq in seqs}
            else:
                colors = {}
        # make msa
        if all([aln_file is None, tree_file is None]):
            aln_file = os.path.join(project_dir, '{}.aln'.format(self.name))
            muscle(seqs, aln_file, as_file=True)
        # make treefile
        if tree_file is None:
            tree_file = os.path.join(project_dir, '{}.nw'.format(self.name))
            tree.fast_tree(aln_file, tree_file, is_aa=aa)
        # make phylogeny
        fig_file = os.path.join(project_dir, '{}.pdf'.format(self.name))
        self._make_tree_figure(tree_file, fig_file, colors, orders,
            show_names=show_names, branch_vert_margin=branch_vert_margin, scale=scale,
            tree_orientation=orientation, fontsize=fontsize, min_order_fraction=min_order_fraction)


    def _get_metadata(self, field):
        _metadata = None
        if self.heavies:
            _metadata = [p.heavy._mongo.get(field, None) for p in self.heavies]
            _metadata = [g for g in _metadata if g is not None]
            if len(list(set(_metadata))) == 1:
                metadata = _metadata[0]
            elif len(list(set(_metadata))) > 1:
                mcount = Counter(_metadata)
                metadata = sorted([s for s in mcount.keys()],
                                key=lambda x: mcount[s],
                                reverse=True)[0]
        if self.lights and metadata is None:
            _metadata = [p.light._mongo.get(field, None) for p in self.lights]
            _metadata = [g for g in _metadata if g is not None]
            if len(list(set(_metadata))) == 1:
                metadata = _metadata[0]
            elif len(list(set(_metadata))) > 1:
                mcount = Counter(_metadata)
                metadata = sorted([s for s in mcount.keys()],
                                key=lambda x: mcount[s],
                                reverse=True)[0]
        return metadata


    def _make_tree_figure(self, tree, fig, colors, orders, scale=None, branch_vert_margin=None,
            fontsize=12, show_names=True, tree_orientation=0, min_order_fraction=0.1):
        if show_names is True:
            show_names = [p.name for p in self.pairs]
        elif show_names is False:
            show_names = []
        t = ete2.Tree(tree)
        t.set_outgroup(t&"root")
        # style the nodes
        for node in t.traverse():
            if orders is not None:
                leaves = node.get_leaf_names()
                order_count = Counter([orders[l] for l in leaves])
                for order in sorted(order_count.keys()):
                    if float(order_count[order]) / len(leaves) >= min_order_fraction:
                        color = colors[order]
                        break
            else:
                color = colors.get(node.name, '#000000')
            style = ete2.NodeStyle()
            style['size'] = 0
            style['vt_line_width'] = 1.0
            style['hz_line_width'] = 1.0
            style['vt_line_color'] = color
            style['hz_line_color'] = color
            style['vt_line_type'] = 0
            style['hz_line_type'] = 0
            if node.name in show_names:
                tf = ete2.TextFace(node.name)
                tf.fsize = fontsize
                node.add_face(tf, column=0)
                style['fgcolor'] = '#000000'
            node.set_style(style)
        t.dist = 0
        ts = ete2.TreeStyle()
        ts.orientation = tree_orientation
        ts.show_leaf_name = False
        if scale:
            ts.scale = int(scale)
        if branch_vert_margin:
            ts.branch_vertical_margin = float(branch_vert_margin)
        ts.show_scale = False
        # ladderize
        t.ladderize()
        # render the tree
        t.render(fig, tree_style=ts)


    def _calculate_uca(self, paired_only=False):
        if paired_only:
            heavies = self.just_pairs
            lights = self.just_pairs
        else:
            heavies = self.heavies
            lights = self.lights
        # heavy chain UCA
        if len(heavies) >= 1:
            lmhc = sorted(heavies, key=lambda x: x.heavy['nt_identity']['v'], reverse=True)[0].heavy
            hc_uca = run_abstar(('UCA', lmhc['vdj_germ_nt']))
        else:
            lmhc = None
        # light chain UCA
        if len(lights) >= 1:
            lmlc = sorted(lights, key=lambda x: x.light['nt_identity']['v'], reverse=True)[0].light
            lc_uca = run_abstar(('UCA', lmlc['vdj_germ_nt']))
        else:
            lmlc = None
        return Pair([hc_uca, lc_uca])


    def _calculate_recalled_memory_precursor(self, paired_only=False):
        if paired_only:
            heavies = [p.heavy for p in self.just_pairs]
            lights = [p.light for p in self.just_pairs]
        else:
            heavies = [p.heavy for p in self.heavies]
            lights = [p.light for p in self.lights]
        if len(heavies) >= 1:
            hc_rmp = self._rmp(heavies + [self.uca.heavy])
        else:
            hc_rmp = None
        if len(lights) >= 1:
            lc_rmp = self._rmp(lights + [self.uca.light])
        else:
            lc_rmp = None
        return Pair([hc_rmp, lc_rmp])


    def _rmp(self, sequences):
        rmp = ''
        seqs = [(s['seq_id'], s['vdj_nt']) for s in sequences]
        aln = muscle(seqs)
        g_aln = [a for a in aln if a.id == 'UCA'][0]
        query_seqs = [str(a.seq) for a in aln if a.id != 'UCA']
        for i, g in enumerate(g_aln):
            qcounts = Counter([q[i] for q in query_seqs])
            qmax = sorted(qcounts.keys(),
                          key=lambda x: qcounts[x],
                          reverse=True)[0]
            qmax_fraction = float(qcounts[qmax]) / sum(qcounts.values())
            qmax_alt_mismatches = sum(qcounts.values()) - qcounts[qmax]
            if any([qmax_fraction >= self.rmp_threshold,
                    qmax_alt_mismatches <= self.rmp_alt_allowed_mismatches]):
                rmp += qmax
            else:
                rmp += g
        return run_abstar(Sequence(rmp.replace('-', ''), id='RMP'))


    def _germline_field_map(self, field):
        fmap = {'fr1_nt': 'fr1_germ_nt',
                'fr2_nt': 'fr2_germ_nt',
                'fr3_nt': 'fr3_germ_nt',
                'fr4_nt': 'fr4_germ_nt',
                'cdr1_nt': 'cdr1_germ_nt',
                'cdr2_nt': 'cdr2_germ_nt',
                'fr1_aa': 'fr1_germ_aa',
                'fr2_aa': 'fr2_germ_aa',
                'fr3_aa': 'fr3_germ_aa',
                'fr4_aa': 'fr4_germ_aa',
                'cdr1_aa': 'cdr1_germ_aa',
                'cdr2_aa': 'cdr2_germ_aa',
                'vdj_nt': 'vdj_germ_nt',
                'vdj_aa': 'vdj_germ_aa'}
        return fmap.get(field.lower(), None)
