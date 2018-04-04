#!/usr/bin/env python
# filename: plots.py


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

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns




def germline_plot(sequences, outfile=None, species='human',
        segment='V', chain='heavy', resolution='gene'):
    from abtools.germlines import germline_names
    resolution = resolution.lower()
    germs = germline_names(species.lower(), segment.upper(), chain=chain.lower(), resolution=resolution)
    if resolution not in ['gene', 'fam']:
        err = "ERROR: level must be either 'gene' or 'fam'"
        raise RuntimeError(err)
    if resolution == 'fam':
        keys = list(set([g.split('-')[0] for g in germs]))
        size = (6, 4)
    elif resolution == 'gene':
        keys = list(set([g.split('*')[0] for g in germs]))
        if segment == 'J':
            size = (6, 4)
        elif segment == 'D':
            size = (8, 4)
        else:
            size = (10, 4)
    gtype = '{}_gene'.format(segment.lower())
    _chain = chain[0].upper()
    data = [s[gtype][resolution].rstrip('D') for s in sequences if gtype in s]
    x, y = _aggregate(data, keys=keys)
    x = [i.replace('IG{}{}'.format(_chain, segment), '{}{}'.format(segment, _chain)) for i in x]
    colors = _get_germline_plot_colors(x, resolution)
    _make_barplot(x, y,
                 colors,
                 outfile,
                 ylabel='Frequency (%)',
                 l=resolution,
                 size=size)


def cdr3_plot(sequences, outfile=None, chain='heavy'):
    max_len = 40 if chain == 'heavy' else 20
    size = (9, 4) if chain == 'heavy' else (6, 4)
    chain = ['kappa', 'lambda'] if chain == 'light' else [chain, ]
    sequences = [s for s in sequences if s['chain'] in chain]
    cdr3s = [s['cdr3_len'] for s in sequences if s['cdr3_len'] > 0 and s['cdr3_len'] <= max_len]
    x, y = _aggregate(cdr3s, keys=list(range(1, max_len + 1)))
    color = sns.hls_palette(7)[4]
    x_title = 'CDR3 Length (AA)'
    y_title = 'Frequency (%)'
    _make_barplot([str(i) for i in x], y,
                 color,
                 outfile,
                 xlabel=x_title,
                 ylabel=y_title,
                 grid=True,
                 size=size,
                 xfontsize=7)


def vj_heatmap(seqs, outfile=None, chain='heavy'):
    # plot_file = '{}_VJheatmap_{}.pdf'.format(collection, args.chain)
    vj_data = _group_by_vj(seqs, chain)
    vj_df = pd.DataFrame(vj_data)
    _make_heatmap(vj_df.transpose(), outfile)


def _group_by_vj(data, chain):
    from abtools.germlines import germline_names
    vs = list(set([g.split('*')[0].replace('D', '') for g in germline_names('human', 'V', chain=chain)]))
    js = list(set([g.split('*')[0].replace('D', '') for g in germline_names('human', 'J', chain=chain)]))
    vj = {}
    for v in vs:
        vj[v] = {}
        for j in js:
            vj[v][j] = 0
    for d in data:
        v = d['v_gene']['gene'].replace('D', '')
        j = d['j_gene']['gene'].rstrip('P')
        if v[:3] != j[:3]:
            continue
        if v not in vj:
            continue
        if j not in vj[v]:
            continue
        vj[v][j] = vj[v][j] + 1 if j in vj[v] else 1
    total = len(data)
    for v in list(vj.keys()):
        for j in list(vj[v].keys()):
            vj[v][j] = 100. * vj[v][j] / total
    return vj


def _get_germline_plot_colors(data, l):
    fams = [d.split('-')[0].replace('/OR', '') for d in data]
    nr_fams = sorted(list(set(fams)))
    num_colors = len(nr_fams)
    rgbs = sns.hls_palette(num_colors)
    rgb_dict = {i[0]: i[1] for i in zip(nr_fams, rgbs)}
    return [rgb_dict[f] for f in fams]


def _make_barplot(x, y, colors, ofile=None, xlabel=None, ylabel=None,
        l=None, grid=False, size=None, xfontsize=None):
    sns.set_style('ticks')
    # set bar locations and width
    ind = np.arange(len(x))
    width = 0.75
    # plot objects
    if size:
        fig = plt.figure(figsize=size)
    else:
        fig = plt.figure()
    ax = fig.add_subplot(111)
    # axis limits and ticks
    ax.set_ylim(0, 1.05 * max(y))
    ax.set_xlim(-width / 2, len(ind))
    ax.set_xticks(ind + width / 2)
    xtick_names = ax.set_xticklabels(x)
    if l == 'gene':
        plt.setp(xtick_names, rotation=90, fontsize=7)
    if grid:
        ax.yaxis.grid(True, alpha=0.5)
    # axis labels
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xfontsize:
        plt.setp(xtick_names, fontsize=xfontsize)
    ax.tick_params(axis='x', which='both', top='off', length=3, pad=1.5)
    ax.tick_params(axis='y', which='both', right='off', length=3, pad=1.5)
    # make the plot
    bar = ax.bar(ind, y, width, color=colors)
    fig.tight_layout()
    if ofile is not None:
        plt.savefig(os.path.join(args.output, ofile))
        plt.close()
    else:
        plt.show()


def _make_heatmap(df, outfile):
    sns.set()
    # set up plot, determine plot size
    h, w = df.shape
    f, ax = plt.subplots(figsize=(w / 1.75, h / 3))
    sns.heatmap(df,
                square=True,
                cmap='Blues',
                cbar=True,
                cbar_kws={'orientation': 'horizontal',
                          'fraction': 0.02,
                          'pad': 0.02,
                          'shrink': 0.675})
    # adjust labels
    ax.xaxis.tick_top()
    plt.xticks(rotation=90)
    f.tight_layout()
    # make the plot
    if outfile is not None:
        plt.savefig(os.path.join(args.output, outfile))
        plt.close()
    else:
        plt.show()


def _aggregate(data, norm=True, sort_by='value', keys=None):
    '''
    Counts the number of occurances of each item in 'data'.

    Inputs
    data: a list of values.
    norm: normalize the resulting counts (as percent)
    sort_by: how to sort the retured data. Options are 'value' and 'count'.

    Output
    a non-redundant list of values (from 'data') and a list of counts.
    '''
    if keys:
        vdict = {k: 0 for k in keys}
        for d in data:
            if d in keys:
                vdict[d] += 1
    else:
        vdict = Counter(data)
    vals = [(k, v) for k, v in vdict.items()]
    if sort_by == 'value':
        vals.sort(key=lambda x: x[0])
    else:
        vals.sort(key=lambda x: x[1])
    xs = [v[0] for v in vals]
    if norm:
        raw_y = [v[1] for v in vals]
        total_y = sum(raw_y)
        ys = [100. * y / total_y for y in raw_y]
    else:
        ys = [v[1] for v in vals]
    return xs, ys
