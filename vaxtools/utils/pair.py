#!/usr/bin/env python
# filename: pair.py


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


class PairHL(object):
	'''
	Holds a pair of sequences, corresponding to HC and LC of a single mAb.

	Input is a list of dicts, with each dict containing sequence information from a single
	chain, formatted as would be returned from a query on a MongoDB database containing
	AbStar output.
	'''
	def __init__(self, seqs, name=None):
		self._seqs = seqs
		self._heavy = None
		self._light = None
		self._name = name
		self._sample = None
		self._group = None
		self._is_pair = None
		self._vrc01_like = None
		self._lineage = None

	@property
	def heavy(self):
		if self._heavy is None:
			h_chains = [s for s in self._seqs if s['chain'] == 'heavy']
			if len(h_chains) > 0:
				self._heavy = h_chains[0]
			else:
				self._heavy = None
		return self._heavy

	@property
	def light(self):
		if self._light is None:
			l_chains = [s for s in self._seqs if s['chain'] in ['kappa', 'lambda']]
			if len(l_chains) > 0:
				self._light = l_chains[0]
			else:
				self._light = None
		return self._light

	@property
	def is_pair(self):
		if all([self.heavy is not None, self.light is not None]):
			return True
		return False

	@property
	def lineage(self):
		if self._lineage is None:
			self._lineage = self.heavy['clonify']['id']
		return self._lineage

	@property
	def vrc01_like(self):
		if self._vrc01_like is None:
			self._vrc01_like = all([self.heavy['v_gene']['gene'] == 'IGHV1-2', self.light['cdr3_len'] == 5])
		return self._vrc01_like

	@property
	def name(self):
		return self._name

	@name.setter
	def name(self, name):
		self._name = name

	@property
	def sample(self):
		if self._sample is None:
			if self.heavy is not None and 'sample' in self.heavy.keys():
				self._sample = self.heavy['sample']
			elif self.light is not None and 'sample' in self.light.keys():
				self._sample = self.light['sample']
		return self._sample

	@sample.setter
	def sample(self, sample):
		self._sample = sample

	@property
	def group(self):
		if self._group is None:
			if self.heavy is not None and 'group' in self.heavy.keys():
				self._group = self.heavy['group']
			elif self.light is not None and 'group' in self.light.keys():
				self._group = self.light['group']
		return self._group

	@group.setter
	def group(self, group):
		self._group = group


def get_pairs(db, collection, sample=None, group=None, name='seq_id',
	delim=None, delim_occurance=1, pairs_only=False):
	'''
	Gets seuqences and assigns them to the appropriate mAb pair, based on the sequence name.

	Inputs:

	::db:: is a pymongo database connection object
	::collection:: is the collection name, as a string
	If ::sample:: is provided, only sequences with a 'sample' field matching ::sample:: will
		be included. ::sample:: can be either a single sample (as a string) or an iterable
		(list or tuple) of sample strings.
	If ::group:: is provided, only sequences with a 'group' field matching ::group:: will
		be included. ::group:: can be either a single group (as a string) or an iterable
		(list or tuple) of group strings.
	::name:: is the dict key of the field to be used to group the sequences into pairs.
		Default is 'seq_id'
	::delim:: is an optional delimiter used to truncate the contents of the ::name:: field.
		Default is None, which results in no name truncation.
	::delim_occurance:: is the occurance of the delimiter at which to trim. Trimming is performed
		as delim.join(name.split(delim)[:delim_occurance]), so setting delim_occurance to -1 will
		trucate after the last occurance of delim. Default is 1.
	::pairs_only:: setting to True results in only truly paired sequences (pair.is_pair == True)
		will be returned. Default is False.

	Returns a list of PairHL objects, one for each mAb pair.
	'''
	match = {}
	if sample is not None:
		if type(sample) in (list, tuple):
			match['sample'] = {'$in': sample}
		elif type(sample) in (str, unicode):
			match['sample'] = sample
	if group is not None:
		if type(group) in (list, tuple):
			match['group'] = {'$in': group}
		elif type(group) in (str, unicode):
			match['group'] = group
	seqs = list(db[collection].find(match))
	return assign_pairs(seqs, name=name, delim=delim,
		delim_occurance=delim_occurance, pairs_only=pairs_only)


def assign_pairs(seqs, name='seq_id', delim=None, delim_occurance=1, pairs_only=False):
	'''
	Assigns sequences to the appropriate mAb pair, based on the sequence name.

	Inputs:

	::seqs:: is a list of dicts, of the format returned by querying a MongoDB containing
		Abstar output.
	::name:: is the dict key of the field to be used to group the sequences into pairs.
		Default is 'seq_id'
	::delim:: is an optional delimiter used to truncate the contents of the ::name:: field.
		Default is None, which results in no name truncation.
	::delim_occurance:: is the occurance of the delimiter at which to trim. Trimming is performed
		as delim.join(name.split(delim)[:delim_occurance]), so setting delim_occurance to -1 will
		trucate after the last occurance of delim. Default is 1.
	::pairs_only:: setting to True results in only truly paired sequences (pair.is_pair == True)
		will be returned. Default is False.

	Returns a list of PairHL objects, one for each mAb pair.
	'''
	pdict = {}
	for s in seqs:
		if delim is not None:
			pname = delim.join(s[name].split(delim)[:delim_occurance])
		else:
			pname = s[name]
		if pname not in pdict:
			pdict[pname] = [s, ]
		else:
			pdict[pname].append(s)
	pairs = [PairHL(pdict[n], name=n) for n in pdict.keys()]
	if pairs_only:
		pairs = [p for p in pairs if p.is_pair]
	return pairs
