#!/usr/bin/python
# filename: samples.py

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


import os

from abtools.utils import mongodb
from abtools.utils.pipeline import list_files


def assign_samples(db, collection, samplemap_dir):
	'''
	Assigns sample names (updates the provided MongoDB db/collection with
	a 'sample' field).

	Inputs:
		::db:: is a pymongo database object, containing the sequences
			to be updated
		::collection:: is a MongoDB collection name, as a string
		::platemap_dir:: is the directory containing one or more samplemap files
			(output from vaxtools.scpcr.platemap)
	'''
	samples = parse_samplemaps(samplemap_dir)
	update_samples(db, collection, samples)


def parse_samplemaps(samplemap_dir):
	samples = {}
	for samplemap in list_files(samplemap_dir):
		plate = os.path.basename(samplemap)
		with open(samplemap) as f:
			for line in f:
				sline = line.strip().split()
				if len(sline) >= 2:
					well = sline[0]
					sample = sline[1]
					platewell = '{}-{}'.format(plate, well)
					if sample not in samples:
						samples[sample] = []
					samples[sample].append(platewell)
	return samples


def update_samples(db, collection, samples):
	mongodb.remove_padding(db, collection)
	for sample in samples:
		seq_ids = samples[sample]
		match = {'seq_id': {'$in': seq_ids}}
		mongodb.update(field='sample', value=sample,
			db=db, collection=collection, match=match)



def assign_groups(db, collection, groupmap):
	'''
	Assigns group names (updates the provided MongoDB db/collection with
	a 'group' field).

	Inputs:
		::db:: is a pymongo database object, containing the sequences
			to be updated
		::collection:: is a MongoDB collection name, as a string
		::groupmap:: is a file containing group assignments, of the format:
				sample_name  group_name
			separated by any whitespace (one sample/group entry per line).
	'''
	groups = parse_groupmap(groupmap)
	update_groups(db, collection, groups)


def parse_groupmap(groupmap):
	groups = {}
	with open(groupmap) as f:
		for line in f:
			sline = line.strip().split()
			if sline:
				sample = sline[0]
				group = sline[1]
				groups[sample] = group
	return groups


def update_groups(db, collection, groups):
	mongodb.remove_padding(db, collection)
	samples = sorted(groups.keys())
	for sample in samples:
		group = groups.get(sample, None)
		if group is None:
			continue
		seq_ids = get_sample_seq_ids(db, collection, sample)
		match = {'seq_id': {'$in': seq_ids}}
		mongodb.update(field='group', value=group,
			db=db, collection=collection, match=match)


def get_sample_seq_ids(db, collection, sample):
	c = db[collection]
	seqs = c.find({'sample': sample}, {'seq_id': 1, '_id': 0})
	return [s['seq_id'] for s in seqs]
