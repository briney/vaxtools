#!/usr/bin/env python
# filename: demultiplex.py


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

import glob
import logging
import math
import os
import sqlite3
import string
from StringIO import StringIO
import subprocess as sp
import sys
import tempfile
import time
import traceback
import urllib
import uuid

from pymongo import MongoClient

import numpy as np
import pandas as pd

from Bio import AlignIO
from Bio.Align import AlignInfo

# import matplotlib
# # Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from vaxtools.utils import pixel

from abtools.utils import log, mongodb
from abtools.utils.alignment import mafft
from abtools.utils.pipeline import make_dir


def parse_args():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-o', '--output', dest='output', required=True,
						help="Output file for demultiplexed FASTAs. Required. \
						If the parent folder doesn't exist, it will be created.")
	parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
						help="Directory for temporary files. Required. \
						If the directory doesn't exist, it will be created.")
	parser.add_argument('-r', '--raw-sequences', dest='raw_sequence_dir', default=None,
						help="Directory for saving the raw sequence data (not just consensus/centroid sequences. \
						If not provided, raw binned sequence data will not be retained. \
						If the directory doesn't exist, it will be created.")
	parser.add_argument('-a', '--alignment-pixel-dir', dest='alignment_pixel_dir', default=None,
						help="Directory for saving the alignment pixel figures. \
						If not provided, figures will not be generated. \
						If the directory doesn't exist, it will be created.")
	parser.add_argument('-l', '--log', dest='log',
						help="Location for the log file. \
						Default is <output_dir>/demultiplex.log.")
	parser.add_argument('-d', '--database', dest='db', required=True,
						help="Name of the MongoDB database containing un-demultiplexed sequences. Required.")
	parser.add_argument('-c', '--collection', dest='collection', default=None,
						help="Name of the MongoDB collection to query. \
						If not provided, all collections in the given database will be processed iteratively.")
	parser.add_argument('--collection-prefix', dest='collection_prefix', default=None,
						help="Only collections starting with <collection-prefix> will be processed. \
						Useful if a single MongoDB database contains collections from multiple experiments \
						and the collections from each experiment should be processed separately. \
						Default is None, which (unless <collection> or <collection-suffix> is provided) \
						processes all collections.")
	parser.add_argument('--collection-suffix', dest='collection_suffix', default=None,
						help="Only collections ending with <collection-suffix> will be processed. \
						Useful if a single MongoDB database contains collections from multiple experiments \
						and the collections from each experiment should be processed separately. \
						Default is None, which (unless <collection> or <collection-prefix>is provided) \
						processes all collections.")
	parser.add_argument('-i', '--ip', dest='ip', default='localhost',
						help="IP address for the MongoDB server.  Defaults to 'localhost'.")
	parser.add_argument('--port', dest='port', default=27017, type=int,
						help="Port for the MongoDB server.  Defaults to '27017'.")
	parser.add_argument('-u', '--user', dest='user', default=None,
						help="Username for the MongoDB server. Not used if not provided.")
	parser.add_argument('-p', '--password', dest='password', default=None,
						help="Password for the MongoDB server. Not used if not provided.")
	parser.add_argument('-I', '--index', dest='index', default=None,
						help="Type of built-in indexes to use. Choices are: 'top-odd', 'top-even', \
						'bottom-odd', 'bottom-even', or '384'. \
						If a custom index file is desired (instead of the built-in indexes), use \
						--index-file. Default is None, which doesn't use built-in indexes.")
	parser.add_argument('--index-file', dest='index_file', default=None,
						help="File containing the indexes, one per line. \
						File must have either 96 or 384 indexes, and indexes should be in columnar order \
						(A01, B01, C01, etc...). Only required if not using built-in indexes (via --index).")
	parser.add_argument('-M', '--plate-map', dest='plate_map', required=True,
						help="Plate map. Can either be provided as a list of well names \
						(one per line, in row order starting at A01) or as a list of wells \
						and well names (one per line, separated by white space) in any order. \
						If providing simply a list of well names and not every well is used, \
						leave blank lines for unused wells.")
	parser.add_argument('--index-position', default='start', choices=['start', 'end'],
						help="Position of the indexes. Choices are 'start', if they're \
						at the start of the raw merged read (start of read 1) or 'end' if \
						they're at the end of the raw merged read (start of read 2).")
	parser.add_argument('--index-length', default=0, type=int,
						help="Length of the index, in nucleotides. \
						Default is to parse the index length from the file of index sequences.")
	parser.add_argument('--index-reverse-complement', default=False, action='store_true',
						help="Set if the indexes in the supplied index file are the reverse \
						complement of the indexes as they will appear in the sequences. \
						Default is False.")
	parser.add_argument('--score-cutoff-heavy', default=200, type=int,
						help="V-gene alignment score cutoff for heavy chains. \
						Alignment score must be equal to or higher than cutoff to be considered for clustering. \
						Default is 200.")
	parser.add_argument('--score-cutoff-light', default=100, type=int,
						help="V-gene alignment score cutoff for kappa/lambda chains. \
						Alignment score must be equal to or higher than cutoff to be considered for clustering. \
						Default is 100.")
	parser.add_argument('--cdhit-threshold', default=0.96, type=float,
						help="Threshold for CD-HIT clustering. \
						Default is 0.96.")
	parser.add_argument('--minimum-well-size', default='relative',
						help="Minimum size of a CD-HIT cluster of sequences. \
						Centroids will not be determined for clusters below this cutoff. \
						If set to 'relative', the cutoff will be a fraction of the total number of reads \
						for the entire plate. \
						Default is 'relative'.")
	parser.add_argument('--minimum-max-well-size', default=250, type=int,
						help="If --minimum-well-size is set to relative, this sets the minimum number of \
						sequences in the largest well on the plate. \
						Without this, a plate with only 1 or 2 sequences per well would be processed, \
						which would produce some strange (and likely wrong) results. \
						Default is 250.")
	parser.add_argument('--minimum-cluster-fraction', default='largest',
						help="Minimum fraction for a CD-HIT cluster (relative to the total \
						number of sequences in a given well) for a centroid to be determined. \
						For example, if set to 0.7, the largest CD-HIT cluster for a well must comprise \
						70 percent of all sequences in that well. \
						If set to 'largest', the largest cluster is selected, regardless of size relative \
						to the total number of sequences in the well. \
						Default is 'largest'.")
	parser.add_argument('--minimum-well-size-denom', default=96, type=int,
						help="If --minimum-well-size is set to 'relative', the cutoff will be a fraction \
						of the total number of sequences in the plate. \
						This is the denominator for that fraction. \
						Default is 96.")
	parser.add_argument('--cluster-cutoff-gradient', default=False, action='store_true',
						help="Allows for a gradiated cluster cutoff, with the percent of high-homology \
						sequences within each well varying based on the number of sequences in the well. \
						Default is False.")
	parser.add_argument('--centroid', dest='consensus', default=True, action='store_false',
						help="If set, will compute the centroid sequence as the representative sequence \
						for each passed cluster. \
						Default is to calculate the consensus.")
	parser.add_argument('--debug', dest='debug', action='store_true', default=False,
						help="If set, will run in debug mode.")
	return parser.parse_args()


class Args(object):
	"""docstring for Args"""
	def __init__(self, output=None, temp=None, raw_sequence_dir=None, alignment_pixel_dir=None, log=None,
		db=None, collection=None, collection_prefix=None, collection_suffix=None, ip='localhost', port=27017,
		user=None, password=None, index=None, index_file=None, plate_map=None, index_position='start', index_length=0,
		index_reverse_complement=False, score_cutoff_heavy=200, score_cutoff_light=100, cdhit_threshold=0.96,
		minimum_well_size='relative', minimum_max_well_size=250, minimum_cluster_fraction='largest',
		minimum_well_size_denom=96, cluster_cutoff_gradient=False, consensus=True, debug=False):
		super(Args, self).__init__()
		if not all([output, temp_dir, db]):
			logger.critical("You must supply an output directory, \
				a temp directory and the name of a MongoDB database.")
			sys.exit(1)
		self.output = output
		self.temp_dir = temp
		self.raw_sequence_dir = raw_sequence_dir
		self.alignment_pixel_dir = alignment_pixel_dir
		self.log = log
		self.db = db
		self.collection = collection
		self.collection_prefix = collection_prefix
		self.collection_suffix = collection_suffix
		self.ip = ip
		self.port = int(port)
		self.user = user
		self.password = password
		self.index = index
		self.index_file = index_file
		self.plate_map = plate_map
		self.index_position = index_position
		self.index_reverse_complement = index_reverse_complement
		self.index_length = int(index_length)
		self.score_cutoff_heavy = int(score_cutoff_heavy)
		self.score_cutoff_light = int(score_cutoff_light)
		self.cdhit_threshold = float(cdhit_threshold)
		self.minimum_well_size = minimum_well_size if minimum_well_size == 'relative' else int(minimum_well_size)
		self.minimum_max_well_size = int(minimum_max_well_size)
		self.minimum_cluster_fraction = minimum_cluster_fraction if minimum_cluster_fraction == 'largest' else float(minimum_cluster_fraction)
		self.minimum_well_size_denom = int(minimum_well_size_denom)
		self.cluster_cutoff_gradient = cluster_cutoff_gradient
		self.consensus = consensus
		self.debug = debug


###############
#   MongoDB
###############


def get_sequences(db, collection, chain, score_cutoff):
	seqs = db[collection].find({'chain': chain, 'prod': 'yes', 'v_gene.score': {'$gte': score_cutoff}},
							   {'seq_id': 1, 'raw_input': 1, 'raw_query': 1, 'vdj_nt': 1})
	return [s for s in seqs]


######################
#   SQLite database
######################


def build_seq_db(seqs, temp_dir):
	sys.stdout.flush()
	db_path = os.path.join(temp_dir, 'seq_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	create_cmd = get_seq_db_creation_cmd()
	insert_cmd = get_seq_db_insert_cmd()
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute(create_cmd)
	c.executemany(insert_cmd, seqs)
	sys.stdout.flush()
	start = time.time()
	c.execute('CREATE INDEX seq_index ON seqs (seq_id)')
	return c


def get_seq_db_creation_cmd():
	return '''CREATE TABLE seqs (seq_id text, vdj_nt text)'''


def get_seq_db_insert_cmd():
	return 'INSERT INTO seqs VALUES (?,?)'


def remove_sqlite_db(temp_dir):
	db_path = os.path.join(temp_dir, 'seq_db')
	os.unlink(db_path)


##############
#   CD-HIT
##############


def cdhit_clustering(seqs, bin_id, plate_name, temp_dir, num_plate_seqs,
	minimum_well_size, minimum_well_size_denom, minimum_cluster_fraction,
	raw_sequence_dir, alignment_pixel_dir, consensus, cdhit_threshold):
	logger.info('clustering...')
	seq_db = build_seq_db(seqs, temp_dir)
	infile = make_cdhit_input(seqs, temp_dir)
	outfile = os.path.join(temp_dir, 'clust')
	logfile = open(os.path.join(temp_dir, 'log'), 'a')
	threshold = 0.9
	do_cdhit(infile.name, outfile, logfile, cdhit_threshold)
	clust_handle = open('{}.clstr'.format(outfile), 'r')
	seq, size, total_count, legit = parse_clusters(clust_handle,
												   seq_db,
												   bin_id,
												   plate_name,
												   num_plate_seqs,
												   minimum_well_size,
												   minimum_well_size_denom,
												   minimum_cluster_fraction,
												   temp_dir,
												   raw_sequence_dir,
												   alignment_pixel_dir,
												   consensus)
	os.unlink(infile.name)
	os.unlink(os.path.join(temp_dir, 'log'))
	os.unlink(outfile)
	os.unlink(outfile + '.clstr')
	remove_sqlite_db(temp_dir)
	if legit:
		logger.info('PASSED')
		return seq
	logger.info('FAILED')
	return None


def check_cluster_size(clust_size, num_well_seqs, num_plate_seqs,
	minimum_well_size, minimum_well_size_denom, minimum_cluster_fraction):
	if minimum_well_size == 'relative':
		rel_size = int(num_plate_seqs / float(minimum_well_size_denom))
		logger.info('Minimum well size: {}'.format(rel_size))
		if num_well_seqs >= rel_size:
			return(check_cluster_fraction(clust_size,
										  num_well_seqs,
										  rel_size,
										  minimum_cluster_fraction))
	elif num_well_seqs >= int(minimum_well_size):
		return(check_cluster_fraction(clust_size,
									  num_well_seqs,
									  minimum_well_size,
									  minimum_cluster_fraction))
	logger.info('Minimum cluster fraction: {}'.format(minimum_cluster_fraction))
	return False


def check_cluster_fraction(clust_size, num_well_seqs,
	minimum_well_size, minimum_cluster_fraction):
	if minimum_cluster_fraction == 'largest':
		return True
	if 1. * clust_size / num_well_seqs >= float(minimum_well_size):
		return True
	return False


def make_cdhit_input(seqs, temp_dir):
	fastas = ['>{}\n{}'.format(s[0], s[1]) for s in seqs]
	infile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
	infile.write('\n'.join(fastas))
	infile.close()
	return infile


def do_cdhit(fasta, clust, log, threshold):
	sys.stdout.flush()
	start_time = time.time()
	cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(fasta,
																		  clust,
																		  threshold)
	cluster = sp.Popen(cdhit_cmd, shell=True, stdout=log)
	cluster.communicate()


def parse_centroids(centroid_handle, sizes=None):
	counter = 0
	centroids = []
	for seq in SeqIO.parse(centroid_handle, 'fasta'):
		if sizes:
			size = sizes[counter]
			centroids.append('>{}_{}\n{}'.format(seq.id, size, str(seq.seq)))
		else:
			centroids.append('>{}\n{}'.format(seq.id, str(seq.seq)))
		counter += 1
	return centroids


def parse_clusters(cluster_handle, seq_db, well, plate, num_plate_seqs,
	minimum_well_size, minimum_well_size_denom, minimum_cluster_fraction,
	temp_dir, raw_sequence_dir, alignment_pixel_dir, consensus):
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	cluster_lengths = [(len(c) - 1, c) for c in clusters]
	cluster_lengths.sort(key=lambda x: x[0], reverse=True)
	total_seqs = sum([l[0] for l in cluster_lengths])
	biggest_cluster_size, biggest_cluster = cluster_lengths[0]
	legit = check_cluster_size(biggest_cluster_size,
							   total_seqs,
							   num_plate_seqs,
							   minimum_well_size,
							   minimum_well_size_denom,
							   minimum_cluster_fraction)
	logger.info('Total number of sequences: {}'.format(total_seqs))
	logger.info('Largest cluster: {}'.format(cluster_lengths[0][0]))
	if legit:
		cluster_seq = get_cluster_seq(biggest_cluster,
									  seq_db,
									  plate,
									  well,
									  temp_dir,
									  raw_sequence_dir,
									  consensus)
	else:
		cluster_seq = None
	if all([alignment_pixel_dir is not None, cluster_seq is not None]):
		logger.info('making alignment pixel...')
		ffile = os.path.join(alignment_pixel_dir, '{}-{}.png'.format(plate, well))
		try:
			pixel.make_pixel(get_all_cluster_seqs(biggest_cluster, seq_db),
							ffile,
							temp_dir,
							consentroid=cluster_seq)
		except:
			logger.info('PIXEL EXCEPTION: {}'.format(traceback.format_exc()))
			pass
	return cluster_seq, biggest_cluster_size, total_seqs, legit


def parse_cluster_sizes(cluster_handle):
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	lengths = []
	for cluster in clusters:
		lengths.append(len(cluster) - 1)
	return lengths


def get_cluster_seq(cluster, seq_db, plate, well, temp_dir, raw_sequence_dir, consensus):
	if raw_sequence_dir:
		seqs = get_all_cluster_seqs(cluster, seq_db)
		ofile = os.path.join(raw_sequence_dir, '{}-{}.fasta'.format(plate, well))
		write_raw_binned_data(seqs, ofile)
	if consensus:
		return get_cluster_consensus(cluster, seq_db, temp_dir)
	return get_cluster_centroid(cluster, seq_db)


###############################
#   Centroids and Consensus
###############################


def get_cluster_centroid(cluster, seq_db):
	logger.info('identifying centroid...')
	centroid_id = None
	for c in cluster[1:]:
		if c:
			c = c.strip().split()
			if c[-1] == '*':
				centroid_id = c[2][1:-3]
				break
	if centroid_id:
		seq_db.execute('''SELECT seqs.seq_id, seqs.vdj_nt
						  FROM seqs
						  WHERE seqs.seq_id LIKE "{}"'''.format(centroid_id))
		return seq_db.fetchone()
	# return None


def get_cluster_consensus(cluster, seq_db, temp_dir):
	logger.info('calculating consensus...')
	cluster_seqs = get_all_cluster_seqs(cluster, seq_db)
	return calculate_consensus(cluster_seqs, temp_dir)


def calculate_consensus(cluster_seqs, temp_dir):
	if len(cluster_seqs) == 1:
		return (cluster_seqs[0])
	aln = mafft(cluster_seqs)
	if aln is None:
		return None
	# aln = muscle(cluster_seqs)
	summary_align = AlignInfo.SummaryInfo(aln)
	consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
	consensus_string = str(consensus).replace('-', '')
	return consensus_string.upper()


def get_all_cluster_seqs(cluster, seq_db):
	'''
	Parses a cluster and returns cluster sequences as
	(seq_id, sequence) tuples.
	'''
	ids = []
	seqs = []
	for c in cluster[1:]:
		if c:
			ids.append(c.split()[2][1:-3])
	for chunk in chunker(ids):
		schunk = seq_db.execute('''SELECT seqs.seq_id, seqs.vdj_nt
								FROM seqs
								WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
		seqs.extend(schunk)
	return seqs


def chunker(l, size=900):
	return (l[pos:pos + size] for pos in xrange(0, len(l), size))


###############
#   Indexes
###############


def parse_indexes(index, index_file, index_length):
	indexes = {}
	index_seqs = []
	# parse index sequences from file
	if index is not None:
		index_file = get_builtin_index_file(index)
	with open(index_file) as f:
		for line in f:
			index_seqs.append(line.strip())
	# define row and column ranges
	if len(index_seqs) == 384:
		rows = [c for c in string.ascii_uppercase[:16]]
		columns = [str(c) if len(str(c)) == 2 else '0{}'.format(c) for c in range(1, 25)]
	elif len(index_seqs) == 96:
		rows = [c for c in string.ascii_uppercase[:8]]
		columns = [str(c) if len(str(c)) == 2 else '0{}'.format(c) for c in range(1, 13)]
	else:
		logger.critical('Index file must contain either 96 or 384 index sequences.')
		sys.exit(1)
	# perform a couple of sanity checks on the index sequences
	index_lengths = list(set([len(i) for i in index_seqs]))
	if len(index_lengths) > 1 and not index_length:
		logger.critical('Indexes must all be the same length, or --index-length must be provided.')
		sys.exit(1)
	if not index_length:
		index_length = index_lengths[0]
	if min(index_lengths) < index_length:
		logger.critical('All indexes must be at least as long as --index-length.')
		sys.exit(1)
	# build a dictionary of well locations and index sequences
	wells = []
	for column in columns:
		for row in rows:
			wells.append('{}{}'.format(row, column))
	for well, index in zip(wells, index_seqs):
		indexes[index] = well
	return indexes


def get_builtin_index_file(index):
	mod_dir = os.path.dirname(os.path.abspath(__file__))
	index_dir = os.path.join(mod_dir, 'indexes')
	index_files = {'top-odd': os.path.join(index_dir, 'topodd.txt'),
				   'top-even': os.path.join(index_dir, 'topeven.txt'),
				   'bottom-odd': os.path.join(index_dir, 'bottomodd.txt'),
				   'bottom-even': os.path.join(index_dir, 'bottomeven.txt'),
				   '384': os.path.join(index_dir, '384.txt')}
	return index_files[index]


def parse_plate_map(platemap_file, wells):
	well_names = []
	well_map = {}
	with open(platemap_file) as f:
		for line in f:
			well_names.append(line.strip().split())
	if all([len(w) <= 1 for w in well_names]):
		well_names = [(w, n[0]) for w, n in zip(wells, well_names) if n]
	for wname in well_names:
		if len(wname) == 1:
			continue
		well_map[wname[0]] = wname[1]
	return well_map


def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(seq)
	rc = reversed([complement.get(b, b) for b in bases])
	return ''.join(rc)


def bin_by_index(sequences, indexes, index_length, index_position, rev_comp):
	bins = {w: [] for w in indexes.values()}
	if not index_length:
		index_length = min(list(set([len(i) for i in indexes.keys()])))
	pos = index_length if index_position == 'start' else -1 * index_length
	for seq in sequences:
		if index_position == 'start':
			index = seq['raw_query'][:pos]
		else:
			index = seq['raw_query'][pos:]
		if rev_comp:
			index = reverse_complement(index)
		if index in indexes:
			well = indexes[index]
			bins[well].append((seq['seq_id'], seq['vdj_nt']))
	return bins


##########################
#   Logging and Output
##########################


# def setup_logging(args):
# 	logfile = args.log if args.log else os.path.join(args.output, '{}.log'.format(args.db))
# 	if args.debug:
# 		logger = logging.basicConfig(filename=logfile,
# 							filemode='w',
# 							format='[%(levelname)s] %(asctime)s %(message)s',
# 							level=logging.DEBUG)
# 	else:
# 		logger = logging.basicConfig(filename=logfile,
# 							filemode='w',
# 							format='[%(levelname)s] %(asctime)s %(message)s',
# 							level=logging.INFO)

# 	# set up a streamHandler so that all logged messages also print to console
# 	formatter = logging.Formatter("%(message)s")
# 	ch = logging.StreamHandler()
# 	ch.setLevel(logging.INFO)
# 	ch.setFormatter(formatter)
# 	logger.addHandler(ch)
# 	logger.info('\n')
# 	logger.info('LOG LOCATION: {}'.format(logfile))
# 	return logger


def print_start_info():
	logger.info('')
	logger.info('')
	logger.info('')
	logger.info('-' * 25)
	logger.info('DEMULTIPLEX')
	logger.info('-' * 25)
	logger.info('')


def log_options(args, logfile=None):
	if logfile is not None:
		logger.info('LOG LOCATION: {}'.format(logfile))
	logger.info('OUTPUT DIRECTORY: {}'.format(args.output))
	logger.info('DATABASE: {}'.format(args.db))
	logger.info('INDEX FILE: {}'.format(args.index_file))
	logger.info('PLATE-MAP FILE: {}'.format(args.plate_map))
	logger.info('HEAVY CHAIN SCORE THRESHOLD: {}'.format(args.score_cutoff_heavy))
	logger.info('LIGHT CHAIN SCORE THRESHOLD: {}'.format(args.score_cutoff_light))
	logger.info('USER-SUPPLIED INDEX LENGTH: {}'.format(
		args.index_length if args.index_length > 0 else 'None'))
	logger.info('INDEX POSITION: {}'.format(args.index_position))
	logger.info('INDEX REVERSE COMPLEMENT: {}'.format(args.index_reverse_complement))
	logger.info('CD-HIT THRESHOLD: {}'.format(args.cdhit_threshold))
	logger.info('MINIMUM CLUSTER SIZE: {}'.format(args.minimum_well_size))
	logger.info('MINIMUM CLUSTER SIZE DENOMINATOR: {}'.format(args.minimum_well_size_denom))
	logger.info('MINIMUM CLUSTER FRACTION: {}'.format(args.minimum_cluster_fraction))
	logger.info('CENTROID OR CONSENSUS: {}'.format('consensus' if args.consensus else 'centroid'))


def log_output(bins, seqs, minimum_well_size):
	num_bins = len([b for b in bins.keys() if len(bins[b]) >= minimum_well_size])
	num_seqs = len(seqs)
	logger.info('')
	logger.info('RESULTS: Of {} wells with at least {} reads, {} passed filter'.format(
		num_bins, minimum_well_size, num_seqs))


def make_directories(args):
	output_dir = os.path.dirname(args.output)
	make_dir(output_dir)
	make_dir(args.temp_dir)
	if args.raw_sequence_dir is not None:
		make_dir(args.raw_sequence_dir)
	if args.alignment_pixel_dir is not None:
		make_dir(args.alignment_pixel_dir)


def write_raw_binned_data(seqs, ofile):
	ohandle = open(ofile, 'w')
	ohandle.write('\n'.join(['>{}\n{}'.format(s[0], s[1]) for s in seqs]))


def write_output(seqs, output_file):
	seq_string = '\n'.join(['>{}\n{}'.format(s[0], s[1]) for s in seqs])
	open(output_file, 'w').write(seq_string)


################
#   Printing
################


def print_plate_info(name, collection):
	logger.info('')
	logger.info('')
	logger.info('COLLECTION: {}'.format(collection))
	logger.info('PLATE NAME: {}'.format(name))


def print_bin_info(b):
	bin_string = '  {}  '.format(b)
	logger.info('')
	logger.info(bin_string)
	logger.info('-' * len(bin_string))


def run(**kwargs):
	args = Args(**kwargs)
	main(args)


def main(args, logfile=None):
	global logger
	logger = log.get_logger('demultiplex')
	print_start_info()
	if all([args.index is None, args.index_file is None]):
		err = 'Indexes must be provided, either using --index or --index-file'
		raise RuntimeError(err)
	log_options(args, logfile=logfile)
	make_directories(args)
	db = mongodb.get_db(args.db, ip=args.ip, port=args.port,
					  user=args.user, password=args.password)
	indexes = parse_indexes(args.index, args.index_file, args.index_length)
	plate_map = parse_plate_map(args.plate_map, sorted(indexes.values()))
	all_seqs = []
	collections = mongodb.get_collections(db, args.collection,
		prefix=args.collection_prefix, suffix=args.collection_suffix)
	for collection in collections:
		if collection not in plate_map:
			logger.info('\n\n{} was not found in the supplied plate map file.'.format(
				collection))
			continue
		plate_name = plate_map[collection]
		print_plate_info(plate_name, collection)
		for chain in ['heavy', 'kappa', 'lambda']:
			plate_seqs = []
			logger.info('')
			logger.info('Querying for {} chain sequences'.format(chain))
			score_cutoff = args.score_cutoff_heavy if chain == 'heavy' else args.score_cutoff_light
			sequences = get_sequences(db,
									  collection,
									  chain,
									  score_cutoff)
			logger.info('QUERY RESULTS: {} {} chain sequences met the quality threshold'.format(
				len(sequences), chain.lower()))
			bins = bin_by_index(sequences,
								indexes,
								args.index_length,
								args.index_position,
								args.index_reverse_complement)
			if args.minimum_well_size == 'relative':
				min_well_size = int(len(sequences) / float(args.minimum_well_size_denom))
			else:
				min_well_size = int(args.minimum_well_size)
			min_max_well_size = max(min_well_size, args.minimum_max_well_size)
			if max([len(b) for b in bins.values()]) < int(min_max_well_size):
				logger.info('The biggest well had fewer than {} sequences, so the plate was not processed'.format(min_max_well_size))
				continue
			for b in sorted(bins.keys()):
				if len(bins[b]) < 25:
					continue
				print_bin_info(b)
				consentroid = cdhit_clustering(bins[b],
											   b,
											   plate_name,
											   args.temp_dir,
											   len(sequences),
											   args.minimum_well_size,
											   args.minimum_well_size_denom,
											   args.minimum_cluster_fraction,
											   args.raw_sequence_dir,
											   args.alignment_pixel_dir,
											   args.consensus,
											   args.cdhit_threshold)
				if consentroid:
					consentroid_name = '{}-{}'.format(plate_name, b)
					plate_seqs.append((consentroid_name, consentroid))
			log_output(bins, plate_seqs, min_well_size)
			all_seqs.extend(plate_seqs)
	write_output(all_seqs, args.output)
	logger.info('')


if __name__ == '__main__':
	args = parse_args()
	logfile = args.log if args.log else os.path.join(os.path.dirname(args.output), 'demultiplex.log')
	log.setup_logging(logfile, print_log_location=False, debug=args.debug)
	main(args, logfile=logfile)
