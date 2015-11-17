#!/usr/bin/env python
# filename: scpcr_demultiplexing.py


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

import argparse
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
import urllib
import uuid

from pymongo import MongoClient

import numpy as np
import pandas as pd

from Bio import AlignIO
from Bio.Align import AlignInfo

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', dest='output', required=True,
					help="Output directory for demultiplexed FASTA files. Required.")
parser.add_argument('-t', '--temp', dest='temp_dir', required=True,
					help="Directory for temporary files. Required.")
parser.add_argument('-r', '--raw-sequences', dest='raw_sequence_dir', default=None,
					help="Directory for saving the raw sequence data (not just consensus/centroid sequences. \
						If not provided, raw binned sequence data will not be retained.")
parser.add_argument('-a', '--alignment-pixel-dir', dest='alignment_pixel_dir', default=None,
					help="Directory for saving the alignment pixel figures. \
						If not provided, figures will not be generated.")
parser.add_argument('-l', '--log', dest='log',
					help="Location for the log file. \
					Default is <output>/<database>.log.")
parser.add_argument('-d', '--database', dest='db', required=True,
					help="Name of the MongoDB database containing un-demultiplexed sequences. Required.")
parser.add_argument('-c', '--collection', dest='collection', default=None,
					help="Name of the MongoDB collection to query. \
					If not provided, all collections in the given database will be processed iteratively.")
parser.add_argument('-i', '--ip', dest='ip', default='localhost',
					help="IP address for the MongoDB server.  Defaults to 'localhost'.")
parser.add_argument('-u', '--user', dest='user', default=None,
					help="Username for the MongoDB server. Not used if not provided.")
parser.add_argument('-p', '--password', dest='password', default=None,
					help="Password for the MongoDB server. Not used if not provided.")
parser.add_argument('-I', '--index-file', dest='index_file', required=True,
					help="File containing the indexes, one per line. \
					File must have either 96 or 384 indexes, and indexes should be in columnar order \
					(A01, B01, C01, etc...).")
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
					which would produce some strange (and likely wrong) results \
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
parser.add_argument('--consensus', default=False, action='store_true',
					help="Calculate the consensus sequence as the representative sequence \
					for each passed cluster. \
					Default is to calculate the centroid, not the consensus.")
# parser.add_argument('--all-raw-sequences', default=False, action='store_true',
# 					help="If set, all sequences from each well will be saved. \
# 					Default is to save raw sequences only for the cluster used to calculate \
# 					centroid or consensus sequences.")
parser.add_argument('--debug', dest='debug', action='store_true', default=False,
					help="If set, will run in debug mode.")
args = parser.parse_args()



###############
#   MongoDB
###############

if args.user and args.password:
	password = urllib.quote_plus(password)
	uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
	conn = MongoClient(uri)
else:
	conn = MongoClient(args.ip, 27017)
db = conn[args.db]


def get_collections():
	if args.collection:
		return [args.collection, ]
	collections = db.collection_names(include_system_collections=False)
	return sorted(collections)


def get_sequences(collection, chain):
	score_cutoff = args.score_cutoff_heavy if chain == 'heavy' else args.score_cutoff_light
	seqs = db[collection].find({'chain': chain, 'prod': 'yes', 'v_gene.score': {'$gte': score_cutoff}},
							   {'seq_id': 1, 'raw_input': 1, 'raw_query': 1, 'vdj_nt': 1})
	return [s for s in seqs]


######################
#   SQLite database
######################


def build_seq_db(seqs):
	sys.stdout.flush()
	db_path = os.path.join(args.temp_dir, 'seq_db')
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


def remove_sqlite_db():
	db_path = os.path.join(args.temp_dir, 'seq_db')
	os.unlink(db_path)


##############
#   CD-HIT
##############


def cdhit_clustering(seqs, bin_id, plate_name, num_plate_seqs):
	print('clustering...')
	# legit = False
	seq_db = build_seq_db(seqs)
	temp_dir = args.temp_dir
	infile = make_cdhit_input(seqs)
	outfile = os.path.join(temp_dir, 'clust')
	logfile = open(os.path.join(temp_dir, 'log'), 'a')
	threshold = 0.9
	do_cdhit(infile.name, outfile, logfile)
	clust_handle = open('{}.clstr'.format(outfile), 'r')
	seq, size, total_count, legit = parse_clusters(
		clust_handle, seq_db, bin_id, plate_name, num_plate_seqs)
	os.unlink(infile.name)
	os.unlink(os.path.join(args.temp_dir, 'log'))
	os.unlink(outfile)
	os.unlink(outfile + '.clstr')
	remove_sqlite_db()
	# legit = check_cluster_size(size, len(seqs), num_plate_seqs)
	# if args.cluster_cutoff_gradient:
	# 	legit = cutoff_gradient(size, len(seqs))
	# if size >= args.minimum_well_size and 1. * size / len(seqs) >= args.minimum_cluster_fraction:
	# 	legit = True
	if legit:
		print('PASS')
		logging.info('{}: PASSED, {} total sequences, {} in largest cluster'.format(
			bin_id, total_count, size))
		return seq
	# if size < args.minimum_well_size:
	# 	print('FAIL (cluster size is too small)')
	# elif 1. * size / len(seqs) < args.minimum_cluster_fraction:
	# 	print('FAIL')
	print('FAIL')
	logging.info('{}: FAILED, {} total sequences, {} in largest cluster'.format(
			bin_id, total_count, size))
	return None


def check_cluster_size(clust_size, num_well_seqs, num_plate_seqs):
	if args.cluster_cutoff_gradient:
		return (cluster_gradient(clust_size, num_well_seqs))
	if args.minimum_well_size == 'relative':
		rel_size = int(num_plate_seqs / float(args.minimum_well_size_denom))
		print('Minimum well size: {}'.format(rel_size))
		logging.info('Minimum well size: {}'.format(rel_size))
		if num_well_seqs >= rel_size:
			return(check_cluster_fraction(clust_size, num_well_seqs))
	elif num_well_seqs >= int(args.minimum_well_size):
		return(check_cluster_fraction(clust_size, num_well_seqs))
	print('Minimum cluster fraction: {}'.format(args.minimum_cluster_fraction))
	logging.info('Minimum cluster fraction: {}'.format(args.minimum_cluster_fraction))
	return False


def check_cluster_fraction(clust_size, num_well_seqs):
	print('Minimum cluster fraction: {}'.format(args.minimum_cluster_fraction))
	logging.info('Minimum cluster fraction: {}'.format(args.minimum_cluster_fraction))
	if args.minimum_cluster_fraction == 'largest':
		return True
	if 1. * clust_size / num_well_seqs >= float(args.minimum_well_size):
		return True
	return False


def cluster_gradient(size, total):
	if size >= 25 and 1. * size / len(seqs) >= 0.8:
		return True
	if size >= 100 and 1. * size / len(seqs) >= 0.6:
		return True
	if size >= 250 and 1. * size / len(seqs) >= 0.5:
		return True
	return False


def make_cdhit_input(seqs):
	fastas = ['>{}\n{}'.format(s[0], s[1]) for s in seqs]
	infile = tempfile.NamedTemporaryFile(dir=args.temp_dir, delete=False)
	infile.write('\n'.join(fastas))
	infile.close()
	return infile


def do_cdhit(fasta, clust, log):
	sys.stdout.flush()
	start_time = time.time()
	cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(fasta, clust, args.cdhit_threshold)
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


def parse_clusters(cluster_handle, seq_db, well, plate, num_plate_seqs):
	sys.stdout.flush()
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	sys.stdout.flush()
	start = time.time()
	# cluster_lengths = []
	# for cluster in clusters:
	# 	length = len(cluster) - 1
	# 	# cluster_id = get_cluster_id(cluster)
	# 	cluster_seq_id, cluster_seq = get_cluster_seq(cluster, seq_db)
	# 	cluster_lengths.append((length, (cluster_seq_id, cluster_seq)))
	cluster_lengths = [(len(c) - 1, c) for c in clusters]
	cluster_lengths.sort(key=lambda x: x[0], reverse=True)
	total_seqs = sum([l[0] for l in cluster_lengths])
	biggest_cluster_size, biggest_cluster = cluster_lengths[0]
	legit = check_cluster_size(biggest_cluster_size, total_seqs, num_plate_seqs)
	sum_string = 'Total number of sequences: {}\nLargest cluster: {}'.format(total_seqs, cluster_lengths[0][0])
	print(sum_string)
	logging.info(sum_string)
	if legit:
		cluster_seq_id, cluster_seq = get_cluster_seq(biggest_cluster, seq_db, plate, well)
	else:
		cluster_seq = None
	if args.alignment_pixel_dir:
		print('making alignment pixel...')
		ffile = os.path.join(args.alignment_pixel_dir, '{}_{}.png'.format(plate, well))
		try:
			make_alignment_pixel(cluster_seq, get_all_cluster_seqs(biggest_cluster, seq_db), ffile)
		except:
			pass
	return cluster_seq, biggest_cluster_size, total_seqs, legit


def parse_cluster_sizes(cluster_handle):
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	lengths = []
	for cluster in clusters:
		lengths.append(len(cluster) - 1)
	return lengths


def get_cluster_seq(cluster, seq_db, plate, well):
	if args.raw_sequence_dir:
		seqs = get_all_cluster_seqs(cluster, seq_db)
		ofile = os.path.join(args.raw_sequence_dir, '{}_{}.fasta'.format(plate, well))
		write_raw_binned_data(seqs, ofile)
	if args.consensus:
		return get_cluster_consensus(cluster, seq_db)
	return get_cluster_centroid(cluster, seq_db)


###############################
#   Centroids and Consensus
###############################


def get_cluster_centroid(cluster, seq_db):
	print('identifying centroid...')
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
	return None


def get_cluster_consensus(cluster, seq_db):
	print('calculating consensus...')
	cluster_seqs = get_all_cluster_seqs(cluster, seq_db)
	return (calculate_consensus(cluster_seqs))


def calculate_consensus(cluster_seqs):
	if len(cluster_seqs) == 1:
		return (cluster_seqs[0])
	fasta_string = '.\n'.join(['>{}\n{}'.format(c[0], c[1]) for c in cluster_seqs])
	if len(cluster_seqs) < 100:
		muscle_cline = 'muscle -clwstrict'
	elif len(cluster_seqs) < 1000:
		muscle_cline = 'muscle -clwstrict -maxiters 2'
	else:
		muscle_cline = 'muscle -clwstrict -maxiters 1 -diags'
	muscle = sp.Popen(str(muscle_cline),
					  stdin=sp.PIPE,
					  stdout=sp.PIPE,
					  stderr=sp.PIPE,
					  universal_newlines=True,
					  shell=True)
	alignment = muscle.communicate(input=fasta_string)[0]
	alignment_string = AlignIO.read(StringIO(alignment), 'clustal')
	summary_align = AlignInfo.SummaryInfo(alignment_string)
	consensus = summary_align.gap_consensus(threshold=0.51, ambiguous='n')
	consensus_id = uuid.uuid4()
	consensus_string = str(consensus).replace('-', '')
	return (consensus_id, consensus_string.upper())


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


def parse_indexes():
	indexes = {}
	index_seqs = []
	# parse index sequences from file
	with open(args.index_file) as f:
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
		print('Index file must contain either 96 or 384 index sequences.')
		sys.exit(1)
	# perform a couple of sanity checks on the index sequences
	index_lengths = list(set([len(i) for i in index_seqs]))
	if len(index_lengths) > 1 and not args.index_length:
		print('Indexes must all be the same length, or --index-length must be provided.')
		sys.exit(1)
	if not args.index_length:
		args.index_length = index_lengths[0]
	if min(index_lengths) < args.index_length:
		print('All indexes must be at least as long as --index-length.')
		sys.exit(1)
	# build a dictionary of well locations and index sequences
	wells = []
	for column in columns:
		for row in rows:
			wells.append('{}{}'.format(row, column))
	for well, index in zip(wells, index_seqs):
		indexes[index] = well
	return indexes


def parse_plate_map(wells):
	well_names = []
	well_map = {}
	with open(args.plate_map) as f:
		for line in f:
			well_names.append(line.strip().split())
	if len(well_names[0]) == 1:
		well_names = [(w, n[0]) for w, n in zip(wells, well_names) if n]
	for (w, n) in well_names:
		well_map[w] = n
	return well_map


def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(seq)
	rc = reversed([complement.get(b, b) for b in bases])
	return ''.join(rc)


def bin_by_index(sequences, indexes):
	bins = {w: [] for w in indexes.values()}
	pos = args.index_length if args.index_position == 'start' else -1 * args.index_length
	for seq in sequences:
		if args.index_position == 'start':
			index = seq['raw_query'][:pos]
		else:
			index = seq['raw_query'][pos:]
		if args.index_reverse_complement:
			index = reverse_complement(index)
		if index in indexes:
			well = indexes[index]
			bins[well].append((seq['seq_id'], seq['vdj_nt']))
	return bins


##########################
#   Logging and Output
##########################


def setup_logging():
	logfile = args.log if args.log else os.path.join(args.output, '{}.log'.format(args.db))
	if args.debug:
		logging.basicConfig(filename=logfile,
							filemode='w',
							format='[%(levelname)s] %(asctime)s %(message)s',
							level=logging.DEBUG)
	else:
		logging.basicConfig(filename=logfile,
							filemode='w',
							format='[%(levelname)s] %(asctime)s %(message)s',
							level=logging.INFO)
	logging.info('LOG LOCATION: {}'.format(logfile))
	log_options()


def log_options():
	logging.info('OUTPUT DIRECTORY: {}'.format(args.output))
	logging.info('DATABASE: {}'.format(args.db))
	logging.info('INDEX FILE: {}'.format(args.index_file))
	logging.info('PLATE-MAP FILE: {}'.format(args.plate_map))
	logging.info('HEAVY CHAIN SCORE THRESHOLD: {}'.format(args.score_cutoff_heavy))
	logging.info('LIGHT CHAIN SCORE THRESHOLD: {}'.format(args.score_cutoff_light))
	logging.info('USER-SUPPLIED INDEX LENGTH: {}'.format(
		args.index_length if args.index_length > 0 else 'None'))
	logging.info('INDEX POSITION: {}'.format(args.index_position))
	logging.info('INDEX REVERSE COMPLEMENT: {}'.format(args.index_reverse_complement))
	logging.info('CD-HIT THRESHOLD: {}'.format(args.cdhit_threshold))
	logging.info('MINIMUM CLUSTER SIZE: {}'.format(args.minimum_well_size))
	logging.info('MINIMUM CLUSTER SIZE DENOMINATOR: {}'.format(args.minimum_well_size_denom))
	logging.info('MINIMUM CLUSTER FRACTION: {}'.format(args.minimum_cluster_fraction))
	logging.info('CENTROID OR CONSENSUS: {}'.format('consensus' if args.consensus else 'centroid'))


def log_output(bins, seqs):
	num_bins = len([b for b in bins.keys() if len(bins[b]) >= args.minimum_well_size])
	num_seqs = len(seqs)
	logging.info('RESULTS: Of {} wells with at least {} reads, {} passed filter'.format(
		num_bins, args.minimum_well_size, num_seqs))


def write_raw_binned_data(seqs, ofile):
	ohandle = open(ofile, 'w')
	ohandle.write('\n'.join(['>{}\n{}'.format(s[0], s[1]) for s in seqs]))


def write_output(seqs, outname):
	outpath = os.path.join(args.output, outname)
	seq_string = '\n'.join(['>{}\n{}'.format(s[0], s[1]) for s in seqs])
	open(outpath, 'w').write(seq_string)


#######################
#   Alignment Pixel
#######################


def make_alignment_pixel(consentroid, seqs, ffile):
	nuc_vals = {'-': 0,
				'A': 1,
				'C': 2,
				'G': 3,
				'T': 4}
	colors = {'-': 'w',
			  'A': '#E12427',
			  'C': '#3B7FB6',
			  'G': '#63BE7A',
			  'T': '#E1E383'}
	cmap = ListedColormap(['w', '#E12427', '#3B7FB6', '#63BE7A', '#E1E383', 'k'])
	data = []
	aligned_seqs = pixel_msa(consentroid, seqs)
	for seq in aligned_seqs:
		seq_data = []
		s_id = seq[0]
		sequence = seq[1]
		if s_id == 'consentroid':
			sequence = 'XXXXXX--' + sequence
		else:
			sequence = '--------' + sequence
		for res in sequence:
			seq_data.append(nuc_vals.get(res.upper(), 5))
		data.append(seq_data)
	mag = (magnitude(len(data[0])) + magnitude(len(data))) / 2
	x_dim = len(data[0]) / 10**mag
	y_dim = len(data) / 10**mag
	plt.figure(figsize=(x_dim, y_dim), dpi=100)
	plt.imshow(data, cmap=cmap, interpolation='none')
	plt.axis('off')
	plt.savefig(ffile, bbox_inches='tight', dpi=400)
	plt.close()


def pixel_msa(consentroid, seqs):
	fasta = '>consentroid\n{}\n'.format(consentroid)
	fasta += '\n'.join(['>{}\n{}'.format(s[0], s[1]) for s in seqs])
	if len(seqs) < 100:
		muscle_cline = 'muscle -clwstrict'
	elif len(seqs) < 1000:
		muscle_cline = 'muscle -clwstrict -maxiters 2'
	else:
		muscle_cline = 'muscle -clwstrict -maxiters 1 -diags'
	muscle = sp.Popen(str(muscle_cline),
					  stdin=sp.PIPE,
					  stdout=sp.PIPE,
					  stderr=sp.PIPE,
					  universal_newlines=True,
					  shell=True)
	alignment = muscle.communicate(input=fasta)[0]
	aln = AlignIO.read(StringIO(alignment), 'clustal')
	aln_seqs = []
	for record in aln:
		aln_seqs.append((record.id, str(record.seq)))
	return aln_seqs


def magnitude(x):
	return int(math.log10(x))


################
#   Printing
################


def print_plate_info(name, collection):
	name_string = '     Processing {} (collection {})     '.format(name, collection)
	print('\n\n')
	print('=' * len(name_string))
	print(name_string)
	print('=' * len(name_string))
	logging.info('')
	logging.info('COLLECTION: {}'.format(collection))
	logging.info('PLATE NAME: {}'.format(name))


def print_bin_info(b):
	bin_string = '  {}  '.format(b)
	print('\n')
	print(bin_string)
	print('-' * len(bin_string))




def main():
	indexes = parse_indexes()
	plate_map = parse_plate_map(sorted(indexes.values()))
	for collection in get_collections():
		if collection not in plate_map:
			print('\n\n{} was not found in the supplied plate map file.'.format(
				collection))
			continue
		plate_name = plate_map[collection]
		print_plate_info(plate_name, collection)
		plate_seqs = []
		for chain in ['heavy', 'kappa', 'lambda']:
			print('\n\nQuerying for {} chain sequences'.format(chain))
			logging.info('{} CHAIN'.format(chain.upper()))
			sequences = get_sequences(collection, chain)
			print('Retrieved {} sequences.\n'.format(len(sequences)))
			logging.info('QUERY RESULTS: {} {} chain sequences met the quality threshold'.format(
				len(sequences), chain.lower()))
			bins = bin_by_index(sequences, indexes)
			min_well_size = args.minimum_max_well_size if args.minimum_well_size == 'relative' else args.minimum_well_size
			if max([len(b) for b in bins.values()]) < int(min_well_size):
				logging.info('The biggest well had fewer than {} sequences, so the plate was not processed'.format(min_well_size))
				continue
			for b in sorted(bins.keys()):
				if len(bins[b]) < 25:
					continue
				print_bin_info(b)
				centroid = cdhit_clustering(bins[b], b, plate_name, len(sequences))
				if centroid:
					plate_seqs.append((b, centroid))
			log_output(bins, plate_seqs)
		write_output(plate_seqs, plate_name)
	print('\n')


if __name__ == '__main__':
	setup_logging()
	main()
