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

from abtools import mongodb
from abtools.pipeline import list_files


def assign_sample_metadata(database, collection, ip='localhost', port=27017,
                           user=None, password=None, subjects=None,
                           groups=None, experiments=None, timepoints=None):
    db = mongodb.get_db(database, ip=ip, port=port)
    if subjects is not None:
        assign_subjects(db, collection, subjects)
    if groups is not None:
        assign_groups(db, collection, groups)
    if experiments is not None:
        assign_experiments(db, collection, experiments)
    if timepoints is not None:
        assign_timepoints(db, collection, timeponts)


def unset_sample_metadata(db, collection, subject=False, group=False,
                          experiment=False, timepoint=False):
    if subject:
        mongodb.unset(db, collection, field='subject')
    if group:
        mongodb.unset(db, collection, field='group')
    if experiment:
        mongodb.unset(db, collection, field='experiment')
    if timepoint:
        mongodb.unset(db, collection, field='timepoint')



def assign_subjects(db, collection, platemap_dir):
    '''
    Assigns subject names (updates the provided MongoDB db/collection with
    a 'subject' field).

    Inputs:
        ::db:: is a pymongo database object, containing the sequences
            to be updated
        ::collection:: is a MongoDB collection name, as a string
        ::platemap_dir:: is the directory containing one or more platemap files
            (output from vaxtools.scpcr.platemap)
    '''
    mongodb.remove_padding(db, collection)
    subjects = parse_samplemaps(platemap_dir)
    update_subjects(db, collection, subjects)


def parse_samplemaps(platemap_dir):
    subjects = {}
    for platemap in list_files(platemap_dir):
        plate = os.path.basename(platemap)
        with open(platemap) as f:
            for line in f:
                sline = line.strip().split()
                if len(sline) >= 2:
                    seq_id = sline[0]
                    subject = sline[1]
                    # platewell = '{}-{}'.format(plate, well)
                    if subject not in subjects:
                        subjects[subject] = []
                    subjects[subject].append(seq_id)
    return subjects


def update_subjects(db, collection, subjects):
    for subject in subjects:
        seq_ids = subjects[subject]
        match = {'seq_id': {'$in': seq_ids}}
        mongodb.update(field='subject', value=subject,
            db=db, collection=collection, match=match)



def assign_groups(db, collection, groupmap, by_sequence=False, subject_field='subject'):
    '''
    Assigns group names (updates the provided MongoDB db/collection with
    a 'group' field).

    Inputs:
        ::db:: is a pymongo database object, containing the sequences
            to be updated
        ::collection:: is a MongoDB collection name, as a string
        ::groupmap:: is a file containing group assignments, of the format:
                subject_name  group_name
            separated by any whitespace (one subject/group entry per line).

    Only the first occurance of each subject_name will be used.

    If the groupmap file contains sequence names instead of subject names,
    set ::by_sequence:: to True.
    '''
    mongodb.remove_padding(db, collection)
    groups = parse_groupmap(groupmap)
    update_groups(db, collection, groups, by_sequence, subject_field)


def parse_groupmap(groupmap):
    groups = {}
    with open(groupmap) as f:
        for line in f:
            sline = line.strip().split()
            if sline:
                sample = sline[0]
                if sample in groups:
                    continue
                group = sline[1]
                groups[sample] = group
    return groups


def update_groups(db, collection, groups, by_sequence=False, subject_field='subject'):
    mongodb.remove_padding(db, collection)
    subjects = sorted(groups.keys())
    for subject in subjects:
        group = groups.get(subject, None)
        if group is None:
            continue
        if by_sequence:
            seq_ids = [subject, ]
        else:
            seq_ids = get_subject_seq_ids(db, collection, subject, subject_field)
        match = {'seq_id': {'$in': seq_ids}}
        mongodb.update(field='group', value=group,
            db=db, collection=collection, match=match)


def get_subject_seq_ids(db, collection, subject, subject_field='subject'):
    c = db[collection]
    seqs = c.find({subject_field: subject}, {'seq_id': 1, '_id': 0})
    return [s['seq_id'] for s in seqs]



def assign_experiments(db, collection, experimentmap, subject_field='subject'):
    '''
    Assigns experiment names (updates the provided MongoDB db/collection with
    a 'experiment' field).

    Inputs:
        ::db:: is a pymongo database object, containing the sequences
            to be updated
        ::collection:: is a MongoDB collection name, as a string
        ::experimentmap:: is a file containing experiment assignments, of the format:
                subject_name  experiment_name
            separated by any whitespace (one sample/experiment entry per line).
    '''
    mongodb.remove_padding(db, collection)
    experiments = parse_experimentmap(experimentmap)
    update_experiments(db, collection, experiments, subject_field)


def parse_experimentmap(experimentmap):
    experiments = {}
    with open(experimentmap) as f:
        for line in f:
            sline = line.strip().split()
            if sline:
                subject = sline[0]
                experiment = sline[1]
                experiments[subject] = experiment
    return experiments


def update_experiments(db, collection, experiments, subject_field):
    subjects = sorted(experiments.keys())
    for subject in subjects:
        experiment = experiments.get(subject, None)
        if experiment is None:
            continue
        seq_ids = get_subject_seq_ids(db, collection, subject, subject_field)
        match = {'seq_id': {'$in': seq_ids}}
        mongodb.update(field='experiment', value=experiment,
            db=db, collection=collection, match=match)



def assign_timepoints(db, collection, timepointmap, by_subject=False, subject_field='subject'):
    '''
    Assigns timepoint names (updates the provided MongoDB db/collection with
    a 'timepoint' field).

    Inputs:
        ::db:: is a pymongo database object, containing the sequences
            to be updated
        ::collection:: is a MongoDB collection name, as a string
        ::timepointmap:: is a file containing timepoint assignments, of the format:
                seq_id  timepoint_name
            separated by any whitespace (one seq_id/timepoint entry per line).
    '''
    mongodb.remove_padding(db, collection)
    timepoints = parse_timepointmap(timepointmap)
    update_timepoints(db, collection, timepoints, by_subject, subject_field)


def parse_timepointmap(timepointmap):
    timepoints = {}
    with open(timepointmap) as f:
        for line in f:
            sline = line.strip().split()
            if sline:
                seq_id = sline[0]
                timepoint = sline[1]
                if timepoint not in timepoints:
                    timepoints[timepoint] = []
                timepoints[timepoint].append(seq_id)
    return timepoints


def update_timepoints(db, collection, timepoints, by_subject=False, subject_field='subject'):
    timepoint_names = sorted(timepoints.keys())
    for timepoint in timepoint_names:
        if by_subject:
            seq_ids = []
            for subject in timepoints.get(timepoint, None):
                if subject is None:
                    continue
                seq_ids.extend(get_subject_seq_ids(db, collection, subject, subject_field))
        else:
            seq_ids = timepoints.get(timepoint, None)
        match = {'seq_id': {'$in': seq_ids}}
        mongodb.update(field='timepoint', value=timepoint,
            db=db, collection=collection, match=match)
