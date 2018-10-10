#!/usr/bin/env python
# filename: input_data.py


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

from abc import ABC, abstractmethod
import json
import os

from abutils.core.sequence import Sequence
from abutils.utils import mongodb
from abutils.utils.decorators import lazy_property



class InputData(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def get_sequences(self):
        pass
    
    def parse_platemap(self, platemap_file):
        platemap = {}
        with open(platemap_file) as f:
            for line in f:
                if len(line.strip().split()) == 2:
                    k, v = line.strip().split()
                    platemap[k] = v
        return platemap


class MongoData(InputData):
    '''
    Docstring for MongoData
    '''
    def __init__(self, db, collection, args):
        super().__init__()
        self.db = db
        self.collection = collection
        self.args = args
        self.raw_name = collection
        self.data_type = 'MongoDB collection'
        self.sequence_retrieval_string = 'Querying for {} chain sequences'

    @lazy_property
    def name(self):
        if self.args.plate_map is None:
            delim = self.args.plate_name_delimiter
            pos = self.args.plate_name_delimiter_position
            if delim is None:
                return self.collection
            else:
                return delim.join(self.collection.split(delim)[:pos])
        else:
            platemap = self.parse_platemap(self.args.plate_map)
            return platemap.get(self.collection, None)
    
    def get_sequences(self, chain, score_cutoff):
        query = {'chain': chain, 'prod': 'yes', 'v_gene.score': {'$gte': score_cutoff}}
        projection = {'seq_id': 1, 'raw_input': 1, 'raw_query': 1, 'oriented_input': 1, 'vdj_nt': 1}
        seqs = self.db[self.collection].find(query, projection)
        return [Sequence(s) for s in seqs]

    
class JsonData(InputData):
    '''
    Docstring for JsonData.
    '''
    def __init__(self, json_file, args):
        super().__init__()
        self.json_file = json_file
        self.args = args
        self.raw_name = json_file
        self.data_type = 'JSON file'
        self.sequence_retrieval_string = 'Reading {} chain sequences from JSON file'

    @lazy_property
    def name(self):
        name = os.path.basename(self.json_file)
        if self.args.plate_map is None:
            delim = self.args.plate_name_delimiter
            pos = self.args.plate_name_delimiter_position
            if delim is None:
                return name.rstrip('.json')
            else:
                return delim.join(name.split(delim)[:pos])
        else:
            platemap = self.parse_platemap(self.args.plate_map)
            return platemap.get(name, None)

    def get_sequences(self, chain, score_cutoff):
        seqs = []
        keys = ['seq_id', 'raw_input', 'raw_query', 'oriented_input', 'vdj_nt']
        with open(self.json_file) as f:
            for line in f:
                if line.strip():
                    j = json.loads(line.strip())
                    if all([j['chain'] == chain, j['prod'] == 'yes', j['v_gene']['score'] >= score_cutoff]):
                        seq = {key: j[key] for key in keys if j.get(key, None) is not None}
                        seqs.append(Sequence(seq))
        return seqs
        

