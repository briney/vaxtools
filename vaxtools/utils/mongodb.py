#!/usr/bin/env python
# filename: mongodb.py


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

import logging
import subprocess as sp

from pymongo import MongoClient


def update(field, value, database, collection,
	match_field=None, match=None, ip='localhost', port=27017):
	'''
	::match:: can be a single value or a list of values.
	If ::match:: is a single value, all documents with ::lookup_field:: matching ::match:: will be updated.
	If ::match:: is a list of values, matches to any values in the list will be updated.
	'''
	conn = MongoClient(ip, int(port))
	db = conn[database]
	c = db[collection]
	if lookup_field and match:
		if type(match) != list:
			match = [match, ]
		query_match = {match_field: {'$in': match}}
	else:
		query_match = {}
	c.update_many(query_match, {field: value}, multi=True)
	conn.close()


def mongoimport(jfile, database, collection, ip='localhost', port=27017, user=None, password=None):
	u = " -u {}".format(user) if user else ""
	p = " -p {}".format(password) if password else ""
	# user_password = "{}{} --authenticationDatabase admin".format(username, password)
	mongo_cmd = "mongoimport --host {}:{}{}{} --db {} --collection {} --file {}".format(
		ip, port, u, p, database, collection, jfile)
	mongo = sp.Popen(mongo_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	stdout, stderr = mongo.communicate()
	return stdout
