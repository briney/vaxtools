#!/usr/bin/env python
# filename: cellario.py


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


class Line(object):
    """
    Holds a single line of Run Order data

    Input is a single Run Order line, as a string.
    Optionally, a delimiter can be provided.
    Default delimiter is ','.
    """
    def __init__(self, line, delim=','):
        super(Line, self).__init__()
        self.raw_line = line
        l = line.rstrip('\n').split(delim)
        self.sample_operation_id = l[0]
        self.sample_id = l[1]
        self.protocol_step_id = l[2]
        self.operation_name = l[3]
        self.resource_name = l[4]
        self.destination = l[5]
        self.sample_op_start = l[6]
        self.sample_op_end = l[7]
        self.notes = l[8]
        self.sample_name = l[9]
        self.order_sample_id = l[10]
        self.parent_sample_id = l[11]
        self.current_resource = l[12]
        self.order_sample_id1 = l[13]
        self.run_order_id = l[14]
        self.plate_protocol_id = l[15]
        self.barcode = l[16]
        self.sample_type = l[17]
        self.sample_position_id = l[18]
        self.batch_id = l[19]


class Plate(object):
    """
    Stores barcode and process information for a single plate (either 96-well
    or 384-well). Main purpose is to be able to easily relate the 384-well
    (parent) plate name/info with the 96-well subplates (children).

    Input is a list of lines that correspond to Run Order information
    for a single plate.
    """
    def __init__(self, lines):
        super(Plate, self).__init__()
        self.lines = lines
        self.sample_id = self._get_sample_id()
        self.barcode = self._get_barcode()
        self._parent = None
        self._is_parent = self._check_if_parent()
        self.subplates = []
        self.quadrants = {}
        self.subplate_barcodes = []
        self.plate_type = self._get_plate_type()

    def __repr__(self):
        string = '\n' + self.barcode
        string += '\n{}\n'.format('-' * len(self.barcode))
        subplates = sorted(self.quadrants.items(), key=lambda x: x[1], reverse=True)
        for bc, loc in subplates:
            string += '{}: {}\n'.format(loc, bc)
        return string


    @property
    def is_parent(self):
        return self._is_parent

    @is_parent.setter
    def is_parent(self, val):
        self._is_parent = val

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        if self.is_parent:
            self.is_parent = False
        self._parent = parent


    def add_subplate(self, child, location):
        valid_locs = ['top_odd', 'top_even',
                      'bottom_odd', 'bottom_even']
        if location not in valid_locs:
            estring = 'Location must be one of the following:\n'
            estring += ' '.join(valid_locs)
            raise ValueError(estring)
        self.quadrants[child.barcode] = location
        self.subplates.append(child)
        self.subplate_barcodes.append(child.barcode)

    def get_quadrant(self):
        quadrants = {'P2': 'top_odd', 'P5': 'top_even',
                     'P8': 'bottom_odd', 'P11': 'bottom_even'}
        tfr_line = [l for l in self.lines if l.operation_name.startswith('LiquidTransfer')][0]
        position = tfr_line.resource_name.strip().split()[-1]
        return quadrants[position]


    def _check_if_parent(self):
        parent_sample_ids = list(set([l.parent_sample_id for l in self.lines]))
        if parent_sample_ids[0] == '-1':
            return True
        self._parent = parent_sample_ids[0]
        return False

    def _get_sample_id(self):
        sample_ids = list(set([l.sample_id for l in self.lines]))
        if len(sample_ids) > 1:
            raise RuntimeError('Multiple Sample IDs found for a single plate')
        return sample_ids[0]

    def _get_barcode(self):
        barcodes = list(set([l.sample_name for l in self.lines]))
        if len(barcodes) > 1:
            raise RuntimeError('Multiple barcodes found for a single plate')
        return barcodes[0]

    def _get_plate_type(self):
        types = {'PCR Plate': 'pcr',
                 'ELISA Plate': 'pcr',
                 'P20  Tips': 'tips'}
        sample_type = self.lines[0].sample_type
        return types[sample_type]


def parse_run_order(run_order_file, plates=None, return_parents_only=True):
    """
    Parses a Cellario Run Order to determine which 96-well template plates were
    pooled into a single 384-well PCR plate.

    Inputs are a Run Order file and, optionally, a list of Plate objects to which the
    current parsed plates should be appended. The plate list option provides a simple
    means to combine the parsing results from multiple Run Order files:

    plates = []
    for f in run_order_file_list:
        plates = parse_pooled_plates(f, plates=plates)

    If ::return_parents_only:: is True (whcih is the default), only the 384-well
    (parent) plates will be returned. Set to False to return all plates (including the
    96-well 'children' plates)

    Returns a list of Plate objects.
    """
    if plates is None:
        plates = []
    _plates = []
    lines = []
    header = False
    prev_line = ''
    with open(run_order_file) as f:
        for line in f:
            # we don't need the file header
            if not header:
                while line[:10] != 'SAMPLEOPER':
                    line = f.next()
                if line[:10] == 'SAMPLEOPER':
                    header = True
                    line = f.next()
            # WAIT lines are truncated, we can skip
            if 'WAIT' in line:
                continue
            l = Line(line)
            if l.sample_name != prev_line:
                prev_line = l.sample_name
                if lines:
                    _plates.append(Plate(lines))
                    lines = []
            lines.append(l)
        # assign the last batch of lines to a plate
        if lines:
            _plates.append(Plate(lines))

    # we don't care about the tip boxes, which are handled by the Run Order
    # file as if they were plates
    pcr_plates = [p for p in _plates if p.plate_type == 'pcr']
    parents = {p.sample_id: p for p in pcr_plates if p.is_parent}
    children = [p for p in pcr_plates if not p.is_parent]
    for c in children:
        q = c.get_quadrant()
        parents[c.parent].add_subplate(c, q)
    if return_parents_only:
        par_plates = plates + parents.values()
        return sorted(par_plates, key=lambda x: x.barcode)
    all_plates = plates + parents.values() + children
    return sorted(all_plates, key=lambda x: x.barcode)
