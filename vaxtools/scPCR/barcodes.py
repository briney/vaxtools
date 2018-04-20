#!/usr/bin/env python
# filename: barcodes.py

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
import itertools
import os
import string

from openpyxl import load_workbook

from abutils.utils import log
from abutils.utils.pipeline import list_files

from vaxtools.utils.containers import Well, Sample, Plate



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True,
                        help="The input file, in Excel format. Required.")
    parser.add_argument('-o', '--out', dest='output', required=True,
                        help="The output file, into which the platemap data will be written. \
                        Required.")
    parser.add_argument('-l', '--log', dest='log', default=None,
                        help="The log file. If not provided, defaults to <output>/platemap.log.")
    parser.add_argument('-p', '--plate-name-prefix', dest='prefix', default=None,
                        help="Prefix for the plate name.")
    parser.add_argument('-s', '--plate-name-suffix', dest='suffix', default=None,
                        help="Suffix for the plate name.")
    parser.add_argument('--plate-delim-prefix', required=True,
                        help='Text at the start of a well marking the boundary between two plates. \
                        Often, the best choice is the well with the plate name (if part of the platemap) \
                        or the sample row. Required.')
    parser.add_argument('--plate-delim-well-number', required=True, type=int,
                        help="The well number (using 1-based indexing) of the plate delimiter well. \
                        Used in conjunction with --plate-delim-row-number. \
                        Required.")
    parser.add_argument('--plate-name-row-number', default=None, type=int,
                        help="The number of rows between the plate delimiter row and the row containing \
                        the plate name. Only required if parsing plate names. \
                        Default is None, which results in plate names not being parsed.")
    parser.add_argument('--plate-name-well-number', default=None, type=int,
                        help="The well number (using 1-based indexing) of the plate name well. \
                        Used in conjunction with --plate-name-row-number. \
                        Default is None, which results in plate names not being parsed.")
    parser.add_argument('--barcode-row-number', default=None, type=int,
                        help="The number of rows between the plate delimiter row and the row containing \
                        the plate name. Only required if parsing plate names. \
                        Default is None, which results in plate names not being parsed.")
    parser.add_argument('--barcode-well-number', default=None, type=int,
                        help="The well number (using 1-based indexing) of the plate name well. \
                        Used in conjunction with --plate-name-row-number. \
                        Default is None, which results in plate names not being parsed.")
    parser.add_argument('-D', '--debug')
    return parser.parse_args()


class Args(object):
    """docstring for Args"""
    def __init__(self, input=None, output=None, log=None,
        prefix=None, suffix=None,
        plate_delim_prefix=None, plate_delim_well_number=None,
        plate_name_row_number=None, plate_name_well_number=None,
        barcode_row_number=None, barcode_well_number=None,
        debug=False):
        super(Args, self).__init__()

        if not all([input, output, plate_delim_prefix, plate_delim_well_number, barcode_row_number, barcode_well_number]):
            err_string = 'The following options are all required:\n'
            err_string += '--input\n--output\n'
            err_string += '--plate-delim-prefix\n--plate-delim-well-number\n'
            err_string += '--sample-row-number\n--sample-well-number'
            raise RuntimeError(err_string)

        self.input = input
        self.output = output
        self.log = log
        self.prefix = prefix
        self.suffix = suffix
        self.plate_delim_prefix = plate_delim_prefix
        self.plate_delim_well_number = int(plate_delim_well_number)
        self.plate_name_row_number = int(plate_name_row_number) if plate_name_row_number is not None else None
        self.plate_name_well_number = int(plate_name_well_number) if plate_name_well_number is not None else None
        self.barcode_row_number = int(barcode_row_number) if barcode_row_number is not None else None
        self.barcode_well_number = int(barcode_well_number) if barcode_well_number is not None else None
        args.debug = debug


# def get_experiment(f, args):
#     if args.experiment:
#         return args.experiment
#     exp = os.path.basename(f).rstrip('.xlsx').rstrip('.xls')
#     return exp


def get_plate_blocks(ws, args):
    num = args.plate_delim_well_number
    prefix = args.plate_delim_prefix
    plate_blocks = []
    curr_plate = []
    for row in ws.rows:
        if not row[num - 1].value:
            curr_plate.append(row)
            continue
        if str(row[num - 1].value).startswith(prefix):
            plate_blocks.append(curr_plate)
            curr_plate = [row, ]
        else:
            curr_plate.append(row)
    if curr_plate:
        plate_blocks.append(curr_plate)
    return plate_blocks


def parse_barcodes(raw_plates, args):
    plate_info = []
    for rp in raw_plates:
        # parse the plate name
        name_row = args.plate_name_row_number
        name_well = args.plate_name_well_number - 1
        plate_name = rp[name_row][name_well].value
        if type(plate_name) in [float, int]:
            plate_name = str(int(rp[name_row][name_well].value))
        if args.prefix is not None:
            plate_name = args.prefix + plate_name
        if args.suffix is not None:
            plate_name += args.suffix
        # parse the plate barcode
        bc_row = args.barcode_row_number
        bc_well = args.barcode_well_number - 1
        barcode = str(rp[bc_row][bc_well].value)
        print('NAME: {}\nBARCODE: {}\n'.format(plate_name, barcode))
        plate_info.append({'name': plate_name, 'barcode': barcode})
    return plate_info



# def parse_plates(raw_plates, args):
#     plates = []
#     for i, rp in enumerate(raw_plates):
#         if args.plate_name_row_number is not None and args.plate_name_well_number is not None:
#             row_num = args.plate_name_row_number
#             well_num = args.plate_name_well_number - 1
#             plate_name = rp[row_num][well_num].value
#         else:
#             num = str(i + args.plate_numbering_start)
#             if len(num) < 2:
#                 num = '0' + num
#             plate_name = '{}{}'.format(args.plate_prefix, num)
#         logger.info('Processing plate: {}'.format(plate_name))
#         samples = parse_samples(rp, args)
#         wells = parse_plate_grid(rp)
#         plates.append(Plate(plate_name, wells, samples))
#     return plates


# def parse_samples(plate, args):
#     samples = []
#     row = plate[args.sample_row_number]
#     start = args.sample_well_number - 1
#     end = start + args.max_sample_number
#     step = args.sample_well_offset
#     for cell in row[start:end:step]:
#         samples.append(Sample(cell))
#     return samples


# def parse_plate_grid(raw_plate):
#     # trim empty rows from the end of the raw plate
#     rp = [list(row) for row in raw_plate]
#     rp_set = list(set([c.value for c in rp[-1]]))
#     while len(rp_set) == 1 and rp_set[0] is None:
#         rp = rp[:-1]
#         rp_set = list(set([c.value for c in rp[-1]]))
#     # find the column containing row labels
#     labelc = 0
#     for i in range(len(rp)):
#         column = [r[i] for r in rp]
#         cvals = [c.value for c in column if c.value is not None]
#         if len([c for c in cvals if str(c) in string.ascii_uppercase[:8]]) == 8:
#             labelc = i
#     gridc = labelc + 1
#     # find the row that contains the column labels
#     labelr = 0
#     for i, row in enumerate(raw_plate):
#         if len([r for r in row[labelc + 1:] if r.value in range(1, 13)]) == 12:
#             labelr = i
#     gridr = labelr + 1
#     # parse the cells from the plate grid
#     plate = []
#     for row in raw_plate[gridr:gridr + 8]:
#         rname = row[labelc].value
#         for cname, cell in enumerate(row[gridc:gridc + 12]):
#             well = Well(cell)
#             well.set_well(str(rname), str(cname + 1))
#             plate.append(well)
#     return plate


def write_output(plates, args):
    output = []
    for p in plates:
        output.append('{}\t{}'.format(p['name'], p['barcode']))
    open(args.output, 'w').write('\n'.join(output))

# def write_output(plates, experiment, args):
#     output = []
#     ohandle = open(args.output, 'w')
#     for i, plate in enumerate(plates):
#         for well in plate.wells:
#             output.append('{}-{}\t{}\t{}\t{}'.format(
#                 plate.name, well.well, well.sample, well.value, experiment))
#     ohandle.write('\n'.join(output))
#     # for i, plate in enumerate(plates):
#     #     ohandle = open(os.path.join(args.output, plate.name), 'w')
#     #     output = []
#     #     for well in plate.wells:
#     #         output.append('{}\t{}\t{}\t{}'.format(
#     #             well.well, well.sample, well.value, experiment))
#     #     ohandle.write('\n'.join(output))


def run(**kwargs):
    args = Args(**kwargs)
    global logger
    logger = log.get_logger('barcodes')
    main(args)


def run_standalone(args):
    global logger
    logger = log.get_logger('barcodes')
    main(args)


def main(args):
    for f in list_files(args.input):
        # experiment = get_experiment(f, args)
        wb = load_workbook(f)
        ws = wb[wb.get_sheet_names()[0]]
        plate_blocks = get_plate_blocks(ws, args)
        plural = '' if len(plate_blocks) <= 2 else 's'
        logger.info('\nFound {} plate{} in the input file'.format(len(plate_blocks) - 1, plural))
        # logger.info('Experiment name: {}\n'.format(experiment))
        # plates = parse_plates(plate_blocks[1:], args)
        plates = parse_barcodes(plate_blocks[1:], args)
        write_output(plates, args)
        logger.info('')


if __name__ == '__main__':
    args = parse_args()
    if args.log is None:
        args.log = os.path.join(args.output, 'barcodes.log')
    log.setup_logging(args.log)
    logger = log.get_logger('barcodes')
    main(args)
