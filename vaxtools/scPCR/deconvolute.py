#!/usr/bin/env python
# filename: deconvolute.py


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


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cdna', dest='cdna_plates', required=True,
                        help='''A file containing subject assignments for each well \
                        in the cDNA plate, with the following fields:
                        Plate  Well  Subject
                        separated by whitespace. One entry per line. \
                        Required.''')
    parser.add_argument('-1', '--pcr1', dest='pcr1_plates', required=True,
                        help='''A file containing PCR1 plate barcodes and the barcode \
                        of the cDNA template plate, in the following format:
                        cDNA_plate_barcode  PCR1_plate_barcode
                        separated by whitespace. One plate per line. \
                        Required.''')
    parser.add_argument('-2', '--pcr2', dest='pcr2_plates', default=None,
                        help='''A file containing PCR2 plate barcodes and the barcode \
                        of the cDNA template plate along with the position of the PCR1 plate \
                        in the PCR2 plate, in the following format:
                        PCR1_plate_barcode  PCR2_plate_barcode  Position
                        separated by whitespace. One plate per line. \
                        If the PCR1 plate is 96-well and the PCR2 plate is 384-well, the following \
                        position options are valid: 'top-odd', 'top-even', 'bottom-odd', 'bottom-even'. \
                        If PCR1 and PCR2 plates are of the same format (both 96-well or both 384-well), \
                        the position should be 'whole-plate'. \
                        Required.''')
    parser.add_argument('-i', '--indexing-PCR', dest='ipcr_plates', default=None,
                        help='''A file containing PCR2 plate barcodes and the barcode \
                        of the cDNA template plate along with the position of the PCR1 plate \
                        in the PCR2 plate, in the following format:
                        PCR2_plate_barcode  indexing_PCR_plate_barcode  Position
                        separated by whitespace. One plate per line. \
                        If the PCR2 plate is 96-well and the indexing PCR plate is 384-well, the following \
                        position options are valid: 'top-odd', 'top-even', 'bottom-odd', 'bottom-even'. \
                        If PCR2 and indexing PCR plates are of the same format (both 96-well or both 384-well), \
                        the position should be 'whole-plate'. \
                        Required.''')
    parser.add_argument('--pcr1-plate-type', dest='pcr1_plate_type', default=96,
                        help="The plate type of the PCR1 plates. \
                        Valid options are '96' and '384'. \
                        Default is '96'.")
    parser.add_argument('--pcr2-plate-type', dest='pcr2_plate_type', default=384,
                        help="The plate type of the PCR2 plates. \
                        Valid options are '96' and '384'. \
                        Default is '384'.")
    parser.add_argument('--indexing-pcr-plate-type', dest='indexing_pcr_plate_type', default=384,
                        help="The plate type of the indexing PCR plates. \
                        Valid options are '96' and '384'. \
                        Default is '384'.")




def run(**kwargs):
    pass



def main(*args):
    pass


if __name__ == '__main__':
    main()