import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg', warn=False)

from scPCR.demultiplex import run as demultiplex
from scPCR.deconvolute import run as deconvolute
