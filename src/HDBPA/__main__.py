
import sys

from . import HDBPA, utils

args = utils._parse_args(sys.argv[1:])

HDBPA.HDBPA(args)
