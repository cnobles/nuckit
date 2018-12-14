__author__ = "Christopher Nobles"
__license__ = "GPL3+"

import os
import re
import sys
import csv

from pathlib import Path
from pkg_resources import get_distribution

from semantic_version import Version

__version__ = str(Version.coerce(get_distribution('nuckit').version))
