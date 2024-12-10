#!/usr/bin/python3
__version__ = "1.0"

import coralme.builder
import coralme.core
import coralme.io

import sys
if sys.platform in ['win32']:
	pass
else:
	import coralme.solver.solver

import coralme.util

from . import _version
__version__ = _version.get_versions()['version']
