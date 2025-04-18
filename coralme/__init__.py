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

def check_installed_packages():
	import importlib.metadata
	# installed with python:
	# ast, collections, functools, importlib, typing, warnings, copy, errno, glob, gzip, io, json, logging, math, operator, os
	imports = ['anyconfig', 'Bio', 'cobra', 'coralme', 'docplex', 'gurobipy', 'importlib_resources', 'jsonschema', 'numpy', 'openpyxl', 'optlang', 'pandas', 'pint', 'pyranges', 'pytest', 'python-libsbml', 'scipy', 'sympy', 'tqdm', 'versioneer', 'xlsxwriter']
	modules = {}
	for x in imports:
		try:
			if x == 'python-libsbml':
				modules[x] = __import__('libsbml')
			else:
				modules[x] = __import__(x)

			if x == 'Bio':
				print('{:s}=={:s}'.format('Biopython', importlib.metadata.version('Biopython')))
			else:
				print('{:s}=={:s}'.format(x, importlib.metadata.version(x)))
		except ImportError:
			print("Error importing {:s}.".format(x))
	return modules
