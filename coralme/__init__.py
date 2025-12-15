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

def citation(format = 'bibtex'):
	bibtex = """
	@article{coralME,
	title   = {Metabolism and gene expression models for the microbiome reveal how diet and metabolic dysbiosis impact disease},
	journal = {Cell Systems},
	pages   = {101451},
	year    = {2025},
	issn    = {2405-4712},
	doi     = {10.1016/j.cels.2025.101451},
	url     = {https://www.sciencedirect.com/science/article/pii/S2405471225002844},
	author  = {Juan D. Tibocha-Bonilla and Rodrigo Santib{\\'a}{\\~n}ez-Palominos and Yuhan Weng and Manish Kumar and Karsten Zengler},
	keywords = {gut microbiome, gut dysbiosis, models of metabolism and gene expression, predictive microbiome modeling, systems biology, genome-scale metabolic modeling},
	abstract = {Summary
	The gut microbiome plays a critical role in human health, spurring extensive research using multi-omic technologies. Although these tools offer valuable insights, they often fall short in capturing the complexity of microbial interactions that associate with disease onset, progression, and treatment. Thus, integration of multi-omics datasets with metabolic models is needed to predict associations between microbial activity and disease. Here, we automated the reconstruction of 495 metabolic and gene expression models (ME-models), overcoming the main limitation preventing the wide use of this approach. We integrated them with multi-omics data from patients with inflammatory bowel disease (IBD), identifying taxa associated with variations in amino acids, short-chain fatty acids, and pH in the gut of IBD patients. In general, this approach provides testable hypotheses of the metabolic activity of the gut microbiota, and the automated pipeline opens the opportunity to study microbial interactions in other biologically relevant settings using ME-models.}
	}
	"""
	if format == 'bibtex':
		print(bibtex.replace('\t', ''))
	else:
		return NotImplementedError
