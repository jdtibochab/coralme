import os
import sys
import copy
import json
import sympy
import jsonschema
from collections import OrderedDict

import logging
log_format = '%(asctime)s %(message)s' #%(clientip)-15s %(user)-8s

import pandas
import cobra
import coralme

try:
	# If cannot import SymbolicParameter, assume using cobrapy versions <= 0.5.11
	from optlang.interface import SymbolicParameter
except ImportError:
	from cobra.io.json import save_json_model
	from cobra.io.dict import _metabolite_from_dict
else:
	from cobra.io.json import _metabolite_from_dict as metabolite_from_dict
	from cobra.io.json import save_json_model

cur_dir = os.path.dirname(os.path.abspath(__file__))

def get_schema():
	"""
	Load JSON schema for ME-model JSON saving/loading

	Returns
	-------
	dict
		JSONSCHEMA

	"""
	with open(os.path.join(cur_dir, 'JSONSCHEMA'), 'r') as f:
		return json.load(f)

def save_json_me_model(model, file_name, sort = True, compress = False):
	"""
	Save a full JSON version of the ME-model. Saving/loading a model in this
	format can then be loaded to return a ME-model identical to the one saved,
	which retains all ME-model functionality.

	Parameters
	----------

	model : :class:`coralme.core.model.MEModel`
		A full ME-model

	file_name : str or file-like object
		Filename of the JSON output or an open json file

	"""

	#should_close = False
	#if isinstance(file_name, str):
		#file_name = open(file_name, 'w')
		#should_close = True

	# set logger
	log = logging.getLogger() # root logger
	for hdlr in log.handlers[:]: # remove all old handlers
		log.removeHandler(hdlr)

	# Old code works in a separate script; but it works if we remove the old handler
	logging.basicConfig(level = logging.WARNING, format = log_format)
	#log.addHandler(logging.StreamHandler(sys.stdout))
	logging.captureWarnings(True)

	model_dict = coralme.io.dict.me_model_to_dict(model)

	# cannot be before me_model_from_dict
	if 'default_parameters' in model.global_info:
		model_dict['global_info']['default_parameters'] = { str(k):v for k,v in model_dict['global_info']['default_parameters'].items() }

	logging.shutdown()

	# Confirm that dictionary representation of model adheres to JSONSCHEMA
	try:
		jsonschema.validate(model_dict, get_schema())
	except jsonschema.ValidationError:
		raise Exception('Must pass valid ME-model json file')

	#json.dump(model_dict, file_name)

	#if should_close:
		#file_name.close()

	def set_default(obj):
		if isinstance(obj, set):
			return list(obj)
		if isinstance(obj, pandas.DataFrame):
			return obj.to_dict()
		#raise TypeError

	json_str = json.dumps(model_dict, indent = 2, sort_keys = sort, default = set_default)
	if compress:
		import gzip
		with gzip.open(file_name + '.gz', 'w') as outfile:
			outfile.write(json_str.encode('utf-8'))
	else:
		with open(file_name, 'w') as outfile:
			outfile.write(json_str)

def load_json_me_model(file_name):
	"""
	Load a full JSON version of the ME-model. Loading a model in this format
	will return a ME-model identical to the one saved, which retains all
	ME-model functionality.

	Parameters
	----------
	file_name : str or file-like object
		Filename of the JSON output or an open json file

	Returns
	-------
	:class:`coralme.core.model.MEModel`
		A full ME-model

	"""
	# set logger
	log = logging.getLogger() # root logger
	for hdlr in log.handlers[:]: # remove all old handlers
		log.removeHandler(hdlr)

	# Old code works in a separate script; but it works if we remove the old handler
	logging.basicConfig(filename = 'MELoader.log', filemode = 'w', level = logging.WARNING, format = log_format)
	log.addHandler(coralme.builder.main.ListHandler([]))
	#log.addHandler(logging.StreamHandler(sys.stdout))
	logging.captureWarnings(True)

	if isinstance(file_name, str):
		with open(file_name, 'r') as f:
			model_dict = json.load(f)
	else:
		model_dict = json.load(file_name)

	try:
		jsonschema.validate(model_dict, get_schema())
	except jsonschema.ValidationError:
		raise Exception('Must pass valid ME-model json file')

	model = coralme.io.dict.me_model_from_dict(model_dict)

	logging.shutdown()

	# We will remove duplicates entries in the log output
	with open('MELoader.log', 'w') as outfile:
		logger = log.handlers[1].log_list

		tmp = pandas.DataFrame(logger)
		for idx, data in tmp.drop_duplicates(subset = 1).iterrows():
			outfile.write('{:s} {:s}\n'.format(data[0], data[1]))

	return model

# -----------------------------------------------------------------------------
# Functions below here facilitate json dumping/loading of reduced ME-models
# without all process_data/reaction info intact.

def save_reduced_json_me_model(me0, file_name):
	"""
	Save a stripped-down JSON version of the ME-model. This will exclude all of
	ME-model information except the reaction stoichiometry information and the
	reaction bounds. Saving/loading a model in this format will thus occur much
	quicker, but limit the ability to edit the model and use most of its
	features.

	Parameters
	----------
	me0 : :class:`coralme.core.model.MEModel`
		A full ME-model

	file_name : str or file-like object
		Filename of the JSON output

	"""
	me = copy.deepcopy(me0)

	for rxn in me.reactions:
		for met in rxn.metabolites:
			s = rxn._metabolites[met]
			if isinstance(s, sympy.Basic):
				rxn._metabolites[met] = str(s)
		if isinstance(rxn.lower_bound, sympy.Basic):
			rxn.lower_bound = str(rxn.lower_bound)
		if isinstance(rxn.upper_bound, sympy.Basic):
			rxn.upper_bound = str(rxn.upper_bound)

	for met in me.metabolites:
		if isinstance(met._bound, sympy.Basic):
			met._bound = str(met._bound)

	save_json_model(me, file_name)

def load_reduced_json_me_model(file_name):
	"""
	Load a stripped-down JSON version of the ME-model. This will exclude all of
	ME-model information except the reaction stoichiometry information and the
	reaction bounds. Saving/loading a model in this format will thus occur much
	quicker, but limit the ability to edit the model and use most of its
	features.

	Parameters
	----------
	file_name : str or file-like object
		Filename of the JSON ME-model

	Returns
	-------
	:class:`cobra.core.model.Model`
		COBRA Model representation of the ME-model. This will not include
		all of the functionality of a :class:`~coralme.core.model.MEModel` but
		will solve identically compared to the full model.
	"""
	if isinstance(file_name, str):
		with open(file_name, 'r') as f:
			obj = json.load(f)
	else:
		obj = json.load(file_name)

	model = cobra.Model()

	# If cannot import SymbolicParameter, assume using cobrapy
	# versions <= 0.5.11. If versions >= 0.8.0 are used, a ME-model interface
	# must be assigned as the solver interface
	try:
		from optlang.interface import SymbolicParameter
	except ImportError:
		pass
	else:
		model.solver = coralme.util.me_model_interface

	default_reactions = [i.id for i in model.reactions]

	for k, v in obj.items():
		if k in {'id', 'name'}:
			setattr(model, k, v)

	def _reaction_from_dict(reaction, model):
		new_reaction = cobra.Reaction()
		for k, v in reaction.items():
			if k in {'objective_coefficient', 'reversibility', 'reaction'}:
				continue
			elif k == 'metabolites':
				new_reaction.add_metabolites(OrderedDict(
					(model.metabolites.get_by_id(str(met)),
					 coralme.io.dict.get_numeric_from_string(coeff))
					for met, coeff in v.items()))
			elif k in {'upper_bound', 'lower_bound'}:
				v = coralme.io.dict.get_numeric_from_string(v)
				setattr(new_reaction, k, v)
			else:
				setattr(new_reaction, k, v)
		return new_reaction

	model.add_metabolites([metabolite_from_dict(metabolite) for metabolite in obj['metabolites']])

	new_reactions = [_reaction_from_dict(reaction, model) for reaction in obj['reactions']]

	model.remove_reactions(default_reactions)
	model.add_reactions(new_reactions)

	return model
