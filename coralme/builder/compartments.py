import coralme
import coralme
import collections

def _return_compartments_of_complexes(model, cplx):
	try:
		data = model.process_data.get_by_id(cplx.id)
	except KeyError:
		data = coralme.builder.helper_functions.get_base_complex_data(model, cplx.id)

	mem_dict = collections.defaultdict(int)
	for s in data.stoichiometry:
		if '_Inner_Membrane' in s:
			mem_dict['im'] += 1
		elif '_Outer_Membrane' in s:
			mem_dict['om'] += 1
		elif '_Periplasm' in s:
			mem_dict['p'] += 1

	# if no membrane associated with membrane subunit, assume complex is cytosolic
	if len(mem_dict) == 0:
		return 'c'
	# if only one membrane is represented in protein subunits, use this membrane for the compartment
	elif len(mem_dict) == 1:
		return mem_dict.popitem()[0]
	# if multiple membrane compartments are represented, use generic "m" for "membrane" for now
	else:
		return 'm'

def add_compartments_to_model(model):
	"""First adds compartments based on suffix of metabolite ID. If metabolite
	is a complex, the protein subunit stoichiometry is used to infer
	compartment. All remaining metabolites without a compartment suffix (RNAs,
	generic metabolites, nonmembrane proteins, etc.) are assumed to be
	cytosolic"""

	for met in model.metabolites:
		if met.compartment:
			continue

		if isinstance(met, coralme.core.component.Constraint):
			met.compartment = 'mc'
		elif '_Inner_Membrane' in met.id:
			met.compartment = 'im'
		elif '_Outer_Membrane' in met.id:
			met.compartment = 'om'
		elif '_Periplasm' in met.id or met.id.endswith('_p'):
			met.compartment = 'p'
		elif met.id.endswith('_e'):
			met.compartment = 'e'
		elif met.id.endswith('_c'):
			met.compartment = 'c'
		elif isinstance(met, coralme.core.component.Complex):
			met.compartment = _return_compartments_of_complexes(model, met)
		else:
			met.compartment = 'c'
