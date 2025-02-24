import copy
import sympy
import coralme
from collections import defaultdict

def stringify(element, number):
	number = int(number)
	return element if number == 1 else (element + str(number).rstrip('.'))

def elements_to_formula(obj, elements):
	#obj.formula = ''.join(stringify(e, n) for e, n in sorted(iteritems(elements)))
	obj.formula = ''.join(stringify(e, n) for e, n in sorted(elements.items()) if n != 0)

def eval_reaction_at_growth_rate(reaction, growth_rate):
	for key, value in reaction._metabolites.items():
		if isinstance(value, sympy.Basic):
			reaction._metabolites[key] = value.subs(coralme.util.mu, growth_rate)
		if isinstance(reaction, coralme.TranslationReaction) and (isinstance(key, coralme.Ribosome) or isinstance(key, coralme.TranscribedGene)):
			reaction._metabolites[key] = 0.
		if isinstance(reaction, coralme.TranscriptionReaction) and isinstance(key, coralme.RNAP):
			reaction._metabolites[key] = 0.
	return reaction

def get_elements_from_process_data(reaction, process_data, elements):
	"""
	If a modification is required to form a functioning macromolecule,
	update the element dictionary accordingly.
	"""

	#for sub, count in iteritems(process_data.subreactions):
	for sub, count in process_data.subreactions.items():
		sub_obj = reaction._model.process_data.get_by_id(sub)
		for e, n in sub_obj.element_contribution.items():
			elements[e] += n * count
	return elements

def check_transcription_mass_balance(reaction):
	transcription_data = reaction.transcription_data
	if len(transcription_data.excised_bases) > 0:
		return reaction.check_me_mass_balance()

	# Store and reset the formula strings for each RNA
	conserve_dict = {}
	for rna in transcription_data.RNA_products:
		rna_obj = reaction._model.metabolites.get_by_id(rna)
		conserve_dict[rna_obj] = rna_obj.formula
		rna_obj.formula = ''
		rna_obj.elements = {}

	elements = defaultdict(int)

	nucleotide_count = {
		#i.replace('tp_c', 'mp_c'): v for i, v in iteritems(transcription_data.nucleotide_count)
		i.replace('tp_c', 'mp_c'): v for i, v in transcription_data.nucleotide_count.items()
		}

	for nuc, value in nucleotide_count.items():
		nuc_obj = reaction._model.metabolites.get_by_id(nuc)
		for e, n in nuc_obj.elements.items():
			elements[e] += value * n

	# Remove -OH for each
	elements['H'] -= len(transcription_data.nucleotide_sequence)
	elements['O'] -= len(transcription_data.nucleotide_sequence)

	elements_to_formula(rna_obj, elements)
	mass_balance = reaction.check_me_mass_balance()

	for rna_obj, formula in conserve_dict.items():
		rna_obj.formula = formula

	return mass_balance

def check_me_model_mass_balance(model0):
	model = copy.deepcopy(model0)

	def should_skip(reaction):
		# Exchanges and demand reactions will not be mass balanced
		if reaction.id.startswith('DM_') or reaction.id.startswith('EX_'):
			return True
		# Generic metabolites do not have a formula and Dummy reaction will not
		# be mass balanced
		elif 'to_generic' in reaction.id or 'dummy_reaction' in reaction.id:
			return True
		# Global constraints do not have a formula (for summary variables)
		# tRNA charging reactions are mass balance checked by confirming
		# mass balance of each modification/subreaction
		elif isinstance(reaction, (coralme.SummaryVariable, coralme.tRNAChargingReaction)):
			return True
		else:
			return False

	# Set k_deg = 0 so model will be feasible at growth rate = 0
	model.global_info['k_deg'] = 0
	for r in model.reactions:
		if isinstance(r, coralme.TranslationReaction):
			r.update()

	output = {}
	for r in model.reactions:
		me_reaction = True
		if should_skip(r):
			continue
		if not isinstance(r, coralme.MEReaction) and r.check_mass_balance():
			output[r.id] = r.check_mass_balance()
			me_reaction = False
		elif not isinstance(r, coralme.MEReaction):
			me_reaction = False

		if me_reaction:
			# Lipoprotein metabolites can have different formulas based on
			# reaction used to synthesized them (TODO: How to handle this?)
			if 'lipid_modification' in r.id:
				r.update()
			# Check mass balance at growth rate = 0
			eval_reaction_at_growth_rate(r, 0.)
			# Transcription reactions must be handled differently to account
			# for multiple RNAs per TU
			if isinstance(r, coralme.TranscriptionReaction) and check_transcription_mass_balance(r):
				output[r.id] = check_transcription_mass_balance(r)
			elif not isinstance(r, coralme.TranscriptionReaction) and r.check_me_mass_balance():
				output[r.id] = r.check_me_mass_balance()

	# Each individual subreaction data should be mass balanced, if necessary
	for data in model.process_data:
		if hasattr(data, 'element_contribution'):
			# Skip trivial modifications (only one metabolite as reactant)
			# These are simple modifications
			if len(data.stoichiometry) == 1 and list(data.stoichiometry.values())[0] < 0:
				continue

			calculated_element_contribution = {}
			for key, value in data.calculate_element_contribution().items():
				if value != 0:
					calculated_element_contribution[key] = value

			set_element_contribution = {}
			for key, value in data._element_contribution.items():
				if value != 0:
					set_element_contribution[key] = value

			if calculated_element_contribution != set_element_contribution:
				txt = 'Calculated element contribution ({:s}) not equal to user defined element contribution({:s})'
				output[data.id] = txt.format(calculated_element_contribution, set_element_contribution)

	return output
