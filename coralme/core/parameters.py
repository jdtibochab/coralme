import logging
import numpy
import pint
import sympy
import tqdm
import coralme

class MEParameters():
	def __init__(self, model):
		self._model = model

		# Create a unit registry
		ureg = pint.UnitRegistry()
		model.unit_registry = ureg

		model.symbols = {}
		for var, unit in [('P', 'grams per gram'), ('R', 'grams per gram'), ('k_t', '1 per hour'), ('r_0', None), ('k^mRNA_deg', '1 per hour'), ('m_rr', 'gram per mmol'), ('m_aa', 'gram per mmol'), ('m_nt', 'gram per mmol'), ('f_rRNA', None), ('f_mRNA', None), ('f_tRNA', None), ('m_tRNA', 'gram per mmol'), ('k^default_cat', '1 per second'), ('temperature', 'K'), ('propensity_scaling', None), ('g_p_gdw_0', 'grams per gram'), ('g_per_gdw_inf', 'grams per gram'), ('b', '1 per hour**{:f}'.format(model.default_parameters[sympy.Symbol('d', positive = True)])), ('d', None)]:
			# WARNING: [b] is nominally per hour**d, but mu**d cannot be calculated if the types of mu and d are pint.Quantity
			if var == 'b':
				model.symbols[var] = ureg.Quantity(sympy.Symbol(var, positive = True))
			elif unit is None:
				model.symbols[var] = ureg.Quantity(sympy.Symbol(var, positive = True))
			else:
				model.symbols[var] = sympy.Symbol(var, positive = True) * ureg.parse_units(unit)

		# set growth rate symbolic variable
		mu = model.global_info['growth_key']
		model._mu = sympy.Symbol(mu, positive = True) * ureg.parse_units('1 per hour')
		# this allows the change of symbolic variables through the ME-model object
		model._mu_old = model._mu

		# derived parameters that are common throughout the ME-model
		# WARNING: The equations are written following O'Brien 2013 paper, no COBRAme documentation
		# https://www.embopress.org/doi/full/10.1038/msb.2013.52#supplementary-materials
		# Empirical relationship between measured ratio of RNA (R) to Protein (P)
		model.symbols['P'] # grams of amino acids per gDW := dimensionless
		model.symbols['R'] # grams of nucleotides per gDW := dimensionless

		# [R/P] = grams of nucleotides per grams of amino acids := dimensionless
		model.symbols['R/P'] = (model._mu / model.symbols['k_t']) + model.symbols['r_0'] # eq 1, page 15
		# [P/R] = grams of amino acids per grams of nucleotides := dimensionless
		model.symbols['P/R'] = 1. / model.symbols['R/P']

		# 70S ribosomes (page 16)
		# this is Ps in the supplementary material; [Ps] = millimoles of average amino acids per gDW per hour
		model.symbols['p_rate'] = model._mu * model.symbols['P'] / model.symbols['m_aa']
		# [R times f_rRNA] = grams of nucleotides in rRNA per gDW
		# this is nr in the supplementary material; [nr] = millimoles of nucleotides in rRNA per gDW
		model.symbols['n_ribo'] = model.symbols['R'] * model.symbols['f_rRNA'] / model.symbols['m_rr']

		# Hyperbolic ribosome catalytic rate
		model.symbols['c_ribo'] = model.symbols['m_rr'] / (model.symbols['m_aa'] * model.symbols['f_rRNA']) # eq 2, page 16
		# [kribo = p_rate / n_ribo] = millimoles of average amino acids per millimoles of nucleotides in rRNA per hour := per hour
		model.symbols['k_ribo'] = model.symbols['c_ribo'] * model._mu / model.symbols['R/P']
		# WARNING: the ribosome coupling coefficient in translation reactions is 'v_ribo' times protein length
		model.symbols['v_ribo'] = 1. / (1. * model.symbols['k_ribo'] / model._mu) # page 17

		# RNA Polymerase
		# WARNING: the RNAP coupling coefficient in transcription reactions is 'v_rnap' times RNA length
		model.symbols['v_rnap'] = 1. / (3. * model.symbols['k_ribo'] / model._mu) # page 17

		# mRNA coupling
		model.symbols['c_mRNA'] = model.symbols['m_nt'] / (model.symbols['f_mRNA'] * model.symbols['m_aa']) # page 19
		# Hyperbolic mRNA catalytic rate
		model.symbols['k_mRNA'] = 3 * model.symbols['c_mRNA'] * model._mu / model.symbols['R/P'] # 3 nt per aa

		# mRNA dilution, degradation, and translation
		model.symbols['alpha_1'] = model._mu / model.symbols['k^mRNA_deg']
		# WARNING: There is an error in O'Brien 2013; corrected in COBRAme docs
		model.symbols['alpha_2'] = model.symbols['R/P'] / (3 * model.symbols['alpha_1'] * model.symbols['c_mRNA'])
		# mRNA dilution, degradation, and translation
		model.symbols['rna_amount'] = model._mu / model.symbols['k_mRNA'] # == alpha_1 * alpha_2
		model.symbols['deg_amount'] = model.symbols['k^mRNA_deg'] / model.symbols['k_mRNA'] # == alpha_2

		# tRNA coupling
		model.symbols['c_tRNA'] = model.symbols['m_tRNA'] / (model.symbols['f_tRNA'] * model.symbols['m_aa']) # page 20
		# Hyperbolic tRNA efficiency
		model.symbols['k_tRNA'] = model.symbols['c_tRNA'] * model._mu / model.symbols['R/P']

		# Remaining Macromolecular Synthesis Machinery
		model.symbols['v^default_enz'] = 1. / (1. * (model.symbols['k^default_cat'].to('1 per hour')) / model._mu) # page 20, k^default_cat in 1/s

		# DNA replication (derivation not in documentation or supplementary material)
		# c = g_per_gdw_inf
		# a = g_p_gdw_0 - g_per_gdw_inf
		# g_p_gdw = (-a * gr ** d) / (b + gr ** d) + a + c, with a + c => g_p_gdw_0 - g_per_gdw_inf + g_per_gdw_inf <=> g_p_gdw_0
		model.symbols['dna_g_per_g'] = ((model.symbols['g_p_gdw_0'] - model.symbols['g_per_gdw_inf']) * model._mu.magnitude**model.symbols['d'] / (model.symbols['b'] + model._mu.magnitude**model.symbols['d'])) + model.symbols['g_p_gdw_0']

	@property
	def fundamental_equations(self):
		return self._model.symbols

	@fundamental_equations.setter
	def fundamental_equations(self, value: dict):
		self._model.symbols.update(value)

	# MetabolicReaction
	@staticmethod
	def coupling_coefficient_enzyme(obj, value):
		if isinstance(obj, coralme.core.reaction.MetabolicReaction):
			self._coupling_coefficient_enzyme = obj._model.mu * value.to('1 per hour')**-1

	# SubreactionData
	@staticmethod
	def coupling_coefficient_subreaction(obj, value):
		if isinstance(obj, coralme.core.processdata.SubreactionData):
			obj._coupling_coefficient_subreaction = obj._model.mu * value.to('1 per hour')**-1

	# tRNAData
	@property
	def _recalculate_all_synthetase_keff(obj, value):
		for obj in tqdm.tqdm(self._model.tRNA_data):
			obj._synthetase_keff = self._model.symbols['k^default_cat']

	@staticmethod
	def synthetase_keff(obj, value):
		if isinstance(obj, coralme.core.processdata.tRNAData):
			obj._synthetase_keff = value

	@property
	def _recalculate_all_coupling_coefficient_trna_keff(self):
		for obj in tqdm.tqdm(self._model.tRNA_data):
			obj._coupling_coefficient_trna_keff = self._model.symbols['k_tRNA']

	@staticmethod
	def coupling_coefficient_trna_keff(obj, value):
		if isinstance(obj, coralme.core.processdata.tRNAData):
			obj._coupling_coefficient_trna_keff = value

	@property
	def _recalculate_all_coupling_coefficient_trna_amount(self):
		for obj in tqdm.tqdm(self._model.tRNA_data):
			obj._coupling_coefficient_trna_amount = self._model.mu * obj.coupling_coefficient_trna_keff**-1

	@staticmethod
	def coupling_coefficient_trna_amount(obj, value):
		if isinstance(obj, coralme.core.processdata.tRNAData):
			obj._coupling_coefficient_trna_amount = value

	@property
	def _recalculate_all_coupling_coefficient_synthetase(self):
		self._coupling_coefficient_synthetase = self._model.mu * self.synthetase_keff.to('1 per hour')**-1 * (1 + self.coupling_coefficient_trna_amount)

	@staticmethod
	def coupling_coefficient_synthetase(obj, value):
		if isinstance(obj, coralme.core.processdata.tRNAData):
			obj._coupling_coefficient_synthetase = value

	# TranscriptionData
	@property
	def _recalculate_all_coupling_coefficient_rnapol(self):
		for obj in tqdm.tqdm(self._model.transcription_data):
			# this sets beta_transcription^RNAP (see overlead, page 8)
			obj._coupling_coefficient_rnapol = len(obj.nucleotide_sequence) * self._model.symbols['v_rnap']

	@staticmethod
	def coupling_coefficient_rnapol(obj, value):
		if isinstance(obj, coralme.core.processdata.TranscriptionData):
			obj._coupling_coefficient_rnapol = value # == len(nucleotide_sequence) * model.symbols['v_rnap']

	# TranslationData
	@property
	def _recalculate_all_coupling_coefficient_ribosome(self):
		# WARNING: k_ribo is the effective ribosomal translation rate (see overleaf, page 1)
		# WARNING: v_ribo is the coupling coefficient := mu/k_ribo (see overleaf, page 2)
		# this is beta_translation^ribosome (see overleaf, page 7)
		for obj in tqdm.tqdm(self._model.translation_data):
			obj._coupling_coefficient_ribosome = len(obj.translation) * self._model.symbols['v_ribo']

	@staticmethod
	def coupling_coefficient_ribosome(obj, value):
		if isinstance(obj, coralme.core.processdata.TranslationData):
			obj._coupling_coefficient_ribosome = value

	@property
	def _recalculate_all_coupling_coefficient_rna_synthesis(self):
		for obj in tqdm.tqdm(self._model.translation_data):
			obj._coupling_coefficient_rna_synthesis = self._model.symbols['rna_amount'] + self._model.symbols['deg_amount']

	@staticmethod
	def coupling_coefficient_rna_synthesis(obj, value):
		if isinstance(obj, coralme.core.processdata.TranslationData):
			obj._coupling_coefficient_rna_synthesis = value

	@property
	def _recalculate_all_coupling_coefficient_hydrolysis(self):
		for obj in tqdm.tqdm(self._model.translation_data):
			obj._coupling_coefficient_hydrolysis = ((len(obj.nucleotide_sequence) - 1) / 4.) * self._model.symbols['deg_amount']

	@staticmethod
	def coupling_coefficient_hydrolysis(obj, value):
		if isinstance(obj, coralme.core.processdata.TranslationData):
			obj._coupling_coefficient_hydrolysis = value
