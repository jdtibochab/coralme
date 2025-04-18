import logging
import numpy
import pint
import sympy
import tqdm
import coralme

"""
| coralme†     | Symbol‡      | Value                      | Unit                          | Definition                                              |
|              |              | (Follows "Symbol" column)  |                               |                                                         |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| me.mu        | μ            | call model.optimize()      | h^-1                          | Specific growth rate = ln(2) / doubling time            |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| R            | R            | See R/P                    | g gDW^-1                      | Total cellular RNA mass per gram of biomass             |
| P            | P            | See R/P                    | g gDW^-1                      | Total cellular protein mass per gram of biomass         |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| k_t          | κ_t          | 4.5                        | h^-1                          | Slope (Scott et al., 2010)                              |
| r_0          | r_0          | 0.087                      | dimensionless                 | Intercept (Scott et al., 2010)                          |
| f_rRNA       | f_rRNA       | 0.86                       | dimensionless                 | Fraction of RNA that is rRNA                            |
| f_mRNA       | f_mRNA       | 0.02                       | dimensionless                 | Fraction of RNA that is mRNA                            |
| f_tRNA       | f_tRNA       | 0.12                       | dimensionless                 | Fraction of RNA that is tRNA                            |
| m_aa         | m_aa         | 0.109                      | g mmol^-1 (:= kDa)            | Molecular weight of average amino acid                  |
| m_nt         | m_nt         | 0.324                      | g mmol^-1 (:= kDa)            | Molecular weight of average mRNA nucleotide             |
| m_tRNA       | m_tRNA       | 25.0                       | g mmol^-1 (:= kDa)            | Molecular weight of average tRNA                        |
| m_rr         | m_rr         | 1453.0                     | g mmol^-1 (:= kDa)            | Molecular weight of rRNA per ribosome                   |
| k^mRNA_deg   | k^mRNA_deg   | 12.0                       | h^-1                          | First-order mRNA degradation constant                   |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| R/P          | R/P          | (μ κ_t^-1) + r_0           | dimensionless                 | RNA to protein ratio (Scott et al., 2010)               |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| n_ribo       | n_r          | R f_rRNA m_rr^-1           | mmol rRNA gDW^-1              | Concentration of rRNA (Assumption #rRNA = #ribosomes)   |
| p_rate       | v_trans_mRNA | μ P m_aa^-1                | mmol gDW^-1 h^-1              | Protein synthesis rate (aka, Ps)                        |
| §            | k'_ribo      |                            | h^-1                          | Average translation rate of an active ribosome          |
| §            | f_r          |                            | dimensionless                 | Fraction of active ribosomes                            |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| k_ribo       | k_ribo       | v_trans_mRNA n_r^-1        | h^-1                          | Effective ribosomal translation rate (aka, k'_ribo f_r) |
| v_ribo       | v_ribo       | μ k_ribo^-1                | dimensionless                 | Ribosome coup coeff per aa in translation rxns          |
| k_rnap       | k_RNAP       | 3 k_ribo                   | h^-1                          | RNA Polymerase transcription rate                       |
| v_rnap       | v_RNAP       | μ (k_RNAP)^-1              | dimensionless                 | RNA Polymerase per nt coup coeff in transcription rxns  |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| [mRNA]       | [mRNA]       | R f_mRNA m_nt^-1           | mmol gDW^-1                   | mmol of ribonucleotides in mRNA per gDW^-1              |
| v_dil_mRNA   | v_dil_mRNA   | μ [mRNA]                   | mmol gDW^-1 h^-1              | Dilution rate for mRNA                                  |
| v_deg_mRNA   | v_deg_mRNA   | k^mRNA_deg [mRNA]          | mmol gDW^-1 h^-1              | Degradation rate for mRNA                               |
| alpha_1      | α_1          | v_dil_mRNA v_deg_mRNA^-1   | dimensionless                 | Dil-to-deg coup coeff for mRNAs in translation rxns     |
| 3*deg_amount | α_2          | v_deg_mRNA v_trans_mRNA^-1 | dimensionless                 | Deg-to-trans coup coeff for mRNAs in translation rxns   |
| k_mRNA/3     | k_mRNA       | v_trans_mRNA [mRNA]^-1     | mmol prot (mmol mRNA)^-1 h^-1 | mRNA catalytic rate                                     |
| rna_amount   | α_1 α_2      | μ k_mRNA^-1                | dimensionless                 | Dil-to-trans coup coeff for mRNAs in translation rxns   |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| [tRNA]       | [tRNA]       | R f_tRNA m_tRNA^-1         | mmol gDW^-1                   | mmol of ribonucleotides in tRNA per gDW^-1              |
| v_dil_tRNA   | v_dil_tRNA   | μ [tRNA]                   | mmol gDW^-1 h^-1              | Dilution rate for tRNAs in translation rxns             |
| v_chrg_tRNA  | v_charg_tRNA | v_trans_mRNA               | mmol gDW^-1 h^-1              | Charging rate for tRNAs (We assume it is v_trans_mRNA)  |
| alpha_3      | α_3          | v_dil_tRNA v_charg_tRNA^-1 | dimensionless                 | Dil-to-charg coup coeff for tRNAs in translation rxns   |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|
| c_ribo       | c_ribo       | m_rr f_rRNA^-1 m_aa^-1     | dimensionless                 | Ribosome catalytic rate                                 |
| c_mRNA       | c_mRNA       | m_nt f_mRNA^-1 m_aa^-1     | dimensionless                 |                                                         |
| c_tRNA       | c_tRNA       | m_tRNA f_tRNA^-1 m_aa^-1   | dimensionless                 |                                                         |
| k_tRNA       | k_tRNA       | μ c_tRNA (R/P)^-1          | mmol prot (mmol tRNA)^-1 h^-1 | tRNA catalytic rate (aka, μ / α_3)                      |
|--------------|--------------|----------------------------|-------------------------------|---------------------------------------------------------|

# See overleaf
beta^rnap_transcription   := len(TU) * v_rnap

beta^ribosome_translation := len(protein) * v_ribo
beta^mRNA_translation     := beta^mRNA_dilution + beta^mRNA_degradation
beta^mRNA_dilution        := 1/3 α_1 α_2
beta^mRNA_degradation     := 1/3 α_2

beta^tRNA_translation     := mu / α_3
beta^tRNA ligase_charging := mu / free_parameter * (1 + beta^tRNA_translation)

Abbreviations:
aa         : amino acid
nt         : (ribo)nucleotide
rnx(s)     : reaction(s)
coup coeff : coupling coefficient
dil        : dilution
deg        : degradation
trans      : translation
charg      : charging of tRNA with an amino acid

† Access to the parameter or expression using model.symbols dictionary, e.g., model.symbols['p_rate']
‡ As in the documentation (overleaf)
§ Defined in documentation (overleaf, cobrame.readthedocs.io, and/or O'Brien et al., 2013), but never used
"""

class DefaultParameters(dict):
	def __init__(self, *args, **kwargs):
		"""Normalize keys on initialization."""
		super().__init__()  # Avoid passing args directly

		# Merge args and kwargs into a single dict
		initial_data = dict(*args, **kwargs)
		for key, value in initial_data.items():
			self[key] = value  # Triggers __setitem__

	def __setitem__(self, key, value):
		"""Ensure all keys are stored as sympy.Symbol with positive=True."""
		if isinstance(key, str):
			key = sympy.Symbol(key, positive=True)  # Convert string to Symbol
		elif isinstance(key, sympy.Symbol):
			key = sympy.Symbol(key.name, positive=True)  # Normalize existing Symbol

		super().__setitem__(key, value)

	def __getitem__(self, key):
		"""Ensure keys are normalized before retrieval."""
		if isinstance(key, str):
			key = sympy.Symbol(key, positive=True)  # Convert string to Symbol
		elif isinstance(key, sympy.Symbol):
			key = sympy.Symbol(key.name, positive=True)  # Normalize existing Symbol

		return super().__getitem__(key)

	def get(self, key, default=None):
		"""Retrieve the value for a given sympy.Symbol key."""
		if isinstance(key, str):  # Allow lookup by string name
			key = sympy.Symbol(key, positive=True)
		elif isinstance(key, sympy.Symbol):
			key = sympy.Symbol(key.name, positive=True)
		return super().get(key, default)

	def update(self, *args, **kwargs):
		"""Override update to ensure all keys are sympy.Symbol with positive=True."""
		new_data = dict(*args, **kwargs)
		converted_data = {
			(sympy.Symbol(k, positive=True) if isinstance(k, str) else sympy.Symbol(k.name, positive=True)): v
			for k, v in new_data.items()
		}
		super().update(converted_data)  # Call original dict update method

class MEParameters():
	def __init__(self, model):
		self._model = model

		# Create a unit registry
		ureg = pint.UnitRegistry()
		model.unit_registry = ureg

		# WARNING: DefaultParameters class ensures keys are sympy's Symbols and positive
		model.global_info['default_parameters'] = DefaultParameters({
			'k_t' : 4.5, # per hour
			'r_0' : 0.087, # dimensionless
			'k^mRNA_deg' : 12.0, # per hour
			'm_rr' : 1453.0, # kDa = g per millimole
			'm_aa' : 0.109, # kDa = g per millimole
			'm_nt' : 0.324, # kDa = g per millimole
			'f_rRNA' : 0.86, # dimensionless, between 0 and 1
			'f_mRNA' : 0.02, # dimensionless, between 0 and 1
			'f_tRNA' : 0.12, # dimensionless, between 0 and 1
			'm_tRNA' : 25.0, # kDa = g per millimole
			'k^default_cat' : 65.0, # per second, internally converted to per hour
			'temperature' : 37.0, # kelvin
			'propensity_scaling' : 0.45, # dimensionless
			# DNA replication; see dna_replication.percent_dna_template_function
			'g_p_gdw_0' : 0.059314110730022594, # dimensionless
			'g_per_gdw_inf' : 0.02087208296776481, # dimensionless
			# WARNING: [b] is nominally per hour**d, but mu**d cannot be calculated if mu and d types are pint.Quantity
			'b' : 0.1168587392731988,
			'd' : 3.903641432780327 # dimensionless
			})

		# set growth rate symbolic variable
		mu = model.global_info['growth_key']
		model._mu = sympy.Symbol(mu, positive = True) * ureg.parse_units('1 per hour')
		# this allows the change of symbolic variables through the ME-model object
		model._mu_old = model._mu

		# set symbols and derived equations
		model.symbols = {}
		for var, unit in [('P', None), ('R', None), ('k_t', '1 per hour'), ('r_0', None), ('k^mRNA_deg', '1 per hour'), ('m_rr', 'gram per mmol'), ('m_aa', 'gram per mmol'), ('m_nt', 'gram per mmol'), ('f_rRNA', None), ('f_mRNA', None), ('f_tRNA', None), ('m_tRNA', 'gram per mmol'), ('k^default_cat', '1 per second'), ('temperature', 'K'), ('propensity_scaling', None), ('g_p_gdw_0', 'grams per gram'), ('g_per_gdw_inf', 'grams per gram'), ('b', '1 per hour**{:f}'.format(model.global_info['default_parameters']['d'])), ('d', None)]:
			# WARNING: [b] is nominally per hour**d, but mu**d cannot be calculated if the types of mu and d are pint.Quantity
			if var == 'b':
				model.symbols[var] = ureg.Quantity(sympy.Symbol(var, positive = True))
			elif unit is None:
				model.symbols[var] = ureg.Quantity(sympy.Symbol(var, positive = True))
			else:
				model.symbols[var] = sympy.Symbol(var, positive = True) * ureg.parse_units(unit)

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
		# this is v_trans_mRNA in the overleaf document
		model.symbols['p_rate'] = model._mu * model.symbols['P'] / model.symbols['m_aa']
		# [R times f_rRNA] = grams of nucleotides in rRNA per gDW
		# this is nr in the supplementary material; [nr] = millimoles of nucleotides in rRNA per gDW
		model.symbols['n_ribo'] = model.symbols['R'] * model.symbols['f_rRNA'] / model.symbols['m_rr']

		# definitions
		model.symbols['[mRNA]'] = model.symbols['R'] * model.symbols['f_mRNA'] / model.symbols['m_nt']
		model.symbols['[tRNA]'] = model.symbols['R'] * model.symbols['f_tRNA'] / model.symbols['m_tRNA']
		model.symbols['[rRNA]'] = model.symbols['R'] * model.symbols['f_rRNA'] / model.symbols['m_rr']
		model.symbols['v_dil_mRNA'] = model._mu * model.symbols['[mRNA]']
		model.symbols['v_deg_mRNA'] = model.symbols['k^mRNA_deg'] * model.symbols['[mRNA]']
		model.symbols['v_dil_tRNA'] = model._mu * model.symbols['[tRNA]']
		model.symbols['v_charg_tRNA'] = model.symbols['p_rate'] # Assumption

		# Hyperbolic ribosome catalytic rate
		model.symbols['c_ribo'] = model.symbols['m_rr'] / (model.symbols['f_rRNA'] * model.symbols['m_aa']) # eq 2, page 16
		# [kribo = p_rate / n_ribo] = millimoles of average amino acids per millimoles of nucleotides in rRNA per hour := per hour
		model.symbols['k_ribo'] = model.symbols['c_ribo'] * model._mu / model.symbols['R/P']
		# WARNING: the ribosome coupling coefficient in translation reactions is 'v_ribo' times protein length
		model.symbols['v_ribo'] = model._mu / model.symbols['k_ribo']  # page 17

		# RNA Polymerase
		model.symbols['k_rnap'] = 3. * model.symbols['k_ribo'] # Assumption
		# WARNING: the RNAP coupling coefficient in transcription reactions is 'v_rnap' times RNA length
		model.symbols['v_rnap'] = model._mu / model.symbols['k_rnap'] # page 17

		# mRNA coupling
		model.symbols['c_mRNA'] = model.symbols['m_nt'] / (model.symbols['f_mRNA'] * model.symbols['m_aa']) # page 19
		# Hyperbolic mRNA catalytic rate
		model.symbols['k_mRNA'] = 3. * model.symbols['c_mRNA'] * model._mu / model.symbols['R/P'] # 3 nt per aa

		# mRNA dilution, degradation, and translation
		model.symbols['alpha_1'] = model._mu / model.symbols['k^mRNA_deg']
		# WARNING: There is an error in O'Brien 2013; corrected in COBRAme docs
		model.symbols['alpha_2'] = model.symbols['R/P'] / (3. * model.symbols['alpha_1'] * model.symbols['c_mRNA'])
		# mRNA dilution, degradation, and translation
		model.symbols['rna_amount'] = model._mu / model.symbols['k_mRNA'] # == alpha_1 * alpha_2
		model.symbols['deg_amount'] = model.symbols['k^mRNA_deg'] / model.symbols['k_mRNA'] # == alpha_2

		# tRNA coupling
		model.symbols['c_tRNA'] = model.symbols['m_tRNA'] / (model.symbols['f_tRNA'] * model.symbols['m_aa']) # page 20
		# Hyperbolic tRNA efficiency
		model.symbols['k_tRNA'] = model.symbols['c_tRNA'] * model._mu / model.symbols['R/P']
		model.symbols['alpha_3'] = model.symbols['v_dil_tRNA'] / model.symbols['v_charg_tRNA']

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

	# Recalculate all
	@property
	def _recalculate_all_symbolic_stoichiometries(self):
		# tRNAData
		self._recalculate_all_synthetase_keff
		self._recalculate_all_coupling_coefficient_trna_keff
		self._recalculate_all_coupling_coefficient_trna_amount
		self._recalculate_all_coupling_coefficient_synthetase
		# TranscriptionData
		self._recalculate_all_coupling_coefficient_rnapol
		# TranslationData
		self._recalculate_all_coupling_coefficient_ribosome
		self._recalculate_all_coupling_coefficient_rna_synthesis
		self._recalculate_all_coupling_coefficient_hydrolysis

		# update
		for rxn in self._model.reactions.query('translation_'):
			rxn.update()
		for rxn in self._model.reactions.query('transcription_'):
			rxn.update()
		for rxn in self._model.reactions.query('charging_'):
			rxn.update()

	# MetabolicReaction
	@staticmethod
	def coupling_coefficient_enzyme(obj, value):
		if isinstance(obj, coralme.core.reaction.MetabolicReaction):
			self._coupling_coefficient_enzyme = obj._model.mu * value.to('1 per hour')**-1 # mu/k_eff

	# SubreactionData
	@staticmethod
	def coupling_coefficient_subreaction(obj, value):
		if isinstance(obj, coralme.core.processdata.SubreactionData):
			obj._coupling_coefficient_subreaction = obj._model.mu * value.to('1 per hour')**-1 # mu/k_eff

	# tRNAData
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
			obj._coupling_coefficient_trna_amount = self._model.mu * obj._coupling_coefficient_trna_keff**-1

	@staticmethod
	def coupling_coefficient_trna_amount(obj, value):
		if isinstance(obj, coralme.core.processdata.tRNAData):
			obj._coupling_coefficient_trna_amount = value

	@property
	def _recalculate_all_synthetase_keff(self):
		for obj in tqdm.tqdm(self._model.tRNA_data):
			obj._synthetase_keff = self._model.symbols['k^default_cat']

	@staticmethod
	def synthetase_keff(obj, value):
		if isinstance(obj, coralme.core.processdata.tRNAData):
			obj._synthetase_keff = value

	@property
	def _recalculate_all_coupling_coefficient_synthetase(self):
		for obj in tqdm.tqdm(self._model.tRNA_data):
			obj._coupling_coefficient_synthetase = self._model.mu * obj._synthetase_keff.to('1 per hour')**-1 * (1 + obj._coupling_coefficient_trna_amount)

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

	@staticmethod
	def check_parameter(value):
		if not value > 0.:
			raise ValueError('The coupling coefficient cannot be negative or zero.')

		if isinstance(value, pint.Quantity):
			return float(value.to('1 per second').magnitude)
		# WARNING: to check for numpy.int or numpy.float types, use numpy.issubdtype per type, i.e., numpy.integer and numpy.floating
		elif isinstance(float(value), float):
			return float(value)
		else:
			raise NotImplementedError
