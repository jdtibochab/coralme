import copy
import re
import pickle
import typing
import collections

import logging
log = logging.getLogger(__name__)

# install by the user
import tqdm
bar_format = '{desc:<75}: {percentage:.1f}%|{bar:10}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
import numpy
import pandas
import pint
import scipy
import sympy
import cobra
import coralme
import sys

# due to a circular import
from coralme.core.component import Metabolite as Metabolite
from coralme.core.reaction import MEReaction as MEReaction

def _update(MEReaction):
	"""updates all component reactions"""
	MEReaction.update()
	return None

import types
class MEModel(cobra.core.object.Object):
	def __init__(self, id_or_model = 'coralME', name = 'coralME', mu = 'mu'):
		#cobra.Model.__init__(self, name)
		# to avoid setting the solver interface to gurobi or any other
		cobra.core.object.Object.__init__(self, id_or_model, name = name)

		# simulation methods in optimization.py
		self.optimize = types.MethodType(coralme.core.optimization.optimize, self)
		self.optimize_windows = types.MethodType(coralme.core.optimization.optimize_windows, self)
		self.feasibility = types.MethodType(coralme.core.optimization.feasibility, self) # qminos
		self.feas_gurobi = types.MethodType(coralme.core.optimization.feas_gurobi, self) # gurobi
		self.feas_cplex = types.MethodType(coralme.core.optimization.feas_cplex, self) # cplex
		self.construct_lp_problem = types.MethodType(coralme.core.optimization.construct_lp_problem, self)
		self.fva = types.MethodType(coralme.core.optimization.fva, self)

		self.model_version = coralme.__version__

		self.global_info = {
			'domain' : 'Prokaryote',
			'growth_key' : mu,
			'ME-Model-ID' : id_or_model,

			'dnapol_id' : 'DNAP',
			'ribosome_id' : 'ribosome',
			'dummy_rxn_id' : 'dummy_reaction',
			'degradosome_id' : 'RNA_degradosome',
			'mg2_per_ribosome' : 171,
			'amino_acid_loader' : 'generic_Tuf',
			'feature_types' : [ 'CDS', 'rRNA', 'tRNA', 'ncRNA', 'tmRNA', 'misc_RNA' ],

			'electron_transfers' : {
				'cytochromes' : [],
				'ferredoxins' : [],
				'thioredoxins': [],
				'glutaredoxins': [],
				'flavodoxins': [],
				'peroxiredoxins': [],
				},

			# analysis
			'add_lipoproteins' : False, #
			'add_translocases' : True, # actually, assign CPLX_dummy to missing enzymes
			'include_pseudo_genes' : False,
			'run_bbh_blast' : True,

			# TODO: We should test if the user set this correctly as { codon : { amino_acid : tRNAs }}
			'genetic_recoding' : {},

			'peptide_release_factors' : {
				'UAG': 'PrfA_mono',
				'UGA': 'PrfB_mono',
				'UAA': 'generic_RF',
				},

			'transcription_subreactions' : {
				'Transcription_normal_rho_independent' : '',
				'Transcription_normal_rho_dependent' : 'atp_hydrolysis_rho',
				'Transcription_stable_rho_independent' : '',
				'Transcription_stable_rho_dependent' : 'atp_hydrolysis_rho',
				},

			'translation_subreactions' : {
				'Translation_initiation_factor_InfA' : '',
				'Translation_initiation_factor_InfC' : '',
				'Translation_initiation_fmet_addition_at_START' : 'FMETTRS',
				'Translation_initiation_gtp_factor_InfB' : 'gtp_hydrolysis',
				'Translation_elongation_FusA_mono' : 'gtp_hydrolysis',
				'Translation_elongation_Tuf_gtp_regeneration' : '',
				'Translation_termination_PrfA_mono_mediated' : '',
				'Translation_termination_PrfB_mono_mediated' : '',
				'Translation_termination_generic_RF_mediated' : '',
				'Translation_termination_peptide_deformylase_processing' : 'DEF',
				'Translation_termination_peptide_chain_release' : 'gtp_hydrolysis',
				'Translation_termination_ribosome_recycler' : '',
				'Protein_processing_GroEL_dependent_folding' : 'atp_hydrolysis_groel',
				'Protein_processing_DnaK_dependent_folding' : 'atp_hydrolysis',
				'Protein_processing_N_terminal_methionine_cleavage' : 'MAP',
				'Ribosome_RbfA_mono_assembly_factor_phase1' : '',
				'Ribosome_RimM_mono_assembly_factor_phase1' : '',
				'Ribosome_gtp_bound_30S_assembly_factor_phase1' : 'gtp_hydrolysis_era'
				},

			'complex_cofactors' : {
				'fes_transfers' : [],
				'biotin_subreactions' : { 'mod_btn_c' : [ 'biotin_ligase' ] },
				'lipoate_subreactions' : { 'mod_lipoyl_c' : [ 'lipoyl_denovo', 'lipoyl_scavenging' ] },
				'fes_chaperones' : {},
				'bmocogdp_chaperones' : {},
				'FeFe/NiFe' : { 'mod_FeFe_cofactor_c' : '', 'mod_NiFe_cofactor_c' : '' }
				},

			'peptide_processing_subreactions' : [
				'Translation_termination_peptide_deformylase_processing',
				'Translation_termination_peptide_chain_release',
				'Translation_termination_ribosome_recycler'
				],

			'translocation_pathway' : {
				'sec' : {
					'abbrev' : 's',
					'keff' : 4.0000,
					'length_dependent_energy' : True,
					'stoichiometry' : 'atp_hydrolysis_sec_pathway'
					},
				'secA' : {
					'abbrev' : 'a',
					'keff' : 4.0000,
					'length_dependent_energy' : True,
					'stoichiometry' : 'atp_hydrolysis_secA'
					},
				'tat' : {
					'abbrev' : 't',
					'keff' : 0.0125,
					'length_dependent_energy' : False,
					'stoichiometry' : ''
					},
				'tat_alt' : {
					'abbrev' : 't',
					'keff' : 0.0125,
					'length_dependent_energy' : False,
					'stoichiometry' : ''
					},
				'yidC' : {
					'abbrev' : 'y',
					'keff' : 20.000,
					'length_dependent_energy' : False,
					'stoichiometry' : 'gtp_hydrolysis'
					},
				'srp' : {
					'abbrev' : 'r',
					'keff' : 20.000,
					'length_dependent_energy' : False,
					'stoichiometry' : 'gtp_hydrolysis_srp_pathway',
					'FtsY' : 'FtsY_MONOMER'
					},
				'srp_yidC' : {
					'abbrev' : 'p',
					'keff' : 20.000,
					'length_dependent_energy' : False,
					'stoichiometry' : 'gtp_hydrolysis'
					},
				'lol' : {
					'abbrev' : 'l',
					'keff' : 0.9000,
					'length_dependent_energy' : False,
					'stoichiometry' : 'atp_hydrolysis'
					},
				'bam' : {
					'abbrev' : 'b',
					'keff' : 0.0270,
					'length_dependent_energy' : False,
					'stoichiometry' : ''
					}
				},

			'excision_machinery' : [
				'rRNA_containing',
				'monocistronic',
				'polycistronic_wout_rRNA'
				],

			'biomass_constraints' : [
				'protein_biomass',
				'mRNA_biomass',
				'tRNA_biomass',
				'rRNA_biomass',
				'ncRNA_biomass',
				'tmRNA_biomass',
				'DNA_biomass',
				'lipid_biomass',
				'constituent_biomass',
				'prosthetic_group_biomass',
				'peptidoglycan_biomass'
				],

			'compartments' : {
				'c' : 'Cytoplasm',
				'e' : 'Extracellular',
				'p' : 'Periplasm',
				'mc': 'ME-model Constraint'
				},

			'START_tRNA' : [],
			'rna_components' : [],
			'knockouts' : [],
			'genome_mods' : {},
			'trna_misacylation' : {},
			'trna_to_codon' : {},
			'trna_to_aa' : {},

			'gam' : 34.98,
			'ngam' : 1.,
			'unmodeled_protein_fraction' : 0.36,

			'braun\'s_lipoproteins' : [],
			'braun\'s_lipid_mod' : 'murein5px4p_p',
			'braun\'s_lpp_flux' : -0.0,
			'braun\'s_murein_flux' : -0.0,

			# active biomass reaction, default value
			'active_biomass_reaction' : 'biomass_constituent_demand'
			}

		# instantiate model parameters as symbols
		# check me.default_parameters and me.symbols
		self.parameters = coralme.core.parameters.MEParameters(self)

		# Create basic M-model structures
		self.reactions = cobra.core.dictlist.DictList()
		self.metabolites = cobra.core.dictlist.DictList()
		self.process_data = cobra.core.dictlist.DictList()
		self._all_genes = cobra.core.dictlist.DictList()

		self._compartments = {}
		self._contexts = []

		# Create the biomass dilution constraint
		self._biomass = coralme.core.component.Constraint('biomass')
		self._biomass_dilution = coralme.core.reaction.SummaryVariable('biomass_dilution')
		self._biomass_dilution.add_metabolites({self._biomass: -1})
		self.add_reactions([self._biomass_dilution])

		# cobra/core/reaction.py:328 Cannot convert expression to float
		# Solved: Check if variable type is sympy.core.symbol.Symbol or float
		# Solved: Removed _populate_solver from reactions -> no need to modify optlang
		self._biomass_dilution.upper_bound = self.mu
		self._biomass_dilution.lower_bound = self.mu

		# Maintenance energy
		self._gam = self.global_info['gam'] # default/user value
		self._ngam = self.global_info['ngam'] # default/user value

		"""
		Unmodeled protein is handled by converting protein_biomass to
		biomass, and requiring production of the appropriate amount of dummy
		protein
		"""
		self._unmodeled_protein_fraction = self.global_info['unmodeled_protein_fraction'] # default/user value

		# troubleshooting flags
		self.troubleshooted = False
		self.troubleshooting = False

		# merging flags
		self.merged_models = {}

		# aliases
		self._aliases = { 'reactions' : {}, 'metabolites' : {} }

	def __getstate__(self):
		state = self.__dict__.copy()

		# Don't pickle unit_registry
		del state["unit_registry"]
		# Don't pickle optimization methods
		del state["optimize"]
		del state["optimize_windows"]
		del state["feasibility"]
		del state["feas_gurobi"]
		del state["feas_cplex"]
		del state["construct_lp_problem"]
		del state["fva"]
		# Don't pickle troubleshooting methods
		if hasattr(self, 'get_solution'):
			del state['get_solution']
			del state['check_feasibility']
		if hasattr(self, 'get_feasibility'):
			del state['get_feasibility']
		return state

	def __setstate__(self, state):
		self.__dict__.update(state)
		# Add unit_registry back since it doesn't exist in the pickle
		self.unit_registry = self.mu._REGISTRY

		# simulation methods in optimization.py
		self.optimize = types.MethodType(coralme.core.optimization.optimize, self)
		self.optimize_windows = types.MethodType(coralme.core.optimization.optimize_windows, self)
		self.feasibility = types.MethodType(coralme.core.optimization.feasibility, self) # qminos
		self.feas_gurobi = types.MethodType(coralme.core.optimization.feas_gurobi, self) # gurobi
		self.feas_cplex = types.MethodType(coralme.core.optimization.feas_cplex, self) # cplex
		self.construct_lp_problem = types.MethodType(coralme.core.optimization.construct_lp_problem, self)
		self.fva = types.MethodType(coralme.core.optimization.fva, self)

	@property
	def active_biomass_reaction(self):
		return self.get(self.global_info['active_biomass_reaction'])

	@active_biomass_reaction.setter
	def active_biomass_reaction(self, name):
		if self.global_info['biomass_reactions'] != ['biomass_constituent_demand']:
			name = 'biomass_constituent_demand_' + name
			biomass_reactions = [ 'biomass_constituent_demand_{:s}'.format(x) for x in self.biomass_reactions ]
		else:
			name = 'biomass_constituent_demand'
			biomass_reactions = self.global_info['biomass_reactions']

		# close reaction bounds
		for rxn in biomass_reactions:
			self.reactions.get_by_id(rxn).bounds = (0., 0.) # close bounds for every biomass reaction

		if 'lipid_demand_per_condition' in self.global_info:
			for cond, rxns in self.global_info['lipid_demand_per_condition'].items():
				for rxn in rxns:
					self.reactions.get_by_id(rxn).bounds = (0., 0.) # close bounds for every lipid composition reaction

		if not self.reactions.has_id(name):
			raise ValueError('ME-model has no biomass reaction \'biomass_constituent_demand_{:s}\''.format(name))

		# open bounds for active biomass reaction
		self.reactions.get_by_id(name).bounds = (self.mu, self.mu)
		if 'lipid_demand_per_condition' in self.global_info:
			for rxn in self.global_info['lipid_demand_per_condition'][name.replace('biomass_constituent_demand_', '')]:
				self.reactions.get_by_id(rxn).bounds = (self.mu, self.mu)
		self.global_info['active_biomass_reaction'] = name

	@property
	def aliases(self):
		return self._aliases

	@aliases.setter
	def aliases(self, args):
		if args == {}:
			self._aliases = { 'reactions' : {}, 'metabolites' : {} }
		else:
			self._aliases['metabolites'].update(args.get('metabolites', {}))
			self._aliases['reactions'].update(args.get('reactions', {}))
			# add new aliases from metabolite aliases
			for key, value in self._aliases['metabolites'].items():
				for reaction in self.reactions.query(value.replace('(', r'\(').replace(')', r'\)').replace('[', r'\[').replace(']', r'\]')):
					self._aliases['reactions'][reaction.id.replace(value, key)] = reaction.id

	def perform_gene_knockouts(self, genes):
		return coralme.util.essentiality.perform_gene_knockouts(self, genes)

	def to_json(self, outfile):
		coralme.io.json.save_json_me_model(self, outfile)

	def to_pickle(self, outfile):
		coralme.io.pickle.save_pickle_me_model(self, outfile)

	def minimize(self, id_or_model = 'copy', name = 'copy', include_original_m_model = False, include_processed_m_model = False, include_processdata = True):
		new_model = coralme.core.model.MEModel(id_or_model = id_or_model, name = name)
		# add_processdata, add_metabolites, and add_reactions take care of
		# new memory addresses for associated data
		new_model.add_processdata([ x.copy() for x in self.process_data ])
		new_model.global_info = copy.deepcopy(self.global_info)
		new_model.metabolites[0].remove_from_model()
		new_model.add_metabolites([ x.copy() for x in self.metabolites ])
		new_model.reactions[0].remove_from_model()
		# reaction copies should be associated to new process data
		# the copy includes the objective coefficient and process data
		new_model.add_reactions([ x.copy() for x in self.reactions ])
		new_model.compartments = self.compartments
		if not include_processdata:
			del new_model.process_data
		if include_original_m_model:
			new_model.gem = self.gem
		if include_processed_m_model:
			new_model.processed_m_model = self.processed_m_model
		return new_model

	@staticmethod
	def from_cobra(model, objective = None):
		if model.notes.get('from cobra', False):
			return model

		def reaction_from_cobra_model(model, reaction):
			new_reaction = MEReaction(reaction.id)
			new_reaction.name = reaction.name
			new_reaction.subsystem = reaction.subsystem
			new_reaction.lower_bound = reaction.lower_bound
			new_reaction.upper_bound = reaction.upper_bound
			new_reaction.gpr = reaction.gpr
			for met, stoichiometry in reaction.metabolites.items():
				new_reaction.add_metabolites({ model.metabolites.get_by_id(met.id): stoichiometry })
			new_reaction.cofactors = reaction.cofactors if hasattr(reaction, 'cofactors') else cobra.core.GPR.from_string('')
			return new_reaction

		def metabolite_from_cobra_model(model, metabolite):
			new_metabolite = Metabolite(metabolite.id)
			new_metabolite.name = metabolite.name
			new_metabolite.formula = metabolite.formula
			new_metabolite.compartment = metabolite.compartment
			new_metabolite.charge = metabolite.charge
			new_metabolite.annotation = metabolite.annotation
			new_metabolite.notes = metabolite.notes
			new_metabolite.functional = True
			return new_metabolite

		new_model = MEModel()
		new_model.metabolites[0].remove_from_model()
		new_model.add_metabolites([ metabolite_from_cobra_model(model, x) for x in model.metabolites ])
		new_model.reactions[0].remove_from_model()
		new_model.add_reactions([ reaction_from_cobra_model(model, x) for x in model.reactions ])
		new_model.all_genes = model.genes

		if objective is not None:
			new_model.reactions.get_by_id(objective).objective_coefficient = +1
		else:
			# bof: defaultdict = { (optlang.gurobi_interface.Variable, coeff) }
			bof = model.objective.expression.as_coefficients_dict()
			for variable, objective_coefficient in bof.items():
				if 'reverse' in variable.name:
					continue
				new_model.reactions.get_by_id(variable.name).objective_coefficient = objective_coefficient
		new_model.gem = copy.deepcopy(model)
		new_model.notes = {
			'from cobra' : True
			}

		# simulation methods in optimization.py
		new_model.optimize = types.MethodType(coralme.core.optimization.optimize, new_model)
		new_model.feasibility = types.MethodType(coralme.core.optimization.feasibility, new_model)
		new_model.construct_lp_problem = types.MethodType(coralme.core.optimization.construct_lp_problem, new_model)
		new_model.fva = types.MethodType(coralme.core.optimization.fva, new_model)

		return new_model

	@property
	def default_parameters(self):
		return self.global_info.get('default_parameters', {})

	@default_parameters.setter
	def default_parameters(self, args):
		"""
		This will only update original MEModel symbols.

		Using an empty dictionary will reset the values.

		Use `me.global_info['default_parameters'].update` to add new symbols and values.

		Use 'kt' instead of 'k_t'
		Use 'r0' instead of 'r_0'
		Use 'k_deg' instead of 'k^mRNA_deg'
		Use 'kcat' instead of 'k^default_cat'
		"""
		self.global_info['default_parameters'].update({
			sympy.Symbol('k_t', positive = True) : args.get('kt', 4.5),
			sympy.Symbol('r_0', positive = True) : args.get('r0', 0.087),
			sympy.Symbol('k^mRNA_deg', positive = True) : args.get('k_deg', 12.0),
			sympy.Symbol('m_rr', positive = True) : args.get('m_rr', 1453.0),
			sympy.Symbol('m_aa', positive = True) : args.get('m_aa', 0.109),
			sympy.Symbol('m_nt', positive = True) : args.get('m_nt', 0.324),
			sympy.Symbol('f_rRNA', positive = True) : args.get('f_rRNA', 0.86),
			sympy.Symbol('f_mRNA', positive = True) : args.get('f_mRNA', 0.02),
			sympy.Symbol('f_tRNA', positive = True) : args.get('f_tRNA', 0.12),
			sympy.Symbol('m_tRNA', positive = True) : args.get('m_tRNA', 25.0),
			sympy.Symbol('k^default_cat', positive = True) : args.get('kcat', 65.0), # not stored in json with coralME v1.0
			sympy.Symbol('temperature', positive = True) : args.get('temperature', 37.0),
			sympy.Symbol('propensity_scaling', positive = True) : args.get('propensity_scaling', 0.45),
			# DNA replication; see dna_replication.percent_dna_template_function
			sympy.Symbol('g_p_gdw_0', positive = True) : args.get('g_p_gdw_0', 0.059314110730022594), # dimensionless
			sympy.Symbol('g_per_gdw_inf', positive = True) : args.get('g_per_gdw_inf', 0.02087208296776481), # dimensionless
			sympy.Symbol('b', positive = True) : args.get('b', 0.1168587392731988), # per hour**d
			sympy.Symbol('d', positive = True) : args.get('c', 3.903641432780327) # dimensionless
			})

	@property
	def mu(self):
		return self._mu

	@mu.setter
	def mu(self, value: str):
		# set growth rate symbolic variable
		self._mu_old = self._mu
		self._mu = sympy.Symbol(value, positive = True) * self.unit_registry.parse_units('1 per hour')

		if self._mu_old == self._mu:
			return # doing nothing because user changed to the current mu

		for rxn in self.reactions:
			if hasattr(rxn.lower_bound, 'subs'):
				rxn.lower_bound = rxn.lower_bound.magnitude.subs({ self._mu_old.magnitude : self._mu.magnitude }) * self.unit_registry.parse_units('1 per hour')
			if hasattr(rxn.upper_bound, 'subs'):
				rxn.upper_bound = rxn.upper_bound.magnitude.subs({ self._mu_old.magnitude : self._mu.magnitude }) * self.unit_registry.parse_units('1 per hour')
			for met, coeff in rxn.metabolites.items():
				if hasattr(coeff, 'subs'):
					rxn._metabolites[met] = coeff.subs({ self._mu_old.magnitude : self._mu.magnitude }) * self.unit_registry.parse_units('dimensionless')

		for symbol, fn in self.symbols.items():
			if hasattr(fn, 'units'):
				if str(fn.units) == 'dimensionless':
					fn.magnitude.subs({ self._mu_old.magnitude : self._mu.magnitude }) * self.unit_registry.parse_units('dimensionless')
				else:
					self.symbols[symbol] = fn.magnitude.subs({ self._mu_old.magnitude : self._mu.magnitude }) * self.unit_registry.parse_units(str(fn.units))
			else:
				self.symbols[symbol] = fn.subs({ self._mu_old.magnitude : self._mu.magnitude })

	# WARNING: FROM COBRAPY WITHOUT MODIFICATIONS
	@property
	def compartments(self) -> typing.Dict:
		"""Return all metabolites' compartments.

		Returns
		-------
		dict
			A dictionary of metabolite compartments, where the keys are the short
			version (one letter version) of the compartmetns, and the values are the
			full names (if they exist).
		"""
		return {
			met.compartment: self._compartments.get(met.compartment, "")
			for met in self.metabolites
			if met.compartment is not None
		}

	# WARNING: FROM COBRAPY WITHOUT MODIFICATIONS
	@compartments.setter
	def compartments(self, value: typing.Dict) -> None:
		"""Get or set the dictionary of current compartment descriptions.

		Assigning a dictionary to this property updates the model's
		dictionary of compartment descriptions with the new values.

		Parameters
		----------
		value : dict
			Dictionary mapping compartments abbreviations to full names.

		Examples
		--------
		>>> from cobra.io import load_model
		>>> model = load_model("textbook")
		>>> model.compartments = {'c': 'the cytosol'}
		>>> model.compartments
		{'c': 'the cytosol', 'e': 'extracellular'}
		"""
		self._compartments.update(value)

	# WARNING: FROM COBRAPY WITHOUT MODIFICATIONS
	@property
	def medium(self):
		"""Get the constraints on the model exchanges.

		`model.medium` returns a dictionary of the bounds for each of the
		boundary reactions, in the form of `{rxn_id: bound}`, where `bound`
		specifies the absolute value of the bound in direction of metabolite
		creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

		Returns
		-------
		Dict[str, float]
			A dictionary with rxn.id (str) as key, bound (float) as value.
		"""

		def is_active(reaction) -> bool:
			"""Determine if boundary reaction permits flux towards creating metabolites.

			Parameters
			----------
			reaction: cobra.Reaction

			Returns
			-------
			bool
				True if reaction produces metaoblites and has upper_bound above 0
				or if reaction consumes metabolites and has lower_bound below 0 (so
				could be reversed).
			"""
			return (bool(reaction.products) and (reaction.upper_bound > 0)) or (
				bool(reaction.reactants) and (reaction.lower_bound < 0)
			)

		def get_active_bound(reaction) -> float:
			"""For an active boundary reaction, return the relevant bound.

			Parameters
			----------
			reaction: cobra.Reaction

			Returns
			-------
			float:
				upper or minus lower bound, depenending if the reaction produces or
				consumes metaoblties.
			"""
			if reaction.reactants:
				return -reaction.lower_bound
			elif reaction.products:
				return reaction.upper_bound

		return {
			rxn.id: get_active_bound(rxn) for rxn in self.get_exchange_reactions if is_active(rxn)
		}

	# WARNING: FROM COBRAPY WITHOUT MODIFICATIONS
	@medium.setter
	def medium(self, medium) -> None:
		"""Set the constraints on the model exchanges.

		`model.medium` returns a dictionary of the bounds for each of the
		boundary reactions, in the form of `{rxn_id: rxn_bound}`, where `rxn_bound`
		specifies the absolute value of the bound in direction of metabolite
		creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

		Parameters
		----------
		medium: dict
			The medium to initialize. medium should be a dictionary defining
			`{rxn_id: bound}` pairs.
		"""

		def set_active_bound(reaction, bound: float) -> None:
			"""Set active bound.

			Parameters
			----------
			reaction: cobra.Reaction
				Reaction to set
			bound: float
				Value to set bound to. The bound is reversed and set as lower bound
				if reaction has reactants (metabolites that are consumed). If reaction
				has reactants, it seems the upper bound won't be set.
			"""
			if reaction.reactants:
				reaction.lower_bound = -bound
			elif reaction.products:
				reaction.upper_bound = bound

		# Set the given media bounds
		media_rxns = []
		exchange_rxns = frozenset(self.get_exchange_reactions)
		for rxn_id, rxn_bound in medium.items():
			rxn = self.reactions.get_by_id(rxn_id)
			if rxn not in exchange_rxns:
				logger.warning(
					f"{rxn.id} does not seem to be an an exchange reaction. "
					f"Applying bounds anyway."
				)
			media_rxns.append(rxn)
			# noinspection PyTypeChecker
			set_active_bound(rxn, rxn_bound)

		frozen_media_rxns = frozenset(media_rxns)

		# Turn off reactions not present in media
		for rxn in exchange_rxns - frozen_media_rxns:
			is_export = rxn.reactants and not rxn.products
			set_active_bound(
				rxn, min(0.0, -rxn.lower_bound if is_export else rxn.upper_bound)
			)

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def copy(self):
		return copy.deepcopy(self)

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def prune_unused_metabolites(self):
		# originally at cobra.manipulation.delete.prune_unused_metabolites, but it requires to make a copy of the model
		inactive_metabolites = [ m for m in self.metabolites if len(m.reactions) == 0 ]
		self.remove_metabolites(inactive_metabolites)
		return inactive_metabolites

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def prune_unused_reactions(self):
		# originally at cobra.manipulation.delete.prune_unused_reactions, but it requires to make a copy of the model
		reactions_to_prune = [ r for r in self.reactions if len(r.metabolites) == 0 ]
		self.remove_reactions(reactions_to_prune)
		return reactions_to_prune

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	# WARNING: NEW IMPLEMENTATION AND VERY EXPERIMENTAL
	@staticmethod
	def merge(models_to_merge = {}, id_or_model = 'merge', name = 'merge'):
		# check if models' mu values are different
		mus = [ v.mu for v in models_to_merge.values() ]
		if not len(mus) == len(set(mus)):
			raise ValueError('')

		# create an empty coralME model, and copy the merging models into merged_models dictionary
		merge = coralme.core.model.MEModel(id_or_model = id_or_model, name = name)
		merge.reactions[0].remove_from_model()
		merge.metabolites[0].remove_from_model()

		for org, me in list(models_to_merge.items()):
			merge.merged_models[org] = me.copy()
			merge.merged_models[org].merging_key = org

			# add tags to merging models to make them unique in the new model
			for data in merge.merged_models[org].process_data:
				data.id = '{:s}_{:s}'.format(org, data.id)
				if hasattr(data, '_stoichiometry'):
					for met in list(data._stoichiometry.keys()):
						data._stoichiometry['{:s}_{:s}'.format(org, met)] = data._stoichiometry.pop(met)

			for data in merge.merged_models[org].metabolites:
				if not data.id.endswith('_e'): # do not modify medium
					data._id = '{:s}_{:s}'.format(org, data.id)

			for data in merge.merged_models[org].reactions:
				if not data.id.startswith('EX_'): # do not modify medium
					data._id = '{:s}_{:s}'.format(org, data.id)

		# add renamed process_data, metabolites, and reactions to merge model
		for org, me in list(merge.merged_models.items()):
			merge.add_processdata(merge.merged_models[org].process_data)
			merge.add_metabolites(merge.merged_models[org].metabolites)
			merge.add_reactions(merge.merged_models[org].reactions)

		return merge

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	@property
	def objective(self):
		# TODO: make it look like cobrapy output?
		return [ (x, x.objective_coefficient) for x in self.reactions if x.objective_coefficient != 0 ]

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	@objective.setter
	def objective(self, dct = { 'dummy_reaction_FWD_SPONT' : +1. } ):
		for rxn in self.reactions:
			rxn.objective_coefficient = 0.

		for rxn, coeff in dct.items():
			if self.reactions.has_id(rxn):
				self.reactions.get_by_id(rxn).objective_coefficient = coeff
			else:
				raise ValueError('Reaction \'{:s}\' does not exist in the ME-model'.format(rxn))

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def add_metabolites(self, metabolite_list):
		"""Add new metabolites to a model.

		Will add a list of metabolites to the model object.

		This function is different from COBRApy and it won't:
			Add new constraints accordingly.
			Revert changes upon exit when using the model as a context.

		Parameters
		----------
		metabolite_list : list or Metabolite.
			A list of `cobra.core.Metabolite` objects. If it isn't an iterable
			container, the metabolite will be placed into a list.

		"""
		if not hasattr(metabolite_list, "__iter__"):
			metabolite_list = [metabolite_list]
		if len(metabolite_list) == 0:
			return None

		# First check whether the metabolites exist in the model
		metabolite_list = [x for x in metabolite_list if x.id not in self.metabolites]

		bad_ids = [
			m for m in metabolite_list if not isinstance(m.id, str) or len(m.id) < 1
		]
		if len(bad_ids) != 0:
			raise ValueError("invalid identifiers in {}".format(repr(bad_ids)))

		for x in metabolite_list:
			x._model = self
		self.metabolites += metabolite_list

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def remove_metabolites(self, metabolite_list, destructive=False):
		"""Remove a list of metabolites from the the object.

		This function is different from COBRApy and it won't:
			Revert changes upon exit when using the model as a context.

		Parameters
		----------
		metabolite_list : list or Metaoblite
			A list of `cobra.core.Metabolite` objects. If it isn't an iterable
			container, the metabolite will be placed into a list.

		destructive : bool, optional
			If False then the metabolite is removed from all
			associated reactions.  If True then all associated
			reactions are removed from the Model (default False).
		"""
		if not hasattr(metabolite_list, "__iter__"):
			metabolite_list = [metabolite_list]
		# Make sure metabolites exist in model
		metabolite_list = [x for x in metabolite_list if x.id in self.metabolites]
		for x in metabolite_list:
			x._model = None

			# remove reference to the metabolite in all groups
			#associated_groups = self.get_associated_groups(x)
			#for group in associated_groups:
				#group.remove_members(x)

			if not destructive:
				for the_reaction in list(x._reaction):  # noqa W0212
					the_coefficient = the_reaction._metabolites[x]  # noqa W0212
					the_reaction.subtract_metabolites({x: the_coefficient})

			else:
				for x2 in list(x._reaction):  # noqa W0212
					x2.remove_from_model()

		self.metabolites -= metabolite_list

	# WARNING: New method based on add_reactions
	def add_processdata(self, processdata_list):
		def existing_filter(data):
			if data.id in self.process_data:
				return False
			return True

		# First check whether the reactions exist in the model.
		pruned = cobra.core.dictlist.DictList(filter(existing_filter, processdata_list))

		#
		for data in pruned:
			data._model = self

		self.process_data += pruned

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def add_reactions(self, reaction_list):
		"""Add reactions to the model.

		Reactions with identifiers identical to a reaction already in the
		model are ignored.

		This function is different from COBRApy and it won't:
			Revert changes upon exit when using the model as a context.

		Parameters
		----------
		reaction_list : list
			A list of `cobra.Reaction` objects
		"""
		def existing_filter(rxn):
			if rxn.id in self.reactions:
				return False
			return True

		# First check whether the reactions exist in the model.
		pruned = cobra.core.dictlist.DictList(filter(existing_filter, reaction_list))

		# Add reactions. Also take care of genes and metabolites in the loop.
		for reaction in pruned:
			reaction._model = self

			# WARNING: DO NOT DELETE
			# Build a `list()` because the dict will be modified in the loop.
			for metabolite in list(reaction.metabolites):
				# TODO: Maybe this can happen with
				#  Reaction.add_metabolites(combine=False)
				# TODO: Should we add a copy of the metabolite instead?
				if metabolite not in self.metabolites:
					self.add_metabolites(metabolite)
				# A copy of the metabolite exists in the model, the reaction
				# needs to point to the metabolite in the model.
				else:
					# FIXME: Modifying 'private' attributes is horrible.
					stoichiometry = reaction._metabolites.pop(metabolite)
					model_metabolite = self.metabolites.get_by_id(metabolite.id)
					reaction._metabolites[model_metabolite] = stoichiometry
					model_metabolite._reaction.add(reaction)

			# WARNING: coralme reactions can have process_data associated to them
			if hasattr(reaction, 'process_data'):
				for key, value in reaction.process_data.items():
					if value is None:
						setattr(reaction, key, value)
					else:
						setattr(reaction, key, self.process_data.get_by_id(value.id))
				delattr(reaction, 'process_data')

			# WARNING: units system is associated to the model
			if isinstance(reaction.lower_bound, (numpy.floating, float, numpy.integer, int, sympy.Symbol)):
				reaction.lower_bound = reaction.lower_bound * reaction._model.unit_registry.parse_units('mmols per gram per hour')
			else:
				reaction.lower_bound = reaction.lower_bound.magnitude * reaction._model.unit_registry.parse_units('mmols per gram per hour')

			if isinstance(reaction.upper_bound, (numpy.floating, float, numpy.integer, int, sympy.Symbol)):
				reaction.upper_bound = reaction.upper_bound * reaction._model.unit_registry.parse_units('mmols per gram per hour')
			else:
				reaction.upper_bound = reaction.upper_bound.magnitude * reaction._model.unit_registry.parse_units('mmols per gram per hour')

		self.reactions += pruned

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def remove_reactions(self, reactions, remove_orphans=False):
		"""Remove reactions from the model.

		This function is different from COBRApy and it won't:
			Revert changes upon exit when using the model as a context.
			Remove orphaned genes

		Parameters
		----------
		reactions : list or reaction or str
			A list with reactions (`cobra.Reaction`), or their id's, to remove.
			Reaction will be placed in a list. Str will be placed in a list and used to
			find the reaction in the model.
		remove_orphans : bool, optional
			Remove orphaned genes and metabolites from the model as
			well (default False).
		"""
		if isinstance(reactions, str) or hasattr(reactions, "id"):
			reactions = [reactions]

		for reaction in reactions:
			# Make sure the reaction is in the model
			try:
				reaction = self.reactions[self.reactions.index(reaction)]
			except ValueError:
				warn(f"{reaction} not in {self}")
			else:
				self.reactions.remove(reaction)
				reaction._model = None

				for met in reaction._metabolites:
					if reaction in met._reaction:
						met._reaction.remove(reaction)
						if remove_orphans and len(met._reaction) == 0:
							self.remove_metabolites(met)

				#for gene in reaction._genes:
					#if reaction in gene._reaction:
						#gene._reaction.remove(reaction)
						#if remove_orphans and len(gene._reaction) == 0:
							#self.genes.remove(gene)

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def add_boundary(
		self,
		metabolite: Metabolite,
		type: str = "exchange",
		reaction_id: typing.Optional[str] = None,
		lb: typing.Optional[float] = None,
		ub: typing.Optional[float] = None,
		sbo_term: typing.Optional[str] = None,
	) -> MEReaction:
		"""
		Add a boundary reaction for a given metabolite.

		There are three different types of pre-defined boundary reactions:
		exchange, demand, and sink reactions.
		An exchange reaction is a reversible, unbalanced reaction that adds
		to or removes an extracellular metabolite from the extracellular
		compartment.
		A demand reaction is an irreversible reaction that consumes an
		intracellular metabolite.
		A sink is similar to an exchange but specifically for intracellular
		metabolites, i.e., a reversible reaction that adds or removes an
		intracellular metabolite.

		If you set the reaction `type` to something else, you must specify the
		desired identifier of the created reaction along with its upper and
		lower bound. The name will be given by the metabolite name and the
		given `type`.

		The change is reverted upon exit when using the model as a context.

		Parameters
		----------
		metabolite : cobra.Metabolite
			Any given metabolite. The compartment is not checked but you are
			encouraged to stick to the definition of exchanges and sinks.
		type : {"exchange", "demand", "sink"}
			Using one of the pre-defined reaction types is easiest. If you
			want to create your own kind of boundary reaction choose
			any other string, e.g., 'my-boundary' (default "exchange").
		reaction_id : str, optional
			The ID of the resulting reaction. This takes precedence over the
			auto-generated identifiers but beware that it might make boundary
			reactions harder to identify afterwards when using `model.boundary`
			or specifically `model.exchanges` etc. (default None).
		lb : float, optional
			The lower bound of the resulting reaction (default None).
		ub : float, optional
			The upper bound of the resulting reaction (default None).
		sbo_term : str, optional
			A correct SBO term is set for the available types. If a custom
			type is chosen, a suitable SBO term should also be set (default None).

		Returns
		-------
		cobra.Reaction
			The created boundary reaction.

		Examples
		--------
		>>> from cobra.io load_model
		>>> model = load_model("textbook")
		>>> demand = model.add_boundary(model.metabolites.atp_c, type="demand")
		>>> demand.id
		'DM_atp_c'
		>>> demand.name
		'ATP demand'
		>>> demand.bounds
		(0, 1000.0)
		>>> demand.build_reaction_string()
		'atp_c --> '

		"""
		ub = cobra.Configuration().upper_bound if ub is None else ub
		lb = cobra.Configuration().lower_bound if lb is None else lb
		types = {
			"exchange": ("EX", lb, ub, cobra.medium.sbo_terms["exchange"]),
			"demand": ("DM", 0, ub, cobra.medium.sbo_terms["demand"]),
			"sink": ("SK", lb, ub, cobra.medium.sbo_terms["sink"]),
		}
		if type == "exchange":
			external = cobra.medium.find_external_compartment(self)
			if metabolite.compartment != external:
				raise ValueError(
					f"The metabolite is not an external metabolite (compartment is "
					f"`{metabolite.compartment}` but should be `{external}`). "
					f"Did you mean to add a demand or sink? If not, either change"
					f" its compartment or rename the model compartments to fix this."
				)
		if type in types:
			prefix, lb, ub, default_term = types[type]
			if reaction_id is None:
				reaction_id = f"{prefix}_{metabolite.id}"
			if sbo_term is None:
				sbo_term = default_term
		if reaction_id is None:
			raise ValueError(
				"Custom types of boundary reactions require a custom "
				"identifier. Please set the `reaction_id`."
			)
		if reaction_id in self.reactions:
			return None
			#raise ValueError(f"Boundary reaction '{reaction_id}' already exists.")
		name = f"{metabolite.name} {type}"
		rxn = MEReaction(id=reaction_id, name=name)
		# WARNING: setting lb and ub through MEReaction definition is not working
		rxn.lower_bound = lb
		rxn.upper_bound = ub
		rxn.add_metabolites({metabolite: -1})
		if sbo_term:
			rxn.annotation["sbo"] = sbo_term
		self.add_reactions([rxn])
		return rxn

	# WARNING: Modified functions from COBRAme and new functions
	@property
	def get_exchange_reactions(self):
		return self.reactions.query('^EX_')

	@property
	def get_sink_reactions(self):
		return self.reactions.query('^SK_')

	@property
	def get_demand_reactions(self):
		return self.reactions.query('^DM_')

	@property
	def get_troubleshooted_reactions(self):
		return self.reactions.query('^TS_')

	def remove_troubleshooted_reactions(self):
		return self.remove_reactions(self.get_troubleshooted_reactions)

	@property
	def get_unbounded_reactions(self):
		return [ x for x in self.reactions if x.bound_violation[0] ]

	@property
	def get_spontaneous_reactions(self):
		return self.reactions.query('_FWD_SPONT$|_REV_SPONT$')

	@property
	def get_null_gpr_metabolic_reactions(self):
		# TODO: remove false positive reactions (aka, reactions with an enzyme that also use CPLX_dummy)
		return [ x for x in self.get('CPLX_dummy').reactions ]

	@property
	def get_mass_unbalanced_reactions(self):
		return [ x for x in self.reactions if isinstance(x.get_me_mass_balance(), dict) and x.get_me_mass_balance() != {} ]

	def add_biomass_constraints_to_model(self, biomass_types):
		for biomass_type in tqdm.tqdm(biomass_types, 'Adding biomass constraint(s) into the ME-model...', bar_format = bar_format):
			if '_biomass' not in biomass_type:
				raise ValueError('Biomass types should be suffixed with \'_biomass\'.')
			constraint_obj = coralme.core.component.Constraint(biomass_type)
			summary_variable_obj = coralme.core.reaction.SummaryVariable('{:s}_to_biomass'.format(biomass_type))
			summary_variable_obj.add_metabolites({constraint_obj: -1, self._biomass: 1})
			self.add_reactions([summary_variable_obj])

	@property
	def get_unmodeled_protein(self):
		return self.metabolites.get_by_id('protein_dummy')

	@property
	def get_unmodeled_protein_biomass(self):
		return self.metabolites.get_by_id('unmodeled_protein_biomass')

	@property
	def unmodeled_protein_fraction(self):
		return self._unmodeled_protein_fraction

	@unmodeled_protein_fraction.setter
	def unmodeled_protein_fraction(self, value):
		if 'protein_biomass_to_biomass' not in self.reactions:
			raise UserWarning(
				'Must add SummaryVariable handling the protein '
				'biomass constraint (via :meth:`add_biomass_constraints_to_model`) '
				'before defining the unmodeled protein fraction'
				)

		# See the Biomass_formulations for an explanation (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006302)
		if 0 <= value < 1.:
			amount = value / (1 - value)
		else:
			raise ValueError('The unmodeled protein fraction cannot be exactly 1 or greater.')

		self.reactions.protein_biomass_to_biomass.add_metabolites({self.get_unmodeled_protein_biomass: -amount}, combine = False)
		self.reactions.protein_biomass_to_biomass.add_metabolites({self._biomass: 1 + amount}, combine = False)
		self._unmodeled_protein_fraction = value

	@property
	def gam(self):
		return self._gam

	@gam.setter
	def gam(self, value):
		if 'GAM' not in self.reactions:
			logging.warning('Adding GAM (ATP requirement for growth) reaction into the ME-model.')
			self.add_reactions([coralme.core.reaction.SummaryVariable('GAM')])
			self.reactions.GAM.lower_bound = self.mu
			self.reactions.GAM.upper_bound = 1000.
		#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1, 'pi_c': 1} # charges: -4, 0 => -3, +1, -2
		atp_hydrolysis = self.process_data.get_by_id('atp_hydrolysis').stoichiometry
		for met, coeff in atp_hydrolysis.items():
			self.reactions.GAM.add_metabolites({met: value * coeff}, combine = False)
		self._gam = value

		# check stoichiometry
		if self.reactions.GAM.check_mass_balance().get('H', False):
			tmp = collections.Counter(self.reactions.GAM.metabolites)
			tmp.update({ self.metabolites.h_c : -1*self.reactions.GAM.check_mass_balance()['H'] })
			tmp = { k:v for k,v in tmp.items() if v != 0. }
			self.reactions.GAM._metabolites = dict(tmp)

	@property
	def ngam(self):
		return self._ngam

	@ngam.setter
	def ngam(self, value):
		if 'ATPM' not in self.reactions:
			logging.warning('Adding ATPM (ATP requirement for maintenance) reaction into the ME-model.')
			#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1, 'pi_c': 1} # charges: -4, 0 => -3, +1, -2
			atp_hydrolysis = self.process_data.get_by_id('atp_hydrolysis').stoichiometry
			self.add_reactions([coralme.core.reaction.SummaryVariable('ATPM')])
			self.reactions.ATPM.add_metabolites(atp_hydrolysis)
		self.reactions.ATPM.lower_bound = value
		self.reactions.ATPM.upper_bound = 1000.
		self._ngam = value

		# check stoichiometry
		if self.reactions.ATPM.check_mass_balance().get('H', False):
			tmp = collections.Counter(self.reactions.ATPM.metabolites)
			tmp.update({ self.metabolites.h_c : -1*self.reactions.ATPM.check_mass_balance()['H'] })
			tmp = { k:v for k,v in tmp.items() if v != 0. }
			self.reactions.ATPM._metabolites = dict(tmp)

	def add_translocation_pathway(self, key = 'new', abbrev: str = 'n', keff: float = 65., length_dependent_energy: bool = False, stoichiometry: str = '', enzymes: dict = {}):
		# check properties of enzymes
		for k,v in enzymes.items():
			if 'fixed_keff' in v and 'length_dependent' in v:
				pass

		self.global_info['translocation_pathway'][key] = {
			'abbrev': abbrev,
			'keff': keff,
			'length_dependent_energy': length_dependent_energy,
			'stoichiometry': stoichiometry if self.reactions.has_id(stoichiometry) and stoichiometry != '' else '',
			'enzymes': enzymes
			}

		return self.global_info['translocation_pathway']

	# data types generators:
	# StoichiometricData, ComplexData, TranslationData, TranscriptionData,
	# GenericData, tRNAData, TranslocationData, PostTranslationData, SubreactionData
	@property
	def stoichiometric_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.StoichiometricData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.StoichiometricData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def complex_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.ComplexData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.ComplexData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def translation_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.TranslationData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.TranslationData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def transcription_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.TranscriptionData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.TranscriptionData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def generic_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.GenericData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.GenericData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def tRNA_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.tRNAData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.tRNAData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def translocation_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.TranslocationData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.TranslocationData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def posttranslation_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.PostTranslationData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.PostTranslationData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def subreaction_data(self):
		#for data in self.process_data:
			#if isinstance(data, coralme.core.processdata.SubreactionData):
				#yield data
		lst = [ x for x in self.process_data if isinstance(x, coralme.core.processdata.SubreactionData)]
		return cobra.core.dictlist.DictList(lst)

	@property
	def genes(self):
		return self._all_genes

	@property
	def all_genes(self):
		if len(self._all_genes) == 0.:
			lst = [ g for g in self.metabolites if isinstance(g, coralme.core.component.TranscribedGene) and "dummy" not in g.id ]
			self._all_genes = cobra.core.dictlist.DictList(lst)
		return self._all_genes

	@all_genes.setter
	def all_genes(self, values):
		if self.notes.get('from cobra', False):
			lst = [ g for g in self.metabolites if isinstance(g, coralme.core.component.TranscribedGene) and "dummy" not in g.id ]
			self._all_genes = cobra.core.dictlist.DictList(lst)
		else:
			self._all_genes = values

	@property
	def mRNA_genes(self):
		lst = [ g for g in self.all_genes if hasattr(g, 'RNA_type') and g.RNA_type == 'mRNA' ]
		return cobra.core.dictlist.DictList(lst)

	@property
	def rRNA_genes(self):
		lst = [ g for g in self.all_genes if hasattr(g, 'RNA_type') and g.RNA_type == 'rRNA' ]
		return cobra.core.dictlist.DictList(lst)

	@property
	def tRNA_genes(self):
		lst = [ g for g in self.all_genes if hasattr(g, 'RNA_type') and g.RNA_type == 'tRNA' ]
		return cobra.core.dictlist.DictList(lst)

	@property
	def pseudo_genes(self):
		lst = [ g.mRNA for g in [ g for g in self.translation_data if g.pseudo ] if not g.id.endswith('dummy') ]
		return lst

	def get_metabolic_flux(self, solution = None):
		"""Extract the flux state for Metabolic reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {r.id: 0 for r in tqdm.tqdm(list(self.stoichiometric_data), 'Building reaction dictionary...', bar_format = bar_format)}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.MetabolicReaction):
				m_reaction_id = reaction.stoichiometric_data.id
				if reaction.reverse:
					flux_dict[m_reaction_id] -= solution.fluxes[reaction.id]
				else:
					flux_dict[m_reaction_id] += solution.fluxes[reaction.id]
			# SummaryVariable in M-model
			elif reaction.id == 'ATPM':
				flux_dict[reaction.id] = solution.fluxes[reaction.id]
			# Exchange, Demand, and Sink reactions
			elif reaction.id.startswith('EX_') or reaction.id.startswith('DM_') or reaction.id.startswith('SK_'):
				flux_dict[reaction.id] = solution.fluxes[reaction.id]
		return flux_dict

	def get_transcription_flux(self, solution = None):
		"""Extract the flux state of Transcription reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.TranscriptionReaction):
				for rna_id in reaction.transcription_data.RNA_products:
					locus_id = rna_id.replace('RNA_', '', 1)
					if locus_id not in flux_dict:
						flux_dict[locus_id] = 0
					flux_dict[locus_id] += solution.fluxes[reaction.id]
		return flux_dict

	def get_translation_flux(self, solution = None):
		"""Extract the flux state of Translation reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {r.id: 0 for r in tqdm.tqdm(list(self.translation_data), 'Building reaction dictionary...', bar_format = bar_format)}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.TranslationReaction):
				protein_id = reaction.translation_data.id
				flux_dict[protein_id] += solution.fluxes[reaction.id]
		return flux_dict

	def prune(self, skip = None):
		"""
		Remove all unused metabolites and reactions
		This should be run after the model is fully built. It will be
		difficult to add new content to the model once this has been run.
		skip: list
			List of complexes/proteins/mRNAs/TUs to remain unpruned from model.
		"""
		if not skip:
			skip = []

		#inactive_reactions = [ x for x in self.reactions if x.lower_bound == 0 and x.upper_bound == 0 ]
		#for r in tqdm.tqdm(inactive_reactions, 'Pruning inactive MetabolicReaction\'s...', bar_format = bar_format):
			#logging.warning('Removing inactive MetabolicReaction {}'.format(r.id))
			#r.remove_from_model(remove_orphans = False)

		complex_data_list = [ i.id for i in self.complex_data if i.id not in skip ]
		for c_d in tqdm.tqdm(complex_data_list, 'Pruning unnecessary ComplexData reactions...', bar_format = bar_format):
			c = self.process_data.get_by_id(c_d)
			cplx = c.complex
			if len(cplx.reactions) == 1:
				list(cplx.reactions)[0].delete(remove_orphans = True)
				logging.warning('Removing unnecessary ComplexData reactions for \'{:s}\''.format(c_d))
				self.process_data.remove(self.process_data.get_by_id(c_d))

		for p in tqdm.tqdm(list(self.metabolites.query('_folded')), 'Pruning unnecessary FoldedProtein reactions...', bar_format = bar_format):
			if 'partially' not in p.id and p.id not in skip:
				delete = True
				for rxn in p.reactions:
					if rxn.metabolites[p] < 0:
						delete = False
						break

				if delete:
					while len(p.reactions) > 0:
						list(p.reactions)[0].delete(remove_orphans = True)
						for data in self.process_data.query(p.id):
							logging.warning('Removing unnecessary FoldedProtein reactions for \'{:s}\''.format(p.id))
							self.process_data.remove(data.id)

		for p in tqdm.tqdm(self.metabolites.query('^protein_'), 'Pruning unnecessary ProcessedProtein reactions...', bar_format = bar_format):
			if isinstance(p, coralme.core.component.ProcessedProtein) and p.id not in skip:
				delete = True
				for rxn in p.reactions:
					if rxn.metabolites[p] < 0:
						delete = False
						break
				if delete:
					for rxn in list(p.reactions):
						logging.warning('Removing unnecessary ProcessedProtein reactions for \'{:s}\''.format(rxn.posttranslation_data.id))
						self.process_data.remove(rxn.posttranslation_data.id)
						rxn.delete(remove_orphans = True)

		for p in tqdm.tqdm(self.metabolites.query('^protein_'), 'Pruning unnecessary TranslatedGene reactions...', bar_format = bar_format):
			if isinstance(p, coralme.core.component.TranslatedGene) and p.id not in skip:
				delete = True
				for rxn in p.reactions:
					if rxn.metabolites[p] < 0 and not rxn.id.startswith('degradation'):
						delete = False
						break
				if delete:
					for rxn in p.reactions:
						p_id = p.id.replace('protein_', '')
						data = self.process_data.get_by_id(p_id)
						self.process_data.remove(data.id)
						logging.warning('Removing unnecessary TranslatedGene reactions for \'{:s}\''.format(p_id))
						rxn.delete(remove_orphans = True)

		removed_rna = set()
		for m in tqdm.tqdm(self.metabolites.query('^RNA_'), 'Pruning unnecessary TranscribedGene reactions...', bar_format = bar_format):
			delete = False if m.id in skip else True
			for rxn in m.reactions:
				if rxn.metabolites[m] < 0 and not rxn.id.startswith('DM_'):
					delete = False
			if delete and self.reactions.has_id('DM_' + m.id):
				#try:
					#WARNING: for some reason, m._model returns None and the try/except fails to catch a KeyError at m.remove_from_model
					#self.reactions.get_by_id('DM_' + m.id).remove_from_model(remove_orphans = True)
					#if m in self.metabolites:
						#Defaults to subtractive when removing reaction
						#m.remove_from_model(destructive = False)
				#except KeyError:
					#pass
				self.reactions.get_by_id('DM_' + m.id).remove_from_model(remove_orphans = True)
				try:
					logging.warning('Removing unnecessary TranscribedGene reactions for \'{:s}\''.format(m.id))
					m.remove_from_model(destructive = False)
				except AttributeError:
					logging.warning('AttributeError for \'{:s}\''.format(m.id))
					pass
				removed_rna.add(m.id)

		for t in tqdm.tqdm(self.reactions.query('transcription_TU'), 'Pruning unnecessary Transcriptional Units...', bar_format = bar_format):
			if t.id in skip:
				delete = False
			else:
				delete = True

			for product in t.products:
				if isinstance(product, coralme.core.component.TranscribedGene):
					delete = False

			t_process_id = t.id.replace('transcription_', '')
			if delete:
				t.remove_from_model(remove_orphans = True)
				logging.warning('Removing the unnecessary \'{:s}\' transcriptional unit.'.format(t_process_id))
				self.process_data.remove(t_process_id)
			else:
				# gets rid of the removed RNA from the products
				self.process_data.get_by_id(t_process_id).RNA_products.difference_update(removed_rna)

			# update the TranscriptionReaction mRNA biomass stoichiometry with new RNA_products
			# WARNING: The deletion of RNA(s) from a TU increases the number of nucleotides that should be degraded using the degradosome
			# WARNING: However, n_cuts and n_excised are not recalculated using coralme.builder.transcription.add_rna_splicing
			if not delete:
				t.update()

	def remove_genes_from_model(self, gene_list):
		for gene in tqdm.tqdm(gene_list, 'Removing gene(s) from ME-model...', bar_format = bar_format):
			# defaults to subtractive when removing model
			self.metabolites.get_by_id('RNA_' + gene).remove_from_model()
			protein = self.metabolites.get_by_id('protein_'+gene)
			for cplx in protein.complexes:
				print('Complex \'{:s}\' removed from ME-model.'.format(cplx.id))
				for rxn in cplx.metabolic_reactions:
					try:
						self.process_data.remove(rxn.id.split('_')[0])
					except ValueError:
						pass
					rxn.remove_from_model()

			protein.remove_from_model(destructive = True)

		# Remove all transcription reactions that now do not form a used transcript
		for tu in tqdm.tqdm(self.reactions.query('transcription_TU'), 'Removing unnecessary Transcriptional Units...', bar_format = bar_format):
			delete = True
			for product in tu.products:
				if isinstance(product, coralme.core.component.TranscribedGene):
					delete = False
			if delete:
				tu.remove_from_model(remove_orphans = True)
				t_process_id = tu.id.replace('transcription_', '')
				self.process_data.remove(t_process_id)

	def set_sasa_keffs(self, median_keff):
		# Get median SASA value considering all complexes in model
		sasa_list = []
		for met in tqdm.tqdm(self.metabolites, 'Processing Complexes...', bar_format = bar_format):
			cplx_sasa = 0.
			if not isinstance(met, coralme.core.component.Complex):
				continue
			cplx_sasa += met.formula_weight ** (3. / 4.)
			sasa_list.append(cplx_sasa)
		median_sasa = numpy.median(numpy.array(sasa_list))

		# redo scaling average SASA to 65.
		for rxn in tqdm.tqdm(self.reactions, 'Processing Reactions...', bar_format = bar_format):
			if hasattr(rxn, 'keff') and rxn.complex_data is not None:
				sasa = rxn.complex_data.complex.formula_weight ** (3. / 4.)
				if sasa == 0:
					raise UserWarning('No SASA for reaction \'{:s}\'.'.format(rxn.id))
				rxn.keff = sasa * median_keff / median_sasa

		for data in tqdm.tqdm(self.process_data, 'Processing ProcessData...', bar_format = bar_format):
			sasa = 0.
			if isinstance(data, coralme.core.processdata.TranslocationData):
				continue
			if hasattr(data, 'keff') and hasattr(data, 'formula_weight') and data.enzyme is not None:
				cplxs = [data.enzyme] if type(data.enzyme) == str else data.enzyme
				for cplx in cplxs:
					sasa += self.metabolites.get_by_id(cplx).formula_weight ** (3. / 4.)
				if sasa == 0:
					raise UserWarning('No SASA for reaction \'{:s}\'.'.format(data.id))
				data.keff = sasa * median_keff / median_sasa

		self.update()

	def update(self):
		new = []
		for r in self.reactions:
			if hasattr(r, 'update'):
				new.append(r)
		for r in tqdm.tqdm(new, 'Updating ME-model Reactions...', bar_format = bar_format):
			_update(r)

	# me.update() cannot be parallelized without considering new constraints being added into the model.
	# New constraints must have a different name, so me.update() fails if two reactions are changed to add the same constraint:
	# ContainerAlreadyContains: Container '<optlang.container.Container object at 0x...>' already contains an object with name 'Name'.
	def _parallel_update(self):
		return NotImplemented

	def get(self, x: typing.Union[cobra.core.object.Object, str]) -> cobra.core.object.Object:
		"""
		Return the element with a matching id from model.reactions or model.metabolites attributes.
		"""
		if isinstance(x, cobra.core.object.Object):
			x = x.id
		if isinstance(x, str):
			if self.metabolites.has_id(x):
				return self.metabolites.get_by_id(x)
			elif self.metabolites.has_id(self.aliases['metabolites'].get(x, None)):
				return self.metabolites.get_by_id(self.aliases['metabolites'][x])
			elif self.reactions.has_id(x):
				return self.reactions.get_by_id(x)
			elif self.reactions.has_id(self.aliases['reactions'].get(x, None)):
				return self.reactions.get_by_id(self.aliases['reactions'][x])
			else:
				raise ValueError('Query not found.')
		else:
			return NotImplemented

	def query(self, queries, filter_out_blocked_reactions = False):
		"""
		Return the elements with a matching substring or substrings (AND logic) from
		model.reactions, model.metabolites, and model.process_data attributes.

		For OR logic, use pipe symbol ('|'), e.g. 'ACP|ac'

		Parenthesis and square brackets are allowed without escape symbol.
		"""
		res = []
		if isinstance(queries, list):
			pass
		else:
			queries = [queries]

		# correct queries
		queries = [ x.replace('(', r'\(').replace(')', r'\)').replace('[', r'\[').replace(']', r'\]') for x in queries ]

		for query in queries:
			res.append(self.metabolites.query(query))
			if filter_out_blocked_reactions:
				res.append([ x for x in self.reactions.query(query) if x.bounds != (0, 0) ])
			else:
				res.append(self.reactions.query(query))
			res.append(self.process_data.query(query))

		# compress
		res = [ x for y in res for x in y ]

		if len(queries) > 1:
			# remove from output (AND logic)
			for query in queries[1:]:
				res = [ x for x in res if query in x.id ]

		return res

	# Originally developed by JDTB@UCSD, 2022
	def relax_bounds(self, copy = False):
		if copy:
			test = self.copy()
		else:
			test = self

		for rxn in test.reactions:
			if rxn.id == 'biomass_dilution':
				continue
			if hasattr(rxn.upper_bound, 'subs') or rxn.upper_bound > 0:
				rxn.upper_bound = 1000
			else:
				rxn.upper_bound = 0

			if hasattr(rxn.lower_bound, 'subs') or rxn.lower_bound > 0: # Is this OK?
				rxn.lower_bound = 0
			elif rxn.lower_bound < 0:
				rxn.lower_bound = -1000

		return test

	# Modified from COBRApy
	def _repr_html_(self) -> str:
		"""Get HTML representation of the model.

		Returns
		-------
		str
			Model representation as HTML string.
		"""

		if hasattr(self, 'solution'):
			if self.notes.get('from cobra', False):
				mu = self.solution.objective_value
				dt = numpy.log(2) / mu
			else:
				mu = self.solution.fluxes['biomass_dilution']
				dt = numpy.log(2) / mu
		else:
			mu = dt = numpy.nan

		if hasattr(self, 'process_data'):
			process_data = len(self.process_data)
		else:
			process_data = numpy.nan

		return f"""
		<table>
			<tr>
				<td><strong>Name</strong></td>
				<td>{self.id}</td>
			</tr><tr>
				<td><strong>Memory address</strong></td>
				<td>{f"{id(self):x}"}</td>
			</tr><tr>
				<td><strong>Growth rate</strong></td>
				<td>{mu:g} per hour</td>
			</tr><tr>
				<td><strong>Doubling time</strong></td>
				<td>{dt:g} hours</td>
			</tr><tr>
				<td><strong>Number of metabolites</strong></td>
				<td>{len(self.metabolites)}</td>
			</tr><tr>
				<td><strong>Number of reactions</strong></td>
				<td>{len(self.reactions)}</td>
			</tr><tr>
				<td><strong>Number of process data</strong></td>
				<td>{process_data}</td>
			</tr><tr>
				<td><strong>Number of genes</strong></td>
				<td>{len(self.all_genes)}</td>
			</tr><tr>
				<td><strong>Number of mRNA genes</strong></td>
				<td>{len(self.mRNA_genes)}</td>
			</tr><tr>
				<td><strong>Number of rRNA genes</strong></td>
				<td>{len(self.rRNA_genes)}</td>
			</tr><tr>
				<td><strong>Number of tRNA genes</strong></td>
				<td>{len(self.tRNA_genes)}</td>
			</tr><tr>
				<td><strong>Number of pseudogenes</strong></td>
				<td>{len(self.pseudo_genes)}</td>
			</tr><tr>
				<td><strong>Objective expression</strong></td>
				<td>{cobra.util.util.format_long_string(" + ".join([ '{:.1f}*{:s}'.format(r[1], r[0].id) for r in self.objective]), 100)}</td>
			</tr><tr>
				<td><strong>Compartments</strong></td>
				<td>{", ".join(v if v else k for k, v in
								self.compartments.items())}</td>
			</tr>
			</table>"""
