import re
import pickle
import typing

import logging
log = logging.getLogger(__name__)

# install by the user
import tqdm
bar_format = '{desc:<75}: {percentage:.1f}%|{bar}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
import numpy
import pandas
import scipy
import sympy
import cobra
import coralme

# due to a circular import
from coralme.core.component import Metabolite as Metabolite
from coralme.core.reaction import MEReaction as MEReaction

def _update(MEReaction):
	"""updates all component reactions"""
	MEReaction.update()
	return None

class MEModel(cobra.core.model.Model):
	def __init__(self, name = 'coralME', mu = 'mu'):
		cobra.Model.__init__(self, name)

		self.model_version = coralme.__version__

		self.global_info = {
			'domain' : 'Prokaryote',

			'dnapol_id' : 'DNAP',
			'ribosome_id' : 'ribosome',
			'dummy_rxn_id' : 'dummy_reaction',
			'degradosome_id' : 'RNA_degradosome',
			'mg2_per_ribosome' : 171,
			'amino_acid_loader' : 'generic_Tuf',
			'feature_types' : [ 'CDS', 'rRNA', 'tRNA', 'ncRNA', 'tmRNA', 'misc_RNA' ],

			# analysis
			'add_lipoproteins' : False, #
			'add_translocases' : True, # actually, assign CPLX_dummy to missing enzymes
			'include_pseudo_genes' : False,
			'run_bbh_blast' : True,

			'translation_stop_dict' : {
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

			'default_parameters' : {
				sympy.Symbol('k_t', positive = True) : 4.5,
				sympy.Symbol('r_0', positive = True) : 0.087,
				sympy.Symbol('k^mRNA_deg', positive = True) : 12.0,
				sympy.Symbol('m_rr', positive = True) : 1453.0,
				sympy.Symbol('m_aa', positive = True) : 0.109,
				sympy.Symbol('m_nt', positive = True) : 0.324,
				sympy.Symbol('f_rRNA', positive = True) : 0.86,
				sympy.Symbol('f_mRNA', positive = True) : 0.02,
				sympy.Symbol('f_tRNA', positive = True) : 0.12,
				sympy.Symbol('m_tRNA', positive = True) : 25.0,
				sympy.Symbol('k^default_cat', positive = True) : 65.0,
				sympy.Symbol('temperature', positive = True) : 37.0,
				sympy.Symbol('propensity_scaling', positive = True) : 0.45
				}
			}

		# model parameters as symbols
		self.symbols = {}
		for var in ['k_t', 'r_0', 'k^mRNA_deg', 'm_rr', 'm_aa', 'm_nt', 'f_rRNA', 'f_mRNA', 'f_tRNA', 'm_tRNA', 'k^default_cat', 'temperature', 'propensity_scaling']:
			self.symbols[var] = sympy.Symbol(var, positive = True)

		self.process_data = cobra.core.dictlist.DictList()
		self.metabolites = cobra.core.dictlist.DictList()

		# set growth rate symbolic variable
		self._mu = sympy.Symbol(mu, positive = True)
		# allows the change of symbolic variables through the ME-model object
		self._mu_old = self.mu

		# derived parameters that are common throughout the ME-model
		# WARNING: The equation are written following O'Brien 2013 paper, no COBRAme documentation
		# Empirical relationship between measured ratio of RNA (R) to Protein (P)
		# https://cobrame.readthedocs.io/en/master/coupling_constraint_derivation.html#Parameters
		self.symbols['R/P'] = (self._mu / self.symbols['k_t']) + self.symbols['r_0'] # eq 1, page 15

		# 70S ribosomes (page 16)
		self.symbols['c_ribo'] = self.symbols['m_rr'] / (self.symbols['m_aa'] * self.symbols['f_rRNA']) # eq 2, page 16
		# Hyperbolic ribosome catalytic rate
		self.symbols['k_ribo'] = self.symbols['c_ribo'] * self._mu / self.symbols['R/P']
		# WARNING: the ribosome coupling coefficient in translation reactions is 'v_ribo' times protein length
		self.symbols['v_ribo'] = 1. / (1. * self.symbols['k_ribo'] / self._mu) # page 17

		# RNA Polymerase
		# WARNING: the RNAP coupling coefficient in transcription reactions is 'v_rnap' times RNA length
		self.symbols['v_rnap'] = 1. / (3. * self.symbols['k_ribo'] / self._mu) # page 17

		# mRNA coupling
		self.symbols['c_mRNA'] = self.symbols['m_nt'] / (self.symbols['f_mRNA'] * self.symbols['m_aa']) # page 19
		# Hyperbolic mRNA catalytic rate
		self.symbols['k_mRNA'] = self.symbols['c_mRNA'] * self._mu / self.symbols['R/P']

		# mRNA dilution, degradation, and translation
		self.symbols['alpha_1'] = self._mu / self.symbols['k^mRNA_deg']
		# WARNING: There is an error in O'Brien 2013; corrected in COBRAme docs
		self.symbols['alpha_2'] = self.symbols['R/P'] / (3 * self.symbols['alpha_1'] * self.symbols['c_mRNA'])
		# mRNA dilution, degradation, and translation
		self.symbols['rna_amount'] = self._mu / self.symbols['k_mRNA'] # == 3 * alpha_1 * alpha_2
		self.symbols['deg_amount'] = self.symbols['k^mRNA_deg'] / self.symbols['k_mRNA'] # == 3 * alpha_2

		# tRNA coupling
		self.symbols['c_tRNA'] = self.symbols['m_tRNA'] / (self.symbols['f_tRNA'] * self.symbols['m_aa']) # page 20
		# Hyperbolic tRNA efficiency
		self.symbols['k_tRNA'] = self.symbols['c_tRNA'] * self._mu / self.symbols['R/P']

		# Remaining Macromolecular Synthesis Machinery
		self.symbols['v^default_enz'] = 1. / (1. * (self.symbols['k^default_cat'] / 3600.) / self._mu) # page 20, k^default_cat in 1/s

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

	@property
	def mu(self):
		return self._mu

	@mu.setter
	def mu(self, value):
		# set growth rate symbolic variable
		self._mu_old = self._mu
		self._mu = sympy.Symbol(value, positive = True)

		if self._mu_old == self._mu:
			return # doing nothing because user changed to the current mu

		for rxn in self.reactions:
			if hasattr(rxn.lower_bound, 'subs'):
				rxn._lower_bound = rxn.lower_bound.subs({ self._mu_old : self.mu })
			if hasattr(rxn.upper_bound, 'subs'):
				rxn._upper_bound = rxn.upper_bound.subs({ self._mu_old : self.mu })
			for met, coeff in rxn.metabolites.items():
				if hasattr(coeff, 'subs'):
					rxn._metabolites[met] = coeff.subs({ self._mu_old : self.mu })

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def copy(self):
		return NotImplemented

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def merge(self, right, prefix_existing=None, inplace=True, objective='left'):
		return NotImplemented

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	@property
	def objective(self):
		return [ x for x in self.reactions if x.objective_coefficient != 0 ]

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

	def add_biomass_constraints_to_model(self, biomass_types):
		for biomass_type in tqdm.tqdm(biomass_types, 'Adding biomass constraint(s) into the ME-model...', bar_format = bar_format):
			if '_biomass' not in biomass_type:
				raise ValueError('Biomass types should be suffixed with \'_biomass\'.')
			constraint_obj = coralme.core.component.Constraint(biomass_type)
			summary_variable_obj = coralme.core.reaction.SummaryVariable('{:s}_to_biomass'.format(biomass_type))
			summary_variable_obj.add_metabolites({constraint_obj: -1, self._biomass: 1})
			self.add_reactions([summary_variable_obj])

	@property
	def unmodeled_protein(self):
		return self.metabolites.get_by_id('protein_dummy')

	@property
	def unmodeled_protein_biomass(self):
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

		self.reactions.protein_biomass_to_biomass.add_metabolites({self.unmodeled_protein_biomass: -amount}, combine = False)
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
		if self.reactions.GAM.check_mass_balance() == {'charge': -1.0, 'H': -1.0}:
			self.reactions.GAM._metabolites.update({self.metabolites.h_c : +1})

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
		if self.reactions.ATPM.check_mass_balance() == {'charge': -1.0, 'H': -1.0}:
			self.reactions.ATPM._metabolites.update({self.metabolites.h_c : +1})

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
	def all_genes(self):
		lst = [ g for g in self.metabolites if isinstance(g, coralme.core.component.TranscribedGene) and "dummy" not in g.id]
		return cobra.core.dictlist.DictList(lst)

	@property
	def find_complex(m):
		if isinstance(m,cobrame.core.component.TranslatedGene):
			cplxs = []
			for r in m.reactions:
				cplxs += find_complex(r)
			return cplxs
		if isinstance(m,cobrame.core.reaction.PostTranslationReaction):
			return find_complex(next(i for i in m.metabolites if isinstance(i,cobrame.core.component.ProcessedProtein)))
		if isinstance(m,cobrame.core.component.ProcessedProtein):
			return find_complex(next(i for i in m.reactions if isinstance(i,cobrame.core.reaction.ComplexFormation)))
		if isinstance(m,cobrame.core.reaction.ComplexFormation):
			return find_complex(next(i for i in m.metabolites if isinstance(i,cobrame.core.component.Complex)))
		if isinstance(m,cobrame.core.component.Complex):
			return [m]
		return []

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

	def construct_s_matrix(self, growth_rate) -> scipy.sparse.dok_matrix:
		"""Build the stoichiometric matrix at a specific growth rate."""
		# initialize to 0
		s_matrix = scipy.sparse.dok_matrix((len(self.metabolites), len(self.reactions)))
		# populate with stoichiometry
		for idx, rxn in tqdm.tqdm(list(enumerate(self.reactions)), 'Constructing stoichiometric matrix', bar_format = bar_format):
			for met, value in rxn.metabolites.items():
				met_index = self.metabolites.index(met)
				if hasattr(value, 'subs'):
					s_matrix[met_index, idx] = float(value.subs(self.mu, growth_rate))
				else:
					s_matrix[met_index, idx] = float(value)
		return s_matrix

	def _construct_attribute_vector(self, attr_name, growth_rate):
		"""
		Build a vector of a reaction attribute at a specific growth rate.
		Mainly used for upper and lower bounds.
		"""
		return numpy.array([
			float(value.subs(self.mu, growth_rate))
			if hasattr(value, 'subs') else float(value)
			for value in tqdm.tqdm(self.reactions.list_attr(attr_name), 'Constructing vector of bounds', bar_format = bar_format)
			])

	def compute_solution_error(self, solution = None):
		errors = {}

		if solution is None:
			solution = self.solution

		s_matrix = self.construct_s_matrix(solution.objective_value)
		lb = self._construct_attribute_vector('lower_bound', solution.objective_value)
		ub = self._construct_attribute_vector('upper_bound', solution.objective_value)
		# old code
		#x = numpy.array(solution.x)
		x = numpy.array(list(solution.fluxes.values()))
		err = abs(s_matrix * x)
		errors['max_error'] = err.max()
		errors['sum_error'] = err.sum()
		ub_err = min(ub - x)
		errors['upper_bound_error'] = abs(ub_err) if ub_err < 0 else 0
		lb_err = min(x - lb)
		errors['lower_bound_error'] = abs(lb_err) if lb_err < 0 else 0
		return errors

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
		if isinstance(x, str) or isinstance(x, cobra.core.object.Object):
			if self.metabolites.has_id(x):
				return self.metabolites.get_by_id(x)
			elif self.reactions.has_id(x):
				return self.reactions.get_by_id(x)
			else:
				return
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
		queries = [ x.replace('(', '\(').replace(')', '\)').replace('[', '\[').replace(']', '\]') for x in queries ]

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

	def construct_lp_problem(self, parameters : dict = {}, lambdify = False) -> tuple:
		"""
		lambdify
		    Returns lambda functions for each symbolic stoichiometric coefficient

		Output:
		    A tuple of 9 elements:
		        Dictionary with numeric stoichiometric coefficients
		        Dictionary with symbolic stoichiometric coefficients
		        List of lower bounds (numeric and symbolic)
		        List of upper bounds (numeric and symbolic)
		        List of metabolic bounds (see metabolites._bound property)
		        List of objectives (see reaction.objective_coefficient property)
		        List of constraint senses (always 'E')
		        Set of atoms (i.e., free symbols in symbolic stoichiometric coefficients)
		        Dictionary of lambda functions for each symbolic stoichiometry coefficient
		"""

		# populate empty dictionaries with stoichiometry
		Sf = dict() # floats
		Se = dict() # expressions
		Lr = [ x.id for x in self.reactions ] # reaction identifiers
		Lm = [ x.id for x in self.metabolites ] # reaction identifiers

		# override default parameters
		if hasattr(self, 'global_info'):
			default = self.global_info['default_parameters']
			default.update(parameters)

		# check how many variables are in the ME-model
		atoms = set()

		for idx, rxn in enumerate(self.reactions):
			for met, value in rxn.metabolites.items():
				met_index = self.metabolites.index(met)
				if hasattr(value, 'subs'):
					#atoms.add(list(value.free_symbols)[0])
					atoms.update(list(value.free_symbols))
					Se[met_index, idx] = value.xreplace(default)
				else:
					Sf[met_index, idx] = value

		lb, ub = zip(*[ rxn.bounds for rxn in self.reactions ])
		b = [ m._bound for m in self.metabolites ] # accumulation
		c = [ r.objective_coefficient for r in self.reactions ]
		# constraint sense eventually will be in the metabolite object
		cs = [ 'E' for m in self.metabolites ]

		if lambdify:
			fn = numpy.vectorize(lambda x: sympy.lambdify(list(atoms), x))
			lb = [ x for x in fn(lb) ]
			ub = [ x for x in fn(ub) ]
			lambdas = { k:v for k,v in zip(Se.keys(), fn(list(Se.values()))) }
		else:
			lambdas = None

		return Sf, Se, list(lb), list(ub), b, c, cs, atoms, lambdas, Lr, Lm

	def rank(self, mu = 0.001):
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = self.construct_lp_problem()
		Sp = scipy.sparse.dok_matrix((len(b), len(c)))

		for idx, idj in Sf.keys():
		    Sp[idx, idj] = Sf[idx, idj]

		for idx, idj in Se.keys():
		    Sp[idx, idj] = float(Se[idx, idj].subs({ self.mu : mu }))

		return numpy.linalg.matrix_rank(Sp.todense())

	def fva(self,
		reaction_list, fraction_of_optimum, mu_fixed = None, objective = 'biomass_dilution',
		max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True,
		tolerance = 1e-6, precision = 'quad', verbose = True):

		"""
		Determine the minimum and maximum flux value for each reaction constrained
		to a fraction of the current growth rate (default = 1.0)

		Parameters
		----------
		reaction_list : list of cobra.Reaction or str, optional
			List of reactions IDs and/or reaction objects
		fraction_of_optimum : float, optional
			Must be <= 1.0. Requires that the objective value is at least the
			fraction times maximum objective value. A value of 0.85 for instance
			means that the objective has to be at least at 85% percent of its
			maximum (default 1.0).
		mu_fixed : float, optional
			Set it to avoid the optimization of a ME-model. The growth rate must
			be feasible. If not, the ME-model will be optimized with the following
			options:

			max_mu : float, optional
				Maximum growth rate for initializing the growth rate binary search (GRBS).
			min_mu : float, optional
				Minimum growth rate for initializing GRBS.
			maxIter : int
				Maximum number of iterations for GRBS.
			lambdify : bool
				If True, returns a dictionary of lambda functions for each symbolic
				stoichiometric coefficient
			tolerance : float
				Tolerance for the convergence of GRBS.
			precision : str, {"quad", "double", "dq", "dqq"}
				Precision (quad or double precision) for the GRBS

		verbose : bool
			If True, allow printing.
		"""

		# max_mu is constrained by the fastest-growing bacterium (14.8 doubling time)
		# https://www.nature.com/articles/s41564-019-0423-8

		# check options
		tolerance = tolerance if tolerance >= 1e-15 else 1e-6
		precision = precision if precision in [ 'quad', 'double', 'dq', 'dqq' ] else 'quad'
		fraction_of_optimum = fraction_of_optimum if fraction_of_optimum <= 1.0 and fraction_of_optimum >= 0.0 else 1.0
		if isinstance(reaction_list, str):
			reaction_list = [reaction_list]

		# populate with stoichiometry, no replacement of mu's
		if hasattr(self, 'construct_lp_problem'):
			# check if the ME-model has a solution
			if mu_fixed is not None and not hasattr(self, 'solution'):
				self.optimize(max_mu = max_mu, min_mu = min_mu, maxIter = maxIter, lambdify = lambdify,
					tolerance = tolerance, precision = precision, verbose = verbose)

			# set mu_fixed for replacement in a ME-model.
			mu_fixed = self.solution.fluxes.get(objective, mu_fixed) * fraction_of_optimum

			# get mathematical representation
			Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = self.construct_lp_problem(lambdify = lambdify)
		else:
			# not a ME-model, and objective bounds usually are (0, 1000)
			if self.reactions.has_id(objective):
				self.reactions.get_by_id(objective).lower_bound = mu_fixed * fraction_of_optimum
				self.reactions.get_by_id(objective).upper_bound = mu_fixed
			else:
				raise ValueError('Objective reaction \'{:s}\' not in the M-model.'.format(objective))

			# get mathematical representation
			Sf, Se, lb, ub, b, c, cs, atoms, lambdas = coralme.core.model.MEModel.construct_lp_problem(self)

		if verbose:
			print('Running FVA for {:d} reactions. Maximum growth rate fixed to {:g}'.format(len(reaction_list), mu_fixed))

		from coralme.solver.solver import ME_NLP
		me_nlp = ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)

		# We need only reaction objects
		rxns_fva = []
		for rxn in reaction_list:
			if isinstance(rxn, str) and self.reactions.has_id(rxn):
				rxns_fva.append(self.reactions.get_by_id(rxn))
			else:
				rxns_fva.append(rxn)

		obj_inds0 = [ self.reactions.index(rxn) for rxn in rxns_fva for j in range(0, 2) ]
		obj_coeffs = [ ci for rxn in rxns_fva for ci in (1.0, -1.0) ]

		# varyME is a specialized method for multiple min/maximization problems
		obj_inds0, nVary, obj_vals = me_nlp.varyme(mu_fixed, obj_inds0, obj_coeffs, basis = None, verbosity = verbose)

		# Return result consistent with cobrapy FVA
		fva_result = {
			(self.reactions[obj_inds0[2*i]].id): {
				'maximum':obj_vals[2*i],
				'minimum':obj_vals[2*i+1]
				} for i in range(0, nVary//2) }

		return pandas.DataFrame(fva_result).T

	def _solver_solution_to_cobrapy_solution(self, muopt, xopt, yopt, zopt, stat, solver = 'qminos'):
		if solver == 'qminos':
			#f = sum([ rxn.objective_coefficient * xopt[idx] for idx, rxn in enumerate(self.reactions) ])
			#x_primal = xopt[ 0:len(self.reactions) ]   # The remainder are the slacks
			x_dict = { rxn.id : xopt[idx] for idx, rxn in enumerate(self.reactions) }
			y_dict = { met.id : yopt[idx] for idx, met in enumerate(self.metabolites) }
			z_dict = { rxn.id : zopt[idx] for idx, rxn in enumerate(self.reactions) }
		if solver == 'cplex':
			#x_primal =
			x_dict = { rxn.id: xopt[rxn.id].solution_value for idx, rxn in enumerate(self.reactions) }
			y_dict = { met.id: yopt[met.id].dual_value for idx, met in enumerate(self.metabolites) }
			z_dict = { rxn.id: zopt[rxn.id].reduced_cost for idx, rxn in enumerate(self.reactions) }
		if solver == 'gurobi':
			#x_primal = gpModel.x[ 0:len(self.reactions) ]
			x_dict = { rxn.id : xopt[idx] for idx, rxn in enumerate(self.reactions) }
			y_dict = { met.id : yopt[idx] for idx, met in enumerate(self.metabolites) }
			z_dict = { rxn.id : zopt[idx] for idx, rxn in enumerate(self.reactions) }

		#self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
		return cobra.core.Solution(
			objective_value = muopt,
			status = stat,
			fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
			reduced_costs = z_dict,
			shadow_prices = y_dict,
			)

	def optimize(self,
		max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True, basis = None,
		tolerance = 1e-6, precision = 'quad', verbose = True, get_reduced_costs = False):

		"""Solves the NLP problem to obtain reaction fluxes for a ME-model.

		Parameters
		----------
		max_mu : float
			Maximum growth rate for initializing the growth rate binary search (GRBS).
		min_mu : float
			Minimum growth rate for initializing GRBS.
		maxIter : int
			Maximum number of iterations for GRBS.
		lambdify : bool
			If True, returns a dictionary of lambda functions for each symbolic
			stoichiometric coefficient.
		tolerance : float
			Tolerance for the convergence of GRBS.
		precision : str, {"quad", "double", "dq", "dqq"}
			Precision (quad or double precision) for the GRBS
		verbose : bool
			If True, allow printing.
		get_reduced_costs : bool
			If True, re-optimize but changing the objective function to 'biomass_dilution'
			and its bounds. New reduced costs and shadow prices will be returned.
		"""

		# max_mu is constrained by the fastest-growing bacterium (14.8 min, doubling time)
		# https://www.nature.com/articles/s41564-019-0423-8

		# check options
		min_mu = min_mu if min_mu >= 0. else 0.
		max_mu = max_mu if max_mu <= 2.8100561374051836 else 2.8100561374051836
		tolerance = tolerance if tolerance >= 1e-15 else 1e-6
		precision = precision if precision in [ 'quad', 'double', 'dq', 'dqq' ] else 'quad'

		if hasattr(self, 'troubleshooting') and not self.troubleshooting or not hasattr(self, 'troubleshooting'):
			print('The MINOS and quad MINOS solvers are a courtesy of Prof Michael A. Saunders. Please cite Ma, D., Yang, L., Fleming, R. et al. Reliable and efficient solution of genome-scale models of Metabolism and macromolecular Expression. Sci Rep 7, 40863 (2017). https://doi.org/10.1038/srep40863\n')

		# populate with stoichiometry, no replacement of mu's
		if hasattr(self, 'construct_lp_problem'):
			Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = self.construct_lp_problem(lambdify = lambdify)
		else:
			Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = coralme.core.model.MEModel.construct_lp_problem(self)
			me_nlp = coralme.solver.solver.ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)
			xopt, yopt, zopt, stat, basis = me_nlp.solvelp(.1, None, 'quad', probname = 'lp')

			if stat == 'optimal':
				muopt = [ x for x,c in zip(xopt, c) if c != 0 ][0]
				self.solution = self._solver_solution_to_cobrapy_solution(self, muopt, xopt, yopt, zopt, stat)
				return True
			else:
				if hasattr(self, 'solution'):
					del self.solution
				return False

		if len(atoms) > 1:
			print('Use `me_model.map_feasibility()` to obtain the boundary of feasible solutions.')
			print('Optimization will proceed replacing all growth keys with the same value.')

		me_nlp = coralme.solver.solver.ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)
		muopt, xopt, yopt, zopt, basis, stat = me_nlp.bisectmu(
				mumax = max_mu,
				mumin = min_mu,
				maxIter = maxIter,
				basis = basis,
				tolerance = tolerance,
				precision = precision,
				verbose = verbose)

		if stat == 'optimal':
			# Adapted from Maxwell Neal, 2024
			if get_reduced_costs:
				rxn_idx =  {rxn.id : idx for idx, rxn in enumerate(self.reactions)}
				# Open biomass dilution bounds
				me_nlp.xl[rxn_idx["biomass_dilution"]] =  lambda mu : 0
				me_nlp.xu[rxn_idx["biomass_dilution"]] = lambda mu : 1000
				# Set new objective coefficient
				me_nlp.c = [1.0 if r=="biomass_dilution" else 0.0 for r in rxn_idx]
				# Solve at muopt
				_xopt, yopt, zopt, _stat, _basis = me_nlp.solvelp(muf = muopt, basis = basis, precision = precision)

			self.solution = coralme.core.model.MEModel._solver_solution_to_cobrapy_solution(self, muopt, xopt, yopt, zopt, stat)
			self.basis = basis
			return True
		else:
			if hasattr(self, 'solution'):
				del self.solution
			if hasattr(self, 'basis'):
				self.basis = None
			return False

	# WARNING: Experimental. We could not compile qminos under WinOS, and qminos has a licence restriction for its source code
	def optimize_windows(self,
		max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True,
		tolerance = 1e-6, precision = 'quad', verbose = True, solver = 'gurobi'):

		"""Solves the NLP problem to obtain reaction fluxes for a ME-model. This
		method is used when setting a solver other than qMINOS. It allows to
		use coralME in other OS than Linux.

		Parameters
		----------
		max_mu : float
			Maximum growth rate for initializing the growth rate binary search (GRBS).
		min_mu : float
			Minimum growth rate for initializing GRBS.
		maxIter : int
			Maximum number of iterations for GRBS.
		lambdify : bool
			If True, returns a dictionary of lambda functions for each symbolic
   			stoichiometric coefficient
		tolerance : float
			Tolerance for the convergence of GRBS.
		precision : str, {"quad", "double", "dq", "dqq"}
			Precision (quad or double precision) for the GRBS
		verbose : bool
			If True, allow printing.
		"""

		# check options
		tolerance = tolerance if tolerance >= 1e-15 else 1e-6
		solver = solver if solver in [ 'gurobi', 'cplex' ] else 'gurobi'

		if solver == 'gurobi':
			self.check_feasibility = self.feas_gurobi
		elif solver == 'cplex':
			self.check_feasibility = self.feas_cplex
		else:
			print('The \'solver\' must be \'gurobi\' or \'cplex\'.')

		# populate with stoichiometry with replacement of mu's (Sf contains Se)
		# for multiple evaluations of the LP problem, replacement in lambdify'ed Se is faster overall
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas = self.construct_lp_problem(lambdify = lambdify)

		# test max_mu
		self.check_feasibility(keys = { self.mu:max_mu }, precision = 'quad', **{ 'lp' : [Sf, Se, lb, ub, b, c, cs, atoms, lambdas] })
		if hasattr(self, 'solution') and self.solution.status == 'optimal':
			return True
		else:
			for idx in range(1, maxIter + 1):
				# Just a sequence of feasibility checks
				muf = (min_mu + max_mu) / 2.
				self.check_feasibility(keys = { self.mu:muf }, precision = 'quad', **{ 'lp' : [Sf, Se, lb, ub, b, c, cs, atoms, lambdas] })

				if hasattr(self, 'solution') and self.solution.status == 'optimal':
					stat_new = 'optimal'
					min_mu = muf
				else:
					stat_new = 1
					max_mu = muf

				if verbose:
					print('{:s}\t{:.16f}\t{:s}'.format(str(idx).rjust(9), muf, 'Not feasible' if stat_new == 1 else stat_new.capitalize()))

				if abs(max_mu - min_mu) <= tolerance and stat_new == 'optimal':
					return True

				if max_mu <= tolerance:
					return False

	def feas_windows(self, solver = 'gurobi'):
		if solver == 'gurobi':
			return self.feas_gurobi
		elif solver == 'cplex':
			return self.feas_cplex
		else:
			print('The \'solver\' must be \'gurobi\' or \'cplex\'.')
			return None

	# WARNING: Experimental. We could not compile qminos under WinOS, and qminos has a licence restriction for its source code
	def feas_cplex(self, keys = { sympy.Symbol('mu', positive = True) : 0.1 }, **kwargs):
		# check options
		for key in list(keys.keys()):
			if isinstance(key, sympy.Symbol):
				pass
			else:
				keys[sympy.Symbol(key, positive = True)] = keys.pop(key)

		# populate with stoichiometry with replacement of mu's (Sf contains Se)
		# for single evaluations of the LP problem, direct replacement is faster than lambdify
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas = kwargs.get('lp', self.construct_lp_problem(lambdify = False))

		if lambdas is None:
			Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, Se, lb, ub, keys, atoms)
		else:
			Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, lambdas, lb, ub, keys, atoms)

		from docplex.mp.model import Model

		# create a cplex model
		mpModel = Model(float_precision = 17)

		# Define decision variables
		x = {}
		for idx, rxn in enumerate(self.reactions):
			x[idx] = mpModel.continuous_var(lb = lb[idx], ub = ub[idx], name = rxn.id)

		# Set objective function
		lst = [ x[idx] for idx, rxn in enumerate(self.reactions) if rxn.objective_coefficient != 0 ]
		mpModel.maximize(mpModel.sum(lst))

		# Add constraints for system of linear equations
		for jdx, met in enumerate(self.metabolites):
			lhs = mpModel.linear_expr()
			for idx, rxn in enumerate(self.reactions):
				if (jdx, idx) in Sf: # Sf is a dictionary
					lhs += Sf[(jdx, idx)] * x[idx]
			mpModel.add_constraint(lhs == 0, ctname = met.id)

		mpModel.solve()

		# output solution
		if mpModel.solve_details.status == 'optimal':
			self.solution = self._solver_solution_to_cobrapy_solution(muopt, mpModel._vars_by_name, mpModel._cts_by_name, mpModel._vars_by_name, stat = 'optimal', solver = 'cplex')
			return True
		else:
			if hasattr(self, 'solution'):
				del self.solution
			return False

	# WARNING: Experimental. We could not compile qminos under WinOS, and qminos has a licence restriction for its source code
	def feas_gurobi(self, keys = { sympy.Symbol('mu', positive = True) : 0.1 }, precision = 'quad', **kwargs):
		# check options
		precision = precision if precision in [ 'quad', None, False ] else 'quad'

		for key in list(keys.keys()):
			if isinstance(key, sympy.Symbol):
				pass
			else:
				keys[sympy.Symbol(key, positive = True)] = keys.pop(key)

		# populate with stoichiometry with replacement of mu's (Sf contains Se)
		# for single evaluations of the LP problem, direct replacement is faster than lambdify
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas = kwargs.get('lp', self.construct_lp_problem(lambdify = False))

		if lambdas is None:
			Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, Se, lb, ub, keys, atoms)
		else:
			Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, lambdas, lb, ub, keys, atoms)

		import gurobipy as gp
		from gurobipy import GRB

		# create a gurobi model
		gpModel = gp.Model()

		# Set params
		gpModel.Params.OutputFlag = 0
		gpModel.Params.Presolve = 0
		if precision == 'quad':
			gpModel.Params.Quad = 1
		gpModel.Params.NumericFocus = 3
		gpModel.Params.FeasibilityTol = 1e-9
		gpModel.Params.IntFeasTol = 1e-9
		gpModel.Params.OptimalityTol = 1e-9
		gpModel.Params.Method = 0
		gpModel.Params.BarQCPConvTol = 1.
		#gpModel.Params.PivotTolG = 10.
		#gpModel.Params.UpdateTol = 10.
		#gpModel.Params.SingularTol = 1e-30
		gpModel.Params.BarHomogeneous = 1.

		# Define decision variables
		x = {}
		for idx, rxn in enumerate(self.reactions):
			x[idx] = gpModel.addVar(lb = lb[idx], ub = ub[idx], name = rxn.id, vtype = GRB.CONTINUOUS)

		# Set objective function
		lst = [ x[idx] for idx, rxn in enumerate(self.reactions) if rxn.objective_coefficient != 0 ]
		gpModel.setObjective(gp.quicksum(lst), gp.GRB.MAXIMIZE)

		# Add constraints for system of linear equations
		for jdx, met in enumerate(self.metabolites):
			lhs = gp.LinExpr()
			for idx, rxn in enumerate(self.reactions):
				if (jdx, idx) in Sf: # Sf is a dictionary
					lhs += Sf[(jdx, idx)] * x[idx]
			gpModel.addConstr(lhs == 0)

		# Optimize the model
		gpModel.optimize()

		# output solution
		if gpModel.status == gp.GRB.OPTIMAL:
			self.solution = self._solver_solution_to_cobrapy_solution(muopt, gpModel.x, gpModel.pi, gpModel.RC, stat = 'optimal', solver = 'gurobi')
			return True
		else:
			if hasattr(self, 'solution'):
				del self.solution
			return False

	def feasibility(self, keys = { sympy.Symbol('mu', positive = True) : 0.001 }, tolerance = 1e-6, precision = 'quad', basis = None, **kwargs):
		# check options
		tolerance = tolerance if tolerance >= 1e-15 else 1e-6
		precision = precision if precision in [ 'quad', 'double', 'dq', 'dqq' ] else 'quad'

		for key in list(keys.keys()):
			if isinstance(key, sympy.Symbol):
				pass
			else:
				keys[sympy.Symbol(key, positive = True)] = keys.pop(key)

		# populate with stoichiometry with replacement of mu's (Sf contains Se)
		# for single evaluations of the LP problem, direct replacement is faster than lambdify
		Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = kwargs.get('lp', self.construct_lp_problem(lambdify = False))

		if lambdas is None:
			Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, Se, lb, ub, keys, atoms)
		else:
			Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, lambdas, lb, ub, keys, atoms)

		from coralme.solver.solver import ME_NLP
		#me_nlp = ME_NLP(me)
		me_nlp = ME_NLP(Sf, dict(), b, c, lb, ub, cs, set(keys.keys()), None)
		muopt, xopt, yopt, zopt, basis, stat = me_nlp.bisectmu(
				mumax = 1., # mu was already replaced and maxIter is one, so a value here doesn't matter
				mumin = 0.,
				maxIter = 1,
				basis = basis,
				tolerance = tolerance,
				precision = precision,
				verbose = False)

		if stat == 'optimal':
			if len(self.reactions) > 1 and len(self.metabolites) > 1:
				self.solution = self._solver_solution_to_cobrapy_solution(list(keys.values())[0], xopt, yopt, zopt, stat)
			else:
				x_primal = xopt[ 0:len(Lr) ]   # The remainder are the slacks
				x_dict = { rxn : xopt[idx] for idx, rxn in enumerate(Lr) }
				y_dict = { met : yopt[idx] for idx, met in enumerate(Lm) }
				z_dict = { rxn : zopt[idx] for idx, rxn in enumerate(Lr) }
				self.solution = cobra.core.Solution(
					objective_value = muopt,
					status = stat,
					fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
					reduced_costs = z_dict,
					shadow_prices = y_dict,
					)

			self.basis = basis
			return True
		else:
			if hasattr(self, 'solution'):
				del self.solution
			if hasattr(self, 'basis'):
				self.basis = None
			return False

	def map_feasibility(self, keys = { sympy.Symbol('mu', positive = True) : 1. }, tolerance = 1e-6, precision = 'quad'):
		return NotImplemented

	# Originally developed by JDTB@UCSD, 2022
	def relax_bounds(self):
		for rxn in self.reactions:
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

	# Modified from COBRApy
	def _repr_html_(self) -> str:
		"""Get HTML representation of the model.

		Returns
		-------
		str
			Model representation as HTML string.
		"""
		return f"""
		<table>
			<tr>
				<td><strong>Name</strong></td>
				<td>{self.id}</td>
			</tr><tr>
				<td><strong>Memory address</strong></td>
				<td>{f"{id(self):x}"}</td>
			</tr><tr>
				<td><strong>Number of metabolites</strong></td>
				<td>{len(self.metabolites)}</td>
			</tr><tr>
				<td><strong>Number of reactions</strong></td>
				<td>{len(self.reactions)}</td>
			</tr><tr>
				<td><strong>Number of genes</strong></td>
				<td>{len(self.all_genes)-1}</td>
			</tr><tr>
				<td><strong>Number of groups</strong></td>
				<td>{len(self.groups)}</td>
			</tr><tr>
				<td><strong>Objective expression</strong></td>
				<td>{cobra.util.util.format_long_string(",".join([r.id for r in self.objective]), 100)}</td>
			</tr><tr>
				<td><strong>Compartments</strong></td>
				<td>{", ".join(v if v else k for k, v in
								self.compartments.items())}</td>
			</tr>
			</table>"""
