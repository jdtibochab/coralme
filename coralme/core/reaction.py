import numpy
import sympy
import cobra
import logging

import coralme
# use this because recursive import leads to a parcial import and an error
from coralme.core.component import Metabolite as Metabolite

from collections import defaultdict
from operator import attrgetter

class MEReaction(cobra.core.reaction.Reaction):
	# TODO set _upper and _lower bounds as a property
	"""
	MEReaction is a general reaction class from which all ME-Model reactions
	will inherit.

	This class contains functionality that can be used by all ME-model
	reactions.

	Parameters
	----------
	id : str
		Identifier of the MEReaction. Should follow best practices of child
		class

	"""
	def __init__(self, id = None, name = ''):
		cobra.core.reaction.Reaction.__init__(self, id, name)
		self._objective_coefficient = 0.

	@property
	def objective_coefficient(self):
		"""
		Get and set objective coefficient of reaction

		Overrides method in parent class in order to enable use of optlang
		interfaces.

		Returns
		-------
		float
			Objective coefficient of reaction
		"""
		return self._objective_coefficient

	@objective_coefficient.setter
	def objective_coefficient(self, value):
		self._objective_coefficient = value

	def check_me_mass_balance(self):
		"""
		Checks the mass balance of ME reaction, ignoring charge balances

		Returns
		-------
		dict
			{element: number_of_elemental_imbalances}

		"""
		mass_balance = self.check_mass_balance()

		# ME-model is not currently charge balanced
		if 'charge' in mass_balance:
			mass_balance.pop('charge')

		return {met: value for met, value in mass_balance.items() if abs(value) > 1e-11}

	def add_subreactions(self, process_data_id, stoichiometry, scale = 1.):
		"""
		Function to add subreaction process data to reaction stoichiometry

		Parameters
		----------
		process_data_id : str
			ID of the process data associated with the metabolic reaction.

			For example, if the modifications are being added to a complex
			formation reaction, the process data id would be the name of the
			complex.

		stoichiometry : dict
			Dictionary of {metabolite_id: float} or
			{metabolite_id: float * (sympy.Symbol)}

		scale : float
		   Some processes (ie. tRNA charging) are reformulated such that other
		   involved metabolites need scaling

		Returns
		-------
		dict
			Stoichiometry dictionary with updated entries
		"""
		process_info = self._model.process_data.get_by_id(process_data_id)
		for subreaction_id, count in process_info.subreactions.items():
			subreaction_data = self._model.process_data.get_by_id(subreaction_id)

			if isinstance(subreaction_data.enzyme, list) or isinstance(subreaction_data.enzyme, set):
				for enzyme in subreaction_data.enzyme:
					stoichiometry[enzyme] -= self._model.mu / subreaction_data.keff / 3600. * count * scale

			elif isinstance(subreaction_data.enzyme, str):
				stoichiometry[subreaction_data.enzyme] -= self._model.mu / subreaction_data.keff / 3600. * count * scale

			for met, stoich in subreaction_data.stoichiometry.items():
				stoichiometry[met] += count * stoich * scale

		return stoichiometry

	def get_components_from_ids(self, id_stoichiometry, default_type = Metabolite, verbose = True):
		"""
		Function to convert stoichiometry dictionary entries from strings to
		cobra objects.

		{metabolite_id: value} to {:class:`coralme.core.component.Metabolite`:
		value}

		Parameters
		----------
		id_stoichiometry: Dict {string: float}
			Input Dict of {metabolite_id: value}

		default_type: String
			The type of cobra.Metabolite to default to if the metabolite is not
			yet present in the model

		verbose: Boolean
			If True, print metabolites added to model if not yet present in
			model

		Returns
		-------
		dict
			{:class:`coralme.core.component.Metabolite`: float}
		"""

		stoic = id_stoichiometry
		object_stoichiometry = {}
		mets = self._model.metabolites
		for key, value in stoic.items():
			try:
				object_stoichiometry[mets.get_by_id(key)] = value
			except KeyError:
				new_met = coralme.core.component.create_component(key, default_type = default_type)
				if verbose:
					#logging.warning('Metabolite created \'{:s}\' in ME-Model \'{:s}\'.'.format(repr(new_met), repr(self)))
					logging.warning('Metabolite \'{:s}\' created in Reaction \'{:s}\''.format(new_met.id, self.id))
				object_stoichiometry[new_met] = value
				self._model.add_metabolites([new_met])
		return object_stoichiometry

	def add_biomass_from_subreactions(self, process_data, biomass = 0.):
		"""
		Account for the biomass of metabolites added to macromolecule (protein,
		complex, etc.) due to a modification such as prosthetic group addition.

		Parameters
		----------
		process_data : :class:`coralme.core.processdata.ProcessData`
			ProcessData that is used to construct MEReaction

		biomass : float
			Initial biomass value in kDa

		Returns
		-------
		float
			Initial biomass value + biomass added from subreactions in kDa

		"""
		for subrxn, count in process_data.subreactions.items():
			subrxn_obj = self._model.process_data.get_by_id(subrxn)
			biomass += subrxn_obj.calculate_biomass_contribution() / 1000. * count
		return biomass  # in kDa

	def clear_metabolites(self):
		"""
		Remove all metabolites from the reaction
		"""
		for metabolite in list(self._metabolites.keys()):
			self.add_metabolites({metabolite: 0}, combine = False)

	# overwrite methods from cobrapy
	def _set_id_with_model(self, value: str) -> None:
		"""Set Reaction id in model, check that it doesn't already exist.

		The function will rebuild the model reaction index.

		Parameters
		----------
		value: str
			A string that represents the id.

		Raises
		------
		ValueError
			If the model already contains a reaction with the id value.
		"""
		if value in self.model.reactions:
			raise ValueError(
				f"The model already contains a reaction with the id: {value}"
			)
		#forward_variable = self.forward_variable
		#reverse_variable = self.reverse_variable
		self._id = value
		self.model.reactions._generate_index()
		#forward_variable.name = self.id
		#reverse_variable.name = self.reverse_id

	def _check_bounds(self, lb, ub):
		#logging.warning('New cobraME \'_check_bounds\' method superseeds \'_check_bounds\' from cobrapy')
		if isinstance(lb, float) and isinstance(ub, float):
			if lb > ub:
				raise ValueError('The lower bound must be less than or equal to the upper bound ({:s} <= {:s}).'.format(lb, ub))

	def update_variable_bounds(self):
		#logging.warning('New cobraME \'update_variable_bounds\' method superseeds \'update_variable_bounds\' from cobrapy')

		if self.model is None:
			return

		# sympy.core.symbol.Symbol > 0 = True
		if isinstance(self.lower_bound, sympy.core.symbol.Symbol) or isinstance(self.lower_bound, sympy.core.mul.Mul) or isinstance(self.lower_bound, sympy.core.add.Add):
			lb = self.lower_bound
		elif isinstance(self.upper_bound, sympy.core.symbol.Symbol) or isinstance(self.upper_bound, sympy.core.mul.Mul) or isinstance(self.upper_bound, sympy.core.add.Add):
			ub = self.upper_bound
		else:
			pass

		# We know that `lb <= ub`.
		if isinstance(self.lower_bound, float) and isinstance(self.upper_bound, float):
			if self.lower_bound > 0:
				self.forward_variable.set_bounds(
					lb = None if numpy.isinf(self._lower_bound) else self.lower_bound,
					ub = None if numpy.isinf(self._upper_bound) else self.upper_bound,
					)
				self.reverse_variable.set_bounds(lb = 0, ub = 0)
			elif self.upper_bound < 0:
				self.forward_variable.set_bounds(lb = 0, ub = 0)
				self.reverse_variable.set_bounds(
					lb = None if numpy.isinf(self.upper_bound) else -self.upper_bound,
					ub = None if numpy.isinf(self.lower_bound) else -self.lower_bound,
					)
			else:
				self.forward_variable.set_bounds(lb = 0, ub = None if numpy.isinf(self.upper_bound) else +self.upper_bound)
				self.reverse_variable.set_bounds(lb = 0, ub = None if numpy.isinf(self.lower_bound) else -self.lower_bound)

	def build_reaction_string(self, use_metabolite_names: bool = False) -> str:
		"""Generate a human readable reaction str.

		Parameters
		----------
		use_metabolite_names: bool
			Whether to use metabolite names (when True) or metabolite ids (when False,
			default).

		Returns
		-------
		str
			A human readable str.
		"""

		def _format(number) -> str:
			return "1.0 " if number == 1 else "{:s} ".format(str(number).rstrip("."))

		id_type = "id"
		if use_metabolite_names:
			id_type = "name"
		reactant_bits = []
		product_bits = []
		for met in sorted(self.metabolites, key = attrgetter("id")):
			coefficient = self.metabolites[met]
			name = str(getattr(met, id_type)) if str(getattr(met, id_type)) != '' else met.id
			if isinstance(coefficient, sympy.core.symbol.Symbol) or isinstance(coefficient, sympy.core.mul.Mul) or isinstance(coefficient, sympy.core.add.Add):
				if coefficient.subs([(self._model.mu, 1e-6)]) >= 0:
					product_bits.append('[{:s}] '.format(_format(coefficient).strip()) + name)
				else:
					reactant_bits.append('[{:s}] '.format(_format(coefficient * -1).strip()) + name)
			elif coefficient >= 0:
				product_bits.append(_format(coefficient) + name)
			else:
				reactant_bits.append(_format(abs(coefficient)) + name)

		reaction_string = " + ".join(reactant_bits)
		try:
			if not self.reversibility:
				if self.lower_bound < 0 and self.upper_bound <= 0:
					reaction_string += " <-- "
				else:
					reaction_string += " --> "
			else:
				reaction_string += " <=> "
		except:
			reaction_string += " --> "
		reaction_string += " + ".join(product_bits)
		return reaction_string

	def _repr_html_(self) -> str:
			"""Generate html representation of reaction.

			Returns
			-------
			str
				HTML representation of the reaction.
			"""
			rxn = cobra.util.util.format_long_string(str(self.id), 500)
			name = cobra.util.util.format_long_string(str(self.name), 500)
			subs = cobra.util.util.format_long_string(self.build_reaction_string(), 1000)
			prod = cobra.util.util.format_long_string(self.build_reaction_string(True), 1000)
			gpr = cobra.util.util.format_long_string(self.gene_reaction_rule, 500)
			lower = self.lower_bound
			upper = self.upper_bound

			return f"""
			<table>
				<tr><td><strong>Reaction identifier</strong></td><td>{rxn}</td></tr>
				<tr><td><strong>Name</strong></td><td>{name}</td></tr>
				<tr><td><strong>Memory address</strong></td><td>{f'{id(self):#x}'}</td></tr>
				<tr><td><strong>Stoichiometry</strong></td><td>
					<p style='text-align:right'>{subs}</p>
					<p style='text-align:right'>{prod}</p>
				</td></tr>
				<tr><td><strong>GPR</strong></td><td>{gpr}</td></tr>
				<tr><td><strong>Lower bound</strong></td><td>{lower}</td></tr>
				<tr><td><strong>Upper bound</strong></td><td>{upper}</td></tr>
				<tr><td><strong>Type</strong></td><td>{f'{repr(self):s}'}</td></tr>
			</table>
		"""

class MetabolicReaction(MEReaction):
	"""Irreversible metabolic reaction including required enzymatic complex

	This reaction class's update function processes the information contained
	in the complex data for the enzyme that catalyzes this reaction as well as
	the stoichiometric data which contains the stoichiometry of the metabolic
	conversion being performed (i.e. the stoichiometry of the M-model reaction
	analog)

	Parameters
	----------
	id : str
		Identifier of the metabolic reaction. As a best practice, this
		ID should use the following template (FWD=forward, REV=reverse):
		'<StoichiometricData.id> + _ + <FWD or REV> + _ + <Complex.id>'

	Attributes
	----------
	keff : float
		The turnover rete (keff) couples enzymatic dilution to metabolic flux
	reverse : boolean
		If True, the reaction corresponds to the reverse direction of the
		reaction. This is necessary since all reversible enzymatic reactions
		in an ME-model are broken into two irreversible reactions

	"""

	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._complex_data = None
		self._stoichiometric_data = None
		self.keff = 65.  # in per second
		self.reverse = False

	@property
	def complex_data(self):
		"""
		Get or set the ComplexData instance that details the enzyme that
		catalyzes the metabolic reaction.  Can be set with instance of
		ComplexData or with its id.

		Returns
		-------
		:class:`coralme.core.processdata.ComplexData`
			Complex data detailing enzyme that catalyzes this reaction
		"""
		return self._complex_data

	@complex_data.setter
	def complex_data(self, process_data):
		if isinstance(process_data, str):
			process_data = self._model.process_data.get_by_id(process_data)
		self._complex_data = process_data
		#if not hasattr(process_data, 'complex_id'):
			#raise TypeError('The \'{:s}\' is not a ComplexData instance.'.format(process_data.id))
		if process_data is not None:
			process_data._parent_reactions.add(self.id)

	@property
	def stoichiometric_data(self):
		"""
		Get or set the StoichiometricData instance that details the metabolic
		conversion of the metabolic reaction.  Can be set with instance of
		StoichiometricData or with its id.

		Returns
		-------
		:class:`coralme.core.processdata.StoichiometricData`
		   Stoichiometric data detailing enzyme that catalyzes this reaction
		"""
		return self._stoichiometric_data

	@stoichiometric_data.setter
	def stoichiometric_data(self, process_data):
		if isinstance(process_data, str):
			process_data = self._model.process_data.get_by_id(process_data)
		self._stoichiometric_data = process_data
		process_data._parent_reactions.add(self.id)

	def update(self, verbose = True):
		"""
		Creates reaction using the associated stoichiometric data and
		complex data.

		This function adds the following components to the reaction
		stoichiometry (using 'data' as shorthand for
		:class:`coralme.core.processdata.StoichiometricData`):

		1) Complex w/ coupling coefficients defined in self.complex_data.id
		   and self.keff

		2) Metabolite stoichiometry defined in data.stoichiometry. Sign is
		   flipped if self.reverse == True

		Also sets the lower and upper bounds based on self.reverse and
		data.upper_bound and data.lower_bound.

		Parameters
		----------
		verbose : bool
			Prints when new metabolites are added to the model when executing
			update()
		"""
		self.clear_metabolites()
		new_stoichiometry = defaultdict(float)
		stoichiometric_data = self.stoichiometric_data

		# Add complex if enzyme catalyzed
		if self.complex_data:
			new_stoichiometry[self.complex_data.complex.id] = -self._model.mu / self.keff / 3600.  # s-1 / (3600 s/h)

		# Update new stoichiometry values
		sign = -1 if self.reverse else 1
		for component, value in stoichiometric_data.stoichiometry.items():
			new_stoichiometry[component] += value * sign

		new_stoichiometry = self.add_subreactions(stoichiometric_data.id, new_stoichiometry)

		# Convert component ids to cobra metabolites
		object_stoichiometry = self.get_components_from_ids(new_stoichiometry, verbose = verbose)

		# Replace old stoichiometry with new one
		try:
			self.add_metabolites(object_stoichiometry)
		except:
			print('core/reaction.py:422 ' + str(object_stoichiometry))

		# Set the bounds
		if self.reverse:
			self.lower_bound = max(0, -self.stoichiometric_data.upper_bound)
			self.upper_bound = max(0, -self.stoichiometric_data.lower_bound)
		else:
			self.lower_bound = max(0, +self.stoichiometric_data.lower_bound)
			self.upper_bound = max(0, +self.stoichiometric_data.upper_bound)

class ComplexFormation(MEReaction):
	"""Formation of a functioning enzyme complex that can act as a catalyst for
	a ME-model reaction.

	This reaction class produces a reaction that combines the protein subunits
	and adds any coenyzmes, prosthetic groups or enzyme modifications to form
	complete enzyme complex.

	Parameters
	----------
	id : str
		Identifier of the complex formation reaction. As a best practice, this
		ID should be prefixed with 'formation + _ + <complex_id>'. If there
		are multiple ways of producing complex, this can be suffixed with
		'_ + alt'

	Attributes
	----------
	_complex_id : str
		Name of the complex being produced by the complex formation reaction

	complex_data_id : str
		Name of ComplexData that defines the subunit stoichiometry or
		subreactions (modfications). This will not always be the same as the
		_complex_id. Sometimes complexes can be modified using different
		processes/enzymes

	"""
	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._complex_id = None
		self.complex_data_id = None

	@property
	def complex(self):
		"""
		Get the metabolite product of the complex formation reaction

		Returns
		-------
		:class:`coralme.core.component.Complex`
			Instance of complex metabolite from self._complex_id
		"""
		return self._model.metabolites.get_by_id(self._complex_id)

	def _add_formula_to_complex(self, complex_data, complex_met):
		"""
		Add chemical formula as sum of all protein and modification components
		detailed in subreaction data.

		Parameters
		----------
		complex_data : :class:`coralme.core.processdata.ComplexData`
			Complex data for complex being formed in the reaction

		complex_met : :class:`coralme.core.processdata.ComplexData`
			Metabolite of complex being formed in the reaction

		"""
		elements = defaultdict(int)
		for component, count in complex_data.stoichiometry.items():
			component_obj = self._model.metabolites.get_by_id(component)
			for e, n in component_obj.elements.items():
				elements[e] += n * count

		elements = coralme.util.massbalance.get_elements_from_process_data(self, complex_data, elements)

		# Convert element dict to formula string and associate it with complex
		coralme.util.massbalance.elements_to_formula(complex_met, elements)

	def update(self, verbose=True):
		"""
		Creates reaction using the associated complex data and adds chemical
		formula to complex metabolite product.

		This function adds the following components to the reaction
		stoichiometry (using 'data' as shorthand for
		:class:`coralme.core.processdata.ComplexData`):

		1) Complex product defined in self._complex_id

		2) Protein subunits with stoichiometry defined in data.stoichiometry

		3) Metabolites and enzymes w/ coupling coefficients defined in
		   data.subreactions. This often includes enzyme complex
		   modifications by coenzymes or prosthetic groups.

		4) Biomass :class:`coralme.core.component.Constraint` corresponding to
		   modifications detailed in data.subreactions, if any

		Parameters
		----------
		verbose : bool
			Prints when new metabolites are added to the model when executing
			update()
		"""
		self.clear_metabolites()
		stoichiometry = defaultdict(float)
		metabolites = self._model.metabolites
		complex_info = self._model.process_data.get_by_id(self.complex_data_id)

		# Find or create complex product and add it to stoichiometry dict
		try:
			complex_met = metabolites.get_by_id(self._complex_id)
		except KeyError:
			complex_met = coralme.core.component.create_component(self._complex_id, default_type = coralme.core.component.Complex)
			self._model.add_metabolites([complex_met])
		stoichiometry[complex_met.id] = 1

		# build the complex itself
		for component_id, value in complex_info.stoichiometry.items():
			stoichiometry[component_id] -= value

		# add in the subreactions and modifications
		stoichiometry = self.add_subreactions(complex_info.id, stoichiometry)

		# convert string stoichiometry representations to coralme metabolites
		object_stoichiometry = self.get_components_from_ids(stoichiometry, default_type = coralme.core.component.Complex, verbose = verbose)

		# Add formula to complex
		self._add_formula_to_complex(complex_info, complex_met)

		# Biomass accounting of protein subunits is handled in translation
		# reactions. Handle cofactors and prosthetic groups here
		biomass = 0.
		biomass = self.add_biomass_from_subreactions(complex_info, biomass)
		if biomass > 0:
			self.add_metabolites({metabolites.prosthetic_group_biomass: biomass})

		self.add_metabolites(object_stoichiometry, combine = False)

class PostTranslationReaction(MEReaction):
	"""
	Reaction class that includes all posttranslational modification reactions
	(translocation, protein folding, modification (for lipoproteins) etc)

	There are often multiple different reactions/enzymes that can accomplish
	the same modification/function. In order to account for these and
	maintain one translation reaction per protein, these processes need to be
	modeled as separate reactions.

	Parameters
	----------
	id : str
		Identifier of the post translation reaction

	"""
	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._posttranslation_data = None

	@property
	def posttranslation_data(self):
		"""
		Get or set PostTranslationData that defines the type of post
		translation modification/process (folding/translocation) that the
		reaction accounts for. Can be set with instance of
		PostTranslationData or with its id.

		Returns
		-------
		:class:`coralme.core.processdata.PostTranslationData`
			The PostTranslationData that defines the PostTranslationReaction

		"""
		return self._posttranslation_data

	@posttranslation_data.setter
	def posttranslation_data(self, process_data):
		if isinstance(process_data, str):
			process_data = self._model.process_data.get_by_id(process_data)
		self._posttranslation_data = process_data
		process_data._parent_reactions.add(self.id)

	def add_translocation_pathways(self, process_data_id, protein_id, stoichiometry = None):
		"""
		Add complexes and metabolites required to translocate the protein into
		cell membranes.

		Parameters
		----------
		process_data_id : str
			ID of translocation data defining post translation reaction

		protein_id : str
			ID of protein being translocated via post translation reaction

		stoichiometry : dict
			Dictionary of {metabolite_id: float} or
			{metabolite_id: float * (sympy.Symbol)}

		Returns
		-------
		dict
			Stoichiometry dictionary with updated entries from translocation
		"""
		if not stoichiometry:
			stoichiometry = defaultdict(float)

		process_info = self._model.process_data.get_by_id(process_data_id)
		protein = self._model.metabolites.get_by_id(protein_id)
		protein_length = len(protein.amino_acid_sequence)

		for translocation in process_info.translocation:
			translocation_data = self._model.process_data.get_by_id(translocation)
			for metabolite, amount in translocation_data.stoichiometry.items():
				if translocation_data.length_dependent_energy:
					stoichiometry[metabolite] += amount * protein_length
				else:
					stoichiometry[metabolite] += amount

			# Requirement of some translocation complexes vary depending
			# on protein being translocated
			multiplier_dict = process_info.translocation_multipliers
			for enzyme, enzyme_info in translocation_data.enzyme_dict.items():
				length_dependent = enzyme_info['length_dependent']
				fixed_keff = enzyme_info['fixed_keff']
				multiplier = multiplier_dict.get(enzyme, 1.)
				length = protein_length if length_dependent else 1.

				# keff = translocation_data.keff
				keff = 65. if fixed_keff else translocation_data.keff / length
				enzyme_stoichiometry = multiplier * self._model.mu / keff / 3600.
				stoichiometry[enzyme] -= enzyme_stoichiometry

		return stoichiometry

	def update(self, verbose = True):
		"""
		Creates reaction using the associated posttranslation data and adds
		chemical formula to processed protein product

		This function adds the following components to the reaction
		stoichiometry (using 'data' as shorthand for
		:class:`coralme.core.processdata.PostTranslationData`):

		1) Processed protein product defined in data.processed_protein_id

		2) Unprocessed protein reactant defined in data.unprocessed_protein_id

		3) Metabolites and enzymes defined in data.subreactions

		4) Translocation pathways defined in data.translocation

		5) Folding mechanism defined in data.folding_mechanims w/ coupling
		   coefficients defined in data.keq_folding, data.k_folding,
		   model.global_info['temperature'], data.aggregation_propensity,
		   and data.propensity_scaling

		6) Surface area constraints defined in data.surface_are

		7) Biomass if a significant chemical modification takes place (i.e.
		   lipid modifications for lipoproteins)

		Parameters
		----------
		verbose : bool
			Prints when new metabolites are added to the model when executing
			update()

		"""
		self.clear_metabolites()
		stoichiometry = defaultdict(float)
		metabolites = self._model.metabolites
		posttranslation_data = self.posttranslation_data
		unprocessed_protein = posttranslation_data.unprocessed_protein_id
		processed_protein = posttranslation_data.processed_protein_id

		# folding properties
		folding_mechanism = posttranslation_data.folding_mechanism
		aggregation_propensity = posttranslation_data.aggregation_propensity
		scaling = posttranslation_data.propensity_scaling
		if folding_mechanism:
			temp = str(self._model.global_info['temperature'])
			keq_folding = posttranslation_data.keq_folding[temp]
			k_folding = posttranslation_data.k_folding[temp] * 3600.  # in hr-1

		# Get or make processed protein metabolite
		try:
			protein_met = metabolites.get_by_id(processed_protein)
		except KeyError:
			protein_met = coralme.core.component.ProcessedProtein(processed_protein, unprocessed_protein)
			self._model.add_metabolites(protein_met)

		# Add subreactions (e.g. lipid modifications for lipoproteins)
		stoichiometry = self.add_subreactions(posttranslation_data.id, stoichiometry)

		# Add translocation pathways, if applicable
		if posttranslation_data.translocation:
			stoichiometry = self.add_translocation_pathways(
				posttranslation_data.id, unprocessed_protein, stoichiometry)

		# Add folding protein coupling coefficients, if applicable
		if folding_mechanism == 'folding_spontaneous':
			dilution = (keq_folding + self._model.mu / k_folding)
			stoichiometry[unprocessed_protein] -= (dilution + 1.)
			stoichiometry[protein_met.id] += 1.

		elif folding_mechanism:
			dilution = aggregation_propensity * scaling * (keq_folding + 1.) + 1.
			stoichiometry[unprocessed_protein] -= (1. / dilution + 1.)
			stoichiometry[protein_met.id] += 1. / dilution
			stoichiometry[protein_met.id.replace('_folded', '')] += (1.)
		else:
			stoichiometry[unprocessed_protein] = -1.
			stoichiometry[protein_met.id] = 1.

		# Add surface area constraints for all translocated proteins, if applicable
		surface_area = posttranslation_data.surface_area
		if surface_area:
			for SA, value in surface_area.items():
				try:
					sa_constraint = metabolites.get_by_id(SA)
				except KeyError:
					logging.warning('Constraint \'{:s}\' added to ME-Model.'.format(SA))
					sa_constraint = coralme.Constraint(SA)
					self._model.add_metabolites([sa_constraint])

				stoichiometry[sa_constraint.id] += value

		# Convert metabolite strings to metabolite objects
		object_stoichiometry = self.get_components_from_ids(stoichiometry, verbose = verbose)

		# Add formula as sum of unprocessed protein and modification components
		elements = defaultdict(int)
		elements.update(metabolites.get_by_id(unprocessed_protein).elements)
		elements = coralme.util.massbalance.get_elements_from_process_data(self, posttranslation_data, elements)

		# Convert element dict to formula string and associate it with protein
		coralme.util.massbalance.elements_to_formula(protein_met, elements)

		# Add biomass from significant modifications (i.e. lipids for lipoproteins)
		biomass = self.add_biomass_from_subreactions(posttranslation_data)
		if biomass > 0 and posttranslation_data.biomass_type:
			self.add_metabolites({metabolites.get_by_id(posttranslation_data.biomass_type): biomass})
		elif biomass > 0 and not posttranslation_data.biomass_type:
			raise ValueError('If SubReactions in PostTranslationData modify the protein, the \'biomass_type\' must be provided.')

		self.add_metabolites(object_stoichiometry, combine = False)

class TranscriptionReaction(MEReaction):
	"""Transcription of a TU to produced TranscribedGene.

	RNA is transcribed on a transcription unit (TU) level. This type of
	reaction produces all of the RNAs contained within a TU, as well as
	accounts for the splicing/excision of RNA between tRNAs and rRNAs.
	The appropriate RNA_biomass constrain is produced based on the molecular
	weight of the RNAs being transcribed

	Parameters
	----------
	id : str
		Identifier of the transcription reaction. As a best practice, this ID
		should be prefixed with 'transcription + _'

	"""

	# TODO double check how initiation is used as well as ATP cost etc.
	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._transcription_data = None

	@property
	def transcription_data(self):
		"""
		Get or set the :class:`coralme.core.processdata.TranscriptionData`
		that defines the transcription unit architecture and the features of
		the RNAs being transcribed.

		"""
		return self._transcription_data

	@transcription_data.setter
	def transcription_data(self, process_data):
		if isinstance(process_data, str):
			process_data = self._model.process_data.get_by_id(process_data)
		self._transcription_data = process_data
		process_data._parent_reactions.add(self.id)

	def _add_formula_to_transcript(self, transcript):
		"""

		Add element formula to transcript based on nucleotide composition.
		1 OH group is removed for each nucleotide to account for polymerization
		of mononucleotides. This was done to instead of considering the 3'
		diphosphate group as a simplification to avoid keeping track of the
		3' nucleotide in cases of transcription unit splicing.

		Parameters
		----------
		transcript : :class:`cobra.core.component.TranscribedGene`
			Instance of gene being transcribed

		"""

		elements = defaultdict(int)

		for nuc, value in transcript.nucleotide_count.items():
			nuc_obj = self._model.metabolites.get_by_id(nuc)
			for e, n in nuc_obj.elements.items():
				elements[e] += value * n

		# Remove -OH for each
		elements['H'] -= len(transcript.nucleotide_sequence)
		elements['O'] -= len(transcript.nucleotide_sequence)

		coralme.util.massbalance.elements_to_formula(transcript, elements)

	def _add_or_update_demand_reaction(self, transcript):
		"""
		This is in case the TU makes multiple products and one needs a sink.
		If the demand reaction is used, it means the RNA biomass doesn't count
		toward the overall biomass constraint

		Parameters
		----------
		transcript : :class:`coralme.core.component.TranscribedGene`
			Instance of gene having its demand reaction updated/added

		"""
		metabolites = self._model.metabolites
		demand_reaction_id = 'DM_' + transcript.id
		if demand_reaction_id not in self._model.reactions:
			demand_reaction = MEReaction(demand_reaction_id)
			self._model.add_reaction(demand_reaction)
			demand_reaction.add_metabolites({transcript.id: -1})
		else:
			demand_reaction = self._model.reactions.get_by_id(demand_reaction_id)

		mass_in_kda = transcript.formula_weight / 1000.
		# Add biomass drain for each demand reaction
		if transcript.RNA_type == 'tRNA':
			demand_reaction.add_metabolites({metabolites.tRNA_biomass: -mass_in_kda}, combine = False)
		elif transcript.RNA_type == 'rRNA':
			demand_reaction.add_metabolites({metabolites.rRNA_biomass: -mass_in_kda}, combine = False)
		elif transcript.RNA_type == 'ncRNA':
			demand_reaction.add_metabolites({metabolites.ncRNA_biomass: -mass_in_kda}, combine = False)
		elif transcript.RNA_type == 'mRNA':
			demand_reaction.add_metabolites({metabolites.mRNA_biomass: -mass_in_kda}, combine = False)
		elif transcript.RNA_type == 'tmRNA':
			demand_reaction.add_metabolites({metabolites.tmRNA_biomass: -mass_in_kda}, combine = False)
		else:
			logging.warning('Gene locus ID has an invalid RNA type (Valid types are mRNA, rRNA, tRNA, ncRNA, and tmRNA).')

	def update(self, verbose = True):
		"""
		Creates reaction using the associated transcription data and adds
		chemical formula to RNA products

		This function adds the following components to the reaction
		stoichiometry (using 'data' as shorthand for
		:class:`coralme.core.processdata.TranscriptionData`):

		1) RNA_polymerase from data.RNA_polymerase w/ coupling
		   coefficient (if present)

		2) RNA products defined in data.RNA_products

		3) Nucleotide reactants defined in data.nucleotide_counts

		4) If tRNA or rRNA contained in data.RNA_types, excised base products

		5) Metabolites + enzymes w/ coupling coefficients defined in
		   data.subreactions (if present)

		6) Biomass :class:`coralme.core.component.Constraint` corresponding to
		   data.RNA_products and their associated masses

		7) Demand reactions for each transcript product of this reaction

		Parameters
		----------
		verbose : bool
			Prints when new metabolites are added to the model when executing
			update()

		"""
		self.clear_metabolites()

		metabolites = self._model.metabolites
		stoichiometry = defaultdict(int)

		tu_id = self.transcription_data.id
		tu_length = len(self.transcription_data.nucleotide_sequence)
		rna_polymerase = self.transcription_data.RNA_polymerase

		# Set Parameters
		kt = self._model.global_info['kt']
		r0 = self._model.global_info['r0']
		m_rr = self._model.global_info['m_rr']
		f_rrna = self._model.global_info['f_rRNA']
		m_aa = self._model.global_info['m_aa']

		c_ribo = m_rr / f_rrna / m_aa

		try:
			rnap = self._model.metabolites.get_by_id(rna_polymerase)
		except KeyError:
			if verbose:
				logging.warning('The \'RNA_polymerase\' Component not found for {:s}.'.format(tu_id))
		else:
			num = self._model.mu * c_ribo * kt
			den = self._model.mu + kt * r0
			k_rnap = (num / (den)) * 3  # (3*k_ribo) hr-1
			coupling = -tu_length * self._model.mu / k_rnap
			stoichiometry[rnap.id] = coupling

		# All genes in TU must be added to model prior to creating
		# transcription reaction
		for transcript_id in self.transcription_data.RNA_products:
			if transcript_id not in metabolites:
				raise UserWarning('Transcript \'{:s}\' not found in the ME-Model.'.format(transcript_id))
			else:
				transcript = self._model.metabolites.get_by_id(transcript_id)
			stoichiometry[transcript.id] += 1

			try:
				self._add_formula_to_transcript(transcript)
			except:
				logging.warning('Problem adding formula to transcript \'{:s}\'.'.format(transcript.id))

		# Add modifications and subreactions to reaction stoichiometry
		stoichiometry = self.add_subreactions(tu_id, stoichiometry)

		base_counts = self.transcription_data.nucleotide_count
		for base, count in base_counts.items():
			stoichiometry[base] -= count

		for base, count in self.transcription_data.excised_bases.items():
			stoichiometry[base] += count
			stoichiometry['h2o_c'] -= count
			stoichiometry['h_c'] += count

		stoichiometry['ppi_c'] += tu_length

		new_stoich = self.get_components_from_ids(stoichiometry, verbose = verbose, default_type = coralme.core.component.TranscribedGene)

		self.add_metabolites(new_stoich, combine = False)

		# add biomass constraints for RNA products
		trna_mass = rrna_mass = ncrna_mass = mrna_mass = tmrna_mass = 0.

		for met, v in new_stoich.items():
			if v < 0 or not hasattr(met, 'RNA_type'):
				continue
			if met.RNA_type == 'tRNA':
				trna_mass += met.formula_weight / 1000.  # kDa
			if met.RNA_type == 'rRNA':
				rrna_mass += met.formula_weight / 1000.  # kDa
			if met.RNA_type == 'ncRNA':
				ncrna_mass += met.formula_weight / 1000.  # kDa
			if met.RNA_type == 'mRNA':
				mrna_mass += met.formula_weight / 1000.  # kDa
			if met.RNA_type == 'tmRNA':
				tmrna_mass += met.formula_weight / 1000.  # kDa

			# Add demand of each transcript
			self._add_or_update_demand_reaction(met)

		# Add the appropriate biomass constraints for each RNA contained in
		# the transcription unit
		if trna_mass > 0:
			self.add_metabolites({metabolites.tRNA_biomass: trna_mass}, combine = False)
		if rrna_mass > 0:
			self.add_metabolites({metabolites.rRNA_biomass: rrna_mass}, combine = False)
		if ncrna_mass > 0:
			self.add_metabolites({metabolites.ncRNA_biomass: ncrna_mass}, combine = False)
		if mrna_mass > 0:
			self.add_metabolites({metabolites.mRNA_biomass: mrna_mass}, combine = False)
		if tmrna_mass > 0:
			self.add_metabolites({metabolites.tmRNA_biomass: tmrna_mass}, combine = False)

class GenericFormationReaction(MEReaction):
	"""
	Some components in an ME-model can perform exactly the same function. To
	handle this, GenericFormationReactions are used to create generic forms
	of these components.

	Parameters
	----------
	id : str
		Identifier of the generic formation reaction. As a best practice, this
		ID should be prefixed with
		'metabolite_id + _to_ + generic_metabolite_id'
	"""

	def __init__(self, id = None):
		MEReaction.__init__(self, id)

	# TODO GenericFormation should have update function to account for changes
	# in component_list of GenericData

class TranslationReaction(MEReaction):
	"""Reaction class for the translation of a TranscribedGene to a
	TranslatedGene

	Parameters
	----------
	id : str
		Identifier of the translation reaction. As a best practice, this ID
		should be prefixed with 'translation + _'

	"""

	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._translation_data = None

	@property
	def translation_data(self):
		"""
		Get and set the :class:`cobra.core.processdata.TranslationData` that
		defines the translation of the gene. Can be set with instance of
		TranslationData or with its id.

		Returns
		-------
		:class:`cobra.core.processdata.TranslationData`

		"""
		return self._translation_data

	@translation_data.setter
	def translation_data(self, process_data):
		if isinstance(process_data, str):
			process_data = self._model.process_data.get_by_id(process_data)
		self._translation_data = process_data
		process_data._parent_reactions.add(self.id)

	def _add_formula_to_protein(self, translation_data, protein):
		"""
		Adds formula to protein based on amino acid sequence and subreactions

		Some subreactions modify the composition of the protein, therefore
		this must be accounted for.

		Water is subtracted from the formula to with a multiplier of
		len(amino_acid_sequence) - 1 to account for the condensation
		reactions that occur during amino acid polymerization.

		Parameters
		----------
		translation_data : :class:`cobra.core.processdata.TranslationData`
			This is required to subtract elements removed/added to protein
			when applying reaction defined in subreaction

		protein : :class:`cobra.core.processdata.TranslationData`
			Protein product that needs a chemical formula

		"""
		elements = defaultdict(int)
		aa_count = self.translation_data.amino_acid_count
		for aa_name, value in aa_count.items():
			aa_obj = self._model.metabolites.get_by_id(aa_name)
			for e, n in aa_obj.elements.items():
				elements[e] += n * value

		elements = coralme.util.massbalance.get_elements_from_process_data(self, translation_data, elements)

		# subtract water from composition
		protein_length = len(translation_data.amino_acid_sequence)
		elements['H'] -= (protein_length - 1) * 2
		elements['O'] -= (protein_length - 1)

		coralme.util.massbalance.elements_to_formula(protein, elements)

	def update(self, verbose = True):
		"""
		Creates reaction using the associated translation data and adds
		chemical formula to protein product

		This function adds the following components to the reaction
		stoichiometry (using 'data' as shorthand for
		:class:`coralme.core.processdata.TranslationData`):

		1) Amino acids defined in data.amino_acid_sequence. Subtracting water
		   to account for condensation reactions during polymerization

		2) Ribosome w/ translation coupling coefficient (if present)

		3) mRNA defined in data.mRNA w/ translation coupling coefficient

		4) mRNA + nucleotides + hydrolysis ATP cost w/ degradation coupling
		   coefficient (if kdeg (defined in model.global_info) > 0)

		5) RNA_degradosome w/ degradation coupling coefficient (if present and
		   kdeg > 0)

		6) Protein product defined in data.protein

		7) Subreactions defined in data.subreactions

		8) protein_biomass :class:`coralme.core.component.Constraint`
		   corresponding to the protein product's mass

		9) Subtract mRNA_biomass :class:`coralme.core.component.Constraint`
		   defined by mRNA degradation coupling coefficinet (if kdeg > 0)

		Parameters
		----------
		verbose : bool
			Prints when new metabolites are added to the model when executing
			update()

		"""
		self.clear_metabolites()

		model = self._model
		metabolites = self._model.metabolites
		new_stoichiometry = defaultdict(int)

		translation_data = self.translation_data
		protein_id = translation_data.protein
		mrna_id = translation_data.mRNA
		protein_length = len(translation_data.amino_acid_sequence)
		nucleotide_sequence = translation_data.nucleotide_sequence

		# Set Parameters
		kt = self._model.global_info['kt']
		k_deg = self._model.global_info['k_deg']
		r0 = self._model.global_info['r0']
		m_rr = self._model.global_info['m_rr']
		f_rrna = self._model.global_info['f_rRNA']
		m_aa = self._model.global_info['m_aa']

		m_nt = self._model.global_info['m_nt']
		f_mrna = self._model.global_info['f_mRNA']

		c_ribo = m_rr / f_rrna / m_aa
		c_mrna = m_nt / f_mrna / m_aa

		ribosome_id = self._model.global_info['ribosome_id']
		degradosome_id = self._model.global_info['degradosome_id']
		trna_misacylation = self._model.global_info['trna_misacylation']

		# -----------------Add Amino Acids----------------------------------
		for aa, value in translation_data.amino_acid_count.items():
			if aa.replace('__L_c', '') in [ x.lower() for x in trna_misacylation.keys() ]:
				new_stoichiometry[aa] = 0
				aa = trna_misacylation[aa.replace('__L_c', '').capitalize()] + '__L_c'
				aa = aa[0].lower() + aa[1:]
				new_stoichiometry[aa] -= value
				new_stoichiometry['h2o_c'] += value
				continue

			new_stoichiometry[aa] -= value
			new_stoichiometry['h2o_c'] += value

		# Length protein - 1 dehydration reactions
		new_stoichiometry['h2o_c'] -= 1.

		# -----------------Add Ribosome Coupling----------------------------
		try:
			ribosome = metabolites.get_by_id(ribosome_id)
		except KeyError:
			if verbose:
				logging.warning('The \'{:s}\' component not found in the ME-Model. Coupling coefficient not added to \'{:s}\'.'.format(ribosome_id, protein_id))
		else:
			num = self._model.mu * c_ribo * kt
			den = self._model.mu + kt * r0
			k_ribo = num / (den)  # in hr-1
			coupling = -protein_length * self._model.mu / k_ribo
			new_stoichiometry[ribosome.id] = coupling

		# -------------------Add mRNA Coupling------------------------------
		try:
			transcript = metabolites.get_by_id(mrna_id)
		except KeyError:
			# If transcript not found add to the model as the mRNA_id
			transcript = coralme.core.component.TranscribedGene(mrna_id, mrna_id, nucleotide_sequence)
			model.add_metabolites(transcript)
			logging.warning('Transcript \'{:s}\' not found in ME-Model. Added into the ME-Model.'.format(mrna_id))

		# Calculate coupling constraints for mRNA and degradation
		num = self._model.mu * c_mrna * kt
		den = self._model.mu + kt * r0
		k_mrna = num / (den) * 3.  # 3 nucleotides per AA
		rna_amount = self._model.mu / k_mrna
		deg_amount = k_deg / k_mrna

		# Add mRNA coupling to stoichiometry
		new_stoichiometry[transcript.id] = -(rna_amount + deg_amount)

		# ---------------Add Degradation Requirements -------------------------
		# Add degraded nucleotides to stoichiometry
		for nucleotide, count in transcript.nucleotide_count.items():
			new_stoichiometry[nucleotide] += count * deg_amount

		# ATP hydrolysis required for cleaving
		nucleotide_length = len(transcript.nucleotide_sequence)

		# .25 ATP required per nucleotide hydrolysis
		hydrolysis_amount = (nucleotide_length - 1) / 4. * deg_amount
		# old code; now set as a global_info
		#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'pi_c': 1, 'h_c': 1}
		atp_hydrolysis = self._model.process_data.get_by_id('atp_hydrolysis').stoichiometry
		for metabolite, value in atp_hydrolysis.items():
			new_stoichiometry[metabolite] += hydrolysis_amount * value

		# Add degradosome coupling, if known
		try:
			rna_degradosome = metabolites.get_by_id(degradosome_id)
		except KeyError:
			if verbose:
				logging.warning('The \'{:s}\' component not found in the ME-Model. Coupling coefficient not added to \'{:s}\'.'.format(degradosome_id, protein_id))
		else:
			deg_coupling = -deg_amount * self._model.mu / 65. / 3600  # keff of degradosome
			new_stoichiometry[rna_degradosome.id] = deg_coupling

		# --------------- Add Protein to Stoichiometry ------------------------
		# Added protein to model if not already included
		try:
			protein = metabolites.get_by_id(protein_id)
		except KeyError:
			protein = coralme.core.component.TranslatedGene(protein_id)
			model.add_metabolites(protein)
		new_stoichiometry[protein.id] = 1

		# ------- Convert ids to metabolites and add to model -----------------
		# add subreactions to stoichiometry
		new_stoichiometry = self.add_subreactions(self.translation_data.id, new_stoichiometry)

		# convert metabolite ids to cobra metabolites
		object_stoichiometry = self.get_components_from_ids(new_stoichiometry, verbose = verbose)
		# add metabolites to reaction
		self.add_metabolites(object_stoichiometry, combine = False)

		# -------------Update Element Dictionary and Formula-------------------
		self._add_formula_to_protein(translation_data, protein)

		# ------------------ Add biomass constraints --------------------------
		# add biomass constraint for protein translated¨
		protein_mass = protein.formula_weight / 1000.  # kDa
		self.add_metabolites({metabolites.protein_biomass: protein_mass}, combine = False)

		# RNA biomass consumed due to degradation
		mrna_mass = transcript.formula_weight / 1000.  # kDa
		self.add_metabolites({metabolites.mRNA_biomass: (-mrna_mass * deg_amount)}, combine = False)

class tRNAChargingReaction(MEReaction):
	"""
	Reaction class for the charging of a tRNA with an amino acid

	Parameters
	----------
	id : str
		Identifier for the charging reaction. As a best practice, ID should
		follow the template 'charging_tRNA + _ + <tRNA_locus> + _ + <codon>'.
		If tRNA initiates translation, <codon> should be replaced with START.

	"""
	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._tRNA_data = None

	@property
	def tRNA_data(self):
		"""
		Get and set the :class:`cobra.core.processdata.tRNAData` that
		defines the translation of the gene. Can be set with instance of
		tRNAData or with its id.

		Returns
		-------
		:class:`cobra.core.processdata.tRNAData`
		"""
		return self._tRNA_data

	@tRNA_data.setter
	def tRNA_data(self, process_data):
		if isinstance(process_data, str):
			process_data = self._model.process_data.get_by_id(process_data)
		self._tRNA_data = process_data
		process_data._parent_reactions.add(self.id)

	def update(self, verbose = True):
		"""
		Creates reaction using the associated tRNA data

		This function adds the following components to the reaction
		stoichiometry (using 'data' as shorthand for
		:class:`coralme.core.processdata.tRNAData`):

		1) Charged tRNA product following template:
		   'generic_tRNA + _ + <data.codon> + _ + <data.amino_acid>'

		2) tRNA metabolite (defined in data.RNA) w/ charging coupling
		   coefficient

		3) Charged amino acid (defined in data.amino_acid) w/ charging
		   coupling coefficient

		5) Synthetase (defined in data.synthetase) w/ synthetase coupling
		   coefficient found, in part, using data.synthetase_keff

		6) Post transcriptional modifications defined in data.subreactions

		Parameters
		----------
		verbose : bool
			Prints when new metabolites are added to the model when executing
			update()

		"""
		self.clear_metabolites()
		new_stoichiometry = defaultdict(float)
		data = self.tRNA_data

		# set tRNA coupling parameters
		m_trna = self._model.global_info['m_tRNA']
		m_aa = self._model.global_info['m_aa']
		f_trna = self._model.global_info['f_tRNA']
		kt = self._model.global_info['kt']  # hr-1
		r0 = self._model.global_info['r0']
		c_trna = m_trna / m_aa / f_trna

		# The meaning of a generic tRNA is described in the
		# TranslationReaction comments
		generic_trna = 'generic_tRNA_' + data.codon + '_' + data.amino_acid
		new_stoichiometry[generic_trna] = 1

		# Compute tRNA (and amino acid) coupling and add to stoichiometry
		num = c_trna * kt * self._model.mu
		den = self._model.mu + r0 * kt
		trna_keff = num / (den)  # per hr
		trna_amount = self._model.mu / trna_keff
		new_stoichiometry[data.RNA] = -trna_amount
		new_stoichiometry[data.amino_acid] = -trna_amount

		# Add synthetase coupling and enzyme, if known
		synthetase_amount = self._model.mu / data.synthetase_keff / 3600. * (1 + trna_amount)
		if data.synthetase is not None:
			new_stoichiometry[data.synthetase] = -synthetase_amount

		# Add tRNA modifications to stoichiometry
		new_stoichiometry = self.add_subreactions(self.tRNA_data.id, new_stoichiometry, scale = trna_amount)

		# Convert component ids to cobra metabolites
		object_stoichiometry = self.get_components_from_ids(new_stoichiometry, verbose = verbose)

		# Replace reaction stoichiometry with updated stoichiometry
		self.add_metabolites(object_stoichiometry)

class SummaryVariable(MEReaction):
	"""
	SummaryVariables are reactions that impose global constraints on the model.

	The primary example of this is the biomass_dilution SummaryVariable which
	forces the rate of biomass production of macromolecules, etc., to be equal
	to the rate of their dilution to daughter cells during growth.

	Parameters
	----------
	id : str
		Identifier of the SummaryVariable

	"""
	def __init__(self, id = None):
		MEReaction.__init__(self, id)
		self._objective_coefficient = 0.

	# WARNING: included to add the DNAPol into the DNA_replication SummaryVariable
	def update(self, verbose = True):
		if self.id == 'DNA_replication':
			model = self._model
			metabolites = self._model.metabolites
			new_stoichiometry = defaultdict(int)

			dnapol_id = self._model.global_info['dnapol_id']

			# -----------------Add DNAP Coupling----------------------------
			try:
				dnap = metabolites.get_by_id(dnapol_id)
			except KeyError:
				if verbose:
					logging.warning('The \'{:s}\' component not found in the ME-Model. Coupling coefficient not added to \'{:s}\'.'.format(dnapol_id, self.id))
			else:
				#num = self._model.mu * c_ribo * kt
				#den = self._model.mu + kt * r0
				#k_ribo = num / (den)  # in hr-1
				#coupling = -protein_length * self._model.mu / k_ribo
				new_stoichiometry[dnap.id] = -1e-6 # must be low

			# Convert component ids to cobra metabolites
			object_stoichiometry = self.get_components_from_ids(new_stoichiometry, verbose = verbose)

			# Replace reaction stoichiometry with updated stoichiometry
			self.add_metabolites(object_stoichiometry, combine = True)