import Bio
import copy
import pandas
import sympy

import cobra
import coralme

import logging
log = logging.getLogger(__name__)

import collections

class ProcessData(object):
	"""Generic class for storing information about a process

	This class essentially acts as a database that contains all of the
	relevant information needed to construct a particular reaction. For
	example, to construct a transcription reaction, following information must
	be accessed in some way:

	 - nucleotide sequence of the transcription unit
	 - RNA_polymerase (w/ sigma factor)
	 - RNAs transcribed from transcription unit
	 - other processes involved in transcription of RNAs (splicing, etc.)

	ME-model reactions are built from information in these objects.

	Parameters
	----------
	id : str
		Identifier of the ProcessData instance.

	model : :class:`coralme.core.model.MEModel`
		ME-model that the ProcessData is associated with

	"""

	def __init__(self, id, model):
		self.id = id
		self._model = model
		# parents need to be updated every time the process is updated
		# a parent must have an update method
		self._parent_reactions = set()
		model.process_data.append(self)

	# WARNING: MODIFIED FUNCTION FROM COBRAPY
	def copy(self) -> "ProcessData":
		"""Copy a ProcessData.
		"""
		# no references to model when copying
		model = self._model
		self._model = None
		# now we can copy
		new_processdata = copy.deepcopy(self)
		# restore the references
		self._model = model
		if not hasattr(self, 'id'):
			print(new_processdata)
		return new_processdata

	@property
	def model(self):
		"""
		Get the ME-model the process data is associated with

		Returns
		-------
		:class:`coralme.core.model.MEModel
			ME-model that uses this process data
		"""
		return self._model

	@property
	def parent_reactions(self):
		"""
		Get reactions that the ProcessData instance is used to construct.

		Returns
		-------
		set
			Parent reactions of ProcessData

		"""
		return {self._model.reactions.get_by_id(i) for i in self._parent_reactions}

	def update_parent_reactions(self):
		"""

		Executes the update() function for all reactions that the ProcessData
		instance is used to construct.

		"""
		reactions = self._model.reactions
		for i in self._parent_reactions:
			reactions.get_by_id(i).update()

	def __repr__(self):
		return '<{:s} {:s} at 0x{:x}>'.format(self.__class__.__name__, self.id, id(self))

class StoichiometricData(ProcessData):
	"""Encodes the stoichiometry for a metabolic reaction.

	StoichiometricData defines the metabolite stoichiometry and upper/lower
	bounds of metabolic reaction

	Parameters
	----------
	id : str
		Identifier of the metabolic reaction. Should be identical to the
		M-model reactions in most cases.

	model : :class:`coralme.core.model.MEModel`
		ME-model that the StoichiometricData is associated with

	Attributes
	----------
	_stoichiometry : dict
		Dictionary of {metabolite_id: stoichiometry} for reaction

	subreactions : :class:`collections.DefaultDict(int)`
		Cases where multiple enzymes (often carriers ie. Acyl Carrier Protein)
		are involved in a metabolic reactions.

	upper_bound : int
		Upper reaction bound of metabolic reaction. Should be identical to the
		M-model reactions in most cases.

	lower_bound : int
		Lower reaction bound of metabolic reaction. Should be identical to the
		M-model reactions in most cases.
	"""
	def __init__(self, id, model):
		ProcessData.__init__(self, id, model)
		self._stoichiometry = {}
		self.subreactions = collections.defaultdict(int)
		self.lower_bound = 0.
		self.upper_bound = 1000.

	@property
	def stoichiometry(self):
		"""
		Get or set metabolite stoichiometry for reaction.

		Returns
		-------
		dict
			Dictionary of {metabolite_id: stoichiometry}
		"""
		return self._stoichiometry

	@stoichiometry.setter
	def stoichiometry(self, value):
		if not isinstance(value, dict):
			raise TypeError('Stoichiometry must be a dict, not ({:s})'.format(type(value)))
		for k in value:
			if not isinstance(k, str):
				raise TypeError('Stoichiometry keys must be strings, not \'{:s}\'.'.format(type(k)))
		self._stoichiometry = value

class SubreactionData(ProcessData):
	"""
	Parameters
	----------
	id : str
		Identifier of the subreaction data. As a best practice, if the
		subreaction data details a modification, the ID should be prefixed
		with 'mod + _'

	model : :class:`coralme.core.model.MEModel`
		ME-model that the SubreactionData is associated with

	Attributes
	----------
	enzyme : list or str or None
		List of :attr:`coralme.core.component.Complex.id` s for enzymes that
		catalyze this process

		or

		String of single :attr:`coralme.core.component.Complex.id` for enzyme
		that catalyzes this process

	keff : float
		Effective turnover rate of enzyme(s) in subreaction process

	_element_contribution : dict
		If subreaction adds a chemical moiety to a macromolecules via a
		modification or other means, net element contribution of the
		modification process should be accounted for. This can be used to
		mass balance check each of the individual processes.

		Dictionary of {element: net_number_of_contributions}

	"""
	def __init__(self, id, model):
		ProcessData.__init__(self, id, model)
		self.stoichiometry = {}
		self.enzyme = None
		self._coupling_coefficient_subreaction = sympy.Mul(self._model.mu, sympy.Rational('1/3600'), model.symbols['k^default_cat']**-1, evaluate = False)
		self._element_contribution = {}

	# Backward compatibility
	@property
	def keff(self):
		"""
		returns the keff value, not the coupling coefficient, in per second
		"""
		value = self._model.mu * sympy.Rational('1/3600') / self._coupling_coefficient_subreaction
		try:
			return float(value)
		except:
			return float(value.xreplace(self._model.global_info['default_parameters']))

	# Backward compatibility
	@keff.setter
	def keff(self, value):
		"""
		value is the keff in per second, not the coupling coefficient in per hour
		this returns the coupling coefficient as growth rate divided by the keff, in per hour
		"""
		self._coupling_coefficient_subreaction = sympy.Mul(self._model.mu, sympy.Rational('1/3600'), value**-1, evaluate = False)

	@property
	def coupling_coefficient_subreaction(self):
		"""
		returns the coupling coefficient, not the keff value, in per hour
		"""
		return self._coupling_coefficient_subreaction

	@coupling_coefficient_subreaction.setter
	def coupling_coefficient_subreaction(self, value):
		"""
		value is the keff in per second, not the coupling coefficient in per hour
		this returns the coupling coefficient as growth rate divided by the keff, in per hour
		"""
		self._coupling_coefficient_subreaction = sympy.Mul(self._model.mu, sympy.Rational('1/3600'), value**-1, evaluate = False)

	@property
	def element_contribution(self):
		"""
		Get net contribution of elements from subreaction process to
		macromolecule

		If subreaction adds a chemical moiety to a macromolecules via a
		modification or other means, net element contribution of the
		modification process should be accounted for. This can be used to
		mass balance check each of the individual processes.

		Returns
		-------
		dict
			Dictionary of {element: net_number_of_contributions}

		"""
		if self._element_contribution:
			return self._element_contribution
		else:
			contribution = { k:v for k,v in self.calculate_element_contribution().items() if v }

		# Return 'trivial' cases (only one modifying metabolite in the
		# reactants and no products) without warning
		if len(self.stoichiometry) == 1 and list(self.stoichiometry.values())[0] < 0:
			return self.calculate_element_contribution()
		elif contribution:
			logging.warning('No element contribution input for SubReaction \'{:s}\'. Calculating based on stoichiometry instead.'.format(self.id))
			return self.calculate_element_contribution()
		else:
			return {}

	@element_contribution.setter
	def element_contribution(self, value):
		if not isinstance(value, dict):
			raise TypeError('Elemental_contribution must be a dictionary, not \'{:s}\'.'.format(type(value)))
		self._element_contribution = value

	def calculate_element_contribution(self):
		"""
		Calculate net contribution of chemical elements based on the
		stoichiometry of the subreaction data

		Returns
		-------
		dict
			Dictionary of {element: net_number_of_contributions}

		"""
		elements = collections.defaultdict(int)
		for met, coefficient in self.stoichiometry.items():
			if self._model.metabolites.has_id(met):
				met_obj = self._model.metabolites.get_by_id(met)
			else:
				logging.warning('The metabolite \'{:s}\' must exist in the ME-model to calculate the element contribution.'.format(met))
				continue

			# elements lost in conversion are added to complex, protein, etc.
			if not met_obj.elements and not isinstance(met_obj, coralme.core.component.GenerictRNA):
				logging.warning('Metabolite \'{:s}\' does not have a formula. If it is a \'Complex\', its formula will be determined from amino acid composition and prosthetic groups stoichiometry. Otherwise, please add it to the M-model.'.format(met_obj.id))

			for e, n in met_obj.elements.items():
				elements[e] -= n * coefficient

		return elements

	def calculate_biomass_contribution(self):
		"""
		Calculate net biomass increase/decrease as a result of the subreaction
		process.

		If subreaction adds a chemical moiety to a macromolecules via a
		modification or other means, the biomass contribution of the
		modification process should be accounted for and ultimately included
		in the reaction it is involved in.

		Returns
		-------
		float
			Mass of moiety transferred to macromolecule by subreaction

		"""
		elements = self.element_contribution

		# Create temporary metabolite for calculating formula weight
		tmp_met = cobra.Metabolite('mass')
		coralme.util.massbalance.elements_to_formula(tmp_met, elements)

		return tmp_met.formula_weight

	def get_complex_data(self):
		"""
		Get the complex data that the subreaction is involved in

		Yields
		------
		:class:`coralme.core.processdata.ComplexData`
			ComplexData that subreaction is involved in
		"""
		for i in self._model.complex_data:
			if self.id in i.subreactions:
				yield i

	def get_all_usages(self):
		"""
		Get all process data that the subreaction is involved in

		Yields
		------
		:class:`coralme.core.processdata.ProcessData`
			ProcessData that subreaction is involved in
		"""
		for i in self._model.process_data:
			if hasattr(i, 'subreactions') and self.id in i.subreactions:
				yield i

class ComplexData(ProcessData):
	"""Contains all information associated with the formation of an
	functional enzyme complex.

	This can include any enzyme complex modifications required for the enzyme
	to become active.

	Parameters
	----------
	id : str
		Identifier of the complex data. As a best practice, this should
		typically use the same ID as the complex being formed. In cases with
		multiple ways to form complex '_ + alt' or similar suffixes can be
		used.

	model : :class:`coralme.core.model.MEModel`
		ME-model that the ComplexData is associated with

	Attributes
	----------

	stoichiometry : :class:`collections.DefaultDict(int)`
		Dictionary containing {protein_id: count} for all protein subunits
		comprising enzyme complex

	subreactions : dict
		Dictionary of {subreaction_data_id: count} for all complex formation
		subreactions/modifications. This can include cofactor/prosthetic group
		binding or enzyme side group addition.

	"""

	def __init__(self, id, model):
		ProcessData.__init__(self, id, model)
		# {Component.id: stoichiometry}
		self.stoichiometry = collections.defaultdict(int)
		# {SubreactionData.id : number}
		# Forming some metacomplexes occur in multiple steps
		self.subreactions = {}
		self._complex_id = None  # assumed to be the same as id if None

	@property
	def formation(self):
		"""Get the formation reaction object

		Returns
		-------
		:class:`coralme.core.reaction.ComplexFormation`
			Complex formation reaction detailed in ComplexData
		"""
		try:
			return self._model.reactions.get_by_id('formation_' + self.id)
		except KeyError:
			return None

	@property
	def complex(self):
		"""
		Get complex metabolite object

		Returns
		-------
		:class:`coralme.core.component.Complex`
			Instance of complex metabolite that ComplexData is used to
			synthesize
		"""
		return self._model.metabolites.get_by_id(self.complex_id)

	@property
	def complex_id(self):
		"""
		Get  and set complex ID for product of complex formation reaction

		There are cases where multiple equivalent processes can result in
		the same final complex. This allows the equivalent final complex
		complex_id to be queried. This only needs set in the above case

		Returns
		-------
		str
			ID of complex that ComplexData is used to synthesize
		"""

		return self.id if self._complex_id is None else self._complex_id

	@complex_id.setter
	def complex_id(self, value):
		self._complex_id = None if value == self.id else value

	def create_complex_formation(self, verbose=True):
		"""creates a complex formation reaction

		This assumes none exists already. Will create a reaction (prefixed by
		'formation') which forms the complex

		Parameters
		----------
		verbose : bool
			If True, print if a metabolite is added to model during update

		"""
		formation_id = 'formation_' + self.id
		if formation_id in self._model.reactions:
			raise ValueError('Reaction \'{:s}\' already in the ME-model.'.format(formation_id))
		formation = coralme.core.reaction.ComplexFormation(formation_id)
		formation.complex_data_id = self.id
		formation._complex_id = self.complex_id
		self._model.add_reactions([formation])
		formation.update(verbose = verbose)

class TranscriptionData(ProcessData):
	"""
	Class for storing information needed to define a transcription reaction

	Parameters
	----------
	id : str
		Identifier of the transcription unit, typically beginning with 'TU'

	model : :class:`coralme.core.model.MEModel`
		ME-model that the TranscriptionData is associated with

	Attributes
	----------

	nucleotide_sequence : str
		String of base pair abbreviations for nucleotides contained in the
		transcription unit

	RNA_products : set
		IDs of :class:`coralme.core.component.TranscribedGene` that the
		transcription unit encodes. Each member should be prefixed with
		'RNA + _'

	RNA_polymerase : str
		ID of the :class:`coralme.core.component.RNAP` that transcribes the
		transcription unit. Different IDs are used for different sigma factors

	subreactions : :class:`collections.DefaultDict(int)`
		Dictionary of
		{:class:`coralme.core.processdata.SubreactionData` ID: num_usages}
		required for the transcription unit to be transcribed

	"""
	def __init__(self, id, model, nucleotide_sequence, rnap, rna_products, organelle):
		ProcessData.__init__(self, id, model)
		self.nucleotide_sequence = nucleotide_sequence
		self.RNA_products = rna_products
		self.original_RNA_products = rna_products
		self.RNA_polymerase = rnap
		self.organelle = organelle

		self._subreactions = collections.defaultdict(int)
		self._coupling_coefficient_rnapol = sympy.Mul(len(nucleotide_sequence), model.symbols['v_rnap'], evaluate = False)

	@property
	def coupling_coefficient_rnapol(self):
		return self._coupling_coefficient_rnapol

	@coupling_coefficient_rnapol.setter
	def coupling_coefficient_rnapol(self, value):
		self._coupling_coefficient_rnapol = value

	@property
	def n_cuts(self):
		# Number of cuts depends on the type of the RNAs in the TU
		return len([ x for x in self.RNA_types if x in ['rRNA', 'tRNA']]) * 2.

	@property
	def n_excised(self):
		# Number of excised bases depends on the type of the RNAs in the TU
		if set(self.RNA_types) == {'mRNA'}:
			return 0
		else:
			return sum(self.excised_bases.values())

	@property
	def n_overlapping(self):
		if self.id == 'RNA_dummy' or len(self.RNA_products) == 0:
			return 0

		import pyranges

		ranges = []
		for rna in self.RNA_products:
			data = self.model.metabolites.get_by_id(rna)
			left_pos = data.left_pos[0].replace('>', '').replace('<', '')
			right_pos = data.right_pos[0].replace('>', '').replace('<', '')
			ranges.append(['X', left_pos, right_pos, data.strand])

		df = pandas.DataFrame(ranges, columns = ['Chromosome', 'Start', 'End', 'Strand'])
		ranges = pyranges.PyRanges(df)

		# add overlapping ranges to df
		res = ranges.intersect(ranges, strandedness = 'same').df

		# remove original ranges
		# WARNING: What does happen if a gene overlaps completely another?
		tmp = pandas.merge(res, ranges.df, how = 'outer', indicator = True)
		tmp = tmp[tmp['_merge'] == 'left_only'].drop_duplicates()

		# return total length of the overlaps
		return abs(tmp['Start'] - tmp['End']).sum().sum()

	@property
	def subreactions(self):
		data = self._subreactions

		# Number of cuts and excised bases depend on the type of the RNAs in the TU
		if 'rRNA' not in set(self.RNA_types) and 'tRNA' not in set(self.RNA_types):
			return data

		n_overlapping = self.n_overlapping
		n_excised = self.n_excised
		n_cuts = self.n_cuts

		# WARNING: Because first 'if', n_cuts cannot be zero
		#if n_excised == 0 or (n_excised + n_overlapping) == 0 or n_cuts == 0:
			#return data
		#if n_excised == 0:
			#n_cuts = 0

		rna_types = list(self.RNA_types)
		n_trna = rna_types.count('tRNA')

		if 'rRNA' in set(rna_types):
			data['rRNA_containing_excision'] = n_cuts
		elif n_trna == 1:
			data['monocistronic_excision'] = n_cuts
		elif n_trna > 1:
			data['polycistronic_wout_rRNA_excision'] = n_cuts
		else: # only applies to rnpB (RNase P catalytic RNA component)
			data['monocistronic_excision'] = n_cuts

		# The non functional RNA segments need degraded back to nucleotides
		# TODO check if RNA_degradation requirement is per nucleotide
		data['RNA_degradation_machine'] = n_cuts
		data['RNA_degradation_atp_requirement'] = n_excised + n_overlapping

		return { k:v for k,v in data.items() if v != 0 }

	@property
	def nucleotide_count(self):
		"""
		Get count of each nucleotide contained in the nucleotide sequence

		Returns
		-------
		dict
			{nucleotide_id: number_of_occurences}

		"""
		#return { coralme.util.dogma.transcription_table[i]: self.nucleotide_sequence.count(i) for i in ['A', 'T', 'G', 'C'] }
		#return { coralme.util.dogma.transcription_table[k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }
		if self.organelle is None:
			if self._model.global_info['domain'].lower() in ['prokaryote', 'bacteria']:
				return { coralme.util.dogma.transcription_table['c'][k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }
			if self._model.global_info['domain'].lower() in ['eukarya', 'eukaryote']:
				return { coralme.util.dogma.transcription_table['n'][k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }
			#return { coralme.util.dogma.transcription_table['n'][k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }
		elif self.organelle.lower() in ['mitochondria', 'mitochondrion']:
			return { coralme.util.dogma.transcription_table['m'][k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }
		elif self.organelle.lower() in ['chloroplast', 'plastid']:
			return { coralme.util.dogma.transcription_table['h'][k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }
		else:
			logging.warning('The \'organelle\' property of the feature \'{:s}\' is not \'mitochondria\' or \'chloroplast\'.'.format(self.id))
			return { coralme.util.dogma.transcription_table['n'][k]:v for k,v in collections.Counter(self.nucleotide_sequence).items() }

	@property
	def RNA_types(self):
		"""
		Get generator consisting of the RNA type for each RNA product

		Yields
		------
		str
			(mRNA, tRNA, rRNA, ncRNA)
		"""
		for rna in self.RNA_products:
			rna_type = self._model.metabolites.get_by_id(rna).RNA_type
			if rna_type:
				yield rna_type

	@property
	def excised_bases(self):
		"""
		Get count of bases that are excised during transcription

		If a stable RNA (e.g. tRNA or rRNA) is coded for in the transcription
		unit, the transcript must be spliced in order for these to function.

		This determines whether the transcription unit requires splicing and,
		if so, returns the count of nucleotides within the transcription unit
		that are not accounted for in the RNA products, thus identifying the
		appropriate introns nucleotides.

		Returns
		-------
		dict
			{nucleotide_monophosphate_id: number_excised}

			i.e. {'amp_c': 10, 'gmp_c': 11, 'ump_c': 9, 'cmp_c': 11}

		"""
		rna_types = set(self.RNA_types)

		# Skip if TU does not have any annotated RNA Products
		if len(rna_types) == 0:
			return {'amp_c': 0, 'gmp_c': 0, 'ump_c': 0, 'cmp_c': 0}

		# Skip if TU only codes for mRNA
		# WARNING: The GenBank can contain other types of RNAs that break the condition of only mRNAs in the TU
		if rna_types == {'mRNA'}:
			return {'amp_c': 0, 'gmp_c': 0, 'ump_c': 0, 'cmp_c': 0}

		# WARNING: Features in the TU can overlap, thus this calculation must be corrected
		# Get dictionary of all nucleotide counts for TU
		seq = self.nucleotide_sequence
		#counts = { i: seq.count(i) for i in ('A', 'T', 'G', 'C') }
		counts = collections.Counter(seq)

		# Subtract bases contained in RNA_product from dictionary
		metabolites = self._model.metabolites
		for product_id in self.RNA_products:
			gene_seq = metabolites.get_by_id(product_id).nucleotide_sequence
			#for b in ('A', 'T', 'G', 'C'):
				#counts[b] -= gene_seq.count(b)
			counts.subtract(collections.Counter(gene_seq)) # inplace

		# First base being a triphosphate will be handled by the reaction
		# producing an extra ppi during transcription. But generally, we add
		# triphosphate bases when transcribing, but excise monophosphate bases.
		#monophosphate_counts = { coralme.util.dogma.transcription_table[k].replace('tp_c', 'mp_c'):v for k,v in counts.items() }
		monophosphate_counts = { coralme.util.dogma.transcription_table['c'][k].replace('tp_c', 'mp_c'):v for k,v in counts.items() }
		return monophosphate_counts

	@property
	def codes_stable_rna(self):
		"""
		Get whether transcription unit codes for a stable RNA

		Returns
		-------
		bool
			True if tRNA or rRNA in RNA products
			False if not

		"""
		has_stable_rna = False
		for RNA in self.RNA_products:
			try:
				gene = self._model.metabolites.get_by_id(RNA)
			except KeyError:
				pass
			else:
				if gene.RNA_type in ['tRNA', 'rRNA', 'ncRNA']:
					has_stable_rna = True
		return has_stable_rna

class GenericData(ProcessData):
	"""
	Class for storing information about generic metabolites

	Parameters
	----------
	id : str
		Identifier of the generic metabolite. As a best practice, this ID
		should be prefixed with 'generic + _'

	model : :class:`coralme.core.model.MEModel`
		ME-model that the GenericData is associated with

	component_list : list
		List of metabolite ids for all metabolites that can provide
		identical functionality
	"""
	def __init__(self, id, model, component_list):
		if not id.startswith('generic_'):
			logging.warning('Best practice for generic id to start with the \'generic_\' prefix.')
		ProcessData.__init__(self, id, model)
		self.component_list = component_list

		# bypass problems with GenericData not having complex and complex_id attributes
		self._complex_id = None  # assumed to be the same as id if None

	# bypass problems with GenericData not having complex and complex_id attributes
	@property
	def complex(self):
		"""
		Get complex metabolite object

		Returns
		-------
		:class:`coralme.core.component.Complex`
			Instance of complex metabolite that ComplexData is used to
			synthesize
		"""
		return self._model.metabolites.get_by_id(self.complex_id)

	@property
	def complex_id(self):
		"""
		Get  and set complex ID for product of complex formation reaction

		There are cases where multiple equivalent processes can result in
		the same final complex. This allows the equivalent final complex
		complex_id to be queried. This only needs set in the above case

		Returns
		-------
		str
			ID of complex that ComplexData is used to synthesize
		"""

		return self.id if self._complex_id is None else self._complex_id

	@complex_id.setter
	def complex_id(self, value):
		self._complex_id = None if value == self.id else value

	def create_reactions(self):
		"""

		Adds reaction with id '<metabolite_id> + _ + to + _ + <generic_id>'
		for each metabolite in self.component_list.

		Creates generic metabolite and generic reaction, if they do not already
		exist.
		"""
		model = self._model
		try:
			generic_metabolite = model.metabolites.get_by_id(self.id)
		except KeyError:
			generic_metabolite = coralme.core.component.GenericComponent(self.id)
			model.add_metabolites([generic_metabolite])
		for c_id in [ x for x in self.component_list if x.replace('RNA_', '') not in model.global_info['knockouts'] ]:
			reaction_id = c_id + '_to_' + self.id
			try:
				reaction = model.reactions.get_by_id(reaction_id)
			except KeyError:
				reaction = coralme.core.reaction.GenericFormationReaction(reaction_id)
				model.add_reactions([reaction])
			stoic = {
				generic_metabolite: 1,
				model.metabolites.get_by_id(c_id): -1
				}
			reaction.add_metabolites(stoic, combine=False)

class TranslationData(ProcessData):
	"""
	Class for storing information about a translation reaction.

	Parameters
	----------
	id : str
		Identifier of the gene being translated, typically the locus tag

	model : :class:`coralme.core.model.MEModel`
		ME-model that the TranslationData is associated with

	mrna : str
		ID of the mRNA that is being translated

	protein : str
		ID of the protein product.

	Attributes
	----------
	mRNA : str
		ID of the mRNA that is being translated

	protein : str
		ID of the protein product.

	subreactions : :class:`collections.DefaultDict(int)`
		Dictionary of
		{:attr:`coralme.core.processdata.SubreactionData.id`: num_usages}
		required for the mRNA to be translated

	nucleotide_sequence : str
		String of base pair abbreviations for nucleotides contained in the gene
		being translated

	"""
	def __init__(self, id, model, mrna, protein, nucleotide_sequence, organelle, translation, transl_table, pseudo):
		ProcessData.__init__(self, id, model)
		self.mRNA = mrna
		self.protein = protein
		self.nucleotide_sequence = nucleotide_sequence
		self.organelle = organelle
		self.translation = translation
		self.transl_table = transl_table
		self.pseudo = pseudo
		self.notes = []

		self.subreactions = collections.defaultdict(int)
		self._coupling_coefficient_ribosome = sympy.Mul(len(translation), model.symbols['v_ribo'], evaluate = False)
		self._coupling_coefficient_rna_synthesis = self._model.symbols['rna_amount'] + self._model.symbols['deg_amount']
		self._coupling_coefficient_hydrolysis = sympy.Mul((len(nucleotide_sequence) - 1) / 4., self._model.symbols['deg_amount'], evaluate = False) # deg_amount
		self._translational_efficiency = 1.

	@property
	def translational_efficiency(self):
		return self._translational_efficiency

	@translational_efficiency.setter
	def translational_efficiency(self, value):
		self._translational_efficiency = value

	@property
	def last_codon(self):
		"""
		Get the last codon contained in the mRNA sequence. This should
		correspond to the stop codon for the gene.

		Returns
		-------
		str
			Last 3 nucleotides comprising the last codon in the mRNA gene
			sequence

		"""
		return self.nucleotide_sequence[-3:].replace('T', 'U')

	@property
	def first_codon(self):
		"""
		Get the first codon contained in the mRNA sequence. This should
		correspond to the start codon for the gene.

		Returns
		-------
		str
			First 3 nucleotides comprising the first codon in the mRNA gene
			sequence

		"""
		return self.nucleotide_sequence[:+3].replace('T', 'U')

	def _itercodons(self):
		yield [i for i in self.codon_count]

	@property
	def sequence_as_codons(self):
		return [ self.nucleotide_sequence[i:i+3] for i in range(0, len(self.nucleotide_sequence), 3) ]

	@property
	def amino_acid_sequence(self):
		"""
		Get amino acid sequence from mRNA's nucleotide sequence

		Returns
		-------
		str
			Amino acid sequence

		"""
		#codons = (self.nucleotide_sequence[i: i + 3] for i in range(0, (len(self.nucleotide_sequence)), 3))
		#amino_acid_sequence = ''.join(coralme.util.dogma.codon_table[i] for i in codons)
		#amino_acid_sequence = str(Bio.Seq.Seq(self.nucleotide_sequence).translate(self._model.global_info['codon_table']))

		codons = self.sequence_as_codons # it includes the stop codon

		# translate first codon
		amino_acid_sequence = Bio.Seq.Seq('M')
		if isinstance(self.transl_table, str): # from JSON
			self.transl_table = Bio.Data.CodonTable.generic_by_id[int(list(self.transl_table)[0])]

		if codons[0] not in self.transl_table.start_codons:
			logging.warning('First codon in \'{:s}\' does not encode a start methionine. A methionine replaces the first amino acid.'.format(self.id))
			self.notes.append('Codon at position 0 does not encode a start methionine.')

		# translate rest of the sequence
		for idx, codon in enumerate(codons[1:-1]): # avoid start and stop codon
			if codon == 'TGA' and self.transl_table.id == 11:
				# Ser-tRNA is the precursor of Sec-tRNA.
				# Reaction Ser-tRNA(Sec) => Sec-tRNA(Sec) is added as a subreaction in translation reactions
				aa = 'S'
				logging.warning('Internal stop codon UGA identified in \'{:s}\' and translated into Selenocysteine tRNA precursor.'.format(self.id))
				self.notes.append('Codon at position {:d} encodes Selenocysteine.'.format(idx+1))
			elif codon in self._model.global_info['genetic_recoding'].keys():
				aa = '_' # placeholder to identify a recoded stop codon in self.amino_acid_sequence when creating TranslationReactions
				logging.warning('Internal stop codon \'{:s}\' identified in \'{:s}\' following user input.'.format(codon, self.id))
				self.notes.append('Codon at position {:d} encodes a recoded stop codon.'.format(idx+1))
			elif codon in self.transl_table.stop_codons:
				aa = '_' # placeholder to identify a recoded stop codon in self.amino_acid_sequence when creating TranslationReactions
				logging.warning('Internal stop codon \'{:s}\' identified in \'{:s}\'. Translation will not proceed. Please check if the gene is a pseudogene.'.format(codon, self.id))
				self.notes.append('Codon at position {:d} encodes an internal stop codon (\'{:s}\') not recoded.'.format(idx+1, codon))
				#break
			else:
				aa = Bio.Seq.Seq(codon).translate(self.transl_table)
			# append the translated nucleotide
			amino_acid_sequence += aa

		# last codon does not need translation unless is not a stop codon
		if codons[-1] not in self.transl_table.stop_codons:
			logging.warning('Last codon in \'{:s}\' does not encode a stop codon.'.format(self.id))
			self.notes.append('Codon at position {:d} does not encode a stop codon.'.format(idx+2))
			amino_acid_sequence += Bio.Seq.Seq(codons[-1]).translate(self.transl_table)

		self.notes = sorted(set(self.notes))

		#amino_acid_sequence = amino_acid_sequence.rstrip('*')

		if self.id != 'dummy':
			if amino_acid_sequence != self.translation.rstrip('*'):
				logging.warning('Protein sequence for \'{:s}\' from the GenBank file differs from the inferred from nucleotide sequence and translation table.'.format(self.id))

		# WARNING: This was replaced by code above
		#if '*' in amino_acid_sequence or 'U' in amino_acid_sequence: # translation of selenocysteine
			##amino_acid_sequence = amino_acid_sequence.replace('*', 'C') # Cysteine?
			#amino_acid_sequence = amino_acid_sequence.replace('*', 'S') # Ser-tRNA is the precursor of Sec-tRNA
			#amino_acid_sequence = amino_acid_sequence.replace('U', 'S') # Ser-tRNA is the precursor of Sec-tRNA

		return amino_acid_sequence

	@property
	def amino_acid_count(self):
		"""Get number of each amino acid in the translated protein

		Returns
		-------
		dict
			{amino_acid_id: number_of_occurrences}
		"""

		#aa_count = collections.defaultdict(int)
		#for i in self.amino_acid_sequence:
			#aa_count[coralme.util.dogma.amino_acids[i]] += 1
		#return aa_count

		# Set compartment
		if self.organelle is None:
			compartment = '_c'
		elif self.organelle.lower() in ['mitochondria', 'mitochondrion']:
			compartment = '_m'
		elif self.organelle.lower() in ['chloroplast', 'plastid']:
			compartment = '_h'

		#return { coralme.util.dogma.amino_acids[k] + compartment:v for k,v in collections.Counter(self.amino_acid_sequence).items() }

		precount = []
		# WARNING: sequence_as_codons includes the stop codon, therefor it is +1 longer that amino_acid_sequence
		for idx, (codon, amino_acid) in enumerate(zip(self.sequence_as_codons, self.amino_acid_sequence)):
			if amino_acid == '_' and codon == '' and self.transl_table.id == 11:
				precount.append('ser__L_c') # Ser in the precursor for Selenocysteine
			elif amino_acid == '_' and codon in self._model.global_info['genetic_recoding']:
				# TODO: set compartment of the recoded stop codon?
				precount.append(list(self._model.global_info['genetic_recoding'][codon].keys())[0])
			elif amino_acid == '_' and codon not in self._model.global_info['genetic_recoding']:
				logging.warning('Internal stop codon at position \'{:d}\' has no recoding alternative. Please set up \'genetic_recoding\' dictionary. '.format(idx))
			else:
				precount.append(coralme.util.dogma.amino_acids[amino_acid] + compartment)
		return collections.Counter(precount)

	@property
	def codon_count(self):
		"""
		Get the number of each codon contained within the gene sequence

		Returns
		-------
		dict
			{codon_sequence: number_of_occurrences}

		"""
		#codons = (self.nucleotide_sequence[i: i+3] for i in range(0, len(self.nucleotide_sequence), 3))
		#codon_count = collections.defaultdict(int)
		#for i in codons:
			#codon_count[i.replace('T', 'U')] += 1

		codons = self.sequence_as_codons
		codons = [ x.replace('T', 'U') for x in codons if len(x) == 3 ]
		codon_count = collections.Counter(codons)

		return codon_count

	@property
	def subreactions_from_sequence(self):
		"""
		Get subreactions associated with each tRNA/AA addition.

		tRNA activity is accounted for as subreactions. This returns the
		subreaction counts associated with each amino acid addition, based
		on the sequence of the mRNA.

		Returns
		-------
		dict
			{:attr:`coralme.core.processdata.SubreactionData.id`: num_usages}
		"""
		subreactions = {}

		#table = self._model.global_info['translation_table']
		# Trip first and last codon. Not translated during elongation
		codon_count = self.codon_count
		codon_count[self.first_codon] -= 1
		codon_count[self.last_codon] -= 1

		for codon, count in codon_count.items():
			if count == 0:
				continue

			codon = codon.replace('U', 'T')
			#if codon == 'TGA' and table == 11:

			#abbreviated_aa = coralme.util.dogma.codon_table[codon]
			#abbreviated_aa = Bio.Seq.Seq(codon).translate(self._model.global_info['codon_table'])
			abbreviated_aa = Bio.Seq.Seq(codon).translate(self.transl_table)
			# Filter out the compartment and stereochemistry from aa id
			aa = coralme.util.dogma.amino_acids.get(abbreviated_aa, 'STOP_CODON').split('_')[0]

			if aa == 'STOP':
				if str(abbreviated_aa) == '*' and codon in self._model.global_info.get('genetic_recoding', {}).keys():
					# perform recoding of internal stop codons
					aa = list(self._model.global_info['genetic_recoding'][codon].keys())[0].replace('__L_c', '')
					logging.warning('Recoded stop codon \'{:s}\' to \'{:s}__L_c\' in \'{:s}\'.'.format(codon, aa, self.id))
				elif str(abbreviated_aa) == '*' and codon == 'TGA' and self.transl_table.id == 11:
					logging.warning('Adding selenocysteine for \'{:s}\', following translation table {:d} (See more https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG{:d}).'.format(self.id, self.transl_table.id, self.transl_table.id))
					aa = 'sec'
				else:
					logging.warning('Internal stop codon \'{:s}\' detected in \'{:s}\'. Please review if the gene is a pseudogene.'.format(codon, self.id))

			if aa == 'STOP':
				break # do not remove break or STOP_addition_at_UAA or similar will be added

			codon = codon.replace('T', 'U')
			subreaction_id = aa + '_addition_at_' + codon
			#try:
				#self._model.process_data.get_by_id(subreaction_id)
			#except KeyError:
				#logging.warning('The tRNA SubReaction \'{:s}\' is not in the ME-model.'.format(subreaction_id))
			if self._model.process_data.has_id(subreaction_id):
				subreactions[subreaction_id] = count
			else:
				logging.warning('The tRNA SubReaction \'{:s}\' is not in the ME-model.'.format(subreaction_id))

		return subreactions

	def add_elongation_subreactions(self, elongation_subreactions=set()):
		"""
		Add all subreactions involved in translation elongation.

		This includes:

		 - tRNA activity subreactions returned with
		   :meth:`subreactions_from_sequence` which is called within this
		   function.

		 - Elongation subreactions passed into this function. These will be
		   added with a value of len(amino_acid_sequence) - 1 as these are
		   involved in each amino acid addition

		Some additional enzymatic processes are required for each amino acid
		addition during translation elongation

		Parameters
		----------
		elongation_subreactions : set
			Subreactions that are required for each amino acid addition

		"""

		for subreaction_id in elongation_subreactions:
			try:
				self._model.process_data.get_by_id(subreaction_id)
			except KeyError:
				logging.warning('Elongation SubReaction \'{:s}\' is not in ME-model. However, it can be added later.'.format(subreaction_id))
			else:
				# No elongation subreactions needed for start codon
				self.subreactions[subreaction_id] = len(self.amino_acid_sequence) - 1.

		for subreaction_id, value in self.subreactions_from_sequence.items():
			self.subreactions[subreaction_id] = value

	def add_initiation_subreactions(self, start_codons=set(), start_subreactions=set()):
		"""
		Add all subreactions involved in translation initiation.

		Parameters
		----------
		start_codons : set, optional
			Start codon sequences for the organism being modeled

		start_subreactions : set, optional
			Subreactions required to initiate translation, including the
			activity by the start tRNA

		"""
		#print(self.mRNA, type(self.mRNA), self.first_codon, type(self.first_codon))
		if self.first_codon not in start_codons:
			logging.warning('\'{:s}\' starts with \'{:s}\', which is not a start codon'.format(self.mRNA, str(self.first_codon)))

		for subreaction_id in start_subreactions:
			try:
				self._model.process_data.get_by_id(subreaction_id)
			except KeyError:
				logging.warning('Initiation SubReaction \'{:s}\' is not in the ME-model. However, it can be added later.'.format(subreaction_id))
			else:
				self.subreactions[subreaction_id] = 1.

	def add_termination_subreactions(self, translation_terminator_dict=None):
		"""
		Add all subreactions involved in translation termination.

		Parameters
		----------
		translation_terminator_dict : dict or None
			{stop_codon : enzyme_id_of_terminator_enzyme}

		"""
		if not translation_terminator_dict:
			translation_terminator_dict = {}
		last_codon = self.last_codon
		term_enzyme = translation_terminator_dict.get(last_codon, None)
		if term_enzyme:
			termination_subreaction_id = last_codon + '_' + term_enzyme + '_mediated_termination_c'
			try:
				self._model.process_data.get_by_id(termination_subreaction_id)
			except KeyError:
				logging.warning('Termination SubReaction \'{:s}\' is not in ME-model. However, it can be added later.'.format(termination_subreaction_id))
			else:
				self.subreactions[termination_subreaction_id] = 1.
		else:
			logging.warning('No termination enzyme for \'{:s}\'. Please review if the gene is a pseudogene.'.format(self.mRNA))

class tRNAData(ProcessData):
	"""
	Class for storing information about a tRNA charging reaction.

	Parameters
	----------
	id : str
		Identifier for tRNA charging process. As best practice, this should
		be follow 'tRNA + _ + <tRNA_locus> + _ + <codon>' template. If tRNA
		initiates translation, <codon> should be replaced with START.

	model : :class:`coralme.core.model.MEModel`
		ME-model that the tRNAData is associated with

	amino_acid : str
		Amino acid that the tRNA transfers to an peptide

	rna : str
		ID of the uncharged tRNA metabolite. As a best practice, this ID should
		be prefixed with 'RNA + _'

	Attributes
	----------
	subreactions : :class:`collections.DefaultDict(int)`
		Dictionary of
		{:attr:`coralme.core.processdata.SubreactionData.id`: num_usages}
		required for the tRNA to be charged

	synthetase : str
		ID of the tRNA synthetase required to charge the tRNA with an amino
		acid

	synthetase_keff : float
		Effective turnover rate of the tRNA synthetase

	"""

	def __init__(self, id, model, amino_acid, rna, codon):
		ProcessData.__init__(self, id, model)
		self.codon = codon
		self.amino_acid = amino_acid
		self.RNA = rna
		self.subreactions = collections.defaultdict(int)
		self.synthetase = None

		self._synthetase_keff = self._model.symbols['k^default_cat']

		self._coupling_coefficient_trna_keff = self._model.symbols['k_tRNA']
		self._coupling_coefficient_trna_amount = sympy.Mul(self._model.mu, self._coupling_coefficient_trna_keff**-1, evaluate = False)
		self._coupling_coefficient_synthetase = sympy.Mul(self._model.mu, sympy.Rational('1/3600'), self._synthetase_keff**-1, (1 + self._coupling_coefficient_trna_amount), evaluate = False)

		self.organelle = None

	@property
	def coupling_coefficient_trna_amount(self):
		return self._coupling_coefficient_trna_amount

	@coupling_coefficient_trna_amount.setter
	def coupling_coefficient_trna_amount(self, value):
		return NotImplemented # user should modify model.default_parameters

	# Backward compatibility
	@property
	def synthetase_keff(self):
		"""
		returns the synthetase keff value, not the coupling coefficient, in per second
		"""
		value = self._synthetase_keff
		try:
			return float(value)
		except:
			return float(value.xreplace(self._model.global_info['default_parameters']))

	# Backward compatibility
	@synthetase_keff.setter
	def synthetase_keff(self, value):
		"""
		value is the synthetase keff in per second, not the coupling coefficient in per hour
		this returns the coupling coefficient as growth rate divided by the keff times (1 + coupling_coefficient_trna_amount), in per hour
		"""
		self._coupling_coefficient_synthetase = sympy.Mul(self._model.mu, sympy.Rational('1/3600'), value**-1, (1 + self._coupling_coefficient_trna_amount), evaluate = False)

	@property
	def coupling_coefficient_synthetase(self):
		"""
		returns the coupling coefficient, not the synthetase keff value, in per hour
		"""
		return self._coupling_coefficient_synthetase

	@coupling_coefficient_synthetase.setter
	def coupling_coefficient_synthetase(self, value):
		"""
		value is the synthetase keff in per second, not the coupling coefficient in per hour
		this returns the coupling coefficient as growth rate divided by the keff times (1 + coupling_coefficient_trna_amount), in per hour
		"""
		self._coupling_coefficient_synthetase = sympy.Mul(self._model.mu, sympy.Rational('1/3600'), value**-1, (1 + self._coupling_coefficient_trna_amount), evaluate = False)

class TranslocationData(ProcessData):
	"""
	Class for storing information about a protein translocation pathway

	Parameters
	----------
	id : str
		Identifier for translocation pathway.

	model : :class:`coralme.core.model.MEModel`
		ME-model that the TranslocationData is associated with

	Attributes
	----------
	keff : float
		Effective turnover rate of the enzymes in the translocation pathway

	enzyme_dict : dict
		Dictionary containing enzyme specific information about the way it is
		coupled to protein translocation

		{enzyme_id: {length_dependent: <True or False>,
		 fixed_keff: <True or False>}}

	length_dependent_energy : bool
		True if the ATP cost of translocation is dependent on the length of
		the protein

	stoichiometry : dict
		Stoichiometry of translocation pathway, typically ATP/GTP hydrolysis

	"""

	def __init__(self, id, model):
		ProcessData.__init__(self, id, model)
		self.keff = 65.
		self.enzyme_dict = {}
		self.length_dependent_energy = False
		self.stoichiometry = {}

class PostTranslationData(ProcessData):
	"""
	Parameters
	----------
	id : str
		Identifier for post translation process.

	model : :class:`coralme.core.model.MEModel`
		ME-model that the PostTranslationData is associated with

	processed_protein : str
		ID of protein following post translational process

	preprocessed_protein : str
		ID of protein before post translational process

	Attributes
	----------
	translocation : set
		Translocation pathways involved in post translation reaction.

		Set of {:attr:`coralme.core.processdata.TranslocationData.id`}

	translocation_multipliers : dict
		Some proteins require different coupling of translocation enzymes.

		Dictionary of
		{:attr:`coralme.core.processdata.TranslocationData.id`: float}

	surface_area : dict
		If protein is translated into the inner or outer membrane, the surface
		area the protein occupies can be accounted for as well.

		Dictionary of {SA_+<inner_membrane or outer_membrane>: float}

	subreactions : :class:`collections.DefaultDict(int)`
		If a protein is modified following translation, this is accounted for
		here

		Dictionary of {subreaction_id: float}

	biomass_type : str
		If the subreactions add biomass to the translated gene, the
		biomass type (:attr:`coralme.core.compontent.Constraint.id`) of the
		modification must be defined.

	folding_mechanism : str
		ID of folding mechanism for post translation reaction

	aggregation_propensity : float
		Aggregation propensity for the protein

	keq_folding : dict
		Temperature dependant keq for folding protein

		Dictionary of {str(temperature): value}

	k_folding : dict
		Temperature dependant rate constant (k) for folding protein

		Dictionary of {str(temperature): value}

	propensity_scaling : float
		Some small peptides are more likely to be folded by certain
		chaperones. This is accounted for using propensity_scaling.

	"""

	def __init__(self, id, model, processed_protein, preprocessed_protein):
		ProcessData.__init__(self, id, model)
		self.processed_protein_id = processed_protein
		self.unprocessed_protein_id = preprocessed_protein

		# For translocation post translation reactions
		self.translocation = set()
		self.translocation_multipliers = {}
		self.surface_area = {}

		# For post translation modifications
		self.subreactions = collections.defaultdict(int)
		self.biomass_type = ''

		# For protein folding reactions (FoldME)
		self.folding_mechanism = ''
		self.aggregation_propensity = 0.
		self.keq_folding = {}
		self.k_folding = {}
		self.propensity_scaling = 1.
