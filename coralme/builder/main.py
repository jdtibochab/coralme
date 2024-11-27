#!/usr/bin/python3
# python imports
import io
import os
import re
import sys
import pickle
import shutil
import pathlib
import warnings
import collections
import subprocess

# third party imports
import tqdm
import sympy
import numpy
import pandas
import anyconfig

import cobra
import coralme

# configuration
log_format = '%(asctime)s %(message)s' #%(clientip)-15s %(user)-8s
bar_format = '{desc:<75}: {percentage:.1f}%|{bar:10}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
try:
	warnings.simplefilter(action = 'ignore', category = pandas.errors.SettingWithCopyWarning)
except:
	warnings.warn("This pandas version does not allow for correct warning handling. Pandas >=1.5.1 is suggested.")

import logging
#https://stackoverflow.com/questions/36408496/python-logging-handler-to-append-to-list
#Here is a naive, non thread-safe implementation:

class ListHandler(logging.Handler): # Inherit from logging.Handler
	"""
	ListHandler class to handle prints and logs.

	"""
	def __init__(self, log_list):
		# run the regular Handler __init__
		logging.Handler.__init__(self)
		# Our custom argument
		self.level = logging.WARNING
		self.formatter = log_format
		self.log_list = log_list
	def emit(self, record):
		# record.message is the log message
		try:
			self.log_list.append((record.asctime, record.message))
		except:
			pass
	def print_and_log(msg):
		print(msg)
		logging.warning(msg)

	def log(msg):
		logging.warning(msg)

class MEBuilder(object):
	"""
	MEBuilder class to coordinate the reconstruction of ME-models.

	Parameters
	----------
	*args:
		Positional arguments are passed as paths to JSON files that
		update the configuration of the parent class.
	**kwargs:
		Further keyword arguments are passed on as dictionaries
		to update the configuration of the parent class.

	"""
	def __init__(self, *args, **kwargs):
		config = {}
		for input_file in args:
			with open(input_file, 'r') as infile:
				config.update(anyconfig.load(infile))

		if kwargs:
			config.update(kwargs)
		self.configuration = config
		self.me_model = coralme.core.model.MEModel(config.get('ME-Model-ID', 'coralME'), config.get('growth_key', 'mu'))
		self.curation_notes = coralme.builder.notes.load_curation_notes(
			self.configuration['out_directory'] + '/curation_notes.json'
		)
		self.logger = {
			'MEBuilder' : coralme.builder.main.ListHandler([]),
			'MEReconstruction-step1' : coralme.builder.main.ListHandler([]),
			'MEReconstruction-step2' : coralme.builder.main.ListHandler([]),
			'METroubleshooter' : coralme.builder.main.ListHandler([])
			}

		data = \
			'code,interpretation,gram\n' \
			'CCI-CW-BAC-POS-GP,Cell_Wall,pos\n' \
			'CCI-OUTER-MEM-GN,Outer_Membrane,neg\n' \
			'CCI-PM-BAC-NEG-GN,Inner_Membrane,neg\n' \
			'CCI-PM-BAC-POS-GP,Plasma_Membrane,pos\n' \
			'CCO-MEMBRANE,Membrane,'

		self.location_interpreter = pandas.read_csv(io.StringIO(data), index_col = 0)

		# check user options
		exists = []
		for filename in [
			"m-model-path",
			"genbank-path",

			"df_TranscriptionalUnits",
			"df_matrix_stoichiometry",
			"df_matrix_subrxn_stoich",
			"df_metadata_orphan_rxns",
			"df_metadata_metabolites",
			"df_reaction_keff_consts",

			"biocyc.genes",
			"biocyc.prots",
			"biocyc.TUs",
			"biocyc.RNAs",
			"biocyc.seqs",
			]:

			if config.get(filename, None) is None:
				pass
			elif config.get(filename, '') == '':
				pass
			else:
				exists.append([config[filename], pathlib.Path(config[filename]).exists()])

		if not all([ y for x,y in exists ]):
			raise FileNotFoundError('Check the path to the {:s} file(s).'.format(', '.join([ x for x,y in exists if y == False ])))

		for option in [ "out_directory", "log_directory" ]:
			if config.get(option, '.') == '':
				self.configuration[option] = config.get('ME-Model-ID', 'coralME')

		# debug data
		self.debug_data = {}

	def generate_files(self, overwrite = True):
		"""Performs the Synchronize and Complement steps of the reconstruction.

		This function will read the Organism and the Reference. It will
		synchronize the input files, complement them, and finally build
		the OSM for the Organism.

		Parameters
		----------
		overwrite : bool
			If True, overwrite the OSM using the defined path in the configuration.
		"""
		config = self.configuration
		model = config.get('ME-Model-ID', 'coralME')
		directory = config.get('log_directory', '.')
		#if overwrite and os.path.exists(directory):
			#shutil.rmtree(directory + '/blast_files_and_results')
			#shutil.rmtree(directory + '/building_data')

		if not os.path.exists(directory):
			os.mkdir(directory)

		log = logging.getLogger() # root logger
		for hdlr in log.handlers[:]: # remove all old handlers
			log.removeHandler(hdlr)

		logging.basicConfig(filename = '{:s}/MEBuilder-{:s}.log'.format(directory, model), filemode = 'w', level = logging.WARNING, format = log_format)
		log.addHandler(self.logger['MEBuilder'])
		#log.addHandler(logging.StreamHandler(sys.stdout))
		logging.captureWarnings(True)

		sep = ''
		ListHandler.print_and_log("{}Initiating file processing...".format(sep))

		# Read organism
		self.org = coralme.builder.organism.Organism(config, is_reference = False)
		self.org.get_organism()
		self.curation_notes = self.org.curation_notes
		# self.org.rpod = ''
		# self.org.get_rna_polymerase(force_RNAP_as='')

		logging.warning("Modifying and preparing M-model")
# 		self.curate()
		self.prepare_model()

		# ## Homology with reference
		# Curation note: Check which locus tag field in the genbank agrees with those in the biocyc files.
		# Make sure you have specified that field in the beginning of this notebook.
		# Reference
		if bool(config.get('dev_reference', False)) or bool(config.get('user_reference', False)):
			logging.warning("Reading reference")

			self.ref = coralme.builder.organism.Organism(config.copy(), is_reference = True)
			self.ref.get_organism()

			folder = self.org.blast_directory
			if bool(config.get('run_bbh_blast', True)):
				blast_threads = config.get('blast_threads', 1)
				ListHandler.print_and_log("~ Running BLAST with {} threads...".format(blast_threads))
				self.org.gb_to_faa('org', element_types = {'CDS'}, outdir = self.org.blast_directory)
				self.ref.gb_to_faa('ref', element_types = {'CDS'}, outdir = self.org.blast_directory)

				def execute(cmd):
					if os.name == 'nt':
						cmd = 'cmd /c {:s}'.format(cmd)
					cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
					out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

				# make blast databases
				execute('makeblastdb -in {:s}/org.faa -dbtype prot -out {:s}/org'.format(folder, folder))
				execute('makeblastdb -in {:s}/ref.faa -dbtype prot -out {:s}/ref'.format(folder, folder))

				# bidirectional blast
				execute('blastp -db {:s}/org -query {:s}/ref.faa -num_threads {} -out {:s}/org_as_db.txt -outfmt 6'.format(folder, folder, blast_threads, folder))
				execute('blastp -db {:s}/ref -query {:s}/org.faa -num_threads {} -out {:s}/ref_as_db.txt -outfmt 6'.format(folder, folder, blast_threads, folder))

				#os.system('{}/auto_blast.sh {}'.format(self.directory,self.org.directory))
				ListHandler.print_and_log('BLAST done.')

			# #### Reciprocal hits
			logging.warning("Getting homologs")

			self.get_homology(evalue = self.org.config.get("e_value_cutoff", 1e-10))
			self.homology.mutual_hits_df.to_csv('{:s}/mutual_hits.txt'.format(folder))
			#self.homology.mutual_hits_df.to_csv(self.org.directory+'{:s}/mutual_hits.txt'.format(folder))

			# #### Get enzyme homology
			self.homology.get_complex_homology()

			# #### Update model info with homology
			logging.warning("Updating from homology")
			self.update_from_homology()

		filename = self.org.config.get('df_TranscriptionalUnits', self.org.directory + "TUs_from_biocyc.txt")
		filename = self.org.directory + "TUs_from_biocyc.txt" if filename == '' else filename

		df = self.org.TU_df
		df = df.sort_index(inplace = False)

		if overwrite:
			with open(filename, 'w') as outfile:
				self.org.TU_df.to_csv(outfile, sep = '\t')
				logging.warning('The BioCyc transcriptional data file was processed and overwritten into the {:s} file.'.format(filename))
		else:
			if pathlib.Path(filename).exists():
				logging.warning('Set \'overwrite = True\' to overwrite the {:s} file.'.format(filename))
			else:
				with open(filename, 'w') as outfile:
					self.org.TU_df.to_csv(outfile, sep = '\t')
					logging.warning('The BioCyc transcriptional data file was saved to the ./{:s} file.'.format(filename))
		self.configuration['df_TranscriptionalUnits'] = filename

		filename = self.org.directory + "subreaction_matrix.txt"
		if overwrite:
			with open(filename, 'w') as outfile:
				self.org.subreaction_matrix.to_csv(outfile, sep = '\t')
				logging.warning('The subreaction data file was processed and overwritten into the {:s} file.'.format(filename))
		else:
			if pathlib.Path(filename).exists():
				logging.warning('Set \'overwrite = True\' to overwrite the {:s} file.'.format(filename))
			else:
				with open(filename, 'w') as outfile:
					self.org.subreaction_matrix.to_csv(outfile, sep = '\t')
					logging.warning('The subreaction data file was saved to the ./{:s} file.'.format(filename))

		filename = self.org.directory + "me_metabolites.txt"
		if overwrite:
			with open(filename, 'w') as outfile:
				self.org.me_mets.to_csv(outfile, sep = '\t')
				logging.warning('The M to ME metabolite mapping file was processed and overwritten into the {:s} file.'.format(filename))
		else:
			if pathlib.Path(filename).exists():
				logging.warning('Set \'overwrite = True\' to overwrite the {:s} file.'.format(filename))
			else:
				with open(filename, 'w') as outfile:
					self.org.me_mets.to_csv(outfile, sep = '\t')
					logging.warning('The M to ME metabolite mapping file was saved to the ./{:s} file.'.format(filename))

		logging.warning("Updating enzyme-reaction association")
		self.update_enzyme_reaction_association()
		logging.warning("Processing RNA modifications")
		self.org.process_rna_modifications()
		logging.warning("Getting tRNA to codon information")
		self.get_trna_to_codon()

		# #### Biomass constituents
		self.org.biomass_constituents = config.get('flux_of_biomass_constituents', {})

		# Fill builder with dummy
		logging.warning("Filling files with CPLX_dummy")
		self.fill()
		# Final checks of builder
		logging.warning("Performing final checks of files")
		self.check()

		# Update manual curation files for user reference
		logging.warning("Generating filled manual curation files")
		self.org.manual_curation.save()

		# Update notes
		logging.warning("Generating curation notes")
		coralme.builder.notes.save_curation_notes(
				self.curation_notes,
				self.configuration['out_directory'] + '/curation_notes.json'
			)
		coralme.builder.notes.publish_curation_notes(
				self.curation_notes,
				self.configuration['out_directory']+ '/curation_notes.txt'
			)

		logging.warning("Saving modified M-model")
		filename = '{:s}/building_data/m_model_modified.json'.format(config.get('out_directory', './'))
		cobra.io.save_json_model(self.org.m_model, filename)
		config['m-model-path'] = filename

		logging.warning("Generating new configuration file")
		self.input_data(self.org.m_model, overwrite)
		ListHandler.print_and_log("{}File processing done.".format(sep))

		logging.shutdown()

		# We will remove duplicates entries in the log output
		with open('{:s}/MEBuilder-{:s}.log'.format(config.get('log_directory', '.'), config.get('ME-Model-ID', 'coralME')), 'w') as outfile:
			logger = self.logger['MEBuilder'].log_list

			tmp = pandas.DataFrame(logger)
			for idx, data in tmp.drop_duplicates(subset = 1).iterrows():
				outfile.write('{:s} {:s}\n'.format(data[0], data[1]))

	def prepare_model(self):
		"""Performs initial preparation of the M-model.

		This function will fix some known issues that M-models can

		Parameters
		----------
		overwrite : bool
			If True, overwrite the OSM using the defined path in the configuration.
		"""
		m_model = self.org.m_model
		target_compartments = {"c": "Cytosol", "e": "Extra-organism", "p": "Periplasm"}
		new_dict = {}
		warn_compartments = []
		for k,v in tqdm.tqdm(m_model.compartments.items(),
						'Gathering M-model compartments...',
						bar_format = bar_format,
						total=len(m_model.compartments)):
			if k in target_compartments:
				new_dict[k] = target_compartments[k]
			else:
				warn_compartments.append(k)
				new_dict[k] = k
		m_model.compartments = new_dict
		gene_list = []
		for m in tqdm.tqdm(m_model.metabolites,
						'Fixing compartments in M-model metabolites...',
						bar_format = bar_format):
			if not m.compartment:
				logging.warning("Fixing compartment for metabolite {}".format(m.id))
				m.compartment = m.id[-1]
		for r in tqdm.tqdm(m_model.reactions,
						'Fixing missing names in M-model reactions...',
						bar_format = bar_format):
			if not r.name or isinstance(r.name, float):
				logging.warning("Fixing name for reaction {}".format(r.id))
				r.name = r.id

		# Solve m_model
		solution = m_model.optimize()
		if not solution or solution.objective_value < 1e-6:
			self.org.curation_notes['prepare_model'].append({
				'msg':'M-model cannot grow',
				'importance':'critical',
				'to_do':'Check that the model you provided works'})

		try:
			self.org.biomass = str(m_model.objective.expression.as_two_terms()[0]).split('*')[1]
			biomass_rxn = m_model.reactions.get_by_id(self.org.biomass)
			logging.warning('{} was identified as the biomass reaction'.format(biomass_rxn.id))
		except:
			self.org.biomass = None
			biomass_rxn = None
			logging.warning('Could not identify biomass reaction')

		self.org.GAM = self.configuration.get('gam',None)
		# Get GAM
		if self.org.GAM is None:
			adp = m_model.metabolites.adp_c
			if biomass_rxn is not None and adp in biomass_rxn.metabolites:
				self.org.GAM = biomass_rxn.metabolites[adp]
				logging.warning('GAM identified with value {}'.format(self.org.GAM))
			else:
				self.org.GAM = self.configuration.get('gam', self.me_model._gam)
				self.org.curation_notes['prepare_model'].append({
					'msg':'GAM could not be identified from biomass reaction, setting a standard value of {:g}. adp_c is not present as a product.'.format(self.me_model._gam),
					'importance':'high',
					'to_do':'Check whether the biomass reaction was read or defined correctly. You can define GAM with builder.configuration[\'gam\'] = GAM_value'})
		self.org.NGAM = self.configuration.get('ngam',None)
		if self.org.NGAM is None:
			# Get NGAM
			NGAMs = ['NGAM','ATPM']
			for r in NGAMs:
				if r in m_model.reactions:
					rxn = m_model.reactions.get_by_id(r)
					if rxn.lower_bound <= 0:
						continue
					self.org.NGAM = rxn.lower_bound
					logging.warning('{} was identified as NGAM with value {}'.format(r,self.org.NGAM))
					break
		if self.org.NGAM is None:
			self.org.NGAM = self.configuration.get('ngam', self.me_model._ngam)
			self.org.curation_notes['prepare_model'].append({
				'msg':'NGAM could not be identified in M-model, setting a standard value of {:g}.'.format(self.me_model._ngam),
				'importance':'high',
				'to_do':'Manually define NGAM with builder.configuration[\'ngam\'] = NGAM_value. Check if a reaction with identifier NGAM or ATPM has a zero or negative lower bound.'})
		elif self.org.NGAM == 0:
			self.org.NGAM = self.configuration.get('ngam', self.me_model._ngam)
			self.org.curation_notes['prepare_model'].append({
				'msg':'NGAM was identified from reaction {}, but its lower bound is 0. NGAM set to a standard value of {:g}.'.format(rxn.id, self.me_model._ngam),
				'importance':'high',
				'to_do':'Manually define NGAM with builder.configuration[\'ngam\'] = NGAM_value'})
		# Warnings
		if warn_compartments:
			self.org.curation_notes['prepare_model'].append({
				'msg':'Some compartments in m_model are unknown',
				'triggered_by':warn_compartments,
				'importance':'medium',
				'to_do':'Check whether the compartment is correct. If not, change it in the reaction ID in the m_model.'})

	def get_homology(self, evalue=1e-10):
		"""Calculates homology between Organism and Reference.

		Parameters
		----------
		evalue : float, default 1e-10
			Sets the E-value cutoff for calling protein
			homologs using BLAST.
		"""
		self.homology = coralme.builder.homology.Homology(self.org, self.ref, evalue = evalue)

	def get_trna_to_codon(self):
		"""Gets tRNA to codon association from the Genome.
		"""
		import Bio
		from Bio import SeqIO, Seq, SeqFeature, SeqUtils

		me_model = self.me_model
		contigs = self.org.contigs

		# Dictionary of tRNA locus ID to the model.metabolite object. It accounts for misacylation
		trna_to_aa = {}

		# Dictionary of tRNA locus ID to the model.metabolite object. It accounts for misacylation
		trna_misacylation = self.configuration.get('trna_misacylation',dict())

		# Dictionary of tRNA locus ID to amino acid, one dict of tRNAs per organelle type
		# aa2trna does not account for misacylation
		aa2trna = { 'c' : {} } # prokaryotes and eukaryotes
		if self.configuration.get('domain','prokaryote').lower() not in ['prokaryote', 'bacteria']:
			aa2trna.update({'m' : {}, 'h' : {}}) # mitochondria and plastids

		# Translation tables, one table per organelle type
		transl_tables = { 'c' : set() } # prokaryotes and eukaryotes
		if self.configuration.get('domain','prokaryote') not in ['prokaryote', 'bacteria']:
			transl_tables.update({'m' : set(), 'h' : set()}) # mitochondria and plastids

		# Messages
		msg1 = 'From the tRNA misacylation dictionary, the {:s} gene [tRNA({:s})] is loaded and converted into {:s}-tRNA({:s}). Make sure a MetabolicReaction to convert a {:s}-tRNA({:s}) into a {:s}-tRNA({:s}) is present in the ME-model.'
		msg2 = 'From the tRNA misacylation dictionary, the {:s} gene [tRNA({:s})] is loaded and converted into {:s}-tRNA({:s}). No further modification needs to take place.'

		canonical_aas = [
			'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
			'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
			]

		for contig in contigs:
			iterator = tqdm.tqdm(contig.features, 'Getting tRNA to codon dictionary from {}'.format(contig.id), bar_format = bar_format) if len(contigs) < 10 else contig.features
			for feature in iterator:
				# Find organelle in source
				if feature.type in ['source', 'region']: # genbanks from gff+fna do not contain a 'source' feature
					organelle = feature.qualifiers.get('organelle', [None])[0]
					continue
				seq = feature.extract(contig).seq.replace('-', '')
				if feature.type == 'CDS':
					# Add the translation table
					prot = feature.qualifiers.get('translation', [''])[0]
					transl_table = feature.qualifiers.get('transl_table', ['1'])[0]

					# Add the translation table per organelle
					if organelle is None:
						transl_tables['c'].add(int(transl_table))
						if me_model.global_info['domain'].lower() not in ['prokaryote', 'bacteria']:
							logging.warning('Contig \'{:s}\' does not report an organelle type.'.format(contig.id))
					elif organelle.lower() in ['mitochondria', 'mitochondrion']:
						transl_tables['m'].add(int(transl_table))
					elif organelle.lower() in ['chloroplast', 'plastid']:
						transl_tables['h'].add(int(transl_table))
					else:
						continue

				if not feature.type == 'tRNA':
					continue

				bnum = self.org._get_feature_locus_tag(feature)
				if bnum is None:
					continue
				aa = feature.qualifiers.get('product', ['tRNA-None'])[0].split('-')[1]
				if aa in canonical_aas + ['Asx', 'Glx', 'fMet', 'Sec']:
					pass
				else:
					logging.warning('The tRNA \'{:s}\' is not associated to a valid product name (tRNA-Amino acid 3 letters code)'.format(bnum))
					continue
				msg = 'The tRNA \'{:s}\' is associated to two amino acids. The \'trna_misacylation\' dictionary was modified to attempt load the correct amino acid.'
				# Special tRNA(Asx) that can be loaded with Asn (EC 6.1.1.22) or Asp (EC 6.1.1.12)
				if aa == 'Asx':
					trna_misacylation['Asx'] = 'Asp'
					logging.warning(msg.format(bnum))
				# Special tRNA(Glx) that can be loaded with Gln (EC 6.1.1.18) or Glu (EC 6.1.1.17)
				if aa == 'Glx':
					trna_misacylation['Glx'] = 'Glu'
					logging.warning(msg.format(bnum))

				if aa in trna_misacylation.keys():
					# misacylation only in mitochondria and chloroplasts
					filter1a = me_model.global_info['domain'].lower() in ['eukarya', 'eukaryote']
					filter1b = str(organelle).lower() in ['mitochondria', 'mitochondrion', 'chloroplast', 'plastid']
					# misacylation in the cytoplasm of Gram-positive eubacteria (and other bacteria such as cyanobacteria)
					filter2 = me_model.global_info['domain'].lower() in ['bacteria', 'prokaryote']

					if filter1a and filter1b or filter2:
						trna_to_aa[bnum] = trna_misacylation[aa]
						if aa.endswith('x'):
							logging.warning(msg2.format(bnum, aa, trna_misacylation[aa], aa))
						else:
							logging.warning(msg1.format(bnum, aa, trna_misacylation[aa], aa, trna_misacylation[aa], aa, aa, aa))
					else:
						# misacylation is not valid in the compartment and domain
						trna_to_aa[bnum] = aa
				else:
					trna_to_aa[bnum] = aa

				if organelle is None:
					aa2trna['c'][bnum] = aa
				elif organelle.lower() in ['mitochondria', 'mitochondrion']:
					aa2trna['m'][bnum] = aa
				elif organelle.lower() in ['chloroplast', 'plastid']:
					aa2trna['h'][bnum] = aa

		# trna_to_codon does not account for misacylation: { 'tRNA ID' : 'Amino acid to load into the tRNA' }
		trna_to_aa = { k:v.replace('fMet', 'Met') for k,v in trna_to_aa.items() }

		for organelle, aa2trna_dct in aa2trna.items():
			aa2trna_dct = { k:v.capitalize().split('_')[0] if 'fMet' not in v else 'fMet' for k,v in aa2trna_dct.items() }
			aa2trna_df = pandas.DataFrame(data = [aa2trna_dct.values(), aa2trna_dct.keys()]).T
			aa2trna_df = aa2trna_df.groupby(0).agg({1: lambda x: x.tolist()})
			if 'fMet' in aa2trna_df.index:
				aa2trna_df.loc['Met'] = aa2trna_df.loc['fMet'] + aa2trna_df.loc['Met']

			aa2trna[organelle] = aa2trna_df

			if not aa2trna_df.empty:
				# assign START tRNAs to every fMet-tRNA (Met-tRNA if not) and check if at least one tRNA was identified
				if 'fMet' in aa2trna_df.index:
					if len(me_model.global_info['START_tRNA']) == 0:
						me_model.global_info['START_tRNA'] = list(aa2trna_df.loc['fMet'])[0]

				if 'Met' in aa2trna_df.index:
					if len(me_model.global_info['START_tRNA']) == 0:
						me_model.global_info['START_tRNA'] = list(aa2trna_df.loc['Met'])[0]

				# final check
				if len(me_model.global_info['START_tRNA']) == 0:
					logging.warning('Unable to identify at least one \'tRNA-Met\' or \'tRNA-fMet\' annotation from the \'Definition\' column in the organism-specific matrix.')
			else:
				logging.warning('No tRNA genes were identified from their locus tags.')

		# DataFrame mapping tRNAs (list) and the encoded amino acid (index), per organelle
		# aa2trna derives from trna_to_aa, so it also accounts for misacylation: { 'organelle ID' : 'DataFrame of amino acid to load into the tRNA' }
		aa2codon = {}
		trna_to_codon = {}
		for organelle, transl_table in transl_tables.items():
			if len(transl_table) == 0:
				continue

			codon_table = Bio.Data.CodonTable.generic_by_id[list(transl_table)[0]]

			dct = { k.replace('T', 'U'):SeqUtils.seq3(v) for k,v in codon_table.forward_table.items() if 'U' not in k }
			aa2codons = pandas.DataFrame(data = [dct.keys(), dct.values()]).T.groupby(1).agg({0: lambda x: x.tolist()})
			# aa2codons derives from the translation table and maps amino acids to the codon
			aa2codon[organelle] = aa2codons
			if list(transl_table)[0] == 11:
				aa2codons.loc['Sec'] = [['UGA']] # an internal UGA encodes selenocysteine

			df = pandas.concat([aa2codons, aa2trna[organelle]], axis = 1).dropna(how = 'any').explode(1)

			# Check amino acids
			warn_trnas = []
			for aa in canonical_aas:
				if aa in aa2trna[organelle].index:
					pass
				else:
					warn_trnas.append(aa)
			if warn_trnas:
				self.org.curation_notes['get_trna_to_codon'].append({
					'msg':'Some tRNAs could be missing in the files.',
					'triggered_by':warn_trnas,
					'importance':'critical',
					'to_do':'Identify the respective missing tRNAs and add them to a file (RNAs.txt or genome.gb)'})

			trna_to_codon_organelle = { k:v + ['START'] if k in me_model.global_info['START_tRNA'] else v for k,v in zip(df[1].values, df[0].values) }
			trna_to_codon[organelle] = trna_to_codon_organelle

		me_model.global_info['trna_to_aa'] = trna_to_aa
		me_model.global_info['trna_to_codon'] = trna_to_codon

		me_model.global_info['trna_misacylation'] = trna_misacylation

		return None

	def update_enzyme_stoichiometry(self):
		complexes_df = self.org.complexes_df
		org_cplx_homolog = self.homology.org_cplx_homolog
		ref_complexes_df = self.ref.complexes_df
		mutual_hits = self.homology.mutual_hits

		cplx_stoich = {}
		for c, rc in org_cplx_homolog.items():
			if "generic" in c:
				continue
			if rc not in ref_complexes_df.index:
				continue
			for rg in ref_complexes_df["genes"][rc].split(" AND "):
				rg_id = re.findall(r'.*(?=\(\d*\))', rg)[0]
				coeff = re.findall(r"(?<=\()[0-9]{1,}", rg)
				if coeff:
					coeff = coeff[0]
				else:
					coeff = ""
				if c not in cplx_stoich:
					cplx_stoich[c] = {}
				cplx_stoich[c][mutual_hits[rg_id]] = coeff

		for c, stoich in cplx_stoich.items():
			complexes_df["genes"][c] = " AND ".join(
				[g + "({})".format(coeff) for g, coeff in stoich.items()]
			)
		self.org.complexes_df = complexes_df

	def update_protein_modification(self):
		cplx_homolog = self.homology.org_cplx_homolog
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		complexes_df = self.org.complexes_df
		ref_mod_df = self.ref.protein_mod.reset_index().set_index("Core_enzyme")
		cplx_cofactor_dict = {}
		protein_mod_dict = {}
		for c, row in complexes_df.iterrows():
			if c in cplx_homolog:
				rc = cplx_homolog[c]
				if rc in ref_mod_df.index:
					if c not in cplx_cofactor_dict:
						cplx_cofactor_dict[c] = {}
					ref_mods = ref_mod_df.loc[[rc]]
					for _,row in ref_mods.iterrows():
						mods = row["Modifications"].split(" AND ")
						cplx = c
						cofs = []
						coeffs = []
						for mod in mods:
							cof = re.findall(r".*(?=\()", mod)[0]
							coeff = re.findall(r"(?<=\()[0-9]{1,}", mod)
							if coeff:
								coeff = coeff[0]
								cplx += "_mod_{}({})".format(cof, coeff)
							else:
								cplx += "_mod_{}".format(cof)
								coeff = ""
							cofs.append(cof)
							coeffs.append(coeff)
						ref_cplx_homolog[row["Modified_enzyme"]] = cplx
						cplx_homolog[cplx] = row["Modified_enzyme"]
						if cplx in self.org.protein_mod.index:
							continue
						protein_mod_dict[cplx] = {}
						protein_mod_dict[cplx]["Core_enzyme"] = c
						protein_mod_dict[cplx]["Modifications"] = " AND ".join(
							"{}({})".format(cof, coeff)
							for cof, coeff in zip(cofs, coeffs)
						)
						protein_mod_dict[cplx]["Source"] = "Homology"

		protein_mod = pandas.DataFrame.from_dict(protein_mod_dict).T
		protein_mod.index.name = "Modified_enzyme"
		self.org.protein_mod = pandas.concat([self.org.protein_mod,protein_mod])

	def update_enzyme_reaction_association(self):
		enz_rxn_assoc_df = self.org.enz_rxn_assoc_df
		protein_mod = self.org.protein_mod
		for r,row in tqdm.tqdm(enz_rxn_assoc_df.iterrows(),
					'Updating enzyme reaction association...',
					bar_format = bar_format,
					total=enz_rxn_assoc_df.shape[0]):
			if r in self.org.manual_curation.enz_rxn_assoc_df.data.index:
				# Only update those not in manual curation
				continue
			cplxs = row['Complexes']
			reaction_cplx_list = []
			for cplx_id in cplxs.split(" OR "):
				cplx_id = cplx_id.split("_mod_")[0]
				if cplx_id in protein_mod["Core_enzyme"].values:
					# Use modifications
					cplx_mods = protein_mod[
						protein_mod["Core_enzyme"].eq(cplx_id)
					].index
					for cplx_id in cplx_mods:
						if "Oxidized" in cplx_id:
							cplx_id = cplx_id.split("_mod_Oxidized")[0]
						reaction_cplx_list.append(cplx_id)
				else:
					reaction_cplx_list.append(cplx_id)
			enz_rxn_assoc_df.at[r,'Complexes'] = " OR ".join(reaction_cplx_list)

	def update_TU_df(self):
		return NotImplemented

	def protein_location_from_homology(self):
		protein_location = self.org.protein_location
		complexes_df = self.org.complexes_df
		org_proteins_df = self.org.proteins_df
		org_cplx_homolog = self.homology.org_cplx_homolog
		mutual_hits = self.homology.mutual_hits
		ref_protein_location = self.ref.protein_location
		if not isinstance(ref_protein_location, pandas.DataFrame):
			return
		for rc,row in tqdm.tqdm(ref_protein_location.iterrows(),
					'Updating protein location from homology...',
					bar_format = bar_format,
					total=ref_protein_location.shape[0]):
			ref_gene = re.findall(r'.*(?=\(.*\))', row['Protein'])[0]
			if ref_gene not in mutual_hits:
								continue
			org_gene = mutual_hits[ref_gene]
			ref_info = ref_protein_location[ref_protein_location['Protein'].str.contains(ref_gene)]
			gene_string = r'{}\('.format(org_gene)
			org_cplxs = complexes_df[complexes_df['genes'].str.contains(gene_string)].index
			for org_cplx in org_cplxs:
				if protein_location.any().any() and \
					org_cplx in protein_location.index and \
					protein_location.loc[[org_cplx]]['Protein'].str.contains(gene_string).any().any():
							# Check if already in protein location, if not add.
					continue
				tmp = pandas.DataFrame.from_dict(
					{ org_cplx: {
						"Complex_compartment": ref_info["Complex_compartment"].values[0],
						"Protein": '{}()'.format(org_gene),
						"Protein_compartment": ref_info["Protein_compartment"].values[0],
						"translocase_pathway": ref_info["translocase_pathway"].values[0],
									}
								}).T
				protein_location = pandas.concat([protein_location, tmp], axis = 0, join = 'outer')
		protein_location.index.name = 'Complex'
		self.org.protein_location = protein_location

	def update_translocation_multipliers(self):
		ref_multipliers = self.ref.translocation_multipliers
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		mutual_hits = self.homology.mutual_hits
		multipliers = {}
		for ref_cplx, genes in tqdm.tqdm(ref_multipliers.items(),
					'Updating translocation multipliers from homology...',
					bar_format = bar_format,
					total=len(ref_multipliers)):
			if ref_cplx not in ref_cplx_homolog:
				continue
			org_cplx = ref_cplx_homolog[ref_cplx]
			multipliers[org_cplx] = {}
			for g, value in genes.items():
				if not value:
					continue
				if g not in mutual_hits:
					continue
				org_g = mutual_hits[g]
				multipliers[org_cplx][org_g] = value
		self.org.translocation_multipliers = multipliers

	def update_lipoprotein_precursors(self):
		ref_lipoprotein_precursors = self.ref.lipoprotein_precursors
		mutual_hits = self.homology.mutual_hits
		lipoprotein_precursors = {}
		for c, g in tqdm.tqdm(ref_lipoprotein_precursors.items(),
					'Updating lipoprotein precursors from homology...',
					bar_format = bar_format,
					total=len(ref_lipoprotein_precursors)):
			if g in mutual_hits:
				og = mutual_hits[g]
				lipoprotein_precursors[c] = og
		self.org.lipoprotein_precursors = lipoprotein_precursors

	def update_cleaved_methionine(self):
		cleaved_methionine = self.org.cleaved_methionine
		ref_cleaved_methionine = self.ref.cleaved_methionine
		mutual_hits = self.homology.mutual_hits
		for g in tqdm.tqdm(ref_cleaved_methionine,
					'Updating cleaved-methionine proteins from homology...',
					bar_format = bar_format):
			if g not in mutual_hits:
				continue
			cleaved_methionine.append(mutual_hits[g])
		self.org.cleaved_methionine = cleaved_methionine

	def update_me_mets(self):
		ref_me_mets = self.ref.me_mets
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		protein_mod = self.org.protein_mod.reset_index().set_index("Core_enzyme")
		me_mets = self.org.me_mets
		m_model = self.org.m_model
		d = {}
		warn_skip = []
		warn_found = []
		for ref_m, row in tqdm.tqdm(ref_me_mets.iterrows(),
					'Mapping M-metabolites to E-metabolites...',
					bar_format = bar_format,
					total=ref_me_mets.shape[0]):
			if ref_m not in m_model.metabolites:
				if "__" in ref_m:
					ref_m = ref_m.replace("__", "_")
					if ref_m not in m_model.metabolites:
						warn_skip.append(ref_m)
						continue
					else:
						warn_found.append(ref_m)
			if ref_m in me_mets.index:
				if me_mets.loc[ref_m]['type'] != 'CURATE':
					continue
				me_mets.drop(ref_m,inplace=True)
			ref_me = row["me_id"]
			ref_changetype = row['type']
			d[ref_m] = {}
			me_id = ''
			if ref_me in ref_cplx_homolog:
				org_me = ref_cplx_homolog[ref_me]
				me_id = org_me
				changetype = "REPLACE"
			elif ref_changetype == "REMOVE":
				changetype = "REMOVE"
			else:
				changetype = "CURATE"
			d[ref_m]["me_id"] = me_id
			d[ref_m]["type"] = changetype
		if d:
			df = pandas.DataFrame.from_dict(d).T
			df.index.name = "id"
			me_mets = pandas.concat([me_mets, df], axis = 0, join = 'outer')
			self.org.me_mets = me_mets.fillna('')
		if warn_skip:
			self.org.curation_notes['update_me_mets'].append({
				'msg':'Some metabolites in me_metabolites.txt are not in m_model, so they were skipped.',
				'triggered_by':warn_skip,
				'importance':'medium',
				'to_do':'Confirm these metabolites are correctly defined in me_metabolites.txt'})
		if warn_found:
			self.org.curation_notes['update_me_mets'].append({
				'msg':'Some metabolites in me_metabolites.txt were found in reference m_model after replacing __ with _',
				'triggered_by':warn_found,
				'importance':'medium',
				'to_do':'Confirm these metabolites are correctly defined in me_metabolites.txt'})

	def update_generics_from_homology(self):
		generic_dict = self.org.generic_dict
		ref_generic_dict = self.ref.generic_dict
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_generic_dict.items(),
					'Updating generics from homology...',
					bar_format = bar_format,
					total=len(ref_generic_dict)):
			if k not in generic_dict:
				generic_dict[k] = {"enzymes":[]}
				#continue
			if generic_dict[k]['enzymes']:
				continue
			ref_cplxs = v['enzymes']
			for i in ref_cplxs:
				if i in ref_cplx_homolog:
					homolog = ref_cplx_homolog[i]
					if homolog not in generic_dict[k]['enzymes']:
						generic_dict[k]['enzymes'].append(homolog)

	def update_folding_dict_from_homology(self):
		org_folding_dict = self.org.folding_dict
		ref_folding_dict = self.ref.folding_dict
		mutual_hits = self.homology.mutual_hits
		for k, v in tqdm.tqdm(ref_folding_dict.items(),
					'Updating folding from homology...',
					bar_format = bar_format,
					total=len(ref_folding_dict)):
			ref_cplxs = v['enzymes']
			for i in ref_cplxs:
				if i in mutual_hits:
					homolog = mutual_hits[i]
					if homolog not in org_folding_dict[k]['enzymes']:
						org_folding_dict[k]['enzymes'].append(homolog)

	def update_ribosome_subreactions_from_homology(self):
		ref_ribosome_subreactions = self.ref.ribosome_subreactions
		org_ribosome_subreactions = self.org.ribosome_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_proteins = []
		for k, v in tqdm.tqdm(ref_ribosome_subreactions.items(),
					'Updating ribosome subreaction machinery from homology...',
					bar_format = bar_format,
					total=len(ref_ribosome_subreactions)):
			ref_cplx = v["enzyme"]
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[v["enzyme"]]
				if k not in org_ribosome_subreactions:
					org_ribosome_subreactions[k] = {"enzyme":""}
				defined_cplx = org_ribosome_subreactions[k]["enzyme"]
				if defined_cplx: continue
# 				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
				org_ribosome_subreactions[k]["enzyme"] = org_cplx
# 				else:
				warn_proteins.append({
					'subreaction':k,
					'defined_complex':defined_cplx,
					'inferred_complex':org_cplx
				})
		# Warnings
		if warn_proteins:
			self.org.curation_notes['update_ribosome_subreactions_from_homology'].append({
				'msg':'Some enzymes defined in me_builder.org.ribosome_subreactions are different from the ones inferred from homology',
				'triggered_by':warn_proteins,
				'importance':'medium',
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.ribosome_subreactions. Curate the inputs in ribosome_subreactions.txt accordingly.'})

# 	def update_rrna_modifications_from_homology(self):
# 		ref_rrna_modifications = self.ref.rrna_modifications
# 		org_rrna_modifications = self.org.rrna_modifications
# 		ref_cplx_homolog = self.homology.ref_cplx_homolog
# 		warn_proteins = []
# 		for k, v in tqdm.tqdm(ref_rrna_modifications.items(),
# 					'Updating rRNA modifications from homology...',
# 					bar_format = bar_format,
# 					total=len(ref_rrna_modifications)):
# 			ref_cplx = v["machine"]
# 			if ref_cplx in ref_cplx_homolog:
# 				org_cplx = ref_cplx_homolog[v["machine"]]
# 				defined_cplx = org_rrna_modifications[k]["machine"]
# 				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
# 					org_rrna_modifications[k]["machine"] = org_cplx
# 				else:
# 					warn_proteins.append({
# 						'subreaction':k,
# 						'defined_complex':defined_cplx,
# 						'inferred_complex':org_cplx
# 					})
# 		# Warnings
# 		if warn_proteins:
# 			self.org.curation_notes['update_rrna_modifications_from_homology'].append({
# 				'msg':'Some enzymes defined in me_builder.org.rrna_modifications are different from the ones inferred from homology',
# 				'triggered_by':warn_proteins,
# 				'importance':'medium',
# 				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.rrna_modifications. Curate the inputs in rrna_modifications.txt accordingly.'})

	def update_amino_acid_trna_synthetases_from_homology(self):
		ref_amino_acid_trna_synthetase = self.ref.amino_acid_trna_synthetase
		org_amino_acid_trna_synthetase = self.org.amino_acid_trna_synthetase
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_proteins = []
		for k, v in tqdm.tqdm(ref_amino_acid_trna_synthetase.items(),
					'Updating tRNA synthetases from homology...',
					bar_format = bar_format,
					total=len(ref_amino_acid_trna_synthetase)):
			ref_cplx = v
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[v]
				defined_cplx = org_amino_acid_trna_synthetase[k]
				if self.configuration.get('user_data', False) and defined_cplx: continue
# 				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
				org_amino_acid_trna_synthetase[k] = org_cplx
# 				else:
				warn_proteins.append({
					'amino_acid':k,
					'defined_ligase':defined_cplx,
					'inferred_ligase':org_cplx
				})
		# Warnings
		if warn_proteins:
			self.org.curation_notes['update_amino_acid_trna_synthetases_from_homology'].append({
				'msg':'Some enzymes defined in me_builder.org.amino_acid_trna_synthetase are different from the ones inferred from homology',
				'triggered_by':warn_proteins,
				'importance':'medium',
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.amino_acid_trna_synthetase. Curate the inputs in amino_acid_trna_synthetase.txt accordingly.'})

	def update_peptide_release_factors_from_homology(self):
		ref_peptide_release_factors = self.ref.peptide_release_factors
		org_peptide_release_factors = self.org.peptide_release_factors
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_proteins = []
		for k, v in tqdm.tqdm(ref_peptide_release_factors.items(),
					'Updating peptide release factors from homology...',
					bar_format = bar_format,
					total=len(ref_peptide_release_factors)):
			ref_cplx = v['enzyme']
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[ref_cplx]
				defined_cplx = org_peptide_release_factors[k]['enzyme']
				if defined_cplx: continue
# 				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
				org_peptide_release_factors[k]['enzyme'] = org_cplx
# 				else:
				warn_proteins.append({
					'subreaction':k,
					'defined_complex':defined_cplx,
					'inferred_complex':org_cplx
				})
		# Warnings
		if warn_proteins:
			self.org.curation_notes['update_peptide_release_factors_from_homology'].append({
				'msg':'Some enzymes defined in me_builder.org.peptide_release_factors are different from the ones inferred from homology',
				'triggered_by':warn_proteins,
				'importance':'medium',
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.peptide_release_factors. Curate the inputs in peptide_release_factors.txt accordingly.'})

	def update_initiation_subreactions_from_homology(self):
		ref_initiation_subreactions = self.ref.initiation_subreactions
		org_initiation_subreactions = self.org.initiation_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_initiation_subreactions.items(),
					'Updating translation initiation subreactions from homology...',
					bar_format = bar_format,
					total=len(ref_initiation_subreactions)):
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_initiation_subreactions[k]["enzymes"]
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i in defined_cplxs:
					continue
				defined_cplxs.append(i)

	def update_elongation_subreactions_from_homology(self):
		ref_elongation_subreactions = self.ref.elongation_subreactions
		org_elongation_subreactions = self.org.elongation_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_elongation_subreactions.items(),
					'Updating translation elongation subreactions from homology...',
					bar_format = bar_format,
					total=len(ref_elongation_subreactions)):
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_elongation_subreactions[k]["enzymes"]
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i in defined_cplxs:
					continue
				defined_cplxs.append(i)

	def update_termination_subreactions_from_homology(self):
		ref_termination_subreactions = self.ref.termination_subreactions
		org_termination_subreactions = self.org.termination_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_termination_subreactions.items(),
					'Updating translation termination subreactions from homology...',
					bar_format = bar_format,
					total=len(ref_termination_subreactions)):
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_termination_subreactions[k]["enzymes"]
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i in defined_cplxs:
					continue
				defined_cplxs.append(i)

	def update_special_trna_subreactions_from_homology(self):
		ref_special_trna_subreactions = self.ref.special_trna_subreactions
		org_special_trna_subreactions = self.org.special_trna_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_special_trna_subreactions.items(),
					'Updating special tRNA subreactions from homology...',
					bar_format = bar_format,
					total=len(ref_special_trna_subreactions)):
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_special_trna_subreactions[k]["enzymes"]
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_rna_degradosome_from_homology(self):
		ref_rna_degradosome = self.ref.rna_degradosome
		org_rna_degradosome = self.org.rna_degradosome
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_rna_degradosome.items(),
					'Updating RNA degradosome composition from homology...',
					bar_format = bar_format,
					total=len(ref_rna_degradosome)):
			ref_cplxs = v['enzymes']
			defined_cplxs = org_rna_degradosome[k]['enzymes']
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i in defined_cplxs:
					continue
				defined_cplxs.append(i)

	def update_excision_machinery_from_homology(self):
		ref_excision_machinery = self.ref.excision_machinery
		org_excision_machinery = self.org.excision_machinery
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_excision_machinery.items(),
					'Updating excision machinery from homology...',
					bar_format = bar_format,
					total=len(ref_excision_machinery)):
			ref_cplxs = v['enzymes']
			defined_cplxs = org_excision_machinery[k]['enzymes']
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if self._is_base_complex_in_list(i,defined_cplxs):
					continue
				defined_cplxs.append(i)

	def update_special_modifications_from_homology(self):
		ref_special_trna_subreactions = self.ref.special_modifications
		org_special_trna_subreactions = self.org.special_modifications
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_special_trna_subreactions.items(),
					'Updating tRNA subreactions from homology...',
					bar_format = bar_format,
					total=len(ref_special_trna_subreactions)):
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_special_trna_subreactions[k]["enzymes"]
			if defined_cplxs: continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if v["stoich"]:
					org_special_trna_subreactions[k]["stoich"] = v["stoich"]
				if self._is_base_complex_in_list(i,defined_cplxs):
					continue
				defined_cplxs.append(i)

	def _is_base_complex_in_list(self,cplx,lst):
		return cplx in set(i.split('_mod_')[0] for i in lst)

	def update_rna_modification_from_homology(self):
		ref_rna_modification = self.ref.rna_modification_df
		org_rna_modification = self.org.rna_modification_df
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for mod,row in ref_rna_modification.iterrows():
			enzymes = set(row['enzymes'].split('AND'))
			positions = row['positions'].split(',')
			if mod in org_rna_modification.index:
				df = org_rna_modification.loc[[mod]]
				mod_type = row['type']
				if df['type'].str.contains(mod_type).any():
					continue
			hits = enzymes.intersection(set(ref_cplx_homolog.keys()))
			if len(hits) == len(enzymes):
				row = row.copy()
				row['enzymes'] = ' AND '.join([ref_cplx_homolog[i] for i in enzymes])
				row['positions'] = ','.join(positions)
				row['source'] = 'Homology'
				df = pandas.DataFrame(row).T
				org_rna_modification = pandas.concat([org_rna_modification,df])
		org_rna_modification.index.name = 'modification'
		self.org.rna_modification_df = org_rna_modification

# 	def update_rna_modification_from_homology(self):
# 		ref_rna_modification = self.ref.rna_modification
# 		org_rna_modification = self.org.rna_modification
# 		ref_cplx_homolog = self.homology.ref_cplx_homolog
# 		defined_mods = set()
# 		for _,v in org_rna_modification.items():
# 			for i in v:
# 				defined_mods.add(i)
# 		for k, v in tqdm.tqdm(ref_rna_modification.items(),
# 					'Updating RNA modification machinery from homology...',
# 					bar_format = bar_format,
# 					total=len(ref_rna_modification)):
# 			if k not in ref_cplx_homolog: continue
# 			org_cplx = ref_cplx_homolog[k]
# 			if self._is_base_complex_in_list(org_cplx,list(org_rna_modification.keys())):
# 				continue
# 			if org_cplx not in org_rna_modification:
# 				org_rna_modification[org_cplx] = []
# 			org_rna_modification[org_cplx] += [i for i in v.copy() if i not in defined_mods]
# 			org_rna_modification[org_cplx] = list(set(org_rna_modification[org_cplx]))
# 		org_rna_modification = {k:v for k,v in org_rna_modification.items() if v}

	def update_lipid_modifications_from_homology(self):
		ref_lipid_modifications = self.ref.lipid_modifications
		org_lipid_modifications = self.org.lipid_modifications
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_lipid_modifications.items(),
					'Updating lipid modification machinery from homology...',
					bar_format = bar_format,
					total=len(ref_lipid_modifications)):
			if k not in org_lipid_modifications:
				org_lipid_modifications[k] = []
			for i in v:
				if i not in ref_cplx_homolog: continue
				org_cplx = ref_cplx_homolog[i]
# 				if org_cplx not in org_lipid_modifications[k]:
				if self._is_base_complex_in_list(org_cplx,org_lipid_modifications[k]):
					continue
				org_lipid_modifications[k].append(org_cplx)

	def update_transcription_subreactions_from_homology(self):
		ref_transcription_subreactions = self.ref.transcription_subreactions
		org_transcription_subreactions = self.org.transcription_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_transcription_subreactions.items(),
					'Updating transcription subreactions machinery from homology...',
					bar_format = bar_format,
					total=len(ref_transcription_subreactions)):
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_transcription_subreactions[k]["enzymes"]
			if defined_cplxs:
				continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if self._is_base_complex_in_list(i,defined_cplxs):
					continue
				defined_cplxs.append(i)

	def update_translocation_pathways_from_homology(self):
		ref_translocation_pathways = self.ref.translocation_pathways
		org_translocation_pathways = self.org.translocation_pathways
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in tqdm.tqdm(ref_translocation_pathways.items(),
					'Updating translocation machinery from homology...',
					bar_format = bar_format,
					total=len(ref_translocation_pathways)):
			ref_cplxs = v["enzymes"]
			if k not in org_translocation_pathways:
				org_translocation_pathways[k] = {"enzymes":[]}#v.copy()
			defined_cplxs = org_translocation_pathways[k]["enzymes"]
			if defined_cplxs:
				continue
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if self._is_base_complex_in_list(i,defined_cplxs):
					continue
				defined_cplxs.append(i)

	def update_m_model(self):
		org_model = self.org.m_model
		ref_model = self.ref.m_model

		for m in tqdm.tqdm(org_model.metabolites,
					'Fixing M-model metabolites with homology...',
					bar_format = bar_format):
			if m.id not in ref_model.metabolites:
				continue
			ref_m = ref_model.metabolites.get_by_id(m.id)
			if not m.formula:
				m.formula = ref_m.formula

	def update_subreaction_matrix(self):
		org_subreaction_matrix = self.org.subreaction_matrix
		if not org_subreaction_matrix.empty:return
		ref_subreaction_matrix = self.ref.subreaction_matrix
		org_model = self.org.m_model
		ref_model = self.ref.m_model
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_mets = []
		warn_cplxs = []
		for subrxn,row in tqdm.tqdm(ref_subreaction_matrix.iterrows(),
					'Updating subreaction matrix with homology...',
					bar_format = bar_format,
					total=ref_subreaction_matrix.shape[0]):
			d = {}
			d[subrxn] = {'Metabolites':'','Stoichiometry':''}
			met = row['Metabolites'].split('_mod_')[0]
			if met in self.ref.complexes_df.index:
				if met in ref_cplx_homolog:
					if 'mod' in row['Metabolites']:
						mods = '_mod_' + '_mod_'.join(row['Metabolites'].split('_mod_')[1:])
					else:
						mods = ''
					d[subrxn]['Metabolites'] = ref_cplx_homolog[met] + mods
				else:
					warn_cplxs = []
					d[subrxn]['Metabolites'] = 'CPLX_dummy'
				d[subrxn]['Stoichiometry'] = row['Stoichiometry']
			else:
				if not org_model.metabolites.has_id(met):
					warn_mets.append(met)
				d[subrxn]['Metabolites'] = met
				d[subrxn]['Stoichiometry'] = row['Stoichiometry']
			org_subreaction_matrix = \
				pandas.concat([org_subreaction_matrix,
							  pandas.DataFrame.from_dict(d).T],
							 axis = 0, join = 'outer')
		self.org.subreaction_matrix = org_subreaction_matrix
		self.org.subreaction_matrix.index.name = 'Reaction'
# 		self.org.subreaction_matrix.to_csv(self.org.directory + 'subreaction_matrix.txt')

		if warn_mets:
			self.curation_notes['update_subreaction_matrix'].append({
				'msg':'Some metabolites in subreaction_matrix were added from reference but are not in M-model',
				'triggered_by':warn_mets,
				'importance':'high',
				'to_do':'Map these metabolites or replace the subreaction'})
		if warn_cplxs:
			self.curation_notes['update_subreaction_matrix'].append({
				'msg':'Some complexes in subreaction_matrix of reference could not be mapped',
				'triggered_by':warn_cplxs,
				'importance':'high',
				'to_do':'Map these complexes or replace the subreaction'})

	def update_from_homology(self):
		self.update_enzyme_stoichiometry()
		self.update_protein_modification()
		self.update_TU_df()
		self.update_translocation_pathways_from_homology()
		self.protein_location_from_homology()
		self.update_translocation_multipliers()
		self.update_lipoprotein_precursors()
		self.update_cleaved_methionine()
		self.update_me_mets()
		self.update_generics_from_homology()
		self.update_folding_dict_from_homology()
		self.update_ribosome_subreactions_from_homology()
# 		self.update_rrna_modifications_from_homology()
		self.update_amino_acid_trna_synthetases_from_homology()
		self.update_peptide_release_factors_from_homology()
		self.update_transcription_subreactions_from_homology()
		self.update_initiation_subreactions_from_homology()
		self.update_elongation_subreactions_from_homology()
		self.update_termination_subreactions_from_homology()
		self.update_special_trna_subreactions_from_homology()
		self.update_rna_degradosome_from_homology()
		self.update_excision_machinery_from_homology()
		self.update_special_modifications_from_homology()
		self.update_rna_modification_from_homology()
		self.update_lipid_modifications_from_homology()
		self.update_m_model()
		self.update_subreaction_matrix()

	def fill(self, fill_with='CPLX_dummy'):
		coralme.builder.troubleshooting.fill_builder(self,fill_with='CPLX_dummy')

	def check(self):
		t_pathways = self.org.protein_location['translocase_pathway'].unique()

		file_t_paths = set()
		for t in t_pathways:
			for i in t:
				file_t_paths.add(i)

		defined_t_paths = set()
		for t in tqdm.tqdm(self.org.translocation_pathways.keys(),
					'Checking defined translocation pathways...',
					bar_format = bar_format):
			defined_t_paths.add(coralme.builder.dictionaries.pathway_to_abbreviation[t])
		missing_pathways = file_t_paths - defined_t_paths
		if missing_pathways:
			self.org.curation_notes['check'].append({
									'msg':'Some translocase pathways in org.protein_location are not defined in org.translocation_pathways.',
									'triggered_by':list(missing_pathways),
									'importance':'high',
									'to_do':'Fill in translocation pathways in org.translocation_pathways or in translocation_pathways.txt'
				})

		me_mets = self.org.me_mets
		warn_mets = list(me_mets[me_mets['type'] == 'CURATE'].index)
		# Warnings
		if warn_mets:
			self.org.curation_notes['check'].append({
				'msg':'Some metabolites in me_metabolites.txt need curation',
				'triggered_by':warn_mets,
				'importance':'medium',
				'to_do':'Map or remove these metabolites in me_metabolites.txt'})

	def load(self, directory):
		with open(directory, "rb") as f:
			tmp = pickle.load(f)
			return tmp

	def save(self, directory=False):
		if not directory:
			directory = self.org.directory + "builder.pickle"
		with open(directory, "wb") as f:
			pickle.dump(self, f)

# 	def load_me(self,filename='me_model.pickle'):
# 		with open(self.org.directory + '/'+filename, "rb") as f:
# 			return pickle.load(f)

	def save_builder_info(self):
		include = [float,int,str,pandas.DataFrame,dict]
		floats = {}
		dataframes = {}
		for i in dir(self.org):
			if i[0] == '_':
				continue
			attr = getattr(self.org,i)
			for c in include:
				if isinstance(attr,c):
					if c == float:
						floats[i] = attr
					elif c == dict:
						dataframes[i] = pandas.DataFrame.from_dict({'value':attr})
					elif c == pandas.DataFrame:
						dataframes[i] = attr
					break
			dataframes['parameters'] = pandas.DataFrame.from_dict({'value':floats})

		directory = self.org.directory + 'builder_info/'
		if not os.path.exists(directory):
			os.mkdir(directory)
		for k,v in dataframes.items():
			v.to_csv(directory + k + '.txt')

	# shortcuts to methods in the MECurator class
	def find_issue(self,query):
		coralme.builder.curation.MECurator(self.org).find_issue_with_query(query)

	# shortcuts to methods in the MEReconstruction and METroubleshooter classes
	def build_me_model(self, update = True, prune = True, overwrite = False, skip = None):
		coralme.builder.main.MEReconstruction(self).build_me_model(update = update, prune = prune, overwrite = overwrite, skip = skip)

	def troubleshoot(self, growth_key_and_value = None, skip = set(), guesses = set(), met_types = set(), platform = None, solver = 'qminos', savefile = None, gapfill_cofactors = False):
		"""
		growth_key_and_value:
			dictionary of Sympy.Symbol and value to replace

		skip:
			set of ME-components to not evaluate during gapfilling

		guesses:
			set of ME-components to try first before any other set of components

		platform:
			'win32' to use gurobi (default) or cplex as solver

		solver:
			'gurobi' (default, if platform is 'win32') or 'cplex'

		savefile:
			file path (absolute or relative) to save the ME-model as a pickle file
		"""

		coralme.builder.main.METroubleshooter(self).troubleshoot(growth_key_and_value, skip = skip, guesses = guesses, platform = platform, solver = solver, savefile = savefile,gapfill_cofactors=gapfill_cofactors)
		coralme.builder.notes.save_curation_notes(
				self.curation_notes,
				self.configuration['out_directory'] + '/curation_notes.json'
			)
		coralme.builder.notes.publish_curation_notes(
				self.curation_notes,
				self.configuration['out_directory']+ '/curation_notes.txt'
			)

	def input_data(self, gem, overwrite):
		tmp1, tmp2 = coralme.builder.main.MEReconstruction(self).input_data(gem, overwrite)
		self.df_tus, self.df_rmsc, self.df_subs, self.df_mets, self.df_keffs = tmp1
		self.df_data, self.df_rxns, self.df_cplxs, self.df_ptms, self.df_enz2rxn, self.df_rna_mods, self.df_protloc, self.df_transpaths = tmp2
		return tmp1, tmp2

class MEReconstruction(MEBuilder):
	"""
	MEReconstruction class for reconstructing a ME-model from user/automated input

	Parameters
	----------
	MEBuilder : coralme.builder.main.MEBuilder

	"""
	def __init__(self, builder):
		# only if builder.generate_files() was run before builder.build_me_model()
		if hasattr(builder, 'org'):
			self.org = builder.org
		if hasattr(builder, 'ref'):
			self.ref = builder.ref
		if hasattr(builder, 'homology'):
			self.homology = builder.homology
		if hasattr(builder, 'df_data'):
			self.df_data = builder.df_data
			self.df_tus = builder.df_tus
			self.df_rmsc = builder.df_rmsc
			self.df_subs = builder.df_subs
			self.df_mets = builder.df_mets
			self.df_keffs = builder.df_keffs
			self.df_rxns = builder.df_rxns
			self.df_cplxs = builder.df_cplxs
			self.df_ptms = builder.df_ptms
			self.df_enz2rxn = builder.df_enz2rxn
			self.df_rna_mods = builder.df_rna_mods
			self.df_protloc = builder.df_protloc
			self.df_transpaths = builder.df_transpaths
		if hasattr(builder, 'debug_data'):
			self.debug_data = builder.debug_data

		self.logger = builder.logger
		self.configuration = builder.configuration
		# reboot me_model
		if len(builder.me_model.reactions) == 1 and len(builder.me_model.metabolites):
			self.me_model = builder.me_model
		else:
			self.me_model = coralme.core.model.MEModel(self.configuration.get('ME-Model-ID', 'coralME'), self.configuration.get('growth_key', 'mu'))
		self.curation_notes = builder.curation_notes

	def input_data(self, m_model, overwrite = False):
		if hasattr(self, 'df_data'):
			return (self.df_tus, self.df_rmsc, self.df_subs, self.df_mets, self.df_keffs), (self.df_data, self.df_rxns, self.df_cplxs, self.df_ptms, self.df_enz2rxn, self.df_rna_mods, self.df_protloc, self.df_transpaths)

		config = self.configuration

		# Inferred information
		if hasattr(self, 'org'):
			config['selenocysteine_enzymes'] = self.org.special_trna_subreactions['sec_addition_at_UGA']['enzymes']
			logging.warning('The selenocysteine complex SelAB was set from homology data.')

			config['pg_pe_160'] = self.org.lipid_modifications.get('pg_pe_160', 'CPLX_dummy')
			logging.warning('The prolipoprotein diacylglyceryl transferase and the signal peptidase homologs were set from homology data.')

			config['other_lipids'] = self.org.lipid_modifications.get('other_lipids', 'CPLX_dummy')
			logging.warning('The apolipoprotein N-acyltransferase homolog was set from homology data.')

			lst = self.org.generic_dict.get('generic_fes_transfers_complex', {'enzymes' : ['CPLX_dummy']})['enzymes'] # ecolime = ['CPLX0-7617', 'CPLX0-7824', 'IscA_tetra']
# 			if len(lst) != 3:
# 				lst = lst + ['CPLX_dummy'] * (3 - len(lst))
			config['complex_cofactors']['fes_transfers'] = lst
			logging.warning('The iron-sulfur cluster insertion homologs were set from homology data.')

			self.org.get_reaction_keffs() # saves a file to <out_directory>/building_data/reaction_median_keffs.txt
# 			config['df_reaction_keff_consts'] = config.get('out_directory', '.') + '/building_data/reaction_median_keffs.txt'

		# include rna_polymerases, lipids and lipoproteins from automated info and save new configuration file
		if config.get('rna_polymerases', None) is None or config.get('rna_polymerases') == {}:
			if hasattr(self, 'org'):
				config['default_sigma_factor'] = self.org.rpod
				config['rna_polymerases'] = self.org.rna_polymerase_id_by_sigma_factor
				logging.warning('The RNA Polymerases (core enzyme and sigma factors) information was set from homology data.')

			## replace IDs
			#for name, rnap in config['rna_polymerases'].items():
				#if hasattr(self, 'homology'):
					#for key, value in rnap.items():
						#config['rna_polymerases'][name][key] = self.homology.org_cplx_homolog.get(value, value)
						#config['rna_polymerases'][name][key] = self.homology.org_cplx_homolog.get(value, value.replace('-MONOMER', '_MONOMER'))

		if config.get('lipid_modifications', None) is None or len(config.get('lipid_modifications')) == 0:
			if hasattr(self, 'org'):
				config['lipid_modifications'] = [ x for x in self.org.lipids if x.endswith('_p') and x.startswith('pg') and not x.startswith('pgp') ]
				logging.warning('The lipid modifications were set from M-model metabolites.')

		if config.get('lipoprotein_precursors', None) is None or len(config.get('lipoprotein_precursors')) == 0:
			if hasattr(self, 'org'):
				config['lipoprotein_precursors'] = self.org.lipoprotein_precursors
				logging.warning('The lipoprotein precursors were set from homology data.')

		if config.get('ngam', None) is None:
			if hasattr(self, 'org'):
				config['ngam'] = self.org.NGAM
				logging.warning('The ATPM value (ATP requirement for maintenance; also NGAM) was set from the M-model.')

		if config.get('gam', None) is None:
			if hasattr(self, 'org'):
				config['gam'] = self.org.GAM
				logging.warning('The GAM value (ATP requirement for growth) was set from the M-model or default value.')

		# modify options
		#config['create_files'] = False
		#config['run_bbh_blast'] = False
		#config['dev_reference'] = False

		if hasattr(self, 'org') and len(config.get('translocation_multipliers', {})) == 0:
			config['translocation_multipliers'] = { k:{ k:v for k,v in v.items() if v != 0 } for k,v in self.org.translocation_multipliers.items() }
			logging.warning('The translocation multipliers for yidC and tat homologs were set from homology data.')

		if hasattr(self, 'org') and len(config.get('amino_acid_trna_synthetase', {})) == 0:
			config['amino_acid_trna_synthetase'] = self.org.amino_acid_trna_synthetase
			logging.warning('The tRNA synthetases were set from homology data.')

		if hasattr(self, 'org') and len(config.get('defer_to_rxn_matrix', [])) == 0:
			config['defer_to_rxn_matrix'] = [self.org.biomass] if self.org.biomass is not None else []
			logging.warning('The biomass reaction will be skipped during the ME reconstruction steps.')

		if hasattr(self, 'org') and len(config.get('peptide_release_factors', [])) == 0:
			config['peptide_release_factors'] = { k:v['enzyme'] for k,v in self.org.peptide_release_factors.items() }
			logging.warning('The peptide release factors were set from homology data.')

		if not 'FMETTRS' in config.get('defer_to_rxn_matrix', []):
			config['defer_to_rxn_matrix'].append('FMETTRS')
			logging.warning('The FMETTRS reaction from the M-model will be replaced by a SubReaction during the ME-model reconstruction steps.')

		if not 'ATPM' in config.get('defer_to_rxn_matrix', []):
			config['defer_to_rxn_matrix'].append('ATPM')
			logging.warning('The ATPM reaction from the M-model will be replaced by a SummaryVariable during the ME-model reconstruction steps.')

		if hasattr(self, 'org') and len(config.get('braun\'s_lipoproteins', [])) == 0:
			lst = [ k.split('_mod_')[0] for k,v in self.org.protein_mod.to_dict()['Modifications'].items() if 'palmitate' in v ]
			config['braun\'s_lipoproteins'] = lst if isinstance(lst, list) else [lst]
			if len(lst) != 0:
				logging.warning('The Braun\'s lipoprotein homologs list was set to \'{:s}\'.'.format(', '.join(lst)))

		def read(filecode, input_type, default, columns = []):
			filename = config.get(filecode, '')
			if pathlib.Path(filename).is_file():
				file_to_read = filename
			else:
				#return pandas.DataFrame(columns = columns)
				if filecode == 'df_reaction_keff_consts':
					return pandas.DataFrame(columns = columns)
				file_to_read = '{:}/building_data/{:s}'.format(config.get('out_directory', '.'), default)

			df = coralme.builder.flat_files.read(file_to_read)
			if set(columns).issubset(set(df.columns)):
				return df
			else:
				#logging.warning('Column names in \'{:s}\' does not comply default values: {:s}.'.format(filename, ', '.join(columns)))
				raise Exception('Column names in \'{:s}\' does not comply default values: {:s}.'.format(filename, ', '.join(columns)))

		# User inputs
		# Transcriptional Units
		cols = ['TU_id', 'replicon', 'genes', 'start', 'stop', 'tss', 'strand', 'rho_dependent', 'rnapol']
		df_tus = read('df_TranscriptionalUnits', 'transcriptional units data', 'TUs_from_biocyc.txt', cols)
		df_tus = df_tus.set_index('TU_id', inplace = False)

		# Reaction Matrix: reactions, metabolites, compartments, stoichiometric coefficients
		cols = ['Reaction', 'Metabolites', 'Stoichiometry']
		df_rmsc = read('df_matrix_stoichiometry', 'reaction stoichiometry data', 'reaction_matrix.txt', cols)

		# SubReaction Matrix: subreactions, metabolites, compartments, stoichiometric coefficients
		cols = ['Reaction', 'Metabolites', 'Stoichiometry']
		df_subs = read('df_matrix_subrxn_stoich', 'subreaction stoichiometry data', 'subreaction_matrix.txt', cols)

		# Orphan and Spontaneous reaction metadata
		cols = ['name', 'description', 'subsystems', 'is_reversible', 'is_spontaneous']
		df_rxns = read('df_metadata_orphan_rxns', 'new reactions metadata', 'orphan_and_spont_reactions.txt', cols)
		df_rxns = df_rxns.set_index('name', inplace = False)

		# Metabolites metadata
		cols = ['id', 'me_id', 'name', 'formula', 'compartment', 'type']
		df_mets = read('df_metadata_metabolites', 'new metabolites metadata', 'me_metabolites.txt', cols)
		df_mets = df_mets.set_index('id', inplace = False)

		# Effective turnover rates
		cols = ['reaction', 'direction', 'complex', 'mods', 'keff']
		df_keffs = read('df_reaction_keff_consts', 'effective turnover rates', 'reaction_median_keffs.txt', cols)

		# set new options in the MEBuilder object
		self.configuration.update(config)

		# detect if the genbank file was modified using biocyc data
		gb = '{:s}/building_data/genome_modified.gb'.format(config.get('out_directory', '.'))
		config['genbank-path'] = gb if pathlib.Path(gb).exists() else config['genbank-path']

		if overwrite:
			new = config.get('new_config_file', 'coralme-config.yaml')
			yaml = new if new.endswith('.yaml') else '{:s}.yaml'.format(new)
			with open('{:s}/{:s}'.format(config.get('out_directory', '.'), yaml), 'w') as outfile:
				anyconfig.dump(config, outfile)
			logging.warning('New configuration file \'{:s}\' was written with inferred options.'.format(yaml))
			#with open('{:s}/{:s}'.format(config.get('out_directory', '.'), new), 'w') as outfile:
				#anyconfig.dump(config, outfile)
			#with open('{:s}/{:s}'.format(config.get('out_directory', '.'), new), 'w') as outfile:
				#anyconfig.dump(config, outfile)

		# Drop-in replacement of input files:
		# step1: reactions, complexes, modification of complexes, enzyme-to-reaction mapping
		# step2a: generics, dnap stoichiometry, ribosome stoichiometry, degradosome stoichiometry, tRNA ligases, RNA modifications
		# step2b: folding pathways (DnaK, GroEL), N-terminal Methionine Processing, translocation pathways
		filename = config.get('df_gene_cplxs_mods_rxns', '')
		if overwrite:
			try:
				pathlib.Path(filename).unlink(missing_ok = True) # python>=3.8
			except:
				if pathlib.Path(filename).exists():
					pathlib.Path(filename).unlink() # python==3.7

		if pathlib.Path(filename).is_file() and filename.endswith('.xlsx'):
			df_data = pandas.read_excel(filename, dtype = str).dropna(how = 'all')
		elif pathlib.Path(filename).is_file() and filename.endswith('.txt'):
			df_data = pandas.read_csv(filename, sep = '\t', header = 0, dtype = str).dropna(how = 'all')
		else:
			# detect if the genbank file was modified using biocyc data
			gb = '{:s}/building_data/genome_modified.gb'.format(config.get('out_directory', '.'))
			gb = gb if pathlib.Path(gb).exists() else config['genbank-path']
			ListHandler.print_and_log('Writting the Organism-Specific Matrix...')
			# generate a minimal dataframe from the genbank and m-model files
			df_data = coralme.builder.preprocess_inputs.generate_organism_specific_matrix(gb, config.get('locus_tag', 'locus_tag'), model = m_model)
			# complete minimal dataframe with automated info from homology
			df_data = coralme.builder.preprocess_inputs.complete_organism_specific_matrix(self, df_data, model = m_model, output = filename)
			ListHandler.print_and_log('Organism-Specific Matrix saved to {:s} file.'.format(filename))

		# All other inputs and remove unnecessary genes from df_data
		return (df_tus, df_rmsc, df_subs, df_mets, df_keffs), coralme.builder.preprocess_inputs.get_df_input_from_excel(df_data, df_rxns)

	def build_me_model(self, update = True, prune = True, overwrite = False, skip = None):
		"""Performs the Build step of the reconstruction.

		This function will read the synchronized and complemented information
		in the OSM and build a ME-model.

		Parameters
		----------
		update : bool
			If True, runs the update method of all reactions after building.
		prune : bool
			If True, prunes unused reactions and metabolites after building.
		overwrite : bool
			If True, overwrites the OSM and other configuration files.

		"""
		config = self.configuration
		model = config.get('ME-Model-ID', 'coralME')
		out_directory = config.get('out_directory', '.')
		if not os.path.exists(out_directory):
			os.mkdir(out_directory)

		log_directory = config.get('log_directory', '.')
		if not os.path.exists(log_directory):
			os.mkdir(log_directory)

		# ## Part 1: Create a minimum solvable ME-model
		# set logger
		log = logging.getLogger() # root logger
		for hdlr in log.handlers[:]: # remove all old handlers
			log.removeHandler(hdlr)

		# Old code works in a separate script; but it works if we remove the old handler
		logging.basicConfig(filename = '{:s}/MEReconstruction-step1-{:s}.log'.format(log_directory, model), filemode = 'w', level = logging.WARNING, format = log_format)
		log.addHandler(self.logger['MEReconstruction-step1'])
		#log.addHandler(logging.StreamHandler(sys.stdout))
		logging.captureWarnings(True)

		ListHandler.print_and_log("Initiating ME-model reconstruction...")

		# This will include the bare minimum representations of that still produce a working ME-model:
		# - Metabolic Reactions
		# - Transcription Reactions
		# - Translation Reactions
		# - Complex Formation Reactions
		#
		# ### 1) Create a MEModel object and populate its global info
		# This includes important parameters that are used to calculate coupling constraints as well as organism-specific information such as compartments, GC fraction, etc.

		me = self.me_model

		# update default options with user-defined values
		me.global_info.update(self.configuration)

		# Define the types of biomass components that will be synthesized in the model
		me.add_biomass_constraints_to_model(me.global_info['biomass_constraints'])

		# Define ME-model compartments
		me._compartments = me.global_info.get('compartments', {})
		if 'mc' not in me._compartments.keys():
			me._compartments['mc'] = 'ME-model Constraint'
			logging.warning('Pseudo-compartment \'mc\' (\'ME-model Constraint\') was added into the ME-model.')

		# Define M-model
		if hasattr(self, 'org'):
			me.gem = self.org.m_model
		else:
			if me.global_info['m-model-path'].endswith('.json'):
				me.gem = cobra.io.load_json_model(me.global_info['m-model-path'])
			else:
				me.gem = cobra.io.read_sbml_model(me.global_info['m-model-path'])

		# remove unused genes, reactions, and metabolites
		cobra.manipulation.delete.remove_genes(me.gem, [ x for x in me.gem.genes if len(x.reactions) == 0 ], remove_reactions = False)
		coralme.core.model.MEModel.prune_unused_reactions(me.gem) # reactions without metabolites
		coralme.core.model.MEModel.prune_unused_metabolites(me.gem) # orphan metabolites

		# update default options with missing, automated-defined values
		me.global_info.update(self.configuration)

		# Read user inputs
		tmp1, tmp2 = coralme.builder.main.MEReconstruction.input_data(self, me.gem, overwrite)
		(df_tus, df_rmsc, df_subs, df_mets, df_keffs), (df_data, df_rxns, df_cplxs, df_ptms, df_enz2rxn, df_rna_mods, df_protloc, df_transpaths) = tmp1, tmp2

		#TODO: Iterate tmp1 and tmp2 and warn empty dataframes

		me.internal_data = {}
		for key in ['df_tus', 'df_rmsc', 'df_subs', 'df_mets', 'df_keffs', 'df_data', 'df_rxns', 'df_cplxs', 'df_ptms', 'df_enz2rxn', 'df_rna_mods', 'df_protloc', 'df_transpaths']:
			exec('me.internal_data[\'{:s}\'] = {:s}'.format(key, key))

		# Remove default ME-model SubReactions from global_info that are not mapped in the organism-specific matrix
		subrxns = set(df_data[df_data['ME-model SubReaction'].notnull()]['ME-model SubReaction'])

		# list of subreactions
		for key in ['peptide_processing_subreactions']:
			me.global_info[key] = list(set(me.global_info[key]).intersection(subrxns))

		# dictionary mapping subreactions and stoichiometry
		for key in ['transcription_subreactions', 'translation_subreactions']:
			me.global_info[key] = { k:v for k,v in me.global_info[key].items() if k in subrxns }

		# ### 2) Load metabolites and build Metabolic reactions
		#
		# It creates a new M-model, then incorporates it into the ME-model using the *add_m_model_content* function. Reactions are added directly and metabolites are added as *StoichiometricData* (*me.stoichiometric_data*).
		#
		# Different metabolite types have different properties in a ME-model, so complexes are added to the model as a *ComplexData*, not as a *Metabolite*. Components in the M-model that are actually *Complexes* are compiled in the *cplx_lst* variable.

		# Modify M-model
		m_model = coralme.builder.flat_files.process_m_model(
			m_model = me.gem,
			rxns_data = df_rxns, # metadata of new reactions
			mets_data = df_mets, # metadata of new metabolites
			cplx_data = df_cplxs, # metadata of prot-prot/prot-rna complexes
			reaction_matrix = df_rmsc, # stoichiometric coefficients of new! reactions

			me_compartments = { v:k for k,v in me._compartments.items() }, repair = True,
			defer_to_rxn_matrix = me.global_info['defer_to_rxn_matrix'])

		me.processed_m_model = m_model

		# Some of the 'metabolites' in the M-model are actually complexes.
		# We pass those in so they get created as complexes, not metabolites.
		cplx_dct = coralme.builder.flat_files.get_complex_subunit_stoichiometry(df_cplxs)
		cplx_lst = [ i.id for i in m_model.metabolites if i.id.split('_mod_')[0] in cplx_dct.keys() ]

		# Complexes in `cplx_lst` are added as a coralme.core.component.Complex object
		# Other `metabolites` are added as a coralme.core.component.Metabolite object
		coralme.util.building.add_m_model_content(me, m_model, complex_metabolite_ids = set(cplx_lst))

		# NEW! Add metabolic subreactions (e.g. 10-Formyltetrahydrofolate:L-methionyl-tRNA N-formyltransferase)
		rxn_to_cplx_dict = coralme.builder.flat_files.get_reaction_to_complex(m_model, df_enz2rxn)

		reaction_matrix_dict = coralme.builder.flat_files.process_reaction_matrix_dict(
			df_subs, cplx_data = df_cplxs, me_compartments = { v:k for k,v in me._compartments.items() })

		for subreaction, stoichiometry in reaction_matrix_dict.items():
			coralme.util.building.add_subreaction_data(me, subreaction, stoichiometry, rxn_to_cplx_dict[subreaction])

		# ### 3) Add Transcription and Translation reactions
		#
		# To construct the bare minimum components of a transcription and translation reactions.
		# For example, transcription reactions at this point include nucleotides and the synthesized RNAs.

		# RNA and protein names are prefixed in the ME-model following then the locus tag
		lst = set(df_data['Gene Locus ID'].str.replace('^protein_', '', regex = True).str.replace('^RNA_', '', regex = True).tolist())

		# detect if the genbank file was modified using biocyc data
		gb = '{:s}/building_data/genome_modified.gb'.format(out_directory)
		gb = gb if pathlib.Path(gb).exists() else config['genbank-path']

		coralme.util.building.build_reactions_from_genbank(
			#me_model = me, gb_filename = me.global_info['genbank-path'],
			me_model = me, gb_filename = gb,
			tu_frame = df_tus, genes_to_add = lst, update = True, verbose = True,
			feature_types = me.global_info['feature_types'],
			trna_misacylation = me.global_info['trna_misacylation'],
			genome_mods = me.global_info['genome_mods'],
			knockouts = me.global_info['knockouts'])

		# ### 4) Add in ComplexFormation reactions without modifications (for now)

		# RNA components different from tRNAs and from 5S, 16S and 23S rRNAs
		rna_components = set(me.global_info['rna_components']) # in order: RNase_P_RNA, SRP_RNA, 6S RNA

		# cplx_dct is a dictionary {complex_name: {locus_tag/generic_name/subcomplex_name: count}}
		cplx_dct = coralme.builder.flat_files.get_complex_subunit_stoichiometry(df_cplxs, rna_components)

		# mods_dct is a dictionary {complex_name_with_mods: {core_enzyme: complex_name, modifications: {cofactor: stoichiometry}}
		mods_dct = coralme.builder.flat_files.get_complex_modifications(
			reaction_matrix = df_rmsc,
			protein_complexes = df_cplxs,
			complex_mods = df_ptms,
			compartments = { v:k for k,v in me._compartments.items() })

		# Add complexes into the ME-model as coralme.core.processdata.ComplexData objects.
		# Modifications are added as SubReactions
		coralme.util.building.add_model_complexes(me, cplx_dct, mods_dct)

		# Remove modifications. They will be added back in later (See Building Step 2, subsection 3).
		for data in tqdm.tqdm(list(me.complex_data), 'Removing SubReactions from ComplexData...', bar_format = bar_format):
			data.subreactions = {}

		# Add ComplexFormation reactions for each of the ComplexData
		for data in tqdm.tqdm(list(me.complex_data), 'Adding ComplexFormation into the ME-model...', bar_format = bar_format):
			formation = data.formation
			if formation:
				formation.update()
			elif data.id != 'CPLX_dummy':
				data.create_complex_formation()
				logging.warning('Added ComplexFormation for \'{:s}\'.'.format(data.id))
			else:
				pass

		# ### 5) Add *GenericData* and its reactions
		# Multiple entities can perform the same role. To prevent a combinatorial explosion, we create "generic" versions of the components, where any of those entities can fill in.
		#
		# This step was originally in the second part of the builder process, and was moved here to be able to use generics in enzymes associated to metabolic reactions.

		generics = coralme.builder.preprocess_inputs.get_generics(df_data)
		for generic, components in tqdm.tqdm(generics, 'Adding Generic(s) into the ME-model...', bar_format = bar_format):
			if 'generic_fes_transfers_complex' in generic:
				continue
			coralme.core.processdata.GenericData(generic, me, components).create_reactions()

		# ### 6) Add dummy reactions to model and the *unmodeled_protein_fraction* constraint
		#
		# This includes the Transcription, Translation, ComplexFormation, and Metabolic reactions for a dummy RNA/protein/complex. Sequence for *dummy RNA* is based on the prevalence of each codon found in the genbank file.

		coralme.util.building.add_dummy_reactions(me, me.global_info['transl_tables']['c'], update = True)

		# The 'dummy protein' is associated to orphan reactions.
		# This ensures that orphan reactions will not become favored to fulfil unmodeled protein fraction requirement.
		rxn = coralme.core.reaction.SummaryVariable('dummy_protein_to_mass')
		me.add_reactions([rxn])
		mass = me.metabolites.protein_dummy.formula_weight / 1000. # in kDa
		met = coralme.core.component.Constraint('unmodeled_protein_biomass')
		rxn.add_metabolites({'protein_biomass': -mass, 'protein_dummy': -1, met: mass})

		# WARNING: DO NOT REMOVE. This is needed if troubleshooter adds sink reactions for 4fe4s/2fe2s
		# Add CPLX_dummy_mod_2fe2s(1) and CPLX_dummy_mod_4fe4s(1) if fes_transfers is {'CPLX_dummy'}
		if 'CPLX_dummy' in me.global_info['complex_cofactors']['fes_transfers']:
			for fes in ['2fe2s', '4fe4s']:
				coralme.util.building.add_complex_to_model(me, 'CPLX_dummy_mod_{:s}(1)'.format(fes), { 'protein_dummy' : 1.0, fes + '_c': 1.0})
			# WARNING: do not change escape characters
			me.complex_data.query(r'CPLX_dummy_mod_2fe2s\(1\)')[0].create_complex_formation()
			me.complex_data.query(r'CPLX_dummy_mod_4fe4s\(1\)')[0].create_complex_formation()

		# ### 7) Associate Complexes to Metabolic reactions and build the ME-model metabolic network

		# Associate a reaction id with the ME-model complex id (including modifications)
		rxn_to_cplx_dict = coralme.builder.flat_files.get_reaction_to_complex(m_model, df_enz2rxn)
		spontaneous_rxns = [me.global_info['dummy_rxn_id']] + list(df_rxns[df_rxns['is_spontaneous'].str.match('1|TRUE|True|true')].index.values)
		# Correct rxn_to_cplx_dict with spontaneous reactions
		for key in spontaneous_rxns:
			if rxn_to_cplx_dict.get(key, None) is None:
				rxn_to_cplx_dict[key] = {None}
			else:
				rxn_to_cplx_dict[key].add(None)

		coralme.util.building.add_reactions_from_stoichiometric_data(
			me, rxn_to_cplx_dict, is_spontaneous = spontaneous_rxns, update = True)

		# correct the enzyme list of a subreaction if it matches the component_list of a generic
		# WARNING: This allows the use of generics in subreactions
		# e.g.: acpP_activation is associated to EG12221-MONOMER OR HOLO-ACP-SYNTH-CPLX_mod_mg2(1) OR HOLO-ACP-SYNTH-CPLX_mod_mn2(1)
		# e.g.: EG12221-MONOMER OR HOLO-ACP-SYNTH-CPLX_mod_mg2(1) OR HOLO-ACP-SYNTH-CPLX_mod_mn2(1) is associated to generic_acp_synthase
		for data in list(me.subreaction_data):
			for generic in list(me.generic_data):
				if data.enzyme == set(generic.component_list):
					data.enzyme = set([generic.id])

		# ### 8) Incorporate remaining biomass constituents
		# ### 1. General Demand Requirements
		# There are leftover components from the biomass equation that either:
		# 1. they have no mechanistic function in the model (e.g., *glycogen*)
		# 2. they are cofactors that are regenerated (e.g., *nad*)
		#
		# Applies demands and coefficients from the biomass objective function from the M-model.

		# set values from configuration as properties of the ME-model
		for key in [ 'gam', 'ngam', 'unmodeled_protein_fraction' ]:
			if key in me.global_info:
				setattr(me, key, me.global_info[key])

		biomass_constituents = me.global_info.get('flux_of_biomass_constituents', {})
		# replace IDs. New metabolites x types are created during the processing of the m_model
		# old ID (M-model) : new ID (ME-model)
		dct = df_mets[df_mets['type'].str.contains('REPLACE')].to_dict()['me_id']
		for key in list(biomass_constituents.keys()):
			if key in dct:
				biomass_constituents[dct[key]] = biomass_constituents.pop(key)
				logging.warning('Metabolite \'{:s}\' was replaced with \'{:s}\' in MetabolicReaction \'{:s}\'.'.format(key, dct[key], 'biomass_constituent_demand'))

		# remove metabolites not in the model or without molecular weight
		biomass_constituents = { k:v for k,v in biomass_constituents.items() if me.metabolites.has_id(k) and me.metabolites.get_by_id(k).formula_weight }

		problems = list(set(me.global_info.get('flux_of_biomass_constituents', {})).difference(biomass_constituents))
		if problems:
			logging.warning('The following biomass constituents are not in the ME-model or have no formula: {:s}.'.format(', '.join(problems)))
			logging.warning('A second attempt to add biomass constituents will be perform after update of formulas.')

		rxn = coralme.core.reaction.SummaryVariable('biomass_constituent_demand')
		me.add_reactions([rxn])
		rxn.add_metabolites({ k:-(abs(v)) for k,v in biomass_constituents.items() })
		rxn.lower_bound = me.mu # coralme.util.mu
		rxn.upper_bound = me.mu # coralme.util.mu
		constituent_mass = sum([me.metabolites.get_by_id(c).formula_weight / 1000. * abs(v) for c,v in biomass_constituents.items()])
		rxn.add_metabolites({me.metabolites.get_by_id('constituent_biomass'): constituent_mass})

		# ### 2. Lipid Demand Requirements
		# Metabolites and coefficients from biomass objective function

		lipid_demand = me.global_info.get('flux_of_lipid_constituents', {})
		#for key, value in me.global_info.get('flux_of_lipid_constituents', {}).items():
			#lipid_demand[key] = abs(value)

		if lipid_demand:
			for met, requirement in lipid_demand.items():
				try:
					component_mass = me.metabolites.get_by_id(met).formula_weight / 1000.
					rxn = coralme.core.reaction.SummaryVariable('DM_' + met)
					me.add_reactions([rxn])
					rxn.add_metabolites({met: -1 * abs(requirement), 'lipid_biomass': component_mass * abs(requirement)})
					rxn.lower_bound = me.mu # coralme.util.mu
					rxn.upper_bound = me.mu # originally 1000.
				except:
					msg = 'Metabolite \'{:s}\' lacks a formula. Please correct it in the M-model or the \'metabolites.txt\' metadata file.'
					logging.warning(msg.format(met))

		# ### 3. DNA Demand Requirements
		# Added based on growth rate dependent DNA levels as in [O'brien EJ et al 2013](https://www.ncbi.nlm.nih.gov/pubmed/24084808) (*E. coli* data)
		if me.global_info.get('flux_of_dna_constituents', False):
			# Fraction of each nucleotide in DNA, based on M-model's BOF
			dna_demand_stoich = {
				'datp_c': -abs(me.global_info['flux_of_dna_constituents']['datp_c']),
				'dctp_c': -abs(me.global_info['flux_of_dna_constituents']['dctp_c']),
				'dgtp_c': -abs(me.global_info['flux_of_dna_constituents']['dgtp_c']),
				'dttp_c': -abs(me.global_info['flux_of_dna_constituents']['dttp_c']),
				}
			dna_demand_stoich['ppi_c'] = abs(me.global_info['flux_of_dna_constituents'].get('ppi_c', sum([ v for k,v in dna_demand_stoich.items() if k != 'ppi_c'])))
		else:
			# Fraction of each nucleotide in DNA, based on gc_fraction
			dna_demand_stoich = {
				'datp_c': -((1 - me.global_info['GC_fraction']) / 2),
				'dctp_c': -(me.global_info['GC_fraction'] / 2),
				'dgtp_c': -(me.global_info['GC_fraction'] / 2),
				'dttp_c': -((1 - me.global_info['GC_fraction']) / 2),
				'ppi_c': 1
				}

		dna_replication = coralme.core.reaction.SummaryVariable('DNA_replication')
		me.add_reactions([dna_replication])
		dna_replication.add_metabolites(dna_demand_stoich)

		dna_mw = 0
		dna_mw_no_ppi = coralme.builder.dna_replication.get_dna_mw_no_ppi_dict(me)
		for met, value in me.reactions.get_by_id('DNA_replication').metabolites.items():
			if met.id != 'ppi_c':
				dna_mw -= value * dna_mw_no_ppi[met.id.replace('_c', '')] / 1000. # in kDa

		dna_biomass = coralme.core.component.Constraint('DNA_biomass')
		dna_replication.add_metabolites({dna_biomass: dna_mw})

		if me.global_info.get('flux_of_dna_constituents', False):
			dna_demand_bound = me.mu
		else:
			dna_demand_bound = coralme.builder.dna_replication.return_gr_dependent_dna_demand(
				me, me.global_info['GC_fraction'], me.global_info['percent_dna_data'], me.global_info['gr_data_doublings_per_hour'])

		dna_replication.lower_bound = dna_demand_bound
		dna_replication.upper_bound = dna_demand_bound

		# ### 9) Save ME-model as a pickle file
		with open('{:s}/MEModel-step1-{:s}.pkl'.format(out_directory, model), 'wb') as outfile:
			pickle.dump(me, outfile)

		ListHandler.print_and_log('ME-model was saved in the {:s} directory as MEModel-step1-{:s}.pkl'.format(out_directory, model))

		# ## Part 2: Add metastructures to solving ME-model
		# set logger
		log = logging.getLogger() # root logger
		for hdlr in log.handlers[:]: # remove all old handlers
			log.removeHandler(hdlr)

		# Old code works in a separate script; but it works if we remove the old handler
		logging.basicConfig(filename = '{:s}/MEReconstruction-step2-{:s}.log'.format(log_directory, model), filemode = 'w', level = logging.WARNING, format = log_format)
		log.addHandler(self.logger['MEReconstruction-step2'])
		#log.addHandler(logging.StreamHandler(sys.stdout))
		logging.captureWarnings(True)

		# ### 1) Add *GenericData* and its reactions
		# Multiple entities can perform the same role. To prevent a combinatorial explosion, we create "generic" versions of the components, where any of those entities can fill in.
		#
		# This step was moved to the first part of the builder process to be able to use generics in enzymes linked to metabolic reactions.

		# ### 2) Add DNA polymerase and the ribosome and its rRNAs modifications
		# ~~This uses the ribosome composition and subreaction definitions in **coralme/ribosome.py**~~

		# WARNING: EXPERIMENTAL, not included in the original cobrame paper.
		# 1) Is the coupling coefficient dependent on the length of the DNA?
		# 2) Do we need to modify the DNA replication to have one reaction per replicon?
		# dnapol_id = me.global_info['dnapol_id']

		# me.add_metabolites([cobrame.core.component.Complex(dnapol_id)])
		# data = cobrame.core.processdata.ComplexData(dnapol_id, me)
		# data.stoichiometry.update(coralme.preprocess_inputs.dnapolymerase_stoichiometry(df_data))
		# data.create_complex_formation(verbose = False)

		ribosome_stoich = coralme.builder.preprocess_inputs.ribosome_stoichiometry(df_data)
		ribosome_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Ribosome')
		df_rrna_mods = df_rna_mods[df_rna_mods['bnum'].isin(['16s_rRNAs', '23s_rRNAs'])]
		coralme.builder.ribosome.add_ribosome(me, ribosome_stoich, ribosome_subreactions, df_rrna_mods, verbose = True)

		# ### 3) Add charged tRNA reactions

		# The tRNA charging reactions were automatically added when loading the genome from the GenBank file. However, the charging reactions still need to be made aware of the tRNA synthetases which are responsible. Generic charged tRNAs are added to translation reactions via *SubreactionData* below.

		# tRNA synthetases per organelle
		if hasattr(self, 'org') and len(me.global_info.get('amino_acid_trna_synthetase', {})) == 0:
			aa_synthetase_dict = me.global_info['amino_acid_trna_synthetase']
		elif len(me.global_info.get('amino_acid_trna_synthetase', {})) != 0:
			aa_synthetase_dict = me.global_info.get('amino_acid_trna_synthetase', {})
		else:
			aa_synthetase_dict = coralme.builder.preprocess_inputs.aa_synthetase_dict(df_data)
			logging.warning('Association of tRNA synthetases and amino acids was inferred from GenBank annotation. It can be incomplete.')

		for data in tqdm.tqdm(list(me.tRNA_data), 'Adding tRNA synthetase(s) information into the ME-model...', bar_format = bar_format):
			data.synthetase = str(aa_synthetase_dict.get(data.amino_acid, 'CPLX_dummy'))

		# Correct 'translation_stop_dict' if PrfA and/or PrfB homologs were not identified
		PrfA_mono = me.global_info['peptide_release_factors'].get('UAG', 'PrfA_mono')
		# In Trans Table 4, UGA encodes Trp, not a stop codon. User might choose not to add 'UGA' in peptide_release_factors in organism.json
		PrfB_mono = me.global_info['peptide_release_factors'].get('UGA', 'PrfB_mono')
		generic_RF = me.global_info['peptide_release_factors'].get('UAA', 'generic_RF')

		if     me.metabolites.has_id(PrfA_mono) and not me.metabolites.has_id(PrfB_mono):
			me.global_info['peptide_release_factors']['UGA'] = PrfA_mono # originally assigned to PrfB_mono
			me.global_info['peptide_release_factors']['UAA'] = PrfA_mono # originally assigned to generic_RF
		if not me.metabolites.has_id(PrfA_mono) and     me.metabolites.has_id(PrfB_mono):
			me.global_info['peptide_release_factors']['UAG'] = PrfB_mono # originally assigned to PrfA_mono
			me.global_info['peptide_release_factors']['UAA'] = PrfB_mono # originally assigned to generic_RF
		if not me.metabolites.has_id(PrfA_mono) and not me.metabolites.has_id(PrfB_mono):
			me.global_info['peptide_release_factors']['UAG'] = 'CPLX_dummy' # originally assigned to PrfA_mono
			me.global_info['peptide_release_factors']['UGA'] = 'CPLX_dummy' # originally assigned to PrfB_mono

		if me.global_info['peptide_release_factors']['UAG'] == 'CPLX_dummy' and me.global_info['peptide_release_factors']['UGA'] == 'CPLX_dummy':
			me.global_info['peptide_release_factors']['UAA'] = 'CPLX_dummy'

		# charged tRNAs
		for organelle, transl_table in me.global_info['transl_tables'].items():
			if len(transl_table) == 0:
				continue

			coralme.builder.translation.add_charged_trna_subreactions(me, organelle, transl_table, translation_stop_dict = me.global_info['peptide_release_factors'], selenocysteine_enzymes = me.global_info.get('selenocysteine_enzymes', []))

		# ### 4) Add tRNA modifications into the ME-model and associate them with tRNA charging reactions

		# Add tRNA modifications to ME-model per type of organism
		if me.global_info['domain'].lower() in ['prokaryote', 'bacteria']:
			df_trna_mods = df_rna_mods[~df_rna_mods['bnum'].isin(['16s_rRNAs', '23s_rRNAs'])]
		elif me.global_info['domain'].lower() in ['eukarya', 'eukaryote']:
			df_trna_mods = df_rna_mods[~df_rna_mods['bnum'].isin(['18S_rRNAs', '25S_rRNAs', '28S_rRNAs'])]
		else:
			logging.warning('The \'domain\' property is not valid. A valid value is \'Prokaryote\', \'Bacteria\', \'Eukarya\', or \'Eukaryote\'.')

		coralme.builder.trna_charging.add_trna_modification_procedures(me, df_trna_mods)
		# Associate tRNA modifications to tRNAs
		trna_modifications = coralme.builder.flat_files.get_trna_modification_targets(df_trna_mods)
		for idx, trna in tqdm.tqdm(list(df_trna_mods.iterrows()), 'Associating tRNA modification enzyme(s) to tRNA(s)...', bar_format = bar_format):
			for data in me.process_data.query('tRNA_' + trna['bnum']):
				data.subreactions = trna_modifications[trna['bnum']]

		# ### 5) Add translation SubReactions into the ME-model and associate them with Translation reactions

		initiation_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_initiation')
		elongation_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_elongation')
		termination_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_termination')
		#termination_subreactions.update(coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_termination_generic_RF_mediated'))
		processing_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Protein_processing')

		for data in tqdm.tqdm(list(me.translation_data), 'Adding SubReactions into TranslationReactions...', bar_format = bar_format):
			data.add_initiation_subreactions(start_codons = me.global_info['start_codons'], start_subreactions = initiation_subreactions)
			data.add_elongation_subreactions(elongation_subreactions = elongation_subreactions)
			data.add_termination_subreactions(translation_terminator_dict = me.global_info['peptide_release_factors'])

		# ### 6) Add Transcription Metacomplexes: RNA Polymerase(s)

		rna_polymerases = me.global_info.get('rna_polymerases', {})

		# Create polymerase "metabolites"
		for rnap, components in tqdm.tqdm(rna_polymerases.items(), 'Adding RNA Polymerase(s) into the ME-model...', bar_format = bar_format):
			if me.metabolites.has_id(components['sigma_factor']) and me.metabolites.has_id(components['polymerase']):
				rnap_obj = coralme.core.component.RNAP(rnap)
				me.add_metabolites(rnap_obj)
				logging.warning('The RNA Polymerase \'{:s}\' was created in the ME-model successfully.'.format(rnap))
			else:
				if not me.metabolites.has_id(components['sigma_factor']):
					logging.warning('The complex ID \'{:s}\' from \'rna_polymerases\' in the configuration does not exist in the organism-specific matrix. Please check if it is the correct behaviour.'.format(components['sigma_factor']))
				if not me.metabolites.has_id(components['polymerase']):
					logging.warning('The complex ID \'{:s}\' from \'rna_polymerases\' in the configuration does not exist in the organism-specific matrix. Please check if it is the correct behaviour.'.format(components['polymerase']))

		# Add polymerase complexes in the model
		coralme.builder.transcription.add_rna_polymerase_complexes(me, rna_polymerases, verbose = False)

		# Associate the correct RNA_polymerase and factors to TUs
		sigma_to_rnap = { v['sigma_factor']:k for k,v in rna_polymerases.items() }
		for tu_id in tqdm.tqdm(df_tus.index, 'Associating a RNA Polymerase to each Transcriptional Unit...', bar_format = bar_format):
			if me.process_data.has_id(tu_id):
				transcription_data = me.process_data.get_by_id(tu_id)
				transcription_data.RNA_polymerase = sigma_to_rnap.get(df_tus['rnapol'][tu_id], None)
			else:
				logging.warning('Transcription Unit \'{:s}\' is missing from ProcessData, likely the associated gene(s) (\'{:s}\') were filtered out before the reconstruction. Check if it is the correct behavior.'.format(tu_id, df_tus['genes'][tu_id]))
				pass

		# WARNING: Without a TUs file, the 'most common' polymerase should be an empty string
		rnap_counter = collections.Counter([ x.RNA_polymerase for x in me.transcription_data ])
		# override 'most_common' polymerase if it is empty
		user_default_sigma = me.global_info.get('default_sigma_factor', '')
		user_default_rnap = sigma_to_rnap.get(user_default_sigma, None)

		if user_default_rnap is None:
			msg = 'Assigning most common RNAP \'{:s}\' to missing polymerase in \'{:s}\''
			most_common = max(rnap_counter, key = rnap_counter.get)
		else:
			msg = 'Assigning user default RNAP \'{:s}\' to missing polymerase in \'{:s}\''
			most_common = user_default_rnap

		for transcription_data in me.transcription_data:
			if transcription_data.RNA_polymerase == '' and most_common != '':
				logging.warning(msg.format(most_common,transcription_data.id))
				transcription_data.RNA_polymerase = most_common

		# ### 7) Add Transcription Metacomplexes: Degradosome (both for RNA degradation and RNA splicing)
		degradosome_id = me.global_info['degradosome_id']

		me.add_metabolites([coralme.core.component.Complex(degradosome_id)])
		data = coralme.core.processdata.ComplexData(degradosome_id, me)
		stoich = coralme.builder.preprocess_inputs.degradosome_stoichiometry(df_data)
		if stoich is not None:
			data.stoichiometry.update(stoich)
		else:
			data.stoichiometry.update({'CPLX_dummy' : -1})
		data.create_complex_formation(verbose = False)

		# Used for RNA splicing
		data = coralme.core.processdata.SubreactionData('RNA_degradation_machine', me)
		data.enzyme = [me.global_info['degradosome_id']]

		# .25 water equivalent for ATP hydrolysis per nucleotide
		data = coralme.core.processdata.SubreactionData('RNA_degradation_atp_requirement', me)
		data.stoichiometry = { 'atp_c': -0.25, 'h2o_c': -0.25, 'adp_c': +0.25, 'h_c': +0.25, 'pi_c': +0.25 }

		for excision_type in me.global_info['excision_machinery']:
			stoichiometry = coralme.builder.preprocess_inputs.excision_machinery_stoichiometry(df_data, excision_type)
			if stoichiometry is not None:
				coralme.builder.transcription.add_rna_excision_machinery(me, excision_type, stoichiometry)
			else:
				coralme.builder.transcription.add_rna_excision_machinery(me, excision_type, {'CPLX_dummy' : +1}) # +1 is correct
				logging.warning('All the components of the excision complex for \'{:s}\' were not identified from homology and it was assigned to the \'CPLX_dummy\' complex.'.format(excision_type))

		# add excision machineries into TranscriptionData
		# WARNING: subreactions is now a property of the TranscriptionData recalculated when accessed
		#coralme.builder.transcription.add_rna_splicing(me)

		# ## Part 3: Add remaining modifications (including iron clusters and lipoate)

		## WARNING: This is calculated above (See Building Step 1, subsection 3)
		## mods_dct is a dictionary {complex_name_with_mods: {core_enzyme: complex_name, modifications: {stoichiometry}}
		#mods_dct = coralme.builder.flat_files.get_complex_modifications(
			#reaction_matrix = df_rmsc,
			#protein_complexes = df_cplxs,
			#complex_mods = df_ptms,
			#compartments = { v:k for k,v in me._compartments.items() })

		for complex_id, info in tqdm.tqdm(mods_dct.items(), 'Processing ComplexData in ME-model...', bar_format = bar_format):
			modifications = {}
			for mod, value in info['modifications'].items():
				# stoichiometry of modification determined in subreaction_data.stoichiometry
				modifications['mod_' + mod] = abs(value)
			me.process_data.get_by_id(complex_id).subreactions = modifications

		# Check if complex modifications are set on any component in the organism-specific matrix
		# Dipyrromethane
		# dpm modification never from the free metabolite
		if me.process_data.has_id('mod_dpm_c'):
			if me.metabolites.has_id('dpm_c'):
				me.remove_metabolites([me.metabolites.dpm_c])

		# WARNING: use a different ID for spontaneous modification vs enzymatic modification (?)
		# biotin from the free metabolite in malonate decarboxylase (EC 7.2.4.4); don't remove biotin from the model (EC 4.1.1.88 is biotin-independent)
		# biotin from the free metabolite in acetyl-CoA carboxylase, but using biotin---[acetyl-CoA-carboxylase] ligase
		if me.process_data.has_id('mod_btn_c'):
			coralme.builder.modifications.add_btn_modifications(me)

		# 2'-(5''-triphosphoribosyl)-3'-dephospho-CoA in CitD catalyzed by CitX
		# 2'-(5''-triphosphoribosyl)-3'-dephospho-CoA loaded from the free metabolite; don't remove it from the model
		if me.process_data.has_id('mod_2tpr3dpcoa_c'):
			coralme.builder.modifications.add_2tpr3dpcoa_modifications(me)

		# activation of glycyl radical enzymes
		# glycyl modification never from the free metabolite
		if me.process_data.has_id('mod_glycyl_c'):
			coralme.builder.modifications.add_glycyl_modifications(me)
			if me.metabolites.has_id('glycyl_c'):
				me.remove_metabolites([me.metabolites.glycyl_c])

		# pap4p in AcpP catalyzed by AcpS
		# pan4p loaded from free CoA, producing pap as byproduct; don't remove it from the model
		if me.process_data.has_id('mod_pan4p_c'):
			coralme.builder.modifications.add_pan4p_modifications(me)

		# https://www.genome.jp/pathway/map00785
		# lipoyl is a pseudo metabolite; it rather comes from free lipoate or from octanoate-ACP
		if me.process_data.has_id('mod_lipoyl_c'):
			coralme.builder.modifications.add_lipoyl_modifications(me)
			if me.metabolites.has_id('lipoyl_c'):
				me.remove_metabolites([me.metabolites.lipoyl_c])

		if me.process_data.has_id('mod_3fe4s_c') or me.process_data.has_id('mod_4fe4s_c') or me.process_data.has_id('mod_2fe2s_c'):
			coralme.builder.modifications.add_iron_sulfur_modifications(me)
		if me.process_data.has_id('mod_FeFe_cofactor_c') or me.process_data.has_id('mod_NiFe_cofactor_c'):
			coralme.builder.modifications.add_FeFe_and_NiFe_modifications(me)
		if me.process_data.has_id('mod_bmocogdp_c'):
			coralme.builder.modifications.add_bmocogdp_chaperones(me)

		# add formation reactions for each of the ComplexData
		for data in tqdm.tqdm(list(me.complex_data), 'Adding ComplexFormation into the ME-model...', bar_format = bar_format):
			formation = data.formation
			if formation:
				formation.update()
			else:
				data.create_complex_formation()
				logging.warning('Added a ComplexFormation reaction for \'{:s}\'.'.format(data.id))

		# ## Part 4: Add remaining subreactions
		# ### 1. Add translation related subreactions

		# get list of processed proteins from df_data
		processing_pathways = [
			x.replace('Protein_processing_', '')
			for x in me.global_info['translation_subreactions'].keys()
			if x.startswith('Protein_processing')
			]
		protein_processing = { k:df_data[~df_data[k].isnull() & df_data[k].str.match('1|TRUE|True|true')]['Gene Locus ID'].tolist() for k in processing_pathways }

		# Add the translation subreaction data objects to the ME-model
		coralme.builder.translation.add_subreactions_to_model(
			me, [initiation_subreactions, elongation_subreactions, termination_subreactions, processing_subreactions])

		for data in tqdm.tqdm(list(me.translation_data), 'Adding SubReactions into TranslationReactions...', bar_format = bar_format):
			for process in protein_processing:
				if data.id in protein_processing[process]:
					data.subreactions['Protein_processing_' + process] = 1.

			# This block was run above, but it should be run again to incorporate any subreactions not added previously
			data.add_initiation_subreactions(start_codons = me.global_info['start_codons'], start_subreactions = initiation_subreactions)
			data.add_elongation_subreactions(elongation_subreactions = elongation_subreactions)
			data.add_termination_subreactions(translation_terminator_dict = me.global_info['peptide_release_factors'])

			# Add organism specific subreactions associated with peptide processing
			for subrxn in me.global_info['peptide_processing_subreactions']:
				data.subreactions[subrxn] = 1.

		# ### 2) Add transcription related subreactions
		transcription_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Transcription')
		coralme.builder.transcription.add_subreactions_to_model(me, [transcription_subreactions])

		for transcription_data in tqdm.tqdm(list(me.transcription_data), 'Adding Transcription SubReactions...', bar_format = bar_format):
			# Assume false if not in tu_df
			is_rho_dependent = df_tus.rho_dependent.get(transcription_data.id, 'False')
			rho_dependency = 'dependent' if is_rho_dependent in ['1', 'TRUE', 'True', 'true'] else 'independent'
			stable = 'stable' if transcription_data.codes_stable_rna else 'normal'
			if 'Transcription_{:s}_rho_{:s}'.format(stable, rho_dependency) in me.global_info['transcription_subreactions']:
				transcription_data._subreactions['Transcription_{:s}_rho_{:s}'.format(stable, rho_dependency)] = 1.
			else:
				logging.warning('The SubReaction \'Transcription_{:s}_rho_{:s}\' is not defined in the organism-specific matrix.'.format(stable, rho_dependency))

		# ## Part 5: Add in Translocation reactions

		v1 = { 'fixed_keff' : False, 'length_dependent' : True } # default
		v2 = { 'fixed_keff' : True,  'length_dependent' : False } # only for FtsY in the SRP pathway
		v3 = { 'fixed_keff' : False, 'length_dependent' : False } # for all the enzymes from the tat, tat_alt, lol and bam pathways

		for key, value in me.global_info['translocation_pathway'].items():
			if 'translocation_pathway_' + key in df_transpaths.index:
				me.global_info['translocation_pathway'][key]['enzymes'] = {
					k:(v2 if k == value.get('FtsY', None) else v1 if (key.lower() not in ['tat', 'tat_alt', 'lol', 'bam']) else v3) \
						for k in df_transpaths.loc['translocation_pathway_' + key].tolist()[0] }

			# TO ADD PATHWAYS WITHOUT HOMOLOGS
			# Check if the user wants to add dummies to the translocation pathways
			elif bool(me.global_info.get('add_translocases', False)) and value.get('enzymes', None) is None:
				me.global_info['translocation_pathway'][key]['enzymes'] = { 'CPLX_dummy':(v2 if value.get('FtsY', None) else v1 if (key.lower() not in ['tat', 'tat_alt', 'lol', 'bam']) else v3) }
				logging.warning('The component \'CPLX_dummy\' was associated to translocation pathways without defined homologs.')

		dct = { k:v['abbrev'] for k,v in me.global_info['translocation_pathway'].items() }
		dct = dict([(v, [k + '_translocation' for k,v1 in dct.items() if v1 == v]) for v in set(dct.values())])
		dct = { k:(v[0] if len(v) == 1 else v) for k,v in dct.items() } # tat pathways should be a list, but not the others

		# Check if the user added the reactions for translocation data
		if not me.process_data.has_id('atp_hydrolysis_sec_pathway'):
			stoichiometry = {'atp_c': -0.04, 'h2o_c': -0.04, 'adp_c': +0.04, 'h_c': +0.04, 'pi_c' : +0.04}
			coralme.util.building.add_subreaction_data(
				me, modification_id = 'atp_hydrolysis_sec_pathway', modification_stoichiometry = stoichiometry, modification_enzyme = 'CPLX_dummy')
		if not me.process_data.has_id('atp_hydrolysis_secA'):
			stoichiometry = {'atp_c': -1/75, 'h2o_c': -1/75, 'adp_c': +1/75, 'h_c': +1/75, 'pi_c': +1/75}
			coralme.util.building.add_subreaction_data(
				me, modification_id = 'atp_hydrolysis_secA', modification_stoichiometry = stoichiometry, modification_enzyme = 'CPLX_dummy')
		if not me.process_data.has_id('gtp_hydrolysis_srp_pathway'):
			stoichiometry = {'gtp_c': -2.0, 'h2o_c': -2.0, 'gdp_c': +2.0, 'h_c': +2.0, 'pi_c': +2.0}
			coralme.util.building.add_subreaction_data(
				me, modification_id = 'gtp_hydrolysis_srp_pathway', modification_stoichiometry = stoichiometry, modification_enzyme = 'CPLX_dummy')

		# for pathway, info in coralme.translocation.pathway.items():
		for pathway, info in me.global_info['translocation_pathway'].items():
			if me.global_info['translocation_pathway'][pathway].get('enzymes', None) is not None:
				transloc_data = coralme.core.processdata.TranslocationData(pathway + '_translocation', me)
				transloc_data.enzyme_dict = info['enzymes']
				transloc_data.keff = info['keff']
				transloc_data.length_dependent_energy = info['length_dependent_energy']
				if me.process_data.has_id(info['stoichiometry']):
					transloc_data.stoichiometry = me.process_data.get_by_id(info['stoichiometry']).stoichiometry
				else:
					transloc_data.stoichiometry = {}

		# Associate data and add translocation reactions
		multipliers = collections.defaultdict(dict)
		for enzyme, value in me.global_info.get('translocation_multipliers', {}).items():
			for bnum in value.keys():
				multipliers['protein_' + bnum][enzyme] = value[bnum]

		coralme.builder.translocation.add_translocation_pathways(me, pathways_df = df_protloc, abbreviation_to_pathway = dct, multipliers = multipliers, membrane_constraints = False)

		# Update stoichiometry of membrane complexes
		new_stoich = collections.defaultdict(dict)

		for idx, row in tqdm.tqdm(list(df_protloc.iterrows()), 'Processing StoichiometricData in SubReactionData...', bar_format = bar_format):
			protID = row['Protein']
			protIDLoc = protID + '_' + row['Protein_compartment']
			# to use dict.update() later; It replaces removing the key using the dict.pop method
			new_stoich[idx]['protein_' + protID] = 0
			new_stoich[idx]['protein_' + protIDLoc] = float(me.process_data.get_by_id(idx).stoichiometry['protein_' + protID])

		for cplx, stoich in new_stoich.items():
			complex_data = me.process_data.get_by_id(cplx)
			complex_data.stoichiometry.update(new_stoich[cplx])
			# remove zeroes from complex_data.stoichiometry
			complex_data.stoichiometry = { k:v for k,v in complex_data.stoichiometry.items() if v != 0 }
			complex_data.formation.update()

			# Complex IDs in protein compartment file don't include modifications
			# Some have multiple alternative modifications so must loop through these
			for complex_data in me.process_data.query('^{:s}_mod_'.format(cplx)):
				# WARNING: FeFe and NiFe cofactors reform the formation reactions as follow:
				# requires a formation -> base_complex + FeFe/NiFe => base_complex_mod_FeFe/_mod_NiFe <- should not have a formation reaction
				# base_complex_mod_FeFe/_mod_NiFe + other cofactors => final modified complex
				lst = [ type(me.metabolites.get_by_id(x)) for x in complex_data.stoichiometry.keys() ]
				if coralme.core.component.Complex in lst:
					continue

				complex_data.stoichiometry.update(new_stoich[cplx])
				# remove zeroes from complex_data.stoichiometry
				complex_data.stoichiometry = { k:v for k,v in complex_data.stoichiometry.items() if v != 0 }
				complex_data.formation.update()

		# ## Part 6: Add Cell Wall Components
		# ### 1. Add lipid modification SubreactionData

		compartment_dict = {}
		for idx, compartment in df_protloc.set_index('Protein').Protein_compartment.to_dict().items():
			compartment_dict[idx] = compartment

		lipid_modifications = me.global_info.get('lipid_modifications')
		lipoprotein_precursors = me.global_info.get('lipoprotein_precursors', [])

		# Step1: assign enzymes to lipid modifications
		for data in me.process_data.query('^mod_1st'):
			data.enzyme = me.global_info.get('pg_pe_160', 'CPLX_dummy')
		for data in me.process_data.query('^mod_2nd'):
			data.enzyme = me.global_info.get('other_lipids', 'CPLX_dummy')

		# Step2: add reactions of lipoprotein formation
		if bool(config.get('add_lipoproteins', False)) and lipoprotein_precursors:
			coralme.builder.translocation.add_lipoprotein_formation(
				me, compartment_dict, lipoprotein_precursors, lipid_modifications, membrane_constraints = False, update = True)

		# ### 2. Correct complex formation IDs if they contain lipoproteins

		#for gene in tqdm.tqdm(coralme.builder.translocation.lipoprotein_precursors.values()):
		if bool(config.get('add_lipoproteins', False)) and lipoprotein_precursors:
			for gene in tqdm.tqdm(lipoprotein_precursors.values(), 'Adding lipid precursors and lipoproteins...', bar_format = bar_format):
				compartment = compartment_dict.get(gene)
				if compartment is None:
					pass
				else:
					for rxn in me.metabolites.get_by_id('protein_' + gene + '_' + compartment).reactions:
						if isinstance(rxn, coralme.core.reaction.ComplexFormation):
							data = me.process_data.get_by_id(rxn.complex_data_id)
							value = data.stoichiometry.pop('protein_' + gene + '_' + compartment)
							data.stoichiometry['protein_' + gene + '_lipoprotein' + '_' + compartment] = value
							rxn.update()

		# ### 3. Braun's lipoprotein demand
		# Metabolites and coefficients as defined in [Liu et al 2014](http://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-014-0110-6)
		brauns_lipid_mod = me.global_info['braun\'s_lipid_mod']

		if len(me.global_info['braun\'s_lipoproteins']) >= 1:
			for brauns_lipoprotein in me.global_info['braun\'s_lipoproteins']:
				# Perform checks before attempt to incorporate the Braun's lipoproteins
				if not me.metabolites.has_id(brauns_lipid_mod):
					logging.warning('The metabolite \'{:s}\' is not present in the ME-model and the Braun\'s lipoprotein demand cannot be set. Please check if it is the correct behavior. See http://bigg.ucsd.edu/universal/metabolites/murein5px4p for more information.'.format(brauns_lipid_mod))
					continue
				if not me.reactions.has_id('formation_{:s}'.format(brauns_lipoprotein)):
					logging.warning('The \'add_lipoproteins\' option is \'False\' or coralme failed to add the \'{:s}\' ComplexFormation reaction and the Braun\'s lipoprotein demand cannot be set. Please check if it is the correct behavior.'.format(brauns_lipoprotein))
					continue

				rxn = coralme.core.reaction.SummaryVariable('core_structural_demand_brauns_{:s}'.format(brauns_lipoprotein))
				murein5px4p = me.metabolites.get_by_id(brauns_lipid_mod)
				murein5px4p_mass = murein5px4p.formula_weight / 1000.
				# ecolime: 1.0 protein_b1677_lipoprotein_Outer_Membrane --> 1.0 EG10544-MONOMER (brauns_lipoprotein ID)
				lipoprotein = me.metabolites.get_by_id(brauns_lipoprotein)
				me.add_reactions([rxn])

				# biomass of lipoprotein accounted for in translation and lipid modification
				rxn.add_metabolites({
					murein5px4p : -abs(me.global_info['braun\'s_murein_flux']),
					lipoprotein : -abs(me.global_info['braun\'s_lpp_flux']),
					me.metabolites.peptidoglycan_biomass : abs(me.global_info['braun\'s_murein_flux']) * murein5px4p_mass
					},
					combine = False)
				rxn.lower_bound = me.mu # coralme.util.mu
				rxn.upper_bound = me.mu # coralme.util.mu
		else:
			logging.warning('No Braun\'s lipoprotein (lpp gene) homolog was set. Please check if it is the correct behavior.')

		# WARNING: Part 7 was originally "set keffs", however, formulas of complexes are corrected later and sasa can be underestimated
		# ## Part 7: Model updates and corrections
		# ### 1. Subsystems

		# Add reaction subsystems from M-model to ME-model
		for rxn in tqdm.tqdm(me.processed_m_model.reactions, 'Adding reaction subsystems from M-model into the ME-model...', bar_format = bar_format):
			if rxn.id in me.process_data:
				data = me.process_data.get_by_id(rxn.id)
			else:
				continue

			for parent_rxn in data.parent_reactions:
				parent_rxn.subsystem = rxn.subsystem

		# ### 2. Add enzymatic coupling for "carriers"
		# These are enzyme complexes that act as metabolites in a metabolic reaction.
		for data in tqdm.tqdm(list(me.stoichiometric_data), 'Processing StoichiometricData in ME-model...', bar_format = bar_format):
			if data.id == 'dummy_reaction':
				continue

			for met, value in data.stoichiometry.items():
				# This will add any complex to the list of enzymes
				if not isinstance(me.metabolites.get_by_id(met), coralme.core.component.Complex):
					continue

				subreaction_id = met + '_carrier_activity'
				if subreaction_id not in me.process_data:
					sub = coralme.core.processdata.SubreactionData(subreaction_id, me)
					sub.enzyme = [met]

				data.subreactions[subreaction_id] = abs(value)

		# ### 3. Update ME-model
		# trick to obtain shadow prices and reduced costs from the optimizer
		me.reactions.dummy_reaction_FWD_SPONT.objective_coefficient = 1.

		if update:
			me.update() # not parallelized

		# ### 4. Add remaining formulas and compartments to the ME-model
		for r in tqdm.tqdm(me.reactions.query('^formation_'), 'Updating all FormationReactions...', bar_format = bar_format):
			r.update()

		modification_formulas = df_mets[df_mets['type'].str.match('COFACTOR|MOD|MODIFICATION')]
		modification_formulas = dict(zip(modification_formulas['me_id'], modification_formulas['formula']))

		# Correct formula of complexes based on their base complex
		# This will add the formula to complexes not formed from a complex formation reaction (e.g. CPLX + na2_c -> CPLX_mod_na2(1))
		#coralme.builder.formulas.add_remaining_complex_formulas(me, modification_formulas)
		for met in [ x for x in me.metabolites if '_mod_' in x.id and isinstance(x, coralme.core.component.Complex)]:
			met.formula = None
			met.elements = {}

			base_complex = met.id.split('_mod_')[0]
			base_complex_elements = collections.Counter(me.metabolites.get_by_id(base_complex).elements)

			for mod in met.id.split('_mod_')[1:]:
				#for num in range(int(mod.rstrip(')').split('(')[1])):
				mod_elements = None
				mod_name = mod.split('(')[0]

				if mod_name in modification_formulas:
					mod_elements = coralme.builder.helper_functions.parse_composition(modification_formulas[mod_name])
					if me.metabolites.has_id(mod_name + '_c') and me.metabolites.get_by_id(mod_name + '_c').formula is None:
						me.metabolites.get_by_id(mod_name + '_c').formula = modification_formulas[mod_name]

				elif me.metabolites.has_id(mod_name + '_c') and me.metabolites.get_by_id(mod_name + '_c').formula is not None:
					mod_elements = me.metabolites.get_by_id(mod_name + '_c').elements

				# WARNING: flavodoxin homologs might have a different base_complex ID compared to ecolime
				# WARNING: Negative elemental contributions cannot be set in the metabolites.txt input file
				elif 'Oxidized(1)' == mod and 'FLAVODOXIN' not in base_complex:
					mod_elements = {'H': -2}
				elif 'Oxidized(2)' == mod and 'FLAVODOXIN' not in base_complex:
					mod_elements = {'H': -4}
				elif 'Oxidized(1)' == mod and 'FLAVODOXIN' in base_complex: # TODO: is the fmn cofactor in flavodoxin neutral?
					mod_elements = {'H': 0}

				elif 'glycyl(1)' == mod:
					mod_elements = {'H': -1}
				elif 'cosh(1)' == mod:
					mod_elements = {'H': +1, 'O': -1, 'S': +1}

				if mod_elements:
					mod_elements = collections.Counter(mod_elements)
					mod_elements = { k:v * int(re.findall(r'\((\d+)\)', mod)[0]) for k,v in mod_elements.items() }
					base_complex_elements.update(mod_elements)
				else:
					logging.warning('Attempt to calculate a corrected formula for \'{:s}\' failed. Please check if it is the correct behaviour, or if the modification \'{:s}_c\' exists as a metabolite in the ME-model or a formula is included in the me_mets.txt file.'.format(met.id, mod_name))
			complex_elements = { k:base_complex_elements[k] for k in sorted(base_complex_elements) if base_complex_elements[k] != 0 }
			met.formula = ''.join([ '{:s}{:d}'.format(k, v) for k,v in complex_elements.items() ])
			met.elements = coralme.builder.helper_functions.parse_composition(met.formula)
			logging.warning('Setting new formula for \'{:s}\' to \'{:s}\' successfully.'.format(met.id, met.formula))

		# Update a second time to incorporate all of the metabolite formulas correctly
		for data in tqdm.tqdm(me.subreaction_data.query(r'(?!^\w\w\w_addition_at_\w\w\w$)'), 'Recalculation of the elemental contribution in SubReactions...', bar_format = bar_format):
			data._element_contribution = data.calculate_element_contribution()

		# Update reactions affected by formula update
		for r in tqdm.tqdm(me.reactions.query('^formation_'), 'Updating all FormationReactions...', bar_format = bar_format):
			r.update()

		for r in tqdm.tqdm(me.reactions.query('_mod_lipoyl'), 'Updating FormationReactions involving a lipoyl prosthetic group...', bar_format = bar_format):
			r.update()

		for r in tqdm.tqdm(me.reactions.query('_mod_glycyl'), 'Updating FormationReactions involving a glycyl radical...', bar_format = bar_format):
			r.update()

		# Update biomass_constituent_demand reaction
		constituent_mass = sum(me.metabolites.get_by_id(c).formula_weight / 1000. * abs(v) for c,v in biomass_constituents.items())
		rxn = me.reactions.get_by_id('biomass_constituent_demand')
		rxn.add_metabolites({ k:-(abs(v)) for k,v in biomass_constituents.items() }, combine = False)
		rxn.add_metabolites({me.metabolites.get_by_id('constituent_biomass'): constituent_mass}, combine = False)

		# ## Part 8: Set keffs
		# Step 1. Determine SASA and median SASA
		if bool(config.get('estimate_keffs', True)):
			sasa_dct = {
				x.id:( x.formula_weight ** (3. / 4.) if x.formula_weight else 0, x.id if not x.formula_weight else False )
				for x in me.metabolites if type(x) == coralme.core.component.Complex
				}

			for met in [ v[1] for k,v in sasa_dct.items() ]:
				if met:
					logging.warning('The complex \'{:s}\' has no valid formula to determine its molecular weight.'.format(met))
					logging.warning('Please, set a value in the keff input file for reactions associated to the \'{:s}\' complex.'.format(met))

			median_sasa = numpy.median([ v[0] for k,v in sasa_dct.items() if v[0] != 0 ])

			me.global_info['median_sasa'] = median_sasa
			me.global_info['sasa_estimation'] = sasa_dct

			# Step 2: Estimate keff for all the reactions in the model
			#if "complex" not in df_keffs.columns: #df_keffs.empty: # The if True avoids the estimation if the user uses an "incomplete" input
			reaction_ids = [
				rxn for rxn in me.reactions if isinstance(rxn, coralme.core.reaction.MetabolicReaction)
				if rxn.id not in [ 'dummy_reaction_FWD_SPONT', 'dummy_reaction_REV_SPONT' ]
				if rxn._complex_data is not None
				]

			#if 'complex' in df_keffs.columns: # user provided a file with keffs # This is always True because w/o a user file, df_keffs is empty
			with open('{:s}/building_data/reaction_median_keffs.txt'.format(me.global_info['out_directory']), 'r') as infile:
				reaction_median_keffs = pandas.read_csv(infile, sep = '\t').set_index('reaction')

			mapped_keffs = {}
			for rxn in tqdm.tqdm(reaction_ids, 'Estimating effective turnover rates for reactions using the SASA method...', bar_format = bar_format):
				logging.warning('Estimating effective turnover rates for reaction \'{:s}\''.format(rxn.id))

				base_id = rxn._stoichiometric_data.id
				cplx_id = me.metabolites.get_by_id(rxn._complex_data.id).id

				if base_id not in reaction_median_keffs.index:
					continue

				median_keff = reaction_median_keffs.T[base_id]['keff']
				sasa = sasa_dct[cplx_id][0]
				keff = sasa * median_keff / median_sasa
				mapped_keffs[rxn] = 3000. if keff > 3000. else 0.01 if keff < 0.01 else keff

			# WARNING: keff won't be deprecated by coupling_coefficient_enzyme or coupling_coefficient_subreaction properties
			# WARNING: keff maps back to the value used to determine coupling_coefficient_enzyme or coupling_coefficient_subreaction
			# dictionary of reaction IDs : coralme.core.reaction objects
			rxns_to_map = { x.id:x for x in me.subreaction_data + me.reactions if hasattr(x, 'keff') }
			# Step 3: Replace user values if they match
			for idx, row in tqdm.tqdm(list(df_keffs.iterrows()), 'Mapping effective turnover rates from user input...', bar_format = bar_format):
				if row['direction'] == '' and row['complex'] == '' and row['mods'] == '':
					# subreactions have ID = reaction_name
					idx = row['reaction']
				else:
					# metabolic reactions have ID = reaction_name + direction + complex
					idx = '{:s}_{:s}_{:s}'.format(row['reaction'], row['direction'], row['complex'])
					if row['mods'] != '':
						idx = '{:s}_mod_{:s}'.format(idx, '_mod_'.join(row['mods'].split(' AND ')))

				if idx in rxns_to_map.keys():
					mapped_keffs[rxns_to_map[idx]] = 3000. if float(row['keff']) > 3000. else 0.01 if float(row['keff']) < 0.01 else row['keff']
					logging.warning('Mapping of the effective turnover rate for \'{:}\' with a user provided value.'.format(idx))
				else:
					logging.warning('Mapping of the effective turnover rate for \'{:}\' reaction failed. Please check if the reaction or subreaction is in the ME-model.'.format(idx))

			# Step 4: Set keffs
			if mapped_keffs:
				for rxn, keff in tqdm.tqdm(sorted(mapped_keffs.items(), key = lambda x: x[0].id), 'Setting the effective turnover rates using user input...', bar_format = bar_format):
					if isinstance(rxn, coralme.core.reaction.MetabolicReaction):
						rxn.coupling_coefficient_enzyme = float(keff)
					if isinstance(rxn, coralme.core.processdata.SubreactionData):
						rxn.coupling_coefficient_subreaction = float(keff)
					# update() will convert the coupling coefficient into mu / value / 3600.
					if hasattr(rxn, 'update'): # subreactions has no update attribute
						rxn.update()
					logging.warning('Setting the effective turnover rate for \'{:s}\' in {:f} successfully.'.format(rxn.id, float(keff)))

		# ### 5. Add metabolite compartments
		coralme.builder.compartments.add_compartments_to_model(me)

		# ### 6. Prune reactions from ME-model
		# WARNING: Do it recursively to reduce further the size of the ME-model.
		if prune:
			rnum = len(me.reactions)
			delta = 1
			while delta > 0:
				me.prune(skip=skip)
				delta = rnum - len(me.reactions)
				rnum = len(me.reactions)

		# Part 9. Save and report
		with open('{:s}/MEModel-step2-{:s}.pkl'.format(out_directory, model), 'wb') as outfile:
			pickle.dump(me, outfile)

		ListHandler.print_and_log('ME-model was saved in the {:s} directory as MEModel-step2-{:s}.pkl'.format(out_directory, model))

		n_mets = len(me.metabolites)
		new_mets = n_mets * 100. / len(me.gem.metabolites) - 100
		n_rxns = len(me.reactions)
		new_rxns = n_rxns * 100. / len(me.gem.reactions) - 100
		n_genes = len(me.metabolites.query(re.compile('^RNA_(?!biomass|dummy|degradosome)')))
		new_genes = n_genes * 100. / len(me.gem.genes) - 100

		ListHandler.print_and_log('ME-model reconstruction is done.')
		ListHandler.print_and_log('Number of metabolites in the ME-model is {:d} (+{:.2f}%, from {:d})'.format(n_mets, new_mets, len(me.gem.metabolites)))
		ListHandler.print_and_log('Number of reactions in the ME-model is {:d} (+{:.2f}%, from {:d})'.format(n_rxns, new_rxns, len(me.gem.reactions)))
		ListHandler.print_and_log('Number of genes in the ME-model is {:d} (+{:.2f}%, from {:d})'.format(n_genes, new_genes, len(me.gem.genes)))
		if hasattr(self, 'ref'):
			ListHandler.print_and_log('Number of missing genes from reconstruction with homology, but no function is {:d}. Check the curation notes for more details.'.format(coralme.builder.helper_functions.check_me_coverage(self)))
		else:
			ListHandler.print_and_log('Number of missing genes from reconstruction cannot be determined.')

		logging.shutdown()

		coralme.builder.notes.save_curation_notes(
				self.curation_notes,
				self.configuration['out_directory'] + '/curation_notes.json'
			)
		coralme.builder.notes.publish_curation_notes(
				self.curation_notes,
				self.configuration['out_directory']+ '/curation_notes.txt'
			)

		with open('{:s}/MEReconstruction-{:s}.log'.format(log_directory, model), 'w') as outfile:
			for filename in [
				'{:s}/MEReconstruction-step1-{:s}.log'.format(log_directory, model),
				'{:s}/MEReconstruction-step2-{:s}.log'.format(log_directory, model)
				]:

				try:
					pathlib.Path(filename).unlink(missing_ok = True) # python>=3.8
				except:
					if pathlib.Path(filename).exists():
						pathlib.Path(filename).unlink() # python==3.7

			logger = self.logger['MEReconstruction-step1'].log_list
			logger += self.logger['MEReconstruction-step2'].log_list

			tmp = pandas.DataFrame(logger)
			for idx, data in tmp.drop_duplicates(subset = 1).iterrows():
				outfile.write('{:s} {:s}\n'.format(data[0], data[1]))

class METroubleshooter(object):
	"""
	METroubleshooter class for troubleshooting growth in a ME-model

	This class contains methods to identify gaps and obtain a feasible ME-model.

	Parameters
	----------
	MEBuilder : coralme.builder.main.MEBuilder
	"""

	def __init__(self, builder):
		self.logger = builder.logger
		self.me_model = builder.me_model
		self.configuration = builder.configuration
		self.curation_notes = builder.curation_notes

	def troubleshoot(self, growth_key_and_value = None, skip = set(),
		guesses = [], met_types = [], platform = None, solver = 'qminos', savefile = None,
		gapfill_cofactors = False):
		"""Performs the Gap-finding step of the reconstruction.

		This function will iterate through different parts of the M-
		and E-matrices, looking for a minimal set of sinks that
		allows for growth.

		Parameters
		----------
		growth_key_and_value : dict
			A dictionary of Sympy.Symbol and value(s) to evaluate. It defines
			the parameters for the feasibility checks in each iteration.
		skip : set
			A set of ME-components to not evaluate
		met_types: list
			Any combination of 'ME-Deadends', 'Cofactors', 'All-Deadends',
			'Metabolite', 'GenerictRNA', 'Complex', 'TranscribedGene',
			'TranslatedGene', 'ProcessedProtein', and/or 'GenericComponent'.
		platform: str
			'win32' or 'darwin' to use gurobi (default) or cplex as solver
		solver: str
			Solver to use. Values: 'gurobi' (default) or 'cplex'
		"""
		types = {
			'M-matrix' : [ 'ME-Deadends', 'Cofactors', 'All-Deadends', 'Metabolite' ],
			'E-matrix' : [ 'GenerictRNA', 'Complex', 'TranscribedGene', 'TranslatedGene', 'ProcessedProtein', 'GenericComponent' ]
			}

		if hasattr(self, 'notes') and self.notes.get('from cobra', False):
			met_types = [ ('M-matrix', 'Metabolite') ]
		elif len(met_types) > 0:
			met_types = [ ('M-matrix', x) if x in types['M-matrix'] else ('E-matrix', x) if x in types['E-matrix'] else None for x in set(met_types) ]
			met_types = [ x for x in met_types if x is not None ]
		else:
			logging.warning('Metabolite types valid values are {:s}. The predefined order of metabolites will be tested.\n'.format(', '.join(types['M-matrix'] + types['E-matrix'])))
			met_types = []
			for x, y in types.items():
				for met in y:
					met_types.append((x, met))

		if not hasattr(self, 'me_model'):
			me = self
			self = coralme.builder.main.MEBuilder(**{'out_directory' : '.'})
			self.me_model = me
			self.configuration['out_directory'] = './'
			self.configuration['log_directory'] = './'

		self.me_model.troubleshooting = True
		if solver in ['gurobi', 'cplex']:
			self.me_model.get_solution = self.me_model.optimize_windows
			self.me_model.check_feasibility = self.me_model.feas_windows(solver = solver)
		elif solver == "qminos":
			self.me_model.get_solution = self.me_model.optimize
			self.me_model.check_feasibility = self.me_model.feasibility
			print('The MINOS and quad MINOS solvers are a courtesy of Prof Michael A. Saunders. Please cite Ma, D., Yang, L., Fleming, R. et al. Reliable and efficient solution of genome-scale models of Metabolism and macromolecular Expression. Sci Rep 7, 40863 (2017). https://doi.org/10.1038/srep40863\n')
		else:
			raise ValueError('Solver not supported.')

		config = self.configuration
		model = config.get('ME-Model-ID', 'coralME')
		out_directory = config.get('out_directory', '.')
		log_directory = config.get('log_directory', '.')

		# clean restart if troubleshooter is killed
		if not hasattr(self.me_model, "troubleshooted") or self.me_model.troubleshooted:
			pass # do not remove TS reactions
		else:
			rxns = self.me_model.reactions.query('^TS_')
			self.me_model.remove_reactions(rxns)

		# set logger
		log = logging.getLogger() # root logger
		for hdlr in log.handlers[:]: # remove all old handlers
			log.removeHandler(hdlr)

		# Old code works in a separate script; but it works if we remove the old handler
		logging.basicConfig(filename = '{:s}/METroubleshooter-{:s}.log'.format(log_directory, model), filemode = 'w', level = logging.WARNING, format = log_format)
		log.addHandler(self.logger['METroubleshooter'])
		log.addHandler(logging.StreamHandler(sys.stdout))
		logging.captureWarnings(True)

		if growth_key_and_value is None:
			growth_key_and_value = { self.me_model.mu : 0.001 }

		growth_key, growth_value = zip(*growth_key_and_value.items())

		logging.warning('~ '*1 + 'Troubleshooting started...')

		if gapfill_cofactors:
			# Ensure cofactors can be produced
			logging.warning('  '*1 + 'Ensuring the ME-model can produce all cofactors')
			cofactors = coralme.builder.troubleshooting.get_cofactors_in_me_model(self.me_model)
			ts_cofactors = coralme.builder.troubleshooting.add_exchange_reactions(self.me_model, cofactors, prefix = 'COFACTOR_TS_')
			for ts in ts_cofactors:
				ts.bounds = (1e-6,1000)

		# Step 1. Test if current ME-model is feasible
		logging.warning('  '*1 + 'Checking if the ME-model can simulate growth without gapfilling reactions...')
		if self.me_model.check_feasibility(keys = growth_key_and_value):
			logging.warning('  '*1 + 'Original ME-model is feasible with a tested growth rate of {:f} 1/h'.format(list(growth_value)[0]))
			works = True
		else:
			logging.warning('  '*1 + 'Original ME-model is not feasible with a tested growth rate of {:f} 1/h'.format(list(growth_value)[0]))
			works = False

		# Step 2. Test different sets of MEComponents
		if len(guesses) > 0:
			guesses = [ x for x in guesses if self.me_model.metabolites.has_id(x) ]
			if len(guesses) > 0:
				met_types.insert(0, (guesses, 'User guesses'))

		e_gaps = []
		history = dict()
		if works == False:
			#logging.warning('~ '*1 + 'Step 3. Attempt gapfilling different groups of E-matrix components.')
			for idx, met_type in enumerate(met_types):
				logging.warning('  '*1 + 'Step {}. Gapfill reactions to provide components of type \'{:s}\' using brute force.'.format(idx + 1, met_type[1]))
				if met_type[0] == 'E-matrix':
					logging.warning('  '*5 + 'Relaxing bounds for E-matrix gap-fill')
					self.me_model.relax_bounds()
					self.me_model.reactions.protein_biomass_to_biomass.lower_bound = growth_value[0]/100 # Needed to enforce protein production
				if met_type[1] == 'User guesses':
					history, output = coralme.builder.troubleshooting.brute_check(self.me_model, growth_key_and_value, met_type, skip = skip, history = history)
				else:
					history, output = coralme.builder.troubleshooting.brute_check(self.me_model, growth_key_and_value, met_type[1], skip = skip, history = history)
				bf_gaps, no_gaps, works = output
				# close sink reactions that are not gaps
				if no_gaps:
					self.me_model.remove_reactions(no_gaps)
				if works:
					e_gaps = bf_gaps
					break

		if works: # Meaning it can grow in any case
			# Save warnings
			if isinstance(e_gaps, list) and e_gaps:
				self.curation_notes['troubleshoot'].append({
					'msg':'Some metabolites are necessary for growth',
					'triggered_by':e_gaps,#[ x for y in e_gaps for x in y ],
					'importance':'critical',
					'to_do':'Fix the gaps by adding reactions or solving other warnings. If some items are from the E-matrix, fix these first!'})

			# delete added sink reactions with lb == 0 and ub == 0
			sinks = []
			for rxn in self.me_model.reactions.query('^TS_'):
				sinks.append(rxn.id)
				#f = self.me_model.solution.fluxes[rxn.id]
				if rxn.lower_bound == 0 and rxn.upper_bound == 0:# or f == 0:
					self.me_model.remove_reactions([rxn])
			if sinks:
				logging.warning('~ '*1 + 'Troubleshooter added the following sinks: {:s}.'.format(', '.join(sinks)))
			logging.warning('~ '*1 + 'Final step. Fully optimizing with precision 1e-6 and save solution into the ME-model...')

			# Delete demand reactions for cofactor gapfilling
			if gapfill_cofactors:
				rxns = self.me_model.reactions.query('^COFACTOR_TS_')
				self.me_model.remove_reactions(rxns)

			# final optimization
			if self.me_model.get_solution(max_mu = 3.0, precision = 1e-6, verbose = False):
				if hasattr(self.me_model, 'gem'):
					logging.warning('  '*1 + 'Gapfilled ME-model is feasible with growth rate {:f} (M-model: {:f}).'.format(self.me_model.solution.objective_value, self.me_model.gem.optimize().objective_value))
				else:
					logging.warning('  '*1 + 'Gapfilled ME-model is feasible with growth rate {:f}.'.format(self.me_model.solution.objective_value))
			else:
				logging.warning('  '*1 + 'Error: Gapfilled ME-model is not feasible ?')

			# save model as a pickle file
			if savefile is None:
				savefile = '{:s}/MEModel-step3-{:s}-TS.pkl'.format(out_directory, self.me_model.id)
				message = 'ME-model was saved in the {:s} directory as MEModel-step3-{:s}-TS.pkl'.format(out_directory, self.me_model.id)
			elif pathlib.Path(savefile).parent.exists():
				message = 'ME-model was saved to {:s}.'.format(savefile)
			else:
				message = False
				logging.warning('Model was not saved. Please do it manually.')

			self.me_model.troubleshooted = True
			if message:
				with open(savefile, 'wb') as outfile:
					pickle.dump(self.me_model, outfile)
				logging.warning(message)

		else:
			logging.warning('~ '*1 + 'METroubleshooter failed to determine a set of problematic metabolites.')
			self.me_model.troubleshooted = False

		logging.shutdown()

		coralme.builder.notes.save_curation_notes(
				self.curation_notes,
				self.configuration['out_directory'] + '/curation_notes.json'
			)
		coralme.builder.notes.publish_curation_notes(
				self.curation_notes,
				self.configuration['out_directory']+ '/curation_notes.txt'
			)

		# We will remove duplicates entries in the log output
		with open('{:s}/METroubleshooter-{:s}.log'.format(log_directory, model), 'w') as outfile:
			logger = self.logger['METroubleshooter'].log_list

			tmp = pandas.DataFrame(logger)
			for idx, data in tmp.drop_duplicates(subset = 1).iterrows():
				outfile.write('{:s} {:s}\n'.format(data[0], data[1]))

		self.me_model.troubleshooting = False
