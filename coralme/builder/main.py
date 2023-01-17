#!/usr/bin/python3
# python imports
import io
import os
import re
import pickle
import pathlib
import logging
import warnings
import collections
import subprocess

# third party imports
import tqdm
import sympy
import pandas
import anyconfig

import cobra
import coralme

# configuration
log_format = '%(asctime)s %(message)s' #%(clientip)-15s %(user)-8s
bar_format = '{desc:<75}: {percentage:.1f}%|{bar}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
try:
	warnings.simplefilter(action = 'ignore', category = pandas.errors.SettingWithCopyWarning)
except:
	warnings.warn("This pandas version does not allow for correct warning handling. Pandas >=1.5.1 is suggested.")

class MEBuilder(object):
	"""
	MEBuilder class to construct input data from protein homology.

	Parameters
	----------

	"""
	def __init__(self, *args, **kwargs):
		config = {}
		for input_file in args:
			with open(input_file, 'r') as infile:
				config.update(anyconfig.load(infile))

		if kwargs:
			config.update(kwargs)

		self.me_model = coralme.core.model.MEModel(config.get('model_id', 'coralME'), config.get('growth_key', 'mu'))
		self.configuration = config
		self.curation_notes = { 'troubleshoot' : [] }

		data = \
			'code,interpretation,gram\n' \
			'CCI-CW-BAC-POS-GP,Cell_Wall,pos\n' \
			'CCI-OUTER-MEM-GN,Outer_Membrane,neg\n' \
			'CCI-PM-BAC-NEG-GN,Inner_Membrane,neg\n' \
			'CCI-PM-BAC-POS-GP,Plasma_Membrane,pos\n' \
			'CCO-MEMBRANE,Membrane,'

		self.location_interpreter = pandas.read_csv(io.StringIO(data))

		return None

	def generate_files(self):
		sep = ""
		config = self.configuration

		print("{}Reading organism{}".format(sep, sep))

		# Read organism
		self.org = coralme.builder.organism.Organism(config, is_reference = False)

		# self.org.rpod = ''
		# self.org.get_rna_polymerase(force_RNAP_as='')

		print("{} Modifying and preparing M-model {}".format(sep, sep))
		self.modify_metabolic_reactions('reaction_corrections.csv')
		self.prepare_model()

		# ## Homology with reference
		# Curation note: Check which locus tag field in the genbank agrees with those in the biocyc files.
		# Make sure you have specified that field in the beginning of this notebook.
		# Reference
		if bool(config.get('dev_reference', False)) or bool(config.get('user_reference', False)):
			print("{} Reading reference {}".format(sep, sep))
			#with open(org_files_dir + ref + '/parameters.txt') as rf:
				#ref_parameters = rf.read()
			#ref_parameters = ast.literal_eval(ref_parameters)
			#self.read_reference(ref, locus_tag=ref_parameters['locus_tag'])

			self.ref = coralme.builder.organism.Organism(config, is_reference = True)

			folder = self.org.blast_directory
			if bool(config.get('run_bbh_blast', False)):
				print("{} Running BLAST {}".format(sep, sep))
				self.org.gb_to_faa('org', element_types = {'CDS'}, outdir = self.org.blast_directory)
				self.ref.gb_to_faa('ref', element_types = {'CDS'}, outdir = self.org.blast_directory)

				def execute(cmd):
					cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
					out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

				# make blast databases
				execute('makeblastdb -in {:s}/org.faa -dbtype prot -out {:s}/org'.format(folder, folder))
				execute('makeblastdb -in {:s}/ref.faa -dbtype prot -out {:s}/ref'.format(folder, folder))

				# bidirectional blast
				execute('blastp -db {:s}/org -query {:s}/ref.faa -num_threads 4 -out {:s}/org_as_db.txt -outfmt 6'.format(folder, folder, folder))
				execute('blastp -db {:s}/ref -query {:s}/org.faa -num_threads 4 -out {:s}/ref_as_db.txt -outfmt 6'.format(folder, folder, folder))

				#os.system('{}/auto_blast.sh {}'.format(self.directory,self.org.directory))

			# #### Reciprocal hits
			print("{} Getting homologs {}".format(sep, sep))

			self.get_homology(evalue = 1e-10)
			self.homology.mutual_hits_df.to_csv('{:s}/mutual_hits.csv'.format(folder))
			#self.homology.mutual_hits_df.to_csv(self.org.directory+'{:s}/mutual_hits.csv'.format(folder))

			# #### Get enzyme homology
			self.homology.get_complex_homology()

			# #### Update model info with homology
			print("{} Updating from homology {}".format(sep, sep))
			self.update_from_homology()

		# #### Manual curation
		print("{} Integrating manual curation of complexes {}".format(sep, sep))
		self.add_manual_complexes()

		# ## enzyme_reaction_association.txt
		print("{} Getting enzyme-reaction association {}".format(sep, sep))
		self.get_enzyme_reaction_association()

		# ## Keffs
		print("{} Setting reaction Keffs {}".format(sep, sep))
		self.org.get_reaction_keffs()

		# ## Files from model
		# #### Metabolites
		print("{} Generating metabolites file {}".format(sep, sep))
		self.org.generate_metabolites_file()

		# #### Reactions
		print("{} Generating reactions file {}".format(sep, sep))
		self.org.generate_reactions_file()

		# #### Reaction matrix
		print("{} Generating reaction matrix file {}".format(sep, sep))
		self.org.generate_reaction_matrix()

		# #### Biomass constituents
		self.org.biomass_constituents = config.get('flux_of_biomass_constituents', {})

		# Fill builder with dummy
		#print("{} Filling files with CPLX_dummy {}".format(sep, sep))
		#self.fill()
		# Final checks of builder
		print("{} Performing final checks of files {}".format(sep, sep))
		self.check()
		# Save builder
		#print("{} Saving builder {}".format(sep, sep))
		#self.save()
		# Save builer info in human readable form
		#print("{} Saving builder info {}".format(sep, sep))
		#self.save_builder_info()
		# Update notes
		print("{} Generating curation notes {}".format(sep, sep))
		self.org.generate_curation_notes()

	def prepare_model(self):
		m_model = self.org.m_model
		target_compartments = {"c": "Cytosol", "e": "Extra-organism", "p": "Periplasm"}
		new_dict = {}
		warn_compartments = []
		for k,v in m_model.compartments.items():
			if k in target_compartments:
				new_dict[k] = target_compartments[k]
			else:
				warn_compartments.append(k)
				new_dict[k] = k
		m_model.compartments = new_dict
		gene_list = []
		for m in m_model.metabolites:
			if not m.compartment:
				print("Fixing compartment for metabolite {}".format(m.id))
				m.compartment = m.id[-1]
		for r in m_model.reactions:
			if not r.name or isinstance(r.name, float):
				print("Fixing name for reaction {}".format(r.id))
				r.name = r.id

		# Solve m_model
		solution = m_model.optimize()
		if not solution or solution.objective_value < 1e-6:
			self.org.curation_notes['prepare_model'].append({
				'msg':'M-model cannot grow',
				'importance':'critical',
				'to_do':'Check that the model you provided works'})

		self.org.biomass = str(m_model.objective.expression.as_two_terms()[0]).split('*')[1]
		biomass_rxn = m_model.reactions.get_by_id(self.org.biomass)
		print('{} was identified as the biomass reaction'.format(biomass_rxn.id))


		adp = m_model.metabolites.adp_c
		# Get GAM
		self.org.GAM = None
		if adp in biomass_rxn.metabolites:
			self.org.GAM = biomass_rxn.metabolites[adp]
			print('GAM identified with value {}'.format(self.org.GAM))
		else:
			self.org.GAM = 45.
			self.org.curation_notes['prepare_model'].append({
				'msg':'GAM could not be identified from biomass reaction, setting a standard value of 45. adp_c is not present as a product.',
				'importance':'high',
				'to_do':'Check whether the biomass reaction was read or defined correctly. You can define GAM with me_builder.org.GAM = GAM_value'})
		# Get NGAM
		NGAMs = ['NGAM','ATPM']
		self.org.NGAM = None
		for r in NGAMs:
			if r in m_model.reactions:
				rxn = m_model.reactions.get_by_id(r)
				self.org.NGAM = rxn.lower_bound
				print('{} was identified as NGAM with value {}'.format(r,self.org.NGAM))
				break
		if self.org.NGAM == None:
			self.org.NGAM = 1.
			self.org.curation_notes['prepare_model'].append({
				'msg':'NGAM could not be identified in M-model, setting a standard value of 1.',
				'importance':'high',
				'to_do':'Manually define NGAM with me_builder.org.NGAM = NGAM_value'})
		elif self.org.NGAM == 0:
			self.org.NGAM = 1.
			self.org.curation_notes['prepare_model'].append({
				'msg':'NGAM was identified from reaction {}, but its lower bound is 0. NGAM set to a standard value of 1.0.'.format(rxn.id),
				'importance':'high',
				'to_do':'Manually define NGAM with me_builder.org.NGAM = NGAM_value'})
		# Warnings
		if warn_compartments:
			self.org.curation_notes['prepare_model'].append({
				'msg':'Some compartments in m_model are unknown',
				'triggered_by':warn_compartments,
				'importance':'medium',
				'to_do':'Check whether the compartment is correct. If not, change it in the reaction ID in the m_model.'})

	def modify_metabolic_reactions(self,
								   filename):
		m_model = self.org.m_model
		filename = self.org.directory + filename
		try:
			new_reactions_dict = (
				pandas.read_csv(filename, index_col=0)
				.fillna({"gene_reaction_rule": "", "notes": "", "reaction": ""})
				.T.to_dict()
			)
		except:
			new_reactions_dict = {}
			self.org.curation_notes['modify_metabolic_reactions'].append({
				'msg':'No {} file, creating one.'.format(filename),
				'importance':'low',
				'to_do':'Fill {}'.format(filename)})
			pandas.DataFrame.from_dict(
				{
					"reaction_id": {},
					"name": {},
					"gene_reaction_rule": {},
					"reaction": {},
					"notes": {},
				}
			).set_index("reaction_id").to_csv(filename)
		for rxn_id, info in new_reactions_dict.items():
			if info["reaction"] == "eliminate":
				m_model.reactions.get_by_id(rxn_id).remove_from_model()
			else:
				if rxn_id not in m_model.reactions:
					rxn = cobra.Reaction(rxn_id)
					m_model.add_reaction(rxn)
				else:
					rxn = m_model.reactions.get_by_id(rxn_id)
				if info["reaction"]:
					rxn.build_reaction_from_string(info["reaction"])
					rxn.name = info["name"]
					rxn.gene_reaction_rule = info["gene_reaction_rule"]
				if info["gene_reaction_rule"]:
					if info["gene_reaction_rule"] == "no_gene":
						rxn.gene_reaction_rule = ""
						rxn.name = str(rxn.name) + ', spontaneous'
					else:
						rxn.gene_reaction_rule = info["gene_reaction_rule"]
				if info["name"]:
					rxn.name = info["name"]

	def get_homology(self, evalue=1e-10):
		self.homology = coralme.builder.homology.Homology(self.org, self.ref, evalue = evalue)

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
				rg_id = re.findall('.*(?=\(\d*\))', rg)[0]
				coeff = re.findall("(?<=\()[0-9]{1,3}", rg)
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
					for idx, (_, row) in enumerate(ref_mods.iterrows()):
						mods = row["Modifications"].split(" AND ")
						cplx = c
						cofs = []
						coeffs = []
						for mod in mods:
							cof = re.findall(".*(?=\()", mod)[0]
							coeff = re.findall("(?<=\()[0-9]{1}", mod)
							if coeff:
								coeff = coeff[0]
								cplx += "_mod_{}({})".format(cof, coeff)
							else:
								cplx += "_mod_{}".format(cof)
								coeff = ""
							cofs.append(cof)
							coeffs.append(coeff)
						protein_mod_dict[cplx] = {}
						protein_mod_dict[cplx]["Core_enzyme"] = c
						protein_mod_dict[cplx]["Modifications"] = " AND ".join(
							"{}({})".format(cof, coeff)
							for cof, coeff in zip(cofs, coeffs)
						)
						protein_mod_dict[cplx]["Source"] = "E_coli_homology"
						ref_cplx_homolog[row["Modified_enzyme"]] = cplx
						cplx_homolog[cplx] = row["Modified_enzyme"]
		protein_mod = pandas.DataFrame.from_dict(protein_mod_dict).T
		protein_mod.index.name = "Modified_enzyme"
		self.org.protein_mod = pandas.concat([self.org.protein_mod,protein_mod])

	def add_manual_complexes(self):
		manual_complexes = self.org.manual_complexes
		complexes_df = self.org.complexes_df
		protein_mod = self.org.protein_mod
		warn_manual_mod = []
		warn_replace = []
		for new_complex, info in manual_complexes.iterrows():
			if info["genes"]:
				if new_complex not in complexes_df:
					complexes_df = complexes_df.append(
						pandas.DataFrame.from_dict(
							{new_complex: {"name": "", "genes": "", "source": "Manual"}}
						).T
					)
				complexes_df.loc[new_complex, "genes"] = info["genes"]
				complexes_df.loc[new_complex, "name"] = str(info["name"])
			if info["mod"]:
				mod_complex = (
					new_complex
					+ "".join(
						[
							"_mod_{}".format(m)
							for m in info['mod'].split(' AND ')
						]
					)
					if info["mod"]
					else new_complex
				)
				if mod_complex in protein_mod.index:
					warn_manual_mod.append(mod_complex)
					continue
				if info["replace"]:
					if info["replace"] in protein_mod.index:
						protein_mod = protein_mod.drop(info["replace"])
					else:
						warn_replace.append(mod_complex)
				protein_mod = protein_mod.append(
					pandas.DataFrame.from_dict(
						{
							mod_complex: {
								"Core_enzyme": new_complex,
								"Modifications": " AND ".join(
									"{}({})".format(cof, coeff)
									for cof, coeff in ast.literal_eval(
										info["mod"]
									).items()
								),
								"Source": "Manual",
							}
						}
					).T
				)
		complexes_df.index.name = "complex"

		self.org.complexes_df = complexes_df
		self.org.protein_mod = protein_mod

		# Warnings
		if warn_manual_mod or warn_replace:
			if warn_manual_mod:
				self.org.curation_notes['add_manual_complexes'].append({
					'msg':'Some modifications in protein_corrections.csv are already in me_builder.org.protein_mod and were skipped.',
					'triggered_by':warn_manual_mod,
					'importance':'low',
					'to_do':'Check whether the protein modification specified in protein_corrections.csv is correct and not duplicated.'})
			if warn_replace:
				self.org.curation_notes['add_manual_complexes'].append({
					'msg':'Some modified proteins marked for replacement in protein_corrections.csv are not in me_builder.org.protein_mod. Did nothing.',
					'triggered_by':warn_replace,
					'importance':'low',
					'to_do':'Check whether the marked modified protein in protein_corrections.csv for replacement is correctly defined.'})

	def get_enzyme_reaction_association(self, gpr_combination_cutoff = 100):
		import numpy as np
		#from draft_cobrame.util.helper_functions import process_rule_dict, find_match

		m_model = self.org.m_model
		org_complexes_df = self.org.complexes_df
		protein_mod = self.org.protein_mod
		gene_dictionary = (
			self.org.gene_dictionary.reset_index()
			.set_index("Accession-1")
		)
		generic_dict = self.org.generic_dict
		enz_rxn_assoc_dict = {}

		for rxn in m_model.reactions:
			unnamed_counter = 0
			rule = str(rxn.gene_reaction_rule)
			if not rule:
				continue
			enz_rxn_assoc_dict[rxn.id] = []
			#rule_list = expand_gpr(listify_gpr(rule)).split(" or ")
			rule_list = coralme.builder.helper_functions.expand_gpr(rule)
			if len(rule_list) <= gpr_combination_cutoff:
				enz_rxn_assoc = []
				reaction_cplx_list = []
				for rule_gene_list in rule_list:
					identified_genes = []
					for i in rule_gene_list:
						identified_genes.append(i)
					cplx_id = coralme.builder.helper_functions.find_match(org_complexes_df["genes"].to_dict(),identified_genes)
					if not cplx_id:
						if len(identified_genes) > 1:
							# New cplx not found in BioCyc files
							cplx_id = "CPLX_{}-{}".format(rxn.id,unnamed_counter)
							unnamed_counter += 1
						else:
							gene = identified_genes[0]
							cplx_id = "{}-MONOMER".format(gene_dictionary.loc[gene]['Gene Name'])
						if cplx_id not in org_complexes_df.index:
							print("Adding {} to complexes from m_model".format(cplx_id))
							#org_complexes_df = org_complexes_df.append(
								#pandas.DataFrame.from_dict(
									#{
										#cplx_id: {
											#"name": str(rxn.name),
											#"genes": " AND ".join(
												#["{}()".format(g) for g in identified_genes]
											#),
											#"source": "{}({})".format(m_model.id, rxn.id),
										#}
									#}
								#).T
							#)
							tmp = pandas.DataFrame.from_dict({
								cplx_id: {
									"name": str(rxn.name),
									"genes": " AND ".join(["{}()".format(g) for g in identified_genes]),
									"source": "{}({})".format(m_model.id, rxn.id),
									}}).T
							org_complexes_df = pandas.concat([org_complexes_df, tmp], axis = 0, join = 'outer')
					if cplx_id in protein_mod["Core_enzyme"].values:
						cplx_id = protein_mod[
							protein_mod["Core_enzyme"].str.contains(cplx_id)
						].index[0]
						if "Oxidized" in cplx_id:
							cplx_id = cplx_id.split("_mod_Oxidized")[0]
					reaction_cplx_list.append(cplx_id)
				enz_rxn_assoc_dict[rxn.id] = " OR ".join(reaction_cplx_list)
			else:
				print('{} contains a GPR rule that has {} possible gene combinations. Generifying it.'.format(rxn.id,len(rule_list)))
				listified_gpr = coralme.builder.helper_functions.listify_gpr(rule)
				n,rule_dict = coralme.builder.helper_functions.generify_gpr(listified_gpr,rxn.id,d={})
				if not rule_dict: # n in gene_dictionary.index:
					product = gene_dictionary.loc[n,'Product']
					rule_dict[product] = n
					n = product
				n,rule_dict = coralme.builder.helper_functions.process_rule_dict(n,rule_dict,org_complexes_df["genes"].to_dict(),protein_mod)
				generified_rule = n
				#print(rxn.id)
				#print(rxn.gene_reaction_rule)
				#print(rule_dict)
				#print(generified_rule)
				#print()
				for cplx,rule in rule_dict.items():
					if 'mod' in cplx:
						cplx_id = cplx.split('_mod_')[0]
					else:
						cplx_id = cplx
					if 'generic' in cplx_id:
						print("Adding {} to generics from m_model".format(cplx_id))
						generic_dict[cplx_id] = {
							'enzymes':[gene_dictionary.loc[i,'Product'] if i in gene_dictionary.index else i for i in rule.split(' or ')]
						}
					elif cplx_id not in org_complexes_df.index:
						# New cplx not found in BioCyc files
						print("Adding {} to complexes from m_model".format(cplx_id))
						#org_complexes_df = org_complexes_df.append(
							#pandas.DataFrame.from_dict(
								#{
									#cplx_id: {
										#"name": str(rxn.name),
										#"genes": " AND ".join(
											#["{}()".format(g) for g in rule.split(' and ')]
										#),
										#"source": "{}({})".format(m_model.id, rxn.id),
									#}
								#}
							#).T
						#)
						tmp = pandas.DataFrame.from_dict({
							cplx_id: {
								"name": str(rxn.name),
								"genes": " AND ".join(["{}()".format(g) for g in rule.split(' and ')]),
								"source": "{}({})".format(m_model.id, rxn.id),
								}}).T
						org_complexes_df = pandas.concat([org_complexes_df, tmp], axis = 0, join = 'outer')
				enz_rxn_assoc_dict[rxn.id] = generified_rule
			enz_rxn_assoc_df = pandas.DataFrame.from_dict({"Complexes": enz_rxn_assoc_dict})
			enz_rxn_assoc_df = enz_rxn_assoc_df.replace(
				"", np.nan
			).dropna()  # Remove empty rules

		enz_rxn_assoc_df.index.name = "Reaction"
		self.org.enz_rxn_assoc_df = enz_rxn_assoc_df
		self.org.complexes_df = org_complexes_df
		self.org.protein_mod = protein_mod

	def update_TU_df(self):
		org_TU_to_genes = self.org.TU_to_genes
		org_TUs = self.org.TUs
		org_sigmas = self.org.sigmas
		org_complexes_df = self.org.complexes_df
		ref_TUs = self.ref.TUs
		ref_TU_df = self.ref.TU_df
		gene_dictionary = self.org.gene_dictionary
		mutual_hits = self.homology.mutual_hits
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		rpod = self.org.rpod
		ref_genes_to_TU = self.ref.genes_to_TU
		ref_sigmas = self.ref.sigmas
		TU_df = self.org.TU_df
		remove_TUs = []
		TU_dict = {}
		for tu_id, row in TU_df.iterrows():
			tu = tu_id.split("_from_")[0]
			rho_dependent = True
			sigma = rpod
			genes = org_TU_to_genes[tu]
			if set(genes).issubset(mutual_hits):
				ref_TU = [
					ref_genes_to_TU[mutual_hits[g]]
					for g in genes
					if mutual_hits[g] in ref_genes_to_TU
				]
				if (
					len(ref_TU) == 1
				):  # All mapped genes are from only one TU. TU identified!
					TU_hit = ref_TU_df[ref_TU_df.index.str.contains(ref_TU[0])]
					if not TU_hit.empty:
						rho_dependent = TU_hit["rho_dependent"].tolist()[0]
						ref_sigma = TU_hit["sigma"].tolist()[0]
						if ref_sigma in ref_cplx_homolog:
							sigma = ref_cplx_homolog[ref_sigma]
							if sigma not in org_sigmas.index:
								org_sigmas = org_sigmas.append(
									pandas.DataFrame.from_dict(
										{
											sigma: {
												"complex": "RNAP_{}".format(sigma),
												"genes": org_complexes_df.loc[sigma]["genes"],
												"name": org_complexes_df.loc[sigma]["name"],
											}
										}
									).T
								)
				tu_name = "{}_from_{}".format(tu, sigma)
				if tu_name not in TU_df.index:
					remove_TUs.append(tu_id)
					TU_df.loc[tu_name] = [0, 0, 0, 0, 0, 0]
					TU_df.loc[tu_name]["strand"] = row["strand"]
					TU_df.loc[tu_name]["start"] = int(row["start"])
					TU_df.loc[tu_name]["stop"] = int(row["stop"])
					TU_df.loc[tu_name]["tss"] = None
				TU_df.loc[tu_name]["rho_dependent"] = rho_dependent
				TU_df.loc[tu_name]["sigma"] = sigma
		self.org.TU_df = TU_df
		org_sigmas.index.name = "sigma"
		self.org.sigmas = org_sigmas

	def protein_location_from_homology(self):
		protein_location = self.org.protein_location
		complexes_df = self.org.complexes_df
		org_proteins_df = self.org.proteins_df
		org_cplx_homolog = self.homology.org_cplx_homolog
		mutual_hits = self.homology.mutual_hits
		ref_protein_location = self.ref.protein_location
		if not isinstance(ref_protein_location, pandas.DataFrame):
			return
		for rc,row in ref_protein_location.iterrows():
			ref_gene = re.findall('.*(?=\(.*\))', row['Protein'])[0]
			if ref_gene not in mutual_hits:
								continue
			org_gene = mutual_hits[ref_gene]
			ref_info = ref_protein_location[ref_protein_location['Protein'].str.contains(ref_gene)]
			gene_string = '{}\('.format(org_gene)
			org_cplxs = org_complexes_df[org_complexes_df['genes'].str.contains(gene_string)].index
			for org_cplx in org_cplxs:
				if protein_location.any().any() and \
					org_cplx in protein_location.index and \
					protein_location.loc[[org_cplx]]['Protein'].str.contains(gene_string).any().any():
					pass
					# Check if already in protein location, if not add.
				else:
					tmp = pandas.DataFrame.from_dict(
						{ org_cplx: {
							"Complex_compartment": ref_info["Complex_compartment"].values[0],
							"Protein": '{}()'.format(org_gene),
							"Protein_compartment": ref_info["Protein_compartment"].values[0],
							"translocase_pathway": ref_info["translocase_pathway"].values[0],
										}
									}).T
					protein_location = pandas.concat([protein_location, tmp], axis = 0, join = 'outer')
		self.org.protein_location = protein_location

	def update_translocation_multipliers(self):
		ref_multipliers = self.ref.translocation_multipliers
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		mutual_hits = self.homology.mutual_hits
		multipliers = {}
		for ref_cplx, genes in ref_multipliers.items():
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
		self.org.translocation_multipliers = pandas.DataFrame.from_dict(multipliers).fillna(
			0.0
		)

	def update_lipoprotein_precursors(self):
		ref_lipoprotein_precursors = self.ref.lipoprotein_precursors
		mutual_hits = self.homology.mutual_hits
		lipoprotein_precursors = {}
		for c, g in ref_lipoprotein_precursors.items():
			if g in mutual_hits:
				og = mutual_hits[g]
				lipoprotein_precursors[c] = og
		self.org.lipoprotein_precursors = lipoprotein_precursors

	def update_cleaved_methionine(self):
		cleaved_methionine = self.org.cleaved_methionine
		ref_cleaved_methionine = self.ref.cleaved_methionine
		mutual_hits = self.homology.mutual_hits
		for g in ref_cleaved_methionine:
			if g not in mutual_hits:
				continue
			cleaved_methionine.append(mutual_hits[g])
		self.org.cleaved_methionine = cleaved_methionine

	def update_m_to_me_mets(self):
		ref_m_to_me_mets = self.ref.m_to_me_mets
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		protein_mod = self.org.protein_mod.reset_index().set_index("Core_enzyme")
		m_to_me_mets = self.org.m_to_me_mets
		m_model = self.org.m_model
		d = {}
		warn_skip = []
		warn_found = []
		warn_skip_2 = []
		for ref_m, row in ref_m_to_me_mets.iterrows():
			if ref_m not in m_model.metabolites:
				if "__" in ref_m:
					ref_m = ref_m.replace("__", "_")
					if ref_m not in m_model.metabolites:
						warn_skip.append(ref_m)
						continue
					else:
						warn_found.append(ref_m)
			if ref_m in m_to_me_mets.index:
				continue
			ref_me = row["me_name"]
			d[ref_m] = {}
			if ref_me in ref_cplx_homolog:
				org_me = ref_cplx_homolog[ref_me]
				me_name = org_me
			elif ref_me == "eliminate":
				me_name = "eliminate"
			else:
				me_name = "curate"
			d[ref_m]["me_name"] = me_name
		if d:
			df = pandas.DataFrame.from_dict(d).T
			df.index.name = "m_name"
			df = df.reset_index()
			m_to_me_mets = m_to_me_mets.reset_index(drop=True)
			#m_to_me_mets = m_to_me_mets.append(df)
			m_to_me_mets = pandas.concat([m_to_me_mets, df], axis = 0, join = 'outer')
			m_to_me_mets = m_to_me_mets.set_index("m_name")
			self.org.m_to_me_mets = m_to_me_mets
		# Warnings
		if warn_skip or warn_found or warn_skip_2:
			if warn_skip:
				self.org.curation_notes['update_m_to_me_mets'].append({
					'msg':'Some metabolites in m_to_me_mets.csv are not in m_model, so they were skipped.',
					'triggered_by':warn_skip,
					'importance':'medium',
					'to_do':'Confirm these metabolites are correctly defined in m_to_me_mets.csv'})
			if warn_found:
				self.org.curation_notes['update_m_to_me_mets'].append({
					'msg':'Some metabolites in m_to_me_mets.csv were found in reference m_model after replacing __ with _',
					'triggered_by':warn_found,
					'importance':'medium',
					'to_do':'Confirm these metabolites are correctly defined in m_to_me_mets.csv'})
	def update_generics_from_homology(self):
		generic_dict = self.org.generic_dict
		ref_generic_dict = self.ref.generic_dict
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_generic_dict.items():
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
		for k, v in ref_folding_dict.items():
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
		for k, v in ref_ribosome_subreactions.items():
			ref_cplx = v["enzyme"]
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[v["enzyme"]]
				defined_cplx = org_ribosome_subreactions[k]["enzyme"]
				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
					org_ribosome_subreactions[k]["enzyme"] = org_cplx
				else:
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
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.ribosome_subreactions. Curate the inputs in ribosome_subreactions.csv accordingly.'})
	def update_rrna_modifications_from_homology(self):
		ref_rrna_modifications = self.ref.rrna_modifications
		org_rrna_modifications = self.org.rrna_modifications
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_proteins = []
		for k, v in ref_rrna_modifications.items():
			ref_cplx = v["machine"]
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[v["machine"]]
				defined_cplx = org_rrna_modifications[k]["machine"]
				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
					org_rrna_modifications[k]["machine"] = org_cplx
				else:
					warn_proteins.append({
						'subreaction':k,
						'defined_complex':defined_cplx,
						'inferred_complex':org_cplx
					})
		# Warnings
		if warn_proteins:
			self.org.curation_notes['update_rrna_modifications_from_homology'].append({
				'msg':'Some enzymes defined in me_builder.org.rrna_modifications are different from the ones inferred from homology',
				'triggered_by':warn_proteins,
				'importance':'medium',
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.rrna_modifications. Curate the inputs in rrna_modifications.csv accordingly.'})
	def update_amino_acid_trna_synthetases_from_homology(self):
		ref_amino_acid_trna_synthetase = self.ref.amino_acid_trna_synthetase
		org_amino_acid_trna_synthetase = self.org.amino_acid_trna_synthetase
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_proteins = []
		for k, v in ref_amino_acid_trna_synthetase.items():
			ref_cplx = v
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[v]
				defined_cplx = org_amino_acid_trna_synthetase[k]
				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
					org_amino_acid_trna_synthetase[k] = org_cplx
				else:
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
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.amino_acid_trna_synthetase. Curate the inputs in amino_acid_trna_synthetase.csv accordingly.'})
	def update_peptide_release_factors_from_homology(self):
		ref_peptide_release_factors = self.ref.peptide_release_factors
		org_peptide_release_factors = self.org.peptide_release_factors
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		warn_proteins = []
		for k, v in ref_peptide_release_factors.items():
			ref_cplx = v['enzyme']
			if ref_cplx in ref_cplx_homolog:
				org_cplx = ref_cplx_homolog[ref_cplx]
				defined_cplx = org_peptide_release_factors[k]['enzyme']
				if not defined_cplx or defined_cplx in org_cplx or 'CPLX_dummy' in defined_cplx:
					org_peptide_release_factors[k]['enzyme'] = org_cplx
				else:
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
				'to_do':'Confirm whether the definitions or homology calls are correct in me_builder.org.peptide_release_factors. Curate the inputs in peptide_release_factors.csv accordingly.'})
	def update_initiation_subreactions_from_homology(self):
		ref_initiation_subreactions = self.ref.initiation_subreactions
		org_initiation_subreactions = self.org.initiation_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_initiation_subreactions.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_initiation_subreactions[k]["enzymes"]
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_elongation_subreactions_from_homology(self):
		ref_elongation_subreactions = self.ref.elongation_subreactions
		org_elongation_subreactions = self.org.elongation_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_elongation_subreactions.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_elongation_subreactions[k]["enzymes"]
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_termination_subreactions_from_homology(self):
		ref_termination_subreactions = self.ref.termination_subreactions
		org_termination_subreactions = self.org.termination_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_termination_subreactions.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_termination_subreactions[k]["enzymes"]
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_special_trna_subreactions_from_homology(self):
		ref_special_trna_subreactions = self.ref.special_trna_subreactions
		org_special_trna_subreactions = self.org.special_trna_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_special_trna_subreactions.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_special_trna_subreactions[k]["enzymes"]
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
		for k, v in ref_rna_degradosome.items():
			ref_cplxs = v['enzymes']
			defined_cplxs = org_rna_degradosome[k]['enzymes']
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_excision_machinery_from_homology(self):
		ref_excision_machinery = self.ref.excision_machinery
		org_excision_machinery = self.org.excision_machinery
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_excision_machinery.items():
			ref_cplxs = v['enzymes']
			defined_cplxs = org_excision_machinery[k]['enzymes']
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_special_modifications_from_homology(self):
		ref_special_trna_subreactions = self.ref.special_modifications
		org_special_trna_subreactions = self.org.special_modifications
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_special_trna_subreactions.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_special_trna_subreactions[k]["enzymes"]
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if v["stoich"]:
					org_special_trna_subreactions[k]["stoich"] = v["stoich"]
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_trna_modification_from_homology(self):
		ref_trna_modification = self.ref.trna_modification
		org_trna_modification = self.org.trna_modification
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_trna_modification.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_trna_modification[k]["enzymes"]
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if v["stoich"]:
					org_trna_modification[k]["stoich"] = v["stoich"]
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_transcription_subreactions_from_homology(self):
		ref_transcription_subreactions = self.ref.transcription_subreactions
		org_transcription_subreactions = self.org.transcription_subreactions
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_transcription_subreactions.items():
			ref_cplxs = v["enzymes"]
			defined_cplxs = org_transcription_subreactions[k]["enzymes"]
			org_cplxs = [
				ref_cplx_homolog[i] for i in ref_cplxs if i in ref_cplx_homolog
			]
			for i in org_cplxs:
				if i not in defined_cplxs:
					defined_cplxs.append(i)

	def update_translocation_pathways_from_homology(self):
		ref_translocation_pathways = self.ref.translocation_pathways
		org_translocation_pathways = self.org.translocation_pathways
		ref_cplx_homolog = self.homology.ref_cplx_homolog
		for k, v in ref_translocation_pathways.items():
			if k not in org_translocation_pathways:
				org_translocation_pathways[k] = v.copy()
			org_translocation_pathways[k]["enzymes"] = {}
			for ref_cplx, ref_dict in v["enzymes"].items():
				if ref_cplx in ref_cplx_homolog:
					org_cplx = ref_cplx_homolog[ref_cplx]
					org_translocation_pathways[k]["enzymes"][org_cplx] = ref_dict.copy()


	def update_m_model(self):
		org_model = self.org.m_model
		ref_model = self.ref.m_model

		for m in org_model.metabolites:
			if m.id not in ref_model.metabolites:
				continue
			ref_m = ref_model.metabolites.get_by_id(m.id)
			if not m.formula:
				m.formula = ref_m.formula

	def update_from_homology(self):
		self.update_enzyme_stoichiometry()
		self.update_protein_modification()
		self.update_TU_df()
		self.update_translocation_pathways_from_homology()
		self.protein_location_from_homology()
		self.update_translocation_multipliers()
		self.update_lipoprotein_precursors()
		self.update_cleaved_methionine()
		self.update_m_to_me_mets()
		self.update_generics_from_homology()
		self.update_folding_dict_from_homology()
		self.update_ribosome_subreactions_from_homology()
		self.update_rrna_modifications_from_homology()
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
		self.update_trna_modification_from_homology()
		self.update_m_model()

	def fill(self,
			 fill_with='CPLX_dummy'):
		fill_builder(self,fill_with='CPLX_dummy')

	def check(self):
		t_pathways = self.org.protein_location['translocase_pathway'].unique()

		file_t_paths = set()
		for t in t_pathways:
			for i in t:
				file_t_paths.add(i)

		defined_t_paths = set()
		for t in self.org.translocation_pathways.keys():
			defined_t_paths.add(coralme.builder.dictionaries.pathway_to_abbreviation[t])
		missing_pathways = file_t_paths - defined_t_paths
		if missing_pathways:
			self.org.curation_notes['check'].append({
									'msg':'Some translocase pathways in org.protein_location are not defined in org.translocation_pathways.',
									'triggered_by':missing_pathways,
									'importance':'high',
									'to_do':'Fill in translocation pathways in org.translocation_pathways or in translocation_pathways.csv'
				})

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
			v.to_csv(directory + k + '.csv')

	def build_me_model(self, overwrite = False):
		coralme.builder.main.MEReconstruction(self).build_me_model(overwrite = overwrite)

	def troubleshoot(self, growth_key_and_value = None):
		coralme.builder.main.METroubleshooter(self).troubleshoot(growth_key_and_value)

class MEReconstruction(object):
	"""
	MEReconstruction class for reconstructing a ME-model from user/automated input

	Parameters
	----------

	"""
	def __init__(self, builder):
		try:
			# only if builder.generate_files() was before builder.build_me_model()
			self.org = builder.org
			self.homology = builder.homology
		except:
			pass

		self.me_model = builder.me_model
		self.configuration = builder.configuration
		self.curation_notes = builder.curation_notes

		return None

	def input_data(self, m_model, overwrite = False):
		config = self.configuration

		# include rna_polymerases, lipids and lipoproteins from automated info and save new configuration file
		if config.get('rna_polymerases', None) is None or config.get('rna_polymerases') == {}:
			config['rna_polymerases'] = self.org.rna_polymerase_id_by_sigma_factor

			# replace IDs
			for name, rnap in config['rna_polymerases'].items():
				for key, value in rnap.items():
					config['rna_polymerases'][name][key] = self.homology.org_cplx_homolog.get(value, value.replace('-MONOMER', '_MONOMER'))

		if config.get('lipid_modifications', None) is None or len(config.get('lipid_modifications')) == 0:
			config['lipid_modifications'] = [ x for x in self.org.lipids if x.endswith('_p') and (x.startswith('pg') or x.startswith('pe')) and not x.startswith('pgp') ]

		if config.get('lipoprotein_precursors', None) is None or len(config.get('lipoprotein_precursors')) == 0:
			config['lipoprotein_precursors'] = self.org.lipoprotein_precursors

		if config.get('ngam', None) is None:
			config['ngam'] = self.org.NGAM

		if config.get('gam', None) is None:
			config['gam'] = self.org.GAM

		# modify options to not run again blast
		config['create_files'] = False
		config['run_bbh_blast'] = False
		config['dev_reference'] = False

		# set new options in the MEBuilder object
		self.configuration.update(config)

		with open('coralme-config.json', 'w') as outfile:
			anyconfig.dump(config, outfile)
		with open('coralme-config.yaml', 'w') as outfile:
			anyconfig.dump(config, outfile)
		#with open('automated.toml', 'w') as outfile:
			#anyconfig.dump(config, outfile)

		def read(filename, columns = []):
			if pathlib.Path(filename).is_file():
				df = coralme.builder.flat_files.read(filename)
				if sorted(df.columns) == sorted(columns):
					return df
				else:
					logging.warning('Column names in \'{:s}\' does not comply default values.'.format(filename))
			else:
				logging.warning('File \'{:s}\' does not exist.'.format(filename))
				return pandas.DataFrame(columns = columns)

		# INPUTS: We capture if the file exists or if it is empty
		# Transcriptional Units
		filename = config.get('df_TranscriptionalUnits', '')
		self.df_tus = read(filename, ['TU_id', 'replicon', 'genes', 'start', 'stop', 'tss', 'strand', 'rho_dependent', 'dnapol']).set_index('TU_id', inplace = False)

		# Reaction Matrix: reactions, metabolites, compartments, stoichiometric coefficientes
		filename = config.get('df_matrix_stoichiometry', '')
		self.df_rmsc = read(filename, ['Reaction', 'Metabolites', 'Compartment', 'Stoichiometry'])
		# SubReaction Matrix: subreactions, metabolites, compartments, stoichiometric coefficientes
		filename = config.get('df_matrix_subrxn_stoich', '')
		self.df_subs = read(filename, ['Reaction', 'Metabolites', 'Compartment', 'Stoichiometry'])
		# Orphan and Spontaneous reaction metadata
		filename = config.get('df_metadata_orphan_rxns', '')
		self.df_rxns = read(filename, ['name', 'description', 'is_reversible', 'is_spontaneous']).set_index('name', inplace = False)
		# Metabolites metadata
		filename = config.get('df_metadata_metabolites', '')
		self.df_mets = read(filename, ['id', 'me_id', 'name', 'formula', 'compartment', 'type']).set_index('id', inplace = False)

		# Drop-in replacement of input files:
		# step1: reactions, complexes, modification of complexes, enzyme-to-reaction mapping
		# step2a: generics, dnap stoichiometry, ribosome stoichiometry, degradosome stoichiometry, tRNA ligases, RNA modifications
		# step2b: folding pathways (DnaK, GroEL), N-terminal Methionine Processing, translocation pathways
		filename = config.get('df_gene_cplxs_mods_rxns', '')
		if pathlib.Path(filename).is_file() or overwrite == False:
			self.df_data = pandas.read_excel(filename).dropna(how = 'all')
		else:
			# generate a minimal dataframe from genbank and m-model
			self.df_data = coralme.builder.preprocess_inputs.generate_organism_specific_matrix(config['genbank-path'], model = m_model)
			# complete minimal dataframe with automated info
			self.df_data = coralme.builder.preprocess_inputs.complete_organism_specific_matrix(self, self.df_data, model = m_model, output = filename)

		# All other inputs and remove unnecessary genes from df_data
		return coralme.builder.preprocess_inputs.get_df_input_from_excel(self.df_data, self.df_rxns)

	def build_me_model(self, overwrite = False):
		config = self.configuration

		model = config.get('model_id', 'coralME')
		directory = config.get('log_directory', './')

		# ## Part 1: Create a minimum solveable ME-model
		logging.basicConfig(filename = '{:s}/step1-{:s}.log'.format(directory, model), filemode = 'w', level = logging.WARNING, format = log_format)
		logging.captureWarnings(True)

		# This will include the bare minimum representations of that still produce a working ME-model:
		# - Metabolic Reactions
		# - Transcription Reactions
		# - Translation Reactions
		# - Complex Formation Reactions
		#
		# ### 1) Create a MEModel object and populate its global info
		# This includes important parameters that are used to calculate coupling constraints as well as organism-specific information such as comparments, GC fraction, etc.

		me = self.me_model

		# update default options with user-defined values
		me.global_info.update(self.configuration)

		# Define the types of biomass components that will be synthesized in the model
		me.add_biomass_constraints_to_model(me.global_info['biomass_constraints'])

		# Define ME-model compartments
		me.compartments = me.global_info['compartments']

		# Define M-Model
		if me.global_info['m-model-path'].endswith('.json'):
			me.gem = cobra.io.load_json_model(me.global_info['m-model-path'])
		else:
			me.gem = cobra.io.read_sbml_model(me.global_info['m-model-path'])

		# Read user inputs
		df_data, df_rxns, df_cplxs, df_ptms, df_enz2rxn, df_rna_mods, df_protloc, df_transpaths = coralme.builder.main.MEReconstruction.input_data(self, me.gem)

		# update default options with missing, automated-defined values
		me.global_info.update(self.configuration)

		df_tus = self.df_tus
		df_rmsc = self.df_rmsc
		df_subs = self.df_subs
		df_mets = self.df_mets

		# Remove default ME-model SubReactions that are not mapped in the organism-specific matrix
		subrxns = set(df_data[df_data['ME-Model SubReaction'].notnull()]['ME-Model SubReaction'])

		# list of subreactions
		for key in ['peptide_processing_subreactions']:
			me.global_info[key] = list(set(me.global_info[key]).intersection(subrxns))

		# dictionary mapping subreactions and stoichiometry
		for key in ['transcription_subreactions', 'translation_subreactions']:
			me.global_info[key] = { k:v for k,v in me.global_info[key].items() if k in subrxns }

		# ### 2) Load metabolites and build Metabolic reactions
		#
		# It creates a new M-Model, then incorporates it into the ME-model using the *add_m_model_content* function. Reactions are added directly and metabolites are added as *StoichiometricData* (*me.stoichiometric_data*).
		#
		# Different metabolite types have different properties in a ME-model, so complexes are added to the model as a *ComplexData*, not as a *Metabolite*. Components in the M-Model that are actually *Complexes* are compiled in the *cplx_lst* variable.

		# Modify M-Model
		m_model = coralme.builder.flat_files.process_m_model(
			m_model = me.gem,
			rxns_data = df_rxns, # metadata of new reactions
			mets_data = df_mets, # metadata of new metabolites
			cplx_data = df_cplxs, # metadata of prot-prot/prot-rna complexes
			reaction_matrix = df_rmsc, # stoichiometric coefficients of new! reactions

			me_compartments = { v:k for k,v in me._compartments.items() }, repair = True,
			defer_to_rxn_matrix = me.global_info['defer_to_rxn_matrix'])

		# Some of the 'metabolites' in the M-Model are actually complexes.
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
		# To construct the bare minimimum components of a transcription and translation reactions.
		# For example, transcription reactions at this point include nucleotides and the synthesized RNAs.

		lst = set(df_data['Gene Locus ID'].str.replace('protein_', '').str.replace('RNA_', '').tolist())

		coralme.util.building.build_reactions_from_genbank(
			me_model = me, gb_filename = me.global_info['genbank-path'],
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

		# Remove modifications. They will be added back in later (See Building Step 2, Part 3).
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
			coralme.core.processdata.GenericData(generic, me, components).create_reactions()

		# ### 6) Add dummy reactions to model and the *unmodeled_protein_fraction* constraint
		#
		# This includes the Transcription, Translation, ComplexFormation, and Metabolic reactions for a dummy RNA/protein/complex. Sequence for *dummy RNA* is based on the prevalance of each codon found in the genbank file.

		coralme.util.building.add_dummy_reactions(me, update = True)

		# The 'dummy protein' is associated to orphan reactions.
		# This ensures that orphan reactions will not become favored to fulfil unmodeled protein fraction requirement.
		rxn = coralme.core.reaction.SummaryVariable('dummy_protein_to_mass')
		me.add_reactions([rxn])
		mass = me.metabolites.protein_dummy.formula_weight / 1000. # in kDa
		met = coralme.core.component.Constraint('unmodeled_protein_biomass')
		rxn.add_metabolites({'protein_biomass': -mass, 'protein_dummy': -1, met: mass})

		# ### 7) Associate Complex(es) to Metabolic reactions and build the ME-model metabolic network

		# Associate a reaction id with the ME-model complex id (including modifications)
		rxn_to_cplx_dict = coralme.builder.flat_files.get_reaction_to_complex(m_model, df_enz2rxn)
		spontaneous_rxns = [me.global_info['dummy_rxn_id']] + list(df_rxns[df_rxns['is_spontaneous'] == True].index.values)

		coralme.util.building.add_reactions_from_stoichiometric_data(
			me, rxn_to_cplx_dict, is_spontaneous = spontaneous_rxns, update = True)

		# ### 8) Incorporate remaining biomass constituents
		# ### 1. General Demand Requirements
		# There are leftover components from the biomass equation that either:
		# 1. have no mechanistic function in the model (e.g., *glycogen*)
		# 2. are cofactors that are regenerated (e.g., *nad*)
		#
		# Applies demands and coefficients from the biomass objective function from the M-Model.

		for key in [ 'gam', 'ngam', 'unmodeled_protein_fraction' ]:
			if key in me.global_info:
				setattr(me, key, me.global_info[key])

		biomass_constituents = me.global_info['flux_of_biomass_constituents']

		rxn = coralme.core.reaction.SummaryVariable('biomass_constituent_demand')
		me.add_reactions([rxn])
		rxn.add_metabolites({ k:-(abs(v)) for k,v in biomass_constituents.items() })
		rxn.lower_bound = me.mu # coralme.util.mu
		rxn.upper_bound = me.mu # coralme.util.mu
		constituent_mass = sum(me.metabolites.get_by_id(c).formula_weight / 1000. * abs(v) for c,v in biomass_constituents.items())
		rxn.add_metabolites({me.metabolites.constituent_biomass: constituent_mass})

		# ### 2. Lipid Demand Requirements
		# Metabolites and coefficients from biomass objective function

		lipid_demand = {}
		for key, value in me.global_info['flux_of_lipid_constituents'].items():
			lipid_demand[key] = abs(value)

		for met, requirement in lipid_demand.items():
			try:
				component_mass = me.metabolites.get_by_id(met).formula_weight / 1000.
				rxn = coralme.core.reaction.SummaryVariable('DM_' + met)
				me.add_reactions([rxn])
				rxn.add_metabolites({met: -1 * requirement, 'lipid_biomass': component_mass * requirement})
				rxn.lower_bound = me.mu # coralme.util.mu
				rxn.upper_bound = 1000. # coralme.util.mu?
			except:
				msg = 'Metabolite \'{:s}\' lacks a formula. Please correct it in the M-Model or the \'metabolites.txt\' metadata file.'
				logging.warning(msg.format(met))

		# ### 3. DNA Demand Requirements
		# Added based on growth rate dependent DNA levels as in [O'brien EJ et al 2013](https://www.ncbi.nlm.nih.gov/pubmed/24084808) (*E. coli* data)

		dna_demand_bound = coralme.builder.dna_replication.return_gr_dependent_dna_demand(
			me, me.global_info['GC_fraction'], me.global_info['percent_dna_data'], me.global_info['gr_data_doublings_per_hour'])

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
		dna_replication.lower_bound = dna_demand_bound
		dna_replication.upper_bound = dna_demand_bound

		# ### 9) Save ME-model as a pickle file

		import pickle
		with open('MEModel-step1-{:s}.pkl'.format(model), 'wb') as outfile:
			pickle.dump(me, outfile)

		print('ME-model was saved as MEModel-step1-{:s}.pkl'.format(model))

		# ## Part 2: Add metastructures to solving ME-model
		logging.basicConfig(filename = '{:s}/step2-{:s}.log'.format(directory, model), filemode = 'w', level = logging.WARNING, format = log_format)
		logging.captureWarnings(True)

		# ### 1) Add *GenericData* and its reactions
		# Multiple entities can perform the same role. To prevent a combinatorial explosion, we create "generic" versions of the components, where any of those entities can fill in.
		#
		# This step was moved to the first part of the builder process to be able to use generics in enzymes linked to metabolic reactions.

		# ### 2) Add DNA polymerase and the ribosome and its rRNAs modifications
		# ~~This uses the ribosome composition and subreaction definitions in **coralme/ribosome.py**~~

		# WARNING: EXPERIMENTAL
		# dnapol_id = me.global_info['dnapol_id']

		# me.add_metabolites([cobrame.core.component.Complex(dnapol_id)])
		# data = cobrame.core.processdata.ComplexData(dnapol_id, me)
		# data.stoichiometry.update(coralme.preprocess_inputs.dnapolymerase_stoichiometry(df_data))
		# data.create_complex_formation(verbose = False)

		ribosome_stoich = coralme.builder.preprocess_inputs.ribosome_stoichiometry(df_data)
		ribosome_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Ribosome')
		df_rrna_mods = df_rna_mods[df_rna_mods['bnum'].isin(['16S_rRNAs', '23S_rRNAs'])]
		coralme.builder.ribosome.add_ribosome(me, ribosome_stoich, ribosome_subreactions, df_rrna_mods, verbose = True)

		# ### 3) Add charged tRNA reactions

		# The tRNA charging reactions were automatically added when loading the genome from the GenBank file. However, the charging reactions still need to be made aware of the tRNA synthetases which are responsible. Generic charged tRNAs are added to translation reactions via *SubreactionData* below.

		aa_synthetase_dict = coralme.builder.preprocess_inputs.aa_synthetase_dict(df_data)
		for data in tqdm.tqdm(list(me.tRNA_data), 'Adding tRNA synthetase(s) information into the ME-model...', bar_format = bar_format):
			data.synthetase = str(aa_synthetase_dict.get(data.amino_acid, 'CPLX_dummy'))

		special_trna_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Special_tRNA')
		coralme.builder.translation.add_charged_trna_subreactions(me, translation_stop_dict = me.global_info['translation_stop_dict'])

		# ### 4) Add tRNA modifications into the ME-model and asocciate them with tRNA charging reactions

		# Add tRNA modifications to ME-model
		df_trna_mods = df_rna_mods[~df_rna_mods['bnum'].isin(['16S_rRNAs', '23S_rRNAs'])]
		coralme.builder.trna_charging.add_trna_modification_procedures(me, df_trna_mods)

		# Associate tRNA modifications to tRNAs
		trna_modifications = coralme.builder.flat_files.get_trna_modification_targets(df_trna_mods)
		for idx, trna in tqdm.tqdm(list(df_trna_mods.iterrows()), 'Associating tRNA modification enzyme(s) to tRNA(s)...', bar_format = bar_format):
			for data in me.process_data.query(trna['bnum']):
				data.subreactions = trna_modifications[trna['bnum']]

		# ### 5) Add translation SubReactions into the ME-model and associate them with Translation reactions

		initiation_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_initiation')
		elongation_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_elongation')
		termination_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Translation_termination')
		processing_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Protein_processing')

		for data in tqdm.tqdm(list(me.translation_data), 'Adding SubReactions into TranslationReactions...', bar_format = bar_format):
			data.add_initiation_subreactions(start_codons = me.global_info['start_codons'], start_subreactions = initiation_subreactions)
			data.add_elongation_subreactions(elongation_subreactions = elongation_subreactions)
			data.add_termination_subreactions(translation_terminator_dict = me.global_info['translation_stop_dict'])

		# ### 6) Add Transcription Metacomplexes: RNA Polymerase(s)

		rna_polymerases = me.global_info.get('rna_polymerases', {})

		# Create polymerase "metabolites"
		for rnap in tqdm.tqdm(rna_polymerases.keys(), 'Adding RNA Polymerase(s) into the ME-model...', bar_format = bar_format):
			rnap_obj = coralme.core.component.RNAP(rnap)
			me.add_metabolites(rnap_obj)

		# Add polymerase complexes in the model
		coralme.builder.transcription.add_rna_polymerase_complexes(me, rna_polymerases, verbose = False)

		# Associate the correct RNA_polymerase and factors to TUs
		for tu_id in tqdm.tqdm(df_tus.index, 'Associating a RNA Polymerase to each Transcriptional Unit...', bar_format = bar_format):
			try:
				transcription_data = me.process_data.get_by_id(tu_id)
				transcription_data.RNA_polymerase = df_tus['dnapol'][tu_id]
			except KeyError as e:
				logging.warning('Transcription Unit \'{:s}\' is missing from ProcessData. Check if it is the correct behavior.'.format(tu_id))
				pass

		# ### 7) Add Transcription Metacomplexes: Degradosome (both for RNA degradation and RNA splicing)

		degradosome_id = me.global_info['degradosome_id']

		me.add_metabolites([coralme.core.component.Complex(degradosome_id)])
		data = coralme.core.processdata.ComplexData(degradosome_id, me)
		data.stoichiometry.update(coralme.builder.preprocess_inputs.degradosome_stoichiometry(df_data))
		data.create_complex_formation(verbose = False)

		# Used for RNA splicing
		data = coralme.core.processdata.SubreactionData('RNA_degradation_machine', me)
		data.enzyme = me.global_info['degradosome_id']

		# .25 water equivalent for atp hydrolysis per nucleotide
		data = coralme.core.processdata.SubreactionData('RNA_degradation_atp_requirement', me)
		data.stoichiometry = { 'atp_c': -0.25, 'h2o_c': -0.25, 'adp_c': 0.25, 'pi_c': 0.25 }

		for excision_type in me.global_info['excision_machinery']:
			stoichiometry = coralme.builder.preprocess_inputs.excision_machinery_stoichiometry(df_data, excision_type)
			if stoichiometry is not None:
				coralme.builder.transcription.add_rna_excision_machinery(me, excision_type, stoichiometry)

		# add excision machineries into TranscriptionData
		coralme.builder.transcription.add_rna_splicing(me)

		# ## Part 3: Add remaining modifications (including iron clusters and lipoate)

		# mods_dct is a dictionary {complex_name_with_mods: {core_enzyme: complex_name, modifications: {stoichiometry}}
		mods_dct = coralme.builder.flat_files.get_complex_modifications(
			reaction_matrix = df_rmsc,
			protein_complexes = df_cplxs,
			complex_mods = df_ptms,
			compartments = { v:k for k,v in me._compartments.items() })

		for complex_id, info in tqdm.tqdm(mods_dct.items(), 'Processing ComplexData in ME-model...', bar_format = bar_format):
			modifications = {}
			for mod, value in info['modifications'].items():
				# stoichiometry of modification determined in subreaction_data.stoichiometry
				modifications['mod_' + mod] = abs(value)
			me.process_data.get_by_id(complex_id).subreactions = modifications

		# Check if the modification is set on any component in the organism-specific matrix
		if me.process_data.has_id('mod_btn_c'):
			coralme.builder.modifications.add_biotin_modifications(me)
		if me.process_data.has_id('mod_lipoyl_c'):
			coralme.builder.modifications.add_lipoate_modifications(me)
		if me.process_data.has_id('mod_3fe4s_c') or me.process_data.has_id('mod_4fe4s_c') or me.process_data.has_id('mod_2fe2s_c'):
			coralme.builder.modifications.add_iron_sulfur_modifications(me)
		if me.process_data.has_id('mod_FeFe_cofactor_c') or me.process_data.has_id('mod_NiFe_cofactor_c'):
			coralme.builder.modifications.add_FeFe_and_NiFe_modifications(me)
		if me.process_data.has_id('mod_bmocogdp_c'):
			coralme.builder.modifications.add_bmocogdp_modifications(me)

		# add formation reactions for each of the ComplexData
		for data in tqdm.tqdm(list(me.complex_data), 'Adding ComplexFormation into the ME-model...', bar_format = bar_format):
			formation = data.formation
			if formation:
				formation.update()
			else:
				data.create_complex_formation()
				logging.warning('Added ComplexFormation for \'{:s}\'.'.format(data.id))

		# ## Part 4: Add remaining subreactions
		# ### 1. Add translation related subreactions

		# get list of processed proteins from df_data
		processing_pathways = [
			x.replace('Protein_processing_', '')
			for x in me.global_info['translation_subreactions'].keys()
			if x.startswith('Protein_processing')
			]
		protein_processing = { k:df_data[~df_data[k].isnull() & df_data[k]]['Gene Locus ID'].tolist() for k in processing_pathways }

		# Add the translation subreaction data objects to the ME-model
		coralme.builder.translation.add_subreactions_to_model(
			me, [initiation_subreactions, elongation_subreactions, termination_subreactions, processing_subreactions])

		for data in tqdm.tqdm(list(me.translation_data), 'Adding SubReactions into TranslationReactions...', bar_format = bar_format):
			for process in protein_processing:
				if data.id in protein_processing[process]:
					data.subreactions['Protein_processing_' + process] = 1

			# This block was run above, but it should be run again to incorporate any subreactions not added previously
			data.add_initiation_subreactions(start_codons = me.global_info['start_codons'], start_subreactions = initiation_subreactions)
			data.add_elongation_subreactions(elongation_subreactions = elongation_subreactions)
			data.add_termination_subreactions(translation_terminator_dict = me.global_info['translation_stop_dict'])

			# Add organism specific subreactions associated with peptide processing
			for subrxn in me.global_info['peptide_processing_subreactions']:
				data.subreactions[subrxn] = 1

		# ### 2) Add transcription related subreactions

		transcription_subreactions = coralme.builder.preprocess_inputs.get_subreactions(df_data, 'Transcription')
		coralme.builder.transcription.add_subreactions_to_model(me, [transcription_subreactions])

		for transcription_data in tqdm.tqdm(list(me.transcription_data), 'Adding Transcription SubReactions...', bar_format = bar_format):
			# Assume false if not in tu_df
			rho_dependent = df_tus.rho_dependent.get(transcription_data.id, False)
			rho = 'dependent' if rho_dependent else 'independent'
			stable = 'stable' if transcription_data.codes_stable_rna else 'normal'
			transcription_data.subreactions['Transcription_{:s}_rho_{:s}'.format(stable, rho)] = 1

		# ## Part 5: Add in Translocation reactions

		v1 = { 'length_dependent' : True, 'fixed_keff' : False } # default
		v2 = { 'length_dependent' : False, 'fixed_keff' : True } # only for FtsY in the SRP pathway
		v3 = { 'length_dependent' : False, 'fixed_keff' : False } # for all the lol and bam pathway enzymes

		for key, value in me.global_info['translocation_pathway'].items():
			if 'translocation_pathway_' + key in df_transpaths.index:
				value['enzymes'] = {
					k:(v2 if k == value.get('FtsY', None) else v1 if (key.lower() not in ['lol', 'bam']) else v3) \
						for k in df_transpaths.loc['translocation_pathway_' + key].tolist()[0] }

			# TO ADD PATHWAYS WITHOUT HOMOLOGS
			# TODO: Check if the user wants to add dummies to the translocation pathways
			elif value.get('enzymes', None) is None:
				value['enzymes'] = { 'CPLX_dummy':(v2 if value.get('FtsY', None) else v1 if (key.lower() not in ['lol', 'bam']) else v3) }

		dct = { k:v['abbrev'] for k,v in me.global_info['translocation_pathway'].items() }
		dct = dict([(v, [k + '_translocation' for k,v1 in dct.items() if v1 == v]) for v in set(dct.values())])
		dct = { k:(v[0] if len(v) == 1 else v) for k,v in dct.items() } # tat pathways should be a list, but not the others

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
		coralme.builder.translocation.add_translocation_pathways(me, pathways_df = df_protloc, abbreviation_to_pathway = dct, membrane_constraints = False)

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
			complex_data.formation.update()

			# Complex IDs in protein compartment file don't include modifications
			# Some have multiple alternative modifications so must loop through these
			for complex_data in me.process_data.query('^{:s}_mod_'.format(cplx)):
				complex_data.stoichiometry.update(new_stoich[cplx])
				complex_data.formation.update()

		# ## Part 6: Add Cell Wall Components
		# ### 1. Add lipid modification SubreactionData

		compartment_dict = {}
		for idx, compartment in df_protloc.set_index('Protein').Protein_compartment.to_dict().items():
			compartment_dict[idx] = compartment

		lipid_modifications = me.global_info.get('lipid_modifications')
		lipoprotein_precursors = me.global_info.get('lipoprotein_precursors')

		# TODO: associate subreactions of lipid modifications with enzyme

		# Step2: add reactions of lipoprotein formation
		coralme.builder.translocation.add_lipoprotein_formation(
			me, compartment_dict, lipoprotein_precursors, lipid_modifications, membrane_constraints = False, update = True)

		# ### 2. Correct complex formation IDs if they contain lipoproteins

		# TODO:
		#for gene in tqdm.tqdm(coralme.builder.translocation.lipoprotein_precursors.values()):
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

		if len(me.global_info['braun\'s_lipoprotein']) >= 1:
			for Braun_lipoprotein in me.global_info['braun\'s_lipoprotein']:
				rxn = coralme.core.reaction.SummaryVariable('core_structural_demand_brauns_{:s}'.format(Braun_lipoprotein))
				murein5px4p = me.metabolites.get_by_id(me.global_info['braun\'s_lipid_mod'])
				murein5px4p_mass = murein5px4p.formula_weight / 1000.
				lipoprotein = me.metabolites.get_by_id('protein_{:s}_lipoprotein_Outer_Membrane'.format(Braun_lipoprotein))
				me.add_reactions([rxn])

				# biomass of lipoprotein accounted for in translation and lipid modification
				rxn.add_metabolites({
					murein5px4p : -abs(me.global_info['braun\'s_murein_flux']),
					lipoprotein : -abs(me.global_info['braun\'s_lpp_flux']),
					me.metabolites.peptidoglycan_biomass : abs(me.global_info['brauns_murein_flux']) * murein5px4p_mass
					},
					combine = False)
				rxn.lower_bound = me.mu # coralme.util.mu
				rxn.upper_bound = me.mu # coralme.util.mu

		else:
			logging.warning('No Braun\'s lipoprotein homolog was set. Please check if it is the correct behavior.')

		# ## Part 7: Set keffs
		# Either entirely based on SASA or using fit keffs from [Ebrahim et al 2016](https://www.ncbi.nlm.nih.gov/pubmed/27782110?dopt=Abstract)

		# ## Part 8: Model updates and corrections
		# ### 1. Subsystems

		# Add reaction subsystems from M-Model to ME-model
		for rxn in tqdm.tqdm(me.gem.reactions, 'Adding reaction subsystems from M-Model into the ME-model...', bar_format = bar_format):
			if rxn.id in me.process_data:
				data = me.process_data.get_by_id(rxn.id)
			else:
				continue

			for parent_rxn in data.parent_reactions:
				rxn.subsystem = parent_rxn.subsystem

		# ### 2. Add enzymatic coupling for "carriers"
		# These are enzyme complexes that act as metabolites in a metabolic reaction.

		for data in tqdm.tqdm(list(me.stoichiometric_data), 'Processing StoichiometricData in ME-model...', bar_format = bar_format):
			if data.id == 'dummy_reaction':
				continue

			for met, value in data.stoichiometry.items():
				if not isinstance(me.metabolites.get_by_id(met), coralme.core.component.Complex) or value > 0:
					continue

				subreaction_id = met + '_carrier_activity'
				if subreaction_id not in me.process_data:
					sub = coralme.core.processdata.SubreactionData(met + '_carrier_activity', me)
					sub.enzyme = met

				data.subreactions[subreaction_id] = abs(value)

		# ### 3. Add remaining complex formulas and compartments to the ME-model

		# Update a second time to incorporate all of the metabolite formulas correctly
		for r in tqdm.tqdm(me.reactions.query('formation_'), 'Updating all FormationReactions...', bar_format = bar_format):
			r.update()

		# Update complex formulas
		modification_formulas = df_mets[df_mets['type'].str.match('MOD')]
		modification_formulas = dict(zip(modification_formulas['me_id'], modification_formulas['formula']))
		# This will add the formula to complexes not formed from a complex formation reaction (e.g. CPLX + na2_c -> CPLX_mod_na2(1))
		coralme.builder.formulas.add_remaining_complex_formulas(me, modification_formulas)

		# Update reactions affected by formula update
		for r in tqdm.tqdm(me.reactions.query('_mod_lipoyl'), 'Updating FormationReactions involving a lipoyl prostethic group...', bar_format = bar_format):
			r.update()

		for r in tqdm.tqdm(me.reactions.query('_mod_glycyl'), 'Updating FormationReactions involving a glycyl prostethic group...', bar_format = bar_format):
			r.update()

		# add metabolite compartments
		coralme.builder.compartments.add_compartments_to_model(me)

		# ## Part 9: Save

		import pickle

		# 26.5 seconds without constraints
		me.update() # ~15 mins without prefilter
		me.prune() # ~3.5 mins without prefilter

		with open('MEModel-step2-{:s}.pkl'.format(model), 'wb') as outfile:
			pickle.dump(me, outfile)

		print('ME-model was saved as MEModel-step2-{:s}.pkl'.format(model))

		n_genes = len(me.metabolites.query(re.compile('^RNA_(?!biomass|dummy|degradosome)')))
		new_genes = n_genes * 100. / len(me.gem.genes) - 100
		print('Done. Number of genes in the ME-model is {:d} (+{:.2f}%, from {:d})'.format(n_genes, new_genes, len(me.gem.genes)))

		self.me_model = me

		return None

class METroubleshooter(object):
	"""
	METroubleshooter class for troubleshooting growth in a ME-model

	This class contains methods to identify gaps and obtain a feasible ME-model.

	Parameters
	----------
	MEBuilder : coralme.builder.main.MEBuilder
	"""

	def __init__(self, builder):
		self.me_model = builder.me_model
		self.configuration = builder.configuration
		self.curation_notes = builder.curation_notes

	def gap_find(self):
		#from draft_coralme.util.helper_functions import find_gaps
		print('  '*5 + 'Finding gaps from the M-Model only...')
		m_gaps = coralme.builder.helper_functions.find_gaps(self.me_model.gem)

		print('  '*5 + 'Finding gaps in the ME-model...')
		me_gaps = coralme.builder.helper_functions.find_gaps(self.me_model, growth_key = self.me_model.mu)

		idx = list(set(me_gaps.index) - set(m_gaps.index))
		new_gaps = me_gaps.loc[idx]

		filt1 = new_gaps['p'] == 1
		filt2 = new_gaps['c'] == 1
		filt3 = new_gaps['u'] == 1

		deadends = list(new_gaps[filt1 | filt2 | filt3].index)
		deadends = sorted([ x for x in deadends if 'biomass' not in x if not x.endswith('_e') ])

		print('  '*5 + '{:d} metabolites were identified as deadends.'.format(len(deadends)))
		for met in deadends:
			name = self.me_model.metabolites.get_by_id(met).name
			print('  '*6 + '{:s}: {:s}'.format(met, 'Missing metabolite in the M-Model.' if name == '' else name))
		return deadends

	def gap_fill(self, deadends = [], growth_key_and_value = { sympy.Symbol('mu', positive = True) : 0.01 }, met_types = 'Metabolite'):
		if len(deadends) != 0:
			print('  '*5 + 'Adding a sink reaction for each identified deadend metabolite...')
			coralme.builder.helper_functions.add_exchange_reactions(self.me_model, deadends)
		else:
			print('  '*5 + 'Empty set of deadends metabolites to test.')
			return None

		print('  '*5 + 'Optimizing gapfilled ME-model...', end = '')

		if self.me_model.feasibility(keys = growth_key_and_value):
			print(' The ME-model is feasible.')
			print('  '*5 + 'Gapfilled ME-model is feasible with growth rate {:g} 1/h.'.format(list(growth_key_and_value.values())[0]))
			return True
		else:
			print(' The ME-model is not feasible.')
			print('  '*5 + 'Provided set of sink reactions for deadend metabolites does not allow growth.')
			return False

	def brute_check(self, growth_key_and_value, met_types = 'Metabolite'):
		if isinstance(met_types, str):
			met_types = [met_types]

		mets = []
		for met_type in met_types:
			for met in self.me_model.metabolites:
				filter1 = type(met) == getattr(coralme.core.component, met_type)
				filter2 = met.id.startswith('trna')
				filter3 = met.id.endswith('trna_c')
				filter4 = met.id.endswith('_e')
				if filter1 and not filter2 and not filter3 and not filter4:
					mets.append(met.id)

# 		if met_types == 'Metabolite':
# 			# remove metabolites that are fed into the model through transport reactions
# 			medium = set([ '{:s}_c'.format(x[3:-2]) for x in self.me_model.gem.medium.keys() ])
# 			mets = set(mets).difference(medium)

# 			# filter out manually
# 			mets = set(mets).difference(set(['ppi_c', 'ACP_c', 'h_c']))
# 			mets = set(mets).difference(set(['ade_c']))
# 			mets = set(mets).difference(set(['adp_c', 'amp_c', 'atp_c']))
# 			mets = set(mets).difference(set(['cdp_c', 'cmp_c', 'ctp_c']))
# 			mets = set(mets).difference(set(['gdp_c', 'gmp_c', 'gtp_c']))
# 			mets = set(mets).difference(set(['udp_c', 'ump_c', 'utp_c']))
# 			mets = set(mets).difference(set(['dadp_c', 'dcdp_c', 'dgdp_c', 'dtdp_c', 'dudp_c']))
# 			mets = set(mets).difference(set(['damp_c', 'dcmp_c', 'dgmp_c', 'dtmp_c', 'dump_c']))
# 			mets = set(mets).difference(set(['datp_c', 'dctp_c', 'dgtp_c', 'dttp_c', 'dutp_c']))
# 			mets = set(mets).difference(set(['nad_c', 'nadh_c', 'nadp_c', 'nadph_c']))
# 			mets = set(mets).difference(set(['5fthf_c', '10fthf_c', '5mthf_c', 'dhf_c', 'methf_c', 'mlthf_c', 'thf_c']))
# 			mets = set(mets).difference(set(['fad_c', 'fadh2_c', 'fmn_c']))
# 			mets = set(mets).difference(set(['coa_c']))

		bf_gaps = coralme.builder.helper_functions.brute_force_check(self.me_model, mets, growth_key_and_value)
		if self.me_model.feasibility(keys = growth_key_and_value):
			return bf_gaps, True
		else:
			return bf_gaps, False

	def troubleshoot(self, growth_key_and_value = None):
		if growth_key_and_value is None:
			growth_key_and_value = { self.me_model.mu : 0.01 }

		growth_key, growth_value = zip(*growth_key_and_value.items())

		print('~ '*1 + 'Troubleshooting started...')
		print('  '*1 + 'Checking if the ME-model can simulate growth without gapfilling reactions...', end = '')
		if self.me_model.feasibility(keys = growth_key_and_value):
			print('  '*5 + '\nOriginal ME-model is feasible with a tested growth rate of {:f} 1/h'.format(list(keys.values)[0]))
			return None
		else:
			print(' FALSE.')
			works = False

		# Step 1. Find topological gaps
		print('~ '*1 + 'Step 1. Find topological gaps in the ME-model.')
		deadends = self.gap_find()

		if len(deadends) != 0:
			self.curation_notes['troubleshoot'].append({
				'msg':'Some deadends were identified',
				'triggered_by':deadends,
				'importance':'high',
				'to_do':'Fix metabolic deadends by adding reactions or solving other warnings.'})

		# Step 2. Test feasibility adding all topological gaps
		if len(deadends) != 0:
			print('~ '*1 + 'Step 2. Solve gap-filled ME-model with all identified deadend metabolites.')
			print('  '*5 + 'Attempt optimization gapfilling the identified metabolites from Step 1')
			works = self.gap_fill(deadends = deadends, growth_key_and_value = growth_key_and_value)

		if len(deadends) == 0 and works == False:
			print('~ '*1 + 'Step 2. Solve gap-filled ME-model with provided sink reactions for deadend metabolites.')
		if works == False:
			met_type = 'Metabolite'
			print('  '*5 + 'Checking reactions that provide components of type \'{:s}\' using brute force...'.format(met_type))
			bf_gaps, works = self.brute_check(growth_key_and_value, met_types = met_type)

			if len(bf_gaps) != 0:
				self.curation_notes['troubleshoot'].append({
					'msg':'Additional deadends were identified',
					'triggered_by':bf_gaps,
					'importance':'high',
					'to_do':'Fix metabolic deadends by adding reactions or solving other warnings.'})

		# Step 3. Test different sets of MEComponents
		e_gaps = []
		if works == False:
			print('~ '*1 + 'Step 3. Attempt gapfilling different groups of E-matrix components.')

			met_types = [ 'Complex', 'GenerictRNA', 'TranscribedGene', 'TranslatedGene', 'ProcessedProtein', 'GenericComponent' ]

			for met_type in met_types:
				print('  '*5 + 'Gapfill reactions to provide components of type \'{:s}\' using brute force.'.format(mt))

				self.me_model.relax_bounds()
				self.me_model.reactions.protein_biomass_to_biomass.lower_bound = growth_value[0]/100 # Needed to enforce protein production

				bf_gaps, works = self.brute_check(growth_key_and_value, met_types = met_type)
				if works:
					e_gaps.append(bf_gaps)
					break

		if works: # Meaning it can grow in any case
			if isinstance(e_gaps, list) and e_gaps:
				self.curation_notes['troubleshoot'].append({
					'msg':'Some metabolites are necessary for growth',
					'triggered_by':[ x for y in e_gaps for x in y ],
					'importance':'critical',
					'to_do':'Fix the gaps by adding reactions or solving other warnings. If some items are from the E-matrix, fix these first!'})

			# delete added sink reactions with lb == 0 and ub == 0
			for rxn in self.me_model.reactions.query('^SK_'):
				if rxn.lower_bound == 0 and rxn.upper_bound == 0:
					self.me_model.remove_reactions([rxn])

			print('~ '*1 + 'Final step. Fully optimizing with precision 1e-6 and save solution into the ME-model...')
			self.me_model.optimize(max_mu = 3.0, precision = 1e-6)
			print('  '*1 + 'Gapfilled ME-model is feasible with growth rate {:f}.'.format(self.me_model.solution.objective_value))

			with open('MEModel-step3-{:s}-TS.pkl'.format(self.me_model.id), 'wb') as outfile:
				pickle.dump(self.me_model, outfile)

			print('\nME-model was saved as MEModel-step3-{:s}-Troubleshooted.pkl'.format(self.me_model.id))
		else:
			print('~ '*1 + 'METroubleshooter failed to determine a set of problematic metabolites.')

		return None