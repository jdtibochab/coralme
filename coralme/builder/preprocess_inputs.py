#!/usr/bin/python3
import numpy
import pathlib
import pandas
import xlsxwriter
import warnings

import logging
log = logging.getLogger(__name__)

from Bio import SeqIO

def _save_to_excel(data, output):
	#if overwrite:
	try:
		pathlib.Path(output).unlink(missing_ok = True) # python>=3.8
	except:
		if pathlib.Path(output).exists():
			pathlib.Path(output).unlink() # python==3.7

	with open(output, 'wb') as outfile:
		writer = pandas.ExcelWriter(outfile, engine = 'xlsxwriter')
		data.to_excel(writer, index = False, freeze_panes = (1, 8))
		(max_row, max_col) = data.shape

		# Get the xlsxwriter workbook and worksheet objects.
		workbook  = writer.book
		worksheet = writer.sheets['Sheet1']

		# Set the autofilter.
		worksheet.autofilter(0, 0, max_row, max_col - 1)

		# Make the columns wider for clarity.
		worksheet.set_column_pixels(0,  max_col - 1, 96)

		# Close the Pandas Excel writer and output the Excel file.
		writer.close()

def generate_organism_specific_matrix(genbank, locus_tag, model):
	contigs = []
	for contig in SeqIO.parse(genbank, 'genbank'):
		contigs.append(contig)

	# get all features
	lst = [ x for y in [ x.features for x in contigs ] for x in y ]
	lst = [ x for x in lst if x.type in [ 'CDS', 'rRNA', 'tRNA', 'ncRNA', 'tmRNA', 'misc_RNA' ] ]

	# create a pandas DataFrame with organism-specific information to be completed with the builder data
	df = pandas.DataFrame(columns = [
		'Gene Locus ID',
		'Definition',
		'Gene Names',
		'Old Locus Tag',
		'BioCyc',
		'Reference BBH',
		'Reference Definition',
		'Feature Type',
		'Complex Name',
		'Complex ID',
		'Cofactors in Modified Complex',
		'Generic Complex ID',
		'MetaComplex ID',
		'ME-model SubReaction',
		'M-model Reaction ID',
		'Reaction Name',
		'Reversibility',
		'GroEL_dependent_folding',
		'DnaK_dependent_folding',
		'N_terminal_methionine_cleavage',
		'RNA mods/enzyme',
		'tRNA-codon association',
		'RNA stability',
		'Complex Location',
		'Subunit Location',
		'Translocation Pathway',
		#'Translocation Multiplier'
	])

	def _get_reaction(x):
		if x is None:
			return []
		else:
			lst = []
			for gene in x.split(';'):
				if model.genes.has_id(gene):
					lst.append([ x.id for x in model.genes.get_by_id(gene).reactions ])

			# (An)A (old)locustag/gene name can be associated to many reactions
			lst = [ x for y in lst for x in y ]
			return lst

	def _get_reaction_name(x):
		if model.reactions.has_id(x):
			return model.reactions.get_by_id(x).name

	def _get_reversibility(x):
		if model.reactions.has_id(x):
			return 'True' if model.reactions.get_by_id(x).reversibility else 'False'

	def _get_spontaneous(x):
		if model.reactions.has_id(x):
			return 'True' if 'spontaneous' in model.reactions.get_by_id(x).name else 'False'

	#df['Gene Locus ID'] = [ x.qualifiers.get('locus_tag', [None])[0] for x in lst ]
	df['Gene Locus ID'] = [ x.qualifiers.get(locus_tag, [None])[0] for x in lst ]
	df['Definition'] = [ x.qualifiers.get('product', [None])[0] for x in lst ]
	df['Feature Type'] = [ x.type if x.qualifiers.get('pseudo') is None else 'pseudo' for x in lst ]

	tmp = [ x.qualifiers.get('gene', None) for x in lst ]
	df['Gene Names'] = [ ';'.join(x) if x is not None else None for x in tmp ]

	tmp = [ x.qualifiers['old_locus_tag'] if x.qualifiers.get('old_locus_tag', None) is not None else None for x in lst ]
	df['Old Locus Tag'] = [ ';'.join(x) if x is not None else None for x in tmp ]

	#df['M-model Reaction ID'] = df['Old Locus Tag'].apply(lambda x: _get_reaction(x))
	#df['M-model Reaction ID'] += df['Gene Names'].apply(lambda x: _get_reaction(x))
	#df['M-model Reaction ID'] += df['Gene Locus ID'].apply(lambda x: _get_reaction(x))
	#df = df.explode('M-model Reaction ID')

	#df['Reaction Name'] = df['M-model Reaction ID'].apply(lambda x: _get_reaction_name(x))
	#df['Reversibility'] = df['M-model Reaction ID'].apply(lambda x: _get_reversibility(x))

	# df.set_index(['Gene Locus ID', 'Definition', 'Feature type'], inplace = True)
	return df.sort_values(['M-model Reaction ID', 'Gene Locus ID'])

def complete_organism_specific_matrix(builder, data, model, output = False):
	if not hasattr(builder, 'org'):
		raise Exception('Please, run MEBuilder(*[configuration file]).generate_files() to generate the Organism-Specific Matrix.')

	# ME-model homology to reference
	def _get_homology(x, dct, keys):
		#tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = x[keys].to_list()
		tags = [ str(x).split(';') for x in tags ]

		lst = []
		for tag in [ x for y in tags for x in y ]:
			lst.append(dct.get(tag, None))
		lst = [ x for x in lst if x is not None ]

		if len(lst) != 0:
			return lst[0]

	def _get_definition(contigs, locus_tag):
		for contig in contigs:
			for feat in contig.features:
				if feat.type == 'CDS' and feat.qualifiers.get('locus_tag', [None])[0] == locus_tag:
					return ', '.join(feat.qualifiers.get('product', [None]))

	# Step 1a: Add reference data
	if hasattr(builder, 'homology'):
		dct = builder.homology.mutual_hits
		data['Reference BBH'] = data.apply(lambda x: _get_homology(x, dct, keys = ['Gene Locus ID', 'Old Locus Tag', 'BioCyc']), axis = 1)
		ref = builder.ref.contigs
		data['Reference Definition'] = data.apply(lambda x: _get_definition(ref, x['Reference BBH']), axis = 1)

	# Step 1b: Add BioCyc data if present
	if builder.configuration.get('biocyc.genes', False):
		# We reuse the _get_homology function, but changing the dictionary of relationships between IDs
		dct = { v['Accession-1']:k for k,v in builder.org.gene_dictionary.iterrows() }
		data['BioCyc'] = data.apply(lambda x: _get_homology(x, dct, keys = ['Gene Locus ID', 'Old Locus Tag', 'Gene Names']), axis = 1)
	else:
		data['BioCyc'] = None

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step1_get_homology'] = data.copy(deep = True)

	# Step 2: Add ME-model complexes
	# Restructure complexes_df to obtain the correct complex stoichiometry from the index
	df = builder.org.complexes_df.copy(deep = True)
	#df = df[~df.index.str.contains('MONOMER')]
	df = df[df['genes'].str.contains(r'\(')]
	df['genes'] = df['genes'].str.split(' AND ')
	df = df.explode('genes')
	df['stoich'] = df['genes'].apply(lambda x: '1' if x.split('(')[1][:-1] == '' else str(x.split('(')[1][:-1]))
	df['genes'] = df['genes'].apply(lambda x: x.split('(')[0].replace('-MONOMER', ''))
	df.index = df.index + ':' + df['stoich']

	def _get_complexes(x, df):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = [ str(x).split(';') for x in tags ]

		cplx = []
		for tag in [ x for y in tags for x in y ]:
			res = df[df['genes'].str.fullmatch(tag)].index.to_list()
			if len(res) != 0:
				cplx.append(res)
		if len(cplx) != 0:
			return cplx[0]

	data['Complex ID'] = data.apply(lambda x: _get_complexes(x, df), axis = 1)
	data = data.explode('Complex ID')

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step2_get_complexes'] = data.copy(deep = True)

	# Step 3: Add ME-model cofactors
	def _get_cofactors(x, dct):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], str(x['Complex ID']).split(':')[0] ]
		tags = [ str(x).split(';') for x in tags ]

		all_mods = []
		for tag in [ x for y in tags for x in y ]:
			all_mods.append([ x.split('_mod_')[1:] for x in dct.get(tag, [])]) # return a list of mods
		all_mods = [ x for y in all_mods for x in y ]
		# if multiple modifications: [['mn2(1)'], ['2fe2s(1)', 'mn2(1)'], ['4fe4s(1)', 'mn2(1)']]
		# into ['mn2(1)', '2fe2s(1) AND mn2(1), '4fe4s(1) AND mn2(1)]

		mod_strings = [None]
		if len(all_mods) != 0:
			for mods in all_mods:
				mods = [ x for x in mods if not x.startswith('SH') ] # Sulfur transfer in [enzyme]-S-sulfanylcysteine are enzymatic reactions
				logging.warning('The modification \'SH\' was removed. Add MetabolicReaction(s) to transfer sulfur from cysteine and to the many acceptors.')
				mods = [ x for x in mods if not x.startswith('Oxidized') ] # metabolic modification in ferredoxin and other proteins
				mods = [ x for x in mods if not x.startswith('palmitate') ] # metabolic modification from 2agpg160 in the lpp gene
				mods = [ x.replace('LI', 'li') for x in mods ]
				#mods = [ x.replace('lipo', 'lipoyl') for x in mods ]
				#logging.warning('The modification \'lipo\' was renamed to \'lipoyl\'.')
				#logging.warning('Add MetabolicReaction\'s to salvage or de novo synthesize lipoyl moieties. See https://www.genome.jp/pathway/map00785 for more information.')

				mods = [ x.replace('NiFeCoCN2', 'NiFe_cofactor') for x in mods ]
				logging.warning('The modification \'NiFeCoCN2\' was renamed to \'NiFe_cofactor\'.')
				logging.warning('Add MetabolicReactions to synthesize and transfer the NiFe_cofactor into the final acceptor enzyme. See https://biocyc.org/ECOLI/NEW-IMAGE?type=PATHWAY&object=PWY-8319 for more information.')

				if len(mods) != 0:
					if 'pan4p(1)' not in mods:
						mod_strings.append(' AND '.join(mods))
					else:
						# This remove all modification from holo-ACP protein. Modifications are metabolic.
						mod_strings.append('pan4p(1)')
						logging.warning('The modification \'pan4p\' was added. A MetabolicReaction will transfer pan4p from coenzyme A into Acyl-carrier protein(s).')

		if len(mod_strings) != 0:
			return mod_strings

	dct = {}
	for k,v in builder.org.protein_mod[['Core_enzyme']].to_dict()['Core_enzyme'].items():
		dct.setdefault(v, []).append(k)
	data['Cofactors in Modified Complex'] = data.apply(lambda x: _get_cofactors(x, dct), axis = 1)
	data = data.explode('Cofactors in Modified Complex')
	data = data.drop_duplicates()

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step3_get_cofactors'] = data.copy(deep = True)

	# Step 4: Add ME-model generics
	def _get_generics_from_gene(x, dct):
		generics = []
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for key, lst in dct.items():
			for tag in [ x for y in tags for x in y ]:
				if tag + '-MONOMER' == key[0]:
					generics.append(lst[0])
				if 'RNA_' + tag == key[0]:
					generics.append(lst[0])

		if len(generics) != 0:
			return [ x for y in generics for x in y ]

	def _get_generics_from_complex(x, dct):
		generics = []
		for key, lst in dct.items():
			if str(x['Complex ID']).split(':')[0] == key[0] and str(x['Cofactors in Modified Complex']).split(':')[0] == key[1]:
				generics.append(lst[0])
		if len(generics) != 0:
			return [ x for y in generics for x in y ]

	dct = pandas.DataFrame.from_dict({ k.replace('generic_', ''):v for k,v in builder.org.generic_dict.items() }).T
	dct = dct.explode('enzymes')
	dct['Complex'] = dct['enzymes'].apply(lambda gpr: [ x.split('_mod_')[0] for x in gpr.split(' OR ') ])
	dct['Cofactors'] = dct['enzymes'].apply(lambda gpr: [ ' AND '.join(x.split('_mod_')[1:]) for x in gpr.split(' OR ') ])
	dct['Generic'] = dct.index

	dct = dct.explode(['Complex', 'Cofactors']).replace('', 'None')
	dct = dct.groupby(['Complex', 'Cofactors']).agg({'Generic' : lambda x: set(x.tolist())})
	# final dictionary: (Complex, Cofactors) : Generic ID
	# we will repeat the Complex if it can form two or more generics
	#dct = { k:list(v['Generic'])[0] for k,v in dct.iterrows() }
	# final dictionary: (Complex, Cofactors) : [[Generic IDs]]
	tmp = {}
	for k,v in dct.iterrows():
		tmp.setdefault(k, []).append(list(v['Generic']))

	data['Generic Complex ID'] = data.apply(lambda x: _get_generics_from_gene(x, tmp), axis = 1)
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['Generic Complex ID'].update(data.apply(lambda x: _get_generics_from_complex(x, tmp), axis = 1))
	data.update({ 'Generic Complex ID' : data.apply(lambda x: _get_generics_from_complex(x, tmp), axis = 1) })
	data = data.explode('Generic Complex ID')

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step4_get_generics'] = data.copy(deep = True)

	# Step 5: Add M-model Reaction IDs, names, and reversibility from builder
	# It considers complex name and cofactors
	dct = builder.org.enz_rxn_assoc_df.copy(deep = True)
	dct['Complex'] = dct['Complexes'].apply(lambda gpr: [ x.split('_mod_')[0] for x in gpr.split(' OR ') ])
	dct['Cofactors'] = dct['Complexes'].apply(lambda gpr: [ ' AND '.join(x.split('_mod_')[1:]) for x in gpr.split(' OR ') ])
	dct['Reaction'] = dct.index

	dct = dct.explode(['Complex', 'Cofactors']).replace('', 'None')
	dct = dct.groupby(['Complex', 'Cofactors']).agg({'Reaction' : lambda x: x.tolist()})

	# intermediate dictionary: (Complex, Cofactors) : Reaction ID
	dct = { k:v['Reaction'] for k,v in dct.iterrows() }
	# correction based on generics -> NEW complex <-> list of reactions
	df = builder.org.complexes_df.copy(deep = True)
	df = { idx:[ x[8:-2] for x in row['genes'].split(' AND ')] for idx,row in df[df['genes'].str.contains('generic')].iterrows() }
	for cplx, generics in df.items():
		for generic in generics:
			dct.update({ (generic, 'None') : dct[(cplx, 'None')] })

	fn = lambda x: dct.get((str(x['Complex ID']).split(':')[0], str(x['Cofactors in Modified Complex'])), [None]) + dct.get((str(x['Generic Complex ID']).split(':')[0], str(x['Cofactors in Modified Complex'])), [None])

	data['M-model Reaction ID'] = data.apply(lambda x: fn(x), axis = 1)
	data = data.explode('M-model Reaction ID')
	data['Reaction Name'] = data['M-model Reaction ID'].apply(lambda x: model.reactions.get_by_id(x).name if model.reactions.has_id(x) else None)
	data['Reversibility'] = data['M-model Reaction ID'].apply(lambda x: str(model.reactions.get_by_id(x).reversibility) if model.reactions.has_id(x) else None)

	# correction based on generics -> NEW complex <-> list of reactions
	fn = lambda x: 'CPLX_{:s}-0:1({:s})'.format(x['M-model Reaction ID'], x['Generic Complex ID']) if x['M-model Reaction ID'] is not None and x['Generic Complex ID'] is not None else x['Complex ID']
	data['Complex ID'] = data.apply(lambda x: fn(x), axis = 1)
	fn = lambda x: None if x['M-model Reaction ID'] is not None else x['Generic Complex ID']
	data['Generic Complex ID'] = data.apply(lambda x: fn(x), axis = 1)
	data = data.drop_duplicates()

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step5_get_enzrxn_associations'] = data.copy(deep = True)

	# Step 6a: Add ME-model metacomplexes (e.g., ribosome)
	def _get_rnapol(x, lst):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			if '{:s}-MONOMER'.format(tag) in lst or tag.split(':')[0] in lst:
				return 'RNAP'

	lst = builder.org.RNAP if isinstance(builder.org.RNAP, list) else [builder.org.RNAP]
	lst += list(builder.org.sigmas.index)
	data['MetaComplex ID'] = data.apply(lambda x: _get_rnapol(x, lst), axis = 1)

	def _get_ribosome(x, dct):
		mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
		#tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
		tags = [ x['Complex ID'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			#if '{:s}-MONOMER'.format(tag) in lst or tag.split(':')[0] in lst or 'generic_{:s}'.format(tag) in lst:
			filter1 = '{:s}-MONOMER'.format(tag) in dct and mods in dct.get('{:s}-MONOMER'.format(tag), [None])
			filter2 = tag.split(':')[0] in dct and mods in dct.get(tag.split(':')[0], [None])
			filter3 = 'generic_{:s}'.format(tag) in dct #and mods in dct.get('generic_{:s}'.format(tag), [None])

			if (filter1 or filter2 or filter3) and mods == '':
				return 'ribosome:1'
			elif (filter1 or filter2 or filter3) and mods == 'acetyl(1)':
				return 'ribosome:2'

	dct = {}
	for x in [ x for y in [builder.org.ribosome_stoich[x]['stoich'].keys() for x in ['30_S_assembly', '50_S_assembly']] for x in y ]:
		dct.setdefault(x.split('_mod_')[0], []).append(' AND '.join(x.split('_mod_')[1:]))
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['MetaComplex ID'].update(data.apply(lambda x: _get_ribosome(x, dct), axis = 1))
	data.update({ 'MetaComplex ID' : data.apply(lambda x: _get_ribosome(x, dct), axis = 1) })

	def _get_degradosome_and_excision(x, dct):
		subrxns = []
		for key, subdct in dct.items():
			mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
			#tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
			tags = [ x['Complex ID'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				#if '{:s}-MONOMER'.format(tag) in lst or tag.split(':')[0] in lst or 'generic_{:s}'.format(tag) in lst:
				filter1 = '{:s}-MONOMER'.format(tag) in subdct and mods in subdct.get('{:s}-MONOMER'.format(tag), [None])
				filter2 = tag.split(':')[0] in subdct and mods in subdct.get(tag.split(':')[0], [None])
				filter3 = 'generic_{:s}'.format(tag) in subdct #and mods in subdct.get('generic_{:s}'.format(tag), [None])

				if filter1 or filter2 or filter3:
					subrxns.append(key + ':1')

		if len(subrxns) != 0:
			return subrxns

	dct = builder.org.excision_machinery
	dct['RNA_degradosome'] = builder.org.rna_degradosome['rna_degradosome']

	dct = { k:[ x for x in v['enzymes'] ] for k,v in dct.items() }
	for key in dct:
		tmp = {}
		for x in dct[key]:
			tmp.setdefault(x.split('_mod_')[0], []).append(' AND '.join(x.split('_mod_')[1:]))
		dct[key] = tmp

	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['MetaComplex ID'].update(data.apply(lambda x: _get_degradosome_and_excision(x, dct), axis = 1))
	data.update({ 'MetaComplex ID' : data.apply(lambda x: _get_degradosome_and_excision(x, dct), axis = 1) })

	# Step 6b: set RNA modifiers (CDS<->mod_at_position). Part1: Add 'RNA_modifier_enzyme' in 'MetaComplex ID' column
	def _get_rna_modifiers(x, dct):
		mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			filter1 = '{:s}-MONOMER'.format(tag) in dct and mods in dct.get('{:s}-MONOMER'.format(tag), [None])
			filter2 = tag.split(':')[0] in dct and mods in dct.get(tag.split(':')[0], [None])
			filter3 = 'generic_{:s}'.format(tag) in dct #and mods in dct.get('generic_{:s}'.format(tag), [None])

			if filter1 or filter2 or filter3:
				return 'RNA_modifier_enzyme'

	dct = { k.split('_mod_')[0]:' AND '.join(k.split('_mod_')[1:]) for k,v in builder.org.rna_modification.items() }
	dct = { k:v for k,v in dct.items() } # if v != '' } # just in case of empty values
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['MetaComplex ID'].update(data.apply(lambda x: _get_rna_modifiers(x, dct), axis = 1))
	data.update({ 'MetaComplex ID' : data.apply(lambda x: _get_rna_modifiers(x, dct), axis = 1) })

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step6_get_metacomplexes'] = data.copy(deep = True)

	# Step 7a: set RNA targets (tRNA<->mod_at_position)
	dct = builder.org.rna_modification_targets.copy(deep = True)
	dct['at'] = dct['modification'] + '_at_' + dct['position'].astype(str)
	dct = dct.groupby('bnum').agg({'at': lambda x: ','.join(x.tolist())}).to_dict()['at']
	data['RNA mods/enzyme'] = data['Gene Locus ID'].apply(lambda x: dct.get(str(x), None))

	# Step 7b: set RNA modifiers (CDS<->mod_at_position). Part2: Add 'modification' list in 'RNA mods/enzyme' column
	def _get_rna_modifications(x, dct):
		mod_at_pos = set()
		mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			filter1 = ('{:s}-MONOMER'.format(tag), mods) in dct
			filter2 = (tag.split(':')[0], mods) in dct
			filter3 = ('generic_{:s}'.format(tag), '') in dct

			if filter1 or filter2 or filter3:
				mod_at_pos.add(dct.get(('{:s}-MONOMER'.format(tag), mods), None))
				mod_at_pos.add(dct.get((tag.split(':')[0], mods), None))
				mod_at_pos.add(dct.get(('generic_{:s}'.format(tag), ''), None))

		lst = list(set([ x for x in mod_at_pos if x is not None ]))
		if len(lst) != 0:
			return ','.join(lst)

	dct = { (k.split('_mod_')[0], ' AND '.join(k.split('_mod_')[1:])):','.join(v) for k,v in builder.org.rna_modification.items() }
	dct = { k:v for k,v in dct.items() } # if v != '' } # just in case of empty values
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['RNA mods/enzyme'].update(data.apply(lambda x: _get_rna_modifications(x, dct), axis = 1))
	data.update({ 'RNA mods/enzyme' : data.apply(lambda x: _get_rna_modifications(x, dct), axis = 1) })

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step7_get_trna_mods_and_modifiers'] = data.copy(deep = True)

	# Step 8: Add translocation pathways
	def _get_transpaths(x, dct):
		pathways = []
		for key, lst in dct.items():
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or tag.split(':')[0] in lst or 'generic_{:s}'.format(tag) in lst:
					pathways.append('translocation_pathway_' + key)
		if len(pathways) != 0:
			return pathways

	dct = { k:v['enzymes'] for k,v in builder.org.translocation_pathways.items() if len(v['enzymes']) != 0 }
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['MetaComplex ID'].update(data.apply(lambda x: _get_transpaths(x, dct), axis = 1))
	data.update({ 'MetaComplex ID' : data.apply(lambda x: _get_transpaths(x, dct), axis = 1) })
	data = data.explode('MetaComplex ID')

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step8_get_translocation_pathways'] = data.copy(deep = True)

	# Step 9: Add ME-model subreactions
	# TODO: generalize to complexes with cofactors
	def _get_ribosome_subrxns(x, dct):
		for key, lst in dct.items():
			mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or tag.split(':')[0] in lst or 'generic_{:s}'.format(tag) in lst:
					return 'Ribosome_' + key if 'Ribosome_' not in key else key

	dct = { k:[ x.split('_mod_')[0] for x in [v['enzyme']] ] for k,v in builder.org.ribosome_subreactions.items() }
	data['ME-model SubReaction'] = data.apply(lambda x: _get_ribosome_subrxns(x, dct), axis = 1)

	def _get_translation_subrxns(x, dct):
		subprocess = []
		for key, subdct in dct.items():
			mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				filter1 = '{:s}-MONOMER'.format(tag) in subdct and mods in subdct.get('{:s}-MONOMER'.format(tag), [None])
				filter2 = tag.split(':')[0] in subdct and mods in subdct.get(tag.split(':')[0], [None])
				filter3 = 'generic_{:s}'.format(tag) in subdct #and mods in subdct.get('generic_{:s}'.format(tag), [None])

				if filter1 or filter2 or filter3:
					if key.endswith('InfA') or key.endswith('InfC'):
						subprocess.append(key)
					elif key.endswith('InfB'):
						subprocess.append('Translation_initiation_gtp_factor_InfB')
					elif key == 'fmet_addition_at_START':
						subprocess.append('Translation_initiation_' + key)
					elif key == 'FusA_mono_elongation':
						subprocess.append('Translation_elongation_FusA_mono')
					elif key == 'Tuf_gtp_regeneration':
						subprocess.append('Translation_elongation_' + key)
					elif key in ['N_terminal_methionine_cleavage', 'DnaK_dependent_folding', 'GroEL_dependent_folding']:
						subprocess.append('Protein_processing_' + key)
					elif key == 'PrfA_mono_mediated_termination':
						subprocess.append('Translation_termination_PrfA_mono_mediated')
					elif key == 'PrfB_mono_mediated_termination':
						subprocess.append('Translation_termination_PrfB_mono_mediated')
					elif key == 'generic_RF_mediated_termination':
						subprocess.append('Translation_termination_generic_RF_mediated')
					else:
						subprocess.append('Translation_termination_' + key)

		if len(subprocess) != 0:
			return subprocess

	dct = { k:[ x for x in v['enzymes'] ] for k,v in builder.org.initiation_subreactions.items() }
	dct.update({ k:[ x for x in v['enzymes'] ] for k,v in builder.org.elongation_subreactions.items() })
	dct.update({ k:[ x for x in v['enzymes'] ] for k,v in builder.org.termination_subreactions.items() })

	for key in dct:
		tmp = {}
		for x in dct[key]:
			tmp.setdefault(x.split('_mod_')[0], []).append(' AND '.join(x.split('_mod_')[1:]))
		dct[key] = tmp

	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['ME-model SubReaction'].update(data.apply(lambda x: _get_translation_subrxns(x, dct), axis = 1))
	data.update({ 'ME-model SubReaction' : data.apply(lambda x: _get_translation_subrxns(x, dct), axis = 1) })
	data = data.explode('ME-model SubReaction')

	def _get_transcription_subrxns(x, dct):
		subrxns = []
		for key, subdct in dct.items():
			mods = '' if x['Cofactors in Modified Complex'] is None else x['Cofactors in Modified Complex']
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Complex ID'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				filter1 = '{:s}-MONOMER'.format(tag) in subdct and mods in subdct.get('{:s}-MONOMER'.format(tag), [None])
				filter2 = tag.split(':')[0] in subdct and mods in subdct.get(tag.split(':')[0], [None])
				filter3 = 'generic_{:s}'.format(tag) in subdct #and mods in subdct.get('generic_{:s}'.format(tag), [None])

				if filter1 or filter2 or filter3:
					subrxns.append(key)
		if len(subrxns) != 0:
			return subrxns

	dct = { k:[ x for x in v['enzymes'] ] for k,v in builder.org.transcription_subreactions.items() }
	for key in dct:
		tmp = {}
		for x in dct[key]:
			tmp.setdefault(x.split('_mod_')[0], []).append(' AND '.join(x.split('_mod_')[1:]))
		dct[key] = tmp

	# Copy-on-Write pandas 3.0 behavior is default to be True
	# data['ME-model SubReaction'].update(data.apply(lambda x: _get_transcription_subrxns(x, dct), axis = 1))
	data.update({ 'ME-model SubReaction' : data.apply(lambda x: _get_transcription_subrxns(x, dct), axis = 1) })
	data = data.explode('ME-model SubReaction')

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step9_get_subreactions'] = data.copy(deep = True)

	# Step 10: Add post-translation processing
	def _get_processing_targets(x, lst):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			if tag in lst:
				return 'True'

	lst = builder.org.folding_dict['GroEL_dependent_folding']['enzymes']
	data['GroEL_dependent_folding'] = data.apply(lambda x: _get_processing_targets(x, lst), axis = 1)
	lst = builder.org.folding_dict['DnaK_dependent_folding']['enzymes']
	data['DnaK_dependent_folding'] = data.apply(lambda x: _get_processing_targets(x, lst), axis = 1)
	lst = builder.org.cleaved_methionine
	data['N_terminal_methionine_cleavage'] = data.apply(lambda x: _get_processing_targets(x, lst), axis = 1)
# 	lst = builder.org.stable_RNAs
# 	data['RNA stability'] = data.apply(lambda x: _get_processing_targets(x, lst), axis = 1)

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step10_get_posttranslation_processing'] = data.copy(deep = True)

	# Step 11: Add tRNA to codon association from GenBank data
	dct = { k:','.join(v) for x in [ v for k,v in builder.me_model.global_info['trna_to_codon'].items() ] for k,v in x.items() }
	data['tRNA-codon association'] = data.apply(lambda x: dct.get(x['Gene Locus ID'], None), axis = 1)

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step11_get_trna2codon_specificity'] = data.copy(deep = True)

	# Step 12: protein location
	def _get_protein_location(x, df):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ] # we convert here None to 'None'
		tags = [ x for y in tags for x in y ]
		# extra modification if compared this function to the others above
		tags = [ r'{:s}(\(\d*\)|\(\d*\:\d*\))'.format(x) for x in tags if x != 'None' ]

		#res = df[df['Protein'].str.match(x)][['Complex_compartment', 'Protein_compartment', 'translocase_pathway']].values
		#query = [ '{:s}()'.format(y) for y in tags ] # sp old_locus_tag's, the query can be a string of tags separated by ';'
		res = df[df['Protein'].str.fullmatch('|'.join(tags))][['Complex_compartment', 'Protein_compartment', 'translocase_pathway']].values

		if len(res) != 0:
			return res[0]
		else:
			return [None, None, None]

	df = builder.org.protein_location.dropna(how = 'all', subset = ['Complex_compartment', 'Protein', 'Protein_compartment'])
	data['Complex Location'], data['Subunit Location'], data['Translocation Pathway'] = zip(*data.apply(lambda x: _get_protein_location(x, df), axis = 1))
	data = data.explode('Complex Location')

	# debug purposes
	if builder.configuration.get('collect_debug_data', False):
		builder.debug_data['complete_organism_specific_matrix_step12_get_protein_location'] = data.copy(deep = True)

	# final sorting
	data = data.sort_values(['M-model Reaction ID', 'Gene Locus ID'])
	data = data.drop_duplicates()

	if output:
		try:
			# Save file as excel or tsv depending on the size
			if not str(output).endswith('.xlsx'):
				output = '.'.join(output.split('.')[:-1]) + '.xlsx'

			# GPRs can expand the model specifications beyond the max size of an excel file (1048576 rows)
			if data.shape[0] > 1048575: # one less to accommodate the header
				# Divide the DataFrame into pieces and save them
				rxn_ids = data.groupby(['M-model Reaction ID'], dropna = True).size() # Anything with 1 or more occurrences

				gene_ids = data[data['M-model Reaction ID'].isin(rxn_ids[rxn_ids < 1000].index)]['Gene Locus ID']
				tmp = pandas.concat([ data[data['Gene Locus ID'].isin(gene_ids)], data[data['M-model Reaction ID'].isna()]], axis = 1)
				_save_to_excel(tmp, '.'.join(output.split('.')[:-1]) + '_{:02d}.xlsx'.format(0))

				if output.endswith('.xlsx'):
					for idx, gene_id in enumerate(rxn_ids[rxn_ids >= 1000].index):
						tmp = data[data['Gene Locus ID'].isin([gene_id])]
						_save_to_excel(tmp, '.'.join(output.split('.')[:-1]) + '_{:03d}.xlsx'.format(idx+1))

				# save the output as a tsv file
				if not output.endswith('.txt'):
					output = '.'.join(output.split('.')[:-2]) + '.txt'

				with open(output, 'w') as outfile:
					data.to_csv(outfile, index = False, sep = '\t')
			else:
				_save_to_excel(data, output)
		except:
			logging.warning('The builder.df_data was not saved to the \'{:s}\' file.'.format(output))

	return data

def _correct_input(df):
	# correct Gene Locus ID to reflect if they are proteins or RNAs
	fn = lambda x: \
		'protein_{:s}'.format(x['Gene Locus ID']) if (x['Feature Type'] == 'CDS' and not x['Gene Locus ID'].startswith('protein_')) \
		else 'RNA_{:s}'.format(x['Gene Locus ID']) if (x['Feature Type'] != 'pseudo' and not x['Gene Locus ID'].startswith('RNA_')) \
		else 'pseudo_{:s}'.format(x['Gene Locus ID'])
	df.loc[:, 'Gene Locus ID'] = df[['Gene Locus ID', 'Feature Type']].apply(fn, axis = 1)

	# correct generics to contain the `generic_` prefix
	fn = lambda x: 'generic_{:s}'.format(x) if isinstance(x, str) and not x.startswith('generic_') else x
	df.loc[:, 'Generic Complex ID'] = df['Generic Complex ID'].apply(fn)
	return df

# 2nd step
def get_generics(df):
	#tmp = df[df['Generic Complex ID'].notna() & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['Generic Complex ID'].notna()].copy() # this was a view
	tmp = _correct_input(tmp)

	fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp.loc[:, 'Modified Complex'] = tmp.apply(fn, axis = 1).values

	# collapse
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
	# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

	tmp.loc[:, 'Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)
	tmp = tmp.groupby(['Generic Complex ID']).agg({'Gene Locus ID': lambda x: sorted(set(x.tolist()))}).to_dict()['Gene Locus ID'].items()
	return tmp

def _get_metacomplex_stoichiometry(df, key):
	#tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.startswith(key) & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.startswith(key)].copy() # this was a view
	if not tmp.empty:
		tmp = _correct_input(tmp)

		tmp.loc[:, 'Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

		fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
			if isinstance(x['Cofactors in Modified Complex'], str) else None #numpy.nan
		tmp.loc[:, 'Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		# collapse
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
		tmp.update({ 'Modified Complex' : tmp['Generic Complex ID'] })
		# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
		tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
		# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
		tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

		tmp.loc[:, 'stoich'] = tmp['MetaComplex ID'].str.split(':', expand = True).iloc[:, 1]
		return tmp

def ribosome_stoichiometry(df):
	tmp = _get_metacomplex_stoichiometry(df, 'ribosome')
	ribosome_stoich = { 'ribosome' : { 'stoich' : { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich'])}}}
	ribosome_stoich['ribosome']['stoich']['gtp_c'] = 1
	return ribosome_stoich

def degradosome_stoichiometry(df):
	tmp = _get_metacomplex_stoichiometry(df, 'RNA_degradosome')
	return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) } if tmp is not None else None

def dnapolymerase_stoichiometry(df):
	tmp = _get_metacomplex_stoichiometry(df, 'DNAP')
	return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) }

def excision_machinery_stoichiometry(df, keys):
	#tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.match(keys) & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.match(keys)].copy() # this was a view
	if not tmp.empty:
		tmp = _correct_input(tmp)

		fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
			if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
		tmp.loc[:, 'Modified Complex'] = tmp.apply(fn, axis = 1)

		# collapse
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
		tmp.update({ 'Modified Complex' : tmp['Generic Complex ID'] })
		# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
		tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
		# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
		tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

		tmp.loc[:, 'stoich'] = tmp['MetaComplex ID'].str.split(':', expand = True).iloc[:, 1]
		return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) }

def aa_synthetase_dict(df):
	#tmp = df[df['Definition'].str.contains('--tRNA ligase|-tRNA synthetase') & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['Definition'].str.contains('--tRNA ligase|-tRNA synthetase') & df['Definition'].notna()]
	tmp = _correct_input(tmp)

	fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp.apply(fn, axis = 1)

	# collapse
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
	tmp.update({ 'Modified Complex' : tmp['Generic Complex ID'] })
	# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
	# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

	tmp = tmp.drop_duplicates(['Gene Locus ID'])

	tmp = { k.split('--tRNA')[0].lower():v.split(':')[0] for k,v in zip(tmp['Definition'], tmp['Gene Locus ID']) }
	tmp = { k.replace('asparagine', 'asn'):v for k,v in tmp.items() }
	tmp = { k.replace('isoleucine', 'ile'):v for k,v in tmp.items() }
	tmp = { k.replace('tryptophan', 'trp'):v for k,v in tmp.items() }
	return { (k[:3]+'__L_c').replace('gly__L_c', 'gly_c'):v for k,v in tmp.items() }

def get_subreactions(df, key: str):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = _correct_input(tmp)
	tmp = tmp[tmp['ME-model SubReaction'].notna()]

	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}_cplx'.format(x))
	tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

	fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
	tmp.update({ 'Modified Complex' : tmp['Generic Complex ID'] })
	# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
	# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

	# Correct Translation_termination_generic_RF_mediated
	if key.startswith('Translation_termination_generic_RF_mediated'):
		cplxs = tmp[tmp['ME-model SubReaction'].str.match('Translation_termination_generic_RF_mediated')]['Complex ID'].tolist()
		cplxs = set([ x.split(':')[0] for x in cplxs ])

	tmp = tmp.groupby(['ME-model SubReaction']).agg({'Gene Locus ID': lambda x: x.tolist()})
	if key.startswith('Translation_termination_generic_RF_mediated'):
		tmp = tmp.apply(lambda x: [['generic_RF']] if set(x['Gene Locus ID']) == cplxs else x, axis = 1)

	return { k : {'enzymes' : list(set(v)), 'stoich': {}, 'element_contribution' : {}, 'keff' : []}
		for k,v in zip(tmp.index, tmp['Gene Locus ID']) if k.startswith(key) }

def _get_df_rxns(df):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['M-model Reaction ID'].notna()]
	tmp = tmp[['M-model Reaction ID', 'Reaction Name', 'Reversibility']]
	tmp.columns = ['name', 'description', 'is_reversible']
	tmp = tmp.dropna().drop_duplicates('name', keep = 'first').set_index('name')

	return tmp

def _get_df_cplxs(df, generics = False):
	tmp = df[df['Feature Type'].isin(['CDS', 'pseudo', 'ncRNA'])].fillna('')

	# get enzymatic complexes with generics subunits
	df_generics = df[df['Feature Type'].isin(['CDS'])].fillna('')
	df_generics = df_generics[df_generics['Complex ID'].str.contains(r'\(')]
	df_generics['Gene Locus ID'] = 'generic_' + df_generics['Complex ID'].str.extract(r'[A-Za-z0-9\_-]+:\d+\(([A-Za-z0-9\_-]+)\)', expand = True)
	df_generics['Complex ID'] = df_generics['Complex ID'].str.extract(r'([A-Za-z0-9\_-]+:\d+)\([A-Za-z0-9\_-]+\)', expand = True)

	tmp = pandas.concat([tmp, df_generics])
	tmp = tmp[~tmp['Complex ID'].str.contains(r'\(')]

	fn = lambda x: '{:s} complex'.format(x['Definition']) if x['Complex Name'] == '' else x['Complex Name']
	tmp['Complex Name'] = tmp[['Definition', 'Complex Name']].apply(fn, axis = 1)

	fn = lambda x: '{:s}({:s})'.format(x['Gene Locus ID'], x['Complex ID'].split(':')[1].split('(')[0]) if ':' in x['Complex ID'] else '{:s}(1)'.format(x['Gene Locus ID'])
	tmp['Gene Locus ID'] = tmp[['Gene Locus ID', 'Complex ID']].apply(fn, axis = 1)

	fn = lambda x: '{:s}-MONOMER'.format(x['Gene Locus ID'].split('(')[0]) if x['Complex ID'] == '' else x['Complex ID'].split('(')[0].split(':')[0]
	tmp['Complex ID'] = tmp[['Gene Locus ID', 'Complex ID']].apply(fn, axis = 1)

	fn = lambda x: x[0] if isinstance(x, list) else x.tolist()[0]
	tmp = tmp.groupby(['Complex ID']).agg({'Gene Locus ID': lambda x: ' AND '.join(sorted(set(x.tolist()))), 'Complex Name': fn}).reset_index()
	tmp = tmp[['Complex ID', 'Complex Name', 'Gene Locus ID']].set_index('Complex ID')

	return tmp

def _get_df_ptms(df):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['Cofactors in Modified Complex'].notna()]

	if not tmp.empty:
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split('(')[0].split(':')[0] if isinstance(x, str) else x)

		fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')])
		tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		# Modified complex names can be replaced using an alias
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Modified Complex'].update(tmp['Complex Name'].str.replace(' ', '_')) # inplace
		tmp.update({ 'Modified Complex' : tmp['Complex Name'].str.replace(' ', '_') })

		tmp = tmp[['Modified Complex', 'Complex ID', 'Cofactors in Modified Complex']]
		tmp.columns = ['Modified_enzyme', 'Core_enzyme', 'Modifications']
		tmp = tmp.drop_duplicates(keep = 'first').set_index('Modified_enzyme')

		return tmp
	else:
		return pandas.DataFrame(columns = ['Modified Complex', 'Complex ID', 'Cofactors in Modified Complex'])

def _get_df_enz2rxn(df, filter_in = set(), generics = False):
	tmp = df.copy(deep = True)
	#tmp = _correct_input(tmp)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['M-model Reaction ID'].notna()]

	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}-MONOMER'.format(x))
	tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

	fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	# Copy-on-Write pandas 3.0 behavior is default to be True
	#tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
	tmp.update({ 'Modified Complex' : tmp['Generic Complex ID'] })
	# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
	# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

	# Modified complex names can be replaced using an alias
	# Copy-on-Write pandas 3.0 behavior is default to be True
	# tmp['Gene Locus ID'].update(tmp['Complex Name'].str.replace(' ', '_')) # inplace
	tmp.update({ 'Gene Locus ID' : tmp['Complex Name'].str.replace(' ', '_') })

	tmp = tmp.groupby(['M-model Reaction ID']).agg({'Gene Locus ID': lambda x: ' OR '.join(sorted(set(x.tolist())))})
	tmp.columns = ['Complexes']

	return tmp

def _get_df_rna_enzs(df, filter_in = set(), generics = False):
	tmp = df[~df['Feature Type'].isin(['tRNA', 'rRNA'])]
	tmp = tmp[tmp['RNA mods/enzyme'].notna()]

	if not tmp.empty:
		tmp = _correct_input(tmp)
		tmp['RNA mods/enzyme'] = tmp['RNA mods/enzyme'].str.split(',')
		tmp = tmp.explode('RNA mods/enzyme')

		tmp['modification'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[0])
		tmp['positions'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[1])

		tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}-MONOMER'.format(x))
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

		fn = lambda x: x['Complex ID'].split(':')[0] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
			if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
		tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		# collapse
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
		tmp.update({ 'Modified Complex' : tmp['Generic Complex ID'] })
		# tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
		tmp.update({ 'Complex ID' : tmp['Modified Complex'] })
		# tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
		tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

		tmp = tmp[['Gene Locus ID', 'modification', 'positions']]
		tmp.columns = ['enzymes', 'modification', 'positions']
		tmp = tmp.drop_duplicates(keep = 'first')

		return tmp.groupby(['modification', 'positions']).agg({'enzymes': lambda x: ' AND '.join(x.to_list())}).reset_index()
	else:
		return pandas.DataFrame(columns = ['enzymes', 'modification', 'positions'])

def _get_df_rna_ptms(df, filter_in = set(), generics = False):
	tmp = df[df['Feature Type'].isin(['tRNA', 'rRNA'])]
	tmp = tmp[tmp['RNA mods/enzyme'].notna()]

	if not tmp.empty:
		tmp['RNA mods/enzyme'] = tmp['RNA mods/enzyme'].str.split(',')
		tmp = tmp.explode('RNA mods/enzyme')

		tmp['modification'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[0])
		tmp['positions'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[1])

		# collapse
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Gene Locus ID'].update(tmp['Generic Complex ID']) # inplace
		tmp.update({ 'Gene Locus ID' : tmp['Generic Complex ID'] })

		tmp = tmp[['Gene Locus ID', 'modification', 'positions']]
		tmp.columns = ['bnum', 'modification', 'positions']
		tmp = tmp.drop_duplicates(keep = 'first')

		return tmp
	else:
		return pandas.DataFrame(columns = ['bnum', 'modification', 'positions'])

def _get_df_protloc(df, filter_in = set(), generics = False):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['Complex Location'].notna()]

	if not tmp.empty:
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)
		tmp['Monomer ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}-MONOMER'.format(x))

		# collapse
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Monomer ID'].update(tmp['Complex ID']) # inplace
		tmp.update({ 'Monomer ID' : tmp['Complex ID'] })

		tmp = tmp[['Monomer ID', 'Complex Location', 'Gene Locus ID', 'Subunit Location', 'Translocation Pathway']]
		tmp.columns = ['Complex', 'Complex_compartment', 'Protein', 'Protein_compartment', 'translocase_pathway']
		tmp = tmp.drop_duplicates(subset = ['Complex', 'Protein', 'Protein_compartment'], keep = 'first').set_index('Complex')

		return tmp
	else:
		#return pandas.DataFrame(columns = ['Complex ID', 'Complex Location', 'Gene Locus ID', 'Subunit Location', 'Translocation Pathway'])
		return pandas.DataFrame(columns = ['Complex', 'Complex_compartment', 'Protein', 'Protein_compartment', 'translocase_pathway'])

def _get_df_transpaths(df, filter_in = set(), generics = False):
	tmp = df[~df['MetaComplex ID'].isnull()]
	tmp = tmp[tmp['MetaComplex ID'].str.match('translocation_pathway')]

	if not tmp.empty:
		tmp = _correct_input(tmp)
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

		# collapse
		# Copy-on-Write pandas 3.0 behavior is default to be True
		# tmp['Complex ID'].update(tmp['Generic Complex ID'])
		tmp.update({ 'Complex ID' : tmp['Generic Complex ID'] })
		# tmp['Gene Locus ID'].update(tmp['Complex ID'])
		tmp.update({ 'Gene Locus ID' : tmp['Complex ID'] })

		return tmp[['Gene Locus ID', 'MetaComplex ID']].groupby('MetaComplex ID').agg(lambda x: list(set(x.tolist())))
	else:
		return pandas.DataFrame(columns = ['Gene Locus ID'], index = ['MetaComplex ID'])

def get_df_input_from_excel(df, df_rxns):
	lst = [
		#'Complex Name',
		#'Complex ID',
		'Cofactors in Modified Complex',
		'Generic Complex ID',
		'MetaComplex ID',
		'ME-model SubReaction',
		'M-model Reaction ID',
		'RNA mods/enzyme'
		]

	tmp1 = df[df['Feature Type'].str.match('CDS')]#.dropna(subset = lst, how = 'all', axis = 0)
	tmp2 = df[df['Feature Type'].isin(['rRNA', 'tRNA', 'ncRNA', 'tmRNA', 'misc_RNA', 'pseudo'])]
	df = pandas.concat([tmp1, tmp2], axis = 0).drop_duplicates()
	df = df.fillna({
		'Gene Locus ID' : '',
		'Reversibility' : numpy.nan,
		#'Spontaneous' : 'False', # see reactions.txt input file
		'GroEL_dependent_folding' : 'False',
		'DnaK_dependent_folding' : 'False',
		'N_terminal_methionine_cleavage' : 'False',
		'RNA stability' : 'False'
		})

	tmp = _get_df_rxns(df)
	for idx, x in df_rxns[df_rxns.index.isin(tmp.index)].iterrows():
		logging.warning('The reaction \'{:s}\' appears in the M-model and in the \'df_metadata_orphan_rxns\' input (default value)'.format(idx))
		logging.warning('If you want to use the M-model \'{:s}\' reaction, delete the ID from \'df_metadata_orphan_rxns\'. Otherwise, add the ID to \'defer_to_rxn_matrix\'.'.format(idx))

	# remove entries that are in df_rxns (user input overrides m-model info)
	tmp = tmp[~tmp.index.isin(df_rxns.index)]
	df_rxns = pandas.concat([tmp, df_rxns]).fillna('False')
	df_cplxs = _get_df_cplxs(df)
	df_ptms = _get_df_ptms(df)
	df_enz2rxn = _get_df_enz2rxn(df)
	df_protloc = _get_df_protloc(df)
	df_transpaths = _get_df_transpaths(df)

	# modifications of tRNAs and rRNAs
	df_rna_enzs = _get_df_rna_enzs(df)
	df_rna_ptms = _get_df_rna_ptms(df)
	cols = ['modification', 'positions']

	if not df_rna_ptms.empty or not df_rna_enzs.empty:
		df_rna_mods = pandas.merge(df_rna_ptms, df_rna_enzs, how = 'left', left_on = cols, right_on = cols).fillna('No_Machine')
	else:
		df_rna_mods = pandas.DataFrame(columns = ['bnum', 'modification', 'positions', 'enzymes'])

	return df, df_rxns, df_cplxs, df_ptms, df_enz2rxn, df_rna_mods, df_protloc, df_transpaths

def filter_out_cofactor(df, cofactors_to_skip):
	df = df.copy(deep = True)
	def filt(old_cofactors, *lst):
		if isinstance(old_cofactors, str):
			new_string = ' AND '.join([ cofactor for cofactor in old_cofactors.split(' AND ') if cofactor.split('(')[0] not in lst ])
			return new_string if new_string != '' else None
		else:
			return old_cofactors

	df['Cofactors in Modified Complex'] = df['Cofactors in Modified Complex'].apply(filt, args = (cofactors_to_skip))
	return df

def filter_in_cofactor(df, cofactors_to_keep):
	df = df.copy(deep = True)
	def filt(old_cofactors, *lst):
		if isinstance(old_cofactors, str):
			new_string = ' AND '.join([ cofactor for cofactor in old_cofactors.split(' AND ') if cofactor.split('(')[0] in lst ])
			return new_string if new_string != '' else None
		else:
			return old_cofactors

	df['Cofactors in Modified Complex'] = df['Cofactors in Modified Complex'].apply(filt, args = (cofactors_to_keep))
	return df
