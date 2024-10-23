import re
import sys
import numpy
import cobra
import pandas
from IPython.display import display, HTML

def FromExcel(infile:str, model_name:str, outfile:str, f_replace:dict = {}, debug:bool = False) -> cobra.core.model.Model:
	code = []
	code.append('import cobra\n')
	#code.append('try:\n')
	#code.append('\tdel model\n')
	#code.append('except:\n')
	#code.append('\tpass\n')
	code.append('model = cobra.Model(\'{:s}\')\n'.format(model_name))

	# metabolites
	mets = pandas.read_excel(infile, sheet_name = 'metabolites').dropna(axis = 0, how = 'all').fillna('').reset_index()

	# apply replace
	if f_replace:
		mets['_id'] = mets['_id'].replace(f_replace)

	base = '_{:s} = cobra.Metabolite(\'{:s}\', formula = \'{:s}\', name = \'{:s}\', compartment = \'{:s}\', charge = {:d})\n' \
		'_{:s}.annotation = {:s}\n' \
		'model.add_metabolites([_{:s}])\n'

	for idx in mets.index:
		#if 'model' in mets.columns:
		if mets.iloc[idx]['model'] == '__REMOVE__' or mets.iloc[idx]['model'] == '__TEST__':
			continue
		else:
			annots = {}
			if 'sbo' in mets.columns:
				annots['sbo'] = mets.iloc[idx]['sbo']
			if 'pubchem.compound' in mets.columns and isinstance(mets.iloc[idx]['pubchem.compound'], float):
				annots['pubchem.compound'] = ['{:.0f}'.format(mets.iloc[idx]['pubchem.compound'])]
			if 'kegg.compound' in mets.columns:
				annots['kegg.compound'] = mets.iloc[idx]['kegg.compound'].split(';') if (mets.iloc[idx]['kegg.compound'] != '') else []
			if 'seed.compound' in mets.columns:
				annots['seed.compound'] = mets.iloc[idx]['seed.compound'].split(';') if (mets.iloc[idx]['seed.compound'] != '') else []
			if 'inchikey' in mets.columns:
				annots['inchikey'] = mets.iloc[idx]['inchikey'].split(';') if (mets.iloc[idx]['inchikey'] != '') else []
			if 'inchi' in mets.columns:
				annots['inchi'] = mets.iloc[idx]['inchi'].split(';') if (mets.iloc[idx]['inchi'] != '') else []
			if 'bigg.metabolite' in mets.columns:
				annots['bigg.metabolite'] = mets.iloc[idx]['bigg.metabolite']
			if 'SMILES' in mets.columns:
				annots['SMILES'] = mets.iloc[idx]['SMILES'].split(';') if (mets.iloc[idx]['SMILES'] != '') else []
			if 'chebi' in mets.columns:
				annots['chebi'] = mets.iloc[idx]['chebi'].split(';') if (mets.iloc[idx]['chebi'] != '') else []
			if 'biocyc' in mets.columns:
				annots['biocyc'] = mets.iloc[idx]['biocyc'].split(';') if (mets.iloc[idx]['biocyc'] != '') else []
			if 'hmdb' in mets.columns:
				annots['hmdb'] = mets.iloc[idx]['hmdb'].split(';') if (mets.iloc[idx]['hmdb'] != '') else []
			if 'lipidmaps' in mets.columns:
				annots['lipidmaps'] = mets.iloc[idx]['lipidmaps'].split(';') if (mets.iloc[idx]['lipidmaps'] != '') else []
			if 'metanetx.chemical' in mets.columns:
				annots['metanetx.chemical'] = mets.iloc[idx]['metanetx.chemical'].split(';') if (mets.iloc[idx]['metanetx.chemical'] != '') else []
			if 'seed.compound' in mets.columns:
				annots['seed.compound'] = mets.iloc[idx]['seed.compound'].split(';') if (mets.iloc[idx]['seed.compound'] != '') else []
			if 'reactome' in mets.columns:
				annots['reactome'] = mets.iloc[idx]['reactome.compound'].split(';') if (mets.iloc[idx]['reactome.compound'] != '') else []
			annots = str(annots).replace('{', '{\n\t').replace('}', '\n\t}').replace(', ', ',\n\t') if annots else '{}'

			tmp = base.format(
				mets.iloc[idx]['_id'].replace('-', '_DASH_'),
				mets.iloc[idx]['_id'],
				mets.iloc[idx]['formula'],
				mets.iloc[idx]['name'].replace('\'', '\\\''),
				mets.iloc[idx]['compartment'],
				int(mets.iloc[idx]['charge']),
				mets.iloc[idx]['_id'].replace('-', '_DASH_'),
				annots,
				mets.iloc[idx]['_id'].replace('-', '_DASH_'))

			code.append(tmp)

	# reactions
	rxns = pandas.read_excel(infile, sheet_name = 'reactions').dropna(axis = 0, how = 'all').fillna('').reset_index()

	base = '\n' \
		'reaction = cobra.Reaction(\'{:s}\')\n' \
		'reaction.name = \'{:s}\'\n' \
		'reaction.subsystem = \'{:s}\'\n' \
		'reaction.lower_bound = {:f}\n' \
		'reaction.upper_bound = {:f}\n' \
		'reaction.add_metabolites({:s})\n' \
		'reaction.annotation = {:s}\n' \
		'reaction.gene_reaction_rule = \'{:s}\'\n' \
		'model.add_reactions([reaction])\n'

	for idx in rxns.index:
		#if 'model' in rxns.columns:
		if rxns.iloc[idx]['model'] == '__REMOVE__' or rxns.iloc[idx]['model'] == '__TEST__':
			continue
		else:
			# reconstruct reaction from string
			regex = re.compile(r'([A-Za-z0-9\.\_\-]+)')
			reaction = rxns.iloc[idx]['_metabolites']

			# apply replace
			for key, value in f_replace.items():
				reaction = reaction.replace(key, value)

			try:
				subs = re.findall(regex, reaction.split('=')[0].strip())
			except Exception as e:
				if rxns.iloc[idx]['model'] != '__BIOMASS__' and debug:
					print('Error in:')
					display(rxns.iloc[idx].to_frame().T)
				#traceback.print_exc()

			try:
				prods = re.findall(regex, reaction.split('=')[1].strip())
			except Exception as e:
				if rxns.iloc[idx]['model'] != '__BIOMASS__' and debug:
					print('Error in:')
					display(rxns.iloc[idx].to_frame().T)
				#traceback.print_exc()

			# correct strings
			#new_subs = []
			#for idj, x in enumerate(subs):
				#if x.replace('.', '').isnumeric():
					#new_subs.append(x)
					#subs.pop(idj)
				#else:
					#new_subs.append('1')

			#new_prods = []
			#for idj, x in enumerate(prods):
				#if x.replace('.', '').isnumeric():
					#new_prods.append(x)
					#prods.pop(idj)
				#else:
					#new_prods.append('1')

			#subs = [ x for y in list(zip(new_subs, subs)) for x in y ]
			#prods = [ x for y in list(zip(new_prods, prods)) for x in y ]

			invert = +1
			if 'model' in rxns.columns:
				if rxns.iloc[idx]['model'] == '__INVERT__':
					invert = -1
			lower = rxns.iloc[idx]['_lower_bound'] if invert == +1 else -1*rxns.iloc[idx]['_upper_bound']
			upper = rxns.iloc[idx]['_upper_bound'] if invert == +1 else -1*rxns.iloc[idx]['_lower_bound']

			mets = {}
			for coeff, substrate in zip(subs[0::2], subs[1::2]):
				substrate = substrate.replace('-', '_DASH_')
				mets['_' + substrate] = -1*float(coeff)*invert

			for coeff, product in zip(prods[0::2], prods[1::2]):
				product = product.replace('-', '_DASH_')
				try:
					mets['_' + product] = mets['_' + product] + 1*float(coeff)*invert
				except:
					mets['_' + product] = +1*float(coeff)*invert
			mets = str(mets).replace('{', '{\n\t').replace('}', '\n\t}').replace(', ', ',\n\t').replace('\'', '')

			annots = {}
			if 'sbo' in rxns.columns:
				annots['sbo'] = rxns.iloc[idx]['sbo']
			if 'bigg.reaction' in rxns.columns:
				annots['bigg.reaction'] = rxns.iloc[idx]['bigg.reaction'] if (rxns.iloc[idx]['bigg.reaction'] != '') else []
			if 'biocyc' in rxns.columns:
				annots['biocyc'] = rxns.iloc[idx]['biocyc'].replace('META:', '').split(';') if (rxns.iloc[idx]['biocyc'] != '') else []
			if 'ec-code' in rxns.columns:
				annots['ec-code'] = rxns.iloc[idx]['ec-code'].split(';') if (rxns.iloc[idx]['ec-code'] != '') else []
			if 'kegg.reaction' in rxns.columns:
				annots['kegg.reaction'] = rxns.iloc[idx]['kegg.reaction'].split(';') if (rxns.iloc[idx]['kegg.reaction'] != '') else []
			annots = str(annots).replace('{', '{\n\t').replace('}', '\n\t}').replace(', ', ',\n\t') if annots else '{}'

			tmp = base.format(
				rxns.iloc[idx]['_id'],
				rxns.iloc[idx]['name'].replace('\'', '\\\''),
				rxns.iloc[idx]['subsystem'],
				lower, upper, mets, annots,
				rxns.iloc[idx]['_gpr'])

			code.append(tmp)

		if 'model' in rxns.columns:
			if rxns.iloc[idx]['model'] == '__OBJ__':
				code.append('\nmodel.objective = \'{:s}\'\n'.format(rxns.iloc[idx]['_id']))

	# gene annotations
	genes = pandas.read_excel(infile, sheet_name = 'genes', dtype = str).dropna(axis = 0, how = 'all').fillna('').reset_index()

	base = 'try:\n' \
		'\tmodel.genes.get_by_id(\'{:s}\').annotation = {:s}\n' \
		'except:\n' \
		'\tprint(\'INFO: gene ID {:s} not associated to any reaction in the M-model.\')\n\n'

	code.append('\n')
	for idx in genes.index:
		annots = {}
		for annot in genes.columns[2:]:
			if isinstance(genes.iloc[idx][annot], float):
				annots[annot] = '{:.0f}'.format(genes.iloc[idx][annot]) if (genes.iloc[idx][annot] != '') else []
			elif ';' in genes.iloc[idx][annot]:
				annots[annot] = genes.iloc[idx][annot].split(';') if (genes.iloc[idx][annot] != '') else []
			else:
				annots[annot] = [genes.iloc[idx][annot]] if (genes.iloc[idx][annot] != '') else []

		annots = str(annots).replace('{', '{\n\t\t').replace('}', '\n\t\t}').replace(', ', ',\n\t\t') if annots else '{}'
		code.append(base.format(genes.iloc[idx]['_id'], annots, genes.iloc[idx]['_id']))

	code.append('print(\'INFO: genes are added from the \\\'reactions\\\' spreadsheet, and gene annotations from the \\\'genes\\\' spreadsheet.\')')

	if outfile:
		with open(outfile, 'w') as fhandle:
			for line in code:
				fhandle.write(line)

	loc = {}
	for line in code:
		try:
			exec(line, globals(), loc)
		except Exception as e:
			print('Error in {:s}'.format(line))
			print(e)
			return None

	return loc['model']
