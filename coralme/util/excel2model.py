import re
import sys
import numpy
import cobra
import pandas
# from IPython.display import display, HTML

def FromExcel(infile:str, model_name:str, outfile:str, f_replace:dict = {}, debug:bool = False, report_missing_mets:bool = False) -> cobra.Model:
	code = []
	code.append('import cobra\n')
	code.append('from cobra.core.gene import GPR\n')
	#code.append('try:\n')
	#code.append('\tdel model\n')
	#code.append('except:\n')
	#code.append('\tpass\n')
	code.append('model = cobra.Model(\'{:s}\')\n'.format(model_name))

	# metabolites
	df_mets = pandas.read_excel(infile, sheet_name = 'metabolites', engine = 'openpyxl').dropna(axis = 0, how = 'all').fillna('').reset_index()

	# apply replace
	if f_replace:
		df_mets['_id'] = df_mets['_id'].replace(f_replace)

	base = '_{:s} = cobra.Metabolite(\'{:s}\', formula = \'{:s}\', name = \'{:s}\', compartment = \'{:s}\', charge = {:d})\n' \
		'_{:s}.annotation = {:s}\n' \
		'model.add_metabolites([_{:s}])\n'

	for idx in df_mets.index:
		#if 'model' in df_mets.columns:
		if 'model' in df_mets.columns and (df_mets.iloc[idx]['model'] == '__REMOVE__' or df_mets.iloc[idx]['model'] == '__TEST__'):
			continue
		else:
			annots = {}
			if 'sbo' in df_mets.columns and df_mets.iloc[idx]['sbo'] != '': # SBML: cannot be empty string
				annots['sbo'] = df_mets.iloc[idx]['sbo']
			if 'pubchem.compound' in df_mets.columns and isinstance(df_mets.iloc[idx]['pubchem.compound'], float):
				annots['pubchem.compound'] = ['{:.0f}'.format(df_mets.iloc[idx]['pubchem.compound'])]
			if 'kegg.compound' in df_mets.columns:
				annots['kegg.compound'] = df_mets.iloc[idx]['kegg.compound'].split(';') if (df_mets.iloc[idx]['kegg.compound'] != '') else []
			if 'seed.compound' in df_mets.columns:
				annots['seed.compound'] = df_mets.iloc[idx]['seed.compound'].split(';') if (df_mets.iloc[idx]['seed.compound'] != '') else []
			if 'inchikey' in df_mets.columns:
				annots['inchikey'] = df_mets.iloc[idx]['inchikey'].split(';') if (df_mets.iloc[idx]['inchikey'] != '') else []
			if 'inchi' in df_mets.columns:
				annots['inchi'] = df_mets.iloc[idx]['inchi'].split(';') if (df_mets.iloc[idx]['inchi'] != '') else []
			if 'bigg.metabolite' in df_mets.columns and df_mets.iloc[idx]['bigg.metabolite'] != '': # SBML: cannot be empty string
				annots['bigg.metabolite'] = df_mets.iloc[idx]['bigg.metabolite']
			if 'SMILES' in df_mets.columns:
				annots['SMILES'] = df_mets.iloc[idx]['SMILES'].split(';') if (df_mets.iloc[idx]['SMILES'] != '') else []
			if 'chebi' in df_mets.columns:
				annots['chebi'] = df_mets.iloc[idx]['chebi'].split(';') if (df_mets.iloc[idx]['chebi'] != '') else []
			if 'biocyc' in df_mets.columns:
				annots['biocyc'] = df_mets.iloc[idx]['biocyc'].split(';') if (df_mets.iloc[idx]['biocyc'] != '') else []
			if 'hmdb' in df_mets.columns:
				annots['hmdb'] = df_mets.iloc[idx]['hmdb'].split(';') if (df_mets.iloc[idx]['hmdb'] != '') else []
			if 'lipidmaps' in df_mets.columns:
				annots['lipidmaps'] = df_mets.iloc[idx]['lipidmaps'].split(';') if (df_mets.iloc[idx]['lipidmaps'] != '') else []
			if 'metanetx.chemical' in df_mets.columns:
				annots['metanetx.chemical'] = df_mets.iloc[idx]['metanetx.chemical'].split(';') if (df_mets.iloc[idx]['metanetx.chemical'] != '') else []
			if 'seed.compound' in df_mets.columns:
				annots['seed.compound'] = df_mets.iloc[idx]['seed.compound'].split(';') if (df_mets.iloc[idx]['seed.compound'] != '') else []
			if 'reactome' in df_mets.columns:
				annots['reactome'] = df_mets.iloc[idx]['reactome.compound'].split(';') if (df_mets.iloc[idx]['reactome.compound'] != '') else []
			annots = str(annots).replace('{', '{\n\t').replace('}', '\n\t}').replace(', ', ',\n\t') if annots else '{}'

			try:
				tmp = base.format(
					# ME-model IDs can have invalid python's variable characters
					df_mets.iloc[idx]['_id'].replace('-', '_DASH_').replace('(', '_LPAR_').replace(')', '_RPAR_'),
					df_mets.iloc[idx]['_id'],
					df_mets.iloc[idx]['formula'] if 'formula' in df_mets.columns else '',
					df_mets.iloc[idx]['name'].replace('\'', '\\\''),
					df_mets.iloc[idx]['compartment'] if 'compartment' in df_mets.columns else '',
					int(df_mets.iloc[idx]['charge'] if 'charge' in df_mets.columns else 0) or 0,
					df_mets.iloc[idx]['_id'].replace('-', '_DASH_').replace('(', '_LPAR_').replace(')', '_RPAR_'),
					annots,
					df_mets.iloc[idx]['_id'].replace('-', '_DASH_').replace('(', '_LPAR_').replace(')', '_RPAR_'))
			except:
				raise ValueError('Incorrect format detected at \'{:s}\' metabolite entry.'.format(df_mets.iloc[idx]['_id']))

			code.append(tmp)

	# reactions
	df_rxns = pandas.read_excel(infile, sheet_name = 'reactions', engine = 'openpyxl').dropna(axis = 0, how = 'all').fillna('').reset_index()
	df_rxns['_metabolites'] = df_rxns['_metabolites'].str.replace(r'<=>|->|<-', '=', regex=True)

	base = '\n' \
		'reaction = cobra.Reaction(\'{:s}\')\n' \
		'reaction.name = \'{:s}\'\n' \
		'reaction.subsystem = \'{:s}\'\n' \
		'reaction.lower_bound = {:f}\n' \
		'reaction.upper_bound = {:f}\n' \
		'reaction.add_metabolites({:s})\n' \
		'reaction.annotation = {:s}\n' \
		'reaction.gene_reaction_rule = \'{:s}\'\n' \
		'reaction.cofactors = GPR.from_string(\'{:s}\')\n' \
		'model.add_reactions([reaction])\n'

	uniq_mets = set()
	for idx in df_rxns.index:
		#if 'model' in df_rxns.columns:
		if 'model' in df_rxns.columns and (df_rxns.iloc[idx]['model'] == '__REMOVE__' or df_rxns.iloc[idx]['model'] == '__TEST__'):
			continue
		else:
			# reconstruct reaction from string
			# regex = re.compile(r'([A-Za-z0-9\.\_\-\(\)]+)')
			# regex = re.compile(r'(?:(\d+)\s+)?([A-Za-z][A-Za-z0-9._\-\(\)]*)')
			regex = re.compile(r'(?:(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+)?([A-Za-z][A-Za-z0-9._\-\(\)]*)')
			reaction = df_rxns.iloc[idx]['_metabolites']

			# apply replace
			for key, value in f_replace.items():
				reaction = reaction.replace(key, value)

			try:
				# reaction strings are in the form of "2 A + B = C + D"
				subs = reaction.split('=')[0].strip().split(' + ')
			except Exception as e:
				if df_rxns.iloc[idx]['model'] != '__BIOMASS__' and debug:
					print('Error in:')
					print(df_rxns.iloc[idx].to_frame().T.to_string())
				#traceback.print_exc()

			try:
				# reaction strings are in the form of "2 A + B = C + D"
				prods = reaction.split('=')[1].strip().split(' + ')
			except Exception as e:
				if df_rxns.iloc[idx]['model'] != '__BIOMASS__' and debug:
					print('Error in:')
					print(df_rxns.iloc[idx].to_frame().T.to_string())
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
			if 'model' in df_rxns.columns:
				if df_rxns.iloc[idx]['model'] == '__INVERT__':
					invert = -1
			lower = df_rxns.iloc[idx]['_lower_bound'] if invert == +1 else -1*df_rxns.iloc[idx]['_upper_bound']
			upper = df_rxns.iloc[idx]['_upper_bound'] if invert == +1 else -1*df_rxns.iloc[idx]['_lower_bound']

			mets = {}
			for substrate in subs: # zip(subs[0::2], subs[1::2]):
				if substrate == '': # exchange and other reactions
					continue
				if ' ' in substrate: # metabolite has a coefficient, e.g., "2 A"
					coeff, substrate = substrate.split(' ', 1)
				else:
					coeff = 1
				uniq_mets.add(substrate)
				substrate = substrate.replace('-', '_DASH_').replace('(', '_LPAR_').replace(')', '_RPAR_')
				mets['_' + substrate] = -1*float(coeff)*invert

			for product in prods: # zip(prods[0::2], prods[1::2]):
				if product == '': # exchange and other reactions
					continue
				if ' ' in product: # metabolite has a coefficient, e.g., "2 A"
					coeff, product = product.split(' ', 1)
				else:
					coeff = 1
				uniq_mets.add(product)
				product = product.replace('-', '_DASH_').replace('(', '_LPAR_').replace(')', '_RPAR_')
				mets['_' + product] = +1*float(coeff)*invert

			mets = str(mets).replace('{', '{\n\t').replace('}', '\n\t}').replace(', ', ',\n\t').replace('\'', '')

			annots = {}
			if 'sbo' in df_rxns.columns:
				annots['sbo'] = df_rxns.iloc[idx]['sbo']
			if 'bigg.reaction' in df_rxns.columns:
				annots['bigg.reaction'] = df_rxns.iloc[idx]['bigg.reaction'] if (df_rxns.iloc[idx]['bigg.reaction'] != '') else []
			if 'biocyc' in df_rxns.columns:
				annots['biocyc'] = df_rxns.iloc[idx]['biocyc'].replace('META:', '').split(';') if (df_rxns.iloc[idx]['biocyc'] != '') else []
			if 'ec-code' in df_rxns.columns:
				annots['ec-code'] = df_rxns.iloc[idx]['ec-code'].split(';') if (df_rxns.iloc[idx]['ec-code'] != '') else []
			if 'kegg.reaction' in df_rxns.columns:
				annots['kegg.reaction'] = df_rxns.iloc[idx]['kegg.reaction'].split(';') if (df_rxns.iloc[idx]['kegg.reaction'] != '') else []
			annots = str(annots).replace('{', '{\n\t').replace('}', '\n\t}').replace(', ', ',\n\t') if annots else '{}'

			try:
				tmp = base.format(
					df_rxns.iloc[idx]['_id'],
					df_rxns.iloc[idx]['name'].replace('\'', '\\\''),
					df_rxns.iloc[idx]['subsystem'] if 'subsystem' in df_rxns.columns else '',
					lower, upper, mets, annots,
					df_rxns.iloc[idx]['_gpr'] if '_gpr' in df_rxns.columns else '',
					df_rxns.iloc[idx]['_cofactors'] if '_cofactors' in df_rxns.columns else '')
			except:
				for col in df_rxns.columns:
					print(col, df_rxns.iloc[idx][col], type(df_rxns.iloc[idx][col]))
				raise ValueError('Incorrect format detected at \'{:s}\' reaction entry.'.format(df_rxns.iloc[idx]['_id']))

			code.append(tmp)

		if 'model' in df_rxns.columns:
			if df_rxns.iloc[idx]['model'] == '__OBJ__':
				code.append('\nmodel.objective = \'{:s}\'\n'.format(df_rxns.iloc[idx]['_id']))

	if report_missing_mets:
		for met in uniq_mets.difference(df_mets['_id']):
			print(met)

	# gene annotations
	df_genes = pandas.read_excel(infile, sheet_name = 'genes', engine = 'openpyxl', dtype = str).dropna(axis = 0, how = 'all').fillna('').reset_index()
	df_genes['_id'] = df_genes['_id'].apply(lambda x: x.replace('(', '').replace(')', '').split(' or '))
	df_genes = df_genes.explode('_id')

	base = 'try:\n' \
		'\tmodel.genes.get_by_id(\'{:s}\').functional = {:s}\n' \
		'\tmodel.genes.get_by_id(\'{:s}\').annotation = {:s}\n' \
		'except:\n' \
		'\tprint(\'INFO: gene ID {:s} not associated to any reaction in the M-model.\')\n\n'

	code.append('\n')
	for idx in df_genes.index:
		annots = {}
		for annot in [ x for x in df_genes.columns if x not in [ 'index', '_id', 'functional' ] ]:
			if isinstance(df_genes.iloc[idx][annot], float):
				annots[annot] = '{:.0f}'.format(df_genes.iloc[idx][annot]) if (df_genes.iloc[idx][annot] != '') else []
			elif ';' in df_genes.iloc[idx][annot]:
				annots[annot] = df_genes.iloc[idx][annot].split(';') if (df_genes.iloc[idx][annot] != '') else []
			elif ',' in df_genes.iloc[idx][annot]:
				annots[annot] = df_genes.iloc[idx][annot].split(',') if (df_genes.iloc[idx][annot] != '') else []
			else:
				annots[annot] = [df_genes.iloc[idx][annot]] if (df_genes.iloc[idx][annot] != '') else []

		annots = str(annots).replace('{', '{\n\t\t').replace('}', '\n\t\t}').replace(', ', ',\n\t\t') if annots else '{}'
		if 'functional' in df_genes.columns:
			code.append(base.format(df_genes.iloc[idx]['_id'], df_genes.iloc[idx]['functional'], df_genes.iloc[idx]['_id'], annots, df_genes.iloc[idx]['_id']))
		else:
			code.append(base.format(df_genes.iloc[idx]['_id'], 'True', df_genes.iloc[idx]['_id'], annots, df_genes.iloc[idx]['_id']))

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
