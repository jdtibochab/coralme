import numpy
import pandas
import sympy

def _mets2df(model):
	#metabolites
	mets = {
		'notes' : [],
		'model' : [],
		'_id' : [],
		'name' : [],
		'formula' : [],
		'compartment' : [],
		'charge' : [],
		'_annotation' : [],
		}

	for metabolite in model.metabolites:
		for key in list(mets.keys()):
			if key in ['model', 'notes']:
				mets[key].append(numpy.nan)
			elif key == 'charge':
				mets[key].append(0) if metabolite.__dict__[key] is None else mets[key].append(metabolite.__dict__[key])
			else:
				mets[key].append(metabolite.__dict__[key])

	mets = pandas.DataFrame.from_dict(mets)
	# formatting annotations
	tmp = pandas.json_normalize(mets['_annotation'])
	tmp = tmp.astype(str)
	tmp = tmp.apply(lambda col: col.str.replace('[\'', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\']', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\', \'', ';', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('nan', '', regex = False))
	# add annotations
	mets = pandas.concat([mets.iloc[:, :-1], tmp], axis = 1)
	return mets

def _rxns2df(model, keys = { sympy.Symbol('mu', positive = True) : 1.0 }):
	#reactions
	rxns = {
		'notes' : [],
		'model' : [],
		'_id' : [],
		'name' : [],
		'_metabolites' : [],
		'_substrates' : [],
		'_products' : [],
		'_lower_bound' : [],
		'_upper_bound' : [],
		'_gpr' : [],
		'_cofactors' : [],
		'subsystem' : [],
		'_annotation' : [],
		}

	for reaction in model.reactions:
		rxn_subs = [ [x._id,y.subs(keys)] if hasattr(y, 'subs') else [x._id,y] for x,y in reaction.metabolites.items() if y < 0]
		rxn_prod = [ [x._id,y.subs(keys)] if hasattr(y, 'subs') else [x._id,y] for x,y in reaction.metabolites.items() if y > 0]

		for key in rxns.keys():
			if key in ['model']:
				if reaction.objective_coefficient != 0.:
					rxns[key].append('__OBJ__')
				else:
					rxns[key].append(numpy.nan)

			elif key in ['notes']:
				if reaction.objective_coefficient != 0.:
					rxns[key].append('__BIOMASS__')
				else:
					rxns[key].append(numpy.nan)

			elif key == '_metabolites':
				tmp_subs = ' + '.join([ '{:.6g} {:s}'.format(y*-1., x) if isinstance(y, (int, float, sympy.core.numbers.Float)) else '[{:s}] {:s}'.format(str(y*-1), x) for x,y in rxn_subs ])
				tmp_prod = ' + '.join([ '{:.6g} {:s}'.format(y*+1., x) if isinstance(y, (int, float, sympy.core.numbers.Float)) else '[{:s}] {:s}'.format(str(y*+1), x) for x,y in rxn_prod ])
				rxns['_metabolites'].append('{:s} = {:s}'.format(tmp_subs, tmp_prod))
			elif key == '_substrates':
				rxns['_substrates'].append('{:s}'.format(tmp_subs))
			elif key == '_products':
				rxns['_products'].append('{:s}'.format(tmp_prod))

			elif key == '_lower_bound':
				bound = reaction._lower_bound
				bound = bound.magnitude if hasattr(bound, 'subs') else bound
				bound = '{:.6g}'.format(bound.subs(model.default_parameters).subs(keys) if hasattr(bound, 'subs') else bound)
				rxns['_lower_bound'].append(bound)

			elif key == '_upper_bound':
				bound = reaction._upper_bound
				bound = bound.magnitude if hasattr(bound, 'subs') else bound
				bound = '{:.6g}'.format(bound.subs(model.default_parameters).subs(keys) if hasattr(bound, 'subs') else bound)
				rxns['_upper_bound'].append(bound)

			elif key == '_gpr':
				rxns['_gpr'].append(reaction._gpr.to_string())

			elif key == '_cofactors':
				if hasattr(reaction, '_cofactors'):
					rxns['_cofactors'].append(reaction._cofactors.to_string())
				else:
					rxns['_cofactors'].append(None)

			else:
				rxns[key].append(reaction.__dict__[key])

	rxns = pandas.DataFrame.from_dict(rxns)
	# formatting annotations
	tmp = pandas.json_normalize(rxns['_annotation'])
	tmp = tmp.astype(str)
	tmp = tmp.apply(lambda col: col.str.replace('[\'', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\']', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\', \'', ';', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('nan', '', regex = False))
	# add annotations
	rxns = pandas.concat([rxns.iloc[:, :-1], tmp], axis = 1)
	return rxns

def _genes2df(model):
	# genes
	genes = {
		'_id' : [],
		# 'name' : [], # it is included in _annotation
		'functional' : [],
		'_annotation' : [],
		}

	for gene in model.genes:
		for key in genes.keys():
			if key == 'functional':
				genes['functional'].append(gene.functional)
			else:
				genes[key].append(gene.__dict__[key])

	genes = pandas.DataFrame.from_dict(genes)
	# formatting annotations
	tmp = pandas.json_normalize(genes['_annotation'])
	tmp = tmp.astype(str)
	tmp = tmp.apply(lambda col: col.str.replace('[\'', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\']', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\', \'', ';', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('nan', '', regex = False))
	# add annotations
	genes = pandas.concat([genes.iloc[:, :-1], tmp], axis = 1)
	return genes

def ToExcel(model, outfile: str, keys : dict = dict()):
	if not outfile.endswith('.xlsx'):
		raise('Filename must end with \'.xlsx\'.')

	mets = _mets2df(model)
	rxns = _rxns2df(model, keys = keys)
	genes = _genes2df(model)

	if outfile.endswith('.xlsx'):
		with open(outfile, 'wb') as outfile:
			writer = pandas.ExcelWriter(outfile, engine = 'xlsxwriter')

			rxns.to_excel(writer, index = False, sheet_name = 'reactions')
			mets.to_excel(writer, index = False, sheet_name = 'metabolites')
			genes.to_excel(writer, index = False, sheet_name = 'genes')

			for data, sheet in zip([ rxns, mets, genes ], [ 'reactions', 'metabolites', 'genes' ]):
				(max_row, max_col) = data.shape

				# Get the xlsxwriter workbook and worksheet objects
				workbook  = writer.book
				worksheet = writer.sheets[sheet]

				# Freeze first row
				worksheet.freeze_panes(1, 0)

				# Set the autofilter
				worksheet.autofilter(0, 0, max_row, max_col - 1)

				# Make the columns wider for clarity
				worksheet.set_column_pixels(0,  max_col - 1, 96)

				# Set zoom level
				worksheet.set_zoom(120)

			# Close the Pandas Excel writer and output the Excel file.
			writer.close()

	return None
