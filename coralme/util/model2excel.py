import numpy
import pandas
import cobra
import coralme

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
			else:
				mets[key].append(metabolite.__dict__[key])

	mets = pandas.DataFrame.from_dict(mets)
	# formating annotations
	tmp = pandas.json_normalize(mets['_annotation'])
	tmp = tmp.astype(str)
	tmp = tmp.apply(lambda col: col.str.replace('[\'', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\']', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\', \'', ';', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('nan', '', regex = False))
	# add annotations
	mets = pandas.concat([mets.iloc[:, :-1], tmp], axis = 1)
	return mets

def _rxns2df(model):
	#reactions
	rxns = {
		'notes' : [],
		'model' : [],
		'_id' : [],
		'name' : [],
		'_metabolites' : [],
		'_lower_bound' : [],
		'_upper_bound' : [],
		'_gpr' : [],
		'_cofactors' : [],
		'subsystem' : [],
		'_annotation' : [],
		}

	for reaction in model.reactions:
		rxn_subs = [ [x,y] for x,y in zip([x._id for x in reaction._metabolites], reaction._metabolites.values()) if y < 0]
		rxn_prod = [ [x,y] for x,y in zip([x._id for x in reaction._metabolites], reaction._metabolites.values()) if y > 0]

		if isinstance(model, coralme.core.model.MEModel):
			rxn_subs = [ [x[0],x[1]] if isinstance(x[1], (int, float)) else [x[0],x[1].subs(model.default_parameters)] for x in rxn_subs ]
			rxn_prod = [ [x[0],x[1]] if isinstance(x[1], (int, float)) else [x[0],x[1].subs(model.default_parameters)] for x in rxn_prod ]

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
				rxns['_metabolites'].append(
					'{:s} = {:s}'.format(
						' + '.join([ '{:.6g} {:s}'.format(y*-1., x) if isinstance(y, (int, float)) else '[{:s}] {:s}'.format(str(y*-1), x) for x,y in rxn_subs ]),
						' + '.join([ '{:.6g} {:s}'.format(y*+1., x) if isinstance(y, (int, float)) else '[{:s}] {:s}'.format(str(y*+1), x) for x,y in rxn_prod ])))

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
	# formating annotations
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
		'name' : [],
		'_annotation' : [],
		}

	for gene in model.genes:
		for key in genes.keys():
			genes[key].append(gene.__dict__[key])

	genes = pandas.DataFrame.from_dict(genes)
	# formating annotations
	tmp = pandas.json_normalize(genes['_annotation'])
	tmp = tmp.astype(str)
	tmp = tmp.apply(lambda col: col.str.replace('[\'', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\']', '', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('\', \'', ';', regex = False))
	tmp = tmp.apply(lambda col: col.str.replace('nan', '', regex = False))
	# add annotations
	genes = pandas.concat([genes.iloc[:, :-1], tmp], axis = 1)
	return genes

def ToExcel(model, outfile: str):
	mets = _mets2df(model)
	rxns = _rxns2df(model)
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
