import numpy
import pandas
import cobra

def ToExcel(model, outfile: str):
	#metabolites
	try:
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

	except:
		print('metabolite dump failed')

	#reactions
	try:
		rxns = {
			'notes' : [],
			'model' : [],
			'_id' : [],
			'name' : [],
			'_metabolites' : [],
			'_lower_bound' : [],
			'_upper_bound' : [],
			'_gpr' : [],
			'subsystem' : [],
			'_annotation' : [],
			}

		for reaction in model.reactions:
			rxn_subs = [ [x,y] for x,y in zip([x._id for x in reaction._metabolites], reaction._metabolites.values()) if y < 0]
			rxn_prod = [ [x,y] for x,y in zip([x._id for x in reaction._metabolites], reaction._metabolites.values()) if y > 0]

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
							' + '.join([ '{:.6g} {:s}'.format(y*-1., x) for x,y in rxn_subs ]),
							' + '.join([ '{:.6g} {:s}'.format(y*+1., x) for x,y in rxn_prod ])))

				elif key == '_gpr':
					rxns['_gpr'].append(reaction._gpr.to_string())

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

	except:
		print('reactions dump failed')

	# genes
	if True:
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

	else:
		print('genes dump failed')

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
