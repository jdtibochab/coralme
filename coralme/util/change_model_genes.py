#!/usr/bin/python3
import cobra
import glob
import pandas
import re
import tqdm
from Bio import SeqIO

# Written originally by Rodrigo Santibanez
def get_orthofinder_results(path, save = False):
	# index: Orthogroups; columns: organisms; data: single string, comma separated
	data = pandas.read_csv(path, sep = '\t', index_col = 0)

	# convert NaN to strings -> split string -> get unique IDs -> convert to list -> collapse columns into one
	fn = lambda x: list(set(str(x).split(', ')))
	tmp = data.map(fn).sum(axis = 1)

	# remove NaN and duplicates from lists
	fn = lambda x: [ y for y in list(set(x)) if y != 'nan' ]
	tmp = tmp.apply(fn)

	# explode lists
	tmp = tmp.reset_index(drop = False).explode(0)

	# save
	tmp = tmp[[0, 'Orthogroup']].sort_values([0, 'Orthogroup']).reset_index(drop = True)
	if save:
		tmp.to_csv('01.proteinID-to-orthogroupID.txt', sep = '\t', index = False, header = None)
	return tmp

def _get_genbank_features_as_dataframe(genbank_path):
	contigs = []

	for contig in SeqIO.parse(genbank_path, 'genbank'):
		contigs.append(contig)

	dct = {}
	for contig in contigs:
		for feature in contig.features:
			if feature.type == 'CDS':
				# if feature.qualifiers.get('protein_id', [None])[0] is not None:
				#     if feature.qualifiers.get('old_locus_tag', None) is not None:
				protID = feature.qualifiers.get('protein_id', [None])[0]
				if protID is not None:
					dct[protID] = feature.qualifiers.get('locus_tag', [None])[0]
				else:
					dct[feature.qualifiers.get('locus_tag', [None])[0]] = feature.qualifiers.get('locus_tag')[0]

	return pandas.DataFrame([dct]).T

def get_features_from_genbanks(reference_genbank_path : str = None, genbanks_target_directory : list = []):
	"""
	This assumes gene IDs in M-model are contained in protein_id/locus_tag only
	"""
	data = [_get_genbank_features_as_dataframe(reference_genbank_path)]

	files = sorted(glob.glob(genbanks_target_directory + '/*'))
	for genbank in files:
		data.append(_get_genbank_features_as_dataframe(genbank))

	tmp = pandas.concat(data, axis = 1)
	tmp.columns = ['Original'] + [ x.split('/')[1].split('.gb')[0] for x in files ]
	return tmp

def _get_homologues(x):
	if any(isinstance(x, list) for x in x.dropna().tolist()):
		return [ x for y in x.dropna().tolist() for x in y ]
	else:
		return x.dropna().tolist()

def map_features_between_assemblies(mapping, genbank_features, save = False):
	mapping = mapping.groupby('Orthogroup').agg(lambda x: x.tolist())

	data = []
	for idx, value in tqdm.tqdm(mapping.iterrows(), total = mapping.shape[0]):
		tmp = genbank_features[genbank_features.index.isin([ x.replace('kb|', '') for x in value.tolist()[0]])].agg(lambda x: _get_homologues(x))

		if type(tmp) is pandas.core.series.Series:
			tmp = tmp.to_frame().T

		# tmp.index = [idx] # why?!
		data.append(tmp)

	tmp = pandas.concat(data)
	if save:
		tmp.to_csv('02.locustagID-to-locustagID.txt', sep = '\t', index = False)
	return tmp

def map_genes_between_models(model, mapping_features, save = False):
	model = cobra.io.load_json_model(model)
	genes = [ x.id for x in model.genes ]

	def fn(x, gene):
		if isinstance(x['Original'], list) and gene in x['Original']:
			return True
		if isinstance(x['Original'], str) and gene == x['Original']:
			return True
		return False

	data = []
	found_genes = []
	missing_genes = []

	for gene in tqdm.tqdm(genes):
		tmp = mapping_features[mapping_features.apply(lambda x: fn(x, gene), axis = 1)]
		data.append(tmp)

		shape = tmp.shape
		if shape[0] == 0:
			missing_genes.append(gene)
		else:
			found_genes.append(gene)

	tmp = pandas.concat(data)#.drop('Original', axis = 1)
	tmp = tmp.explode('Original')
	tmp = tmp.map(lambda x: ' or '.join(x) if type(x) == list else x).drop_duplicates()
	if save:
		tmp.to_csv('03.mapping-refmodel-to-others.txt', sep = '\t')

	#tmp5 = tmp2[~tmp2['Original'].astype(str).str.contains('|'.join(genes)) & tmp2['Original'].notna()]
	#tmp5 = tmp5.drop('Original', axis = 1)
	#tmp5 = tmp5.map(lambda x: ' or '.join(x) if type(x) == list else x).drop_duplicates()
	#tmp5.to_csv('05.mapping-of-genes-not-in-Mmodel.txt', sep = '\t')

	return tmp, found_genes, missing_genes

def convert_gene_IDs(target, ref_model, mapping_genes, do_not_convert_gprs = ['rtranscription', 'dreplication', 'pbiosynthesis']):
	deleted = []
	model = cobra.io.load_json_model(ref_model)
	genes = [ x.id for x in model.genes ]

	dct = { re.compile(r'\b{:s}\b'.format(row['Original'])):row[target] for idx, row in mapping_genes[['Original', target]].iterrows() }
	dct

	## Modify IDs for consistency with genbank
	for rxn in list(model.reactions):
		if rxn.id in do_not_convert_gprs:
			continue
		if rxn.gene_reaction_rule != '':
			for key, value in dct.items():
				if value == '':
					value = 'NULL'
				if rxn.gene_reaction_rule == value:
					continue
				if len(re.findall(key, rxn.gene_reaction_rule)) == 0:
					continue
				try:
					new_gpr = re.sub(key, value, rxn.gene_reaction_rule).replace('(NULL or ', '(').replace('NULL or ', '').replace(' or NULL', '')
					rxn.gene_reaction_rule = new_gpr
				except TypeError:
					print(rxn.id, new_gpr)

	for rxn in list(model.reactions):
		rxn.update_genes_from_gpr()

	# remove reactions with "and NULL" <- pseudogenes associated to the model
	for rxn in list(model.reactions):
		if rxn.id in do_not_convert_gprs:
			continue
		if 'NULL' in rxn.gene_reaction_rule:
			# rxn.remove_from_model()
			deleted.append((rxn.id, rxn.gene_reaction_rule))
			try:
				new_gpr = rxn.gene_reaction_rule.replace('(NULL and', '(').replace(' and NULL', '').replace('NULL and ', '').replace('NULL', '')
				rxn.gene_reaction_rule = new_gpr
			except TypeError:
				print(rxn.id, new_gpr)

	# update_genes_from_gpr adds genes to model.genes
	cobra.manipulation.remove_genes(model, [g for g in model.genes if not g.reactions])
	cobra.io.save_json_model(model, 'New-M-models/m_model_{:s}.json'.format(target))
	del model
	return (target, deleted)

def run_all(orthogroups, ref_model, new_model_ids, ref_genbank, genbanks_directory):
	mapping_orthologues = get_orthofinder_results(orthogroups)
	genbank_features = get_features_from_genbanks(reference_genbank_path = ref_genbank, genbanks_target_directory = genbanks_directory)
	mapping_features = map_features_between_assemblies(mapping_orthologues, genbank_features)
	mapping_genes, genes, missing = map_genes_between_models(ref_model, mapping_features)

	res = []
	for idx in new_model_ids:
		if idx in mapping_genes.columns[1:]: # 'Original' should be position 0
			res.append(convert_gene_IDs(idx, ref_model, mapping_genes))
	return res

if __name__ == '__main__':
	# to be completed
	run_all()
