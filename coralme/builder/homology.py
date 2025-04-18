import re
import pandas
import logging

class Homology(object):
	"""
	Homology class for storing information about homology of the
	main and reference organisms.

	This class contains methods to predict and process homology
	of the main and reference organisms. Homology is inferred from
	the reciprocal best hits of a BLAST. The results are used to
	update and complement the attributes of the class Organism.

	Parameters
	----------
	org : str
		Identifier of the main organism. Has to be the same as its
		containing folder name.

	ref : str
		Identifier of the reference organism. Has to be the same as
		its containing folder name.

	evalue : float
		E-value cutoff to call enzyme homologs from the BLAST. Two
		reciprocal best hits are considered homologs if their
		E-value is less than this parameter.
	"""

	def __init__(self, org, ref, evalue = False, verbose = False):
		self.org = org
		self.ref = ref

		column_names = 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'.split('\t')
		self.org_df = pandas.read_csv(self.org.blast_directory + '/org_as_db.txt', sep = '\t', index_col = 0, names = column_names)
		self.ref_df = pandas.read_csv(self.org.blast_directory + '/ref_as_db.txt', sep = '\t', index_col = 0, names = column_names)

		self.get_mutual_hits(evalue = evalue)

	def get_mutual_hits(self, evalue = False, verbose = False):
		# to not import from coralme.builder.homology
		def get_top_hits(df, evalue = False):
			top_hits = {}
			for g in df.index.unique():
				hits = df.loc[g]
				if isinstance(hits, pandas.DataFrame):
					top_hit = hits.sort_values(by = 'evalue').iloc[0]
				else:
					top_hit = hits
				if evalue and top_hit['evalue'] > evalue:
					continue
				top_hits[g] = top_hit
			return top_hits

		org_df = self.org_df
		ref_df = self.ref_df

		if verbose:
			logging.warning('Getting top hits...')

		org_top_hits = get_top_hits(org_df, evalue = evalue)
		ref_top_hits = get_top_hits(ref_df, evalue = evalue)

		if verbose:
			logging.warning('Getting reciprocal top hits...')

		mutual_hits = {}
		for g, hit in org_top_hits.items():
			hit_g = hit['sseqid']

			if hit_g not in ref_top_hits:
				continue

			if g == ref_top_hits[hit_g]['sseqid']:
				mutual_hits[g] = {}
				mutual_hits[hit_g] = {}

				mutual_hits[g]['query'] = hit_g
				mutual_hits[g]['evalue'] = hit['evalue']

				mutual_hits[hit_g]['query'] = g
				mutual_hits[hit_g]['evalue'] = ref_top_hits[hit_g]['evalue']

		if verbose:
			logging.warning('Done')

		mutual_hits_df = pandas.DataFrame.from_dict(mutual_hits).T.sort_index()
		self.mutual_hits_df = mutual_hits_df
		self.mutual_hits = mutual_hits_df['query'].to_dict()

	def get_complex_homology(self):
		org_complexes_df = self.org.complexes_df
		ref_complexes_df = self.ref.complexes_df
		mutual_hits = self.mutual_hits

		# Get ref homology for complexes already annotated in BioCyc
		org_cplx_homolog = {}
		ref_cplx_homolog = {}
		not_annotated_candidates = set()
		warn_candidates = []

		for c, row in org_complexes_df.iterrows():
			if isinstance(row['genes'], float) or not row['genes']:
				continue

			genes = [re.findall(r'.*(?=\(\d*\))', g)[0] for g in row['genes'].split(' AND ')]

			if not genes or not set(genes).issubset(set(mutual_hits.keys())):
				continue  # All org genes must have a hit

			for g in genes:
				if c in org_cplx_homolog:
					break

				ref_gene = mutual_hits[g]
				ref_complexes = ref_complexes_df[ref_complexes_df['genes'].str.contains(ref_gene)]
				for rc, rrow in ref_complexes.iterrows():
					rgenes = [re.findall(r'.*(?=\(\d*\))', g)[0] for g in rrow['genes'].split(' AND ')]
					if not set(rgenes).issubset(set(mutual_hits.keys())):
						continue  # All ref genes must have a hit
					ogenes = [mutual_hits[og] for og in genes]
					if set(rgenes) == set(ogenes):  # Complex identified
						org_cplx_homolog[c] = rc
						ref_cplx_homolog[rc] = c
					else:
						if rc not in not_annotated_candidates:
							warn_candidates.append({
								'complex': c,
								'reference_complex' : rc
								})
						not_annotated_candidates.add(rc)

		not_annotated_candidates = not_annotated_candidates.difference(set(ref_cplx_homolog.keys()))
		logging.warning('{} complexes were mapped successfully'.format(len(org_cplx_homolog)))

		if warn_candidates:
			self.org.curation_notes['org.get_complex_homology'].append({
				'msg':'Some complexes were partial hits in the BLAST',
				'triggered_by':warn_candidates,
				'importance':'medium',
				'to_do':'Curate these manually in protein_corrections.txt'})
		for i in self.org.generic_dict.keys():
			org_cplx_homolog[i] = i
			ref_cplx_homolog[i] = i
		self.org_cplx_homolog = org_cplx_homolog
		self.ref_cplx_homolog = ref_cplx_homolog
		self.not_annotated_candidates = not_annotated_candidates
