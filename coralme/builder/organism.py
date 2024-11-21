#!/usr/bin/python3
import os
import re
import random
import io
import anyconfig
import numpy
import collections

import Bio
import cobra
import pandas
import tqdm
bar_format = '{desc:<75}: {percentage:.1f}%|{bar:10}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
import coralme

import warnings
try:
    warnings.simplefilter(action = 'ignore', category = Bio.BiopythonWarning)
except:
    warnings.warn("This biopython version does not allow for correct warning handling. Biopython >=1.80 is suggested.")

import logging
log = logging.getLogger(__name__)

element_types = {'CDS', 'rRNA','tRNA', 'ncRNA','misc_RNA','RNA','tmRNA'}

class Organism(object):
    """Organism class for storing information about an organism

    This class acts as a database containing all necessary information
    to reconstruct a ME-model. It is used to retrieve and store
    information of the main (org) and the reference (ref) organisms.
    Information in Organism is read and manipulated by methods in
    the MEBuilder class.

    Parameters
    ----------
    config : dict
        Dictionary containing configuration and settings.

    is_reference : bool
        If True, process as reference organism.
    """

    def __init__(self, config, is_reference):
        if is_reference:
            # If Organism is reference, set some default values for future processing
            if bool(config.get('dev_reference', False)) and not bool(config.get('user_reference', False)):
                # If no reference set, use the default iJL1678b E. coli ME-model files as reference
                self.id = 'iJL1678b'
                config["locus_tag"] = "locus_tag"
            elif not bool(config.get('dev_reference', False)) and bool(config.get('user_reference', False)):
                # If User set a reference, use manually set values
                # Set ID
                self.id = config['user_reference']
                config = config.copy()
                # Update the configuration using the configuration from the User's defined reference
                for input_file in [config['user_reference'] + "/organism.json", \
                                    config['user_reference'] + "/input.json"]:
                    with open(input_file, 'r') as infile:
                        config.update(anyconfig.load(infile))
            elif bool(config.get('dev_reference', False)) and bool(config.get('user_reference', False)):
                # Reference configuration is contradictory. Using default reference instead.
                logging.warning('The \'dev_reference\' and \'user-reference\' options are mutually exclusive. Using default iJL1678b as reference.')
                self.id = 'iJL1678b'
                config["locus_tag"] = "locus_tag"
            else:
                # In any other case, use the default reference instead
                logging.warning('No reference was set. Using default iJL1678b as reference.')
                self.id = 'iJL1678b'
                config["locus_tag"] = "locus_tag"
        else:
            # If it is not the reference, use the organism files
            self.id = config['ME-Model-ID']

        # Set initial properties
        self.is_reference = is_reference
        self.curation_notes = collections.defaultdict(list)
        self.config = config

        #if self.is_reference:
            ## If it is the reference, use the gene ID conventions of the reference
            #self.locus_tag = config.get('reference_tag','locus_tag')
        #else:
            # If it is not the reference, use the gene ID conventions of the organism
        self.locus_tag = config.get('locus_tag','locus_tag')

        # Create the protein localization interpreter
        # TODO: Move this to a more standard location, builder/?
        data = \
            'code,interpretation,gram\n' \
            'CCI-CW-BAC-POS-GP,Cell_Wall,pos\n' \
            'CCI-OUTER-MEM-GN,Outer_Membrane,neg\n' \
            'CCI-PM-BAC-NEG-GN,Inner_Membrane,neg\n' \
            'CCI-PERI-BAC-GN,Periplasm,neg\n' \
            'CCI-PM-BAC-POS-GP,Plasma_Membrane,pos\n' \
            'CCI-EXTRACELLULAR-GP,Extracellular_Space,pos\n' \
            'CCO-MEMBRANE,Membrane,'
        self.location_interpreter = pandas.read_csv(io.StringIO(data), index_col=0)

        # Initialize the curation fields
        self.manual_curation = coralme.builder.curation.CurationList()

    @property
    def directory(self):
        """Working file directory with reconstruction files"""
        if self.is_reference and self.id == 'iJL1678b':
            try:
                from importlib.resources import files
            except ImportError:
                from importlib_resources import files
            return str(files("coralme") / self.id) + "-ME/building_data/"
        else:
            return self.config.get('out_directory', self.id) + "/building_data/"
        #return self.id + "/building_data/"

    @property
    def blast_directory(self):
        """Directory containing the blast files and results"""
        if not self.is_reference:
            return self.config.get('out_directory', self.id) + "/blast_files_and_results/"
        # If reference organism, no blast information
        return None

    @property
    def _complexes_df(self):
        """
        Returns the dataframe containing the protein complex information from the file.
        If the file does not exist, it initializes it.
        """
        filename = self.directory + "protein_complexes.txt"
        if os.path.isfile(filename):
            # If protein_complexes.txt exists in the files directory,
            # use it as the final complexes dataframe.
            return pandas.read_csv(
                filename, index_col=0, sep="\t",comment='#'
            ).fillna('')
        # If it does not exist, generate it from raw files
        return self.generate_complexes_df()

    @property
    def _protein_mod(self):
        """
        Returns the dataframe containing the protein complex modifications from the file.
        If the file does not exist, it initializes it.
        """
        filename = self.directory + "protein_modification.txt"
        if os.path.isfile(filename):
            # If protein_modification.txt exists in the files directory,
            # use it as the final protein modification dataframe.
            return pandas.read_csv(
                filename, index_col=0, sep="\t",comment='#'
            )
        # If it does not exist, generate it from raw files
        return pandas.DataFrame.from_dict(
                {
                    "Modified_enzyme": {},
                    "Core_enzyme": {},
                    "Modifications": {},
                    "Source": {},
                }
            ).set_index("Modified_enzyme")

    @property
    def _TU_df(self):
        """
        Returns the dataframe containing the TU information from the file.
        If the file does not exist, it initializes it.
        """
        if self.is_reference:
            # Reference file should contain TUs_from_biocyc.txt
            filename = self.directory + "TUs_from_biocyc.txt"
        else:
            # Use TUs dataframe file path as defined
            filename = self.config.get('df_TranscriptionalUnits', self.directory + "TUs_from_biocyc.txt")

        if os.path.isfile(filename):
            # If TUs dataframe exists, read it
            tmp = pandas.read_csv(filename, index_col = 0, sep = "\t")
            tmp = tmp.dropna(subset=['start', 'stop', 'genes'], how = 'any')
            return tmp
        else:
            # If it does not exist, generate it from raw files
            return self.get_TU_df()

    @property
    def _m_model(self):
        """
        Returns the M-model.
        """
        if self.id == 'iJL1678b':
            # If reference organism is iJL1678b, read it from m_model.json
            model = self.directory + 'm_model.json'
        else:
            # Read M-model from the configuration file path
            model = self.config['m-model-path']

        if model.endswith('.json'):
            # Read from JSON model
            return cobra.io.load_json_model(model)
        elif model.endswith('.xml'):
            # Read from SBML model
            return cobra.io.read_sbml_model(model)
        else:
            raise ValueError('M-model input file must be json or xml format.')

    @property
    def rna_components(self):
        """
        Returns the RNAs from the identified products.
        """
        product_types = self.product_types
        # Return all products whose type is RNA
        return set(g for g,t in product_types.items() if 'RNA' in t)

    def get_organism(self):
        """ Processes input files, and creates an instance of
        Organism.
        """
        sep = '~ '*1
        print("{}Processing files for {}...".format(sep,self.id))
        if not self.is_reference:
            logging.warning('Checking folder')
            self.check_folder()
        logging.warning("Loading M-model")
        self.m_model = self._m_model
        logging.warning("Checking M-model")
        self.check_m_model()
        logging.warning("Loading genbank file")
        self.get_genbank_contigs()
        logging.warning("Loading optional files")
        self.load_optional_files()
        logging.warning("Checking gene overlap")
        self.check_gene_overlap()
        logging.warning("Generating complexes dataframe")
        self.complexes_df = self._complexes_df
        logging.warning("Syncing files")
        self.sync_files()
        logging.warning("Looking for duplicates in provided files")
        self.check_for_duplicates()
        logging.warning('Pruning genbank from unwanted feature types')
        self.prune_genbank()
        logging.warning('Completing genbank with provided files')
        self.update_genbank_from_files()

        logging.warning("Updating genes and complexes from genbank")
        self.update_complexes_genes_with_genbank()
        logging.warning("Generating protein modifications dataframe")
        self.protein_mod = self._protein_mod
        logging.warning("Purging genes in optional files")
        self.purge_genes_in_file()

        logging.warning("Loading manual curation")
        self.load_manual_curation()

        logging.warning("Integrating manual metabolites")
        self.modify_metabolites()

        logging.warning("Integrating manual metabolic reactions")
        self.modify_metabolic_reactions()

        logging.warning("Integrating manual complexes")
        self.add_manual_complexes()

        logging.warning("Getting sigma factors from BioCyc")
        self.get_sigma_factors()
        self.get_rpod()
        logging.warning("Getting RNA polymerase from BioCyc")
        self.get_rna_polymerase()

        logging.warning("Updating generics with genbank")
        self.get_generics_from_genbank()

        logging.warning("Generating transcription units dataframe")
        self.TU_df = self._TU_df
        self.get_TU_genes()
        logging.warning("Updating ribosomal proteins with BioCyc")
        self.update_ribosome_stoich()
        logging.warning("Updating protein location with BioCyc")
        self.get_protein_location()
        logging.warning("Updating tRNA synthetases with BioCyc")
        self.get_trna_synthetase()
        logging.warning("Getting lipids")
        self.lipids = self.get_lipids()
        logging.warning("Getting phospholipids")
        self.phospholipids = self.get_phospholipids()
        logging.warning("Updating peptide release factors with BioCyc")
        self.get_peptide_release_factors()
        logging.warning("Complementing non-metabolic metabolites in M-model")
        self.get_nonmetabolic()

        logging.warning("Purging genes in M-model")
        self.purge_genes_in_model()

        logging.warning("Getting enzyme-reaction association")
        self.get_enzyme_reaction_association()

        print("Reading {} done.".format(self.id))

    def get_genbank_contigs(self):
        """ Reads GenBank file as a list of contigs.
        """
        if self.id == 'iJL1678b':
            # If default reference iJL1678b, read from genome.gb
            gb_it = Bio.SeqIO.parse(self.directory + "genome.gb", "gb")
        else:
            # Read from configuration file path
            gb_it = Bio.SeqIO.parse(self.config['genbank-path'], "gb")
        self.contigs = [ i for i in gb_it ]

    def check_folder(self):
        """ Checks that the necessary directories are present.
        """
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)
            logging.warning("{} directory was created.".format(self.directory))
        if not os.path.isdir(self.blast_directory):
            os.makedirs(self.blast_directory)
            logging.warning("{} directory was created.".format(self.blast_directory))

    def check_m_model(self):
        """ Performs checks on the M-model
        """
        m_model = self.m_model

        # Metabolites
        RNA_mets = []
        formula_mets = []
        formulaweight_mets = []
        deadend_mets = []

        for m in tqdm.tqdm(m_model.metabolites,
                           'Checking M-model metabolites...',
                           bar_format = bar_format):
            if m.id.startswith("RNA"):
                RNA_mets.append(m)
            if not m.formula:
                formula_mets.append(m.id)
            try:
                float(m.formula_weight)
            except:
                formulaweight_mets.append(m.id)
            if len(m.reactions) == 0:
                deadend_mets.append(m.id)

        unused_genes = []
        for g in tqdm.tqdm(m_model.genes,
                           'Checking M-model genes...',
                           bar_format = bar_format):
            if not g.reactions:
                unused_genes.append(g.id)

        # Reactions
        subsystem_RXNS = []
        for r in tqdm.tqdm(m_model.reactions,
                           'Checking M-model reactions...',
                           bar_format = bar_format):
            if not r.subsystem:
                subsystem_RXNS.append(r.id)

        # Warnings
        if RNA_mets:
            self.curation_notes['org.check_m_model'].append({
                'msg':"Metabolites starting with the prefix RNA were removed.",
                'triggered_by':[i.id for i in RNA_mets],
                'importance':'high',
                'to_do':'Check whether the model needs the removed metabolites, and change their name. RNA as a prefix is used for TranscribedGene objects.'})
            m_model.remove_metabolites(RNA_mets)
        if formula_mets:
            self.curation_notes['org.check_m_model'].append({
                'msg':"Some metabolites are missing their formula",
                'triggered_by':formula_mets,
                'importance':'critical',
                'to_do':'Correct the formulas of the listed metabolites. Some metabolite formulas are critical for the completion of this pipeline. If homology is ON, this pipeline will try to fill in the formulas from the reference.'})
        if formulaweight_mets:
            self.curation_notes['org.check_m_model'].append({
                'msg':"Some metabolites have a problematic formula. If these metabolites are used in protein modifications, or other subreactions, it will cause an error.",
                'triggered_by':formulaweight_mets,
                'importance':'critical',
                'to_do':'Correct the formulas of the listed metabolites. Some metabolite formulas are critical for the completion of this pipeline. If homology is ON, this pipeline will try to fill in the formulas from the reference.'})
        if subsystem_RXNS:
            self.curation_notes['org.check_m_model'].append({
                'msg':"Some reactions are missing their subsystem",
                'triggered_by':subsystem_RXNS,
                'importance':'high',
                'to_do':'Make sure the subsystems of these reactions are correct'})
        if deadend_mets:
            self.curation_notes['org.check_m_model'].append({
                'msg':"Some metabolites have no reactions associated",
                'triggered_by':deadend_mets,
                'importance':'critical',
                'to_do':'Make sure these metabolites are removed or connected properly'})
        if unused_genes:
            self.curation_notes['org.check_m_model'].append({
                'msg':"Some genes have no reactions associated",
                'triggered_by':unused_genes,
                'importance':'critical',
                'to_do':'Make sure these genes are removed or associated properly'})

    def load_optional_files(self):
        """ Loads optional files.
        """
        logging.warning("Loading gene dictionary")
        self.gene_dictionary = self.read_gene_dictionary(
            self.config.get('biocyc.genes', self.directory + "genes.txt")
        )
        self.gene_sequences = self.read_gene_sequences(
            self.config.get('biocyc.seqs', self.directory + "sequences.fasta")
        )
        logging.warning("Getting proteins from BioCyc")
        self.proteins_df = self.read_proteins_df(
            self.config.get('biocyc.prots', self.directory + "proteins.txt")
        )
        logging.warning("Getting RNAs from BioCyc")
        self.RNA_df = self.read_RNA_df(
            self.config.get('biocyc.RNAs', self.directory + "RNAs.txt")
        )
        logging.warning("Getting transcription units from BioCyc")
        self.TUs = self.read_TU_df(
            self.config.get('biocyc.TUs', self.directory + "TUs.txt")
        )

    def load_manual_curation(self):
        """ Loads manual curation to Organism instance
        """
        coralme.builder.curation.MEManualCuration(self).load_manual_curation()

    def purge_genes_in_file(self):
        """ Checks genes in files and purges problematic ones.
        """
        if self.is_reference:
            return

        all_genenames_in_gb = []
        for i in self.all_genes_in_gb:
            genenames = self.gene_dictionary[self.gene_dictionary['Accession-1'].eq(i)].index
            all_genenames_in_gb += list(genenames)
        warn_products = set(self.gene_dictionary[self.gene_dictionary["Product"] == ''].index)
        warn_replicons = set(self.gene_dictionary[self.gene_dictionary["replicon"] == ''].index)
        warn_rightpos = set(self.gene_dictionary[self.gene_dictionary["Left-End-Position"] == ''].index)
        warn_leftpos = set(self.gene_dictionary[self.gene_dictionary["Right-End-Position"] == ''].index)
        warn_sequences = set(self.gene_dictionary.index) - set(self.gene_sequences.keys()) - set(all_genenames_in_gb)
        warn_genenames = set(self.gene_dictionary[self.gene_dictionary.index == ''].index)

        self.gene_dictionary.drop(list(
            warn_products|warn_genenames|warn_replicons|warn_sequences|warn_rightpos|warn_leftpos
        ),inplace=True)

    def _get_product_type(self,
                         gene_name,
                         gene_dictionary = None,
                         complexes_df = None,
                         RNA_df = None,
                         warn_genes = []):
        """Gets the product type of a gene (RNA, CDS, ...) from the files"""
        if gene_dictionary is None:
            gene_dictionary = self.gene_dictionary
        if complexes_df is None:
            complexes_df = self.complexes_df
        if RNA_df is None:
            RNA_df = self.RNA_df
        row = gene_dictionary.loc[gene_name]
        gene_id = row['Accession-1']
        product = row['Product'].split(' // ')[0]
        ### Try to get product type from gene id of type LOCUST_TAG-RNA
        product_type = ''
        if '-' in product and ' ' not in product:
            product_type = re.findall('[a-zA-Z]+',product.split('-')[-1])
            if product_type:product_type = product_type[0]
        ### Set product type to RNA if it is in ID
        if ' ' in product or ('RNA' not in product and 'MONOMER' not in product) or not product_type:
            if 'RNA' in gene_id \
                    or RNA_df['Gene'].str.match(gene_name).any() \
                    or product in RNA_df.index:
                product_type = 'RNA'
            elif 'MONOMER' in gene_id \
                    or complexes_df['genes'].str.contains('{}\(\d*\)'.format(gene_id),regex=True).any() \
                    or product in complexes_df.index:
                product_type = 'MONOMER'
            else:
                warn_genes.append(gene_id)
                return gene_id,None,None
        return gene_id,product,product_type

    def _correct_product(self,
                        gene_name,
                        product_type,
                        gene_dictionary = None):
        """
        Product identifier contains wrong characters or is missing critical information
        and must be corrected for parsing
        """
        if gene_dictionary is None:
            gene_dictionary = self.gene_dictionary
        ## Correct product. Likely product is a description and not an actual
        ## product ID like GENE-MONOMER or GENE-tRNA
        product = '{}-{}'.format(gene_name,product_type)
        gene_dictionary.at[gene_name,'Product'] = product
        return product

    def _add_entry_to_df(self,
                         df,
                         tmp):
        """Append one entry to the end of a dataframe"""
        indexname = df.index.name
        df = pandas.concat([df,
                              pandas.DataFrame.from_dict(tmp).T],
                             axis = 0, join = 'outer')
        df.index.name = indexname
        return df

    def _add_entry_to_rna(self,
                         gene_id,
                         name,
                         product,
                         RNA_df,
                         source):
        """Append one entry to the RNA dataframe"""
        logging.warning('Adding {} ({}) to RNAs from {}'.format(gene_id,product,source))
        tmp = {product : {"Common-Name": name,
                          "Gene": gene_id}}
        return self._add_entry_to_df(RNA_df,tmp)

    def _add_entry_to_complexes(self,
                               gene_id,
                               name,
                               product,
                               complexes_df,
                               source):
        """Append one entry to the complexes dataframe"""
        if product in complexes_df.index:
            logging.warning('Could not add {} ({}) to complexes from {}. Already in complexes_df'.format(gene_id,product,source))
            return complexes_df
        if isinstance(gene_id,str):
            gene_id = '{}()'.format(gene_id)
        elif isinstance(gene_id,list):
            gene_id = ' AND '.join(['{}()'.format(g) for g in gene_id])
        elif isinstance(gene_id,dict):
            gene_id = ' AND '.join(['{}({})'.format(k,v) for k,v in gene_id.items()])
        else:
            raise TypeError("Unsupported entry to add to complexes of type " + type(gene_id))
        logging.warning('Adding {} ({}) to complexes from {}'.format(product,gene_id,source))
        tmp = {product: {
                "name": name,
                "genes": gene_id,
                "source": source,
                }}
        return self._add_entry_to_df(complexes_df,tmp)

    def _add_entry_to_protein_mod(self,
                                  protein_mod,
                                  mod_complex,
                                  core_enzyme,
                                  mods,
                                  source):
        """Append one entry to the protein modification dataframe"""
        logging.warning('Adding {} to protein_mod from {}'.format(mod_complex, source))
        tmp = {mod_complex: {
                "Core_enzyme": core_enzyme,
                "Modifications": mods,
                "Source": source,
                }}
        return self._add_entry_to_df(protein_mod,tmp)

    def sync_files(self):
        """
        Syncs provided files. Cross checks information between
        optional files and the genome file
        """
        if self.is_reference:
            return

        gene_dictionary = self.gene_dictionary
        RNA_df = self.RNA_df
        complexes_df = self.complexes_df
        product_types = {}
        warn_genes = []
        for gene_name,row in tqdm.tqdm(gene_dictionary.iterrows(),
                           'Syncing optional genes file...',
                           bar_format = bar_format,
                           total=gene_dictionary.shape[0]):
            gene_id = row['Accession-1']
            if not gene_name or isinstance(gene_name,float):
                warn_genes.append(gene_id)
                continue

            gene_id,product,product_type = \
                self._get_product_type(
                         gene_name,
                         gene_dictionary=gene_dictionary,
                         complexes_df=complexes_df,
                         RNA_df=RNA_df,
                         warn_genes=warn_genes)

            if product is None:
                warn_genes.append(gene_id)
                continue
            if ' ' in product or ('RNA' not in product and 'MONOMER' not in product):
                product = \
                    self._correct_product(
                        gene_name,
                        product_type,
                        gene_dictionary)

            product_types[gene_id] = product_type

            ## Sync files
            if 'RNA' in product_type and product not in RNA_df.index:
                RNA_df = \
                    self._add_entry_to_rna(gene_id,
                                           product,
                                           product,
                                           RNA_df,
                                           "Sync")

            elif product_type == 'MONOMER' and product not in complexes_df.index:
                complexes_df = \
                    self._add_entry_to_complexes(gene_id,
                                                 product,
                                                 product,
                                                 complexes_df,
                                                 "Sync")

        self.gene_dictionary = gene_dictionary[pandas.notnull(gene_dictionary.index)]
        self.RNA_df = RNA_df
        self.complexes_df = complexes_df
        self.product_types = product_types

        # Warnings
        if warn_genes:
            self.curation_notes['org.sync_files'].append({
                                'msg':'The types of some genes (e.g. CDS, RNA...) could not be identified. Is Product or Gene Name missing?',
                                'triggered_by':warn_genes,
                                'importance':'medium',
                                'to_do':'Manually fill the products (with types) of these genes in genes.txt'
            })

    def _create_genbank_contig(self,
                               contig_id,
                               seq,
                               name,
                               description,
                               source):
        """Create a Biopython contig"""
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, ExactPosition, SimpleLocation
        new_contig = SeqRecord(seq=seq,
                              id = contig_id,
                              name = name,
                              description = description,
                              annotations = {
                                  'molecule_type' : 'DNA'
                              })

        new_contig.features = [SeqFeature(SimpleLocation(ExactPosition(0),ExactPosition(len(seq))),
                              type='source',
                              id = contig_id,
                              qualifiers = {'note':'Added from {}'.format(source)})]
        return new_contig

    def _create_contig_feature(self,
                                gene_id,
                                seq,
                                strand,
                                feature_type,
                                product_name):
        """Create a Biopython feature to add to a contig"""
        from Bio.SeqFeature import SeqFeature, ExactPosition, SimpleLocation
        return SeqFeature(SimpleLocation(ExactPosition(0),ExactPosition(len(seq)),strand),
                      type=feature_type,
                      id = gene_id,
                      qualifiers = {
                          self.locus_tag:[gene_id],
                          'product':[product_name]
                      })

    def _add_entry_to_genbank(self,
                             gene_id,
                             gene_name,
                             product_type,
                             product_name,
                             row,
                             contigs,
                             gene_sequences):
        """Append one entry to the GenBank"""
        if self.duplicated_genes is not None and gene_id in self.duplicated_genes:
            gene_id = '{};{}'.format(gene_id,gene_name)
        logging.warning('Adding {} to genbank file as {}'.format(gene_id,product_type))
        from Bio.SeqFeature import SeqFeature, CompoundLocation, ExactPosition, FeatureLocation, SimpleLocation
        gene_seq = gene_sequences[gene_name]
        gene_left = int(row['Left-End-Position'])
        gene_right = int(row['Right-End-Position'])

        new_contig = self._create_genbank_contig('{}'.format(gene_id),
                                                 gene_seq.seq,
                                                 gene_seq.name,
                                                 gene_seq.description,
                                                 'BioCyc')
        feature = self._create_contig_feature(gene_id,
                                               gene_seq.seq,
                                               1 if gene_left < gene_right else -1,
                                               product_type,
                                               product_name)
        feature.qualifiers['transl_table'] = [self.transl_table]
        new_contig.features += [feature]
        contigs.append(new_contig)

    def _get_product_name_if_present(self,
                         gene_id,
                         product,
                         product_type,
                         query_types,
                         dfs,
                         warns,
                         columns):
        """Get the name of a gene product from the files"""
        for qt,df,warn,col in zip(query_types,dfs,warns,columns):
            if qt in product_type:
                if product not in df.index:
                    warn.append(gene_id)
                    return None
                product_name = df.loc[product][col]
                return product_name if product_name else None
        return product if product else None

    def _read_product_type(self,
                           gene_id,
                           product_types):
        """Get the type of a gene product from the files"""
        product_type = product_types[gene_id] \
                        if gene_id in product_types else 'gene'
        if 'MONOMER' in product_type:
            return 'CDS'
        else:
            return product_type

    def update_genbank_from_files(self):
        """ Complements GenBank file from optional files."""
        if self.is_reference:
            return
        contigs = self.contigs
        gene_sequences = self.gene_sequences
        gene_dictionary = self.gene_dictionary
        RNA_df = self.RNA_df
        complexes_df = self.complexes_df
        product_types = self.product_types
        all_genes_in_gb = self.all_genes_in_gb

        warn_rnas = []
        warn_proteins = []
        warn_position = []
        warn_sequence = []

        # Add new genes
        for gene_name,row in tqdm.tqdm(gene_dictionary.iterrows(),
                           'Updating Genbank file with optional files...',
                           bar_format = bar_format,
                           total=gene_dictionary.shape[0]):
            gene_id = row['Accession-1']

            if gene_id not in all_genes_in_gb:
                product = row['Product'].split(' // ')[0]
                ### Try to get product type from gene id of type LOCUST_TAG-RNA
                product_type = self._read_product_type(gene_id,
                                                       product_types)
                ### Retrieve values to sync with genbank
                product_name = \
                    self._get_product_name_if_present(gene_id,
                                                 product,
                                                 product_type,
                                                 ['RNA','MONOMER'],
                                                 [RNA_df,complexes_df],
                                                 [warn_rnas,warn_proteins],
                                                 ['Common-Name','name'])
                if product_name is None:
                    continue
                if not row['Left-End-Position'] or not row['Right-End-Position']:
                    warn_position.append(gene_id)
                    continue
                if gene_name not in gene_sequences:
                    warn_sequence.append(gene_id)
                    continue
                self._add_entry_to_genbank(
                     gene_id,
                     gene_name,
                     product_type,
                     product_name,
                     row,
                     contigs,
                     gene_sequences)

        # Ensure translation is in qualifiers
        warn_translation = []
        for record in self.contigs:
            for feature in record.features:
                if self.locus_tag not in feature.qualifiers:
                    continue
                if feature.type != "CDS":
                    continue
                if "translation" in feature.qualifiers:
                    continue
                warn_translation.append(feature.qualifiers[self.locus_tag][0])
                seq = feature.extract(record).seq
                feature.qualifiers["translation"] = [seq.translate(self.transl_table)]
        with open(self.directory + 'genome_modified.gb', 'w') as outfile:
            for contig in self.contigs:
                Bio.SeqIO.write(contig, outfile, 'genbank')
        # Warnings
        if warn_rnas:
            self.curation_notes['org.update_genbank_from_files'].append({
                                'msg':'Some genes were identified as RNA from their locus_tags, but they are not present in RNAs.txt',
                                'triggered_by':warn_rnas,
                                'importance':'medium',
                                'to_do':'Check whether you should add these genes to RNAs.txt or fix its product value in genes.txt'
            })
        if warn_proteins:
            self.curation_notes['org.update_genbank_from_files'].append({
                                'msg':'Some genes were identified as CDS from their locus_tags, but they are not present in proteins.txt',
                                'triggered_by':warn_proteins,
                                'importance':'medium',
                                'to_do':'Check whether you should add these genes to proteins.txt or fix its product value in genes.txt'
            })
        if warn_position:
            self.curation_notes['org.update_genbank_from_files'].append({
                                'msg':'Could not add some genes in genes.txt to genbank.gb since they lack position information',
                                'triggered_by':warn_position,
                                'importance':'medium',
                                'to_do':'Fill in position information in genes.txt'
            })
        if warn_sequence:
            self.curation_notes['org.update_genbank_from_files'].append({
                                'msg':'Could not add some genes in genes.txt to genbank.gb since they lack sequence information. Are your BioCyc files from the same database version?',
                                'triggered_by':warn_sequence,
                                'importance':'medium',
                                'to_do':'Add gene sequence in sequences.fasta. Check whether you downloaded the database files from the same BioCyc version.'
            })
        if warn_translation:
            self.curation_notes['org.update_genbank_from_files'].append({
                                'msg':'Some feature in genbank are CDS but have no translation qualifier. Translated sequences from Biopython were filled in instead',
                                'triggered_by':warn_translation,
                                'importance':'high',
                                'to_do':'Check whether the genbank was downloaded or constructed correctly.'
            })

    def _create_complexes_entry(self,
                                row,
                                genes,
                                stoich):
        """Create an entry to append to the complexes dataframe"""
        return {"name" : row["Common-Name"],
                "genes" : " AND ".join(
                        [
                            self.gene_dictionary["Accession-1"][g] + "({})".format(stoich)
                            for g in genes
                        ]
                    ),
                "source" : "BioCyc"}

    def generate_complexes_df(self):
        """
        Creates a DataFrame containing complex composition
        information from the provided files.
        """
        proteins_df = self.proteins_df

        if proteins_df.empty:
            return pandas.DataFrame(
                columns = [
                    'complex',
                    'name',
                    'genes',
                    'source'
                ]
            ).set_index('complex')

        gene_dictionary = self.gene_dictionary
        complexes = {}
        warn_proteins = []
        for p, row in tqdm.tqdm(proteins_df.iterrows(),
                           'Generating complexes dataframe from optional proteins file...',
                           bar_format = bar_format,
                           total=proteins_df.shape[0]):
            stoich = "" if "dimer" not in str(row["Common-Name"]) else "2"
            genes = row["Genes of polypeptide, complex, or RNA"]
            if not genes:
                warn_proteins.append(p)
                continue
            genes = [
                g for g in genes.split(" // ") if g in gene_dictionary["Accession-1"]
            ]

            complexes[p] = self._create_complexes_entry(row,
                                                        genes,
                                                        stoich)

        complexes_df = pandas.DataFrame.from_dict(complexes).T[["name", "genes", "source"]]
        complexes_df.index.name = "complex"
        # Warnings
        if warn_proteins:
            self.curation_notes['org.generate_complexes_df'].append({
                        'msg':'Some proteins have no genes in proteins.txt',
                        'triggered_by':warn_proteins,
                        'importance':'medium',
                        'to_do':'Fill genes in proteins.txt'})
        return complexes_df.fillna({"name": ""})

    def read_optional_file(self,filetype,filename,columns):
        """ Method for reading an optional file."""
        if os.path.isfile(filename):
            file = pandas.read_csv(filename, sep="\t",index_col=0)
        else:
            self.curation_notes['org.read_optional_file'].append({
                            'msg':'No {} file was found. Initializing an empty one.'.format(filetype),
                            'importance':'high',
                            'to_do':'Download {} from BioCyc if available'.format(filetype)})
            file = pandas.DataFrame(columns=columns).set_index(columns[0],inplace=False)
        return file.fillna('')

    def read_gene_dictionary(self,filename):
        """ Loads the genes file."""
        gene_dictionary = self.read_optional_file(
            'genes',
            filename,
            columns=[
                'Gene Name',
                'Accession-1',
                'Left-End-Position',
                'Right-End-Position',
                'Product'
            ]).reset_index().set_index('Gene Name')
        gene_dictionary['replicon'] = ''
        warn_genes = []
        if not self.is_reference:
            warn_start = list(gene_dictionary[gene_dictionary['Left-End-Position'].isna()].index)
            warn_end = list(gene_dictionary[gene_dictionary['Right-End-Position'].isna()].index)
            if warn_start:
                self.curation_notes['org.read_gene_dictionary'].append({
                            'msg':'Some genes are missing start positions in genes.txt',
                            'triggered_by':warn_start,
                            'importance':'medium',
                            'to_do':'Complete start positions in genes.txt if those genes are important.'})
            if warn_end:
                self.curation_notes['org.read_gene_dictionary'].append({
                            'msg':'Some genes are missing end positions in genes.txt',
                            'triggered_by':warn_end,
                            'importance':'medium',
                            'to_do':'Complete end positions in genes.txt if those genes are important.'})

            for g, row in gene_dictionary.iterrows():
                if not row["Accession-1"] or isinstance(row["Accession-1"],float):
                    gene_dictionary.at[g, "Accession-1"] = g
                    warn_genes.append(g)

        if warn_genes:
            self.curation_notes['org.read_gene_dictionary'].append({
                        'msg':'Some genes are missing Accession-1 IDs in genes.txt',
                        'triggered_by':warn_genes,
                        'importance':'medium',
                        'to_do':'Complete Accession-1 IDs in genes.txt if those genes are important.'})
        return gene_dictionary

    def read_proteins_df(self,filename):
        """ Loads the proteins file."""
        return self.read_optional_file(
            'proteins',
            filename,
            columns = [
                'Proteins',
                'Common-Name',
                'Genes of polypeptide, complex, or RNA',
                'Locations'
            ]
        )

    def read_gene_sequences(self,filename):
        """ Loads the gene sequences file."""
        if os.path.isfile(filename):
            d = {}
            for i in Bio.SeqIO.parse(filename,'fasta'):
                for g in i.id.split('|'):
                    d[g] = i
            return d
        return {}

    def read_RNA_df(self,filename):
        """ Loads the RNAs file."""
        return self.read_optional_file(
            'RNAs',
            filename,
            columns = [
                '(All-tRNAs RNAs Misc-RNAs rRNAs)',
                'Common-Name',
                'Gene'
            ]
        )

    def read_TU_df(self,filename):
        """ Loads the TUs file."""
        return self.read_optional_file(
            'TUs',
            filename,
            columns = [
                'Transcription-Units',
                'Genes of transcription unit',
                'Direction'
            ]
        )

    def check_gene_overlap(self):
        """ Assesses gene identifier overlap between files."""
        if self.is_reference:
            return

        def get_severity(o):
            if o < 50 :
                return 'critical'
            elif o < 60 :
                return 'high'
            elif o < 70 :
                return 'medium'
            elif o < 80 :
                return 'low'
            else:
                return 0

        m_model_genes = set([g.id for g in self.m_model.genes])
        file_genes = set(self.gene_dictionary['Accession-1'].values)
        all_genes_in_gb = []
        transl_table = []
        warn_table = []
        for record in self.contigs:
            for feature in record.features:
                if self.locus_tag not in feature.qualifiers:
                    continue
                all_genes_in_gb.append(feature.qualifiers[self.locus_tag][0])
                transl_table+=(feature.qualifiers.get('transl_table',[None]))
#         self.all_genes_in_gb = all_genes_in_gb
        transl_table = set(i for i in set(transl_table) if i is not None)
        if len(transl_table) > 1:
            warn_table = transl_table
        elif not transl_table:
            transl_table = ['11']
        self.transl_table = list(transl_table)[0]

        genbank_genes = set(all_genes_in_gb)

        # Overlaps
        file_overlap = int((len(file_genes & m_model_genes) / len(m_model_genes))*100)
        gb_overlap = int((len(genbank_genes & m_model_genes) / len(m_model_genes))*100)

        logging.warning('Gene overlap between M-model and Genbank : {}%'.format(gb_overlap))
        logging.warning('Gene overlap between M-model and optional files : {}%'.format(file_overlap))

        if gb_overlap < 1:
            raise ValueError('Overlap of M-model genes with genbank is too low ({}%)'.format(gb_overlap))

        fs = get_severity(file_overlap)
        gs = get_severity(gb_overlap)

        if fs:
            self.curation_notes['org.check_gene_overlap'].append({
                'msg':'M-model has a {} gene overlap with optional files (BioCyc)'.format(file_overlap),
                'importance':fs,
                'to_do':'Check whether optional files where downloaded correctly.'})
        if gs:
            self.curation_notes['org.check_gene_overlap'].append({
                'msg':'M-model has a {} gene overlap with Genbank'.format(gb_overlap),
                'importance':gs,
                'to_do':'Check whether genbank was downloaded correctly.'})
        if warn_table:
            self.curation_notes['org.check_gene_overlap'].append({
                'msg':'Provided GenBank file contains more than one translation table. Is this correct?',
                'triggered_by':warn_table,
                'importance':'high',
                'to_do':'Check if translation tables are correct.'})

    def update_ribosome_stoich(self):
        """ Updated ribosome composition from files."""
        if self.is_reference:
            return
        complexes_df = self.complexes_df
        protein_mod = self.protein_mod.reset_index().set_index('Core_enzyme')
        ribo_df = complexes_df.loc[
            complexes_df["name"].str.contains("ribosomal.*(?:subunit)?.* protein", regex=True)
        ]
        self.ribosomal_proteins = ribo_df
        ribosome_stoich = self.ribosome_stoich
        ribo_30S = ribosome_stoich["30_S_assembly"]["stoich"]
        if [i for i in ribo_30S if "generic" not in i]:
            update_30S = False
        else:
            update_30S = True
        ribo_50S = ribosome_stoich["50_S_assembly"]["stoich"]
        if [i for i in ribo_50S if "generic" not in i]:
            update_50S = False
        else:
            update_50S = True
        if not(update_30S or update_50S):
            # Only update if it has not been user-defined
            return
        trigger_factor = list(complexes_df[complexes_df['name'].str.contains('[T,t]rigger factor',regex=True)].index)
        if trigger_factor:
            ribo_50S[trigger_factor[0]] = 1
        warn_proteins = []

        for p, row in tqdm.tqdm(ribo_df.iterrows(),
                           'Gathering ribosome stoichiometry...',
                           bar_format = bar_format,
                           total=ribo_df.shape[0]):
            p_mod_list = []
            if p in protein_mod.index:
                p_mod_list = protein_mod.loc[[p]]['Modified_enzyme'].values
            if re.search("30S|small.*subunit",row["name"],re.IGNORECASE) and update_30S:
                if set(p_mod_list) & set(ribo_30S.keys()):
                    continue
                ribo_30S[p] = 1
            elif re.search("50S|large.*subunit",row["name"],re.IGNORECASE) and update_50S:
                if set(p_mod_list) & set(ribo_50S.keys()):
                    continue
                ribo_50S[p] = 1
            else:
                if set(p_mod_list) & set(ribo_50S.keys()):
                    continue
                ribo_50S[p] = 1 # Add it to stoichiometry but warn it might not be a good mapping
                warn_proteins.append(p)
        if warn_proteins:
            self.curation_notes['org.update_ribosome_stoich'].append({
                'msg':'Some ribosomal proteins do not contain subunit information (30S, 50S). Check whether they are ribosomal proteins!',
                'triggered_by':warn_proteins,
                'importance':'high',
                'to_do':'Curate them in ribosomal_proteins.txt'})

    def _add_entry_to_gene_dictionary(self,
                                gene_dictionary,
                                gene_id,
                                feature,
                                left_end,
                                right_end):
        """Append one entry to the gene information dataframe"""
        logging.warning("Adding {} to genes from genbank".format(gene_id))
        feature_type = feature.type
        if feature_type == 'CDS':
            feature_type = 'MONOMER'
        tmp = {gene_id: {
                        "Accession-1": gene_id,
                        "Left-End-Position": left_end,
                        "Right-End-Position": right_end,
                        "Product": "{}-{}".format(gene_id,feature_type)
                }}
        return self._add_entry_to_df(gene_dictionary, tmp)

    def _add_entry_to_complexes_or_rna(self,
                                       complexes_df,
                                       RNA_df,
                                       gene_name,
                                       gene_id,
                                       feature,
                                      ):
        """Append one entry to the Complexes of RNA dataframe according to its type"""
        name_annotation = feature.qualifiers["product"][0] if 'product' in feature.qualifiers \
                else gene_name
        if feature.type == 'CDS':
            product = gene_name + '-MONOMER'
            if not complexes_df["genes"].str.contains(gene_id).any():
                complexes_df = \
                    self._add_entry_to_complexes(gene_id,
                                                 name_annotation,
                                                 product,
                                                 complexes_df,
                                                 "GenBank")
        else: # It's not CDS, but an RNA
            product = "{}-{}".format(gene_name,feature.type)
            if not RNA_df["Gene"].str.contains(gene_name).any():
                RNA_df = \
                    self._add_entry_to_rna(gene_id,
                                           name_annotation,
                                           product,
                                           RNA_df,
                                           "GenBank")
        return complexes_df,RNA_df,product

    def _add_entries_to_optional_files(self,
                                       gene_dictionary,
                                       complexes_df,
                                       RNA_df,
                                       feature,
                                       record,
                                       product_types):
        """Add entries to the optional files from the GenBank"""

#         gene_id = feature.qualifiers[self.locus_tag][0]
        gene_id = self._get_feature_locus_tag(feature)
        if ';' in gene_id:
            gene_id = gene_id.split(';')[0]
        left_end = min([i.start for i in feature.location.parts])
        right_end = max([i.end for i in feature.location.parts])
        if not gene_dictionary["Accession-1"].eq(gene_id).any():
            gene_dictionary = \
                self._add_entry_to_gene_dictionary(
                        gene_dictionary,
                        gene_id,
                        feature,
                        left_end,
                        right_end)

        gene_names = gene_dictionary[gene_dictionary["Accession-1"].eq(gene_id)].index
        warn_genes = []
        for gene_name in gene_names:
            complexes_df,RNA_df,product = \
                self._add_entry_to_complexes_or_rna(
                                   complexes_df,
                                   RNA_df,
                                   gene_name,
                                   gene_id,
                                   feature,
                                  )
            gene_dictionary.at[gene_name,'Product'] = product # Ensuring product is the same.
            gene_dictionary.at[gene_name,"Left-End-Position"] = left_end
            gene_dictionary.at[gene_name,"Right-End-Position"] = right_end
            gene_dictionary.at[gene_name,"replicon"] = record.id

            # Update product types
            gid,product,product_type = \
                self._get_product_type(
                         gene_name,
                         gene_dictionary = gene_dictionary,
                         complexes_df = complexes_df,
                         RNA_df = RNA_df,
                         warn_genes = warn_genes)
            if product is None:
                warn_genes.append(gid)
                continue
            if ' ' in product or ('RNA' not in product and 'MONOMER' not in product):
                product = \
                    self._correct_product(
                        gene_name,
                        product_type)

            product_types[gene_id] = product_type

        return gene_dictionary,complexes_df,RNA_df

    def update_complexes_genes_with_genbank(self):
        """ Complements complexes and genes with genome"""
        if self.is_reference:
            return

        # In some genbanks, CDS are duplicated with gene features. See staph or pputida
        complexes_df = self.complexes_df
        gene_dictionary = self.gene_dictionary
        RNA_df = self.RNA_df
        product_types = self.product_types
        warn_locus = []
        for record in tqdm.tqdm(self.contigs,
                           'Syncing optional files with genbank contigs...',
                           bar_format = bar_format):
            for feature in record.features:
                if feature.type not in element_types:
                    continue
                if self._get_feature_locus_tag(feature) is None:
                    continue
                gene_dictionary,complexes_df,RNA_df = \
                    self._add_entries_to_optional_files(
                                       gene_dictionary,
                                       complexes_df,
                                       RNA_df,
                                       feature,
                                       record,
                                       product_types)
        self.complexes_df = complexes_df
        gene_dictionary.index.name = "Gene Name"
        self.gene_dictionary = gene_dictionary
        self.RNA_df = RNA_df

        if warn_locus:
            self.curation_notes['org.update_complexes_genes_with_genbank'].append({
                'msg':'Some genbank features do not have locus_tag.',
                'triggered_by':warn_locus,
                'importance':'high',
                'to_do':'Check whether these features are necessary, and correct their locus_tag. If they have been completed from other provided files, ignore.'})

    def purge_genes_in_model(self):
        """ Purges problematic genes in the M-model"""
        m_model = self.m_model
        gene_dictionary = self.gene_dictionary
        RNA_df = self.RNA_df
        gene_list = []
        wrong_assoc = []
        for g in tqdm.tqdm(m_model.genes,
                           'Purging M-model genes...',
                           bar_format = bar_format):
            if g.id not in gene_dictionary['Accession-1'].values:
                gene_list.append(g)
            else:
                product = gene_dictionary[self.gene_dictionary['Accession-1'].eq(g.id)]['Product'].values[0]
                if product in RNA_df.index:
                    wrong_assoc.append(g)

        self.skip_genes = [g.id for g in gene_list + wrong_assoc]
#         cobra.manipulation.delete.remove_genes(m_model,
#                      gene_list + wrong_assoc,
#                      remove_reactions=False)
        # Warnings
        if gene_list:
            self.curation_notes['org.purge_genes_in_model'].append({
                'msg':'Some genes in M-model were not found in genes.txt or genome.gb. These genes will be skipped in reconstruction.',
                'triggered_by':[g.id for g in gene_list],
                'importance':'high',
                'to_do':'Confirm the gene is correct in the m_model. If so, add it to genes.txt'})
        if wrong_assoc:
            self.curation_notes['org.purge_genes_in_model'].append({
                'msg':'Some genes in M-model are RNAs. These genes will be skipped in reconstruction.',
                'triggered_by':[g.id for g in wrong_assoc],
                'importance':'high',
                'to_do':'Confirm the gene is correct in the m_model. If so, then annotation from GenBank or BioCyc marked them as a different type'})

    def _get_ligases_from_regex(self,
                                complexes_df):
        """Call ligases from provided files using regex"""
        return self._get_slice_from_regex(
            complexes_df,
            "[-]{,2}tRNA (?:synthetase|ligase)(?!.*subunit.*)")

    def _get_ligases_subunits_from_regex(self,
                                complexes_df):
        """Call ligases subunits from provided files using regex"""
        return self._get_slice_from_regex(
            complexes_df,
            "[-]{,2}tRNA (?:synthetase|ligase)(?=.*subunit.*)")
    def _extract_trna_string(self,
                             trna_string):
        """Srip the tRNA ligase information to call the amino acid"""
        t = re.findall(".*[-]{,2}tRNA (?:synthetase|ligase)",trna_string)
        return t[0] if t else None

    def _is_base_complex_in_list(self,cplx,lst):
        """Check if a complex (base) is already in a list of complexes"""
        return cplx in set(i.split('_mod_')[0] for i in lst)

    def _get_genes_of_cplx(self,cplx):
        """Get the genes associated with a complex from the files"""
        d = {}
        for i in self.complexes_df.loc[cplx]['genes'].split(' AND '):
            gene = re.findall('.*(?=\(\d*\))',i)[0]
            coeff = re.findall('(?<=\().*(?=\))',i)[0]
            d[gene] = coeff
        return d
#         return [re.findall('.*(?=\(\d*\))',i)[0] \
#                 for i in self.complexes_df.loc[cplx]['genes'].split(' AND ')]

    def get_trna_synthetase(self):
        """ Gets tRNA synthetases from files."""
        if self.is_reference:
            return

        def find_aminoacid(trna_string):
            """Call the aminoacid from a ligase description string"""
            trna_string = trna_string.lower()
            for aa, rx in coralme.builder.dictionaries.amino_acid_regex.items():
                if re.search(rx, trna_string):
                    return aa
            return None

        org_amino_acid_trna_synthetase = self.amino_acid_trna_synthetase
        manually_curated_aa = [k for k,v in org_amino_acid_trna_synthetase.items() if v]
        generic_dict = self.generic_dict
        complexes_df = self.complexes_df
        warn_generic = []
        d = collections.defaultdict(set)
        for k,v in org_amino_acid_trna_synthetase.copy().items():
            if isinstance(v,list):
                d[k] = set(v)
            elif isinstance(v,str):
                if not v:
                    d[k] = set()
                    continue
                if 'generic' in v:
                    if v not in generic_dict:
                        warn_generic.append(v)
                        d[k] = set()
                        continue
#                     d[k] = set(generic_dict[v]['enzymes'])
#                     continue
                d[k] = set([v])
        trna_ligases = self._get_ligases_from_regex(complexes_df).to_dict()['name']
        for cplx, trna_string in trna_ligases.items():
            aa = find_aminoacid(trna_string)
            if aa is None:continue
            if aa in manually_curated_aa: continue
#             if aa not in d: d[aa] = set()
            if self._is_base_complex_in_list(cplx,d[aa]): continue
            d[aa].add(cplx)
        trna_ligases_from_subunits = self._get_ligases_subunits_from_regex(complexes_df).to_dict()['name']
#         new_cplxs = {k:dict() for k in d.copy()}
        new_cplxs = collections.defaultdict(dict)
        for cplx,trna_string in trna_ligases_from_subunits.items():
            trna_string = self._extract_trna_string(trna_string)
            aa = find_aminoacid(trna_string)
            if aa is None:continue
#             if aa not in d: d[aa] = set()
            if d[aa]: continue
            cplx_genes = self._get_genes_of_cplx(cplx)
            for k,v in cplx_genes.items():
#                 if aa not in new_cplxs: new_cplxs[aa] = dict()
                new_cplxs[aa][k] = v
#             new_cplxs[aa].add(cplx)

        for k,v in new_cplxs.items():
            if not v: continue
            cplx_id = "CPLX-tRNA-{}-LIGASE".format(k.upper()[:3])
            complexes_df = self._add_entry_to_complexes(
                               v,
                               "tRNA-{} ligase".format(k[0].upper() + k[1:3]),
                               cplx_id,
                               complexes_df,
                               "Inferred from subunits")
            d[k] = set([cplx_id])

        warn_ligases = []
        for aa,c_set in d.items():
            c_list = list(c_set)
            if len(c_list) == 1:
                d[aa] = c_list[0]
            elif len(c_list) > 1:
                generic = 'generic_{}_ligase'.format(aa.split('__L_c')[0])
                generic_dict[generic] = {'enzymes':c_list}
                d[aa] = generic
            else:
                d[aa] = 'CPLX_dummy'
                warn_ligases.append(aa)
        self.amino_acid_trna_synthetase = dict(d)
        self.complexes_df = complexes_df

        # Warnings
        if warn_ligases:
            self.curation_notes['org.get_trna_synthetase'].append({
                'msg':'No tRNA ligases were found for some amino acids. Assigned CPLX_dummy.',
                'triggered_by':warn_ligases,
                'importance':'high',
                'to_do':'Check whether your organism should have a ligase for these amino acids, or if you need to add a reaction to get it (e.g. tRNA amidotransferases)'})

        if warn_generic:
            self.curation_notes['org.get_trna_synthetase'].append({
                'msg':'A generic tRNA ligase was defined in amino_acid_trna_synthetase, but it is not defined in generic_dict.',
                'triggered_by':warn_generic,
                'importance':'high',
                'to_do':'Fix the definition in generic_dict'})

    def get_peptide_release_factors(self):
        """ Gets peptide release factors from files."""
        if self.is_reference:
            return

        proteins_df = self.proteins_df["Common-Name"].dropna()
        peptide_release_factors = self.peptide_release_factors
        generics = self.generic_dict
        rf = proteins_df[proteins_df.str.contains("peptide.*release.*factor")]
        if rf.size:
            if not peptide_release_factors["UAG"]['enzyme'] and rf.str.contains("1").any():
                peptide_release_factors["UAG"]['enzyme'] = rf[rf.str.contains("1")].index[0]
                generics["generic_RF"]['enzymes'].append(peptide_release_factors["UAG"]['enzyme'])
            if not peptide_release_factors["UGA"]['enzyme'] and rf.str.contains("2").any():
                peptide_release_factors["UGA"]['enzyme'] = rf[rf.str.contains("2")].index[0]
                generics["generic_RF"]['enzymes'].append(peptide_release_factors["UGA"]['enzyme'])

    def get_nonmetabolic(self):
        """ Gets nonmetabolic metabolites in M-model."""
        m_model = self.m_model
        queries = ['ACP','trna']
        for m in m_model.metabolites.query('|'.join(queries)):
            if m.id in self.me_mets.index:
                continue
            tmp = {m.id:{
                'me_id':'',
                'name':'',
                'formula':'',
                'compartment':'',
                'type':'CURATE'
            }}
            self.me_mets = self._add_entry_to_df(self.me_mets,tmp)
        return None

    def gb_to_faa(self, org_id, outdir = False, element_types = {"CDS"}):
        """ Generates a protein FASTA from genome for BLAST"""

        # Create FASTA file with AA sequences for BLAST
        contigs = self.contigs

        if not outdir:
            outdir = self.blast_directory

        #outdir += "blast_files_and_results/"
        FASTA_file = outdir + "{}.faa".format(org_id)
#         FASTA_file = "{}.faa".format(org_id)

        file = open(FASTA_file, "w")
        for contig in tqdm.tqdm(contigs,
                           'Converting Genbank contigs to FASTA for BLAST...',
                           bar_format = bar_format):
            for feature in contig.features:
                if feature.type not in element_types \
                    or "translation" not in feature.qualifiers \
                    or self.locus_tag not in feature.qualifiers:
                    continue
                file.write(
                    ">{}\n".format(feature.qualifiers[self.locus_tag][0])
                )  # Some way to identify which qualifier meets regular expression?
                file.write("{}\n".format(feature.qualifiers["translation"][0]))
        file.close()

    def _process_sigma_name(self,name, row):
        """Generate a sigma factor identifier from its description"""
        name = name.split("RNA polymerase")[-1]
        replace_list = ["sigma", "factor", "sup"]
        for r in replace_list:
            name = name.replace(r, "")
        name = "".join(re.findall("[a-zA-Z0-9]{1,}", name))
        if not name:
            name = "_".join(row["genes"])
        return "RNAP_" + name

    def get_sigma_factors(self):
        """ Gets sigma factors from files."""
        complexes_df = self.complexes_df

        sigma_df = complexes_df.loc[
            complexes_df["name"].str.contains(
                "RNA polymerase.*sigma.*[factor.*]?.*|sigma.*[factor.*]?.*RNA polymerase", regex=True
            )
        ]
        if not sigma_df.shape[0]:
            self.curation_notes['org.get_sigma_factors'].append({
                'msg':"No sigma factors could be identified from proteins.txt",
                'importance':'critical',
                'to_do':'Manually define sigmas in sigma_factors.txt'})
            random_cplx = random.choice(complexes_df.index)
            sigma_df = complexes_df.loc[[random_cplx]]
        ## Get sigmas automatically
        # Find RpoD to add as default sigma
        for s, row in tqdm.tqdm(sigma_df.iterrows(),
                           'Getting sigma factors...',
                           bar_format = bar_format,
                           total=sigma_df.shape[0]):

            tmp = {
                s : {
                    "complex" : self._process_sigma_name(s, row),
                    "genes" : row["genes"],
                    "name" : row["name"]
                }
            }
            self.sigmas = self._add_entry_to_df(self.sigmas,tmp)

    def get_rpod(self):
        """ Gets RpoD from files."""
        sigma_df = self.sigmas
        rpod = sigma_df[sigma_df["name"].str.contains("RpoD")].index.to_list()
        if not rpod:
            rpod_re = "|".join(["70", "sigma-A", "sigA", "SigA","Sigma-A"])
            rpod = sigma_df[sigma_df["name"].str.contains(rpod_re)].index.to_list()
        if rpod:
            rpod = rpod[0]
            # Warnings
            self.curation_notes['org.get_rpod'].append({
                'msg':"{} was identified as RpoD. If this is not true, define RpoD!".format(rpod),
                'importance':'high',
                'to_do':'Check whether you need to correct RpoD by running me_builder.org.rpod = correct_rpod'})
        else:
            rpod = random.choice(sigma_df.index)
            # Warnings
            self.curation_notes['org.get_rpod'].append({
                'msg':"RpoD randomly assigned to {}".format(rpod),
                'importance':'critical',
                'to_do':'genome.gb does not have a valid annotation for RpoD. A random identified sigma factor in me_builder.org.sigmas was set as RpoD so that the builder can continue running. Set the correct RpoD by running me_builder.org.rpod = correct_rpod'})
        self.rpod = rpod

    def _get_slice_from_regex(self,
                        df,
                        regex,
                       ):
        """The a slice of a dataframe from a regex"""
        return df[df["name"].str.contains(regex,regex=True)]

    def _get_complex_from_regex(self,
                               complexes_df,
                               cplx_regex,
                               subunit_regex=None):
        """Get complex as one entry from complexes"""
        cplx = self._get_slice_from_regex(complexes_df,cplx_regex)
        if cplx.shape[0] == 1:
            return cplx.iloc[[0],:],'cplx'
        if subunit_regex is None:
            return None,None
        # Get complex as composed from subunits
        subunits = self._get_slice_from_regex(complexes_df,subunit_regex)
        if subunits.empty:
            return None,None
        if subunits.shape[0] == 1:
            return subunits,'cplx'
        return subunits,'subunits'

    def _get_rna_polymerase_from_regex(self,
                                        complexes_df):
        """Call the RNAP from regex"""
        cplx,flag = self._get_complex_from_regex(
            complexes_df,
            "(?:RNA polymerase.*core enzyme|DNA.*directed.*RNA polymerase)(?!.*subunit.*|.*chain.*)",
            subunit_regex = "(?:RNA polymerase.*core enzyme|DNA.*directed.*RNA polymerase)(?=.*subunit.*|.*chain.*)")

        if cplx is not None:
            return cplx,flag

        cplx,flag = self._get_complex_from_regex(
                complexes_df,
                "(?:^RNA polymerase$)",
                subunit_regex = "(?:^RNA polymerase)(?=.*subunit.*|.*chain.*)")
        return cplx,flag

    def _add_rna_polymerase_to_complexes(self,
                                        complexes_df,
                                        RNAP_genes):
        """Add the RNAP entry to the complexes dataframe"""
        tmp = pandas.DataFrame.from_dict({
            "RNAP-CPLX": {
                "name": "DNA-directed RNA polymerase",
                "genes": " AND ".join(["{}()".format(g) for g in RNAP_genes]),
                "source": "GenBank",
                }}).T
        return pandas.concat([complexes_df, tmp], axis = 0, join = 'outer')

    def _is_beta_prime_in_RNAP(self,RNAP,complexes_df):
        """Checks if beta prime subunit is in RNAP"""
        genes = [re.findall('.*(?=\(\d*\))',i)[0] for i in complexes_df.loc[RNAP]['genes'].split(' AND ')]
        df = complexes_df[complexes_df['genes'].str.contains('|'.join(genes))]
        return df['name'].str.contains("beta(?:\'|.*prime)|rpoc|RNA polymerase.*(?:subunit|chain).*beta",regex=True,case=False).any()

    def get_rna_polymerase(self, force_RNAP_as=""):
        """Call the RNAP from files"""
        # TODO: Allow user to define RNAP, skip inferring?
        complexes_df = self.complexes_df
        protein_mod = self.protein_mod
        RNAP = ""
        if force_RNAP_as:
            RNAP = force_RNAP_as
        else:
            RNAP,flag = self._get_rna_polymerase_from_regex(complexes_df)
            if RNAP is None:
                RNAP = random.choice(complexes_df.index)
                self.curation_notes['org.get_rna_polymerase'].append({
                    'msg':"Could not identify RNA polymerase".format(RNAP),
                    'importance':'critical',
                    'to_do':'Find correct RNAP complex and run me_builder.org.get_rna_polymerase(force_RNAP_as=correct_RNAP)'})
            elif flag == 'cplx':
                RNAP = RNAP.index[0]
                # Warnings
                self.curation_notes['org.get_rna_polymerase'].append({
                    'msg':"{} was identified as RNA polymerase".format(RNAP),
                    'importance':'high',
                    'to_do':'Check whether you need to correct RNAP by running me_builder.org.get_rna_polymerase(force_RNAP_as=correct_RNAP)'})
            elif flag == 'subunits':
                RNAP_genes = [g.split("-MONOMER")[0] for g in RNAP.index if "-MONOMER" in g]
                RNAP_genes = [self.gene_dictionary.loc[g]['Accession-1'] for g in RNAP_genes]
                RNAP = 'RNAP-CPLX'
                complexes_df = self._add_rna_polymerase_to_complexes(complexes_df,
                                                                    RNAP_genes)
                self.curation_notes['org.get_rna_polymerase'].append({
                    'msg':"RNAP was identified with subunits {}".format(
                        ", ".join(RNAP_genes)
                    ),
                    'importance':'medium',
                    'to_do':'Check whether the correct proteins were called as subunits of RNAP. If not find correct RNAP complex and run me_builder.org.get_rna_polymerase(force_RNAP_as=correct_RNAP)'})

        # Identify if beta prime in RNAP, if so, add zn2 and mg2. https://pubmed.ncbi.nlm.nih.gov/15351641/
        if self._is_beta_prime_in_RNAP(RNAP,complexes_df):
            RNAP_mod = RNAP + '_mod_zn2(1)_mod_mg2(2)'
            protein_mod = \
                self._add_entry_to_protein_mod(protein_mod,
                                             RNAP_mod,
                                             RNAP,
                                             "zn2(1) AND mg2(2)",
                                             "RNA_Polymerase")
            RNAP = RNAP_mod
        self.RNAP = RNAP
        self.complexes_df = complexes_df
        self.protein_mod = protein_mod
        self.sigma_factor_complex_to_rna_polymerase_dict = self.sigmas[
            "complex"
        ].to_dict()
        self.rna_polymerase_id_by_sigma_factor = {}
        for k, v in self.sigma_factor_complex_to_rna_polymerase_dict.items():
            self.rna_polymerase_id_by_sigma_factor[v] = {
                "sigma_factor": k,
                "polymerase": self.RNAP,
            }
        self.rna_polymerases = list(self.rna_polymerase_id_by_sigma_factor.keys())

    def get_TU_genes(self):
        """ Gets TU composition from files."""
        TUs = self.TUs
        gene_dictionary = self.gene_dictionary
        genes_to_TU = {}
        TU_to_genes = {}
        for tu, row in tqdm.tqdm(TUs.iterrows(),
                           'Getting TU-gene associations from optional TUs file...',
                           bar_format = bar_format,
                           total=TUs.shape[0]):
            genes = row["Genes of transcription unit"]
            if not genes:
                continue
            for g in genes.split(" // "):
                if g not in gene_dictionary.index:
                    continue
                genes_to_TU[gene_dictionary["Accession-1"][g]] = tu

                if tu not in TU_to_genes:
                    TU_to_genes[tu] = []
                TU_to_genes[tu].append(g)
        self.genes_to_TU = genes_to_TU
        self.TU_to_genes = TU_to_genes

    def get_TU_df(self):
        """ Generates TUs_from_biocyc."""
        TUs = self.TUs
        gene_dictionary = self.gene_dictionary
        rpod = self.rpod
        rho_independent = self.rho_independent
        TU_dict = {}
        warn_genes = []
        warn_tus = []
        if TUs.empty:
            return pandas.DataFrame(
                columns = [
                    'TU_id',
                    'replicon',
                    'genes',
                    'start',
                    'stop',
                    'tss',
                    'strand',
                    'rho_dependent',
                    'rnapol'
                ]
            ).set_index('TU_id')
        for tu, row in tqdm.tqdm(TUs.iterrows(),
                           'Processing optional TUs file...',
                           bar_format = bar_format,
                           total=TUs.shape[0]):
            sites = []
            start = []
            stop = []
            genes = []
            replicons = []
            for g in row["Genes of transcription unit"].split(" // "):
                if g not in gene_dictionary.index:
                    warn_genes.append(g)
                    continue
                genes.append(gene_dictionary["Accession-1"][g])
                sites.append(int(gene_dictionary["Left-End-Position"][g]))
                sites.append(int(gene_dictionary["Right-End-Position"][g]))
#                 start.append(int(gene_dictionary["Left-End-Position"][g]))
#                 stop.append(int(gene_dictionary["Right-End-Position"][g]))
                replicons.append(gene_dictionary["replicon"][g])
            if not genes:
                warn_tus.append(tu)
                continue
            sigma = rpod  # Default RpoD
            tu_name = "{}_from_{}".format(tu, sigma)
            TU_dict[tu_name] = {}
            TU_dict[tu_name]["genes"] =  ','.join(genes)
            TU_dict[tu_name]["rho_dependent"] = False if tu in rho_independent else True
            TU_dict[tu_name]["rnapol"] = sigma
            TU_dict[tu_name]["tss"] = None
            TU_dict[tu_name]["strand"] = row["Direction"] if row["Direction"] else '+'
            TU_dict[tu_name]["start"] = int(min(sites))+1
#             TU_dict[tu_name]["start"] = ','.join([ str(x+1) for x in start ])
            TU_dict[tu_name]["stop"] = int(max(sites))
#             TU_dict[tu_name]["stop"] = ','.join([ str(x) for x in stop ])
            TU_dict[tu_name]["replicon"] = ','.join(replicons) if set(replicons) != {''} else None
        df = pandas.DataFrame.from_dict(TU_dict).T[
            #["start", "stop", "tss", "strand", "rho_dependent", "rnapol","replicon"]
            ['replicon', 'genes', 'start', 'stop', 'tss', 'strand', 'rho_dependent', 'rnapol']
        ]
        df.index.name = "TU_id"

        # Warnings
        if warn_genes or warn_tus:
            if warn_genes:
                self.curation_notes['org.get_TU_df'].append({
                        'msg':"Some genes appear in TUs.txt but not in genes.txt",
                        'triggered_by':warn_genes,
                        'importance':'medium',
                        'to_do':'If those genes are supposed to be in the model, fill them in genes.txt'})
            if warn_tus:
                self.curation_notes['org.get_TU_df'].append({
                        'msg':"Some TUs are either empty (no genes in TUs.txt) or contain genes that are not in genes.txt",
                        'triggered_by':warn_tus,
                        'importance':'medium',
                        'to_do':'If those TUs contain genes that are supposed to be in the model, fill them in TUs.txt and genes.txt'})
        return df

    def _is_cytosol_in_locations(self,locs):
        """Checks if cytosol is in the locations list of a protein"""
        for i in locs:
            if 'CCI-CYTOSOL' in i:
                return True
        return False
    def _process_location_dict(self,
                               location,
                               location_interpreter):
        new_location = {}
        for k, v in location.items():
            if not v or isinstance(v, float):
                continue
            locations = []
            for loc in v.split(" // "):
                if loc not in location_interpreter.index:
                    continue
                locations.append(location_interpreter["interpretation"][loc])
            if self._is_cytosol_in_locations(locations):
                continue
            if locations:
                new_location[k] = locations[0]
        return new_location

    def _add_entry_to_protein_location(self,
                                       c,
                                       c_loc,
                                       gene_string,
                                       gene_dictionary,
                                      protein_location,
                                      gene_location):
        """Append one entry to the protein locations dataframe"""
        gene = re.findall('.*(?=\(\d*\))', gene_string)[0]
        if gene not in gene_dictionary.index:
            return protein_location
        gene = gene_dictionary.loc[[gene]]["Gene Name"]
        for gene_ in gene: # In case of duplicates
            if gene_ in gene_location:
                tmp = {
                    c: {
                        "Complex_compartment": c_loc,
                        "Protein": gene_string,
                        "Protein_compartment": gene_location[gene_],
                        "translocase_pathway": "s",
                        }}
                protein_location = self._add_entry_to_df(protein_location,tmp)
        return protein_location

    def get_protein_location(self):
        """ Gets protein location from files."""
        complexes_df = self.complexes_df
        proteins_df = self.proteins_df
        gene_dictionary = self.gene_dictionary
        location_interpreter = self.location_interpreter
        protein_location = self.protein_location
        gene_location = self._process_location_dict(
            proteins_df.set_index("Genes of polypeptide, complex, or RNA")["Locations"]
            .dropna()
            .to_dict(),
            location_interpreter,
        )
        cplx_location = self._process_location_dict(
            proteins_df["Locations"].dropna().to_dict(), location_interpreter
        )
        gene_dictionary = gene_dictionary.reset_index().set_index("Accession-1")
        for c, row in tqdm.tqdm(complexes_df.iterrows(),
                           'Adding protein location...',
                           bar_format = bar_format,
                           total=complexes_df.shape[0]):

            if c not in cplx_location:
                continue
            c_loc = cplx_location[c]
            for gene_string in complexes_df["genes"][c].split(' AND '):
                protein_location = self._add_entry_to_protein_location(
                                                    c,
                                                    c_loc,
                                                    gene_string,
                                                    gene_dictionary,
                                                    protein_location,
                                                    gene_location)
        self.protein_location = protein_location

        #TODO: Add warnings for empty location

    # TODO: New format of keffs file
    def get_reaction_keffs(self):
        """ Gets reaction Keffs from files."""
        if self.is_reference:
            return None
        # Keff estimator from https://pubs.acs.org/doi/10.1021/bi2002289
        m_model = self.m_model
        subsystem_classification = self.subsystem_classification
        enz_rxn_assoc_df = self.enz_rxn_assoc_df
        rxn_keff_dict = {}
        reaction_dirs = ["FWD", "REV"]
        keffs = {
            "central_CE": 79,
            "central_AFN": 18,
            "intermediate": 5.2,
            "secondary": 2.5,
            "other": 65,
        }
        for reaction, row in tqdm.tqdm(enz_rxn_assoc_df.iterrows(),
                           'Getting reaction Keffs...',
                           bar_format = bar_format,
                           total=enz_rxn_assoc_df.shape[0]):
            if reaction not in m_model.reactions:
                #TODO: Change this so that Keffs of new reactions in reaction_matrix can be estimated
                logging.warning('Tried setting Keffs for {}, but the reaction is not in model.'.format(reaction))
                continue
            r = m_model.reactions.get_by_id(reaction)
            subsystem = r.subsystem
            if subsystem in subsystem_classification:
                classification = subsystem_classification[subsystem]
            else:
                classification = "other"
            keff = keffs[classification]
            cplxs = row["Complexes"].split(" OR ")
            for c in cplxs:
                for d in reaction_dirs:
                    r = "{}".format(reaction)
                    rxn_keff_dict[r] = {}
                    rxn_keff_dict[r]["keff"] = keff
        self.reaction_median_keffs = pandas.DataFrame.from_dict(rxn_keff_dict).T
        self.reaction_median_keffs.index.name = "reaction"
        self.reaction_median_keffs.to_csv(self.directory + 'reaction_median_keffs.txt', sep='\t')
        return self.reaction_median_keffs['keff'].to_dict()

    def get_phospholipids(self):
        """ Gets phospholipids from M-model."""
        m_model = self.m_model
        return [
            str(m.id) for m in m_model.metabolites.query(re.compile("^pg[0-9]{2,3}_.$"))
        ]
    def get_lipids(self):
        """ Gets lipids from M-model."""
        m_model = self.m_model
        return [
            str(m.id) for m in m_model.metabolites.query(re.compile("^[a-z]*[0-9]{2,3}_.$"))
        ]

    def _get_feature_locus_tag(self,
                               feature):
        """Get the locus tag identifier from a Biopython feature"""
        lt = feature.qualifiers.get(self.locus_tag,None)
        if lt is not None:
            return lt[0]
        lt = feature.qualifiers.get('locus_tag',None)
        return lt[0] if lt is not None else None

    def _map_to_a_generic(self,
                          feature,
                          generic_dict):
        """Associate a Biopython feature to a generic component"""
        gene_id = self._get_feature_locus_tag(feature)
        if gene_id is None:
            logging.warning('Could not get {} of a feature at location {}'.format(feature.location))
            return
        gene = "RNA_" + gene_id
        if any("5S" in i for i in feature.qualifiers["product"]):
            cat = "generic_5s_rRNAs"
        elif any("16S" in i for i in feature.qualifiers["product"]):
            cat = "generic_16s_rRNAs"
        elif any("23S" in i for i in feature.qualifiers["product"]):
            cat = "generic_23s_rRNAs"
        else:
            cat = 0
        if cat:
            logging.warning("Setting {} to {}".format(gene, cat))
            generic_dict[cat]['enzymes'].append(gene)

    def get_generics_from_genbank(self):
        """ Gets generics from genome."""
        if self.is_reference:
            return None
        contigs = self.contigs
        generic_dict = self.generic_dict
        warn_generics = []
        for contig in tqdm.tqdm(contigs,
                           'Getting generics from Genbank contigs...',
                           bar_format = bar_format):
            for feature in contig.features:
                if "rRNA" in feature.type:
                    self._map_to_a_generic(
                          feature,
                          generic_dict)
        if self.duplicated_genes is not None:
            for d in self.duplicated_genes:
                if not d: continue
                dups = self.gene_dictionary[self.gene_dictionary['Accession-1'].eq(d)]
                generic_dict['generic_{}'.format(d)] = {"enzymes":[i for i in dups['Product'].values if i]}
        for k, v in generic_dict.items():
            if not v:
                warn_generics.append(k)
        if warn_generics:
            self.curation_notes['org.get_generics_from_genbank'].append({
                'msg':'Some generics in me_builder.org.generic_dict are empty.',
                'triggered_by':warn_generics,
                'importance':'high',
                'to_do':'Curate and fill generics in generics.txt or directly in me_builder.org.generic_dict'})

    def _modify_rna_modification_from_load(self,df):
        """Read and preprocess the RNA modification file"""
        d = {}
        for idx,row in df.iterrows():
            mods = ['{}_at_{}'.format(idx,i) for i in row['positions'].split(',')]
            for enz in row['enzymes'].split(' AND '):
                #if enz not in d: d[enz] = []
                #d[enz] = mods
                for mod in mods:
                    d.setdefault(enz, []).append(mod)
        return d

    def _get_rrna_genes(self):
        """Get the genes associated with rRNAs"""
        rrnas = ['generic_5s_rRNAs','generic_16s_rRNAs','generic_23s_rRNAs']
        generic_dict = self.generic_dict
        d = {}
        for key in rrnas:
            rrnaid = key.split('_')[1].upper()
            d[rrnaid] = [i[4:] for i in generic_dict[key]['enzymes']]
        return d

    def process_rna_modifications(self):
        """ Processes RNA modification information."""
        rna_mods = self.rna_modification_df
        self.rna_modification = self._modify_rna_modification_from_load(rna_mods)

        mod_targets = self.rna_modification_targets
        for subunit,genes in self._get_rrna_genes().items():
            mod_df = rna_mods[rna_mods['type'].eq(subunit)].T
            if mod_df.empty:continue
            for mod,row in mod_df.items():
                positions = row['positions'].split(',')
                for g in genes:
                    gene_mods = pandas.DataFrame()
                    if g in mod_targets.index:
                        gene_mods = mod_targets.loc[[g]]
                    d = {}
                    for p in positions:
                        if not gene_mods.empty and gene_mods['position'].eq(p).any():
                            continue
                        d[g] = {'modification':mod,'position':p}
                        df = pandas.DataFrame.from_dict(d).T
                        mod_targets = pandas.concat([mod_targets,df])
        mod_targets.index.name = 'bnum'
        self.rna_modification_targets = mod_targets

    def _check_for_duplicates_within_datasets(self,
                                             info):
        """Checks for duplicates within optinal files"""
        import collections
        warn_dups = {}
        for k,v in tqdm.tqdm(info.items(),
                           'Looking for duplicates within datasets...',
                           bar_format = bar_format,
                           total = len(info)):
            if len(v) != len(set(v)):
                warn_dups[k] = [item for item, count in collections.Counter(v).items() if count > 1]

        if warn_dups:
            self.curation_notes['org.check_for_duplicates'].append({
                'msg':'Some datasets contain duplicate indices or Accession IDs.',
                'triggered_by' : [warn_dups],
                'importance':'critical',
                'to_do':'Remove or fix duplicates. If duplicates are in Accession-1, they are processed as different possibilities to get the same enzyme, so they are added as generic complexes. Check!'})
            if 'Accession-1' in warn_dups and not self.is_reference:
                self.duplicated_genes = warn_dups['Accession-1']
                return
        self.duplicated_genes = None

    def _check_for_duplicates_between_datasets(self,
                                               info):
        """Checks for duplicates between optinal files"""
        cplxs = set(info['complexes_df'])
        rnas = set(info['RNA_df'])
        genes = set(info['gene_dictionary'])
        rxns = set(info['reactions'])
        occ = {}
        for i in tqdm.tqdm(cplxs|rnas|genes|rxns,
                           'Gathering ID occurrences across datasets...',
                           bar_format = bar_format):
            occ[i] = {'proteins':0,'RNAs':0,'genes':0,'reactions':0}
            if i in cplxs:
                occ[i]['proteins'] += 1
            if i in rnas:
                occ[i]['RNAs'] += 1
            if i in genes:
                occ[i]['genes'] += 1
            if i in rxns:
                occ[i]['reactions'] += 1
        df = pandas.DataFrame.from_dict(occ).T
        df = df[df['reactions'] == 1]
        dup_df = df[df.sum(1)>1]
        return dup_df

    def _solve_duplicates_between_datasets(self,
                                           dup_df):
        """Solves found duplicates according to hierarchy"""
        from coralme.builder.helper_functions import change_reaction_id
        for c,row in tqdm.tqdm(dup_df.iterrows(),
                           'Solving duplicates across datasets...',
                           bar_format = bar_format,
                           total=dup_df.shape[0]):
            if row['reactions']:
                change_reaction_id(self.m_model,c,c+'_rxn')
                logging.warning('Changed reaction ID from {} to {} to prevent the conflict between: {}'.format(c,c+'_rxn',' and '.join([j for j,k in row.items() if k])))
            else:
                raise ValueError('The identifier {} is duplicated in {}. Please fix!'.format(c,' and '.join([j for j,k in row.items() if k])))

    def check_for_duplicates(self):
        """ Checks for problematic duplicates in provided files."""
        # Duplicates within datasets
        info = {
            'complexes_df' : list(self.complexes_df.index),
            'RNA_df' : list(self.RNA_df.index),
            'gene_dictionary' : list(self.gene_dictionary.index),
            'reactions' : list([i.id for i in self.m_model.reactions]),
            'Accession-1' : list(self.gene_dictionary['Accession-1'].values)
        }
        self._check_for_duplicates_within_datasets(info)
        dup_df = self._check_for_duplicates_between_datasets(info)
        self._solve_duplicates_between_datasets(dup_df)

    def _is_sequence_divisible_by_three(self,
                                        contig,
                                        f):
        """Checks whether a gene sequence is divisible by three"""
        if f.type == 'source' or 'RNA' in f.type:
            return True
        return not bool(len(f.extract(contig).seq.replace('-', '')) % 3)

    def prune_genbank(self):
        """ Prunes and cleans genome file"""
        if self.is_reference:
            return
        contigs = self.contigs
        exclude_prune_types = list(element_types) #+ ['source','gene']
        new_contigs = []
        warn_sequence = []
        self.all_genes_in_gb = []

        for contig in tqdm.tqdm(contigs,
                           'Pruning GenBank...',
                           bar_format = bar_format):
            new_contig = \
                self._create_genbank_contig(
                               contig.id,
                               contig.seq,
                               contig.name,
                               contig.description,
                               'GenBank')
            new_contig.features = []
            for feature in contig.features:
                if feature.type in ['source', 'region']: # genbanks from gff+fna do not contain a 'source' feature
                    new_contig.features.append(feature)
                    continue
                if feature.type not in exclude_prune_types:
                    continue
                if 'pseudo' in feature.qualifiers:
                    if 'RNA' in feature.type:
                        del feature.qualifiers['pseudo'] # pseudo RNA?
                    elif not self.config.get('include_pseudo_genes', False):
                        continue
                if 'transl_table' not in feature.qualifiers and 'RNA' not in feature.type:
                    feature.qualifiers['transl_table'] = [self.transl_table]
                if not self._is_sequence_divisible_by_three(new_contig,
                                                        feature):
                    d = feature.qualifiers.copy()
                    d['location'] = str(feature.location)
                    warn_sequence.append(d)

                locus_tag = self._get_feature_locus_tag(feature)
                if locus_tag is not None:
                    feature.qualifiers[self.locus_tag] = [locus_tag]
                    self.all_genes_in_gb.append(locus_tag)
                new_contig.features.append(feature)
            if len(new_contig.features) <= 1:
                # only source, no feature
                continue
            new_contigs.append(new_contig)
        self.contigs = new_contigs
        if warn_sequence:
            self.curation_notes['org.prune_genbank'].append({
                    'msg':'Some features contain a sequence that is not divisible by 3.',
                    'triggered_by':warn_sequence,
                    'importance':'critical',
                    'to_do':'Check whether any of these genes are translated in your final ME-model. If so, fix the positions of the gene in genome_modified.gb'})

    def modify_metabolic_reactions(self):
        """ Modifies metabolic reactions from manual input"""
        if self.is_reference:
            return
        m_model = self.m_model
        new_reactions_dict = self.reaction_corrections

        for rxn_id, info in tqdm.tqdm(new_reactions_dict.items(),
                            'Modifying metabolic reactions with manual curation...',
                            bar_format = bar_format,
                            total=len(new_reactions_dict)):
            if info["reaction"] == "eliminate":
                m_model.reactions.get_by_id(rxn_id).remove_from_model()
            else:
                if rxn_id not in m_model.reactions:
                    rxn = cobra.Reaction(rxn_id)
                    m_model.add_reactions([rxn])
                else:
                    rxn = m_model.reactions.get_by_id(rxn_id)
                if info["reaction"]:
                    rxn.build_reaction_from_string(info["reaction"])
                if info["gene_reaction_rule"]:
                    if info["gene_reaction_rule"] == "no_gene":
                        rxn.gene_reaction_rule = ""
                    else:
                        rxn.gene_reaction_rule = info["gene_reaction_rule"]
                if info["name"]:
                    rxn.name = info["name"]

    def modify_metabolites(self):
        """ Modifies metabolites from manual input"""
        if self.is_reference:
            return
        m_model = self.m_model
        new_metabolites_dict = self.metabolite_corrections

        for met_id, info in tqdm.tqdm(new_metabolites_dict.items(),
                            'Modifying metabolites with manual curation...',
                            bar_format = bar_format,
                            total=len(new_metabolites_dict)):
            if info["name"] == "eliminate":
                m_model.metabolites.get_by_id(met_id).remove_from_model()
            else:
                if met_id not in m_model.metabolites:
                    met = cobra.Metabolite(met_id)
                    m_model.add_metabolites([met])
                else:
                    met = m_model.metabolites.get_by_id(met_id)
                if info["name"]:
                    met.name = info["name"]
                if info["formula"]:
                    met.formula = info["formula"]

    def add_manual_complexes(self):
        """ Modifies complexes from manual input"""
        if self.is_reference:
            return
        manual_complexes = self.manual_complexes
        complexes_df = self.complexes_df
        protein_mod = self.protein_mod
        warn_manual_mod = []
        warn_replace = []
        for new_complex, info in tqdm.tqdm(manual_complexes.iterrows(),
                    'Adding manual curation of complexes...',
                    bar_format = bar_format,
                    total=manual_complexes.shape[0]):
            if info["genes"]:
                if new_complex not in complexes_df:
                    complexes_df = \
                    self._add_entry_to_complexes(
                               "",
                               "",
                               new_complex,
                               complexes_df,
                               "Manual")
                complexes_df.loc[new_complex, "genes"] = info["genes"]
                complexes_df.loc[new_complex, "name"] = str(info["name"])

            if info["replace"]:
                if info["replace"] in protein_mod.index:
                    protein_mod = protein_mod.drop(info["replace"])
                else:
                    warn_replace.append(mod_complex)
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

                protein_mod = self._add_entry_to_protein_mod(protein_mod,
                                                             mod_complex,
                                                             new_complex,
                                                             info["mod"],
                                                             "Manual")
        complexes_df.index.name = "complex"

        self.complexes_df = complexes_df
        self.protein_mod = protein_mod

        # Warnings
        if warn_manual_mod or warn_replace:
            if warn_manual_mod:
                self.curation_notes['org.add_manual_complexes'].append({
                    'msg':'Some modifications in protein_corrections.txt are already in me_builder.org.protein_mod and were skipped.',
                    'triggered_by':warn_manual_mod,
                    'importance':'low',
                    'to_do':'Check whether the protein modification specified in protein_corrections.txt is correct and not duplicated.'})
            if warn_replace:
                self.curation_notes['org.add_manual_complexes'].append({
                    'msg':'Some modified proteins marked for replacement in protein_corrections.txt are not in me_builder.org.protein_mod. Did nothing.',
                    'triggered_by':warn_replace,
                    'importance':'low',
                    'to_do':'Check whether the marked modified protein in protein_corrections.txt for replacement is correctly defined.'})

    def get_enzyme_reaction_association(self):
        """
        Maps the M-model GPRs to the provided files and links reactions
        to enzymatic complexes
        """
        if self.is_reference:
            return
        m_model = self.m_model
        org_complexes_df = self.complexes_df
        protein_mod = self.protein_mod
        gene_dictionary = (
            self.gene_dictionary.reset_index()
            .set_index("Accession-1")
        )
        generic_dict = self.generic_dict
        enz_rxn_assoc_dict = {}
        new_generics = {}
        gpr_combination_cutoff = self.config.get("gpr_combination_cutoff",100)

        for rxn in tqdm.tqdm(m_model.reactions,
                    'Getting enzyme-reaction associations...',
                    bar_format = bar_format):
            if rxn.id in self.manual_curation.enz_rxn_assoc_df.data.index:
                # Only complete those not in manual curation
                continue
            unnamed_counter = 0
            rule = str(rxn.gene_reaction_rule)
            if not rule:
                continue
            enz_rxn_assoc_dict[rxn.id] = []
            #rule_list = expand_gpr(listify_gpr(rule)).split(" or ")
            rule_list = coralme.builder.gpr.expand_gpr(rule,threshold=gpr_combination_cutoff)
            #if rule_list != "STOP" and len(rule_list) <= gpr_combination_cutoff:
            if rule_list != "STOP":
                enz_rxn_assoc = []
                reaction_cplx_list = []
                for rule_gene_list in rule_list:
                    identified_genes = [i for i in rule_gene_list if i not in self.skip_genes]
                    if not identified_genes:
                        continue
                    cplx_id = coralme.builder.gpr.find_match(org_complexes_df["genes"].to_dict(),identified_genes)
                    if not cplx_id:
                        if len(identified_genes) > 1:
                            # New cplx not found in BioCyc files
                            cplx_id = "CPLX_{}-{}".format(rxn.id,unnamed_counter)
                            unnamed_counter += 1
                        else:
                            gene = identified_genes[0]
                            cplx_id = "{}-MONOMER".format(gene_dictionary.loc[gene]['Gene Name'])
                        if cplx_id not in org_complexes_df.index:
                            logging.warning("Adding {} to complexes from m_model".format(cplx_id))
                            tmp = pandas.DataFrame.from_dict({
                                cplx_id: {
                                    "name": str(rxn.name),
                                    "genes": " AND ".join(["{}()".format(g) for g in identified_genes]),
                                    "source": "{}({})".format(m_model.id, rxn.id),
                                    }}).T
                            org_complexes_df = pandas.concat([org_complexes_df, tmp], axis = 0, join = 'outer')
                    if cplx_id in protein_mod["Core_enzyme"].values:
                        # Use modifications
                        cplx_mods = protein_mod[
                            protein_mod["Core_enzyme"].eq(cplx_id)
                        ].index
                        for cplx_id in cplx_mods:
                            if "Oxidized" in cplx_id:
                                reaction_cplx_list.append(cplx_id.split("_mod_Oxidized")[0])
                            else:
                                reaction_cplx_list.append(cplx_id)
                    else:
                        # Use base complex
                        reaction_cplx_list.append(cplx_id)
                enz_rxn_assoc_dict[rxn.id] = " OR ".join(reaction_cplx_list)
            else:
                logging.warning('{} contains a GPR rule that has more gene combinations than the specified cutoff ({}). Generifying it.'.format(rxn.id,gpr_combination_cutoff))
                listified_gpr = coralme.builder.gpr.listify_gpr(rule)
                n,rule_dict = coralme.builder.gpr.generify_gpr(listified_gpr,rxn.id,d={},generic_gene_dict=new_generics)
                if not rule_dict: # n in gene_dictionary.index:
                    product = gene_dictionary.loc[n,'Product']
                    rule_dict[product] = n
                    n = product
                n,rule_dict = coralme.builder.gpr.process_rule_dict(n,rule_dict,org_complexes_df["genes"].to_dict(),protein_mod)
                generified_rule = n
                for cplx,rule in rule_dict.items():
                    if 'mod' in cplx:
                        cplx_id = cplx.split('_mod_')[0]
                    else:
                        cplx_id = cplx
                    if 'generic' in cplx_id and cplx_id not in generic_dict:
                        logging.warning("Adding {} to generics from m_model".format(cplx_id))
                        new_generics[cplx_id] = rule.split(' or ')
                        generic_dict[cplx_id] = {
                            'enzymes':[gene_dictionary.loc[i,'Product'] if i in gene_dictionary.index else i for i in rule.split(' or ')]
                        }
                    elif 'generic' not in cplx_id and cplx_id not in org_complexes_df.index:
                        # New cplx not found in BioCyc files
                        logging.warning("Adding {} to complexes from m_model".format(cplx_id))
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
            "", numpy.nan
        ).dropna()  # Remove empty rules

        if not enz_rxn_assoc_df.empty: # Only if it inferred any new GPRs
            self.enz_rxn_assoc_df = pandas.concat([enz_rxn_assoc_df, self.enz_rxn_assoc_df], axis = 0, join = 'outer')
        else:
            logging.warning('No new GPR was inferred. If you provided all GPRs in enzyme_reaction_association.txt, no further action is needed.')
        self.enz_rxn_assoc_df.index.name = "Reaction"
        self.complexes_df = org_complexes_df
        self.protein_mod = protein_mod
