from coralme.builder.organism import Organism
import anyconfig


def input_nobiocyc(temp_dir):
  temp_dir = str(temp_dir)
  return {
    # Inputs
    "m-model-path": "./tests/data/m_model.json", # Path to model file
    "genbank-path": "./tests/data/genome.gb", # Path to genome genbank file
    # Outputs
    "df_gene_cplxs_mods_rxns": "building_data/OSM.xlsx", # Desired output path of OSM
    "out_directory": "{:s}/base_model-only_organism".format(temp_dir), # Output directory
    "log_directory": "{:s}/base_model-only_organism".format(temp_dir), # Log directory
    "locus_tag": "locus_tag", # What IDs were used in the M-model? e.g. locus_tag, old_locus_tag
    "ME-Model-ID" : "EXAMPLE-ME" # Name of the ME-model
  }

organism = {}
with open("./tests/data/organism.json", 'r') as infile:
    organism.update(anyconfig.load(infile))

def test_organism_init(temp_dir):
    config = organism.copy()
    config.update(input_nobiocyc(temp_dir))
    org = Organism(config,is_reference=False)
    ref = Organism(config,is_reference=True,available_reference_models={'iJL1678b':'locus_tag'})

def test_get_organism_nobiocyc(temp_dir):
    config = organism.copy()
    config.update(input_nobiocyc(temp_dir))
    Organism(config,is_reference=False).get_organism()
    Organism(config,is_reference=True,available_reference_models={'iJL1678b':'locus_tag'}).get_organism()

def test_get_organism_withbiocyc(temp_dir):
    config = organism.copy()
    config.update(input_nobiocyc(temp_dir))
    config.update({
      "biocyc.genes": "./tests/data/genes.txt",
      "biocyc.prots": "./tests/data/proteins.txt",
      "biocyc.TUs": "./tests/data/TUs.txt",
      "biocyc.RNAs": "./tests/data/RNAs.txt",
      "biocyc.seqs": "./tests/data/sequences.fasta",
    })
    Organism(config,is_reference=False).get_organism()
    Organism(config,is_reference=True,available_reference_models={'iJL1678b':'locus_tag'}).get_organism()
