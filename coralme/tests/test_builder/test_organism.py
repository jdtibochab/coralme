from coralme.builder.organism import Organism
from importlib.resources import files
import anyconfig

dir = str(files("coralme"))

input_nobiocyc = {
  # Inputs
  "m-model-path": "{}/tests/data/m_model.json".format(dir), # Path to model file
  "genbank-path": "{}/tests/data/genome.gb".format(dir), # Path to genome genbank file
  # Outputs
  "df_gene_cplxs_mods_rxns": "{}/tests/data/building_data/OSM.xlsx".format(dir), # Desired output path of OSM
  "out_directory": "{}/tests/data/".format(dir), # Output directory
  "log_directory": "{}/tests/data/".format(dir), # Log directory
  "locus_tag": "locus_tag", # What IDs were used in the M-model? e.g. locus_tag, old_locus_tag
  "ME-Model-ID" : "EXAMPLE-ME" # Name of the ME-model
}
organism = {}
with open("{}/tests/data/organism.json".format(dir), 'r') as infile:
    organism.update(anyconfig.load(infile))

def test_organism_init():
    config = organism.copy()
    config.update(input_nobiocyc)
    org = Organism(config,is_reference=False)
    ref = Organism(config,is_reference=True,available_reference_models={'iJL1678b':'locus_tag'})

def test_get_organism_nobiocyc():
    config = organism.copy()
    config.update(input_nobiocyc)
    Organism(config,is_reference=False).get_organism()
    Organism(config,is_reference=True,available_reference_models={'iJL1678b':'locus_tag'}).get_organism()

def test_get_organism_withbiocyc():
    config = organism.copy()
    config.update(input_nobiocyc)
    config.update({
      "biocyc.genes": "{}/tests/data/genes.txt".format(dir),
      "biocyc.prots": "{}/tests/data/proteins.txt".format(dir),
      "biocyc.TUs": "{}/tests/data/TUs.txt".format(dir),
      "biocyc.RNAs": "{}/tests/data/RNAs.txt".format(dir),
      "biocyc.seqs": "{}/tests/data/sequences.fasta".format(dir),
    })
    Organism(config,is_reference=False).get_organism()
    Organism(config,is_reference=True,available_reference_models={'iJL1678b':'locus_tag'}).get_organism()
