from coralme.builder.main import MEBuilder
from importlib.resources import files
import anyconfig
import pytest
# dir = str(files("coralme"))

def create_builder(troubleshoot = False, reference = True):
  config = {
    # Inputs
    "m-model-path": "./tests/data/m_model.json", # Path to model file
    "genbank-path": "./tests/data/genome.gb", # Path to genome genbank file
    # Outputs
    "df_gene_cplxs_mods_rxns": "./tests/data/base_model/building_data/OSM.xlsx", # Desired output path of OSM
    "out_directory": "./tests/data/base_model/", # Output directory
    "log_directory": "./tests/data/base_model/", # Log directory
    "locus_tag": "locus_tag", # What IDs were used in the M-model? e.g. locus_tag, old_locus_tag
    "run_bbh_blast" : True,
    "dev_reference" : reference,
    "blast_threads" : 4,
    "ME-Model-ID" : "EXAMPLE-ME" # Name of the ME-model
  }
  with open("./tests/data/organism.json", 'r') as infile:
      config.update(anyconfig.load(infile))
  config.update({
    "biocyc.genes": "./tests/data/genes.txt",
    "biocyc.prots": "./tests/data/proteins.txt",
    "biocyc.TUs": "./tests/data/TUs.txt",
    "biocyc.RNAs": "./tests/data/RNAs.txt",
    "biocyc.seqs": "./tests/data/sequences.fasta",
  })
  builder = MEBuilder(**config)
  builder.generate_files(overwrite=True)
  builder.save_builder_info()
  builder.build_me_model(overwrite=False)
  if troubleshoot:
    guesses = ["dUTPase_c", "pg_c", "apoACP_c", "protein_c", "biomass_c", "zn2_c", "actp_c", "f6p_c","ACP_R_c","fdp_c","fe2_c", "mn2_c"]
    builder.troubleshoot(growth_key_and_value = { builder.me_model.mu.magnitude : 0.001 },guesses=guesses)
  assert builder.me_model.id is not None
  return builder

@pytest.fixture(scope="session")
def shared_builder():
  pytest.shared_builder = create_builder(troubleshoot = False)

@pytest.fixture(scope="session")
def shared_builder_bsub_reference():
  pytest.shared_builder_bsub_reference = create_builder(troubleshoot = False, reference = 'iJT964')

@pytest.fixture(scope="session")
def shared_builder_troubleshooted():
  pytest.shared_builder_troubleshooted = create_builder(troubleshoot = True)

@pytest.fixture(scope="session")
def shared_builder_troubleshooted_bsub_reference():
  pytest.shared_builder_troubleshooted_bsub_reference = create_builder(troubleshoot = True, reference = 'iJT964')

# TODO: fix this test
# @pytest.fixture(scope="session")
# def shared_generification_builder(shared_builder):
#     builder = MEBuilder(**shared_builder.configuration.copy())
#     # Generate files
#     builder.configuration["m-model-path"] = "./tests/data/m_model-generification.json" # Path to model file
#     builder.configuration["out_directory"] = "./tests/data/generified_model/"
#     builder.configuration["log_directory"] = "./tests/data/generified_model/"
#     builder.configuration["df_gene_cplxs_mods_rxns"] = "./tests/data/generified_model/building_data/OSM.xlsx" # Desired output path of OSM
#     builder.generate_files(overwrite=True)
#     # Build
#     builder.save_builder_info()
#     builder.build_me_model(overwrite=False)
#     return builder

