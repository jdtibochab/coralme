from coralme.builder.main import MEBuilder
# from importlib.resources import files
import anyconfig
import pytest
# dir = str(files("coralme"))

def create_builder(temp_dir, troubleshoot = False, reference = True, builder_type = None):
  temp_dir = str(temp_dir)
  print(temp_dir)
  config = {
    # Inputs
    "m-model-path": "./tests/data/m_model.json", # Path to model file
    "genbank-path": "./tests/data/genome.gb", # Path to genome genbank file
    # Outputs
    "df_gene_cplxs_mods_rxns": "building_data/OSM.xlsx", # Desired output path of OSM := out_directory + df_gene_cplxs_mods_rxns
    "out_directory": "{:s}/base_model_{:s}/".format(temp_dir, builder_type), # Output directory
    "log_directory": "{:s}/base_model_{:s}/".format(temp_dir, builder_type), # Log directory
    "locus_tag": "locus_tag", # What IDs were used in the M-model? e.g. locus_tag, old_locus_tag, protein_id
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
    builder.troubleshoot(growth_key_and_value = { builder.me_model.mu.magnitude : 0.01 },guesses=guesses)
  assert builder.me_model.id is not None
  return builder

@pytest.fixture(scope="session")
def temp_dir(tmp_path_factory):
  path = tmp_path_factory.mktemp("session")
  yield path

@pytest.fixture(scope="session")
def shared_builder(temp_dir):
  pytest.shared_builder = create_builder(temp_dir, troubleshoot = False, reference = 'iJL1678b', builder_type = 'no-troubleshooted')

@pytest.fixture(scope="session")
def shared_builder_bsub_reference(temp_dir):
  pytest.shared_builder_bsub_reference = create_builder(temp_dir, troubleshoot = False, reference = 'iJT964', builder_type = 'no-troubleshooted-bsub-reference')

@pytest.fixture(scope="session")
def shared_builder_troubleshooted(temp_dir):
  pytest.shared_builder_troubleshooted = create_builder(temp_dir, troubleshoot = True, reference = 'iJL1678b', builder_type = 'troubleshooted')

@pytest.fixture(scope="session")
def shared_builder_troubleshooted_bsub_reference(temp_dir):
  pytest.shared_builder_troubleshooted_bsub_reference = create_builder(temp_dir, troubleshoot = True, reference = 'iJT964', builder_type = 'troubleshooted-bsub-reference')

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

