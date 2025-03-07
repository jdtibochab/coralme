from coralme.builder.main import MEBuilder
from importlib.resources import files
import anyconfig
import pytest
dir = str(files("coralme"))

def create_builder():
    config = {
      # Inputs
      "m-model-path": "{}/tests/data/m_model.json".format(dir), # Path to model file
      "genbank-path": "{}/tests/data/genome.gb".format(dir), # Path to genome genbank file
      # Outputs
      "df_gene_cplxs_mods_rxns": "{}/tests/data/base_model/building_data/OSM.xlsx".format(dir), # Desired output path of OSM
      "out_directory": "{}/tests/data/base_model/".format(dir), # Output directory
      "log_directory": "{}/tests/data/base_model/".format(dir), # Log directory
      "locus_tag": "locus_tag", # What IDs were used in the M-model? e.g. locus_tag, old_locus_tag
      "run_bbh_blast" : True,
      "dev_reference" : True,
      "blast_threads" : 4,
      "ME-Model-ID" : "EXAMPLE-ME" # Name of the ME-model
    }
    with open("{}/tests/data/organism.json".format(dir), 'r') as infile:
        config.update(anyconfig.load(infile))
    config.update({
      "biocyc.genes": "{}/tests/data/genes.txt".format(dir),
      "biocyc.prots": "{}/tests/data/proteins.txt".format(dir),
      "biocyc.TUs": "{}/tests/data/TUs.txt".format(dir),
      "biocyc.RNAs": "{}/tests/data/RNAs.txt".format(dir),
      "biocyc.seqs": "{}/tests/data/sequences.fasta".format(dir),
    })
    builder = MEBuilder(**config)
    builder.generate_files(overwrite=True)
    builder.save_builder_info()
    builder.build_me_model(overwrite=False)
    assert builder.me_model.id is not None
    return builder

@pytest.fixture(scope="session")
def shared_builder():
    return create_builder()

# TODO: fix this test
# @pytest.fixture(scope="session")
# def shared_generification_builder(shared_builder):
#     builder = MEBuilder(**shared_builder.configuration.copy())
#     # Generate files
#     builder.configuration["m-model-path"] = "{}/tests/data/m_model-generification.json".format(dir) # Path to model file
#     builder.configuration["out_directory"] = "{}/tests/data/generified_model/".format(dir)
#     builder.configuration["log_directory"] = "{}/tests/data/generified_model/".format(dir)
#     builder.configuration["df_gene_cplxs_mods_rxns"] = "{}/tests/data/generified_model/building_data/OSM.xlsx".format(dir) # Desired output path of OSM
#     builder.generate_files(overwrite=True)
#     # Build
#     builder.save_builder_info()
#     builder.build_me_model(overwrite=False)
#     return builder

