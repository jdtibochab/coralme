from typing import Literal
import coralme
import pytest

from coralme.builder.main import MEBuilder

def test_get_enzyme_reaction_association(shared_builder: MEBuilder):
    enz_rxn_assoc = shared_builder.org.enz_rxn_assoc_df
    assert "ARTIFICIAL_RXN_0" in enz_rxn_assoc.index
    assert "ARTIFICIAL_COMPLEX_0" in enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"]
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"].split(" OR ")) == set("ARTIFICIAL_COMPLEX_0 OR ARTIFICIAL_GENE_1-MONOMER".split(" OR ")), "GPR processing went wrong for ARTIFICIAL_RXN_0"
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_2"].split(" OR ")) == set("ARTIFICIAL_COMPLEX_0 OR ARTIFICIAL_GENE_1-MONOMER".split(" OR ")), "GPR processing went wrong for ARTIFICIAL_RXN_2"

def test_troubleshoot(shared_builder: MEBuilder):
    guesses = ["dUTPase_c", "pg_c", "apoACP_c", "protein_c", "biomass_c", "zn2_c", "actp_c", "f6p_c","ACP_R_c","fdp_c","fe2_c", "mn2_c"]
    shared_builder.troubleshoot(growth_key_and_value = { shared_builder.me_model.mu : 0.001 },guesses=guesses)
    assert shared_builder.me_model.solution is not None

@pytest.mark.parametrize("builder_fixture,reaction_id",
                         [
                             ("shared_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_COMPLEX_0"),
                             ("shared_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                             ("shared_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_COMPLEX_0"),
                             ("shared_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                             
                            # TODO: fix this test
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_COMPLEX_0"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_1_FWD_CPLX_ARTIFICIAL_RXN_1-0"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_COMPLEX_0"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                             ])
def test_gpr_build(builder_fixture: Literal['shared_builder'],reaction_id: Literal['ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_COMPLEX_0'] | Literal['ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_GENE_1-MONOMER'] | Literal['ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_COMPLEX_0'] | Literal['ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_GENE_1-MONOMER'],request: pytest.FixtureRequest):
    builder = request.getfixturevalue(builder_fixture)  # Retrieve the actual fixture
    model = builder.me_model
    assert model.reactions.has_id(reaction_id), "GPR build for {} went wrong".format(reaction_id)

def test_generification(shared_generification_builder: MEBuilder):
    # Check generification
    enz_rxn_assoc = shared_generification_builder.org.enz_rxn_assoc_df
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_1"].split(" OR ")) == set(["CPLX_ARTIFICIAL_RXN_1-0"]), "GPR processing went wrong for ARTIFICIAL_RXN_1"
