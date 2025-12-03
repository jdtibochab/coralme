import coralme
import pytest

def test_get_enzyme_reaction_association(shared_builder):
    builder1 = pytest.shared_builder
    enz_rxn_assoc = builder1.org.enz_rxn_assoc_df
    assert "ARTIFICIAL_RXN_0" in enz_rxn_assoc.index
    assert "ARTIFICIAL_COMPLEX_0" in enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"]
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"].split(" OR ")) == set("ARTIFICIAL_COMPLEX_0 OR ARTIFICIAL_GENE_1-MONOMER".split(" OR ")), "GPR processing went wrong for ARTIFICIAL_RXN_0"
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_2"].split(" OR ")) == set("ARTIFICIAL_COMPLEX_0 OR ARTIFICIAL_GENE_1-MONOMER".split(" OR ")), "GPR processing went wrong for ARTIFICIAL_RXN_2"

# def test_troubleshoot(shared_builder):
#     builder1 = pytest.shared_builder
#     guesses = ["dUTPase_c", "pg_c", "apoACP_c", "protein_c", "biomass_c", "zn2_c", "actp_c", "f6p_c","ACP_R_c","fdp_c","fe2_c", "mn2_c"]
#     builder1.troubleshoot(growth_key_and_value = { builder1.me_model.mu.magnitude : 0.001 },guesses=guesses)
#     assert builder1.me_model.solution is not None
#
# def test_troubleshoot_bsub_reference(shared_builder_bsub_reference):
#     builder1 = pytest.shared_builder_bsub_reference
#     guesses = ["dUTPase_c", "pg_c", "apoACP_c", "protein_c", "biomass_c", "zn2_c", "actp_c", "f6p_c","ACP_R_c","fdp_c","fe2_c", "mn2_c"]
#     builder1.troubleshoot(growth_key_and_value = { builder1.me_model.mu.magnitude : 0.001 },guesses=guesses)
#     assert builder1.me_model.solution is not None

@pytest.mark.parametrize("builder_fixture,reaction_id",
                         [
                             ("shared_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_COMPLEX_0"),
                             ("shared_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                             ("shared_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_COMPLEX_0"),
                             ("shared_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                            # TODO: fix this tests
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_COMPLEX_0"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_1_FWD_CPLX_ARTIFICIAL_RXN_1-0"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_COMPLEX_0"),
                            #  ("shared_generification_builder","ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_GENE_1-MONOMER"),
                             ])
def test_gpr_build(builder_fixture,reaction_id,request):
    builder = pytest.shared_builder # request.getfixturevalue(builder_fixture)  # Retrieve the actual fixture
    model = builder.me_model
    assert model.reactions.has_id(reaction_id), "GPR build for {} went wrong".format(reaction_id)

# TODO: fix this test
# def test_generification(shared_generification_builder):
#     # Check generification
#     enz_rxn_assoc = shared_generification_builder.org.enz_rxn_assoc_df
#     assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_1"].split(" OR ")) == set(["CPLX_ARTIFICIAL_RXN_1-0"]), "GPR processing went wrong for ARTIFICIAL_RXN_1"
