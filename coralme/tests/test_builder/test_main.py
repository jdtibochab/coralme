import coralme

def test_get_enzyme_reaction_association(shared_builder):
    enz_rxn_assoc = shared_builder.org.enz_rxn_assoc_df
    assert "ARTIFICIAL_RXN_0" in enz_rxn_assoc.index
    assert "ARTIFICIAL_COMPLEX_0" in enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"]
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"].split(" OR ")) == set("ARTIFICIAL_COMPLEX_0 OR ARTIFICIAL_GENE_1-MONOMER".split(" OR ")), "GPR processing went wrong for ARTIFICIAL_RXN_0"
    # assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_1"].split(" OR ")) == set(["CPLX_ARTIFICIAL_RXN_1-0"]), "GPR processing went wrong for ARTIFICIAL_RXN_1"
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_2"].split(" OR ")) == set("ARTIFICIAL_COMPLEX_0 OR ARTIFICIAL_GENE_1-MONOMER".split(" OR ")), "GPR processing went wrong for ARTIFICIAL_RXN_2"

def test_reconstruct(shared_builder):
    shared_builder.build_me_model(overwrite=False)
    assert shared_builder.me_model.id is not None

def test_gpr_builds(shared_builder):
    model = shared_builder.me_model
    rxn0 = [r for r in model.reactions.query("ARTIFICIAL_RXN_0") if isinstance(r,coralme.core.reaction.MetabolicReaction)]
    rxn2 = [r for r in model.reactions.query("ARTIFICIAL_RXN_2") if isinstance(r,coralme.core.reaction.MetabolicReaction)]
    assert model.reactions.has_id("ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_COMPLEX_0"), "GPR build for ARTIFICIAL_RXN_0 went wrong"
    assert model.reactions.has_id("ARTIFICIAL_RXN_0_FWD_ARTIFICIAL_GENE_1-MONOMER"), "GPR build for ARTIFICIAL_RXN_0 went wrong"
    assert model.reactions.has_id("ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_COMPLEX_0"), "GPR build for ARTIFICIAL_RXN_2 went wrong"
    assert model.reactions.has_id("ARTIFICIAL_RXN_2_FWD_ARTIFICIAL_GENE_1-MONOMER"), "GPR build for ARTIFICIAL_RXN_0 went wrong"
    
# def test_generification(shared_builder):
#     model = shared_builder.me_model
#     rxn1 = [r for r in model.reactions.query("ARTIFICIAL_RXN_1") if isinstance(r,coralme.core.reaction.MetabolicReaction)]
#     assert len(rxn1) == 2, "Generification for ARTIFICIAL_RXN_1 went wrong"

def test_troubleshoot(shared_builder):
    guesses = ["dUTPase_c", "pg_c", "apoACP_c", "protein_c", "biomass_c", "zn2_c", "actp_c", "f6p_c","ACP_R_c","fdp_c","fe2_c", "mn2_c"]
    shared_builder.troubleshoot(growth_key_and_value = { shared_builder.me_model.mu : 0.001 },guesses=guesses)
    assert shared_builder.me_model.solution is not None
