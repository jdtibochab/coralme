import coralme

# guesses = ["dUTPase_c", "zn2_c", "protein_c", "biomass_c", "pg_c", "apoACP_c", "acgam6p_c", "actp_c"]
guesses = ["dUTPase_c", "pg_c", "apoACP_c", "protein_c", "biomass_c", "zn2_c", "actp_c", "f6p_c","ACP_R_c","fdp_c","fe2_c", "mn2_c"]

def test_get_enzyme_reaction_association(shared_builder):
    enz_rxn_assoc = shared_builder.org.enz_rxn_assoc_df
    assert "ARTIFICIAL_RXN_0" in enz_rxn_assoc.index
    assert "ARTIFICIAL_COMPLEX_0" in enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"]
    assert set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_0"].split(" OR ")) == set(enz_rxn_assoc["Complexes"]["ARTIFICIAL_RXN_2"].split(" OR "))
    assert not enz_rxn_assoc["Complexes"].str.contains("MISSING_GENE_0").any(), "Enzyme-reaction association went wrong for ARTIFICIAL_RXN_0 and ARTIFICIAL_RXN_2"

def test_reconstruct(shared_builder):
    shared_builder.build_me_model(overwrite=False)
    assert shared_builder.me_model.id is not None

def test_generification(shared_builder):
    model = shared_builder.me_model
    rxn0 = [r for r in model.reactions.query("ARTIFICIAL_RXN_0") if isinstance(r,coralme.core.reaction.MetabolicReaction)]
    rxn1 = [r for r in model.reactions.query("ARTIFICIAL_RXN_1") if isinstance(r,coralme.core.reaction.MetabolicReaction)]
    rxn2 = [r for r in model.reactions.query("ARTIFICIAL_RXN_2") if isinstance(r,coralme.core.reaction.MetabolicReaction)]

    assert len(rxn1) == 2, "Generification for ARTIFICIAL_RXN_1 went wrong"
    assert len(rxn0) == len(rxn2) == 4, "Generification for ARTIFICIAL_RXN_0 and ARTIFICIAL_RXN_2 went wrong"

def test_troubleshoot(shared_builder):
    shared_builder.troubleshoot(growth_key_and_value = { shared_builder.me_model.mu : 0.001 },guesses=guesses)
