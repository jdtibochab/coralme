#!/usr/bin/python3
import pandas
import pickle
import coralme

_ = coralme.check_installed_packages()

with open('lp.pkl', 'rb') as infile:
    lp = pickle.load(infile)

nlp = coralme.solver.solver.ME_NLP(**lp)

# test qwarmLP.qwarmlp
x, pi, rc, inform, hs = nlp.solvelp(muf = 0.1, precision = 'quad')
assert inform == 'optimal'

# test warmLP.warmlp
x, pi, rc, inform, hs = nlp.solvelp(muf = 0.1, precision = 'double')
assert inform == 'optimal'

# test qvaryME.qvaryme
obj_inds0 = [ idx for idx, rxn in enumerate(lp['Lr']) for j in range(0, 2) if rxn in ['EX_ac_e'] ]
obj_coeffs = [ ci for rxn in lp['Lr'] for ci in (1.0, -1.0) if rxn in ['EX_ac_e'] ]

# varyME is a specialized method for multiple min/maximization problems
obj_inds0, nVary, obj_vals = nlp.varyme(mu_fixed = 0.1, obj_inds0 = obj_inds0, obj_coeffs = obj_coeffs, basis = None, verbosity = False)

# Return result consistent with cobrapy FVA
fva_result = {(lp['Lr'][obj_inds0[2*i]]) : { 'maximum':obj_vals[2*i], 'minimum':obj_vals[2*i+1] } for i in range(0, nVary//2) }
assert pandas.DataFrame(fva_result).T.empty is False
