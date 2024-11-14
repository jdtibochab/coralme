from coralme.solver.solver import ME_NLP
import cobra

def get_nlp(model):
    # Call construct LP problem function from model to get precursor objects.
    # lamdify = True -> Creates lambdify functions to calculate bounds as a function of mu
    # per_position = True -> LB and UB bounds as list of lambdify instead of a lambdify to 
    # be able to change individual values
    Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = model.construct_lp_problem(lambdify=True,per_position=True)

    # Construct NLP object from precursor objects
    me_nlp = ME_NLP(Sf, Se,b, c, lb, ub,  cs, atoms, lambdas)
    return me_nlp

def get_feasibility(me_nlp, basis=None):
    x_new,y_new,z_new,stat_new,hs_new = me_nlp.solvelp(0.001,basis,'quad')
    return (True,hs_new) if stat_new=="optimal" else (False,hs_new)

def optimize(rxn_index_dct,met_index_dct,me_nlp,max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100,
		tolerance = 1e-6, precision = 'quad', verbose = True,basis=None):
    muopt, xopt, yopt, zopt, basis, stat = me_nlp.bisectmu(
				mumax = max_mu,
				mumin = min_mu,
				maxIter = maxIter,
				tolerance = tolerance,
				precision = precision,
				verbose = verbose,
                basis=basis)

    if stat == 'optimal':
        #f = sum([ rxn.objective_coefficient * xopt[idx] for idx, rxn in enumerate(self.reactions) ])
        x_primal = xopt[ 0:len(rxn_index_dct) ]   # The remainder are the slacks
        x_dict = { rxn : xopt[idx] for rxn,idx in rxn_index_dct.items() }
        #y = pi
        # J = [S; c]
        y_dict = { met : yopt[idx] for met,idx in met_index_dct.items() }
        z_dict = { rxn : zopt[idx] for rxn,idx in rxn_index_dct.items() }
        #y_dict['linear_objective'] = y[len(y)-1]

        #self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
        return cobra.core.Solution(
            objective_value = muopt,
            status = stat,
            fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
            reduced_costs = z_dict,
            shadow_prices = y_dict,
            ),basis
    else:
        return None,None