#!/usr/bin/python3

import sys
import numpy
import scipy
import sympy

from coralme.solver import qwarmLP, warmLP, qvaryME

def makeME_VA(S, b, c, xl, xu, csense, obj_inds, obj_coeffs):
    """
    Creates ME_LP data for qvaryME, solved using qMINOS with warm-start
    obj_inds: obj column indices
    obj_coeffs: explicitly state -1 or 1 for min or max
    Thus, obj_inds & obj_coeffs have 2*n elements
    [LY] 11 Aug 2015: first version
    [RS] 13 Jul 2023: Adapted to new structure of arguments
    """

    # qMINOS requires an extra row holding the objective
    # Put ANY non-zero for all columns that will be min/maxed
    Sm, Sn = S.shape
    c = [ 0. for j in range(0, Sn)]
    for j, v in zip(obj_inds, obj_coeffs):
        c[j] = 1.0

    # Just change to 1-based indexing for Fortran
    obj_indsf = [ i+1 for i in obj_inds ]

    J, ne, P, I, V, bl, bu = makeME_LP(S, b, c, xl, xu, csense)

    return J, ne, P, I, V, bl, bu, obj_indsf

def makeME_LP(S, b, c, xl, xu, csense):
    """
    Create simple LP for qMINOS and MINOS
    Inputs:
    nlp_compat  Make matrices compatible with NLP so that basis can
                be used to warm start NLP by setting
    12 Aug 2015: first version
    13 Jul 2023: Adapted to new structure of arguments by Rodrigo Santibanez
    """

    # c is added as a free (unbounded slacks) row,
    # so that MINOS treats problem as an LP - Ding Ma
    J = scipy.sparse.vstack((S, [c]), dtype = float).tocsc()
    J.sort_indices()

    b2 = b + [0.0]
    m, n = J.shape
    ne = J.nnz
    # Finally, make the P, I, J, V, as well
    # Row indices: recall fortran is 1-based indexing
    I = [ i+1 for i in J.indices ]
    V = J.data
    # Pointers to start of each column
    # Just change to 1-based indexing for Fortran
    P = [ pi+1 for pi in J.indptr ]

    # Make primal and slack bounds
    bigbnd = 1e+40
    # For csense==E rows (equality)
    sl = [ bi for bi in b2 ]
    su = [ bi for bi in b2 ]

    # It can be avoided since csense contains always 'E'
    for row, csen in enumerate(csense):
        if csen == 'L':
            sl[row] = -bigbnd
        elif csen == 'G':
            su[row] = +bigbnd

    # Objective row has free bounds
    sl[m - 1] = -bigbnd
    su[m - 1] = +bigbnd

    bl = scipy.sparse.vstack([ numpy.matrix(xl).transpose(), numpy.matrix(sl).transpose() ], dtype = float)
    bu = scipy.sparse.vstack([ numpy.matrix(xu).transpose(), numpy.matrix(su).transpose() ], dtype = float)

    return J, ne, P, I, V, bl, bu

# Modified from solvemepy.me2
class ME_NLP:
    """
    Contains the data matrices needed for solving ME as an NLP using qMINOS
    """
    #def __init__(self, me, growth_rxn = 'biomass_dilution'):
    def __init__(self, Sf, Se, b, c, xl, xu, cs, mu = sympy.Symbol('mu', positive = True), lambdas = None, Lr = list(), Lm = list()):
        # The ME model object
        #self.me = me

        # Reformulation of ME to NLP
        self.Sf = Sf # dict of float coefficients
        self.Se = Se # dict of symbolic coefficients
        self.b  = b  # accumulation
        self.c  = c  # objective
        self.xl = xl # lower bounds
        self.xu = xu # upper bounds
        self.cs = cs # constraint sense
        self.mu = mu # symbol to evaluate
        self.fn = lambdas
        self.Lr = Lr
        self.Lm = Lm

        # Initialize solver options
        self.init_solver_opts()

    def init_solver_opts(self):
        #----------------------------------------------------
        # Solver options
        self.opt_strwidth  = {}
        self.opt_realwidth = {}
        self.opt_intwidth  = {}
        self.opt_strlist   = {}
        self.opt_intdict   = {}
        self.opt_realdict  = {}
        self.opt_stropts   = {}
        self.opt_intopts   = {}
        self.opt_intvals   = {}
        self.opt_realopts  = {}
        self.opt_realvals  = {}
        #----------------------------------------------------
        # NLP solver options
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['nlp'] = 72
        self.opt_realwidth['nlp'] = 55
        self.opt_intwidth['nlp'] = 55
        self.opt_strlist['nlp'] = [
            'Maximize',     # Default obj sense is to maximize
            'Completion full',
            'Print level (jflxb) 00001',
            'Solution No'
            ]
        self.opt_intdict['nlp'] = {
            'Major iterations': 1000,
            'Superbasics limit': 40,
            'Verify level': 0,
            'Scale option': 2,
            'Partial price': 1,
            'Iterations': 10000,
            'Print frequency': 100000,
            'Summary level': 0,
            'Summary frequency': 100,
            'Solution file': 9,
            'New basis file': 11,
            'Save frequency': 500000
            }
        self.opt_realdict['nlp'] = {
            'Penalty parameter':100.0,
            'LU factor tol': 1.1,
            'LU update tol': 1.1,
            'LU singularity tol': 1e-30,
            'Feasibility tol': 1e-15,
            'Optimality tol': 1e-15,
            'Unbounded step size': 1e+30
            }

        #----------------------------------------------------
        # LP options
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['lp'] = 72
        self.opt_realwidth['lp'] = 55
        self.opt_intwidth['lp'] = 55
        self.opt_strlist['lp'] = [
            'Maximize',     # Default obj sense is to maximize
            'Solution No'
            ]
        self.opt_intdict['lp'] = {
            'New basis file': 11,
            'Save frequency': 500000,
            'Print level': 0,
            'Print frequency': 100000,
            'Scale option': 2,
            'Iteration limit': 2000000,
            'Expand frequency': 100000
            }
        self.opt_realdict['lp'] = {
            'Penalty parameter':100.0,
            'LU factor tol': 10.0,
            'LU update tol': 10.0,
            'LU singularity tol': 1e-30,
            'Feasibility tol': 1e-20,
            'Optimality tol': 1e-20,
            'Unbounded step size': 1e+30
            }

        #----------------------------------------------------
        # LP options: double-precision (can't set as strict tols)
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['lp_d'] = 72
        self.opt_realwidth['lp_d'] = 55
        self.opt_intwidth['lp_d'] = 55
        self.opt_strlist['lp_d'] = [
            'Maximize',     # Default obj sense is to maximize
            'Solution No'
            ]
        self.opt_intdict['lp_d'] = {
            'New basis file': 11,
            'Save frequency': 500000,
            'Print level': 0,
            'Print frequency': 100000,
            'Scale option': 2,
            'Iteration limit': 2000000,
            'Expand frequency': 100000
            }
        self.opt_realdict['lp_d'] = {
            'Penalty parameter':100.0,
            'LU factor tol': 1.9,
            'LU update tol': 1.9,
            'LU singularity tol': 1e-12,
            'Feasibility tol': 1e-7,
            'Optimality tol': 1e-7,
            'Unbounded step size': 1e+18
            }

    def get_solver_opts(self, prob = 'lp'):
        """
        Return options that will be passed as arguments to minoss
        """
        if sys.version_info >= (3, 11): # qMINOS was compiled using numpy 2.0
            stropts = [ s.ljust(self.opt_strwidth[prob]).encode() for s in self.opt_strlist[prob] ]
            intopts = [ s.ljust(self.opt_intwidth[prob]).encode() for s in self.opt_intdict[prob].keys() ]
            realopts = [ s.ljust(self.opt_realwidth[prob]).encode() for s in self.opt_realdict[prob].keys() ]
        else:
            lst = [ s.ljust(self.opt_strwidth[prob]) for s in self.opt_strlist[prob] ]
            stropts = numpy.array(numpy.array(lst, dtype = 'c').T)
            lst = [ s.ljust(self.opt_intwidth[prob]) for s in self.opt_intdict[prob].keys() ]
            intopts = numpy.array(numpy.array(lst, dtype = 'c').T)
            lst = [ s.ljust(self.opt_realwidth[prob]) for s in self.opt_realdict[prob].keys() ]
            realopts = numpy.array(numpy.array(lst, dtype = 'c').T)

        intvals = numpy.array([ self.opt_intdict[prob][k] for k in self.opt_intdict[prob].keys() ], dtype = 'i4')
        realvals = numpy.array([ self.opt_realdict[prob][k] for k in self.opt_realdict[prob].keys() ], dtype = 'd')

        self.opt_stropts[prob] = stropts
        self.opt_intopts[prob] = intopts
        self.opt_realopts[prob] = realopts
        self.opt_intvals[prob] = intvals
        self.opt_realvals[prob]= realvals

        nStrOpts = len(self.opt_strlist[prob])
        nIntOpts = len(self.opt_intdict[prob].keys())
        nRealOpts = len(self.opt_realdict[prob].keys())

        return stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts

    def make_lp(self, muf, obj_inds0 = None, obj_coeffs = None):
        """
        Construct LP problem for qMINOS or MINOS.
        """

        # self.mu is a set of sympy Symbols
        dct = { symbol:muf for symbol in self.mu }

        if self.fn is None:
            # If lambdify == False
            xl = [ float(x.xreplace(dct)) if hasattr(x, 'subs') else x for x in self.xl ]
            xu = [ float(x.xreplace(dct)) if hasattr(x, 'subs') else x for x in self.xu ]
            Se = { k:float(x.xreplace(dct)) if hasattr(x, 'subs') else x for k,x in self.Se.items() }
            self.Sf.update(Se)
        else:
            # If lambdify == True
            if isinstance(self.xl,list):
                # If per_position == True
                xl = [ fn(*[muf]*len(dct)) for fn in self.xl ]
                xu = [ fn(*[muf]*len(dct)) for fn in self.xu ]
                Se = { k:fn(*[muf]*len(dct)) for k,fn in self.fn.items() }
            else:
                # If per_position == False
                xl = self.xl(*[muf]*len(dct))
                xu = self.xu(*[muf]*len(dct))
                Se = { k:v for k,v in zip(self.fn[0], self.fn[1](*[muf]*len(dct))) }
            self.Sf.update(Se)

        Sp = scipy.sparse.dok_matrix((len(self.b), len(self.c)), dtype = float)
        for idx, idj in self.Sf.keys():
            Sp[idx, idj] = self.Sf[idx, idj]

        if obj_inds0 and obj_coeffs:
            J, ne, P, I, V, bl, bu, obj_inds = makeME_VA(Sp, self.b, self.c, xl, xu, self.cs, obj_inds0, obj_coeffs)
        else:
            J, ne, P, I, V, bl, bu           = makeME_LP(Sp, self.b, self.c, xl, xu, self.cs)
            obj_inds = None

        # Solve a single LP
        m, n = J.shape
        ha = I
        ka = P
        ad = V
        bld = [ bi for bi in bl.toarray() ]
        bud = [ bi for bi in bu.toarray() ]
        nb = m + n
        hs = numpy.zeros(nb, numpy.dtype('i4'))

        return m, n, ha, ka, ad, bld, bud, hs, obj_inds

    def solvelp(self, muf, basis = None, precision = 'quad', probname = 'me_lp', output_as_dict = False):
        """
        x, status, hs = solvelp(self, muf, basis = None, precision = 'quad')

        Solve LP at mu using qMINOS or MINOS.
        Pass the basis (hs) back and forth with Fortran for warm-start.

        Inputs:
        muf: fixed growth rate
        basis: basis vector

        Outputs:
        x: primal solution
        status: solver status
        hs: basis
        """
        #me = self.me

        m, n, ha, ka, ad, bld, bud, hs0, _ = self.make_lp(muf)

        hs = basis
        if hs is None:
            warm = False
            hs = hs0
        else:
            warm = True

        inform = numpy.array(0)
        precision = precision.lower()

        if precision == 'quad':
            optimizer = qwarmLP.qwarmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

        elif precision == 'double':
            optimizer = warmLP.warmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp_d')

        elif precision == 'dq' or precision == 'dqq':
            # D
            self.opt_intdict['lp_d']['Scale option'] = 2
            optimizer = warmLP.warmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp_d')

            # Q1: pass optimal basis hs and scale = 2
            warm = True
            self.opt_intdict['lp']['Scale option'] = 2
            optimizer = qwarmLP.qwarmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

            # Last Q2 if requested: pass optimal basis hs and scale = 0
            if precision == 'dqq':
                self.opt_intdict['lp']['Scale option'] = 0
                optimizer = qwarmLP.qwarmlp
                stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

                # Kindly reset scale option to default
                self.opt_intdict['lp']['Scale option'] = 2

        x, pi, rc = optimizer(
            inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
            stropts, intopts, realopts, intvals, realvals,
            nstropts = nStrOpts, nintopts = nIntOpts, nrealopts = nRealOpts
            )

        #self.inform = inform
        #self.hs = hs
        #self.x = x
        #self.y = pi # shadow prices
        #self.z = rc # reduced costs

        #status = self.inform
        if int(inform) == 0:
            inform = 'optimal'

        if output_as_dict:
            return { 'fluxes' : x, 'shadow_prices' : pi, 'reduced_costs' : rc, 'status' : inform, 'basis' : hs }
        else:
            return x, pi, rc, inform, hs

    def bisectmu(
        self, mumin = 0.0, mumax = 2.0, maxIter = 100, basis = None,
        tolerance = 1e-6, precision = 'quad', verbose = False
        ):
        """
        muopt, hs, xopt, cache = bisectmu(
            self, mumin = 0.0, mumax = 2.0, maxIter = 100,
            tolerance = 1e-6, precision = 'quad', verbose = False
            )

        Bisection to maximize the growth rate using qMINOS.
        Sequence of feasibility problems.
        """

        # basis (hs) intent(inout). Will generate first basis from Cold-start if hs=None
        #basis = None

        if verbose:
            #print('Iteration\t       Growth Rate\t Solution to check\tSolver Status')
            #print('---------\t------------------\t------------------\t-------------')
            print('Iteration\t Solution to check\tSolver Status')
            print('---------\t------------------\t-------------')

        # test mumax
        x_new, y_new, z_new, stat_new, hs_new = self.solvelp(mumax, basis, precision)

        def get_stat_msg(stat):
            if stat == "optimal":
                return "Optimal"
            if stat == 1:
                return "Not feasible"
            return str(stat)

        if stat_new == 'optimal' and verbose:
            print('{:s}\t{:.16f}\t{:s}'.format(
                        str(0).rjust(9), mumax, get_stat_msg(stat_new)))
            return mumax, x_new, y_new, z_new, basis, stat_new

        else:
            res = []
            for idx in range(1, maxIter + 1):
                # Just a sequence of feasibility checks
                muf = (mumin + mumax) / 2.
                # Retrieve evaluation from cache if it exists: golden section advantage
                #xopt, yopt, basis, stat = self.checkmu(muf, basis, precision)
                x_new, y_new, z_new, stat_new, hs_new = self.solvelp(muf, basis, precision)
                res.append([muf, x_new, y_new, z_new, hs_new, stat_new])

                if stat_new == 'optimal':
                    basis = hs_new
                    mumin = muf
                else:
                    mumax = muf

                if verbose:
                    #print('{:s}\t{:.16f}\t{:.16f}\t{:s}'.format(
                    print('{:s}\t{:.16f}\t{:s}'.format(
                        str(idx).rjust(9), muf, get_stat_msg(stat_new)))

                if abs(mumax - mumin) <= tolerance:# and stat_new == 'optimal':
                    valid_solutions = [ x for x in res if x[-1] == 'optimal' ]
                    if len(valid_solutions) != 0:
                        return valid_solutions[-1]
                    else:
                        return res[-1] # stat_new is not optimal

                if mumax <= tolerance:
                    return res[-1] # stat_new is not optimal

            # # Save feasible basis
            # self.feas_basis = basis

            return muf, x_new, y_new, z_new, basis, stat_new

    def varyme(self, mu_fixed, obj_inds0, obj_coeffs, basis = None, verbosity = False):
        """
        fva_result, fva_stats = varyme(self, mu_fixed, obj_inds0, obj_coeffs)

        rxns_fva0:  list of reactions to be varied (Reaction objects or ID)

        High-level interface for qvaryME (quad-prec FVA)
        12 Aug 2015: first version. Must fix bugs.
        13 Jul 2023: Adapted to new structure of arguments by Rodrigo Santibanez
        """

        m, n, ha, ka, ad, bld, bud, hs0, obj_inds = self.make_lp(mu_fixed, obj_inds0, obj_coeffs)

        hs = basis
        if hs is None:
            warm = False
            hs = hs0
        else:
            warm = True

        # Get MINOS options
        stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

        nVary = len(obj_inds)
        obj_vals = numpy.zeros(nVary)
        fva_stats = numpy.zeros(nVary, numpy.dtype('i4'))
        probname = 'varyme'

        qvaryME.qvaryme(
            fva_stats, probname, m, ha, ka, ad, bld, bud, hs, warm,
            obj_inds, obj_coeffs, obj_vals,
            stropts, intopts, realopts, intvals, realvals,
            nstropts = nStrOpts, nintopts = nIntOpts, nrealopts = nRealOpts
            )

        #return fva_result, fva_stats
        return obj_inds0, nVary, obj_vals
