#!/usr/bin/python3

# ME model troubleshooting
# Originally developed by JDTB@UCSD, 2022
# Modified by RSP@UCSD, 2022

import coralme
import cobra
import sympy
import pandas
import logging
import sys
log = logging.getLogger(__name__)

def process_model(model, growth_key = sympy.Symbol('mu', positive = True)):
	"""
	Get a dictionary containing information on whether a metabolite has
	producing or consuming reactions. This is used to find gaps.
	"""
	dct = {}
	for met in model.metabolites:
		filter1 = type(met) == cobra.core.metabolite.Metabolite or type(met) == coralme.core.component.Metabolite
		filter2 = met.id.startswith('trna')
		filter3 = met.id.endswith('trna_c')

		if filter1 and not filter2 and not filter3:
			t = { 'c' : set(), 'p' : set() }
			#seen = [] #?
			for rxn in met.reactions:
				if rxn.id.startswith('BIOMASS_'):
					continue

				lb, ub = rxn.lower_bound, rxn.upper_bound

				# Replace 'growth_key' if model is a ME-model
				if hasattr(lb, 'subs'):
					lb = lb.subs(growth_key, 1.)
				if hasattr(ub, 'subs'):
					ub = ub.subs(growth_key, 1.)

				coeff = rxn.metabolites[met]
				if hasattr(coeff, 'subs'):
					coeff = coeff.subs(growth_key, 1.)

				pos = 1 if coeff > 0 else -1
				rev = 1 if lb < 0 else 0
				fwd = 1 if ub > 0 else 0
				if pos*fwd == -1 or pos*rev == +1:
					t['c'].add(rxn.id)
				if pos*fwd == +1 or pos*rev == -1:
					t['p'].add(rxn.id)
			dct[met.id] = t
	return dct

def add_exchange_reactions(me, metabolites, prefix = 'SK_'):
	"""Add exchange/sink reaction to the model"""
	rxns = []
	for met in metabolites:
		rxn_id = prefix + met
		if rxn_id not in me.reactions:
			r = coralme.core.reaction.MEReaction(rxn_id)
			me.add_reactions([r])
			r.add_metabolites({ met: -1 })
		else:
			r = me.reactions.get_by_id(rxn_id)
		r.bounds = (-10, 1000)
		rxns.append(r)
		#print(r.id,r.lower_bound,r.upper_bound,r.reaction)
	return rxns

def find_gaps(model, growth_key = sympy.Symbol('mu', positive = True)):
	"""Find gaps in the model"""
	g = {}
	dct = process_model(model, growth_key = growth_key)
	for met, t in dct.items():
		# not producing, not consuming, not uerever
		g[met] = { 'p' : 0, 'c' : 0, 'u' : 0 }
		if not t['c']:
			g[met]['c'] = 1
		if not t['p']:
			g[met]['p'] = 1
		if len(t['c']) == 1 and t['c'] == t['p']:
			g[met]['u'] = 1
	df = pandas.DataFrame.from_dict(g).T
	df = df[df.any(axis = 1)]
	df = df.sort_index()
	return df

def find_issue(query,d,msg = ''):
	"""Retrieve any warning message in curation notes associated with an ID"""
	if isinstance(d,dict):
		if 'msg' in d:
			msg = d['msg']
			if 'triggered_by' in d:
				trigger = d['triggered_by']
				find_issue(query,trigger,msg=msg)
		else:
			for k,v in d.items():
				find_issue(query,v,msg=msg)
	elif isinstance(d,list):
		for i in d:
			find_issue(query,i,msg=msg)
	elif isinstance(d,str):
		if query == d:
			print(msg)
	else:
		raise TypeError("unsupported type  " + type(d))

def fill_builder(b,fill_with='CPLX_dummy',key=None,d=None,fieldname=None,warnings=None):
	"""Fill empty fields in builder with CPLX dummy"""
	if isinstance(b,coralme.builder.main.MEBuilder):
		for i in dir(b.org):
			if i[0] == '_':
				continue
			attr = getattr(b.org,i)
			if not isinstance(attr,dict):
				continue
			fill_builder(attr,fill_with=fill_with,fieldname=i,warnings = warnings)
	elif isinstance(b,dict):
		for k,v in b.items():
			fill_builder(v,key=k,d=b,fill_with=fill_with,fieldname=fieldname,warnings=warnings)
	elif isinstance(b,list):
		include_keys = ['enzymes','proteins','enzyme','protein','machine']
		for ik in include_keys:
			if key in ik:
				if not b:
					d[key] = ['CPLX_dummy']
	elif isinstance(b,str):
		include_keys = ['enzymes','proteins','enzyme','protein','machine']
		for ik in include_keys:
			if key in ik:
				if not b:
					d[key] = 'CPLX_dummy'
	else:
		pass

def gap_find(me_model,de_type = None):
	"""Find and classify gaps in the model"""

	logging.warning('  '*5 + 'Finding gaps in the ME-model...')
	me_gaps = coralme.builder.troubleshooting.find_gaps(me_model, growth_key = me_model.mu)

	if de_type == 'me_only':
		logging.warning('  '*5 + 'Finding gaps from the M-model only...')
		m_gaps = coralme.builder.troubleshooting.find_gaps(me_model.gem)
		idx = list(set(me_gaps.index) - set(m_gaps.index))
	else:
		idx = list(set(me_gaps.index))
	new_gaps = me_gaps.loc[idx]

	filt1 = new_gaps['p'] == 1
	filt2 = new_gaps['c'] == 1
	filt3 = new_gaps['u'] == 1

	deadends = list(new_gaps[filt1 | filt2 | filt3].index)
	deadends = sorted([ x for x in deadends if 'biomass' not in x if not x.endswith('_e') ])

	logging.warning('  '*5 + '{:d} metabolites were identified as deadends.'.format(len(deadends)))
	for met in deadends:
		name = me_model.metabolites.get_by_id(met).name
		logging.warning('  '*6 + '{:s}: {:s}'.format(met, 'Missing metabolite in the M-model.' if name == '' else name))
	return deadends

def gap_fill(me_model, deadends = [], growth_key_and_value = { sympy.Symbol('mu', positive = True) : 0.1 }, met_types = 'Metabolite'):
	"""Add sink reactions of gap metabolites to the model"""
	if sys.platform in ['win32', 'darwin']:
		self.me_model.get_solution = self.me_model.optimize_windows
		self.me_model.get_feasibility = self.me_model.feas_gurobi
	else:
		self.me_model.get_solution = self.me_model.optimize
		self.me_model.get_feasibility = self.me_model.feasibility

	if len(deadends) != 0:
		logging.warning('  '*5 + 'Adding a sink reaction for each identified deadend metabolite...')
		coralme.builder.troubleshooting.add_exchange_reactions(me_model, deadends, prefix='TS_')
	else:
		logging.warning('  '*5 + 'Empty set of deadends metabolites to test.')
		return None

	logging.warning('  '*5 + 'Optimizing gapfilled ME-model...')

	if me_model.get_feasibility(keys = growth_key_and_value):
		#logging.warning('  '*5 + 'The ME-model is feasible.')
		logging.warning('  '*5 + 'Gapfilled ME-model is feasible with growth rate {:g} 1/h.'.format(list(growth_key_and_value.values())[0]))
		return True
	else:
		#logging.warning('  '*5 + 'The ME-model is not feasible.')
		logging.warning('  '*5 + 'Provided set of sink reactions for deadend metabolites does not allow growth.')
		return False

def brute_force_check(me_model, metabolites_to_add, growth_key_and_value):
	"""
	Iteratively search for minimal set of metabolites that are needed as
	sinks to allow for growth. This function searches by batches of
	different types of metabolites.
	"""
	if sys.platform == 'win32':
		me_model.get_solution = me_model.opt_gurobi
		me_model.get_feasibility = me_model.feas_gurobi
	else:
		me_model.get_solution = me_model.optimize
		me_model.get_feasibility = me_model.feasibility

	logging.warning('  '*5 + 'Adding sink reactions for {:d} metabolites...'.format(len(metabolites_to_add)))
# 	existing_sinks = [r.id for r in me_model.reactions.query('^TS_')]
	sk_rxns = coralme.builder.troubleshooting.add_exchange_reactions(me_model, metabolites_to_add, prefix='TS_')

	if me_model.get_feasibility(keys = growth_key_and_value):
		pass
	else:
		logging.warning('  '*5 + 'Provided metabolites through sink reactions cannot recover growth. Proceeding to next set of metabolites.')
		return metabolites_to_add, [], False

	rxns = []
	rxns_to_drop = []
# 	rxns_to_append = []
# 	for idx, flux in me_model.solution.fluxes.items():
	for r in sk_rxns:
		idx = r.id
		flux = me_model.solution.fluxes[idx]
		if idx.startswith('TS_') and idx.split('TS_')[1] in metabolites_to_add:
# 			if r.id in existing_sinks:
# 				rxns_to_append.append(idx)
# 				continue
			if abs(flux) > 0:
				rxns.append(idx)
			else:
				#logging.warning('Closing {}'.format(idx))
				rxns_to_drop.append(idx)
				me_model.reactions.get_by_id(idx).bounds = (0, 0)

	logging.warning('  '*6 + 'Sink reactions shortlisted to {:d} metabolites.'.format(len(rxns)))

	# reaction_id:position in the model.reactions DictList object
# 	rxns = rxns + rxns_to_append# Try present SKs the last.
# 	logging.warning('  '*6 + 'Will try a total of {:d} metabolites including previous iterations:'.format(len(rxns)))
	ridx = []
	for r in rxns:
		ridx.append((r,me_model.reactions._dict[r]))
# 	ridx = { k:v for k,v in me_model.reactions._dict.items() if k in rxns }

	# populate with stoichiometry
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas, Lr, Lm = me_model.construct_lp_problem()

	if lambdas is None:
		Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, Se, lb, ub, growth_key_and_value, atoms)
	else:
		Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, lambdas, lb, ub, growth_key_and_value, atoms)

	res = []
	msg = 'Processed: {:s}/{:d}, Gaps: {:d}. The ME-model is {:s}feasible if {:s} is closed.'
	for idx, (rxn, pos) in enumerate(ridx):
		lb[pos] = 0
		ub[pos] = 0
		if me_model.get_feasibility(keys = growth_key_and_value, **{'lp' : [Sf, dict(), lb, ub, b, c, cs, set(), lambdas, Lr, Lm]}):
			res.append(False)
			logging.warning('{:s} {:s}'.format('  '*6, msg.format(str(idx+1).rjust(len(str(len(ridx)))), len(ridx), len([ x for x in res if x ]), '', rxn)))
		else:
			lb[pos] = -1000
			ub[pos] = +1000
			res.append(True)
			logging.warning('{:s} {:s}'.format('  '*6, msg.format(str(idx+1).rjust(len(str(len(ridx)))), len(ridx), len([ x for x in res if x ]), 'not ', rxn)))

	bf_gaps = [ y for x,y in zip(res, rxns) if x ] # True
	no_gaps = [ y for x,y in zip(res, rxns) if not x ] + rxns_to_drop

	return bf_gaps, no_gaps, True

def get_mets_from_type(me_model,met_type):
	"""Get metabolites by type relevant to gap filling the model"""
	if met_type[1] == 'User guesses':
		return set(met_type[0])
	elif met_type == 'ME-Deadends':
		return set(coralme.builder.troubleshooting.gap_find(me_model,de_type='me_only'))
	elif met_type == 'All-Deadends':
		return set(coralme.builder.troubleshooting.gap_find(me_model))
	elif met_type == 'Cofactors':
		return set(get_cofactors_in_me_model(me_model))
	else:
		mets = set()
		for met in me_model.metabolites:
			filter1 = type(met) == getattr(coralme.core.component, met_type)
			filter2 = met.id.startswith('trna')
			filter3 = met.id.endswith('trna_c')
			filter4 = met.id.endswith('_e')
			if filter1 and not filter2 and not filter3 and not filter4:
				mets.add(met.id)
		return mets

def _append_metabolites(mets,new_mets):
	"""Merge metabolite lists"""
	return mets + [m for m in new_mets if m not in mets]

def brute_check(me_model, growth_key_and_value, met_type, skip = set(), history = dict()):
	"""Remove metabolites from our heuristics and call the brute force search algorithm"""
	mets = get_mets_from_type(me_model,met_type)
	if met_type == 'Metabolite':
		#remove from the metabolites to test that are fed into the model through transport reactions
		medium = set([ '{:s}_c'.format(x[3:-2]) for x in me_model.gem.medium.keys() ])
		mets = set(mets).difference(medium)
		# filter out manually
		mets = set(mets).difference(set(['ppi_c', 'ACP_c', 'h_c']))
		mets = set(mets).difference(set(['adp_c', 'amp_c', 'atp_c']))
		mets = set(mets).difference(set(['cdp_c', 'cmp_c', 'ctp_c']))
		mets = set(mets).difference(set(['gdp_c', 'gmp_c', 'gtp_c']))
		mets = set(mets).difference(set(['udp_c', 'ump_c', 'utp_c']))
		mets = set(mets).difference(set(['dadp_c', 'dcdp_c', 'dgdp_c', 'dtdp_c', 'dudp_c']))
		mets = set(mets).difference(set(['damp_c', 'dcmp_c', 'dgmp_c', 'dtmp_c', 'dump_c']))
		mets = set(mets).difference(set(['datp_c', 'dctp_c', 'dgtp_c', 'dttp_c', 'dutp_c']))
		mets = set(mets).difference(set(['nad_c', 'nadh_c', 'nadp_c', 'nadph_c']))
		mets = set(mets).difference(set(['5fthf_c', '10fthf_c', '5mthf_c', 'dhf_c', 'methf_c', 'mlthf_c', 'thf_c']))
		mets = set(mets).difference(set(['fad_c', 'fadh2_c', 'fmn_c']))
		mets = set(mets).difference(set(['coa_c']))
	mets = set(mets).difference(skip)
	if met_type[1] == 'User guesses':
		history['User guesses'] = mets
	else:
		history[met_type] = mets

	mets_to_check = []
	for k,v in history.items():
		mets_to_check = _append_metabolites(mets_to_check,v)
	return history,coralme.builder.troubleshooting.brute_force_check(me_model,
															  mets_to_check[::-1],
															  growth_key_and_value)

def get_cofactors_in_me_model(me):
	"""Get metabolites that work as cofactors in the model"""
	cofactors = set()
	for i in me.process_data.query('^mod_'):
		for k,v in i.stoichiometry.items():
			if not me.metabolites.has_id(k):
				continue
			if v < 0:
				cofactors.add(k)
	return list(cofactors)
