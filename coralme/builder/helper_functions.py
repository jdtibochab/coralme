#!/usr/bin/python3
import re
import copy
import sympy
import pandas

from ast import parse as ast_parse, Name, And, Or, BitOr, BitAnd, BoolOp, Expression, NodeTransformer

import cobra
import coralme

# from cobrame without changes
def get_base_complex_data(model, complex_id):
	"""
	If a complex is modified in a metabolic reaction it will not
	have a formation reaction associated with it. This function returns
	the complex data of the "base" complex, which will have the subunit
	stoichiometry of that complex
	"""

	# First try unmodified complex id
	try_1 = complex_id.split('_mod_')[0]
	if try_1 in model.process_data:
		return model.process_data.get_by_id(try_1)

	count = 0
	try_2 = complex_id.split('_mod_')[0] + '_'
	for i in model.process_data.query(try_2):
		if isinstance(i, coralme.core.processdata.ComplexData):
			count += 1
			data = i

	if count == 0:
		raise UserWarning('No base complex found for \'{:s}\'.'.format(complex_id))

	elif count > 1:
		raise UserWarning('More than one possible base complex found for \'{:s}\'.'.format(complex_id))

	return data

# Originally developed by JDTB, UCSD (2022)
def close_sink_and_solve(rxn_id):
	global _model
	global _muf
	model = copy.deepcopy(_model)
	model.reactions.get_by_id(rxn_id).bounds = (0, 0)
	model.optimize(max_mu = _muf, maxIter = 1)
	if not model.solution:
		return (rxn_id, False)
	else:
		return (rxn_id, True)

def change_reaction_id(model,old_id,new_id):
	import cobra
	old_rxn = model.reactions.get_by_id(old_id)
	rxn = cobra.Reaction(new_id)
	model.add_reactions([rxn])
	for k,v in old_rxn.metabolites.items():
		rxn.add_metabolites({k.id:v})
	rxn.bounds = old_rxn.bounds
	rxn.name = old_rxn.name
	rxn.subsystem = old_rxn.subsystem
	rxn.notes = old_rxn.notes
	model.remove_reactions([old_rxn])

def listify_gpr(expr,level = 0):
	import cobra

	if level == 0:
		return listify_gpr(cobra.core.gene.parse_gpr(str(expr))[0], level = 1)
	if isinstance(expr, cobra.core.gene.GPR):
		return listify_gpr(expr.body, level = 1) if hasattr(expr, "body") else ""
	elif isinstance(expr, Name):
		return expr.id
	elif isinstance(expr, BoolOp):
		op = expr.op
		if isinstance(op, Or):
			str_exp = list([listify_gpr(i, level = 1) for i in expr.values])
		elif isinstance(op, And):
			str_exp = tuple([listify_gpr(i, level = 1) for i in expr.values])
		return str_exp
	elif expr is None:
		return ""
	else:
		raise TypeError("unsupported operation  " + repr(expr))

def get_combinations(l_gpr):
	l_gpr = [[i] if isinstance(i,str) else i for i in l_gpr]
	return list(product(*l_gpr))

def get_chain(l_gpr):
	l_gpr = [[i] if isinstance(i,str) else i for i in l_gpr]
	return list(chain(*l_gpr))

def _expand_gpr(l_gpr):
	if isinstance(l_gpr,str):
		return l_gpr
	elif isinstance(l_gpr,list):
		return ' or '.join([expand_gpr(i) for i in l_gpr])
	elif isinstance(l_gpr,tuple):
		if len(l_gpr) == 1:
			return expand_gpr(l_gpr[0])
		elif any(isinstance(i,list) for i in l_gpr):
			return ' or '.join(chain([expand_gpr(i) for i in get_combinations(l_gpr)]))
		else:
			return ' and '.join(get_chain(l_gpr))

def generify_gpr(l_gpr,rxn_id,d={}):
	if isinstance(l_gpr,str):
		name = l_gpr
		return name,d
	elif isinstance(l_gpr,list):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d)
			l.append(n)
		base_name = 'generic_{}'.format(rxn_id)
		name = '{}_{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' or '.join(l)
		return name,d
	elif isinstance(l_gpr,tuple):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d)
			l.append(n)
		base_name = 'CPLX_{}'.format(rxn_id)
		name = '{}-{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' and '.join(l)
		return name,d

def print_check(i,l_gpr,T):
	print(i)
	print(l_gpr)
	print(T)
	print()

def get_tree(l_gpr,T={}):
	if isinstance(l_gpr,str):
		return l_gpr
	else:
		if isinstance(l_gpr,list):
			op = 'or'
		elif isinstance(l_gpr,tuple):
			op = 'and'
		T[op] = []
		for idx,i in enumerate(l_gpr):
			d = {}
			T[op].append(get_tree(i,T=d))
		return T

def append_graph(G,g):
	if G == '$':
		return g.copy()
	if isinstance(G,dict):
		for k,v in G.items():
			G[k] = append_graph(v,g)
		return G
def concatenate_graphs(L,r=[]):
	if r:
		for i in r:
			L = append_graph(L,i)
		return L
	elif isinstance(L,list):
		if len(L) == 1:
			return L[0]
		else:
			b = L[0]
			r = L[1:]
			L = concatenate_graphs(b,r)
		return L

def get_graph(T,G={}):
	if isinstance(T,str):
		G[T] = '$'
		return G
	elif isinstance(T,dict):
		if 'and' in T:
			l = []
			for i in T['and']:
				d = {}
				l.append(get_graph(i,d))
			d = concatenate_graphs(l)
			for k,v in d.items():
				G[k] = v
			return G
		elif 'or' in T:
			for i in T['or']:
				G = get_graph(i,G)
		return G

def traverse_graph(G,L = [], C = []):
	if G == '$':
		C.append(L)
		return L,C
	if isinstance(G,dict):
		for k,v in G.items():
			l = L + [k]
			l,C = traverse_graph(v,l,C)
		return L,C

def expand_gpr(rule):
	l = listify_gpr(rule)
	T = get_tree(l,T={})
	G = get_graph(T,G={})
	return traverse_graph(G,L=[],C=[])[1]

def generify_gpr(l_gpr,rxn_id,d={}):
	if isinstance(l_gpr,str):
		name = l_gpr
		return name,d
	elif isinstance(l_gpr,list):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d)
			l.append(n)
		base_name = 'generic_{}'.format(rxn_id)
		name = '{}_{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' or '.join(l)
		return name,d
	elif isinstance(l_gpr,tuple):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d)
			l.append(n)
		base_name = 'CPLX_{}'.format(rxn_id)
		name = '{}-{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' and '.join(l)
		return name,d

def process_rule_dict(n,rule_dict,gene_dict,protein_mod):
	corrected_ids = {}
	for cplx,rule in rule_dict.items():
		cplx_id = 0
		if 'CPLX' in cplx:
			rule_gene_list = rule.split(" and ")
			identified_genes = rule_gene_list
			cplx_id = find_match(gene_dict,identified_genes)
		if not cplx_id:
			cplx_id = cplx
		corrected_ids[cplx] = cplx_id
	corrected_rule_dict = {}

	for cplx,rule in rule_dict.items():
		if cplx in corrected_ids:
			cplx_id = corrected_ids[cplx]
		else:
			cplx_id = cplx
		# if cplx_id in protein_mod
		if cplx_id in protein_mod["Core_enzyme"].values:
			cplx_mod_id = protein_mod[
				protein_mod["Core_enzyme"].str.contains(cplx_id)
			].index[0]
			if "Oxidized" in cplx_mod_id:
				cplx_mod_id = cplx_mod_id.split("_mod_Oxidized")[0]
			if corrected_ids[n] == cplx_id:
				rule = corrected_ids.pop(n)
				corrected_ids[n] = cplx_mod_id
			cplx_id = cplx_mod_id
		for c,cid in corrected_ids.items():
			regex = '{}(?!\d)'
			corrected_rule_dict[cplx_id] = re.sub(regex.format(c), cid, rule)
			rule = corrected_rule_dict[cplx_id]
	return corrected_ids[n],corrected_rule_dict

def find_match(d,items):
	for c, cg in d.items():
		if not cg: continue
		cg = [re.findall('.*(?=\(\d*\))', g)[0] for g in cg.split(' AND ')]
		if set(cg) == set(items):
			return c
	return 0

# Originally developed by JDTB@UCSD, 2022
# Modified by RSP@UCSD, 2022
def process_model(model, growth_key = sympy.Symbol('mu', positive = True)):
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

				# Replace 'growth_key' if model is a ME-Model
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

def find_gaps(model, growth_key = sympy.Symbol('mu', positive = True)):
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

def add_exchange_reactions(me, metabolites, prefix = 'SK_'):
	for met in metabolites:
		rxn_id = prefix + met
		if rxn_id not in me.reactions:
			r = coralme.core.reaction.MEReaction(rxn_id)
			me.add_reaction(r)
			r.add_metabolites({ met: -1 })
		else:
			r = me.reactions.get_by_id(rxn_id)
		r.bounds = (-10, 1000)
		#print(r.id,r.lower_bound,r.upper_bound,r.reaction)
	return me

def brute_force_check(me, metabolites_to_add, growth_key_and_value):
	print('  '*6 + 'Adding sink reactions for {:d} metabolites'.format(len(metabolites_to_add)))
	add_exchange_reactions(me, metabolites_to_add)

	if me.feasibility(keys = growth_key_and_value):
		pass
	else:
		return False

	rxns = []
	for idx, flux in me.solution.fluxes.items():
		if idx.startswith('SK_') and idx.split('SK_')[1] in metabolites_to_add:
			if abs(flux) > 0:
				rxns.append(idx)
			else:
				#print('Closing {}'.format(idx))
				me.reactions.get_by_id(idx).bounds = (0, 0)

	print('  '*6 + 'Sink reactions shortlisted to {:d} metabolites:'.format(len(rxns)))

	# reaction ID : position in the model.reactions DictList object
	ridx = { k:v for k,v in me.reactions._dict.items() if k in rxns }

	# populate with stoichiometry
	Sf, Se, lb, ub, b, c, cs, atoms = me.construct_lp_problem(keys = growth_key_and_value)

	res = []
	msg = 'Processed: {:s}/{:d}, Gaps: {:d}. The ME-Model is {:s}feasible if {:s} is closed.'
	for idx, (rxn, pos) in enumerate(ridx.items()):
		lb[pos] = 0
		ub[pos] = 0
		if me.feasibility(keys = growth_key_and_value, **{'lp' : [Sf, dict(), lb, ub, b, c, cs, set()]}):
			res.append(False)
			print('  '*6, msg.format(str(idx+1).rjust(len(str(len(ridx))), ' '), len(ridx), len([ x for x in res if x ]), '', rxn))
		else:
			lb[pos] = -1000
			ub[pos] = +1000
			res.append(True)
			print('  '*6, msg.format(str(idx+1).rjust(len(str(len(ridx))), ' '), len(ridx), len([ x for x in res if x ]), 'not ', rxn))

	return [ y for x,y in zip(res, rxns) if x ]