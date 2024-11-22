#!/usr/bin/python3
import cobra
from ast import parse as ast_parse, Name, And, Or, BitOr, BitAnd, BoolOp, Expression, NodeTransformer
import re

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
def get_size(G):
	return len(re.findall("\$",str(G)))
def get_graph(T,G={},length=1,threshold=100):
	#print(1, G,length)
	if G == "STOP":
		return "STOP",length
	if isinstance(T,str):
		if T in G:
			T = T + '_REPETITIONMARK_' + str(len(G))
		G[T] = '$'
		if length > threshold:
			G = "STOP"
		return G,1
	elif isinstance(T,dict):
		if 'and' in T:
			l = []
			for i in T['and']:
				d = {}
				g,_length = get_graph(i,d,threshold=threshold,length=length)
				if g == "STOP":
					return "STOP",length
				#print(f"{length} * {_length} = {length*_length}")
				length = length*_length
				l.append(g)
			d = concatenate_graphs(l)
			for k,v in d.items():
				if k in G:
					k = k + '_REPETITIONMARK_' + str(len(G))
				G[k] = v
			if length > threshold:
				G = "STOP"
			return G,length
		elif 'or' in T:
			for i in T['or']:
				G,_length = get_graph(i,G,threshold=threshold,length=length)
				if length > threshold:
					G = "STOP"
		#print(get_size(G))
		if G == "STOP":
			return G,length
		length = get_size(G)
		if length > threshold:
			G = "STOP"
		return G,length

def traverse_graph(G,L = [], C = []):
	if G == '$':
		C.append(L)
		return L,C
	if isinstance(G,dict):
		for k,v in G.items():
			k = k.split('_REPETITIONMARK_')[0]
			l = L + [k]
			l,C = traverse_graph(v,l,C)
		return L,C

def expand_gpr(rule,threshold=100):
	l = listify_gpr(rule)
	T = get_tree(l,T={})
	G,_ = get_graph(T,G={},threshold=threshold)
	if G == "STOP":
		return G
	return traverse_graph(G,L=[],C=[])[1]

def generify_gpr(l_gpr,rxn_id,d={},generic_gene_dict={}):
	if isinstance(l_gpr,str):
		name = l_gpr
		return name,d
	elif isinstance(l_gpr,list):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d,generic_gene_dict=generic_gene_dict)
			l.append(n)
		existing_generic = find_match(generic_gene_dict,l)
		if existing_generic:
			name = existing_generic
		else:
			base_name = 'generic_{}'.format(rxn_id)
			name = '{}_{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' or '.join(l)
		return name,d
	elif isinstance(l_gpr,tuple):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d,generic_gene_dict=generic_gene_dict)
			l.append(n)
		base_name = 'CPLX_{}'.format(rxn_id)
		name = '{}-{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' and '.join(l)
		return name,d

def listify_gpr(expr,level = 0,length=0):
	"""
	Modified from COBRApy
	"""
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
        if isinstance(cg,str):
            cg = [re.findall('.*(?=\(\d*\))', g)[0] for g in cg.split(' AND ')]
        if set(cg) == set(items):
            return c
    return 0
