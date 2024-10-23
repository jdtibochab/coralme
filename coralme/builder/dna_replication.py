import math
import numpy
import scipy

import coralme

# Experimental Data
#gr_data_doublings_per_hour = [0, 0.6, 1.0, 1.5, 2.0, 2.5]

# First point is mass 1 genome in 0.4 um^3 at density
#percent_dna_data = [0.0592, 0.0512, 0.0330, 0.0252, 0.0222, 0.0208]

def percent_dna_template_function(params, gr):
	[g_p_gdw_0, g_per_gdw_inf, b, d] = params
	c = g_per_gdw_inf
	a = g_p_gdw_0 - g_per_gdw_inf
	g_p_gdw = (-a * gr ** d) / (b + gr ** d) + a + c
	return g_p_gdw

def optimize_dna_function(gr, percent_dna):
	params = numpy.array([0.9, 0.3, 0.2, 1.0])

	def _minimization_function(params, gr, percent_dna):
		return percent_dna_template_function(params, gr) - percent_dna
	a = scipy.optimize.leastsq(_minimization_function, params, args = (gr, percent_dna))
	return a[0]

def get_dna_mw_no_ppi_dict(model):
	"""
	Return the molecular weight of dna component with the diphosphate group
	removed.
	"""
	default = {
		'dctp': 467.156923,
		'dgtp': 507.181023,
		'datp': 491.181623,
		'dttp': 482.168263,
		'ppi' : 177.97508200000001
		}

	dna_mw_no_ppi = {}
	ppi_mw = model.metabolites.ppi_c.formula_weight
	if ppi_mw == 0:
		ppi_mw = default['ppi']

	for dna in ['dctp', 'dgtp', 'datp', 'dttp']:
		dna_mw = model.metabolites.get_by_id(dna + '_c').formula_weight
		if dna_mw == 0:
			dna_mw = default[dna]

		dna_mw_no_ppi[dna] = dna_mw - ppi_mw

	return dna_mw_no_ppi

def return_gr_dependent_dna_demand(model, gc_fraction, percent_dna_data, gr_data_doublings_per_hour):
	"""
	Returns dNTP coefficients and lower/upper bounds of DNA_replication
	reaction
	"""
	if len(percent_dna_data) == len(gr_data_doublings_per_hour):
		pass
	else:
		raise Exception('The \'percent DNA data\' and \'growth data have different\' lengths.')

	gr_data = [m * math.log(2) for m in gr_data_doublings_per_hour]
	fit_params = optimize_dna_function(gr_data, percent_dna_data)
	dna_g_per_g = percent_dna_template_function(fit_params, model.mu)  # gDNA / gDW

	# average dinucleotide molecular weight
	dna_mw = get_dna_mw_no_ppi_dict(model)
	dntp_mw = (gc_fraction * (dna_mw['dctp'] + dna_mw['dgtp'])) / 2
	dntp_mw += ((1 - gc_fraction) * (dna_mw['datp'] + dna_mw['dttp'])) / 2

	# 1 / (gDNA / mol) * (1000 mmol / 1 mol)
	mmol_dntps_per_gram_dna = 1 / dntp_mw * 1000

	# (mmol / gDNA) * (gDNA / gDW)
	mmol_dntp_per_gdw = dna_g_per_g * mmol_dntps_per_gram_dna

	# lower and upper bound
	mmol_dntp_per_gdw_per_hr = mmol_dntp_per_gdw * model.mu

	return mmol_dntp_per_gdw_per_hr
