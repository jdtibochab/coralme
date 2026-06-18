import copy
import sympy
import coralme
import collections

def stringify(element, number):
	number = int(number)
	return element if number == 1 else (element + str(number).rstrip('.'))

def elements_to_formula(obj, elements):
	#obj.formula = ''.join(stringify(e, n) for e, n in sorted(iteritems(elements)))
	obj.formula = ''.join(stringify(e, n) for e, n in sorted(elements.items()) if n != 0)

def get_charge_from_process_data(reaction, process_data):
	"""
	If a modification is required to form a functioning macromolecule,
	update the element dictionary accordingly.
	"""

	#for sub, count in iteritems(process_data.subreactions):
	charge = 0
	for sub, count in process_data.subreactions.items():
		sub_obj = reaction._model.process_data.get_by_id(sub)
		charge += sub_obj.charge_contribution * count
	return charge

def get_elements_from_process_data(reaction, process_data, elements):
	"""
	If a modification is required to form a functioning macromolecule,
	update the element dictionary accordingly.
	"""

	#for sub, count in iteritems(process_data.subreactions):
	for sub, count in process_data.subreactions.items():
		sub_obj = reaction._model.process_data.get_by_id(sub)
		for e, n in sub_obj.element_contribution.items():
			elements[e] += n * count
	return elements