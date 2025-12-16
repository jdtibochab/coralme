import collections
import functools
import logging
import pint
import pandas
import sympy

# configuration
log_format = '%(asctime)s %(funcName)s %(filename)s:%(lineno)d %(message)s' #%(clientip)-15s %(user)-8s

#https://stackoverflow.com/questions/36408496/python-logging-handler-to-append-to-list
#Here is a naive, non thread-safe implementation:
class ListHandler(logging.Handler): # Inherit from logging.Handler
	"""
	ListHandler class to handle prints and logs.

	"""
	def __init__(self, log_list, debug = False):
		# run the regular Handler __init__
		logging.Handler.__init__(self)
		# Our custom argument
		self.level = logging.WARNING
		self.formatter = log_format
		self.log_list = log_list
		self.debug = debug

	def emit(self, record):
		# WARNING: should match log_format
		try:
			self.log_list.append((record.asctime, record.funcName, record.filename, str(record.lineno), record.message))
		except:
			pass

	def print_and_log(msg):
		print(msg)
		logging.warning(msg)

	def log(msg):
		logging.warning(msg)

	def log_to_file(logger, debug, outfile):
		tmp = pandas.DataFrame(logger, columns = ['asctime', 'funcName', 'filename', 'lineno', 'message'])
		for idx, data in tmp.drop_duplicates(subset = 'message').iterrows():
			if debug:
				outfile.write('{:s} {:s} {:s}:{:s} {:s}\n'.format(data.iloc[0], data.iloc[1], data.iloc[2], str(data.iloc[3]), data.iloc[4]))
			else:
				outfile.write('{:s} {:s}\n'.format(data.iloc[0], data.iloc[4]))

class ExtendedQuantity(pint.Quantity):
	"""
	An extension of pint.Quantity that supports symbolic substitution
	with sympy.

	Provides:
	- subs() for replacing sympy symbols or symbolic quantities
	- Unit-aware arithmetic with substituted values
	- Rich comparison operators against quantities or scalars
	- Clean string representation

	Example:
		>>> t = sympy.Symbol("t")
		>>> q = ExtendedQuantity(t, "s")
		>>> q.subs(t, 5)
		5 second
	"""
	def __new__(cls, value, units=None):
		if isinstance(value, pint.Quantity):
			return super().__new__(cls, value.magnitude, value.units)
		return super().__new__(cls, value, units)

	def __init__(self, value, units=None):
		pass

	def subs(self, *args):
		"""Substitute symbolic variables in the attached expression.

		Example:
			>>> t = sympy.Symbol("t")
			>>> q = ExtendedQuantity(t, "s")
			>>> q.subs(t, 10)
			10 second
		"""
		if len(args) == 1 and isinstance(args[0], dict):
			substitutions = args[0]
		elif len(args) % 2 == 0:
			substitutions = {args[i]: args[i + 1] for i in range(0, len(args), 2)}
		else:
			raise TypeError("subs() takes either (dict) or (old, new) arguments.")

		magnitude_subs = {}
		resulting_units = self.units

		for key, value in substitutions.items():
			if isinstance(key, sympy.Symbol):
				symbol = key
				symbol_units = 1
			elif isinstance(key, pint.Quantity) and isinstance(key.magnitude, sympy.Symbol):
				symbol = key.magnitude
				symbol_units = key.units
			else:
				raise TypeError("subs() expects a dict or even number of arguments as (old, new) pairs.")

			if isinstance(value, pint.Quantity):
				magnitude_subs[symbol] = value.to(symbol_units).magnitude
				resulting_units *= symbol_units * value.units / symbol_units
			elif isinstance(float(value), float):
				magnitude_subs[symbol] = value
			else:
				raise TypeError("Substitution values must be float, int, or pint.Quantity.")

		new_expr = self.magnitude.subs(magnitude_subs)
		return ExtendedQuantity(new_expr, resulting_units)

	def _compare(self, other, op):
		"""Helper for comparisons."""
		if isinstance(other, ExtendedQuantity):
			return op(self, other.magnitude)
		elif isinstance(float(other), float):
			return op(self.magnitude, other)
		else:
			return NotImplemented

	def __lt__(self, other):
		"""Less-than comparison."""
		return self._compare(other, lambda x, y: x < y)

	def __le__(self, other):
		"""Less-than-or-equal comparison."""
		return self._compare(other, lambda x, y: x <= y)

	def __gt__(self, other):
		"""Greater-than comparison."""
		return self._compare(other, lambda x, y: x > y)

	def __ge__(self, other):
		"""Greater-than-or-equal comparison."""
		return self._compare(other, lambda x, y: x >= y)

	def __eq__(self, other):
		"""Equality comparison."""
		return self._compare(other, lambda x, y: x == y)

	def __ne__(self, other):
		"""Inequality comparison."""
		return self._compare(other, lambda x, y: x != y)

	def __repr__(self):
		"""Return string with magnitude and units."""
		return "{} {}".format(self.magnitude, self.units)

class DefaultParameters(dict):
	"""
	A dict subclass that normalizes all keys into sympy.Symbol objects
	with positive=True.

	Supports flexible initialization, access, and updates by string or
	symbol.

	Example:
		>>> params = DefaultParameters({"x": 5}, y=10)
		>>> params["x"]
		5
		>>> params["y"]
		10
	"""
	def __init__(self, *args, **kwargs):
		super().__init__()
		initial_data = dict(*args, **kwargs)
		for key, value in initial_data.items():
			self[key] = value

	def __setitem__(self, key, value):
		"""Ensure all keys are sympy.Symbol with positive=True."""
		if isinstance(key, str):
			key = sympy.Symbol(key, positive=True)
		elif isinstance(key, sympy.Symbol):
			key = sympy.Symbol(key.name, positive=True)
		super().__setitem__(key, value)

	def __getitem__(self, key):
		"""Normalize keys before retrieval."""
		if isinstance(key, str):
			key = sympy.Symbol(key, positive=True)
		elif isinstance(key, sympy.Symbol):
			key = sympy.Symbol(key.name, positive=True)
		return super().__getitem__(key)

	def get(self, key, default=None):
		"""Retrieve value for key (string or symbol)."""
		if isinstance(key, str):
			key = sympy.Symbol(key, positive=True)
		elif isinstance(key, sympy.Symbol):
			key = sympy.Symbol(key.name, positive=True)
		return super().get(key, default)

	def update(self, *args, **kwargs):
		"""Update dict with normalized keys."""
		new_data = dict(*args, **kwargs)
		converted_data = {
			(sympy.Symbol(k, positive=True) if isinstance(k, str) else sympy.Symbol(k.name, positive=True)): v
			for k, v in new_data.items()
		}
		super().update(converted_data)

	def query(self, substring):
		"""
		Return a new DefaultParameters instance containing only the
		entries whose key string contains `substring`.

		Example:
			>>> params = DefaultParameters({
			...     "Vmax_glc": 10,
			...     "Km_glc": 0.5,
			...     "Vmax_lac": 3
			... })
			>>> params.query("glc")
			DefaultParameters({Vmax_glc: 10, Km_glc: 0.5})
		"""
		substring = str(substring)
		filtered = {
			key: value
			for key, value in self.items()
			if substring in str(key)
		}
		return DefaultParameters(filtered)

class ScalableKeyDict(dict):
	"""
	A dict subclass supporting arithmetic on keys and values.

	Features:
	- Multiply/divide all values by a scalar
	- Shift keys by adding a scalar
	- Merge with other dicts, summing overlapping keys

	Example:
		>>> d = ScalableKeyDict({1: 10, 2: 20})
		>>> 2 * d
		{1: 20, 2: 40}
		>>> d + 1
		{2: 10, 3: 20}
	"""
	def __mul__(self, scalar):
		if not isinstance(scalar, (int, float)):
			raise TypeError("Can only multiply by int or float")
		return ScalableKeyDict({key: value * scalar for key, value in self.items()})

	def __rmul__(self, scalar):
		return self.__mul__(scalar)

	def __truediv__(self, scalar):
		return self.__mul__(1. / scalar)

	def __add__(self, other):
		if isinstance(other, (int, float)):
			return ScalableKeyDict({key + other: value for key, value in self.items()})
		elif isinstance(other, dict):
			result = ScalableKeyDict(self)
			for key, value in other.items():
				if key in result:
					result[key] += value
				else:
					result[key] = value
			return result
		else:
			raise TypeError("Can only add int, float, or dict")

	def __radd__(self, other):
		return self.__add__(other)

class MappableList(list):
	"""
	A list subclass with functional and mapping utilities.

	Features:
	- map(), filter(), reduce()
	- Attribute/key access with dot notation
	- Preserves MappableList type for operations

	Example:
		>>> lst = MappableList([1, 2, 3])
		>>> lst.map(lambda x: x*2)
		MappableList([2, 4, 6])
	"""
	def map(self, func):
		"""Apply a function to each element."""
		return MappableList([func(x) for x in self])

	def filter(self, func):
		"""Keep elements where func(element) is True."""
		return MappableList([x for x in self if func(x)])

	def reduce(self, func, initial=None):
		"""Reduce the list to a single value."""
		return functools.reduce(func, self, initial) if initial is not None else reduce(func, self)

	def __getattr__(self, attr):
		"""Access attributes or dict keys across elements."""
		values = []
		for x in self:
			if isinstance(x, dict) and attr in x:
				values.append(x[attr])
			elif hasattr(x, attr):
				values.append(getattr(x, attr))
			else:
				raise AttributeError(f"Element {x!r} has no attribute or key '{attr}'")
		return MappableList(values)

	def get(self, *attrs):
		"""Get one or more attributes/keys (supports dot notation)."""
		values = []
		for x in self:
			row = []
			for attr in attrs:
				val = self._get_nested(x, attr)
				row.append(val)
			values.append(tuple(row) if len(row) > 1 else row[0])
		return MappableList(values)

	def _get_nested(self, obj, attr):
		"""Recursively get nested attributes or keys using dot notation."""
		for part in attr.split('.'):
			if isinstance(obj, dict) and part in obj:
				obj = obj[part]
			elif hasattr(obj, part):
				obj = getattr(obj, part)
			else:
				raise AttributeError(f"Object {obj!r} has no attribute or key '{part}'")
		return obj

	def to_list(self):
		"""Convert to a plain Python list."""
		return list(self)

	def __getitem__(self, key):
		"""Preserve MappableList type for slices."""
		result = super().__getitem__(key)
		return MappableList(result) if isinstance(key, slice) else result

	def __repr__(self):
		"""String representation."""
		return f"MappableList({list(self)})"

class MCounter(collections.Counter):
	"""
	Extension of collections.Counter that supports:

	- Multiplication by int or float
	- Addition with other Counters
	- Sorted output by keys
	- Expansion into repeated lists

	Example:
		>>> c = MCounter({"a": 2, "b": 1})
		>>> 2 * c
		MCounter({'a': 4, 'b': 2})
		>>> c.expand()
		['a', 'a', 'b']
	"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def __mul__(self, other):
		return MCounter({k: other * v for k, v in self.items()})

	def __rmul__(self, other):
		return self * other

	def __add__(self, other):
		return MCounter(super().__add__(other))

	def sorted(self, reverse=False):
		"""Return a new MCounter sorted by keys."""
		return MCounter(dict(sorted(self.items(), key=lambda kv: kv[0], reverse=reverse)))

	def expand(self):
		"""Return a list with keys repeated by counts (rounded)."""
		return [k for k, v in self.items() for _ in range(round(v))]

	@classmethod
	def from_keys(cls, keys):
		"""
		Create an MCounter where each key appears with initial count 0.

		Example:
			>>> MCounter.from_keys(["a", "b", "c"])
			MCounter({'a': 0, 'b': 0, 'c': 0})
		"""
		return cls({key : 0 for key in keys})

class MEModelWrapper:
	"""
	A wrapper around a ME model that allows direct attribute-style access
	to metabolites, reactions, and other model elements.

	Example:
		>>> wrapper = MEModelWrapper(me_model)
		>>> wrapper.glc__D_c       # metabolite
		<Metabolite object>
		>>> wrapper.rRNA_genes     # model attribute
		{...}
	"""
	def __init__(self, me_model):
		"""Initialize with a ME model instance."""
		self._me_model = me_model

	def __getattr__(self, name):
		"""
		Provide attribute-style access. First checks wrapper attributes,
		then attempts to retrieve the attribute from the ME model.
		"""
		# Avoid overriding existing attributes
		if name in self.__dict__:
			return self.__dict__[name]

		# First, try to access as a normal attribute
		try:
			return getattr(self._me_model, name)
		except AttributeError:
			pass

		# If that fails, try me_model.get(name) (for metabolites/reactions)
		try:
			return self._me_model.get(name)
		except Exception as e:
			raise AttributeError(f"{name!r} not found in me_model") from e

	def __getitem__(self, key):
		"""Optional dict-style access to elements."""
		return self._me_model.get(key)

	def __dir__(self):
		"""Tab-completion: list wrapper attributes plus metabolite and reaction IDs."""
		base_dir = set(super().__dir__()) | set(dir(self._me_model))
		try:
			metabolite_ids = set(self._me_model.metabolites.list_attr("id"))
		except Exception:
			metabolite_ids = set()
		try:
			reaction_ids = set(self._me_model.reactions.list_attr("id"))
		except Exception:
			reaction_ids = set()
		return sorted(base_dir | metabolite_ids | reaction_ids)
