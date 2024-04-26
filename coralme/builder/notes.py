#!/usr/bin/python3

import collections
import json
import os

def dict_to_defaultdict(dct):
	"""Cast a dict to a defaultdict"""
	return collections.defaultdict(lambda: [], dct)

def save_curation_notes(curation_notes,filepath):
	"""Save curation notes as JSON"""
	file = open(filepath,'w')
	file.write(json.dumps(curation_notes, indent=4))
	file.close()

def load_curation_notes(filepath):
	"""Load curation notes from a JSON"""
	if not os.path.isfile(filepath):
		return collections.defaultdict(list)
	with open(filepath) as json_file:
		return json.load(json_file,object_hook=dict_to_defaultdict)

def publish_curation_notes(curation_notes,filepath):
	"""Print the curation notes to a txt file"""
	file = open(filepath,'w')
	for k,v in curation_notes.items():
		file.write('\n')
		for w in v:
			file.write('{} {}@{} {}\n'.format('#'*20,w['importance'],k,'#'*20))
			file.write('{} {}\n'.format('*'*10,w['msg']))
			if 'triggered_by' in w:
				file.write('The following items triggered the warning:\n')
				for i in w['triggered_by']:
					if isinstance(i,dict):
						file.write('\n')
						file.write(json.dumps(i))
						file.write('\n')
					else:
						file.write(i + '\n')
			file.write('\n{}Solution:\n{}\n\n'.format('*'*10,w['to_do']))
		file.write('\n\n')
	file.close()
