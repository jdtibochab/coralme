amino_acids = {
	'A': 'ala__L',
	'R': 'arg__L',
	'N': 'asn__L',
	'D': 'asp__L',
	'C': 'cys__L',
	'E': 'glu__L',
	'Q': 'gln__L',
	'G': 'gly',
	'H': 'his__L',
	'I': 'ile__L',
	'L': 'leu__L',
	'K': 'lys__L',
	'M': 'met__L',
	'F': 'phe__L',
	'P': 'pro__L',
	'S': 'ser__L',
	'T': 'thr__L',
	'W': 'trp__L',
	'Y': 'tyr__L',
	'V': 'val__L',
	'U': 'ser__L',  # precursor to selenocysteine
	'X': 'gly'
	}

amino_acids_3letters = {
	'A': 'Ala',
	'R': 'Arg',
	'N': 'Asn',
	'D': 'Asp',
	'C': 'Cys',
	'E': 'Glu',
	'Q': 'Gln',
	'G': 'Gly',
	'H': 'His',
	'I': 'Ile',
	'L': 'Leu',
	'K': 'Lys',
	'M': 'Met',
	'F': 'Phe',
	'P': 'Pro',
	'S': 'Ser',
	'T': 'Thr',
	'W': 'Trp',
	'Y': 'Tyr',
	'V': 'Val',
	'U': 'Sec',
	}

amino_acids_fullname = {
	'A': 'Alanine',
	'R': 'Arginine',
	'N': 'Asparagine',
	'D': 'Aspartate',
	'C': 'Cysteine',
	'G': 'Glycine', # replace first to avoid replace it again
	'E': 'Glutamate',
	'Q': 'Glutamine',
	'H': 'Histidine',
	'I': 'Isoleucine',
	'L': 'Leucine',
	'K': 'Lysine',
	'M': 'Methionine',
	'P': 'Proline', # replace first to avoid replace it again
	'F': 'Phenylalanine',
	'S': 'Serine',
	'T': 'Threonine',
	'W': 'Tryptophan',
	'Y': 'Tyrosine',
	'V': 'Valine',
	'U': 'Selenocysteine'
	}

# We use now Biopython to get the codon table specific for the organism
#codon_table = {
	#'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	#'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	#'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
	#'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
	#'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	#'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	#'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	#'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	#'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	#'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	#'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	#'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	#'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	#'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	#'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	#'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
	#}

#codon_table_3letters = {
	#k:amino_acids[v].split('_')[0].capitalize()
	#if v != '*' else 'STOP'
	#for k,v in codon_table.items()
	#}

#transcription_table = {
	#'A': 'utp_c',
	#'T': 'atp_c',
	#'C': 'gtp_c',
	#'G': 'ctp_c'
	#}

#import random
#ambiguous_nucleotides = {
	#'N': transcription_table[random.choice(['A', 'T', 'G', 'C'])],
	#'R': transcription_table[random.choice(['A', 'G'])],
	#'Y': transcription_table[random.choice(['T', 'C'])],
	#'K': transcription_table[random.choice(['G', 'T'])],
	#'M': transcription_table[random.choice(['T', 'C'])],
	#'S': transcription_table[random.choice(['C', 'G'])],
	#'W': transcription_table[random.choice(['A', 'T'])],
	#'B': transcription_table[random.choice(['C', 'G', 'T'])],
	#'D': transcription_table[random.choice(['A', 'G', 'T'])],
	#'H': transcription_table[random.choice(['A', 'C', 'T'])],
	#'V': transcription_table[random.choice(['A', 'C', 'G'])]
	#}

#transcription_table.update(ambiguous_nucleotides)

transcription_table = {
	'c' : { 'A': 'atp_c', 'T': 'utp_c', 'C': 'ctp_c', 'G': 'gtp_c' }, # cytosol
	'n' : { 'A': 'atp_n', 'T': 'utp_n', 'C': 'ctp_n', 'G': 'gtp_n' }, # nucleus
	'm' : { 'A': 'atp_m', 'T': 'utp_m', 'C': 'ctp_m', 'G': 'gtp_m' }, # mitochondria
	'h' : { 'A': 'atp_h', 'T': 'utp_h', 'C': 'ctp_h', 'G': 'gtp_h' }, # chloroplast
	}

import random
def ambiguous_nucleotides(compartment):
	return {
		'N': transcription_table[compartment][random.choice(['A', 'T', 'G', 'C'])],
		'R': transcription_table[compartment][random.choice(['A', 'G'])],
		'Y': transcription_table[compartment][random.choice(['T', 'C'])],
		'K': transcription_table[compartment][random.choice(['G', 'T'])],
		'M': transcription_table[compartment][random.choice(['T', 'C'])],
		'S': transcription_table[compartment][random.choice(['C', 'G'])],
		'W': transcription_table[compartment][random.choice(['A', 'T'])],
		'B': transcription_table[compartment][random.choice(['C', 'G', 'T'])],
		'D': transcription_table[compartment][random.choice(['A', 'G', 'T'])],
		'H': transcription_table[compartment][random.choice(['A', 'C', 'T'])],
		'V': transcription_table[compartment][random.choice(['A', 'C', 'G'])]
		}

for key, value in transcription_table.items():
	transcription_table[key].update(ambiguous_nucleotides(key))

base_pairs = {
	'A': 'T',
	'T': 'A',
	'G': 'C',
	'C': 'G'
	}

def reverse_transcribe(seq):
	return ''.join(base_pairs[i] for i in reversed(seq))

def return_frameshift_sequence(full_seq, frameshift_string):
	# Subtract 1 from start position to account for 0 indexing
	seq = ''
	for x in frameshift_string.split(','):
		left_pos, right_pos = x.split(':')
		seq += full_seq[int(left_pos)-1:int(right_pos)]
	return seq

def extract_sequence(full_seq, left_pos, right_pos, strand):
	#seq = full_seq[left_pos:right_pos] # old code
	if len(left_pos.split(',')) == len(right_pos.split(',')):
		left_pos = left_pos.split(',')
		right_pos = right_pos.split(',')
		seq = '' # it can be Bio.Seq.Seq()
		for idx in range(len(left_pos)):
			seq += full_seq[int(left_pos[idx]):int(right_pos[idx])]
	else:
		raise ValueError('Gene coordinates are ill-defined.')

	if strand == '+':
		return seq
	elif strand == '-':
		return reverse_transcribe(seq)
	else:
		raise ValueError('strand must be either \'+\' or \'-\'.')

def get_amino_acid_sequence_from_dna(dna_seq):
	if len(dna_seq) % 3 != 0:
		raise ValueError('Gene nucleotide sequence is not a valid length.')

	codons = (dna_seq[i: i + 3] for i in range(0, (len(dna_seq)), 3))
	amino_acid_sequence = ''.join(codon_table[i] for i in codons)
	amino_acid_sequence = amino_acid_sequence.rstrip('*')
	if '*' in amino_acid_sequence:
		# one letter code for selenocysteine, not C(ysteine) (and wrong because precursor is serine)
		amino_acid_sequence = amino_acid_sequence.replace('*', 'U')
	return amino_acid_sequence
