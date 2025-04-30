import difflib
import glob
import pathlib

def main(path1, path2):
	files = [ '/'.join(x.split('/')[1:]) for x in sorted(glob.glob('{:s}/**/*'.format(path1), recursive=True)) ]
	files = [ x for x in files if not x.endswith(('.log', '.pkl', '.json', '.yaml', '.xlsx')) ]
	files = [ x for x in files if not x.startswith(('blast_files_and_results')) ]
	# files

	def sort(data):
		if len(data) > 1:
			data = [data[0]] + sorted(data[1:])
		return data

	def compare_files(file1_path, file2_path):
		"""Compares two files and prints the differences using difflib."""
		with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
			file1_lines = file1.readlines()
			file2_lines = file2.readlines()

		file1_lines = sort(file1_lines)
		file2_lines = sort(file2_lines)

		diff = difflib.unified_diff(file1_lines, file2_lines, fromfile=file1_path, tofile=file2_path, n=0)

		for line in diff:
			print(line, end='')  # Print the diff output
		return file1_lines, file2_lines

	for file in files:
		if file.endswith(('m_model_modified.json', 'genome_modified.gb', 'curation_notes.txt', 'config.txt')):
			continue
		if pathlib.Path(f'{path1}/{file}').is_dir():
			continue

		f1, f2 = compare_files(f'{path1}/{file}', f'{path2}/{file}')
		if f1 != f2:
			print()

if __name__ == '__main__':
	main(path1, path2)
