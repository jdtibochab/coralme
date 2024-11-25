# -*- coding: utf-8 -*-

# Always prefer setuptools over distutils
#from ez_setup import use_setuptools
#use_setuptools()
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path, walk
import versioneer

def main():
	# additional files
	data_files = []
	for dirpath, dirnames, filenames in walk('coralme/iJL1678b-ME'):
		tmp = []
		for filename in filenames:
			tmp.append(path.join(dirpath, filename))
		data_files.append(('coralme', tmp))

	# Get the long description from the README file
	here = path.abspath(path.dirname(__file__))
	with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
		long_description = f.read()

	setup(
		name='coralME',
		license='MIT',
		#version='1.0', # use tags

		version=versioneer.get_version(), # from versioneer
		cmdclass=versioneer.get_cmdclass(), # from versioneer

		description='Comprehensive Reconstruction Algorithm for ME-models (coralME)',
		long_description=long_description,
		long_description_content_type='text/markdown',
		url='https://github.com/jdtibochab/coralme',
		author='Juan D. Tibocha-Bonilla and Rodrigo Santibanez-Palominos',
		author_email='jdtibochab@users.noreply.github.com',
		classifiers=[
			#'Development Status :: 1 - Planning',
			#'Development Status :: 2 - Pre-Alpha',
			#'Development Status :: 3 - Alpha',
			#'Development Status :: 4 - Beta',
			'Development Status :: 5 - Production/Stable',
			#'Development Status :: 6 - Mature',
			#'Development Status :: 7 - Inactive',

			# Indicate who your project is intended for
			'Intended Audience :: Science/Research',
			'Topic :: Scientific/Engineering :: Bio-Informatics',

			# Pick your license as you wish
			'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

			# Specify the Python versions you support here. In particular, ensure
			# that you indicate whether you support Python 2, Python 3 or both.
			#'Programming Language :: Python :: 2',
			#'Programming Language :: Python :: 2.7',
			'Programming Language :: Python :: 3',
			#'Programming Language :: Python :: 3.4',
			#'Programming Language :: Python :: 3.5',
			#'Programming Language :: Python :: 3.6',
			'Programming Language :: Python :: 3.7',
			'Programming Language :: Python :: 3.8',
			'Programming Language :: Python :: 3.9',
			'Programming Language :: Python :: 3.10',

			# other classifiers added by author
			'Environment :: Console',
			'Operating System :: Unix',
			'Operating System :: Microsoft :: Windows',
			'Operating System :: MacOS',
		],

		python_requires='~=3.0',
		keywords=[
			'metabolism',
			'biology',
			'constraint-based',
			'linear programming',
			'mixed-integer',
			'optimization',
			'flux-balance analysis',
			'reconstruction'
			],
		install_requires=[
			'cobra<=0.29.0',
			'python-libsbml<=5.20.1',
			'Biopython==1.80',
			'anyconfig<=0.13.0',
			'pyranges<=0.0.129',
			'XlsxWriter<=3.1.2',
			'openpyxl<=3.1.2',
			'numpy<=1.26.4',
			'scipy<=1.11.1',
			'sympy<=1.12',
			'pandas<=1.5.1',
			'jsonschema<=4.21.1',
			'tqdm<=4.62.3',
			],

		# WARNING: seems to be bdist_wheel only
		packages=find_packages(exclude=('contrib', 'docs', 'tests', 'tutorials')),
		# using the MANIFEST.in file to exclude same folders from sdist
		include_package_data=False,

		# WARNING: use this way to install in the package folder and
		# to have access to files using importlib_resources or importlib.resources
		# bdist_wheel only
		package_data = {
			'coralme' : [
				'iJL1678b-ME/*',
				'iJL1678b-ME/building_data/*',
				'solver/*.so',
				'io/JSONSCHEMA',
				'builder/column_format.json',
				'licenses/*'
				]
			},

		# WARNING: do not use data_files as the installation path is hard to determine
		# e.g.: ubuntu 18.04 installs to /usr/local/installation_path or /$USER/.local/installation_path
		#data_files=[('installation_path', ['filename'])],
		#data_files=data_files,

		# others
		project_urls={
			'Manual': 'https://coralme.readthedocs.io',
			'Bug Reports': 'https://github.com/jdtibochab/coralme/issues',
			'Source': 'https://github.com/jdtibochab/coralme',
		},
	)

if __name__ == '__main__':
    main()
