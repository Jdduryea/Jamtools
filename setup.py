# BETA
import os

from setuptools import setup, find_packages

# Convert README.md to README.rst automagically
# http://stackoverflow.com/a/23265673
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()

here = os.path.abspath(os.path.dirname(__file__))

README = open(os.path.join(here, 'README.md')).read()
CHANGES = '' # open(os.path.join(here, 'CHANGES.md')).read()

setup(
  name='jamtools',
  version='1.0',
  description='methtuple',
  long_description=read_md('README.md'),
  license='MIT',
  classifiers=[
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.2",
    "Programming Language :: Python :: 3.4",
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    "License :: OSI Approved :: MIT License",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
  author='Jack Duryea',
  author_email='duryea@bcm.edu',
  url='https://github.com/Jdduryea/Jamtools/',
  keywords='bisulfite sequencing methylation bismark bioinformatics samtools',
  packages=find_packages(),
  include_package_data=True,
  zip_safe=False,
  install_requires=['pysam >= 0.8.4', 'matplotlib >= 2.0.0'],
  test_suite="methtuple.tests",
  scripts = [
    'methtuple/scripts/methtuple'
  ]
)
