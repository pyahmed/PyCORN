#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    print("setuptools not found, falling back to distutils")
    from distutils.core import setup

setup(
    name='pycorn',
    version='0.18',
    author='Yasar L. Ahmed',
    packages=['pycorn'],
    extras_require = {'plotting':  ["matplotlib"], 'xlsx-output': ['xlsxwriter']},
    scripts=['examplescripts/pycorn-bin.py'],
    platforms=['Linux', 'Windows', 'MacOSX'],
    zip_safe=False,
    classifiers=["License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
                 "Environment :: Console",
                 "Intended Audience :: Science/Research",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 2.7",
                 "Programming Language :: Python :: 3.4",],
    package_data={'pycorn': ['docs/*.*']},
    license='GNU General Public License v2 (GPLv2)',
    description='A script to extract data from UNICORN result (res) files',
    long_description=open('README.rst').read(),
    url='https://github.com/pyahmed/PyCORN',
)
