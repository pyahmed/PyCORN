#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    print("setuptools not found, falling back to distutils")
    from distutils.core import setup

setup(
    name='pycorn',
    version='0.14',
    author='Yasar L. Ahmed',
    packages=['pycorn',],
    extras_require = {'plotting':  ["matplotlib"], 'xlsx-output': ['xlsxwriter']},
    scripts=['examplescripts/pycorn-bin.py'],
    license='GNU General Public License v2 (GPLv2)',
    description='A script to extract data from UNICORN result (res) files',
    long_description=open('README.rst').read(),
    url='https://github.com/pyahmed/PyCORN',
)
