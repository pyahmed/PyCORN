try:
    from setuptools import setup
except ImportError:
    print("setuptools not found, falling back to distutils")
    from distutils.core import setup

setup(
    name='pycorn',
    version='0.1dev0',
    packages=['pycorn',],
    license='GNU General Public License v2 (GPLv2)',
    long_description=open('README.md').read(),
)
