from setuptools import setup, find_packages

setup(
    name='octopus',
    version='0.1-beta',
    packages=find_packages(),
    install_requires=['numpy','scipy','thermo'],
    include_package_data=True,
    license='GNU General Public License v3.0',
    author_email='ec765@cam.ac.uk',
    description='Utility for analysis of 2-phase compressible flow through an injector'
)
