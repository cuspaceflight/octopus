from distutils.core import setup

setup(
    name='octopus',
    version='0.1',
    packages=['octopus'],
    install_requires=['numpy', 'scipy', 'matplotlib', 'CoolProp', 'h5py'],
    package_data={'octopus': ['data/*.json']},
    license='GNU General Public License v3.0',
    author_email='ec765@cam.ac.uk',
    description='Utility for analysis of 2-phase compressible flow through an injector'
)
