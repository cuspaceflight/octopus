from distutils.core import setup

setup(
    name='octopus',
    version='0.1',
    packages=['octopus'],
    install_requires=['numpy', 'scipy', 'thermo'],
    data_files=[('lib/site-packages/octopus/data', ['octopus/10024-97-2.json'])],
    license='GNU General Public License v3.0',
    author_email='ec765@cam.ac.uk',
    description='Utility for analysis of 2-phase compressible flow through an injector'
)
