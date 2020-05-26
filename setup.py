from distutils.core import setup
import sys

script_args = sys.argv[1:]

setup(name='assess',
      version='0.2',
      author='Stephan Kramer',
      author_email='s.kramer@imperial.ac.uk',
      description='Assess is a python package that implements a number of '
      'analytical solutions to the Stokes equations in cylindrical and '
      'spherical shell domains.',
      long_description=open('README.rst').read(),
      url='https://github.com/stephankramer/assess',
      packages=['assess'],
      install_requires=['numpy', 'scipy'],
      keywords=['Stokes equations', 'analytical solutions',
                'mantle flow', 'spherical shell'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering'],
      script_args=script_args,
      ext_package='assess',
      )
