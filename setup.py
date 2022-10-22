#!/usr/bin/env python
from setuptools import setup
    
setup(name='mujpy',
      version='2.0.1',
      description='A Python MuSR data analysis program designed for jupyter.',
      author='Roberto De Renzi, Pietro Bonfa',
      author_email='roberto.derenzi@unipr.it',
      url='https://github.com/RDeRenzi/mujpy',
      packages=['mujpy',
                'mujpy.aux',
                'mujpy.logo',
                'mujpy.mucomponents',
                ],
      include_package_data=True,
      package_dir={'mujpy': 'mujpy' },
      install_requires=[
            'numpy >= 1.6',
            'ipywidgets >= 7.5.1',
            'iminuit >= 1.2',
            'matplotlib >= 2.0'
            'musr2py >= 0.0.1'
      ],
      long_description='A Python MuSR data analysis program designed for jupyter.'
      license = 'MIT',
      classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Users',
            'Topic :: Physics',
            'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3',
],
      install_requires=['iminuit,ipywidgets,jupiter,numpy,scipy,matplotlib'],
)
