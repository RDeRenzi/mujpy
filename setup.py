#!/usr/bin/env python
from setuptools import setup
    
setup(name='mujpy',
      version='2.2.7',
      description='A Python MuSR data analysis program designed for Jupyterlab. ',
      author='Roberto De Renzi, Pietro Bonfa',
      author_email='roberto.derenzi@unipr.it',
      url='https://github.com/RDeRenzi/mujpy',
      packages=['mujpy',
                'mujpy.tools',
                'mujpy.logo',
                'mujpy.mucomponents',
                ],
      include_package_data=True,
      package_dir={'mujpy': 'mujpy' },
      install_requires=[
            'numpy',
            'scipy',
            'jupyterlab'
            'ipywidgets ',
            'iminuit >= 2.17',
            'matplotlib >= 2.0',
            'NeXus>=4.4.1',
            'nexusformat>= 0.5.3',
            'python-nexus>= 2.0.1',
            'musr2py>=0.0.2',
            'munexus>=0.1'
      ],
      long_description='A Python MuSR data analysis program designed for Jupyterlab. User-friendly, does single, sequential and global fits, nice animation display for multi-run fits. Works on Linux, Windows, MacOS',
      license = 'MIT',
      keywords = 'musr fit',
      classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Physics',
            'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3',
]
)
