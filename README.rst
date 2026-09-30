.. encoding: utf-8
.. |mSR| replace:: :math:`\mu`\ SR

=====
mujpy
=====

A Python |mSR| data analysis based on classes, with a graphical interface
designed for jupyter, released under the GPL-3 licence. 

It aims at the power of musrfit
with the user-friendly appearance of mulab.

Version 3.1. Main technical features: 

- a model built on two-letter bricks: mg, for
  Gaussian-damped cosine, ml for Lorentzian-damped cosine etc. 
- sequential fits by the same model, driven by a run list and a list of
  grouping dictionaries, for asymmetry definition 
- global fits by user-defined parameters assigned to model parameters in a json file 
- a mulab-like gui interface in jupyterlab that allows fit model and parameter editing,
  hopefully with a gentler learning curve than musrfit
- direct standalone gui web interface by ```voila```, included

**Try mujpy!** 

on linux, or on Windows by lightweight WSL2 (~10GB of ubuntu-in-Windows). Follow simple `installation instructions <https://mujpy.readthedocs.io/latest/Installation.html>`_. Both os come with python, you only need to:  
create a venv ``python -m venv ~/.mujpy-env``, activate it ``source ~/.mujpy-env/bin/activate``, and invoke ``pip install mujpy``.
Now try your freshly installed ``mujpy`` from command line, to demonstrate its capabilities. Just type ``python -m mujpy.tests.tests`` and 16 fully automatic cases will popup, in order of growing complexity (find a description in the `Tutorial <https://mujpy.readthedocs.io/latest/Tutorial.html>`_).

There is also a GUI and its demos, see basic instructions on `ReaTheDocs <https://mujpy.readthedocs.io/latest/Tutorial.html>`_. Choice of jupyter notebook ``%matplotlib qt`` for popup plots, or ``%matplotlib widget`` for Log-tab side-by-side log+plot. 

Please email bugs to the author.



