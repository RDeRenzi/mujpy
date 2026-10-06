.. encoding: utf-8

=====
mujpy
=====

A Python μSR data analysis package, for both command line and a graphical interface
running under Jupyter and Voilà, released under the GPL-3 licence. 

It aims at the power of musrfit with the user-friendly appearance of mulab.

From version 3.1.3 on, the main technical features are: 

- a model built on two-letter bricks: mg, for
  Gaussian-damped cosine, ml for Lorentzian-damped cosine etc. 
- sequential fits by the same model, driven by a run list and a list of
  grouping dictionaries, for asymmetry definition 
- global fits by user-defined parameters assigned to model parameters in a json file 
- a mulab-like gui interface in jupyterlab that allows fit model and parameter editing,
  hopefully with a gentler learning curve than musrfit
- direct standalone gui web interface by ```voila```, included

**Try mujpy!** 

on linux, or on Windows by wsl (requires ~10GB for ubuntu-in-Windows). Follow simple `installation instructions <https://mujpy.readthedocs.io/latest/Installation.html>`_. Both os come with python (wsl requires a further step), you only need to:  

* create a venv ``python -m venv ~/.mujpy-env``, 
* activate it ``source ~/.mujpy-env/bin/activate``, 
* invoke ``python -m pip install --upgrade mujpy`` 

Now try your freshly installed ``mujpy`` from command line, to demonstrate its capabilities. Type ``python -m mujpy.tests.tests``: 16 fully automated command line executions will start and popup their graphics (find the description in the `Tutorial <https://mujpy.readthedocs.io/latest/Tutorial.html>`_).

For normal use a GUI (with its onw demos) is available (see basic instructions on `ReaTheDocs <https://mujpy.readthedocs.io/latest/Tutorial.html>`_) with choice of popuop plots - jupyter notebook ``%matplotlib qt`` - or Log-tab side-by-side log+plot - ``%matplotlib widget``. Dircet browser access with ``voila``

Please email bugs to the author.



