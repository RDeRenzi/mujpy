mujpy
=====

A Python MuSR data analysis based on classes, with a graphical interface
designed for jupyter, released under the MIT licence. It aims at the
useri-friendly appearance of mulab and the power of musrfit.

Version 3.0 refactoring. See Changelog for the most recent version
details.

Main features: - a model built on two-letter bricks: mg, for
Gaussian-damped cosine, ml for Lorentzian-damped cosine etc. -
sequential fits by the same model driven by a run list and a list of
grouping dicts for asymmetry definition - global fits by user-defined
parameters assigned to model parameters in a json file - a mulab-like
interface in Jupyter-lab that allow fit model and parameter editing,
hopefully with a gentler learning curve than musrfit

To try mujpy: 

  1. Install python (python3 is assumed) 
  2. Install git and clone `mujpy <https://github.com/RDeRenzi/mujpy.git>`_

The scripts collected in test.py show how to drive mujpy in command line
mode. To see some fit capabilities try:

::

   cd example
   python test.py 

Modify these scripts for your purposes. The clumsy part is to modify by
hand the models fit/… .json files. To avoid that

   3.a Install jupiterlab to run the gui editor notebook Mudashed.ipynb

or 

   3.b ```pip nstall voila``` and run ```voila Mudashed.ipynb``` (with `snap-installed firefox bug workaround? <https://github.com/voila-dashboards/voila/issues/1508>`_)

Docs and installation instructions still are work-in-progress
(`ReadTheDocs <http://mujpy.readthedocs.io/en/latest/>`_ obsolete)

--------------

Old installation instructions (v 1.1)
-------------------------------------

-  Make sure you have python, standard on linux, and jupyter lab.
   Otherwise install them (see
   https://docs.python.org/3/using/windows.html,
   https://docs.python.org/3/using/mac.html, jupyter.readthedoc.io).

-  Install mujpy. Clone or download from
   https://github.com/RDeRenzi/mujpy, unzip into the directory of your
   choice:

   ::

      cd mujpy/mujpy/musr2py
      make
      sudo make install

-  Start jupyter lab:

   ::

      jupyter-lab

