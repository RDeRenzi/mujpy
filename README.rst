.. encoding: utf-8
.. |mSR| replace:: :math:`\mu`\ SR

=====
mujpy
=====

A Python |mSR| data analysis based on classes, with a graphical interface
designed for jupyter, released under the GPL-3 licence. 

It aims at the power of musrfit
with the user-friendly appearance of mulab.

Version 3.0 refactoring. See Changelog for the most recent version
details. Main technical features: 

- a model built on two-letter bricks: mg, for
  Gaussian-damped cosine, ml for Lorentzian-damped cosine etc. 
- sequential fits by the same model, driven by a run list and a list of
  grouping dictionaries, for asymmetry definition 
- global fits by user-defined parameters assigned to model parameters in a json file 
- a mulab-like interface in jupyterlab that allow fit model and parameter editing,
  hopefully with a gentler learning curve than musrfit
- web interface by ```voila```, included

**Try mujpy!** 

(for the time being linux-only, Windows coming soon). Any linux has python (or has it?), you only need to:  [#1]_
create a venv ``python -m venv ~/.mujpy-env``, activate it ``source ~/.mujpy-env/bin/activate``, and invoke ``pip install mujpy``
Now try your freshly installed ``mujpy`` from command line to demonstrate its capabilities  (16 fully automatic casesin order of growing complexity). Just type ``python -m mujpy.tests.tests``. Kill each graphic window after inspection, see ReadTheDocs `Tutorial https://mujpy.readthedocs.io/latest/Tutorial.html<>`_ for more info.

This is the command line demos. There is also a GUI and its demos, coming instructions on ReaTheDocs.


.. rubric:: Footnotes

.. [#1]  in principle ``mujpy`` works on all OS, but for the moment it is tested only on linux
