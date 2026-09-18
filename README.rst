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

Now try your freshly installed ``mujpy`` from command line to demonstrate its capabilities  (fully automatic). Just type ``python -m unittest mujpy.tests/test``, e.g. from ``/tmp``.


.. code-block::

     %matplotlib qt  

and

.. code-block::

        from mujpy.mudashed import dashed as mudash
        the_dash = mudash()

Otherwise, download ```Mudashed.ipynb``` and ```Mudashed-demos.ipynb``` from the ```tests/``` folder of this repository and 
test the gui: 

   a. create project directory and cd to it
   b. copy ``Mudashed-demos.ipynb`` into it
   d. type ``voila Mudashed-demos.ipynb`` in a terminal (``voila`` comes with ``mujpy``) [#2]_
   e. Group0 already reads ``3-4`` (Up-Down in GPS), leave it
   f. press DL and choose the first datafile (e.g. a transverse field run, ``822``)
   g. insert the run number in run list and hit ``Enter`` (or press RL)
   h. press LF and double-click on ```almgml.822.3-4.1_fit.json``` in teh pop-up
   i. press Fit, check the result
   j. change run list to ``822,834`` or ``822,827:834:-1`` and press LF again, check the result

For more complex fits see the Introduction 
in `ReadTheDocs <http://mujpy.readthedocs.io/en/latest/>`_
For a first demo try out [#3]_ the script ``test.py``: choose an empty folder, edit a file ``savetests.py`` and copy-paste the following code:

.. code-block::

   from mujpy.tools.tools import savetests
   savetests()

Execute ``python savetests.py``. This allows you to run ``python test.py``, to get a look and feel of the 16 different types of fit that mujpy provides. 

Modify the script for your purposes. The clumsy part here is to modify by
hand the models in the ``fit/… .json`` files. Use the gui instead.


.. rubric:: Footnotes

.. [#1]  in principle ``mujpy`` works on all OS, but for the moment it is tested only on linux
.. [#2]  snap-installed firefox has a known bug with ``voila``, follow ` <https://github.com/voila-dashboards/voila/issues/1508>`_ for a simple workaround
.. [#3]  if ``pip install mujpy`` was successful
