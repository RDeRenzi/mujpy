.. _installation:

Installation
============

Mujpy is `python3` native. 

`Linux`_ `Windows`_ `Macos`_

Linux
-----

Python comes with all distributions, but ``pip`` does not. Run ``python3 -m pip install --upgrade pip`` and then [#1]_

.. code::

    python3 -m pip install mujpy

This provides all the required dependencies. 

If you want a lean mujpy python environment it is strongly suggested that you create a virtual environment: ``python3 -m venv ~/.mujpy-venv`` (see the `official notes <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>`_ for more details) and lauch it by ``source ~/.mujpy-venv/bin/activate`` (``deactivate`` will exit venv). 

The terminal prompt has now ``(.mujpy-venv)`` prepended. Install ``mujpy`` once and, from now on any python or jupyter command sees mujpy each time you activate the venv.

To upgrade to the newest distribution
 git clone https://github.com/RDeRenzi/mujpy/ 
 
Windows
-------

Download the python `installer <https://www.python.org/downloads/>`_ and run it. See if you have already pip by ``py -m pip --version``, if you do not ``py -m ensurepip --default-pip``. Create a vitrual environment: ``py -m venv mujpy-env``, ``mujpy-env/Script/activate``, choose where you want to install mujpy  and ``pip install mujpy``. If it still does not work for windows you can install git.

Download `git <https://git-scm.com/install/windowsi>`, run it, accepting all defaults and including Git Bash.  From the Git Bash term ``git clone https://github.com/RDeRenzi/mujpy.git``, ``cd mujpy`` and ``pip install .``


Macos
-----

Should work as in linux


.. rubric:: Footnotes

.. [#1] after starting the venv, if you opt for this
