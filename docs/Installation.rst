.. _installation:

Installation
============

Mujpy is `python3` native. 

`Linux`_ `Windows`_ `Macos`_

Linux
-----

Python comes with all distributions, but ``pip`` does not. Check typing ``pip --version``.  If you get `pip: command not found` run ``python3 -m pip install --upgrade pip``. 
It is best to create a virtual environment: run ``python3 -m venv ~/.mujpy-venv`` (see the `official notes <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>`_ for more details) and lauch it by ``source ~/.mujpy-venv/bin/activate`` (``deactivate`` will exit venv). 
The terminal prompt has now ``(.mujpy-venv)`` prepended. Install ``mujpy`` once. Normally, after ``activate``, ``cd`` to a suitable directory and use command line ``python -m mujpy.tests.tests`` to see a full demo, or ``jupyter lab`` to start the GUI interface by luaunching a new python3 notebook with either ``%matplotlib qt`` (popup plots) or ``%matplotlib widget`` (notebook plots) in the first cell, and ``from mujpy.mudashed import dashed``, ``the_dash=dasehd(test='GPS')`` in the second cell. Run it all.

.. code::

    python3 -m pip install mujpy

This provides all the required dependencies and, from now on, each time you activate the venv any python or jupyter command knows mujpy. 

To upgrade to the newest distribution, ``python3 -m install mujpy --upgrade``. If you are impatient you can also

.. code::

        git clone https://github.com/RDeRenzi/mujpy/ 
 
Windows
-------

As of 09/26, `ipywidgets` breaks on Windows. If you can afford 10 GB to install Ubuntu (`wsl --install -d Ubuntu` in cmd) you can then run it in an Ubuntu WSL2 window, install `mujpy` as in linux. Command line works and notebook works with ``%matplotlib widget`` first cell (plot to the right of log) (does it?). After `sudo apt update && sudo apt upgrade -y && sudo apt install python3.14-venv`, follow linux instructions, replacing ``python3`` for ``python`` and chooseing ``widget`` plots.

Macos
-----

Should work as in linux, but never checked


To check your mujpy installation see :doc:`Tutorial`

.. rubric:: Footnotes

.. [#1] after starting the venv, if you opt for this
