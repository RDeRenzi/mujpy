.. _installation:

Installation
============

Mujpy is `python3` native. 

`Linux`_ `Windows`_ `Macos`_

Linux
-----

Python comes with all distributions, but ``pip`` does not. Check typing ``pip --version``.  If you get `pip: command not found` run ``python3 -m pip install --upgrade pip``. 
It is best to create a virtual environment: run ``python3 -m venv ~/.mujpy-venv`` (see the `official notes <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>`_ for more details) and lauch it by ``source ~/.mujpy-venv/bin/activate`` (``deactivate`` will exit venv). 
The terminal prompt has now ``(.mujpy-venv)`` prepended. Install ``mujpy`` by running

.. code::

    python3 -m pip install mujpy

This provides all the required dependencies and, from now on, each time you activate the venv any python or jupyter command knows mujpy. 

In the future, to upgrade to the newest distribution, ``python3 -m install mujpy --upgrade``. If you are impatient you can also

.. code::

        git clone https://github.com/RDeRenzi/mujpy/ 
 
Windows
-------

Install Ubuntu under WSL2 (`wsl --install -d Ubuntu` in cmd Admin, it requires ~10 GB of disk space).
After restarting Windows, start an Ubuntu app (terminal), cut&paste the following:

.. code::

      sudo apt update
      sudo apt install python3-pip python3-venv python-is-python3 \
      libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0 \
      libxcb-shape0 libxcb-xinerama0 libxcb-xkb1 libxkbcommon-x11-0  

wait for completion and then follow the linux instructions to install `mujpy` under a venv. 

Macos
-----

Should work as in linux, but never checked


To check your mujpy installation see :doc:`Tutorial`

.. rubric:: Footnotes

.. [#1] after starting the venv, if you opt for this
