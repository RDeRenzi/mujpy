.. _installation:

Installation
============

Mujpy is `python3` native. 

`Linux`_ `Windows`_ `Macos`_

Linux
-----

Python comes with all distributions, but ``pip`` does not. Check typing ``pip --version``.  If you get `pip: command not found` run ``python3 -m pip install --upgrade pip``. 
It is best to create a virtual environment: run ``python3 -m venv ~/.mujpy-venv`` (see the `official notes <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>`_ for more details) and lauch it by ``source ~/.mujpy-venv/bin/activate`` (``deactivate`` will exit venv). 
The terminal prompt has now ``(.mujpy-venv)`` prepended. Install ``mujpy`` once 

.. code::

    python3 -m pip install mujpy

This provides all the required dependencies and, from now on, each time you activate the venv any python or jupyter command knows mujpy. 

To upgrade to the newest distribution, ``python3 -m install mujpy --upgrade``. If you are impatient you can also

.. code::

        git clone https://github.com/RDeRenzi/mujpy/ 
 
Windows
-------

As of 09/26, `mujpy` is still broken on Windows. Just in case from `link <https://www.python.org/downloads/windows>`_ click on the link ``Latest Python install manager`` and then on the link ``using the Microsoft Store app``. Install it. You now have also ``pip``, check by opening ``Start PowerShell`` (not ``Start cmd``!!) and type  ``py -m pip --version``. 

Always using PowerShell, create a virtual environment with e.g. ``py -m venv .pip-mujpy``, and in the same folder run ``.pip-mujpy/Script/activate``. Your PowerShell terminal now has `(.git-mujpy)` prepended to the prompt. It guarantees that you are inside the venv.  Now run ``pip install mujpy``. 


Macos
-----

Should work as in linux, but never checked


To check your mujpy installation see :doc:`Tutorial`

.. rubric:: Footnotes

.. [#1] after starting the venv, if you opt for this
