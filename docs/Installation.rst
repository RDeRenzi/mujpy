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

As of 09/26, from `i <https://www.python.org/downloads/windows>`_ click on the link``Latest Python install manager`` and then on the link ``using the Microsoft Store app``. Install it. You now have ``pip``, check by opening ``Start PowerShell`` (not ``Start cmd``!!) and type  ``py -m pip --version``. 

Always using PowerShell, create a virtual environment with ``py -m venv .pip-mujpy``, and in the same folder run ``.pip-mujpy/Script/activate``. Your PowerShell terminal now has (.git-mujpy) prepended to the prompt. It means that you are in the venv.  Now run ``pip install mujpy``. 

To check your mujpy installation ``cp 
Or ``cp .git-mujpy\lib\site-packages\mujpy\test\Mudashed-demos.ipynb .`` and launch ``voila Mudashed-d
Macos
-----

Should work as in linux


.. rubric:: Footnotes

.. [#1] after starting the venv, if you opt for this
