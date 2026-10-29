.. _installation:

Installation
============

Mujpy is `python3` native! 

`Linux`_ `Windows`_ `Macos`_

Linux
-----

Python comes with all distributions. [1]_
It is best to create a virtual environment: run ``python3 -m venv ~/.mujpy-venv`` and lauch it by ``source ~/.mujpy-venv/bin/activate`` (``deactivate`` will exit venv). 
The terminal prompt has now ``(.mujpy-venv)`` prepended. Install ``mujpy`` by running

.. code::

    python3 -m pip install --upgrade mujpy

This provides all the required dependencies for the jupyter widget/ipympl matplotlib backend, and, from now on, each time you activate the venv any python or jupyter command knows mujpy. If you have Qt5 and want to enable also popup Qt graphics, run instead 

.. code::

    python -m pip install "mujpy[qt]"  # also enable Qt popup

For those who want to have the source, the repository can by obtained creating a git folder and running from there 

.. code::

        git clone https://github.com/RDeRenzi/mujpy/ 
 
Windows
-------

Install Ubuntu under WSL2 (it requires ~10 GB of disk space, from cmd Admin run `wsl --install -d Ubuntu`).
Restart Windows, start an Ubuntu app (terminal), and cut&paste [2]_ the following:

.. code::

      sudo apt update
      sudo apt install python3-pip python3-venv python-is-python3 \
      libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0 \
      libxcb-shape0 libxcb-xinerama0 libxcb-xkb1 libxkbcommon-x11-0  

Wait for completion and follow the linux instructions to install `mujpy` under a venv. 

Macos
-----

Should work as in linux (never checked, though)


To check your mujpy installation see :doc:`Tutorial`

From git
--------

mujpy is not released very often to PyPI.  If you need the latest developments use git.
Both Linux and WSL2 on Windows have it natively. From `$HOME`

.. code::

   mkdir git
   cd git
   git clone https://github.com/RDeRenzi/mujpy.git

then activate your mujpy venv, as explained in `Linux`_, make sure you are in `git/mujpy` and run

.. code::

   pyton -m pip install --upgrade .

and continue with the Tutorial instructions

.. rubric:: Footnotes

.. [1] but ``pip`` does not. Check typing ``pip --version``.  If you get ``pip: command not found`` run ``python3 -m pip install --upgrade pip``. Just use ``python`` if you installed ``python-is-python3`` (in Ubuntu).
.. [2] hover on the green box and click the copy widget.
