.. _installation:

Installation
============

Mujpy is `python3` native. 

`Linux`_ `Windows`_ `Macos`_

Linux
-----

Python comes with most distributions. To install jupyterlab

::

    python3 -m pip install --upgrade pip
    python3 -m pip install jupyterlab

::

    python3 -m pip install mujpy

To upgrade to the newest distribution
 git clone https://github.com/RDeRenzi/mujpy/ 
 
Windows
-------

download Anaconda from https://anaconda.com_
execute Anaconda.exe 
launch Anaconda Prompt (Anaconda 3)
from the Anaconda Prompt run

::

    conda install -c conda-force jupyterlab
    
download https://bootstrap.pypa.io/get-pip.py_
find the Download path (download_path), and run the following in a jupyterlab cell 
launch jupyterlab

::

    %cd download_path
    %run get-pip.py
    pip install mujpy


Macos
-----

to come
