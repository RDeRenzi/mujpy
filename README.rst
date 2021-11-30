*****
MuJPy
*****

A Python MuSR data analysis graphical interface, based on classes, designed for jupyter.

Released under the MIT licence.

v. 2.0.alpha major refactoring

for the time being works only via jupyter-lab notebooks
e.g. once downloaded, make suure you have afs access to PSI data and launch the Tst-xxx.ipynb notebook

Performs 
* TF calibration (A1-calib), 
* single run fit (A1), 
* sequential single run fits (B1),
* single run multi-group global fit (A2)
with plots, static for single fit, animated plots for multiple fits (both sequential or global) 

mugui is broken 
  
Linux installation instructions, Valid on WIN10 that has a linux shell, to come.
Docs not updated yet

----
Old installation instructions for v 1.0
* Make sure you have python, standard on linux, and jupyter lab. Otherwise install them (see https://docs.python.org/3/using/windows.html, https://docs.python.org/3/using/mac.html, jupyter.readthedoc.io).
* Install mujpy. Download from https://github.com/RDeRenzi/mujpy, unzip into the directory of your choice::

   cd mujpy/mujpy/musr2py
   make
   sudo make install
   cd ../..
   # preferably use python 3
   python[3] setup.py install

* Check dependencies, see requirements.txt. When the first mujpy release is on Pipy, pip will sort this out.

* Start jupyter lab::

   jupyter-lab

* Use mu_test_damg_lab.ipnb. Alternatively the first cell type::

  >>>%matplotlib # "%matplotlib tk"  avoids a line of output if you are using the tk backend

* Add a cell and write::

   from mujpy.mugui import mugui
   MuJPy = mugui()

* Run all the cells (if ipywidget does not work at notebook start select Kernel/Restart & Run all)

* install jq https://stedolan.github.io/jq/download/
Documentatation on the GUI usage at http://mujpy.readthedocs.io/en/latest/
