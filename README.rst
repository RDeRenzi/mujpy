*****
MuJPy
*****

A Python MuSR data analysis based on classes, with a graphical interface designed for jupyter, released under the MIT licence.

v. 2.0.alpha major refactoring. See _version.py for more recent version details. For the time being works only via python3 scripts (see example test.py) and jupyter-lab notebooks 
(mudash is broken). Once cloned, run test.py, or::

  jupyter-lab

and launch the example/Tst_xxx.ipynb notebooks (check path to suitable data files for the tested type of fit)

or run, e.g. for a single run single group fit:: 
 
  python mgml.822.3-4.1_fit.py

Each notebook shows how to perform
* TF calibration, both single group (A1_calib) and multigroup, either sequential (A20_calib) or global (A21_calib) 
* correspondinmg single run fit (A1,A20,A21), 
* sequential single run fits (B1),
with plots: static for single fit, animated plots for multiple fits (both sequential or global) and optional fft

Modify these scripts to your purposes. The clumsy part is to modify by hand the models fit/... .json
Linux installation instructions, Valid on WIN10 that has a linux shell, to come.
Docs not updated yet

----
Old installation instructions (v 1.1)
* Make sure you have python, standard on linux, and jupyter lab. Otherwise install them (see https://docs.python.org/3/using/windows.html, https://docs.python.org/3/using/mac.html, jupyter.readthedoc.io).
* Install mujpy. Clone or download from https://github.com/RDeRenzi/mujpy, unzip into the directory of your choice::

   cd mujpy/mujpy/musr2py
   make
   sudo make install

* Check dependencies, see requirements.txt. When the first mujpy release is on Pipy, pip will sort this out.

* Start jupyter lab::

   jupyter-lab

* Copy example/Tst_xxx_ipynb to a path of your choice and modify path /home/roberto.derenzi/mujpy/log/ accordingly

Documentatation on the GUI usage at http://mujpy.readthedocs.io/en/latest/ obsolete. DOcs to come
