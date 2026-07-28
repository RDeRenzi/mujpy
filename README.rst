*****
MuJPy
*****

A Python MuSR data analysis based on classes, with a graphical interface designed for jupyter, released under the MIT licence.

v. 2.9.beta refactoring. See _version.py for the most recent version details. 
For the time being works only via python3 scripts (test.py ties them all)
(mudash, next step, is still broken). i

Once cloned, run test.pyr,
or run the individual tests, 
e.g. for a single run single group fit
:: 
 
  python mgml.822.3-4.1_fit.py

Modify these scripts to your purposes. The clumsy part is to modify by hand the models fit/... .json files.
mudash will soon come to make this user-friendly
Linux installation instructionsi. In principle valid under windows with possible os related problems. 
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
