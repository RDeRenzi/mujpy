MuJPy
=====

A Python MuSR data analysis based on classes, with a graphical interface designed for jupyter, released under the MIT licence. 
It aims at the useri-friendly appearance of mulab and the power of musrfit.

v. 3.0 refactoring. See _version.py for the most recent version details. 

Main features: 
- a model built on two-letter bricks: mg, for Gaussian-damped cosine, ml for Lorentzian-damped cosine etc.
- sequential fits by the same model driven by a run list and a list of grouping dicts for asymmetry definition
- global fits by user-defined parameters assigned to model parameters in a json file
- a mulab-like interface in Jupyter-lab that allow fit model and parameter editing, hopefully with a gentler learning curve than musrfit   

To try mujpy:
1. Install python (python3 is assumed)
2. Install git and clone https://github.com/RDeRenzi/mujpy.git

The scripts collected in test.py show how to drive mujpy in command line mode.
To see some fit capabilities try

```
cd example
python test.py 
```

Modify these scripts for your purposes. The clumsy part is to modify by hand the models fit/... .json files. To avoid that

3. Install jupiterlab for the gui editor class mudashed

Docs and installation instructions still are work-in-progress


Old installation instructions (v 1.1)
-------------------------------------
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
