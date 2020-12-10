*****
MuJPy
*****

A Python MuSR data analysis graphical interface, based on classes, designed for jupyter.

Released under the MIT licence.

Linux instructions, for the time being. 
Valid on WIN10 that has a linux shell.

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

Documentatation on the GUI usage at http://mujpy.readthedocs.io/en/latest/
