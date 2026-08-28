.. |mSR| replace:: :math:`\mu`\ SR
.. |a| replace:: :math:`\alpha`

Introduction
============
Welcome to **mujpy**  |mSR| data analysis, v.3.0, a set of a few python classes that can be called by command line. Here we assume that you are familiar with the basic technique, with asymmetry and its |a| ratio. This docs provides also an ultrashort, *nuts and bolts* introduction to the physics.
The full aims of this suite of programs are 
        * easy and efficient Minuit fit of muon asymmetries from any facility (limited to PSI and ISIS at present)
        * an intuitive and flexible GUI model editor, both for the simple single dataset fit and for more complex global fits across many datasets and different detector groupings
        * a powerful fit display, based on animations, with optional distinct packing of early and late times, to inspect fast precessions and slow decays at the same time.    
        * simple logs and a csv output for easy plotting of fit parameters
        * saved json fit input for reproducible results

In a sentence, mujpy joins the power of `musrfit <http://lmu.web.psi.ch/musrfit/technical/main.html>`_ by Andreas Suter, a golden standard of |mSR| data analyisis, albeit with a notorious *steep learning curve*, with the ease of the old Matlab `mulab <http://www.fis.unipr.it/~derenzi/dispense/pmwiki.php?n=MuSR.Mulab>`_ interface by R. De Renzi. The main virtue of the latter were

    * an easy way to combine elementary fit components, labelled by two letters: damped cosines (``mg``, ``ml``), Bessel functions (``jg``, ``js``), simple relaxations (``bg``, ``bl``), Kubo-Toyabe functions (``kg``, ``kl``), FMuF functions (``fm``), etc. The model is defined by an acronym that is a combination of two-letter components, e.g. ``mgml`` is the sum of one Gaussian-damped and one Lorentzian-damped cosines.
    * an easy GUI input for sharing or combining parameters among these components 
    * an efficient plot with residues to immediately grasp goodness of fit 


The Mujpy GUI is based on the `jupyterlab <https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html>`_  notebook interface and `ipywidgets <https://ipywidgets.readthedocs.io/en/latest/>`_  

.. sidebar:: Global fit

   A global fit is the minimization of a single cost function, the sum of a cost function for each asymmetry. 

Fits come in eight types, increasing in complexity, conventionally labelled *A1, A20, A21, B1, B20, B21, C1, C2*. *A* fits are single run, *B* are sequential fits on a set of runs, *C1* is a **global** fit of a set of runs. *X20* (with *X* = *A,B*\ ) are sequential fits on two or more detector groupings (e.g. 2-1 and 3-4 on PSI GPS). *X21* are **global** fits of the same set of groupings.

Finally *C2* is a global fit of a set of runs and groupings. 
In the simplest *A21* case the global parameters are first defined by the user, and then assigned to the standard model parameters.
Notice that *C1, C2* can define local virtual parameters, called *hash* parameters, that give rise to a separate replica for each run. This nomenclature is irrelevant in the GUI editor, but it is good to know it anyway.

.. note:: Each of these eight fit types comes in two versions: one, say ``mgml``, with fixed values of the |a| ratios that define a detector grouping asymmetry, and one with the |a| values as  Minuit fit parameters for an automatic calibration, when the data allow it. The acronym of the second version must be ``almgml``, the |a| parameter being formally treated as the first fit component (al), although it is not an additive one. 

Mujpy comes with a test script showing all these 16 subtypes on a standard TF |mSR| data set (a temperature scan approaching a magnetic transition from the paramagnetic side). Try it! Just run ```python test.py``` from subfolder ```example```. Each script in the sequence can be run independently and used as templates for the corresponding fit.

Basic GUI editor usage
----------------------

Using the ```.py``` templates is efficient but has one big drawback: the fit is defined in its json file, automatically stored in the ``fit`` subfolder. It is pretty transparent, but its modification is very error prone. Furthermore, the logic of a global fit is not so easy to plan when concentrating on the json syntax. The solution is the mudashed gui class, to be invoked inside a jupyter notebook. The command sequence is straightforward. If you already installed python, mujpy and jupyter lab, choose your project folder, download the data in subfolder ```data```, launch jupyter lab, launch a new empty notebook andtype 

.. code-block::

    %matplotlib ipympl

in the first cell to activate the correct graphic backend. Then type

.. code-block::

        from mujpy.mudashed import dashed as mudash
        mudash()
  
In the second cell (this is actually the content of Mudashed.ipynb, just start that one if you are lazy). Three rows of widgets appear 

.. image:: Mudashed-0.png

For a single run, single grouping fit, check that Group 0 is OK and press DL to select your data file. Type the run number (no leading zeros) in the run list box and press Enter. A new row o widgets appears. Run title and other info have now appeared in the fist two widget rows. 

.. image:: Mudashed-1.png

Type your chosen model acronym in the box, e.g. ```mgml``` and press ```Enter```. The model editor appears. 

.. image:: Mudashed-2.png 

Just select reasonable guess values and press fit to see a minimization.

If you are a |mSR| beginner refer to a muon primer, such as [Blundell]_, also at `arXiv <https://arxiv.org/abs/cond-mat/0207699>`_ or to a textbook, like [BDLP]_, [Yaouanc]_ or [Amato]_ for this purpose.

You will find more detailed description of the different mujpy methods in the :ref:`reference`. 
Caution: tutorial for v. 3.0 is *work-in-progress*, together with :ref:`Reference`. 

References
----------

.. [BDLP] S.J. Blundell, R. De Renzi, T. Lancaster and F. Pratt, 
   Muon spin spectroscopy, OUP 2021
.. [Blundell] S.J. Blundell, 
   Contemporary Physics 40, 175-192 (1999)
.. [Yaouanc] A. Yaouanc, P. Dalmas De Reotier, 
   Muon spin rotation relaxation and resonance, Oxford University Press, 2011
.. [Amato] A. Amato, E. Morenzoni, 
   Introduction to muon spin spectroscopy, Springer, 2024

