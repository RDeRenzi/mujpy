.. _intro:

Introduction
============
Welcome to mujpy v.2.7.3, the `jupyterlab <https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html>`_  notebook interface for muon spin spectroscopy analysis. 
This suite is based loosely on the general layout of the Matlab `mulab <http://www.fis.unipr.it/~derenzi/dispense/pmwiki.php?n=MuSR.Mulab>`_ interface and on contamination of ideas with `musrfit <http://lmu.web.psi.ch/musrfit/technical/main.html>`_ by Andreas Suter.  


Mujpy has just been refactored in version 2.0. As of v. 2.7.3 it comes with:

    * a number of "Tst_mrun..." jupyter notebook examples that can be used as templates for existing categories of fits
        * TF :math:`\alpha` calibration on single run
        * fit on a single run
        * sequential fits on a suite of runs
        * sequential fits on a single run, distinct groups of detectors
        * global fit on a single run, distinct groups of detectors
        * sequential fits on multiple runs
        * global fit on multiple runs, single group of detectors
        * yet to come global fit on multiple runs, multiple detectors  

It comes with a Graphical User Interface (gui), or Dashboard, mimicking and extending mulab, meant to be intuitive for the experienced MuSR user. 

    * a smaller number of MuDash jupyter notebook examples to start the gui
    * the Dashboard features unobtrusive tiptools: hovering with the mouse overt the descriptive text or the buttons additional instructions appear.  

.. image:: setup.png

You must be familiar with standard MuSR concepts such as the :math:`\alpha` calibration factor, or the muon asymmetry. 

Refer to a muon primer, such as [Blundell]_, also at `arXiv <https://arxiv.org/abs/cond-mat/0207699>`_ or to a textbook, like [BDLP]_ or [Yaouanc]_ for this purpose.

Caution: tutorial for v. 2.7.3 is *work-in-progress*. The example notebooks themselved clarify the way to use mujpy and it is briefly commented in :ref:`Reference`. 
 
You can find more detailed description of the different mujpy methods in the :ref:`reference`. 

.. [BDLP] S.J. Blundell, R. De Renzi, T. Lancaster and F. Pratt, Muon spin spectroscopy, OUP 2021
.. [Blundell] S.J. Blundell, Contemporary Physics 40, 175-192 (1999)
.. [Yaouanc] A. Yaouanc and P. Dalmas de Reotier, MUON SPIN ROTATION, RELAXATION and RESONANCE, Oxford University Press, 2011


