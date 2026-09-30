.. |mSR| replace:: :math:`\mu`\ SR
.. |a| replace:: :math:`\alpha`

Introduction
============
Welcome to **mujpy**  |mSR| data analysis, v.3.1, based on a small set of command-line python classes, with an independent gui interface. Basic  |mSR| concepts are recalled in :doc:`NutsNdBoltsMuSR` for the aim of this package.

The classes provide 
        * simultaneous access to sets of asymmetries on the same sample: list of runs and/or list of detector groups.
        * Minuit fit of muon asymmetries from PSI and ISIS (add your own)
        * an intuitive GUI model editor, for simple single dataset fit and more complex global fits of many datasets and detector groups
        * a fit display based on animations, optional distinct packing of early and late times, to show fast and slow evolution.    
        * logs and csv output for plotting of fit parameters
        * reproducible by saved json fit input

`mujpy` merges some of the best of `musrfit <http://lmu.web.psi.ch/musrfit/technical/main.html>`_ by Andreas Suter, the *steep learning curve* standard, with the old intuitive `mulab <http://www.fis.unipr.it/~derenzi/dispense/pmwiki.php?n=MuSR.Mulab>`_\ (Matlab) approach by R. De Renzi. The ``mulab`` concepts inherited by ``mujpy`` are

    * elementary fit components, labelled by two letters: damped cosines (``mg``, ``ml``), Bessel functions (``jg``, ``js``), simple relaxations (``bg``, ``bl``), Kubo-Toyabe functions (``kg``, ``kl``), FMuF functions (``fm``), etc. 
    * model selected  by acronym, combining the two-letter components: ``mgml`` is the sum of a Gaussian- and a Lorentzian-damped cosine.
    * GUI model editor for combining parameters among components, based on jupyter notebooks interface and `ipywidgets <https://ipywidgets.readthedocs.io/en/latest/>`_  
    * plot with residues to visualize goodness of fit 


.. sidebar:: Global fit

   A global fit is the minimization of a single cost function, the sum of a cost function for each asymmetry. 

Conceptually, labels *A1, A20, A21, B1, B20, B21, C1, C2* distinguish 8 types of fit: *A* fits are single run, *B* are sequential fits on a set of asymmetries, *C1* is a **global** fit of a set of runs. *X20* fits (with *X* = *A,B*\ ) are also group-sequential, on two or more detector groups (e.g. 2-1 and 3-4 on PSI GPS). *X21* are **global** fits of the same groups. Finally, *C2* is a global fit across both runs and groups.  

Luckily, the GUI just selects sequential or global fits, but the code is aware of this nomenclature.
As in ``musrfit``, global parameters are defined in advance. They are the sole Minuit fit parameters. Each model parameter is assigned to one of them. 
Furthermore, *C1, C2* fit mark a few global parameters as *hash* (virtual). They give rise to a distict Minuit replica for each run. 

.. note:: Each of these 8 fit types comes in two versions: one, say ``mgml``, with fixed values of the |a| ratio that defines a detector grouping asymmetry, and one with the |a| values as Minuit fit parameters, for automatic calibration, when the data allow it. The acronym of the second version must be ``almgml``, the |a| parameter being formally treated as the first fit component ``al``, although it is not an additive one. 

The package includes a ```test.py``` script demonstrating all these 16 subtypes on a standard TF |mSR| data set (a temperature scan approaching a magnetic transition from the paramagnetic side). Try it! Just run ```python test.py``` from subfolder ```example```. Each script in the sequence can be run independently and used as templates for the corresponding fit type.

The ```.py``` templates have one BIG! drawback: the fit is defined in a json file, automatically stored in the ``fit`` subfolder. It is transparent, but its modification is very error prone. Furthermore, planning global fits dows not agree with ckecking the json syntax. The `dashed` editor, is a jupyter notebook GUI that makes it much easier, if you proceed to the Tutorial.

If you are a |mSR| beginner, :doc:`NutsNdBoltsMuSR` is not enough, please refer to a muon primer, such as [Blundell]_, also at `arXiv <https://arxiv.org/abs/cond-mat/0207699>`_ or to a textbook, like [BDLP]_, [Yaouanc]_ or [Amato]_ for this purpose.

lease find details on the `mujpy` methods in the :ref:`Reference` (work in progress). 

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

