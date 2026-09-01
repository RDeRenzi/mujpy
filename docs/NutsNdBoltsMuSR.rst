.. |mSR| replace:: :math:`\mu`\ SR
.. |a| replace:: :math:`\alpha`

MuSR, a nuts & bolts approach
=============================

Spin polarized muons decay asymmetrically, 30% more often in the spin direction then opposite to it. 
This is captured by the asymmetry of the decay in two [groups of] detectors placed symmetrically along a given axis, one **forward** (spin pointing to it) and one **backward** (opposite). 
The angle between the initial muon spin and the detector axis determines the decay asymmetry between the two detectors (e.g. if the axis is at 90 deg from the spin, the asymmetry vanishes).
This information is passed to :code:`mujpy` in the **Group** dictionary.

Asymmetry
---------
The experimental asymmetry :math:`A_e(t)` coincides with the normalized difference in forward, :math:`N_f(t)`\ ,  and backward, :math:`N_b(t)` count rates:

.. math::
   :label: experimental

        A_e(t_i) = \frac{N_f(t_i)-\alpha N_b(t_i)}{N_f(t_i)+\alpha N_b(t_i)}
   
Look closely:

   * discretized time :math:`t_i` implies an array of time bins, of fixed bin-width, with a start time :math:`t_i=0`, when the muon stops in the sample and starts to precess  
   * :math:`N_{f,b}(t_i)` imply a *histogram* per detector, an array with average rate at time :math:`t_i` in each bin. 
   * The data file - the **run** - contains the bin-width and the histograms. So the first operation is to load datafiles into :code:`mujpy`, after downloading them according to facility instructions.
   * :math:`\alpha` is a normalizer, correcting for forward and backward detector different experimental geometries/efficiency. It must be **calibrated** with a standard experiment when the sample is mounted, and goes in the :code:`mujpy` **Group** dictionary.
   * Detector rates are decreasing exponentially with time due to the finite lifetime :math:`\tau_\mu` of the muon, but this exponential and the muon lifetime miraculously disappear from :eq:`experimental`, as they must, since the asymmetry is determined uniquely by the muon spin direction at the time of parity-violation weak muon decay. 
     
Fitting the asymmetry
---------------------

The time dependence of the asymmetry :math:`A(t)` during a run is thus due to the ensemble spin evolution of the implanted muons, 100% spin polarized at :math:`t=0`. The evolution provides chemical and physical information on your sample, the implanted muon environment.
In many instances a classical description is sufficient: the average muon spin coherently *precesses* in the local magnetic field :math:`B_\mu` at the muon site (see below) and *relaxes* toward thermodynamic equilibrium (zero spin polarization). So typically the asymmetry contains 
   * damped oscillating functions for coherent precessions, typically :math:`\cos(\gamma_\mu B_\mu t+\phi) g(t)`, where :math:`g(t)` is a relaxation function
   * pure relaxation functions, the simplest being Lorentzian :math:`\exp(-\lambda t)` or Gaussian :math:`\exp(-\sigma^2t^2)`
   * further more complex intermediate or quantum cases

In all cases the model asymmetry may then be written as 

.. math::
   :label: asymmetry 
   
       {\cal A}(t) = A_0 \sum_i f_i G_i(t) 

where

    * :math:`A_0` is the maximum asymmetry, due to sample geometry and other experimental conditions; it must be **calibrated** after mounting each sample.
    *  :math:`f_i` is the fraction of the muon ensemble experiencing the same chemical, physical and geometrical environment (call it improperly **muon site**), i.e. providing the same :math:`G_i(t)` 

Basic implementation
--------------------

:code:`mujpy` produces model functions that correspond exactly to :eq:`asymmetry`

    * components are additive and correspond to damped oscillations, pure relaxations etc.
    * they have a unique two letter acronym and a fixed set of uniquely named parameters
    * an amplitude :math:`A` is always the first parameter of the component


Minimization is performed by :code:`iminuit` on a cost function, typically the :math:`\chi^2`.


.. math::
       
        \chi^2 = \sum_i \left(\frac{A_e(t_i)-{\cal A}(t_i)}{\sigma(t_i)}\right)^2

where :math:`\sigma(t_i)` is the standard deviation of the experimental asymmetry at time :math:`t_i`, assuming Poisson statistics. The model function can be further refined:

    * parameters can be shared among components. E.g. often different oscillating fractions share the same initial phase :math:`\phi`, and this is allowed by model parameter **flag** and **function**, one component defines the parameter and the second *shares* it.

Finally, :code:`mujpy` performs global fits on multiple asymmetries, minimizing the sum of their cost functions. This makes sense if the different asymmetries have common parameters. For example:

    * two detector groups, :math:`a,b` generally do not share the same initial asymmetry :math:`A_a,A_b`, but they may share the same muon fraction :math:`f`. A global fit may define :math:`A_a,A-b,f` as global parameters and assign them through **flag** and **function** to the model component amplitudes as :math:`A_1 f, A_2 f`.
    * This is convenient for two components or more, the amplitudes of the second being :math:`A_1(1-f),A_2(1-f)`: three parameters instead of four, here.
    * an obvious example is a global fit for a **list** of runs, i.e. asymmetries on the same sample and geometry, while scanning, say, temperature. The model, the maximum asymmetry, the phase and a few other parameters  may be the same across the suite, while other interesting parameters depend on the individual run.  

This nuts and bolts description illustrates some of the capabilities of :code:`mujpy` and defines the key concepts for its practical usage: 

.. code::

   datafile, group, alpha, run list, bin-width, bin range, 
   model components, model parameters, flag, function, global parameters

 


 

