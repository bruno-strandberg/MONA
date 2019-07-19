Example macros
===============

Prerequisities
---------------

What is required to run the macros?

1) MONA needs to be compiled and set-up, see https://bstrandberg.pages.km3net.de/MONA/index.html, sections `Prerequisities` and `Setup in lyon and elsewhere`.
2) The MONA MC files at iRODS location /in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_20m_9m/v5.0/postprocessing/MONA_190222_v1.0/dataprocessing_ECAP190222_20m need to be downloaded and fed to the scripts.

How to run
-----------

~~~
root
.x llh_scan.C+("ORCA_MCsummary_SEv2_ORCA115_20x9m_ECAP190222.root", "EffMass_ORCA115_20x9m_ECAP190222.root")
.x sim_fit.C+("ORCA_MCsummary_SEv2_ORCA115_20x9m_ECAP190222.root", "EffMass_ORCA115_20x9m_ECAP190222.root")
~~~

Short description
-----------------

* `llh_scan.C`: this macro performs a likelihood scan in theta-23
* `llh_scan.py`: this macro does exactly the same as `llh_scan.C`, but in python
* `sim_fit.C`: this macro demonstrates how to perform simultaneous fits in track and shower channel