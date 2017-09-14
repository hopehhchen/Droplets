# Droplets
Repo for an upcoming paper analyzing gravitationally unbound, coherent structures with significant velocity gradients in nearby star forming molecular clouds.  The repo is prepared to work with [binder](http://mybinder.org) and will be shared alongside the paper on Authorea.  By working between Github and other data/code/article sharing services, we hope to present what a journal article could look like in an era of open data and reproducible/reusable codes.

Please contact Hope Chen at hhchen@cfa.harvard.edu if you have any questions.  Any suggestions are also welcome via [the issue tracker of this repo](https://github.com/hopehhchen/Droplets/issues).

## Abstract
Using data from the [*GBT Ammonia Survey* (GAS) Data Release 1](https://dataverse.harvard.edu/dataverse/GAS_Project) ([Friesen and Pineda et al., 2017](https://ui.adsabs.harvard.edu/#abs/2017ApJ...843...63F/abstract)), we look into the physical properties of structures that are coherent (previous examples in [Goodman et al., 1998](https://ui.adsabs.harvard.edu/#abs/1998ApJ...504..223G/abstract) and in [Pineda et al, 2010](https://ui.adsabs.harvard.edu/#abs/2010ApJ...712L.116P/abstract)) and show significant velocity gradients (previous examples in [Goodman et al., 1993](https://ui.adsabs.harvard.edu/#abs/1993ApJ...406..528G/abstract)).  With a much improved physical resolution of ~4000 AU (compared to ~0.1 pc in the 1990s) at the distance of Ophiuchus and Taurus, one goal of the analysis is to provide updated numbers for properties of (potentially) rotational motions within 0.1 pc, which has been essential in simulations and analytic models alike in the past two decades.  In our first analysis, these cores seem to be gravitationally unbound and are termed "droplets."  In the paper, we hope to provide a sensible guess to what role these structures could play in the process of star formation.

## Data from GAS DR1 and Herschel
The main dataset used in the analyses are from the GAS DR1 ([Friesen and Pineda et al., 2017](https://ui.adsabs.harvard.edu/#abs/2017ApJ...843...63F/abstract)) and is hosted on [Dataverse](https://dataverse.harvard.edu/dataverse/GAS_Project).  The *Herschel* column density maps are derived by Ayushi Singh and Peter Martin at University of Toronto, using archival data from the *Herschel Space Observatory*.  The data in this repo are copies of the data hosted on Dataverse, without any modification.

Due to the Github policy, data files larger than 100 MB are not hosted here.  These include the raw data cubes and the position-position-velocity cubes from velocity decomposition based on Gaussian line fitting to the ammonia emission lines.  These large files are not necessary to run codes in this repo.  Please look in the GBT DR1 dataverse for the files.

## Status of the Project and the Github Repo
A manuscript is being prepared on Authorea.  Meanwhile, results of analyses are actively shared in this repo, together with the code used for these analyses.  Currently, the codes are *not* yet stable enough to be applied on other datasets.  Please approach with caution.

## Collaboration
The project is lead by [Hope Chen](https://github.com/hopehhchen) at Harvard-Smithsonian Center for Astrophysics, and is a collaboration with [Jaime Pineda](https://github.com/jpinedaf), Alyssa Goodman, and Andreas Burkert.
