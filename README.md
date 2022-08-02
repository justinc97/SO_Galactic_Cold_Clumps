Galactic Cold Clump map generation and frequency scale functions
=======

Galactic Cold Clump map generation and frequency scale script for Simons Observatory using the Planck Galactic Cold Clumps (PGCC) catalogue. Generate HEALPix maps at given Nside with the PGCC objects embedded at their observed galactic latitude and longitude with flux density, scale, shape and orientation given in the PGCC catalogue. 

Source flux can be scaled from default measurement at Planck 353 GHz to frequencies necessary for SO and other CMB observatories.

The `SPoCC.ipynb` workbook is a step by step example of taking hte Planck Galactic Cold Clump Catalogue, observing the distribution of the key parameters and simulating new cold clumps from this information to populate maps at generalised frequencies in I, Q and U for injecting these sources into PySM3 Galactic dust simulations.

The `Test_Notebook.ipynb` workbook guides through the functions available in the module including generation of maps with PGCC sources, adding sources to dust maps to observe their position, shape, scale, orientationm etc., comparing various dust maps such as GNILC, Commander and other PySM3 template dust sky models. The notebook also shows how we can measure the polarisation fraction around these sources.
Finally it demonstrates scaling the flux of sources to new frequencies such as those to be utilised by SO.

## Contains
--------
- `TEST_Notebook.ipynb` reads Planck GNILC maps located in NERSC directory (NERSC access required unless these files are downloaded from https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/foregrounds.html) 
- `PGCC.csv` Obtained from https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd: See https://arxiv.org/pdf/1502.01599.pdf for more information.

If you need to use of anything in this work or have any questions, please contact me.