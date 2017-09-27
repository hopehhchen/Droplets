import os

#
import numpy as np

#
import astropy.units as u
import astropy.constants as c
from astropy.coordinates import SkyCoord, FK4, FK5

#
import pandas as pd


#### conversion is done below

# Read the original tables.
direcData = os.getcwd()+'/data/Goodman93/'
table1 = pd.read_csv(direcData+'table1.csv',
                     index_col = False)
table2 = pd.read_csv(direcData+'table2.csv',
                     index_col = False)

# Convert Table 1: #####
table1_updated = table1.copy()

# distance
## The references are in the jupyter notebook.
table1_updated['Distance'] = pd.Series([260., 260., 315., 135., 135., 135., 135., 135., 135., 135.,
                                        170., 135., 135., 170., 135., 135., 135., 135., 135., 420.,
                                        420., 110., 137.3, 125., 125., 125., 125., 125., 125., 125.,
                                        200., 200., 325., 400., 400., 400., 288., 288., 350., 900.,
                                        286., 286., 286.])
## uncertainty in D (ref in the notebook)
table1_updated['eDistance'] = pd.Series([26., 26., 32., 20., 20., 20., 20., 20., 20., 20.,
                                         50., 20., 20., 50., 20., 20., 20., 20., 20., 42.,
                                         42., 10., 6., 45., 45., 45., 45., 45., 45., 45.,
                                         np.nan, np.nan, 13., 50., 50., 50., 25., 25., np.nan, np.nan,
                                         20., 20., 20.])

# coordinates
## RA
table1_updated['RA'] = pd.Series([SkyCoord(ra, dec,
                                           frame = FK4,
                                           equinox = 'B1950').transform_to(FK5(equinox = 'J2000')).ra.deg
                                  for ra, dec in zip(table1['RA'], table1['Dec'])])
## Dec
table1_updated['Dec'] = pd.Series([SkyCoord(ra, dec,
                                            frame = FK4,
                                            equinox = 'B1950').transform_to(FK5(equinox = 'J2000')).dec.deg
                                   for ra, dec in zip(table1['RA'], table1['Dec'])])

# Convert FWHMTot to SigmaTot.
table1_updated['SigmaTot'] = table1['FWHMTot']/(2.*np.sqrt(2.*np.log(2.)))


# mass
table1_updated['M'] = table1['M']*(2.37/2.33)*(table1_updated['Distance']/table1['Distance'])**3.

# sizes
## major FWHM
table1_updated['major'] = table1['major']*(table1_updated['Distance']/table1['Distance'])
## minor FWHM
table1_updated['minor'] = table1['minor']*(table1_updated['Distance']/table1['Distance'])
## R
table1_updated['R'] = np.sqrt(table1_updated['major']*table1_updated['minor'])


# Tkin and SigmaNH3 from Benson & Myers 1989 and Ladd et al. 1994
## Tkin
table1_updated['Tkin'] = pd.Series([11., 11., 10.8, 9.5, 10., np.nan, 10., np.nan, np.nan, np.nan,
                                    10.2, 15., 8.8, 10.8, 10., 8.0, 11.1, 8.9, 11.1, 14.7,
                                    16.6, 9.0, 15., np.nan, 8.5, 8.8, 11.6, 14.5, 10.3, 9.7,
                                    np.nan, 11.2, 10.4, 9.2, 10., 10., 10., 14.2, 10., 12.5,
                                    10.5, np.nan, 10.3])
## SigmaNH3
table1_updated['SigmaNH3'] = pd.Series([.39, .31, .34, .28, .21, .23, .29, np.nan, np.nan, np.nan,
                                        .30, .34, .46, .26, .32, .34, .28, .24, .23, .44,
                                        .82, .32, .26, np.nan, .32, .21, .29, .30, .23, .26,
                                        np.nan, .40, .33, .36, .40, .27, .74, .43, .70, .94,
                                        .37, np.nan, .37])/(2.*np.sqrt(2.*np.log(2.)))
table1_updated['SigmaNH3'] = np.around(table1_updated['SigmaNH3'],
                                       decimals = 2)


# Reorder the columns.
col_list = ['ID', 'Distance', 'eDistance', 'RA', 'Dec',
            'Vlsr', 'SigmaTot', 'M', 'major', 'minor',
            'R', 'PA', 'AspectRatio', 'Tkin', 'SigmaNH3']
table1_updated = table1_updated[col_list]


# Save to the folder.
direcData = os.getcwd()+'/data/Goodman93/'
table1_updated.to_csv(direcData+'table1_updated.csv',
                      index = False)


# Convert Table 2: #####
table2_updated = table2.copy()

# gradient
table2_updated['Gradient'] = table2['Gradient']*(table1_updated['Distance']/table1['Distance'])**-1.
## uncertainty in G
table2_updated['eGradient'] = table2['eGradient']*(table1_updated['Distance']/table1['Distance'])**-1.

# beta
table2_updated['beta'] = table2['beta']*(2.37/2.33)**-1.*(table1_updated['Distance']/table1['Distance'])**-2.

# JoverM
## Note that J/M was calculated from ~ pGR^2, and does not depend on the mass.
table2_updated['JoverM'] = table2['JoverM']*(table1_updated['Distance']/table1['Distance'])


# Save to the folder.
direcData = os.getcwd()+'/data/Goodman93/'
table2_updated.to_csv(direcData+'table2_updated.csv',
                      index = False)
