import sys
import warnings

#
import numpy as np

#
from astropy.io import fits
import astropy.wcs as wcs
import astropy.units as u
import astropy.constants as c
import astropy.modeling as modeling

#
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#
import pandas as pd

#
from constants import *
import styles



def centroidMask(mask):
    '''
    The function used to calculate the centroid of the boolean mask.

    Input
    ------
    mask: a 2D-array, either of type boolean, or any type that can be converted
    to the boolean.

    Output
    ------
    xcent, ycent: (x, y) [axis-1 and axis-0]
    '''

    # Create the coordinate grids.
    mask = mask.astype(bool)
    xgrid, ygrid = np.meshgrid(np.arange(mask.shape[1], dtype = float),
                               np.arange(mask.shape[0], dtype = float))

    # Calculate the centroid based on the boolean mask.
    xcent = np.average(xgrid[mask],
                       weights = mask[mask].astype(float))
    ycent = np.average(ygrid[mask],
                       weights = mask[mask].astype(float))



    return xcent, ycent


def frameMask(mask):
    '''
    The function used to determine the frame in which the structure is centered.
    This is mainly for plotting.

    Input
    ------
    mask: a 2D-array, either of type boolean, or any type that can be converted
    to the boolean.

    Output
    ------
    xcorner, ycorder: bottom left (x, y) [axis-1 and axis-0] of the frame
    width, height: of the frame
    '''

    # Create the coordinate grids.
    mask = mask.astype(bool)
    xgrid, ygrid = np.meshgrid(np.arange(mask.shape[1], dtype = float),
                               np.arange(mask.shape[0], dtype = float))

    # Calculate the extent of the mask.
    xmin, xmax = np.min(xgrid[mask]), np.max(xgrid[mask])
    ymin, ymax = np.min(ygrid[mask]), np.max(ygrid[mask])

    # Calculate the frame based on the extent
    xcorner = xmin - .5*(xmax-xmin)
    ycorner = ymin - .5*(ymax-ymin)
    width = 2.*(xmax-xmin)
    height = 2.*(ymax-ymin)



    return xcorner, ycorner, width, height




def fitGradient(mask, Vlsr, eVlsr):
    '''
    The function used to fit a 1st-degree 2D polynomial to the Vlsr field.

    Input
    ------
    mask:

    Vlsr:

    Output
    ------
    gradfit: an astropy.modeling.model object

    fitter: an astropy.modeling.fitting object

    Vlsr_predicted:
    '''

    # Create the coordinate grid; shift acoording to the centroid; read Vlsr and eVlsr.
    mask = mask.astype(bool)
    xcent, ycent = centroidMask(mask)
    xgrid, ygrid = np.meshgrid(np.arange(mask.shape[1], dtype = float),
                               np.arange(mask.shape[0], dtype = float))
    xgrid -= xcent
    ygrid -= ycent
    ## Read Vlsr and eVlsr.
    zgrid = Vlsr
    wgrid = 1./eVlsr**2. ## Weight fitting by the reciprocal of uncertainty squared.

    # Fit using `astropy.modeling`.
    gradfit = modeling.polynomial.Polynomial2D(1)
    fitter = modeling.fitting.LevMarLSQFitter()
    gradfit = fitter(gradfit,
                     xgrid[mask&np.isfinite(Vlsr)],
                     ygrid[mask&np.isfinite(Vlsr)],
                     Vlsr[mask&np.isfinite(Vlsr)],
                     weights = wgrid[mask&np.isfinite(Vlsr)])

    # Generate a map of predicted Vlsr.
    Vlsr_predicted = gradfit(xgrid, ygrid)



    return gradfit, fitter, Vlsr_predicted

def convertAngle(angle):
    '''
    The function that converts the numpy angles to PA [E of N] and remove
    redundancies.

    Input
    ------
    angle: in degrees.
    '''

    # Rotate by 90 degrees counter-clockwise.
    angle = angle - 90.

    # Remove degeneracy around the 180 degree mark.
    if angle <= 0.:
        angle = angle + 180.

    # Convert angles in the SE quadrant to positive values.
    if angle <= 0.:
        angle = angle + 180.



    return angle


def readGradient(gradfit, fitter, header, reg):
    '''
    The function that converts the `astropy.modeling` object to physical values.

    Input
    ------
    gradfit:

    Output
    ------
    GradMag: The gradient magnitude, in km/s/pc. [`astropy.units` object]

    eGradMag: The uncertainty in the gradient magnitude measurement.

    GradPA: The position angle of the fitted velocity graidnet, in degrees.
            [E of N]

    eGradPA: The uncertainty in the position angle measurement.
    '''

    # Calculate the pixel scale (corresponding physical scale at the region distance).
    pixscale = np.radians(abs(header['CDELT1']))*distances[reg]

    # Convert the gradients to physical units (km/s/pc).
    gradx = gradfit.parameters[1]*u.km/u.s/pixscale
    grady = gradfit.parameters[2]*u.km/u.s/pixscale

    # Calculate the magnitude and the PA based on the converted x- and y-components.
    GradMag = (np.sqrt(gradx**2.+grady**2.)).to(u.km/u.s/u.pc).value  ## in km/s/pc
    GradPA = convertAngle(np.degrees(np.arctan2(grady.value, gradx.value)))  ## in degrees

    # Estimate the uncertainty from the covariant matrix.
    ## the raw parameters
    x, y = gradfit.parameters[1], gradfit.parameters[2]
    ## uncertainties in x and y
    sigx = np.sqrt(fitter.fit_info['param_cov'][1, 1])
    sigy = np.sqrt(fitter.fit_info['param_cov'][2, 2])
    ## propagation to GradMag
    eGradMag = np.sqrt((sigx*x/np.sqrt(x**2.+y**2.))**2.+(sigy*y/np.sqrt(x**2.+y**2.))**2.)
    eGradMag *= u.km/u.s/pixscale  ## Convert to physical units.
    eGradMag = eGradMag.to(u.km/u.s/u.pc).value  ## in km/s/pc
    ## propagation to GradPA
    eGradPA = np.sqrt((sigx*(1./(1.+(y/x)**2.))*(-y/x**2.))**2.+\
                      (sigy*(1./(1.+(y/x)**2.))*(1./x))**2.)
    eGradPA = np.degrees(eGradPA)  ## Convert to degrees.




    return GradMag, eGradMag, GradPA, eGradPA
