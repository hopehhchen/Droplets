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
import astropy.visualization.wcsaxes as wcsaxes

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
    xcorner = xmin - .75*(xmax-xmin)
    ycorner = ymin - .75*(ymax-ymin)
    width = 2.5*(xmax-xmin)
    height = 2.5*(ymax-ymin)



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


# a plotting function
def plotDroplet(reg, core, list_dictionaries, annotate = True):
    '''
    The function used to plot the droplet for examining the boundary definition.

    Input
    ------
    reg: 'L1688' or 'B18'

    core: The core number.  From 1 to 12 for 'L1688', and from 1 to 6 for 'B18'.
          L1688 has an 'extra' core.

    Output
    ------
    fig: matplotlib figure instantce
    '''

    #reg, core = 'L1688', 6
    dict_data, dict_masks, dict_YSOs, dict_Vlsr_predicted = list_dictionaries
    mask = dict_masks[reg][core]
    header = dict_data[reg]['header_GAS']
    wcs_GAS = wcs.WCS(header)
    frame = frameMask(mask)
    list_images = [dict_data[reg]['colden'],
                   dict_data[reg]['temp'],
                   dict_data[reg]['Tpeak'],
                   dict_data[reg]['Sigma'],
                   dict_data[reg]['Vlsr'],
                   dict_Vlsr_predicted[reg][core]]
    list_names = [r'$N_{H_2}$',
                  r'$T_{dust}$',
                  r'$T_{peak}$',
                  r'$\sigma_{{NH}_3}$',
                  r'$V_{LSR}$',
                  r'$Pred. V_{LSR}$']
    list_norms = [colors.LogNorm(np.nanmedian(list_images[0][mask])/5.,
                                 np.nanmedian(list_images[0][mask])*2.),
                  colors.Normalize(np.nanmedian(list_images[1][mask])-3.,
                                   np.nanmedian(list_images[1][mask])+3.),
                  colors.LogNorm(np.nanmedian(list_images[2][mask])/5.,
                                 np.nanmedian(list_images[2][mask])*2.),
                  colors.Normalize(.05, .45),
                  colors.Normalize(np.nanmedian(list_images[4][mask])-
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[5][mask])),
                                   np.nanmedian(list_images[4][mask])+
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[5][mask]))),
                  colors.Normalize(np.nanmedian(list_images[4][mask])-
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[5][mask])),
                                   np.nanmedian(list_images[4][mask])+
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[5][mask])))]
    list_cmaps = ['Greys',
                  'YlOrRd_r',
                  'Greys',
                  'YlGnBu',
                  'RdYlBu_r',
                  'RdYlBu_r']

    nrows, ncols = 3, 2
    fig = plt.figure(figsize = (20., 20./ncols*nrows/frame[2]*frame[3]))
    figLeft, figRight, figBottom, figTop = .085, .985, .07, .98
    gapHorizontal, gapVertical = .005, .01
    subplotWidth = (figRight-figLeft - gapHorizontal*(ncols-1.))/ncols
    subplotHeight = (figTop-figBottom - gapVertical*(nrows-1.))/nrows

    scalebar = np.array([.02, .05, .1, .2])
    pixscale = (distances[reg]*np.radians(abs(header['CDELT1']))).to(u.pc).value
    scalebar_pix = scalebar/pixscale
    scalebar = scalebar[np.argmin(abs(scalebar_pix-.25*frame[2]))]
    scalebar_pix = scalebar_pix[np.argmin(abs(scalebar_pix-.25*frame[2]))]

    fig.text(.5, .0015, 'R.A.[J2000]',
             color = 'k',
             weight = 'black',
             verticalalignment = 'bottom',
             horizontalalignment = 'center')
    fig.text(.005, .5, 'Dec.[J2000]',
             rotation = 90.,
             color = 'k',
             weight = 'black',
             verticalalignment = 'center',
             horizontalalignment = 'left')

    for i in range(len(list_images)):
        icol, irow = i%ncols, i//ncols
        axis = fig.add_axes([figLeft+icol*(subplotWidth+gapHorizontal),
                            figBottom+(nrows-irow-1)*(subplotHeight+gapVertical),
                            subplotWidth, subplotHeight],
                            projection = wcs_GAS)

        image = list_images[i]
        axis.imshow(image,
                    cmap = list_cmaps[i],
                    norm = list_norms[i])

        count_image = np.sum(np.isfinite(image[int(frame[1]):int(frame[1]+frame[3]),
                                               int(frame[0]):int(frame[0]+frame[2])]))
        count_mask = np.sum(mask)
        if count_image < 3.*count_mask:
            axis.contour(mask,
                         levels = [.5],
                         colors = 'k',
                         linewidths = 8.)
        else:
            axis.contour(mask,
                         levels = [.5],
                         colors = 'w',
                         linewidths = 8.)
        axis.plot(dict_YSOs[reg][:, 0], dict_YSOs[reg][:, 1],
                  color = 'orange',
                  marker = '*',
                  markersize = 48.,
                  markeredgecolor = 'w',
                  markeredgewidth = 2.,
                  linestyle = 'none')


        if i == 3:
            Tkin_median = np.nanmedian(dict_data[reg]['Tkin'][mask])
            NT_sonic = np.sqrt(c.k_B*Tkin_median*u.K/mass['NH3']
                               +c.k_B*Tkin_median*u.K/mass['average'])
            NT_sonic = NT_sonic.to(u.km/u.s).value
            axis.contour((list_images[i] < NT_sonic),
                         levels = [.5],
                         colors = 'r',
                         linewidths = 3.)


        axis.plot([frame[0]+7./8.*frame[2], frame[0]+7./8.*frame[2]-scalebar_pix],
                  [frame[1]+6./7.*frame[3], frame[1]+6./7.*frame[3]],
                  color = 'k',
                  linewidth = 5.)
        axis.text(frame[0]+7./8.*frame[2]-.5*scalebar_pix, frame[1]+4./5.*frame[3], '%.2f pc'%scalebar,
                  color = 'k',
                  weight = 'bold',
                  verticalalignment = 'center',
                  horizontalalignment = 'center')
        axis.fill_between([frame[0]+7./8.*frame[2]+.1*scalebar_pix, frame[0]+7./8.*frame[2]-1.1*scalebar_pix],
                          frame[1]+(6./7.+4./5.)/2.*frame[3] - 1./14.*frame[3],
                          frame[1]+(6./7.+4./5.)/2.*frame[3] + 1./14.*frame[3],
                          color = 'w',
                          linewidth = 0.,
                          alpha = .4)
        '''
        axis.fill_between([frame[0]+7./8.*frame[2]+.1*scalebar_pix, frame[0]+7./8.*frame[2]-1.1*scalebar_pix],
                          frame[1]+(6./7.+4./5.)/2.*frame[3] - 1./14.*frame[3],
                          frame[1]+(6./7.+4./5.)/2.*frame[3] + 1./14.*frame[3],
                          color = 'none',
                          edgecolor = 'k',
                          linewidth = 1.)
        '''
        if annotate:
            corner_max = np.max(image[int(frame[1]):int(frame[1]+.9*frame[3]),
                                      int(frame[0]):int(frame[0]+.1*frame[2])])
            corner_min = np.min(image[int(frame[1]):int(frame[1]+.9*frame[3]),
                                      int(frame[0]):int(frame[0]+.1*frame[2])])

            if (corner_max >= list_norms[i].vmax):
                axis.text(frame[0]+.1*frame[2], frame[1]+.9*frame[3], list_names[i],
                          color = 'w',
                          weight = 'bold',
                          horizontalalignment = 'center',
                          verticalalignment = 'center')
            elif (i in [4, 5]) and (corner_min < list_norms[i].vmin):
                axis.text(frame[0]+.1*frame[2], frame[1]+.9*frame[3], list_names[i],
                          color = 'w',
                          weight = 'bold',
                          horizontalalignment = 'center',
                          verticalalignment = 'center')
            else:
                axis.text(frame[0]+.1*frame[2], frame[1]+.9*frame[3], list_names[i],
                          color = 'k',
                          weight = 'bold',
                          horizontalalignment = 'center',
                          verticalalignment = 'center')

        beam_center = tuple(wcs_GAS.wcs_pix2world([[frame[0]+1./6.*frame[2], frame[1]+1./6.*frame[3]]], 0)[0]*u.deg)
        beam_size = header['BMAJ']/2. * u.degree
        beam = wcsaxes.SphericalCircle(beam_center,
                                       beam_size,
                                       edgecolor = 'w',
                                       facecolor = 'k',
                                       linewidth = 2.,
                                       zorder = 999,
                                       transform = axis.get_transform('fk5'))
        axis.add_patch(beam)

        axis.set_xlim(frame[0], frame[0]+frame[2])
        axis.set_ylim(frame[1], frame[1]+frame[3])
        axis.coords[0].set_major_formatter('hh:mm:ss')
        if irow != (nrows-1):
            axis.coords[0].set_ticklabel_visible(False)
        if icol != 0:
            axis.coords[1].set_ticklabel_visible(False)


    return fig
