import sys
import warnings

#
import numpy as np
import scipy

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
import matplotlib.ticker as ticker
from matplotlib import rcParams

#
import pandas as pd

#
from statBasic2D import *
from constants import *
from ssk_colors import *
from plot_tools import *
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


#### blow are plotting functions; consider moving to a separate script for orga
# a plotting function
def plotDroplet(reg, core, list_dictionaries, annotate = True):
    '''
    The function used to plot the droplet for examining the boundary definition.

    Input
    ------
    reg: 'L1688' or 'B18'

    core: The core number.  From 1 to 12 for 'L1688', and from 1 to 6 for 'B18'.
          L1688 has an 'extra' core.

    list_dictionaries: list of data dictionaries in the order of dict_data,
                       dict_masks, dict_YSOs, and dict_Vlsr_predicted.

    Output
    ------
    fig: matplotlib figure instantce
    '''


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
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[4][mask])),
                                   np.nanmedian(list_images[4][mask])+
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[4][mask]))),
                  colors.Normalize(np.nanmedian(list_images[4][mask])-
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[4][mask])),
                                   np.nanmedian(list_images[4][mask])+
                                   (np.nanmax(list_images[4][mask])-np.nanmin(list_images[4][mask])))]
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

    if annotate:
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
            Tkin_median = np.nanmedian(dict_data[reg]['Tkin'][int(frame[1]):int(frame[1]+frame[3]),
                                                              int(frame[0]):int(frame[0]+frame[2])])
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
        axis.coords[0].set_ticks(size = 8.)
        axis.coords[1].set_ticks(size = 8.)
        if irow != (nrows-1):
            axis.coords[0].set_ticklabel_visible(False)
        if icol != 0:
            axis.coords[1].set_ticklabel_visible(False)


    return fig


# a plotting function
def plotRegion(reg, list_dictionaries, chooseStructure = None, annotate = True):
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

    if reg == 'L1688' and chooseStructure in list(range(1, 13))+['extra']:
        core = chooseStructure
        plotStructure = True
    elif reg == 'B18' and chooseStructure in range(1, 7):
        core = chooseStructure
        plotStructure = True
    elif chooseStructure is None:
        core = -1
        plotStructure = False
    else:
        raise ValueError('"chooseStructure" is an integer. 1-12 for L1688; 1-6 for B18.')


    ####
    dict_data, dict_masks, dict_YSOs, dict_Vlsr_predicted = list_dictionaries
    header = dict_data[reg]['header_GAS']
    wcs_GAS = wcs.WCS(header)
    if reg == 'L1688':
        frame = (7., 0., 261., 196.2)
    elif reg == 'B18':
        frame = (26., -10., 725., 290.)

    if plotStructure:
        mask = dict_masks[reg][core]
        frameCore = frameMask(mask)

    list_images = [dict_data[reg]['colden'],
                   dict_data[reg]['temp'],
                   dict_data[reg]['Tpeak'],
                   dict_data[reg]['Tkin'],
                   dict_data[reg]['Vlsr'],
                   dict_data[reg]['Sigma']]
    list_names = [r'$N_{H_2}$',
                  r'$T_{dust}$',
                  r'$T_{peak}$',
                  r'$T_{kin}$',
                  r'$V_{LSR}$',
                  r'$\sigma_{{NH}_3}$']
    norm_Vlsr = colors.Normalize(2.5, 4.5) if reg == 'L1688'\
                else colors.Normalize(5.5, 7.)
    list_norms = [colors.LogNorm(1e21, 1e23),
                  colors.Normalize(0., 30.),
                  colors.LogNorm(.5, 30.),
                  colors.Normalize(0., 30.),
                  norm_Vlsr,
                  colors.Normalize(.05, .45)]
    list_cmaps = ['Greys',
                  'YlOrRd_r',
                  'Greys',
                  'YlOrRd_r',
                  'RdYlBu_r',
                  'YlGnBu']

    nrows, ncols = 3, 2
    if reg == 'L1688':
        fig = plt.figure(figsize = (16., 18.))
        figLeft, figRight, figBottom, figTop = .085, .985, .07, .98
        listStructures = list(range(1, 13))+['extra']
        markersizeYSOs = 8.
    elif reg == 'B18':
        fig = plt.figure(figsize = (20., 12.))
        figLeft, figRight, figBottom, figTop = .085, .985, .09, .99
        listStructures = range(1, 7)
        markersizeYSOs = 11.

    gapHorizontal, gapVertical = .005, .01
    subplotWidth = (figRight-figLeft - gapHorizontal*(ncols-1.))/ncols
    subplotHeight = (figTop-figBottom - gapVertical*(nrows-1.))/nrows

    scalebar = .5 if reg == 'L1688' else 1.  ## pc
    pixscale = (distances[reg]*np.radians(abs(header['CDELT1']))).to(u.pc).value
    scalebar_pix = scalebar/pixscale

    if annotate:
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

        for j, structure in enumerate(listStructures):
            axis.contour(dict_masks[reg][structure],
                         levels = [.5],
                         colors = 'w',
                         linewidths = 3.)
            axis.contour(dict_masks[reg][structure],
                         levels = [.5],
                         colors = ssk_colors[j],
                         linewidths = 2.)

        '''
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
        '''

        if i in [1, 3]:
            axis.plot(dict_YSOs[reg][:, 0], dict_YSOs[reg][:, 1],
                      color = 'orange',
                      marker = '*',
                      markersize = markersizeYSOs,
                      markeredgecolor = 'k',
                      linestyle = 'none')
        else:
            axis.plot(dict_YSOs[reg][:, 0], dict_YSOs[reg][:, 1],
                      color = 'orange',
                      marker = '*',
                      markersize = markersizeYSOs,
                      markeredgecolor = 'w',
                      linestyle = 'none')

        if plotStructure:
            if i in [1, 3, 4]:
                axis.fill_between([frameCore[0], frameCore[0]+frameCore[2]],
                                  frameCore[1], frameCore[1]+frameCore[3],
                                  edgecolor = 'w',
                                  color = 'none',
                                  linewidth = 3.)
                axis.fill_between([frameCore[0], frameCore[0]+frameCore[2]],
                                  frameCore[1], frameCore[1]+frameCore[3],
                                  edgecolor = 'k',
                                  color = 'none',
                                  linewidth = 2.)
            else:
                axis.fill_between([frameCore[0], frameCore[0]+frameCore[2]],
                                  frameCore[1], frameCore[1]+frameCore[3],
                                  edgecolor = 'w',
                                  color = 'none',
                                  linewidth = 3.)
                axis.fill_between([frameCore[0], frameCore[0]+frameCore[2]],
                                  frameCore[1], frameCore[1]+frameCore[3],
                                  edgecolor = 'k',
                                  color = 'none',
                                  linewidth = 2.)

        '''
        if i == 5:
            Tkin_median = np.nanmedian(dict_data[reg]['Tkin'][mask])
            NT_sonic = np.sqrt(c.k_B*Tkin_median*u.K/mass['NH3']
                               +c.k_B*Tkin_median*u.K/mass['average'])
            NT_sonic = NT_sonic.to(u.km/u.s).value
            axis.contour((list_images[i] < NT_sonic),
                         levels = [.5],
                         colors = 'r',
                         linewidths = 2.)
        '''


        axis.plot([frame[0]+7./8.*frame[2], frame[0]+7./8.*frame[2]-scalebar_pix],
                  [frame[1]+1./7.*frame[3], frame[1]+1./7.*frame[3]],
                  color = 'k',
                  linewidth = 5.)
        axis.text(frame[0]+7./8.*frame[2]-.5*scalebar_pix, frame[1]+1./5.*frame[3], '%.1f pc'%scalebar,
                  color = 'k',
                  weight = 'bold',
                  verticalalignment = 'center',
                  horizontalalignment = 'center')
        axis.fill_between([frame[0]+7./8.*frame[2]+.1*scalebar_pix, frame[0]+7./8.*frame[2]-1.1*scalebar_pix],
                          frame[1]+(1./7.+1./5.)/2.*frame[3] - 1./14.*frame[3],
                          frame[1]+(1./7.+1./5.)/2.*frame[3] + 1./12.*frame[3],
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
                                       zorder = 999,
                                       transform = axis.get_transform('fk5'))
        axis.add_patch(beam)

        if annotate:
            axis.text(frame[0]+1./6.*frame[2], frame[1]+1./6.*frame[3]+20., 'beam',
                      style = 'italic',
                      weight = 'bold',
                      horizontalalignment = 'center',
                      verticalalignment = 'center')

        axis.set_xlim(frame[0], frame[0]+frame[2])
        axis.set_ylim(frame[1], frame[1]+frame[3])
        axis.coords[0].set_major_formatter('hh:mm')
        axis.coords[1].set_major_formatter('dd:mm')
        axis.coords[0].set_ticks(size = 8.)
        axis.coords[1].set_ticks(size = 8.)
        if irow != (nrows-1):
            axis.coords[0].set_ticklabel_visible(False)
        if icol != 0:
            axis.coords[1].set_ticklabel_visible(False)


    return fig


# a plotting function
def plotTpeakSigma(list_dictionaries, xscale = 'log'):

    ###
    rcParams['figure.subplot.left'] = .09
    rcParams['figure.subplot.right'] = .97
    rcParams['figure.subplot.bottom'] = .12
    rcParams['figure.subplot.top'] = .96
    rcParams['font.size'] = 30
    if xscale == 'log':
        rcParams['xtick.major.pad'] = 11
    elif xscale == 'linear':
        rcParams['xtick.major.pad'] = 5
    ###
    dict_data, dict_masks, dict_YSOs, dict_Vlsr_predicted = list_dictionaries

    fig, ax = plt.subplots(figsize = (14., 7.),
                           ncols = 2)
    if xscale == 'log':
        xmin, xmax = .2, 8.5
    elif xscale == 'linear':
        xmin, xmax = .03, 5.
    ymin, ymax = .02, .92

    #
    Sigma_SonicNT = (np.sqrt(c.k_B*10.*u.K/(17.031*u.u)+c.k_B*10.*u.K/(2.37*u.u))).to(u.km/u.s).value
    Sigma_halfSonicNT = (np.sqrt((c.k_B*10.*u.K/(17.031*u.u)+.5**2.*c.k_B*10.*u.K/(2.37*u.u)))).to(u.km/u.s).value

    for i, reg in enumerate(['L1688', 'B18']):

        # list of structures
        if reg == 'L1688':
            listStructures = list(range(1, 13))+['extra']
        elif reg == 'B18':
            listStructures = range(1, 7)

        ## load the data
        mapX = dict_data[reg]['Tpeak']
        mapY = dict_data[reg]['Sigma']
        ### masking
        maskData = np.isfinite(mapX) & np.isfinite(mapY)
        ### all points and contours of their distribution
        if xscale == 'log':
            hist2D = np.histogram2d(mapX[maskData], mapY[maskData],
                                    [np.logspace(np.log10(xmin), np.log10(xmax), 20),
                                     np.linspace(ymin, ymax, 20)])
            xBinCent = 10.**(np.log10(hist2D[1])[:-1]+.5*np.diff(np.log10(hist2D[1])))
            yBinCent = hist2D[2][:-1]+.5*np.diff(hist2D[2])
        elif xscale == 'linear':
            hist2D = np.histogram2d(mapX[maskData], mapY[maskData],
                                    [np.linspace(xmin, xmax, 20),
                                     np.linspace(ymin, ymax, 20)])
            xBinCent = hist2D[1][:-1]+.5*np.diff(hist2D[1])
            yBinCent = hist2D[2][:-1]+.5*np.diff(hist2D[2])

        height2D = hist2D[0]/np.sum(hist2D[0])
        n = 5000
        ## accumulative contours
        t = np.linspace(0., height2D.max(), n)
        integral = ((height2D >= t[:, None, None]) * height2D).sum(axis=(1,2))
        f = scipy.interpolate.interp1d(integral, t)
        contourLevels = f(np.array([.95, .75, .5, .25]))


        ## plot points within the cores
        axis = ax[i]
        mask_all = np.zeros(maskData.shape, dtype = bool)
        for j, structure in enumerate(listStructures[::-1]):
            ### mask of the core
            mask_core = maskData & dict_masks[reg][structure]
            ### plotting
            alphaStructures = .5
            axis.plot(mapX[mask_core], mapY[mask_core],
                      linestyle = 'none',
                      marker = 'o',
                      color = ssk_colors[j],
                      markeredgecolor = 'none',
                      markersize = 5.,
                      alpha = alphaStructures)

            ### record the masks
            mask_all = mask_all | mask_core

        ### plotting the distribution of all points in contours
        CS = axis.contour(xBinCent, yBinCent, height2D.T,
                          levels = contourLevels,
                          colors = 'k',
                          linewidths = [2., .5, .5, .5],
                          zorder = 1)
        #### labeling inline
        fmt = {}
        strs = ['95%', '75%', '50%', '25%']
        #strs = ['95', '75', '50', '25']
        for l, s in zip(CS.levels, strs):
            fmt[l] = s
        axis.clabel(CS, CS.levels,
                    inline = True,
                    inline_spacing = 1.5,
                    fmt = fmt,
                    fontsize = 12,
                    use_clabeltext = True)

        ### plotting the rest of the points
        axis.plot(mapX[~mask_all], mapY[~mask_all],
                  linestyle = 'none',
                  marker = 'o',
                  color = 'none',
                  markeredgecolor = 'k',
                  markersize = 3.,
                  alpha = .1,
                  zorder = 0)

        ### Plot the expected line widths
        axis.hlines([Sigma_SonicNT, Sigma_halfSonicNT], xmin, xmax,
                    linestyles = ['--', ':'],
                    colors = 'k')

        ### adjust the plot
        #### limits
        axis.set_xlim(xmin, xmax)
        axis.set_xscale(xscale)
        axis.set_ylim(ymin, ymax)
        #### ticks
        #axis.set_xticks([0., 5., 10., 15., 20.])
        axis.set_yticks([.2, .4, .6, .8])
        if i != 0:
            axis.set_yticklabels([])

        if xscale == 'log':
            axis.text(7., .85, reg,
                      weight = 'black',
                      horizontalalignment = 'right',
                      verticalalignment = 'top')
            axis.text(8., Sigma_SonicNT+.01, '$\sigma_{NT}=c_{s, ave}$',
                      size = 14.,
                      horizontalalignment = 'right',
                      verticalalignment = 'bottom')
            axis.text(8., Sigma_halfSonicNT+.01, '$\sigma_{NT}=0.5c_{s, ave}$',
                      size = 14.,
                      horizontalalignment = 'right',
                      verticalalignment = 'bottom')
        elif xscale == 'linear':
            axis.text(4.7, .85, reg,
                      weight = 'black',
                      horizontalalignment = 'right',
                      verticalalignment = 'top')
            axis.text(4.95, Sigma_SonicNT+.01, '$\sigma_{NT}=c_{s, ave}$',
                      size = 14.,
                      horizontalalignment = 'right',
                      verticalalignment = 'bottom')
            axis.text(4.95, Sigma_halfSonicNT+.01, '$\sigma_{NT}=0.5c_{s, ave}$',
                      size = 14.,
                      horizontalalignment = 'right',
                      verticalalignment = 'bottom')

        axis.yaxis.set_minor_locator(ticker.AutoMinorLocator(n = 4))
        axis.xaxis.set_minor_formatter(FuncFormatter2(ticks_format, interval = 2))
        axis.tick_params(axis='x', which='minor', labelsize=11)

        #### axis labels
        #axis.set_xlabel(r'Peak T$_{A^*}$ [K]',
        #                labelpad = -5.)
        #if i == 0:
        #    axis.set_ylabel('Line Width [km s$^{-1}$]')

    fig.text(.5, .035, '$T_{peak}$ [main-beam; K]',
             weight = 'bold',
             family = 'StixGeneral',
             horizontalalignment = 'center',
             verticalalignment = 'center')
    fig.text(.03, .5, '$\sigma_{{NH}_3}$ [km s$^{-1}$]',
             rotation = 90,
             weight = 'bold',
             horizontalalignment = 'center',
             verticalalignment = 'center')

    import styles

    return fig

def plotSigmas(list_dictionaries, plotSigma = 'sigma', plotRfromA = False):



    ###
    rcParams['figure.subplot.left'] = .09
    rcParams['figure.subplot.right'] = .97
    rcParams['figure.subplot.bottom'] = .12
    rcParams['figure.subplot.top'] = .96
    rcParams['font.size'] = 30
    ###
    dict_data, dict_masks, dict_YSOs, dict_Vlsr_predicted = list_dictionaries

    #ncols, nrows = 5, 4
    fig = plt.figure(figsize = (18., 18.))

    #Dmax = .145
    #Dmax = 1.3*0.10796408847 ####
    Dmax = .13
    rmin, rmax = 0., Dmax
    rbins = np.linspace(rmin, rmax, 12)  ####
    ymin, ymax = .02, .65

    #
    SigmaNT_Sonic = (np.sqrt(c.k_B*10.*u.K/(2.37*u.u))).to(u.km/u.s).value
    SigmaNT_halfSonic = (np.sqrt(.5**2.*c.k_B*10.*u.K/(2.37*u.u))).to(u.km/u.s).value
    #
    Sigma_SonicNT = (np.sqrt(c.k_B*10.*u.K/(17.031*u.u)+c.k_B*10.*u.K/(2.37*u.u))).to(u.km/u.s).value
    Sigma_halfSonicNT = (np.sqrt((c.k_B*10.*u.K/(17.031*u.u)+.5**2.*c.k_B*10.*u.K/(2.37*u.u)))).to(u.km/u.s).value

    #
    figLeft, figRight, figBottom, figTop = .06, .99, .055, .9
    gapHorizontal, gapVertical = .005, .005
    gapReg = .055 ## horizontal
    gapExtra = .055 ## Vertical
    frameWidth = (figRight-figLeft-gapReg-3.*gapHorizontal)/5.
    frameHeight = (figTop-figBottom-4.*gapVertical)/5.
    dict_frames = {'L1688': [], 'B18': []}
    ## L1688
    for i in range(12):
        irow, icol = i//3, i%3
        frame = fig.add_axes([figLeft+gapReg+gapHorizontal+2.*frameWidth+icol*(frameWidth+gapHorizontal),
                              figBottom+gapExtra+frameHeight+(4-irow-1)*(frameHeight+gapVertical),
                              frameWidth, frameHeight])
        dict_frames['L1688'].append(frame)
    frame = fig.add_axes([figLeft+gapReg+gapHorizontal+2.*frameWidth,
                          figBottom,
                          frameWidth, frameHeight])
    dict_frames['L1688'].append(frame)
    ## B18
    for i in range(6):
        irow, icol = i//2, i%2
        frame = fig.add_axes([figLeft+icol*(frameWidth+gapHorizontal),
                              figBottom+gapExtra+frameHeight+(4-irow-1)*(frameHeight+gapVertical),
                              frameWidth, frameHeight])
        dict_frames['B18'].append(frame)

    #
    lineSpacing = .026

    for i in range(19):

        # list of structures
        if i < 13:
            reg = 'L1688'
            listStructures = list(range(1, 13))+['extra']
            structure = listStructures[i]
            axis = dict_frames[reg][i]
            #j = i
        else:
            reg = 'B18'
            listStructures = range(1, 7)
            structure = listStructures[i-13]
            axis = dict_frames[reg][i-13]
            #j = i-13


        hdr = dict_data[reg]['header_GAS']
        mapNT = dict_data[reg]['SigmaNT']
        mapT = dict_data[reg]['SigmaT']
        mapSigma = dict_data[reg]['Sigma']
        mask = dict_masks[reg][structure]
        #maskFinite = np.isfinite(mapNT)&np.isfinite(mapT)
        distance = distances[reg]

        # deriving the profile for pixels within Dmax
        ## centroid
        meshx, meshy = np.meshgrid(np.arange(mask.shape[1]), np.arange(mask.shape[0]))
        cenx, ceny = np.mean(meshx[mask]), np.mean(meshy[mask])
        ## distance
        meshrPix = np.hypot(meshx-cenx, meshy - ceny)
        meshr = (meshrPix*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value
        ## effective radius
        stat = statBasic2D(mask.astype(float)[mask], (meshy[mask], meshx[mask]))
        stat.calculate()
        Reff = (stat.radius.value*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value
        Reff *= (2.*np.sqrt(2.*np.log(2.))) #FWHM as in Goodman+ 93
        Reff2 = (np.sqrt(stat.area_exact.value/np.pi)*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value
        Reff_low = np.arange(1., 30.)[np.array([np.sum(meshrPix<r)!=np.sum(mask&(meshrPix<r))
                                                for r in np.arange(1., 30.)])][0]
        Reff_low = (Reff_low*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value
        Reff_low = min(Reff_low, Reff-(1.*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value)
        Reff_high = np.arange(1., 30.)[np.array([np.sum(mask&(meshrPix<r))==np.sum(mask)
                                                 for r in np.arange(1., 30.)])][0]
        Reff_high = (Reff_high*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value
        Reff_high = max(Reff_low, Reff+(1.*np.radians(abs(hdr['CDELT1']))*distance).to(u.pc).value)

        if plotSigma == 'components':
            ## points within the mask
            mask_in = mask
            ptr_in, ptNT_in, ptT_in = meshr[mask_in], mapNT[mask_in], mapT[mask_in]
            ## all points within Dmax
            mask_all = (meshr < Dmax)
            ptr_all, ptNT_all, ptT_all = meshr[mask_all], mapNT[mask_all], mapT[mask_all]
            ### bin the points
            binCent_all = rbins[:-1] + .5*np.diff(rbins)
            binNT_all = np.array([np.nanmedian(ptNT_all[(ptr_all >= rbins[k])&(ptr_all < rbins[k+1])])\
                                  for k in range(len(rbins)-1)])
            binNTstd_all = np.array([np.nanstd(ptNT_all[(ptr_all >= rbins[k])&(ptr_all < rbins[k+1])])\
                                     for k in range(len(rbins)-1)])
            binT_all = np.array([np.nanmedian(ptT_all[(ptr_all >= rbins[k])&(ptr_all < rbins[k+1])])\
                                 for k in range(len(rbins)-1)])
            binTstd_all = np.array([np.nanstd(ptT_all[(ptr_all >= rbins[k])&(ptr_all < rbins[k+1])])\
                                    for k in range(len(rbins)-1)])

            # plotting
            # plot
            #axis = ax[i//ncols, i%ncols]
            ## points within the mask
            ### property A
            axis.plot(ptr_in, ptNT_in,
                      linestyle = 'none',
                      linewidth = 0.,
                      marker = 'o',
                      color = ssk_colors[2],
                      markeredgecolor = 'none',
                      markersize = 4.)
            ### property B
            axis.plot(ptr_in, ptT_in,
                      linestyle = 'none',
                      linewidth = 0.,
                      marker = 'o',
                      color = ssk_colors[0],
                      markeredgecolor = 'none',
                      markersize = 4.)
            ### bins of all
            axis.fill_between(binCent_all, binNT_all-.5*binNTstd_all, binNT_all+.5*binNTstd_all,
                              color = ssk_colors[2],
                              linewidth = 0.,
                              edgecolor = 'none',
                              alpha = .15)

            axis.fill_between(binCent_all, binT_all-.5*binTstd_all, binT_all+.5*binTstd_all,
                              color = ssk_colors[0],
                              linewidth = 0.,
                              edgecolor = 'none',
                              alpha = .15)
        elif plotSigma == 'sigma':
            ## points within the mask
            mask_in = mask
            ptr_in, ptSigma_in = meshr[mask_in], mapSigma[mask_in]
            ## all points within Dmax
            mask_all = (meshr < Dmax)
            ptr_all, ptSigma_all = meshr[mask_all], mapSigma[mask_all]
            ### bin the points
            binCent_all = rbins[:-1] + .5*np.diff(rbins)
            binSigma_all = np.array([np.nanmedian(ptSigma_all[(ptr_all >= rbins[k])&(ptr_all < rbins[k+1])])\
                                     for k in range(len(rbins)-1)])
            binSigmastd_all = np.array([np.nanstd(ptSigma_all[(ptr_all >= rbins[k])&(ptr_all < rbins[k+1])])\
                                        for k in range(len(rbins)-1)])

            # plotting
            # plot
            #axis = ax[i//ncols, i%ncols]
            ## points within the mask
            ### property A
            axis.plot(ptr_in, ptSigma_in,
                      linestyle = 'none',
                      linewidth = 0.,
                      marker = 'o',
                      color = ssk_colors[1],
                      markeredgecolor = 'none',
                      markersize = 4.)
            ### bins of all
            axis.fill_between(binCent_all, binSigma_all-.5*binSigmastd_all, binSigma_all+.5*binSigmastd_all,
                              color = ssk_colors[1],
                              linewidth = 0.,
                              edgecolor = 'none',
                              alpha = .15)



        ## effective radius used to calculate the physical properties
        '''
        axis.fill_between([.045, .055], ymin, ymax,
                           edgecolor = 'k',
                           facecolor = 'none',
                           alpha = .55,
                           hatch = '//////',
                           linewidth = .01)
        '''

        axis.fill_between([Reff_low, Reff_high], ymin, ymax,
                           facecolor = 'k',
                           edgecolor = 'none',
                           alpha = .1)
        axis.vlines(Reff, ymin, ymax,
                    linestyle = '-',
                    color = 'k')
        if plotRfromA:
            axis.vlines(Reff2, ymin, ymax,
                        linestyle = '--',
                        color = 'k')


        ### Plot the expected line widths
        if plotSigma == 'components':
            axis.hlines([SigmaNT_Sonic, SigmaNT_halfSonic], rmin, rmax,
                        linestyles = ['--', ':'],
                        colors = 'k')
            axis.text(.125, SigmaNT_Sonic+.02, '$c_{s, ave}$',
                      size = 20,
                      horizontalalignment = 'right')
            axis.text(.125, SigmaNT_halfSonic+.02, '$0.5c_{s, ave}$',
                      size = 20,
                      horizontalalignment = 'right')
        elif plotSigma == 'sigma':
            axis.hlines([Sigma_SonicNT, Sigma_halfSonicNT], rmin, rmax,
                        linestyles = ['--', ':'],
                        colors = 'k')
            axis.text(.125, Sigma_SonicNT+.02, '$\sigma_{NT}=c_{s, ave}$',
                      size = 20,
                      horizontalalignment = 'right')
            axis.text(.125, Sigma_halfSonicNT+.02, '$\sigma_{NT}=0.5c_{s, ave}$',
                      size = 20,
                      horizontalalignment = 'right')

        # annotation
        axis.text(.12, .6, str(structure).capitalize(),
                  size = 26,
                  horizontalalignment = 'right',
                  verticalalignment = 'top')

        ## adjust the subplot
        ### tick labels
        if i not in [0, 3, 6, 9, 12, 13, 15, 17]:
            axis.set_yticks([.3, .6])
            axis.set_yticklabels([])
        else:
            axis.set_yticks([.3, .6])
            axis.tick_params(axis='y', which='major', labelsize=24)
            axis.yaxis.labelpad = -1.5

        if i not in [9, 10, 11, 12, 17, 18]:
            axis.set_xticks([.05, .1])
            axis.set_xticklabels([])
        else:
            axis.set_xticks([.05, .1])
            labels = axis.get_xticklabels()
            axis.tick_params(axis='x', which='major', labelsize=24)
            #axis.xaxis.labelpad = -5.
            #plt.setp(labels, rotation=330)

        if i in [9, 12, 17]:
            axis.set_xlabel('$R_{eff}$ [pc]')
            axis.set_ylabel('$\sigma_{NT}$/$\sigma_{T}$ [km s$^{-1}$]')
        ### axis labels
        #if (i//ncols == (nrows-1)) and (i%ncols == 0):
        #    axis.set_xlabel('Distance [pc]')
        #    axis.set_ylabel(r'$\sigma$ [km s$^{-1}$]')
        ### limits
        axis.set_xlim(rmin, rmax)
        axis.set_ylim(ymin, ymax)

    fig.text(.075, .97, 'B18',
             size = 34,
             weight = 'black',
             horizontalalignment = 'left',
             verticalalignment = 'center')
    fig.text(figLeft+gapHorizontal+gapReg+2.*frameWidth+.015, .97, 'L1688',
             size = 34,
             weight = 'black',
             horizontalalignment = 'left',
             verticalalignment = 'center')

    # legend
    if plotSigma == 'components':
        #
        fig.text(figLeft+gapHorizontal, figBottom+4.*lineSpacing,
                 '$\sigma_{T}$ of Pixels Inside',
                 size = 30,
                 color = ssk_colors[0],
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')
        fig.text(figLeft+gapHorizontal, figBottom+3.*lineSpacing,
                 '$\sigma_{T}$ of All Pixels (binned; 1-$\sigma$)',
                 size = 30,
                 color = ssk_colors[0],
                 alpha = .45,
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')
        fig.text(figLeft+gapHorizontal, figBottom+2.*lineSpacing,
                 '$\sigma_{NT}$ of Pixels Inside',
                 size = 30,
                 color = ssk_colors[2],
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')
        fig.text(figLeft+gapHorizontal, figBottom+1.*lineSpacing,
                 '$\sigma_{NT}$ of All Pixels (binned; 1-$\sigma$)',
                 size = 30,
                 color = ssk_colors[2],
                 alpha = .45,
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')
        fig.text(figLeft+gapHorizontal, figBottom,
                 '$R_{eff}$',
                 size = 30,
                 color = 'k',
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')

        axfig = fig.add_axes([0., 0., 1., 1.])
        axfig.plot([figLeft+gapHorizontal+2.*frameWidth, figRight],
                   [figBottom+frameHeight+1.5*gapVertical, figBottom+frameHeight+1.5*gapVertical],
                   'k-')
        axfig.fill_between([figLeft, figLeft+2.*frameWidth],
                           figBottom-.2*lineSpacing,
                           figBottom+5.2*lineSpacing,
                           color = 'gray',
                           linewidth = 0.,
                           alpha = .15,
                           zorder = 999)
    elif plotSigma == 'sigma':

        #
        fig.text(figLeft+gapHorizontal, figBottom+2.*lineSpacing,
                 '$\sigma_{{NH}_3}$ of Pixels Inside',
                 size = 30,
                 color = ssk_colors[1],
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')
        fig.text(figLeft+gapHorizontal, figBottom+1.*lineSpacing,
                 '$\sigma_{{NH}_3}$ of All Pixels (binned; 1-$\sigma$)',
                 size = 30,
                 color = ssk_colors[1],
                 alpha = .45,
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')
        fig.text(figLeft+gapHorizontal, figBottom,
                 '$R_{eff}$',
                 size = 30,
                 color = 'k',
                 horizontalalignment = 'left',
                 verticalalignment = 'bottom')

        axfig = fig.add_axes([0., 0., 1., 1.])
        axfig.plot([figLeft+gapHorizontal+2.*frameWidth, figRight],
                   [figBottom+frameHeight+1.5*gapVertical, figBottom+frameHeight+1.5*gapVertical],
                   'k-')
        axfig.fill_between([figLeft, figLeft+2.*frameWidth],
                           figBottom-.2*lineSpacing,
                           figBottom+3.2*lineSpacing,
                           color = 'gray',
                           linewidth = 0.,
                           alpha = .15,
                           zorder = 999)

    axfig.set_xlim(0., 1.)
    axfig.set_ylim(0., 1.)
    axfig.set_xticks([])
    axfig.set_yticks([])
    for sp in ['left', 'right', 'bottom', 'top']:
        axfig.spines[sp].set_visible(False)

    import styles

    return fig
