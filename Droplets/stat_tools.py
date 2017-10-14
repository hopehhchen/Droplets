import sys
import warnings

#
import numpy as np

#
import astropy.units as u
import astropy.modeling as modeling

#
from constants import *
from ssk_colors import *
import styles


#### the codes are directly modified from astrodendro (PPStatistics and PPVStatistics)
def mom2_along(selfmom2, direction):
    w = np.atleast_2d(direction).astype(np.float)
    for row in w:
        row /= np.linalg.norm(row)

    result = np.dot(np.dot(w, selfmom2), w.T)
    if result.size == 1:
        result = np.asscalar(result)

    return result

def projected_paxes(selfmom2, axes):
    axes = tuple(axes)
    mom2 = mom2_along(selfmom2, axes)
    w, v = np.linalg.eig(mom2)
    order = np.argsort(w)

    return tuple(v[:, o] for o in order[::-1])

## basic object 2D
class statBasic2D(object):

    def __init__(self, values, indices,
                 metadata = {'spatial_scale': None}):


        self.values = values
        self.indices = indices
        self.spatial_scale = metadata['spatial_scale']

        # mom0/mom1
        self.mom0 = np.nansum(self.values)
        self.mom1 = [np.nansum(i * self.values) / self.mom0 for i in self.indices]

        # mom2 (covariance matrix)
        v = self.values / self.mom0
        nd = len(self.indices)
        zyx = tuple(i - m for i, m in zip(self.indices, self.mom1))
        result = np.zeros((nd, nd))

        for i in range(nd):
            result[i, i] = np.nansum(v * zyx[i] ** 2)
            for j in range(i + 1, nd):
                result[i, j] = result[j, i] = np.nansum(v * zyx[i] * zyx[j])
        self.mom2 = result

        # principal axes
        w, v = np.linalg.eig(self.mom2)
        order = np.argsort(w)
        self.paxes = tuple(v[:, o] for o in order[::-1])

    def calculate(self):
        ## _sky_paxes
        ax = [(1, 0), (0, 1)]
        a, b = projected_paxes(self.mom2, tuple(ax))
        a = list(a) ## was in numpy.array
        b = list(b)
        self._sky_paxes = tuple(a), tuple(b)

        ## major_sigma/minor_sigma/radius
        dx = self.spatial_scale or u.pixel
        a, b = self._sky_paxes
        # We need to multiply the second moment by two to get the major axis
        # rather than the half-major axis.
        self.major_sigma = dx * np.sqrt(mom2_along(self.mom2, tuple(a)))
        self.minor_sigma = dx * np.sqrt(mom2_along(self.mom2, tuple(b)))
        self.radius = self.major_sigma.unit * np.sqrt(self.major_sigma.value * self.minor_sigma.value)
        ## position_angle
        a = list(a)
        self.position_angle = np.degrees(np.arctan2(a[0], a[1])) * u.degree

        ## area_exact
        #dx = self.spatial_scale or u.pixel ## assigned above
        indices = zip(*tuple(self.indices[i] for i in range(2)))
        self.area_exact = len(set(indices)) * dx ** 2


class statGrad(object):

    def __init__(self, dict_data, dict_masks, region, structureID):

        self.Vlsr = dict_data[region]['Vlsr']
        self.eVlsr = dict_data[region]['eVlsr']
        self.mask = dict_masks[region][structureID]

        xmesh, ymesh = np.meshgrid(np.arange(self.Vlsr.shape[1]),
                                   np.arange(self.Vlsr.shape[0]))

        z = self.Vlsr[self.mask & np.isfinite(self.Vlsr)]
        x, y = xmesh[self.mask & np.isfinite(self.Vlsr)],\
               ymesh[self.mask & np.isfinite(self.Vlsr)]
        w = (1./self.eVlsr**2.)[self.mask & np.isfinite(self.Vlsr)]

        # Fit the data using astropy.modeling
        model = modeling.polynomial.Polynomial2D(degree = 1)
        fitter = modeling.fitting.LevMarLSQFitter()
        fitted = fitter(model, x, y, z,
                        weights = w)
        ## fitted parameters
        x1 = fitted.parameters[1]  ## x
        x2 = fitted.parameters[2]  ## y
        e1 = np.sqrt((fitter.fit_info['param_cov'])[1, 1])
        e2 = np.sqrt((fitter.fit_info['param_cov'])[2, 2])
        ## Store the fits
        self._fit = {'model': model,
                     'fitted': fitted,
                     'fitter': fitter,
                     'x': x1,
                     'y': x2,
                     'ex': e1,
                     'ey': e2}


        # result
        ## Grad in km/s/pix
        Grad = np.sqrt(x1**2.+x2**2.)
        eGrad = (.5*1./np.sqrt(x1**2.+x2**2.))**2.\
                *((2.*x1*e1)**2.+(2.*x2*e2)**2.)
        eGrad = np.sqrt(eGrad)
        self.Grad, self.eGrad = Grad, eGrad
        ## PAGrad in degrees
        PAGrad = np.degrees(np.arctan2(x2, x1))
        if (PAGrad >= -180.) and (PAGrad <= -90.): ## convert to E of N
            PAGrad = PAGrad + 270.
        else:
            PAGrad = PAGrad - 90.
        ePAGrad = (1./(1.+(x1/x2)**2.)*(e1/x2))**2.\
                  +(1./(1.+(x1/x2)**2.)*(x1*e2/x2**2.))**2.
        ePAGrad = np.sqrt(ePAGrad)
        self.PAGrad, self.ePAGrad = PAGrad, ePAGrad
