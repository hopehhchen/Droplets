import sys
import warnings

#
import numpy as np

#
import astropy.units as u

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
