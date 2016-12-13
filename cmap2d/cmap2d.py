import abc
import itertools as itools

import numpy as np
from numpy import linalg as la

from . import util
from .simplex import simplex_coordinates1
from .util import __SLACK__

class ColorMapBase(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, coords, colors, verbose=False, **kargs):
        super().__init__(**kargs)
        self._coords = coords
        self._colors = colors
        self.verbose = verbose

        self.set_oob((0.0, 0.0, 0.0))

    def _print(self, *args):
        if self.verbose:
            print(*args)

    def set_oob(self, color):
        # If out of range, default to this color.
        self._oob_color = color

    @abc.abstractmethod
    def __call__(self, point):
        pass

class CompositeColorMapBase(ColorMapBase):
    def __init__(self, *args):
        super().__init__(*args)

        self._gen_sub_bcs()

    def _gen_sub_bcs(self):
        self._bcs = []
        # generate a ternary cmap for each triple of coords
        for indices in itools.combinations(range(len(self._coords)), 3):
            tri_coords = [self._coords[i] for i in indices]
            if util.collinear(*tri_coords):
                continue
            tri_colors = [self._colors[i] for i in indices]
            self._bcs.append([util.TriangleBCS(tri_coords), tri_colors])


class SimplexColorMap(ColorMapBase):
    # Simple idea: map input coords + colors onto simplex with same number of vertices
    # This mapping is obtained from the least squares projection between
    # the two in homogeneous coordinates.
    # Then use same projection on any sample points, compute color based on
    # average of vertex colors (weighted by distance on simplex).

    def __init__(self, *args):
        super().__init__(*args)
        # The dimension of the space we are operating in.
        self._dim = len(self._coords)-1

        a = util.homogenize(self._coords,self._dim)
        self._shape_coords = b = simplex_coordinates1(self._dim).T
        self._x = la.lstsq(a,b)[0]
        #self._max_dist = la.norm(b[0]-b[1])
        self._max_dist = la.norm(util.project(-b[0], b[1]-b[0]))

    def __call__(self, point):
        a = util.homogenize([point],self._dim)
        nps = np.dot(a,self._x)
        #distances = [[max(max_dist-la.norm(p-c),0) for c in b] for p in nps]
        distances = [[max(self._max_dist-la.norm(util.project(-c, p-c)), 0) for c in self._shape_coords] for p in nps]
        norm_dist = [np.array(k)/sum(k) for k in distances]
        out = np.dot(norm_dist, self._colors).clip(min=0.0, max=1.0)
        return out.tolist()[0]

class AvgColorMap(ColorMapBase):
    """Simple class, maps input to color based on weighted average of inverse of
    distance to all color coordinates."""

    def __init__(self, *args):
        super().__init__(*args)
        self._max_dist = 10000

    def __call__(self, point):
        distances = la.norm(self._coords-np.array(point), axis=1, ord=2)
        distances = self._max_dist/distances
        if np.max(distances) < __SLACK__:
            return self._oob_color
        color = np.dot(distances, self._colors)
        if np.max(color) > 1:
            color /= np.max(color)
        np.clip(color, 0, 1, color)
        return color

class TernaryColorMap(ColorMapBase):
    # 2D colormap, as ternary plot
    # Barycentric coordinate approach (https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Interpolation_on_a_triangular_unstructured_grid)
    # aka treat it as a ternary plot https://en.wikipedia.org/wiki/Ternary_plot
    def __init__(self, *args):
        super().__init__(*args)
        self._bcs = util.TriangleBCS(self._coords, False)

    def __call__(self, point):
        bpoints, _ = self._bcs.to(point)
        color = np.dot(bpoints.T, self._colors)
        np.clip(color, 0, 1, color)
        return color[0].tolist()

class ProjectionColorMap(ColorMapBase):
    # Dead-simple - just compute projection matrix from input coords
    # to RGB color space.
    def __init__(self, *args):
        super().__init__(*args)
        self._x = la.lstsq(util.homogenize(self._coords),self._colors)[0]

    def __call__(self, point):
        color = np.dot(util.homogenize([point]),self._x)
        np.clip(color, 0, 1, color)
        return color[0].tolist()

class ColorMapCrude(CompositeColorMapBase):
    def __call__(self, point):
        outs = []
        weights = []
        for bcs,colors in self._bcs:
            bary_coords, location = bcs.to(point)
            self._print("barycentric coordinates:", location, "\n", bary_coords)
            if location == bcs.OUTSIDE:
                color = self._oob_color
            else:
                color = np.dot(bary_coords.T,colors)[0]

            if location <= bcs.VERTEX:
                outs.append(color)
                # INSIDE is obviously weight = 1
                # VERTEX means that all triangles will either have
                # this as a vertex or OUTSIDE, so weight doesn't matter.
                weights.append(1)
            elif location == bcs.EDGE:
                outs.append(color)
                # for each point on an edge of any triangle,
                # we expect n-2 other triangles to share that edge, and
                # produce the same value. To avoid double counting, we
                # give them all an appropriate weight.
                weights.append(1/(len(self._coords)-2))
        self._print("outs:", len(outs), outs)
        if len(outs) == 0:
            color = self._oob_color
        else:
            color = np.average(outs,axis=0, weights=weights)
        if color is None or np.any(np.isnan(color)):
            print("oob!", point, outs)
            color = self._oob_color
        if np.max(color) > 1+__SLACK__ or np.min(color) < -__SLACK__:
            print("oob2!", point, color)
            color = self._oob_color
        color = np.array(color)
        np.clip(color, 0, 1, color)
        color = color.tolist()
        return color

class ColorMap(CompositeColorMapBase):
    def __call__(self, point):
        outs = []
        weights = []
        for bcs,colors in self._bcs:
            bary_coords, location = bcs.to(point)
            # we'll let each cmap vote on points up to an additional
            # unit beyond it's boundaries.
            np.clip(bary_coords, -1, 2, out=bary_coords)
            self._print("barycentric coordinates:", location, "\n", bary_coords)

            outs.extend(colors)
            weights.extend(bary_coords.T[0,0:].tolist())
        # convert all weights to be symmetric about 1,
        # and weights@1=1, decreasing everywhere else
        weights = [1-abs(x-1) for x in weights]
        # TODO: should further normalize all to [0,1] range?

        if len(outs) == 0:
            color = self._oob_color
        elif sum(weights) == 0:
            # weighted average compuation complains, so just hardcode it.
            color = self._oob_color
        else:
            color = np.average(outs,axis=0, weights=weights)
#             color = np.dot(np.array([weights]),outs)[0]
            if np.max(color) > 1:
                color /= np.max(color)
            np.clip(color, 0, 1, color)

            color = color.tolist()
        self._print("outs:", point, color, len(outs), outs)
        if color is None or np.any(np.isnan(color)):
            print("oob!", point, outs)
            color = self._oob_color
        return color

class ColorMap2(ColorMapBase):
    def __call__(self, point):
        # Just treat every vertex equally, and compute the best
        # fit of weights to derive input vertex.
        self._print("point:", point)
        a = util.homogenize(self._coords).T
        b = util.homogenize([point]).T
        weights = la.lstsq(a,b)[0]
        weights = [1-abs(x-1) for x in weights.T[0,0:].tolist()]
        self._print("weights:", weights)
        if sum(weights) == 0:
            color = self._oob_color
        else:
            color = np.dot(np.array([weights]),self._colors)[0]
            m = 1
            if np.max(color) > 1:
                m = np.max(color)
            if np.min(color) < 1:
                m += -np.min(color)
            color /= m
            np.clip(color, 0, 1, color)
            color = color.tolist()
        if color is None or np.any(np.isnan(color)):
            print("oob!", point, outs)
            color = self._oob_color
        self._print("color:", color)
        return color
