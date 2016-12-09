import itertools as itools
import math

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as la
from scipy.spatial import Delaunay

# Tolerance
__SLACK__ = 1e-08

class TriangleBCS(object):
    INSIDE = 0
    VERTEX = 1
    EDGE = 2
    OUTSIDE = 3

    def __init__(self, coords, strict=True):
        if strict:
            assert len(coords)==3, "Must provide 3 Cartesian coordinates for Barycentric conversions."
            assert not collinear(*coords), "Cartesian coordinates must not be collinear."
        self._strict = strict
        self._a = homogenize(coords).T

    def to(self, point):
        # Barycentric coordinate approach
        # (https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Interpolation_on_a_triangular_unstructured_grid)
        b = homogenize([point]).T
        if self._strict:
            bary_coords = la.solve(self._a,b)
        else:
            bary_coords = la.lstsq(self._a,b)[0]

        location = None
        # if any coordinate is negative, point is outside triangle.
        if any(bary_coords < -__SLACK__):
            location = self.OUTSIDE
        # if any coordinate is 1 at this point, means exactly one
        # coord is 1 and must be on vertex
        elif (any(bary_coords > 1-__SLACK__)):
            location = self.VERTEX
        # if any coordinate is near 0 at this point, means on edge.
        elif any(bary_coords < __SLACK__):
            location = self.EDGE
        # otherwise, must be inside
        else:
            location = self.INSIDE

        return [bary_coords, location]

def project(a, b):
    """Project vector b onto line given by vector a"""
    a_ = np.array(a)
    b_ = np.array(b)
    return np.dot(a_.T,b_)/np.dot(a_.T,a_)*a_

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

def collinear(p1,p2,p3):
    x1,y1 = p1
    x2,y2 = p2
    x3,y3 = p3
    m = np.array([[x2-x1, y2-y1], [x3-x1, y3-y1]])
    return abs(la.det(m)) < __SLACK__

def plot_cmap(cmap, show=True, ax=plt, scale=4, bounds=None, constrain=False,
              *args):
    constrained = constrain and bounds is not None
    if constrained:
        hull = Delaunay(bounds[:])
        test = lambda x,y: hull.find_simplex([[x,y]])>=0
    else:
        test = lambda x,y: True

    colors = eval_cmap(cmap=cmap, filter=test, scale=scale, *args)
    ax.imshow(colors, origin='lower')

    if not constrained:
        plot_bounds(bounds ,ax, scale=scale)
    if show:
        plt.show()

def eval_cmap(cmap, vmax=100, filter=lambda x,y: x+y <= 100, scale=1):
    cmap.verbose = False
    n = vmax//scale + 1
    xs = np.linspace(0, vmax, n)
    xv, yv = np.meshgrid(xs, xs)
    colors = np.zeros([n,n]).tolist()
    for i in range(n):
        for j in range(n):
            x,y = (xv[j,i], yv[j,i])
            if filter(x,y):
                res = cmap((x,y))
                if res is None or np.any(np.isnan(res)):
                    res = (0,0,0)
                if np.any(np.isnan(res)):
                    print((x,y), res)
                    colors[j][i] = res
            else:
                colors[j][i] = (0,0,0)
    return colors

def plot_bounds(corners, ax=None, scale=1):
    n = len(corners)
    for i in range(n):
        x1,y1 = corners[i]
        x2,y2 = corners[(i+1)%n]
        if ax is None:
            ax = plt.gca()
        ax.plot([x1/scale,x2/scale], [y1/scale,y2/scale], 'k', linewidth=1)

def plot_color_triangles(colors):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    n = len(colors)
    for indices in itools.combinations(range(n), 3):
        local_coords = [colors[i] for i in indices]
        c = np.average(local_coords, axis=1)
        for i in range(3):
            x1,y1,z1 = local_coords[i]
            x2,y2,z2 = local_coords[(i+1)%3]
            ax.plot([x1,x2], [y1,y2], [z1, z2], color=tuple(c.tolist()))
    for color in colors:
        ax.scatter(color[0],color[1], color[2], color=color, s=20)
    plt.show()

# generate corners of regular n-poly centered at 0,0
def poly_corners(n, scale=1):
    step = 2*math.pi/n
    return np.array([[math.cos(step*i)*scale, math.sin(step*i)*scale] for i in range(n)])

# Convert points to homegenous coordinates for the given dimension.
# Homogeneous coords in 2D space have a 3rd coordinate, etc.
# Note, that input points may be of lower dimensionality than desired
# (e.g. 2D points provided when 3D homo coords desired). Missing
# coordinates will be zeroed out before adding homogeneous coord=1.
def homogenize(points, dim=2):
    n = len(points)
    res = np.array(points)
    in_dim = len(points[0])
    if in_dim < dim:
        res = np.concatenate((res, np.zeros((n,dim-in_dim))), axis=1)
    return np.concatenate((res,np.ones((n,1))), axis=1)

def make_tests(radius):
    return [
        [[[radius*x+50 for x in p] for p in poly_corners(3)], [(0,0,0), (0,0,1), (1,0,0)]],
        [[[radius*x+50 for x in p] for p in poly_corners(4)], [(0,0,0), (0,0,1), (1,1,1), (1,0,0)]],
        [[[radius*x+50 for x in p] for p in poly_corners(5)], [(0,0,0), (0,0,1), (1,0,0), (0,1,0), (1,1,1)]],
        [[[radius*x+50 for x in p] for p in poly_corners(6)], [(0,0,1), (0,1,1), (0,1,0), (1,1,0), (1,0,0), (1,0,1)]]
    ]

def run_tests(classes, shapes, scale=4, constrain=True):
    if type(classes) != list:
        classes = [classes]
    y = len(classes)
    x = len(shapes)
    fig, axes = plt.subplots(y,x, figsize=(x*2, y*2),
                             subplot_kw={'xticks': [], 'yticks': []})
    if y==1:
        axes = [axes]
    if x==1:
        axes = [[a] for a in axes]
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    for i,cls in enumerate(classes):
        for j,shape in enumerate(shapes):
            print("plotting", cls.__name__, 'shape', j)
            coords, vert_colors = shape[:]
            cmap = cls(*shape)
            ax = axes[i][j]

            plot_cmap(cmap, ax=ax, scale=scale, constrain=constrain, bounds=coords,show=False)

            for coord, color in zip(coords, vert_colors):
                ax.scatter(coord[0]/scale, coord[1]/scale, color=color, edgecolors=(.5,.5,.5,.5), s=np.pi*5**2)
                ax.set_xlim([0, 100/scale])
                ax.set_ylim([0, 100/scale])
            if j==0:
                ax.set_ylabel(cls.__name__)

    plt.show()
