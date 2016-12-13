# cmap2d: RGB Colormaps for 2D

This is a python package, cmap2d, which provides tools for generating colormaps
for 2 dimensions. In other words, while typical colormaps map from a single (1D)
value to a single color value, a 2D colormap maps from 2 values (2D) to a single
color value. These 2D colormaps are useful for visualizing 2 variables at once.

The package contains several different ColorMap implementations. A ColorMap is
initialized with a set of 2D points (`coords`) and associated color values
(`colors`). The `coords` define the domain of your input, while the `colors`
defines the range of the desired output colors. Then, when invoked with a sample
point from the 2D data, the initialized ColorMap will return a color
interpolated from the colorspace defined by `colors`. Each ColorMap class in the
package corresponds to a different approach for computing the interpolation.

The package also contains several utility functions for visualizing the
colormaps. These utilities use [matplotlib](http://matplotlib.org/).

You can see the package in use and get a better feel for the different
classes and how to use them [here](https://github.com/oderby/cmap2d/blob/master/demo.ipynb).

## Limitations and Notes

* Input coordinates should be the vertices of a convex polygon, defined in
  clockwise or counterclockwise order. I haven't tested providing random or
  un-ordered points at all, and do not expect it to work.
* Only tested with RGB values. I think the ColorMaps should work with any color
  representation, but the plotting utilities only seem to work with RGB
  (limitation in matplotlib?)
* Only tested with 2D input values. Conceptually, there's no reason these
  shouldn't work with higher-dimensional inputs, but I just haven't spent
  any time testing it.
* All ColorMaps 3+ colors to map. For exactly 3 colors, the mapping is unique,
  and so most ColorMaps will produce the same or nearly the same mapping. For
  4+, there are many possible mappings, and that is why there are several
  different ColorMaps. If you want to map 4+ colors, I recommend testing several
  maps and seeing which you like the most.
* One of the benefits of not assuming constraints on the input domain (many 1D
  colormaps map from the range [0,1) to colors) is being able to potentially
  gracefully handle inputs outside the specified range. Many of the ColorMaps
  included in cmap2d support this.
* Documentation is currently poor - best reference is the worked example in [this
  jupyter notebook](https://github.com/oderby/cmap2d/blob/master/demo.ipynb)
* The utility functions for plotting and comparing ColorMaps are somewhat
  brittle, and might not work well if your colormap isn't defined over at least
  [0,100] range in all dimensions.

# Meta Comment
This package is something I built quickly to satisfy a need as part of a
separate exploration. I had a problem (mapping US Presidential election
results), for which I was exploring better ways to visualize the outcomes. I
built this little package as an exercise, and have decided to open source it,
in case anyone else finds it useful.

In other words, this is definitely provided AS IS, and not something I plan on
spending much further time on. That said, please enjoy and do let me know if you
have any questions or feedback!

# Future Directions

* Are there existing, quantitative metrics for testing/scoring discernability
  of colormaps? Intuitively, I want to score a generated colormap for
  (A) how well it separates the given extrema while (B) remaining
  continuous and (C) evenly spreading out the change (derivative is
  near constant).
* Extend to using CMYK and possibly other, higher dimensional,
  representations of colors, to better handle maps of 4+ colors.
* Test with N dimensional input?
* Could extend with a composite mapper, which combines results from
  several colormaps?
  Alternatively, could deconstruct input params into non-overlapping
  triangles and delegate to TernaryColorMap
* Support other color values (RGBA, HSV, etc.)
* Support mapping arrays of points in one call
