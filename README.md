

* Are there existing, quantitative metrics for testing discernability
  of colormaps? Intuitively, I want to score a generated colormap for
  (A) how well it separates the given extrema while (B) remaining
  continuous and (C) evenly spreading out the change (derivative is
  near constant).
* Extend to using CMYK and possibly other, higher dimensional,
  representations of colors, to better handle maps of 4+ colors.
* Could extend with a composite mapper, which combines results from
  several colormaps?
  Alternatively, could deconstruct input params into non-overlapping
  triangles and delegate to TernaryColorMap
* Support RGBA values
* Support mapping arrays of points in one call