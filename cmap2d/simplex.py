import numpy as np

def simplex_coordinates1 ( m ):

#*****************************************************************************80
#
## SIMPLEX_COORDINATES1 computes the Cartesian coordinates of simplex vertices.
#
#  Discussion:
#
#    The simplex will have its centroid at 0
#
#    The sum of the vertices will be zero.
#
#    The distance of each vertex from the origin will be 1.
#
#    The length of each edge will be constant.
#
#    The dot product of the vectors defining any two vertices will be - 1 / M.
#    This also means the angle subtended by the vectors from the origin
#    to any two distinct vertices will be arccos ( - 1 / M ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    28 June 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer M, the spatial dimension.
#
#    Output, real X(M,M+1), the coordinates of the vertices
#    of a simplex in M dimensions.
#

  x = np.zeros ( [ m, m + 1 ] )

  for k in range ( 0, m ):
#
#  Set X(K,K) so that sum ( X(1:K,K)^2 ) = 1.
#
    s = 0.0
    for i in range ( 0, k ):
      s = s + x[i,k] ** 2

    x[k,k] = np.sqrt ( 1.0 - s )
#
#  Set X(K,J) for J = K+1 to M+1 by using the fact that XK dot XJ = - 1 / M.
#
    for j in range ( k + 1, m + 1 ):
      s = 0.0
      for i in range ( 0, k ):
        s = s + x[i,k] * x[i,j]

      x[k,j] = ( - 1.0 / float ( m ) - s ) / x[k,k]

  return x
