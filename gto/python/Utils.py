################################################################################
#                                                                              #
# Copyright (C) 2024, 2024 Sylwester Wysocki (sw143@wp.pl)                     #
#                                                                              #
# This program is free software: you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# This program is distributed in the hope that it will be useful,              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program. If not, see <http://www.gnu.org/licenses/>          #
#                                                                              #
################################################################################

from math   import gamma
from mpmath import gammainc

# ------------------------------------------------------------------------------
#                               Public API
# ------------------------------------------------------------------------------

def DistSquared(r1, r2):
  distX = r1[0] - r2[0]
  distY = r1[1] - r2[1]
  distZ = r1[2] - r2[2]
  return distX**2 + distY**2 + distZ**2

def Dist(r1, r2):
  return sqrt(DistSquared(r1, r2))

# Calculate approximation of Boys function of order n at point x:
#
#         /1  2n         2
# Fn(x) = |  t  * exp(-xt ) dt
#         /0

def Boys(order, x):
  rv = 0.0
  if x > 0.0:
    g  = gamma(order + 0.5)
    gi = 1.0 - gammainc(order + 0.5, x, regularized=True)
    rv = 0.5 * g * gi * x**(-order - 0.5)
  else:
    rv = 1.0 / (2 * order + 1)
  return rv