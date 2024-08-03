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

import Utils

from math import fabs, exp, pi, sqrt

#-------------------------------------------------------------------------------
#                                 Constants
#-------------------------------------------------------------------------------

#                       v------- total sum of all componens
# q = [ xyz xyz xyz xyz s ]      and all functions
#       ^   ^    ^    ^
#       A    B    C    D
#       func func func func
#       tion tion tion tion
#
# Total length of q is 12 + 1 elements.
# Extra total sum is calculated internally only.
# Caller doesn't need to pass it to public API.

IDX_Q_SUM = 12

IDX_LEFT_A   = 18
IDX_RIGHT_A  = 19
IDX_LEFT_AA  = 20
IDX_RIGHT_AA = 21
IDX_MIXED    = 22

# We use hard-coded Obara-Saika scheme to reduce [ab|cd] to [a-1b|cd].
# If we need to reduce angular momentum at different position (b,c,d),
# then we use one of ERI symmetry to move reduced function back to a.
#
# Below ERIs with real functions are equivalent:
# [ab|cd] = [ba|cd] = [ba|dc] = [ab|dc] = [cd|ab] = [cd|ba] = [dc|ba] = [dc|ab]
#
# We have 2 symmetries:
# - swap functions within the same bra/ket side: [ab|->[ba| and |cd]->|dc],
# - swap the whole pair between bra/ket sides: [ab|cd]->[cd|ab]
#
# Permutate original [ab|cd] functions depending on selected one.
# Below map shows new functions order within integral.
MUTATE_MAP = [
# a, b, c, d    New [a,b,c,d] if...
  0, 1, 2, 3, # ... function no. 0 chosen (old a)
  1, 0, 2, 3, # ... function no. 1 chosen (old b)
  2, 3, 0, 1, # ... function no. 2 chosen (old c)
  3, 2, 0, 1  # ... function no. 3 chosen (old d)
# ^
# function to reduce
# always go to new
# a
]

#-------------------------------------------------------------------------------
#                    Inernal implementation (private)
#-------------------------------------------------------------------------------

#
# Find expansion coefficient for 2-electron repulsion (ERI) integral
# using Obara-Saika scheme.
#
# How does it work:
# -----------------
#   Ref. Obara, S.; Saika, A. J Chem Phys 1986, 84, 3963.
#   https://doi.org/10.1063/1.450106
#
#   At each call we reduce total angular momentum by 1.
#   Eg. [200 000 000 110] -> [100 000 000 110]
#
#   Each call generates up to 8 child calls.
#   We stop recursion if we reduce integral to basic [ss|ss].
#
# How to use results:
# -------------------
# Final integral is:
#   [ab|cd] = C0 * F(0,t) + C1*F1(1,t) + ... + Cn*Fn(n,t)
#
#   - Ci = coefficients calculated by this function (see rv[] param)
#   - F(i,t) = Boys function of i-th order
#
#   - t is the same for all orders (can be precomputed):
#
#        (za + zb) * (zc + zd)       2
#   t =  --------------------- * |PQ|
#        za + zb + zc + zd
#
#   |PQ| - distance between P and Q centers.
#    P   - bi-center on "left" side (AB),
#    Q   - bi-ecnter of "right" side (CD)
#
# Parameters:
# -----------
#   - rv[] - calculated coefficient grouped by Boys() function
#            order (OUT)
#
#   - q[]  - angular momentum of each functions
#            [LAx,LAy,LAy, LBx,LBy,LBz, LCx,LCy,LCz, LDx,LDy,LDz] (IN)
#
#   - currentTerm  - currently accumulated term passed between calls (IN)
#
#   - currentOrder - current accumulated order of Boys() function (IN)
#
#   - params - constants calcualted once, needed to calculate
#              coefficient (IN)
#
# Returns:
# --------
#   After finished rv[] contains coefficients needed to be multiplicated
#   by related Boys(s,t) functions.
#   See ObaraSaika_ERI() for more.
#
def _ObaraSaika_ERI_Expansion(rv, q, currentTerm, currentOrder, params):
  if q[IDX_Q_SUM] == 0:
    # Already reduced to [000 000 | 000 000] = [ss|ss]
    # Accumulate terms per Boys function order.
    rv[currentOrder] += currentTerm

  else:
    # Non-zero total angular momentum.
    # Find first non-zero index to reduce.
    # We can choose any path to get correct result.
    #
    # Example:
    # --------
    #  component: xyz xyz xyz xyz
    #  function:  000 111 222 333
    #             --- --- --- ---
    #  Angular:   000 010 200 002
    #  Momentum        ^
    #  (L)             we reduce here
    #                  function no. 1
    #                  component y (1)
    #                  L1.y (-1 here)
    #
    #   We reduce:
    #   [ 000 010 200 002 ] -> [ 000 000 200 002 ]
    #          ^                      ^
    #         -1 here                 -1 here
    #

    reduceIdx = 0
    while q[reduceIdx] == 0:
      reduceIdx += 1

    reduceFun = reduceIdx // 3 # Function number (0,1,2,3)
    reduceXyz = reduceIdx %  3 # Component of angular momentum (x,y,z)

    # Check does selected function belong to <bra| (left) or |ket> (right) side.
    # [01|23] = <bra| 1/r2 |ket>
    isKet = 1 if reduceFun > 1 else 0

    #
    # Fetch parameters for function being reduced (A,B,C,D) and
    # selected angular momentum component (i=x,y,z)
    #

    # Xi = Ai, Bi, Ci, Di
    # Possible improvement: Remove magic numbers?
    Xi = params[reduceXyz * 6 + reduceFun]

    # leftX/Y/Z or rightX/Y/Z
    # Possible improvement: Remove magic numbers?
    leftOrRightXyz = params[reduceXyz * 6 + isKet + 4]

    # Global - no coord (x,y,z), but depends on all functions together.
    # Related strictly to one side (bra/ket, left/right)
    leftOrRightA   = params[IDX_LEFT_A  + isKet]
    leftOrRightAA  = params[IDX_LEFT_AA + isKet]

    # Appears for terms with mixed bra-ket terms.
    mixed = params[IDX_MIXED]

    #
    # Permutate original [ab|cd] integral to keep function, which we're going
    # to reduce at "a" position.
    # Due to this we can apply the same hard-coded Obara-Saika scheme for
    # each integrals.
    # This trick works due to ERIs symmetries e.g. [ab|cd] = [ba|cd] etc.
    #

    a = MUTATE_MAP[reduceFun * 4 + 0] # <- function to reduce goes here
    b = MUTATE_MAP[reduceFun * 4 + 1]
    c = MUTATE_MAP[reduceFun * 4 + 2]
    d = MUTATE_MAP[reduceFun * 4 + 3]

    #
    # Helper function to make one child call with reduced integral.
    #

    def _goDeeper(factor, extraReduceFun, orderDelta):
      # Prepare new state for child node.
      newQ     = q[:]
      newTerm  = currentTerm  * factor
      newOrder = currentOrder + orderDelta

      # Reduce selected qi on each child.
      newQ[reduceIdx] -= 1
      newQ[IDX_Q_SUM] -= 1

      goOn = True

      # Extra reduction on child node if needed.
      if extraReduceFun != None:
        extraReduceIdx = extraReduceFun * 3 + reduceXyz

        if newQ[extraReduceIdx] == 0:
          # Already reduced to 0.
          # Nothing to do.
          goOn = False

        else:
          # Child angular momentum is still greater than 0.
          # Scale child term by reduced angular momentum.
          newTerm *= newQ[extraReduceIdx]

          # Reduce angular momentum by one and go on deeper.
          newQ[extraReduceIdx] -= 1
          newQ[IDX_Q_SUM]      -= 1

      # Go on deeper recursively if new integral exists (non-negative
      # angular momentum).
      if goOn:
        _ObaraSaika_ERI_Expansion(rv, newQ, newTerm, newOrder, params)

    # Apply Obara-Saika scheme.
    # Each reduction generate up to 8 child nodes.
    # Note, that some of these functions doesn't make sense,
    # becase after subtraction we may get negative angular momentum.
    # These "fake" functions will be skipped.

    # Direct implementation of eq.39 from:
    #   Obara, S.; Saika, A. J Chem Phys 1986, 84, 3963.
    #   https://doi.org/10.1063/1.450106
    _goDeeper(Xi             , None ,  0)
    _goDeeper(leftOrRightXyz , None , +1)
    _goDeeper(leftOrRightA   ,    a ,  0)
    _goDeeper(leftOrRightAA  ,    a , +1)
    _goDeeper(leftOrRightA   ,    b ,  0)
    _goDeeper(leftOrRightAA  ,    b , +1)
    _goDeeper(mixed          ,    c , +1)
    _goDeeper(mixed          ,    d , +1)

#-------------------------------------------------------------------------------
#                              Public API
#-------------------------------------------------------------------------------

#
# Calculate 2-electron repulsion (ERI) integral using Obara-Saika scheme.
#           //                 1
# [ab|cd] = || a(r1)*b(r1) * ----- * c(r2) * d(r2) dr1 dr2
#           //                r12
#
# r1 = [x1, y1, z1] = position of 1-st electron
# r2 = [x2, y2, z2] = position of 2-nd electron
#
# a,b,c,d - primitive Gaussian functions:
#             Lx        Ly        Lz                   2
#   a = (x-Ax)   * (y-Ay)  * (z-Az)  * exp(zetaA * |r-A| )
#   etc.
#
# Parameters:
# -----------
#   za,zb,zc,zd - exponential coefficient of A,B,C,D functions ("zeta") (IN)
#   ra,rb,rc,rd - centers of A,B,C,D functions (IN)
#   angularMomentum - angular momentum of each functions (IN)
#
# Angular momentum contains xyz components for each function.
# Each angular momentum has 3 components: lx, ly, lz.
#
# Example:
# --------
#   za = zb = zc = zd = 1
#   ra = rb = rc = rd = [0,0,0]
#   angularMomentum   = [1,0,0, 0,0,0, 0,1,1, 2,0,0]
#
#   All GTOs are centered on (0,0,0).
#   All GTOs have zeta (exponent coefficient) set to 1.
#   We calcualte integral [px s | dyz dxx] = [100 000 | 011 200].
#
# Returns:
# --------
#   Value of [ab|cd] integral (ERI).
#
def ObaraSaika_ERI(za, zb, zc, zd, ra, rb, rc, rd, angularMomentum):

  if len(angularMomentum) == IDX_Q_SUM:
    # Precalculate sum of all elements to detect [ss|ss] integrals easily.
    # We need to update this sum when one of component is changed.
    totalAngularMomentum = sum(angularMomentum)
    angularMomentum = angularMomentum[:] + [totalAngularMomentum]
  else:
    # Sum already available.
    # Just use it.
    totalAngularMomentum = angularMomentum[IDX_Q_SUM]

  #
  # Precalculate constants needed during calculating expansion terms.
  # We can calculate these paremeters once before starting to work.
  #

  zetaLeft  = za + zb
  zetaRight = zc + zd
  zetaTotal = zetaLeft + zetaRight

  rho = zetaLeft * zetaRight / zetaTotal

  # Convert 4-center integral to 2-center using Gaussian product rule.
  # We can do it, because A,B depends on r1, and C,D depends on r2.
  # P = new center for A*B (left side, bra),
  rp = [
    (za * ra[0] + zb * rb[0]) / zetaLeft, # x
    (za * ra[1] + zb * rb[1]) / zetaLeft, # y
    (za * ra[2] + zb * rb[2]) / zetaLeft  # z
  ]

  # Q = new center for C*D (right side, ket).
  rq = [
    (zc * rc[0] + zd * rd[0]) / zetaRight, # x
    (zc * rc[1] + zd * rd[1]) / zetaRight, # y
    (zc * rc[2] + zd * rd[2]) / zetaRight  # z
  ]

  # W = new center for P*Q.
  rw = [
    (zetaLeft * rp[0] + zetaRight * rq[0]) / zetaTotal, # x
    (zetaLeft * rp[1] + zetaRight * rq[1]) / zetaTotal, # y
    (zetaLeft * rp[2] + zetaRight * rq[2]) / zetaTotal  # z
  ]

  # X-coord related.
  Ax = rp[0] - ra[0]
  Bx = rp[0] - rb[0]
  Cx = rq[0] - rc[0]
  Dx = rq[0] - rd[0]

  leftX  = rw[0] - rp[0]
  rightX = rw[0] - rq[0]

  # Y-coord related.
  Ay = rp[1] - ra[1]
  By = rp[1] - rb[1]
  Cy = rq[1] - rc[1]
  Dy = rq[1] - rd[1]

  leftY  = rw[1] - rp[1]
  rightY = rw[1] - rq[1]

  # Z-coord related.
  Az = rp[2] - ra[2]
  Bz = rp[2] - rb[2]
  Cz = rq[2] - rc[2]
  Dz = rq[2] - rd[2]

  leftZ  = rw[2] - rp[2]
  rightZ = rw[2] - rq[2]

  # Not related to specific coords, but depends on exponential
  # coefficients ("zetas").
  leftA  = 0.5 / zetaLeft  # bra side (AB)
  rightA = 0.5 / zetaRight # ket side (CD)

  leftAA  = -0.5 * rho / (zetaLeft  * zetaLeft)  # bra side (AB)
  rightAA = -0.5 * rho / (zetaRight * zetaRight) # ket side (CD)

  # Global - no coord (x,y,z), but depends on all functions together.
  # Appears in terms mixed between bra-ket sides (left/right, AB-CD).
  mixed = 0.5 / (zetaLeft + zetaRight)

  # Pack precalcuated parametrers in Context object.
  # This context will be passed internally during all Obara-Saika calls.
  params = [
    Ax, Bx, Cx, Dx, leftX, rightX,
    Ay, By, Cy, Dy, leftY, rightY,
    Az, Bz, Cz, Dz, leftZ, rightZ,

    leftA , rightA,
    leftAA, rightAA,
    mixed
  ]

  # Calculate expansion coefficient before each Boys(n) functions,
  # by applying Obara-Saika recursive scheme.
  expansion = [0] * (totalAngularMomentum + 1)
  _ObaraSaika_ERI_Expansion(expansion, angularMomentum, 1, 0, params)

  # Calculate argument passed to all Boys(order, x) functions.
  X = (rp[0] - rq[0]) * (rp[0] - rq[0]) + \
      (rp[1] - rq[1]) * (rp[1] - rq[1]) + \
      (rp[2] - rq[2]) * (rp[2] - rq[2])

  boysArg = rho * X

  # Calculate final integral value by adding up expansion terms:
  # [ab|cd] = C0 * Boys(0,x) + C1 + Boys(1,x) + ... + Cn * Boys(n,x)
  rv = 0
  for order in range(totalAngularMomentum + 1):
    rv += Utils.Boys(order, boysArg) * expansion[order]

  # Calculate extra coeficient for the whole integral.
  Rab2 = (ra[0] - rb[0])**2 + (ra[1] - rb[1])**2 + (ra[2] - rb[2])**2
  Rcd2 = (rc[0] - rd[0])**2 + (rc[1] - rd[1])**2 + (rc[2] - rd[2])**2

  Sab = exp(- za*zb / zetaLeft  * Rab2) / zetaLeft
  Scd = exp(- zc*zd / zetaRight * Rcd2) / zetaRight
  aux = 2 * pi**(5.0 / 2.0) * Sab * Scd / sqrt(zetaTotal)

  return aux * rv
