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

from math import pi, exp

#-------------------------------------------------------------------------------
#                                 Constants
#-------------------------------------------------------------------------------

#                v------- total sum of all componens
# q = [ xyz xyz  s ]      and all functions
#       ^   ^
#       A    B
#       func func
#       tion tion
#
# Total length of q is 6 + 1 elements.
# Extra total sum is calculated internally only.
# Caller doesn't need to pass it to public API.

IDX_Q_SUM = 6

#-------------------------------------------------------------------------------
#                   Internal implementation (private)
#-------------------------------------------------------------------------------

def _ObaraSaika_OneElectronInternal(rv, q, currentTerm,
                                        currentOrder, schemeId, params):
  if q[IDX_Q_SUM] == 0:
    # Already reduced to <000|...|000> = <s|...|s>
    match schemeId:
      case 'S': rv[0] += currentTerm
      case 'T': rv[1] += currentTerm
      case 'N': rv[currentOrder] += currentTerm

  else:
    # Non-zero total angular momentum.
    # Find first non-zero index to reduce.
    # We can choose any path to get correct result.
    #
    # Example:
    # --------
    #  component: xyz xyz
    #  function:  000 111
    #             --- ---
    #  Angular:   000 010
    #  Momentum        ^
    #  (L)             we reduce here
    #                  function no. 1
    #                  component y (1)
    #                  L1.y (-1 here)
    #
    #   We reduce:
    #   [ 000 010 ] -> [ 000 000 ]
    #          ^              ^
    #         -1 here         -1 here
    #

    reduceIdx = 0
    while q[reduceIdx] == 0:
      reduceIdx += 1

    reduceFun = reduceIdx // 3 # Function number (0,1)
    reduceXyz = reduceIdx %  3 # Component of angular momentum (x,y,z)

    # We have hard-coded Obara-Saika schemes to reduce left side ("a"
    # function). To keep it simple, we just swap a<->b to get function
    # being reduced always at "a" position (on left0.
    # This trick works due to integral symmetry with real functions (GTO):
    # <a|b> = <b|a>
    a = reduceFun
    b = 1 if (a == 0) else 0

    # Fetch parameters.
    # Possible improvement: Don't use magic numbers?
    Xi = params[reduceFun * 3 + reduceXyz]
    Ci = params[6 + reduceXyz]

    zetaCenter = params[9]

    #
    # Helper function to make one child call with reduced integral.
    # Possible improvement: Reuse code from ERI?
    #

    def _goDeeper(factor, extraReduceFun, orderDelta, newschemeId):
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
        _ObaraSaika_OneElectronInternal(rv, newQ, newTerm,
                                            newOrder, newschemeId, params)

    #
    # Apply one of Obara-Saika recursion schemes.
    # At each call we reduce angular momentum by one and generate
    # up to 3 (S), 5 (T) or 6 (N) terms recurisively.
    #

    match schemeId:
      case 'S':
        _goDeeper(Xi         , None , 0, 'S')
        _goDeeper(zetaCenter , a    , 0, 'S')
        _goDeeper(zetaCenter , b    , 0, 'S')

      case 'T':
        # Possible improvement: Don't use magic numbers?
        zetaCenterRho   = params[10]
        zetaLeftOrRight = params[11 + reduceFun]

        _goDeeper(Xi              , None , 0, 'T')
        _goDeeper(zetaCenter      , a    , 0, 'T')
        _goDeeper(zetaCenter      , b    , 0, 'T')
        _goDeeper(zetaLeftOrRight , a    , 0, 'S')

        newTerm = currentTerm * zetaCenterRho
        _ObaraSaika_OneElectronInternal(rv, q[:], newTerm, 0, 'S', params)

      case 'N':
        _goDeeper(Xi          , None ,  0, 'N')
        _goDeeper(Ci          , None , +1, 'N')
        _goDeeper(zetaCenter  , a    ,  0, 'N')
        _goDeeper(-zetaCenter , a    , +1, 'N')
        _goDeeper(zetaCenter  , b    ,  0, 'N')
        _goDeeper(-zetaCenter , b    , +1, 'N')

#-------------------------------------------------------------------------------

#
# conventional call - just set up environment and pass to underlying
# core function.
#

def _ObaraSaika_OneElectron(za, zb, ra, rb, rc, q, schemeId):

  if len(q) == IDX_Q_SUM:
    # Precalculate sum of all elements to detect [ss|ss] integrals easily.
    # We need to update this sum when one of component is changed.
    totalAngularMomentum = sum(q)
    q = q[:] + [totalAngularMomentum]
  else:
    # Sum already available.
    # Just use it.
    totalAngularMomentum = q[IDX_Q_SUM]

  # Set up parameters.
  params = _ObaraSaika_CreateOneElectronParams(za, zb, ra, rb, rc, q)

  # Allocate output array.
  # Possible improvement: Possibility to use caller array?
  match schemeId:
    case 'S': rv = [0]
    case 'T': rv = [0, 0]
    case 'N': rv = [0] * (totalAngularMomentum + 1)

  # Find expansion coefficient.
  _ObaraSaika_OneElectronInternal(rv, q, 1, 0, schemeId, params)

  # Pass just filled up array to the caller.
  return rv

#-------------------------------------------------------------------------------

#
# Helper function to precalculate all parameters needed during
# integral expansion and store them in one Context like buffer.
# This context is passed thgrought _ObaraSaika_xxx() calls.
#

def _ObaraSaika_CreateOneElectronParams(za, zb, ra, rb, rc, q):
  # Convert 2-center integral to 1-center using Gaussian product rule.
  # P = new center for A*B.
  rp = Utils.GaussianProduct(za, zb, ra, rb)

  # Transform original A,B (and optionally C) centers into new
  # coordinate system relative to new RP center.
  # Left center (bra).
  Ax = rp[0] - ra[0]
  Ay = rp[1] - ra[1]
  Az = rp[2] - ra[2]

  # Right center (ket).
  Bx = rp[0] - rb[0]
  By = rp[1] - rb[1]
  Bz = rp[2] - rb[2]

  # Optional third center (e.g. for nuclear repulsion integral).
  Cx = Cy = Cz = 0.0
  if rc != None:
    Cx = rc[0] - rp[0]
    Cy = rc[1] - rp[1]
    Cz = rc[2] - rp[2]

  # No-coords, but depends on zeta parameters on both sides
  # (bra/ket, left/right, AB).
  zetaSum        = za + zb
  zetaCenter     = 0.5 / zetaSum
  zetaCenterRho  = za * zb / zetaSum
  zetaCenter2Rho = 2.0 * zetaCenterRho

  # No-coords, but are specific to one side (bra/ket, left/right, A/B).
  zetaLeft  = -zetaCenterRho / za
  zetaRight = -zetaCenterRho / zb

  # Fill up common params struct.
  # This context will be passed throught ObaraSaika_Xxx() calls.
  params = [
    Ax, Ay, Az,
    Bx, By, Bz,
    Cx, Cy, Cz,
    zetaCenter,
    zetaCenter2Rho,
    zetaLeft,
    zetaRight
  ]

  return params

#-------------------------------------------------------------------------------
#                              Public API
#-------------------------------------------------------------------------------

#
# Calculate overlap integral with two S-type GTO orbitals.
#         /               /
# <s|s> = |a(r)*b(r) dr = |exp(-za r) * exp(-zb r) dr
#         /               /
#
# Parameters:
# -----------
#   ra[], rb[] - x,y,z coordinates of A and B functions,
#   za, zb     - exponential coefficients for A and B.
#
# Returns:
# --------
#   Value of <s|s> integral.
#

def ObaraSaika_Overlap_SS(za, zb, ra, rb):
  zetaSum = za + zb
  zetaRho = za * zb / zetaSum
  ra2     = Utils.DistSquared(ra, rb)
  rv      = exp(-zetaRho * ra2) * (pi / zetaSum)**1.5
  return rv

# ------------------------------------------------------------------------------

#
# Calculate overlap integral over two GTO orbitals with arbitrary
# angular momentums.
#
#         /
# <a|b> = |a(r)*b(r) dr
#         /
#
# Parameters:
# -----------
#   ra[], rb[] - x,y,z coordinates of A and B functions,
#   za, zb     - exponential coefficients for A and B.
#   q[]        - angular momentums [ Lax, Lay, Laz, Lbx, Lby, Lbz ]
#
# Returns:
# --------
#   Value of <a|b> integral.
#

def ObaraSaika_Overlap(za, zb, ra, rb, q):
  aux       = ObaraSaika_Overlap_SS(za, zb, ra, rb)
  expansion = _ObaraSaika_OneElectron(za, zb, ra, rb, None, q, 'S')
  return aux * expansion[0]

# ------------------------------------------------------------------------------

#
# Calculate kinetric energy integral (KEI) over two GTO orbitals with arbitrary
# angular momentums.
#
#  /        ( 1      1       1   )
#  | a(r) * ( ---- + ---- + ---- ) * b(r) dr
#  /        ( dx^2   dy^2   dz^2 )
#
# Parameters:
# -----------
#   ra[], rb[] - x,y,z coordinates of A and B functions,
#   za, zb     - exponential coefficients for A and B.
#   q[]        - angular momentums [ Lax, Lay, Laz, Lbx, Lby, Lbz ]
#
# Returns:
# --------
#   Value of kinetic energy intergal.
#

def ObaraSaika_Kinetic(za, zb, ra, rb, q):

  aux     = ObaraSaika_Overlap_SS(za, zb, ra, rb)
  zetaRho = Utils.HarmonicMean(za, zb)
  Rab2    = Utils.DistSquared(ra, rb)

  expansion = _ObaraSaika_OneElectron(za, zb, ra, rb, None, q, 'T')

  rv = expansion[0] + \
     + 0.5 * zetaRho * (3 - zetaRho * Rab2) * expansion[1]

  return aux * rv

# ------------------------------------------------------------------------------

#
# Calculate nuclear-electron attarction (NEI) integral over two GTO orbitals
# with arbitrary angular momentums.
#
#  /         1
#  | a(r) * ---- * b(r) dr
#  /        |r-C|
#
# Parameters:
# -----------
#   ra[], rb[] - x,y,z coordinates of A and B functions,
#   rc[]       - x,y,z coordinates of nuclei,
#   za, zb     - exponential coefficients for A and B.
#   q[]        - angular momentums [ Lax, Lay, Laz, Lbx, Lby, Lbz ]
#
# Returns:
# --------
#   Value of nuclear interaction integral.
#

def ObaraSaika_Nuclear(za, zb, ra, rb, rc, q):
  # Calculate expansion coefficient before each Boys(n) functions,
  # by applying Obara-Saika recursive scheme.
  expansion = _ObaraSaika_OneElectron(za, zb, ra, rb, rc, q, 'N')

  # Calculate argument passed to all Boys(order, x) functions.
  rp      = Utils.GaussianProduct(za, zb, ra, rb)
  zetaSum = za + zb
  boysArg = zetaSum * Utils.DistSquared(rp, rc)

  # Calculate final integral value by adding up expansion terms:
  # <a|1/rc|d> = C0 * Boys(0,x) + C1 + Boys(1,x) + ... + Cn * Boys(n,x)
  rv = 0.0
  for order in range(sum(q) + 1):
    rv += expansion[order] * Utils.Boys(order, boysArg)

  # Calculate extra coeficient for the whole integral.
  S00 = ObaraSaika_Overlap_SS(za, zb, ra, rb)
  aux = -2.0 * (zetaSum / pi)**0.5 * S00

  return aux * rv
