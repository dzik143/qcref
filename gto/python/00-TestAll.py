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

from math import fabs
from time import time

from ObaraSaika_1E  import ObaraSaika_Overlap, ObaraSaika_Nuclear, ObaraSaika_Kinetic
from ObaraSaika_ERI import ObaraSaika_ERI

# ------------------------------------------------------------------------------
#                                 Constants
# ------------------------------------------------------------------------------

ra = [0.1, 0.5, 1.2]
rb = [0.2, 0.7, 1.5]
rc = [0.3, 0.9, 1.8]
rd = [0.4, 1.1, 2.1]

za = 0.6
zb = 0.8
zc = 1.0
zd = 1.3

TOLERANCE = 4.5e-13

# ------------------------------------------------------------------------------
#                            Helper functions
# ------------------------------------------------------------------------------

def _overlap(q):
  return ObaraSaika_Overlap(za, zb, ra, rb, q)

def _kinetic(q):
  return ObaraSaika_Kinetic(za, zb, ra, rb, q)

def _nuclear(q):
  return ObaraSaika_Nuclear(za, zb, ra, rb, rc, q)

def _eri(q):
  return ObaraSaika_ERI(za, zb, zc, zd, ra, rb, rc, rd, q)

def _runTests(fname, fct):
  print('-------------------------------------------------------')
  print('Running tests "%s"...' % fname)
  print('-------------------------------------------------------')

  with open(fname) as f:
    for oneLine in f:
      # Remove end of line character.
      oneLine = oneLine.strip()

      if (len(oneLine) == 0) or (oneLine[0] == '#'):
        # Empty line or comment.
        print(oneLine)

      else:
        # Read patern and reference value to compare.
        tokens = oneLine.split(' ')

        q = []
        for i in range(len(tokens) - 1):
          q.append(int(tokens[i]))

        refValue = float(tokens[-1])

        # Run tested function over reference input parameters.
        t0 = time()
        currentValue = fct(q)
        timeElapsedMs = (time() - t0) * 1000

        # Compare result with reference value.
        print(q, '-> %+.15e... ' % refValue, end='')

        if fabs(currentValue - refValue) < TOLERANCE:
          if timeElapsedMs > 0.001:
            print('OK', '(%.1f ms)' % timeElapsedMs)
          else:
            print('OK')

        else:
          print('FAIL!')
          print('Test file      :', fname)
          print('Reference value:', refValue)
          print('Current value  :', currentValue)
          exit(-1)

# ------------------------------------------------------------------------------
#                               Entry point
# ------------------------------------------------------------------------------

_runTests('../data/ref-overlap.txt', _overlap)
_runTests('../data/ref-kinetic.txt', _kinetic)
_runTests('../data/ref-nuclear.txt', _nuclear)
_runTests('../data/ref-eri.txt'    , _eri)
