/*******************************************************************************
*                                                                              *
* Copyright (C) 2024, 2024 Sylwester Wysocki (sw143@wp.pl)                     *
*                                                                              *
* This program is free software: you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License as published by         *
* the Free Software Foundation, either version 3 of the License, or            *
* (at your option) any later version.                                          *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program. If not, see <http://www.gnu.org/licenses/>          *
*                                                                              *
*******************************************************************************/

const gamma    = require('@stdlib/math-base-special-gamma')
const gammainc = require('@stdlib/math-base-special-gammainc')

// -----------------------------------------------------------------------------
//                               Public API
// -----------------------------------------------------------------------------

function distSquared(r1, r2) {
  distX = r1[0] - r2[0]
  distY = r1[1] - r2[1]
  distZ = r1[2] - r2[2]
  return distX**2 + distY**2 + distZ**2
}

function dist(r1, r2) {
  return Math.sqrt(DistSquared(r1, r2))
}

// Multiply two GTOs with centers A,B gives another GTO with different center.
// Function returns coords of this new center.
function gaussianProduct(za, zb, ra, rb) {
  return [
    (za * ra[0] + zb * rb[0]) / (za + zb), // x
    (za * ra[1] + zb * rb[1]) / (za + zb), // y
    (za * ra[2] + zb * rb[2]) / (za + zb)  // z
  ];
}

// https://en.wikipedia.org/wiki/Harmonic_mean#Harmonic_mean_of_two_or_three_numbers
function harmonicMean(x, y) {
  return 2 * x*y / (x + y);
}

// Calculate approximation of Boys function of order n at point x:
//
//         /1  2n         2
// Fn(x) = |  t  * exp(-xt ) dt
//         /0

function boys(order, x) {
  let rv = 0.0

  if (x > 0.0) {
    const g  = gamma(order + 0.5);
    const gi = 1.0 - gammainc(x, order + 0.5, true, true);
    rv = 0.5 * g * gi * x**(-order - 0.5);

  } else {
    rv = 1.0 / (2 * order + 1);
  }

  return rv;
}

function sum(tab) {
  let rv = 0;
  tab.forEach((v) => {
    rv += v;
  });
  return rv;
}

module.exports = {
  distSquared,
  dist,
  gaussianProduct,
  harmonicMean,
  boys,
  sum
};
