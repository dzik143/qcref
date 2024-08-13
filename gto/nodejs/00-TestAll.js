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

const fs = require('fs');

const {
  ObaraSaika_Overlap,
  ObaraSaika_Nuclear,
  ObaraSaika_Kinetic
} = require('./ObaraSaika_1E.js');

const { ObaraSaika_ERI } = require('./ObaraSaika_ERI.js');

const { abs } = Math;

// -----------------------------------------------------------------------------
//                                Constants
// -----------------------------------------------------------------------------

const ra = [0.1, 0.5, 1.2];
const rb = [0.2, 0.7, 1.5];
const rc = [0.3, 0.9, 1.8];
const rd = [0.4, 1.1, 2.1];

const za = 0.6;
const zb = 0.8;
const zc = 1.0;
const zd = 1.3;

const TOLERANCE = 1.4e-9;

// -----------------------------------------------------------------------------
//                            Helper functions
// -----------------------------------------------------------------------------

function _overlap(q) {
  return ObaraSaika_Overlap(za, zb, ra, rb, q);
}

function _kinetic(q) {
  return ObaraSaika_Kinetic(za, zb, ra, rb, q);
}

function _nuclear(q) {
  return ObaraSaika_Nuclear(za, zb, ra, rb, rc, q);
}

function _eri(q) {
  return ObaraSaika_ERI(za, zb, zc, zd, ra, rb, rc, rd, q);
}

function _runTests(fname, fct) {
  console.log('-------------------------------------------------------')
  console.log(`Running tests "${fname}...`);
  console.log('-------------------------------------------------------')

  const lines = fs.readFileSync(fname).toString().split(/\r?\n/);

  lines.forEach((oneLine) => {
    // Remove end of line character.
    oneLine = oneLine.replace(/\r?\n|\r/g, '').trim();

    if ((oneLine.length == 0) || (oneLine[0] == '#')) {
      // Empty line or comment.
      console.log(oneLine);

    } else {
      // Read patern and reference value to compare.
      const tokens = oneLine.split(' ');

      const q = [];

      for (let idx = 0; idx < tokens.length - 1; idx++) {
        q.push(parseInt(tokens[idx]));
      }

      const refValue = parseFloat(tokens[tokens.length - 1]);

      // Run tested function over reference input parameters.
      const t0 = Date.now();
      const currentValue = fct(q);
      const timeElapsedMs = Date.now() - t0;

      // Compare result with reference value.
      err = abs(currentValue - refValue);

      if (err < TOLERANCE) {
        if (timeElapsedMs > 0) {
          console.log(oneLine, '... OK (', timeElapsedMs, 'ms )');
        } else {
          console.log(oneLine, '... OK');
        }

      } else {
        console.log(oneLine, ' ... FAIL!');
        console.log('Test file      :', fname);
        console.log('Reference value:', refValue);
        console.log('Current value  :', currentValue);
        console.log('Absolute error :', err);
        process.exit(-1);
      }
    }
  });
}

// -----------------------------------------------------------------------------
//                               Entry point
// -----------------------------------------------------------------------------

_runTests('../data/ref-overlap.txt', _overlap);
_runTests('../data/ref-kinetic.txt', _kinetic);
_runTests('../data/ref-nuclear.txt', _nuclear);
_runTests('../data/ref-eri.txt'    , _eri);
