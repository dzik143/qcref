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

import Utils from './Utils.js';

const { exp, PI, floor } = Math;

// -----------------------------------------------------------------------------
//                                Constants
// -----------------------------------------------------------------------------

//                v------- total sum of all componens
// q = [ xyz xyz  s ]      and all functions
//       ^   ^
//       A    B
//       func func
//       tion tion
//
// Total length of q is 6 + 1 elements.
// Extra total sum is calculated internally only.
// Caller doesn't need to pass it to public API.

const IDX_Q_SUM = 6;

// -----------------------------------------------------------------------------
//                   Internal implementation (private)
// -----------------------------------------------------------------------------

function _ObaraSaika_OneElectronInternal(rv, q, currentTerm,
                                         currentOrder, schemeId, params) {
  if (q[IDX_Q_SUM] == 0) {
    // Already reduced to <000|...|000> = <s|...|s>
    switch (schemeId) {
      case 'S': rv[0] += currentTerm; break;
      case 'T': rv[1] += currentTerm; break;
      case 'N': rv[currentOrder] += currentTerm; break;
    }

  } else {
    // Non-zero total angular momentum.
    // Find first non-zero index to reduce.
    // We can choose any path to get correct result.
    //
    // Example:
    // --------
    //  component: xyz xyz
    //  function:  000 111
    //             --- ---
    //  Angular:   000 010
    //  Momentum        ^
    //  (L)             we reduce here
    //                  function no. 1
    //                  component y (1)
    //                  L1.y (-1 here)
    //
    //   We reduce:
    //   [ 000 010 ] -> [ 000 000 ]
    //          ^              ^
    //         -1 here         -1 here
    //

    let reduceIdx = 0;
    while (q[reduceIdx] == 0) {
      reduceIdx += 1;
    }

    const reduceFun = floor(reduceIdx / 3); // Function number (0,1)
    const reduceXyz = reduceIdx % 3;        // Component of angular momentum (x,y,z)

    // We have hard-coded Obara-Saika schemes to reduce left side ("a"
    // function). To keep it simple, we just swap a<->b to get function
    // being reduced always at "a" position (on left0.
    // This trick works due to integral symmetry with real functions (GTO):
    // <a|b> = <b|a>
    const a = reduceFun;
    const b = (a == 0) ? 1 : 0;

    // Fetch parameters.
    // Possible improvement: Don't use magic numbers?
    const Xi = params[reduceFun * 3 + reduceXyz];
    const Ci = params[6 + reduceXyz];

    const zetaCenter = params[9];

    //
    // Helper function to make one child call with reduced integral.
    // Possible improvement: Reuse code from ERI?
    //

    const _goDeeper = (factor, extraReduceFun, orderDelta, newschemeId) => {
      // Prepare new state for child node.
      const newQ     = [...q];
      let   newTerm  = currentTerm  * factor;
      const newOrder = currentOrder + orderDelta;

      // Reduce selected qi on each child.
      newQ[reduceIdx] -= 1;
      newQ[IDX_Q_SUM] -= 1;

      let goOn = true;

      // Extra reduction on child node if needed.
      if (extraReduceFun != null) {
        const extraReduceIdx = extraReduceFun * 3 + reduceXyz;

        if (newQ[extraReduceIdx] == 0) {
          // Already reduced to 0.
          // Nothing to do.
          goOn = false;

        } else {
          // Child angular momentum is still greater than 0.
          // Scale child term by reduced angular momentum.
          newTerm *= newQ[extraReduceIdx];

          // Reduce angular momentum by one and go on deeper.
          newQ[extraReduceIdx] -= 1;
          newQ[IDX_Q_SUM]      -= 1;
        }
      }

      // Go on deeper recursively if new integral exists (non-negative
      // angular momentum).
      if (goOn) {
        _ObaraSaika_OneElectronInternal(rv, newQ, newTerm,
                                            newOrder, newschemeId, params);
      }
    };

    //
    // Apply one of Obara-Saika recursion schemes.
    // At each call we reduce angular momentum by one and generate
    // up to 3 (S), 5 (T) or 6 (N) terms recurisively.
    //

    switch (schemeId) {
      case 'S': {
        _goDeeper(Xi         , null , 0, 'S');
        _goDeeper(zetaCenter , a    , 0, 'S');
        _goDeeper(zetaCenter , b    , 0, 'S');
        break;
      }

      case 'T': {
        // Possible improvement: Don't use magic numbers?
        const zetaCenterRho   = params[10];
        const zetaLeftOrRight = params[11 + reduceFun];

        _goDeeper(Xi              , null , 0, 'T');
        _goDeeper(zetaCenter      , a    , 0, 'T');
        _goDeeper(zetaCenter      , b    , 0, 'T');
        _goDeeper(zetaLeftOrRight , a    , 0, 'S');

        const newTerm = currentTerm * zetaCenterRho;
        _ObaraSaika_OneElectronInternal(rv, [...q], newTerm, 0, 'S', params);
        break;
      }

      case 'N': {
        _goDeeper(Xi          , null ,  0, 'N');
        _goDeeper(Ci          , null , +1, 'N');
        _goDeeper(zetaCenter  , a    ,  0, 'N');
        _goDeeper(-zetaCenter , a    , +1, 'N');
        _goDeeper(zetaCenter  , b    ,  0, 'N');
        _goDeeper(-zetaCenter , b    , +1, 'N');
        break;
      }
    }
  }
}

//
// conventional call - just set up environment and pass to underlying
// core function.
//

function _ObaraSaika_OneElectron(za, zb, ra, rb, rc, q, schemeId) {

  let totalAngularMomentum = 0;

  if (q.length == IDX_Q_SUM) {
    // Precalculate sum of all elements to detect [ss|ss] integrals easily.
    // We need to update this sum when one of component is changed.
    totalAngularMomentum = Utils.sum(q);
    q = [...q, totalAngularMomentum];

  } else {
    // Sum already available.
    // Just use it.
    totalAngularMomentum = q[IDX_Q_SUM];
  }

  // Set up parameters.
  const params = _ObaraSaika_CreateOneElectronParams(za, zb, ra, rb, rc);

  // Allocate output array.
  // Possible improvement: Possibility to use caller array?
  let rv = null;

  switch (schemeId) {
    case 'S': rv = [0]; break;
    case 'T': rv = [0, 0]; break;
    case 'N': rv = Array(totalAngularMomentum + 1).fill(0); break;
  }

  // Find expansion coefficient.
  _ObaraSaika_OneElectronInternal(rv, q, 1, 0, schemeId, params);

  // Pass just filled up array to the caller.
  return rv;
}

//
// Helper function to precalculate all parameters needed during
// integral expansion and store them in one Context like buffer.
// This context is passed thgrought _ObaraSaika_xxx() calls.
//

function _ObaraSaika_CreateOneElectronParams(za, zb, ra, rb, rc) {
  // Convert 2-center integral to 1-center using Gaussian product rule.
  // P = new center for A*B.
  const rp = Utils.gaussianProduct(za, zb, ra, rb);

  // Transform original A,B (and optionally C) centers into new
  // coordinate system relative to new RP center.
  // Left center (bra).
  const Ax = rp[0] - ra[0];
  const Ay = rp[1] - ra[1];
  const Az = rp[2] - ra[2];

  // Right center (ket).
  const Bx = rp[0] - rb[0];
  const By = rp[1] - rb[1];
  const Bz = rp[2] - rb[2];

  // Optional third center (e.g. for nuclear repulsion integral).
  let Cx = 0.0;
  let Cy = 0.0;
  let Cz = 0.0;

  if (rc != null) {
    Cx = rc[0] - rp[0];
    Cy = rc[1] - rp[1];
    Cz = rc[2] - rp[2];
  }

  // No-coords, but depends on zeta parameters on both sides
  // (bra/ket, left/right, AB).
  const zetaSum        = za + zb;
  const zetaCenter     = 0.5 / zetaSum;
  const zetaCenterRho  = za * zb / zetaSum;
  const zetaCenter2Rho = 2.0 * zetaCenterRho;

  // No-coords, but are specific to one side (bra/ket, left/right, A/B).
  const zetaLeft  = -zetaCenterRho / za;
  const zetaRight = -zetaCenterRho / zb;

  // Fill up common params struct.
  // This context will be passed throught ObaraSaika_Xxx() calls.
  const params = [
    Ax, Ay, Az,
    Bx, By, Bz,
    Cx, Cy, Cz,
    zetaCenter,
    zetaCenter2Rho,
    zetaLeft,
    zetaRight
  ];

  return params;
}

// -----------------------------------------------------------------------------
//                              Public API
// -----------------------------------------------------------------------------

//
// Calculate overlap integral with two S-type GTO orbitals.
//         /               /
// <s|s> = |a(r)*b(r) dr = |exp(-za r) * exp(-zb r) dr
//         /               /
//
// Parameters:
// -----------
//   ra[], rb[] - x,y,z coordinates of A and B functions,
//   za, zb     - exponential coefficients for A and B.
//
// Returns:
// --------
//   Value of <s|s> integral.
//

function ObaraSaika_Overlap_SS(za, zb, ra, rb) {
  const zetaSum = za + zb;
  const zetaRho = za * zb / zetaSum;
  const ra2     = Utils.distSquared(ra, rb);
  const rv      = exp(-zetaRho * ra2) * (PI / zetaSum)**1.5;
  return rv;
};

//
// Calculate overlap integral over two GTO orbitals with arbitrary
// angular momentums.
//
//         /
// <a|b> = |a(r)*b(r) dr
//         /
//
// Parameters:
// -----------
//   ra[], rb[] - x,y,z coordinates of A and B functions,
//   za, zb     - exponential coefficients for A and B.
//   q[]        - angular momentums [ Lax, Lay, Laz, Lbx, Lby, Lbz ]
//
// Returns:
// --------
//   Value of <a|b> integral.
//

export function ObaraSaika_Overlap(za, zb, ra, rb, q) {
  const aux       = ObaraSaika_Overlap_SS(za, zb, ra, rb);
  const expansion = _ObaraSaika_OneElectron(za, zb, ra, rb, null, q, 'S');
  return aux * expansion[0];
}

//
// Calculate kinetric energy integral (KEI) over two GTO orbitals with arbitrary
// angular momentums.
//
//  /        ( 1      1       1   )
//  | a(r) * ( ---- + ---- + ---- ) * b(r) dr
//  /        ( dx^2   dy^2   dz^2 )
//
// Parameters:
// -----------
//   ra[], rb[] - x,y,z coordinates of A and B functions,
//   za, zb     - exponential coefficients for A and B.
//   q[]        - angular momentums [ Lax, Lay, Laz, Lbx, Lby, Lbz ]
//
// Returns:
// --------
//   Value of kinetic energy intergal.
//

export function ObaraSaika_Kinetic(za, zb, ra, rb, q) {
  const aux     = ObaraSaika_Overlap_SS(za, zb, ra, rb);
  const zetaRho = Utils.harmonicMean(za, zb);
  const Rab2    = Utils.distSquared(ra, rb);

  const expansion = _ObaraSaika_OneElectron(za, zb, ra, rb, null, q, 'T');

  const rv = expansion[0] + 0.5 * zetaRho * (3 - zetaRho * Rab2) * expansion[1];

  return aux * rv;
}

//
// Calculate nuclear-electron attarction (NEI) integral over two GTO orbitals
// with arbitrary angular momentums.
//
//  /         1
//  | a(r) * ---- * b(r) dr
//  /        |r-C|
//
// Parameters:
// -----------
//   ra[], rb[] - x,y,z coordinates of A and B functions,
//   rc[]       - x,y,z coordinates of nuclei,
//   za, zb     - exponential coefficients for A and B.
//   q[]        - angular momentums [ Lax, Lay, Laz, Lbx, Lby, Lbz ]
//
// Returns:
// --------
//   Value of nuclear interaction integral.
//

export function ObaraSaika_Nuclear(za, zb, ra, rb, rc, q) {
  // Calculate expansion coefficient before each Boys(n) functions,
  // by applying Obara-Saika recursive scheme.
  const expansion = _ObaraSaika_OneElectron(za, zb, ra, rb, rc, q, 'N');

  // Calculate argument passed to all Boys(order, x) functions.
  const rp      = Utils.gaussianProduct(za, zb, ra, rb);
  const zetaSum = za + zb;
  const boysArg = zetaSum * Utils.distSquared(rp, rc);

  // Calculate final integral value by adding up expansion terms:
  // <a|1/rc|d> = C0 * Boys(0,x) + C1 + Boys(1,x) + ... + Cn * Boys(n,x)
  let rv = 0.0;
  expansion.forEach((Ci, order) => {
    rv += Ci * Utils.boys(order, boysArg);
  });

  // Calculate extra coeficient for the whole integral.
  const S00 = ObaraSaika_Overlap_SS(za, zb, ra, rb);
  const aux = -2.0 * (zetaSum / PI)**0.5 * S00;

  return aux * rv;
}
