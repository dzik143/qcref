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

const { exp, PI, floor, sqrt } = Math;

// -----------------------------------------------------------------------------
//                                 Constants
// -----------------------------------------------------------------------------

//                       v------- total sum of all componens
// q = [ xyz xyz xyz xyz s ]      and all functions
//       ^   ^    ^    ^
//       A    B    C    D
//       func func func func
//       tion tion tion tion
//
// Total length of q is 12 + 1 elements.
// Extra total sum is calculated internally only.
// Caller doesn't need to pass it to public API.

const IDX_Q_SUM = 12;

const IDX_LEFT_A   = 18;
const IDX_RIGHT_A  = 19;
const IDX_LEFT_AA  = 20;
const IDX_RIGHT_AA = 21;
const IDX_MIXED    = 22;

// We use hard-coded Obara-Saika scheme to reduce [ab|cd] to [a-1b|cd].
// If we need to reduce angular momentum at different position (b,c,d),
// then we use one of ERI symmetry to move reduced function back to a.
//
// Below ERIs with real functions are equivalent:
// [ab|cd] = [ba|cd] = [ba|dc] = [ab|dc] = [cd|ab] = [cd|ba] = [dc|ba] = [dc|ab]
//
// We have 2 symmetries:
// - swap functions within the same bra/ket side: [ab|->[ba| and |cd]->|dc],
// - swap the whole pair between bra/ket sides: [ab|cd]->[cd|ab].
//
// Permutate original [ab|cd] functions depending on selected one.
// Below map shows new functions order within integral.
const MUTATE_MAP = [
// a, b, c, d    New [a,b,c,d] if...
  0, 1, 2, 3, // ... function no. 0 chosen (old a)
  1, 0, 2, 3, // ... function no. 1 chosen (old b)
  2, 3, 0, 1, // ... function no. 2 chosen (old c)
  3, 2, 0, 1  // ... function no. 3 chosen (old d)
// ^
// function to reduce
// always go to new
// a
];

// -----------------------------------------------------------------------------
//                    Inernal implementation (private)
// -----------------------------------------------------------------------------

//
// Find expansion coefficient for 2-electron repulsion (ERI) integral
// using Obara-Saika scheme.
//
// How does it work:
// -----------------
//   Ref. Obara, S.; Saika, A. J Chem Phys 1986, 84, 3963.
//   https://doi.org/10.1063/1.450106
//
//   At each call we reduce total angular momentum by 1.
//   Eg. [200 000 000 110] -> [100 000 000 110]
//
//   Each call generates up to 8 child calls.
//   We stop recursion if we reduce integral to basic [ss|ss].
//
// How to use results:
// -------------------
// Final integral is:
//   [ab|cd] = C0 * F(0,t) + C1*F1(1,t) + ... + Cn*Fn(n,t)
//
//   - Ci = coefficients calculated by this function (see rv[] param)
//   - F(i,t) = Boys function of i-th order
//
//   - t is the same for all orders (can be precomputed):
//
//        (za + zb) * (zc + zd)       2
//   t =  --------------------- * |PQ|
//        za + zb + zc + zd
//
//   |PQ| - distance between P and Q centers.
//    P   - bi-center on "left" side (AB),
//    Q   - bi-ecnter of "right" side (CD)
//
// Parameters:
// -----------
//   - rv[] - calculated coefficient grouped by Boys() function
//            order (OUT)
//
//   - q[]  - angular momentum of each functions
//            [LAx,LAy,LAy, LBx,LBy,LBz, LCx,LCy,LCz, LDx,LDy,LDz] (IN)
//
//   - currentTerm  - currently accumulated term passed between calls (IN)
//
//   - currentOrder - current accumulated order of Boys() function (IN)
//
//   - params - constants calcualted once, needed to calculate
//              coefficient (IN)
//
// Returns:
// --------
//   After finished rv[] contains coefficients needed to be multiplicated
//   by related Boys(s,t) functions.
//   See ObaraSaika_ERI() for more.
//
function _ObaraSaika_ERI_Expansion(rv, q, currentTerm, currentOrder, params) {
  if (q[IDX_Q_SUM] == 0) {
    // Already reduced to [000 000 | 000 000] = [ss|ss]
    // Accumulate terms per Boys function order.
    rv[currentOrder] += currentTerm;

  } else {
    // Non-zero total angular momentum.
    // Find first non-zero index to reduce.
    // We can choose any path to get correct result.
    //
    // Example:
    // --------
    //  component: xyz xyz xyz xyz
    //  function:  000 111 222 333
    //             --- --- --- ---
    //  Angular:   000 010 200 002
    //  Momentum        ^
    //  (L)             we reduce here
    //                  function no. 1
    //                  component y (1)
    //                  L1.y (-1 here)
    //
    //   We reduce:
    //   [ 000 010 200 002 ] -> [ 000 000 200 002 ]
    //          ^                      ^
    //         -1 here                 -1 here
    //

    let reduceIdx = 0;
    while (q[reduceIdx] == 0) {
      reduceIdx += 1;
    }

    const reduceFun = floor(reduceIdx / 3); // Function number (0,1,2,3)
    const reduceXyz = reduceIdx %  3;       // Component of angular momentum (x,y,z)

    // Check does selected function belong to <bra| (left) or |ket> (right) side.
    // [01|23] = <bra| 1/r2 |ket>
    const isKet = (reduceFun > 1) ? 1 : 0;

    //
    // Fetch parameters for function being reduced (A,B,C,D) and
    // selected angular momentum component (i=x,y,z)
    //

    // Xi = Ai, Bi, Ci, Di
    // Possible improvement: Remove magic numbers?
    const Xi = params[reduceXyz * 6 + reduceFun];

    // leftX/Y/Z or rightX/Y/Z
    // Possible improvement: Remove magic numbers?
    const leftOrRightXyz = params[reduceXyz * 6 + isKet + 4];

    // Global - no coord (x,y,z), but depends on all functions together.
    // Related strictly to one side (bra/ket, left/right)
    const leftOrRightA   = params[IDX_LEFT_A  + isKet];
    const leftOrRightAA  = params[IDX_LEFT_AA + isKet];

    // Appears for terms with mixed bra-ket terms.
    const mixed = params[IDX_MIXED];

    //
    // Permutate original [ab|cd] integral to keep function, which we're going
    // to reduce at "a" position.
    // Due to this we can apply the same hard-coded Obara-Saika scheme for
    // each integrals.
    // This trick works due to ERIs symmetries e.g. [ab|cd] = [ba|cd] etc.
    //

    const a = MUTATE_MAP[reduceFun * 4 + 0]; // <- function to reduce goes here
    const b = MUTATE_MAP[reduceFun * 4 + 1];
    const c = MUTATE_MAP[reduceFun * 4 + 2];
    const d = MUTATE_MAP[reduceFun * 4 + 3];

    //
    // Helper function to make one child call with reduced integral.
    //

    const _goDeeper = (factor, extraReduceFun, orderDelta) => {
      // Prepare new state for child node.
      const newQ     = [...q];
      let   newTerm  = currentTerm  * factor;
      const newOrder = currentOrder + orderDelta;

      // Reduce selected qi on each child.
      newQ[reduceIdx] -= 1;
      newQ[IDX_Q_SUM] -= 1;

      let goOn = true;

      // Extra reduction on child node if (needed.
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

      // Go on deeper recursively if (new integral exists (non-negative
      // angular momentum).
      if (goOn) {
        _ObaraSaika_ERI_Expansion(rv, newQ, newTerm, newOrder, params);
      }
    };

    // Apply Obara-Saika scheme.
    // Each reduction generate up to 8 child nodes.
    // Note, that some of these functions doesn't make sense,
    // becase after subtraction we may get negative angular momentum.
    // These "fake" functions will be skipped.

    // Direct implementation of eq.39 from:
    //   Obara, S.; Saika, A. J Chem Phys 1986, 84, 3963.
    //   https://doi.org/10.1063/1.450106
    _goDeeper(Xi             , null ,  0);
    _goDeeper(leftOrRightXyz , null , +1);
    _goDeeper(leftOrRightA   ,    a ,  0);
    _goDeeper(leftOrRightAA  ,    a , +1);
    _goDeeper(leftOrRightA   ,    b ,  0);
    _goDeeper(leftOrRightAA  ,    b , +1);
    _goDeeper(mixed          ,    c , +1);
    _goDeeper(mixed          ,    d , +1);
  }
}

// -----------------------------------------------------------------------------
//                              Public API
// -----------------------------------------------------------------------------

//
// Calculate 2-electron repulsion (ERI) integral using Obara-Saika scheme.
//           //                 1
// [ab|cd] = || a(r1)*b(r1) * ----- * c(r2) * d(r2) dr1 dr2
//           //                r12
//
// r1 = [x1, y1, z1] = position of 1-st electron
// r2 = [x2, y2, z2] = position of 2-nd electron
//
// a,b,c,d - primitive Gaussian functions:
//             Lx        Ly        Lz                   2
//   a = (x-Ax)   * (y-Ay)  * (z-Az)  * exp(zetaA * |r-A| )
//   etc.
//
// Parameters:
// -----------
//   za,zb,zc,zd - exponential coefficient of A,B,C,D functions ("zeta") (IN)
//   ra,rb,rc,rd - centers of A,B,C,D functions (IN)
//   angularMomentum - angular momentum of each functions (IN)
//
// Angular momentum contains xyz components for each function.
// Each angular momentum has 3 components: lx, ly, lz.
//
// Example:
// --------
//   za = zb = zc = zd = 1
//   ra = rb = rc = rd = [0,0,0]
//   angularMomentum   = [1,0,0, 0,0,0, 0,1,1, 2,0,0]
//
//   All GTOs are centered on (0,0,0).
//   All GTOs have zeta (exponent coefficient) set to 1.
//   We calcualte integral [px s | dyz dxx] = [100 000 | 011 200].
//
// Returns:
// --------
//   Value of [ab|cd] integral (ERI).
//

export function ObaraSaika_ERI(za, zb, zc, zd, ra, rb, rc, rd, angularMomentum) {

  let totalAngularMomentum = 0;

  if (angularMomentum.length == IDX_Q_SUM) {
    // Precalculate sum of all elements to detect [ss|ss] integrals easily.
    // We need to update this sum when one of component is changed.
    totalAngularMomentum = Utils.sum(angularMomentum);
    angularMomentum = [...angularMomentum, totalAngularMomentum];

  } else {
    // Sum already available.
    // Just use it.
    totalAngularMomentum = angularMomentum[IDX_Q_SUM];
  }

  //
  // Precalculate constants needed during finding out expansion terms.
  // We can calculate these paremeters once before starting to work.
  //

  const zetaLeft  = za + zb;
  const zetaRight = zc + zd;
  const zetaTotal = zetaLeft + zetaRight;

  const rho = zetaLeft * zetaRight / zetaTotal;

  // Convert 4-center integral to 2-center using Gaussian product rule.
  // We can do it, because A,B depends on r1, and C,D depends on r2.
  // P = new center for A*B (left side, bra),
  const rp = Utils.gaussianProduct(za, zb, ra, rb);

  // Q = new center for C*D (right side, ket).
  const rq = Utils.gaussianProduct(zc, zd, rc, rd);

  // W = new center for P*Q.
  const rw = Utils.gaussianProduct(zetaLeft, zetaRight, rp, rq);

  // X-coord related.
  const Ax = rp[0] - ra[0];
  const Bx = rp[0] - rb[0];
  const Cx = rq[0] - rc[0];
  const Dx = rq[0] - rd[0];

  const leftX  = rw[0] - rp[0];
  const rightX = rw[0] - rq[0];

  // Y-coord related.
  const Ay = rp[1] - ra[1];
  const By = rp[1] - rb[1];
  const Cy = rq[1] - rc[1];
  const Dy = rq[1] - rd[1];

  const leftY  = rw[1] - rp[1];
  const rightY = rw[1] - rq[1];

  // Z-coord related.
  const Az = rp[2] - ra[2];
  const Bz = rp[2] - rb[2];
  const Cz = rq[2] - rc[2];
  const Dz = rq[2] - rd[2];

  const leftZ  = rw[2] - rp[2];
  const rightZ = rw[2] - rq[2];

  // Not related to specific coords, but depends on exponential
  // coefficients ("zetas").
  const leftA  = 0.5 / zetaLeft;  // bra side (AB)
  const rightA = 0.5 / zetaRight; // ket side (CD)

  const leftAA  = -0.5 * rho / (zetaLeft  * zetaLeft);  // bra side (AB)
  const rightAA = -0.5 * rho / (zetaRight * zetaRight); // ket side (CD)

  // Global - no coord (x,y,z), but depends on all functions together.
  // Appears in terms mixed between bra-ket sides (left/right, AB-CD).
  const mixed = 0.5 / (zetaLeft + zetaRight);

  // Pack precalcuated parametrers in Context object.
  // This context will be passed internally during all Obara-Saika calls.
  const params = [
    Ax, Bx, Cx, Dx, leftX, rightX,
    Ay, By, Cy, Dy, leftY, rightY,
    Az, Bz, Cz, Dz, leftZ, rightZ,

    leftA , rightA,
    leftAA, rightAA,
    mixed
  ];

  // Calculate expansion coefficient before each Boys(n) functions,
  // by applying Obara-Saika recursive scheme.
  const expansion = Array(totalAngularMomentum + 1).fill(0);
  _ObaraSaika_ERI_Expansion(expansion, angularMomentum, 1, 0, params);

  // Calculate argument passed to all Boys(order, x) functions.
  const boysArg = rho * Utils.distSquared(rp, rq);

  // Calculate final integral value by adding up expansion terms:
  // [ab|cd] = C0 * Boys(0,x) + C1 + Boys(1,x) + ... + Cn * Boys(n,x)
  let rv = 0;
  expansion.forEach((Ci, order) => {
    rv += Ci * Utils.boys(order, boysArg);
  });

  // Calculate extra coeficient for the whole integral.
  const Rab2 = Utils.distSquared(ra, rb);
  const Rcd2 = Utils.distSquared(rc, rd);

  const Sab = exp(-za*zb / zetaLeft  * Rab2) / zetaLeft;
  const Scd = exp(-zc*zd / zetaRight * Rcd2) / zetaRight;

  const aux = 2 * PI**(5/2) * Sab * Scd / sqrt(zetaTotal);

  return aux * rv;
}
