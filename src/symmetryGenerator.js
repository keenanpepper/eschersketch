//--------------------------------------------------------------------------------------------------
//
// Eschersketch - A drawing program for exploring symmetrical designs
//
// Geometry Functions and Routines
//
// Copyright (c) 2017 Anselm Levskaya (http://anselmlevskaya.com)
// Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php) license.
//
// Points
// Matrix Transforms
// Affine Transforms
//
//--------------------------------------------------------------------------------------------------


// import core math to local namespace
//--------------------------------------------------------------------------------------------------
import { _ } from 'underscore';

const { min, max, abs, sqrt, floor, round, sin, cos, tan, acos, asin, atan, pow, PI } = Math;
const sign = // sign func -> -1, 1
    function(x) { if (x < 0) { return -1; } else { return 1; } };
const map =  // linear map of i-range onto o-range
    (value, istart, istop, ostart, ostop) =>
    ostart + (((ostop - ostart) * (value - istart)) / (istop - istart));


// 2d Point Class
//--------------------------------------------------------------------------------------------------

class Point2d {
  constructor(x, y) {
    this.x = x || 0;
    this.y = y || 0;
  }

  fromArray(ar) {
    return new Point2(ar[0], ar[1]);
  }

  toArray() {
    return [ x, y ];
  }

  plus(pIn) {
    const pOut = new Point2();
    pOut.x = this.x + pIn.x;
    pOut.y = this.y + pIn.y;
    return pOut;
  }

  _plus(pIn) {
    this.x = this.x + pIn.x;
    this.y = this.y + pIn.y;
    return this;
  }

  minus(pIn) {
    const pOut = new Point2();
    pOut.x = this.x - pIn.x;
    pOut.y = this.y - pIn.y;
    return pOut;
  }

  _minus(pIn) {
    this.x = this.x - pIn.x;
    this.y = this.y - pIn.y;
    return this;
  }

  times(a) {
    const pOut = new Point2();
    pOut.x = a * this.x;
    pOut.y = a * this.y;
    return pOut;
  }

  _times(pIn) {
    this.x = a * this.x;
    this.y = a * this.y;
    return this;
  }

  dot(pIn) {
    return (this.x * pIn.x) + (this.y * pIn.y);
  }

  norm() {
    return sqrt((this.x * this.x) + (this.y * this.y));
  }

  angleWith(pIn) {
    const thisnorm = sqrt((this.x * this.x) + (this.y * this.y));
    const pInnorm = sqrt((pIn.x * pIn.x) + (pIn.y * pIn.y));
    return acos(((this.x * pIn.x) + (this.y * pIn.y)) / thisnorm / pInnorm);
  }

  cross(pIn) {
    return (this.x * pIn.y) - (this.y * pIn.x);
  }
}

const radialPoint =
    (r, theta) => new Point2d(r * cos(theta), r * sin(theta));



// 2D Linear Transform Class
//--------------------------------------------------------------------------------------------------

class Matrix2 {
    constructor(a, b, c, d) {
      this.a = a || 1;
      this.b = b || 0;
      this.c = c || 0;
      this.d = d || 1;
    }

    fromArray(ar) {
      return new Matrix2(ar[0][0], ar[0][1], ar[1][0], ar[1][1]);
    }

    add(matIn) {
      const matOut = new Matrix2();
      matOut.a = this.a + matIn.a;
      matOut.b = this.a + matIn.b;
      matOut.c = this.a + matIn.c;
      matOut.d = this.a + matIn.d;
      return matOut;
    }

    multiply(matIn) {
      const matOut = new Matrix2();
      matOut.a = (this.a * matIn.a) + (this.b * matIn.c);
      matOut.b = (this.a * matIn.b) + (this.b * matIn.d);
      matOut.c = (this.c * matIn.a) + (this.d * matIn.c);
      matOut.d = (this.c * matIn.b) + (this.d * matIn.d);
      return matOut;
    }

    inverse() {
      const matOut = new Matrix2();
      const det = (this.a * this.d) - (this.b * this.c);
      matOut.a = this.d / det;
      matOut.b = (-1 * this.b) / det;
      matOut.c = (-1 * this.c) / det;
      matOut.d = this.a / det;
      return matOut;
    }
}


// 2D Affine Transform Class
//--------------------------------------------------------------------------------------------------

class AffineTransform {
  constructor(a, b, c, d, x, y) {
    this.a = a || 1;
    this.b = b || 0;
    this.c = c || 0;
    this.d = d || 1;
    this.x = x || 0;
    this.y = y || 0;
  }

  // A.multiply(B) = A*B
  multiply(Afin) {
    const Afout = new AffineTransform();
    Afout.a = (this.a * Afin.a) + (this.b * Afin.c);
    Afout.b = (this.a * Afin.b) + (this.b * Afin.d);
    Afout.c = (this.c * Afin.a) + (this.d * Afin.c);
    Afout.d = (this.c * Afin.b) + (this.d * Afin.d);
    Afout.x = this.x + (this.a * Afin.x) + (this.b * Afin.y);
    Afout.y = this.y + (this.c * Afin.x) + (this.d * Afin.y);
    return Afout;
  }

  // A.Lmultiply(B) = B*A
  Lmultiply(Afin) {
    const Afout = new AffineTransform();
    Afout.a = (Afin.a * this.a) + (Afin.b * this.c);
    Afout.b = (Afin.a * this.b) + (Afin.b * this.d);
    Afout.c = (Afin.c * this.a) + (Afin.d * this.c);
    Afout.d = (Afin.c * this.b) + (Afin.d * this.d);
    Afout.x = Afin.x + (Afin.a * this.x) + (Afin.b * this.y);
    Afout.y = Afin.y + (Afin.c * this.x) + (Afin.d * this.y);
    return Afout;
  }

  inverse() {
    const Afout = new AffineTransform();
    const det = (this.a * this.d) - (this.b * this.c);
    Afout.a = this.d / det;
    Afout.b = (-1 * this.b) / det;
    Afout.c = (-1 * this.c) / det;
    Afout.d = this.a / det;
    Afout.x = ((-1 * this.d) + this.x + (this.b * this.y)) / det;
    Afout.y = ((-1 * this.a) + this.y + (this.c * this.x)) / det;
    return Afout;
  }

  on(x, y) {
    const nx = (this.a * x) + (this.b * y) + this.x;
    const ny = (this.c * x) + (this.d * y) + this.y;
    return [ nx, ny ];
  }

  onVec(x) {
    const nx = (this.a * x[0]) + (this.b * x[1]) + this.x;
    const ny = (this.c * x[0]) + (this.d * x[1]) + this.y;
    return [ nx, ny ];
  }

  toList() {
    return [ this.a, this.b, this.c, this.d, this.x, this.y ];
  }

  // approximate comparison function
  sameAs(Afin, tol) {
    tol = tol || 1e-8;
    let sum = 0;
    sum += abs(this.a - Afin.a);
    sum += abs(this.b - Afin.b);
    sum += abs(this.c - Afin.c);
    sum += abs(this.d - Afin.d);
    sum += abs(this.x - Afin.x);
    sum += abs(this.y - Afin.y);
    if (sum < tol) {
      return true;
    } else {
      return false;
    }
  }
}


// Common Affine Transforms
//--------------------------------------------------------------------------------------------------

export const IdentityTransform =
    () => new AffineTransform(1, 0, 0, 1, 0, 0);

export const RotationTransform =
    angle => new AffineTransform(Math.cos(angle), Math.sin(angle), -1 * Math.sin(angle), Math.cos(angle), 0, 0);

export const TranslationTransform =
    (dx, dy) => new AffineTransform(1, 0, 0, 1, dx, dy);

export const ReflectionTransform =
    angle => new AffineTransform(Math.cos(2 * angle), Math.sin(2 * angle), Math.sin(2 * angle), -1 * Math.cos(2 * angle), 0, 0);

export const RotationAbout =
    (angle, px, py) => TranslationTransform(px, py).multiply(RotationTransform(angle)).multiply(TranslationTransform(-px, -py));

export const ReflectionAbout =
    (angle, px, py) => TranslationTransform(px, py).multiply(ReflectionTransform(angle)).multiply(TranslationTransform(-px, -py));

export const GlideTransform =
    (angle, distance, px, py) => ReflectionAbout(angle, px, py).multiply(TranslationTransform(distance * cos(angle), distance * sin(angle)));


// Functions for manipulating Sets of Affine Transforms
//--------------------------------------------------------------------------------------------------

const setProduct = function(X, Y, prodfunc) {
  prodfunc = prodfunc || ((x, y) => [ x, y ]);
  return _.reduce(X, ((memo, x) =>
    _.map(Y, y => prodfunc(x, y)).concat(memo)
  ), []);
};

const affinesetproduct =
    (Afset1, Afset2) => setProduct(Afset1, Afset2, (x, y) => x.multiply(y));

// generates unique subset of array ar using equivalency function eqfunc
// if sortkey is supplied, it must satisfy the following:
// if sortkey(A) and sortkey(B) differ by at least tol, then eqfunc(A,B) must be false 
const uniques = function(ar, eqfunc, discard_func, sortkey, tol) {
  eqfunc = eqfunc || ((x, y) => x === y);
  discard_func = discard_func || (() => false);
  if (sortkey) {
    ar.sort((a,b) => sortkey(a) - sortkey(b));
  }
  const newar = [];
  for (let i = 0; i < ar.length; ++i) {
    if (i == ar.length - 1) {
      newar.push(ar[i]);
      break;
    }
    let sameQ = false;
    for (let j = i+1; j < ar.length && (!sortkey || (sortkey(ar[j]) - sortkey(ar[i]) < tol)); ++j) {
      if (eqfunc(ar[i], ar[j])) {
        sameQ = true;
        break;
      }
    }
    if (!sameQ && !discard_func(ar[i])) { newar.push(ar[i]); }
  }
  return newar;
};

const uniqueaffineset =
    (Afset, discard_func, sortkey_func, tol) => uniques(Afset, (x, y) => x.sameAs(y, tol), discard_func, sortkey_func, tol);

const findclosure = function(Afset, recursion_limit, discard_func, sortkey_func, tol) {
  let uniqueset;
  recursion_limit = recursion_limit || 3;
  let oldset = Afset;
  let i = 0;
  const use_mag2 = !!(Afset[0].mag2);
  while (i < recursion_limit) {
    let smallset = Afset.slice();
    if (use_mag2) {
      smallset.sort((a,b) => a.mag2() - b.mag2());
      smallset = smallset.slice(0,40);
    }
    const setprod = affinesetproduct(Afset, smallset).concat(Afset);
    uniqueset = uniqueaffineset(setprod, discard_func, sortkey_func, tol);
    if (oldset.length === uniqueset.length) { break; }
    Afset = uniqueset;
    oldset = uniqueset;
    i++;
  }
  //console.log "/uniqueset() length: " + uniqueset.length
  return uniqueset;
};

// function maskfilter(Afset,positionfunc){
// }


// Affine Set Generators
//--------------------------------------------------------------------------------------------------

const makegrid = function(nx, ny, d) {
    const Afs = [];
    for (var i = 0; i < nx; i++) {
    for (var j = 0; j < ny; j++) {
        Afs.push(TranslationTransform(i * d, j * d));
    }
    }
    return Afs;
};

const rotateRosette = function(n, x, y) {
  const Afs = [];
  for (var i = 0; i < n; i++) {
    Afs.push(RotationAbout((i * 2*PI)/n, x, y));
  }
  return Afs;
};

const reflectRosette = function(n, x, y, offsetangle) {
  offsetangle = offsetangle || 0;
  const Afs = [];
  Afs.push(IdentityTransform());

  for (var i = 0; i< n; i++) {
    Afs.push(RotationAbout(offsetangle, x, y).multiply(ReflectionAbout((i * PI) / n, x, y)));
  }

  return findclosure(Afs);
};

export const IdentitySet = function(){
  const Afs = [];
  Afs.push(IdentityTransform());
  return Afs;
}

//XXX: rotation is badly broken (but visually interesting!!)
export const RosetteGroup = function(n1, n2, x, y, offsetangle) {
  offsetangle = offsetangle || 0;
  const Af1 = rotateRosette(n1, x, y);
  const Af2 = reflectRosette(n2, x, y);
  return findclosure(Af1.concat(affinesetproduct([ RotationAbout(offsetangle, x, y) ], Af2)));
};
/*
const multiRosette = function(n1, n2, x1, y1, x2, y2) {
  const Af1 = rotateRosette(n1, x1, y1);
  const Af2 = reflectRosette(n2, x2, y2);
  return affinesetproduct(Af1, Af2);
};

const multiRosette2 = function(n1, n2, x1, y1, x2, y2) {
  const Af1 = rotateRosette(n1, x1, y1);
  const Af2 = reflectRosette(n2, x2, y2);
  return affinesetproduct(Af2, Af1);
};

const multiRosette3 = function(n1, n2, n3, a, d, skew, x, y) {
  // n1 - core rosette degree,
  // n2 - rotation around '0,0' degree
  // n3 - number of scaled repeats outward to do
  // a  - scaling factor
  // skew - rotational skew factor for each scale jump
  // x,y - center coords of all this

  const Af1 = rotateRosette(n1, x, y);
  let afset = [];
  for (var i = 0; i < n2; i++) { // rotation loop
      for (let j = 1; j < n3; j++) {// scale loop
  const afop = RotationAbout(j*PI*skew, x, y).multiply( //skew
          TranslationTransform(j*d*cos((i*2*PI)/n2), j*d*sin((i*2*PI)/n2)).multiply(//translate
            RotationAbout((i * 2 * PI) / n2, x, y).multiply(//rotate to maintain rot sym.
              ScalingAbout(pow(a, j), x, y) )
          )
      ); //scale

      //add to the heap of fun
      afset = afset.concat(affinesetproduct([ afop ], Af1));
    }
  }

  return afset;
};
*/

// Generate Lattice
//--------------------------------------------------------------------------------------------------

export const generateLattice = function(spec, nx, ny, d, phi, x, y) {
  const transset = [];
  const { vec0 } = spec;
  const { vec1 } = spec;
  // Build set of translations
  for (var i = -floor(nx/2); i <= nx/2; i++) {
    for (var j = -floor(ny/2); j <= ny/2; j++) {
      transset.push(TranslationTransform(((i*vec0[0]) + (j*vec1[0]))*d,
 ((i*vec0[1]) + (j*vec1[1]))*d));
    }
  }
  // Lattice Rotation
  // - this is broken, gives a slight but noticeable pixel shift
  //   even w. phi=0, 0,0  wtf, roundoff error in similarity trafo?
  //const globalRot = RotationAbout(phi, 0, 0);
  //return transformAffineSet(globalRot, transset);
  return transset;
};


// Master Routine for making Wallpaper Group Sets
//--------------------------------------------------------------------------------------------------

export const generateTiling = function(spec, nx, ny, d, phi, x, y) {
  let rotset = [];
  let refset = [];
  let glideset = [];
  let Afset = [];
  let transset = [];
  const { rots } = spec;
  const { refs } = spec;
  const { vec0 } = spec;
  const { vec1 } = spec;

  rotset.push(IdentityTransform());
  // Add specified rotational symmetries
  if (spec.rots) {
    _.each(rots, r => rotset.push(RotationAbout(r[0], x + (d * r[1]), y + (d * r[2]))));
  //close the group if asked
    if (spec.closerot && (spec.closerot === true)) { rotset = findclosure(rotset); }
  }

  // Add specified reflection symmetries
  if (spec.refs) {
    _.each(refs, r => refset.push(ReflectionAbout(r[0], x + (d * r[1]), y + (d * r[2]))));
  //close the group if asked
    if (spec.closeref && (spec.closeref === true)) { refset = findclosure(refset); }
  }

  // HACK: for p4g: merge rotation and reflection ...there should be a better way
  if (spec.refrot && (spec.refrot === true)) {
    Afset = uniqueaffineset(affinesetproduct(rotset, refset).concat(rotset));
  } else {
    Afset = uniqueaffineset(rotset.concat(refset));
  }

  // Add specified glide symmetries
  if (spec.glides) {
    _.each(spec.glides, r => glideset.push(GlideTransform(r[0], d*r[1], x + (d*r[2]), y + (d*r[3]))));
  }
  Afset = Afset.concat(glideset);

  // Build set of translations
  for (var i = -floor(nx/2); i <= nx/2; i++){
      for (var j = -floor(ny/2); j <= ny/2; j++) {
        transset.push(TranslationTransform(((i*vec0[0]) + (j*vec1[0]))*d,
                                           ((i*vec0[1]) + (j*vec1[1]))*d));
    }
  }

  // Cartesian product of compositions
  const wholeset = affinesetproduct(transset, Afset);

  return wholeset;
};


// Wallpaper Symmetry Specification
//--------------------------------------------------------------------------------------------------

export const planarSymmetries = {
  squaregrid: {
    rots: [],
    refs: [],
    vec0: [ 0, 1 ],
    vec1: [ 1, 0 ],
    tile: [ 1, 1 ]
  },

  diagonalgrid: {
    rots: [],
    refs: [],
    vec0: [ sqrt(2)/2, sqrt(2)/2 ],
    vec1: [ sqrt(2)/2, -sqrt(2)/2 ],
    tile: [ sqrt(2), sqrt(2) ]
  },

  hexgrid: {
    rots: [],
    refs: [],
    vec0: [ sqrt(3), 0.0],
    vec1: [ sqrt(3)/2, -1.5 ],
    tile: [ sqrt(3), 3 ]
  },

  //  Wallpaper Groups ----------------------------
  //  rotation-free groups
  p1: {
    rots: [],
    refs: [],
    //vec0: [ sqrt(3)/2, 1.5 ]
    //vec1: [ sqrt(3), 0.0]
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  pm: {
    rots: [],
    refs: [ [ PI/2, 0, 0 ], [ PI/2, 0, -1/2 ] ],
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  cm: {
    rots: [],
    refs: [ [PI/2, 0, 0] ],
    //does PI/4 introduce unwanted additional symmetry?
    vec0: [ sin(PI/4), cos(PI/4) ],
    vec1: [ -sin(PI/4), cos(PI/4) ],
    tile: [ sqrt(2), sqrt(2) ]
  },

  pg: {
    rots: [],
    refs: [],
    //glides: [ [ PI/2, sqrt(3)/2, 0.0, 0.0] ],
    glides: [ [ PI/2, 1/2, 0.0, 0.0] ],
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  // 180deg rotation containing groups
  p2: {
    rots: [ [ PI, 0, 0 ] ],
    refs: [],
    //vec0: [ sqrt(3)/2, 1.5 ]
    //vec1: [ sqrt(3), 0.0]
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  pmg: {
    rots: [ [ PI, 1.0, 0.25 ] ],
    refs: [ [ 0.0, 0.0, 0.0] ],
    glides: [ [ PI / 2, 1/2.0, 1.0, 0.0] ],
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  pgg: {
    rots: [ [ PI, 0.0, 0.0] ],
    refs: [],
    glides: [ [ 0.0, 0.5, 0.0, 0.25 ], [ PI/2, 1/2, 1/4, 0.0 ] ],
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  pmm: {
    rots: [ [ PI, 0, 0 ] ],
    refs: [ [ 0, 0, 0 ], [ PI / 2, 0, 0 ] ],
    vec0: [ 1, 0 ],
    //vec1: [ 1.61803399, 0 ]
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  cmm: {
    rots: [ [PI, 0, 0] ],
    refs: [ [PI/2, 0, 0], [0, 0, 0] ],
    //does PI/4 introduce unwanted additional symmetry?
    vec0: [ sin(PI/4), cos(PI/4) ],
    vec1: [ -sin(PI/4), cos(PI/4) ],
    tile: [ sqrt(2), sqrt(2) ]
  },

  // Square-ish Groups
  p4: {
    rots: [ [ PI/2, 0, 0 ], [ PI, 0, 0 ], [ (3 * PI)/2, 0, 0 ] ],
    refs: [],
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  p4g: {
    rots: [ [ PI/2, 0, 0 ], [ PI, 0, 0 ], [ (3 * PI)/2, 0, 0 ] ],
    refs: [ [ -PI / 4, 1/2.0, 0.0] ],
    refrot: true,
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  p4m: {
    rots: [ [ PI/2, 0, 0 ], [ PI, 0, 0 ], [ (3 * PI)/2, 0, 0 ] ],
    refs: [ [ -PI / 4, 0.0, 0.0], [ PI / 4, 0.0, 0.0], [ 0.0, 0.0, 0.0] ],
    closeref: true,
    vec0: [ 1, 0 ],
    vec1: [ 0, -1 ],
    tile: [ 1, 1 ]
  },

  // Hex-ish Groups
  p3: {
    rots: [ [ (2 * PI) / 3, sqrt(3)/2, -0.5 ], [ (4 * PI) / 3, sqrt(3)/2, -0.5 ] ],
    refs: [],
    vec0: [ sqrt(3), 0.0],
    vec1: [ sqrt(3)/2, -1.5 ],
    tile: [ sqrt(3), 3 ]
  },

  p6: {
    rots: [ [ (2 * PI) / 3, sqrt(3)/2, -0.5 ], [ (4 * PI) / 3, sqrt(3)/2, -0.5 ],
            [ PI / 3.0, 0.0, 0.0], [ -PI / 3.0, 0.0, 0.0], [ (3 * PI) / 3.0, 0.0, 0.0] ],
    refs: [],
    vec0: [ sqrt(3), 0.0],
    vec1: [ sqrt(3)/2, -1.5 ],
    tile: [ sqrt(3), 3 ]
  },

  p31m: {
    rots: [ [ (2 * PI)/3, sqrt(3)/2, -0.5 ], [ (4*PI)/3, sqrt(3)/2, -0.5 ] ],
    refs: [ [ PI/3, 0.0, 0.0], [ -PI/3, 0.0, 0.0], [ 0.0, 0.0, 0.0] ],
    vec0: [ sqrt(3), 0.0],
    vec1: [ sqrt(3)/2, -1.5 ],
    tile: [ sqrt(3), 3 ]
  },

  p3m1: {
    rots: [ [ (2 * PI) / 3, 0, 0 ], [ (4 * PI) / 3, 0, 0 ] ],
    refs: [ [ -PI/2, 0, 0 ], [ ((-2 * PI)/3) - (PI/2), 0, 0 ], [ ((2 * PI)/3) - (PI/2), 0, 0 ] ],
    vec0: [ sqrt(3), 0.0],
    vec1: [ sqrt(3)/2, -1.5 ],
    tile: [ sqrt(3), 3 ]
  },

  p6m: {
    rots: [],
    refs: [ [ PI / 6, 0, 0 ], [ (2 * PI) / 6, 0, 0 ], [ (3 * PI) / 6, 0, 0 ],
            [ (4 * PI) / 6, 0, 0 ], [ (5 * PI) / 6, 0, 0 ], [ (6 * PI) / 6, 0, 0 ] ],
    closeref: true,
    vec0: [ sqrt(3), 0.0],
    vec1: [ sqrt(3)/2, -1.5 ],
    tile: [ sqrt(3), 3 ]
  }
};

// Isometry of the Poincare disk as a Mobius transform e^(i phi) (z + b)/(bbar z + 1)
// where b = x + i y
// --------------------------------------------------------------------------------------------------

class MobiusTransform {
  constructor(x, y, phi, flipfirst) {
    this.x = x || 0;
    this.y = y || 0;
    this.phi = phi || 0;
    this.flipfirst = flipfirst || false;
  }

  matrix() {
    // A = e^(i phi), B = e^(i phi) b, C = bbar, D = 1
    let A_re = Math.cos(this.phi);
    let A_im = Math.sin(this.phi);
    let B_re = A_re * this.x - A_im * this.y;
    let B_im = A_im * this.x + A_re * this.y;
    let C_re = this.x;
    let C_im = -this.y;
    let D_re = 1;
    let D_im = 0;
    return {A_re: A_re, A_im: A_im, B_re: B_re, B_im: B_im,
            C_re: C_re, C_im: C_im, D_re: D_re, D_im: D_im};
  }

  conjmatrix() {
    const m = this.matrix();
    return {A_re: m.A_re, A_im: -m.A_im, B_re: m.B_re, B_im: -m.B_im,
            C_re: m.C_re, C_im: -m.C_im, D_re: m.D_re, D_im: -m.D_im};
  }

  multiply(Min) {
    const m1 = Min.flipfirst ? this.matrix() : this.conjmatrix();
    const m2 = Min.matrix();
    const m3 = {A_re: m1.A_re * m2.A_re - m1.A_im * m2.A_im + m1.C_re * m2.B_re - m1.C_im * m2.B_im,
                A_im: m1.A_re * m2.A_im + m1.A_im * m2.A_re + m1.C_re * m2.B_im + m1.C_im * m2.B_re,
                B_re: m1.B_re * m2.A_re - m1.B_im * m2.A_im + m1.D_re * m2.B_re - m1.D_im * m2.B_im,
                B_im: m1.B_re * m2.A_im + m1.B_im * m2.A_re + m1.D_re * m2.B_im + m1.D_im * m2.B_re,
                C_re: m1.A_re * m2.C_re - m1.A_im * m2.C_im + m1.C_re * m2.D_re - m1.C_im * m2.D_im,
                C_im: m1.A_re * m2.C_im + m1.A_im * m2.C_re + m1.C_re * m2.D_im + m1.C_im * m2.D_re,
                D_re: m1.B_re * m2.C_re - m1.B_im * m2.C_im + m1.D_re * m2.D_re - m1.D_im * m2.D_im,
                D_im: m1.B_re * m2.C_im + m1.B_im * m2.C_re + m1.D_re * m2.D_im + m1.D_im * m2.D_re};
    const Dmag = m3.D_re * m3.D_re + m3.D_im * m3.D_im;
    const m3_norm = {A_re: (m3.A_re * m3.D_re + m3.A_im * m3.D_im) / Dmag,
                     A_im: (m3.A_im * m3.D_re - m3.A_re * m3.D_im) / Dmag,
                     B_re: (m3.B_re * m3.D_re + m3.B_im * m3.D_im) / Dmag,
                     B_im: (m3.B_im * m3.D_re - m3.B_re * m3.D_im) / Dmag,
                     C_re: (m3.C_re * m3.D_re + m3.C_im * m3.D_im) / Dmag,
                     C_im: (m3.C_im * m3.D_re - m3.C_re * m3.D_im) / Dmag};
    const Mout = new MobiusTransform();
    Mout.phi = Math.atan2(m3_norm.A_im, m3_norm.A_re);
    const Amag = m3_norm.A_re * m3_norm.A_re + m3_norm.A_im * m3_norm.A_im;
    Mout.x = (m3_norm.B_re * m3_norm.A_re + m3_norm.B_im * m3_norm.A_im) / Amag;
    Mout.y = (m3_norm.B_im * m3_norm.A_re - m3_norm.B_re * m3_norm.A_im) / Amag;
    Mout.flipfirst = this.flipfirst ? !Min.flipfirst : Min.flipfirst;
    return Mout;
  }

  Lmultiply(Min) {
    return Min.multiply(this);
  }

  on(x, y) {
    if (this.flipfirst) {
      y = -y;
    }
    const num_re = x + this.x;
    const num_im = y + this.y;
    const denom_re = this.x * x + this.y * y + 1;
    const denom_im = this.x * y - this.y * x;
    const denom_mag = denom_re * denom_re + denom_im * denom_im;
    const quot_re = (num_re * denom_re + num_im * denom_im) / denom_mag;
    const quot_im = (num_im * denom_re - num_re * denom_im) / denom_mag;
    const cosphi = Math.cos(this.phi);
    const sinphi = Math.sin(this.phi);
    const out_re = cosphi * quot_re - sinphi * quot_im;
    const out_im = cosphi * quot_im + sinphi * quot_re;
    return [ out_re, out_im ];
  }

  onVec(x) {
    return this.on(x[0], x[1]);
  }

  // approximate comparison function
  sameAs(Min, tol) {
    if (this.flipfirst !== Min.flipfirst) {
      return false;
    }
    tol = tol || 1e-8;
    let sum = 0;
    let phidiff = abs(this.phi - Min.phi);
    while (phidiff > Math.PI) { phidiff -= 2*Math.PI; }
    sum += abs(phidiff);
    sum += abs(this.x - Min.x);
    sum += abs(this.y - Min.y);
    if (sum < tol) {
      return true;
    } else {
      return false;
    }
  }

  mag2() {
    return this.x*this.x + this.y*this.y;
  }
}

//Mobius transform operating on screen pixel coordinates rather than absolute unit-disk coordinates
//--------------------------------------------------------------------------------------------------

class PixelMobiusTransform {
  constructor(M, scale, xoff, yoff) {
    this.M = M;
    this.scale = scale;
    this.xoff = xoff;
    this.yoff = yoff;
  }

  multiply(PMin) {
    const PMout = new PixelMobiusTransform();
    PMout.M = this.M.multiply(PMin.M);
    PMout.scale = this.scale; // TODO what if this and PMin have different values? can it happen?
    PMout.xoff = this.xoff;
    PMout.yoff = this.yoff;
    return PMout;
  }

  Lmultiply(Min) {
    const PMout = new PixelMobiusTransform();
    PMout.M = PMin.M.multiply(Min.M);
    PMout.scale = this.scale;
    PMout.xoff = this.xoff;
    PMout.yoff = this.yoff;
    return PMout;
  }

  on(px, py) {
    const x = (px - this.xoff)/this.scale;
    const y = (py - this.yoff)/this.scale;
    const vec = this.M.on(x, y);
    return [vec[0]*this.scale + this.xoff, vec[1]*this.scale + this.yoff];
  }

  onVec(x) {
    return this.on(x[0], x[1]);
  }

  // approximate comparison function
  sameAs(PMin, tol) {
    tol = tol || 1e-8;
    return this.M.sameAs(PMin.M, tol); // TODO what if scale, xoff, yoff different
  }

  mag2() {
    return this.M.mag2();
  }
}

//Common Hyperbolic Isometries
//--------------------------------------------------------------------------------------------------

export const HIdentityTransform =
  () => new MobiusTransform(0, 0, 0);

export const HRotationTransform =
  angle => new MobiusTransform(angle, 0, 0);

export const HTranslationTransform =
  (dx, dy) => new MobiusTransform(0, dx, dy);

export const HReflectionTransform =
  angle => new AffineTransform(Math.cos(2 * angle), Math.sin(2 * angle), Math.sin(2 * angle), -1 * Math.cos(2 * angle), 0, 0);

export const HRotationAbout =
  (angle, px, py) => TranslationTransform(px, py).multiply(RotationTransform(angle)).multiply(TranslationTransform(-px, -py));

export const HReflectionAbout =
  (angle, px, py) => TranslationTransform(px, py).multiply(ReflectionTransform(angle)).multiply(TranslationTransform(-px, -py));

export const HGlideTransform =
  (angle, distance, px, py) => ReflectionAbout(angle, px, py).multiply(TranslationTransform(distance * cos(angle), distance * sin(angle)));

//Master Routine for making hyperbolic symmetry groups
//--------------------------------------------------------------------------------------------------

export const generateHyperbolic = function(spec, nx, ny, d, phi, x, y) {
  const basisset = [new PixelMobiusTransform(new MobiusTransform(0, 0, 0, false), d, x, y)];
  spec.transforms.forEach(t => basisset.push(new PixelMobiusTransform(new MobiusTransform(t.x, t.y, t.phi, t.flipfirst), d, x, y)));
  const myset = findclosure(basisset, 10, (pmt => pmt.M.x * pmt.M.x + pmt.M.y * pmt.M.y > 0.98), pmt => pmt.M.x, 1e-8);
  return myset;
};

//Wallpaper Symmetry Specification
//--------------------------------------------------------------------------------------------------

export const hyperbolicSymmetries = {
  '*237': {
    transforms: [
      {x: 0, y: 0, phi: 2*PI/7, flipfirst: false},
      {x: 0.49697042539518089622118039244518853801751504740469543592, y: 0, phi: PI/7, flipfirst: false},
      {x: 0, y: 0, phi: 0, flipfirst: true}
    ]
  },

  '*238': {
    transforms: [
      {x: 0, y: 0, phi: PI/4, flipfirst: false},
      {x: 0.64359425290558262473544343741820980892420274244400765115, y: 0, phi: 0, flipfirst: false},
      {x: 0, y: 0, phi: 0, flipfirst: true}
    ]
  },

  '237': {
    transforms: [
      {x: 0, y: 0, phi: 2*PI/7, flipfirst: false},
      {x: 0.49697042539518089622118039244518853801751504740469543592, y: 0, phi: PI/7, flipfirst: false}
    ]
  },

  '238': {
    transforms: [
      {x: 0, y: 0, phi: PI/4, flipfirst: false},
      {x: 0.64359425290558262473544343741820980892420274244400765115, y: 0, phi: 0, flipfirst: false}
    ]
  },
};
