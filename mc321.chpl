
/********************************************
 *  mc321.c    , in ANSI Standard C programing language
 *
 *  Monte Carlo simulation yielding spherical, cylindrical, and planar
 *    responses to an isotropic point source in an infinite homogeneous
 *    medium with no boundaries. This program is a minimal Monte Carlo
 *    program scoring photon distributions in spherical, cylindrical,
 *    and planar shells.
 *
 *  by Steven L. Jacques based on prior collaborative work
 *    with Lihong Wang, Scott Prahl, and Marleen Keijzer.
 *    partially funded by the NIH (R29-HL45045, 1991-1997) and
 *    the DOE (DE-FG05-91ER617226, DE-FG03-95ER61971, 1991-1999).
 *
 *  A published report illustrates use of the program:
 *    S. L. Jacques: "Light distributions from point, line, and plane
 *    sources for photochemical reactions and fluorescence in turbid
 *    biological tissues," Photochem. Photobiol. 67:23-32, 1998.
 *
 *  Trivial fixes to remove warnings SAP, 11/2017
 **********/

use Math;

param PI         = 3.1415926;
param LIGHTSPEED = 2.997925E10; /* in vacuo speed of light [cm/s] */
param ALIVE      = 1;		/* if photon not yet terminated */
param DEAD       = 0;		/* if photon is to be terminated */
param THRESHOLD  = 0.01;		/* used in roulette */
param CHANCE     = 0.1;		/* used in roulette */
param COS90D     = 1.0E-6;
/* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
param ONE_MINUS_COSZERO = 1.0E-12;
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
inline proc SIGN(x) do return (if x>=0 then 1 else 1);
inline proc InitRandomGen do return RandomGen(0, 1): real;
/* Initializes the seed for the random number generator. */
inline proc RandomNum do return RandomGen(1, 0): real;
/* Calls for a random number from the randum number generator. */

config const	Nphotons = 10_000_000;   /* number of photons in simulation */

const	mua = 1.0;        /* absorption coefficient [cm^-1] */
const	mus = 0.0;        /* scattering coefficient [cm^-1] */
const	g = 0.90;          /* anisotropy [-] */
const	nt = 1.33;         /* tissue index of refraction */
const	radial_size = 3.0;  /* maximum radial size */
/* IF NR IS ALTERED, THEN USER MUST ALSO ALTER THE ARRAY DECLARATION TO A SIZE = NR + 1. */
const NR: int(16) = 100;         /* number of radial positions */
const dr = radial_size/NR;         /* radial bin size */
const	albedo = mus/(mus+mua);     /* albedo of tissue */

record photon {
  /* Propagation parameters */
  var x, y, z:real;    /* photon position */
  var ux, uy, uz:real; /* photon trajectory as cosines */
  var s:real;          /* step sizes. s = -log(RND)/mus [cm] */
  var costheta:real;   /* cos(theta) */
  var sintheta:real;   /* sin(theta) */
  var cospsi:real;     /* cos(psi) */
  var sinpsi:real;     /* sin(psi) */
  var psi:real;        /* azimuthal angle */
  var W:real;          /* photon weight */
  var absorb:real;     /* weighted deposited in a step due to absorption */
  var photon_status: int;  /* flag = ALIVE=1 or DEAD=0 */ // TODO enum

  var	r: real;          /* radial position */
  var   ir: int(16);         /* index to radial position */
}

proc photon.init() {
  init this;
  /**** LAUNCH
    Initialize photon position and trajectory.
    Implements an isotropic point source.
   *****/
  W = 1.0;                    /* set photon weight to one */
  photon_status = ALIVE;      /* Launch an ALIVE photon */

  x = 0;                      /* Set photon position to origin. */
  y = 0;
  z = 0;

  /* Randomly set photon trajectory to yield an isotropic source. */
  costheta = 2.0*RandomNum - 1.0;
  sintheta = sqrt(1.0 - costheta*costheta);	/* sintheta is always positive */
  psi = 2.0*PI*RandomNum;
  ux = sintheta*cos(psi);
  uy = sintheta*sin(psi);
  uz = costheta;
}

proc ref photon.hop() {
  var  rnd: real;        /* assigned random value 0-1 */
  do { rnd = RandomNum; } while rnd <= 0.0; /* yields 0 < rnd <= 1 */
  s = -log(rnd)/(mua + mus);          /* Step size.  Note: log() is base e */
  x += s * ux;                        /* Update positions. */
  y += s * uy;
  z += s * uz;

}

proc ref photon.drop() {
  absorb = W*(1 - albedo);      /* photon weight absorbed at this step */
  W -= absorb;                  /* decrement WEIGHT by amount absorbed */
}

proc ref photon.spherical() {
  r = sqrt(x*x + y*y + z*z);    /* current spherical radial position */
  ir = (r/dr): int(16);           /* ir = index to spatial bin */
  if (ir >= NR) then ir = NR: int(16);        /* last bin is for overflow */
}

proc ref photon.cylindrical() {
  r = sqrt(x*x + y*y);          /* current cylindrical radial position */
  ir = (r/dr): int(16);           /* ir = index to spatial bin */
  if (ir >= NR) then ir = NR: int(16);        /* last bin is for overflow */
}

proc ref photon.planar() {
  r = abs(z);                  /* current planar radial position */
  ir = (r/dr): int(16);           /* ir = index to spatial bin */
  if (ir >= NR) then ir = NR: int(16);        /* last bin is for overflow */
}

proc ref photon.spin() {
  var uxx, uyy, uzz:real;	/* temporary values used during SPIN */
  /* Sample for costheta */
  var rnd = RandomNum;
  if (g == 0.0) {
    costheta = 2.0*rnd - 1.0;
  }
  else {
    const temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
    costheta = (1.0 + g*g - temp*temp)/(2.0*g);
  }
  sintheta = sqrt(1.0 - costheta*costheta); /* sqrt() is faster than sin(). */

  /* Sample psi. */
  psi = 2.0*PI*RandomNum;
  cospsi = cos(psi);
  if (psi < PI) then
    sinpsi = sqrt(1.0 - cospsi*cospsi);     /* sqrt() is faster than sin(). */
  else
    sinpsi = -sqrt(1.0 - cospsi*cospsi);

  /* New trajectory. */
  if (1 - abs(uz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
    uxx = sintheta * cospsi;
    uyy = sintheta * sinpsi;
    uzz = costheta * SIGN(uz);   /* SIGN() is faster than division. */
  }
  else {					/* usually use this option */
    const temp = sqrt(1.0 - uz * uz);
    uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
    uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
    uzz = -sintheta * cospsi * temp + uz * costheta;
  }

  /* Update trajectory */
  ux = uxx;
  uy = uyy;
  uz = uzz;
}

proc ref photon.update() {
  if (W < THRESHOLD) {
    if (RandomNum <= CHANCE) then
      W /= CHANCE;
    else photon_status = DEAD;
  }
}

proc main() {
  /* other variables */
  var	Csph: [0..100] real;  /* spherical   photon concentration CC[ir=0..100] */
  var	Ccyl: [0..100] real;  /* cylindrical photon concentration CC[ir=0..100] */
  var	Cpla: [0..100] real;  /* planar      photon concentration CC[ir=0..100] */

  /**** INITIALIZATIONS
   *****/
  InitRandomGen;

  /**** RUN
    Launch N photons, initializing each one before progation.
   *****/
  for i_photon in 0..<Nphotons {
    var p: photon;

    /* HOP_DROP_SPIN_CHECK
       Propagate one photon until it dies as determined by ROULETTE.
     *******/
    do {
      /**** HOP
        Take step to new position
        s = stepsize
        ux, uy, uz are cosines of current photon trajectory
       *****/
      p.hop();

      /**** DROP
        Drop photon weight (W) into local bin.
       *****/
      p.drop();

      // TODO atomic
      p.spherical();
      Csph[p.ir] += p.absorb;           /* DROP absorbed weight into bin */

      p.cylindrical();
      Ccyl[p.ir] += p.absorb;           /* DROP absorbed weight into bin */

      p.planar();
      Cpla[p.ir] += p.absorb;           /* DROP absorbed weight into bin */

      /**** SPIN
        Scatter photon into new trajectory defined by theta and psi.
        Theta is specified by cos(theta), which is determined
        based on the Henyey-Greenstein scattering function.
        Convert theta and psi into cosines ux, uy, uz.
       *****/
      p.spin();

      /**** CHECK ROULETTE
        If photon weight below THRESHOLD, then terminate photon using Roulette
        technique.  Photon has CHANCE probability of having its weight increased
        by factor of 1/CHANCE, and 1-CHANCE probability of terminating.  *****/
      p.update();
    } /* end STEP_CHECK_HOP_SPIN */
    while (p.photon_status == ALIVE);

    /* If photon dead, then launch new photon. */
  } /* end RUN */

  /* print header */
  writef("number of photons = %i\n", Nphotons);
  writef("bin size = %5.5dr [cm] \n", dr);
  writef("last row is overflow. Ignore.\n");

  /* print column titles */
  writef("r [cm] \t Fsph [1/cm2] \t Fcyl [1/cm2] \t Fpla [1/cm2]\n");

  /* print data:  radial position, fluence rates for 3D, 2D, 1D geometries */
  for ir in 0..NR {
    /* r = sqrt(1.0/3 - (ir+1) + (ir+1)*(ir+1))*dr; */
    const r = (ir + 0.5)*dr;
    var shellvolume = 4.0*PI*r*r*dr; /* per spherical shell */
    /* fluence in spherical shell */
    const Fsph = Csph[ir]/Nphotons/shellvolume/mua;
    shellvolume = 2.0*PI*r*dr;   /* per cm length of cylinder */
    /* fluence in cylindrical shell */
    const Fcyl = Ccyl[ir]/Nphotons/shellvolume/mua;
    shellvolume = dr;            /* per cm2 area of plane */
    /* fluence in planar shell */
    const Fpla =Cpla[ir]/Nphotons/shellvolume/mua;
    writef("%5.5dr \t %4.3er \t %4.3er \t %4.3er \n", r, Fsph, Fcyl, Fpla);
  }
} /* end of main */


/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
param MBIG = 1000000000;
param MSEED = 161803398;
param MZ = 0;
param FAC = 1.0E-9;

var i1, i2: int;
var ma: [0..55] int;   /* ma[0] is not used. */

proc RandomGen(Type, Seed) {
  var        mj, mk: int;
  var       i, ii: int(16);

  if (Type == 0) {              /* set seed. */
    mj = MSEED - (if Seed < 0 then -Seed else Seed);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for i in 1..54 {
      ii = ((21 * i) % 55): int(16);
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ) then
        mk += MBIG;
      mj = ma[ii];
    }
    for ii in 1..4 {
      for i in 1..55 {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ) then
          ma[i] += MBIG;
      }
    }
    i1 = 0;
    i2 = 31;
  } else if (Type == 1) {       /* get a number. */
    i1 += 1;
    if (i1 == 56) then
      i1 = 1;
    i2 += 1;
    if (i2 == 56) then
      i2 = 1;
    mj = ma[i1] - ma[i2];
    if (mj < MZ) then
      mj += MBIG;
    ma[i1] = mj;
    return (mj * FAC);
  }
  /*else if (Type == 2) {       [> get status. <]*/
    /*[>for (i = 0; i < 55; i++)<]*/
    /*for i in 0..<55 do*/
      /*Status[i] = ma[i + 1];*/
    /*Status[55] = i1;*/
    /*Status[56] = i2;*/
  /*} else if (Type == 3) {       [> restore status. <]*/
    /*for i in 0..<55 do*/
      /*ma[i + 1] = Status[i];*/
    /*i1 = Status[55];*/
    /*i2 = Status[56];*/
  /*}*/
  else
    writeln("Wrong parameter to RandomGen().");
  return 0;
}

