
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

use Math, GPU, Time;
use RandomNumberGenerators;

config param useCudaRng = false;

config const Nphotons = 10_000_000;   /* number of photons in simulation */
config const useGpu = true;
config const numGpuThreads = 10_000;

type RNG = if useCudaRng then CudaRng else Mc321Rng;
const numGpus = Locales.size*here.gpus.size;
const NphotonsPerGpuThread = Nphotons/numGpus/numGpuThreads;

param PI         = 3.1415926;
param LIGHTSPEED = 2.997925E10; /* in vacuo speed of light [cm/s] */
param ALIVE      = 1;		/* if photon not yet terminated */
param DEAD       = 0;		/* if photon is to be terminated */
param THRESHOLD  = 0.01;	/* used in roulette */
param CHANCE     = 0.1;		/* used in roulette */
param COS90D     = 1.0E-6;
param ONE_MINUS_COSZERO = 1.0E-12;
param mua = 1.0;           /* absorption coefficient [cm^-1] */
param mus = 0.0;           /* scattering coefficient [cm^-1] */
param g = 0.90;            /* anisotropy [-] */
param nt = 1.33;           /* tissue index of refraction */
param radial_size = 3.0;   /* maximum radial size */
param NR = 100;                    /* number of radial positions */
param dr = radial_size/NR;         /* radial bin size */
param albedo = mus/(mus+mua);    /* albedo of tissue */

record photon {
  /* Propagation parameters */
  var x, y, z:real;    /* photon position */
  var ux, uy, uz:real; /* photon trajectory as cosines */
  var W:real;          /* photon weight */
  var absorb:real;     /* weighted deposited in a step due to absorption */
  var photon_status: int(8);  /* flag = ALIVE=1 or DEAD=0 */ // TODO enum

  var x2Plusy2: real;
}

proc photon.init(ref rng) {
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
  const costheta = 2.0*rng.next() - 1.0;
  const sintheta = sqrt(1.0 - costheta*costheta);	/* sintheta is always positive */
  const psi = 2.0*PI*rng.next();
  ux = sintheta*cos(psi);
  uy = sintheta*sin(psi);
  uz = costheta;
}

inline proc ref photon.hop(ref rng) {
  var  rnd: real;        /* assigned random value 0-1 */
  do { rnd = rng.next(); } while rnd <= 0.0; /* yields 0 < rnd <= 1 */
  const s = -log(rnd)/(mua + mus);          /* Step size.  Note: log() is base e */
  x += s * ux;                        /* Update positions. */
  y += s * uy;
  z += s * uz;

  x2Plusy2 = x*x + y*y;
}

inline proc ref photon.drop() {
  absorb = W*(1 - albedo);      /* photon weight absorbed at this step */
  W -= absorb;                  /* decrement WEIGHT by amount absorbed */
}

inline proc ref photon.spherical() {
  const r = sqrt(x2Plusy2 + z*z);    /* current spherical radial position */
  const ir = min((r/dr): int, NR);
  return ir;
}

inline proc ref photon.cylindrical() {
  const r = sqrt(x2Plusy2);          /* current cylindrical radial position */
  const ir = min((r/dr): int, NR);
  return ir;
}

inline proc ref photon.planar() {
  const r = abs(z);                  /* current planar radial position */
  const ir = min((r/dr): int, NR);
  return ir;
}

inline proc ref photon.spin(ref rng) {
  var uxx, uyy, uzz:real;	/* temporary values used during SPIN */
  /* Sample for costheta */
  var rnd = rng.next();
  const costheta;
  if (g == 0.0) {
    costheta = 2.0*rnd - 1.0;
  }
  else {
    const temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
    costheta = (1.0 + g*g - temp*temp)/(2.0*g);
  }
  const sintheta = sqrt(1.0 - costheta*costheta); /* sqrt() is faster than sin(). */

  /* Sample psi. */
  const psi = 2.0*PI*rng.next();
  const cospsi = cos(psi);
  const sinpsi;
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

proc ref photon.update(ref rng) {
  if (W < THRESHOLD) {
    if (rng.next() <= CHANCE) then
      W /= CHANCE;
    else photon_status = DEAD;
  }
}
/* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
inline proc SIGN(x) do return (if x>=0 then 1 else -1);

proc main() {
  var t: stopwatch;

  var	Csph: [0..NR] real;  /* spherical   photon concentration CC[ir=0..100] */
  var	Ccyl: [0..NR] real;  /* cylindrical photon concentration CC[ir=0..100] */
  var	Cpla: [0..NR] real;  /* planar      photon concentration CC[ir=0..100] */

  t.start();

  coforall loc in Locales with (+ reduce Csph,
                                + reduce Ccyl,
                                + reduce Cpla) do on loc {
    coforall gpu in here.gpus with (+ reduce Csph,
                                    + reduce Ccyl,
                                    + reduce Cpla) do on gpu {
      @gpu.assertEligible
      foreach thread in 0..<numGpuThreads {
        var rng = new RNG(thread);

        for i_photon in 0..<NphotonsPerGpuThread {
          var p = new photon(rng);

          do {
            p.hop(rng);

            p.drop();

            /* DROP absorbed weight into bin */
            gpuAtomicAdd(Csph[p.spherical()], p.absorb);

            /* DROP absorbed weight into bin */
            gpuAtomicAdd(Ccyl[p.cylindrical()], p.absorb);

            /* DROP absorbed weight into bin */
            gpuAtomicAdd(Cpla[p.planar()], p.absorb);

            p.spin(rng);

            p.update(rng);
          }
          while (p.photon_status == ALIVE);

        } /* end RUN */
      }
    }
  }
  t.stop();

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

  writeln("Number of photons : ", Nphotons);
  writeln("MPhotons/s : ", Nphotons/t.elapsed()/1_000_000);
  writeln("Elapsed time(s) : ", t.elapsed());
} /* end of main */
