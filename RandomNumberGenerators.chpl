module RandomNumberGenerators {
  require "curand_kernel.h";
  require "kernel_rng.h";
  extern type curandState_t;

  pragma "codegen for CPU and GPU"
  extern proc rng_init(seed, idx, ref state: curandState_t): void;

  pragma "codegen for CPU and GPU"
  extern proc rng_get(ref state: curandState_t): real;

  record Mc321Rng {
    var i1, i2: int;
    var ma: 56*int;

    proc init(seed) {
      init this;
      RandomGen(0, seed, i1, i2, ma);
    }

    inline proc ref next() {
      return RandomGen(1, 0, i1, i2, ma);
    }

  }

  record CudaRng {
    var rng: curandState_t;

    proc init(idx) {
      init this;
      rng_init(1, idx, rng);
    }

    inline proc ref next() {
      return rng_get(rng);
    }
  }


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
  private param MBIG = 1000000000;
  private param MSEED = 161803398;
  private param MZ = 0;
  private param FAC = 1.0E-9;

  // TODO this could be refactored into helpers to avoid Type
  proc RandomGen(param Type, Seed, ref i1: int, ref i2: int, ref ma: 56*int) {
    if Type != 0 && Type != 1 then compilerError("Wrong type");
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
    return 0;
  }

}
