#ifndef LMPMCT_H
#define LMPMCT_H

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"

NTL_CLIENT


// MACRO to enable CPU cycle counters. Use ATOMIC_CPUCOUNT for threaded functions.
// Used counters should be registered in cpuperf.cpp.
// DO NOT ENABLE THIS MANUALLY, use 'rebuild.sh -c'
/****  PERFORMANCE_COUNTING *****/
#ifdef PERFORMANCE_COUNTING
#include "cpuperformance.hpp"
  #define CPUCOUNT(s) cpu::update_performance_counter _cpuperf##s(perfcounters[s]);
  #define ATOMIC_CPUCOUNT(s) cpu::update_atomic_performance_counter _atomic_cpuperf##s(atomic_perfcounters[s]);
  extern std::vector<uint64_t> perfcounters;
  extern std::vector<std::atomic<uint64_t>> atomic_perfcounters;
  extern cpu::performance_counter_manager perfmanager;
#else
  #define CPUCOUNT(s)
  #define ATOMIC_CPUCOUNT(s)
#endif
void show_cpu_stats();



/***************************************************************************\
 *  Typedef for a divisor with alpha and beta
\***************************************************************************/

class DivAndTrack {
  public:
    ZZ_pJac2 D;
    ZZ alpha;
    ZZ beta;
    int b, hash;
    /* b=1 --> Wild Pt
     * b=0 --> Tame Pt
     * b=-1 --> Wild/Tame undefined
     */

    DivAndTrack() { }
    DivAndTrack(const DivAndTrack& dd) {
      D = dd.D;  alpha = dd.alpha;  beta = dd.beta;  b = -1; hash = 0;
    }

    ~DivAndTrack() { }

    inline DivAndTrack& operator=(const DivAndTrack& dd) {
      this->D = dd.D;
      this->alpha = dd.alpha;
      this->beta = dd.beta;
      this->b = dd.b;
      this->hash = dd.hash;
      return *this;
    }

};

/*typedef struct {
  DivAndTrack Pt[100];

} dic_DivAndTrack;*/

typedef Vec<DivAndTrack> vec_DivAndTrack;

/***************************************************************************\
 *  Typedef for the random walk data
\***************************************************************************/

typedef struct {
  double l1;
  double l2;
  int r;
  DivAndTrack **O;
  ZZ K;
  ZZ_pJac2 KBaseDivisor;
  ZZ_pJac2 mBaseDivisor;
  ZZ_pJac2 pp1mBaseDivisor;
  ZZ_pJac2 BaseDivisor;
  int pD;
  ZZ B1min, B1max;
  ZZ B2min, B2max;
  ZZ p, m, tot_jumps;
  DivAndTrack distP1, distP2;
  int b;
} RWData;

inline int indexFromKKpBR(int k, int kp, int b, int r) {
  return (k + r*kp + r*r*b);
}

/***************************************************************************\
 *   Hash function.
\***************************************************************************/
unsigned int hashvalue(const ZZ_pJac2& D) {
  //ATOMIC_CPUCOUNT(102);
  ZZ h = rep(D.u1) ^ rep(D.v0);
  ZZ hh = rep(D.u0) ^ rep(D.v1);
  h ^= hh << NumBits(h);
  unsigned int x;
  unsigned char *str = (unsigned char *) (&x);
  BytesFromZZ(str, h, 4);
  return x;
}

unsigned int hashvalue2(const ZZ_pJac2& D) {
  //ATOMIC_CPUCOUNT(103);
  ZZ h = rep(D.u0) ^ rep(D.u1);
  ZZ hh = rep(D.v0) ^ rep(D.v1);
  h ^= hh << NumBits(h);
  unsigned int x;
  unsigned char *str = (unsigned char *) (&x);
  BytesFromZZ(str, h, 4);
  return x;
}

/***************************************************************************\
 *   Distinguished property
\***************************************************************************/

// proba gives the number of bits that have to be 0 for being
// distinguished
inline int IsDistinguished(const ZZ_pJac2& D, const int proba) {
  return (hashvalue2(D) & (((1UL)<<proba) - 1)) == 0;
}

// Starting point of a new thread.  Bit b indicates Wild or Tame
// kangaroo.
void SelectRandomStartingPoint(DivAndTrack& D, const int b,
    const RWData &rwdata) {
  ATOMIC_CPUCOUNT(101);
  D.alpha = rwdata.B1min + RandomBnd(rwdata.B1max - rwdata.B1min);
  D.beta = rwdata.B2min + RandomBnd(rwdata.B2max - rwdata.B2min);
  D.D = (-D.alpha)*rwdata.pp1mBaseDivisor + D.beta*rwdata.mBaseDivisor;

  if (b)
    D.D += rwdata.KBaseDivisor;
}

void JumpOneStep(DivAndTrack& Dp, const DivAndTrack& D, const RWData &rwdata) {
  unsigned int hash = hashvalue(D.D);
  int r = rwdata.r;
  int k, kp, b;

  unsigned int il1 = 0;
  if (rwdata.l1 <= 2) {
    il1 = (unsigned int)(1/rwdata.l1);
  }

  // Select k, kp, b according to hash.
  //      hash >>= (rwdata.pD+1);  this is now useless
  b = hash & 1; hash >>= 1;
  kp = hash % r; hash >>= 4;

  if (rwdata.l1 <= 2) {
    if ( (hash % il1) == 0 ) {
      hash >>= 3;
      k = (hash % (rwdata.r - 1)) + 1;
    } else
      k = 0;
  } else {
    k = hash % r;
  }

//  cout << k << " " << kp << " " << b<< " "  << r << endl;
//  cout << indexFromKKpBR(k, kp, b, r) << endl;
  Dp.D = D.D + rwdata.O[indexFromKKpBR(k, kp, b, r)]->D;
  Dp.alpha = D.alpha + rwdata.O[indexFromKKpBR(k, kp, b, r)]->alpha;
  Dp.beta = D.beta + rwdata.O[indexFromKKpBR(k, kp, b, r)]->beta;
}

#endif
