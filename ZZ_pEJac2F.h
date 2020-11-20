#ifndef NTL_ZZ_pEJAC2F__H
#define NTL_ZZ_pEJAC2F__H

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>

NTL_OPEN_NNS

class ZZ_pEJac2FInfoT {
  private:
    ZZ_pEJac2FInfoT() {}
  public:
    ZZ_pEX f;
    ZZ_pE f0, f1, f2, f3, f4;

    ZZ_pEJac2FInfoT(const ZZ_pEX& f);
    ZZ_pEJac2FInfoT(const ZZ_pE& f0, const ZZ_pE& f1, const ZZ_pE& f2,
	const ZZ_pE& f3, const ZZ_pE& f4);
};

extern ZZ_pEJac2FInfoT *ZZ_pEJac2FInfo;  // global variable for the curve

class ZZ_pEJac2F {
  private:
  public:
    // coordinates of the divisor
    //   D = < u2*x^2 + u1*x + u0, \sqrt{Sq}*(v1*x + v0) >
    ZZ_pE u0, u1, u2, v0, v1;
    ZZ_pE Sq; // the formal parameter
    
    // initialize the curve
    static void init(const ZZ_pEX& f);
    static void init(const ZZ_pE& f0, const ZZ_pE& f1, const ZZ_pE& f2,
	const ZZ_pE& f3, const ZZ_pE& f4);

    // constructors and destructor
    ZZ_pEJac2F() { u0 = 1; } 
    ZZ_pEJac2F(const ZZ_pE& sq) { u0 = 1; Sq = sq; }
    ZZ_pEJac2F(const ZZ_pE& U0, const ZZ_pE& U1, const ZZ_pE& U2,
	const ZZ_pE& V0, const ZZ_pE& V1, const ZZ_pE& sq) {
      u0 = U0; u1 = U1; u2 = U2;
      v0 = V0; v1 = V1;
      Sq = sq;
    }
    ~ZZ_pEJac2F() { }
    
    // copy
    ZZ_pEJac2F& operator=(const ZZ_pEJac2F& D) {
      u0 = D.u0; u1 = D.u1; u2 = D.u2;
      v0 = D.v0; v1 = D.v1;
      Sq = D.Sq;
      return *this;
    }
    
    // read-only access to the curve
    static const ZZ_pE& f0() { return ZZ_pEJac2FInfo->f0; }
    static const ZZ_pE& f1() { return ZZ_pEJac2FInfo->f1; }
    static const ZZ_pE& f2() { return ZZ_pEJac2FInfo->f2; }
    static const ZZ_pE& f3() { return ZZ_pEJac2FInfo->f3; }
    static const ZZ_pE& f4() { return ZZ_pEJac2FInfo->f4; }
    static const ZZ_pEX& f() { return ZZ_pEJac2FInfo->f; }
    
    // U and V polys
    ZZ_pEX u() const;
    ZZ_pEX v() const;
};

// basic
inline void clear(ZZ_pEJac2F& D, const ZZ_pE& Sq) {
  D.u0 = 1; D.u1 = 0; D.u2 = 0; D.v0 = 0; D.v1 = 0; D.Sq = Sq;
}

// arithmetic
void add(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
void duplicate(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1);
void sub(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
void mul(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1, const ZZ& k);
inline void mul(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1, const long k) {
  ZZ K;
  K = to_ZZ(k);
  mul(E, D1, K);
}
inline void mul(ZZ_pEJac2F& E, const long k, const ZZ_pEJac2F& D1){
  mul(E, D1, k);
}
inline void mul(ZZ_pEJac2F& E, const ZZ& k, const ZZ_pEJac2F& D1){
  mul(E, D1, k);
}
void negate(ZZ_pEJac2F& E, const ZZ_pEJac2F& D);
ZZ_pEJac2F operator+(const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
ZZ_pEJac2F operator-(const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
ZZ_pEJac2F& operator+=(ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
ZZ_pEJac2F& operator-=(ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
ZZ_pEJac2F operator-(const ZZ_pEJac2F& D);

// tests
int IsZero(const ZZ_pEJac2F& D);
int operator==(const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);
int operator!=(const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2);

// random
void random(ZZ_pEJac2F& D, const ZZ_pE& Sq);
ZZ_pEJac2F random_ZZ_pEJac2F(const ZZ_pE& Sq);

// I/O
void print(const ZZ_pEJac2F& D);


// Pseudo-divHandler
extern long ZZ_pE_FoundFactorOfModulus;
extern ZZ_pX* ZZ_pE_FactorOfModulus;


NTL_CLOSE_NNS

#endif
