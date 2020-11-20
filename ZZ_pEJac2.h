#ifndef NTL_ZZ_pEJAC2__H
#define NTL_ZZ_pEJAC2__H

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>

NTL_OPEN_NNS

class ZZ_pEJac2InfoT {
  private:
    ZZ_pEJac2InfoT() {}
  public:
    ZZ_pEX f;
    ZZ_pE f0, f1, f2, f3, f4;

    ZZ_pEJac2InfoT(const ZZ_pEX& f);
    ZZ_pEJac2InfoT(const ZZ_pE& f0, const ZZ_pE& f1, const ZZ_pE& f2,
	const ZZ_pE& f3, const ZZ_pE& f4);
};

extern ZZ_pEJac2InfoT *ZZ_pEJac2Info;  // global variable for the curve

class ZZ_pEJac2 {
  private:
    
  public:
    // coordinates of the divisor
    ZZ_pE u0, u1, u2, v0, v1;
    
    // initialize the curve
    static void init(const ZZ_pEX& f);
    static void init(const ZZ_pE& f0, const ZZ_pE& f1, const ZZ_pE& f2,
	const ZZ_pE& f3, const ZZ_pE& f4);

    // constructors and destructor
    ZZ_pEJac2() { u0 = 1; }
    ZZ_pEJac2(const ZZ_pE& U0, const ZZ_pE& U1, const ZZ_pE& U2,
	const ZZ_pE& V0, const ZZ_pE& V1) {
      u0 = U0; u1 = U1; u2 = U2;
      v0 = V0; v1 = V1;
    }
    ~ZZ_pEJac2() { }
    
    // copy
    ZZ_pEJac2& operator=(const ZZ_pEJac2& D) {
      u0 = D.u0; u1 = D.u1; u2 = D.u2;
      v0 = D.v0; v1 = D.v1;
      return *this;
    }
    
    // read-only access to the curve
    static const ZZ_pE& f0() { return ZZ_pEJac2Info->f0; }
    static const ZZ_pE& f1() { return ZZ_pEJac2Info->f1; }
    static const ZZ_pE& f2() { return ZZ_pEJac2Info->f2; }
    static const ZZ_pE& f3() { return ZZ_pEJac2Info->f3; }
    static const ZZ_pE& f4() { return ZZ_pEJac2Info->f4; }
    static const ZZ_pEX& f() { return ZZ_pEJac2Info->f; }
    
    // U and V polys
    ZZ_pEX u() const;
    ZZ_pEX v() const;
};

// basic
inline void clear(ZZ_pEJac2& D) {
  D.u0 = 1; D.u1 = 0; D.u2 = 0; D.v0 = 0; D.v1 = 0;
}

// arithmetic
void add(ZZ_pEJac2& E, const ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
void duplicate(ZZ_pEJac2& E, const ZZ_pEJac2& D1);
void sub(ZZ_pEJac2& E, const ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
void mul(ZZ_pEJac2& E, const ZZ_pEJac2& D1, const ZZ& k);
inline void mul(ZZ_pEJac2& E, const ZZ_pEJac2& D1, const long k) {
  ZZ K;
  K = to_ZZ(k);
  mul(E, D1, K);
}
inline void mul(ZZ_pEJac2& E, const long k, const ZZ_pEJac2& D1){
  mul(E, D1, k);
}
inline void mul(ZZ_pEJac2& E, const ZZ& k, const ZZ_pEJac2& D1){
  mul(E, D1, k);
}
void negate(ZZ_pEJac2& E, const ZZ_pEJac2& D);
ZZ_pEJac2 operator+(const ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
ZZ_pEJac2 operator-(const ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
ZZ_pEJac2& operator+=(ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
ZZ_pEJac2& operator-=(ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
ZZ_pEJac2 operator-(const ZZ_pEJac2& D);

// tests
int IsZero(const ZZ_pEJac2& D);
int operator==(const ZZ_pEJac2& D1, const ZZ_pEJac2& D2);
int operator!=(const ZZ_pEJac2& D1, const ZZ_pEJac2& D2);

// random
void random(ZZ_pEJac2& D);
ZZ_pEJac2 random_ZZ_pEJac2();

// I/O
void print(const ZZ_pEJac2& D);

NTL_CLOSE_NNS

#endif
