#ifndef NTL_ZZ_pJAC2F__H
#define NTL_ZZ_pJAC2F__H

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

NTL_OPEN_NNS

class ZZ_pJac2FInfoT {
  private:
    ZZ_pJac2FInfoT() {}
  public:
    ZZ_pX f;
    ZZ_p f0, f1, f2, f3, f4;

    ZZ_pJac2FInfoT(const ZZ_pX& f);
    ZZ_pJac2FInfoT(const ZZ_p& f0, const ZZ_p& f1, const ZZ_p& f2,
	const ZZ_p& f3, const ZZ_p& f4);
};

extern ZZ_pJac2FInfoT *ZZ_pJac2FInfo;  // global variable for the curve

class ZZ_pJac2F {
  private:
  public:
    // coordinates of the divisor
    //   D = < u2*x^2 + u1*x + u0, \sqrt{Sq}*(v1*x + v0) >
    ZZ_p u0, u1, u2, v0, v1;
    ZZ_p Sq; // the formal parameter
    
    // initialize the curve
    static void init(const ZZ_pX& f);
    static void init(const ZZ_p& f0, const ZZ_p& f1, const ZZ_p& f2,
	const ZZ_p& f3, const ZZ_p& f4);

    // constructors and destructor
    ZZ_pJac2F() { u0 = 1; } 
    ZZ_pJac2F(const ZZ_p& sq) { u0 = 1; Sq = sq; }
    ZZ_pJac2F(const ZZ_p& U0, const ZZ_p& U1, const ZZ_p& U2,
	const ZZ_p& V0, const ZZ_p& V1, const ZZ_p& sq) {
      u0 = U0; u1 = U1; u2 = U2;
      v0 = V0; v1 = V1;
      Sq = sq;
    }
    ~ZZ_pJac2F() { }
    
    // copy
    ZZ_pJac2F& operator=(const ZZ_pJac2F& D) {
      u0 = D.u0; u1 = D.u1; u2 = D.u2;
      v0 = D.v0; v1 = D.v1;
      Sq = D.Sq;
      return *this;
    }
    
    // read-only access to the curve
    static const ZZ_p& f0() { return ZZ_pJac2FInfo->f0; }
    static const ZZ_p& f1() { return ZZ_pJac2FInfo->f1; }
    static const ZZ_p& f2() { return ZZ_pJac2FInfo->f2; }
    static const ZZ_p& f3() { return ZZ_pJac2FInfo->f3; }
    static const ZZ_p& f4() { return ZZ_pJac2FInfo->f4; }
    static const ZZ_pX& f() { return ZZ_pJac2FInfo->f; }
    
    // U and V polys
    ZZ_pX u() const;
    ZZ_pX v() const;
};

// basic
inline void clear(ZZ_pJac2F& D, const ZZ_p& Sq) {
  D.u0 = 1; D.u1 = 0; D.u2 = 0; D.v0 = 0; D.v1 = 0; D.Sq = Sq;
}

// arithmetic
void add(ZZ_pJac2F& E, const ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
void duplicate(ZZ_pJac2F& E, const ZZ_pJac2F& D1);
void sub(ZZ_pJac2F& E, const ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
void mul(ZZ_pJac2F& E, const ZZ_pJac2F& D1, const ZZ& k);
inline void mul(ZZ_pJac2F& E, const ZZ_pJac2F& D1, const long k) {
  ZZ K;
  K = to_ZZ(k);
  mul(E, D1, K);
}
inline void mul(ZZ_pJac2F& E, const long k, const ZZ_pJac2F& D1){
  mul(E, D1, k);
}
inline void mul(ZZ_pJac2F& E, const ZZ& k, const ZZ_pJac2F& D1){
  mul(E, D1, k);
}
void negate(ZZ_pJac2F& E, const ZZ_pJac2F& D);
ZZ_pJac2F operator+(const ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
ZZ_pJac2F operator-(const ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
ZZ_pJac2F& operator+=(ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
ZZ_pJac2F& operator-=(ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
ZZ_pJac2F operator-(const ZZ_pJac2F& D);

// tests
int IsZero(const ZZ_pJac2F& D);
int operator==(const ZZ_pJac2F& D1, const ZZ_pJac2F& D2);
int operator!=(const ZZ_pJac2F& D1, const ZZ_pJac2F& D2);

// random
void random(ZZ_pJac2F& D, const ZZ_p& Sq);
ZZ_pJac2F random_ZZ_pJac2F(const ZZ_p& Sq);

// I/O
void print(const ZZ_pJac2F& D);

NTL_CLOSE_NNS

#endif
