#ifndef NTL_ZZ_pJAC2__H
#define NTL_ZZ_pJAC2__H

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#define CANTOR 0
#define LANGE 1
#define LAUTER 2
#define ADD_ALGORITHM LANGE

NTL_OPEN_NNS

class ZZ_pJac2InfoT {
  private:
    ZZ_pJac2InfoT() {}
  public:
    ZZ_pX f;
    ZZ_p f0, f1, f2, f3, f4;

    ZZ_pJac2InfoT(const ZZ_pX& f);
    ZZ_pJac2InfoT(const ZZ_p& f0, const ZZ_p& f1, const ZZ_p& f2,
	const ZZ_p& f3, const ZZ_p& f4);
};

extern ZZ_pJac2InfoT *ZZ_pJac2Info;  // global variable for the curve

class ZZ_pJac2 {
  public:
    // coordinates of the divisor
    ZZ_p u0, u1, u2, v0, v1;

    // initialize the curve
    static void init(const ZZ_pX& f);
    static void init(const ZZ_p& f0, const ZZ_p& f1, const ZZ_p& f2,
	const ZZ_p& f3, const ZZ_p& f4);

    // constructors and destructor
    ZZ_pJac2() { u0 = 1; }
    ZZ_pJac2(const ZZ_p& U0, const ZZ_p& U1, const ZZ_p& U2,
	const ZZ_p& V0, const ZZ_p& V1) {
      u0 = U0; u1 = U1; u2 = U2;
      v0 = V0; v1 = V1;
    }
    ~ZZ_pJac2() { }

    // copy
    ZZ_pJac2& operator=(const ZZ_pJac2& D) {
      u0 = D.u0; u1 = D.u1; u2 = D.u2;
      v0 = D.v0; v1 = D.v1;
      return *this;
    }

    // read-only access to the curve
    static const ZZ_p& f0() { return ZZ_pJac2Info->f0; }
    static const ZZ_p& f1() { return ZZ_pJac2Info->f1; }
    static const ZZ_p& f2() { return ZZ_pJac2Info->f2; }
    static const ZZ_p& f3() { return ZZ_pJac2Info->f3; }
    static const ZZ_p& f4() { return ZZ_pJac2Info->f4; }
    static const ZZ_pX& f() { return ZZ_pJac2Info->f; }

    // U and V polys
    ZZ_pX u() const;
    ZZ_pX v() const;
};

// basic
inline void clear(ZZ_pJac2& D) {
  D.u0 = 1; D.u1 = 0; D.u2 = 0; D.v0 = 0; D.v1 = 0;
}

// arithmetic
void add(ZZ_pJac2& E, const ZZ_pJac2& D1, const ZZ_pJac2& D2);
void addCantor(ZZ_pJac2& F, const ZZ_pJac2& D, const ZZ_pJac2& E);
void addLange(ZZ_pJac2& F, const ZZ_pJac2& D, const ZZ_pJac2& E);
void addLauter(ZZ_pJac2& F, const ZZ_pJac2& D, const ZZ_pJac2& E);
void duplicate(ZZ_pJac2& E, const ZZ_pJac2& D1);
void sub(ZZ_pJac2& E, const ZZ_pJac2& D1, const ZZ_pJac2& D2);
void mul(ZZ_pJac2& E, const ZZ_pJac2& D1, const ZZ& k);
inline void mul(ZZ_pJac2& E, const ZZ_pJac2& D1, const long k) {
  ZZ K;
  K = to_ZZ(k);
  mul(E, D1, K);
}
inline void mul(ZZ_pJac2& E, const long k, const ZZ_pJac2& D1){
  mul(E, D1, k);
}
inline void mul(ZZ_pJac2& E, const ZZ& k, const ZZ_pJac2& D1){
  mul(E, D1, k);
}
void negate(ZZ_pJac2& E, const ZZ_pJac2& D);
ZZ_pJac2 operator+(const ZZ_pJac2& D1, const ZZ_pJac2& D2);
ZZ_pJac2 operator-(const ZZ_pJac2& D1, const ZZ_pJac2& D2);
ZZ_pJac2& operator+=(ZZ_pJac2& D1, const ZZ_pJac2& D2);
ZZ_pJac2& operator-=(ZZ_pJac2& D1, const ZZ_pJac2& D2);
ZZ_pJac2 operator-(const ZZ_pJac2& D);
ZZ_pJac2 operator*(const ZZ_pJac2& D1, const ZZ& k);
ZZ_pJac2 operator*(const ZZ& k, const ZZ_pJac2& D1);
ZZ_pJac2 operator*(const long k, const ZZ_pJac2& D1);
ZZ_pJac2 operator*(const ZZ_pJac2& D1, const long k);



// tests
int IsZero(const ZZ_pJac2& D);
int operator==(const ZZ_pJac2& D1, const ZZ_pJac2& D2);
int operator!=(const ZZ_pJac2& D1, const ZZ_pJac2& D2);

// random
void random(ZZ_pJac2& D);
ZZ_pJac2 random_ZZ_pJac2();

// I/O
void print(const ZZ_pJac2& D);

NTL_CLOSE_NNS

#endif
