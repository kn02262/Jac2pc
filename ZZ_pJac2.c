#include "ZZ_pJac2.h"
#include <NTL/new.h>
#include <assert.h>

using namespace std;

NTL_START_IMPL

ZZ_pJac2InfoT *ZZ_pJac2Info = NULL;

ZZ_pJac2InfoT::ZZ_pJac2InfoT(const ZZ_pX& f) {
  assert ((deg(f) == 5) && (LeadCoeff(f) == 1));
  this->f0 = coeff(f, 0);
  this->f1 = coeff(f, 1);
  this->f2 = coeff(f, 2);
  this->f3 = coeff(f, 3);
  this->f4 = coeff(f, 4);
  this->f = f;
}

ZZ_pJac2InfoT::ZZ_pJac2InfoT(const ZZ_p& f0, const ZZ_p& f1, const ZZ_p& f2,
    const ZZ_p& f3, const ZZ_p& f4) {
  this->f0 = f0;
  this->f1 = f1;
  this->f2 = f2;
  this->f3 = f3;
  this->f4 = f4;
  SetCoeff(this->f, 5);
  SetCoeff(this->f, 0, f0);
  SetCoeff(this->f, 1, f1);
  SetCoeff(this->f, 2, f2);
  SetCoeff(this->f, 3, f3);
  SetCoeff(this->f, 4, f4);
}

void ZZ_pJac2::init(const ZZ_pX& f) {
  ZZ_pJac2Info = NTL_NEW_OP ZZ_pJac2InfoT(f);
}

void ZZ_pJac2::init(const ZZ_p& f0, const ZZ_p& f1, const ZZ_p& f2,
    const ZZ_p& f3, const ZZ_p& f4) {
  ZZ_pJac2Info = NTL_NEW_OP ZZ_pJac2InfoT(f0, f1, f2, f3, f4);
}

ZZ_pX ZZ_pJac2::u() const {
  ZZ_pX P;
  SetCoeff(P, 0, u0); SetCoeff(P, 1, u1); SetCoeff(P, 2, u2);
  return P;
}

ZZ_pX ZZ_pJac2::v() const {
  ZZ_pX P;
  SetCoeff(P, 0, v0); SetCoeff(P, 1, v1);
  return P;
}

void negate(ZZ_pJac2& E, const ZZ_pJac2& D) {
  E.u0 = D.u0; E.u1 = D.u1; E.u2 = D.u2;
  E.v0 = -D.v0; E.v1 = -D.v1;
}

ZZ_pJac2 operator+(const ZZ_pJac2& D1, const ZZ_pJac2& D2) {
  ZZ_pJac2 E;
  add(E, D1, D2);
  return E;
}

ZZ_pJac2 operator-(const ZZ_pJac2& D1, const ZZ_pJac2& D2) {
  ZZ_pJac2 E;
  sub(E, D1, D2);
  return E;
}

ZZ_pJac2& operator+=(ZZ_pJac2& D1, const ZZ_pJac2& D2) {
  add(D1, D1, D2);
  return D1;
}

ZZ_pJac2& operator-=(ZZ_pJac2& D1, const ZZ_pJac2& D2) {
  sub(D1, D1, D2);
  return D1;
}

ZZ_pJac2 operator-(const ZZ_pJac2& D) {
  ZZ_pJac2 E;
  negate(E, D);
  return E;
}

ZZ_pJac2 operator*(const ZZ_pJac2& D1, const ZZ& k) {
  ZZ_pJac2 E;
  mul(E, k, D1);
  return E;
}

ZZ_pJac2 operator*(const ZZ& k, const ZZ_pJac2& D1) {
  ZZ_pJac2 E;
  mul(E, k, D1);
  return E;
}

ZZ_pJac2 operator*(const long k, const ZZ_pJac2& D1) {
  ZZ_pJac2 E;
  mul(E, k, D1);
  return E;
}

ZZ_pJac2 operator*(const ZZ_pJac2& D1, const long k) {
  ZZ_pJac2 E;
  mul(E, k, D1);
  return E;
}


int IsZero(const ZZ_pJac2& D) {
  return ((D.u2 == 0) && (D.u1 == 0) && (D.u0 == 1));
}

int operator==(const ZZ_pJac2& D, const ZZ_pJac2& E) {
  return ((D.u0 == E.u0) && (D.u1 == E.u1) && (D.u2 == E.u2) &&
      (D.v0 == E.v0) && (D.v1 == E.v1));
}
int operator!=(const ZZ_pJac2& D, const ZZ_pJac2& E) {
  return ((D.u0 != E.u0) || (D.u1 != E.u1) || (D.u2 != E.u2) ||
      (D.v0 != E.v0) || (D.v1 != E.v1));
}

void print(const ZZ_pJac2& D) {
  cout << "[ " << D.u2 << "*x^2 + " << D.u1 << "*x + " << D.u0 << ",\n";
  cout << D.v1 << "*x + " << D.v0 << " ]\n";
}

ZZ_pJac2 random_ZZ_pJac2() {
  ZZ_pJac2 D;
  random(D);
  return D;
}

void random(ZZ_pJac2& D) {
  ZZ_p x1, x2, y1, y2;
  ZZ yy1, yy2;
  ZZ p = ZZ_p::modulus();

	  do {
		x1 = random_ZZ_p();
		yy1 = rep(eval(ZZ_pJac2::f(), x1));
	  } while (Jacobi(yy1, p) != 1);
	  y1 = to_ZZ_p(SqrRootMod(yy1, p));

	  do {
		x2 = random_ZZ_p();
		yy2 = rep(eval(ZZ_pJac2::f(), x2));
	  } while (Jacobi(yy2, p) != 1);
	  y2 = to_ZZ_p(SqrRootMod(yy2, p));

  D.u2 = 1;
  D.u1 = -(x1+x2);
  D.u0 = x1*x2;

  D.v1 = (y1-y2)/(x1-x2);
  D.v0 = y1 - D.v1*x1;
}

void mul(ZZ_pJac2& E, const ZZ_pJac2& D, const ZZ& k) {
  if (k == 0) {
    clear(E);
    return;
  }

  ZZ_pJac2 D1, D2;
  ZZ m;

  if (k < 0) {
    negate(D1, D);
    m = -k;
  } else {
    D1 = D;
    m = k;
  }

  while (m > 0) {
    if (IsOdd(m))
      add(D2, D1, D2);
    m >>= 1;
    duplicate(D1, D1);
  }
  E = D2;
}


void addCantor(ZZ_pJac2& F, const ZZ_pJac2& D, const ZZ_pJac2& E) {
  ZZ_pX a1, a2, b1, b2, f;

  a1 = D.u(); b1 = D.v();
  a2 = E.u(); b2 = E.v();
  f = ZZ_pJac2::f();

  // Composition

  ZZ_pX d0, h1, h2;
  ZZ_pX d, l, h3;
  ZZ_pX a, b;

  XGCD(d0, h1, h2, a1, a2);
  if (deg(d0) > 0) {
    ZZ_pX s = b1 + b2;
    if (s == 0) {
      d = d0; l = 1; h3 = 0;
    } else
      XGCD(d, l, h3, d0, s);

    h1 *= l;
    h2 *= l;
    a = (a1*a2) / sqr(d);
    if (deg(a) < 1) {
      clear(F);
      return;
    }
    b = ( ((h1*a1*b2+h2*a2*b1+h3*(b1*b2+f)) / d) % a );
  } else {
    a = a1*a2;
    if (deg(a) < 1) {
      clear(F);
      return;
    }
    b = (h1*a1*b2+h2*a2*b1) % a;
  }

  // Reduction

  while (deg(a) > 2) {
    a = (f - sqr(b)) / a;
    b = (-b) % a;
  }
  if (deg(a) < 1) {
    clear(F);
    return;
  }
  MakeMonic(a);

  F.u0 = coeff(a, 0);
  F.u1 = coeff(a, 1);
  F.u2 = coeff(a, 2);
  F.v0 = coeff(b, 0);
  F.v1 = coeff(b, 1);
}

void addLange(ZZ_pJac2& F, const ZZ_pJac2& D, const ZZ_pJac2& E) {
  if ((D.u2 == 0) || (E.u2 == 0)) {
    addCantor(F, D, E);
    return;
  }

  ZZ_p u11, u10, u21, u20, v11, v10, v21, v20;
  u11 = D.u1; u10 = D.u0;
  u21 = E.u1; u20 = E.u0;
  v11 = D.v1; v10 = D.v0;
  v21 = E.v1; v20 = E.v0;

  ZZ_p z1, z2, z3, r, inv1, inv0;
  ZZ_p w0, w1, w2, w3, w4, w5, sp1, sp0, spp0;
  ZZ_p lp2, lp1, lp0;
  ZZ_p up0, up1, vp0, vp1;
  ZZ_p aux;

  // Step 1
  sub(z1, u11, u21); sub(z2, u20, u10);
  mul(z3, z1, u11); add(z3, z3, z2);
  mul(r, z2, z3); sqr(aux, z1); mul(aux, aux, u10);
  add(r, r, aux);
  if (r == 0) {
    addCantor(F, D, E);
    return;
  }

  // Step 2
  inv1 = z1; inv0 = z3;

  // Step 3
  sub(w0, v10, v20); sub(w1, v11, v21);
  mul(w2, inv0, w0);
  mul(w3, inv1, w1);
  add(aux, w0, w1); add(sp1, inv0, inv1); mul(sp1, sp1, aux);
  add(aux, u11, 1); mul(aux, aux, w3);
  add(aux, aux, w2);
  sub(sp1, sp1, aux);
  mul(sp0, u10, w3);
  sub(sp0, w2, sp0);

  if (sp1 != 0) {
    // Step 4
    mul(w1, r, sp1);
    inv(w1, w1);
    mul(w2, r, w1);
    sqr(w3, sp1); mul(w3, w3, w1);
    mul(w4, r, w2);
    sqr(w5, w4);
    mul(spp0, sp0, w2);

    // Step 5
    add(lp2, u21, spp0);
    mul(lp1, u21, spp0); add(lp1, lp1, u20);
    mul(lp0, u20, spp0);

    // Step 6
    sub(aux, spp0, u11); sub(up0, spp0, z1); mul(up0, up0, aux);
    sub(aux, lp1, u10); add(up0, up0, aux);
    mul(aux, v21, w4); add(aux, aux, aux); add(up0, up0, aux);
    add(aux, u21, u21); add(aux, z1, aux); sub(aux, aux, ZZ_pJac2::f4());
    mul(aux, aux, w5);
    add(up0, up0, aux);
    add(up1, spp0, spp0); add(aux, z1, w5); sub(up1, up1, aux);

    // Step 7
    sub(w1, lp2, up1);
    mul(w2, up1, w1); sub(aux, up0, lp1); add(w2, w2, aux);
    mul(vp1, w2, w3); sub(vp1, vp1, v21);
    mul(w2, up0, w1); sub(w2, w2, lp0);
    mul(vp0, w2, w3); sub(vp0, vp0, v20);
  } else {
    addCantor(F, D, E);
    return;
  }

  F.u2 = 1;
  F.u1 = up1;
  F.u0 = up0;
  F.v1 = vp1;
  F.v0 = vp0;
}

void addLauter(ZZ_pJac2& F, const ZZ_pJac2& D, const ZZ_pJac2& E) { // Warning: requires a curve with f4=0

  if ((D.u2 == 0) || (E.u2 == 0)) {
    addCantor(F, D, E);
    return;
  }

  ZZ_p u1,u0,v1,v0,u1s,u01,u1d,u0d,v1d,v0d,u1ds,u01d;
  u1 = D.u1; u0 = D.u0;
  u1d = E.u1; u0d = E.u0;
  v1 = D.v1; v0 = D.v0;
  v1d = E.v1; v0d = E.v0;

  sqr(u1s, u1);
  mul(u01, u1, u0);
  sqr(u1ds, u1d);
  mul(u01d, u0d, u1d);

  ZZ_p uS, v0D, v1D, M1, M2, M3, M4;
  ZZ_p t1, t2, t3, t4, l2, l3, d, aux1, aux2;
  ZZ_p A, B, C, Dc, Ec, Cs;
  ZZ_p l3sq, Csq, Dcsq;
  ZZ_p u1dd, u0dd, uu1dd, uu0dd, v1dd, v0dd;

  // Todo: rewrite in ASM-style:
  //add(uS, u1, u1d); sub(v0D, v0, v0d); sub(v1D, v1, v1d); sub(M1, u1s, u0); sub(M1, M1, u1ds); sub(M1, M1, u0d); sub(M2, u01d, u01);
  //sub(M3, u1, u1d); sub(M4, u0d, u0); sub(aux1, M2, v0D); sub(aux2, v1D, M1); mul(t1, aux1, aux2); t2=(-v0D-M2)*(v1D+M1);

  uS=u1+u1d; v0D=v0-v0d; v1D=v1-v1d; M1=u1s-u0-u1ds+u0d; M2=u01d-u01;
  M3=u1-u1d; M4=u0d-u0; t1=(M2-v0D)*(v1D-M1); t2=(-v0D-M2)*(v1D+M1); //2M
  t3=(-v0D+M4)*(v1D-M3); t4=(-v0D-M4)*(v1D+M3); //2M
  l2=t1-t2; l3=t3-t4; d=t3+t4-t1-t2-2*(M2-M4)*(M1+M3); //1M
  inv(A, d*l3); mul(B, d, A); mul(C, d, B); mul(Dc, l2, B); //I + 4M
  sqr(l3sq, l3);
  sqr(Csq, C);
  sqr(Dcsq, Dc);
  Ec =l3sq*A; Cs =Csq; u1dd = 2*Dc-Cs-uS; //1M + 2S
  u0dd = Dcsq + C*(v1+v1d) -((u1dd-Cs)*uS+(u1s+u1ds))/2; //2M + 1S
  sqr(uu1dd, u1dd); mul(uu0dd, u1dd, u0dd);  v1dd = Dc*(u1-u1dd)+ uu1dd-u0dd-u1s+u0; //2M + 1S
  v0dd = Dc*(u0-u0dd) + uu0dd - u01; v1dd = -(Ec*v1dd + v1); v0dd = -(Ec*v0dd + v0); //3M

  F.u2 = 1;
  F.u1 = u1dd;
  F.u0 = u0dd;
  F.v1 = v1dd;
  F.v0 = v0dd;
}

static void duplicateLange(ZZ_pJac2& E, const ZZ_pJac2& D) {
  if ((D.u2 == 0)) {
    addCantor(E, D, D);
    return;
  }

  ZZ_p u1, u0, v1, v0;
  u1 = D.u1; u0 = D.u0;
  v1 = D.v1; v0 = D.v0;

  ZZ_p vt1, vt0, r;
  ZZ_p w0, w1, w2, w3, w4, w5;
  ZZ_p inv0, inv1;
  ZZ_p kp0, kp1, sp0, sp1, spp0;
  ZZ_p lp2, lp1, lp0;
  ZZ_p up0, up1, vp0, vp1;
  ZZ_p f4u1;
  ZZ_p aux;

  // Step 1
  add(vt1, v1, v1); add(vt0, v0, v0);

  // Step 2
  sqr(w0, v1); sqr(w1, u1);
  add(aux, w0, w0); add(w2, aux, aux);
  mul(w3, u1, vt1);
  mul(r, u0, w2); sub(aux, vt0, w3); mul(aux, aux, vt0);
  add(r, r, aux);
  if (r == 0) {
    addCantor(E, D, D);
    return;
  }

  // Step 3
  negate(inv1, vt1);
  sub(inv0, vt0, w3);

  // Step 4
  add(w3, ZZ_pJac2::f3(), w1);
  add(w4, u0, u0);
  mul(f4u1, ZZ_pJac2::f4(), u1); sub(aux, w1, f4u1); add(kp1, aux, aux);
  sub(aux, w3, w4); add(kp1, kp1, aux);
  sub(aux, f4u1, w3); add(aux, aux, w4); add(aux, aux, w4);
  mul(kp0, u1, aux);
  mul(aux, ZZ_pJac2::f4(), u0); add(aux, aux, aux); sub(kp0, kp0, aux);
  sub(aux, ZZ_pJac2::f2(), w0); add(kp0, kp0, aux);

  // Step 5
  mul(w0, kp0, inv0);
  mul(w1, kp1, inv1);
  add(aux, kp0, kp1); add(sp1, inv0, inv1); mul(sp1, sp1, aux);
  add(aux, u1, 1); mul(aux, aux, w1); add(aux, aux, w0);
  sub(sp1, sp1, aux);
  mul(aux, u0, w1); sub(sp0, w0, aux);

  if (sp1 != 0) {
    // Step 6
    mul(w1, r, sp1); inv(w1, w1);
    mul(w2, w1, r);
    sqr(w3, sp1); mul(w3, w3, w1);
    mul(w4, r, w2); sqr(w5, w4);
    mul(spp0, sp0, w2);

    // Step 7
    add(lp2, u1, spp0);
    mul(lp1, u1, spp0); add(lp1, lp1, u0);
    mul(lp0, u0, spp0);

    // Step 8
    sqr(up0, spp0); mul(aux, w4, v1); add(aux, aux, aux);
    add(up0, up0, aux);
    sub(aux, u1, ZZ_pJac2::f4()); add(aux, aux, u1); mul(aux, aux, w5);
    add(up0, up0, aux);
    add(up1, spp0, spp0); sub(up1, up1, w5);

    // Step 9
    sub(w1, lp2, up1);
    mul(w2, up1, w1); sub(aux, up0, lp1); add(w2, w2, aux);
    mul(vp1, w2, w3); sub(vp1, vp1, v1);
    mul(w2, up0, w1); sub(w2, w2, lp0);
    mul(vp0, w2, w3); sub(vp0, vp0, v0);
  } else {
    addCantor(E, D, D);
    return;
  }

  E.u2 = 1;
  E.u1 = up1;
  E.u0 = up0;
  E.v1 = vp1;
  E.v0 = vp0;
}


void add(ZZ_pJac2& E, const ZZ_pJac2& D1, const ZZ_pJac2& D2) {
  #if ADD_ALGORITHM == CANTOR
	addCantor(E, D1, D2);
  #elif ADD_ALGORITHM == LANGE
    addLange(E, D1, D2);
  #elif ADD_ALGORITHM == LAUTER
    addLauter(E, D1, D2);
  #endif

  /* DEBUGGING: */
  //cerr << "DEBUGGING START" << endl;
  //cerr << "Divisor D1: " << "x^2+" << D1.u1 << "*x+" << D1.u0 << ", " << D1.v1 << "*x + " << D1.v0 << endl;
  //cerr << "Divisor D2: " << "x^2+" << D2.u1 << "*x+" << D2.u0 << ", " << D2.v1 << "*x + " << D2.v0 << endl;
  // ZZ_pJac2 D1v, D2v, Ev;
  // // Copy divisors D1, D2
  // D1v.u2 = D1.u2; D1v.u1 = D1.u1; D1v.u0 = D1.u0; D1v.v1 = D1.v1; D1v.v0 = D1.v0;
  // D2v.u2 = D2.u2; D2v.u1 = D2.u1; D2v.u0 = D2.u0; D2v.v1 = D2.v1; D2v.v0 = D2.v0;
  // addLange(E, D1, D2);
  //cerr << "SUM1: " << "x^2+" << E.u1 << "*x+" << E.u0 << " " << E.v1 << "*x + " << E.v0 << endl;
  //addLauter(Ev, D1v, D2v);
  //cerr << "SUM2: " << "x^2+" << Ev.u1 << "*x+" << Ev.u0 << " " << Ev.v1 << "*x + " << Ev.v0 << endl;
  //cerr << endl;
  //assert(E.u0 == Ev.u0);
  //assert(E.u1 == Ev.u1);
  //assert(E.v0 == Ev.v0);
  //assert(E.v1 == Ev.v1);
}

void sub(ZZ_pJac2& E, const ZZ_pJac2& D1, const ZZ_pJac2& D2) {
  negate(E, D2);
  add(E, D1, E);
}

void duplicate(ZZ_pJac2& E, const ZZ_pJac2& D1) {
  duplicateLange(E, D1);
}


NTL_END_IMPL


#ifdef TESTJAC2
#warning "Testing ZZ_pJac2"

NTL_CLIENT

int main(int argc, char** argv) {
  ZZ p = to_ZZ("9001");
  ZZ_p::init(p);

  ZZ_p f0, f1, f2, f3, f4;
  f0 = to_ZZ_p(to_ZZ("123"));
  f1 = to_ZZ_p(to_ZZ("627"));
  f2 = to_ZZ_p(to_ZZ("136"));
  f3 = to_ZZ_p(to_ZZ("74"));
  f4 = to_ZZ_p(to_ZZ("976"));

  ZZ_pJac2::init(f0, f1, f2, f3, f4);

  ZZ_pJac2 D, E;
  random(D);
  cout << "D = ";
  print(D);

  mul(E, D, 81645302);
  cout << "E = ";
  print(E);

  assert (IsZero(E));


  long k = 100000;

  for (int i = 0; i < k; ++i) {
    E += D;
  }

  mul(D, k, D);
  assert (E == D);

  return 0;
}

#endif
