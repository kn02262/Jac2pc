#include "ZZ_pEJac2F.h"
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/new.h>
#include <assert.h>

using namespace std;

NTL_START_IMPL

ZZ_pEJac2FInfoT *ZZ_pEJac2FInfo = NULL;

ZZ_pEJac2FInfoT::ZZ_pEJac2FInfoT(const ZZ_pEX& f) {
  assert ((deg(f) == 5) && (LeadCoeff(f) == 1));
  this->f0 = coeff(f, 0);
  this->f1 = coeff(f, 1);
  this->f2 = coeff(f, 2);
  this->f3 = coeff(f, 3);
  this->f4 = coeff(f, 4);
  this->f = f;
}

ZZ_pEJac2FInfoT::ZZ_pEJac2FInfoT(const ZZ_pE& f0, const ZZ_pE& f1, const ZZ_pE& f2,
    const ZZ_pE& f3, const ZZ_pE& f4) {
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

void ZZ_pEJac2F::init(const ZZ_pEX& f) {
  ZZ_pEJac2FInfo = NTL_NEW_OP ZZ_pEJac2FInfoT(f);
}

void ZZ_pEJac2F::init(const ZZ_pE& f0, const ZZ_pE& f1, const ZZ_pE& f2,
    const ZZ_pE& f3, const ZZ_pE& f4) {
  ZZ_pEJac2FInfo = NTL_NEW_OP ZZ_pEJac2FInfoT(f0, f1, f2, f3, f4);
}

ZZ_pEX ZZ_pEJac2F::u() const {
  ZZ_pEX P;
  SetCoeff(P, 0, u0); SetCoeff(P, 1, u1); SetCoeff(P, 2, u2);
  return P;
}
    
ZZ_pEX ZZ_pEJac2F::v() const {
  ZZ_pEX P;
  SetCoeff(P, 0, v0); SetCoeff(P, 1, v1); 
  return P;
}

void negate(ZZ_pEJac2F& E, const ZZ_pEJac2F& D) {
  E.u0 = D.u0; E.u1 = D.u1; E.u2 = D.u2;
  E.v0 = -D.v0; E.v1 = -D.v1;
}

ZZ_pEJac2F operator+(const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2) {
  ZZ_pEJac2F E;
  add(E, D1, D2);
  return E;
}

ZZ_pEJac2F operator-(const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2) {
  ZZ_pEJac2F E;
  sub(E, D1, D2);
  return E;
}

ZZ_pEJac2F& operator+=(ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2) {
  add(D1, D1, D2);
  return D1;
}

ZZ_pEJac2F& operator-=(ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2) {
  sub(D1, D1, D2);
  return D1;
}

ZZ_pEJac2F operator-(const ZZ_pEJac2F& D) {
  ZZ_pEJac2F E;
  negate(E, D);
  return E;
}

int IsZero(const ZZ_pEJac2F& D) {
  return ((D.u2 == 0) && (D.u1 == 0) && (D.u0 == 1));
}

int operator==(const ZZ_pEJac2F& D, const ZZ_pEJac2F& E) {
  return ((D.u0 == E.u0) && (D.u1 == E.u1) && (D.u2 == E.u2) &&
      (D.v0 == E.v0) && (D.v1 == E.v1) && (D.Sq == E.Sq));
}
int operator!=(const ZZ_pEJac2F& D, const ZZ_pEJac2F& E) {
  return ((D.u0 != E.u0) || (D.u1 != E.u1) || (D.u2 != E.u2) ||
      (D.v0 != E.v0) || (D.v1 != E.v1) || (D.Sq != E.Sq));
}

void print(const ZZ_pEJac2F& D) {
  cout << "[ " << D.u2 << "*x^2 + " << D.u1 << "*x + " << D.u0 << ", ";
  cout << D.v1 << "*x + " << D.v0 << " ];  Sq = " << D.Sq << "\n";
}

ZZ_pEJac2F random_ZZ_pEJac2F(const ZZ_pE& Sq) {
  ZZ_pEJac2F D;
  random(D, Sq);
  return D;
}

static int hasRoot(ZZ_pE& root, const ZZ_pEX& f) {
  vec_pair_ZZ_pEX_long factors;
  factors = CanZass(f);

  for (int i = 0; i < factors.length(); ++i) {
    if (deg(factors[i].a) == 1) {
      assert (LeadCoeff(factors[i].a) == 1);
      root = -coeff(factors[i].a, 0);
      return 1;
    }
  }
  clear(root);
  return 0;
}

static int hasSqrt(ZZ_pE& sq, const ZZ_pE& x) {
  ZZ_pEX f;
  SetCoeff(f, 2);
  SetCoeff(f, 1, 0);
  SetCoeff(f, 0, -x);
  return hasRoot(sq, f);
}

void random(ZZ_pEJac2F& D, const ZZ_pE& Sq) {
  ZZ_pE x1, x2, y1, y2;
   
  do {
    x1 = random_ZZ_pE();
    y1 = eval(ZZ_pEJac2F::f(), x1);
  } while (!hasSqrt(y1, y1));

  do {
    x2 = random_ZZ_pE();
    y2 = eval(ZZ_pEJac2F::f(), x2);
  } while (!hasSqrt(y2, y2));

  D.u2 = 1;
  D.u1 = -(x1+x2);
  D.u0 = x1*x2;

  D.v1 = (y1-y2)/(x1-x2);
  D.v0 = y1 - D.v1*x1;

  ZZ_pE sqrtSq;
  int tt = hasSqrt(sqrtSq, Sq);
  assert (tt);
  
  D.v0 /= sqrtSq;
  D.v1 /= sqrtSq;
  D.Sq = Sq; 
}

/******************************************
 *  Gasp!!! There is no DivHandler for ZZ_pE in NTL.... 
 *  
 *  This is a dirty patch to work around this.
 *  It would be much simpler to patch NTL sources, but then we have
 *  to run a specific version of NTL for this program...
 *****************************************/

/****
 * If an inversion is impossible, due to the non-irreducibility of 
 * ZZ_pE::modulus(), then myinv returns 1 otherwise 0 is returned.
 *
 * If 1 is returned, then a factor of ZZ_pE::modulus() has been found.
 * The global variables (beurk)
 *   long ZZ_pE_FoundFactorOfModulus   is set to 1
 * and
 *   ZZ_pX ZZ_pE_FactorOfModulus   is set to the factor
 ****/

long ZZ_pE_FoundFactorOfModulus = 0;
ZZ_pX* ZZ_pE_FactorOfModulus = 0;

static int myinv(ZZ_pE& y, const ZZ_pE& x) {
  long test;

  if (x == 0)
    Error("myinv: division by 0\n");
  
  test = InvModStatus(y._ZZ_pE__rep, x._ZZ_pE__rep, ZZ_pE::modulus());
  if (test) {
    ZZ_pE_FoundFactorOfModulus = 1;
    if (ZZ_pE_FactorOfModulus == 0)
      ZZ_pE_FactorOfModulus = new ZZ_pX();
    *ZZ_pE_FactorOfModulus = y._ZZ_pE__rep;
    cerr << "Found a factor of degree " << deg(*ZZ_pE_FactorOfModulus) << endl;
    return 1;
  } else 
    return 0;
}

static void SetSize(vec_ZZ_pX& x, long n, long m)
{
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      x[i].rep.SetMaxLength(m);
}

static int myDivRem(ZZ_pEX& q, ZZ_pEX& r, const ZZ_pEX& a, const ZZ_pEX& b)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ_pE *bp;
   ZZ_pE *qp;
   ZZ_pX *xp;


   ZZ_pE LCInv, t;
   ZZ_pX s;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZ_pEX: division by zero");

   if (da < db) {
      r = a;
      clear(q);
      return 0;
   }

   ZZ_pEX lb;

   if (&q == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      long test = myinv(LCInv, bp[db]);
      if (test) return 1;
   }

   vec_ZZ_pX x;

   SetSize(x, da+1, 2*ZZ_pE::degree());

   for (i = 0; i <= da; i++) 
      x[i] = rep(a.rep[i]);

   xp = x.elts();

   dq = da - db;
   q.rep.SetLength(dq+1);
   qp = q.rep.elts();

   for (i = dq; i >= 0; i--) {
      conv(t, xp[i+db]);
      if (!LCIsOne)
         mul(t, t, LCInv);
      qp[i] = t;
      negate(t, t);

      for (j = db-1; j >= 0; j--) {
         mul(s, rep(t), rep(bp[j]));
         add(xp[i+j], xp[i+j], s);
      }
   }

   r.rep.SetLength(db);
   for (i = 0; i < db; i++)
      conv(r.rep[i], xp[i]);
   r.normalize();
   return 0;
}

static int myXGCD(ZZ_pEX& d, ZZ_pEX& s, ZZ_pEX& t, const ZZ_pEX& a, const ZZ_pEX& b)
{
   ZZ_pE z;


   if (IsZero(b)) {
      set(s);
      clear(t);
      d = a;
   }
   else if (IsZero(a)) {
      clear(s);
      set(t);
      d = b;
   }
   else {
      long e = max(deg(a), deg(b)) + 1;

      ZZ_pEX temp(INIT_SIZE, e), u(INIT_SIZE, e), v(INIT_SIZE, e), 
            u0(INIT_SIZE, e), v0(INIT_SIZE, e), 
            u1(INIT_SIZE, e), v1(INIT_SIZE, e), 
            u2(INIT_SIZE, e), v2(INIT_SIZE, e), q(INIT_SIZE, e);


      set(u1); clear(v1);
      clear(u2); set(v2);
      u = a; v = b;

      do {
         long test;
	 test = myDivRem(q, u, u, v);
	 if (test) return 1;
         swap(u, v);
         u0 = u2;
         v0 = v2;
         mul(temp, q, u2);
         sub(u2, u1, temp);
         mul(temp, q, v2);
         sub(v2, v1, temp);
         u1 = u0;
         v1 = v0;
      } while (!IsZero(v));

      d = u;
      s = u1;
      t = v1;
   }

   if (IsZero(d)) return 0;
   if (IsOne(LeadCoeff(d))) return 0;

   /* make gcd monic */

   long test = myinv(z, LeadCoeff(d));
   if (test) return 1;
   mul(d, d, z);
   mul(s, s, z);
   mul(t, t, z);
   return 0;
}


void mul(ZZ_pEJac2F& E, const ZZ_pEJac2F& D, const ZZ& k) {
  if (k == 0) {
    clear(E, D.Sq);
    return;
  } 
  
  ZZ_pEJac2F D1, D2;
  ZZ m;

  if (k < 0) {
    negate(D1, D);
    m = -k;
  } else {
    D1 = D;
    m = k;
  }
  clear(D2, D.Sq);

  while (m > 0) {
    if (IsOdd(m)) 
      add(D2, D1, D2);
    m >>= 1;
    duplicate(D1, D1);
  }
  E = D2;
}


static void addCantor(ZZ_pEJac2F& F, const ZZ_pEJac2F& D, const ZZ_pEJac2F& E) {
  ZZ_pEX a1, a2, b1, b2, f;
  ZZ_pE Sq;
  long test;

  a1 = D.u(); b1 = D.v();
  a2 = E.u(); b2 = E.v();
  Sq = D.Sq;  
  f = ZZ_pEJac2F::f();
  
  // Composition
  
  ZZ_pEX d0, h1, h2;
  ZZ_pEX d, l, h3;
  ZZ_pEX a, b;
  
  test = myXGCD(d0, h1, h2, a1, a2);
  if (test) {
    clear(F, Sq);
    return;
  }
  if (deg(d0) > 0) {
    ZZ_pEX s = b1 + b2;
    if (s == 0) {
      d = d0; l = 1; h3 = 0;
    } else {
      test = myXGCD(d, l, h3, d0, s);
      if (test) {
	clear(F, Sq);
	return;
      }
    }
    
    h1 *= l;
    h2 *= l;
    a = (a1*a2) / sqr(d);
    if (deg(a) < 1) {
      clear(F, Sq);
      return;
    }
    assert ((LeadCoeff(a) == 1) && (LeadCoeff(d) == 1));
    b = ( ((h1*a1*b2+h2*a2*b1+h3*(b1*b2*Sq+f)/Sq) / d) % a );
  } else {
    a = a1*a2;
    if (deg(a) < 1) {
      clear(F, Sq);
      return;
    }
    b = (h1*a1*b2+h2*a2*b1) % a;
  }

  // Reduction

  ZZ_pE aux;
  while (deg(a) > 2) {
    a = (f - Sq*sqr(b)) / a;
    test = myinv(aux, LeadCoeff(a));
    if (test) {
      clear(F, Sq);
      return;
    }
    a *= aux;
    b = (-b) % a;
  }
  if (deg(a) < 1) {
    clear(F, Sq);
    return;
  }
  
  F.u0 = coeff(a, 0);
  F.u1 = coeff(a, 1);
  F.u2 = coeff(a, 2);
  F.v0 = coeff(b, 0);
  F.v1 = coeff(b, 1);
  F.Sq = Sq;
}

static void addLange(ZZ_pEJac2F& F, const ZZ_pEJac2F& D, const ZZ_pEJac2F& E) {
  if ((D.u2 == 0) || (E.u2 == 0)) {
    addCantor(F, D, E);
    return;
  }

  ZZ_pE u11, u10, u21, u20, v11, v10, v21, v20;
  ZZ_pE Sq;
  Sq = D.Sq;
  u11 = D.u1; u10 = D.u0;
  u21 = E.u1; u20 = E.u0;
  v11 = D.v1; v10 = D.v0;
  v21 = E.v1; v20 = E.v0;
  
  ZZ_pE z1, z2, z3, r, inv1, inv0;
  ZZ_pE w0, w1, w2, w3, w4, w5, sp1, sp0, spp0;
  ZZ_pE lp2, lp1, lp0;
  ZZ_pE up0, up1, vp0, vp1;
  ZZ_pE aux;
  
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
    mul(w1, r, sp1); mul(w1, w1, Sq);
    long test = myinv(w1, w1); 
    if (test) {
      clear(F, D.Sq);
      return;
    }
    mul(w2, r, w1);
    sqr(w3, sp1); mul(w3, w3, w1); mul(w3, w3, Sq);
    mul(w4, r, w2);
    sqr(w5, w4); mul(w5, w5, Sq);
    mul(spp0, sp0, w2); mul(spp0, spp0, Sq);

    // Step 5
    add(lp2, u21, spp0);
    mul(lp1, u21, spp0); add(lp1, lp1, u20);
    mul(lp0, u20, spp0);

    // Step 6
    sub(aux, spp0, u11); sub(up0, spp0, z1); mul(up0, up0, aux);
    sub(aux, lp1, u10); add(up0, up0, aux);
    mul(aux, v21, w4); add(aux, aux, aux); mul(aux, aux, Sq);
    add(up0, up0, aux);
    add(aux, u21, u21); add(aux, z1, aux); sub(aux, aux, ZZ_pEJac2F::f4());
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
  F.Sq = Sq;
}

static void duplicateLange(ZZ_pEJac2F& E, const ZZ_pEJac2F& D) {
  if ((D.u2 == 0)) {
    addCantor(E, D, D);
    return;
  }

  ZZ_pE u1, u0, v1, v0, Sq;
  u1 = D.u1; u0 = D.u0;
  v1 = D.v1; v0 = D.v0;
  Sq = D.Sq;
  
  ZZ_pE vt1, vt0, r;
  ZZ_pE w0, w1, w2, w3, w4, w5;
  ZZ_pE inv0, inv1;
  ZZ_pE kp0, kp1, sp0, sp1, spp0;
  ZZ_pE lp2, lp1, lp0;
  ZZ_pE up0, up1, vp0, vp1;
  ZZ_pE f4u1;
  ZZ_pE aux;
 
  // Step 1
  add(vt1, v1, v1); add(vt0, v0, v0);

  // Step 2
  sqr(w0, v1); mul(w0, w0, Sq);
  sqr(w1, u1);
  add(aux, w0, w0); add(w2, aux, aux); 
  mul(w3, u1, vt1);
  mul(r, u0, w2); sub(aux, vt0, w3); mul(aux, aux, vt0);
  mul(aux, aux, Sq); add(r, r, aux);
  if (r == 0) {
    addCantor(E, D, D); 
    return;
  }
  
  // Step 3
  negate(inv1, vt1);
  sub(inv0, vt0, w3);

  // Step 4
  add(w3, ZZ_pEJac2F::f3(), w1);
  add(w4, u0, u0);
  mul(f4u1, ZZ_pEJac2F::f4(), u1); sub(aux, w1, f4u1); add(kp1, aux, aux);
  sub(aux, w3, w4); add(kp1, kp1, aux);
  sub(aux, f4u1, w3); add(aux, aux, w4); add(aux, aux, w4);
  mul(kp0, u1, aux);
  mul(aux, ZZ_pEJac2F::f4(), u0); add(aux, aux, aux); sub(kp0, kp0, aux);
  sub(aux, ZZ_pEJac2F::f2(), w0); add(kp0, kp0, aux);

  // Step 5
  mul(w0, kp0, inv0);
  mul(w1, kp1, inv1);
  add(aux, kp0, kp1); add(sp1, inv0, inv1); mul(sp1, sp1, aux);
  add(aux, u1, 1); mul(aux, aux, w1); add(aux, aux, w0);
  sub(sp1, sp1, aux);
  mul(aux, u0, w1); sub(sp0, w0, aux);

  if (sp1 != 0) {
    // Step 6
    mul(w1, r, sp1); mul(w1, w1, Sq);
    long test = myinv(w1, w1);
    if (test) {
      clear(E, D.Sq);
      return;
    }
    mul(w2, w1, r);
    sqr(w3, sp1); mul(w3, w3, w1); mul(w3, w3, Sq);
    mul(w4, r, w2); sqr(w5, w4); mul(w5, w5, Sq);
    mul(spp0, sp0, w2); mul(spp0, spp0, Sq);

    // Step 7
    add(lp2, u1, spp0);
    mul(lp1, u1, spp0); add(lp1, lp1, u0);
    mul(lp0, u0, spp0);
		
    // Step 8
    sqr(up0, spp0); mul(aux, w4, v1); add(aux, aux, aux);
    mul(aux, aux, Sq); add(up0, up0, aux);
    sub(aux, u1, ZZ_pEJac2F::f4()); add(aux, aux, u1); mul(aux, aux, w5);
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
  E.Sq = Sq;
}


void add(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2) {
  assert (D1.Sq == D2.Sq);
  if (D1 == D2)
    duplicate(E, D1);
  else
    addLange(E, D1, D2);
}
void sub(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1, const ZZ_pEJac2F& D2) {
  assert (D1.Sq == D2.Sq);
  negate(E, D2);
  add(E, D1, E);
}

void duplicate(ZZ_pEJac2F& E, const ZZ_pEJac2F& D1) {
  duplicateLange(E, D1);
}


NTL_END_IMPL


#ifdef TESTEJAC2F
#warning "Testing ZZ_pEJac2F"

NTL_CLIENT

int main(int argc, char** argv) {
  ZZ p = to_ZZ("9001");
  ZZ_p::init(p);

  ZZ_pX Pdef;
  SetCoeff(Pdef, 3);
  SetCoeff(Pdef, 2, to_ZZ_p(to_ZZ("8052")));
  SetCoeff(Pdef, 1, to_ZZ_p(to_ZZ("5140")));
  SetCoeff(Pdef, 0, to_ZZ_p(to_ZZ("5819")));

  ZZ_pE::init(Pdef);
  
  ZZ_pE f0, f1, f2, f3, f4;
  random(f0); random(f1); random(f2);
  random(f3); random(f4);
  cout << "f0 = " << f0 << endl;
  cout << "f1 = " << f1 << endl;
  cout << "f2 = " << f2 << endl;
  cout << "f3 = " << f3 << endl;
  cout << "f4 = " << f4 << endl;

  ZZ_pEJac2F::init(f0, f1, f2, f3, f4);

  ZZ_pE Sq;
  random(Sq); sqr(Sq, Sq);

  ZZ_pEJac2F D, E;
  random(D, Sq);
  cout << "D = ";
  print(D);

  mul(E, D, 81645302);
  cout << "E = ";
  print(E);

//  assert (IsZero(E));
  clear(E, Sq);

  long k = 10000;
  
  for (int i = 0; i < k; ++i) {
    E += D;
  }
  
  mul(D, k, D);
  assert (E == D);
  
  return 0;
}

#endif


