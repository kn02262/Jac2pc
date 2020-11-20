#include "elemQuo.h"

NTL_CLIENT

/*******************************************************************
 * Tools To work in the following structure:
 *   ( GF(q) )[X] / (X^2 - s*X + p)
 * where s and p are elements of GF(q)==ZZ_pE
 *******************************************************************/

ZZ_pE * EelemQUO::s;
ZZ_pE * EelemQUO::p;

void EinitQUO() {
  EelemQUO::s = new ZZ_pE;
  EelemQUO::p = new ZZ_pE;
  *(EelemQUO::s) = 0;
  *(EelemQUO::p) = 0;
}

void EsetSandP(const ZZ_pE& s, const ZZ_pE& p) {
  *(EelemQUO::s) = s;
  *(EelemQUO::p) = p;
}

ostream& operator<<(ostream& o, const EelemQUO& X) {
  o << X.y0 << " " << X.y1;
  return o;
}

istream& operator>>(istream& i, EelemQUO& X) {
  i >> X.y0;
  i >> X.y1;
  return i;
}
  

EelemQUO add(const EelemQUO& Y, const EelemQUO& Z) {
  EelemQUO T;
  T.y0 = Y.y0 + Z.y0;
  T.y1 = Y.y1 + Z.y1;
  return T;
}

EelemQUO subtract(const EelemQUO& Y, const EelemQUO& Z) {
  EelemQUO T;
  T.y0 = Y.y0 - Z.y0;
  T.y1 = Y.y1 - Z.y1;
  return T;
}

// TODO: try to precompute a multiplier for p
EelemQUO mul(const EelemQUO& Y, const EelemQUO& Z) {
  EelemQUO T;
  ZZ_pE t0, t1, t2;

  t0 = Y.y0*Z.y0;
  t2 = Y.y1*Z.y1;
  t1 = (Y.y0 + Y.y1)*(Z.y0 + Z.y1) - t0 - t2;

  T.y0 = t0 - (t2 * (*EelemQUO::p));
  T.y1 = t1 + (t2 * (*EelemQUO::s));
  return T;
}

EelemQUO inv(const EelemQUO& Y) {
  EelemQUO T;
  ZZ_pE denom;

  denom = (sqr(Y.y1)*(*EelemQUO::p) + Y.y0*Y.y1*(*EelemQUO::s) + sqr(Y.y0));
  if (denom == 0)
    Error("EelemQUO is not invertible\n");
  inv(denom, denom);

  NTL_NNS negate(T.y1, Y.y1*denom);
  T.y0 = (Y.y1*(*EelemQUO::s) + Y.y0)*denom;
  return T;
}

EelemQUO operator+(const EelemQUO& Y, const EelemQUO& Z) { return add(Y, Z); }
EelemQUO operator-(const EelemQUO& Y, const EelemQUO& Z) { return subtract(Y, Z); }
EelemQUO operator*(const EelemQUO& Y, const EelemQUO& Z) { return mul(Y, Z); }

EelemQUO operator*(const EelemQUO& Y, const ZZ_p& x) {
  EelemQUO Z = Y;
  Z.y0 *= x;
  Z.y1 *= x;
  return Z;
}
EelemQUO operator*(const ZZ_p& x, const EelemQUO& Y) {
    EelemQUO Z = Y;
    Z.y0 *= x;
    Z.y1 *= x;
    return Z;
}


void evalAtX1(EelemQUO& Z, const ZZ_pX& f) {
  EelemQUO x;
  x.y0 = 0;
  x.y1 = 1;

  Z.y0 = to_ZZ_pE(coeff(f, deg(f)));
  Z.y1 = 0;
  for(int i = deg(f) - 1; i >= 0; --i) {
    Z = Z*x;
    Z.y0 = Z.y0 + coeff(f, i);
  }
}

void evalAtX2(EelemQUO& Z, const ZZ_pX& f) {
  EelemQUO x;
  x.y0 = *EelemQUO::s;
  x.y1 = -1;

  Z.y0 = to_ZZ_pE(coeff(f, deg(f)));
  Z.y1 = 0;
  for(int i = deg(f) - 1; i >= 0; --i) {
    Z = Z*x;
    Z.y0 = Z.y0 + coeff(f, i);
  }
}

EelemQUO to_EelemQUO(const ZZ_pX& x) {
  EelemQUO z;
  z.y0 = to_ZZ_pE(x);
  z.y1 = 0;
  return z;
}

/****************************************************************
 ***  END of EelemQUO
 ****************************************************************/


/*******************************************************************
 * Tools To work in the following structure:
 *   ( GF(q) [p] )[X] / (X^2 - s*X + p)
 * where s is an element of GF(q).
 *******************************************************************/



ZZ_p * elemQUO::s;
  
void initQUO() {
  elemQUO::s = new ZZ_p;
  *(elemQUO::s) = 0;
}
    
void setS(const ZZ_p& s) {
  *(elemQUO::s) = s;
}   

ostream& operator<<(ostream& o, const elemQUO& X) {
  o << "(" << X.y1 << ")*X + (" << X.y0 << ")";
  return o;
} 
  

elemQUO add(const elemQUO& Y, const elemQUO& Z) {
  elemQUO T;
  T.y0 = Y.y0 + Z.y0;
  T.y1 = Y.y1 + Z.y1;
  return T; 
} 

elemQUO subtract(const elemQUO& Y, const elemQUO& Z) {
  elemQUO T;
  T.y0 = Y.y0 - Z.y0;
  T.y1 = Y.y1 - Z.y1;
  return T;
}
  
  
elemQUO mul(const elemQUO& Y, const elemQUO& Z) {
  elemQUO T;
  ZZ_pX t0, t1, t2;
  
  t0 = Y.y0*Z.y0;
  t2 = Y.y1*Z.y1;
  t1 = (Y.y0 + Y.y1)*(Z.y0 + Z.y1) - t0 - t2;
  
  T.y0 = t0 - (t2 << 1);
  T.y1 = t1 + t2*(*elemQUO::s);
  return T;
}

elemQUO operator+(const elemQUO& Y, const elemQUO& Z) { return add(Y, Z); }
elemQUO operator-(const elemQUO& Y, const elemQUO& Z) { return subtract(Y, Z); }
elemQUO operator*(const elemQUO& Y, const elemQUO& Z) { return mul(Y, Z); }

elemQUO operator*(const elemQUO& Y, const ZZ_p& x) {
  elemQUO Z = Y;
  Z.y0 *= x;
  Z.y1 *= x;
  return Z;
}
elemQUO operator*(const ZZ_p& x, const elemQUO& Y) {
    elemQUO Z = Y;
    Z.y0 *= x;
    Z.y1 *= x;
    return Z;
}

/****************************************************************
 ***  END of elemQUO
 ****************************************************************/
