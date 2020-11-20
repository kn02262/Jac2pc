#include <NTL/ZZ_pX.h>
#include <math.h>
#include <assert.h>

NTL_CLIENT

/* compile with:

g++ sommes.c -I../../tmp/ntl-5.3/include ../../tmp/ntl-5.3/src/ntl.a -lgmp -lm

g++ sommes.c -lntl -lgmp -lm

*/


////////////////////////////////////////////////////////
//  Compute x = Log(1 + a), where a is a polynomial modulo X^n.
//  a must have a zero constant term.
//  
//  Use Log(1 + a) = int ( 1 / (1+a) )
//
///////////////////////////////////////////////////////

void LogTrunc(ZZ_pX& x, const ZZ_pX& a, const long n) {
  // Make sense only if char > n  and val(a) > 0
  if (ZZ_p::modulus() < n)
    Error("Characteristic is too small to compute a Log");
  if (ConstTerm(a) != 0)
    Error("LogTrunc : non-zero constant term");
  
  x = diff(a) * InvTrunc(1 + a, n);
  for (int i = 0; i < n; ++i) 
    SetCoeff(x, i, coeff(x, i)/(i+1));
  x <<= 1;
  trunc(x, x, n);
}

ZZ_pX LogTrunc(const ZZ_pX& a, const long n) {
  ZZ_pX x;
  LogTrunc(x, a, n);
  return x;
}


////////////////////////////////////////////////////////
//  Compute x = Exp(a) mod X^n
//
//  Newton algorithm (see Brent Kung)
///////////////////////////////////////////////////////

ZZ_pX ExpTrunc(const ZZ_pX& a, const long n) {
  // Make sense only if char > n  and val(a) > 0
  if (ZZ_p::modulus() < n)
    Error("Characteristic is too small to compute a Log");
  if (ConstTerm(a) != 0)
    Error("LogTrunc : non-zero constant term");

  long prec = 1;
  ZZ_pX R;
  ZZ_pX tmp;
  R = 1;
  do {
    prec <<= 1;
    tmp = 1 + a - LogTrunc(R - 1, prec);
    R = MulTrunc(R, tmp, prec);
  } while (prec < n);

  return trunc(R, n);
}

void ExpTrunc(ZZ_pX& x, const ZZ_pX& a, const long n) {
  x = ExpTrunc(a, n);
}


//-------------------------------------------
// Computes the product of [T-(q_i + r_j)], 
// with q_i in Roots(Q), r_j in Roots(R)
// TODO : check that Q,R non-zero and non-constant
//-------------------------------------------

void ResultantSumRoots(ZZ_pX& Res, const ZZ_pX& QQ, const ZZ_pX& RR) {
  long degr = deg(QQ)*deg(RR);
  ZZ_pX Q, R;  
  Q = QQ;  MakeMonic(Q);
  R = RR;  MakeMonic(R);

  // Everything is done with reciprocal polynomials
  reverse(Q, Q);
  reverse(R, R);

  long prec = degr + 3;

  // Rouillier's series
  ZZ_pX Qser, Rser;
  Qser = MulTrunc(diff(Q), InvTrunc(Q, prec), prec);
  Rser = MulTrunc(diff(R), InvTrunc(R, prec), prec);

  // Hadamard's product with Exp(u)
  ZZ_p n;
  n = 2;
  for (int i = 1; i < prec; ++i) {
    ZZ_p coe;
    coe = coeff(Qser, i);  SetCoeff(Qser, i, coe/n);
    coe = coeff(Rser, i);  SetCoeff(Rser, i, coe/n);
    n *= (i+2);
  }
  Qser <<= 1;
  Rser <<= 1;
 
  // Careful on the first term !
  SetCoeff(Qser, 0, -deg(QQ));
  SetCoeff(Rser, 0, -deg(RR));
  
  // Plain product
  ZZ_pX Sser;
  Sser = MulTrunc(Qser, Rser, prec);
 
  // Hadamard's "division" by Exp(u)
  n = 1;
  for (int i = 1; i < prec; ++i) {
    SetCoeff(Sser, i, -coeff(Sser, i)*n);
    n *= i;
  }
  SetCoeff(Sser, 0, 0);
  
  // Exponentiate and reverse
  Sser = ExpTrunc(Sser, prec);
  reverse(Sser, Sser, degr);
  
  Res = Sser;
}

#if 0

int main(int argc, char **argv) {
  ZZ p = to_ZZ(1009);
  ZZ_p::init(p);

  ZZ_pX f, g, h;
  
  random(f, 5);
  random(g, 5);
  cout << "f = " << f << endl;
  cout << "g = " << g << endl;

  ResultantSumRoots(h, f, g);

  cout << "h = " << h << endl;

  return 0;
  
}

#endif

