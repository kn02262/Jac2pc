#include <NTL/ZZ_pX.h>
#include <assert.h>

NTL_CLIENT

ZZ_pX PsiDivPol(const ZZ_pX& f, const long n) {
  static int firsttime = 1;
  static vec_ZZ_pX* PsiTable;
  static ZZ_pX* SaveF;

  if (firsttime) {
    PsiTable = new vec_ZZ_pX();
    SaveF = new ZZ_pX();
    *SaveF = f;
    firsttime = 0;
    append(*PsiTable, to_ZZ_pX(0));
    append(*PsiTable, to_ZZ_pX(0));
  }
  
  if (f != *SaveF) {
    cerr << "Warning: in PsiDivPol, f has changed. Erase PsiTable\n";
    *SaveF = f;
    PsiTable->SetLength(0);
    append(*PsiTable, to_ZZ_pX(0));
    append(*PsiTable, to_ZZ_pX(0));
  }
  
  assert (n >= 0);
  if ((n == 0) || (n == 1)) {
    return to_ZZ_pX(0);
  }
  ZZ_pX Psi;
  if ((PsiTable->length() > n) && ((*PsiTable)[n] != 0)) {
    Psi = (*PsiTable)[n];
    return Psi;
  }
  ZZ_p f1, f2, f3, f4, f5;
  f1 = coeff(f,4);
  f2 = coeff(f,3);
  f3 = coeff(f,2);
  f4 = coeff(f,1);
  f5 = coeff(f,0);
  ZZ_pX x, x2, x3, x4, p1, p2, p3, p4, p5, fp, fp2, fp3, fp4;
  SetX(x); sqr(x2, x); mul(x3, x, x2); sqr(x4, x2);
  
  if (n == 2) {
    Psi = 1;
  } else if (n == 3) {
    Psi = 4*f;
  } else if (n == 4) {
    p1 = -f2-4*f1*x-10*x2;
    p2 = (-3*f2*x2-5*x4-4*f1*x3-2*f3*x-f4)*(3*f2*x+f3+6*f1*x2+10*x3);
    p3 = power((-3*f2*x2-5*x4-4*f1*x3-2*f3*x-f4), 3);
    Psi = -2*(p3-4*f*(p2-2*p1*f));
  } else if (n == 5) {
    fp = -diff(f);
    p1 = f1+5*x;
    p2 = -f2-4*f1*x-10*x2;
    p3 = sqr(3*f2*x+f3+6*f1*x2+10*x3);
    p4 = (3*f2*x+f3+6*f1*x2+10*x3);
    Psi = -4*f*(-5*power(fp, 4)+8*f*(3*p4*sqr(fp)+2*f*(-p3-2*fp*p2+4*p1*f)));
  } else if (n == 6) {
    p1 = f;
    p2 = -diff(f);
    p3 = 3*f2*x+f3+6*f1*x2+10*x3;
    p4 = -(10*x2+4*f1*x+f2);
    p5 = f1+5*x;
    Psi = -4*(8*sqr(f)*p4-4*p2*f*p3+power(p2,3))*
      (128*power(f,4)+64*p3*p4*power(f,3)+64*power(f,3)*p2*p5-
	48*sqr(f*p2)*p4-48*sqr(f)*p2*sqr(p3)+40*power(p2,3)*p3*f-7*power(p2,5))
      -sqr(-32*sqr(f)*p4*p2-5*power(p2,4)+24*f*sqr(p2)*p3-16*sqr(p3*f)+64*p5*power(f,3));
  } else if (n == 7) {
    fp = -3*f2*x2-5*x4-4*f1*x3-2*f3*x-f4;
    fp2 = 3*f2*x+f3+6*f1*x2+10*x3;
    fp3 = -f2-4*f1*x-10*x2;
    fp4 = f1+5*x;
    ZZ_pX f2, f3, f4;
    sqr(f2, f); mul(f3, f, f2); sqr(f4, f2);
    Psi = -8*f * (256*f4*fp4*fp2-256*f4*fp+128*f4*sqr(fp3)-192*f3*fp4*sqr(fp)
	-384*fp2*fp3*f3*fp-64*f3*power(fp2,3)+160*fp3*f2*power(fp,3)+240*sqr(fp2)*f2*sqr(fp)
	-140*fp2*f*power(fp,4)+21*power(fp,6))
      * (64*f3*fp4-5*power(fp,4)-32*f2*fp*fp3+24*fp2*f*sqr(fp)-16*f2*sqr(fp2))
      -16*f * sqr(64*f3*fp4*fp-7*power(fp,5)+40*f*power(fp,3)*fp2-48*f2*sqr(fp)*fp3
	  -48*f2*sqr(fp2)*fp+128*f4+64*f3*fp2*fp3);
  } else {
    int r, s;
    r = n/2;
    s = n-r;
    if (r == s) {
      r--;
      s++;
    }
    if (r == s-1) {
      r--;
      s++;
    }
    
    ZZ_pX m11, m12, m13, m21, m22, m23, m31, m32, m33;

    m11 = PsiDivPol(f, s  ) * PsiDivPol(f, r-2);
    m12 = PsiDivPol(f, s+1) * PsiDivPol(f, r-1);
    m13 = PsiDivPol(f, s+2) * PsiDivPol(f, r  );

    m21 = PsiDivPol(f, s-1) * PsiDivPol(f, r-1);
    m22 = PsiDivPol(f, s  ) * PsiDivPol(f, r  );
    m23 = PsiDivPol(f, s+1) * PsiDivPol(f, r+1);

    m31 = PsiDivPol(f, s-2) * PsiDivPol(f, r  );
    m32 = PsiDivPol(f, s-1) * PsiDivPol(f, r+1);
    m33 = PsiDivPol(f, s  ) * PsiDivPol(f, r+2);

    ZZ_pX detM;
    detM = m11*(m22*m33 - m23*m32) + 
      m12*(-m21*m33 + m23*m31) + 
      m13*(m21*m32 - m22*m31);
    
    ZZ_pX Denom;
    Denom = PsiDivPol(f, s)*PsiDivPol(f, r)*PsiDivPol(f, s-r);
    Psi = -detM / Denom;
  }

  while (PsiTable->length() <= n)
    append(*PsiTable, to_ZZ_pX(0));
//  if (PsiTable->length() <= n)
//    (*PsiTable).SetLength(n+1);
  
  (*PsiTable)[n] = Psi;
  return Psi;
}

ZZ_pX GammaDivPol(const ZZ_pX& f, const int n) {
  static int firsttime = 1;
  static vec_ZZ_pX* GammaTable;
  static ZZ_pX* SaveF;

  if (firsttime) {
    GammaTable = new vec_ZZ_pX();
    SaveF = new ZZ_pX();
    *SaveF = f;
    firsttime = 0;
    append(*GammaTable, to_ZZ_pX(1));
    append(*GammaTable, to_ZZ_pX(1));
  }
  
  if (f != *SaveF) {
    cerr << "Warning: in PsiDivPol, f has changed. Erase PsiTable\n";
    *SaveF = f;
    GammaTable->SetLength(0);
    append(*GammaTable, to_ZZ_pX(1));
    append(*GammaTable, to_ZZ_pX(1));
  }
 
  assert (n >= 0);
  if ((n == 0) || (n == 1)) {
    return to_ZZ_pX(1);
  }
  ZZ_pX Gamma;
  if ((GammaTable->length() > n) && ((*GammaTable)[n] != 0)) {
    Gamma = (*GammaTable)[n];
    return Gamma;
  }

  ZZ_p f1, f2, f3, f4, f5;
  f1 = coeff(f,4);
  f2 = coeff(f,3);
  f3 = coeff(f,2);
  f4 = coeff(f,1);
  f5 = coeff(f,0);
  ZZ_pX x, x2, x3, x4, p1, p2, p3, p4, p5, fp, fp2, fp3, fp4;
  SetX(x); sqr(x2, x); mul(x3, x, x2); sqr(x4, x2);

  fp = -3*f2*x2-5*x4-4*f1*x3-2*f3*x-f4;
  fp2 = 3*f2*x+f3+6*f1*x2+10*x3;
  fp3 = -f2-4*f1*x-10*x2;
  fp4 = f1+5*x;
						   
  if (n == 2) {
    Gamma = 4*f;
  } else if (n == 3) {
    Gamma = -64*fp4*power(f, 3)+32*fp*fp3*power(f, 2)+16*power(fp2, 2)*power(f, 2)
      -24*power(fp, 2)*fp2*f+5*power(fp, 4);

  } else if (n == 4) {
    Gamma = 1024*power(f, 5)+512*power(f, 4)*fp*fp4+512*power(f, 4)*fp2*fp3
      -384*power(f, 3)*power(fp, 2)*fp3-384*power(f, 3)*fp*power(fp2, 2)+320*power(f, 2)*power(fp, 3)*fp2-56*f*power(fp, 5);

  } else if (n == 5) {
    Gamma =-4096*power(f, 6)*power(fp2, 2)+2048*power(f, 5)*power(fp, 2)*fp2+16384*power(f, 7)*fp4
      -256*power(f, 4)*power(fp, 4)-4096*power(f, 6)*power(fp3, 3)-14*power(fp, 9)-2560*power(f, 4)*power(fp3, 2)*power(fp, 3)
      -384*power(f, 2)*fp3*power(fp, 6)+512*power(f, 4)*fp*power(fp2, 4)+512*power(f, 3)*power(fp, 3)*power(fp2, 3)
      -576*power(f, 2)*power(fp, 5)*power(fp2, 2)+160*f*power(fp, 7)*fp2-768*power(f, 3)*power(fp, 5)*fp4
      +8192*power(f, 6)*power(fp4, 2)*fp-4096*power(f, 5)*fp3*power(fp, 2)*fp4+10240*power(f, 5)*power(fp3, 2)*fp*fp2
      -6144*power(f, 4)*fp3*power(fp, 2)*power(fp2, 2)+3072*power(f, 3)*fp3*power(fp, 4)*fp2
      -4096*power(f, 5)*fp*power(fp2, 2)*fp4+4096*power(f, 4)*power(fp, 3)*fp2*fp4;

  } else if (n == 6) {
    Gamma = 524288*power(f, 9)*fp-262144*power(f, 9)*power(fp3, 2)-40960*power(f, 5)*power(fp, 6)
      -288*f*power(fp, 11)-327680*power(f, 7)*power(fp, 3)*fp3-393216*power(f, 7)*power(fp, 2)*power(fp2, 2)
      +262144*power(f, 8)*power(fp, 2)*fp4+262144*power(f, 6)*power(fp, 4)*fp2-131072*power(f, 7)*power(fp4, 2)*power(fp, 3)
      -524288*power(f, 9)*power(fp4, 2)*fp3+16384*power(f, 4)*fp4*power(fp, 7)-36864*power(f, 5)*power(fp, 5)*power(fp3, 2)
      -98304*power(f, 7)*power(fp, 2)*power(fp3, 3)-5248*power(f, 3)*power(fp, 8)*fp3+36864*power(f, 4)*power(fp2, 3)*power(fp, 5)
      -40960*power(f, 5)*power(fp2, 4)*power(fp, 3)+32768*power(f, 6)*power(fp2, 5)*fp-32768*power(f, 7)*power(fp2, 4)*fp3
      -17408*power(f, 3)*power(fp2, 2)*power(fp, 7)+3712*power(f, 2)*power(fp, 9)*fp2-131072*power(f, 8)*fp2*power(fp3, 3)
      +786432*power(f, 8)*fp*fp2*fp3-147456*power(f, 5)*fp4*power(fp, 5)*fp2
      +180224*power(f, 6)*fp4*power(fp, 4)*fp3+393216*power(f, 6)*fp4*power(fp, 3)*power(fp2, 2)
      -786432*power(f, 7)*fp4*power(fp, 2)*fp2*fp3-262144*power(f, 7)*fp4*power(fp2, 3)*fp
      +524288*power(f, 8)*power(fp4, 2)*fp*fp2+524288*power(f, 8)*fp4*power(fp3, 2)*fp
      +262144*power(f, 8)*fp4*fp3*power(fp2, 2)+38912*power(f, 4)*power(fp, 6)*fp3*fp2
      -61440*power(f, 5)*power(fp, 4)*fp3*power(fp2, 2)+98304*power(f, 6)*power(fp, 3)*power(fp3, 2)*fp2
      -32768*power(f, 6)*power(fp, 2)*fp3*power(fp2, 3)+196608*power(f, 7)*fp*power(fp3, 2)*power(fp2, 2);

  } else {
    ZZ_pX P;
    P = PsiDivPol(f, n+1)*GammaDivPol(f, n-1)
      + PsiDivPol(f, n-1)*PsiDivPol(f, n+2);
    Gamma = P / PsiDivPol(f, n);
  }

//  if (GammaTable->length() <= n)
//    (*GammaTable).SetLength(n+1);
  while (GammaTable->length() <= n)
    append(*GammaTable, to_ZZ_pX(0));
  
  (*GammaTable)[n] = Gamma;
  return Gamma;
}

ZZ_pX AlphaDivPol(const ZZ_pX& f, const int n) {
  static int firsttime = 1;
  static vec_ZZ_pX* AlphaTable;
  static ZZ_pX* SaveF;

  if (firsttime) {
    AlphaTable = new vec_ZZ_pX();
    SaveF = new ZZ_pX();
    *SaveF = f;
    firsttime = 0;
    append(*AlphaTable, to_ZZ_pX(-2));
    append(*AlphaTable, to_ZZ_pX(-2));
  }
  
  if (f != *SaveF) {
    cerr << "Warning: in PsiDivPol, f has changed. Erase PsiTable\n";
    *SaveF = f;
    AlphaTable->SetLength(0);
    append(*AlphaTable, to_ZZ_pX(-2));
    append(*AlphaTable, to_ZZ_pX(-2));
  }
 
  assert (n >= 0);
  if ((n == 0) || (n == 1)) {
    return to_ZZ_pX(-2);
  }
  ZZ_pX Alpha;
  if ((AlphaTable->length() > n) && ((*AlphaTable)[n] != 0)) {
    Alpha = (*AlphaTable)[n];
    return Alpha;
  }

  ZZ_p f1, f2, f3, f4, f5;
  f1 = coeff(f,4);
  f2 = coeff(f,3);
  f3 = coeff(f,2);
  f4 = coeff(f,1);
  f5 = coeff(f,0);
  ZZ_pX x, x2, x3, x4, p1, p2, p3, p4, p5, fp, fp2, fp3, fp4;
  SetX(x); sqr(x2, x); mul(x3, x, x2); sqr(x4, x2);

  fp = -3*f2*x2-5*x4-4*f1*x3-2*f3*x-f4;
  fp2 = 3*f2*x+f3+6*f1*x2+10*x3;
  fp3 = -f2-4*f1*x-10*x2;
  fp4 = f1+5*x;
  if (n == 2) {
    Alpha = 0;
  } else if (n == 3) {
    Alpha = -2*fp;
  } else if (n == 4) {
    Alpha = -8*f*fp;
  } else if (n == 5) {
    Alpha = 64*fp*fp3*power(f, 2)-40*power(fp, 2)*fp2*f+9*power(fp, 4)-64*fp4*power(f, 3)
      +16*power(fp2, 2)*power(f, 2);
  } else if (n == 6) {
    Alpha = 1024*power(f, 4)*fp*fp4-640*power(f, 3)*power(fp, 2)*fp3-512*power(f, 3)*fp*power(fp2, 2)
      +512*power(f, 2)*power(fp, 3)*fp2-96*f*power(fp, 5)+1024*power(f, 5)+512*power(f, 4)*fp2*fp3;
  } else {
    ZZ_pX P;
    P = PsiDivPol(f, n-1)*AlphaDivPol(f, n-1)
      + PsiDivPol(f, n)*PsiDivPol(f, n-3);
    Alpha = P / PsiDivPol(f, n-2);
  }

//  if (AlphaTable->length() <= n)
//    (*AlphaTable).SetLength(n+1);
  while (AlphaTable->length() <= n)
    append(*AlphaTable, to_ZZ_pX(0));

  (*AlphaTable)[n] = Alpha;
  return Alpha;
}

void UnnormalizedDeltaDivPol(ZZ_pX& d0, ZZ_pX& d1, ZZ_pX& d2,
    const ZZ_pX& f, const int n) {
  d0 = -PsiDivPol(f, n-1)*PsiDivPol(f, n+1);
  d1 = -PsiDivPol(f, n-1)*GammaDivPol(f, n)
    + PsiDivPol(f, n+1)*AlphaDivPol(f, n);
  d2 = -16*sqr(f*PsiDivPol(f, n));
}

void UnnormalizedEpsilonDivPol(ZZ_pX& eps0, ZZ_pX& eps1, ZZ_pX& denom,
    const ZZ_pX& f, const int n) {
  ZZ_pX d0, d1, d2;
  ZZ_pX d0_prec, d1_prec, d2_prec;
  ZZ_pX d0_next, d1_next, d2_next;

  UnnormalizedDeltaDivPol(d0, d1, d2, f, n);
  UnnormalizedDeltaDivPol(d0_prec, d1_prec, d2_prec, f, n-1);
  UnnormalizedDeltaDivPol(d0_next, d1_next, d2_next, f, n+1);
  
  denom = -PsiDivPol(f, n+1)*PsiDivPol(f, n-1)*sqr(PsiDivPol(f, n)*d2);
  ZZ_pX Pp2, Pn2;
  Pp2 = sqr(PsiDivPol(f, n-1));
  Pn2 = sqr(PsiDivPol(f, n+1));

  eps0 = (-d2*Pp2*d1_next+d2*Pn2*d1_prec+d1*Pp2*d2_next-d1*Pn2*d2_prec)*d0;
  eps1 = sqr(d2)*Pp2*d0_next-sqr(d2)*Pn2*d0_prec-d2*d0*Pp2*d2_next+
    d2*d0*Pn2*d2_prec-d1*d2*Pp2*d1_next+d1*d2*Pn2*d1_prec+
    sqr(d1)*Pp2*d2_next-sqr(d1)*Pn2*d2_prec;
  ZZ_pX Thegcd;
  Thegcd = GCD(eps0, denom);
  Thegcd = GCD(eps1, Thegcd);
  if (deg(Thegcd) >= 1) {
    denom = denom / Thegcd;
    eps0 = eps0 / Thegcd;
    eps1 = eps1 / Thegcd;
  }
}

void DeltaDivPol(ZZ_pX& d0, ZZ_pX& d1, ZZ_pX& d2,
    const ZZ_pX& f, const int n) {
  UnnormalizedDeltaDivPol(d0, d1, d2, f, n);
  
  ZZ_pX D1, D0, x;
  SetX(x);
  d2 = d2 / f;
  D1 = 2*d1 + x*d2;
  D0 = sqr(x)*d2+4*x*d1+16*d0*f;
  d0 = D0; d1 = D1;
}

void EpsilonDivPol(ZZ_pX& eps0, ZZ_pX& eps1, ZZ_pX& denom,
    const ZZ_pX& f, const int n) {
  UnnormalizedEpsilonDivPol(eps0, eps1, denom, f, n);
  ZZ_pX x;
  SetX(x);
  eps0 = 4*f*eps0+x*eps1;
  eps1 = -eps1;
  denom = 4*f*denom;
}	  

#if 0
int main(int argc, char** argv) {
  ZZ p;
  p = to_ZZ("9001");
  ZZ_p::init(p);
  ZZ_pX f;
  SetCoeff(f, 5);
  SetCoeff(f, 4, to_ZZ_p(to_ZZ("123")));
  SetCoeff(f, 3, to_ZZ_p(to_ZZ("54")));
  SetCoeff(f, 2, to_ZZ_p(to_ZZ("6423")));
  SetCoeff(f, 1, to_ZZ_p(to_ZZ("625")));
  SetCoeff(f, 0, to_ZZ_p(to_ZZ("467")));
  
//  for (int i = 0; i < 10; ++i)
//    cout << AlphaDivPol(f, i) << endl;
#if 1
  ZZ_pX pp, qq, rr;
  EpsilonDivPol(pp, qq, rr, f, 23);
  cout << pp << endl;
  cout << deg(pp) << endl;
  cout << qq << endl;
  cout << deg(qq) << endl;
  cout << rr << endl;
  cout << deg(rr) << endl;
#endif
}
#endif
