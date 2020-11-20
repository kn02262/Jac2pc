#include "elemQuo.h"
#include "resultantsumroots.h"
#include "fastinterp.h"
#include "ZZ_pXCantorpoly.h"
#include "ZZ_pXResultant.h"
#include <assert.h>

NTL_CLIENT

/****************************************************************
 * Fonctions for reducing the degree of Res
 * and computing the ordinates.
 *****************************************************************/

void evaluateEpsi(EelemQUO& res0, EelemQUO& res1,
    const ZZ_pX& Eps0, const ZZ_pX& Eps1) {
  EelemQUO x, xk;
  double tm;
  
  x.y0 = 0;
  x.y1 = 1;
  
  int Deg = max(deg(Eps0), deg(Eps1));
  int sqDeg = (int)sqrt((double)Deg);
  
  tm = GetTime();
  // Precompute 1, x, x^2, x^3, ... x^sqDeg
  EelemQUO xi[sqDeg+1];
  xi[0].y0 = 1; xi[0].y1 = 0;
  xi[1] = x; 
  for (int i = 2; i <= sqDeg; ++i) {
    xi[i] = xi[i-1]*x;  // TODO: improve this! (NB: ca coute rien...)
  } 
  cerr << "precomputation of 1, x, x^2, x^3, ... x^sqDeg done in "
    << GetTime() - tm << endl;
    
  // Evaluate di's
  res0.y0 = 0;  res0.y1 = 0;
  res1.y0 = 0;  res1.y1 = 0;
  for (int i = 0; i < sqDeg; ++i) {
    res0 = res0 + xi[i]*coeff(Eps0, i);
    res1 = res1 + xi[i]*coeff(Eps1, i);
  } 
  xk = xi[sqDeg];
  EelemQUO tmp0, tmp1;
  for (int k = 1; k*sqDeg <= Deg; ++k) {
    tmp0.y0 = 0;  tmp0.y1 = 0;
    tmp1.y0 = 0;  tmp1.y1 = 0;
    for (int i = 0; i < sqDeg; ++i) {
      tmp0 = tmp0 + coeff(Eps0, i+k*sqDeg)*xi[i];
      tmp1 = tmp1 + coeff(Eps1, i+k*sqDeg)*xi[i];
    }
    res0 = res0 + tmp0*xk;
    res1 = res1 + tmp1*xk;
    xk = xk*xi[sqDeg];
  } 
} 


void evaluateEps0DenomAtX1X2(EelemQUO& evEps0X1, EelemQUO& evEps0X2,
    EelemQUO& evDenomX1, EelemQUO& evDenomX2,
    const ZZ_pX& Eps0, const ZZ_pX& Denom){
  EelemQUO x, xk;
  double tm;
  
  // Evaluation at X1

  x.y0 = 0;
  x.y1 = 1;

  int Deg = max(deg(Eps0), deg(Denom));
  int sqDeg = (int)sqrt((double)Deg);

  tm = GetTime();
  // Precompute 1, x, x^2, x^3, ... x^sqDeg
  EelemQUO xi[sqDeg+1];
  xi[0].y0 = 1; xi[0].y1 = 0;
  xi[1] = x;
  for (int i = 2; i <= sqDeg; ++i) {
    xi[i] = xi[i-1]*x;  // TODO: improve this! (NB: ca coute rien...)
  }
  cerr << "precomputation of 1, x, x^2, x^3, ... x^sqDeg done in "
    << GetTime() - tm << endl;

  // Evaluate di's
  evEps0X1.y0 = 0;  evEps0X1.y1 = 0;
  evDenomX1.y0 = 0;  evDenomX1.y1 = 0;
  for (int i = 0; i < sqDeg; ++i) {
    evEps0X1 = evEps0X1 + xi[i]*coeff(Eps0, i);
    evDenomX1 = evDenomX1 + xi[i]*coeff(Denom, i);
  }
  xk = xi[sqDeg];
  EelemQUO tmp0, tmp1;
  for (int k = 1; k*sqDeg <= Deg; ++k) {
    tmp0.y0 = 0;  tmp0.y1 = 0;
    tmp1.y0 = 0;  tmp1.y1 = 0;
    for (int i = 0; i < sqDeg; ++i) {
      tmp0 = tmp0 + coeff(Eps0, i+k*sqDeg)*xi[i];
      tmp1 = tmp1 + coeff(Denom, i+k*sqDeg)*xi[i];
    }
    evEps0X1 = evEps0X1 + tmp0*xk;
    evDenomX1 = evDenomX1 + tmp1*xk;
    xk = xk*xi[sqDeg];
  }

  // Evaluation at X2 = s - X1

  x.y0 = (*EelemQUO::s);
  x.y1 = -1;
  tm = GetTime();
  // Precompute 1, x, x^2, x^3, ... x^sqDeg
  xi[0].y0 = 1; xi[0].y1 = 0;
  xi[1] = x;
  for (int i = 2; i <= sqDeg; ++i) {
    xi[i] = xi[i-1]*x;  // TODO: improve this! (NB: ca coute rien...)
  }
  cerr << "precomputation of 1, x, x^2, x^3, ... x^sqDeg done in "
    << GetTime() - tm << endl;

  // Evaluate di's
  evEps0X2.y0 = 0;  evEps0X2.y1 = 0;
  evDenomX2.y0 = 0;  evDenomX2.y1 = 0;
  for (int i = 0; i < sqDeg; ++i) {
    evEps0X2 = evEps0X2 + xi[i]*coeff(Eps0, i);
    evDenomX2 = evDenomX2 + xi[i]*coeff(Denom, i);
  }
  xk = xi[sqDeg];
  for (int k = 1; k*sqDeg <= Deg; ++k) {
    tmp0.y0 = 0;  tmp0.y1 = 0;
    tmp1.y0 = 0;  tmp1.y1 = 0;
    for (int i = 0; i < sqDeg; ++i) {
      tmp0 = tmp0 + coeff(Eps0, i+k*sqDeg)*xi[i];
      tmp1 = tmp1 + coeff(Denom, i+k*sqDeg)*xi[i];
    }
    evEps0X2 = evEps0X2 + tmp0*xk;
    evDenomX2 = evDenomX2 + tmp1*xk;
    xk = xk*xi[sqDeg];
  }
}



void ComputeE3ModRes(ZZ_pX& repE3, const ZZ_pX& Res, const ZZ_pX& param,
    const ZZ_pX& Eps0, const ZZ_pX& Eps1) {
  ZZ_pE::init(Res);
  EinitQUO();

  ZZ_pX tmp;
  ZZ_pE Etmp, Etmp2;
  SetX(tmp);
  Etmp = to_ZZ_pE(tmp);
  Etmp2 = to_ZZ_pE(param);
  EsetSandP(Etmp, Etmp2);

  EelemQUO evEps0, evEps1;
  double tm = GetTime();
  evaluateEpsi(evEps0, evEps1, Eps0, Eps1);
  cerr << "evaluation of Eps0 and Eps1 done in "
    << GetTime()-tm << endl;

  ZZ_pE E3;
  E3 = evEps0.y0*evEps1.y1 - evEps0.y1*evEps1.y0;
  repE3 = rep(E3);
}

void ComputeOrdinateModRes(ZZ_pX& V12, ZZ_pX& V0overV1,
    const ZZ_pX& Res, const ZZ_pX& param,
    const ZZ_pX& Eps0, const ZZ_pX& Denom, const ZZ_pX& f) {
  ZZ_pX Y1Y2, Fx1Fx2, DivDiffF, aux;
  ZZ_pE::init(Res);
  EinitQUO();

  {
    ZZ_pX tmp;
    ZZ_pE Etmp, Etmp2;
    SetX(tmp);
    Etmp = to_ZZ_pE(tmp);
    Etmp2 = to_ZZ_pE(param);
    EsetSandP(Etmp, Etmp2);
  }

  EelemQUO evEps0X1, evEps0X2, evDenomX1, evDenomX2;
  evaluateEps0DenomAtX1X2(evEps0X1, evEps0X2, evDenomX1, evDenomX2,
      Eps0, Denom);
  EelemQUO tmp, tmp2, tmp3, fx1, fx2;
  tmp = evEps0X1*evDenomX2;
  tmp = inv(tmp);
  tmp2 = evEps0X2*evDenomX1;
  tmp = tmp*tmp2;
  evalAtX2(fx2, f);
  tmp = tmp*fx2;

  assert (tmp.y1 == 0);
  Y1Y2 = -rep(tmp.y0);

  evalAtX1(fx1, f);
  tmp = fx1 + fx2;
  assert (tmp.y1 == 0);
  Fx1Fx2 = rep(tmp.y0);

  SetX(aux);
  aux = sqr(aux) - 4*param;
  InvMod(aux, aux, Res);
  MulMod(V12, Fx1Fx2 - 2*Y1Y2, aux, Res);

  tmp.y0 = 0; tmp.y1 = 1; // tmp := x1
  tmp = tmp*fx2;
  tmp2.y0 = *EelemQUO::s; tmp2.y1 = -1; // tmp := x2
  tmp3 = to_EelemQUO(Y1Y2);
  tmp2 = tmp2*tmp3;
  tmp = tmp - tmp2;
  tmp2 = tmp3 - fx2;
  tmp2 = inv(tmp2);
  tmp = tmp*tmp2;

  assert (tmp.y1 == 0);
  V0overV1 = rep(tmp.y0);
}

/**************************************************************
 *  fonctions for computing the resultant and the subresultant
 **************************************************************/


void evalAtXplusSover2(ZZ_pX& Q, const ZZ_pX& P) {
  ZZ_p a;
  a = (*elemQUO::s)/2;
  evalAtXPlusA(Q, P, a);
}

void evalAtS2over4MinusX(ZZ_pX& Q, const ZZ_pX& P) {
  ZZ_p a;
  Q = P;
  for (int i = 1; i <= deg(Q); i+= 2) {
    a = - coeff(Q, i);
    SetCoeff(Q, i, a);
  }
  a = -sqr(*elemQUO::s)/4;
  evalAtXPlusA(Q, Q, a);
}

void splitInOddEven(ZZ_pX& Podd, ZZ_pX& Peven, const ZZ_pX& P) {
  Podd = 0;
  Peven = 0;
  for (int i = 0; i <= deg(P)/2; ++i) {
    SetCoeff(Peven, i, coeff(P, 2*i));
    SetCoeff(Podd, i, coeff(P, 2*i+1));
  }
}


void evaluateDi(elemQUO& res0, elemQUO& res1, elemQUO& res2,
    const ZZ_pX& d0, const ZZ_pX& d1, const ZZ_pX& d2) {

  ZZ_pX B, Bodd, Beven;

  evalAtXplusSover2(B, d0);
  splitInOddEven(Bodd, Beven, B);
  evalAtS2over4MinusX(Bodd, Bodd);
  evalAtS2over4MinusX(Beven, Beven);
  res0.y0 = Beven - (*elemQUO::s)/2*Bodd;
  res0.y1 = Bodd;

  evalAtXplusSover2(B, d1);
  splitInOddEven(Bodd, Beven, B);
  evalAtS2over4MinusX(Bodd, Bodd);
  evalAtS2over4MinusX(Beven, Beven);
  res1.y0 = Beven - (*elemQUO::s)/2*Bodd;
  res1.y1 = Bodd;

  evalAtXplusSover2(B, d2);
  splitInOddEven(Bodd, Beven, B);
  evalAtS2over4MinusX(Bodd, Bodd);
  evalAtS2over4MinusX(Beven, Beven);
  res2.y0 = Beven - (*elemQUO::s)/2*Bodd;
  res2.y1 = Bodd;
}

void computeResultants(vec_ZZ_p& res, vec_ZZ_p& subres0, vec_ZZ_p& subres1,
    const ZZ_pX& d0, const ZZ_pX& d1, const ZZ_pX& d2, const vec_ZZ_p si,
    const vec_ZZ_p& parasites, const vec_ZZ_p& parasites_sub, const int n) {

  initQUO();
  elemQUO d0x, d1x, d2x;
  ZZ_pX poly1, poly2;
  double tm = GetTime();

  res.SetLength(n);
  subres0.SetLength(n);
  subres1.SetLength(n);

  ZZ_pX tmp;

  for (int i = 0; i < n; ++i) {
    setS(si[i]);
    evaluateDi(d0x, d1x, d2x, d0, d1, d2);
    poly1 = d0x.y1 * d2x.y0 - d0x.y0 * d2x.y1;
    poly2 = d1x.y1 * d2x.y0 - d1x.y0 * d2x.y1;

    resultantWithSubRes(res[i], tmp, poly2, poly1);
    res[i] /= parasites[i];
    subres0[i] = coeff(tmp, 0) / parasites_sub[i];
    subres1[i] = coeff(tmp, 1) / parasites_sub[i];

    if (((i+1)%1000)==0) {
      cerr << "evaluated " << (i+1) << " resultants in " <<
        GetTime() -tm << " sec\n";
    }
  }
}


// Compute g(x) = f(x/2)
void evaluateAtHalfX(ZZ_pX& g, const ZZ_pX& f) {
  g = f;
  ZZ_p c;
  c = 1;
  for (int i = 0; i <= deg(f); ++i) {
    SetCoeff(g, i, coeff(g, i)/c);
    c *= 2;
  }
}



void computeParasites(vec_ZZ_p& parasites, vec_ZZ_p& parasites_sub,
    const ZZ_pX& d0, const ZZ_pX& d1, const ZZ_pX& d2, const ZZ_pX& f,
    const vec_ZZ_p& si, const int n) {
  ZZ_pX factorD2, rem;

  DivRem(factorD2, rem, d2, sqr(f)*f);
  if (rem != 0) Error("computeParasites: d2 is not divisible by f^3");
  factorD2 = GCD(diff(factorD2), factorD2);
  // Now, we have: d2 = f^3*factorD2^2

  ZZ_pX doubleD2, doubleF, sumRootsD2, sumRootsF, cross;
  evaluateAtHalfX(doubleD2, factorD2);
  MakeMonic(doubleD2);
  ResultantSumRoots(sumRootsD2, factorD2, factorD2);
  sumRootsD2 /= doubleD2;
  sumRootsD2 = GCD(diff(sumRootsD2), sumRootsD2);

  // same stuff for f
  evaluateAtHalfX(doubleF, f);
  MakeMonic(doubleF);
  ResultantSumRoots(sumRootsF, f, f);
  sumRootsF /= doubleF;
  sumRootsF = GCD(diff(sumRootsF), sumRootsF);

  ResultantSumRoots(cross, factorD2, f);

  // multiple evaluate
  // NB:  Parasite     == Cross^6*DoubleD2*SumRootsD2^4*DoubleF^3*SumRootsF^9;
  //      Parasite_Sub == Cross^2*SumRootsD2*SumRootsF^4;

  // precompute the "FanIn" structure
  vec_ZZ_pX subprd = buildSubProducts(si);

  vec_ZZ_p evalTmp;
  ZZ_p tmp;
  evalTmp.SetLength(n);
  parasites.SetLength(n);
  parasites_sub.SetLength(n);

  // cross
  multievalRec(parasites_sub, cross, subprd);
  for (int i = 0; i < n; ++i) {
    sqr(parasites_sub[i], parasites_sub[i]);
    power(parasites[i], parasites_sub[i], 3);
  }
  // doubleD2
  multievalRec(evalTmp, doubleD2, subprd);
  for (int i = 0; i < n; ++i)
    parasites[i] *= evalTmp[i];
  // sumRootsD2
  multievalRec(evalTmp, sumRootsD2, subprd);
  for (int i = 0; i < n; ++i) {
    parasites_sub[i] *= evalTmp[i];
    power(evalTmp[i], evalTmp[i], 4);
    parasites[i] *= evalTmp[i];
  }
  // doubleF  
  multievalRec(evalTmp, doubleF, subprd);
  for (int i = 0; i < n; ++i) {
    power(evalTmp[i], evalTmp[i], 3);
    parasites[i] *= evalTmp[i];
  }
  // sumRootsF
  multievalRec(evalTmp, sumRootsF, subprd);
  for (int i = 0; i < n; ++i) {
    sqr(tmp, evalTmp[i]); sqr(tmp, tmp);
    parasites_sub[i] *= tmp;
    sqr(tmp, tmp); mul(tmp, tmp, evalTmp[i]);
    parasites[i] *= tmp;
  }
}

void evaluatedResultants(vec_ZZ_p& res, vec_ZZ_p& subres0, vec_ZZ_p& subres1,
    vec_ZZ_p& parasites, vec_ZZ_p& parasites_sub, const ZZ_pX &f,
    const ZZ_pX& d0, const ZZ_pX& d1, const ZZ_pX& d2, const vec_ZZ_p& si) {
  int n = si.length();
  
  double tm = GetTime();
  computeParasites(parasites, parasites_sub, d0, d1, d2, f, si, n);
  cerr << "parasites computed in " << (GetTime()-tm) << " sec\n";

  tm = GetTime();
  computeResultants(res, subres0, subres1, d0, d1, d2, si,
      parasites, parasites_sub, n);
  cerr << "all resultants computed in " << (GetTime()-tm) << " sec\n";
}


void ResultantAndSubResultants(ZZ_pX& Res, ZZ_pX& Subres0, ZZ_pX& Subres1,
    const ZZ_pX &f, const ZZ_pX& d0, const ZZ_pX& d1, const ZZ_pX& d2,
    const int ell) {
  
  int n;
  vec_ZZ_p si;
  double tm;
  
  n = ( (ell*ell*(7*ell*ell - 33)) / 2 ) + 22 + 5;
  si.SetLength(n);
  for (int i = 0; i < n; ++i)
    si[i] = i+1;
  
  vec_ZZ_p parasites, parasites_sub, res, subres0, subres1;
  evaluatedResultants(res, subres0, subres1, parasites, parasites_sub,
      f, d0, d1, d2, si);
  
  // Interpolate
  tm = GetTime();
  fastinterpolate(Res, si, res);
  MakeMonic(Res);
  fastinterpolate(Subres0, si, subres0);
  fastinterpolate(Subres1, si, subres1);
  cerr << "interpolation computed in " << (GetTime()-tm) << " sec\n";
  cerr << "  Res      is of degree  " << deg(Res) << endl;
  cerr << "  Subres0  is of degree  " << deg(Subres0) << endl;
  cerr << "  Subres1  is of degree  " << deg(Subres1) << endl;
}


void EllTorsionIdeal(ZZ_pX& Res, ZZ_pX& param, ZZ_pX& V12, ZZ_pX& V0overV1,
    const ZZ_pX& f, const int ell) {
  ZZ_pX d0, d1, d2, Eps0, Eps1, Denom;

  DeltaDivPol(d0, d1, d2, f, ell);
  EpsilonDivPol(Eps0, Eps1, Denom, f, ell);
  cerr << "delta and epsilon computed.\n";
	
  ZZ_pX Subres0, Subres1;
  ResultantAndSubResultants(Res, Subres0, Subres1, f, d0, d1, d2, ell);
  

  // recover the parametrization
  double tm = GetTime();
  rem(Subres1, Subres1, Res);
  rem(Subres0, Subres0, Res);
  InvMod(Subres1, Subres1, Res);
  MulMod(param, Subres1, Subres0, Res);
  NTL_NNS negate(param, param);
  cerr << "parametrization reconstructed in " << GetTime()-tm << endl;

  clear(Subres0);
  clear(Subres1);
  
       
  // Compute E3 modulo this parametrization
  ZZ_pX E3;
  ComputeE3ModRes(E3, Res, param, Eps0, Eps1);

  // Take the GCD with Res and update the parametrization
  GCD(Res, Res, E3);
  cerr << "deg(Res) = " << deg(Res) << endl;

  rem(param, param, Res);

  // Build the ordinates
  tm = GetTime();
  ComputeOrdinateModRes(V12, V0overV1, Res, param, Eps0, Denom, f);
  cerr << "ordinate constructed in " << GetTime()-tm << endl;
}  
