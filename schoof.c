#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include "ZZ_pEJac2F.h"
#include "vec_pair_long_long.h"
#include <assert.h>
#include <iostream>

NTL_CLIENT

// NB: Modular Composition is useless here.
// (except maybe for small ells and large q)
ZZ_pEJac2F Frobenius(const ZZ_pEJac2F& D) {
  ZZ_pEJac2F Dq(D.Sq);

  ZZ q;
  q = ZZ_p::modulus();
  
  assert (D.u2 == 1);
  Dq.u2 = 1;
  power(Dq.u1, D.u1, q);
  power(Dq.u0, D.u0, q);
  power(Dq.v1, D.v1, q);
  power(Dq.v0, D.v0, q);
  
  ZZ_pE Fact;
  power(Fact, D.Sq, (q-1)/2);
  mul(Dq.v1, Dq.v1, Fact);
  mul(Dq.v0, Dq.v0, Fact);
  return Dq;
}   

// Compute Fr(D), Fr(Fr(D)), Fr^3(D), Fr^4(D).
void Frobeniuses(ZZ_pEJac2F& FD, ZZ_pEJac2F& FFD, ZZ_pEJac2F& FFFD,
    ZZ_pEJac2F& FFFFD, const ZZ_pEJac2F& D) {
  ZZ q;
  q = ZZ_p::modulus();
  assert (D.u2 == 1);
  FD.u2 = 1; FFD.u2 = 1; FFFD.u2 = 1; FFFFD.u2 = 1;
  FD.Sq = D.Sq; FFD.Sq = D.Sq; FFFD.Sq = D.Sq; FFFFD.Sq = D.Sq;
  
  ZZ_pE Fact;
  power(Fact, D.Sq, (q-1)/2);
  
  power(FD.u1, D.u1, q); power(FD.u0, D.u0, q);
  power(FD.v1, D.v1, q); power(FD.v0, D.v0, q);
  mul(FD.v1, FD.v1, Fact); mul(FD.v0, FD.v0, Fact);
    
  power(FFD.u1, FD.u1, q); power(FFD.u0, FD.u0, q);
  power(FFD.v1, FD.v1, q); power(FFD.v0, FD.v0, q);
  mul(FFD.v1, FFD.v1, Fact); mul(FFD.v0, FFD.v0, Fact);

  power(FFFD.u1, FFD.u1, q); power(FFFD.u0, FFD.u0, q);
  power(FFFD.v1, FFD.v1, q); power(FFFD.v0, FFD.v0, q);
  mul(FFFD.v1, FFFD.v1, Fact); mul(FFFD.v0, FFFD.v0, Fact);
  
  power(FFFFD.u1, FFFD.u1, q); power(FFFFD.u0, FFFD.u0, q);
  power(FFFFD.v1, FFFD.v1, q); power(FFFFD.v0, FFFD.v0, q);
  mul(FFFFD.v1, FFFFD.v1, Fact); mul(FFFFD.v0, FFFFD.v0, Fact);
} 





typedef struct {
  ZZ_pX u0, u1, u2, v0, v1, Sq;
} Saved_div;
  
void SaveDiv(Saved_div& SD, const ZZ_pEJac2F& D) {
  SD.u0 = rep(D.u0); SD.u1 = rep(D.u1); SD.u2 = rep(D.u2);
  SD.v0 = rep(D.v0); SD.v1 = rep(D.v1); SD.Sq = rep(D.Sq);
}

void RestoreDiv(ZZ_pEJac2F& D, Saved_div& SD) {
  D.u0 = to_ZZ_pE(SD.u0); D.u1 = to_ZZ_pE(SD.u1); D.u2 = to_ZZ_pE(SD.u2);
  D.v0 = to_ZZ_pE(SD.v0); D.v1 = to_ZZ_pE(SD.v1); D.Sq = to_ZZ_pE(SD.Sq);
}


void HandleFactor(vec_ZZ_pX& ListFactorRes) {
  ZZ_pE_FoundFactorOfModulus = 0;
  // make room in the list for new factors
  int len = ListFactorRes.length();
  ListFactorRes.SetLength(len+2);
  for (int i = len - 1; i >= 0; --i)
    ListFactorRes[i+2] = ListFactorRes[i];

  // add the 2 new factors
  ListFactorRes[0] = *ZZ_pE_FactorOfModulus;
  ListFactorRes[1] = ZZ_pE::modulus() / (*ZZ_pE_FactorOfModulus);
}


void PlainSchoof(vec_pair_long_long& Candidates, const ZZ_pEJac2F& D,
    const int ell) {
  double tm = GetTime(); 
 
  // Compute Frobeniuses of D
  ZZ_pEJac2F FD, FFD, FFFD, FFFFD;
  tm = GetTime(); 
  Frobeniuses(FD, FFD, FFFD, FFFFD, D);
  cerr << "Frobeniuses computed in " << GetTime() - tm << endl;

  // Deduce some factors of Res == ZZ_pEModulus()
  // Remember:
  //   Res(X) is Monic
  //   D.u1 = -X mod Res(X)
  //   FD.u1 = -X^q mod Res(X)
  //   ...
  tm = GetTime();
  const ZZ_pX& Res = ZZ_pE::modulus().val();
  assert (LeadCoeff(Res) == 1);
  vec_ZZ_pX ListFactorRes;
  ZZ_pX gg, X, aux;
  SetX(X);

  // linear factors
  gg = GCD(Res, -rep(FD.u1)-X);
  if (gg == 1) {
    append(ListFactorRes, Res);
  } else {
    cerr << "  found some linear factors of Res.\n";
    append(ListFactorRes, gg);
    append(ListFactorRes, Res/gg);
  }
  // factors of degree 2
  aux = ListFactorRes[ListFactorRes.length()-1];
  gg = GCD(aux, -rep(FFD.u1)-X);
  if (gg != 1) {
    cerr << "  found some factors of deg 2.\n";
    ListFactorRes[ListFactorRes.length()-1] = gg;
    append(ListFactorRes, aux / gg);
  }
  // factors of degree 3
  aux = ListFactorRes[ListFactorRes.length()-1];
  gg = GCD(aux, -rep(FFFD.u1)-X);
  if (gg != 1) {
    cerr << "  found some factors of deg 3.\n";
    ListFactorRes[ListFactorRes.length()-1] = gg;
    append(ListFactorRes, aux / gg);
  }
  // factors of degree 4
  aux = ListFactorRes[ListFactorRes.length()-1];
  gg = GCD(aux, -rep(FFFFD.u1)-X);
  if (gg != 1) {
    cerr << "  found some factors of deg 4.\n";
    ListFactorRes[ListFactorRes.length()-1] = gg;
    append(ListFactorRes, aux / gg);
  }
  cerr << "checked easy factorization in " 
    << GetTime() - tm << endl;
  
  while (ListFactorRes[ListFactorRes.length()-1] == 1)
    ListFactorRes.SetLength(ListFactorRes.length()-1);
  
  
  {
    ZZ_pX Grr;
    Grr = 1;
    for (int k = 0; k < ListFactorRes.length(); ++k) {
      Grr *= ListFactorRes[k];
    }
    assert (Grr == Res);
  }
  
  
  tm = GetTime();
  // Do the Plain Schoof step:
  ZZ q, qmodl, q2modl;
  q = ZZ_p::modulus();
  qmodl = q % ell; q2modl = sqr(q) % ell;


  Saved_div SD, SFD, SFFD, SFFFD, SFFFFD;
  ZZ_pX f0, f1, f2, f3, f4;
  f0 = rep(ZZ_pEJac2F::f0()); f1 = rep(ZZ_pEJac2F::f1());
  f2 = rep(ZZ_pEJac2F::f2()); f3 = rep(ZZ_pEJac2F::f3());
  f4 = rep(ZZ_pEJac2F::f4());
  SaveDiv(SD, D);
  SaveDiv(SFD, FD);
  SaveDiv(SFFD, FFD);
  SaveDiv(SFFFD, FFFD);
  SaveDiv(SFFFFD, FFFFD);
  
  // TODO: should sort by increasing degree, here.
  Candidates.SetLength(0);


  while (ListFactorRes.length() != 0) {
    long test = 0;
    cerr << "computing modulo a factor of degree " 
      << deg(ListFactorRes[0]) << "..." << endl;
    // create objects modulo a factor of Res
    ZZ_pE::init(ListFactorRes[0]);
    ZZ_pEJac2F::init(to_ZZ_pE(f0), to_ZZ_pE(f1), to_ZZ_pE(f2),
	to_ZZ_pE(f3), to_ZZ_pE(f4));
    ZZ_pEJac2F DD, FDD, FFDD, FFFDD, FFFFDD;
    RestoreDiv(DD, SD);
    RestoreDiv(FDD, SFD);
    RestoreDiv(FFDD, SFFD);
    RestoreDiv(FFFDD,  SFFFD);
    RestoreDiv(FFFFDD, SFFFFD);

    // purge ListFactorRes
    int len;
    len = ListFactorRes.length();
    for (int i = 0; i < len-1; ++i)
      ListFactorRes[i] = ListFactorRes[i+1];
    ListFactorRes.SetLength(len-1);
    
    // test all values of (s1, s2) mod ell
    ZZ_pEJac2F LHS, RHS, qD1D3;
    mul(LHS, DD, q2modl);
    if (ZZ_pE_FoundFactorOfModulus) {
      HandleFactor(ListFactorRes);
      continue;
    }
    add(LHS, LHS, FFFFDD);
    if (ZZ_pE_FoundFactorOfModulus) {
      HandleFactor(ListFactorRes);
      continue;
    }

    ZZ_pEJac2F* TableLHS;
    TableLHS = new ZZ_pEJac2F[ell];
    TableLHS[0] = LHS;
    for (int i = 1; i < ell; ++i) {
      add(TableLHS[i], TableLHS[i-1], FFDD);
      if (ZZ_pE_FoundFactorOfModulus) {
	HandleFactor(ListFactorRes);
	test = 1;
	break; 
      }
    }
    if (test) continue;
    
    clear(RHS, LHS.Sq);
    mul(qD1D3, FDD, qmodl);
    if (ZZ_pE_FoundFactorOfModulus) {
      HandleFactor(ListFactorRes);
      continue;
    }
    add(qD1D3, qD1D3, FFFDD);
    if (ZZ_pE_FoundFactorOfModulus) {
      HandleFactor(ListFactorRes);
      continue;
    }
    pair_long_long cand;
    vec_pair_long_long newCandidates;
    for (int i = 0; i < ell; ++i) {
      if (RHS == TableLHS[i]) {
	cerr << "Found a match: 0, " << i << endl;
	cand.a = 0; cand.b = i;
	append(newCandidates, cand);
      }
      if (qD1D3 == TableLHS[i]) {
	cerr << "Found a match: 1, " << i << endl;
	cand.a = 1; cand.b = i;
	append(newCandidates, cand);
      }
    }
    RHS = qD1D3;
    for (int s1 = 2; s1 < ell; ++s1) {
      add(RHS, RHS, qD1D3);
      if (ZZ_pE_FoundFactorOfModulus) {
	HandleFactor(ListFactorRes);
	test = 1;
	break;
      }
      for (int i = 0; i < ell; ++i) {
	if (RHS == TableLHS[i]) {
	  cerr << "Found a match: " << s1 << ", " << i << endl;
	  cand.a = s1; cand.b = i;
	  append(newCandidates, cand);
	}
      }
    }
    delete [] TableLHS;
    if (test) continue;
    
    // Merge lists of candidates
    if (Candidates.length() == 0)
      Candidates = newCandidates;
    else {
      vec_pair_long_long cc;
      for (int i = 0; i < Candidates.length(); ++i) {
	for (int j = 0; j < newCandidates.length(); ++j) {
	  if (Candidates[i] == newCandidates[j]) {
	    append(cc, Candidates[i]);
	    break;
	  }
	}
      }
      Candidates=cc;
    }

    if (Candidates.length() == 0) {
      Error("<torsion> No more candidates (s1,s2) mod ell...\n");
    }
    if (Candidates.length() == 1)
      break;
  }
  
  cerr << "Full Schoof computed in " << GetTime() - tm << endl;
}


void Schoof(vec_pair_long_long& Candidates, 
    const ZZ_pX& Res, const ZZ_pX& param,
    const ZZ_pX& V12, const ZZ_pX& V0overV1,
    const ZZ_pX& f, const int ell) {

  {
    ZZ_pE::init(Res);
    ZZ_pEX fE;
    fE = to_ZZ_pEX(f);
    ZZ_pEJac2F::init(fE);
  }           
  ZZ_pEJac2F D;

  clear(D, to_ZZ_pE(V12));
  D.u2 = 1;
  {
    ZZ_pX grr;
    SetX(grr);
    D.u1 = -to_ZZ_pE(grr);
  }
  D.u0 = to_ZZ_pE(param);
  D.v1 = 1;
  D.v0 = to_ZZ_pE(V0overV1);

  // bench ZZ_pE
  if (0)
  {
    ZZ_pE x, y, z;
    random(x); random(y);
    double tm = GetTime();
    for (int i = 0; i < 10; ++i)
      mul(z, x, y);
    cerr << "<torsion> 10 mul mod Res done in " << GetTime() - tm << endl;
    tm = GetTime();
    inv(z, x);
    cerr << "<torsion> 1 inv mod Res done in " << GetTime() - tm << endl;
  }

  PlainSchoof(Candidates, D, ell);
}

