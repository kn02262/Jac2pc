#include <NTL/ZZ_pX.h>

NTL_USE_NNS

void UnnormalizedDeltaDivPol(ZZ_pX& d0, ZZ_pX& d1, ZZ_pX& d2,
    const ZZ_pX& f, const int n); 
void UnnormalizedEpsilonDivPol(ZZ_pX& eps0, ZZ_pX& eps1, ZZ_pX& denom,
    const ZZ_pX& f, const int n);
void DeltaDivPol(ZZ_pX& d0, ZZ_pX& d1, ZZ_pX& d2,
    const ZZ_pX& f, const int n);
void EpsilonDivPol(ZZ_pX& eps0, ZZ_pX& eps1, ZZ_pX& denom,
    const ZZ_pX& f, const int n);
