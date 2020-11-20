#include "NTL/ZZ_pXFactoring.h"
#include "NTL/vec_ZZ_p.h"

NTL_USE_NNS

vec_ZZ_pX buildSubProducts(const vec_ZZ_p& xi);
void multievalRec(vec_ZZ_p& yi, const ZZ_pX& P, const vec_ZZ_pX& heap);
void multieval(vec_ZZ_p& yi, const ZZ_pX& P, const vec_ZZ_p& xi);
void fastinterpolateRec(ZZ_pX& P, const vec_ZZ_p& xi, const vec_ZZ_p& yi,
    const vec_ZZ_pX& heap);
void fastinterpolate(ZZ_pX& P, const vec_ZZ_p& xi, const vec_ZZ_p& yi);
void evalAtXPlusA(ZZ_pX& Q, const ZZ_pX& P, const ZZ_p& a);

