#include "vec_pair_long_long.h"
#include <NTL/ZZ_pX.h>

NTL_USE_NNS

void Schoof(vec_pair_long_long& Candidates,
    const ZZ_pX& Res, const ZZ_pX& param,
    const ZZ_pX& V12, const ZZ_pX& V0overV1,
    const ZZ_pX& f, const int ell);
  
