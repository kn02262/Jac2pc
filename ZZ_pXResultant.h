#ifndef NTL_ZZ_pXResultant__H
#define NTL_ZZ_pXResultant__H

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

NTL_OPEN_NNS

void resultantWithSubRes(ZZ_p& rres, ZZ_pX& subres, 
    const ZZ_pX& u, const ZZ_pX& v);

NTL_CLOSE_NNS

#endif


