#ifndef _elemQuo__H_
#define _elemQuo__H_

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEX.h>
#include <iostream>

NTL_CLIENT


/*******************************************************************
 * Tools To work in the following structure:
 *   ( GF(q) )[X] / (X^2 - s*X + p)
 * where s and p are elements of GF(q)==ZZ_pE
 *******************************************************************/
  
class EelemQUO { 
  public:
    static ZZ_pE *s, *p;  // the s and p parameters are static
    ZZ_pE y0, y1;         // represent an element Y = y0 + X*y1;
  
    EelemQUO(const ZZ_pE& Y0, const ZZ_pE& Y1) {
      y0 = Y0;
      y1 = Y1;
    }
  
    EelemQUO() {
      y0 = 0;
      y1 = 0;
    }
     
    EelemQUO(const EelemQUO& a) {
      y0 = a.y0;
      y1 = a.y1;
    }
   
    EelemQUO& operator=(const EelemQUO& a) {
      this->y0 = a.y0;
      this->y1 = a.y1;
      return *this;
    }

    ~EelemQUO() { }
};


void EinitQUO();
void EsetSandP(const ZZ_pE& s, const ZZ_pE& p);
ostream& operator<<(ostream& o, const EelemQUO& X);
istream& operator>>(istream& i, EelemQUO& X);
EelemQUO add(const EelemQUO& Y, const EelemQUO& Z);
EelemQUO subtract(const EelemQUO& Y, const EelemQUO& Z);
EelemQUO mul(const EelemQUO& Y, const EelemQUO& Z);
EelemQUO inv(const EelemQUO& Y);
EelemQUO operator+(const EelemQUO& Y, const EelemQUO& Z);
EelemQUO operator-(const EelemQUO& Y, const EelemQUO& Z);
EelemQUO operator*(const EelemQUO& Y, const EelemQUO& Z);
EelemQUO operator*(const EelemQUO& Y, const ZZ_p& x);
EelemQUO operator*(const ZZ_p& x, const EelemQUO& Y);
void evalAtX1(EelemQUO& Z, const ZZ_pX& f);
void evalAtX2(EelemQUO& Z, const ZZ_pX& f);
EelemQUO to_EelemQUO(const ZZ_pX& x);
inline int IsZero(const EelemQUO& X) { return ((X.y0 == 0) && (X.y1 == 0)); }

/****************************************************************
 ***  END of EelemQUO
 ****************************************************************/



/*******************************************************************
 * Tools To work in the following structure:
 *   ( GF(q) [p] )[X] / (X^2 - s*X + p)
 * where s is an element of GF(q).
 *******************************************************************/

class elemQUO {
  public:
    static ZZ_p *s;  // the s parameter is static
    ZZ_pX y0, y1;    // represent an element Y = y0 + X*y1;

    elemQUO(const ZZ_pX& Y0, const ZZ_pX& Y1) {
      y0 = Y0;
      y1 = Y1;
    }

    elemQUO() {
      y0 = 0;
      y1 = 0;
    }

    elemQUO(const elemQUO& a) {
      y0 = a.y0;
      y1 = a.y1;
    }

    elemQUO& operator=(const elemQUO& a) {
      this->y0 = a.y0;
      this->y1 = a.y1;
      return *this;
    }


    ~elemQUO() { }
};

void initQUO();
void setS(const ZZ_p& s);
ostream& operator<<(ostream& o, const elemQUO& X);
elemQUO add(const elemQUO& Y, const elemQUO& Z);
elemQUO subtract(const elemQUO& Y, const elemQUO& Z);
elemQUO mul(const elemQUO& Y, const elemQUO& Z);

elemQUO operator+(const elemQUO& Y, const elemQUO& Z);
elemQUO operator-(const elemQUO& Y, const elemQUO& Z);
elemQUO operator*(const elemQUO& Y, const elemQUO& Z);

elemQUO operator*(const elemQUO& Y, const ZZ_p& x);
elemQUO operator*(const ZZ_p& x, const elemQUO& Y);



#endif
