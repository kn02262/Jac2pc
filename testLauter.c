#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"
#include <assert.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

#define NTESTS 1000

NTL_CLIENT

void randomizeF2(ZZ_pX& f) {
  while(true) {
    ZZ_pX h;
    random(h, 6);
    if (deg(h) != 5) {
      continue;
    }
    MakeMonic(h);
    SetCoeff(h, 4, 0); // Force f4=0

    ZZ_p d = resultant(h, diff(h));
    if (d != 0) {
      SetCoeff(f, 5, coeff(h, 5));
      SetCoeff(f, 4, coeff(h, 4));
      SetCoeff(f, 3, coeff(h, 3));
      SetCoeff(f, 2, coeff(h, 2));
      SetCoeff(f, 1, coeff(h, 1));
      SetCoeff(f, 0, coeff(h, 0));
      break;
    }
  }
}

int main(int argc, char **argv) {
	int nprimes=7;
	int i,j;
	ZZ primes[nprimes] = {conv<ZZ>("1000033"), conv<ZZ>("1000003"), conv<ZZ>("175899439"), conv<ZZ>("10992336559"), conv<ZZ>("1000993647881"), conv<ZZ>("10022432993647901"), conv<ZZ>("1000022432993647879")};
	ZZ p;
	ZZ_pX f;
	ZZ_pJac2 D1, D2, D1v, D2v, D1c, D2c, sum1, sum2, sum3;
	clock_t t, tm_Cantor, tm_Lange, tm_Lauter;
	
	for(i=0; i<nprimes; i++){
		p = primes[i];
		ZZ_p::init(p);
		// Gen a random curve with f4=0;
		randomizeF2(f);
		cerr << "---- Starting " << NTESTS << " tests with p=" << p << " ----" << endl;
		cerr << "f=" << f << endl;
		ZZ_pJac2::init(f);
		tm_Cantor=0;
		tm_Lange=0;
		tm_Lauter=0;
		for(j=0; j<NTESTS; j++){
			
			random(D1);
			random(D2);
			// Copy D1, D2 to D1v, D2v
			D1v.u2 = D1.u2; D1v.u1 = D1.u1; D1v.u0 = D1.u0; D1v.v1 = D1.v1; D1v.v0 = D1.v0;
			D2v.u2 = D2.u2; D2v.u1 = D2.u1; D2v.u0 = D2.u0; D2v.v1 = D2.v1; D2v.v0 = D2.v0;
			// Copy D1, D2 to D1c, D2c
			D1c.u2 = D1.u2; D1c.u1 = D1.u1; D1c.u0 = D1.u0; D1c.v1 = D1.v1; D1c.v0 = D1.v0;
			D2c.u2 = D2.u2; D2c.u1 = D2.u1; D2c.u0 = D2.u0; D2c.v1 = D2.v1; D2c.v0 = D2.v0;
			
			t=clock();
			addCantor(sum1, D1, D2);
			tm_Cantor += clock()-t;
			
			t=clock();
			addLange(sum2, D1v, D2v);
			tm_Lange += clock()-t;
			
			t=clock();
			addLauter(sum3, D1c, D2c);
			tm_Lauter += clock()-t;
		
			assert(sum2 == sum3);
		}
		cerr << "Cantor time per " << NTESTS << " tests: " << tm_Cantor << " CYCLES = " << ((float)tm_Cantor)/CLOCKS_PER_SEC << "Sec." << endl;
		cerr << "Lange time per " << NTESTS << " tests: " << tm_Lange << " CYCLES = " << ((float)tm_Lange)/CLOCKS_PER_SEC << "Sec." << endl;
		cerr << "Lauter time per " << NTESTS << " tests: " << tm_Lauter << " CYCLES = " << ((float)tm_Lauter)/CLOCKS_PER_SEC << "Sec." << endl;
		cerr << endl;
	}
	cerr << "SUCCESSFULLY" << endl;
}
