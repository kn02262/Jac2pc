#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"
#include "elltorsion.h"
#include "schoof.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include "optionparser.h"
#include "RW2dim.h"

NTL_CLIENT

// Convert a string to a polynomial.
// The strings contains the coeffs in the magma format:
//   [ a0, a1, ... , an ]
ZZ_pX PolyMagmaToNTL(char *str) {
  size_t n = strlen(str);
  char *str2 = (char *)malloc(sizeof(char)*(n+1));

  char *ptr1, *ptr2;
  ptr1 = str;
  ptr2 = str2;

  for (unsigned int i = 0; i < n; ++i) {
    if (*ptr1 != ',')
      *ptr2++ = *ptr1;
    ++ptr1;
  }
  *ptr2 = '\0';

  string tmp(str2); istringstream myin(tmp);
  ZZ_pX f;
  myin >> f;
  free(str2);
  return f;
}

void randomizeF2(ZZ_pX& f) {
	while(true) {
		ZZ_pX h;
		random(h, 6);
		if (deg(h) != 5) {
			continue;
		}
		MakeMonic(h);

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

void initSeedTime() {
  pid_t pid = getpid();
  long int tm = time(NULL);

  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  tm += tv.tv_usec;
  SetSeed(to_ZZ(pid+tm));
}

#define MAX_STR_LENGTH 1000000    /* should be enough for l=23 */

int main(int argc, char** argv) {

  assert (sizeof(int)==4);

  cxxopts::Options options("GS-MCT", "Helper: generates a list of random curves with s1,s2 modulo m. Call with: ./Random_curve_generator > inputs/filename.txt");
  options.add_options()

  ("p,prime", "A finite field characteristic, a prime number", cxxopts::value<std::string>())
  ("m", "A small modulus parameter", cxxopts::value<int>()->default_value("0"))
  ("n", "Number of curves to be generated", cxxopts::value<int>()->default_value("0"))
  ("h,help", "Show this help")
  ;

  auto result = options.parse(argc, argv);
  
  // Check parameters
  if (result.count("help")){
    std::cerr << options.help() << std::endl;
    return 0;
  }

  if(result.count("prime") != 1){
    std::cerr << "Wrong input parameter, or not specified: p" << std::endl;
    std::cerr << "for help please type ./Random_curve_generator --help" << std::endl;
    return 0;
  }

  ZZ p = to_ZZ( result["prime"].as<std::string>().c_str() );
  ZZ s1p, s2p;
  int m;
  ZZ_pX f;
  ZZ_p::init(p);

  m = result["m"].as<int>() ;
  int n = result["n"].as<int>();

  for(int i=0; i<n; i++){
    randomizeF2(f);
    if(m == 0){
      cout << "-f " << coeff(f,0) << "," << coeff(f,1) << "," << coeff(f,2) << "," << coeff(f,3) << "," << coeff(f,4) << endl;
      continue;
    }

    // Starting polynomial algorithm
    ZZ_pX Res, param, V12, V0overV1;
    EllTorsionIdeal(Res, param, V12, V0overV1, f, m);

    vec_pair_long_long Candidates;
    Schoof(Candidates, Res, param, V12, V0overV1, f, m);

    if (Candidates.length() == 1) {
      s1p = Candidates[0].a;
      s2p = Candidates[0].b;
      cout << "-f " << coeff(f,0) << "," << coeff(f,1) << "," << coeff(f,2) << "," << coeff(f,3) << "," << coeff(f,4);
      cout << " --s1p " << s1p << " --s2p " << s2p << endl;
    } else {
      continue;
    }

  }

return 0;
}

