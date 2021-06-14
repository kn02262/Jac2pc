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

  cxxopts::Options options("GS-MCT", "Calculates a number of points on Jacobian of Hypereliptic curve over a finite field using a (1) Schoof-Pila polynomial algorithm, (2) Gaudry-Schost exponential algorithm");
  options.add_options()

  ("p,prime", "A finite field characteristic, a prime number", cxxopts::value<std::string>())
  ("r,rand", "Take a random curve")
  ("f,fpoly", "Coefficients of the polynomial, defining a curve y^2=f(x), namely f0,f1,f2,f3,f4 (no whitespaces)", cxxopts::value<std::vector<std::string>>())
  ("m", "A small modulus parameter", cxxopts::value<std::string>()->default_value("1"))
  ("nthreads", "number of threads", cxxopts::value<int>()->default_value("1"))
  ("T", "Cardinality of short lists in multithreaded version, (JUST FOR TESTING, DELETE THIS PARAM!)", cxxopts::value<size_t>()->default_value("1"))
  ("h,help", "Show this help")
  ;

  auto result = options.parse(argc, argv);
  std::chrono::duration<double> ctm = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
	
  // Check parameters
  if (result.count("help")){
    std::cerr << options.help() << std::endl;
    return 0;
  }

  if(result.count("prime") != 1){
    std::cerr << "Wrong input parameter, or not specified: p" << std::endl;
    std::cerr << "for help please type ./LMPMCT --help" << std::endl;
    return 0;
  }

  if(result.count("fpoly") + result.count("rand") != 1){
    std::cerr << "Wrong input parameter, or not specified: f" << std::endl;
    std::cerr << "for help please type ./LMPMCT --help" << std::endl;
    return 0;
  }
	
  // Convert parameters
  ZZ p = to_ZZ( result["prime"].as<std::string>().c_str() );
  ZZ m, s1p, s2p;
  int nthreads;
  ZZ_p::init(p);

  // ToDo: remove f, rewrite polynomial algorithm to take ZZ_pJac2 curve as input
  ZZ_pX f;
  ZZ_pJac2 curve;

  if(result["rand"].as<bool>()){ // Take a random curve if -r flag is set
	cerr << "Random curve!" << endl;
	randomizeF2(f);
    curve.init(f);
  } else { // Take a curve from the input
    std::vector<std::string> fpoly = result["fpoly"].as<std::vector<std::string>>();
    if (fpoly.size() != 5){
      cerr << "Wrong number of coefficients f, expecting f0,f1,f2,f3,f4";
      return 0;
    }
    ZZ_p f_coeffs[5];
    for(int i=0; i<=4; i++){
      f_coeffs[i] = to_ZZ_p(to_ZZ(fpoly[i].c_str()));
    }
    curve.init(f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3], f_coeffs[4]);
  }
  
  f = curve.f();
  m = to_ZZ( result["m"].as<std::string>().c_str() );
  nthreads = result["nthreads"].as<int>();

  size_t T = result["T"].as<size_t>();
  
  // Debugging
  cerr << " ********** INPUT PARAMETERS ********** " << endl;
  cerr << "Modulus p=" << ZZ_p::modulus() << endl;
  cerr << "Curve's equation f=" << curve.f() << endl;
  cerr << "m=" << m << endl;
  cerr << "nthreads=" << nthreads << endl;
  
  cerr << endl << " ********** STARTING POLYNOMIAL ALG ********** " << endl;
  int ell;
  conv(ell, m);
  ZZ_pX Res, param, V12, V0overV1;
  EllTorsionIdeal(Res, param, V12, V0overV1, f, ell);

  vec_pair_long_long Candidates;
  Schoof(Candidates, Res, param, V12, V0overV1, f, ell);

  if (Candidates.length() == 1) {
	  s1p = Candidates[0].a;
	  s2p = Candidates[0].b;
	  cerr << "Polynomial part finished with:" << endl;
	  cerr << "s1m=" << s1p << endl;
	  cerr << "s2m=" << s2p << endl;
    //cout << "(s1, s2) mod ell = " << Candidates[0].a 
    //  << ", " << Candidates[0].b << endl;
  } else {
    //cout << "there are " << Candidates.length() << " candidates remaining\n";
    // ToDo: Restart polynomial part with next_prime(ell)
    cerr << "Polynomial part has not concluded :((((";
    return 0;
  }
  
  cerr << endl << " ********** STARTING EXPONENTIAL ALG ********** " << endl;
  // choose a seed
  initSeedTime();

  int r = 15; // kept for comparisions with old runs, TODO: make as input

  RandomWalk2dim testwalk(p,m,s1p,s2p,curve);
  testwalk.InitializeParameters(r, T, nthreads);

  auto tm_start=std::chrono::high_resolution_clock::now();

  if(nthreads == 1){
	testwalk.start_rw_onethread();
  } else {
    testwalk.start_rw_multithread();
  }

  auto tm_end=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = tm_end - tm_start;

  std::cout << "N: " << testwalk.N << std::endl;
  std::cout << " elapsed time = " << diff.count() << std::endl;
  if (!IsZero(testwalk.N*testwalk.rwparams.BaseDivisor)) {
  //if (!IsZero(testwalk.N*distP1.D) || !IsZero(testwalk.N*distP2.D)) {
   std::cerr << "Order check failed for N = " << testwalk.N << std::endl;
  }
  //std::cout << "s1=" << testwalk.s1 << std::endl;
  //std::cout << "s2=" << testwalk.s2 << std::endl;
  
  
  return 0;
}

