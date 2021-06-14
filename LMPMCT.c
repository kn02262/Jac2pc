
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"
#include <assert.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "optionparser.h"

#include "RW2dim.h"

#include <iostream>


NTL_CLIENT


void initSeedTime() {
  pid_t pid = getpid();
  long int tm = time(NULL);

  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  tm += tv.tv_usec;
  SetSeed(to_ZZ(pid+tm));
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

int main(int argc, char **argv) {

  assert (sizeof(int)==4);

  cxxopts::Options options("GS-MCT", "Calculates a number of points on Jacobian of Hypereliptic curve over a finite field using a Gaudry-Schost exponential algorithm");
  options.add_options()

  ("p,prime", "A finite field characteristic, a prime number", cxxopts::value<std::string>())
  ("r,rand", "Take a random curve")
  ("s,sp", "Choose starting point in pants")
  ("f,fpoly", "Coefficients of the polynomial, defining a curve y^2=f(x), namely f0,f1,f2,f3,f4 (no whitespaces)", cxxopts::value<std::vector<std::string>>())
  ("m", "A small modulus parameter", cxxopts::value<std::string>()->default_value("1"))
  ("s1p", "s1 mod m coefficient", cxxopts::value<std::string>()->default_value("0"))
  ("s2p", "s2 mod m coefficient", cxxopts::value<std::string>()->default_value("0"))
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


  ZZ_pJac2 curve;

  if(result["rand"].as<bool>()){ // Take a random curve if -r flag is set
	  std::cerr << " !!! RANDOM CURVE !!!";
    ZZ_pX f;
	  randomizeF2(f);
    //ZZ_pJac2::init(f);
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
    //ZZ_pJac2::init( f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3], f_coeffs[4] );
    curve.init(f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3], f_coeffs[4]);
  }

  #if ADD_ALGORITHM == LAUTER // Lauter mode only acceptss a curve with f4=0
    if(ZZ_pJac2::f4() != 0){
      cerr << "Expected a curve with f4=0 for Lauter's addition formulas";
      return 0;
    }
  #endif


  m = to_ZZ( result["m"].as<std::string>().c_str() );
  s1p = to_ZZ( result["s1p"].as<std::string>().c_str() );
  s2p = to_ZZ( result["s2p"].as<std::string>().c_str() );
  nthreads = result["nthreads"].as<int>();

  size_t T = result["T"].as<size_t>();


  /* DEBUGGING */
  cerr << "Modulus p=" << ZZ_p::modulus() << endl;
  cerr << "Curve's equation f=" << ZZ_pJac2::f() << endl;
  cerr << "m=" << m << endl;
  cerr << "s1p=" << s1p << endl;
  cerr << "s2p=" << s2p << endl;
  cerr << "nthreads=" << nthreads << endl;
  cerr << "--------------" << endl;

  // choose a seed
  initSeedTime();

  int r = 15; // kept for comparisions with old runs, TODO: make as input

  RandomWalk2dim testwalk(p,m,s1p,s2p,curve);
  testwalk.InitializeParameters(r, T, nthreads);

  if(result["sp"].as<bool>()){
    testwalk.stpoint = 1;
  } else {
    testwalk.stpoint = 0;
  }

  auto tm_start=std::chrono::high_resolution_clock::now();

  if(nthreads == 1){
	testwalk.start_rw_onethread();
  } else {
    testwalk.start_rw_multithread();
  }

  auto tm_end=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = tm_end - tm_start;

  std::cout << "N: " << testwalk.N << std::endl;
  std::cout << "elapsed time = " << diff.count() << std::endl;
  if (!IsZero(testwalk.N*testwalk.rwparams.BaseDivisor)) {
  //if (!IsZero(testwalk.N*distP1.D) || !IsZero(testwalk.N*distP2.D)) {
   std::cerr << "Order check failed for N = " << testwalk.N << std::endl;
  }
  std::cout << "s1=" << testwalk.s1 << std::endl;
  std::cout << "s2=" << testwalk.s2 << std::endl;

  return 0;

}
