////////////////////////////////////////////////////////////
////  Algorithm by Matsuo Chao Tsujii  (Ants5)
////     (Low memory, parallelizable version)
//////////////////////////////////////////////////////////////

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"
#include <assert.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include "pthread.h"
#include "optionparser.h"
#include <thread>

#define NbOfThreads 1
#define ncoeffs 5

NTL_CLIENT


/***************************************************************************\
 *  Typedef for a divisor with alpha and beta
\***************************************************************************/

class DivAndTrack {
  public:
    ZZ_pJac2 D;
    ZZ alpha;
    ZZ beta;
    
    DivAndTrack() { }
    DivAndTrack(const DivAndTrack& dd) {
      D = dd.D;  alpha = dd.alpha;  beta = dd.beta;
    }
    
    ~DivAndTrack() { }
    
    inline DivAndTrack& operator=(const DivAndTrack& dd) {
      this->D = dd.D;
      this->alpha = dd.alpha;
      this->beta = dd.beta;
      return *this;
    }
    
};


typedef Vec<DivAndTrack> vec_DivAndTrack;



/***************************************************************************\
 *  Typedef for the random walk data
\***************************************************************************/

typedef struct {
  double l1;
  double l2;
  int r;
  DivAndTrack **O;
  ZZ K;
  ZZ_pJac2 KBaseDivisor;
  ZZ_pJac2 mBaseDivisor;
  ZZ_pJac2 pp1mBaseDivisor;
  ZZ_pJac2 BaseDivisor;
  int pD;
  ZZ B1min, B1max;
  ZZ B2min, B2max;
  ZZ p;
} RWData;

inline int indexFromKKpBR(int k, int kp, int b, int r) {
  return (k + r*kp + r*r*b);
}

typedef struct {
  ZZ nb_jump;
  RWData* rwdata;
  int b;
  DivAndTrack* distD;
  int finished;
} arg_struct; // Specific thread data

/***************************************************************************\
 *   Hash function.
\***************************************************************************/
unsigned int hashvalue(const ZZ_pJac2& D) {
  ZZ h = rep(D.u1) ^ rep(D.v0);
  ZZ hh = rep(D.u0) ^ rep(D.v1);
  h ^= hh << NumBits(h);
  unsigned int x;
  unsigned char *str = (unsigned char *) (&x);
  BytesFromZZ(str, h, 4);
  return x;
}

unsigned int hashvalue2(const ZZ_pJac2& D) {
  ZZ h = rep(D.u0) ^ rep(D.u1);
  ZZ hh = rep(D.v0) ^ rep(D.v1);
  h ^= hh << NumBits(h);
  unsigned int x;
  unsigned char *str = (unsigned char *) (&x);
  BytesFromZZ(str, h, 4);
  return x;
}

/***************************************************************************\
 *   Distinguished property
\***************************************************************************/

// proba gives the number of bits that have to be 0 for being
// distinguished
inline int IsDistinguished(const ZZ_pJac2& D, const int proba) {
  return (hashvalue2(D) & (((1UL)<<proba) - 1)) == 0;
}
  
/***************************************************************************\
 *   Random walk functions
\***************************************************************************/

// note: if l1 is too small, then k = 0 gives alpha_k = 0
void InitializeRandomWalk(RWData &rwdata, int r) {
  rwdata.r = r;
  rwdata.O = (DivAndTrack **)malloc(2*r*r*sizeof(DivAndTrack *));
  ZZ ak, bkp;
  ZZ_pJac2 akmP, bkpmP;
  for (int k = 0; k < r; ++k) {
    // select a random alpha_k of average value l1
    // If l1 is too small, then alpha_k is 0, 1, or 2.
    if (rwdata.l1 <= 2) {
      if (k == 0) 
    ak = 0;
      else
        // if (RandomWord() & 1UL)
        if (1)
      ak = 1;
    else
      ak = 2;
    } else {
      ak = RandomBnd(to_ZZ(2*rwdata.l1));
    }
    // compute alpha_k*(p+1)*m*P
    akmP = ak*rwdata.pp1mBaseDivisor;
    for (int kp = 0; kp < r; kp++) {
      // select a random betakp of average value l2
      bkp = RandomBnd(to_ZZ(2*rwdata.l2))+1;
      bkpmP = bkp*rwdata.mBaseDivisor;
 
      // Write the corresponding value in the rwdata.O table
      rwdata.O[indexFromKKpBR(k, kp, 0, r)] = new DivAndTrack;
      rwdata.O[indexFromKKpBR(k, kp, 0, r)]->D = akmP + bkpmP;
      rwdata.O[indexFromKKpBR(k, kp, 0, r)]->alpha = -ak;
      rwdata.O[indexFromKKpBR(k, kp, 0, r)]->beta = bkp;
      
      rwdata.O[indexFromKKpBR(k, kp, 1, r)] = new DivAndTrack;
      rwdata.O[indexFromKKpBR(k, kp, 1, r)]->D = bkpmP - akmP;
      rwdata.O[indexFromKKpBR(k, kp, 1, r)]->alpha = ak;
      rwdata.O[indexFromKKpBR(k, kp, 1, r)]->beta = bkp;
    }
  }
}

void JumpOneStep(DivAndTrack& Dp, const DivAndTrack& D, const RWData &rwdata) {
  unsigned int hash = hashvalue(D.D);
  int r = rwdata.r;
  int k, kp, b;

  unsigned int il1 = 0;
  if (rwdata.l1 <= 2) {
    il1 = (unsigned int)(1/rwdata.l1);
  }

  // Select k, kp, b according to hash.
  //      hash >>= (rwdata.pD+1);  this is now useless
  b = hash & 1; hash >>= 1;
  kp = hash % r; hash >>= 4;
  
  if (rwdata.l1 <= 2) {
    if ( (hash % il1) == 0 ) {
      hash >>= 3;
      k = (hash % (rwdata.r - 1)) + 1;
    } else
      k = 0;
  } else {
    k = hash % r;
  }
  
//  cout << k << " " << kp << " " << b<< " "  << r << endl;
//  cout << indexFromKKpBR(k, kp, b, r) << endl;
  Dp.D = D.D + rwdata.O[indexFromKKpBR(k, kp, b, r)]->D;
  Dp.alpha = D.alpha + rwdata.O[indexFromKKpBR(k, kp, b, r)]->alpha;
  Dp.beta = D.beta + rwdata.O[indexFromKKpBR(k, kp, b, r)]->beta;
}

/***************************************************************************\
 *   One thread
\***************************************************************************/

// Starting point of a new thread.  Bit b indicates Wild or Tame
// kangaroo.
void SelectRandomStartingPoint(DivAndTrack& D, const int b,
    const RWData &rwdata) {
  D.alpha = rwdata.B1min + RandomBnd(rwdata.B1max - rwdata.B1min);
  D.beta = rwdata.B2min + RandomBnd(rwdata.B2max - rwdata.B2min);
  D.D = (-D.alpha)*rwdata.pp1mBaseDivisor + D.beta*rwdata.mBaseDivisor;

  if (b)
    D.D += rwdata.KBaseDivisor;
}

void RunThread(DivAndTrack& distD, int b, const RWData &rwdata, ZZ& nb_jump) {
  DivAndTrack D;
  ZZ no_cycle_bound = to_ZZ((1UL)<<rwdata.pD)*10;
  ZZ cpt;
  
  do {  
    SelectRandomStartingPoint(D, b, rwdata);
    cpt = 0;
    while (!IsDistinguished(D.D, rwdata.pD) && (cpt < no_cycle_bound)) {
      JumpOneStep(D, D, rwdata);
      cpt ++;
    }
  } while (cpt == no_cycle_bound);
  //  cerr << "Distinguished point hit after " << cpt << " jumps.\n";
  distD = D;
  nb_jump = cpt;
}

/***************************************************************************\
 *   Multi-thread procedures
\***************************************************************************/

void thread_init_vars(arg_struct& args, DivAndTrack& distD, int b, RWData &rwdata){ // Load local variables for a single thread
    args.nb_jump = 0;
    args.finished = 0;
    args.rwdata = &rwdata;
    args.b = b;
    args.distD = &distD;
    return;
}

void process (arg_struct* args){ // A procedure of a single thread
  ZZ_p::init(args->rwdata->p);
  DivAndTrack D;
  ZZ no_cycle_bound = to_ZZ((1UL)<<args->rwdata->pD)*10;
  ZZ cpt;
  
  do {  
    SelectRandomStartingPoint(D, args->b, *args->rwdata);
    cpt = 0;
    while (!IsDistinguished(D.D, args->rwdata->pD) && (cpt < no_cycle_bound)) {
      JumpOneStep(D, D, *args->rwdata);
      cpt ++;
    }
  } while (cpt == no_cycle_bound);
  //  cerr << "Distinguished point hit after " << cpt << " jumps.\n";
  args->distD->D = D.D;
  args->distD->alpha = D.alpha;
  args->distD->beta = D.beta;
  args->nb_jump = cpt;
  args->finished = 1;
    return;
}

/***************************************************************************\
 *   Initialize parameters
\***************************************************************************/

void InitializeParameters(RWData &rwdata,
    const ZZ& m, const ZZ& s1p, const ZZ& s2p) {
  ZZ p = ZZ_p::modulus();
  
  // Select a random base divisor
  random(rwdata.BaseDivisor);
  
  // B[12]min/max
  rwdata.B1max = (5*SqrRoot(p))/(2*m);
  rwdata.B1min = - rwdata.B1max;
  rwdata.B2min = - ((2*p)/m);
  rwdata.B2max = (3*p)/m;

  // compute K and K*BaseDivisor and m*BaseDivisor
  rwdata.K = p*p+1 - s1p*(p+1) + s2p + m*((rwdata.B2max + rwdata.B2min)/2);
  rwdata.KBaseDivisor = rwdata.K*rwdata.BaseDivisor;
  rwdata.mBaseDivisor = m*rwdata.BaseDivisor;
  rwdata.pp1mBaseDivisor = (p+1)*rwdata.mBaseDivisor;

  // proba of being distinguished
  // NB: the integer pD means that the proba is 2^-pD
  ZZ aux;
  aux = 3*SqrRoot( (rwdata.B1max-rwdata.B1min) * (rwdata.B2max-rwdata.B2min) );
  aux /= 1000;
  rwdata.pD = NumBits(aux) - 2;

  // l1 and l2
  rwdata.l1 = to_double(rwdata.B1max-rwdata.B1min)/9.0
    / sqrt(to_double(1UL<<rwdata.pD));
  rwdata.l2 = to_double(rwdata.B2max-rwdata.B2min)/10.0
    / to_double(1UL<<rwdata.pD);

  // in the case where l1 << 1, the random walk is different, and we have
  // to choose another value for l1
  if (rwdata.l1 <= 2) {
    rwdata.l1 = to_double(rwdata.B1max-rwdata.B1min)/9.0
      / (to_double(1UL<<rwdata.pD));
    assert (rwdata.l1 <= 2);
  }
  
}

/***************************************************************************\
 ***************************************************************************
 *   MAIN
 ***************************************************************************
\***************************************************************************/

void initSeedTime() {
  pid_t pid = getpid();
  long int tm = time(NULL);
  
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  tm += tv.tv_usec;
  SetSeed(to_ZZ(pid+tm));
}


int main(int argc, char **argv) {
  assert (sizeof(int)==4);
  
  cxxopts::Options options("GS-MCT", "Calculates a number of points on Jacobian of Hypereliptic curve over a finite field using a Gaudry-Schost exponential algorithm");
  options.add_options()
  
  ("p,prime", "A finite field characteristic, a prime number", cxxopts::value<std::string>())
  ("f,fpoly", "Coefficients of the polynomial, defining a curve y^2=f(x), namely f0,f1,f2,f3,f4 (no whitespaces)", cxxopts::value<std::vector<std::string>>())
  ("m", "A small modulus parameter", cxxopts::value<std::string>()->default_value("1"))
  ("s1p", "s1 mod m coefficient", cxxopts::value<std::string>()->default_value("0"))
  ("s2p", "s2 mod m coefficient", cxxopts::value<std::string>()->default_value("0"))
  ("nthreads", "number of threads", cxxopts::value<int>()->default_value("1"))
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
    std::cerr << "for help please type ./LMPMCT --help" << std::endl;
    return 0;
  }
    
  if(result.count("fpoly") != 1){
    std::cerr << "Wrong input parameter, or not specified: f" << std::endl;
    std::cerr << "for help please type ./LMPMCT --help" << std::endl;
    return 0;
  }

  // Convert parameters
  ZZ p = to_ZZ( result["prime"].as<std::string>().c_str() );
  std::vector<std::string> fpoly = result["fpoly"].as<std::vector<std::string>>();
  ZZ m, s1p, s2p;
  int nthreads;
  ZZ_p::init(p);
  if (fpoly.size() != 5){
    cerr << "Wrong number of coefficients f, expecting f0,f1,f2,f3,f4";
    return 0;
  }
  ZZ_p f_coeffs[5];
  for(int i=0; i<=4; i++){
    f_coeffs[i] = to_ZZ_p(to_ZZ(fpoly[i].c_str()));
  }
  #if ADD_ALGORITHM == LAUTER // Lauter mode only acceptss a curve with f4=0
    if(ZZ_pJac2::f4() != 0){
      cerr << "Expected a curve with f4=0 for Lauter's addition formulas";
      return 0;
    }
  #endif
  ZZ_pJac2::init( f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3], f_coeffs[4] );
  m = to_ZZ( result["m"].as<std::string>().c_str() );
  s1p = to_ZZ( result["s1p"].as<std::string>().c_str() );
  s2p = to_ZZ( result["s2p"].as<std::string>().c_str() );
  nthreads = result["nthreads"].as<int>();
  
  /* DEBUGGING */
  cerr << "Modulus p=" << ZZ_p::modulus() << endl;
  cerr << "A curve equation f=" << ZZ_pJac2::f() << endl;
  cerr << "m=" << m << endl;
  cerr << "s1p=" << s1p << endl;
  cerr << "s2p=" << s2p << endl;
  cerr << "nthreads=" << nthreads << endl;
  cerr << "--------------" << endl;
  
  // choose a seed
  initSeedTime();
  
  // initialize the parameters
  cerr << "Initializing random walking parameters..." << endl;
  RWData rwdata;
  InitializeParameters(rwdata, m, s1p, s2p);
  rwdata.p = p;
  InitializeRandomWalk(rwdata, 15);

  cerr << "pD = " << rwdata.pD << endl;
  cerr << "l1 = " << rwdata.l1 << endl;
  cerr << "l2 = " << rwdata.l2 << endl;
  cerr << "range B1 = " << rwdata.B1max - rwdata.B1min << endl;
  cerr << "range B2 = " << rwdata.B2max - rwdata.B2min << endl;

  // receive dist points from slaves 
  vec_DivAndTrack ListWild;
  vec_DivAndTrack ListTame;
  DivAndTrack distP[nthreads], distP1, distP2;
  int found = 0;
  int b;
  ZZ nb_jump, tot_jumps;
  tot_jumps = 0;
  b = (RandomWord() & 1UL);
  
  auto tm_start=std::chrono::high_resolution_clock::now();
  if(nthreads == 1){
  // ------------
  // One thread version (without pthread)
  // ------------
      
  while (!found) {
    b = !b;
    RunThread(distP1, b, rwdata, nb_jump);
    //cerr << "Dist point hit after " << nb_jump << " jumps" << endl;
    //cerr << "Coordinates: "
    //    << to_double(distP1.alpha - rwdata.B1min) / to_double(rwdata.B1max-rwdata.B1min) 
    //    << " "
    //    << to_double(distP1.beta - rwdata.B2min) / to_double(rwdata.B2max-rwdata.B2min) << endl;
    tot_jumps += nb_jump;
    if (b) {
      for (int i = 0; i < ListTame.length(); ++i) {
        if (distP1.D == ListTame[i].D) {
          found = 1;
          distP2 = ListTame[i];
          break;
        }
      }
      append(ListWild, distP1);
    } else {
      for (int i = 0; i < ListWild.length(); ++i) {
        if (distP1.D == ListWild[i].D) {
          found = 1;
          distP2 = ListWild[i];
          break;
        }
      }
      append(ListTame, distP1);
    }
    // cerr << "We have now " << ListWild.length() << " wilds and " 
    //  << ListTame.length() << " tames\n";
  }
  
  if (!b) {
    DivAndTrack tmp;
    tmp = distP1;
    distP1 = distP2;
    distP2 = tmp;
  }
 
 
  } else {
  // ------------
  // Multi-thread version (with pthread)
  // ------------
  arg_struct args[nthreads];
  std::thread threads[nthreads];
  for (int i = 0; i < nthreads; ++i) { // Start all threads
    thread_init_vars(args[i], distP[i], b, rwdata);
    threads[i] = std::thread(process, &args[i]);
    b = !b;
  }
  while (!found) {
    for (int i = 0; i < nthreads; ++i) 
    if(args[i].finished == 1){ // Collect computations and give a new task
      threads[i].join();
      //cerr << "Thread " << i << " hit dist point after " << args[i].nb_jump << " jumps" << endl;
      tot_jumps += args[i].nb_jump;
            
      if (args[i].b) {
        for (int j = 0; j < ListTame.length(); ++j) {
          if (distP[i].D == ListTame[j].D) {
            found = 1;
            distP2 = ListTame[j];
            distP1 = distP[i];
            break;
          }
        }
        append(ListWild, distP[i]);
      } else {
        for (int j = 0; j < ListWild.length(); ++j) {
          if (distP[i].D == ListWild[j].D) {
            found = 1;
            distP1 = ListWild[j];
            distP2 = distP[i];
            break;
          }
        }
        append(ListTame, distP[i]);
      }
            
      if(found == 0){ // Give a new task
          thread_init_vars(args[i], distP[i], b, rwdata);
          threads[i] = std::thread(process, &args[i]);
          b = !b;
      }
    } 
  }
  //cerr << "Overall we have " << ListWild.length() << " wilds and " << ListTame.length() << " tames" << endl;  
  for (int i = 0; i < nthreads; ++i) { // Join all threads
    if(threads[i].joinable()){
      threads[i].join();
    }
  }  
  
      
  } // End if(nthreads == 1)
  auto tm_end=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = tm_end - tm_start;

  cerr << "Total nb of jumps: " << tot_jumps << endl;
  double ratio;
  ratio = to_double(p);
  ratio = pow(ratio, 0.75) / to_double(m);
  ratio = to_double(tot_jumps) / ratio;
  cerr << "Ratio: " << ratio << endl;
  

  ZZ N = rwdata.K - (distP1.alpha - distP2.alpha)*(p+1)*m + 
    (distP1.beta - distP2.beta)*m;

  cout << N << " " << diff.count() << endl;
  
  return 0;
}
