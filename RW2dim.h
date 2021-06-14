#ifndef RW2DIM_H
#define RW2DIM_H

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"

#include <list>

#include "statistics.hpp" //collecting statistics
#include "thread_pool.hpp" // managing threads

NTL_CLIENT


// MACRO to enable CPU cycle counters. Use ATOMIC_CPUCOUNT for threaded functions.
// Used counters should be registered in cpuperf.cpp.
// DO NOT ENABLE THIS MANUALLY, use 'rebuild.sh -c'
/****  PERFORMANCE_COUNTING *****/
#ifdef PERFORMANCE_COUNTING
#include "cpuperformance.hpp"
 #define CPUCOUNT(s) cpu::update_performance_counter _cpuperf##s(perfcounters[s]);
 #define ATOMIC_CPUCOUNT(s) cpu::update_atomic_performance_counter _atomic_cpuperf##s(atomic_perfcounters[s]);
 extern std::vector<uint64_t> perfcounters;
 extern std::vector<std::atomic<uint64_t>> atomic_perfcounters;
 extern cpu::performance_counter_manager perfmanager;
#else
 #define CPUCOUNT(s)
 #define ATOMIC_CPUCOUNT(s)
#endif
void show_cpu_stats();

//#define XSTR(x) STR(x)
//#define STR(x) #x
//#pragma message "The value of PERFORMANCE_COUNTING: " XSTR(PERFORMANCE_COUNTING)

/***************************************************************************\
 *  Typedef for a divisor with alpha and beta
\***************************************************************************/

class DivAndTrack {
  public:
    ZZ_pJac2 D; // Mumford coordinates of the divisor
    ZZ alpha; // dlog wrt the base divisor on the horizontal search-line
    ZZ beta;  // dlog wrt the base divisor on the vertical search-line
    
    unsigned hash; // Divisor ordering hash value
    int b; 
    /* b=1 --> Wild Pt
     * b=0 --> Tame Pt
     * b=-1 --> Wild/Tame undefined
     */

    //TODO: Make move constructor
    DivAndTrack() { }
    DivAndTrack(const DivAndTrack& dd) {
      D = dd.D;  alpha = dd.alpha;  beta = dd.beta;  b = -1; hash = 0;
    }

    ~DivAndTrack() { }

    inline DivAndTrack& operator=(const DivAndTrack& dd) {
      this->D = dd.D;
      this->alpha = dd.alpha;
      this->beta = dd.beta;
      this->b = dd.b;
      this->hash = dd.hash;
      return *this;
    }

};

typedef Vec<DivAndTrack> vec_DivAndTrack;

typedef struct {
  double l1;// av. length of horizontal offsets
  double l2; // av. length of vertical offsets
  int r; // total number of offsets is 2*r^2
  DivAndTrack **O;// offsets
  ZZ K; // p^2 + 1 - s1m*(p+1)+m*(TODO)
  ZZ_pJac2 KBaseDivisor; //K*BaseDivisor
  ZZ_pJac2 mBaseDivisor; // m*BaseDivisor
  ZZ_pJac2 pp1mBaseDivisor; //(p+1)*m*BaseDivisor
  ZZ_pJac2 BaseDivisor; // random divisor D
  int pD; // invers probability of distringuished point
  ZZ no_cycle_bound; // max length of one walk
  ZZ B1min, B1max; // search for s1 \in [B1min, B1max]
  ZZ B2min, B2max; // search for s2 \in [B2min, B2max]
  int threads; // number of threads
} RWParams;

typedef struct {
  ZZ nb_jump;
  int b;
} arg_struct; // Thread-specific data



class RandomWalk2dim {
public:

  RandomWalk2dim(const ZZ &p_, const ZZ &m_, const ZZ &s1m_, const ZZ &s2m_, const ZZ_pJac2 &curve_, unsigned long int seed = 0):
               p(p_), m(m_), s1m(s1m_), s2m(s2m_), curve(curve_)
  #ifdef PERFORMANCE_COUNTING
    , _totalcpu(perfcounters[0])
  #endif
    {
        //std::cout << "inside RandomWalk2dim constructor" << std::endl;
        //p = p_;
        // m = m_;
        //s1m = s1m_;
        //s2m = s1m_;
        //curve = curve_;
    };


  void InitializeParameters(int r, size_t T, int threads); //implemented in RW2dim.c
  inline int indexFromKKpBR(int k, int kp, int b, int r) { return (k + r*kp + r*r*b); }

  void InitializeOffsets(int r); //implemented in RW2dim.c

  void start_rw_multithread(); //implemented in LMPMCT_multithread.c
  void start_rw_multithread_task(int b); //implemented in LMPMCT_multithread.c

  void start_rw_onethread(); //implemented in LMPMCT_onethread.c

public:

  ZZ p; // Modulus
  ZZ m; // we know s1, s2 mod m from prepocessing phase
  ZZ s1m; // s1m = s1 mod m
  ZZ s2m; // s2m = s2 mod m
  ZZ_pJac2 curve; // curve's coefficients
  ZZ N; // Output - #Jac
  ZZ s1, s2; // Output - s1 and s2
  int stpoint = 0; // ToDo: Delete this param, sets a starting point choose method, 0=rectangular, 1=pants

  RWParams rwparams;


private:

  list<DivAndTrack*> SortPts = {};
  DivAndTrack ListPts[100000];
  DivAndTrack distP1, distP2; // two colliding distinguished points
  
  int ListPtsLen = 0;
  size_t T;
  std::mutex mtx;
  bool no_collision_found = true;

  thread_pool::thread_pool threadpool;

  CACHELINE_VARIABLE(MCTStatistics, statistics);

  unsigned int hashvalue(const ZZ_pJac2& D); //implemented in RW2dim.c
  unsigned int hashvalue2(const ZZ_pJac2& D); //implemented in RW2dim.c
  inline int IsDistinguished(const ZZ_pJac2& D) { return (hashvalue2(D) & (((1UL)<<rwparams.pD) - 1)) == 0; }

  void SelectRandomStartingPoint(DivAndTrack& D, const int b); //implemented in RW2dim.c
  void JumpOneStep(DivAndTrack& D); //implemented in RW2dim.c
  void calc_N(); //implemented in RW2dim.c


};


#endif
