#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"
#include <assert.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <thread>
#include <mutex>


#include "compat.hpp"
#include "statistics.hpp" // for statistic collection

#include "RW2dim.h"

NTL_CLIENT


// global variables

//DivAndTrack ListPts[1000000];
//list<DivAndTrack*> SortPts = {};
//int ListPtsLen = 0;
//int T=-1;
//std::mutex mtx;
//int no_found = 1;

//CACHELINE_VARIABLE(MCTStatistics, statistics);


/***************************************************************************\
 *   Multi-thread procedures
\***************************************************************************/

void RandomWalk2dim::start_rw_multithread_task(int b){ // Procedure of a single thread
  ZZ_p::init(p); // Segfaults without it. WHY???
  DivAndTrack D;
  DivAndTrack SListPts[this->T];
  size_t SListLen;
  ZZ cpt;

  //std::cout << "rwparams.no_cycle_bound:" << rwparams.no_cycle_bound << std::endl;
  //std::cout << "T = " << this->T << std::endl;
  //std::cout << "rwparams.pD = " << rwparams.pD << std::endl;

  // create variable to collect stats about overjumps
  #if COLLECT_STATISTICS_OVERJUMPS
    auto &&local_stat_overjumps_B1min = merge_on_exit<unsigned long long>([](unsigned long long val)
    {
      statistics.inc_stats_overjumps_B1min(val);
    } );
    auto &&local_stat_overjumps_B1max = merge_on_exit<unsigned long long>([](unsigned long long val)
    {
      statistics.inc_stats_overjumps_B1max(val);
    } );
    auto &&local_stat_overjumps_B2max = merge_on_exit<unsigned long long>([](unsigned long long val)
    {
      statistics.inc_stats_overjumps_B2max(val);
    } );
    auto &&local_stat_overjumps_B2min = merge_on_exit<unsigned long long>([](unsigned long long val)
    {
      statistics.inc_stats_overjumps_B2min(val);
    } );
  #endif

  while(no_collision_found){ // Continue filling ListPts condition

  for(SListLen=0; SListLen<T; SListLen++){
	  
  // --- Find one more point D
  do {
    SelectRandomStartingPoint(D, b);
    for(cpt=0; cpt < rwparams.no_cycle_bound; cpt++){
		JumpOneStep(D);
		// collect stats about overjumps
		#if COLLECT_STATISTICS_OVERJUMPS
			if(D.alpha <= rwparams.B1min) ++local_stat_overjumps_B1min;
			if(D.alpha >= rwparams.B1max) ++local_stat_overjumps_B1max;
			if(D.beta <= rwparams.B2min) ++local_stat_overjumps_B2min;
			if(D.beta >= rwparams.B2max) ++local_stat_overjumps_B2max;
		#endif
		if(IsDistinguished(D.D)){
			// ToDo: Add collision check
			break;
		}
	}
  } while(cpt == rwparams.no_cycle_bound);
    // -- One more point D is found, add to SListPts
      D.b = b;
      D.hash = hashvalue2(D.D); // Divisor ordering hash value
      SListPts[SListLen] = D;

  }
    
    
    //assert(false);
    for(size_t i=0; i<T; i++){ // quadratic sort SListPts[] in non-descending order; FIXME: better sorting!
  	  for(size_t j=0; j<T; j++){
  		  if(SListPts[i].hash < SListPts[j].hash){
  			  D = SListPts[i];
  			  SListPts[i] = SListPts[j];
  			  SListPts[j] = D;
  		  }
  	  }
    }

    // Merge SListPts into ListPts
    mtx.lock();
    list<DivAndTrack*>::iterator it = SortPts.begin();

    //std::cout << "start merging... " <<  std::endl;

    for(size_t i=0; i<T; i++){
  	  ListPts[ListPtsLen+i] = SListPts[i];
  	  while(it != SortPts.end() && (**it).hash < SListPts[i].hash){
  		  it++;
  	  }
  	  while(it != SortPts.end() && (**it).hash == SListPts[i].hash){ // Found a collision

  		  if((**it).D == SListPts[i].D && (**it).b != SListPts[i].b){ // Confirmed collision
  			  no_collision_found = 0;

  			  if((**it).b){
  			  distP1 = **it;
  			  distP2 = SListPts[i];
  			  } else {
  			  distP2 = **it;
  			  distP1 = SListPts[i];
  			  }
  		  }
  		  it++;
  	  }

  	  SortPts.insert(it, &ListPts[ListPtsLen+i]);
    }
    ListPtsLen += T;

    mtx.unlock();

    b = !b;
  } // outer while

  return;
}


/***************************************************************************\
 *   Multi thread
\***************************************************************************/

void RandomWalk2dim::start_rw_multithread(){
  int b = (RandomWord() & 1UL);
  list<DivAndTrack*>::iterator it;

  statistics.clear_statistics();
  int nthreads = rwparams.threads;

  /*
  std::cout << "rwparams.pD: " << rwparams.pD << std::endl;
  for (int k = 0; k < rwparams.r; ++k) {
    for (int kp = 0; kp < rwparams.r; kp++) {
      std::cout<< this->rwparams.O[indexFromKKpBR(k, kp, 1, rwparams.r)]->beta << std::endl;
    }
  }
  */

  // ------------
  // Multi-thread version (with pthread)
  // ------------

  auto task = &RandomWalk2dim::start_rw_multithread_task;

  for (int i = 0; i < nthreads; ++i) { // Start all threads
    threadpool.push([this,b,task](){((*this).*task)(b);});
    b = !b;
  }

  threadpool.wait_work();
  this->calc_N();
  
  //cerr << "ListPtsLen=" << ListPtsLen << endl;
  #if COLLECT_STATISTICS
    statistics.call_print_statistics();
  #endif

  // DEBUGGING
  /*// Fake collision analyzer
  DivAndTrack D;
  cerr << "Detailed info on collisions:" << endl;
  int start = (**SortPts.begin()).hash;
  it=SortPts.begin(); it++;
  for ( ; it!=SortPts.end(); ++it){
	  if((**it).hash == start){
		  cerr << "***** COLLISION *****" << endl;
		  D=**it;
		  cerr << "alpha=" << (**it).alpha << " beta=" << (**it).beta;
		  if((**it).b){cerr << " WILD";} else {cerr << " TAME";}
		  cerr << endl;
		  it--;
		  cerr << "alpha=" << (**it).alpha << " beta=" << (**it).beta;
		  if((**it).b){cerr << " WILD";} else {cerr << " TAME";}
		  cerr << endl;
		  if(D.D != (**it).D){
			  cerr << "BAD COLLISION (due to bad ordering hash function)" << endl;
		  }
		  it++;
	  }
	  if((**it).hash < start){
		  cerr << "SORTING ERROR" << endl;
	  }
	  start = (**it).hash;
  }*/

	return;
}
