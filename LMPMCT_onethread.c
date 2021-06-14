#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "ZZ_pJac2.h"
#include <assert.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

#include "RW2dim.h"

NTL_CLIENT

/***************************************************************************\
 *   One thread
\***************************************************************************/

void RandomWalk2dim::start_rw_onethread(){
DivAndTrack D;
int b = (RandomWord() & 1UL);
ZZ cpt;

while(no_collision_found){ // Continue filling ListPts condition
// Find one more point D
	b = !b;
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
				break;
			}
		}
	} while(cpt == rwparams.no_cycle_bound);

	// -- One more point D is found, add to ListPts
	D.b = b;
	D.hash = hashvalue2(D.D); // Divisor ordering hash value

	ListPts[++ListPtsLen] = D;

	list<DivAndTrack*>::iterator it = SortPts.begin();
	while(it != SortPts.end() && (**it).hash < D.hash){
		it++;
	}
	while(it != SortPts.end() && (**it).hash == D.hash){ // Found a collision
		if((**it).D == D.D && (**it).b != D.b){ // Confirmed collision
			no_collision_found = 0;
			if((**it).b){
				distP1 = **it;
				distP2 = D;
			} else {
				distP2 = **it;
				distP1 = D;
			}
		}
		it++;
	}
	SortPts.insert(it, &ListPts[ListPtsLen]);

} // end while(no_collision_found)

this->calc_N();
return;

}
