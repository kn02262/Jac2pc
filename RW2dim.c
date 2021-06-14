#include "RW2dim.h"
#include <cassert>

void RandomWalk2dim::InitializeParameters(int r, size_t T, int threads = 1)
{
  // Select a random base divisor
  RWParams rwparams_;
  random(rwparams_.BaseDivisor);

  // B[12]min/max
  rwparams_.B1max = (5*SqrRoot(p))/(2*m);
  rwparams_.B1min = - rwparams_.B1max;
  rwparams_.B2min = - ((2*p)/m);
  rwparams_.B2max = (3*p)/m;


  // compute K and K*BaseDivisor and m*BaseDivisor
  rwparams_.K = p*p+1 - s1m*(p+1) + s2m + m*((rwparams_.B2max + rwparams_.B2min)/2);
  rwparams_.KBaseDivisor = rwparams_.K* rwparams_.BaseDivisor;
  rwparams_.mBaseDivisor = m*rwparams_.BaseDivisor;
  rwparams_.pp1mBaseDivisor = (p+1)*rwparams_.mBaseDivisor;


  /*
  std::cout << "m = " << m << std::endl;
  std::cout << "rwparams_.BaseDivisor" << std::endl;
  print(rwparams_.BaseDivisor);

  std::cout << "rwparams_.mBaseDivisor" << std::endl;
  print(rwparams_.mBaseDivisor);
  */

  // NB: the integer pD means that the proba is 2^-pD

  ZZ aux;
  aux = 3*SqrRoot( (rwparams_.B1max-rwparams_.B1min)
                    * (rwparams_.B2max-rwparams_.B2min) );
  aux /= 1000;
  rwparams_.pD = NumBits(aux) - 2;
  rwparams_.no_cycle_bound = to_ZZ((1UL)<<rwparams_.pD)*10;

  // l1 and l2
  rwparams_.l1 = to_double(rwparams_.B1max-rwparams_.B1min)/9.0
    / sqrt(to_double(1UL<<rwparams_.pD));
  rwparams_.l2 = to_double(rwparams_.B2max-rwparams_.B2min)/10.0
    / to_double(1UL<<rwparams_.pD);

  // in the case where l1 << 1, the random walk is different, and we have
  // to choose another value for l1
  if (rwparams_.l1 <= 2) {
    rwparams_.l1 = to_double(rwparams_.B1max-rwparams_.B1min)/9.0
      / (to_double(1UL<<rwparams_.pD));
    assert (rwparams_.l1 <= 2);
  }

  rwparams_.threads = threads;
  rwparams_.no_cycle_bound = to_ZZ((1UL)<<rwparams_.pD)*10; // TODO: check

  this->rwparams = rwparams_;

  this->T = T;

  this->threadpool.resize(threads);

  InitializeOffsets(r);
}


// TODO: remove this->, pass rwparams_ to this->rwparams at the end of InitializeOffsets
void RandomWalk2dim::InitializeOffsets(int r) {
  this->rwparams.r = r;
  //rwdata.b = (RandomWord() & 1UL);
  this->rwparams.O = (DivAndTrack **)malloc(2*r*r*sizeof(DivAndTrack *));
  ZZ ak, bkp;
  ZZ_pJac2 akmP, bkpmP;
  for (int k = 0; k < r; ++k) {
    // select a random alpha_k of average value l1
    // If l1 is too small, then alpha_k is 0, 1, or 2.
    if (this->rwparams.l1 <= 2) {
      if (k == 0) ak = 0;
      else
        if (RandomWord() & 1UL) ak = 1;
        else ak = 2;
    } else {
      ak = RandomBnd(to_ZZ(2*this->rwparams.l1));
    }
    // compute alpha_k*(p+1)*m*P
    akmP = ak*this->rwparams.pp1mBaseDivisor;
    for (int kp = 0; kp < r; kp++) {
      // select a random betakp of average value l2
      bkp = RandomBnd(to_ZZ(2*this->rwparams.l2))+1;
      bkpmP = bkp*this->rwparams.mBaseDivisor;

      // Write the corresponding value in the rwdata.O table
      this->rwparams.O[indexFromKKpBR(k, kp, 0, r)] = new DivAndTrack;
      this->rwparams.O[indexFromKKpBR(k, kp, 0, r)]->D = akmP + bkpmP;
      this->rwparams.O[indexFromKKpBR(k, kp, 0, r)]->alpha = -ak;
      this->rwparams.O[indexFromKKpBR(k, kp, 0, r)]->beta = bkp;

      this->rwparams.O[indexFromKKpBR(k, kp, 1, r)] = new DivAndTrack;
      this->rwparams.O[indexFromKKpBR(k, kp, 1, r)]->D = bkpmP - akmP;
      this->rwparams.O[indexFromKKpBR(k, kp, 1, r)]->alpha = ak;
      this->rwparams.O[indexFromKKpBR(k, kp, 1, r)]->beta = bkp;
    }
  }
}

unsigned int RandomWalk2dim::hashvalue(const ZZ_pJac2& D) {
  //ATOMIC_CPUCOUNT(102);
  ZZ h = rep(D.u1) ^ rep(D.v0);
  ZZ hh = rep(D.u0) ^ rep(D.v1);
  h ^= hh << NumBits(h);
  unsigned int x;
  unsigned char *str = (unsigned char *) (&x);
  BytesFromZZ(str, h, 4);
  return x;
}

unsigned int RandomWalk2dim::hashvalue2(const ZZ_pJac2& D) {
  //ATOMIC_CPUCOUNT(103);
  ZZ h = rep(D.u0) ^ rep(D.u1);
  ZZ hh = rep(D.v0) ^ rep(D.v1);
  h ^= hh << NumBits(h);
  unsigned int x;
  unsigned char *str = (unsigned char *) (&x);
  BytesFromZZ(str, h, 4);
  return x;
}

// Starting point of a new walk.  Bit b indicates Wild or Tame
// kangaroo.
void RandomWalk2dim::SelectRandomStartingPoint(DivAndTrack& D, const int b) {

  if(stpoint == 1){
    D.beta = rwparams.B2min + RandomBnd(rwparams.B2max - rwparams.B2min);
    if(D.beta > 2*p){
    ZZ nB1min = D.beta/SqrRoot(p) - (2*SqrRoot(p))/m;
    ZZ nB1max = D.beta/(2*SqrRoot(p)) + (SqrRoot(p))/m;
    D.alpha = nB1min + RandomBnd(nB1max - nB1min);
    if (RandomWord() & 1UL) D.alpha = -D.alpha;
    } else {
    ZZ nB1max = D.beta/(2*SqrRoot(p)) + (SqrRoot(p))/m;
    ZZ nB1min = -nB1max;
    D.alpha = nB1min + RandomBnd(nB1max - nB1min);
    }
    D.D = (-D.alpha)*rwparams.pp1mBaseDivisor + D.beta*rwparams.mBaseDivisor;
  } else {
    // Starting point in a rectangle [B2min..B2max], [B1min..B1max]
    D.alpha = rwparams.B1min + RandomBnd(rwparams.B1max - rwparams.B1min);
    D.beta = rwparams.B2min + RandomBnd(rwparams.B2max - rwparams.B2min);
    D.D = (-D.alpha)*rwparams.pp1mBaseDivisor + D.beta*rwparams.mBaseDivisor;
  }


  

  /*
  std::cout << "pp1mBaseDivisor" << std::endl;
  print(rwparams.pp1mBaseDivisor);

  std::cout << "mBaseDivisor" << std::endl;
  print(rwparams.mBaseDivisor);
  */

  if (b)
    D.D += rwparams.KBaseDivisor;
}


void RandomWalk2dim::JumpOneStep(DivAndTrack& D) {
  unsigned int hash = hashvalue(D.D); //D.hash;
  int r = rwparams.r;
  int k, kp, b;


  unsigned int il1 = 0;
  if (rwparams.l1 <= 2) {
    il1 = (unsigned int)(1/rwparams.l1);
  }

  // Select k, kp, b according to hash.
  //      hash >>= (rwdata.pD+1);  this is now useless
  b = hash & 1; hash >>= 1;
  kp = hash % r; hash >>= 4;

  if (rwparams.l1 <= 2) {
    if ( (hash % il1) == 0 ) {
      hash >>= 3;
      k = (hash % (rwparams.r - 1)) + 1;
    } else
      k = 0;
  } else {
    k = hash % r;
  }

  D.D += rwparams.O[indexFromKKpBR(k, kp, b, r)]->D;
  D.alpha += rwparams.O[indexFromKKpBR(k, kp, b, r)]->alpha;
  D.beta += rwparams.O[indexFromKKpBR(k, kp, b, r)]->beta;
}

void RandomWalk2dim::calc_N(){
  this->N = this->rwparams.K - (this->distP1.alpha - this->distP2.alpha)*(this->p+1)*this->m +
   (this->distP1.beta - this->distP2.beta)*this->m;
  this->s1 = this->s1m + this->m*(this->distP1.alpha - this->distP2.alpha);
  this->s2 = this->s2m + this->m*((this->rwparams.B2max + this->rwparams.B2min)/2) + this->m*(this->distP1.beta - this->distP2.beta);
	return;
}
