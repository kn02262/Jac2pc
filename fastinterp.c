#include <NTL/ZZ_pX.h>
#include <math.h>
#include <assert.h>

NTL_CLIENT

/* compile with:

g++ fastinterp.c -I../../tmp/ntl-5.3/include ../../tmp/ntl-5.3/src/ntl.a -lgmp -lm

*/






////////////////////////////////////////////////////
// cf Algo 10.3, page 281 of Modern Computer Algebra
// 
//////////////////////////////////////////////////////////////
//// Do it with a heap.

vec_ZZ_pX buildSubProducts(const vec_ZZ_p& xi) {
  vec_ZZ_pX heap;

  int n = xi.length();
  // compute k = Ceiling(Log_2(n));
  int k = 0;
  while ((1 << k) < n) 
    ++k;
  
  // allocate and initialize the heap
  int m = 2*(1 << k);
  heap.SetLength(m);
  // fill in with the (X-xi)
  m >>= 1;
  ZZ_pX f;
  SetCoeff(f, 1, 1);
  for (int i = 0; i < n; ++i) {
    SetCoeff(f, 0, -xi[i]);
    heap[m + i] = f;
  }
  for (int i = n; i < m; ++i) 
    heap[m + i] = 1;
  
  // build the heap
  // NB: the element heap[0] is not used.
  m >>= 1;
  while (m > 0) {
    for (int i = 0; i < m; ++i)
      heap[m + i] = heap[2*(m+i)] * heap[2*(m+i)+1];
    m >>= 1;
  }

  return heap;
}

///////////////////////////////////////////////////
//   Multipoint evaluation.
//
// Recursive algorithm, see algo 10.5, page 282 of Modern Computer
// Algebra
//
///////////////////////////////////////////////////

void multievalRec(vec_ZZ_p& yi, const ZZ_pX& P, const vec_ZZ_pX& heap) {
  int m = heap.length();
  vec_ZZ_pX v, w;

  append(v, P % heap[1]);
  
  for(int mm = 2; mm < m; mm <<= 1) {
    for (int i = 0; 2*i < mm; i++) {  // TODO: do it in place by doing i--
      append(w, v[i] % heap[mm + 2*i]);
      append(w, v[i] % heap[mm + 2*i + 1]);
    }
    v = w;
    w.kill();
  }

  m >>= 1;
  for(int i = 0; i < yi.length(); ++i) {
    if (deg(v[i]) != 0) {
      cout << "ERROR\n";
      exit(1);
    }
    yi[i] = ConstTerm(v[i]);
  }
}

void multieval(vec_ZZ_p& yi, const ZZ_pX& P, const vec_ZZ_p& xi) {
  vec_ZZ_pX SubPrd = buildSubProducts(xi);
  yi.SetLength(xi.length());
  
  multievalRec(yi, P, SubPrd);
}


///////////////////////////////////////////////////
//   Fast interpolation.
//
// Recursive algorithm, see algo 10.11 of Modern Computer Algebra
//
///////////////////////////////////////////////////

// TODO: save memory: free elements of v as soon as possible.

void fastinterpolateRec(ZZ_pX& P, const vec_ZZ_p& xi, const vec_ZZ_p& yi,
    const vec_ZZ_pX& SubPrd) {
  long n = xi.length();
  if (yi.length() != n) Error("interpolate: vector length mismatch");

  if (n == 0) {
    clear(P);
    return;
  }
  
  vec_ZZ_p si;
  si.SetLength(n);
  ZZ_pX pol = diff(SubPrd[1]);
  multievalRec(si, pol, SubPrd);
  for (int i = 0; i < n; ++i) 
    si[i] = yi[i]/si[i];

  vec_ZZ_pX v;
  long m = SubPrd.length() >> 1;
  v.SetLength(m);
  for (int i = 0; i < n; ++i)
    v[i] = si[i];
  m >>= 1;
  
  while (m > 0) {
    for (int k = 0; k < m; ++k)
      v[k] = v[2*k]*SubPrd[2*m + 2*k + 1] + v[2*k + 1]*SubPrd[2*m + 2*k];
    m >>= 1;
  }
  
  P = v[0];
}


void fastinterpolate(ZZ_pX& P, const vec_ZZ_p& xi, const vec_ZZ_p& yi) {
  vec_ZZ_pX SubPrd = buildSubProducts(xi);
  fastinterpolateRec(P, xi, yi, SubPrd);
}


///////////////////////////////////////////////////
//   Evaluate a polynomial at   x + a   (a is a scalar)
//
// Algorithm with Hadamard products/divides by Exp.
// see Bini-Pan: Polynomial and matrix computations
//     page 15, problem 2.6 (VAR-SHIFT)
//
///////////////////////////////////////////////////


ZZ_p& invCoeffOfExpX(const long n) {
  static vec_ZZ_p* SaveC = NULL;
  /* First time: allocate the table of remember */
  if (SaveC == NULL) {
    SaveC = new vec_ZZ_p;
    append(*SaveC, to_ZZ_p(to_ZZ(1)));
    append(*SaveC, to_ZZ_p(to_ZZ(1)));
  }
    
  if (SaveC->length() <= n) {
    long l;
    ZZ_p c;
    l = SaveC->length();
    while (l <= n) {
      c = (*SaveC)[l-1] * l;
      append(*SaveC, c);
      ++l;
    }
  }
  return ((*SaveC)[n]);
}
      
ZZ_p& coeffOfExpX(const long n) {
  static vec_ZZ_p* SaveC = NULL;
  /* First time: allocate the table of remember */
  if (SaveC == NULL) {
    SaveC = new vec_ZZ_p;
    append(*SaveC, to_ZZ_p(to_ZZ(1)));
    append(*SaveC, to_ZZ_p(to_ZZ(1)));
  }
    
  if (SaveC->length() <= n) {
    long l;
    ZZ_p c;
    l = SaveC->length();
    while (l <= n) {
      c = (*SaveC)[l-1] / l;
      append(*SaveC, c);
      ++l;
    }
  }
  return ((*SaveC)[n]);
}
      

void evalAtXPlusA(ZZ_pX& Q, const ZZ_pX& P, const ZZ_p& a) {
  // Reverse of Hadamard division of P with Exp(x)
  ZZ_pX R;
  ZZ_p coe;
  R = P;
  for (int i = 2; i <= deg(R); ++i) {
    coe = invCoeffOfExpX(i);
    SetCoeff(R, i, coeff(R, i)*coe);
  }
  reverse(R, R, deg(P));
 
  // Multiplication with Exp(a*x)
  ZZ_pX Eax;
  ZZ_p pow_a;
  Eax = 1;
  pow_a = a;
  for (int i = 1; i <= deg(P); ++i) {
    coe = coeffOfExpX(i)*pow_a;
    SetCoeff(Eax, i, coe);
    pow_a *= a;
  }
  MulTrunc(R, R, Eax, deg(P)+1);
  
  // Hadamard product with Exp(x)
  reverse(R, R, deg(P));
  coe = 2;
  for (int i = 2; i <= deg(P); ++i) {
    coe = coeffOfExpX(i);
    SetCoeff(R, i, coeff(R, i)*coe);
  }
   
  Q = R;
}




#if 0

int main(int argc, char **argv) {
  ZZ p = to_ZZ("10000000000000000000000013");
  ZZ_p::init(p);

  ZZ_pX f, g;
  ZZ_p a;

  random(f, 5);
  random(a);
  evalAtXPlusA(g, f, a);
  cout << "f = " << f << endl;
  cout << "a = " << a << endl;
  cout << "g = " << g << endl;
  
}

#endif
#if 0
  
  ZZ_pX f, g;
  vec_ZZ_p xi, yi, zi;

  int n = 1000;

  for (int i = 0; i < n; ++i) {
    append(xi, random_ZZ_p());
    append(yi, random_ZZ_p());
  }

  double tps = GetTime();
  fastinterpolate(f, xi, yi);
  cout << "Fast interp done in " << GetTime()-tps << " sec\n";
  tps = GetTime();
  interpolate(g, xi, yi);
  cout << "Naive interp done in " << GetTime()-tps << " sec\n";
    
  assert(f == g);

  return 0;
}

#endif

