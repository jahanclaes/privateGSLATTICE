#ifndef MATRIXOPS_H
#define MATRIXOPS_H
#include <complex>
#include <iostream>
#include "Blitz.h"

#define F77_DGETRF F77_FUNC(dgetrf,DGETRF)

#define F77_ZGETRF F77_FUNC(zgetrf,ZGETRF)


// #define F77_DGEMM  F77_FUNC(dgemm,DGEMM)

extern "C" int // void 
zgetrf_(int *m, int *n, std::complex<double> A[], int *lda, int ipiv[], int *info);

extern "C" int // void 
dgetrf_(int *m, int *n, double A[], int *lda, int ipiv[], int *info);


extern "C" void
zgetri_ (int *N, std::complex<double> A[], int *lda, int ipiv[],
            std::complex<double> work[], int *lwork, int *info);


/* extern "C" void */
/* dgemm_ (char *transA, char *transB, int *m, int *n, int *k, */
/*            double *alpha, const double *A, int *lda, const double *B, int *ldb, */
/*            double *beta,  double *C, int *ldc); */

extern "C"{
  void dgemm_(const char&, const char&,
             const int&, const int&, const int&,
             const double&, const double*, const int&, const double*, 
	     const int& ,const double&, double*, const int&);

}

extern "C"{
  void zgemm_(const char&, const char&,
             const int&, const int&, const int&,
             const std::complex<double>&, const std::complex<double>*, const int&, const std::complex<double>*, 
	      const int& ,const std::complex<double>&, std::complex<double>*, const int&);

}


struct MatrixOps {
// Adapted from Numerical Recipes in C
static void LUdecomp (blitz::Array<double,2> &A, blitz::Array<int,1> &perm,
               double &sign)
{
  int i, imax, j, k;
  int n = A.rows();
  double big, dum, sum, temp;
  blitz::Array<double,1> vv(n);
  sign = 1.0;
  
  perm.resize(n);
  
  for (i=0; i<n; i++) {
    big = 0.0;
    for (int j=0; j<n; j++)
      if (fabs(A(i,j)) > big) big = fabs(A(i,j));
    if (big == 0.0) {
      std::cerr << "Singularity in LUdecomp.\n";
      std::cerr << "Should abort!.\n";
      ///      abort();
    }
    vv(i) = 1.0/big;
  }
  for (j=0; j<n; j++) {
    for (i=0; i<j; i++) {
      sum=A(i,j);
      for(k=0; k<i; k++)
        sum -= A(i,k)*A(k,j);
      A(i,j) = sum;
    }
    big=0.0;
    for (i=j; i<n; i++) {
      sum = A(i,j);
      for (k=0; k<j; k++)
        sum-=A(i,k)*A(k,j);
      A(i,j) = sum;
      if ((dum=vv(i)*fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k=0; k<n; k++) {
        dum=A(imax,k);
        A(imax,k) = A(j,k);
        A(j,k) = dum;
      }
      sign = -sign;
      vv(imax) = vv(j);
    }
    perm(j)=imax;
    if (A(j,j) == 0.0)
      A(j,j) = 1.0e-200;
    if (j != (n-1)) {
      dum=1.0/A(j,j);
      for (i=j+1; i<n; i++)
        A(i,j) *= dum;
    }
  }
}
  
static void MatrixMultiply(blitz::Array<double,2> &A,blitz::Array<double,2> &B,
		    blitz::Array<double,2> &C)
{
  assert(A.extent(1)==B.extent(0));
  assert(A.extent(0)==C.extent(0));
  assert(B.extent(1)==C.extent(1));
  for (int i=0;i<C.extent(0);i++){
    for (int j=0;j<C.extent(1);j++){
      C(i,j)=0.0;
      for (int k=0;k<B.extent(0);k++)
	C(i,j)+=A(i,k)*B(k,j);
    }
  }
}

  
static void MatrixMultiply(blitz::Array<std::complex<double>,2> &A,blitz::Array<std::complex<double>,2> &B,
		    blitz::Array<std::complex<double>,2> &C)
{
  assert(A.extent(1)==B.extent(0));
  assert(A.extent(0)==C.extent(0));
  assert(B.extent(1)==C.extent(1));
  for (int i=0;i<C.extent(0);i++){
    for (int j=0;j<C.extent(1);j++){
      C(i,j)=0.0;
      for (int k=0;k<B.extent(0);k++)
	C(i,j)+=A(i,k)*B(k,j);
    }
  }
}


  inline static void product(const blitz::Array<double,2>& A,
			     const blitz::Array<double,2>& B, blitz::Array<double,2>& C) {
    assert(1==2);
    const char transa = 'N';
    const char transb = 'N';
    const double one=1.0;
    const double zero=0.0;
    dgemm_(transa, transb, B.cols(), A.rows(), B.rows(),
          one, B.data(), B.cols(), A.data(), A.cols(),
          zero, C.data(), C.cols());
  }


  inline static void product(const blitz::Array<std::complex<double>,2>& A,
			     const blitz::Array<std::complex<double>,2>& B, blitz::Array<std::complex<double>,2>& C) {
    const char transa = 'N';
    const char transb = 'N';
    const double one=1.0;
    const double zero=0.0;
    zgemm_(transa, transb, B.cols(), A.rows(), B.rows(),
          one, B.data(), B.cols(), A.data(), A.cols(),
          zero, C.data(), C.cols());
  }

  static void MatMult (const blitz::Array<double,2> &A, const blitz::Array<double,2> &B,
              blitz::Array<double,2> &C)
{
/*   int m = A.rows(); */
/*   int n = B.cols(); */
/*   int k = A.cols(); */
/*   assert (B.rows() == k); */
/*   // We use "transpose" operation because we have C ordering, which fortran */
/*   // thinks is transposed. */
/*   char transA = 'T'; */
/*   char transB = 'T'; */
/*   double alpha = 1.0; */
/*   double beta = 0.0; */
/*   Generalblitz::ArrayStorage<2> colMajor; */
/*   colMajor.ordering() = firstDim, secondDim; */




/*   dgemm_ (&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, */
/*              B.data(), &n, &beta, C.data(), &m); */
}

/*   static void product(blitz::Array<double,2>& A, */
/* 		      blitz::Array<double,2>& B, blitz::Array<double,2>& C) { */
/*     const char transa = 'N'; */
/*     const char transb = 'N'; */
/*     const double one=1.0; */
/*     const double zero=0.0; */
/*     dgemm(transa, transb, B.cols(), A.rows(), B.rows(), */
/*           one, B.data(), B.cols(), A.data(), A.cols(), */
/*           zero, C.data(), C.cols()); */
/*   } */

static void LUsolve (blitz::Array<double,2> &LU, blitz::Array<int,1> &perm,
              blitz::Array<double,1> &b)
{
  int i, ii=-1,ip,j;
  double sum;
  int n = LU.rows();

  for (i=0; i<n; i++) {
    ip = perm(i);
    sum = b(ip);
    b(ip) = b(i);
    if (ii>=0)
      for (j=ii; j<i; j++)
        sum -= LU(i,j)*b(j);
    else if (sum)
      ii = i;
    b(i) = sum;
  }
  for (i=n-1; i>=0; i--) {
    sum = b(i);
    for (j=i+1; j<n; j++)
      sum -= LU(i,j)*b(j);
    b(i) = sum/LU(i,i);
  }
}

static blitz::Array<double,2> Inverse(blitz::Array<double,2> &A)
{
  blitz::Array<double,2> LU(A.rows(), A.cols()), Ainv(A.rows(), A.cols());
  blitz::Array<double,1> col(A.rows());
  blitz::Array<int,1> perm;
  double sign;

  LU = A;
  LUdecomp (LU, perm, sign);
  for (int j=0; j<A.rows(); j++) {
    for (int i=0; i<A.rows(); i++)
      col(i) = 0.0;
    col(j) = 1.0;
    LUsolve (LU, perm, col);
    for (int i=0; i<A.rows(); i++)
      Ainv(i,j) = col(i);
  }
  return (Ainv);
}



/// This function returns the determinant of A and replaces A with its                                
/// cofactors.                                                                                        
static std::complex<double>
ComplexDetCofactors (blitz::Array<std::complex<double>,2> &A,
                     blitz::Array<std::complex<double>,1> &work)
{
  const int maxN = 2000;
  int ipiv[maxN];
  int N = A.rows();
  int M = A.cols();
  assert (N == M);
  assert (N <= maxN);
  // First, transpose A for fortran ordering                                                          
  //   for (int i=0; i<N; i++)                                                                          
  //     for (int j=0; j<i; j++) {                                                                      
  //       double tmp = A(i,j);                                                                         
  //       A(i,j) = A(j,i);                                                                             
  //       A(j,i) = tmp;                                                                                
  //     }                                                                                              
  Transpose(A);
  int info;
  // Do LU factorization                                                                              
  zgetrf_ (&N, &M, A.data(), &N, ipiv, &info);
  std::complex<double> det = 1.0;
  int numPerm = 0;
  for (int i=0; i<N; i++) {
    det *= A(i,i);
    numPerm += (ipiv[i] != (i+1));
  }
  if (numPerm & 1)
    det = -det;

  int lwork = work.size();
  // Now, do inverse                                                                                  
  zgetri_ (&N, A.data(), &N, ipiv, work.data(), &lwork, &info);

  // Now, we have the transpose of Ainv.  Now, just multiply by det:                                  
  //  A = det * A;
  // And we're done!                                                                                  
  return det;
}


static int
ComplexDetCofactorsWorkSize(int N)
{
  std::complex<double> work;
  std::complex<double> dummy;
  int info;
  int ipiv;
  int lwork = -1;

  zgetri_(&N, &dummy, &N, &ipiv, &work, &lwork, &info);

  return ((int)ceil(work.real()));
}


static blitz::Array<std::complex<double>, 2 > Inverse (blitz::Array<std::complex<double >,2 > &A)
{
  //  std::cerr<<"DOING THIS INVERSE"<<endl;
  int size =ComplexDetCofactorsWorkSize(A.extent(0));
  blitz::Array<std::complex<double>,2> myInverse(A.extent(0),A.extent(1));
  blitz::Array<std::complex<double>,1> work(size);
  myInverse(Range::all(),Range::all())=A(Range::all(),Range::all());
  ComplexDetCofactors(myInverse,work);
  Transpose(myInverse);
  return myInverse;
}

static inline void Transpose (blitz::Array<std::complex<double>,2> &A)
{
  int m = A.rows();
  int n = A.cols();
  if (m != n)
    OutOfPlaceTranspose (A);
  else {
    for (int i=0; i<m; i++)
      for (int j=i+1; j<m; j++)
        swap (A(i,j), A(j,i));
  }
}

static inline void OutOfPlaceTranspose (blitz::Array<std::complex<double>,2> &A)
{
  int m = A.rows();
  int n = A.cols();
  blitz::Array<std::complex<double>,2> Atrans (n,m);
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      Atrans(j,i) = A(i,j);
  A.resize(n,m);
  A = Atrans;
}


static double Determinant (const blitz::Array<double,2> &A)
{

  double logDet=0.0;
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take a determinant of a non-square
		    // matrix
  if (A.rows() == 1)
    return (A(0,0));
  if (A.rows() == 2) 
    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
  else {
    blitz::Array<double,2> LU(m,m);
    blitz::Array<int,1> ipiv(m);
    int info;
    LU = A;
    // Do LU factorization
    //    F77_DGETRF (&m, &n, LU.data(), &m, ipiv.data(), &info);
    dgetrf_(&m, &n, LU.data(), &m, ipiv.data(), &info);
    double det = 1.0;
    int numPerm = 0;
    for (int i=0; i<m; i++) {
      det *= LU(i,i);
      logDet+=log(abs(LU(i,i)));
      numPerm += (ipiv(i) != (i+1));
    }
    if (numPerm & 1)
      det *= -1.0;
    return det;
  }
}








static std::complex<double> Determinant (const blitz::Array<std::complex<double> ,2> &A)
{
  std::complex<double> logDet=0.0;
  int m = A.rows();
  int n = A.cols();
  assert (m == n);  // Cannot take a determinant of a non-square
		    // matrix
  if (A.rows() == 1)
    return (A(0,0));
  if (A.rows() == 2) 
    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
  else {
    blitz::Array<std::complex<double> ,2> LU(m,m);
    blitz::Array<int,1> ipiv(m);
    int info;
    LU = A;
    // Do LU factorization
    //    F77_DGETRF (&m, &n, LU.data(), &m, ipiv.data(), &info);
    zgetrf_(&m, &n, LU.data(), &m, ipiv.data(), &info);
    std::complex<double> det = 1.0;
    int numPerm = 0;
    for (int i=0; i<m; i++) {
      logDet+=log(LU(i,i));
      det *= LU(i,i);
      numPerm += (ipiv(i) != (i+1));
    }
    if (numPerm & 1)
      det *= -1.0;
    return det;
  }
}



static int CalcParity(blitz::Array<int,1> &permutation)
{
  blitz::Array<int,1> newp(permutation);
  int countSwaps=0;
  int i=0;
  while (i<newp.size()){
    while (permutation(i)!=i){
      int newSpot=permutation(i);
      swap(permutation(newSpot),permutation(i));
      countSwaps++;
    }
    i++;
  }
  std::cerr<<"Permutation swaps is "<<i<<endl;
  return i %2;

}



static double mag(std::complex<double> a)
{
  return sqrt(a.real()*a.real()+a.imag()*a.imag());
}

};



#endif


